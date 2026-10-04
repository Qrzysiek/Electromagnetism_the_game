//! Beams: many particles in one ODE system, interacting with each other (PHYSICS.md §3.3).
//!
//! All particles share one state vector and one step size. The interaction is exact
//! classical electrodynamics between point charges:
//!
//! - For `c = ∞` it is the pairwise Coulomb force, which is then the whole interaction:
//!   the magnetic and retarded parts, and every interaction of magnetic moments, scale as
//!   `1/c²` and vanish.
//! - For finite `c`, by default, the quasi-static interaction: each particle feels the
//!   exact Liénard–Wiechert fields of the others, with their recent past continued from
//!   their present state along the motion they would have in the uniform fields they
//!   feel now (`FieldMotion`, `Continuation::Fields`; the fields from a first pass with
//!   the fields of uniform motion, `heaviside_fields`), blending into the tapered constant
//!   acceleration (`accelerated_fields`) where that motion swings far, and with the
//!   tapered acceleration alone for a particle with a magnetic moment or out of the
//!   plane. That is exact in the velocities (so it contains the magnetic attraction that
//!   reduces the space charge of a relativistic beam by 1/γ²) and to first order in the
//!   accelerations (their near and radiation fields); it leaves out the jerk terms the
//!   continuation misses, estimated per particle (`BeamRun::neglected_retardation`). It
//!   needs no record and costs a few times the Coulomb force.
//! - For finite `c` with `BeamScenario::retarded`: each particle feels the Liénard–Wiechert
//!   fields (`lienard.rs`) of the others at their retarded times, taken from the recorded
//!   motion (the dense output of every accepted step). Retarded times inside the current
//!   step (pairs closer than `c h`) come from the previous step's polynomial,
//!   extrapolated (the standard treatment of delays shorter than a step), so the step is
//!   limited by that reach, not by the light time between the closest pair. Before launch
//!   the particles move with their launch acceleration. A removed particle acts for as
//!   long as its field from
//!   before the removal is on its way: its charge is taken to be drained at the moment of
//!   removal, as it vanishes for `c = ∞`. The moments' interaction (O(1/c²)) is not
//!   modelled. Each particle's own radiation reaction is optional (Landau–Lifshitz, in
//!   the total field: external plus the others' retarded fields); without it the
//!   radiated energy is estimated per particle (Liénard formula) so that neglecting it
//!   can be checked.
//!
//! Particles with a radius collide as rigid spheres (PHYSICS.md §3.3): when two touch,
//! the step ends there and each gets an impulse along the line of centres that conserves
//! momentum and kinetic energy exactly (relativistically: the positive root of the energy
//! balance). Point particles (radius 0) never touch. The radiation of the impact itself
//! depends on the spheres' structure and is not modelled; for finite `c` an estimate,
//! `q² |Δv|² / (3 a c²)` for a velocity change over the light crossing time of the
//! sphere (radius a), is added to the particle's radiated energy, so it is checked like
//! the rest of the neglected radiation.
//!
//! Each particle has the events and gates of a single flight (obstacles, bounds, detector,
//! gates in order). The earliest event of any particle ends the step there; that
//! particle's outcome is recorded, and the integration restarts at that time without it
//! as a source. It is still integrated for a while as a *ghost* (pushed by the others,
//! pushing none) until its penetration depth into the boundary it crossed is known,
//! exactly as the single-particle runner follows the continued trajectory
//! (`trajectory.rs`): margins are then properties of the trajectories, not of where the
//! steps ended, and the existing verification (`verify::classify`) applies to every
//! particle.

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::{Mutex, RwLock};

use glam::DVec3;

use crate::dynamics::{Kinematics, Particle, landau_lifshitz};
use crate::events::first_crossing;
use crate::field::FieldSolver;
use crate::field_motion::FieldMotion;
use crate::geometry::{Aabb, Region, Shape};
use crate::integrator::OdeSystem;
use crate::integrator::dop853::{Dense, Dop853, Settings, Stats};
use crate::lienard::fields_from;
use crate::spectrum::Emission;
use crate::trajectory::{
    Acceptance, Gate, GateTracker, MARGIN_SAFE, Margins, Outcome, RunSettings, Sample, Trajectory,
    minimize_on, minimum_at_end,
};

/// One particle of a beam with its launch state and its own detector (particles of
/// different species may be meant for different detectors).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BeamParticle {
    pub particle: Particle,
    pub x0: DVec3,
    pub p0: DVec3,
    pub detector: Option<Region>,
    pub acceptance: Option<Acceptance>,
}

/// What happens to a particle at a boundary (PHYSICS.md §3.3). Charge is conserved: an
/// absorbed charge either stays where it was absorbed or is carried away.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fate {
    /// Absorbed and carried away (grounded metal, e.g. a Faraday cup, which screens the
    /// charge's field from outside): its field disappears, at once for `c = ∞`, where the
    /// light cone of the absorption passes otherwise.
    Drain,
    /// Absorbed by a detector that is the mouth of a deep grounded Faraday cup (the
    /// default of levels, PHYSICS.md §3.3): the particle flies on into the cup at its
    /// entry velocity while the cup screens its charge, seen from outside, as
    /// `q e^{−k depth}` (`k = π/w`, w the mouth's width across the entry direction: the
    /// pipe's lowest evanescent mode), at the retarded time; its field fades instead of
    /// vanishing at once. Other boundaries treat it as `Drain`.
    Cup,
    /// Absorbed where it hits, its charge staying there at rest (an insulating body): it
    /// keeps acting as a charge at rest.
    Stop,
    /// Not a physical boundary (the edge of the arena, which is only the view): the
    /// particle flies on and keeps acting; it only counts as lost.
    Pass,
}

/// The beam's energy budget at one time (recorded with `RunSettings::record`), in the
/// units of the fields. Charges present are the particles still flying (also beyond the
/// arena) and the charges of absorbed particles that stayed where they stopped.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergySample {
    pub t: f64,
    /// Kinetic energy of the moving particles.
    pub kinetic: f64,
    /// Potential energy of the charges present in the level's fields (external potential,
    /// induced charges of metal, magnetic moments).
    pub potential: f64,
    /// Coulomb interaction energy of the charges present with each other (for finite `c`
    /// only the electric part: the energy of the magnetic and radiation fields between
    /// them is not counted).
    pub interaction: f64,
    /// Energy given to the bodies and the detector by absorbed particles: the drop of the
    /// energy above at each absorption (kinetic energy; for a drained particle also its
    /// potential and interaction energy).
    pub absorbed: f64,
    /// Energy radiated so far (Liénard formula, all particles).
    pub radiated: f64,
}

/// Fates at the three kinds of boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fates {
    pub detector: Fate,
    pub obstacles: Fate,
    pub bounds: Fate,
}

impl Default for Fates {
    fn default() -> Self {
        Self {
            detector: Fate::Drain,
            obstacles: Fate::Stop,
            bounds: Fate::Pass,
        }
    }
}

/// A beam flight: the shared scene and the particles, all launched at t = 0.
#[derive(Clone, Debug)]
pub struct BeamScenario<F> {
    pub field: F,
    pub obstacles: Vec<Shape>,
    pub particles: Vec<BeamParticle>,
    pub c: f64,
    pub bounds: Option<Aabb>,
    pub t_max: f64,
    /// Interaction between the particles: Coulomb for `c = ∞`; for finite `c` quasi-static
    /// (Liénard–Wiechert fields of the present state continued back with constant
    /// acceleration), or with `retarded` the exact Liénard–Wiechert fields.
    pub interact: bool,
    /// Exact retarded interaction at finite `c` (costly: steps shorter than the light
    /// time between the closest pair); otherwise quasi-static.
    pub retarded: bool,
    /// Gates every particle must pass, in order, before its detector counts (as for single
    /// flights, PHYSICS.md §6.2).
    pub gates: Vec<Gate>,
    /// Each particle's radiation reaction (Landau–Lifshitz, PHYSICS.md §3.1) in the total
    /// field; no effect for `c = ∞`.
    pub radiation_reaction: bool,
    /// What happens to a particle at each kind of boundary.
    pub fates: Fates,
}

/// How the quasi-static interaction continues a source's past from its present state.
#[derive(Clone, Copy, Debug)]
enum Continuation {
    /// Along the motion it would have in the fields it feels now, taken as uniform: exact
    /// and relativistic (`FieldMotion`; the rotation of its acceleration in a magnetic
    /// field included, and never faster than light).
    Fields(FieldMotion),
    /// With its present acceleration, tapered (`tapered`): a particle with a magnetic
    /// moment (its gradient force is not a Lorentz force), or out of the plane.
    Accelerated(DVec3),
}

/// Where the motion in the fields stops being trusted: from an excursion
/// `|U(s) − U₀|/c` of `SWING[0]` at the retarded time it blends smoothly into the tapered
/// constant acceleration, which alone is used from `SWING[1]` on. A charge in an enormous
/// field (plunging into a charge) would otherwise be continued back through turns and
/// reversals whose fields sweep over the others (the steps collapsed); in the shipped
/// levels the excursions stay below 0.35 (a charge circling at 0.4c whose partner's light
/// time is 0.8 rad of its turn).
const SWING: [f64; 2] = [0.5, 1.0];

/// The weight of the motion in the fields at the excursion `e` (`SWING`): 1, then a
/// quintic smoothstep to 0 (continuous to the second derivative).
fn swing_weight(e: f64) -> f64 {
    let u = ((e - SWING[0]) / (SWING[1] - SWING[0])).clamp(0.0, 1.0);
    1.0 - u * u * u * (10.0 - 15.0 * u + 6.0 * u * u)
}

impl Continuation {
    /// The fields `(E, B)` at `x`, time `dt` after the source's present, of the source `q`
    /// now at `r` with velocity `v`, and the weight of the motion in the fields in them (0
    /// for the constant acceleration).
    fn fields(&self, q: f64, c: f64, x: DVec3, dt: f64, r: DVec3, v: DVec3) -> (DVec3, DVec3, f64) {
        match self {
            Continuation::Fields(m) => {
                // No retarded point (beyond the horizon of a charge accelerated forever):
                // the constant acceleration alone.
                let ret = m.retarded(x, dt);
                let w = ret.map_or(0.0, |r| swing_weight(r.excursion));
                let (mut e, mut b) = match ret {
                    Some(r) if w > 0.0 => {
                        let (e, b) = m.fields(q, x, &r);
                        (e * w, b * w)
                    }
                    _ => (DVec3::ZERO, DVec3::ZERO),
                };
                if w < 1.0 {
                    let (ea, ba) = accelerated_fields(q, c, x, dt, r, v, m.present().1);
                    e += ea * (1.0 - w);
                    b += ba * (1.0 - w);
                }
                (e, b, w)
            }
            Continuation::Accelerated(a) => {
                let (e, b) = accelerated_fields(q, c, x, dt, r, v, *a);
                (e, b, 0.0)
            }
        }
    }

    /// The jerk of the continued past at the present (0 for constant acceleration).
    fn jerk(&self) -> DVec3 {
        match self {
            Continuation::Fields(m) => m.jerk(),
            Continuation::Accelerated(_) => DVec3::ZERO,
        }
    }
}

/// A particle absorbed by a screening cup (`Fate::Cup`): it flies on into the cup at its
/// entry velocity, and seen from outside its charge fades as `e^{−k v_n (t − t_off)}`
/// (`k = π/w`, `v_n` its speed into the cup), at the retarded time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fade {
    /// When and where it entered the cup.
    pub t_off: f64,
    pub x_off: DVec3,
    /// Its velocity from then on.
    pub v: DVec3,
    /// `k v_n`: the charge's decay rate.
    pub rate: f64,
}

impl Fade {
    /// Its position inside the cup at `t` (after `t_off`).
    pub fn position(&self, t: f64) -> DVec3 {
        self.x_off + self.v * (t - self.t_off)
    }

    /// The quasi-static interaction's field `(E, B)` of the fading charge `q` at `x`, `t`:
    /// the field of its uniform motion (`heaviside_fields`), with the charge at the
    /// retarded time of that motion (closed form).
    pub fn quasi_static_fields(&self, q: f64, c: f64, x: DVec3, t: f64) -> (DVec3, DVec3) {
        let r_now = self.position(t);
        let d = x - r_now;
        let (dv, v2) = (d.dot(self.v), self.v.length_squared());
        let tau = (dv + (dv * dv + (c * c - v2) * d.length_squared()).sqrt()) / (c * c - v2);
        let q_eff = q * self.factor(t - tau);
        if q_eff == 0.0 {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        heaviside_fields(q_eff, c, x, r_now, self.v)
    }

    /// The factor of its charge seen by the field emitted at the retarded time `t_r`.
    pub fn factor(&self, t_r: f64) -> f64 {
        if t_r <= self.t_off {
            1.0
        } else {
            libm::exp(-self.rate * (t_r - self.t_off))
        }
    }

    /// The fade of a particle entering `region` at `x` with velocity `v` at `t_off`.
    fn at_entry(region: &Region, x: DVec3, v: DVec3, t_off: f64) -> Fade {
        let (width, normal) = match region {
            Region::Box(b) => {
                // The entry face: the one nearest to x; the width: the box's in-plane extent
                // along the other in-plane axis.
                let faces = [
                    ((x.x - b.min.x).abs(), DVec3::X, b.max.y - b.min.y),
                    ((b.max.x - x.x).abs(), -DVec3::X, b.max.y - b.min.y),
                    ((x.y - b.min.y).abs(), DVec3::Y, b.max.x - b.min.x),
                    ((b.max.y - x.y).abs(), -DVec3::Y, b.max.x - b.min.x),
                ];
                let f = faces
                    .iter()
                    .min_by(|a, b| a.0.total_cmp(&b.0))
                    .expect("four faces");
                (f.2, f.1)
            }
            Region::Sphere(sph) => (2.0 * sph.radius, -(x - sph.center).normalize_or_zero()),
        };
        let vn = v.dot(normal).max(1e-3 * v.length());
        Fade {
            t_off,
            x_off: x,
            v,
            rate: std::f64::consts::PI / width.max(1e-9) * vn,
        }
    }
}

/// Result of a beam flight: one trajectory per particle (in launch order; their `stats`
/// are those of the whole system), and diagnostics of the whole system.
#[derive(Clone, Debug)]
pub struct BeamRun {
    pub trajectories: Vec<Trajectory>,
    pub stats: Stats,
    /// Largest relative change of the total energy of the particles still flying
    /// (kinetic, external potential, moments, pair interaction), between removals. NaN for
    /// interacting particles at finite `c`, where the particles alone do not conserve
    /// energy (the field carries energy, and radiates it).
    pub energy_max_rel_error: f64,
    /// Integrator restarts (one per removal and per finished ghost).
    pub restarts: usize,
    /// The energy budget at launch and at every accepted step (empty unless recorded).
    pub energy: Vec<EnergySample>,
    /// Quasi-static interaction: per particle, the estimated relative error of the
    /// interaction's impulse: the time integral of the neglected fields of the others (from
    /// their jerk, `Σ_j |q_j| γ_j² |ȧ_j| / (c³ κ)` with the Doppler factor κ = 1 − n·β
    /// for the longer light delay ahead of a source, plus `|q_j| γ_j² |a_j| / (c² R)` for
    /// pairs so far apart that the source's continued past has left its constant
    /// acceleration, `|a_j| R / (c κ) > 0.7 Δv`, see `tapered`), over that of the
    /// fields kept, `Σ_j |E_j|`; plus, at each drained particle's removal, the impulse
    /// `|q| / (R c)` of its field lingering for the light time R/c (the quasi-static
    /// interaction drops it at once)
    /// (sampled at the step ends; 0 without quasi-static interaction). The trajectory
    /// error follows the impulse error: a short plunge of a neighbour matters little.
    pub neglected_retardation: Vec<f64>,
    /// Per particle: its fade in a screening cup, if it entered one (`Fate::Cup`).
    pub fades: Vec<Option<Fade>>,
}

/// Event functions of one particle, in priority order.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Event {
    Obstacle(usize),
    Bounds,
    Detector,
}

impl Event {
    fn fate(self, fates: &Fates) -> Fate {
        match self {
            Event::Obstacle(_) => fates.obstacles,
            Event::Bounds => fates.bounds,
            Event::Detector => fates.detector,
        }
    }

    fn outcome(self) -> Outcome {
        match self {
            Event::Obstacle(i) => Outcome::Collided(i),
            Event::Bounds => Outcome::LeftBounds,
            Event::Detector => Outcome::Arrived,
        }
    }
}

/// The right-hand side for the particles in `members` (indices into the scenario), of
/// which those with `source` set act on the others. Scaled state per member:
/// `(x, p / p_ref)`.
struct BeamOde<'a, F> {
    scn: &'a BeamScenario<F>,
    members: Vec<usize>,
    source: Vec<bool>,
    kin: Vec<Kinematics>,
    p_ref: f64,
    /// Retarded interaction (finite `c`): the recorded motion of every particle.
    past: Option<&'a Past>,
    /// Charges of absorbed particles at rest where they stopped (`Fate::Stop`), acting as
    /// static sources (in the Coulomb and quasi-static interactions; the retarded one
    /// has them in its record).
    stopped: Vec<(usize, f64, DVec3)>,
    /// Drained charges fading in their cups (`Fade`, experimental).
    fading: Vec<(usize, f64, Fade)>,
}

/// Recorded motion of the beam for the retarded interaction: the dense output of every
/// accepted step (valid up to `t_stop`, where an event may have cut it short), with the
/// particles it contains, and for every particle the time its charge was removed.
struct Past {
    /// Recorded steps, oldest first; steps that no retarded time can reach any more are
    /// dropped (`prune`).
    segments: RwLock<Vec<Segment>>,
    /// Removal time and position of each particle (infinite time while it flies).
    t_off: RwLock<Vec<f64>>,
    x_off: RwLock<Vec<DVec3>>,
    /// Whether a removed particle stayed where it stopped (at rest from then on) rather
    /// than being drained.
    stays: RwLock<Vec<bool>>,
    /// Drained particles fading in their cups (`Fade`, experimental).
    fades: RwLock<Vec<Option<Fade>>>,
    kin: Vec<Kinematics>,
    p_ref: f64,
    /// Acceleration of each particle at launch (external fields and the others'
    /// quasi-static fields): the motion before launch continues it.
    a0: Vec<DVec3>,
    /// The last retarded time found for each (receiver, source) pair: row `i` for
    /// receiver `i` (NaN: none yet). Retarded times move smoothly, so it starts the next
    /// search. One row per receiver: the receivers' forces are computed in parallel, each
    /// touching only its own row (deterministic).
    guess: Vec<Mutex<Vec<f64>>>,
    /// Diagnostics (`EM_BEAM_LOG`): retarded-time solves, position evaluations in them,
    /// Newton iterations, and solves that fell back to bracketing by doubling.
    solves: AtomicU64,
    /// Nanoseconds spent in the right-hand sides (the parallel part).
    rhs_nanos: AtomicU64,
    evals: AtomicU64,
    newton: AtomicU64,
    fallbacks: AtomicU64,
}

/// Seconds from a nanosecond counter (diagnostics).
#[allow(clippy::cast_precision_loss)]
fn rhs_seconds(nanos: &AtomicU64) -> f64 {
    nanos.load(Relaxed) as f64 * 1e-9
}

/// A read-only view of the recorded motion (taken once per field evaluation, shared by
/// the threads computing the receivers' forces).
struct View<'a> {
    past: &'a Past,
    segments: &'a [Segment],
    t_off: &'a [f64],
    x_off: &'a [DVec3],
    stays: &'a [bool],
    fades: &'a [Option<Fade>],
}

struct Segment {
    dense: Dense,
    t_start: f64,
    t_stop: f64,
    /// Particle indices (sorted): member `k` is the state's block `k`.
    members: Vec<usize>,
}

impl Past {
    /// A read-only view of the record.
    fn view(&self) -> ViewGuards<'_> {
        ViewGuards {
            past: self,
            segments: self.segments.read().expect("not poisoned"),
            t_off: self.t_off.read().expect("not poisoned"),
            x_off: self.x_off.read().expect("not poisoned"),
            stays: self.stays.read().expect("not poisoned"),
            fades: self.fades.read().expect("not poisoned"),
        }
    }

    /// Drops the recorded steps that no retarded time can reach any more. A retarded time
    /// never decreases along a world line, so the earliest one needed from now on is at
    /// least `t − D / (c (1 − β))`, with `D` the extent of all current positions and
    /// removal points and `β` the largest speed so far: with a factor 2 to spare, older
    /// steps go.
    fn prune<F: FieldSolver>(
        &self,
        members: &[usize],
        y: &[f64],
        ode: &BeamOde<'_, F>,
        c: f64,
        beta_max: &mut f64,
        t: f64,
    ) {
        let (mut lo, mut hi) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
        for (k, _) in members.iter().enumerate() {
            let x = BeamOde::<F>::x(y, k);
            (lo, hi) = (lo.min(x), hi.max(x));
            *beta_max = beta_max.max(ode.kin[k].velocity(ode.p(y, k)).length() / c);
        }
        for (j, &x) in self.x_off.read().expect("not poisoned").iter().enumerate() {
            if !members.contains(&j) {
                (lo, hi) = (lo.min(x), hi.max(x));
            }
        }
        let extent = (hi - lo).length();
        let cutoff = t - 2.0 * extent / (c * (1.0 - beta_max.min(0.999_999))) - 1e-3;
        let mut segments = self.segments.write().expect("not poisoned");
        let old = segments.partition_point(|s| s.t_stop < cutoff);
        // Drop in batches (and always keep the newest step).
        if old >= 256 && old < segments.len() {
            segments.drain(..old);
        }
    }
}

/// The read guards behind a `View`.
struct ViewGuards<'a> {
    past: &'a Past,
    segments: std::sync::RwLockReadGuard<'a, Vec<Segment>>,
    t_off: std::sync::RwLockReadGuard<'a, Vec<f64>>,
    x_off: std::sync::RwLockReadGuard<'a, Vec<DVec3>>,
    stays: std::sync::RwLockReadGuard<'a, Vec<bool>>,
    fades: std::sync::RwLockReadGuard<'a, Vec<Option<Fade>>>,
}

impl ViewGuards<'_> {
    fn view(&self) -> View<'_> {
        View {
            past: self.past,
            segments: &self.segments,
            t_off: &self.t_off,
            x_off: &self.x_off,
            stays: &self.stays,
            fades: &self.fades,
        }
    }
}

impl View<'_> {
    /// Before launch each particle moves with its launch acceleration (it was already
    /// flying in the fields), tapered smoothly so that its velocity changes by less than
    /// `taper_reach` (at most 0.1c, and never to c: `tapered`). A sudden start from
    /// uniform motion would send a kink in every acceleration field to every other
    /// particle, each forcing tiny steps when it arrives.
    fn before_launch<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3, DVec3) {
        let b = &scn.particles[j];
        let v0 = self.past.kin[j].velocity(b.p0);
        tapered(b.x0, v0, self.past.a0[j], t, scn.c)
    }

    /// Position, velocity and acceleration of particle `j` at time `t`: `before_launch`,
    /// then the recorded motion, extrapolated a little past its end (see `EXTRAPOLATION`).
    fn state<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3, DVec3) {
        // Flying on into its cup.
        if let Some(f) = self.fades[j]
            && t >= f.t_off
        {
            return (f.position(t), f.v, DVec3::ZERO);
        }
        // At rest where it stopped.
        if self.stays[j] && t >= self.t_off[j] {
            return (self.x_off[j], DVec3::ZERO, DVec3::ZERO);
        }
        let kin = &self.past.kin[j];
        // Last segment starting before t that contains j.
        let segments = self.segments;
        let upto = segments.partition_point(|s| s.t_start < t);
        let found = segments[..upto]
            .iter()
            .rev()
            .find_map(|s| s.members.binary_search(&j).ok().map(|k| (s, k)));
        let Some((seg, k)) = found else {
            // Before the record: the motion before launch, or, where old steps were
            // dropped, uniform motion continued back from the oldest step kept. (Only
            // reached while bracketing a retarded time: every retarded time needed lies
            // in the kept record, and the continuation keeps g monotonic.)
            if t > 0.0
                && let Some((s0, k0)) = segments
                    .iter()
                    .find_map(|s| s.members.binary_search(&j).ok().map(|k| (s, k)))
            {
                let d = &s0.dense;
                let comp = |i: usize| d.eval_component(6 * k0 + i, s0.t_start);
                let x = DVec3::new(comp(0), comp(1), comp(2));
                let v = kin.velocity(DVec3::new(comp(3), comp(4), comp(5)) * self.past.p_ref);
                return (x + v * (t - s0.t_start), v, DVec3::ZERO);
            }
            return self.before_launch(scn, j, t);
        };
        let d = &seg.dense;
        // Past the end of the record (only while bracketing a retarded time, or for the
        // tiny offsets of the radiation reaction's field derivative): the step's own
        // polynomial for up to one step length, uniform motion beyond.
        let reach = seg.t_stop + EXTRAPOLATION * (seg.t_stop - seg.t_start);
        let tt = t.min(reach);
        let comp = |i: usize| d.eval_component(6 * k + i, tt);
        let x = DVec3::new(comp(0), comp(1), comp(2));
        let p = DVec3::new(comp(3), comp(4), comp(5)) * self.past.p_ref;
        let v = kin.velocity(p);
        if t > reach {
            return (x + v * (t - reach), v, DVec3::ZERO);
        }
        let dp = DVec3::new(
            d.eval_derivative_component(6 * k + 3, tt),
            d.eval_derivative_component(6 * k + 4, tt),
            d.eval_derivative_component(6 * k + 5, tt),
        ) * self.past.p_ref;
        (x, v, kin.acceleration(p, dp))
    }

    /// Retarded fields `(E, B)` of particle `j` at `x` and time `t` (`guess_row`: the
    /// receiver's last retarded times); zero once the field of its removal has arrived.
    fn fields_of<F>(
        &self,
        scn: &BeamScenario<F>,
        guess_row: &mut [f64],
        j: usize,
        x: DVec3,
        t: f64,
    ) -> (DVec3, DVec3) {
        let c = scn.c;
        let qj = scn.particles[j].particle.charge;
        let g = |tr: f64| {
            self.past.evals.fetch_add(1, Relaxed);
            c * (t - tr) - (x - self.position(scn, j, tr)).length()
        };
        self.past.solves.fetch_add(1, Relaxed);
        // The retarded time is the root of g, which decreases strictly (|v| < c). The
        // field of a removed particle is gone once the light cone has passed its removal
        // (decided from the exact removal point).
        let fade = self.fades[j];
        let t_off = if fade.is_some() {
            f64::INFINITY
        } else {
            self.t_off[j]
        };
        let stays = self.stays[j];
        if !stays && t_off.is_finite() && c * (t - t_off) >= (x - self.x_off[j]).length() {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        let mut hi = if stays { t } else { t.min(t_off) };
        // While the source flies, g(t) = −|x − r(t)| < 0: the retarded time exists and
        // needs no check; after its removal it may not (the light cone has not reached x).
        let removed = !stays && t_off.is_finite();
        let mut g_hi = if removed { g(hi) } else { f64::NAN };
        if g_hi >= 0.0 {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        let guess = guess_row[j];
        let tol = 4.0 * f64::EPSILON * t.abs().max(1.0);
        // Newton's method from the last retarded time of this pair (retarded times move
        // smoothly: it is close, and Newton converges in a few iterations), kept below
        // `hi`; if it does not converge, the bracketed search below.
        if guess.is_finite() && guess < hi {
            let mut tr = guess;
            for _ in 0..8 {
                self.past.newton.fetch_add(1, Relaxed);
                let (r, vr) = self.position_velocity(scn, j, tr);
                let d = x - r;
                let dist = d.length();
                let gv = c * (t - tr) - dist;
                let slope = -c + d.dot(vr) / dist.max(1e-300);
                let mut next = tr - gv / slope;
                if next >= hi {
                    next = 0.5 * (tr + hi);
                }
                // Quadratic convergence: after a step below 1e-8 (relative) the error is
                // at rounding, so that step is the last one.
                if (next - tr).abs() <= tol.max(1e-8 * t.abs().max(1.0)) {
                    guess_row[j] = next;
                    let (r, vr, ar) = self.state(scn, j, next);
                    let q_eff = qj * fade.map_or(1.0, |fd| fd.factor(next));
                    let f = fields_from(q_eff, c, x, next, r, vr, ar);
                    return (f.e(), f.b);
                }
                tr = next;
            }
        }
        // A bracket [lo, hi] with g(lo) > 0 ≥ g(hi), by doubling steps back from `hi`.
        if g_hi.is_nan() {
            g_hi = g(hi);
        }
        self.past.fallbacks.fetch_add(1, Relaxed);
        // Never a zero step (g_hi < 0 here, but keep the doubling alive regardless).
        let mut step = (-g_hi / c).max(4.0 * f64::EPSILON * t.abs().max(1.0));
        let mut lo = hi - step;
        let mut g_lo = g(lo);
        while g_lo <= 0.0 {
            (hi, g_hi) = (lo, g_lo);
            step *= 2.0;
            lo = hi - step;
            g_lo = g(lo);
        }
        // Newton's method, kept inside the bracket (bisection when it would leave it).
        let mut tr = lo + (hi - lo) * g_lo / (g_lo - g_hi);
        for _ in 0..100 {
            self.past.newton.fetch_add(1, Relaxed);
            let (r, vr) = self.position_velocity(scn, j, tr);
            let d = x - r;
            let dist = d.length();
            let gv = c * (t - tr) - dist;
            if gv > 0.0 {
                lo = tr;
            } else {
                hi = tr;
            }
            let slope = -c + d.dot(vr) / dist;
            let mut next = tr - gv / slope;
            if !(next > lo && next < hi) {
                next = 0.5 * (lo + hi);
            }
            let done = (next - tr).abs() <= tol || hi - lo <= tol;
            tr = next;
            if done {
                break;
            }
        }
        guess_row[j] = tr;
        let (r, vr, ar) = self.state(scn, j, tr);
        let q_eff = qj * fade.map_or(1.0, |fd| fd.factor(tr));
        let f = fields_from(q_eff, c, x, tr, r, vr, ar);
        (f.e(), f.b)
    }

    /// Position of particle `j` at time `t` (as `state`, without the rest).
    fn position<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> DVec3 {
        self.position_velocity(scn, j, t).0
    }

    /// Position and velocity of particle `j` at time `t` (as `state`, without the
    /// acceleration, which needs the polynomial's derivative).
    fn position_velocity<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3) {
        if let Some(f) = self.fades[j]
            && t >= f.t_off
        {
            return (f.position(t), f.v);
        }
        if self.stays[j] && t >= self.t_off[j] {
            return (self.x_off[j], DVec3::ZERO);
        }
        let segments = self.segments;
        let upto = segments.partition_point(|s| s.t_start < t);
        let found = segments[..upto]
            .iter()
            .rev()
            .find_map(|s| s.members.binary_search(&j).ok().map(|k| (s, k)));
        let Some((seg, k)) = found else {
            let (x, v, _) = self.state(scn, j, t);
            return (x, v);
        };
        let d = &seg.dense;
        let reach = seg.t_stop + EXTRAPOLATION * (seg.t_stop - seg.t_start);
        let tt = t.min(reach);
        let comp = |i: usize| d.eval_component(6 * k + i, tt);
        let x = DVec3::new(comp(0), comp(1), comp(2));
        let v = self.past.kin[j].velocity(DVec3::new(comp(3), comp(4), comp(5)) * self.past.p_ref);
        if t > reach {
            return (x + v * (t - reach), v);
        }
        (x, v)
    }
}

impl<F: FieldSolver> BeamOde<'_, F> {
    fn x(y: &[f64], k: usize) -> DVec3 {
        DVec3::new(y[6 * k], y[6 * k + 1], y[6 * k + 2])
    }

    fn p(&self, y: &[f64], k: usize) -> DVec3 {
        DVec3::new(y[6 * k + 3], y[6 * k + 4], y[6 * k + 5]) * self.p_ref
    }

    /// Force on member `k` from the external fields and the source members.
    fn force(&self, y: &[f64], k: usize, t: f64) -> DVec3 {
        let acc = self.continuations(y, t);
        self.force_with(y, k, t, acc.as_deref(), true)
    }

    /// Quasi-static interaction: how every member's past is continued (first pass): along
    /// its motion in the fields it feels in the state `y`, the external fields and the
    /// others' fields of uniform motion.
    fn continuations(&self, y: &[f64], t: f64) -> Option<Vec<Continuation>> {
        self.quasi_static().then(|| {
            (0..self.members.len())
                .map(|j| {
                    let part = self.scn.particles[self.members[j]].particle;
                    let x = Self::x(y, j);
                    let p = self.p(y, j);
                    if part.charge != 0.0 && part.moment == 0.0 {
                        let f = self.scn.field.sample(x, t);
                        let e_self = self.scn.field.self_field(x, part.charge).0;
                        let (e, b) = self.quasi_static_fields(y, t, j, x, t, None);
                        if let Some(m) = FieldMotion::new(
                            part.charge,
                            part.mass,
                            self.scn.c,
                            t,
                            x,
                            self.kin[j].velocity(p),
                            f.e + e_self + e,
                            f.b + b,
                        ) {
                            return Continuation::Fields(m);
                        }
                    }
                    let f = self.force_with(y, j, t, None, false);
                    Continuation::Accelerated(self.kin[j].acceleration(p, f))
                })
                .collect()
        })
    }

    /// Force on member `k`; for the quasi-static interaction `acc` are the members'
    /// continued pasts (None: fields of uniform motion), `rr` adds radiation reaction.
    fn force_with(
        &self,
        y: &[f64],
        k: usize,
        t: f64,
        acc: Option<&[Continuation]>,
        rr: bool,
    ) -> DVec3 {
        self.force_in(y, k, t, acc, rr, None)
    }

    /// `force_with`, reading the recorded motion through `view` when given (taken once per
    /// right-hand side and shared by the threads, so that they take no locks).
    fn force_in(
        &self,
        y: &[f64],
        k: usize,
        t: f64,
        acc: Option<&[Continuation]>,
        rr: bool,
        view: Option<&View<'_>>,
    ) -> DVec3 {
        let part = self.scn.particles[self.members[k]].particle;
        let x = Self::x(y, k);
        let v = self.kin[k].velocity(self.p(y, k));
        let f = self.scn.field.sample(x, t);
        let e_self = self.scn.field.self_field(x, part.charge).0;
        let mut force = (f.e + e_self + v.cross(f.b)) * part.charge;
        if part.moment != 0.0 {
            force += self.scn.field.grad_bz(x, t) * part.moment;
        }
        if self.past.is_some() {
            if part.charge != 0.0 {
                let (e, b) = self.retarded_fields(k, x, t, view);
                force += (e + v.cross(b)) * part.charge;
            }
        } else if self.quasi_static() {
            if part.charge != 0.0 {
                let (e, b) = self.quasi_static_fields(y, t, k, x, t, acc);
                force += (e + v.cross(b)) * part.charge;
            }
        } else if self.scn.interact && part.charge != 0.0 {
            for &(j, qs, xs) in &self.stopped {
                // A particle's ghost does not feel its own stopped charge.
                if j != self.members[k] {
                    let d = x - xs;
                    let r2 = d.length_squared();
                    force += d * (part.charge * qs / (r2 * r2.sqrt()));
                }
            }
            for &(j, qf, fd) in &self.fading {
                if j != self.members[k] {
                    let d = x - fd.position(t);
                    let r2 = d.length_squared();
                    force += d * (part.charge * qf * fd.factor(t) / (r2 * r2.sqrt()));
                }
            }
            for (j, &src) in self.source.iter().enumerate() {
                if j == k || !src {
                    continue;
                }
                let qj = self.scn.particles[self.members[j]].particle.charge;
                if qj == 0.0 {
                    continue;
                }
                let d = x - Self::x(y, j);
                let r2 = d.length_squared();
                force += d * (part.charge * qj / (r2 * r2.sqrt()));
            }
        }
        if rr && self.reacts(k) {
            force += self.radiation_reaction_force(y, k, t, view);
        }
        force
    }

    /// Whether member `k` feels its radiation reaction.
    fn reacts(&self, k: usize) -> bool {
        self.scn.radiation_reaction
            && self.scn.c.is_finite()
            && self.scn.particles[self.members[k]].particle.charge != 0.0
    }

    /// Whether the interaction is quasi-static (finite `c`, not retarded).
    fn quasi_static(&self) -> bool {
        self.scn.interact && self.scn.c.is_finite() && self.past.is_none()
    }

    /// Quasi-static fields at `x`, time `t`, of the flying members other than `k`, in the
    /// state `y` of time `t_y`: from their pasts continued as `acc` says, or (None) as
    /// uniform motion.
    fn quasi_static_fields(
        &self,
        y: &[f64],
        t_y: f64,
        k: usize,
        x: DVec3,
        t: f64,
        acc: Option<&[Continuation]>,
    ) -> (DVec3, DVec3) {
        let (mut e, mut b) = (DVec3::ZERO, DVec3::ZERO);
        for (j, &src) in self.source.iter().enumerate() {
            if j == k || !src {
                continue;
            }
            let qj = self.scn.particles[self.members[j]].particle.charge;
            if qj == 0.0 {
                continue;
            }
            let (ej, bj) = match acc.map(|a| a[j]) {
                Some(cont) => {
                    let (e, b, _) = cont.fields(
                        qj,
                        self.scn.c,
                        x,
                        t - t_y,
                        Self::x(y, j),
                        self.kin[j].velocity(self.p(y, j)),
                    );
                    (e, b)
                }
                None => {
                    let vj = self.kin[j].velocity(self.p(y, j));
                    let rj = Self::x(y, j) + vj * (t - t_y);
                    heaviside_fields(qj, self.scn.c, x, rj, vj)
                }
            };
            e += ej;
            b += bj;
        }
        // Drained charges fading in their cups: the field of their uniform motion, with the
        // charge at its retarded time on that motion.
        for &(j, qf, fd) in &self.fading {
            if j == self.members[k] {
                continue;
            }
            let (ej, bj) = fd.quasi_static_fields(qf, self.scn.c, x, t);
            e += ej;
            b += bj;
        }
        // Charges at rest where absorbed particles stopped: their Coulomb fields.
        for &(j, qs, xs) in &self.stopped {
            if j != self.members[k] {
                let d = x - xs;
                let r2 = d.length_squared();
                e += d * (qs / (r2 * r2.sqrt()));
            }
        }
        (e, b)
    }

    /// Electric field at `x`, `t` (member `k` there) of the charges already absorbed that
    /// still act: those fading in a screening cup, and with the retarded interaction also
    /// a drained one whose removal has not yet reached `x` (not those stopped on bodies,
    /// which stay in the energy books as static charges). For the energy budget: their
    /// work on the particles present (only `E` does work).
    fn absorbed_fields(&self, k: usize, x: DVec3, t: f64) -> DVec3 {
        let me = self.members[k];
        let mut e = DVec3::ZERO;
        if let Some(past) = self.past {
            let guards = past.view();
            let view = guards.view();
            // A scratch copy of the receiver's retarded-time guesses: the dynamics' own
            // searches start from theirs, which these evaluations must not move.
            let mut row = past.guess[me].lock().expect("not poisoned").clone();
            for j in 0..self.scn.particles.len() {
                if j != me
                    && view.t_off[j].is_finite()
                    && !view.stays[j]
                    && self.scn.particles[j].particle.charge != 0.0
                {
                    e += view.fields_of(self.scn, &mut row, j, x, t).0;
                }
            }
        } else if self.quasi_static() {
            for &(j, qf, fd) in &self.fading {
                if j != me {
                    e += fd.quasi_static_fields(qf, self.scn.c, x, t).0;
                }
            }
        } else if self.scn.interact {
            for &(j, qf, fd) in &self.fading {
                if j != me {
                    let d = x - fd.position(t);
                    let r2 = d.length_squared();
                    e += d * (qf * fd.factor(t) / (r2 * r2.sqrt()));
                }
            }
        }
        e
    }

    /// Sum of the retarded fields of the other particles at `x`, `t` (member `k`).
    fn retarded_fields(
        &self,
        k: usize,
        x: DVec3,
        t: f64,
        view: Option<&View<'_>>,
    ) -> (DVec3, DVec3) {
        let (mut e, mut b) = (DVec3::ZERO, DVec3::ZERO);
        if let Some(past) = self.past {
            let i = self.members[k];
            let guards;
            let own;
            let view = if let Some(v) = view {
                v
            } else {
                guards = past.view();
                own = guards.view();
                &own
            };
            let mut row = past.guess[i].lock().expect("not poisoned");
            for j in 0..self.scn.particles.len() {
                if j != i && self.scn.particles[j].particle.charge != 0.0 {
                    let (ej, bj) = view.fields_of(self.scn, &mut row, j, x, t);
                    e += ej;
                    b += bj;
                }
            }
        }
        (e, b)
    }

    /// Landau–Lifshitz force on member `k` in the total field: external (with its image in
    /// the metal, which moves with it) plus the other particles' retarded fields.
    fn radiation_reaction_force(
        &self,
        y: &[f64],
        k: usize,
        t: f64,
        view: Option<&View<'_>>,
    ) -> DVec3 {
        let q = self.scn.particles[self.members[k]].particle.charge;
        let t_y = t;
        let fields = |x: DVec3, t: f64| {
            let f = self.scn.field.sample(x, t);
            let e_self = self.scn.field.self_field(x, q).0;
            let (e, b) = if self.quasi_static() {
                self.quasi_static_fields(y, t_y, k, x, t, None)
            } else {
                self.retarded_fields(k, x, t, view)
            };
            (f.e + e_self + e, f.b + b)
        };
        landau_lifshitz(q, &self.kin[k], self.p(y, k), fields, Self::x(y, k), t)
    }

    /// Total energy of the source members: kinetic, external potential, moments, and
    /// their pair interaction.
    fn energy(&self, y: &[f64], t: f64) -> f64 {
        let mut w = 0.0;
        for k in 0..self.members.len() {
            if !self.source[k] {
                continue;
            }
            let part = self.scn.particles[self.members[k]].particle;
            let x = Self::x(y, k);
            let f = self.scn.field.sample(x, t);
            w += self.kin[k].kinetic_energy(self.p(y, k))
                + part.charge * f.phi
                + 0.5 * part.charge * self.scn.field.self_field(x, part.charge).1
                - part.moment * f.b.z;
            if self.scn.interact && self.past.is_none() {
                for j in k + 1..self.members.len() {
                    if self.source[j] {
                        let qj = self.scn.particles[self.members[j]].particle.charge;
                        w += part.charge * qj / (x - Self::x(y, j)).length();
                    }
                }
                for &(j, qs, xs) in &self.stopped {
                    if j != self.members[k] {
                        w += part.charge * qs / (x - xs).length();
                    }
                }
            }
        }
        w
    }
}

impl<F: FieldSolver> OdeSystem for BeamOde<'_, F> {
    fn dim(&self) -> usize {
        6 * self.members.len()
    }

    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        let acc = self.continuations(y, t);
        let n = self.members.len();
        // With the retarded interaction every receiver's force is a sum of retarded
        // fields (the expensive part): computed in parallel, each receiver's own sum in a
        // fixed order and with its own guess row, so the result is deterministic.
        let forces: Vec<DVec3> = if let (Some(past), true) = (self.past, n > 1) {
            use rayon::prelude::*;
            let start = std::time::Instant::now();
            let guards = past.view();
            let view = guards.view();
            let f: Vec<DVec3> = (0..n)
                .into_par_iter()
                .map(|k| self.force_in(y, k, t, acc.as_deref(), true, Some(&view)))
                .collect();
            #[allow(clippy::cast_possible_truncation)]
            past.rhs_nanos
                .fetch_add(start.elapsed().as_nanos() as u64, Relaxed);
            f
        } else {
            (0..n)
                .map(|k| self.force_with(y, k, t, acc.as_deref(), true))
                .collect()
        };
        for (k, f) in forces.into_iter().enumerate() {
            let v = self.kin[k].velocity(self.p(y, k));
            let dp = f / self.p_ref;
            dy[6 * k..6 * k + 6].copy_from_slice(&[v.x, v.y, v.z, dp.x, dp.y, dp.z]);
        }
    }
}

/// State of one particle during the run.
#[derive(Clone, Debug)]
enum Phase {
    Flying,
    /// Finished; its continued trajectory is followed to find the penetration depth
    /// into boundary `event` (index into the particle's events), `depth` so far; it ends
    /// when the depth stops decreasing or reaches `GHOST_DEPTH`. With `real` it is the
    /// particle itself flying on (`Fate::Pass`): a source, and `Free` afterwards.
    Ghost {
        event: usize,
        depth: f64,
        real: bool,
    },
    /// Flew on past a boundary that is not physical (`Fate::Pass`): integrated and acting
    /// on the others, without events; its outcome is decided.
    Free,
    Done,
}

struct Track<'a> {
    traj: Trajectory,
    phase: Phase,
    margin: Vec<f64>,
    g_prev: Vec<f64>,
    gates: GateTracker<'a>,
    /// Largest radiation-reaction force and largest other force so far
    /// (`Trajectory::reaction_ratio_max`).
    rr_max: f64,
    rest_max: f64,
}

/// Runs a beam flight.
pub fn run_beam<F: FieldSolver>(scn: &BeamScenario<F>, rs: &RunSettings) -> BeamRun {
    run_beam_observed(scn, rs, |_, _, _| {})
}

/// As `run_beam`, calling `observer(dense, members, p_ref)` after every accepted step
/// (member `k` of the state is particle `members[k]`; its position is component
/// `6k..6k+3`, its momentum `p_ref ×` components `6k+3..6k+6`). Ghosts are members too:
/// their states after their end time are not part of their trajectories.
pub fn run_beam_observed<F: FieldSolver>(
    scn: &BeamScenario<F>,
    rs: &RunSettings,
    mut observer: impl FnMut(&Dense, &[usize], f64),
) -> BeamRun {
    run_beam_cancellable(scn, rs, |d, m, p| {
        observer(d, m, p);
        true
    })
    .expect("never cancelled")
}

/// As `run_beam_observed`, but the observer returns whether to go on; `None` if it
/// stopped the run.
pub fn run_beam_cancellable<F: FieldSolver>(
    scn: &BeamScenario<F>,
    rs: &RunSettings,
    mut observer: impl FnMut(&Dense, &[usize], f64) -> bool,
) -> Option<BeamRun> {
    let n = scn.particles.len();
    let mut events: Vec<Event> = (0..scn.obstacles.len()).map(Event::Obstacle).collect();
    if scn.bounds.is_some() {
        events.push(Event::Bounds);
    }
    // One detector event for all; a particle without a detector never reaches it.
    if scn.particles.iter().any(|b| b.detector.is_some()) {
        events.push(Event::Detector);
    }
    let event_value = |ev: Event, i: usize, x: DVec3| match ev {
        Event::Obstacle(o) => {
            scn.obstacles[o].signed_distance(x) - scn.particles[i].particle.radius
        }
        Event::Bounds => -scn.bounds.expect("bounds").signed_distance(x),
        Event::Detector => scn.particles[i]
            .detector
            .map_or(f64::INFINITY, |d| d.signed_distance(x)),
    };
    #[allow(clippy::cast_precision_loss)]
    let p_ref = {
        let s: f64 = scn.particles.iter().map(|b| b.p0.length()).sum::<f64>() / n.max(1) as f64;
        if s > 0.0 { s } else { 1.0 }
    };

    let mut tracks: Vec<Track> = scn
        .particles
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let kin = Kinematics::new(b.particle.mass, scn.c);
            let start = Sample {
                t: 0.0,
                x: b.x0,
                p: b.p0,
            };
            let g: Vec<f64> = events.iter().map(|&e| event_value(e, i, b.x0)).collect();
            Track {
                traj: Trajectory {
                    outcome: Outcome::Timeout,
                    end: start,
                    samples: if rs.record { vec![start] } else { Vec::new() },
                    stats: Stats::default(),
                    energy_initial: f64::NAN,
                    kinetic_initial: kin.kinetic_energy(b.p0),
                    energy_max_abs_error: f64::NAN,
                    radiated_energy: 0.0,
                    radiation_work: 0.0,
                    reaction_ratio_max: 0.0,
                    margins: None,
                    radiation: None,
                    emission: Vec::new(),
                },
                phase: Phase::Flying,
                margin: g.clone(),
                g_prev: g,
                gates: GateTracker::new(&scn.gates, b.x0),
                rr_max: 0.0,
                rest_max: 0.0,
            }
        })
        .collect();
    // Particles launched inside a boundary end at once (an arrival is judged below, once
    // the record exists); a boundary whose distance is not a number fails the particle.
    for t in &mut tracks {
        if t.g_prev.iter().any(|g| g.is_nan()) {
            t.traj.outcome = Outcome::Failed(crate::integrator::Error::NonFiniteEvent);
            t.phase = Phase::Done;
        } else if let Some(k) = t.g_prev.iter().position(|&g| g <= 0.0) {
            t.traj.outcome = events[k].outcome();
            t.phase = Phase::Done;
        }
    }

    let settings = Settings {
        rtol: rs.rtol,
        atol: rs.atol,
        max_steps: rs.max_steps,
        dense: true,
        ..Settings::default()
    };
    let mut stats = Stats::default();
    let mut restarts = 0;
    let mut energy_err: f64 = 0.0;
    let mut t_now = 0.0;
    // Current state of every particle that is flying or a ghost.
    let mut states: Vec<(DVec3, DVec3)> = scn.particles.iter().map(|b| (b.x0, b.p0)).collect();
    let mut failed: Option<crate::integrator::dop853::Error> = None;
    let retarded = scn.interact && scn.c.is_finite() && scn.retarded;
    let quasi_static = scn.interact && scn.c.is_finite() && !scn.retarded;
    // Per particle: time integrals of the neglected and of the kept interaction fields.
    let mut neglected = vec![(0.0f64, 0.0f64); n];
    // Largest speed (in units of c) of any particle so far, for pruning the record.
    let mut beta_max: f64 = 0.0;
    // Where absorbed particles stopped, with their charge staying (`Fate::Stop`).
    let mut stopped: Vec<Option<DVec3>> = vec![None; n];
    // Dense samples of every charge for radiation goals (PHYSICS.md §3.4): a receiver sees
    // the fields of all particles together, so a goal measures them all, coherently.
    let record_emission = scn.c.is_finite()
        && scn
            .particles
            .iter()
            .any(|b| b.acceptance.is_some_and(|a| a.radiation.is_some()));
    let mut emission: Vec<Vec<Emission>> = vec![Vec::new(); n];
    let mut fades: Vec<Option<Fade>> = vec![None; n];
    // Energy budget (recorded runs): samples, and the energy absorbed so far.
    let mut energy: Vec<EnergySample> = Vec::new();
    let mut absorbed = 0.0;
    let past = Past {
        segments: RwLock::new(Vec::new()),
        stays: RwLock::new(vec![false; n]),
        fades: RwLock::new(vec![None; n]),
        t_off: RwLock::new(vec![f64::INFINITY; n]),
        x_off: RwLock::new(scn.particles.iter().map(|b| b.x0).collect()),
        kin: scn
            .particles
            .iter()
            .map(|b| Kinematics::new(b.particle.mass, scn.c))
            .collect(),
        p_ref,
        a0: if retarded {
            launch_accelerations(scn)
        } else {
            Vec::new()
        },
        guess: (0..n).map(|_| Mutex::new(vec![f64::NAN; n])).collect(),
        solves: AtomicU64::new(0),
        rhs_nanos: AtomicU64::new(0),
        evals: AtomicU64::new(0),
        newton: AtomicU64::new(0),
        fallbacks: AtomicU64::new(0),
    };
    let log = std::env::var_os("EM_BEAM_LOG").is_some();
    let log_start = std::time::Instant::now();
    let mut log_steps = 0u64;
    // Particles that end at launch: drained, fading in a screening cup, stopped there, or
    // flying on. One launched inside its detector has arrived at t = 0 and is judged by
    // the gates and the detector's conditions like any entry.
    for (i, t) in tracks.iter_mut().enumerate() {
        if !matches!(t.phase, Phase::Done) {
            continue;
        }
        let Some(k) = t.g_prev.iter().position(|&g| g <= 0.0) else {
            // Failed (NaN): gone at once.
            past.t_off.write().expect("not poisoned")[i] = 0.0;
            continue;
        };
        match events[k].fate(&scn.fates) {
            Fate::Drain => past.t_off.write().expect("not poisoned")[i] = 0.0,
            Fate::Cup => {
                past.t_off.write().expect("not poisoned")[i] = 0.0;
                if events[k] == Event::Detector
                    && let Some(region) = scn.particles[i].detector.as_ref()
                {
                    let b = &scn.particles[i];
                    let v0 = Kinematics::new(b.particle.mass, scn.c).velocity(b.p0);
                    let fd = Fade::at_entry(region, b.x0, v0, 0.0);
                    fades[i] = Some(fd);
                    past.fades.write().expect("not poisoned")[i] = Some(fd);
                }
            }
            Fate::Stop => {
                past.t_off.write().expect("not poisoned")[i] = 0.0;
                past.stays.write().expect("not poisoned")[i] = true;
                stopped[i] = Some(scn.particles[i].x0);
            }
            Fate::Pass => t.phase = Phase::Free,
        }
        if t.traj.outcome == Outcome::Arrived {
            let kin = Kinematics::new(scn.particles[i].particle.mass, scn.c);
            judge_beam_arrival(scn, t, i, &kin, &emission);
        }
    }

    if rs.record {
        let present: Vec<Present> = scn
            .particles
            .iter()
            .enumerate()
            .filter(|(i, _)| !matches!(tracks[*i].phase, Phase::Done) || stopped[*i].is_some())
            .map(|(i, b)| {
                let tk = if stopped[i].is_some() {
                    0.0
                } else {
                    tracks[i].traj.kinetic_initial
                };
                (b.particle.charge, b.particle.moment, b.x0, tk)
            })
            .collect();
        let (kinetic, potential, interaction) = budget(scn, 0.0, &present);
        energy.push(EnergySample {
            t: 0.0,
            kinetic,
            potential,
            interaction,
            ..EnergySample::default()
        });
    }

    'segments: loop {
        let members: Vec<usize> = (0..n)
            .filter(|&i| !matches!(tracks[i].phase, Phase::Done))
            .collect();
        // The flight is over when no particle is left in play (only particles flying on
        // outside the arena, which nothing is left to feel).
        if !members
            .iter()
            .any(|&i| !matches!(tracks[i].phase, Phase::Free))
        {
            break;
        }
        let source: Vec<bool> = members
            .iter()
            .map(|&i| {
                matches!(
                    tracks[i].phase,
                    Phase::Flying | Phase::Free | Phase::Ghost { real: true, .. }
                )
            })
            .collect();
        let ode = BeamOde {
            scn,
            members: members.clone(),
            source: source.clone(),
            kin: members
                .iter()
                .map(|&i| Kinematics::new(scn.particles[i].particle.mass, scn.c))
                .collect(),
            p_ref,
            past: retarded.then_some(&past),
            stopped: stopped
                .iter()
                .enumerate()
                .filter_map(|(i, s)| {
                    let q = scn.particles[i].particle.charge;
                    s.filter(|_| q != 0.0).map(|x| (i, q, x))
                })
                .collect(),
            fading: fades
                .iter()
                .enumerate()
                .filter_map(|(i, f)| {
                    let q = scn.particles[i].particle.charge;
                    f.filter(|_| q != 0.0).map(|fd| (i, q, fd))
                })
                .collect(),
        };
        let mut y0 = Vec::with_capacity(6 * members.len());
        for &i in &members {
            let (x, p) = states[i];
            let s = p / p_ref;
            y0.extend_from_slice(&[x.x, x.y, x.z, s.x, s.y, s.z]);
        }
        let mut energy0 = ode.energy(&y0, t_now);
        // One step budget for the whole flight, over all segments.
        if stats.n_step >= rs.max_steps {
            failed = Some(crate::integrator::dop853::Error::MaxStepsReached);
            break 'segments;
        }
        let segment = Settings {
            max_steps: rs.max_steps - stats.n_step,
            ..settings
        };
        let mut int = Dop853::new(&ode, t_now, &y0, segment);
        let mut y = vec![0.0; ode.dim()];
        loop {
            let t_a = int.t();
            // Retarded interaction: the step is capped (`h_cap` below) so that every
            // retarded time a stage needs lies in the record or at most `EXTRAPOLATION` of
            // the previous step's lengths past its end.
            if retarded {
                past.prune(&members, int.y(), &ode, scn.c, &mut beta_max, t_a);
            }
            let h_cap = if retarded {
                // Closest pair and fastest closing speed (upper bound |v_i − v_j|).
                let (mut r_min, mut v_close) = (f64::INFINITY, 0.0f64);
                for (k, &i) in members.iter().enumerate() {
                    if scn.particles[i].particle.charge == 0.0 {
                        continue;
                    }
                    let xi = BeamOde::<F>::x(int.y(), k);
                    let vi = ode.kin[k].velocity(ode.p(int.y(), k));
                    for (l, &j) in members.iter().enumerate() {
                        if j != i && source[l] && scn.particles[j].particle.charge != 0.0 {
                            r_min = r_min.min((xi - BeamOde::<F>::x(int.y(), l)).length());
                            let vj = ode.kin[l].velocity(ode.p(int.y(), l));
                            v_close = v_close.max((vi - vj).length());
                        }
                    }
                }
                // Retarded times inside the step come from the previous step's polynomial,
                // extrapolated at most EXTRAPOLATION of its lengths L ahead. A stage at
                // t_a + s (s ≤ h) needs t_r = t_a + s − R/c with R ≥ r_min − w h (closing
                // speed w, taken as twice the fastest at the step's start, never above
                // 2c): t_r ≤ t_a + E L holds for h ≤ (E L + r_min/c) / (1 + w/c). The first
                // step (no record yet) keeps every retarded time in the past.
                let previous = past
                    .segments
                    .read()
                    .expect("not poisoned")
                    .last()
                    .map_or(0.0, |s| s.t_stop - s.t_start);
                let w = (2.0 * v_close).min(2.0 * scn.c);
                (EXTRAPOLATION * previous + r_min / scn.c) / (1.0 + w / scn.c)
            } else {
                f64::INFINITY
            };
            if log {
                log_steps += 1;
                if log_steps.is_multiple_of(50) {
                    let n_s = past.solves.load(Relaxed).max(1);
                    #[allow(clippy::cast_precision_loss)]
                    let per = |c: &AtomicU64| c.load(Relaxed) as f64 / n_s as f64;
                    eprintln!(
                        "beam log: step {log_steps}, wall {:.3} s of which right-hand sides {:.3} s, t = {:.6}, last h = {:.3e}, cap {h_cap:.3e}, flying {}, solves {n_s}, per solve: {:.2} evals + {:.2} Newton, fallbacks {:.3}, record {} steps",
                        log_start.elapsed().as_secs_f64(),
                        rhs_seconds(&past.rhs_nanos),
                        int.t(),
                        int.h_next(),
                        members.len(),
                        per(&past.evals),
                        per(&past.newton),
                        per(&past.fallbacks),
                        past.segments.read().expect("not poisoned").len()
                    );
                }
            }
            let reached_end = match int.step(&ode, scn.t_max, h_cap) {
                Ok(done) => done,
                Err(e) => {
                    failed = Some(e);
                    break 'segments;
                }
            };
            let t_step = int.t();
            let dense = int.dense();
            if !observer(dense, &members, p_ref) {
                return None;
            }
            let pos = |k: usize, t: f64| {
                DVec3::new(
                    dense.eval_component(6 * k, t),
                    dense.eval_component(6 * k + 1, t),
                    dense.eval_component(6 * k + 2, t),
                )
            };
            let mom = |k: usize, t: f64| {
                DVec3::new(
                    dense.eval_component(6 * k + 3, t),
                    dense.eval_component(6 * k + 4, t),
                    dense.eval_component(6 * k + 5, t),
                ) * p_ref
            };

            // Contacts between flying particles with a radius (rigid spheres): the earliest
            // one ends the step there. Pairs already touching or separating are skipped
            // (just after a collision they move apart).
            let mut contact: Option<(f64, usize, usize)> = None;
            let speed = |k: usize, t: f64| ode.kin[k].velocity(mom(k, t)).length();
            for (k, &i) in members.iter().enumerate() {
                let ri = scn.particles[i].particle.radius;
                if !matches!(tracks[i].phase, Phase::Flying) {
                    continue;
                }
                for (l, &j) in members.iter().enumerate().skip(k + 1) {
                    let rsum = ri + scn.particles[j].particle.radius;
                    if rsum <= 0.0 || !matches!(tracks[j].phase, Phase::Flying) {
                        continue;
                    }
                    let mut g = |t: f64| (pos(k, t) - pos(l, t)).length() - rsum;
                    let (ga, gb) = (g(t_a), g(t_step));
                    if ga <= 0.0 {
                        continue;
                    }
                    let v_rel = ((0..=4)
                        .map(|s| {
                            let t = t_a + (t_step - t_a) * f64::from(s) * 0.25;
                            speed(k, t) + speed(l, t)
                        })
                        .fold(0.0, f64::max)
                        * 1.25)
                        .min(2.0 * scn.c);
                    if let Some(t) = first_crossing(&mut g, t_a, t_step, ga, gb, v_rel).time
                        && contact.is_none_or(|(tc, _, _)| t < tc)
                    {
                        contact = Some((t, k, l));
                    }
                }
            }
            // The step, cut at the contact if there is one.
            let t_b = contact.map_or(t_step, |(t, _, _)| t);
            let reached_end = reached_end && contact.is_none();
            let end_position = |k: usize| {
                if contact.is_some() {
                    pos(k, t_b)
                } else {
                    BeamOde::<F>::x(int.y(), k)
                }
            };

            // Event values at the step end, and the earliest event of a flying particle.
            let mut g_next: Vec<Vec<f64>> = vec![Vec::new(); members.len()];
            let mut v_max = vec![0.0; members.len()];
            let mut first: Option<(f64, usize, usize)> = None; // (t, member, event)
            let mut nan = false;
            for (k, &i) in members.iter().enumerate() {
                v_max[k] = ((0..=4)
                    .map(|s| {
                        let t = t_a + (t_b - t_a) * f64::from(s) * 0.25;
                        ode.kin[k].velocity(mom(k, t)).length()
                    })
                    .fold(0.0, f64::max)
                    * 1.25)
                    .min(scn.c);
                if !matches!(tracks[i].phase, Phase::Flying) {
                    continue;
                }
                let x_b = end_position(k);
                g_next[k] = events.iter().map(|&ev| event_value(ev, i, x_b)).collect();
                for (e, &ev) in events.iter().enumerate() {
                    let ga = tracks[i].g_prev[e];
                    // At the boundary already at the step's start: it reached it at the
                    // same moment as the event that restarted the integration (two
                    // particles entering their detectors together), so its event is now.
                    let time = if ga <= 0.0 {
                        Some(t_a)
                    } else {
                        let mut g = |t: f64| event_value(ev, i, pos(k, t));
                        let r = first_crossing(&mut g, t_a, t_b, ga, g_next[k][e], v_max[k]);
                        nan |= r.nan;
                        r.time
                    };
                    if let Some(t) = time
                        && first.is_none_or(|(tf, _, _)| t < tf)
                    {
                        first = Some((t, k, e));
                    }
                }
            }
            if nan {
                // A boundary's distance is not a number: the step cannot be certified.
                failed = Some(crate::integrator::Error::NonFiniteEvent);
                stats = add_stats(stats, int.stats());
                break 'segments;
            }
            let t_end = first.map_or(t_b, |(t, _, _)| t);
            if retarded {
                past.segments.write().expect("not poisoned").push(Segment {
                    dense: dense.clone(),
                    t_start: t_a,
                    t_stop: t_end,
                    members: members.clone(),
                });
            }
            // Quasi-static interaction: estimated neglected fields (the jerk the sources'
            // continued past misses: the change of their acceleration over the step less the
            // continuation's own), relative to the fields kept, at the step's end.
            if quasi_static && t_end > t_a {
                let mut y_end = vec![0.0; ode.dim()];
                dense.eval(t_end, &mut y_end);
                let accel = |k: usize, t: f64| {
                    let dp = DVec3::new(
                        dense.eval_derivative_component(6 * k + 3, t),
                        dense.eval_derivative_component(6 * k + 4, t),
                        dense.eval_derivative_component(6 * k + 5, t),
                    ) * p_ref;
                    ode.kin[k].acceleration(mom(k, t), dp)
                };
                let continued = ode.continuations(&y_end, t_end).unwrap_or_default();
                let jerk: Vec<DVec3> = (0..members.len())
                    .map(|k| (accel(k, t_end) - accel(k, t_a)) / (t_end - t_a))
                    .collect();
                for (k, &i) in members.iter().enumerate() {
                    let qi = scn.particles[i].particle.charge;
                    if !matches!(tracks[i].phase, Phase::Flying) || qi == 0.0 {
                        continue;
                    }
                    let xi = BeamOde::<F>::x(&y_end, k);
                    let (mut missing, mut present) = (0.0, 0.0);
                    for (l, &j) in members.iter().enumerate() {
                        let qj = scn.particles[j].particle.charge;
                        if l != k && source[l] && qj != 0.0 {
                            let (xj, pj) = (BeamOde::<F>::x(&y_end, l), ode.p(&y_end, l));
                            let g = ode.kin[l].gamma(pj);
                            let vj = ode.kin[l].velocity(pj);
                            // The light delay grows ahead of a moving source: about
                            // R / (c κ) with κ = 1 − n·β, the Doppler factor.
                            let d = xi - xj;
                            let r = d.length();
                            let kappa = (1.0 - d.dot(vj) / (r * scn.c)).max(1e-6);
                            let delay = r / (scn.c * kappa);
                            // The jerk the continuation misses: for the motion in the
                            // fields (weight w for this pair) the difference from its own.
                            let (w, own) = match continued.get(l) {
                                Some(cont @ Continuation::Fields(_)) => {
                                    let (_, _, w) = cont.fields(qj, scn.c, xi, 0.0, xj, vj);
                                    (w, cont.jerk())
                                }
                                _ => (0.0, DVec3::ZERO),
                            };
                            missing += qj.abs() * g * g * (jerk[l] - own * w).length()
                                / (scn.c.powi(3) * kappa);
                            // Beyond the constant-acceleration range of a past continued
                            // with its acceleration (`tapered`: |a τ| > 0.7 Δv) that
                            // acceleration fades: first-order error there.
                            let a_now = accel(l, t_end).length();
                            if w < 1.0 && a_now * delay > TAPER_START * taper_reach(vj, scn.c) {
                                missing +=
                                    (1.0 - w) * qj.abs() * g * g * a_now / (scn.c * scn.c * r);
                            }
                            present += heaviside_fields(qj, scn.c, xi, xj, vj).0.length();
                        }
                    }
                    // Weighted with the step length: integrals over the flight.
                    let dt = t_end - t_a;
                    neglected[i].0 += missing * dt;
                    neglected[i].1 += present * dt;
                }
            }
            // Largest radiation-reaction force over the largest other force (the
            // Landau–Lifshitz treatment needs it small), at the steps' ends.
            if scn.radiation_reaction && scn.c.is_finite() {
                let mut y_end = vec![0.0; ode.dim()];
                dense.eval(t_end, &mut y_end);
                for (k, &i) in members.iter().enumerate() {
                    if matches!(tracks[i].phase, Phase::Flying) && ode.reacts(k) {
                        let rr = ode.radiation_reaction_force(&y_end, k, t_end, None);
                        let rest = (ode.force(&y_end, k, t_end) - rr).length();
                        let tr = &mut tracks[i];
                        tr.rr_max = tr.rr_max.max(rr.length());
                        tr.rest_max = tr.rest_max.max(rest);
                        if tr.rest_max > 0.0 {
                            tr.traj.reaction_ratio_max = tr.rr_max / tr.rest_max;
                        }
                    }
                }
            }
            // Radiated energy of the flying particles (Liénard power, three-point
            // Gauss–Legendre on the step), from the force given by the derivative of the
            // dense output (the single-particle runner's trapezoidal rule was 0.24 % off
            // on a peaked flight: test R5).
            if scn.c.is_finite() {
                for (k, &i) in members.iter().enumerate() {
                    if matches!(tracks[i].phase, Phase::Flying) {
                        let q = scn.particles[i].particle.charge;
                        let power = |t: f64| {
                            let dp = DVec3::new(
                                dense.eval_derivative_component(6 * k + 3, t),
                                dense.eval_derivative_component(6 * k + 4, t),
                                dense.eval_derivative_component(6 * k + 5, t),
                            ) * p_ref;
                            larmor_power(&ode.kin[k], q, mom(k, t), dp)
                        };
                        let h = t_end - t_a;
                        let r = 0.5 * (0.6f64).sqrt();
                        tracks[i].traj.radiated_energy += h
                            * (5.0 / 18.0 * power(t_a + (0.5 - r) * h)
                                + 8.0 / 18.0 * power(t_a + 0.5 * h)
                                + 5.0 / 18.0 * power(t_a + (0.5 + r) * h));
                    }
                }
            }

            // Emission samples of the charges in flight (also beyond the arena), up to
            // `t_end`: 8 pieces per step, an even number, as in the single-particle runner;
            // the acceleration from the derivative of the dense output.
            if record_emission {
                let h = t_end - t_a;
                for (k, &i) in members.iter().enumerate() {
                    let flies = matches!(
                        tracks[i].phase,
                        Phase::Flying | Phase::Free | Phase::Ghost { real: true, .. }
                    );
                    if !flies || scn.particles[i].particle.charge == 0.0 {
                        continue;
                    }
                    for piece in 0..=8 {
                        let t = t_a + h * (f64::from(piece) / 8.0);
                        if emission[i].last().is_some_and(|e| t <= e.0) {
                            continue;
                        }
                        let p = mom(k, t);
                        let dp = DVec3::new(
                            dense.eval_derivative_component(6 * k + 3, t),
                            dense.eval_derivative_component(6 * k + 4, t),
                            dense.eval_derivative_component(6 * k + 5, t),
                        ) * p_ref;
                        emission[i].push((
                            t,
                            pos(k, t),
                            ode.kin[k].velocity(p),
                            ode.kin[k].acceleration(p, dp),
                        ));
                    }
                }
            }

            // Gates passed by the flying particles up to `t_end`.
            for (k, &i) in members.iter().enumerate() {
                if matches!(tracks[i].phase, Phase::Flying) {
                    tracks[i].gates.advance(
                        |t| (pos(k, t), mom(k, t)),
                        &ode.kin[k],
                        t_a,
                        t_end,
                        v_max[k],
                        rs.margins,
                    );
                }
            }

            // Margins of flying particles (up to the event), depths of ghosts.
            let mut ghosts_done = false;
            for (k, &i) in members.iter().enumerate() {
                match tracks[i].phase.clone() {
                    Phase::Flying => {
                        for (e, &ev) in events.iter().enumerate() {
                            let g = |t: f64| event_value(ev, i, pos(k, t));
                            let gb = g_next[k][e];
                            if first.is_some_and(|(_, kk, ee)| kk == k && ee == e) {
                                // Penetration depth, followed past this step as a ghost.
                                let (t_min, m) = minimize_on(&g, t_a, t_b);
                                tracks[i].margin[e] = tracks[i].margin[e].min(m.max(-GHOST_DEPTH));
                                let unfinished = minimum_at_end(t_min, t_a, t_b);
                                let real = ev.fate(&scn.fates) == Fate::Pass;
                                tracks[i].phase = Phase::Ghost {
                                    event: e,
                                    depth: m,
                                    real,
                                };
                                if !rs.margins || !unfinished || m <= -GHOST_DEPTH {
                                    finish_ghost(&mut tracks[i]);
                                }
                            } else if rs.margins {
                                let m = if first.is_some() {
                                    minimize_on(&g, t_a, t_end).1
                                } else {
                                    let gp = tracks[i].g_prev[e];
                                    let lower = 0.5 * (gp + gb - v_max[k] * (t_b - t_a));
                                    if lower < 2.0 * MARGIN_SAFE {
                                        minimize_on(&g, t_a, t_b).1
                                    } else {
                                        gp.min(gb)
                                    }
                                };
                                tracks[i].margin[e] = tracks[i].margin[e].min(m);
                            }
                        }
                        // The next step starts at `t_end` (earlier than `t_b` if an event
                        // cut this one short, and the integration restarts there).
                        if first.is_some() {
                            let x = pos(k, t_end);
                            tracks[i].g_prev =
                                events.iter().map(|&ev| event_value(ev, i, x)).collect();
                        } else {
                            tracks[i].g_prev.clone_from(&g_next[k]);
                        }
                    }
                    Phase::Ghost { event, depth, real } => {
                        let g = |t: f64| event_value(events[event], i, pos(k, t));
                        let (t_min, m) = minimize_on(&g, t_a, t_b);
                        let depth = depth.min(m);
                        let unfinished = minimum_at_end(t_min, t_a, t_b);
                        tracks[i].phase = Phase::Ghost { event, depth, real };
                        if !unfinished || depth <= -GHOST_DEPTH {
                            finish_ghost(&mut tracks[i]);
                            ghosts_done = true;
                        }
                    }
                    Phase::Free | Phase::Done => {}
                }
            }

            // States at the end of this step (or at the event); samples of the particles
            // that were flying until then.
            dense.eval(t_end, &mut y);
            for (k, &i) in members.iter().enumerate() {
                let (x, p) = (BeamOde::<F>::x(&y, k), ode.p(&y, k));
                states[i] = (x, p);
                let was_flying = matches!(tracks[i].phase, Phase::Flying)
                    || first.is_some_and(|(_, kk, _)| kk == k);
                if was_flying {
                    let sample = Sample { t: t_end, x, p };
                    if rs.record {
                        tracks[i].traj.samples.push(sample);
                    }
                    tracks[i].traj.end = sample;
                }
            }
            let e = ode.energy(&y, t_end);
            // While a drained charge fades in its cup the flying particles' energy is not
            // conserved (the fading charge's field does work, and energy goes into the cup):
            // the drift is measured only from where the fades have become negligible.
            if ode.fading.iter().any(|(_, _, fd)| fd.factor(t_end) > 1e-15) {
                energy0 = e;
            } else {
                energy_err = energy_err.max((e - energy0).abs() / energy0.abs().max(1e-300));
            }

            // The work done on the particles present by charges already absorbed that
            // still act on them (fading in a screening cup; with the retarded interaction
            // also the field of a drained one still on its way): energy that their
            // absorption, booked at entry with the whole pair interaction, did not take
            // out. Three-point Gauss–Legendre on the step, from the dense output.
            if rs.record
                && t_end > t_a
                && (!ode.fading.is_empty()
                    || (retarded
                        && past
                            .t_off
                            .read()
                            .expect("not poisoned")
                            .iter()
                            .any(|t| t.is_finite())))
            {
                let power = |t: f64| -> f64 {
                    let mut p = 0.0;
                    for (k, &i) in members.iter().enumerate() {
                        let q = scn.particles[i].particle.charge;
                        let real = matches!(
                            tracks[i].phase,
                            Phase::Flying | Phase::Free | Phase::Ghost { real: true, .. }
                        ) || first.is_some_and(|(_, kk, _)| kk == k);
                        if real && q != 0.0 {
                            let v = ode.kin[k].velocity(mom(k, t));
                            p += q * v.dot(ode.absorbed_fields(k, pos(k, t), t));
                        }
                    }
                    p
                };
                let h = t_end - t_a;
                let r = 0.5 * (0.6f64).sqrt();
                absorbed -= h
                    * (5.0 / 18.0 * power(t_a + (0.5 - r) * h)
                        + 8.0 / 18.0 * power(t_a + 0.5 * h)
                        + 5.0 / 18.0 * power(t_a + (0.5 + r) * h));
            }

            // Energy budget at t_end (the event's particle still counted as present; its
            // absorption is booked below).
            let present_at = |tracks: &[Track], stopped: &[Option<DVec3>], skip: Option<usize>| {
                let mut out: Vec<Present> = Vec::new();
                for (k, &i) in members.iter().enumerate() {
                    let real = matches!(
                        tracks[i].phase,
                        Phase::Flying | Phase::Free | Phase::Ghost { real: true, .. }
                    ) || first.is_some_and(|(_, kk, _)| kk == k);
                    if real && Some(i) != skip {
                        let part = scn.particles[i].particle;
                        let (x, p) = (BeamOde::<F>::x(&y, k), ode.p(&y, k));
                        out.push((part.charge, part.moment, x, ode.kin[k].kinetic_energy(p)));
                    }
                }
                for (i, st) in stopped.iter().enumerate() {
                    if let Some(x) = st
                        && Some(i) != skip
                    {
                        out.push((scn.particles[i].particle.charge, 0.0, *x, 0.0));
                    }
                }
                out
            };
            if rs.record {
                let present = present_at(&tracks, &stopped, None);
                let (kinetic, potential, interaction) = budget(scn, t_end, &present);
                energy.push(EnergySample {
                    t: t_end,
                    kinetic,
                    potential,
                    interaction,
                    absorbed,
                    radiated: tracks.iter().map(|t| t.traj.radiated_energy).sum(),
                });
            }

            if let Some((_, k, e)) = first {
                let i = members[k];
                tracks[i].traj.outcome = events[e].outcome();
                let fate = events[e].fate(&scn.fates);
                if rs.record && fate != Fate::Pass {
                    // The energy the absorption takes out: before minus after.
                    let before = present_at(&tracks, &stopped, None);
                    let mut after = present_at(&tracks, &stopped, Some(i));
                    if fate == Fate::Stop {
                        let part = scn.particles[i].particle;
                        after.push((part.charge, 0.0, tracks[i].traj.end.x, 0.0));
                    }
                    let sum = |v: &[Present]| {
                        let (a, b, c) = budget(scn, t_end, v);
                        a + b + c
                    };
                    absorbed += sum(&before) - sum(&after);
                    if let Some(last) = energy.last_mut() {
                        let (kinetic, potential, interaction) = budget(scn, t_end, &after);
                        *last = EnergySample {
                            kinetic,
                            potential,
                            interaction,
                            absorbed,
                            ..*last
                        };
                    }
                }
                match fate {
                    Fate::Drain | Fate::Cup => {
                        past.t_off.write().expect("not poisoned")[i] = t_end;
                        past.x_off.write().expect("not poisoned")[i] = tracks[i].traj.end.x;
                        // A cup is the detector's fate; other boundaries treat it as a drain.
                        let cup = fate == Fate::Cup && events[e] == Event::Detector;
                        if cup && let Some(region) = scn.particles[i].detector.as_ref() {
                            let kin = Kinematics::new(scn.particles[i].particle.mass, scn.c);
                            let v_end = kin.velocity(tracks[i].traj.end.p);
                            let fd = Fade::at_entry(region, tracks[i].traj.end.x, v_end, t_end);
                            fades[i] = Some(fd);
                            past.fades.write().expect("not poisoned")[i] = Some(fd);
                        }
                        // A drained charge: the quasi-static interaction drops its field at
                        // once; really it lingers at each other particle for the light
                        // time R/c: a neglected impulse of about |q| / (R c) on each. (In
                        // a cup the fade reaches the others with the light delay of its
                        // uniform motion, `Fade::quasi_static_fields`: nothing dropped.)
                        let qi = scn.particles[i].particle.charge;
                        if quasi_static && !cup && qi != 0.0 {
                            let xe = tracks[i].traj.end.x;
                            for (l, &j) in members.iter().enumerate() {
                                let flying = matches!(
                                    tracks[j].phase,
                                    Phase::Flying | Phase::Free | Phase::Ghost { real: true, .. }
                                );
                                if j != i && flying && scn.particles[j].particle.charge != 0.0 {
                                    let r = (BeamOde::<F>::x(&y, l) - xe).length();
                                    neglected[j].0 += qi.abs() / (r * scn.c);
                                }
                            }
                        }
                    }
                    Fate::Stop => {
                        past.t_off.write().expect("not poisoned")[i] = t_end;
                        past.x_off.write().expect("not poisoned")[i] = tracks[i].traj.end.x;
                        past.stays.write().expect("not poisoned")[i] = true;
                        stopped[i] = Some(tracks[i].traj.end.x);
                    }
                    Fate::Pass => {}
                }
                judge_beam_arrival(scn, &mut tracks[i], i, &ode.kin[k], &emission);
                t_now = t_end;
                stats = add_stats(stats, int.stats());
                restarts += 1;
                continue 'segments;
            }
            if let Some((_, k, l)) = contact {
                let (i, j) = (members[k], members[l]);
                let (x1, x2) = (states[i].0, states[j].0);
                let n = (x2 - x1).normalize();
                if let Some(jn) =
                    elastic_impulse(&ode.kin[k], states[i].1, &ode.kin[l], states[j].1, n)
                {
                    for (idx, kk, sign) in [(i, k, -1.0), (j, l, 1.0)] {
                        let p_old = states[idx].1;
                        let p_new = p_old + n * (sign * jn);
                        states[idx].1 = p_new;
                        // The impact's radiation, not modelled: estimated for a velocity
                        // change over the light crossing time of the sphere.
                        let part = scn.particles[idx].particle;
                        if scn.c.is_finite() && part.charge != 0.0 && part.radius > 0.0 {
                            let dv = ode.kin[kk].velocity(p_new) - ode.kin[kk].velocity(p_old);
                            tracks[idx].traj.radiated_energy +=
                                part.charge * part.charge * dv.length_squared()
                                    / (3.0 * part.radius * scn.c * scn.c);
                        }
                        let sample = Sample {
                            t: t_end,
                            x: states[idx].0,
                            p: p_new,
                        };
                        if rs.record {
                            tracks[idx].traj.samples.push(sample);
                        }
                        tracks[idx].traj.end = sample;
                    }
                }
                t_now = t_end;
                stats = add_stats(stats, int.stats());
                restarts += 1;
                continue 'segments;
            }
            if ghosts_done {
                t_now = t_b;
                stats = add_stats(stats, int.stats());
                restarts += 1;
                continue 'segments;
            }
            if reached_end {
                for &i in &members {
                    match tracks[i].phase {
                        Phase::Flying => {
                            tracks[i].traj.outcome = Outcome::Timeout;
                            tracks[i].phase = Phase::Done;
                        }
                        Phase::Ghost { .. } => finish_ghost(&mut tracks[i]),
                        Phase::Free | Phase::Done => {}
                    }
                }
                stats = add_stats(stats, int.stats());
                break 'segments;
            }
        }
    }
    if let Some(e) = failed {
        for t in &mut tracks {
            if matches!(t.phase, Phase::Flying) {
                t.traj.outcome = Outcome::Failed(e);
            }
        }
    }
    let mut trajectories: Vec<Trajectory> = tracks
        .into_iter()
        .zip(&scn.particles)
        .map(|(mut t, b)| {
            t.traj.stats = stats;
            if rs.margins {
                let acceptance = t.traj.margins.as_ref().and_then(|m| m.acceptance);
                let mut m = Margins {
                    obstacles: Vec::new(),
                    bounds: None,
                    detector: None,
                    acceptance,
                    gates: Vec::new(),
                    gate_acceptance: Vec::new(),
                };
                (m.gates, m.gate_acceptance) = t.gates.into_margins();
                for (&ev, &v) in events.iter().zip(&t.margin) {
                    match ev {
                        Event::Obstacle(_) => m.obstacles.push(v),
                        Event::Bounds => m.bounds = Some(v),
                        Event::Detector => m.detector = b.detector.map(|_| v),
                    }
                }
                t.traj.margins = Some(m);
            }
            t.traj
        })
        .collect();
    for (traj, e) in trajectories.iter_mut().zip(emission) {
        traj.emission = e;
    }
    Some(BeamRun {
        trajectories,
        stats,
        neglected_retardation: neglected
            .iter()
            .map(|&(m, p)| if p > 0.0 { m / p } else { 0.0 })
            .collect(),
        energy_max_rel_error: if retarded
            || quasi_static
            || (scn.radiation_reaction && scn.c.is_finite())
        {
            f64::NAN
        } else {
            energy_err
        },
        restarts,
        energy,
        fades,
    })
}

/// Elastic collision of two rigid spheres touching along `n` (unit, from particle 1 to
/// particle 2): the impulse `J n` on 2 and `−J n` on 1 with `J > 0` that conserves the
/// kinetic energy, the positive root of `T₁(p₁ − J n) + T₂(p₂ + J n) = T₁(p₁) + T₂(p₂)`
/// (convex in J, zero at J = 0 with a negative slope while they approach). For `c = ∞`
/// it is `2 μ (v₁ − v₂)·n`. `None` if they are not approaching.
pub fn elastic_impulse(
    k1: &Kinematics,
    p1: DVec3,
    k2: &Kinematics,
    p2: DVec3,
    n: DVec3,
) -> Option<f64> {
    let closing = (k1.velocity(p1) - k2.velocity(p2)).dot(n);
    if closing <= 0.0 {
        return None;
    }
    let before = k1.kinetic_energy(p1) + k2.kinetic_energy(p2);
    let f = |j: f64| k1.kinetic_energy(p1 - n * j) + k2.kinetic_energy(p2 + n * j) - before;
    // Bracket the root: the Newtonian value is exact for c = inf; relativistically the
    // root lies within a few times it.
    let (m1, m2) = (k1.mass, k2.mass);
    let newtonian = 2.0 * m1 * m2 / (m1 + m2) * closing;
    if !k1.c.is_finite() && !k2.c.is_finite() {
        return Some(newtonian);
    }
    let (mut lo, mut hi) = (0.5 * newtonian, 2.0 * newtonian);
    while f(lo) > 0.0 {
        lo *= 0.5;
    }
    while f(hi) < 0.0 {
        hi *= 2.0;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// The gates and the detector's conditions of particle `i` arriving now (at
/// `track.traj.end`, the moment of entry): `SkippedGate` while a gate is still missing,
/// `Rejected` outside the acceptance (the acceptance margin is recorded). A radiation
/// goal measures every charge's radiation up to now, their far fields added (the receiver
/// cannot tell them apart).
fn judge_beam_arrival<F>(
    scn: &BeamScenario<F>,
    track: &mut Track<'_>,
    i: usize,
    kin: &Kinematics,
    emission: &[Vec<Emission>],
) {
    if track.traj.outcome == Outcome::Arrived
        && let Some(g) = track.gates.missing()
    {
        track.traj.outcome = Outcome::SkippedGate(g);
    }
    if track.traj.outcome != Outcome::Arrived {
        return;
    }
    let Some(acc) = scn.particles[i].acceptance else {
        return;
    };
    let p = track.traj.end.p;
    let mut m = acc.margin(kin.velocity(p), kin.kinetic_energy(p));
    if let Some(w) = acc.radiation {
        let sources: Vec<(f64, &[Emission])> = scn
            .particles
            .iter()
            .zip(emission)
            .map(|(b, e)| (b.particle.charge, e.as_slice()))
            .collect();
        let e = w.measure_system(&sources, scn.c);
        track.traj.radiation = Some(e);
        m = m.min(w.margin(e));
    }
    if m < 0.0 {
        track.traj.outcome = Outcome::Rejected;
    }
    track.traj.margins = Some(Margins {
        obstacles: Vec::new(),
        bounds: None,
        detector: None,
        acceptance: Some(m),
        gates: Vec::new(),
        gate_acceptance: Vec::new(),
    });
}

/// A charge present for the energy budget: charge, moment, position, kinetic energy.
type Present = (f64, f64, DVec3, f64);

/// Kinetic, potential and interaction energy of the charges present (`EnergySample`).
fn budget<F: FieldSolver>(scn: &BeamScenario<F>, t: f64, present: &[Present]) -> (f64, f64, f64) {
    let (mut kinetic, mut potential, mut interaction) = (0.0, 0.0, 0.0);
    for (k, &(q, m, x, tk)) in present.iter().enumerate() {
        let f = scn.field.sample(x, t);
        kinetic += tk;
        potential += q * f.phi + 0.5 * q * scn.field.self_field(x, q).1 - m * f.b.z;
        if scn.interact {
            for &(q2, _, x2, _) in &present[k + 1..] {
                interaction += q * q2 / (x - x2).length();
            }
        }
    }
    (kinetic, potential, interaction)
}

/// Depth to which ghosts are followed into the boundary they crossed, and at which their
/// penetration margin is capped. A ghost shares the step size of the whole beam: followed
/// deep into a point charge (as far as `MARGIN_SAFE`, 0.05 cells from its centre for a
/// 0.3-cell charge) its huge force forced tiny steps on every particle, so that whole
/// beams failed (`StepSizeTooSmall`). A depth reached is recorded as at most this much,
/// the same in every run, so it verifies (it is far above the numerical error).
pub const GHOST_DEPTH: f64 = 0.01;

fn finish_ghost(t: &mut Track) {
    let mut real = false;
    if let Phase::Ghost {
        event,
        depth,
        real: r,
    } = t.phase
    {
        t.margin[event] = t.margin[event].min(depth.max(-GHOST_DEPTH));
        real = r;
    }
    t.phase = if real { Phase::Free } else { Phase::Done };
}

/// Fields `(E, B)` at `x` of a charge `q` moving uniformly with velocity `v` that is now
/// at `r` (the boosted Coulomb field; Jackson §11.10):
/// `E = q (1 − β²) R / (R³ (1 − β² + (β·R̂)²)^{3/2})`, `B = v × E / c²`, `R = x − r`.
pub fn heaviside_fields(q: f64, c: f64, x: DVec3, r: DVec3, v: DVec3) -> (DVec3, DVec3) {
    let d = x - r;
    let dist = d.length();
    let beta = v / c;
    let b2 = beta.length_squared();
    let rb = beta.dot(d) / dist;
    let s = 1.0 - b2 + rb * rb;
    let e = d * (q * (1.0 - b2) / (dist * dist * dist * s * s.sqrt()));
    (e, v.cross(e) / (c * c))
}

/// Acceleration of every particle at launch: external fields and the others' fields of
/// uniform motion (the quasi-static interaction) at t = 0.
fn launch_accelerations<F: FieldSolver>(scn: &BeamScenario<F>) -> Vec<DVec3> {
    let kin: Vec<Kinematics> = scn
        .particles
        .iter()
        .map(|b| Kinematics::new(b.particle.mass, scn.c))
        .collect();
    (0..scn.particles.len())
        .map(|i| {
            let b = &scn.particles[i];
            let q = b.particle.charge;
            let v = kin[i].velocity(b.p0);
            let f = scn.field.sample(b.x0, 0.0);
            let (mut e, mut bf) = (f.e, f.b);
            for (j, s) in scn.particles.iter().enumerate() {
                if j != i && s.particle.charge != 0.0 {
                    let vj = kin[j].velocity(s.p0);
                    let (ej, bj) = heaviside_fields(s.particle.charge, scn.c, b.x0, s.x0, vj);
                    e += ej;
                    bf += bj;
                }
            }
            let mut force = (e + v.cross(bf)) * q;
            if b.particle.moment != 0.0 {
                force += scn.field.grad_bz(b.x0, 0.0) * b.particle.moment;
            }
            kin[i].acceleration(b.p0, force)
        })
        .collect()
}

/// How far past the end of the record (in lengths of its last step) the last step's
/// polynomial is used: retarded times inside the current step (pairs closer than `c h`)
/// come from this extrapolation, the standard treatment of delays shorter than a step.
const EXTRAPOLATION: f64 = 2.0;

/// The quasi-static interaction's fields `(E, B)` at `x` of a source `q` (mass `m`) now at
/// `r` with velocity `v` and acceleration `a`, which feels the fields `e`, `b` (the
/// external ones and the others'): its past continued along its motion in those fields,
/// blending into the tapered constant acceleration where that swings far (`Continuation`),
/// or with the acceleration alone out of the plane. For the game's field views, which
/// show what this interaction leaves out.
#[allow(clippy::too_many_arguments)]
pub fn continued_fields(
    q: f64,
    m: f64,
    c: f64,
    x: DVec3,
    r: DVec3,
    v: DVec3,
    a: DVec3,
    fields: (DVec3, DVec3),
) -> (DVec3, DVec3) {
    let cont = FieldMotion::new(q, m, c, 0.0, r, v, fields.0, fields.1)
        .map_or(Continuation::Accelerated(a), Continuation::Fields);
    let (e, b, _) = cont.fields(q, c, x, 0.0, r, v);
    (e, b)
}

/// The weight of the motion in the fields in the quasi-static interaction's field at the
/// excursion `e` (`SWING`), for the game's GPU field view.
pub fn continuation_weight(e: f64) -> f64 {
    swing_weight(e)
}

/// Where the continued past (`tapered`) stops having exactly the constant acceleration: at
/// `|a τ| = TAPER_START · Δv` (`taper_reach`); beyond, the acceleration fades as `sech²`
/// over the rest.
const TAPER_START: f64 = 0.7;

/// The largest change of velocity of the continued past (`tapered`) of a charge moving at
/// `v`: `Δv = min(0.1 c, 0.9 (c − |v|))`. The velocity changes along the fixed direction
/// of `a`, so its speed stays below `|v| + Δv ≤ 0.9 c + 0.1 |v| < c`. (A cap of 0.1c
/// alone let a past braked from `|v| > 0.9 c` exceed c; every shipped beam is slower than
/// 0.89 c, where the cap is the 0.1c.)
pub fn taper_reach(v: DVec3, c: f64) -> f64 {
    (0.1 * c).min(0.9 * (c - v.length()).max(0.0))
}

/// Motion continued from position `r`, velocity `v` and acceleration `a` by `τ` (the
/// quasi-static interaction's and the pre-launch past; the game's field views continue
/// the world lines before launch with it too). With `T = Δv/|a|` (`Δv = taper_reach`)
/// and `u = τ/T`: the constant acceleration while `|u| ≤ u₁ = TAPER_START`; beyond, the
/// acceleration `a sech²(s/w)`, `s = |u| − u₁`, `w = 1 − u₁`, so that the velocity
/// changes by at most `Δv`, never reaching c, and the acceleration and its derivative stay
/// continuous. The first versions switched to uniform motion at `|u| = 1`: a jump of the
/// acceleration, whose field arriving at the neighbours made the steps collapse to 1e-8
/// for a beam launched beside strong charges (the preview took 50 s, the exact
/// verification did not finish). A taper `a sech²(u)` from the start cured that but
/// departed from the constant acceleration early and made the relativistic beam's
/// quasi-static flight 1e-3 cells worse; with `u₁ = 0.7` its difference from the exact
/// one is 1.79e-3 cells (kinked: 1.72e-3).
pub fn tapered(r: DVec3, v: DVec3, a: DVec3, tau: f64, c: f64) -> (DVec3, DVec3, DVec3) {
    let t_scale = taper_reach(v, c).max(1e-300) / a.length().max(1e-300);
    let u1 = TAPER_START;
    let w = 1.0 - u1;
    let u = tau / t_scale;
    if u.abs() <= u1 {
        return (r + v * tau + a * (0.5 * tau * tau), v + a * tau, a);
    }
    let s = (u.abs() - u1) / w;
    let ln_cosh = s + libm::log1p(libm::exp(-2.0 * s)) - std::f64::consts::LN_2;
    let th = libm::tanh(s);
    let sign = u.signum();
    (
        r + v * tau + a * (t_scale * t_scale * (0.5 * u1 * u1 + u1 * w * s + w * w * ln_cosh)),
        v + a * (t_scale * sign * (u1 + w * th)),
        a * (1.0 - th * th),
    )
}

/// Liénard–Wiechert fields `(E, B)` at `x`, a time `dt` after the present, of a charge
/// `q` now at `r` with velocity `v` and acceleration `a`, whose past is taken as
/// `r + v τ + a τ²/2` (velocity `v + a τ`, acceleration `a`) while `|a τ| ≤ 0.7 Δv`, its
/// acceleration fading smoothly before that (`tapered`: the velocity changes by at most
/// `Δv = taper_reach`, so the past never becomes superluminal; far points in a strong
/// field see the source's older past, which constant acceleration would not describe
/// anyway). The retarded time is the root of `g(τ) = c (dt − τ) − |x − r(τ)|`, strictly
/// decreasing on this subluminal curve, with `g(dt) ≤ 0`: Newton's method, kept inside
/// the bracket of the signs seen so far (bisection where it would leave it). The error is
/// of the order of the jerk, `|ȧ| τ³` in position with τ = R/c.
pub fn accelerated_fields(
    q: f64,
    c: f64,
    x: DVec3,
    dt: f64,
    r: DVec3,
    v: DVec3,
    a: DVec3,
) -> (DVec3, DVec3) {
    let at = |tau: f64| {
        let (rc, vc, _) = tapered(r, v, a, tau, c);
        (rc, vc)
    };
    // From the guess of a source at rest.
    let mut tau = dt - (x - r).length() / c;
    let (mut lo, mut hi) = (f64::NEG_INFINITY, dt);
    for _ in 0..100 {
        let (rp, vp) = at(tau);
        let d = x - rp;
        let dist = d.length();
        let g = c * (dt - tau) - dist;
        if g > 0.0 {
            lo = lo.max(tau);
        } else {
            hi = hi.min(tau);
        }
        let slope = -c + d.dot(vp) / dist;
        let mut next = tau - g / slope;
        if !(next < hi && (next > lo || lo == f64::NEG_INFINITY)) || !next.is_finite() {
            next = if lo.is_finite() {
                0.5 * (lo + hi)
            } else {
                hi - 2.0 * (hi - tau).max(dist / c)
            };
        }
        let done = (next - tau).abs() <= 4.0 * f64::EPSILON * (dt.abs() + dist / c).max(1e-300);
        tau = next;
        if done {
            break;
        }
    }
    let (rp, vp, ap) = tapered(r, v, a, tau, c);
    let f = fields_from(q, c, x, tau, rp, vp, ap);
    (f.e(), f.b)
}

/// Liénard power `(2/3) q² γ⁶ (a² − |v×a|²/c²) / c³` of a particle with momentum `p` under
/// the force `dp` (units `k = 1`).
fn larmor_power(kin: &Kinematics, q: f64, p: DVec3, dp: DVec3) -> f64 {
    let c = kin.c;
    let gamma = kin.gamma(p);
    let v = kin.velocity(p);
    let a = kin.acceleration(p, dp);
    (2.0 / 3.0)
        * q
        * q
        * gamma.powi(6)
        * (a.length_squared() - v.cross(a).length_squared() / (c * c))
        / (c * c * c)
}

fn add_stats(a: Stats, b: Stats) -> Stats {
    Stats {
        n_fcn: a.n_fcn + b.n_fcn,
        n_step: a.n_step + b.n_step,
        n_accept: a.n_accept + b.n_accept,
        n_reject: a.n_reject + b.n_reject,
    }
}
