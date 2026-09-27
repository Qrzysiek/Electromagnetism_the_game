//! Beams: many particles in one ODE system, interacting with each other (PHYSICS.md §3.3).
//!
//! All particles share one state vector and one step size. The interaction is exact
//! classical electrodynamics between point charges:
//!
//! - For `c = ∞` it is the pairwise Coulomb force, which is then the whole interaction:
//!   the magnetic and retarded parts, and every interaction of magnetic moments, scale as
//!   `1/c²` and vanish.
//! - For finite `c`, by default, the quasi-static interaction: each particle feels the
//!   exact Liénard–Wiechert fields of the others, with their recent past taken from
//!   their present state with constant acceleration (`accelerated_fields`; the
//!   accelerations from a first pass with the fields of uniform motion,
//!   `heaviside_fields`). That is exact in the velocities (so it contains the magnetic
//!   attraction that reduces the space charge of a relativistic beam by 1/γ²) and to
//!   first order in the accelerations (their near and radiation fields); it leaves out
//!   the jerk terms, estimated per particle (`BeamRun::neglected_retardation`). It needs
//!   no record and costs a few times the Coulomb force.
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
//! Each particle has the events and gates of a single flight (obstacles, bounds, detector,
//! gates in order). The earliest event of any particle ends the step there; that
//! particle's outcome is recorded, and the integration restarts at that time without it
//! as a source. It is still integrated for a while as a *ghost* (pushed by the others,
//! pushing none) until its penetration depth into the boundary it crossed is known,
//! exactly as the single-particle runner follows the continued trajectory
//! (`trajectory.rs`): margins are then properties of the trajectories, not of where the
//! steps ended, and the existing verification (`verify::classify`) applies to every
//! particle.

use std::cell::RefCell;

use glam::DVec3;

use crate::dynamics::{Kinematics, Particle, landau_lifshitz};
use crate::events::first_crossing;
use crate::field::FieldSolver;
use crate::geometry::{Aabb, Region, Shape};
use crate::integrator::OdeSystem;
use crate::integrator::dop853::{Dense, Dop853, Settings, Stats};
use crate::lienard::fields_from;
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
    /// pairs so far apart that the source's past is taken as uniform), over that of the
    /// fields kept, `Σ_j |E_j|`
    /// (sampled at the step ends; 0 without quasi-static interaction). The trajectory
    /// error follows the impulse error: a short plunge of a neighbour matters little.
    pub neglected_retardation: Vec<f64>,
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
}

/// Recorded motion of the beam for the retarded interaction: the dense output of every
/// accepted step (valid up to `t_stop`, where an event may have cut it short), with the
/// particles it contains, and for every particle the time its charge was removed.
struct Past {
    /// Recorded steps, oldest first; steps that no retarded time can reach any more are
    /// dropped (`prune`).
    segments: RefCell<Vec<Segment>>,
    /// Removal time and position of each particle (infinite time while it flies).
    t_off: RefCell<Vec<f64>>,
    x_off: RefCell<Vec<DVec3>>,
    /// Whether a removed particle stayed where it stopped (at rest from then on) rather
    /// than being drained.
    stays: RefCell<Vec<bool>>,
    kin: Vec<Kinematics>,
    p_ref: f64,
    /// Acceleration of each particle at launch (external fields and the others'
    /// quasi-static fields): the motion before launch continues it.
    a0: Vec<DVec3>,
}

struct Segment {
    dense: Dense,
    t_start: f64,
    t_stop: f64,
    /// Particle indices (sorted): member `k` is the state's block `k`.
    members: Vec<usize>,
}

impl Past {
    /// Before launch each particle moves with its launch acceleration (it was already
    /// flying in the fields): `x0 + v0 t + a0 t²/2`, `v0 + a0 t`. A sudden start from
    /// uniform motion would send a kink in every acceleration field to every other
    /// particle, each forcing tiny steps when it arrives. Continued uniformly where
    /// |a0 t| would exceed 0.1c (reached only while bracketing).
    fn before_launch<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3, DVec3) {
        let b = &scn.particles[j];
        let v0 = self.kin[j].velocity(b.p0);
        let a0 = self.a0[j];
        let t_lim = 0.1 * scn.c / a0.length().max(1e-300);
        let tt = t.max(-t_lim);
        let x = b.x0 + v0 * tt + a0 * (0.5 * tt * tt);
        let v = v0 + a0 * tt;
        if t < tt {
            (x + v * (t - tt), v, DVec3::ZERO)
        } else {
            (x, v, a0)
        }
    }

    /// Position, velocity and acceleration of particle `j` at time `t`: `before_launch`,
    /// then the recorded motion, extrapolated a little past its end (see `EXTRAPOLATION`).
    fn state<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3, DVec3) {
        // At rest where it stopped.
        if self.stays.borrow()[j] && t >= self.t_off.borrow()[j] {
            return (self.x_off.borrow()[j], DVec3::ZERO, DVec3::ZERO);
        }
        let kin = &self.kin[j];
        // Last segment starting before t that contains j.
        let segments = self.segments.borrow();
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
                let v = kin.velocity(DVec3::new(comp(3), comp(4), comp(5)) * self.p_ref);
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
        let p = DVec3::new(comp(3), comp(4), comp(5)) * self.p_ref;
        let v = kin.velocity(p);
        if t > reach {
            return (x + v * (t - reach), v, DVec3::ZERO);
        }
        let dp = DVec3::new(
            d.eval_derivative_component(6 * k + 3, tt),
            d.eval_derivative_component(6 * k + 4, tt),
            d.eval_derivative_component(6 * k + 5, tt),
        ) * self.p_ref;
        (x, v, kin.acceleration(p, dp))
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
        for (j, &x) in self.x_off.borrow().iter().enumerate() {
            if !members.contains(&j) {
                (lo, hi) = (lo.min(x), hi.max(x));
            }
        }
        let extent = (hi - lo).length();
        let cutoff = t - 2.0 * extent / (c * (1.0 - beta_max.min(0.999_999))) - 1e-3;
        let mut segments = self.segments.borrow_mut();
        let old = segments.partition_point(|s| s.t_stop < cutoff);
        // Drop in batches (and always keep the newest step).
        if old >= 256 && old < segments.len() {
            segments.drain(..old);
        }
    }

    /// Retarded fields `(E, B)` of particle `j` at `x` and time `t`; zero once the field
    /// of its removal has arrived.
    fn fields_of<F>(&self, scn: &BeamScenario<F>, j: usize, x: DVec3, t: f64) -> (DVec3, DVec3) {
        let c = scn.c;
        let qj = scn.particles[j].particle.charge;
        let g = |tr: f64| c * (t - tr) - (x - self.state(scn, j, tr).0).length();
        // The retarded time is the root of g, which decreases strictly (|v| < c). The
        // field of a removed particle is gone once the light cone has passed its removal
        // (decided from the exact removal point).
        let t_off = self.t_off.borrow()[j];
        let stays = self.stays.borrow()[j];
        if !stays && t_off.is_finite() && c * (t - t_off) >= (x - self.x_off.borrow()[j]).length() {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        let mut hi = if stays { t } else { t.min(t_off) };
        let mut g_hi = g(hi);
        if g_hi >= 0.0 {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        let mut step = -g_hi / c;
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
            let (r, vr, _) = self.state(scn, j, tr);
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
            let tol = 4.0 * f64::EPSILON * t.abs().max(1.0);
            let done = (next - tr).abs() <= tol || hi - lo <= tol;
            tr = next;
            if done {
                break;
            }
        }
        let (r, vr, ar) = self.state(scn, j, tr);
        let f = fields_from(qj, c, x, tr, r, vr, ar);
        (f.e(), f.b)
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
        let acc = self.accelerations(y, t);
        self.force_with(y, k, t, acc.as_deref(), true)
    }

    /// Quasi-static interaction: the accelerations of all members in the state `y`, from
    /// the external fields and the others' fields of uniform motion (first pass).
    fn accelerations(&self, y: &[f64], t: f64) -> Option<Vec<DVec3>> {
        self.quasi_static().then(|| {
            (0..self.members.len())
                .map(|j| {
                    let f = self.force_with(y, j, t, None, false);
                    self.kin[j].acceleration(self.p(y, j), f)
                })
                .collect()
        })
    }

    /// Force on member `k`; for the quasi-static interaction `acc` are the members'
    /// accelerations (None: fields of uniform motion), `rr` adds radiation reaction.
    fn force_with(&self, y: &[f64], k: usize, t: f64, acc: Option<&[DVec3]>, rr: bool) -> DVec3 {
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
                let (e, b) = self.retarded_fields(k, x, t);
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
            force += self.radiation_reaction_force(y, k, t);
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
    /// state `y` of time `t_y` (each continued uniformly to `t`).
    fn quasi_static_fields(
        &self,
        y: &[f64],
        t_y: f64,
        k: usize,
        x: DVec3,
        t: f64,
        acc: Option<&[DVec3]>,
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
            let (ej, bj) = if let Some(acc) = acc {
                accelerated_fields(
                    qj,
                    self.scn.c,
                    x,
                    t - t_y,
                    Self::x(y, j),
                    self.kin[j].velocity(self.p(y, j)),
                    acc[j],
                )
            } else {
                let vj = self.kin[j].velocity(self.p(y, j));
                let rj = Self::x(y, j) + vj * (t - t_y);
                heaviside_fields(qj, self.scn.c, x, rj, vj)
            };
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

    /// Sum of the retarded fields of the other particles at `x`, `t` (member `k`).
    fn retarded_fields(&self, k: usize, x: DVec3, t: f64) -> (DVec3, DVec3) {
        let (mut e, mut b) = (DVec3::ZERO, DVec3::ZERO);
        if let Some(past) = self.past {
            let i = self.members[k];
            for j in 0..self.scn.particles.len() {
                if j != i && self.scn.particles[j].particle.charge != 0.0 {
                    let (ej, bj) = past.fields_of(self.scn, j, x, t);
                    e += ej;
                    b += bj;
                }
            }
        }
        (e, b)
    }

    /// Landau–Lifshitz force on member `k` in the total field: external plus the other
    /// particles' retarded fields.
    fn radiation_reaction_force(&self, y: &[f64], k: usize, t: f64) -> DVec3 {
        let q = self.scn.particles[self.members[k]].particle.charge;
        let t_y = t;
        let fields = |x: DVec3, t: f64| {
            let f = self.scn.field.sample(x, t);
            let (e, b) = if self.quasi_static() {
                self.quasi_static_fields(y, t_y, k, x, t, None)
            } else {
                self.retarded_fields(k, x, t)
            };
            (f.e + e, f.b + b)
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
        let acc = self.accelerations(y, t);
        for k in 0..self.members.len() {
            let v = self.kin[k].velocity(self.p(y, k));
            let dp = self.force_with(y, k, t, acc.as_deref(), true) / self.p_ref;
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
                },
                phase: Phase::Flying,
                margin: g.clone(),
                g_prev: g,
                gates: GateTracker::new(&scn.gates, b.x0),
            }
        })
        .collect();
    // Particles launched inside a boundary end at once.
    for t in &mut tracks {
        if let Some(k) = t.g_prev.iter().position(|&g| g <= 0.0) {
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
    // Energy budget (recorded runs): samples, and the energy absorbed so far.
    let mut energy: Vec<EnergySample> = Vec::new();
    let mut absorbed = 0.0;
    let past = Past {
        segments: RefCell::new(Vec::new()),
        stays: RefCell::new(vec![false; n]),
        t_off: RefCell::new(vec![f64::INFINITY; n]),
        x_off: RefCell::new(scn.particles.iter().map(|b| b.x0).collect()),
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
    };
    // Particles that end at launch: drained, stopped there, or flying on.
    for (i, t) in tracks.iter_mut().enumerate() {
        if matches!(t.phase, Phase::Done) {
            let k = t
                .g_prev
                .iter()
                .position(|&g| g <= 0.0)
                .expect("ended at launch");
            match events[k].fate(&scn.fates) {
                Fate::Drain => past.t_off.borrow_mut()[i] = 0.0,
                Fate::Stop => {
                    past.t_off.borrow_mut()[i] = 0.0;
                    past.stays.borrow_mut()[i] = true;
                    stopped[i] = Some(scn.particles[i].x0);
                }
                Fate::Pass => t.phase = Phase::Free,
            }
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
        };
        let mut y0 = Vec::with_capacity(6 * members.len());
        for &i in &members {
            let (x, p) = states[i];
            let s = p / p_ref;
            y0.extend_from_slice(&[x.x, x.y, x.z, s.x, s.y, s.z]);
        }
        let energy0 = ode.energy(&y0, t_now);
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
            // Retarded interaction: every step shorter than a third of the light travel
            // time between the closest pair of a charged member and a flying source. In
            // a step of length h the distance shrinks by less than 2ch, so it stays above
            // ch, and every retarded time of the step lies before its start: in the record.
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
                    .borrow()
                    .last()
                    .map_or(0.0, |s| s.t_stop - s.t_start);
                let w = (2.0 * v_close).min(2.0 * scn.c);
                (EXTRAPOLATION * previous + r_min / scn.c) / (1.0 + w / scn.c)
            } else {
                f64::INFINITY
            };
            let reached_end = match int.step(&ode, scn.t_max, h_cap) {
                Ok(done) => done,
                Err(e) => {
                    failed = Some(e);
                    break 'segments;
                }
            };
            let t_b = int.t();
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

            // Event values at the step end, and the earliest event of a flying particle.
            let mut g_next: Vec<Vec<f64>> = vec![Vec::new(); members.len()];
            let mut v_max = vec![0.0; members.len()];
            let mut first: Option<(f64, usize, usize)> = None; // (t, member, event)
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
                let x_b = BeamOde::<F>::x(int.y(), k);
                g_next[k] = events.iter().map(|&ev| event_value(ev, i, x_b)).collect();
                for (e, &ev) in events.iter().enumerate() {
                    let mut g = |t: f64| event_value(ev, i, pos(k, t));
                    let r = first_crossing(
                        &mut g,
                        t_a,
                        t_b,
                        tracks[i].g_prev[e],
                        g_next[k][e],
                        v_max[k],
                    );
                    if let Some(t) = r.time
                        && first.is_none_or(|(tf, _, _)| t < tf)
                    {
                        first = Some((t, k, e));
                    }
                }
            }
            let t_end = first.map_or(t_b, |(t, _, _)| t);
            if retarded {
                past.segments.borrow_mut().push(Segment {
                    dense: dense.clone(),
                    t_start: t_a,
                    t_stop: t_end,
                    members: members.clone(),
                });
            }
            // Quasi-static interaction: estimated neglected fields (the sources' jerk: their
            // past is taken with constant acceleration), relative to the fields kept, at the
            // step's end. The jerk is the change of the acceleration over the step.
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
                            missing +=
                                qj.abs() * g * g * jerk[l].length() / (scn.c.powi(3) * kappa);
                            // Beyond the constant-acceleration range (|a τ| > 0.1 c) the
                            // past is continued uniformly: first-order error there.
                            let a_now = accel(l, t_end).length();
                            if a_now * delay > 0.1 * scn.c {
                                missing += qj.abs() * g * g * a_now / (scn.c * scn.c * r);
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
            // Largest ratio of the radiation-reaction force to the rest (the Landau–Lifshitz
            // treatment needs it small), at the step's end.
            if scn.radiation_reaction && scn.c.is_finite() {
                let mut y_end = vec![0.0; ode.dim()];
                dense.eval(t_end, &mut y_end);
                for (k, &i) in members.iter().enumerate() {
                    if matches!(tracks[i].phase, Phase::Flying) && ode.reacts(k) {
                        let rr = ode.radiation_reaction_force(&y_end, k, t_end);
                        let rest = (ode.force(&y_end, k, t_end) - rr).length();
                        if rest > 0.0 {
                            let t = &mut tracks[i].traj;
                            t.reaction_ratio_max = t.reaction_ratio_max.max(rr.length() / rest);
                        }
                    }
                }
            }
            // Radiated energy of the flying particles (Liénard power, trapezoidal rule),
            // from the force given by the derivative of the dense output.
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
                        tracks[i].traj.radiated_energy +=
                            0.5 * (power(t_a) + power(t_end)) * (t_end - t_a);
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
            energy_err = energy_err.max((e - energy0).abs() / energy0.abs().max(1e-300));

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
                    Fate::Drain => {
                        past.t_off.borrow_mut()[i] = t_end;
                        past.x_off.borrow_mut()[i] = tracks[i].traj.end.x;
                    }
                    Fate::Stop => {
                        past.t_off.borrow_mut()[i] = t_end;
                        past.x_off.borrow_mut()[i] = tracks[i].traj.end.x;
                        past.stays.borrow_mut()[i] = true;
                        stopped[i] = Some(tracks[i].traj.end.x);
                    }
                    Fate::Pass => {}
                }
                if tracks[i].traj.outcome == Outcome::Arrived
                    && let Some(g) = tracks[i].gates.missing()
                {
                    tracks[i].traj.outcome = Outcome::SkippedGate(g);
                }
                // Acceptance, decided at the moment of entry.
                if tracks[i].traj.outcome == Outcome::Arrived
                    && let Some(acc) = scn.particles[i].acceptance
                {
                    let p = tracks[i].traj.end.p;
                    let m = acc.margin(ode.kin[k].velocity(p), ode.kin[k].kinetic_energy(p));
                    if m < 0.0 {
                        tracks[i].traj.outcome = Outcome::Rejected;
                    }
                    tracks[i].traj.margins = Some(Margins {
                        obstacles: Vec::new(),
                        bounds: None,
                        detector: None,
                        acceptance: Some(m),
                        gates: Vec::new(),
                        gate_acceptance: Vec::new(),
                    });
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
    let trajectories = tracks
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
    })
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

/// Liénard–Wiechert fields `(E, B)` at `x`, a time `dt` after the present, of a charge
/// `q` now at `r` with velocity `v` and acceleration `a`, whose past is taken as
/// `r + v τ + a τ²/2` (velocity `v + a τ`, acceleration `a`) while `|a τ| ≤ 0.1 c`, and
/// uniform motion before that (so that the past never becomes superluminal: far points
/// in a strong field see the source's older past, which constant acceleration would not
/// describe anyway). The retarded time is solved on this curve by Newton's method. The
/// error is of the order of the jerk, `|ȧ| τ³` in position with τ = R/c.
pub fn accelerated_fields(
    q: f64,
    c: f64,
    x: DVec3,
    dt: f64,
    r: DVec3,
    v: DVec3,
    a: DVec3,
) -> (DVec3, DVec3) {
    let tau_lim = 0.1 * c / a.length().max(1e-300);
    let at = |tau: f64| {
        let tc = tau.max(-tau_lim);
        let (rc, vc) = (r + v * tc + a * (0.5 * tc * tc), v + a * tc);
        (rc + vc * (tau - tc), vc)
    };
    // Retarded τ: the root of g(τ) = c (dt − τ) − |x − r(τ)|, from the uniform-motion guess.
    let mut tau = dt - (x - r).length() / c;
    for _ in 0..50 {
        let (rp, vp) = at(tau);
        let d = x - rp;
        let dist = d.length();
        let g = c * (dt - tau) - dist;
        let slope = -c + d.dot(vp) / dist;
        let next = tau - g / slope;
        let done = (next - tau).abs() <= 4.0 * f64::EPSILON * (dt.abs() + dist / c).max(1e-300);
        tau = next;
        if done {
            break;
        }
    }
    let (rp, vp) = at(tau);
    let ap = if tau >= -tau_lim { a } else { DVec3::ZERO };
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
