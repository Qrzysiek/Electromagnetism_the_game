//! Beams: many particles in one ODE system, interacting with each other (PHYSICS.md §3.3).
//!
//! All particles share one state vector and one step size. The interaction is exact
//! classical electrodynamics between point charges:
//!
//! - For `c = ∞` it is the pairwise Coulomb force, which is then the whole interaction:
//!   the magnetic and retarded parts, and every interaction of magnetic moments, scale as
//!   `1/c²` and vanish.
//! - For finite `c` each particle feels the Liénard–Wiechert fields (`lienard.rs`) of the
//!   others at their retarded times, taken from the recorded motion (the dense output of
//!   every accepted step). Every step is shorter than a third of the light travel time
//!   between the closest pair, so that every retarded time falls into the recorded past
//!   and the system stays an explicit ODE. Before launch the particles move uniformly
//!   (with their launch velocity). A removed particle acts for as long as its field from
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

/// A beam flight: the shared scene and the particles, all launched at t = 0.
#[derive(Clone, Debug)]
pub struct BeamScenario<F> {
    pub field: F,
    pub obstacles: Vec<Shape>,
    pub particles: Vec<BeamParticle>,
    pub c: f64,
    pub bounds: Option<Aabb>,
    pub t_max: f64,
    /// Interaction between the particles: Coulomb for `c = ∞`, retarded Liénard–Wiechert
    /// fields for finite `c`.
    pub interact: bool,
    /// Gates every particle must pass, in order, before its detector counts (as for single
    /// flights, PHYSICS.md §6.2).
    pub gates: Vec<Gate>,
    /// Each particle's radiation reaction (Landau–Lifshitz, PHYSICS.md §3.1) in the total
    /// field; no effect for `c = ∞`.
    pub radiation_reaction: bool,
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
}

/// Event functions of one particle, in priority order.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Event {
    Obstacle(usize),
    Bounds,
    Detector,
}

impl Event {
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
    kin: Vec<Kinematics>,
    p_ref: f64,
}

struct Segment {
    dense: Dense,
    t_start: f64,
    t_stop: f64,
    /// Particle indices (sorted): member `k` is the state's block `k`.
    members: Vec<usize>,
}

impl Past {
    /// Position, velocity and acceleration of particle `j` at time `t`: uniform motion
    /// before launch, the recorded motion after, continued uniformly beyond what is
    /// recorded (reached only while bracketing a retarded time).
    fn state<F>(&self, scn: &BeamScenario<F>, j: usize, t: f64) -> (DVec3, DVec3, DVec3) {
        let b = &scn.particles[j];
        let kin = &self.kin[j];
        // Last segment starting before t that contains j.
        let segments = self.segments.borrow();
        let upto = segments.partition_point(|s| s.t_start < t);
        let found = segments[..upto]
            .iter()
            .rev()
            .find_map(|s| s.members.binary_search(&j).ok().map(|k| (s, k)));
        let Some((seg, k)) = found else {
            // Before the record: uniform motion before launch, or, where old steps were
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
            let v0 = kin.velocity(b.p0);
            return (b.x0 + v0 * t, v0, DVec3::ZERO);
        };
        let d = &seg.dense;
        // Past the end of the record (only while bracketing a retarded time, or for the
        // tiny offsets of the radiation reaction's field derivative): the step's own
        // polynomial for up to one step length, uniform motion beyond.
        let reach = seg.t_stop + (seg.t_stop - seg.t_start);
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
        if t_off.is_finite() && c * (t - t_off) >= (x - self.x_off.borrow()[j]).length() {
            return (DVec3::ZERO, DVec3::ZERO);
        }
        let mut hi = t.min(t_off);
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
        } else if self.scn.interact && part.charge != 0.0 {
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
        if self.reacts(k) {
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
        let fields = |x: DVec3, t: f64| {
            let f = self.scn.field.sample(x, t);
            let (e, b) = self.retarded_fields(k, x, t);
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
        for k in 0..self.members.len() {
            let v = self.kin[k].velocity(self.p(y, k));
            let dp = self.force(y, k, t) / self.p_ref;
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
    /// when the depth stops decreasing or reaches `GHOST_DEPTH`.
    Ghost {
        event: usize,
        depth: f64,
    },
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
    let retarded = scn.interact && scn.c.is_finite();
    // Largest speed (in units of c) of any particle so far, for pruning the record.
    let mut beta_max: f64 = 0.0;
    let past = Past {
        segments: RefCell::new(Vec::new()),
        t_off: RefCell::new(vec![f64::INFINITY; n]),
        x_off: RefCell::new(scn.particles.iter().map(|b| b.x0).collect()),
        kin: scn
            .particles
            .iter()
            .map(|b| Kinematics::new(b.particle.mass, scn.c))
            .collect(),
        p_ref,
    };
    // Particles that end at launch never act.
    for (i, t) in tracks.iter().enumerate() {
        if matches!(t.phase, Phase::Done) {
            past.t_off.borrow_mut()[i] = 0.0;
        }
    }

    'segments: loop {
        let members: Vec<usize> = (0..n)
            .filter(|&i| !matches!(tracks[i].phase, Phase::Done))
            .collect();
        if members.is_empty() {
            break;
        }
        let source: Vec<bool> = members
            .iter()
            .map(|&i| matches!(tracks[i].phase, Phase::Flying))
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
                let mut r_min = f64::INFINITY;
                for (k, &i) in members.iter().enumerate() {
                    if scn.particles[i].particle.charge == 0.0 {
                        continue;
                    }
                    let xi = BeamOde::<F>::x(int.y(), k);
                    for (l, &j) in members.iter().enumerate() {
                        if j != i && source[l] && scn.particles[j].particle.charge != 0.0 {
                            r_min = r_min.min((xi - BeamOde::<F>::x(int.y(), l)).length());
                        }
                    }
                }
                r_min / (3.0 * scn.c)
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
                                tracks[i].phase = Phase::Ghost { event: e, depth: m };
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
                    Phase::Ghost { event, depth, .. } => {
                        let g = |t: f64| event_value(events[event], i, pos(k, t));
                        let (t_min, m) = minimize_on(&g, t_a, t_b);
                        let depth = depth.min(m);
                        let unfinished = minimum_at_end(t_min, t_a, t_b);
                        tracks[i].phase = Phase::Ghost { event, depth };
                        if !unfinished || depth <= -GHOST_DEPTH {
                            finish_ghost(&mut tracks[i]);
                            ghosts_done = true;
                        }
                    }
                    Phase::Done => {}
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

            if let Some((_, k, e)) = first {
                let i = members[k];
                tracks[i].traj.outcome = events[e].outcome();
                past.t_off.borrow_mut()[i] = t_end;
                past.x_off.borrow_mut()[i] = tracks[i].traj.end.x;
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
                        Phase::Done => {}
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
        energy_max_rel_error: if retarded || (scn.radiation_reaction && scn.c.is_finite()) {
            f64::NAN
        } else {
            energy_err
        },
        restarts,
    })
}

/// Depth to which ghosts are followed into the boundary they crossed, and at which their
/// penetration margin is capped. A ghost shares the step size of the whole beam: followed
/// deep into a point charge (as far as `MARGIN_SAFE`, 0.05 cells from its centre for a
/// 0.3-cell charge) its huge force forced tiny steps on every particle, so that whole
/// beams failed (`StepSizeTooSmall`). A depth reached is recorded as at most this much,
/// the same in every run, so it verifies (it is far above the numerical error).
pub const GHOST_DEPTH: f64 = 0.01;

fn finish_ghost(t: &mut Track) {
    if let Phase::Ghost { event, depth, .. } = t.phase {
        t.margin[event] = t.margin[event].min(depth.max(-GHOST_DEPTH));
    }
    t.phase = Phase::Done;
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
