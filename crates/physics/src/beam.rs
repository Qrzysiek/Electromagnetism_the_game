//! Beams: many particles in one ODE system, interacting with each other (PHYSICS.md §3.3).
//!
//! All particles share one state vector and one step size. The interaction is the exact
//! pairwise Coulomb force for `c = ∞`, where it is the whole interaction: the magnetic
//! and retarded parts, and every interaction of magnetic moments, scale as `1/c²` and
//! vanish. For finite `c` the interaction is not modelled here (`BeamScenario::interact`
//! must be off, or the level must be Newtonian).
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

use glam::DVec3;

use crate::dynamics::{Kinematics, Particle};
use crate::events::first_crossing;
use crate::field::FieldSolver;
use crate::geometry::{Aabb, Region, Shape};
use crate::integrator::OdeSystem;
use crate::integrator::dop853::{Dense, Dop853, Settings, Stats};
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
    /// Coulomb interaction between the particles (requires `c = ∞`).
    pub interact: bool,
    /// Gates every particle must pass, in order, before its detector counts (as for single
    /// flights, PHYSICS.md §6.2).
    pub gates: Vec<Gate>,
}

/// Result of a beam flight: one trajectory per particle (in launch order; their `stats`
/// are those of the whole system), and diagnostics of the whole system.
#[derive(Clone, Debug)]
pub struct BeamRun {
    pub trajectories: Vec<Trajectory>,
    pub stats: Stats,
    /// Largest relative change of the total energy of the particles still flying
    /// (kinetic, external potential, moments, pair interaction), between removals.
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
        if self.scn.interact && part.charge != 0.0 {
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
        force
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
            if self.scn.interact {
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
    assert!(
        !scn.interact || scn.c.is_infinite(),
        "beam interaction is only modelled for c = ∞"
    );
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
            let reached_end = match int.step(&ode, scn.t_max, f64::INFINITY) {
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
        energy_max_rel_error: energy_err,
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

fn add_stats(a: Stats, b: Stats) -> Stats {
    Stats {
        n_fcn: a.n_fcn + b.n_fcn,
        n_step: a.n_step + b.n_step,
        n_accept: a.n_accept + b.n_accept,
        n_reject: a.n_reject + b.n_reject,
    }
}
