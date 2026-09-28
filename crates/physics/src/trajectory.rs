//! Integration of a single test-particle trajectory with events, margins and diagnostics
//! (PHYSICS.md §3–§8).

use glam::DVec3;

use crate::dynamics::{Kinematics, Particle, ParticleOde};
use crate::events::first_crossing;
use crate::field::FieldSolver;
use crate::geometry::{Aabb, Region, Shape};
use crate::integrator::{self, Dense, Dop853, OdeSystem, Settings, Stats};
use crate::spectrum::{Emission, PHASE_STEP, RadiationWindow};

/// Margins at or above this value (grid units) are not refined further: they are far too
/// large for numerical error to matter. Refinement starts below twice this value.
pub const MARGIN_SAFE: f64 = 0.25;

/// Everything that defines one flight.
#[derive(Clone, Debug)]
pub struct Scenario<F> {
    pub field: F,
    /// Solid obstacles (charges, magnets, coil wires); touching one loses the particle.
    pub obstacles: Vec<Shape>,
    pub particle: Particle,
    /// Speed of light in internal units (`f64::INFINITY` for Newtonian mechanics).
    pub c: f64,
    pub x0: DVec3,
    pub p0: DVec3,
    /// Detector B: the flight succeeds when the particle centre enters it.
    pub detector: Option<Region>,
    /// World bounds: leaving them loses the particle (game rule).
    pub bounds: Option<Aabb>,
    /// Maximum flight time (game rule).
    pub t_max: f64,
    /// Include the particle's radiation reaction (Landau–Lifshitz, PHYSICS.md §3.1).
    pub radiation_reaction: bool,
    /// Optional conditions on arrival at the detector (direction, kinetic energy). A
    /// particle entering the detector outside them is absorbed but not counted
    /// (`Outcome::Rejected`).
    pub acceptance: Option<Acceptance>,
    /// Gates the particle must pass, in order, before the detector counts (multi-stage
    /// instruments, PHYSICS.md §6.2). Entering the detector with a gate still missing
    /// ends the flight as `Outcome::SkippedGate`.
    pub gates: Vec<Gate>,
}

/// A pass-through region with optional conditions on the particle entering it
/// (PHYSICS.md §6.2). Gates must not overlap each other or contain the launch point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gate {
    pub region: Region,
    pub acceptance: Option<Acceptance>,
}

/// Conditions on the particle when it enters the detector (PHYSICS.md §6.1).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Acceptance {
    /// Allowed direction of motion: unit axis and half-angle (radians).
    pub direction: Option<(DVec3, f64)>,
    /// Allowed kinetic energy `[min, max]`.
    pub kinetic: Option<(f64, f64)>,
    /// Required radiation of the whole flight into an arc of directions (PHYSICS.md
    /// §3.4), measured when the particle arrives.
    pub radiation: Option<RadiationWindow>,
}

impl Acceptance {
    /// Signed margin of a state entering the detector: the smallest of the angular
    /// margin `half-angle − deviation` (radians) and the relative energy margin
    /// `min(T − min, max − T) / max`. Negative: outside the acceptance.
    pub fn margin(&self, velocity: DVec3, kinetic: f64) -> f64 {
        let mut m = f64::INFINITY;
        if let Some((axis, half)) = self.direction {
            let dev = velocity.cross(axis).length().atan2(velocity.dot(axis));
            m = m.min(half - dev);
        }
        if let Some((lo, hi)) = self.kinetic {
            m = m.min((kinetic - lo).min(hi - kinetic) / hi.abs().max(1e-300));
        }
        m
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RunSettings {
    pub rtol: f64,
    pub atol: f64,
    pub max_steps: u64,
    /// Record the state at every accepted step.
    pub record: bool,
    /// Compute accurate margins to all event boundaries (needed for verification).
    pub margins: bool,
}

impl RunSettings {
    pub fn with_tolerance(tol: f64) -> Self {
        Self {
            rtol: tol,
            atol: tol,
            max_steps: 1_000_000,
            record: true,
            margins: true,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outcome {
    Arrived,
    /// Entered the detector outside its acceptance (direction or energy).
    Rejected,
    /// Entered the detector without having passed gate `i` (and the ones after it).
    SkippedGate(usize),
    /// Collision with obstacle `i`.
    Collided(usize),
    LeftBounds,
    Timeout,
    Failed(integrator::Error),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    pub t: f64,
    pub x: DVec3,
    pub p: DVec3,
}

/// Minimum of each event function over the flight (PHYSICS.md §7).
///
/// For a boundary that was not crossed this is the closest approach (positive). For the
/// boundary that ended the flight it is the minimum over the step containing the event,
/// including the part after the event (negative: the penetration depth the trajectory
/// would have reached). Values `≥ MARGIN_SAFE` are lower-bounded but not refined.
#[derive(Clone, Debug, PartialEq)]
pub struct Margins {
    pub obstacles: Vec<f64>,
    pub bounds: Option<f64>,
    pub detector: Option<f64>,
    /// Acceptance margin at the detector (see `Acceptance::margin`), if the flight
    /// entered a detector with conditions.
    pub acceptance: Option<f64>,
    /// Per gate: the smallest signed distance to it after the previous gate was passed
    /// (negative: the depth reached inside; positive: the closest approach). Infinite if
    /// the previous gate was never passed.
    pub gates: Vec<f64>,
    /// Per gate with conditions: the acceptance margin at the last entry.
    pub gate_acceptance: Vec<Option<f64>>,
}

impl Margins {
    /// All margins, obstacles first.
    pub fn all(&self) -> impl Iterator<Item = f64> + '_ {
        self.obstacles
            .iter()
            .copied()
            .chain(self.bounds)
            .chain(self.detector)
            .chain(self.acceptance)
            .chain(self.gates.iter().copied())
            .chain(self.gate_acceptance.iter().flatten().copied())
    }
}

#[derive(Clone, Debug)]
pub struct Trajectory {
    pub outcome: Outcome,
    /// Final state (at the event, or at `t_max`).
    pub end: Sample,
    /// States at accepted steps (if recorded), ending with `end`.
    pub samples: Vec<Sample>,
    pub stats: Stats,
    /// Total energy `W − mc² = (γ−1)mc² + qφ` at the start.
    pub energy_initial: f64,
    /// Kinetic energy at the start.
    pub kinetic_initial: f64,
    /// Largest `|W(t) − W(0)|` over the accepted steps (static fields only).
    pub energy_max_abs_error: f64,
    /// Energy radiated according to the Liénard formula (trapezoidal rule over steps).
    pub radiated_energy: f64,
    /// Work done on the particle by the radiation-reaction force (≤ 0 over a flight
    /// between force-free states; 0 without radiation reaction). Integrated as part of
    /// the ODE state.
    pub radiation_work: f64,
    /// Largest radiation-reaction force over the largest Lorentz force, both over the
    /// accepted steps' ends (PHYSICS.md §3.1): the Landau–Lifshitz treatment requires it to
    /// be small. (A ratio of the two at the same instant would diverge wherever the
    /// Lorentz force passes through zero, e.g. between alternating magnets, where the
    /// reaction force is no problem.)
    pub reaction_ratio_max: f64,
    /// Present when `RunSettings::margins` is set.
    pub margins: Option<Margins>,
    /// Energy per steradian radiated into the acceptance's radiation window (averaged
    /// over its arc), when it has one and the particle arrived.
    pub radiation: Option<f64>,
    /// With a radiation goal: the dense samples `(t, x, v, a)` of the flight it was
    /// measured on (for displaying the spectrum).
    pub emission: Vec<Emission>,
}

/// View of one accepted step, for observers.
pub struct StepView<'a, F> {
    pub ode: &'a ParticleOde<&'a F>,
    pub dense: &'a Dense,
}

impl<F: FieldSolver> StepView<'_, F> {
    pub fn t_start(&self) -> f64 {
        self.dense.t_start()
    }

    pub fn t_end(&self) -> f64 {
        self.dense.t_end()
    }

    /// Interpolated work done by the radiation-reaction force since launch.
    pub fn radiation_work(&self, t: f64) -> f64 {
        if self.ode.radiation_reaction {
            self.ode.radiation_work(&[
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                self.dense.eval_component(6, t),
            ])
        } else {
            0.0
        }
    }

    /// Interpolated position and momentum.
    pub fn state(&self, t: f64) -> (DVec3, DVec3) {
        let mut y = [0.0; 6];
        self.dense.eval(t, &mut y);
        (ParticleOde::<&F>::position(&y), self.ode.momentum(&y))
    }
}

pub fn run<F: FieldSolver>(scn: &Scenario<F>, rs: &RunSettings) -> Trajectory {
    run_observed(scn, rs, |_| {})
}

/// Event functions, in priority order for simultaneous events.
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

fn event_value<F>(scn: &Scenario<F>, ev: Event, x: DVec3) -> f64 {
    match ev {
        Event::Obstacle(i) => scn.obstacles[i].signed_distance(x) - scn.particle.radius,
        Event::Bounds => -scn
            .bounds
            .expect("bounds event without bounds")
            .signed_distance(x),
        Event::Detector => scn
            .detector
            .expect("detector event without detector")
            .signed_distance(x),
    }
}

/// Runs the flight and calls `observer` after every accepted step (before events are
/// applied, so the last step may extend past the event time).
pub fn run_observed<F: FieldSolver>(
    scn: &Scenario<F>,
    rs: &RunSettings,
    mut observer: impl FnMut(&StepView<'_, F>),
) -> Trajectory {
    run_cancellable(scn, rs, |step| {
        observer(step);
        true
    })
    .expect("never cancelled")
}

/// As `run_observed`, but the observer returns whether to go on; `None` if it stopped the
/// flight (e.g. because the setup changed and the result is no longer wanted).
pub fn run_cancellable<F: FieldSolver>(
    scn: &Scenario<F>,
    rs: &RunSettings,
    mut observer: impl FnMut(&StepView<'_, F>) -> bool,
) -> Option<Trajectory> {
    let p_norm = scn.p0.length();
    let p_ref = if p_norm > 0.0 {
        p_norm
    } else {
        scn.particle.mass
    };
    let e_ref = {
        let k = Kinematics::new(scn.particle.mass, scn.c).kinetic_energy(scn.p0);
        if k > 0.0 { k } else { scn.particle.mass }
    };
    let ode = ParticleOde::new(&scn.field, &scn.particle, scn.c, p_ref)
        .with_radiation_reaction(scn.radiation_reaction, e_ref);
    let q = scn.particle.charge;

    let mut events: Vec<Event> = (0..scn.obstacles.len()).map(Event::Obstacle).collect();
    if scn.bounds.is_some() {
        events.push(Event::Bounds);
    }
    if scn.detector.is_some() {
        events.push(Event::Detector);
    }

    // With conductors, the particle's interaction with its own induced charges adds
    // ½ q φ_self (PHYSICS.md §2.6).
    // A magnetic moment adds −m B_z (PHYSICS.md §3.2).
    let energy = |x: DVec3, p: DVec3, t: f64| {
        let f = scn.field.sample(x, t);
        ode.kin.kinetic_energy(p) + q * f.phi + 0.5 * q * scn.field.self_field(x, q).1
            - scn.particle.moment * f.b.z
    };
    let power = |x: DVec3, p: DVec3, t: f64| radiated_power(&ode, x, p, t);

    // In time-dependent fields the energy is not conserved: no diagnostic (NaN).
    let is_static = scn.field.is_static();
    let kinetic_initial = ode.kin.kinetic_energy(scn.p0);
    let energy_initial = energy(scn.x0, scn.p0, 0.0);
    let start = Sample {
        t: 0.0,
        x: scn.x0,
        p: scn.p0,
    };
    let mut traj = Trajectory {
        outcome: Outcome::Timeout,
        end: start,
        samples: if rs.record { vec![start] } else { Vec::new() },
        stats: Stats::default(),
        energy_initial,
        kinetic_initial,
        energy_max_abs_error: if is_static { 0.0 } else { f64::NAN },
        radiated_energy: 0.0,
        radiation_work: 0.0,
        reaction_ratio_max: 0.0,
        margins: None,
        radiation: None,
        emission: Vec::new(),
    };

    // Dense samples of the flight for its radiation (only with a radiation goal).
    let window = scn.acceptance.and_then(|a| a.radiation);
    let emission_at = |x: DVec3, p: DVec3, t: f64| -> Emission {
        let mut f = ode.force(x, p, t);
        if ode.radiation_reaction {
            f += ode.radiation_reaction_force(x, p, t);
        }
        let v = ode.kin.velocity(p);
        let a = (f - v * (v.dot(f) / (scn.c * scn.c))) / (ode.kin.gamma(p) * scn.particle.mass);
        (t, x, v, a)
    };
    let mut emission: Vec<Emission> = Vec::new();
    if window.is_some() {
        emission.push(emission_at(scn.x0, scn.p0, 0.0));
    }

    let mut gates = GateTracker::new(&scn.gates, scn.x0);
    // Largest Lorentz and radiation-reaction forces so far (`reaction_ratio_max`).
    let (mut lorentz_max, mut rr_max) = (0.0_f64, 0.0_f64);

    let mut g_prev: Vec<f64> = events
        .iter()
        .map(|&e| event_value(scn, e, scn.x0))
        .collect();
    let mut margin: Vec<f64> = g_prev.clone();
    if let Some(k) = g_prev.iter().position(|&g| g <= 0.0) {
        traj.outcome = events[k].outcome();
        traj.margins = rs.margins.then(|| collect_margins(&events, &margin));
        return Some(traj);
    }

    let settings = Settings {
        rtol: rs.rtol,
        atol: rs.atol,
        max_steps: rs.max_steps,
        dense: true,
        ..Settings::default()
    };
    let mut int = Dop853::new(&ode, 0.0, &ode.pack(scn.x0, scn.p0), settings);
    let mut p_prev = power(scn.x0, scn.p0, 0.0);
    let mut g_next = vec![0.0; events.len()];
    // Triggered event whose penetration depth must be followed past the event step:
    // (event index, depth so far, minimum not yet reached).
    let mut pending_depth: Option<(usize, f64, bool)> = None;

    loop {
        let t_a = int.t();
        let reached_end = match int.step(&ode, scn.t_max, f64::INFINITY) {
            Ok(done) => done,
            Err(e) => {
                traj.outcome = Outcome::Failed(e);
                break;
            }
        };
        let t_b = int.t();
        let dense = int.dense();
        let view = StepView { ode: &ode, dense };
        if !observer(&view) {
            return None;
        }

        // Speed bound on the step: largest sampled speed with a safety margin, never more
        // than c (PHYSICS.md §6).
        let v_max = (0..=4)
            .map(|i| {
                let t = t_a + (t_b - t_a) * f64::from(i) * 0.25;
                ode.kin.velocity(view.state(t).1).length()
            })
            .fold(0.0, f64::max)
            * 1.25;
        let v_max = v_max.min(scn.c);

        let x_b = ParticleOde::<&F>::position(int.y());
        for (k, &ev) in events.iter().enumerate() {
            g_next[k] = event_value(scn, ev, x_b);
        }

        let mut first: Option<(f64, Event)> = None;
        for (k, &ev) in events.iter().enumerate() {
            let mut g = |t: f64| event_value(scn, ev, view.state(t).0);
            let r = first_crossing(&mut g, t_a, t_b, g_prev[k], g_next[k], v_max);
            if let Some(t) = r.time
                && first.is_none_or(|(tf, _)| t < tf)
            {
                first = Some((t, ev));
            }
        }

        let (t_end, event) = match first {
            Some((t, ev)) => (t, Some(ev)),
            None => (t_b, None),
        };

        if rs.margins {
            for (k, &ev) in events.iter().enumerate() {
                let g = |t: f64| event_value(scn, ev, view.state(t).0);
                let m = if event == Some(ev) {
                    // Penetration depth along the continued trajectory; followed past
                    // this step below if the minimum is not reached within it.
                    let (t_min, m) = minimize_on(&g, t_a, t_b);
                    pending_depth = Some((k, m, minimum_at_end(t_min, t_a, t_b)));
                    m
                } else if event.is_some() {
                    // Closest approach up to the event time only.
                    minimize_on(&g, t_a, t_end).1
                } else {
                    let lower_bound = 0.5 * (g_prev[k] + g_next[k] - v_max * (t_b - t_a));
                    if lower_bound < 2.0 * MARGIN_SAFE {
                        minimize_on(&g, t_a, t_b).1
                    } else {
                        g_prev[k].min(g_next[k])
                    }
                };
                margin[k] = margin[k].min(m);
            }
        }

        gates.advance(|t| view.state(t), &ode.kin, t_a, t_end, v_max, rs.margins);

        let (x, p) = if event.is_some() {
            view.state(t_end)
        } else {
            (x_b, ode.momentum(int.y()))
        };
        let sample = Sample { t: t_end, x, p };
        if let Some(w) = window {
            // At least 8 pieces per step, and short enough that the phase ω(t − n·r/c)
            // advances by at most PHASE_STEP at the band's top: its rate 1 − n·v/c,
            // largest over the arc and 5 times in the step, with a margin of 1.5 (much
            // below 2 for directions near the motion, where the radiation is beamed). An
            // even number, so that Simpson's pairs (`lienard_energy`) stay in a step.
            let h = t_end - t_a;
            let rate = (0..=4)
                .map(|i| {
                    let t = t_a + h * f64::from(i) * 0.25;
                    w.phase_rate(ode.kin.velocity(view.state(t).1), scn.c)
                })
                .fold(0.0, f64::max)
                * 1.5;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let pieces = ((w.top_frequency() * rate * h / PHASE_STEP).ceil() as usize).max(8);
            let pieces = pieces + pieces % 2;
            for i in 1..=pieces {
                #[allow(clippy::cast_precision_loss)]
                let t = t_a + h * (i as f64 / pieces as f64);
                let (x, p) = view.state(t);
                emission.push(emission_at(x, p, t));
            }
        }

        let work = if event.is_some() {
            let mut y = vec![0.0; ode.dim()];
            view.dense.eval(t_end, &mut y);
            ode.radiation_work(&y)
        } else {
            ode.radiation_work(int.y())
        };
        traj.radiation_work = work;
        if ode.radiation_reaction {
            lorentz_max = lorentz_max.max(ode.force(x, p, t_end).length());
            rr_max = rr_max.max(ode.radiation_reaction_force(x, p, t_end).length());
            if lorentz_max > 0.0 {
                traj.reaction_ratio_max = rr_max / lorentz_max;
            }
        }
        if is_static {
            // With radiation reaction, W(t) − W(0) equals the work of that force.
            let dw = (energy(x, p, t_end) - energy_initial - work).abs();
            traj.energy_max_abs_error = traj.energy_max_abs_error.max(dw);
        }
        let p_now = power(x, p, t_end);
        traj.radiated_energy += 0.5 * (p_prev + p_now) * (t_end - t_a);
        p_prev = p_now;
        if rs.record {
            traj.samples.push(sample);
        }
        traj.end = sample;

        if let Some(ev) = event {
            traj.outcome = ev.outcome();
            break;
        }
        if reached_end {
            traj.outcome = Outcome::Timeout;
            break;
        }
        std::mem::swap(&mut g_prev, &mut g_next);
    }
    traj.stats = int.stats();
    // Acceptance of the detector: decided at the moment of entry.
    let mut acceptance_margin = None;
    if traj.outcome == Outcome::Arrived
        && let Some(k) = gates.missing()
    {
        traj.outcome = Outcome::SkippedGate(k);
    }
    if traj.outcome == Outcome::Arrived
        && let Some(acc) = scn.acceptance
    {
        let v = ode.kin.velocity(traj.end.p);
        let mut m = acc.margin(v, ode.kin.kinetic_energy(traj.end.p));
        if let Some(w) = acc.radiation {
            let e = w.measure(&emission, q, scn.c);
            traj.radiation = Some(e);
            m = m.min(w.margin(e));
        }
        acceptance_margin = Some(m);
        if m < 0.0 {
            traj.outcome = Outcome::Rejected;
        }
    }

    // Follow the continued trajectory (as if the boundary were not there) until the
    // event function reaches its minimum or the depth is clearly safe. The depth is then
    // a property of the trajectory, independent of where the steps happened to end.
    if let Some((k, mut depth, mut unfinished)) = pending_depth {
        let ev = events[k];
        let mut steps = 0;
        while unfinished && depth > -MARGIN_SAFE && steps < 10_000 {
            steps += 1;
            let t_a = int.t();
            if int.step(&ode, f64::MAX, f64::INFINITY).is_err() {
                // Integration broke down inside the obstacle (e.g. at a point charge):
                // the penetration is certainly deep.
                depth = f64::NEG_INFINITY;
                break;
            }
            let view = StepView {
                ode: &ode,
                dense: int.dense(),
            };
            let g = |t: f64| event_value(scn, ev, view.state(t).0);
            let (t_min, m) = minimize_on(&g, t_a, int.t());
            depth = depth.min(m);
            unfinished = minimum_at_end(t_min, t_a, int.t());
        }
        margin[k] = margin[k].min(depth);
    }
    traj.emission = emission;
    traj.margins = rs.margins.then(|| {
        let mut m = collect_margins(&events, &margin);
        m.acceptance = acceptance_margin;
        (m.gates, m.gate_acceptance) = gates.into_margins();
        m
    });
    Some(traj)
}

/// Progress of one flight through its gates (PHYSICS.md §6.2): the next gate to pass, the
/// time from which each gate's margin is tracked (when the previous one was passed), the
/// margins and the acceptance margins. Shared by the single-particle and beam runners.
pub(crate) struct GateTracker<'a> {
    gates: &'a [Gate],
    next: usize,
    from: Vec<f64>,
    margin: Vec<f64>,
    acc: Vec<Option<f64>>,
    /// Signed distance to the next gate at the end of the last step.
    prev: f64,
}

impl<'a> GateTracker<'a> {
    pub(crate) fn new(gates: &'a [Gate], x0: DVec3) -> Self {
        let n = gates.len();
        Self {
            gates,
            next: 0,
            from: vec![0.0; n],
            margin: vec![f64::INFINITY; n],
            acc: vec![None; n],
            prev: gates
                .first()
                .map_or(f64::INFINITY, |g| g.region.signed_distance(x0)),
        }
    }

    fn value(&self, k: usize, x: DVec3) -> f64 {
        self.gates[k].region.signed_distance(x)
    }

    /// Gates crossed on `[t_a, t_end]` of a step (up to a terminal event), in order; then
    /// the margins of the gates whose tracking has started. `state(t)` is the dense
    /// output `(x, p)`, `v_max` the step's speed bound.
    pub(crate) fn advance(
        &mut self,
        state: impl Fn(f64) -> (DVec3, DVec3),
        kin: &Kinematics,
        t_a: f64,
        t_end: f64,
        v_max: f64,
        margins: bool,
    ) {
        let n = self.gates.len();
        if n == 0 {
            return;
        }
        let mut t_from = t_a;
        while self.next < n {
            let k = self.next;
            let g_end = self.value(k, state(t_end).0);
            if self.prev <= 0.0 {
                // Inside after an entry outside its conditions: wait until it leaves.
                self.prev = g_end;
                break;
            }
            let mut g = |t: f64| self.value(k, state(t).0);
            match first_crossing(&mut g, t_from, t_end, self.prev, g_end, v_max).time {
                Some(tc) => {
                    let (_, p) = state(tc);
                    let ok = match self.gates[k].acceptance {
                        Some(acc) => {
                            let m = acc.margin(kin.velocity(p), kin.kinetic_energy(p));
                            self.acc[k] = Some(m);
                            m >= 0.0
                        }
                        None => true,
                    };
                    if !ok {
                        self.prev = g_end;
                        break;
                    }
                    self.next += 1;
                    t_from = tc;
                    if self.next < n {
                        self.from[self.next] = tc;
                        self.prev = self.value(self.next, state(tc).0);
                    }
                }
                None => {
                    self.prev = g_end;
                    break;
                }
            }
        }
        if margins {
            for k in 0..n.min(self.next + 1) {
                let from = self.from[k].max(t_a);
                if from < t_end {
                    let g = |t: f64| self.value(k, state(t).0);
                    self.margin[k] = self.margin[k].min(minimize_on(&g, from, t_end).1);
                }
            }
        }
    }

    /// The first gate not passed, if any.
    pub(crate) fn missing(&self) -> Option<usize> {
        (self.next < self.gates.len()).then_some(self.next)
    }

    /// Gate margins and acceptance margins (`Margins::gates`, `Margins::gate_acceptance`).
    pub(crate) fn into_margins(self) -> (Vec<f64>, Vec<Option<f64>>) {
        (self.margin, self.acc)
    }
}

/// Whether a minimum found on `[a, b]` lies at the right end, i.e. the function is still
/// decreasing there. Decided from the location of the minimum, not by comparing values:
/// the integrator state and the dense output at `b` differ by rounding.
pub(crate) fn minimum_at_end(t_min: f64, a: f64, b: f64) -> bool {
    t_min >= b - 1e-6 * (b - a)
}

fn collect_margins(events: &[Event], margin: &[f64]) -> Margins {
    let mut m = Margins {
        obstacles: Vec::new(),
        bounds: None,
        detector: None,
        acceptance: None,
        gates: Vec::new(),
        gate_acceptance: Vec::new(),
    };
    for (&ev, &v) in events.iter().zip(margin) {
        match ev {
            Event::Obstacle(_) => m.obstacles.push(v),
            Event::Bounds => m.bounds = Some(v),
            Event::Detector => m.detector = Some(v),
        }
    }
    m
}

/// Minimum `(t, g(t))` of a smooth function on `[a, b]`: 32 samples, then golden-section
/// search in the bracket around the smallest sample.
pub(crate) fn minimize_on(g: &impl Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    const N: u32 = 32;
    let t_at = |i: u32| a + (b - a) * f64::from(i) / f64::from(N);
    let (mut i_min, mut g_min) = (0, g(a));
    let mut t_min = a;
    for i in 1..=N {
        let v = g(t_at(i));
        if v < g_min {
            i_min = i;
            g_min = v;
            t_min = t_at(i);
        }
    }
    let (mut lo, mut hi) = (t_at(i_min.saturating_sub(1)), t_at((i_min + 1).min(N)));
    let r = 0.5 * (5f64.sqrt() - 1.0);
    let mut x1 = hi - r * (hi - lo);
    let mut x2 = lo + r * (hi - lo);
    let (mut f1, mut f2) = (g(x1), g(x2));
    for _ in 0..80 {
        if f1 < f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - r * (hi - lo);
            f1 = g(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + r * (hi - lo);
            f2 = g(x2);
        }
        if hi - lo <= f64::EPSILON * hi.abs().max(lo.abs()) {
            break;
        }
    }
    [(x1, f1), (x2, f2)].into_iter().fold(
        (t_min, g_min),
        |best, cand| if cand.1 < best.1 { cand } else { best },
    )
}

/// Liénard power `P = (2/3) q² γ⁶ (a² − |v×a|²/c²) / c³` (internal units, `k = 1`).
fn radiated_power<F: FieldSolver>(ode: &ParticleOde<F>, x: DVec3, p: DVec3, t: f64) -> f64 {
    let c = ode.kin.c;
    if c.is_infinite() {
        return 0.0;
    }
    let gamma = ode.kin.gamma(p);
    let v = ode.kin.velocity(p);
    let force = ode.force(x, p, t);
    // dv/dt = (F − v (v·F)/c²) / (γ m)
    let a = (force - v * (v.dot(force) / (c * c))) / (gamma * ode.kin.mass);
    let q = ode.charge;
    (2.0 / 3.0)
        * q
        * q
        * gamma.powi(6)
        * (a.length_squared() - v.cross(a).length_squared() / (c * c))
        / (c * c * c)
}
