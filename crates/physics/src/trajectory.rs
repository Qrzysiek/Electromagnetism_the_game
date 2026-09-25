//! Integration of a single test-particle trajectory with events and diagnostics
//! (PHYSICS.md §3–§6, §8).

use glam::DVec3;

use crate::dynamics::{Particle, ParticleOde};
use crate::events::first_crossing;
use crate::field::FieldSolver;
use crate::geometry::{Aabb, Region, Sphere};
use crate::integrator::{self, Dense, Dop853, Settings, Stats};

/// Everything that defines one flight.
#[derive(Clone, Debug)]
pub struct Scenario<F> {
    pub field: F,
    /// Solid spheres (fixed charges); touching one loses the particle.
    pub obstacles: Vec<Sphere>,
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
}

#[derive(Clone, Copy, Debug)]
pub struct RunSettings {
    pub rtol: f64,
    pub atol: f64,
    pub max_steps: u64,
    /// Record the state at every accepted step.
    pub record: bool,
}

impl RunSettings {
    pub fn with_tolerance(tol: f64) -> Self {
        Self {
            rtol: tol,
            atol: tol,
            max_steps: 1_000_000,
            record: true,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outcome {
    Arrived,
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

#[derive(Clone, Debug)]
pub struct Trajectory {
    pub outcome: Outcome,
    /// Final state (at the event, or at `t_max`).
    pub end: Sample,
    /// States at accepted steps (if recorded), ending with `end`.
    pub samples: Vec<Sample>,
    pub stats: Stats,
    /// Total energy `W = (γ−1)mc² + qφ` at the start.
    pub energy_initial: f64,
    /// Kinetic energy at the start.
    pub kinetic_initial: f64,
    /// Largest `|W(t) − W(0)|` over the accepted steps (static fields only).
    pub energy_max_abs_error: f64,
    /// Energy radiated according to the Liénard formula (trapezoidal rule over steps).
    pub radiated_energy: f64,
    /// Per obstacle: smallest surface distance among evaluated points (a sampled upper
    /// bound of the true closest approach).
    pub closest_sampled: Vec<f64>,
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
    let p_norm = scn.p0.length();
    let p_ref = if p_norm > 0.0 {
        p_norm
    } else {
        scn.particle.mass
    };
    let ode = ParticleOde::new(&scn.field, &scn.particle, scn.c, p_ref);
    let q = scn.particle.charge;

    let mut events: Vec<Event> = (0..scn.obstacles.len()).map(Event::Obstacle).collect();
    if scn.bounds.is_some() {
        events.push(Event::Bounds);
    }
    if scn.detector.is_some() {
        events.push(Event::Detector);
    }

    let energy =
        |x: DVec3, p: DVec3, t: f64| ode.kin.kinetic_energy(p) + q * scn.field.sample(x, t).phi;
    let power = |x: DVec3, p: DVec3, t: f64| radiated_power(&ode, x, p, t);

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
        energy_max_abs_error: 0.0,
        radiated_energy: 0.0,
        closest_sampled: vec![f64::INFINITY; scn.obstacles.len()],
    };

    let mut g_prev: Vec<f64> = events
        .iter()
        .map(|&e| event_value(scn, e, scn.x0))
        .collect();
    record_closest(&mut traj.closest_sampled, &events, &g_prev);
    if let Some(k) = g_prev.iter().position(|&g| g <= 0.0) {
        traj.outcome = events[k].outcome();
        return traj;
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
        observer(&view);

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
            if let Event::Obstacle(i) = ev {
                traj.closest_sampled[i] = traj.closest_sampled[i].min(r.min_sampled);
            }
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
        let (x, p) = if event.is_some() {
            view.state(t_end)
        } else {
            (x_b, ode.momentum(int.y()))
        };
        let sample = Sample { t: t_end, x, p };

        let dw = (energy(x, p, t_end) - energy_initial).abs();
        traj.energy_max_abs_error = traj.energy_max_abs_error.max(dw);
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
    traj
}

fn record_closest(closest: &mut [f64], events: &[Event], g: &[f64]) {
    for (&ev, &gv) in events.iter().zip(g) {
        if let Event::Obstacle(i) = ev {
            closest[i] = closest[i].min(gv);
        }
    }
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
