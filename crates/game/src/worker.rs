//! Background physics thread. The editor submits the latest setup (one scenario per
//! shot); stale requests are dropped. For each request the worker sends every shot's
//! preview trajectory first, then every shot's verification verdict (SPEC §2.3), so
//! rendering never waits for physics.

use std::sync::Mutex;
use std::sync::mpsc::{Receiver, Sender, channel};

use physics::DVec3;
use physics::field::{FieldSolver, StaticField};
use physics::trajectory::{Outcome, RunSettings, Scenario, Trajectory, run, run_observed};
use physics::verify::{Status, Tolerances, classify};

pub struct Request {
    pub revision: u64,
    pub scenarios: Vec<Scenario<StaticField>>,
    pub tolerances: Tolerances,
}

#[derive(Clone, Copy, Debug)]
pub struct PathPoint {
    pub t: f64,
    pub x: DVec3,
    pub p: DVec3,
    pub kinetic: f64,
    /// Potential energy relative to the launch point, `q(φ(x) − φ(A))`.
    pub potential: f64,
    pub force: DVec3,
    pub speed_over_c: f64,
}

#[derive(Clone, Debug)]
pub struct Preview {
    pub path: Vec<PathPoint>,
    pub outcome: Outcome,
    pub flight_time: f64,
    pub energy_rel_error: f64,
    pub radiated_fraction: f64,
    pub max_speed_over_c: f64,
}

#[derive(Clone, Debug)]
pub enum Response {
    Preview {
        revision: u64,
        shot: usize,
        preview: Preview,
    },
    Verified {
        revision: u64,
        shot: usize,
        status: Status,
        outcome: Outcome,
    },
}

pub struct Worker {
    tx: Sender<Request>,
    rx: Mutex<Receiver<Response>>,
}

impl Worker {
    pub fn spawn() -> Self {
        let (req_tx, req_rx) = channel::<Request>();
        let (resp_tx, resp_rx) = channel::<Response>();
        std::thread::Builder::new()
            .name("physics".into())
            .spawn(move || worker_loop(&req_rx, &resp_tx))
            .expect("spawn physics thread");
        Self {
            tx: req_tx,
            rx: Mutex::new(resp_rx),
        }
    }

    pub fn submit(&self, req: Request) {
        let _ = self.tx.send(req);
    }

    pub fn poll(&self) -> Vec<Response> {
        self.rx.lock().unwrap().try_iter().collect()
    }
}

fn latest(rx: &Receiver<Request>, mut req: Request) -> Request {
    while let Ok(newer) = rx.try_recv() {
        req = newer;
    }
    req
}

fn worker_loop(rx: &Receiver<Request>, tx: &Sender<Response>) {
    let mut next: Option<Request> = None;
    'requests: loop {
        let req = match next.take() {
            Some(r) => r,
            None => match rx.recv() {
                Ok(r) => r,
                Err(_) => return,
            },
        };
        let req = latest(rx, req);

        let mut previews: Vec<Trajectory> = Vec::new();
        for (shot, scn) in req.scenarios.iter().enumerate() {
            let (traj, preview) = preview_shot(scn, req.tolerances.preview);
            previews.push(traj);
            let msg = Response::Preview {
                revision: req.revision,
                shot,
                preview,
            };
            if tx.send(msg).is_err() {
                return;
            }
            if let Ok(newer) = rx.try_recv() {
                next = Some(newer);
                continue 'requests;
            }
        }
        for (shot, scn) in req.scenarios.iter().enumerate() {
            // Skip verification if the setup has already changed.
            if let Ok(newer) = rx.try_recv() {
                next = Some(newer);
                continue 'requests;
            }
            let verified = run(scn, &RunSettings::with_tolerance(req.tolerances.verify));
            let status = classify(&previews[shot], &verified, scn.t_max);
            let msg = Response::Verified {
                revision: req.revision,
                shot,
                status,
                outcome: verified.outcome,
            };
            if tx.send(msg).is_err() {
                return;
            }
        }
    }
}

fn preview_shot(scn: &Scenario<StaticField>, tol: f64) -> (Trajectory, Preview) {
    let mut dense_points: Vec<(f64, DVec3, DVec3)> = Vec::new();
    let traj = run_observed(scn, &RunSettings::with_tolerance(tol), |step| {
        let (a, b) = (step.t_start(), step.t_end());
        for i in 1..=8 {
            let t = a + (b - a) * f64::from(i) / 8.0;
            let (x, p) = step.state(t);
            dense_points.push((t, x, p));
        }
    });
    let path = build_path(scn, &traj.samples[0], &dense_points, &traj.end);
    let max_speed_over_c = path.iter().map(|p| p.speed_over_c).fold(0.0, f64::max);
    let preview = Preview {
        path,
        outcome: traj.outcome,
        flight_time: traj.end.t,
        energy_rel_error: traj.energy_max_abs_error / traj.kinetic_initial,
        radiated_fraction: traj.radiated_energy / traj.kinetic_initial,
        max_speed_over_c,
    };
    (traj, preview)
}

fn build_path(
    scn: &Scenario<StaticField>,
    start: &physics::trajectory::Sample,
    dense: &[(f64, DVec3, DVec3)],
    end: &physics::trajectory::Sample,
) -> Vec<PathPoint> {
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let q = scn.particle.charge;
    let phi_a = scn.field.sample(scn.x0, 0.0).phi;
    let point = |t: f64, x: DVec3, p: DVec3| {
        let f = scn.field.sample(x, t);
        let v = kin.velocity(p);
        PathPoint {
            t,
            x,
            p,
            kinetic: kin.kinetic_energy(p),
            potential: q * (f.phi - phi_a),
            force: (f.e + v.cross(f.b)) * q,
            speed_over_c: if scn.c.is_finite() {
                v.length() / scn.c
            } else {
                0.0
            },
        }
    };
    let mut path = vec![point(start.t, start.x, start.p)];
    path.extend(
        dense
            .iter()
            .filter(|(t, _, _)| *t < end.t)
            .map(|&(t, x, p)| point(t, x, p)),
    );
    path.push(point(end.t, end.x, end.p));
    path
}
