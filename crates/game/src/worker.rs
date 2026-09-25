//! Background physics thread. The editor submits the latest setup; stale requests are
//! dropped. For each request the worker sends the preview trajectory first, then the
//! verification verdict (SPEC §2.3), so rendering never waits for physics.

use std::sync::Mutex;
use std::sync::mpsc::{Receiver, Sender, channel};

use physics::DVec3;
use physics::field::{Coulomb, FieldSolver};
use physics::trajectory::{Outcome, RunSettings, Scenario, run, run_observed};
use physics::verify::{Status, Tolerances, classify};

pub struct Request {
    pub revision: u64,
    pub scenario: Scenario<Coulomb>,
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
    pub revision: u64,
    pub path: Vec<PathPoint>,
    pub outcome: Outcome,
    pub flight_time: f64,
    pub energy_rel_error: f64,
    pub radiated_fraction: f64,
    pub max_speed_over_c: f64,
}

#[derive(Clone, Debug)]
pub enum Response {
    Preview(Preview),
    Verified {
        revision: u64,
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
    loop {
        let req = match next.take() {
            Some(r) => r,
            None => match rx.recv() {
                Ok(r) => r,
                Err(_) => return,
            },
        };
        let req = latest(rx, req);
        let scn = &req.scenario;

        let preview_settings = RunSettings::with_tolerance(req.tolerances.preview);
        let mut dense_points: Vec<(f64, DVec3, DVec3)> = Vec::new();
        let preview = run_observed(scn, &preview_settings, |step| {
            let (a, b) = (step.t_start(), step.t_end());
            for i in 1..=8 {
                let t = a + (b - a) * f64::from(i) / 8.0;
                let (x, p) = step.state(t);
                dense_points.push((t, x, p));
            }
        });
        let path = build_path(scn, &preview.samples[0], &dense_points, &preview.end);
        let max_speed_over_c = path.iter().map(|p| p.speed_over_c).fold(0.0, f64::max);
        if tx
            .send(Response::Preview(Preview {
                revision: req.revision,
                path,
                outcome: preview.outcome,
                flight_time: preview.end.t,
                energy_rel_error: preview.energy_max_abs_error / preview.kinetic_initial,
                radiated_fraction: preview.radiated_energy / preview.kinetic_initial,
                max_speed_over_c,
            }))
            .is_err()
        {
            return;
        }

        // Skip verification if the setup has already changed.
        if let Ok(newer) = rx.try_recv() {
            next = Some(newer);
            continue;
        }
        let verified = run(scn, &RunSettings::with_tolerance(req.tolerances.verify));
        let status = classify(&preview, &verified, scn.t_max);
        if tx
            .send(Response::Verified {
                revision: req.revision,
                status,
                outcome: verified.outcome,
            })
            .is_err()
        {
            return;
        }
    }
}

fn build_path(
    scn: &Scenario<Coulomb>,
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
