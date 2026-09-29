//! Background physics thread. The editor submits the latest setup (level and player
//! elements); stale requests are dropped. The worker builds the scenarios (metal spheres
//! need a one-off factorization per geometry, which must not stall rendering), sends
//! every flight's preview trajectory first, then every flight's verification verdict
//! (SPEC §2.3). Verification uses the field at verification resolution (it differs from
//! the preview's only with metal spheres, PHYSICS.md §2.6). Every flight is timed for the
//! sandbox's resource meters (`level::cost`).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use physics::DVec3;
use physics::conductor::Resolution;
use physics::field::{FieldSolver, LevelField};
use physics::trajectory::{Outcome, RunSettings, Scenario, Trajectory, run_cancellable};
use physics::verify::{Status, Tolerances, classify};

pub struct Request {
    pub revision: u64,
    pub level: level::Level,
    pub placement: Vec<level::Element>,
}

#[derive(Clone, Copy, Debug)]
pub struct PathPoint {
    pub t: f64,
    pub x: DVec3,
    pub p: DVec3,
    pub kinetic: f64,
    /// Potential energy relative to the launch point, `q(φ(x) − φ(A)) − m(B_z(x) − B_z(A))`.
    pub potential: f64,
    /// Energy lost to radiation so far (minus the work of the radiation-reaction force;
    /// 0 when radiation reaction is off).
    pub radiated: f64,
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
    /// Radiation reaction included in the dynamics.
    pub radiation_reaction: bool,
    /// Energy lost to radiation (with radiation reaction), relative to T₀.
    pub radiation_loss_fraction: f64,
    /// Largest |F_RR| / |F_Lorentz| (validity of the Landau–Lifshitz treatment).
    pub reaction_ratio_max: f64,
    /// Bound on the electrodes' neglected image force along the path, relative to the
    /// force that matters (`level::electrode_image_force_bound`); 0 without electrodes.
    pub image_force_bound: f64,
    /// With a radiation goal, on arrival: the measured energy per steradian.
    pub radiation: Option<f64>,
    /// With a radiation goal: the spectrum over the goal's directions, `(ω, d²I/dω dΩ)`,
    /// from 0 to `spectrum_top` (twice the band's top; without a band, three times the
    /// flight's largest critical frequency).
    pub spectrum: Vec<(f64, f64)>,
    pub spectrum_top: f64,
}

/// Preview of a beam flight (`level::beam`): every particle's path `(t, x)`, its shot
/// and its outcome.
#[derive(Clone, Debug)]
pub struct BeamPreview {
    pub paths: Vec<Vec<(f64, DVec3)>>,
    pub shots: Vec<usize>,
    pub outcomes: Vec<Outcome>,
    /// Energy drift of the whole system (between removals), relative; NaN at finite c
    /// with interaction (the field carries energy).
    pub energy_rel_error: f64,
    /// Largest radiated energy of a particle (Liénard), relative to its launch energy.
    pub radiated_max: f64,
    /// Quasi-static interaction: the largest estimated relative error of a particle's
    /// interaction impulse (`BeamRun::neglected_retardation`).
    pub retardation_max: f64,
    /// Per particle: world-line samples `(t, x, v, a)` (for the field views; a particle
    /// that left the arena flies on, and its samples go on), and where and when it was
    /// absorbed (None if it flew to the time limit or left the arena).
    pub worldlines: Vec<Vec<(f64, DVec3, DVec3, DVec3)>>,
    pub ends: Vec<Option<(f64, DVec3)>>,
    /// Per particle: its fade in a screening cup, if it entered one (its world line then
    /// goes on uniformly into the cup).
    pub fades: Vec<Option<physics::beam::Fade>>,
    /// The beam's energy budget along the flight.
    pub energy: Vec<physics::beam::EnergySample>,
    /// Flown with the exact retarded interaction.
    pub retarded: bool,
}

#[derive(Clone, Debug)]
pub enum Response {
    Preview {
        revision: u64,
        /// Flight index (`Level::flight_of`).
        flight: usize,
        preview: Preview,
    },
    Verified {
        revision: u64,
        /// Flight index (`Level::flight_of`).
        flight: usize,
        status: Status,
        outcome: Outcome,
    },
    /// Beam levels: the preview of flight `flight` (one per disturbance).
    BeamPreview {
        revision: u64,
        flight: usize,
        preview: BeamPreview,
    },
    /// Beam levels: each particle's verification verdict and outcome, and the verdict's
    /// own flight when its model differs from the preview's (the exact retarded
    /// interaction at finite c, metal at verification resolution): the views show it
    /// once it arrives.
    BeamVerified {
        revision: u64,
        flight: usize,
        results: Vec<(Status, Outcome)>,
        run: Option<BeamPreview>,
    },
    /// Measured cost so far (after every flight).
    Cost {
        revision: u64,
        cost: level::cost::Cost,
    },
}

pub struct Worker {
    tx: Sender<Request>,
    rx: Mutex<Receiver<Response>>,
    /// Revision of the newest submitted setup: flights of older ones are abandoned
    /// mid-flight, so a stale (possibly slow) computation never delays the current one.
    latest: Arc<AtomicU64>,
}

impl Worker {
    pub fn spawn() -> Self {
        let (req_tx, req_rx) = channel::<Request>();
        let (resp_tx, resp_rx) = channel::<Response>();
        let latest = Arc::new(AtomicU64::new(0));
        let seen = latest.clone();
        std::thread::Builder::new()
            .name("physics".into())
            .spawn(move || worker_loop(&req_rx, &resp_tx, &seen))
            .expect("spawn physics thread");
        Self {
            tx: req_tx,
            rx: Mutex::new(resp_rx),
            latest,
        }
    }

    pub fn submit(&self, req: Request) {
        self.latest.store(req.revision, Ordering::Release);
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

/// Sends responses to the game and records them, so that a finished computation can be
/// replayed for an identical setup.
struct Sink<'a> {
    tx: &'a Sender<Response>,
    record: std::cell::RefCell<Vec<Response>>,
}

impl Sink<'_> {
    /// `Err` if the game has gone away.
    fn send(&self, r: Response) -> Result<(), ()> {
        self.record.borrow_mut().push(r.clone());
        self.tx.send(r).map_err(|_| ())
    }
}

/// Results of recent setups, keyed by their content (level and placement): editing a
/// setup back to one computed before (undo, moving an element back, switching levels
/// back and forth) replays its results instead of computing them again. Only complete
/// computations are kept, the most recent first.
struct ResultCache {
    entries: std::collections::VecDeque<(u64, Vec<Response>)>,
}

/// How many setups the cache keeps.
const CACHE_SETUPS: usize = 24;

impl ResultCache {
    fn key(req: &Request) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        serde_json::to_string(&(&req.level, &req.placement))
            .unwrap_or_default()
            .hash(&mut h);
        h.finish()
    }

    fn get(&mut self, key: u64) -> Option<Vec<Response>> {
        let i = self.entries.iter().position(|e| e.0 == key)?;
        let entry = self.entries.remove(i)?;
        let out = entry.1.clone();
        self.entries.push_front(entry);
        Some(out)
    }

    fn put(&mut self, key: u64, responses: Vec<Response>) {
        self.entries.retain(|e| e.0 != key);
        self.entries.push_front((key, responses));
        self.entries.truncate(CACHE_SETUPS);
    }
}

/// A recorded response for another revision (the same setup).
fn retag(r: Response, revision: u64) -> Response {
    match r {
        Response::Preview {
            flight, preview, ..
        } => Response::Preview {
            revision,
            flight,
            preview,
        },
        Response::Verified {
            flight,
            status,
            outcome,
            ..
        } => Response::Verified {
            revision,
            flight,
            status,
            outcome,
        },
        Response::BeamPreview {
            flight, preview, ..
        } => Response::BeamPreview {
            revision,
            flight,
            preview,
        },
        Response::BeamVerified {
            flight,
            results,
            run,
            ..
        } => Response::BeamVerified {
            revision,
            flight,
            results,
            run,
        },
        Response::Cost { cost, .. } => Response::Cost { revision, cost },
    }
}

/// Runs a metal-field build that a newer setup (another level, an edit) abandons at its
/// next checkpoint (`physics::cancel`); `None` then.
fn unless_stale<T>(newest: &Arc<AtomicU64>, revision: u64, f: impl FnOnce() -> T) -> Option<T> {
    let wanted = newest.clone();
    physics::cancel::cancellable(move || wanted.load(Ordering::Acquire) == revision, f)
}

fn worker_loop(rx: &Receiver<Request>, out: &Sender<Response>, newest: &Arc<AtomicU64>) {
    let mut next: Option<Request> = None;
    let mut cache = ResultCache {
        entries: std::collections::VecDeque::new(),
    };
    'requests: loop {
        let req = match next.take() {
            Some(r) => r,
            None => match rx.recv() {
                Ok(r) => r,
                Err(_) => return,
            },
        };
        let req = latest(rx, req);
        let key = ResultCache::key(&req);
        if let Some(responses) = cache.get(key) {
            for r in responses {
                if out.send(retag(r, req.revision)).is_err() {
                    return;
                }
            }
            continue 'requests;
        }
        let sink = Sink {
            tx: out,
            record: std::cell::RefCell::new(Vec::new()),
        };
        let tx = &sink;
        let current = |r: u64| newest.load(Ordering::Acquire) == r;
        if req.level.has_beams() {
            if !beam_request(&req, tx, newest) {
                return;
            }
            if current(req.revision) {
                cache.put(key, sink.record.into_inner());
            }
            continue 'requests;
        }
        let metal = req.level.has_metal(&req.placement);
        let mut cost = level::cost::Cost::new(&req.level);
        let Some(scenarios) =
            unless_stale(newest, req.revision, || req.level.scenarios(&req.placement))
        else {
            continue 'requests;
        };
        let tolerances: Tolerances = req.level.tolerances();
        if metal {
            // Also builds (and caches) the display systems off the render thread.
            let Some(display) = unless_stale(newest, req.revision, || {
                req.level.field_at(&req.placement, Resolution::Display).0
            }) else {
                continue 'requests;
            };
            cost.add_field(&display);
            cost.add_field(&scenarios[0].field);
        }
        let send_cost = |cost: &level::cost::Cost| {
            tx.send(Response::Cost {
                revision: req.revision,
                cost: *cost,
            })
            .is_ok()
        };
        if !send_cost(&cost) {
            return;
        }

        let mut previews: Vec<Trajectory> = Vec::new();
        for (shot, scn) in scenarios.iter().enumerate() {
            let start = Instant::now();
            let Some((traj, preview)) =
                preview_shot(scn, tolerances.preview, || current(req.revision))
            else {
                continue 'requests;
            };
            cost.add_preview(start.elapsed().as_secs_f64(), &traj);
            if !send_cost(&cost) {
                return;
            }
            previews.push(traj);
            let msg = Response::Preview {
                revision: req.revision,
                flight: shot,
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
        if !current(req.revision) {
            continue 'requests;
        }
        let fine = if metal {
            let Some(fine) = unless_stale(newest, req.revision, || {
                req.level.scenarios_at(&req.placement, Resolution::Verify)
            }) else {
                continue 'requests;
            };
            cost.add_verification_field(&fine[0].field);
            fine
        } else {
            scenarios
        };
        for (shot, scn) in fine.iter().enumerate() {
            // Skip verification if the setup has already changed.
            if let Ok(newer) = rx.try_recv() {
                next = Some(newer);
                continue 'requests;
            }
            let start = Instant::now();
            let Some(verified) =
                run_cancellable(scn, &RunSettings::with_tolerance(tolerances.verify), |_| {
                    current(req.revision)
                })
            else {
                continue 'requests;
            };
            cost.add_verification(start.elapsed().as_secs_f64(), &verified);
            if !send_cost(&cost) {
                return;
            }
            let status = classify(&previews[shot], &verified, scn.t_max);
            let msg = Response::Verified {
                revision: req.revision,
                flight: shot,
                status,
                outcome: verified.outcome,
            };
            if tx.send(msg).is_err() {
                return;
            }
        }
        // Complete (a newer request would have abandoned it): keep for replay.
        if current(req.revision) {
            cache.put(key, sink.record.into_inner());
        }
    }
}

/// A beam level: per flight the preview (every particle's path), then the per-particle
/// verification, with the measured cost. Returns `false` if the game has gone away.
fn beam_request(req: &Request, tx: &Sink<'_>, newest: &Arc<AtomicU64>) -> bool {
    let current = || newest.load(Ordering::Acquire) == req.revision;
    let tol = req.level.tolerances();
    let mut cost = level::cost::Cost::new(&req.level);
    let send_cost = |cost: &level::cost::Cost| {
        tx.send(Response::Cost {
            revision: req.revision,
            cost: *cost,
        })
        .is_ok()
    };
    let shots = req.level.beam_shots(&req.placement);
    let Some(preview_scns) = unless_stale(newest, req.revision, || {
        req.level
            .beam_scenarios(&req.placement, Resolution::Preview)
    }) else {
        return true;
    };
    let mut previews = Vec::new();
    for (flight, scn) in preview_scns.iter().enumerate() {
        let start = Instant::now();
        let Some((run, preview)) = fly_beam(scn, tol.preview, 4, &shots, current) else {
            return true;
        };
        if let Some(t) = run.trajectories.first() {
            cost.add_preview(start.elapsed().as_secs_f64(), t);
        }
        let msg = Response::BeamPreview {
            revision: req.revision,
            flight,
            preview,
        };
        if !send_cost(&cost) || tx.send(msg).is_err() {
            return false;
        }
        previews.push(run);
    }
    // The verdict: metal at verification resolution, and at finite c the exact retarded
    // interaction (the preview above is quasi-static). Its flight is then recorded for
    // the views too (they show it once it arrives).
    let exact = req.level.physics.beam_interaction && req.level.physics.c.is_some();
    let distinct = req.level.has_metal(&req.placement) || exact;
    let fine = if distinct {
        let Some(fine) = unless_stale(newest, req.revision, || {
            req.level.verification_beam_scenarios(&req.placement)
        }) else {
            return true;
        };
        if req.level.has_metal(&req.placement) {
            cost.add_verification_field(&fine[0].field);
        }
        fine
    } else {
        preview_scns
    };
    for (flight, scn) in fine.iter().enumerate() {
        let start = Instant::now();
        let (verified, run) = if distinct {
            // The verification's steps are shorter: fewer points per step.
            let Some((verified, view)) = fly_beam(scn, tol.verify, 2, &shots, current) else {
                return true;
            };
            (verified, Some(view))
        } else {
            let Some(verified) = physics::beam::run_beam_cancellable(
                scn,
                &RunSettings::with_tolerance(tol.verify),
                |_, _, _| current(),
            ) else {
                return true;
            };
            (verified, None)
        };
        if let Some(t) = verified.trajectories.first() {
            cost.add_verification(start.elapsed().as_secs_f64(), t);
        }
        let results = previews[flight]
            .trajectories
            .iter()
            .zip(&verified.trajectories)
            .map(|(a, b)| (classify(a, b, scn.t_max), b.outcome))
            .collect();
        let msg = Response::BeamVerified {
            revision: req.revision,
            flight,
            results,
            run,
        };
        if !send_cost(&cost) || tx.send(msg).is_err() {
            return false;
        }
    }
    true
}

/// A beam flight recorded for the views: every particle's path, world line (position,
/// velocity, and acceleration from the derivative of the dense output, at `per_step + 1`
/// points of every step), outcome, end and fade, and the beam's energy budget. `None` if
/// `go_on` stopped it.
fn fly_beam(
    scn: &physics::beam::BeamScenario<LevelField>,
    tol: f64,
    per_step: u32,
    shots: &[usize],
    go_on: impl Fn() -> bool,
) -> Option<(physics::beam::BeamRun, BeamPreview)> {
    let mut paths: Vec<Vec<(f64, DVec3)>> =
        scn.particles.iter().map(|b| vec![(0.0, b.x0)]).collect();
    let kins: Vec<physics::dynamics::Kinematics> = scn
        .particles
        .iter()
        .map(|b| physics::dynamics::Kinematics::new(b.particle.mass, scn.c))
        .collect();
    let mut lines: Vec<Vec<(f64, DVec3, DVec3, DVec3)>> = vec![Vec::new(); scn.particles.len()];
    let run = physics::beam::run_beam_cancellable(
        scn,
        &RunSettings::with_tolerance(tol),
        |dense, members, p_ref| {
            let (a, b) = (dense.t_start(), dense.t_end());
            let at = |s: u32| a + (b - a) * f64::from(s) / f64::from(per_step);
            for (k, &i) in members.iter().enumerate() {
                for s in 0..=per_step {
                    let t = at(s);
                    if lines[i].last().is_some_and(|l| t <= l.0) {
                        continue;
                    }
                    let comp = |c: usize| dense.eval_component(6 * k + c, t);
                    let der = |c: usize| dense.eval_derivative_component(6 * k + c, t);
                    let x = DVec3::new(comp(0), comp(1), comp(2));
                    let p = DVec3::new(comp(3), comp(4), comp(5)) * p_ref;
                    let dp = DVec3::new(der(3), der(4), der(5)) * p_ref;
                    lines[i].push((t, x, kins[i].velocity(p), kins[i].acceleration(p, dp)));
                }
                for s in 1..=per_step {
                    let t = at(s);
                    let x = DVec3::new(
                        dense.eval_component(6 * k, t),
                        dense.eval_component(6 * k + 1, t),
                        dense.eval_component(6 * k + 2, t),
                    );
                    paths[i].push((t, x));
                }
            }
            go_on()
        },
    )?;
    // Up to each particle's end (ghosts are followed a little further).
    for (path, traj) in paths.iter_mut().zip(&run.trajectories) {
        path.retain(|(t, _)| *t < traj.end.t);
        path.push((traj.end.t, traj.end.x));
    }
    for ((line, traj), fade) in lines.iter_mut().zip(&run.trajectories).zip(&run.fades) {
        if traj.outcome != Outcome::LeftBounds {
            line.retain(|l| l.0 <= traj.end.t);
        }
        // Into a screening cup it flies on uniformly from where it entered (the world
        // line's continuation after its last sample).
        if let Some(f) = fade
            && line.last().is_some_and(|l| l.0 < f.t_off)
        {
            let a = line.last().map_or(DVec3::ZERO, |l| l.3);
            line.push((f.t_off, f.x_off, f.v, a));
        }
    }
    let ends = run
        .trajectories
        .iter()
        .map(|t| {
            (!matches!(t.outcome, Outcome::Timeout | Outcome::LeftBounds))
                .then_some((t.end.t, t.end.x))
        })
        .collect();
    let preview = BeamPreview {
        paths,
        shots: shots.to_vec(),
        outcomes: run.trajectories.iter().map(|t| t.outcome).collect(),
        energy_rel_error: run.energy_max_rel_error,
        worldlines: lines,
        ends,
        fades: run.fades.clone(),
        energy: run.energy.clone(),
        retardation_max: run
            .neglected_retardation
            .iter()
            .fold(0.0, |m: f64, &x| m.max(x)),
        radiated_max: run
            .trajectories
            .iter()
            .map(|t| t.radiated_energy / t.kinetic_initial)
            .fold(0.0, f64::max),
        retarded: scn.retarded && scn.c.is_finite(),
    };
    Some((run, preview))
}

/// Preview flight with its dense path; `None` if `go_on` stopped it.
fn preview_shot(
    scn: &Scenario<LevelField>,
    tol: f64,
    go_on: impl Fn() -> bool,
) -> Option<(Trajectory, Preview)> {
    let mut dense_points: Vec<(f64, DVec3, DVec3, f64)> = Vec::new();
    let traj = run_cancellable(scn, &RunSettings::with_tolerance(tol), |step| {
        let (a, b) = (step.t_start(), step.t_end());
        for i in 1..=8 {
            let t = a + (b - a) * f64::from(i) / 8.0;
            let (x, p) = step.state(t);
            dense_points.push((t, x, p, -step.radiation_work(t)));
        }
        go_on()
    })?;
    let path = build_path(
        scn,
        &traj.samples[0],
        &dense_points,
        &traj.end,
        -traj.radiation_work,
    );
    let max_speed_over_c = path.iter().map(|p| p.speed_over_c).fold(0.0, f64::max);
    let image_force_bound = level::electrode_image_force_bound(
        scn,
        traj.kinetic_initial,
        path.iter().map(|p| (p.x, p.t)),
    );
    // Charge radiation plus the magnetic moment's (estimate), both neglected.
    let moment_radiation =
        level::moment_radiation_estimate(scn, path.iter().map(|p| (p.x, p.p, p.t)));
    let spectrum_top = match scn.acceptance.and_then(|a| a.radiation) {
        Some(w) if w.band.is_some() => 2.0 * w.top_frequency(),
        Some(_) => physics::spectrum::display_range(&traj.emission, scn.c),
        None => 0.0,
    };
    let preview = Preview {
        path,
        outcome: traj.outcome,
        flight_time: traj.end.t,
        energy_rel_error: traj.energy_max_abs_error / traj.kinetic_initial,
        radiated_fraction: (traj.radiated_energy + moment_radiation) / traj.kinetic_initial,
        max_speed_over_c,
        radiation_reaction: scn.radiation_reaction && scn.c.is_finite(),
        radiation_loss_fraction: -traj.radiation_work / traj.kinetic_initial,
        reaction_ratio_max: traj.reaction_ratio_max,
        image_force_bound,
        // Measured for the verdict on arrival; otherwise for display, on the flight so far.
        radiation: traj.radiation.or_else(|| {
            let w = scn.acceptance.and_then(|a| a.radiation)?;
            (!traj.emission.is_empty())
                .then(|| w.measure(&traj.emission, scn.particle.charge, scn.c))
        }),
        spectrum: match scn.acceptance.and_then(|a| a.radiation) {
            Some(w) if !traj.emission.is_empty() => physics::spectrum::arc_spectrum(
                &w,
                &traj.emission,
                scn.particle.charge,
                scn.c,
                spectrum_top,
                160,
            ),
            _ => Vec::new(),
        },
        spectrum_top,
    };
    Some((traj, preview))
}

fn build_path(
    scn: &Scenario<LevelField>,
    start: &physics::trajectory::Sample,
    dense: &[(f64, DVec3, DVec3, f64)],
    end: &physics::trajectory::Sample,
    radiated_end: f64,
) -> Vec<PathPoint> {
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let (q, moment) = (scn.particle.charge, scn.particle.moment);
    let at_a = scn.field.sample(scn.x0, 0.0);
    let (phi_a, bz_a) = (at_a.phi, at_a.b.z);
    let point = |t: f64, x: DVec3, p: DVec3, radiated: f64| {
        let f = scn.field.sample(x, t);
        let v = kin.velocity(p);
        PathPoint {
            t,
            x,
            p,
            kinetic: kin.kinetic_energy(p),
            potential: q * (f.phi - phi_a) - moment * (f.b.z - bz_a),
            radiated,
            force: (f.e + v.cross(f.b)) * q + scn.field.grad_bz(x, t) * moment,
            speed_over_c: if scn.c.is_finite() {
                v.length() / scn.c
            } else {
                0.0
            },
        }
    };
    let mut path = vec![point(start.t, start.x, start.p, 0.0)];
    path.extend(
        dense
            .iter()
            .filter(|(t, _, _, _)| *t < end.t)
            .map(|&(t, x, p, r)| point(t, x, p, r)),
    );
    path.push(point(end.t, end.x, end.p, radiated_end));
    path
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cost(revision: u64) -> Response {
        Response::BeamVerified {
            revision,
            flight: 0,
            results: Vec::new(),
            run: None,
        }
    }

    /// At finite c the verdict of an interacting beam flies the exact retarded
    /// interaction, and the worker sends that flight for the views: flown retarded (the
    /// preview quasi-static), with the verdict's outcomes, a world line for every particle
    /// and a screening cup's fade for every one that entered its detector.
    #[test]
    fn the_verdicts_exact_flight_is_sent_for_the_views() {
        let path = crate::levels_dir().join("55_soft_landing_current.json");
        let text = std::fs::read_to_string(path).expect("level 55");
        let level = level::Level::from_json(&text).expect("valid level");
        let (tx, rx) = channel();
        let sink = Sink {
            tx: &tx,
            record: std::cell::RefCell::new(Vec::new()),
        };
        let newest = Arc::new(AtomicU64::new(1));
        let req = Request {
            revision: 1,
            placement: level.reference_solution.clone(),
            level,
        };
        assert!(beam_request(&req, &sink, &newest));
        let responses: Vec<Response> = rx.try_iter().collect();
        let preview = responses
            .iter()
            .find_map(|r| match r {
                Response::BeamPreview { preview, .. } => Some(preview),
                _ => None,
            })
            .expect("a preview");
        let (results, run) = responses
            .iter()
            .find_map(|r| match r {
                Response::BeamVerified { results, run, .. } => Some((results, run)),
                _ => None,
            })
            .expect("a verdict");
        let run = run.as_ref().expect("the verdict's flight");
        assert!(run.retarded && !preview.retarded);
        assert_eq!(run.outcomes.len(), results.len());
        for (i, (_, outcome)) in results.iter().enumerate() {
            assert_eq!(run.outcomes[i], *outcome, "particle {i}");
            assert!(!run.worldlines[i].is_empty(), "particle {i}");
            assert_eq!(
                run.fades[i].is_some(),
                *outcome == Outcome::Arrived,
                "particle {i}"
            );
        }
    }

    /// The cache returns what was put under a key, most recent first, forgets the
    /// oldest beyond its size, and replayed responses carry the new revision.
    #[test]
    fn result_cache_replays_and_evicts() {
        let mut c = ResultCache {
            entries: std::collections::VecDeque::new(),
        };
        for k in 0..(CACHE_SETUPS as u64 + 3) {
            c.put(k, vec![cost(k)]);
        }
        assert!(c.get(0).is_none(), "oldest evicted");
        let r = c.get(5).expect("kept");
        assert!(matches!(
            retag(r[0].clone(), 99),
            Response::BeamVerified { revision: 99, .. }
        ));
        assert_eq!(c.entries.front().map(|e| e.0), Some(5), "moved to front");
    }
}
