//! Offline level tools: solving and checking levels. (Random level generation: M5.)

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use level::Level;
use level::solve as search;

#[derive(Parser)]
#[command(about = "Level tools for Electromagnetism the game")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search for verified solutions of a level.
    Solve {
        path: PathBuf,
        /// Annealing restarts per charge count.
        #[arg(long, default_value_t = 64)]
        restarts: u64,
        /// Iterations per restart.
        #[arg(long, default_value_t = 400)]
        iterations: u32,
        /// Store the first solution with the fewest charges as the reference solution.
        #[arg(long)]
        write: bool,
    },
    /// Verify the reference solution of a level.
    Check { path: PathBuf },
    /// Measure difficulty: configuration space, random-guess rate, heuristic search.
    Analyze {
        paths: Vec<PathBuf>,
        /// Random placements sampled.
        #[arg(long, default_value_t = 2000)]
        samples: usize,
        /// Heuristic search runs.
        #[arg(long, default_value_t = 32)]
        runs: usize,
        /// Objective evaluations per search run.
        #[arg(long, default_value_t = 400)]
        budget: u32,
        /// Also report the fewest elements that solve the level: the reference solution's
        /// count, unless a search with fewer elements succeeds (a heuristic: a failed
        /// search is not a proof).
        #[arg(long)]
        fewest: bool,
    },
    /// Rewrite level files in the current format (older formats are migrated on load).
    Normalize { paths: Vec<PathBuf> },
    /// For each power supply and plate of the reference solution: scan its potential over
    /// the slider's range (the others as in the reference) and print the intervals that
    /// solve the level (preview tolerance; each interval's middle verified).
    Window {
        path: PathBuf,
        /// Potentials scanned across the range.
        #[arg(long, default_value_t = 240)]
        steps: u32,
    },
    /// Tube levels: for each placement in a JSON file (a list of placements, each a
    /// list of elements), the preview's currents into every electrode over the goal's
    /// window (per emitter, positive for arriving carriers), what left the board, and the
    /// projectiles' arrivals. For designing tube levels.
    TubeScan {
        path: PathBuf,
        placements: PathBuf,
        /// The resolution: 1 the preview, 2 the verification, higher for convergence
        /// checks (segments, step and macroparticle weight divided by it).
        #[arg(long, default_value_t = 1)]
        refine: u32,
        /// Convergence studies: this many times as many macroparticles, alone.
        #[arg(long, default_value_t = 1.0)]
        particles: f64,
        /// Convergence studies: this many times the particles' softening, alone.
        #[arg(long, default_value_t = 1.0)]
        softening: f64,
        /// Convergence studies: the emitters' segments this many times finer, alone.
        #[arg(long, default_value_t = 1.0)]
        emitter: f64,
        /// Also print this electrode's (1-based) potential every step: `trace t V`
        /// (a tank's oscillation, for designing).
        #[arg(long)]
        trace: Option<usize>,
    },
    /// Tube levels: the electromagnetic field of the reference's preview run
    /// (`level::tube_em`), every `every` time units to `--end`: the power radiated out of
    /// the board (over the last 2 time units) and, for a level with an output, its load's
    /// `R⟨I²⟩` over the same span.
    TubeField {
        path: PathBuf,
        #[arg(long, default_value_t = 10.0)]
        every: f64,
        /// The run's end (default: the level's display end).
        #[arg(long)]
        end: Option<f64>,
    },
    /// Print the reference flight of one flight index: t, x, y, |p|, radiated energy.
    Trace {
        path: PathBuf,
        #[arg(long, default_value_t = 0)]
        flight: usize,
        /// Print every n-th accepted step.
        #[arg(long, default_value_t = 20)]
        every: usize,
    },
    /// For each placement in a JSON file (a list of placements, each a list of elements;
    /// without a file, the reference solution): the verdict of every flight, beam and tube
    /// goal, and each flight's end, as JSON at full precision. The game's numbers for the
    /// level documents (`scripts/level_docs.py game`, docs/levels/README.md).
    Scan {
        path: PathBuf,
        placements: Option<PathBuf>,
        /// Also list every accepted step of each verified flight (t, x, p).
        #[arg(long)]
        samples: bool,
    },
}

fn load(path: &PathBuf) -> Level {
    let s = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    Level::from_json(&s).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The most elements a player may place: charges, magnets, antennas and plates within
/// the limits, and a power supply on every tunable electrode.
fn most_elements(level: &Level) -> usize {
    let l = &level.limits;
    let supplies = level.electrodes.iter().filter(|e| e.tunable).count();
    (l.max_charges + l.max_magnets + l.max_antennas + l.max_plates) as usize + supplies
}

/// `generator scan`: the verdicts of several placements, as JSON.
fn scan(level: &Level, placements: Option<&PathBuf>, samples: bool) {
    let list: Vec<Vec<level::Element>> = placements.map_or_else(
        || vec![level.reference_solution.clone()],
        |p| {
            let text = std::fs::read_to_string(p).expect("placements file");
            serde_json::from_str(&text).expect("placements")
        },
    );
    let out: Vec<serde_json::Value> = list
        .iter()
        .map(|placement| scan_one(level, placement, samples))
        .collect();
    println!("{}", serde_json::to_string_pretty(&out).expect("JSON"));
}

/// One placement's verdicts (`scan`).
fn scan_one(level: &Level, placement: &[level::Element], samples: bool) -> serde_json::Value {
    use serde_json::json;
    let mut refused = level.setup_issues(placement);
    if let Err(e) = level.check_placement(placement) {
        refused.push(format!("{e:?}"));
    }
    if !refused.is_empty() {
        return json!({ "placement": placement, "refused": refused });
    }
    let state = |s: &physics::trajectory::Sample| json!({ "t": s.t, "x": [s.x.x, s.x.y, s.x.z], "p": [s.p.x, s.p.y, s.p.z] });
    let flight = |t: &physics::trajectory::Trajectory| {
        let mut f = json!({
            "outcome": format!("{:?}", t.outcome),
            "end": state(&t.end),
            "kinetic_initial": t.kinetic_initial,
            "energy_max_abs_error": t.energy_max_abs_error,
            "radiated_energy": t.radiated_energy,
            "steps": t.stats.n_step,
        });
        if let Some(r) = t.radiation {
            f["radiation_per_sr"] = json!(r);
        }
        if samples {
            f["samples"] = t.samples.iter().map(state).collect();
        }
        f
    };
    let mut out = json!({ "placement": placement });
    if let Some(spec) = &level.tube {
        let v = level.tube_verdict(placement);
        let goal = |c: &level::tube::CurrentVerdict, min: f64, max: f64| {
            json!({ "value": c.fine, "error": c.error, "preview": c.preview,
                    "min": min, "max": max, "status": format!("{:?}", c.status) })
        };
        out["tube"] = json!({
            "currents": v.currents.iter().zip(&spec.goals)
                .map(|(c, g)| goal(c, g.min, g.max)).collect::<Vec<_>>(),
            "potentials": v.potential.iter().zip(&spec.potential_goals)
                .map(|(c, g)| goal(c, g.min, g.max)).collect::<Vec<_>>(),
            "arrivals": v.arrivals.iter().map(|(j, a, b, s)| json!({
                "projectile": j + 1, "preview": a, "verification": b,
                "status": format!("{s:?}") })).collect::<Vec<_>>(),
            "capped": v.capped,
            "status": format!("{:?}", v.status),
        });
    }
    if level.has_beams() {
        let beams: Vec<serde_json::Value> = level
            .verify_beams(placement)
            .iter()
            .enumerate()
            .map(|(d, v)| {
                let shots: Vec<serde_json::Value> = (0..level.shots.len())
                    .map(|s| {
                        let (ok, n) = v.transmitted(s);
                        json!({ "shot": s + 1, "arrived_verified": ok, "particles": n })
                    })
                    .collect();
                let particles: Vec<serde_json::Value> = v
                    .verified
                    .trajectories
                    .iter()
                    .zip(&v.status)
                    .zip(&v.shot)
                    .map(|((t, st), s)| {
                        let mut f = flight(t);
                        f["shot"] = json!(s + 1);
                        f["status"] = json!(format!("{st:?}"));
                        f
                    })
                    .collect();
                json!({ "disturbance": d + 1, "shots": shots, "particles": particles,
                        "energy_max_rel_error": v.verified.energy_max_rel_error })
            })
            .collect();
        out["beams"] = json!(beams);
    } else if !level.shots.is_empty() {
        let flights: Vec<serde_json::Value> = level
            .verify_flights(placement)
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let (shot, d) = level.flight_of(i);
                json!({
                    "shot": shot + 1,
                    "disturbance": d + 1,
                    "outcome": format!("{:?}", v.outcome()),
                    "status": format!("{:?}", v.status),
                    "verified": flight(&v.verified),
                    "preview": flight(&v.preview),
                })
            })
            .collect();
        out["flights"] = json!(flights);
    }
    out
}

/// `generator tube-scan`: currents into every electrode for several placements.
fn tube_scan(
    level: &Level,
    placements: &PathBuf,
    res: level::tube::TubeResolution,
    trace: Option<usize>,
) {
    let text = std::fs::read_to_string(placements).expect("placements file");
    let list: Vec<Vec<level::Element>> = serde_json::from_str(&text).expect("placements");
    let Some(spec) = &level.tube else {
        eprintln!("not a tube level");
        return;
    };
    // The electrodes' currents over the first current goal's window.
    let (start, end) = spec
        .goals
        .first()
        .map_or((0.0, level.tube_verify_end()), |g| (g.start, g.end));
    for (i, placement) in list.iter().enumerate() {
        let issues = level.setup_issues(placement);
        if !issues.is_empty() {
            println!("{i}: refused: {}", issues.join(" "));
            continue;
        }
        let mut sim = level.tube_sim_at(
            placement,
            level::tube::TubeResolution {
                refine: res.refine.max(1),
                ..res
            },
        );
        let end = end.max(level.tube_verify_end());
        let step = spec.step / f64::from(res.refine.max(1));
        let advance = |sim: &mut level::tube::TubeSim, to: f64| match trace {
            Some(e) => {
                while sim.t() < to - 0.5 * step {
                    sim.advance_to(sim.t() + step);
                    println!("trace {:.6} {:.9}", sim.t(), sim.potentials()[e - 1]);
                }
            }
            None => sim.advance_to(to),
        };
        advance(&mut sim, start);
        let (t0, q0, lost0) = (sim.t(), sim.collected().to_vec(), sim.lost());
        advance(&mut sim, end);
        let dt = (sim.t() - t0).max(1e-300);
        let signs = [spec.charge_per_mass.signum(), 1.0];
        let currents: Vec<String> = sim
            .collected()
            .iter()
            .zip(&q0)
            .enumerate()
            .map(|(e, (q, q0))| {
                let per: Vec<String> = q
                    .iter()
                    .zip(q0)
                    .zip(signs)
                    .map(|((a, b), s)| format!("{:.4}", s * (a - b) / dt))
                    .collect();
                format!("e{}={}", e + 1, per.join("/"))
            })
            .collect();
        let fates: Vec<String> = sim
            .fates()
            .iter()
            .map(|f| match f {
                Some(physics::tube::Fate::Absorbed { electrode, t, x }) => {
                    format!("e{} t{t:.1} ({:.2},{:.2})", electrode + 1, x.x, x.y)
                }
                Some(physics::tube::Fate::Lost { t, x }) => {
                    format!("out t{t:.1} ({:.2},{:.2})", x.x, x.y)
                }
                None => "flying".into(),
            })
            .collect();
        println!(
            "{i}: {} lost {:.4} goals {:?} potential {:?} arrivals {:?} fates {} in flight {} weights {:?}",
            currents.join(" "),
            (sim.lost() - lost0) / dt,
            sim.currents(),
            sim.results().potential,
            sim.arrivals(),
            fates.join("; "),
            sim.in_flight(),
            sim.weights()
        );
    }
}

/// `generator tube-field`: the radiated power of the reference's run along it.
fn tube_field(level: &Level, every: f64, end: Option<f64>) {
    let end = end.unwrap_or_else(|| level.tube_display_end());
    let started = std::time::Instant::now();
    let run = level.tube_run(&level.reference_solution, 1, end);
    let ran = started.elapsed().as_secs_f64();
    let mut em = level::tube_em::TubeEm::new(level, &run);
    let g = em.fdtd.grid;
    println!(
        "run {ran:.1} s; field grid {} x {} (h {}), c {}, dt {:.5}, shown {:.3} later \
         (a closed box's smoothing)",
        g.nx,
        g.ny,
        g.h,
        em.fdtd.c,
        em.fdtd.dt,
        em.delay()
    );
    let spec = level.tube.as_ref().expect("a tube level");
    let load = spec.output.and_then(|o| {
        spec.circuits
            .get(o.circuit)
            .and_then(|c| c.resistance)
            .map(|r| (o.circuit, r))
    });
    let started = std::time::Instant::now();
    let mut t = every;
    while t <= end + 1e-9 {
        em.show(&run, t);
        let absorbed = load.map_or(0.0, |(k, r)| {
            let i2: Vec<f64> = run
                .frames
                .iter()
                .filter(|f| f.t <= t && f.t >= t - 2.0)
                .filter_map(|f| f.currents.get(k).map(|i| f64::from(*i).powi(2)))
                .collect();
            #[allow(clippy::cast_precision_loss)]
            let mean = i2.iter().sum::<f64>() / i2.len().max(1) as f64;
            r * mean
        });
        // The field's energy and its magnetic part, the field's own and the shown one's (a
        // closed box keeps everything radiated in it; its smoothing shows the run's
        // currents' field).
        let magnetic = |bz: &[f64]| {
            bz.iter().map(|b| b * b).sum::<f64>() * g.h * g.h / (8.0 * std::f64::consts::PI)
        };
        println!(
            "t {:7.2}: radiated {:.5e}, load R<I^2> {:.5e}, field energy {:.5e} (magnetic \
             {:.5e}; shown {:.5e})",
            em.shown().t,
            em.mean_outflow(2.0),
            absorbed,
            em.fdtd.energy(),
            magnetic(&em.fdtd.bz),
            magnetic(em.shown().bz)
        );
        t += every;
    }
    println!("field {:.1} s", started.elapsed().as_secs_f64());
}

/// `generator window`: the solving intervals of each slider-set potential.
#[allow(clippy::cast_precision_loss)]
fn window(level: &Level, steps: u32) {
    use level::ElementKind;
    let reference = &level.reference_solution;
    for (i, e) in reference.iter().enumerate() {
        // A supply's own slider: an electrode's potentials or a coil's ramp rates.
        let targets = level.supply_targets();
        let list: &[f64] = match e.kind {
            ElementKind::Supply => targets
                .iter()
                .find(|(c, _)| *c == e.node)
                .map_or(&[][..], |&(_, t)| level.supply_list(t)),
            ElementKind::Plate => &level.limits.plate_voltages,
            _ => continue,
        };
        let Some((lo, hi)) = level::value_range(list) else {
            continue;
        };
        let mut trial = reference.clone();
        let mut solves = |v: f64| {
            trial[i].value = v;
            search::objective(level, &trial).1 == physics::trajectory::Outcome::Arrived
        };
        let at = |k: u32| lo + (hi - lo) * f64::from(k) / f64::from(steps);
        let mut intervals: Vec<(f64, f64)> = Vec::new();
        let mut start: Option<f64> = None;
        for k in 0..=steps {
            let v = at(k);
            match (solves(v), start) {
                (true, None) => start = Some(v),
                (false, Some(s)) => {
                    intervals.push((s, at(k - 1)));
                    start = None;
                }
                _ => {}
            }
        }
        if let Some(s) = start {
            intervals.push((s, hi));
        }
        println!(
            "{}: {:?} at {:?}, reference {:.4e}, range [{lo:.4e}, {hi:.4e}]",
            level.name, e.kind, e.node, e.value
        );
        for (a, b) in intervals {
            let mut mid = reference.clone();
            mid[i].value = 0.5 * (a + b);
            println!(
                "  solves on [{a:.4e}, {b:.4e}] ({:.1} % of the range; middle verified: {})",
                100.0 * (b - a) / (hi - lo),
                search::is_verified_solution(level, &mid)
            );
        }
    }
}

/// The fewest player elements that solve the level: the reference solution's count, unless
/// the search finds a verified solution with fewer (exhaustive for one element, annealing
/// with 64 restarts of 600 iterations for more).
fn fewest_elements(level: &Level) -> usize {
    let reference = level.reference_solution.len();
    if reference > 1 && !search::single_element_solutions(level).is_empty() {
        return 1;
    }
    (2..reference)
        .find(|&k| !search::anneal(level, k, 64, 600, 0x5EED + k as u64).is_empty())
        .unwrap_or(reference)
}

fn main() {
    match Cli::parse().command {
        Command::Solve {
            path,
            restarts,
            iterations,
            write,
        } => {
            let mut level = load(&path);
            println!("level '{}'", level.name);
            let (no_charges, outcome) = search::objective(&level, &[]);
            println!(
                "  without player charges: {outcome:?} (closest approach to detector {no_charges:.3})"
            );

            let mut best: Option<Vec<level::Element>> = None;
            let singles = search::single_element_solutions(&level);
            println!("  verified 1-charge solutions: {}", singles.len());
            if singles.len() <= 10 {
                for c in &singles {
                    println!("    {:?} {:?} {:e}", c.node, c.kind, c.value);
                }
            }
            if let Some(c) = singles.first() {
                best = Some(vec![*c]);
            }
            for k in 2..=most_elements(&level) {
                let found = search::anneal(&level, k, restarts, iterations, 0x5EED + k as u64);
                println!(
                    "  verified {k}-charge solutions found by annealing: {}",
                    found.len()
                );
                if best.is_none() {
                    best = found.into_iter().next();
                }
            }
            match best {
                Some(sol) if write => {
                    println!("  reference solution: {sol:?}");
                    level.reference_solution = sol;
                    std::fs::write(&path, level.to_json() + "\n").expect("write level");
                }
                Some(sol) => println!("  best solution: {sol:?}"),
                None => println!("  no solution found"),
            }
        }
        Command::Analyze {
            paths,
            samples,
            runs,
            budget,
            fewest,
        } => {
            println!(
                "| level | log10 configs | random solve rate | search success | mean evals | expected effort: search / guessing | smoothness |{}",
                if fewest {
                    " fewest elements found |"
                } else {
                    ""
                }
            );
            println!(
                "|---|---|---|---|---|---|---|{}",
                if fewest { "---|" } else { "" }
            );
            for path in paths {
                let level = load(&path);
                let a = level::analysis::analyze(&level, samples, runs, budget, 0xD1FF);
                let rate = if a.random_solutions == 0 {
                    format!("< {:.1e}", a.random_rate())
                } else {
                    format!("{:.1e}", a.random_rate())
                };
                let fewest = if fewest {
                    format!(" {} |", fewest_elements(&level))
                } else {
                    String::new()
                };
                println!(
                    "| {} | {:.1} | {rate} | {}/{} | {:.0} | {:.0} / {:.0} | {:.2} |{fewest}",
                    path.file_stem().unwrap().to_string_lossy(),
                    a.config_space_log10,
                    a.search_successes,
                    a.search_runs,
                    a.search_mean_evaluations,
                    a.expected_search_effort(),
                    a.expected_guess_effort(),
                    a.smoothness
                );
            }
        }
        Command::Normalize { paths } => {
            for path in paths {
                let level = load(&path);
                std::fs::write(
                    &path,
                    level.to_json()
                        + "
",
                )
                .expect("write level");
                println!("{}: format {}", path.display(), level.format_version);
            }
        }
        Command::Window { path, steps } => window(&load(&path), steps),
        Command::TubeScan {
            path,
            placements,
            refine,
            particles,
            softening,
            emitter,
            trace,
        } => tube_scan(
            &load(&path),
            &placements,
            level::tube::TubeResolution {
                refine,
                particles,
                softening,
                emitter,
            },
            trace,
        ),
        Command::TubeField { path, every, end } => tube_field(&load(&path), every, end),
        Command::Scan {
            path,
            placements,
            samples,
        } => scan(&load(&path), placements.as_ref(), samples),
        Command::Trace {
            path,
            flight,
            every,
        } => {
            let level = load(&path);
            let scn = &level.scenarios(&level.reference_solution)[flight];
            let tr = physics::trajectory::run(
                scn,
                &physics::trajectory::RunSettings::with_tolerance(level.physics.tolerances.preview),
            );
            for s in tr.samples.iter().step_by(every.max(1)) {
                println!(
                    "{:10.3} {:9.4} {:9.4} {:.6e}",
                    s.t,
                    s.x.x,
                    s.x.y,
                    s.p.length()
                );
            }
            println!(
                "end {:?} t = {:.3}, radiated (LL work) {:.4e}, max |F_RR|/|F_L| {:.2e}",
                tr.outcome, tr.end.t, -tr.radiation_work, tr.reaction_ratio_max
            );
        }
        Command::Check { path } => {
            let level = load(&path);
            let placement = &level.reference_solution;
            println!(
                "{}: placement {:?}",
                level.name,
                level.check_placement(placement)
            );
            if let Some(spec) = &level.tube {
                // Tube levels: the verified current against the goal and the projectiles'
                // arrivals (printed like shots, so the tools check them too).
                let v = level.tube_verdict(placement);
                let verdict = |s: level::tube::TubeStatus| {
                    if s == level::tube::TubeStatus::Met {
                        "outcome Arrived, status Verified"
                    } else {
                        "outcome Rejected, status NotVerified"
                    }
                };
                for (c, g) in v.currents.iter().zip(&spec.goals) {
                    println!(
                        "  shot tube: {:?} current {:.6} ± {:.6} (preview {:.6}), goal {} to {}: {:?}: {}",
                        g.carriers,
                        c.fine,
                        c.error,
                        c.preview,
                        g.min,
                        g.max,
                        c.status,
                        verdict(c.status)
                    );
                }
                for (c, g) in v.potential.iter().zip(&spec.potential_goals) {
                    println!(
                        "  shot potential: {:?} {:.6} ± {:.6} (preview {:.6}), goal {} to {}: {:?}: {}",
                        g.measure,
                        c.fine,
                        c.error,
                        c.preview,
                        g.min,
                        g.max,
                        c.status,
                        verdict(c.status)
                    );
                }
                for (j, a, b, s) in &v.arrivals {
                    println!(
                        "  shot projectile {}: arrives at {a:?} (preview), {b:?} (verification): {s:?}: {}",
                        j + 1,
                        verdict(*s)
                    );
                }
            }
            if level.has_beams() {
                // Beam levels: the verified transmission of every beam shot per flight.
                for (d, v) in level.verify_beams(placement).iter().enumerate() {
                    for (s, shot) in level.shots.iter().enumerate() {
                        let (ok, n) = v.transmitted(s);
                        let need = shot.beam.map_or(1.0, |b| b.transmission);
                        #[allow(clippy::cast_precision_loss)]
                        let enough = ok as f64 >= need * n as f64 - 1e-9;
                        println!(
                            "  shot {} disturbance {}: beam {ok}/{n} arrive verified (need {:.0} %): {}",
                            s + 1,
                            d + 1,
                            need * 100.0,
                            if enough {
                                "outcome Arrived, status Verified"
                            } else {
                                "outcome Short, status NotVerified"
                            }
                        );
                    }
                    // Every particle's end (EM_BEAM_ENDS=1, for comparing models).
                    if std::env::var_os("EM_BEAM_ENDS").is_some() {
                        for (i, t) in v.verified.trajectories.iter().enumerate() {
                            println!(
                                "  end {} disturbance {}: {:?} {:.12} {:.12} t {:.9}",
                                i,
                                d + 1,
                                t.outcome,
                                t.end.x.x,
                                t.end.x.y,
                                t.end.t
                            );
                        }
                        for (i, t) in v.preview.trajectories.iter().enumerate() {
                            println!(
                                "  preview end {} disturbance {}: {:?} {:.12} {:.12} t {:.9}",
                                i,
                                d + 1,
                                t.outcome,
                                t.end.x.x,
                                t.end.x.y,
                                t.end.t
                            );
                        }
                    }
                    // Goal particles (the level's free particles with a detector): each
                    // must arrive, verified (printed like a shot, so tools check them too).
                    for (k, &i) in level.goal_particles().iter().enumerate() {
                        let t = &v.verified.trajectories[i];
                        println!(
                            "  shot goal-{} disturbance {}: free particle, outcome {:?}, status {:?}, ends at ({:.4}, {:.4})",
                            k + 1,
                            d + 1,
                            t.outcome,
                            v.status[i],
                            t.end.x.x,
                            t.end.x.y
                        );
                    }
                    // Radiation goals: every charge's radiation together (PHYSICS.md §3.4).
                    for (i, t) in v.verified.trajectories.iter().enumerate() {
                        if let Some(e) = t.radiation {
                            println!(
                                "  particle {i} disturbance {}: radiation into its goal {e:.6e} per sr (preview {:.6e})",
                                d + 1,
                                v.preview.trajectories[i].radiation.unwrap_or(f64::NAN)
                            );
                        }
                    }
                    println!(
                        "  beam energy drift {:.1e}, {} steps",
                        v.verified.energy_max_rel_error, v.verified.stats.n_step
                    );
                    // Where each shot's particles end (for designing detectors).
                    let shots = level.beam_shots(&level.reference_solution);
                    for s in 0..level.shots.len() {
                        let ends: Vec<_> = v
                            .verified
                            .trajectories
                            .iter()
                            .zip(&shots)
                            .filter(|(_, k)| **k == s)
                            .map(|(t, _)| t.end.x)
                            .collect();
                        let lo = ends.iter().fold(ends[0], |a, b| a.min(*b));
                        let hi = ends.iter().fold(ends[0], |a, b| a.max(*b));
                        println!(
                            "  beam {} ends within x {:.2}..{:.2}, y {:.2}..{:.2}",
                            s + 1,
                            lo.x,
                            hi.x,
                            lo.y,
                            hi.y
                        );
                    }
                }
            }
            for (i, v) in level
                .verify_flights(placement)
                .iter()
                .enumerate()
                .filter(|_| !level.has_beams())
            {
                let end = v.verified.end.x;
                let (shot, d) = level.flight_of(i);
                let flight = if level.disturbances.is_empty() {
                    String::new()
                } else {
                    format!(" disturbance {}", d + 1)
                };
                println!(
                    "  shot {}{flight}: outcome {:?}, status {:?}, flight time {:.3}, ends at ({:.4}, {:.4})",
                    shot + 1,
                    v.outcome(),
                    v.status,
                    v.verified.end.t,
                    end.x,
                    end.y
                );
                if let Some(e) = v.verified.radiation {
                    println!(
                        "  flight {}{flight} radiation into its goal: {e:.6e} per sr (preview {:.6e})",
                        shot + 1,
                        v.preview.radiation.unwrap_or(f64::NAN)
                    );
                }
            }
            let cost = level.measure_cost(placement);
            let meters: Vec<String> = cost
                .meters()
                .iter()
                .map(|m| format!("{} {} ({:?})", m.label, m.text(), m.load()))
                .collect();
            println!("  cost: {}", meters.join(", "));
        }
    }
}
