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
    /// Print the reference flight of one flight index: t, x, y, |p|, radiated energy.
    Trace {
        path: PathBuf,
        #[arg(long, default_value_t = 0)]
        flight: usize,
        /// Print every n-th accepted step.
        #[arg(long, default_value_t = 20)]
        every: usize,
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
            if let Some(spec) = level.tube {
                // Tube levels: the verified current against the goal (printed like a
                // shot, so the tools check it too).
                let v = level.tube_verdict(placement);
                let g = spec.goal;
                println!(
                    "  shot tube: current {:.6} ± {:.6} (preview {:.6}), goal {} to {}: {:?}: {}",
                    v.fine,
                    v.error,
                    v.preview,
                    g.min,
                    g.max,
                    v.status,
                    if v.status == level::tube::TubeStatus::Met {
                        "outcome Arrived, status Verified"
                    } else {
                        "outcome Rejected, status NotVerified"
                    }
                );
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
