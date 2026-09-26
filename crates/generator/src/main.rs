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
    },
    /// Rewrite level files in the current format (older formats are migrated on load).
    Normalize { paths: Vec<PathBuf> },
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
            for k in 2..=(level.limits.max_charges + level.limits.max_magnets) as usize {
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
        } => {
            println!(
                "| level | log10 configs | random solve rate | search success | mean evals | expected effort: search / guessing | smoothness |"
            );
            println!("|---|---|---|---|---|---|---|");
            for path in paths {
                let level = load(&path);
                let a = level::analysis::analyze(&level, samples, runs, budget, 0xD1FF);
                let rate = if a.random_solutions == 0 {
                    format!("< {:.1e}", a.random_rate())
                } else {
                    format!("{:.1e}", a.random_rate())
                };
                println!(
                    "| {} | {:.1} | {rate} | {}/{} | {:.0} | {:.0} / {:.0} | {:.2} |",
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
            for (i, v) in level.verify_flights(placement).iter().enumerate() {
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
