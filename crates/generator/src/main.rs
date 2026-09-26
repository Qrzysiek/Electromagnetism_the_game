//! Offline level tools: solving and checking levels. (Random level generation: M5.)

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use level::Level;
use level::solve as search;
use physics::verify::verify;

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
    /// Rewrite level files in the current format (older formats are migrated on load).
    Normalize { paths: Vec<PathBuf> },
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
        Command::Check { path } => {
            let level = load(&path);
            let placement = &level.reference_solution;
            println!(
                "{}: placement {:?}",
                level.name,
                level.check_placement(placement)
            );
            for (i, scn) in level.scenarios(placement).iter().enumerate() {
                let v = verify(scn, level.tolerances());
                let end = v.verified.end.x;
                println!(
                    "  shot {}: outcome {:?}, status {:?}, flight time {:.3}, ends at ({:.4}, {:.4})",
                    i + 1,
                    v.outcome(),
                    v.status,
                    v.verified.end.t,
                    end.x,
                    end.y
                );
            }
        }
    }
}
