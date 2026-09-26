//! Checks on the shipped levels.

use std::path::PathBuf;

use level::Level;
use physics::trajectory::{Outcome, RunSettings, run};

fn shipped_levels() -> Vec<(String, Level)> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut out: Vec<(String, Level)> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| {
            p.extension().is_some_and(|e| e == "json") && !p.ends_with("golden_hashes.json")
        })
        .map(|p| {
            let name = p.file_stem().unwrap().to_string_lossy().into_owned();
            (
                name,
                Level::from_json(&std::fs::read_to_string(&p).unwrap()).unwrap(),
            )
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[test]
fn reference_solutions_are_verified() {
    for (name, level) in shipped_levels() {
        assert_eq!(
            level.check_placement(&level.reference_solution),
            Ok(()),
            "{name}"
        );
        for (i, v) in level
            .verify_flights(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let (shot, d) = level.flight_of(i);
            let flight = format!("{name} shot {} disturbance {}", shot + 1, d + 1);
            assert_eq!(v.outcome(), Outcome::Arrived, "{flight}");
            assert!(v.status.is_verified(), "{flight}: {:?}", v.status);
        }
    }
}

/// The model neglects radiation (PHYSICS.md §8). A level is only physically consistent if
/// the energy radiated along its reference flight, by the Liénard formula, is below the
/// numerical accuracy. The radiated fraction per close pass is about r_cl / r with
/// r_cl = q²/(mc²), which is why levels use weakly charged particles and strongly charged
/// electrodes.
#[test]
fn neglected_radiation_is_below_numerical_accuracy() {
    for (name, level) in shipped_levels() {
        for (i, scn) in level
            .scenarios(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let tr = run(
                scn,
                &RunSettings::with_tolerance(level.physics.tolerances.verify),
            );
            let fraction = tr.radiated_energy / tr.kinetic_initial;
            if level.physics.radiation_reaction {
                // Radiation is part of the model; the Landau–Lifshitz treatment must be
                // valid instead (radiation reaction small against the Lorentz force).
                assert!(
                    tr.reaction_ratio_max < 0.05,
                    "{name}: |F_RR|/|F_L| = {:.3e}",
                    tr.reaction_ratio_max
                );
                continue;
            }
            let (shot, d) = level.flight_of(i);
            println!(
                "{name} shot {} disturbance {}: radiated / T0 = {fraction:.2e}",
                shot + 1,
                d + 1
            );
            assert!(fraction < 1e-10, "{name}: radiated fraction {fraction:.3e}");
        }
    }
}

/// Levels that include radiation reaction must need it: with radiation switched off, the
/// reference solution must not solve them (otherwise the radiation is decoration).
#[test]
fn radiation_levels_need_radiation() {
    for (name, mut level) in shipped_levels() {
        if !level.physics.radiation_reaction {
            continue;
        }
        level.physics.radiation_reaction = false;
        let all_arrive = level
            .verify_flights(&level.reference_solution)
            .iter()
            .all(|v| v.outcome() == Outcome::Arrived);
        assert!(!all_arrive, "{name} is solved without radiation reaction");
    }
}

/// Levels with metal spheres: a consistent model (no unsupported combinations), and the
/// conductor model accurate at verification resolution for the reference placement
/// (boundary residual < 1e-10, PHYSICS.md §2.6).
#[test]
fn metal_levels_are_accurate_and_consistent() {
    for (name, level) in shipped_levels() {
        assert_eq!(level.model_issues(), Vec::<String>::new(), "{name}");
        if level.conductors.is_empty() {
            continue;
        }
        let (field, _) = level.field_at(
            &level.reference_solution,
            physics::conductor::Resolution::Verify,
        );
        let sources: Vec<_> = field.coulomb.charges().collect();
        let residual = field.conductors.boundary_residual(&sources);
        println!("{name}: conductor boundary residual {residual:.1e}");
        assert!(residual < 1e-10, "{name}: residual {residual:.3e}");
    }
}
