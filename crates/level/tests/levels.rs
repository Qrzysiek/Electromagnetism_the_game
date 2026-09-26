//! Checks on the shipped levels.

use std::path::PathBuf;

use level::Level;
use physics::trajectory::{Outcome, RunSettings, run};
use physics::verify::verify;

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
        for (i, scn) in level
            .scenarios(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let v = verify(scn, level.tolerances());
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
