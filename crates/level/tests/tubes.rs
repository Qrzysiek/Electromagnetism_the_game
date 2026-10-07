//! Tube levels (`level::tube`, PHYSICS.md §2.11): runs of a planar diode at the game's
//! scale. Run with `cargo test --release -p level --test tubes -- --nocapture`.

#![allow(clippy::disallowed_methods)] // reference values; the engine uses libm (clippy.toml)

use level::tube::{CurrentGoal, TubeSpec, TubeVerdict};
use level::{ConductorBias, Electrode, Level};

/// A planar diode on a 24 × 16 board: cathode and anode plates 10 long, 1 thick, their
/// faces 7 apart, the anode at `v`.
pub fn diode(v: f64) -> Level {
    let mut level =
        Level::from_json(include_str!("../../../levels/01_first_bend.json")).expect("a level file");
    level.shots.clear();
    level.elements.clear();
    level.reference_solution.clear();
    // Tubes are Newtonian (`tube_issues`).
    level.physics.c = None;
    level.grid.nx = 24;
    level.grid.ny = 16;
    let plate = |x: i64, bias| Electrode {
        center: [x, 8, 0],
        length: 10.0,
        thickness: 1.0,
        height: 10.0,
        angle_deg: 90.0,
        bias,
        tunable: false,
        drive: None,
    };
    level.electrodes = vec![
        plate(8, ConductorBias::Grounded),
        plate(16, ConductorBias::Potential(v)),
    ];
    level.tube = Some(TubeSpec {
        cathode: 0,
        charge_per_mass: -1.0,
        emit_toward_deg: None,
        goal: CurrentGoal {
            electrode: 1,
            min: 0.0,
            max: 1.0,
            start: 10.0,
            end: 20.0,
        },
    });
    level
}

/// The planar diode's preview and verification at two voltages: the setup is accepted,
/// and the currents scale as V^{3/2} (exact in any geometry: the steady space-charge
/// equations are invariant under V → λV, J → λ^{3/2} J) within the verdicts' errors
/// (measured 7.938 for 8, errors 5 % and 8 %). Printed for comparison: the 1D Child–
/// Langmuir law per unit length of a 10-wide strip (the plates' finite width adds edge
/// emission, so it is not a reference), and the run times (0.4 s and 3.1 s).
#[test]
fn planar_diode_runs() {
    let mut currents = Vec::new();
    for v in [10.0, 40.0] {
        let level = diode(v);
        assert!(
            level.setup_issues(&[]).is_empty(),
            "{:?}",
            level.setup_issues(&[])
        );
        let start = std::time::Instant::now();
        let coarse = level.tube_run(&[], 1, |_| {});
        let t1 = start.elapsed().as_secs_f64();
        let fine = level.tube_run(&[], 2, |_| {});
        let t2 = start.elapsed().as_secs_f64() - t1;
        let verdict = TubeVerdict::new(coarse.current, fine.current, &level.tube.unwrap().goal);
        // 1D Child–Langmuir (k = 1): J = (√2/9π) √|q/m| V^{3/2}/d², times the width.
        let cl = 2f64.sqrt() / (9.0 * std::f64::consts::PI) * v.powf(1.5) / 49.0 * 10.0;
        println!(
            "V = {v}: preview {:.5} ({} in flight, {t1:.1} s), fine {:.5} ({} in flight, {t2:.1} s), verdict {:.5} ± {:.5}; 1D CL × width {cl:.5}",
            coarse.current,
            coarse.in_flight,
            fine.current,
            fine.in_flight,
            verdict.fine,
            verdict.error
        );
        currents.push((verdict.fine, verdict.error));
    }
    let ratio = currents[1].0 / currents[0].0;
    let rel = currents[0].1 / currents[0].0 + currents[1].1 / currents[1].0;
    println!("ratio {ratio:.3} (V^3/2: 8, ± {:.3})", 8.0 * rel);
    assert!((ratio - 8.0).abs() <= 8.0 * rel, "ratio {ratio}");
}

/// Calibration of `ERROR_FACTOR` at the game's scale: the planar diode at V = 40 with
/// refine 1, 2 and 4 (measured 2.68930, 2.55052, 2.53288: the refine-2 run is 0.7 % from
/// refine 4, 0.13 of its difference to refine 1; first-order extrapolation of (1, 2)
/// overshoots by 4 %).
/// Slow (about a minute): run on demand with `--ignored`.
#[test]
#[ignore = "calibration, about a minute"]
fn planar_diode_convergence() {
    let level = diode(40.0);
    let runs: Vec<f64> = [1, 2, 4]
        .iter()
        .map(|&r| {
            let start = std::time::Instant::now();
            let c = level.tube_run(&[], r, |_| {}).current;
            println!(
                "refine {r}: {c:.5} ({:.1} s)",
                start.elapsed().as_secs_f64()
            );
            c
        })
        .collect();
    println!(
        "refine 2 off refine 4 by {:.2e}: {:.2} of |I1 − I2|; first-order extrapolation off by {:.2e}",
        runs[1] / runs[2] - 1.0,
        (runs[1] - runs[2]).abs() / (runs[0] - runs[1]).abs(),
        (2.0 * runs[1] - runs[0]) / runs[2] - 1.0
    );
}

/// What breaks z-invariance is refused: a placed charge, a shot, and a goal on the
/// cathode itself.
#[test]
fn tube_levels_refuse_what_breaks_the_symmetry() {
    let level = diode(10.0);
    let charge = level::Element::charge([12, 8, 0], 1.0);
    assert!(
        level
            .setup_issues(&[charge])
            .iter()
            .any(|s| s.contains("placed elements"))
    );
    let mut with_shot = level.clone();
    with_shot.shots = diode_with_shot().shots;
    assert!(
        with_shot
            .setup_issues(&[])
            .iter()
            .any(|s| s.contains("shots"))
    );
    let mut bad_goal = level.clone();
    bad_goal.tube.as_mut().unwrap().goal.electrode = 0;
    assert!(
        bad_goal
            .setup_issues(&[])
            .iter()
            .any(|s| s.contains("goal"))
    );
}

fn diode_with_shot() -> Level {
    Level::from_json(include_str!("../../../levels/01_first_bend.json")).expect("a level file")
}
