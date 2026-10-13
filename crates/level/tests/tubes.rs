//! Tube levels (`level::tube`, PHYSICS.md §2.11): runs of a planar diode at the game's
//! scale. Run with `cargo test --release -p level --test tubes -- --nocapture`.

#![allow(clippy::disallowed_methods)] // reference values; the engine uses libm (clippy.toml)

use level::tube::{Carriers, CurrentGoal, CurrentVerdict, TubeSpec};
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
        shape: level::ElectrodeShape::Box,
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
        saturation: None,
        ions: None,
        projectiles: Vec::new(),
        b_z: 0.0,
        solenoid: None,
        particle_rate: level::tube::TL_RATE,
        step: level::tube::DEFAULT_STEP,
        goals: vec![CurrentGoal {
            electrode: 1,
            min: 0.0,
            max: 1.0,
            start: 10.0,
            end: 20.0,
            carriers: Carriers::Cathode,
        }],
        circuits: Vec::new(),
        potential_goals: Vec::new(),
    });
    level
}

/// The planar diode's preview and verification at two voltages: the setup is accepted,
/// and the currents scale as V^{3/2} (exact in any geometry: the steady space-charge
/// equations are invariant under V → λV, J → λ^{3/2} J) within the verdicts' errors
/// (measured 8.045 for 8, errors 4 % and 4 %). Printed for comparison: the 1D Child–
/// Langmuir law per unit length of a 10-wide strip (the plates' finite width adds edge
/// emission, so it is not a reference), and the run times (0.2 s and 1.0 s).
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
        let goal = level.tube.as_ref().unwrap().goals[0];
        let end = goal.end;
        let coarse = level.tube_run(&[], 1, end);
        let t1 = start.elapsed().as_secs_f64();
        let fine = level.tube_run(&[], 2, end);
        let t2 = start.elapsed().as_secs_f64() - t1;
        let (ic, i_f) = (coarse.currents[0].unwrap(), fine.currents[0].unwrap());
        let verdict = CurrentVerdict::new(ic, i_f, goal.min, goal.max);
        // 1D Child–Langmuir (k = 1): J = (√2/9π) √|q/m| V^{3/2}/d², times the width.
        let cl = 2f64.sqrt() / (9.0 * std::f64::consts::PI) * v.powf(1.5) / 49.0 * 10.0;
        println!(
            "V = {v}: preview {:.5} ({} in flight, {t1:.1} s), fine {:.5} ({} in flight, {t2:.1} s), verdict {:.5} ± {:.5}; 1D CL × width {cl:.5}",
            ic, coarse.in_flight, i_f, fine.in_flight, verdict.fine, verdict.error
        );
        currents.push((verdict.fine, verdict.error));
    }
    let ratio = currents[1].0 / currents[0].0;
    let rel = currents[0].1 / currents[0].0 + currents[1].1 / currents[1].0;
    println!("ratio {ratio:.3} (V^3/2: 8, ± {:.3})", 8.0 * rel);
    assert!((ratio - 8.0).abs() <= 8.0 * rel, "ratio {ratio}");
}

/// Calibration of `ERROR_FACTOR` at the game's scale: the planar diode at V = 40 with
/// refine 1, 2 and 4 (measured 2.20803, 2.14708, 2.11007: the refine-2 run is 1.7 % from
/// refine 4, 0.61 of its difference to refine 1; PHYSICS.md §2.11).
/// Slow (about a minute): run on demand with `--ignored`.
#[test]
#[ignore = "calibration, about a minute"]
fn planar_diode_convergence() {
    let level = diode(40.0);
    let runs: Vec<f64> = [1, 2, 4]
        .iter()
        .map(|&r| {
            let start = std::time::Instant::now();
            let c = level.tube_results(&[], r).currents[0].unwrap();
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
    bad_goal.tube.as_mut().unwrap().goals[0].electrode = 0;
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

/// Per-carrier current goals: an older level file's single `goal` loads as a list of one;
/// an ion goal needs an ion source and may not be on the ion source itself; with ions
/// from the anode, the ion goal on the cathode counts the ions arriving there (and the
/// electron goal is unchanged by the second goal).
#[test]
fn current_goals_per_carrier() {
    let old = r#"{"cathode": 0, "goal": {"electrode": 1, "min": 0.0, "max": 1.0, "start": 10.0, "end": 20.0}}"#;
    let spec: level::tube::TubeSpec = serde_json::from_str(old).unwrap();
    assert_eq!(spec.goals.len(), 1);
    assert_eq!(spec.goals[0].carriers, Carriers::Cathode);

    let mut level = diode(20.0);
    let ion_goal = CurrentGoal {
        electrode: 0,
        min: 0.0,
        max: 1.0,
        start: 10.0,
        end: 20.0,
        carriers: Carriers::Ions,
    };
    level.tube.as_mut().unwrap().goals.push(ion_goal);
    assert!(
        level
            .setup_issues(&[])
            .iter()
            .any(|s| s.contains("ion source"))
    );
    let spec = level.tube.as_mut().unwrap();
    spec.ions = Some(level::tube::IonSource {
        electrode: 1,
        charge_per_mass: 0.04,
        emit_toward_deg: None,
    });
    spec.goals[1].electrode = 1;
    assert!(
        level
            .setup_issues(&[])
            .iter()
            .any(|s| s.contains("emitter"))
    );
    level.tube.as_mut().unwrap().goals[1].electrode = 0;
    assert!(
        level.setup_issues(&[]).is_empty(),
        "{:?}",
        level.setup_issues(&[])
    );
    let mut sim = level.tube_sim(&[], 1);
    sim.advance_to(20.0);
    let (electrons, ions) = (sim.currents()[0].unwrap(), sim.currents()[1].unwrap());
    println!(
        "electron current into the anode {electrons:.5}, ion current into the cathode {ions:.5}"
    );
    assert!(ions > 0.0 && electrons > 0.0);
    let mut single = level.clone();
    single.tube.as_mut().unwrap().goals.truncate(1);
    let mut alone = single.tube_sim(&[], 1);
    alone.advance_to(20.0);
    assert_eq!(alone.currents()[0], Some(electrons));
}

/// Coils in tube circuits: a coil needs a finite resistance (an open branch carries
/// nothing), a branch may not return to its own electrode, and a tank (the anode through R
/// and L to its supply) is accepted and runs.
#[test]
fn coils_in_tube_circuits() {
    let mut level = diode(10.0);
    let tank = level::tube::TubeCircuit {
        electrode: 1,
        source: level::Source::Dc { value: 10.0 },
        resistance: Some(1.0),
        capacitance: 1.0,
        supply_offset: false,
        linked: Vec::new(),
        inductance: 2.0,
        return_to: None,
    };
    let refused = |level: &Level, c: level::tube::TubeCircuit| {
        let mut l = level.clone();
        l.tube.as_mut().unwrap().circuits = vec![c];
        l.setup_issues(&[]).iter().any(|s| s.contains("circuit 1"))
    };
    assert!(refused(
        &level,
        level::tube::TubeCircuit {
            resistance: None,
            ..tank.clone()
        }
    ));
    assert!(refused(
        &level,
        level::tube::TubeCircuit {
            return_to: Some(1),
            ..tank.clone()
        }
    ));
    assert!(!refused(&level, tank.clone()));
    level.tube.as_mut().unwrap().circuits = vec![tank];
    let run = level.tube_run(&[], 1, 5.0);
    assert!(
        run.frames
            .iter()
            .all(|f| f.potentials.iter().all(|v| v.is_finite()))
    );
}

/// A coaxial tube: a rod (diameter 2) inside a ring (mid-diameter 6, wall 1), on a 24 × 16
/// board.
fn coax_level() -> Level {
    let mut level = diode(10.0);
    let round = |shape, length: f64, thickness: f64, bias| Electrode {
        center: [12, 8, 0],
        length,
        thickness,
        height: 10.0,
        angle_deg: 0.0,
        shape,
        bias,
        tunable: false,
        drive: None,
    };
    level.electrodes = vec![
        round(
            level::ElectrodeShape::Rod,
            2.0,
            2.0,
            ConductorBias::Potential(1.0),
        ),
        round(
            level::ElectrodeShape::Ring,
            6.0,
            1.0,
            ConductorBias::Grounded,
        ),
    ];
    level
}

/// Rods and rings (tube levels): the level's shapes reach the boundary elements as circles
/// and annuli, so the rod's capacitance per unit length inside the ring is the coaxial
/// `1/(2 ln(b/a))` with a = 1 and b = 2.5 (the ring's inner radius; the grounded envelope
/// outside the ring does not matter): required to 1e-3 with segments of 0.1 (V2: 2e-4 at
/// 0.05). The setup checks use the shapes: the rod in the ring's hole is clear of it, a
/// rod reaching into the wall touches it, a ring beyond the board is refused, and a 3D
/// level refuses rods.
#[test]
fn rods_and_rings() {
    let level = coax_level();
    assert!(
        level.setup_issues(&[]).is_empty(),
        "{:?}",
        level.setup_issues(&[])
    );
    let prisms = physics::zinv::Electrodes::new(level.tube_prisms(&[]), &[], 0.1);
    let c = prisms.capacitance[0][0];
    let exact = 1.0 / (2.0 * (2.5_f64).ln());
    println!(
        "coaxial rod in a ring: C = {c:.6}, exact {exact:.6} ({:.2e})",
        c / exact - 1.0
    );
    assert!((c / exact - 1.0).abs() < 1e-3, "C {c} vs {exact}");

    let mut fat = level.clone();
    fat.electrodes[0].length = 5.0;
    assert!(
        fat.setup_issues(&[])
            .iter()
            .any(|s| s.contains("touch or overlap"))
    );
    let mut wide = level.clone();
    wide.electrodes[1].length = 20.0;
    assert!(
        wide.setup_issues(&[])
            .iter()
            .any(|s| s.contains("beyond the board"))
    );
    let mut flat = level::shipped("first_bend");
    flat.electrodes = level.electrodes.clone();
    assert!(
        flat.setup_issues(&[])
            .iter()
            .any(|s| s.contains("rods and rings"))
    );
    // A free charge in the ring's hole starts in vacuum, one in its wall in metal.
    let free = |x: i64| level::Element {
        kind: level::ElementKind::Free,
        speed: Some(0.0),
        ..level::Element::charge([x, 8, 0], -0.01)
    };
    assert!(
        !level
            .setup_issues(&[free(14)])
            .iter()
            .any(|s| s.contains("starts in metal"))
    );
    assert!(
        level
            .setup_issues(&[free(15)])
            .iter()
            .any(|s| s.contains("starts in metal"))
    );
}

/// The klystron (level Klystron) against small-signal ballistic bunching theory, from the
/// level's own fields: the buncher (its two plates at 1, everything else at 0) gives an
/// electron crossing it at v0 the energy `V1 Re(G e^{iωt0})`, `G = ∫ E_x e^{iωx/v0} dx`
/// along the beam's axis (y = 10), so the speed is modulated by `Δv/v0 = V1 |G| / (2 V0)`
/// and after the drift L between the buncher's and the catcher's centres (both
/// symmetric: 21) the fundamental of the beam current is `2 I0 J1(X)`,
/// `X = (Δv/v0) ω L / v0` (Webster), which the catcher's tank turns into its amplitude
/// (a constant factor). Predicted: the optimum drive `V1 = 1.841 · 2 V0 v0 / (|G| ω L)`.
/// Measured: the preview's catcher amplitude at V0 = 50 over drives around it. Required
/// (fixed before measuring): the measured optimum within 15 % of the prediction (at the
/// optimum the speeds are modulated by ~14 %, where the small-signal theory's
/// corrections are of that order), and the rising side's shape, the amplitude at half
/// the measured optimum's drive relative to the largest, within 0.05 of
/// J1(0.92)/J1(1.841) = 0.711. A first version compared the shape at fractions of the
/// predicted optimum, which counts the optimum's allowed 15 % twice (it failed by 0.09
/// with the optimum 13 % off). Past the optimum the thin-gap theory does not hold here
/// (a buncher of two gaps 7 apart, between which the electrons already drift, and a
/// catcher 4 long): printed, not required (PHYSICS.md §2.11).
#[test]
#[allow(clippy::many_single_char_names)]
fn klystron_against_bunching_theory() {
    let level = level::shipped("klystron");
    let (v0_cathode, omega, l_drift) = (50.0_f64, 2.0 * std::f64::consts::PI, 21.0);
    let v0 = (2.0 * v0_cathode).sqrt();
    // The buncher's field per volt along the axis.
    let mut prisms = level.tube_prisms(&[]);
    for (i, p) in prisms.iter_mut().enumerate() {
        p.bias = if i == 5 || i == 6 {
            physics::conductor::Bias::Potential(1.0)
        } else {
            physics::conductor::Bias::Grounded
        };
    }
    let field = physics::zinv::Electrodes::new(prisms, &[], 0.1);
    let (dx, mut g_re, mut g_im) = (0.02, 0.0, 0.0);
    let mut x = 7.0;
    while x < 30.0 {
        let e = field.sample(physics::DVec3::new(x, 10.0, 0.0)).1;
        // Energy gained by an electron (charge -1) moving +x: -E_x dx per volt.
        let k = omega * x / v0;
        g_re -= e.x * k.cos() * dx;
        g_im -= e.x * k.sin() * dx;
        x += dx;
    }
    let g = (g_re * g_re + g_im * g_im).sqrt();
    let theta0 = omega * l_drift / v0;
    let predicted = 1.841 * 2.0 * v0_cathode / (g * theta0);
    // The measured drive curve (preview).
    let amplitude = |v1: f64| {
        let placement = [
            level::Element {
                kind: level::ElementKind::Supply,
                ..level::Element::charge([3, 10, 0], -v0_cathode)
            },
            level::Element {
                kind: level::ElementKind::Supply,
                ..level::Element::charge([17, 13, 0], v1)
            },
        ];
        level.tube_results(&placement, 1).potential[0].unwrap()
    };
    let drives: Vec<f64> = [0.7, 0.85, 1.0, 1.15, 1.3]
        .iter()
        .map(|f| f * predicted)
        .collect();
    let curve: Vec<f64> = drives.iter().map(|&v| amplitude(v)).collect();
    // The peak of the parabola through the largest sample and its neighbours.
    let k = (1..curve.len() - 1)
        .max_by(|&a, &b| curve[a].total_cmp(&curve[b]))
        .unwrap();
    let (a, b, c) = (curve[k - 1], curve[k], curve[k + 1]);
    let h = drives[1] - drives[0];
    let peak_at = drives[k] + 0.5 * h * (a - c) / (a - 2.0 * b + c);
    let peak = b - 0.125 * (a - c) * (a - c) / (a - 2.0 * b + c);
    let half = amplitude(0.5 * peak_at) / peak;
    let over = amplitude(1.5 * peak_at) / peak;
    let j1 = |x: f64| libm::j1(x) / libm::j1(1.841);
    println!(
        "klystron: |G| = {g:.4}, theta0 = {theta0:.3}, predicted optimum V1 = {predicted:.3}, \
         measured {peak_at:.3} ({:+.1} %); amplitudes {curve:.4?} at {drives:.3?}; at half the \
         optimum's drive {half:.3} (J1: {:.3}); at 1.5 times {over:.3} (thin-gap J1: {:.3})",
        100.0 * (peak_at / predicted - 1.0),
        j1(0.5 * 1.841),
        j1(1.5 * 1.841)
    );
    assert!(
        (peak_at / predicted - 1.0).abs() < 0.15,
        "optimum {peak_at} vs {predicted}"
    );
    assert!((half - j1(0.5 * 1.841)).abs() < 0.05, "half drive {half}");
}
