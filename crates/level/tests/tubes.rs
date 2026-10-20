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
        solenoid_drive: None,
        particle_rate: level::tube::TL_RATE,
        particles: 1.0,
        emitter_mesh: 1.0,
        output: None,
        metals: Vec::new(),
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
        heat_goals: Vec::new(),
    });
    level
}

/// The planar diode's preview and verification at two voltages: the setup is accepted,
/// and the currents scale as V^{3/2} (exact in any geometry: the steady space-charge
/// equations are invariant under V → λV, J → λ^{3/2} J) within the verdicts' errors
/// (measured 8.015 for 8, errors 5 % and 4 %). Printed for comparison: the 1D Child–
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
/// refine 1, 2 and 4 (measured 2.19471, 2.14177, 2.11472: the refine-2 run is 1.3 % from
/// refine 4, 0.51 of its difference to refine 1; PHYSICS.md §2.11).
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

/// The electromagnetic field of a run (`level::tube_em`, PHYSICS.md §2.11) on Bipolar diode's
/// preview: electrons and ions emitted, absorbed, and the electrodes' circuit currents on
/// their leads. Gauss's law must hold at every vacuum node at every 50th step, to 1e-9 of
/// the largest charge density (the deposition is exact; rounding accumulates over the
/// steps).
#[test]
fn em_field_keeps_gauss_law() {
    let level = level::shipped("bipolar_diode");
    let player: Vec<level::Element> = level.reference_solution.clone();
    let run = level.tube_run(&player, 1, 8.0);
    let mut em = level::tube_em::TubeEm::new(&level, &run);
    assert!(em.leads.iter().all(Option::is_some), "a lead per electrode");
    let (nx, ny) = (em.fdtd.grid.nx, em.fdtd.grid.ny);
    let mut worst = 0.0_f64;
    let mut checks = 0;
    let mut t = run.frames[0].t;
    while t < 8.0 {
        t += 50.0 * em.fdtd.dt;
        em.advance_to(&run, t);
        let rho = em.charge_density();
        let peak = rho.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
        if peak == 0.0 {
            continue;
        }
        for j in 2..ny - 1 {
            for i in 2..nx - 1 {
                if em.in_vacuum(i, j) {
                    let r = rho[i + (nx + 1) * j];
                    worst = worst.max((em.fdtd.gauss(i, j) - r).abs() / peak);
                }
            }
        }
        checks += 1;
    }
    println!(
        "EM: Gauss's law to {worst:.2e} of the peak density over {checks} checks, c = {:.1}",
        em.fdtd.c
    );
    assert!(checks > 5 && worst < 1e-9, "EM Gauss {worst:e}");
}

/// A closed box's magnetic map shows the run's currents (PHYSICS.md §2.11, the smoothing):
/// on Vacuum diode's reference preview, the shown `B_z` sampled every 0.05 over its
/// steady current (t = 12–20) must fluctuate about its mean by less than 1 % of the mean's
/// energy (10 % in amplitude, below what the map's logarithmic colours show). The field's
/// own `B_z`, which keeps every wave the macroparticles radiate in the lossless box, is
/// measured beside it (printed).
#[test]
fn closed_box_map_shows_the_circuit() {
    let level = level::shipped("vacuum_diode");
    let run = level.tube_run(&level.reference_solution, 1, 22.0);
    let mut em = level::tube_em::TubeEm::new(&level, &run);
    let (mut shown, mut own) = (Vec::new(), Vec::new());
    for k in 0..=160 {
        em.show(&run, 12.0 + 0.05 * f64::from(k));
        shown.push(em.shown().bz.to_vec());
        own.push(em.fdtd.bz.clone());
    }
    // The energies of the mean and of the fluctuation about it.
    let split = |samples: &[Vec<f64>]| {
        #[allow(clippy::cast_precision_loss)]
        let n = samples.len() as f64;
        let (mut mean, mut fluctuation) = (0.0, 0.0);
        for c in 0..samples[0].len() {
            let m = samples.iter().map(|s| s[c]).sum::<f64>() / n;
            mean += m * m;
            fluctuation += samples.iter().map(|s| (s[c] - m).powi(2)).sum::<f64>() / n;
        }
        (mean, fluctuation)
    };
    let (m, f) = split(&shown);
    let (om, of) = split(&own);
    println!(
        "closed box: shown B_z's fluctuation {:.2e} of its mean's energy; the field's own \
         {:.2e} (c = {:.1}, shown {:.3} later)",
        f / m,
        of / om,
        em.fdtd.c,
        em.delay()
    );
    assert!(f < 0.01 * m, "fluctuation {f:e} against {m:e}");
}

/// A lead's midpoint and the unit normal along which `B_z` rises across it.
type Probe = ((f64, f64), (f64, f64));

/// The maps' magnetic field obeys Ampère's law (PHYSICS.md §2.11): across a lead carrying a
/// current I the shown `B_z` jumps by `4π I/c`, I the lead's own current from the run's
/// charges (`TubeEm::lead_currents`, independent of the field), the jump between the means
/// over squares of side 1 one unit either side of the lead's midpoint. Over a window where
/// the currents are steady (Vacuum diode's beam over t = 12–20, Hall probe's metal over
/// 20–40) or over whole periods of a level's signal (Triode, t = 40–80) the displacement
/// current through the probe averages out (measured: below 0.3 % of the lead's). Every lead
/// carrying a fifth of its level's largest current and with no other lead within 2 of its
/// midpoint (the grid's wires' leads lie closer together); required 3 % (exact in the
/// scheme for a steady current: what remains is the run's noise and the window's ends).
#[test]
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn em_field_obeys_amperes_law() {
    let cases = [
        ("vacuum_diode", 12.0, 20.0),
        ("hall_probe", 20.0, 40.0),
        ("triode", 40.0, 80.0),
    ];
    let mut checked = 0;
    for (slug, t0, t1) in cases {
        let level = level::shipped(slug);
        let run = level.tube_run(&level.reference_solution, 1, t1 + 3.0);
        let mut em = level::tube_em::TubeEm::new(&level, &run);
        let g = em.fdtd.grid;
        let leads = em.leads.clone();
        // B_z averaged over a square of side 1 about a point.
        let patch = |bz: &[f64], x: f64, y: f64| {
            let (mut s, mut n) = (0.0_f64, 0.0_f64);
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let cx = g.origin.0 + (i as f64 + 0.5) * g.h;
                    let cy = g.origin.1 + (j as f64 + 0.5) * g.h;
                    if (cx - x).abs() < 0.5 && (cy - y).abs() < 0.5 {
                        s += bz[i + g.nx * j];
                        n += 1.0;
                    }
                }
            }
            s / n.max(1.0)
        };
        // Each lead's midpoint and the unit normal `ẑ × direction` (B_z rises across the
        // lead along it by 4π I/c).
        let probes: Vec<Option<Probe>> = leads
            .iter()
            .map(|l| {
                l.map(|(a, b)| {
                    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                    let len = libm::hypot(dx, dy);
                    (
                        (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1)),
                        (-dy / len, dx / len),
                    )
                })
            })
            .collect();
        let n_el = leads.len();
        let (mut current, mut jump, mut samples) = (vec![0.0; n_el], vec![0.0; n_el], 0.0);
        let mut t = t0;
        while t <= t1 + 1e-9 {
            em.show(&run, t);
            let now = em.lead_currents();
            let bz = em.shown().bz;
            for e in 0..n_el {
                let Some((m, n)) = probes[e] else { continue };
                current[e] += now[e];
                jump[e] += patch(bz, m.0 + n.0, m.1 + n.1) - patch(bz, m.0 - n.0, m.1 - n.1);
            }
            samples += 1.0;
            t += 0.05;
        }
        let largest = current.iter().fold(0.0_f64, |a, c| a.max(c.abs())) / samples;
        // The distance from a point to a lead.
        let distance = |p: (f64, f64), l: (f64, f64, f64, f64)| {
            let (ax, ay, bx, by) = l;
            let (dx, dy) = (bx - ax, by - ay);
            let s = (((p.0 - ax) * dx + (p.1 - ay) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
            libm::hypot(p.0 - ax - s * dx, p.1 - ay - s * dy)
        };
        for e in 0..n_el {
            let Some((m, _)) = probes[e] else { continue };
            let i = current[e] / samples;
            let isolated = leads.iter().enumerate().all(|(o, l)| {
                o == e || l.is_none_or(|(a, b)| distance(m, (a.0, a.1, b.0, b.1)) > 2.0)
            });
            if i.abs() < 0.2 * largest || !isolated {
                continue;
            }
            let expected = 4.0 * std::f64::consts::PI * i / em.fdtd.c;
            let measured = jump[e] / samples;
            println!(
                "Ampère: {slug} electrode {}: <I> {i:.4}, B_z's jump {measured:.5} against 4 pi I/c \
                 {expected:.5} ({:+.2} %)",
                e + 1,
                100.0 * (measured / expected - 1.0)
            );
            assert!(
                (measured / expected - 1.0).abs() < 0.03,
                "{slug} electrode {}: {measured} vs {expected}",
                e + 1
            );
            checked += 1;
        }
    }
    assert!(checked >= 6, "only {checked} leads checked");
}

/// The maps show no waves at the run's frame rate (PHYSICS.md §2.11, the smoothing): the
/// run's currents step from frame to frame (Magnetron's every 0.05), and the steps radiated
/// stripes along its antenna's feed, alternating every 1.5–2 (ω ≈ 126), strongest where the
/// feed's own field passes zero. Across the feed (x = 62, y = 44–52, every 0.2, the feed's
/// line left out), the shown `B_z`'s largest second difference over t = 109–110 must stay
/// below 1 % of the largest field beside the feed (the map's logarithmic colours show
/// little below that); it was 4.5 % unsmoothed.
#[test]
fn magnetron_feed_field_is_smooth() {
    let level = level::shipped("magnetron");
    let run = level.tube_run(&level.reference_solution, 1, 112.0);
    let mut em = level::tube_em::TubeEm::new(&level, &run);
    let (mut peak, mut wiggle) = (0.0_f64, 0.0_f64);
    for k in 0..=8 {
        em.show(&run, 109.0 + 0.125 * f64::from(k));
        let mut f = em.fdtd.clone();
        f.bz.copy_from_slice(em.shown().bz);
        let profile: Vec<f64> = (0..=40)
            .map(|j| f.bz_at(62.0, 44.0 + 0.2 * f64::from(j)))
            .collect();
        for (j, b) in profile.iter().enumerate() {
            if j != 20 {
                peak = peak.max(b.abs());
            }
        }
        // Each side alone (the feed at j = 20 is a jump).
        for j in (1..19).chain(22..40) {
            wiggle = wiggle.max((profile[j - 1] - 2.0 * profile[j] + profile[j + 1]).abs());
        }
    }
    println!(
        "feed: largest second difference {wiggle:.3e} of a field up to {peak:.3e} ({:.2} %)",
        100.0 * wiggle / peak
    );
    assert!(wiggle < 0.01 * peak, "{wiggle} vs {peak}");
}

/// The Hall effect's forces as the game draws them (`TubeRun::hall_forces`, PHYSICS.md
/// §2.11): on Hall probe's reference run, and on both strips of Positive carriers' (the
/// lower one's carriers positive), once the charge has piled up (t = 20–60, every 2), the
/// electric push across each strip must oppose the magnetic push on the carriers' drift
/// at every sample, and its size summed over them be between 0.85 and 1.5 of the magnetic
/// one's (the model's Hall factor lies between 1 in a strong field and 4/π in a weak one,
/// D4, whatever the carriers' sign; the pushes' scatter over a time unit, up to ~40 % on
/// Positive carriers' weakly driven lower strip, a few per cent over the 21 samples; a
/// ratio of sums, as a mean of ratios leans high when the divisor scatters). The charge
/// piles up within the metal's dielectric relaxation time, `1/(4πσ)` ≈ 0.2: earlier than
/// the drift can be measured (over a time unit), so the build-up itself is not tested.
#[test]
fn hall_forces_balance() {
    for (slug, metals) in [("hall_probe", 1), ("positive_carriers", 2)] {
        let level = level::shipped(slug);
        let run = level.tube_run(&level.reference_solution, 1, 61.0);
        for metal in 0..metals {
            let (mut electric, mut magnetic) = (0.0, 0.0);
            for k in 0..=20 {
                let t = 20.0 + 2.0 * f64::from(k);
                let now = run.frame_between(t).expect("a frame");
                let before = run.frame_between(t - 1.0).expect("a frame");
                let h = run
                    .hall_forces(&now, &before, metal)
                    .expect("the Hall forces");
                // Across the strip (the electric push also drives the carriers along it).
                let (e, m) = (h.electric.dot(h.normal), h.magnetic.dot(h.normal));
                assert!(
                    e * m < 0.0,
                    "{slug} {metal}: the pushes must oppose at t = {t}"
                );
                electric += e.abs();
                magnetic += m.abs();
            }
            let steady = electric / magnetic;
            println!("{slug} {metal}: electric/magnetic {steady:.3} once piled up");
            assert!(
                (0.85..1.5).contains(&steady),
                "{slug} {metal}: steady {steady}"
            );
        }
    }
}

/// The leads of the field's display (`level::tube_em`, PHYSICS.md §2.11) cross no metal: on
/// every shipped level with metals, each lead, sampled at 200 points, stays out of every
/// metal's interior (a wire through a conductor; Positive carriers' inner Hall contacts
/// were wired through the other strip until the leads' rule took the metals in).
#[test]
fn leads_cross_no_metal() {
    let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut checked = 0;
    for path in std::fs::read_dir(dir)
        .expect("levels")
        .map(|e| e.expect("entry").path())
    {
        if !level::is_level_file(&path) {
            continue;
        }
        let level =
            Level::from_json(&std::fs::read_to_string(&path).expect("file")).expect("a level");
        if level.tube.as_ref().is_none_or(|t| t.metals.is_empty()) {
            continue;
        }
        let run = level.tube_run(&level.reference_solution, 1, 2.0);
        let em = level::tube_em::TubeEm::new(&level, &run);
        for (e, lead) in em.leads.iter().enumerate() {
            let Some(((ax, ay), (bx, by))) = *lead else {
                continue;
            };
            for m in &run.metals {
                let inside = (1..200).any(|k| {
                    let u = f64::from(k) / 200.0;
                    let (x, y) = (ax + (bx - ax) * u, ay + (by - ay) * u);
                    x > m.lo.x && x < m.hi.x && y > m.lo.y && y < m.hi.y
                });
                assert!(
                    !inside,
                    "{}: electrode {}'s lead crosses a metal",
                    path.display(),
                    e + 1
                );
            }
        }
        checked += 1;
    }
    assert!(checked >= 5, "{checked} levels with metals");
}

/// The metals' chevrons move with the carriers' drift (`level::tube::metal_steps`, PHYSICS.md
/// §2.11): on Ohm's law's reference over the goal's window (t = 10–40) the carriers' mean
/// displacement along the strip per time must be the current over the background's charge
/// per area times the width, `v = I/(ρW)`, within 10 % (the carriers a contact takes during
/// a step are left out of that step's mean; the noise ~1 %), and across it below 5 % of
/// that.
#[test]
fn metal_drift_carries_the_current() {
    let level = level::shipped("ohms_law");
    let run = level.tube_run(&level.reference_solution, 1, 41.0);
    let current = run.currents[0].expect("the goal's current");
    let metal = &run.metals[0];
    let (density, width) = (metal.background.abs(), metal.hi.y - metal.lo.y);
    let expected = current / (density * width);
    let k0 = run.frames.partition_point(|f| f.t < 10.0);
    let k1 = run.frames.partition_point(|f| f.t <= 40.0) - 1;
    let mut shift = [0.0; 2];
    for k in k0 + 1..=k1 {
        let step = level::tube::metal_steps(&run.metals, &run.frames[k - 1], &run.frames[k]);
        shift[0] += step[0][0];
        shift[1] += step[0][1];
    }
    let span = run.frames[k1].t - run.frames[k0].t;
    let (along, across) = (shift[0] / span, shift[1] / span);
    println!(
        "drift {along:.4} along, {across:.4} across; I/(ρW) = {current:.3}/({density} × {width}) \
         = {expected:.4}"
    );
    assert!(
        (along / expected - 1.0).abs() < 0.10,
        "drift {along} against {expected}"
    );
    assert!(across.abs() < 0.05 * expected, "across {across}");
}

/// The magnetron's antenna radiates its load (PHYSICS.md §2.11): on the reference's preview
/// run, the electromagnetic field's power out of the board over t = 98–100 (the field
/// solver's) against the tank load's R⟨I²⟩ over the same span (the circuit's): two
/// independent calculations of the same energy. Required 10 % (measured within ~3 % at
/// t = 90–130 when the load was set to the antenna's radiation resistance; +4.6 % now).
#[test]
fn magnetron_antenna_radiates_its_load() {
    let level = level::shipped("magnetron");
    let run = level.tube_run(&level.reference_solution, 1, 100.0);
    let mut em = level::tube_em::TubeEm::new(&level, &run);
    em.advance_to(&run, 100.0);
    let radiated = em.mean_outflow(2.0);
    let spec = level.tube.as_ref().unwrap();
    let o = spec.output.unwrap();
    let r = spec.circuits[o.circuit].resistance.unwrap();
    let i2: Vec<f64> = run
        .frames
        .iter()
        .filter(|f| f.t >= em.fdtd.t - 2.0 && f.t <= em.fdtd.t)
        .map(|f| f64::from(f.currents[o.circuit]).powi(2))
        .collect();
    #[allow(clippy::cast_precision_loss)]
    let absorbed = r * i2.iter().sum::<f64>() / i2.len() as f64;
    println!(
        "magnetron: radiated {radiated:.4} against the load's {absorbed:.4} ({:+.1} %)",
        100.0 * (radiated / absorbed - 1.0)
    );
    assert!(
        (radiated / absorbed - 1.0).abs() < 0.1,
        "{radiated} vs {absorbed}"
    );
}

/// A cut (`ElementKind::Cut`) through a node of Magnetoresistance's metals.
fn cut(x: i64, y: i64, angle_deg: f64) -> level::Element {
    level::Element {
        node: [x, y, 0],
        kind: level::ElementKind::Cut,
        value: 0.0,
        angle_deg,
        omega: None,
        speed: None,
    }
}

/// The cuts' placement rules (`Level::check_cuts`, `metal_pieces`) on Magnetoresistance's
/// two metals (a strip 16 × 4 from (4, 13), one 4 × 10 from (10, 1); mean free path 1):
/// a cut along y at x = 12 splits the long strip in two 8 × 4, one along x at y = 6 the short
/// one in two 4 × 5; refused: a node outside the metals or on an edge, an angle not along x
/// or y, a piece narrower than four mean free paths, more cuts than the level allows.
#[test]
fn cuts_follow_their_rules() {
    use level::PlacementError;
    let mut level = level::shipped("magnetoresistance");
    level.limits.max_cuts = 2;
    let pieces = level.metal_pieces(&[cut(12, 15, 90.0), cut(12, 6, 0.0)]);
    let rects: Vec<(usize, [f64; 4])> = pieces
        .iter()
        .map(|(k, lo, hi)| (*k, [lo.x, lo.y, hi.x, hi.y]))
        .collect();
    assert_eq!(
        rects,
        vec![
            (0, [4.0, 13.0, 12.0, 17.0]),
            (0, [12.0, 13.0, 20.0, 17.0]),
            (1, [10.0, 1.0, 14.0, 6.0]),
            (1, [10.0, 6.0, 14.0, 11.0]),
        ]
    );
    assert_eq!(level.metal_pieces(&[]).len(), 2);
    let ok = |p: &[level::Element]| level.check_placement(p);
    assert!(ok(&[cut(12, 15, 90.0), cut(12, 6, 0.0)]).is_ok());
    assert!(ok(&[cut(12, 15, 270.0)]).is_ok());
    assert_eq!(
        ok(&[cut(2, 2, 0.0)]),
        Err(PlacementError::CutNotInMetal([2, 2, 0]))
    );
    assert_eq!(
        ok(&[cut(4, 15, 90.0)]),
        Err(PlacementError::CutNotInMetal([4, 15, 0]))
    );
    assert_eq!(
        ok(&[cut(12, 15, 45.0)]),
        Err(PlacementError::CutNotInMetal([12, 15, 0]))
    );
    assert_eq!(
        ok(&[cut(6, 15, 90.0)]),
        Err(PlacementError::PieceTooNarrow([6, 15, 0]))
    );
    assert_eq!(
        ok(&[cut(12, 15, 0.0)]),
        Err(PlacementError::PieceTooNarrow([12, 15, 0]))
    );
    assert_eq!(
        ok(&[cut(8, 15, 90.0), cut(16, 15, 90.0), cut(12, 6, 0.0)]),
        Err(PlacementError::TooManyCuts)
    );
}

/// A cut insulates: on Magnetoresistance at B = 0 (its supply at 7, no solenoid), over
/// its goals' window t = 10–80, a cut across the long strip (along y at x = 12: no path
/// left from contact to contact) takes its current below 3 % of the uncut strip's, and a
/// cut along the short wide strip (along x at y = 6: two strips 4 × 5 side by side between
/// the same contacts) leaves its current within 10 % (each piece's new wall carries a
/// conductor's surface charge, a few per cent of a narrow strip's conductance, PHYSICS.md
/// §2.11).
#[test]
fn cuts_insulate() {
    let mut level = level::shipped("magnetoresistance");
    level.limits.max_cuts = 2;
    let supply = level::Element::supply([21, 15, 0], 7.0);
    let currents = |placement: &[level::Element]| {
        let run = level.tube_run(placement, 1, 81.0);
        [0, 1].map(|k| run.currents[k].expect("the goal's window has passed"))
    };
    let [long, short] = currents(&[supply]);
    let [long_cut, short_cut] = currents(&[supply, cut(12, 15, 90.0), cut(12, 6, 0.0)]);
    println!(
        "long strip {long:.4} -> {long_cut:.4} (cut across), short {short:.4} -> {short_cut:.4} \
         (cut along: {:+.2} %)",
        100.0 * (short_cut / short - 1.0)
    );
    assert!(long > 0.5, "the uncut long strip conducts: {long}");
    assert!(
        long_cut.abs() < 0.03 * long,
        "a cut across: {long_cut} of {long}"
    );
    assert!(
        (short_cut / short - 1.0).abs() < 0.10,
        "a cut along: {short_cut} of {short}"
    );
}

/// An isolated metal's heat (`TubeState::induced_work`, PHYSICS.md §2.11) on Eddy currents'
/// plate at the preview's resolution over its goal's window (t = 25–125): exactly 0 without
/// a drive (the time step's numerical heating through the electrostatic field, 2.57 per time
/// on this plate, is left out), and quadratic in the drive (linear response): at B₀ = 0.4
/// and 0.8 the heats' ratio within 10 % of 4 (the preview's scatter between seeds is a few
/// per cent on this plate; the carriers' warming at 0.8, `qEℓ/kT` ≈ 0.1 at the edges'
/// middles, a few per cent too).
#[test]
fn isolated_metal_heat_is_the_induced_work() {
    let level = level::shipped("eddy_currents");
    let heat = |placement: &[level::Element]| {
        let run = level.tube_run(placement, 1, 126.0);
        run.heat[0].expect("the goal's window has passed")
    };
    let drive = |b0: f64| level::Element::supply([23, 19, 0], b0);
    let still = heat(&[]);
    let (low, high) = (heat(&[drive(0.4)]), heat(&[drive(0.8)]));
    println!(
        "heat without a drive {still}, at B0 = 0.4 {low:.4}, at 0.8 {high:.4} (ratio {:.3})",
        high / low
    );
    assert!(still == 0.0, "without a drive: {still}");
    assert!(
        (high / low / 4.0 - 1.0).abs() < 0.10,
        "ratio {}",
        high / low
    );
}
