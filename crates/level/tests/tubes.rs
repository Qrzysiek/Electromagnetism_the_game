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
        label: String::new(),
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
        label: String::new(),
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

/// A cut (`ElementKind::Cut`) from node `a` to node `b`, as the player draws it.
fn cut(level: &Level, a: [i64; 2], b: [i64; 2]) -> level::Element {
    level.cut_between([a[0], a[1], 0], [b[0], b[1], 0])
}

/// The cuts' rules (`Level::check_cuts`, `cut_metal`) on Magnetoresistance's two metals (a
/// strip 16 × 4 from (4, 13), one 4 × 10 from (10, 1)): any length and angle within one
/// metal, its ends on its edges too; refused: an end outside, a cut along an edge (its
/// middle not inside), across two metals, without a length, more than the level allows.
#[test]
fn cuts_follow_their_rules() {
    use level::PlacementError;
    let mut level = level::shipped("magnetoresistance");
    level.limits.max_cuts = 2;
    let ok = |p: &[level::Element]| level.check_placement(p);
    for (a, b) in [
        ([12, 13], [12, 17]),
        ([10, 13], [14, 17]),
        ([8, 14], [16, 16]),
        ([10, 6], [14, 6]),
    ] {
        assert!(ok(&[cut(&level, a, b)]).is_ok(), "{a:?} to {b:?}");
        assert_eq!(
            level.cut_metal(&cut(&level, a, b)),
            Some(usize::from(a[1] < 12))
        );
    }
    let (p, q) = level
        .cut_segment(&cut(&level, [10, 13], [14, 17]))
        .expect("a cut");
    assert!((p.x - 10.0).abs() < 1e-12 && (p.y - 13.0).abs() < 1e-12);
    assert!((q.x - 14.0).abs() < 1e-12 && (q.y - 17.0).abs() < 1e-12);
    for (a, b) in [
        ([12, 15], [12, 19]),
        ([4, 13], [20, 13]),
        ([12, 11], [12, 13]),
        ([12, 15], [12, 15]),
    ] {
        assert_eq!(
            ok(&[cut(&level, a, b)]),
            Err(PlacementError::CutNotInMetal([a[0], a[1], 0])),
            "{a:?} to {b:?}"
        );
    }
    assert_eq!(
        ok(&[
            cut(&level, [8, 13], [8, 17]),
            cut(&level, [16, 13], [16, 17]),
            cut(&level, [10, 6], [14, 6])
        ]),
        Err(PlacementError::TooManyCuts)
    );
}

/// A cut insulates (its carriers reflected off it as off the metal's sides): on
/// Magnetoresistance at B = 0 (its supply at 7, no solenoid), over its goals' window
/// t = 10–80, a cut from edge to edge across the long strip, straight (x = 12) or oblique
/// ((10, 13) to (14, 17)), takes its current below 3 % of the uncut strip's (no path left;
/// one that ends on an edge reaches a little past it); one along the short wide strip (y = 6:
/// two strips 4 × 5 between the same contacts) leaves its current within 10 % (each new wall
/// carries a conductor's surface charge, PHYSICS.md §2.11); one across the long strip that
/// stops a unit short of its top leaves a path: the current stays above 5 %.
#[test]
fn cuts_insulate() {
    let mut level = level::shipped("magnetoresistance");
    level.limits.max_cuts = 2;
    let supply = level::Element::supply([21, 15, 0], 7.0);
    let currents = |cuts: &[level::Element]| {
        let mut placement = vec![supply];
        placement.extend_from_slice(cuts);
        assert!(level.check_placement(&placement).is_ok());
        let run = level.tube_run(&placement, 1, 81.0);
        [0, 1].map(|k| run.currents[k].expect("the goal's window has passed"))
    };
    let [long, short] = currents(&[]);
    let [across, along] = currents(&[
        cut(&level, [12, 13], [12, 17]),
        cut(&level, [10, 6], [14, 6]),
    ]);
    let [oblique, _] = currents(&[cut(&level, [10, 13], [14, 17])]);
    let [partial, _] = currents(&[cut(&level, [12, 13], [12, 16])]);
    println!(
        "long strip {long:.4}: cut across {across:.4}, oblique {oblique:.4}, a unit short \
         {partial:.4}; short strip {short:.4}, cut along {along:.4} ({:+.2} %)",
        100.0 * (along / short - 1.0)
    );
    assert!(long > 0.5, "the uncut long strip conducts: {long}");
    assert!(
        across.abs() < 0.03 * long,
        "a cut across: {across} of {long}"
    );
    assert!(
        oblique.abs() < 0.03 * long,
        "an oblique cut: {oblique} of {long}"
    );
    assert!(
        partial > 0.05 * long,
        "a cut a unit short: {partial} of {long}"
    );
    assert!(
        (along / short - 1.0).abs() < 0.10,
        "a cut along: {along} of {short}"
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

/// Y10, Onsager's reciprocity in the Drude metal (PHYSICS.md §9, Symmetries): Hall probe at
/// its reference supply (30) with the solenoid at +3 and at −3, over t = 10–80: the current
/// along the strip is even in B (the same both ways, within 4 %: the runs' scatter is ~1.3 %)
/// and the Hall voltage odd (its two readings cancel to within 15 % of its size: the
/// floating contacts carry ~0.2 of thermal noise).
#[test]
fn y10_onsager_in_the_metal() {
    let level = level::shipped("hall_probe");
    let run = |b: f64| {
        let placement = [
            level::Element::supply([21, 8, 0], 30.0),
            level::Element::supply([23, 15, 0], b),
        ];
        let r = level.tube_run(&placement, 1, 81.0);
        let current = r.currents[0].expect("the window has passed");
        let hall = r.potential[0].expect("the window has passed");
        (current, hall)
    };
    let ((i_up, v_up), (i_down, v_down)) = (run(3.0), run(-3.0));
    println!(
        "Y10 Onsager: current {i_up:.4} at B = +3, {i_down:.4} at −3 ({:+.2} %); Hall voltage \
         {v_up:.4} and {v_down:.4} (sum {:.4})",
        100.0 * (i_down / i_up - 1.0),
        v_up + v_down
    );
    assert!(
        (i_down / i_up - 1.0).abs() < 0.04,
        "current {i_up} {i_down}"
    );
    assert!(
        (v_up + v_down).abs() < 0.15 * 0.5 * (v_up - v_down).abs(),
        "Hall voltage {v_up} {v_down}"
    );
}

/// The quarter turn (x, y) → (ny − y, x) of a tube level: its board, electrodes (their
/// angles turned too), metals, solenoid and reference.
fn turned(level: &Level) -> Level {
    let mut l = level.clone();
    let ny = i64::from(level.grid.ny);
    let turn = |n: [i64; 3]| [ny - n[1], n[0], n[2]];
    (l.grid.nx, l.grid.ny) = (level.grid.ny, level.grid.nx);
    for e in &mut l.electrodes {
        e.center = turn(e.center);
        e.angle_deg += 90.0;
    }
    for e in &mut l.reference_solution {
        e.node = turn(e.node);
    }
    if let Some(t) = &mut l.tube {
        t.solenoid = t.solenoid.map(turn);
        for m in &mut t.metals {
            let (a, b) = (turn(m.corners[0]), turn(m.corners[1]));
            m.corners = [
                [a[0].min(b[0]), a[1].min(b[1]), 0],
                [a[0].max(b[0]), a[1].max(b[1]), 0],
            ];
        }
    }
    l
}

/// Y11, the isotropy of the tube engine with a Drude metal: Ohm's law's strip turned a
/// quarter (its contacts, metal and supply with it; its cores fall elsewhere, from the same
/// seed in the turned rectangle) carries the same current over t = 10–40, within 4 % (the
/// runs' scatter, ~1.3 %): no axis of the board is preferred.
#[test]
fn y11_isotropy_of_the_metal() {
    let level = level::shipped("ohms_law");
    let t = turned(&level);
    assert!(t.check_placement(&t.reference_solution).is_ok());
    assert!(
        t.setup_issues(&t.reference_solution).is_empty(),
        "{:?}",
        t.setup_issues(&t.reference_solution)
    );
    let current = |l: &Level| {
        let r = l.tube_run(&l.reference_solution, 1, 41.0);
        r.currents[0].expect("the window has passed")
    };
    let (along_x, along_y) = (current(&level), current(&t));
    println!(
        "Y11 isotropy: current {along_x:.4} along x, {along_y:.4} along y ({:+.2} %)",
        100.0 * (along_y / along_x - 1.0)
    );
    assert!(
        (along_y / along_x - 1.0).abs() < 0.04,
        "{along_x} {along_y}"
    );
}

/// The mirror x → nx − x of a tube level: its electrodes (their angles mirrored), metals,
/// emission direction, solenoid, reference, and its field reversed (B is axial: the mirror
/// reverses its z).
fn mirrored(level: &Level) -> Level {
    let mut l = level.clone();
    let nx = i64::from(level.grid.nx);
    let flip = |n: [i64; 3]| [nx - n[0], n[1], n[2]];
    for e in &mut l.electrodes {
        e.center = flip(e.center);
        e.angle_deg = 180.0 - e.angle_deg;
    }
    for e in &mut l.reference_solution {
        e.node = flip(e.node);
    }
    if let Some(t) = &mut l.tube {
        t.solenoid = t.solenoid.map(flip);
        t.emit_toward_deg = t.emit_toward_deg.map(|a| 180.0 - a);
        t.b_z = -t.b_z;
        for m in &mut t.metals {
            let (a, b) = (flip(m.corners[0]), flip(m.corners[1]));
            m.corners = [[b[0], a[1], 0], [a[0], b[1], 0]];
        }
    }
    l
}

/// Y18, the vacuum tube engine's mirror symmetry (PHYSICS.md §9, Symmetries): Vacuum diode
/// mirrored left to right (the cathode emitting the other way) delivers the same current to
/// its anode over t = 10–20, within 2 % (the emission points along a face, a low-discrepancy
/// sequence, are not mirror images of each other: the macroparticles differ, ~1 %).
#[test]
fn y18_vacuum_tube_mirror() {
    let level = level::shipped("vacuum_diode");
    let m = mirrored(&level);
    assert!(m.check_placement(&m.reference_solution).is_ok());
    let current = |l: &Level| {
        let r = l.tube_run(&l.reference_solution, 1, 21.0);
        r.currents[0].expect("the window has passed")
    };
    let (here, mirror) = (current(&level), current(&m));
    println!(
        "Y18 vacuum diode: current {here:.5}, mirrored {mirror:.5} ({:+.3} %)",
        100.0 * (mirror / here - 1.0)
    );
    assert!((mirror / here - 1.0).abs() < 0.02, "{here} {mirror}");
}

/// The metal maps (`TubeRun::metal_rates`, PHYSICS.md §2.11, the Drude metal's maps). (a) The
/// heat map adds up to the heat goal: Eddy currents' reference, the plate's cells' power over
/// the goal's window against the goal's measure, within 1 % (the frames' accounts are
/// interpolated between frames 0.28 apart at the window's ends: ~0.3 % of its work). (b) The
/// current map carries the strip's current: Ohm's law's reference (a strip between two
/// contacts alone), the mean `j_x` over its middle half times its width against the current
/// into its end over the goal's window (t = 10–40), within 3 % (the cells' thermal noise
/// averages to ~1 %), flowing from the supply's end; and no current into its walls: in
/// each row along the strip's middle half the mean `j_y` within 15 % of `|j_x|` (the noise of
/// a row of 8 cells over 30 time units ~4 %; booked to the cell where each step began, the
/// charge carried had given the rows along the walls ±70 %). (Not Hall probe: its floating side
/// contacts at mid-length carry part of the current past the metal under them, and charge
/// through the strip until t ≈ 60: PHYSICS.md.) (c) The eddy currents circulate against the field's change as Ohm's law has them:
/// Ohmic, the stream function is `ψ = −σḂu` (`u` the plate's torsion function), so the
/// circulation `L = ∫(r × j)_z dA = 2∫ψ dA = −2σḂ∫u dA` and the heat `P = σḂ²∫u dA`: over
/// t = 40–60 (Ḃ > 0 throughout), `L/P = −2⟨Ḃ⟩/⟨Ḃ²⟩` whatever σ and the plate's shape, within
/// 10 % (the heat came 5–7 % below Ohm's law: the walls). (d) The maps as shown
/// (`TubeRun::metal_maps`) under the swinging field: over the two periods before t = 125 the
/// lock-in current is the rms of the eddy currents, `j_rms = −σ Ḃ_rms ∇ × (u ẑ)`, so its
/// circulation against the power is `−2/Ḃ_rms`, within 10 %. (e) Their error bars are
/// honest: on Ohm's law's uniform middle at t = 40 (4 time units in five batches) the
/// scatter of `j_x` across the cells against the standard errors they report, a ratio
/// within 0.6–1.6 (five batches estimate each error to ~35 %).
#[test]
#[allow(clippy::cast_precision_loss)]
fn metal_maps_add_up() {
    // (a) and (c): Eddy currents.
    let eddy = level::shipped("eddy_currents");
    let (start, end) = {
        let g = &eddy.tube.as_ref().expect("a tube level").heat_goals[0];
        (g.start, g.end)
    };
    let r = eddy.tube_run(&eddy.reference_solution, 1, end + 1.0);
    let m = &r.metals[0];
    let area = (m.hi.x - m.lo.x) * (m.hi.y - m.lo.y) / (m.cells.0 * m.cells.1) as f64;
    let rates = r.metal_rates(start, end, false).expect("frames");
    let mapped: f64 = rates[0].power.iter().sum::<f64>() * area;
    let heat = r.heat[0].expect("the window has passed");
    let (t0, t1) = (40.0, 60.0);
    let rates = r.metal_rates(t0, t1, false).expect("frames");
    let centre = [0.5 * (m.lo.x + m.hi.x), 0.5 * (m.lo.y + m.hi.y)];
    let circulation: f64 = rates[0]
        .j
        .iter()
        .enumerate()
        .map(|(k, j)| {
            let x = m.cell_centre(k);
            ((x[0] - centre[0]) * j[1] - (x[1] - centre[1]) * j[0]) * area
        })
        .sum();
    let power: f64 = rates[0].power.iter().sum::<f64>() * area;
    let n = 4000;
    let h = (t1 - t0) / f64::from(n);
    let bdot = |t: f64| (r.b_at(t + 1e-4) - r.b_at(t - 1e-4)) / 2e-4;
    let mean_bdot = (r.b_at(t1) - r.b_at(t0)) / (t1 - t0);
    let mean_bdot2 = (0..n)
        .map(|k| bdot(t0 + (f64::from(k) + 0.5) * h).powi(2))
        .sum::<f64>()
        / f64::from(n);
    let ohmic = -2.0 * mean_bdot / mean_bdot2;
    // (d): the maps as shown, by lock-in.
    let maps = r.metal_maps(125.0, 4.0).expect("maps");
    assert!(maps.lock_in, "a swinging field's maps lock in");
    let lock_l: f64 = maps.values[0]
        .j
        .iter()
        .enumerate()
        .map(|(k, j)| {
            let x = m.cell_centre(k);
            ((x[0] - centre[0]) * j[1] - (x[1] - centre[1]) * j[0]) * area
        })
        .sum();
    let lock_p: f64 = maps.values[0].power.iter().sum::<f64>() * area;
    let (f0, f1) = (maps.from, maps.to);
    let rms = ((0..n)
        .map(|k| bdot(f0 + (f64::from(k) + 0.5) * (f1 - f0) / f64::from(n)).powi(2))
        .sum::<f64>()
        / f64::from(n))
    .sqrt();
    let lock_ohmic = -2.0 / rms;
    // (b): Ohm's law.
    let ohm = level::shipped("ohms_law");
    let rh = ohm.tube_run(&ohm.reference_solution, 1, 41.0);
    let current = rh.currents[0].expect("the window has passed");
    let mh = &rh.metals[0];
    let rates = rh.metal_rates(10.0, 40.0, false).expect("frames");
    let length = mh.hi.x - mh.lo.x;
    let middle: Vec<f64> = (0..rates[0].j.len())
        .filter(|k| {
            let x = mh.cell_centre(*k)[0];
            x > mh.lo.x + 0.25 * length && x < mh.hi.x - 0.25 * length
        })
        .map(|k| rates[0].j[k][0])
        .collect();
    let jx = middle.iter().sum::<f64>() / middle.len() as f64;
    let carried = -jx * (mh.hi.y - mh.lo.y);
    // (e): the error bars on the uniform middle.
    let shown = rh.metal_maps(40.0, 4.0).expect("maps");
    let middle_cells: Vec<usize> = (0..shown.values[0].j.len())
        .filter(|k| {
            let x = mh.cell_centre(*k)[0];
            x > mh.lo.x + 0.25 * length && x < mh.hi.x - 0.25 * length
        })
        .collect();
    let values: Vec<f64> = middle_cells
        .iter()
        .map(|k| shown.values[0].j[*k][0])
        .collect();
    let mean_value = values.iter().sum::<f64>() / values.len() as f64;
    let scatter = (values.iter().map(|v| (v - mean_value).powi(2)).sum::<f64>()
        / (values.len() - 1) as f64)
        .sqrt();
    let reported = (middle_cells
        .iter()
        .map(|k| shown.errors[0].j[*k][0].powi(2))
        .sum::<f64>()
        / middle_cells.len() as f64)
        .sqrt();
    let (nx, ny) = mh.cells;
    let rows: Vec<f64> = (0..ny)
        .map(|row| {
            let cells: Vec<f64> = (0..nx)
                .filter(|i| {
                    let x = mh.cell_centre(i + nx * row)[0];
                    x > mh.lo.x + 0.25 * length && x < mh.hi.x - 0.25 * length
                })
                .map(|i| rates[0].j[i + nx * row][1])
                .collect();
            cells.iter().sum::<f64>() / cells.len() as f64
        })
        .collect();
    println!(
        "metal maps: (a) heat map {mapped:.5} against the goal's {heat:.5} ({:+.3} %); (b) the \
         strip's middle carries {carried:.4} (j_x {jx:.4}) against {current:.4} into its end \
         ({:+.2} %), j_y by row {rows:.4?}; (c) circulation {circulation:.3}, power {power:.4}: L/P {:.4} against \
         Ohm's −2⟨Ḃ⟩/⟨Ḃ²⟩ = {ohmic:.4} ({:+.2} %); (d) lock-in over t = {f0:.1}–{f1:.1}: \
         circulation {lock_l:.3}, power {lock_p:.4}: L/P {:.4} against −2/Ḃ_rms = {lock_ohmic:.4} \
         ({:+.2} %); (e) j_x {mean_value:.4}: scatter across cells {scatter:.4}, reported error \
         {reported:.4} (ratio {:.2})",
        100.0 * (mapped / heat - 1.0),
        100.0 * (carried / current - 1.0),
        circulation / power,
        100.0 * (circulation / power / ohmic - 1.0),
        lock_l / lock_p,
        100.0 * (lock_l / lock_p / lock_ohmic - 1.0),
        scatter / reported
    );
    assert!((mapped / heat - 1.0).abs() < 0.01, "(a) {mapped} {heat}");
    assert!(
        (carried / current - 1.0).abs() < 0.03,
        "(b) {carried} {current}"
    );
    assert!(
        rows.iter().all(|j| j.abs() < 0.15 * jx.abs()),
        "(b) j_y by row {rows:?} against j_x {jx}"
    );
    assert!(
        (circulation / power / ohmic - 1.0).abs() < 0.10,
        "(c) {circulation} {power} {ohmic}"
    );
    assert!(
        (lock_l / lock_p / lock_ohmic - 1.0).abs() < 0.10,
        "(d) {lock_l} {lock_p} {lock_ohmic}"
    );
    assert!(
        (0.6..1.6).contains(&(scatter / reported)),
        "(e) {scatter} {reported}"
    );
}

/// The energy a circuit gives a metal (PHYSICS.md §2.11; §9, energy): on Ohm's law's
/// reference over t = 10–40 the supply's power `IV` against the field's work on the carriers
/// per time (`TubeRun::work`), within 5 % (set from the work's thermal noise over the window,
/// ~1.5 %); the lattice's heat, the rest of it (D10: the carriers' kinetic energy carried into
/// the contacts), printed. Measured −4.4 %. The account closed (physics test
/// `d11_circuit_energy_account`): `IV = W + S − H`, the contacts' potential steps `S`
/// (`TubeSim::book_potential_steps`) and the step's own remainder `H`, which must not depend
/// on the current: at the reference and with the supply at 0, over t = 10–40 in 6 blocks of
/// 5, within 3 standard errors (their scatter). Both are mostly the carriers' exchange with
/// the contacts' layers, so the gap `W − IV` is the field's work with no current at all.
#[test]
fn metal_energy_balance() {
    let level = level::shipped("ohms_law");
    let r = level.tube_run(&level.reference_solution, 1, 41.0);
    let volts = level.reference_solution[0].value;
    let current = r.currents[0].expect("the window has passed");
    let at = |series: &[(f64, f64)], t: f64| {
        let k = series
            .partition_point(|p| p.0 <= t)
            .max(1)
            .min(series.len() - 1);
        let (a, b) = (series[k - 1], series[k]);
        a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0)
    };
    let work = (at(&r.work[0], 40.0) - at(&r.work[0], 10.0)) / 30.0;
    let heat_series: Vec<(f64, f64)> = r.frames.iter().map(|f| (f.t, f.heat)).collect();
    let heat = (at(&heat_series, 40.0) - at(&heat_series, 10.0)) / 30.0;
    let supplied = volts * current;
    println!(
        "metal energy: the supply gives IV = {volts} × {current:.4} = {supplied:.3}; the field \
         works on the carriers at {work:.3} ({:+.2} %); the lattice takes {heat:.3} ({:.1} %)",
        100.0 * (work / supplied - 1.0),
        100.0 * heat / supplied
    );
    assert!((work / supplied - 1.0).abs() < 0.05, "{work} {supplied}");
    // The account closed (test `d11_circuit_energy_account`): `IV = W + S − H`, with the
    // contacts' potential steps `S` and the step's own remainder `H`, which does not depend
    // on the current: at the reference and with the supply at 0, over t = 10–40 in 6 blocks.
    let account = |v: f64| {
        let mut player = level.reference_solution.clone();
        player[0].value = v;
        let mut sim = level.tube_sim(&player, 1);
        sim.book_potential_steps();
        sim.advance_to(41.0);
        let r = sim.take_chunk();
        let rate = |s: &[(f64, f64)], t0: f64, t1: f64| (at(s, t1) - at(s, t0)) / (t1 - t0);
        let blocks: Vec<f64> = (0..6)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let t0 = 10.0 + 5.0 * k as f64;
                let t1 = t0 + 5.0;
                rate(&r.work[0], t0, t1) + rate(&r.potential_steps, t0, t1)
                    - v * rate(&r.collected[0], t0, t1)
            })
            .collect();
        let rest = blocks.iter().sum::<f64>() / 6.0;
        let spread = (blocks.iter().map(|b| (b - rest).powi(2)).sum::<f64>() / 5.0).sqrt();
        (
            rest,
            spread / 6.0_f64.sqrt(),
            rate(&r.potential_steps, 10.0, 40.0),
        )
    };
    let ((rest, err, steps), (idle, idle_err, idle_steps)) = std::thread::scope(|scope| {
        let on = scope.spawn(|| account(volts));
        let off = account(0.0);
        (on.join().expect("the run"), off)
    });
    println!(
        "metal energy: the contacts' potential steps take {steps:.3} ({:.1} % of IV), the          remainder W + S − IV is {rest:+.3} ± {err:.3} ({:.1} %); with the supply at 0, steps          {idle_steps:+.3}, remainder {idle:+.3} ± {idle_err:.3}: {:.1} standard errors apart",
        100.0 * steps / supplied,
        100.0 * rest / supplied,
        (rest - idle).abs() / libm::hypot(err, idle_err)
    );
    assert!(
        (rest - idle).abs() < 3.0 * libm::hypot(err, idle_err),
        "{rest} {err} {idle} {idle_err}"
    );
}
