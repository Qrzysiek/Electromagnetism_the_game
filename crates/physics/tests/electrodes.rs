//! Validation tests E1–E4 for box electrodes (BEM, PHYSICS.md §2.7). Run with
//! `cargo test --release -p physics --test electrodes -- --nocapture --test-threads=1`.

mod common;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::bem::{Bias, BoxElectrode, Electrodes};
use physics::dynamics::Particle;
use physics::field::{FieldSolver, LevelField};
use physics::trajectory::{RunSettings, Scenario, run};

fn plate(
    x: f64,
    y: f64,
    angle: f64,
    len: f64,
    thick: f64,
    height: f64,
    bias: Bias,
) -> BoxElectrode {
    BoxElectrode {
        center: DVec3::new(x, y, 0.0),
        angle,
        half_length: len / 2.0,
        half_thickness: thick / 2.0,
        half_height: height / 2.0,
        bias,
    }
}

// --- E1: capacitance of the unit cube ------------------------------------------------------

/// The capacitance of a unit cube is C = 0.66067815 · 4πε₀ a (Hwang & Douglas 2004;
/// Read 1997: 0.660678; Mascagni & Simonov 2004: 0.6606785). In units k = 1 its charge at
/// potential 1 is 0.66067815. The BEM value must converge to it as the panels shrink.
#[test]
fn e1_unit_cube_capacitance() {
    let exact = 0.660_678_15;
    let mut previous = f64::INFINITY;
    let mut rows = Vec::new();
    for size in [0.25, 0.125, 0.0625] {
        let cube_e = plate(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, Bias::Potential(1.0));
        let e = Electrodes::with_panel_size(vec![cube_e], &[], size);
        let q = e.charge(0);
        let err = (q - exact).abs() / exact;
        println!(
            "E1 panel size {size}: {} panels, Q = {q:.8}, rel. error {err:.2e}",
            e.sigma.len()
        );
        assert!(err < previous, "no convergence");
        previous = err;
        rows.push(err);
    }
    let rate = (rows[1] / rows[2]).log2();
    println!("E1 observed order {rate:.2}");
    assert!(rows[2] < 2e-3, "finest error {:.3e}", rows[2]);
}

// --- E2: boundary condition between collocation points ----------------------------------

/// A charge near a grounded plate: the surface potential, sampled at points that are not
/// collocation points (panel vertices shifted inwards), is 0 up to the discretization
/// error, which must shrink with the panel size.
#[test]
fn e2_surface_residual_shrinks() {
    let src = [(DVec3::new(0.0, 2.0, 0.0), 1.0)];
    let mut previous = f64::INFINITY;
    for size in [0.5, 0.25, 0.125] {
        let p = plate(0.0, 0.0, 0.0, 4.0, 0.5, 3.0, Bias::Grounded);
        let e = Electrodes::with_panel_size(vec![p], &src, size);
        let mut worst: f64 = 0.0;
        for (t, _) in e.panels() {
            // A point inside the panel, away from the centroid: 1/2 of the way to a vertex.
            let x = (t.centroid() + t.a) * 0.5;
            let phi =
                e.sample(x, 0.0).phi + src.iter().map(|&(q, c)| c / (x - q).length()).sum::<f64>();
            worst = worst.max(phi.abs());
        }
        // Scale: the source potential on the plate (~1/2).
        let rel = worst / 0.5;
        println!("E2 panel size {size}: max |φ| off collocation points {rel:.2e} (relative)");
        assert!(rel < previous);
        previous = rel;
    }
    assert!(previous < 2e-2);
}

// --- E3: energy conservation ---------------------------------------------------------------

/// Two plates at ±V deflect a charged particle: W = T + qφ conserved (the field of a
/// fixed panel set is the exact field of a real charge distribution).
#[test]
fn e3_energy_conservation() {
    for c in [f64::INFINITY, 5.0] {
        let e = Electrodes::with_panel_size(
            vec![
                plate(0.0, 2.0, 0.0, 6.0, 0.4, 4.0, Bias::Potential(0.3)),
                plate(0.0, -2.0, 0.0, 6.0, 0.4, 4.0, Bias::Potential(-0.3)),
            ],
            &[],
            0.5,
        );
        let obstacles = e.obstacles(0.02);
        let field = LevelField {
            electrodes: e,
            ..LevelField::default()
        };
        let scn = Scenario {
            field,
            obstacles,
            particle: UNIT_PARTICLE,
            c,
            x0: DVec3::new(-8.0, 0.4, 0.0),
            p0: DVec3::new(1.0, 0.0, 0.0),
            detector: None,
            bounds: Some(cube(20.0)),
            t_max: 40.0,
            radiation_reaction: false,
            acceptance: None,
        };
        let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
        let rel = tr.energy_max_abs_error / tr.kinetic_initial;
        println!(
            "E3 c = {c}: {:?} after {} steps, max |ΔW|/T0 = {rel:.2e}",
            tr.outcome, tr.stats.n_accept
        );
        assert!(tr.stats.n_accept > 20);
        assert!(rel < 1e-10, "{rel:.3e}");
    }
}

// --- E4: the 2D slice ---------------------------------------------------------------------

/// Electrodes symmetric about z = 0: in the plane E_z = 0 exactly, and a particle starting
/// in the plane stays there bit for bit.
#[test]
fn e4_plane_symmetry() {
    let e = Electrodes::with_panel_size(
        vec![
            plate(1.0, 2.0, 0.4, 5.0, 0.4, 3.0, Bias::Potential(0.5)),
            plate(-2.0, -1.5, -0.3, 3.0, 0.6, 2.0, Bias::Charge(-0.4)),
        ],
        &[(DVec3::new(4.0, -3.0, 0.0), 0.7)],
        0.5,
    );
    for x in [DVec3::new(0.3, 0.1, 0.0), DVec3::new(5.0, -2.0, 0.0)] {
        assert_eq!(e.sample(x, 0.0).e.z.to_bits(), 0.0f64.to_bits());
    }
    let obstacles = e.obstacles(0.02);
    let field = LevelField {
        electrodes: e,
        ..LevelField::default()
    };
    let scn = Scenario {
        field,
        obstacles,
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
        },
        c: 3.0,
        x0: DVec3::new(-8.0, 1.0, 0.0),
        p0: DVec3::new(1.0, -0.2, 0.0),
        detector: None,
        bounds: Some(cube(20.0)),
        t_max: 40.0,
        radiation_reaction: false,
        acceptance: None,
    };
    let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
    println!("E4: {:?} after {} steps", tr.outcome, tr.stats.n_accept);
    assert!(tr.samples.iter().all(|s| s.x.z == 0.0 && s.p.z == 0.0));
}

/// Timing (information only): a deflector pair at preview resolution, a flight through it.
#[test]
fn e5_timing() {
    let e = Electrodes::with_panel_size(
        vec![
            plate(0.0, 3.0, 0.0, 8.0, 0.4, 4.0, Bias::Potential(0.05)),
            plate(0.0, -3.0, 0.0, 8.0, 0.4, 4.0, Bias::Potential(-0.05)),
        ],
        &[],
        0.5,
    );
    let panels = e.sigma.len();
    let obstacles = e.obstacles(0.02);
    let scn = Scenario {
        field: LevelField {
            electrodes: e,
            ..LevelField::default()
        },
        obstacles,
        particle: UNIT_PARTICLE,
        c: f64::INFINITY,
        x0: DVec3::new(-12.0, 0.0, 0.0),
        p0: DVec3::new(1.0, 0.0, 0.0),
        detector: None,
        bounds: Some(cube(15.0)),
        t_max: 40.0,
        radiation_reaction: false,
        acceptance: None,
    };
    let start = std::time::Instant::now();
    let tr = run(&scn, &RunSettings::with_tolerance(1e-10));
    println!(
        "E5: {panels} panels, {:?}, {} steps, {} rhs, {:.1} ms",
        tr.outcome,
        tr.stats.n_accept,
        tr.stats.n_fcn,
        start.elapsed().as_secs_f64() * 1e3
    );
}
