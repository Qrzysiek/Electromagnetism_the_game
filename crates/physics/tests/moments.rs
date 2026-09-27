//! Particles with a magnetic moment (PHYSICS.md §3.2): gradients of B_z, the force
//! m ∇B_z against an analytic impulse, and the conservation laws of the exact
//! (relativistic) dynamics.
//!
//! `cargo test --release -p physics --test moments -- --nocapture`

use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{FieldSolver, LevelField};
use physics::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};
use physics::trajectory::{RunSettings, Scenario, StepView, run, run_observed};

const TOL: f64 = 1e-12;

/// Distance from `x` to the nearest source (wire or magnet centre), for the step of the
/// finite difference (whose error grows like (h / d)⁴ near a source).
fn source_distance(field: &LevelField, x: DVec3) -> f64 {
    let seg = |p: DVec3, a: DVec3, b: DVec3| {
        let ab = b - a;
        let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
        (p - (a + ab * t)).length()
    };
    let mut d = f64::INFINITY;
    for m in &field.dipoles {
        d = d.min((x - m.position).length());
    }
    for l in &field.loops {
        d = d.min(((x - l.center).length() - l.radius).abs());
    }
    for p in &field.polygons {
        let n = p.vertices.len();
        for i in 0..n {
            d = d.min(seg(x, p.vertices[i], p.vertices[(i + 1) % n]));
        }
    }
    d
}

/// 4th-order central difference of B_z along `dir`.
fn fd_bz(field: &LevelField, x: DVec3, dir: DVec3, h: f64) -> f64 {
    let bz = |d: f64| field.sample(x + dir * d, 0.0).b.z;
    (8.0 * (bz(h) - bz(-h)) - (bz(2.0 * h) - bz(-2.0 * h))) / (12.0 * h)
}

fn sources() -> Vec<(&'static str, LevelField)> {
    vec![
        (
            "magnet",
            LevelField {
                dipoles: vec![MagneticDipole {
                    position: DVec3::new(1.0, -0.5, 0.0),
                    moment: DVec3::new(0.0, 0.0, 2.0),
                    radius: 0.3,
                }],
                ..Default::default()
            },
        ),
        (
            "circular coil",
            LevelField {
                loops: vec![CircularLoop {
                    center: DVec3::ZERO,
                    normal: DVec3::Z,
                    radius: 3.0,
                    kappa: 0.7,
                    wire_radius: 0.1,
                }],
                ..Default::default()
            },
        ),
        (
            "polygon coil",
            LevelField {
                polygons: vec![PolygonCoil {
                    vertices: vec![
                        DVec3::new(-2.0, -1.0, 0.0),
                        DVec3::new(3.0, -1.5, 0.0),
                        DVec3::new(2.0, 2.5, 0.0),
                        DVec3::new(-1.5, 2.0, 0.0),
                    ],
                    kappa: -0.4,
                    wire_radius: 0.1,
                }],
                ..Default::default()
            },
        ),
    ]
}

/// G1: the exact gradients (curl-free identity for coils, direct for the dipole) agree
/// with a 4th-order finite difference of B_z to the difference's own accuracy, inside
/// and outside the coils, near and far (down to 0.03 cells from a wire; the difference
/// step scales with the distance to the nearest source), and on the coil axis.
#[test]
fn g1_gradients_of_bz_match_finite_differences() {
    let points = [
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1e-3, 0.0, 0.0),
        DVec3::new(0.4, 0.3, 0.0),
        DVec3::new(-1.7, 0.9, 0.0),
        DVec3::new(2.5, 0.2, 0.0),
        DVec3::new(4.0, -3.0, 0.0),
        DVec3::new(-9.0, 7.0, 0.0),
    ];
    for (name, field) in sources() {
        let mut worst: f64 = 0.0;
        for x in points {
            // Skip points inside a magnet (not physical there).
            if field
                .dipoles
                .iter()
                .any(|d| (x - d.position).length() < 0.5)
            {
                continue;
            }
            let g = field.grad_bz(x, 0.0);
            let h = 1e-3 * source_distance(&field, x).min(1.0);
            let fd = DVec3::new(
                fd_bz(&field, x, DVec3::X, h),
                fd_bz(&field, x, DVec3::Y, h),
                0.0,
            );
            let scale = g.length().max(fd.length()).max(1e-300);
            let err = (g - fd).length() / scale;
            if std::env::var("G1_DEBUG").is_ok() {
                println!("  {name} at {x}: exact {g}, fd {fd}, rel {err:.1e}");
            }
            worst = worst.max(err);
            assert!(g.z.abs() <= 1e-15 * scale, "{name}: ∂B_z/∂z ≠ 0 at {x}");
        }
        println!("G1 {name}: max relative difference from 4th-order FD {worst:.1e}");
        assert!(worst < 1e-9, "{name}: {worst:.3e}");
    }
}

fn dipole_field(mu: f64) -> LevelField {
    LevelField {
        dipoles: vec![MagneticDipole {
            position: DVec3::ZERO,
            moment: DVec3::new(0.0, 0.0, mu),
            radius: 0.3,
        }],
        ..Default::default()
    }
}

fn neutral(moment: f64) -> Particle {
    Particle {
        charge: 0.0,
        mass: 1.0,
        radius: 0.0,
        moment,
    }
}

/// S1: a neutral moment m flying straight past a magnet μ at impact parameter b picks up
/// the transverse momentum ∫ m ∂_y B_z dt = 4 m μ / (v b³) (impulse approximation, exact
/// to first order in m). The odd part (Δp(m) − Δp(−m)) / 2 cancels the second-order
/// term, so it matches to O(m²); spin up and spin down are kicked in opposite directions.
#[test]
fn s1_stern_gerlach_impulse() {
    let (mu, b, v, m) = (1.0, 2.0, 1.0, 1e-4);
    let kick = |m: f64| {
        let scn = Scenario {
            field: dipole_field(mu),
            obstacles: vec![],
            particle: neutral(m),
            c: f64::INFINITY,
            x0: DVec3::new(-300.0, b, 0.0),
            p0: DVec3::new(v, 0.0, 0.0),
            detector: None,
            bounds: None,
            t_max: 600.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        run(&scn, &RunSettings::with_tolerance(TOL)).end.p.y
    };
    let (up, down) = (kick(m), kick(-m));
    let expected = 4.0 * m * mu / (v * b.powi(3));
    let odd = 0.5 * (up - down);
    let rel = (odd - expected).abs() / expected;
    println!(
        "S1: Δp_y(+m) = {up:.9e}, Δp_y(−m) = {down:.9e}, odd part {odd:.9e}, impulse approximation {expected:.9e}, relative difference {rel:.1e}"
    );
    assert!(up > 0.0 && down < 0.0, "opposite kicks");
    // Corrections: O(m²) and the far tails (≈ 1e-9); integration error ≈ 1e-12 / Δp.
    assert!(rel < 1e-6, "{rel:.3e}");
}

/// Energy `(γ−1)mc² − m B_z` and angular momentum `|x × p|` about the magnet, along a
/// scattering orbit in its central potential U = m μ / r³ (attractive for m μ < 0).
fn central_run(c: f64, m: f64, speed_over_c: f64, b: f64) -> (f64, f64) {
    let mu = 1.0;
    let kin = Kinematics::new(1.0, c);
    let v = if c.is_finite() { speed_over_c * c } else { 1.0 };
    let gamma = if c.is_finite() {
        1.0 / (1.0 - speed_over_c * speed_over_c).sqrt()
    } else {
        1.0
    };
    let scn = Scenario {
        field: dipole_field(mu),
        obstacles: vec![],
        particle: neutral(m),
        c,
        x0: DVec3::new(-20.0, b, 0.0),
        p0: DVec3::new(gamma * v, 0.0, 0.0),
        detector: None,
        bounds: None,
        t_max: 40.0 / v,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let energy = |x: DVec3, p: DVec3| kin.kinetic_energy(p) - m * scn.field.sample(x, 0.0).b.z;
    let (e0, l0) = (energy(scn.x0, scn.p0), scn.x0.cross(scn.p0).z);
    let (mut de, mut dl): (f64, f64) = (0.0, 0.0);
    let mut closest = f64::INFINITY;
    run_observed(
        &scn,
        &RunSettings::with_tolerance(TOL),
        |step: &StepView<'_, LevelField>| {
            let (x, p) = step.state(step.t_end());
            de = de.max((energy(x, p) - e0).abs() / e0);
            dl = dl.max((x.cross(p).z - l0).abs() / l0.abs());
            closest = closest.min(x.length());
        },
    );
    println!(
        "  c = {c}, m μ = {m:+}: closest approach {closest:.3}, energy drift {de:.1e}, \
         angular momentum drift {dl:.1e}"
    );
    (de, dl)
}

/// S2: Newtonian central problem, attractive and repulsive (strong enough to deflect
/// the orbit substantially; the impact parameter keeps a centrifugal barrier against
/// the attractive 1/r³ potential).
#[test]
fn s2_newtonian_conservation() {
    println!("S2");
    for m in [0.5, -0.5] {
        let (de, dl) = central_run(f64::INFINITY, m, 0.0, 2.0);
        assert!(de < 1e-10 && dl < 1e-10);
    }
}

/// S3: relativistic central problem at v = 0.8 c. For E = 0 the moment's interaction
/// Lagrangian m B′_z/γ = m B_z is velocity independent (PHYSICS.md §3.2), so γ m c² − m B_z
/// and γ m r² φ̇ are exactly conserved.
#[test]
fn s3_relativistic_conservation() {
    println!("S3");
    for m in [30.0, -30.0] {
        let (de, dl) = central_run(5.0, m, 0.8, 4.0);
        assert!(de < 1e-10 && dl < 1e-10);
    }
}

/// S4: a charged particle with a moment in a coil's field: the magnetic Lorentz force does
/// no work, so `(γ−1)mc² − m B_z` is conserved, through many gyrations.
#[test]
fn s4_charged_moment_in_coil_field() {
    let field = LevelField {
        loops: vec![CircularLoop {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 6.0,
            kappa: 2.0,
            wire_radius: 0.1,
        }],
        ..Default::default()
    };
    let particle = Particle {
        charge: 1.0,
        mass: 1.0,
        radius: 0.0,
        moment: 0.3,
    };
    let c = 5.0;
    let kin = Kinematics::new(1.0, c);
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle,
        c,
        x0: DVec3::new(1.0, 0.5, 0.0),
        p0: DVec3::new(0.0, 1.5, 0.0),
        detector: None,
        bounds: None,
        t_max: 60.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let energy = |x: DVec3, p: DVec3| kin.kinetic_energy(p) - 0.3 * scn.field.sample(x, 0.0).b.z;
    let e0 = energy(scn.x0, scn.p0);
    let mut de: f64 = 0.0;
    let tr = run_observed(
        &scn,
        &RunSettings::with_tolerance(TOL),
        |step: &StepView<'_, LevelField>| {
            let (x, p) = step.state(step.t_end());
            de = de.max((energy(x, p) - e0).abs() / e0);
        },
    );
    // The trajectory's own energy diagnostic includes −m B_z too.
    let own = tr.energy_max_abs_error / e0.abs();
    println!(
        "S4: energy drift {de:.1e} (trajectory diagnostic {own:.1e}), {} steps",
        tr.stats.n_step
    );
    assert!(de < 1e-10 && own < 1e-10);
}
