//! Validation tests N1–N3 for dielectric boxes (BEM, PHYSICS.md §2.7). Run with
//! `cargo test --release -p physics --test dielectrics -- --nocapture --test-threads=1`.
//!
//! The references: the polarizability of a cube, `α = p / (|V| ε₀ E₀)`, from Helsing &
//! Perfekt, "On the polarizability and capacitance of the cube", Appl. Comput. Harmon.
//! Anal. 34 (2013), arXiv:1203.5997, Table 1: 3.644305190268 for a conducting cube
//! (ε → ∞, relative error 1e-11) and −1.638415712936517 for ε = 0 (1e-14), the two ends
//! of the dielectric equation's range λ = (ε − 1)/(ε + 1) ∈ [−1, 1].

#![allow(clippy::disallowed_methods)] // references; the engine uses libm (clippy.toml)

use physics::DVec3;
use physics::bem::{Bias, BodyShape, BoxElectrode, Dielectric, Electrodes};

const ALPHA_CONDUCTOR: f64 = 3.644_305_190_268;
const ALPHA_ZERO: f64 = -1.638_415_712_936_517;
/// Panel sizes: preview (0.5) and verify (0.35) resolutions and a finer one.
const SIZES: [f64; 3] = [0.5, 0.35, 0.25];

/// A cube of side 2 at the origin.
fn cube() -> BoxElectrode {
    BoxElectrode {
        center: DVec3::ZERO,
        angle: 0.0,
        half_length: 1.0,
        half_thickness: 1.0,
        half_height: 1.0,
        bias: Bias::Charge(0.0),
    }
}

/// The nearly uniform field of charges ±Q at ∓D on the x axis (its curvature over the
/// cube is of order (1/D)² = 1e-6), and its strength at the origin.
fn uniform() -> ([(DVec3, f64); 2], f64) {
    let (q, d) = (1e6, 1000.0);
    (
        [(DVec3::new(-d, 0.0, 0.0), q), (DVec3::new(d, 0.0, 0.0), -q)],
        2.0 * q / (d * d),
    )
}

fn dielectric_cube(eps: f64, sources: &[(DVec3, f64)], size: f64) -> Electrodes {
    Electrodes::with_dielectrics(
        Vec::new(),
        vec![Dielectric {
            shape: BodyShape::Box(cube()),
            permittivity: eps,
        }],
        sources,
        size,
    )
}

/// The cube's polarizability from its panels' dipole moment (each constant-density
/// triangle's moment is exactly σ A × centroid; both halves).
fn polarizability(e: &Electrodes, field: f64) -> f64 {
    let p: f64 = e
        .panels()
        .map(|(t, s)| 2.0 * s * t.area() * t.centroid().x)
        .sum();
    p / (8.0 * field / (4.0 * std::f64::consts::PI))
}

/// Observed order of convergence from errors at the sizes `SIZES[1..]`.
fn order(e1: f64, e2: f64) -> f64 {
    (e1.abs() / e2.abs()).ln() / (SIZES[1] / SIZES[2]).ln()
}

/// N1: ε = 1 is vacuum: no bound charge at all (exactly). For any ε each body's bound
/// charge sums to zero: required within 1e-12 of a nearby unit charge (the bordered
/// system's constraint; without it −0.6 % at ε = 2 up to −48 % at ε = 1e4).
#[test]
fn n1_vacuum_and_neutrality() {
    let q = [(DVec3::new(0.3, 2.5, 0.4), 1.0)];
    let vacuum = dielectric_cube(1.0, &q, 0.5);
    assert!(vacuum.sigma.iter().all(|&s| s == 0.0), "N1 ε = 1");
    let mut worst: f64 = 0.0;
    for eps in [0.0, 2.0, 4.0, 80.0, 1e12] {
        let e = dielectric_cube(eps, &q, 0.35);
        let net: f64 = e.panels().map(|(t, s)| 2.0 * s * t.area()).sum();
        println!("N1 ε = {eps:e}: net bound charge {net:.2e}");
        worst = worst.max(net.abs());
    }
    assert!(worst < 1e-12, "N1 {worst:.3e}");
}

/// N2: the cube's polarizability at both ends of the dielectric equation, ε → ∞ (1e12)
/// and ε = 0, against Helsing & Perfekt; and the conducting cube by the electrodes'
/// own method (first kind) for comparison. Required: the error shrinks with every
/// refinement, at an observed order of at least 1 (the electrodes' E1 shows ~2).
/// Printed: the errors at preview and verify resolution, the dielectric's measured
/// accuracy (PHYSICS.md §2.7).
#[test]
fn n2_cube_polarizability_at_the_limits() {
    let (sources, field) = uniform();
    let mut errors: Vec<[f64; 3]> = Vec::new();
    for size in SIZES {
        let metal = Electrodes::with_panel_size(vec![cube()], &sources, size);
        let e = [
            polarizability(&metal, field) / ALPHA_CONDUCTOR - 1.0,
            polarizability(&dielectric_cube(1e12, &sources, size), field) / ALPHA_CONDUCTOR - 1.0,
            polarizability(&dielectric_cube(0.0, &sources, size), field) / ALPHA_ZERO - 1.0,
        ];
        println!(
            "N2 panel size {size}: conductor (first kind) {:.2e}, dielectric ε → ∞ {:.2e}, ε = 0 {:.2e}",
            e[0], e[1], e[2]
        );
        errors.push(e);
    }
    for k in 0..3 {
        let o = order(errors[1][k], errors[2][k]);
        println!("N2 case {k}: observed order {o:.2}");
        assert!(
            errors[0][k].abs() > errors[1][k].abs() && errors[1][k].abs() > errors[2][k].abs(),
            "N2 case {k} not converging: {:?}",
            errors.iter().map(|e| e[k]).collect::<Vec<_>>()
        );
        assert!(o >= 1.0, "N2 case {k} order {o:.2}");
    }
}

/// N3: a cube of ε = 4 (glass), between the limits: the polarizability converges (the
/// differences of successive refinements shrink at an observed order of at least 1);
/// printed, the Richardson extrapolation and the error at verify resolution against it.
#[test]
fn n3_cube_of_glass_converges() {
    let (sources, field) = uniform();
    let a: Vec<f64> = SIZES
        .iter()
        .map(|&s| polarizability(&dielectric_cube(4.0, &sources, s), field))
        .collect();
    let (d1, d2) = (a[1] - a[0], a[2] - a[1]);
    let o = (d1.abs() / d2.abs()).ln() / (SIZES[1] / SIZES[2]).ln();
    let limit = a[2] + d2 / ((SIZES[1] / SIZES[2]).powf(o) - 1.0);
    println!(
        "N3 ε = 4: α = {:.8}, {:.8}, {:.8}; order {o:.2}, extrapolated {limit:.6}: verify \
         resolution {:.2e} off",
        a[0],
        a[1],
        a[2],
        a[1] / limit - 1.0
    );
    assert!(d1.abs() > d2.abs() && o >= 1.0, "N3 {a:?}");
}

fn dielectric_sphere(eps: f64, radius: f64, sources: &[(DVec3, f64)], size: f64) -> Electrodes {
    Electrodes::with_dielectrics(
        Vec::new(),
        vec![Dielectric {
            shape: BodyShape::Sphere {
                center: DVec3::ZERO,
                radius,
            },
            permittivity: eps,
        }],
        sources,
        size,
    )
}

/// N4: a dielectric sphere (radius 1) in a uniform field (as in N2): its polarizability
/// `p/(|V|ε₀E₀) = 3(ε − 1)/(ε + 2)` exactly (Jackson §4.4), at ε = 4, ε → ∞ and ε = 0.
/// Required as N2: the error shrinks with every refinement at an observed order of at
/// least 1. The sphere has no edges, so the method should do better than on the cube.
#[test]
fn n4_sphere_in_a_uniform_field() {
    let (sources, field) = uniform();
    let volume = 4.0 / 3.0 * std::f64::consts::PI;
    for eps in [4.0, 1e12, 0.0] {
        let exact = 3.0 * (eps - 1.0) / (eps + 2.0);
        let errors: Vec<f64> = SIZES
            .iter()
            .map(|&size| {
                let e = dielectric_sphere(eps, 1.0, &sources, size);
                let p: f64 = e
                    .panels()
                    .map(|(t, s)| 2.0 * s * t.area() * t.centroid().x)
                    .sum();
                let alpha = p / (volume * field / (4.0 * std::f64::consts::PI));
                alpha / exact - 1.0
            })
            .collect();
        let o = order(errors[1], errors[2]);
        println!(
            "N4 sphere ε = {eps:e}: errors {:.2e}, {:.2e}, {:.2e} (exact {exact:.6}); order {o:.2}",
            errors[0], errors[1], errors[2]
        );
        assert!(
            errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && o >= 1.0,
            "N4 ε = {eps}: {errors:?}"
        );
    }
}

/// N5: a point charge q = 1 at distance d = 1.6 from the centre of a dielectric sphere
/// (radius 1, ε = 4): the field of the sphere's bound charge at the charge against the
/// exact Legendre series (the boundary conditions give the outside coefficients
/// `B_l = −(ε − 1) l q a^{2l+1} / ((εl + l + 1) d^{l+1})`, so
/// `E_r(d) = −q Σ l(l + 1)(ε − 1) a^{2l+1} / ((εl + l + 1) d^{2l+3})`), summed to 1e-16.
/// Required as N2: converging with every refinement at an observed order of at least 1.
#[test]
fn n5_charge_near_a_dielectric_sphere() {
    let (eps, a, d) = (4.0_f64, 1.0_f64, 1.6_f64);
    let mut exact = 0.0;
    for l in 1..400 {
        let lf = f64::from(l);
        let term = lf * (lf + 1.0) * (eps - 1.0) * a.powi(2 * l + 1)
            / ((eps * lf + lf + 1.0) * d.powi(2 * l + 3));
        exact -= term;
        if term < 1e-17 {
            break;
        }
    }
    let x = DVec3::new(d, 0.0, 0.0);
    let errors: Vec<f64> = SIZES
        .iter()
        .map(|&size| {
            use physics::field::FieldSolver;
            let e = dielectric_sphere(eps, a, &[(x, 1.0)], size);
            e.sample(x, 0.0).e.x / exact - 1.0
        })
        .collect();
    let o = order(errors[1], errors[2]);
    println!(
        "N5 charge 0.6 from a glass sphere: E = {exact:.10} (series); errors {:.2e}, {:.2e}, {:.2e}; order {o:.2}",
        errors[0], errors[1], errors[2]
    );
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && o >= 1.0,
        "N5 {errors:?}"
    );
}
