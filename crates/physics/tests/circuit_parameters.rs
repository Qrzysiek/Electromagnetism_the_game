//! Circuit parameters from the fields (PHYSICS.md §2.8; Jackson §1.11–1.13, §5.17): the
//! groundwork for coupling lumped circuits to the electrodes and coils. Run with
//! `cargo test --release -p physics --test circuit_parameters -- --nocapture`.
//!
//! Inductances are in units of `μ₀/4π` (the game's `1/c²`): the Neumann integral
//! `N = ∮∮ dl₁·dl₂/|x₁ − x₂|`, which a coil's vector potential per unit strength
//! (`κ = μ₀I/4π = 1`) gives as `∮ A·dl` around the other coil. Circular coils through the
//! engine's own functions (`physics::inductance`), which driven coils use.

#![allow(clippy::disallowed_methods)] // references; the engine's fields use libm (clippy.toml)
#![allow(clippy::cast_precision_loss)] // small loop counters

mod common;

use std::f64::consts::PI;

use physics::DVec3;
use physics::inductance::{around_circle, mutual_circles, ring_self};
use physics::magnetic::{CircularLoop, PolygonCoil};

/// A coil of the game's slice: in the plane z = 0, current counter-clockwise.
fn ring(center: DVec3, radius: f64, wire_radius: f64) -> CircularLoop {
    CircularLoop {
        center,
        normal: DVec3::Z,
        radius,
        kappa: 1.0,
        wire_radius,
        rate: 0.0,
    }
}

/// `∮ A·dl` along a closed polygon (Gauss–Legendre on each side).
fn around_polygon(a: impl Fn(DVec3) -> DVec3, vertices: &[DVec3], n: usize) -> f64 {
    let rule = common::gauss_legendre(n);
    let mut sum = 0.0;
    for (i, &p) in vertices.iter().enumerate() {
        let q = vertices[(i + 1) % vertices.len()];
        let d = q - p;
        for &(x, w) in &rule {
            sum += 0.5 * w * a(p + d * (0.5 * (1.0 + x))).dot(d);
        }
    }
    sum
}

/// Z2, Jackson §5.17 (5.155–5.156): the mutual inductance of two coils is the flux of one
/// linked by the other, `∮ A₁·dl₂`, Neumann's symmetric double integral. Coplanar coils,
/// as in the game's slice: (a) concentric loops (radii 3, 5): Jackson Pr. 5.28's closed
/// form with d = 0; (b) identical loops side by side (radius 1, centres 2.5, 5, 10 apart),
/// with Pr. 5.34(c)'s series `−(μ₀πa/4)[(a/R)³ + (9/4)(a/R)⁵ + (375/64)(a/R)⁷ + …]`
/// printed; (c) a circle and a square coil, each linking the other's flux. References:
/// the Neumann double integral taken directly (`scripts/wolfram/z2_inductances.wls`,
/// Wolfram Engine 14.2, 30 digits; for (a) it equals (5.28) to 20 digits). Required:
/// each way within 1e-12 of the reference, and in (c) the two ways within 1e-12 of each
/// other (the formulas are exact).
#[test]
fn z2_mutual_inductance_of_coplanar_coils() {
    let mut worst: f64 = 0.0;
    // (a) Concentric, radii 3 and 5.
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const CONCENTRIC: f64 = 41.804596452269273977;
    let (inner, outer) = (ring(DVec3::ZERO, 3.0, 0.1), ring(DVec3::ZERO, 5.0, 0.1));
    let n21 = mutual_circles(&inner, &outer);
    let n12 = mutual_circles(&outer, &inner);
    for (name, v) in [
        ("inner's flux in the outer", n21),
        ("outer's in the inner", n12),
    ] {
        let e = v / CONCENTRIC - 1.0;
        println!(
            "Z2 (a) concentric 3, 5, {name}: {v:.15e} (Jackson (5.28): {CONCENTRIC:.15e}), {e:.1e}"
        );
        worst = worst.max(e.abs());
    }
    // (b) Identical, side by side.
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const SIDE: [(f64, f64); 3] = [
        (2.5, -1.0411745088500855552),
        (5.0, -0.086899585705737451942),
        (10.0, -0.010097623953420449727),
    ];
    for (r, reference) in SIDE {
        let v = mutual_circles(
            &ring(DVec3::ZERO, 1.0, 0.1),
            &ring(DVec3::new(r, 0.0, 0.0), 1.0, 0.1),
        );
        let e = v / reference - 1.0;
        let u = 1.0 / r;
        let series = -PI * PI * u.powi(3) * (1.0 + 2.25 * u * u + 375.0 / 64.0 * u.powi(4));
        println!(
            "Z2 (b) identical, R/a = {r}: {v:.15e} (Neumann {reference:.15e}), {e:.1e}; \
             Pr. 5.34(c) series {:.1e} off",
            series / reference - 1.0
        );
        worst = worst.max(e.abs());
    }
    // (c) A circle (radius 1, at the origin) and a square coil (side 3, centred at (3.5, 0)).
    let circle = ring(DVec3::ZERO, 1.0, 0.1);
    let square = PolygonCoil {
        vertices: vec![
            DVec3::new(2.0, -1.5, 0.0),
            DVec3::new(5.0, -1.5, 0.0),
            DVec3::new(5.0, 1.5, 0.0),
            DVec3::new(2.0, 1.5, 0.0),
        ],
        kappa: 1.0,
        wire_radius: 0.1,
        rate: 0.0,
    };
    let in_square = around_polygon(|x| circle.unit_vector_potential(x), &square.vertices, 48);
    let in_circle = around_circle(
        |x| square.unit_vector_potential(x),
        DVec3::ZERO,
        DVec3::Z,
        1.0,
        256,
    );
    let e = in_square / in_circle - 1.0;
    println!(
        "Z2 (c) circle and square: the circle's flux in the square {in_square:.15e}, the \
         square's in the circle {in_circle:.15e}: {e:.1e}"
    );
    worst = worst.max(e.abs());
    assert!(worst < 1e-12, "Z2 {worst:.3e}");
}

/// Z3, Jackson §5.17 B and Pr. 5.32: the self-inductance of a ring (radius a) of thin wire
/// (radius b). A current spread uniformly over the wire's surface links the ring's own
/// flux through the circles on that surface, `(ρ, z) = (a + b cos ψ, b sin ψ)`, averaged
/// over ψ. Jackson: `L = μ₀a[ln(8a/b) − 2]` for a surface current, `− 7/4` for a uniform
/// one (the interior adds μ₀a/4); "are the corrections of order b/a or (b/a)²?" The
/// engine's flux against the same average taken by Wolfram Engine 14.2 at 30 digits
/// (`scripts/wolfram/z2_inductances.wls`) for b/a = 1/10 … 1/1000: within 1e-12 (exact
/// formulas). Printed: the corrections to Jackson's, which are of order (b/a)²: their
/// coefficient tends to `(1/8) ln(8a/b) + 1/16` (exact: the remainder is O((b/a)⁴ ln),
/// the Wolfram script shows).
#[test]
fn z3_self_inductance_of_a_ring() {
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const FLUX: [(f64, f64); 5] = [
        (10.0, 30.010140533393935053),
        (30.0, 43.749437269349276977),
        (100.0, 58.869695714544086944),
        (300.0, 72.674280904734077995),
        (1000.0, 87.803719706364956854),
    ];
    let a = 1.0;
    let mut worst: f64 = 0.0;
    for (ratio, reference) in FLUX {
        let b = a / ratio;
        let flux = ring_self(&ring(DVec3::ZERO, a, b));
        let e = flux / reference - 1.0;
        let jackson = 4.0 * PI * a * ((8.0 * a / b).ln() - 2.0);
        let coefficient = (flux - jackson) / (4.0 * PI * a * (b / a).powi(2));
        println!(
            "Z3 b/a = 1/{ratio}: L/(μ₀/4π) = {flux:.15e} (Wolfram {reference:.15e}, {e:.1e}); \
             Jackson's 4πa[ln(8a/b) − 2] = {jackson:.10e}, the rest / (4πa (b/a)²) = \
             {coefficient:.6} ((1/8) ln(8a/b) + 1/16 = {:.6})",
            (8.0 * a / b).ln() / 8.0 + 1.0 / 16.0
        );
        worst = worst.max(e.abs());
    }
    assert!(worst < 1e-12, "Z3 {worst:.3e}");
}

/// Z1, Jackson §1.11 and Problems 1.12–1.13: Green's reciprocation for the electrodes (the
/// boundary-element method of §2.7, at the game's preview and verify resolutions). Two
/// unequal plates, facing across a gap of 1 (an 8 × 8 one, 0.25 thick, and a 6 × 6 one,
/// 0.5 thick, off-centre by 0.5), both grounded, and a charge q: (a) the charge it
/// induces on each, `Q_k = −q φ_k(x)`, φ_k the potential with plate k at 1 and the other
/// grounded (the Shockley–Ramo theorem: a moving charge drives the current
/// `q v·E_k` into plate k), at three places (between the plates, beside them, behind
/// one); (b) between the plates, Problem 1.13: `−q` times the fractional distance from
/// the other plate (infinite planes; the finite plates leak ~1e-4); (c) the coefficients
/// of capacitance are symmetric, `C_AB = C_BA`. Required, set from the method's measured
/// discretization error (E1, E2: 2.6e-3 to 6e-3; preview against verify 6.9e-4): each
/// within 5e-3 at verify, and smaller at verify than at preview.
#[test]
fn z1_green_reciprocation_for_electrodes() {
    use physics::bem::{Bias, BoxElectrode, Electrodes, Resolution};
    use physics::field::FieldSolver;
    let plate = |x: f64, y: f64, len: f64, thick: f64, height: f64, bias: Bias| BoxElectrode {
        center: DVec3::new(x, y, 0.0),
        angle: 0.0,
        half_length: len / 2.0,
        half_thickness: thick / 2.0,
        half_height: height / 2.0,
        bias,
    };
    let pair = |a: Bias, b: Bias| {
        vec![
            plate(0.0, -0.625, 8.0, 0.25, 8.0, a),
            plate(0.5, 0.75, 6.0, 0.5, 6.0, b),
        ]
    };
    let (q, s) = (1.0, 0.3);
    // Between the plates at fractional distance s from A (the inner faces at y = ±0.5).
    let between = DVec3::new(0.0, -0.5 + s, 0.0);
    let places = [
        ("between", between),
        ("beside", DVec3::new(6.0, 0.0, 0.0)),
        ("behind A", DVec3::new(-1.0, -2.0, 0.0)),
    ];
    let mut defects = Vec::new();
    for resolution in [Resolution::Preview, Resolution::Verify] {
        let unit_a = Electrodes::new(pair(Bias::Potential(1.0), Bias::Grounded), &[], resolution);
        let unit_b = Electrodes::new(pair(Bias::Grounded, Bias::Potential(1.0)), &[], resolution);
        let mut worst: f64 = 0.0;
        for (name, x) in places {
            let induced =
                Electrodes::new(pair(Bias::Grounded, Bias::Grounded), &[(x, q)], resolution);
            for (k, unit) in [(0, &unit_a), (1, &unit_b)] {
                let got = induced.charge(k);
                let expected = -q * unit.sample(x, 0.0).phi;
                let e = (got - expected).abs() / q;
                println!(
                    "Z1 (a) {resolution:?}, {name}: Q_{} = {got:.6e}, −q φ = {expected:.6e}: {e:.1e}",
                    ["A", "B"][k]
                );
                worst = worst.max(e);
            }
        }
        let induced = Electrodes::new(
            pair(Bias::Grounded, Bias::Grounded),
            &[(between, q)],
            resolution,
        );
        let pr113 = (induced.charge(0) + q * (1.0 - s)).abs() / q;
        let pr113_b = (induced.charge(1) + q * s).abs() / q;
        let (c_ab, c_ba) = (unit_b.charge(0), unit_a.charge(1));
        let sym = (c_ab - c_ba).abs() / c_ab.abs();
        println!(
            "Z1 {resolution:?}: (a) worst {worst:.1e}; (b) Pr. 1.13: Q_A + q(1 − s) = {:.1e} q, \
             Q_B + q s = {:.1e} q; (c) C_AB = {c_ab:.6e}, C_BA = {c_ba:.6e}: {sym:.1e} \
             (C_AA = {:.4}, C_BB = {:.4}; {} panels)",
            induced.charge(0) + q * (1.0 - s),
            induced.charge(1) + q * s,
            unit_a.charge(0),
            unit_b.charge(1),
            unit_a.sigma.len()
        );
        defects.push([worst, pr113.max(pr113_b), sym]);
    }
    let (preview, verify) = (defects[0], defects[1]);
    for (i, name) in ["(a)", "(b)", "(c)"].iter().enumerate() {
        assert!(
            verify[i] < 5e-3 && verify[i] < preview[i],
            "Z1 {name}: preview {:.2e}, verify {:.2e}",
            preview[i],
            verify[i]
        );
    }
}
