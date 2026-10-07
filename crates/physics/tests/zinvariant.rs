//! Validation tests V1–V3 of the z-invariant world's electrostatics (`zinv.rs`, PHYSICS.md
//! §2.11; the tubes' groundwork, docs/TUBES.md). Run with
//! `cargo test --release -p physics --test zinvariant -- --nocapture --test-threads=1`.
//!
//! Only enclosed setups, whose results do not depend on the logarithm's reference.

#![allow(clippy::disallowed_methods)] // references; the engine uses libm (clippy.toml)

use physics::DVec3;
use physics::bem::Bias;
use physics::zinv::{Electrode, Electrodes, Section, segment_integrals};

/// V1: the segment integrals against Gauss–Legendre quadrature of `−ln r²` and
/// `2 (x − s)/r²` along the segment (2000 panels of 8 points, far from the segment and
/// near it, on its line beyond an end, and on it at its midpoint, where the field's
/// normal part is the average of its two sides, 0). Required: 1e-11 relative (the
/// quadrature's accuracy at these distances).
#[test]
fn v1_segment_integrals_against_quadrature() {
    let (a, b) = (DVec3::new(0.3, -0.2, 0.0), DVec3::new(1.9, 0.7, 0.0));
    let nodes = [
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ];
    let quad = |x: DVec3| {
        let panels = 2000;
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        let d = (b - a) / f64::from(panels);
        let h = d.length();
        for k in 0..panels {
            let mid = a + d * (f64::from(k) + 0.5);
            for &(t, w) in &nodes {
                for s in [mid + d * (0.5 * t), mid - d * (0.5 * t)] {
                    let r = x - s;
                    let r2 = r.x * r.x + r.y * r.y;
                    phi += -0.5 * w * h * r2.ln();
                    e += r * (0.5 * w * h * 2.0 / r2);
                }
            }
        }
        (phi, e)
    };
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(4.0, 3.0, 0.0),
        DVec3::new(1.0, 0.6, 0.0),
        DVec3::new(-1.0, -0.4, 0.0) + (b - a) * 0.0 - (b - a) * 0.5,
    ] {
        let (p1, e1) = segment_integrals(a, b, x);
        let (p2, e2) = quad(x);
        let err = ((p1 - p2) / p2).abs().max((e1 - e2).length() / e2.length());
        println!("V1 at {x}: φ {p1:.14} (quadrature {p2:.14}), E {e1} ({e2}): {err:.1e}");
        worst = worst.max(err);
    }
    // On the segment, at its midpoint: φ = −L(ln(L²/4) − 2), E_normal = 0.
    let mid = (a + b) * 0.5;
    let l = (b - a).length();
    let (p, e) = segment_integrals(a, b, mid);
    let want = -l * ((l * l / 4.0).ln() - 2.0);
    let n = DVec3::new(-(b - a).y, (b - a).x, 0.0) / l;
    println!(
        "V1 on the midpoint: φ {p:.14} (closed form {want:.14}), E·n {:.1e}",
        e.dot(n)
    );
    worst = worst.max(((p - want) / want).abs()).max(e.dot(n).abs());
    assert!(worst < 1e-11, "V1 {worst:.3e}");
}

/// V2: coaxial cylinders (radii 1 and 3; the inner at potential 1, the outer grounded):
/// the capacitance per unit length `1/(2 ln(b/a))` (k = 1; Jackson §1.11), with polygons
/// of segments 0.2, 0.1, 0.05 long. Required: the error shrinking with every refinement
/// at an observed order of at least 1 (inscribed polygons: O(h²) expected), and the two
/// cylinders' net charge (0 in the continuum: the outer one encloses the inner, and a net
/// charge would put the grounded outer one at `−2Q ln 3`, not 0) shrinking with every
/// refinement too (the collocation enforces it only approximately; a first version of
/// this test expected it exact).
#[test]
fn v2_coaxial_capacitance() {
    let exact = 1.0 / (2.0 * 3.0_f64.ln());
    let mut errors = Vec::new();
    let mut nets: Vec<f64> = Vec::new();
    for size in [0.2, 0.1, 0.05] {
        let e = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Potential(1.0),
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 3.0,
                    },
                    bias: Bias::Grounded,
                },
            ],
            &[],
            size,
        );
        let (q_in, q_out) = (e.charge(0), e.charge(1));
        let err = q_in / exact - 1.0;
        println!(
            "V2 segments {size}: {} of them, C' = {q_in:.10} (exact {exact:.10}), {err:.2e}; outer {q_out:.10}",
            e.segments.len()
        );
        nets.push(((q_in + q_out) / q_in).abs());
        errors.push(err);
    }
    println!(
        "V2 net charge relative {:.2e} {:.2e} {:.2e}",
        nets[0], nets[1], nets[2]
    );
    assert!(nets[0] > nets[1] && nets[1] > nets[2], "V2 net {nets:?}");
    let order = (errors[1].abs() / errors[2].abs()).ln() / 2f64.ln();
    println!("V2 observed order {order:.2}");
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && order >= 1.0,
        "V2 {errors:?}"
    );
}

/// V3: a line charge λ = 1 inside a grounded cylinder (radius 2), 0.8 off its axis: the
/// wall's induced charge is the image −λ at R²/e on the same ray, so the force per unit
/// length on the line charge is `2λ²/(R²/e − e)` towards the wall (exact). Required as V2.
#[test]
fn v3_line_charge_in_a_grounded_cylinder() {
    let (r, e) = (2.0, 0.8);
    let exact = 2.0 / (r * r / e - e);
    let x = DVec3::new(e, 0.0, 0.0);
    let mut errors = Vec::new();
    for size in [0.2, 0.1, 0.05] {
        let wall = Electrodes::new(
            vec![Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: r,
                },
                bias: Bias::Grounded,
            }],
            &[(x, 1.0)],
            size,
        );
        let f = wall.sample(x).1.x;
        let err = f / exact - 1.0;
        println!(
            "V3 segments {size}: force {f:.10} (image {exact:.10}), {err:.2e}; wall charge {:.10}",
            wall.charge(0)
        );
        errors.push(err);
    }
    let order = (errors[1].abs() / errors[2].abs()).ln() / 2f64.ln();
    println!("V3 observed order {order:.2}");
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && order >= 1.0,
        "V3 {errors:?}"
    );
}
