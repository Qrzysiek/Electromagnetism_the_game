//! Validation tests F1–F3 for ferrite bodies (high-μ insulators, BEM with a density odd
//! under the mirror z → −z; PHYSICS.md §2.7). Run with
//! `cargo test --release -p physics --test ferrites -- --nocapture --test-threads=1`.
//!
//! The magnetic problem is the dielectric one with B for E and μ for ε, so the exact
//! references of the dielectric tests apply: the cube's polarizability (Helsing &
//! Perfekt 2013, arXiv:1203.5997, Table 1: 3.644305190268 for μ → ∞, −1.638415712936517
//! for μ = 0), which is isotropic (a cube's polarizability tensor is a multiple of the
//! identity), so it holds for B along z, the direction the game's slice has; and the
//! sphere's `3(μ − 1)/(μ + 2)`. Here the field is along z and the bound charge odd under
//! the mirror: the code path of the slice, not the dielectrics' even one.

#![allow(clippy::disallowed_methods)] // references; the engine uses libm (clippy.toml)

use physics::DVec3;
use physics::bem::{Bias, BodyShape, BoxElectrode, Ferrite, Ferrites};

const ALPHA_CONDUCTOR: f64 = 3.644_305_190_268;
const ALPHA_ZERO: f64 = -1.638_415_712_936_517;
const SIZES: [f64; 3] = [0.5, 0.35, 0.25];
const B0: f64 = 1.0;

fn cube() -> BodyShape {
    BodyShape::Box(BoxElectrode {
        center: DVec3::ZERO,
        angle: 0.0,
        half_length: 1.0,
        half_thickness: 1.0,
        half_height: 1.0,
        bias: Bias::Charge(0.0),
    })
}

fn sphere() -> BodyShape {
    BodyShape::Sphere {
        center: DVec3::ZERO,
        radius: 1.0,
    }
}

fn magnetized(shape: BodyShape, mu: f64, size: f64) -> Ferrites {
    Ferrites::new(
        vec![Ferrite {
            shape,
            permeability: mu,
        }],
        |_| DVec3::new(0.0, 0.0, B0),
        size,
    )
}

/// The polarizability `m/(|V| B₀/4π)` from the magnetic dipole moment of the bound
/// charge: each upper panel σ A at its centroid, its mirror −σ A at the mirrored one,
/// so `m_z = 2 Σ σ A z`.
fn polarizability(f: &Ferrites, volume: f64) -> f64 {
    let m: f64 = f
        .panels()
        .map(|(t, s)| 2.0 * s * t.area() * t.centroid().z)
        .sum();
    m / (volume * B0 / (4.0 * std::f64::consts::PI))
}

fn order(e1: f64, e2: f64) -> f64 {
    (e1.abs() / e2.abs()).ln() / (SIZES[1] / SIZES[2]).ln()
}

fn converging(name: &str, errors: &[f64]) {
    let o = order(errors[1], errors[2]);
    println!(
        "{name}: errors {:.2e}, {:.2e}, {:.2e}; order {o:.2}",
        errors[0], errors[1], errors[2]
    );
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && o >= 1.0,
        "{name}: {errors:?}"
    );
}

/// F1: a ferrite cube (side 2) in a uniform B along z, μ → ∞ (1e12) and μ = 0, against
/// Helsing & Perfekt. Required as the dielectric N2: the error shrinking with every
/// refinement at an observed order of at least 1.
#[test]
fn f1_cube_along_z() {
    for (mu, exact) in [(1e12, ALPHA_CONDUCTOR), (0.0, ALPHA_ZERO)] {
        let errors: Vec<f64> = SIZES
            .iter()
            .map(|&s| polarizability(&magnetized(cube(), mu, s), 8.0) / exact - 1.0)
            .collect();
        converging(&format!("F1 cube μ = {mu:e}"), &errors);
    }
}

/// F2: a ferrite sphere (radius 1) in a uniform B along z: `3(μ − 1)/(μ + 2)` at μ = 4,
/// 1e12 and 0. Required as F1.
#[test]
fn f2_sphere_along_z() {
    let volume = 4.0 / 3.0 * std::f64::consts::PI;
    for mu in [4.0, 1e12, 0.0] {
        let exact = 3.0 * (mu - 1.0) / (mu + 2.0);
        let errors: Vec<f64> = SIZES
            .iter()
            .map(|&s| polarizability(&magnetized(sphere(), mu, s), volume) / exact - 1.0)
            .collect();
        converging(&format!("F2 sphere μ = {mu:e}"), &errors);
    }
}

/// F3: symmetry and vacuum. In the plane z = 0 the bound charge's field has no in-plane
/// component (exactly: the mirror's cancels it), so particles in the plane stay in it;
/// μ = 1 magnetizes nothing (exactly zero).
#[test]
fn f3_symmetry_and_vacuum() {
    let f = magnetized(cube(), 50.0, 0.5);
    for x in [
        DVec3::new(2.5, 0.3, 0.0),
        DVec3::new(-1.7, 2.2, 0.0),
        DVec3::new(0.4, -3.0, 0.0),
    ] {
        let b = f.field(x);
        assert!(b.x == 0.0 && b.y == 0.0 && b.z != 0.0, "F3 {x}: {b}");
    }
    let vacuum = magnetized(cube(), 1.0, 0.5);
    assert!(vacuum.sigma.iter().all(|&s| s == 0.0), "F3 μ = 1");
    println!("F3 in-plane field exactly along z; μ = 1 magnetizes nothing");
}
