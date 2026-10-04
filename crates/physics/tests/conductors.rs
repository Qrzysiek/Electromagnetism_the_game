//! Validation tests K1–K4 for conducting spheres (PHYSICS.md §2.6). Run with
//! `cargo test -p physics --test conductors -- --nocapture --test-threads=1`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use std::f64::consts::PI;

use common::cube;
use physics::DVec3;
use physics::conductor::{Bias, Conductors, Resolution, SphereConductor, kelvin_images};
use physics::dynamics::Particle;
use physics::field::{Coulomb, FieldSolver, FixedCharge, LevelField};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{RunSettings, Scenario, run};

fn sphere(x: f64, y: f64, r: f64, bias: Bias) -> SphereConductor {
    SphereConductor {
        center: DVec3::new(x, y, 0.0),
        radius: r,
        bias,
    }
}

/// Potential of the fixed sources plus the induced charges.
fn potential(sources: &[(DVec3, f64)], c: &Conductors, x: DVec3) -> f64 {
    sources
        .iter()
        .map(|&(p, q)| q / (x - p).length())
        .sum::<f64>()
        + c.induced.sample(x, 0.0).phi
}

fn field(sources: &[(DVec3, f64)], c: &Conductors, x: DVec3) -> DVec3 {
    let mut e = c.induced.sample(x, 0.0).e;
    for &(p, q) in sources {
        let d = x - p;
        e += d * (q / d.length().powi(3));
    }
    e
}

/// Points on a sphere (Fibonacci lattice).
fn surface(s: &SphereConductor, n: usize) -> Vec<DVec3> {
    let golden = PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let (fi, fn_) = (i as f64, n as f64);
            let z = 1.0 - 2.0 * (fi + 0.5) / fn_;
            let r = (1.0 - z * z).sqrt();
            let phi = golden * fi;
            s.center + DVec3::new(r * phi.cos(), r * phi.sin(), z) * s.radius
        })
        .collect()
}

/// Outward flux of E through a sphere of radius `r` around `c`, by Gauss–Legendre in
/// cos θ (n nodes) and the trapezoidal rule in φ.
fn flux(e: &dyn Fn(DVec3) -> DVec3, c: DVec3, r: f64) -> f64 {
    // 32-point Gauss–Legendre via Newton on P_n.
    let n = 32;
    let mut total = 0.0;
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let mut u = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, u);
            for k in 2..=n {
                #[allow(clippy::cast_precision_loss)]
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * u * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            #[allow(clippy::cast_precision_loss)]
            let nf = n as f64;
            dp = nf * (u * p1 - p0) / (u * u - 1.0);
            let du = p1 / dp;
            u -= du;
            if du.abs() < 1e-16 {
                break;
            }
        }
        let w = 2.0 / ((1.0 - u * u) * dp * dp);
        let s = (1.0 - u * u).sqrt();
        let m = 64;
        for j in 0..m {
            let phi = 2.0 * PI * f64::from(j) / f64::from(m);
            let nrm = DVec3::new(s * phi.cos(), s * phi.sin(), u);
            total += e(c + nrm * r).dot(nrm) * r * r * w * 2.0 * PI / f64::from(m);
        }
    }
    total
}

// --- K1: a single sphere: boundary condition and image forces ---------------------------

/// Grounded sphere: potential 0 on the surface; the image force on a charge q at
/// distance d is −q² a d / (d² − a²)² (towards the sphere). Isolated neutral sphere:
/// −q² a³ (2d² − a²) / (d³ (d² − a²)²) (Jackson §2.3, Smythe §5.08).
#[test]
fn k1_single_sphere_boundary_condition_and_image_forces() {
    let a = 1.3;
    let src = [(DVec3::new(3.0, 0.7, 0.0), 2.0)];
    let grounded = Conductors::new(
        vec![sphere(0.0, 0.0, a, Bias::Grounded)],
        &src,
        Resolution::Verify,
    );
    let worst = surface(&grounded.spheres[0], 400)
        .into_iter()
        .map(|x| potential(&src, &grounded, x).abs())
        .fold(0.0, f64::max);
    let scale = 2.0 / (src[0].0.length() - a);
    println!(
        "K1 grounded: max |φ| on the surface / scale = {:.2e}",
        worst / scale
    );
    assert!(worst / scale < 1e-14);

    let q = 0.7;
    let mut worst_force: f64 = 0.0;
    for d in [1.5, 2.0, 5.0] {
        let x = DVec3::new(d, 0.0, 0.0);
        let g = Conductors::new(
            vec![sphere(0.0, 0.0, a, Bias::Grounded)],
            &[],
            Resolution::Verify,
        );
        let f = g.self_field(x, q).0 * q;
        let exact = -q * q * a * d / (d * d - a * a).powi(2);
        let neutral = Conductors::new(
            vec![sphere(0.0, 0.0, a, Bias::Charge(0.0))],
            &[],
            Resolution::Verify,
        );
        let fn_ = neutral.self_field(x, q).0 * q;
        let exact_n =
            -q * q * a.powi(3) * (2.0 * d * d - a * a) / (d.powi(3) * (d * d - a * a).powi(2));
        let e1 = (f.x - exact).abs() / exact.abs() + f.y.abs() + f.z.abs();
        let e2 = (fn_.x - exact_n).abs() / exact_n.abs() + fn_.y.abs() + fn_.z.abs();
        println!("K1 d = {d}: grounded force rel. error {e1:.2e}, neutral {e2:.2e}");
        worst_force = worst_force.max(e1).max(e2);
    }
    assert!(worst_force < 1e-13);
}

// --- K2: several spheres with every kind of bias ---------------------------------------

/// Three spheres (grounded, floating with charge 1.5, held at potential 0.8) near two
/// fixed charges. The electrostatic solution is unique given each surface's potential
/// (grounded, fixed) or equipotentiality plus net charge (floating); both are checked
/// with the computed field: surface potentials on 400 points each, and the floating
/// sphere's net charge as the Gauss flux just outside it.
#[test]
fn k2_several_spheres_boundary_conditions() {
    let spheres = vec![
        sphere(0.0, 0.0, 1.0, Bias::Grounded),
        sphere(3.2, 0.5, 0.8, Bias::Charge(1.5)),
        sphere(-1.0, 3.0, 1.2, Bias::Potential(0.8)),
    ];
    let src = [
        (DVec3::new(2.0, -2.0, 0.0), 3.0),
        (DVec3::new(-3.0, -1.0, 0.0), -2.0),
    ];
    let c = Conductors::new(spheres.clone(), &src, Resolution::Verify);
    println!(
        "K2: {} induced charges, boundary residual (independent points) {:.1e}",
        c.induced.charges().count(),
        c.boundary_residual(&src)
    );
    assert!(c.boundary_residual(&src) < 1e-10);
    let mut floating_potential = Vec::new();
    for (i, s) in spheres.iter().enumerate() {
        let v: Vec<f64> = surface(s, 400)
            .into_iter()
            .map(|x| potential(&src, &c, x))
            .collect();
        let mean = v.iter().sum::<f64>() / 400.0;
        let spread = v.iter().map(|p| (p - mean).abs()).fold(0.0, f64::max);
        println!("K2 sphere {i}: potential {mean:.12}, max deviation {spread:.2e}");
        // Requirement (PHYSICS.md §2.6): 1e-10 of the sources' potential on the surfaces (~3).
        let tol = 1e-10 * 3.0;
        assert!(spread < tol, "sphere {i} not equipotential: {spread:.3e}");
        match s.bias {
            Bias::Grounded => assert!(mean.abs() < tol),
            Bias::Potential(v0) => assert!((mean - v0).abs() < tol),
            Bias::Charge(_) => floating_potential.push(mean),
        }
    }
    let s = spheres[1];
    let q_flux = flux(&|x| field(&src, &c, x), s.center, s.radius * 1.05) / (4.0 * PI);
    println!("K2 floating sphere: net charge from Gauss flux {q_flux:.12} (set 1.5)");
    assert!((q_flux - 1.5).abs() < 1e-9);
}

// --- K3: the induced charges of the particle obey the same conditions ----------------

/// The charges induced by a moving particle (its image tree to depth D, plus floating
/// corrections) keep grounded and fixed-potential spheres at their potentials and
/// floating spheres at their net charge, up to the stated truncation bound ρ^D.
#[test]
fn k3_particle_images_obey_the_boundary_conditions() {
    let spheres = vec![
        sphere(0.0, 0.0, 1.0, Bias::Charge(0.0)),
        sphere(3.0, 0.0, 1.0, Bias::Potential(2.0)),
    ];
    let c = Conductors::new(spheres.clone(), &[], Resolution::Verify);
    let (x, q) = (DVec3::new(1.4, 1.1, 0.0), 0.9);
    let induced = c.induced_by(x, q);
    let pot = |p: DVec3| {
        q / (p - x).length()
            + induced
                .iter()
                .map(|&(y, cq)| cq / (p - y).length())
                .sum::<f64>()
    };
    // Scale: the particle's own potential on the nearest surface.
    let scale = q / 0.4;
    let bound = 10.0 * c.self_truncation * scale;
    println!("K3 truncation bound ρ^D = {:.2e}", c.self_truncation);
    for (i, s) in spheres.iter().enumerate() {
        let v: Vec<f64> = surface(s, 300).into_iter().map(pot).collect();
        let mean = v.iter().sum::<f64>() / 300.0;
        let spread = v.iter().map(|p| (p - mean).abs()).fold(0.0, f64::max);
        println!(
            "K3 sphere {i}: particle's contribution {:.2e}, deviation {:.2e} (relative to its potential)",
            mean / scale,
            spread / scale
        );
        assert!(spread < bound);
        if i == 1 {
            assert!(
                mean.abs() < bound,
                "fixed-potential sphere moved: {mean:.3e}"
            );
        }
    }
    let net: f64 = induced
        .iter()
        .filter(|&&(y, _)| (y - spheres[0].center).length() < 1.0)
        .map(|&(_, cq)| cq)
        .sum();
    println!("K3 floating sphere: net induced charge {net:.2e}");
    // With the reciprocity-based floating correction the net charge is kept up to the
    // truncation of the image tree.
    assert!(net.abs() < c.self_truncation * q);
}

// --- K2b: an independent reference: two grounded spheres by the full image series -------

/// For two spheres the image series does not branch, so it can be summed to machine
/// precision (200 generations). The induced field of the image + MFS solver must agree.
#[test]
fn k2b_two_spheres_match_the_full_image_series() {
    let spheres = vec![
        sphere(0.0, 0.0, 1.0, Bias::Grounded),
        sphere(2.6, 0.4, 0.9, Bias::Grounded),
    ];
    let src = [
        (DVec3::new(1.2, 2.0, 0.0), 1.0),
        (DVec3::new(-1.5, -1.0, 0.0), -0.7),
    ];
    let c = Conductors::new(spheres.clone(), &src, Resolution::Verify);
    let s: Vec<(DVec3, f64, Option<usize>)> = src.iter().map(|&(p, q)| (p, q, None)).collect();
    let (exact, _) = kelvin_images(&spheres, &s, 200);
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(1.3, -0.2, 0.0),
        DVec3::new(-2.0, 2.0, 0.5),
        DVec3::new(4.0, 1.0, -0.3),
        DVec3::new(0.2, 1.4, 0.0),
    ] {
        let a = c.induced.sample(x, 0.0);
        let phi_exact: f64 = exact.iter().map(|&(p, q)| q / (x - p).length()).sum();
        let e_exact: DVec3 = exact
            .iter()
            .map(|&(p, q)| (x - p) * (q / (x - p).length().powi(3)))
            .sum();
        worst = worst
            .max((a.phi - phi_exact).abs() / phi_exact.abs())
            .max((a.e - e_exact).length() / e_exact.length());
    }
    println!("K2b: largest relative difference to the full image series {worst:.2e}");
    assert!(worst < 1e-9);
}

// --- K4: energy conservation with image forces ---------------------------------------

/// A strongly charged particle flying past a grounded and a floating sphere and a fixed
/// charge: `W = (γ−1)mc² + qφ + ½ q φ_self` is conserved (the image force is
/// conservative with potential ½ q φ_self for linear conductors).
#[test]
fn k4_energy_conservation_with_image_forces() {
    for c in [f64::INFINITY, 5.0] {
        let spheres = vec![
            sphere(0.0, 0.0, 1.0, Bias::Grounded),
            sphere(4.0, 2.0, 1.5, Bias::Charge(0.5)),
        ];
        let fixed = FixedCharge {
            position: DVec3::new(2.0, -3.0, 0.0),
            charge: 1.0,
            radius: 0.2,
        };
        let conductors = Conductors::new(
            spheres.clone(),
            &[(fixed.position, fixed.charge)],
            Resolution::Verify,
        );
        let field = LevelField {
            coulomb: Coulomb::new(&[fixed]),
            conductors,
            ..LevelField::default()
        };
        let obstacles = spheres
            .iter()
            .map(|s| {
                Shape::Sphere(Sphere {
                    center: s.center,
                    radius: s.radius,
                })
            })
            .collect();
        let scn = Scenario {
            field,
            obstacles,
            particle: Particle {
                charge: 1.0,
                mass: 1.0,
                radius: 0.05,
                moment: 0.0,
            },
            c,
            x0: DVec3::new(-6.0, 1.2, 0.0),
            p0: DVec3::new(0.8, 0.0, 0.0),
            detector: None,
            bounds: Some(cube(20.0)),
            t_max: 60.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
        let rel = tr.energy_max_abs_error / tr.kinetic_initial;
        println!(
            "K4 c = {c}: {:?} after {} steps, max |ΔW|/T0 = {rel:.2e}",
            tr.outcome, tr.stats.n_accept
        );
        assert!(tr.stats.n_accept > 20);
        assert!(rel < 1e-10, "energy error {rel:.3e}");
    }
}

// --- K5: Jackson Problem 2.4, like charges attract near an isolated sphere --------------

/// A point charge q at distance d from the centre of an isolated conducting sphere
/// (radius R) carrying the charge Q of the same sign is repelled far away and attracted
/// close to the surface. The force vanishes where
/// `Q/d² = q R³ (2d² − R²) / (d³ (d² − R²)²)` (the net charge at the centre against the
/// image pair of a neutral sphere). Jackson's answers (Pr. 2.4): d/R − 1 = 0.6178 for
/// Q = q (the golden ratio), 0.4276 for Q = 2q, 0.8823 for Q = q/2. The engine's force
/// (the sphere's field plus the particle's own images) must change sign at the exact
/// root, found by bisection on both, and agree with the book's rounded values.
#[test]
fn k5_like_charges_attract_near_an_isolated_sphere() {
    let (r, q) = (2.0, 0.8);
    let bisect = |f: &dyn Fn(f64) -> f64| {
        let (mut lo, mut hi) = (r * 1.01, r * 3.0);
        assert!(f(lo) < 0.0 && f(hi) > 0.0, "no sign change");
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    for (ratio, book) in [(1.0, 0.6178), (2.0, 0.4276), (0.5, 0.8823)] {
        let big_q = ratio * q;
        let c = Conductors::new(
            vec![sphere(0.0, 0.0, r, Bias::Charge(big_q))],
            &[],
            Resolution::Verify,
        );
        // Radial force on the particle (positive: repelled).
        let engine = |d: f64| {
            let x = DVec3::new(d, 0.0, 0.0);
            q * (c.induced.sample(x, 0.0).e.x + c.self_field(x, q).0.x)
        };
        let exact = |d: f64| {
            q * big_q / (d * d)
                - q * q * r.powi(3) * (2.0 * d * d - r * r) / (d.powi(3) * (d * d - r * r).powi(2))
        };
        let (d_engine, d_exact) = (bisect(&engine), bisect(&exact));
        println!(
            "K5 Q/q = {ratio}: d/R − 1 = {:.6} (engine), {:.6} (exact), {book} (Jackson)",
            d_engine / r - 1.0,
            d_exact / r - 1.0
        );
        assert!((d_engine - d_exact).abs() / r < 1e-10);
        assert!((d_exact / r - 1.0 - book).abs() < 5e-4);
    }
}
