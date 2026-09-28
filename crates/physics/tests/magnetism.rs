//! Validation of magnetic fields and motion in them (PHYSICS.md §9, M1–M6). Run with
//! `cargo test -p physics --test magnetism -- --nocapture --test-threads=1`.

#![allow(clippy::cast_precision_loss)] // small loop counters

mod common;

use std::f64::consts::{PI, TAU};

use common::{Rng, cube};
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{Coulomb, FixedCharge, LevelField, UniformFields};
use physics::geometry::{Capsule, Shape, Sphere, Torus};
use physics::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};
use physics::trajectory::{
    Outcome, RunSettings, Scenario, StepView, Trajectory, run, run_observed,
};

const TOL: f64 = 1e-12;

fn rel(a: DVec3, b: DVec3) -> f64 {
    (a - b).length() / b.length()
}

/// Biot–Savart for a circular loop by the trapezoidal rule (spectrally accurate for this
/// smooth periodic integrand away from the wire).
fn loop_quadrature(l: &CircularLoop, x: DVec3, n: u32) -> DVec3 {
    let (e1, e2) = l.normal.any_orthonormal_pair();
    let mut b = DVec3::ZERO;
    for i in 0..n {
        let t = TAU * f64::from(i) / f64::from(n);
        let pos = l.center + (e1 * t.cos() + e2 * t.sin()) * l.radius;
        let dl = (-e1 * t.sin() + e2 * t.cos()) * (l.radius * TAU / f64::from(n));
        let r = x - pos;
        b += dl.cross(r) / r.length().powi(3);
    }
    b * l.kappa
}

/// Biot–Savart for a straight segment by composite 8-point Gauss–Legendre quadrature.
fn segment_quadrature(a: DVec3, b: DVec3, x: DVec3, kappa: f64) -> DVec3 {
    const NODES: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const WEIGHTS: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let pieces = 200;
    let dl = (b - a) / f64::from(pieces);
    let mut sum = DVec3::ZERO;
    for p in 0..pieces {
        let mid = a + dl * (f64::from(p) + 0.5);
        for (n, w) in NODES.iter().zip(WEIGHTS) {
            for s in [-1.0, 1.0] {
                let pos = mid + dl * (0.5 * s * n);
                let r = x - pos;
                sum += dl.cross(r) * (0.5 * w) / r.length().powi(3);
            }
        }
    }
    sum * kappa
}

// --- M1: field formulas -------------------------------------------------------------------

fn test_loop() -> CircularLoop {
    CircularLoop {
        center: DVec3::new(0.3, -0.2, 0.1),
        normal: DVec3::new(0.2, 0.3, 1.0).normalize(),
        radius: 2.0,
        kappa: 1.7,
        wire_radius: 0.05,
        rate: 0.0,
    }
}

#[test]
fn m1_circular_loop_matches_axis_formula_and_quadrature() {
    // On the axis: B = 2πκ a² / (a² + z²)^(3/2) along the normal.
    let l = test_loop();
    for z in [0.0, 0.5, 3.0, -7.0] {
        let x = l.center + l.normal * z;
        let exact =
            l.normal * (TAU * l.kappa * l.radius.powi(2) / (l.radius.powi(2) + z * z).powf(1.5));
        let e = rel(l.field(x), exact);
        println!("M1 loop on axis z = {z}: relative error {e:.1e}");
        assert!(e < 1e-14);
    }
    // Off the axis, including the plane of the loop, inside and outside.
    let (e1, e2) = l.normal.any_orthonormal_pair();
    let mut worst: f64 = 0.0;
    for (u, v, w) in [
        (0.5, 0.0, 0.0),
        (1.5, 0.3, 0.0),
        (2.6, 0.0, 0.0),
        (1.0, 1.0, 0.7),
        (3.0, -2.0, -1.5),
        (0.1, 0.2, 5.0),
    ] {
        let x = l.center + e1 * u + e2 * v + l.normal * w;
        let e = rel(l.field(x), loop_quadrature(&l, x, 20_000));
        worst = worst.max(e);
    }
    println!("M1 loop off axis vs Biot–Savart quadrature: worst relative error {worst:.1e}");
    assert!(worst < 1e-11);
}

#[test]
fn m1_loop_far_field_is_a_dipole() {
    // Moment of a loop: μ = κ π a² along the normal.
    let l = test_loop();
    let d = MagneticDipole {
        position: l.center,
        moment: l.normal * (l.kappa * PI * l.radius.powi(2)),
        radius: 0.0,
    };
    for r in [200.0, 2000.0] {
        let x = l.center + DVec3::new(0.6, -0.8, 0.3).normalize() * r;
        let e = rel(l.field(x), d.field(x));
        // Next multipole: relative correction O((a/r)²).
        let bound = 3.0 * (l.radius / r).powi(2);
        println!(
            "M1 loop far field at r = {r}: relative difference to dipole {e:.1e} (bound {bound:.1e})"
        );
        assert!(e < bound);
    }
}

#[test]
fn m1_polygon_coil_matches_quadrature_and_converges_to_circle() {
    let square = PolygonCoil {
        vertices: vec![
            DVec3::new(-1.0, -1.0, 0.0),
            DVec3::new(2.0, -1.0, 0.0),
            DVec3::new(2.0, 1.5, 0.0),
            DVec3::new(-1.0, 1.5, 0.0),
        ],
        kappa: 0.8,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(0.3, 0.2, 0.0),
        DVec3::new(4.0, 3.0, 0.0),
        DVec3::new(0.5, 0.5, 1.2),
    ] {
        let n = square.vertices.len();
        let q = (0..n).fold(DVec3::ZERO, |acc, i| {
            acc + segment_quadrature(
                square.vertices[i],
                square.vertices[(i + 1) % n],
                x,
                square.kappa,
            )
        });
        worst = worst.max(rel(square.field(x), q));
    }
    println!("M1 rectangular coil vs Gauss–Legendre quadrature: worst relative error {worst:.1e}");
    assert!(worst < 1e-12);

    // Regular N-gon inscribed in a circle: error O(1/N²).
    let circle = CircularLoop {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        radius: 3.0,
        kappa: 1.0,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let x = DVec3::new(1.1, 0.4, 0.5);
    let err = |n: u32| {
        let poly = PolygonCoil {
            vertices: (0..n)
                .map(|i| {
                    let t = TAU * f64::from(i) / f64::from(n);
                    DVec3::new(t.cos(), t.sin(), 0.0) * 3.0
                })
                .collect(),
            kappa: 1.0,
            wire_radius: 0.05,
            rate: 0.0,
        };
        rel(poly.field(x), circle.field(x))
    };
    let (e1, e2) = (err(500), err(1000));
    let order = (e1 / e2).log2();
    println!("M1 N-gon → circle: error(500) {e1:.2e}, error(1000) {e2:.2e}, order {order:.2}");
    assert!((1.9..2.1).contains(&order));
}

#[test]
fn m1_dipole_is_divergence_and_curl_free() {
    let d = MagneticDipole {
        position: DVec3::new(1.0, 2.0, -0.5),
        moment: DVec3::new(0.3, -0.4, 1.2),
        radius: 0.2,
    };
    let h = 1e-4;
    for x in [DVec3::new(2.0, 1.0, 0.5), DVec3::new(-1.0, 3.5, 0.0)] {
        let g = |axis: DVec3| (d.field(x + axis * h) - d.field(x - axis * h)) / (2.0 * h);
        let (gx, gy, gz) = (g(DVec3::X), g(DVec3::Y), g(DVec3::Z));
        let div = gx.x + gy.y + gz.z;
        let curl = DVec3::new(gy.z - gz.y, gz.x - gx.z, gx.y - gy.x);
        let scale = d.field(x).length() / (x - d.position).length();
        println!(
            "M1 dipole: |div B| / scale = {:.1e}, |curl B| / scale = {:.1e}",
            div.abs() / scale,
            curl.length() / scale
        );
        assert!(div.abs() < 1e-7 * scale && curl.length() < 1e-7 * scale);
    }
}

#[test]
fn m1_in_plane_sources_give_exactly_perpendicular_field_in_the_plane() {
    let field = LevelField {
        coulomb: Coulomb::default(),
        dipoles: vec![MagneticDipole {
            position: DVec3::new(3.0, 1.0, 0.0),
            moment: DVec3::new(0.0, 0.0, -2.0),
            radius: 0.3,
        }],
        loops: vec![CircularLoop {
            center: DVec3::new(-2.0, 0.5, 0.0),
            normal: DVec3::Z,
            radius: 4.0,
            kappa: 1.3,
            wire_radius: 0.05,
            rate: 0.0,
        }],
        polygons: vec![PolygonCoil {
            vertices: vec![
                DVec3::new(5.0, 5.0, 0.0),
                DVec3::new(9.0, 5.0, 0.0),
                DVec3::new(7.0, 8.0, 0.0),
            ],
            kappa: -0.7,
            wire_radius: 0.05,
            rate: 0.0,
        }],
        ..LevelField::default()
    };
    let mut rng = Rng::new(9);
    for _ in 0..200 {
        let x = DVec3::new(rng.range(-8.0, 12.0), rng.range(-6.0, 10.0), 0.0);
        let b = field.magnetic(x);
        assert_eq!(b.x.abs().to_bits(), 0, "{b:?} at {x:?}");
        assert_eq!(b.y.abs().to_bits(), 0, "{b:?} at {x:?}");
    }
}

// --- M2: relativistic cyclotron motion ----------------------------------------------------

#[test]
fn m2_relativistic_cyclotron_radius_and_period() {
    for (c, p) in [(f64::INFINITY, 1.0), (1.0, 0.5), (1.0, 3.0)] {
        let (q, m, b) = (1.0, 1.0, 2.0);
        let kin = Kinematics::new(m, c);
        let p0 = DVec3::new(0.0, p, 0.0);
        let gamma = kin.gamma(p0);
        let radius = p / (q * b);
        let period = TAU * gamma * m / (q * b);
        let turns = 10.0;
        let scn = Scenario {
            field: UniformFields {
                e: DVec3::ZERO,
                b: DVec3::new(0.0, 0.0, b),
            },
            obstacles: vec![],
            particle: Particle {
                charge: q,
                mass: m,
                radius: 0.0,
                moment: 0.0,
            },
            c,
            // q v×B at the start points along +x (ŷ × ẑ = x̂), so starting at x = −r puts
            // the centre of the orbit at the origin.
            x0: DVec3::new(-radius, 0.0, 0.0),
            p0,
            detector: None,
            bounds: None,
            t_max: turns * period,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let mut radius_err: f64 = 0.0;
        let mut p_err: f64 = 0.0;
        let tr = run_observed(
            &scn,
            &RunSettings::with_tolerance(TOL),
            |step: &StepView<'_, UniformFields>| {
                let (x, pp) = step.state(step.t_end());
                radius_err = radius_err.max((x.length() - radius).abs() / radius);
                p_err = p_err.max((pp.length() - p).abs() / p);
            },
        );
        let closure = (tr.end.x - scn.x0).length() / radius;
        println!(
            "M2 c = {c}, p = {p} (γ = {gamma:.3}): radius error {radius_err:.1e}, |p| drift {p_err:.1e}, \
             closure after {turns} periods {closure:.1e}"
        );
        // |p| is conserved exactly by the physics, but DOP853 does not preserve it, so it
        // drifts at the level of the integration tolerance (same criterion as T1).
        assert!(radius_err < 1e-10 && p_err < 1e-10 && closure < 1e-9);
    }
}

// --- M3: E×B drift ------------------------------------------------------------------------

#[test]
fn m3_exb_cycloid_newtonian() {
    // q = m = 1, B = B ẑ, E = E ŷ, from rest: x = (E/Bω)(ωt − sin ωt), y = (E/Bω)(1 − cos ωt).
    let (e, b) = (0.7, 2.0);
    let w = b;
    let t_end = 7.3;
    let scn = Scenario {
        field: UniformFields {
            e: DVec3::new(0.0, e, 0.0),
            b: DVec3::new(0.0, 0.0, b),
        },
        obstacles: vec![],
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c: f64::INFINITY,
        x0: DVec3::ZERO,
        p0: DVec3::ZERO,
        detector: None,
        bounds: None,
        t_max: t_end,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    let a = e / (b * w);
    let exact = DVec3::new(
        a * (w * t_end - (w * t_end).sin()),
        a * (1.0 - (w * t_end).cos()),
        0.0,
    );
    let err = (tr.end.x - exact).length() / exact.length();
    println!("M3 E×B cycloid (Newtonian): relative position error {err:.1e}");
    assert!(err < 1e-11);
}

#[test]
fn m3_exb_drift_relativistic() {
    // Crossed fields with E < cB, particle from rest. In the frame moving with
    // v_d = E/B the field is purely magnetic, B' = B/γ_d, and the particle circles with
    // speed v_d, so in the lab it returns to rest periodically. Between rests it moves
    // Δx = v_d T, with T = γ_d T' = 2π m γ_d³ / (qB).
    let (c, b): (f64, f64) = (1.0, 1.0);
    for ratio in [0.3, 0.6, 0.9] {
        let e = ratio * c * b;
        let v_d = e / b;
        let gamma_d = 1.0 / (1.0 - (v_d / c).powi(2)).sqrt();
        let period = TAU * gamma_d.powi(3) / b;
        let scn = Scenario {
            field: UniformFields {
                e: DVec3::new(0.0, e, 0.0),
                b: DVec3::new(0.0, 0.0, b),
            },
            obstacles: vec![],
            particle: Particle {
                charge: 1.0,
                mass: 1.0,
                radius: 0.0,
                moment: 0.0,
            },
            c,
            x0: DVec3::ZERO,
            p0: DVec3::ZERO,
            detector: None,
            bounds: None,
            t_max: 3.2 * period,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        // Rests: p_y changes sign from negative to positive (the cusps of the trajectory).
        let mut rests: Vec<(f64, DVec3)> = Vec::new();
        run_observed(
            &scn,
            &RunSettings::with_tolerance(TOL),
            |step: &StepView<'_, UniformFields>| {
                let (ta, tb) = (step.t_start(), step.t_end());
                let py = |t: f64| step.state(t).1.y;
                let n = 16;
                for i in 0..n {
                    let (mut lo, mut hi) = (
                        ta + (tb - ta) * f64::from(i) / f64::from(n),
                        ta + (tb - ta) * f64::from(i + 1) / f64::from(n),
                    );
                    if py(lo) < 0.0 && py(hi) >= 0.0 {
                        for _ in 0..200 {
                            let mid = 0.5 * (lo + hi);
                            if mid <= lo || mid >= hi {
                                break;
                            }
                            if py(mid) < 0.0 { lo = mid } else { hi = mid }
                        }
                        rests.push((hi, step.state(hi).0));
                    }
                }
            },
        );
        assert!(rests.len() >= 3, "found {} rests", rests.len());
        let (t3, x3) = rests[2];
        let period_err = (t3 / 3.0 - period).abs() / period;
        let drift_err = (x3.x / t3 - v_d).abs() / v_d;
        println!(
            "M3 relativistic E×B, E/cB = {ratio} (γ_d = {gamma_d:.3}): period error {period_err:.1e}, \
             drift speed error {drift_err:.1e}"
        );
        assert!(period_err < 1e-9 && drift_err < 1e-9);
    }
}

// --- M4: energy conservation --------------------------------------------------------------

fn mixed_field(seed: u64) -> (LevelField, Vec<Shape>) {
    let mut rng = Rng::new(seed);
    let charges: Vec<FixedCharge> = (0..8)
        .map(|_| FixedCharge {
            position: DVec3::new(rng.range(0.0, 20.0), rng.range(-8.0, 8.0), 0.0),
            charge: rng.range(-3.0, 3.0),
            radius: 0.3,
        })
        .collect();
    let dipoles: Vec<MagneticDipole> = (0..4)
        .map(|_| MagneticDipole {
            position: DVec3::new(rng.range(0.0, 20.0), rng.range(-8.0, 8.0), 0.0),
            moment: DVec3::new(0.0, 0.0, rng.range(-30.0, 30.0)),
            radius: 0.3,
        })
        .collect();
    let coil = CircularLoop {
        center: DVec3::new(10.0, 0.0, 0.0),
        normal: DVec3::Z,
        radius: 14.0,
        kappa: 2.0,
        wire_radius: 0.1,
        rate: 0.0,
    };
    let mut obstacles: Vec<Shape> = charges
        .iter()
        .map(|c| {
            Shape::Sphere(Sphere {
                center: c.position,
                radius: c.radius,
            })
        })
        .collect();
    obstacles.extend(dipoles.iter().map(|d| {
        Shape::Sphere(Sphere {
            center: d.position,
            radius: d.radius,
        })
    }));
    obstacles.push(Shape::Torus(Torus {
        center: coil.center,
        normal: coil.normal,
        major: coil.radius,
        minor: coil.wire_radius,
    }));
    let field = LevelField {
        coulomb: Coulomb::new(&charges),
        dipoles,
        loops: vec![coil],
        polygons: vec![],
        ..LevelField::default()
    };
    (field, obstacles)
}

#[test]
fn m4_energy_conservation_with_magnets() {
    for c in [f64::INFINITY, 5.0] {
        for seed in 1..=3 {
            let (field, obstacles) = mixed_field(seed);
            let scn = Scenario {
                field,
                obstacles,
                particle: Particle {
                    charge: 1.0,
                    mass: 1.0,
                    radius: 0.0,
                    moment: 0.0,
                },
                c,
                x0: DVec3::new(-2.0, 0.5, 0.0),
                p0: DVec3::new(1.5, 0.2, 0.0),
                detector: None,
                bounds: Some(cube(30.0)),
                t_max: 200.0,
                radiation_reaction: false,
                acceptance: None,
                gates: Vec::new(),
            };
            let tr = run(&scn, &RunSettings::with_tolerance(TOL));
            let rel = tr.energy_max_abs_error / tr.kinetic_initial;
            println!(
                "M4 c = {c}, seed {seed}: {:?} after {} steps, max |ΔW|/T0 = {rel:.1e}",
                tr.outcome, tr.stats.n_accept
            );
            assert!(tr.stats.n_accept >= 20);
            assert!(rel < 1e-10, "energy error {rel:.3e}");
        }
    }
}

// --- M6: wire collisions ------------------------------------------------------------------

#[test]
fn m6_particle_hits_coil_wire() {
    // No field; the particle flies straight at the wire of an in-plane ring and at a
    // straight wire; the collision time is known exactly.
    let ring = Torus {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        major: 5.0,
        minor: 0.1,
    };
    let wire = Capsule {
        a: DVec3::new(8.0, -3.0, 0.0),
        b: DVec3::new(8.0, 3.0, 0.0),
        radius: 0.1,
    };
    for (shape, x0, t_hit) in [
        (Shape::Torus(ring), DVec3::new(0.0, 0.0, 0.0), 4.9),
        (Shape::Capsule(wire), DVec3::new(6.0, 0.5, 0.0), 1.9),
    ] {
        let scn = Scenario {
            field: UniformFields {
                e: DVec3::ZERO,
                b: DVec3::ZERO,
            },
            obstacles: vec![shape],
            particle: Particle {
                charge: 1.0,
                mass: 1.0,
                radius: 0.0,
                moment: 0.0,
            },
            c: f64::INFINITY,
            x0,
            p0: DVec3::X,
            detector: None,
            bounds: Some(cube(20.0)),
            t_max: 100.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        println!("M6 {:?}: {:?} at t = {:.15}", shape, tr.outcome, tr.end.t);
        assert_eq!(tr.outcome, Outcome::Collided(0));
        assert!((tr.end.t - t_hit).abs() < 1e-12);
    }
}

// --- M7: Jackson Problem 12.9, equatorial drift around a dipole ------------------------

/// A charge gyrating in the equatorial plane of a dipole (gyroradius a at the guiding
/// centre's radius R, a ≪ R) drifts in longitude at `|dφ/dt| = (3/2) (a/R)² ω_B`
/// (Jackson Pr. 12.9b; the gradient drift of §12.4). With a/R = 0.01 the corrections are
/// O(a/R): the drift rate fitted over 3 rad of longitude (about 3000 gyrations) must
/// agree within 2 %.
#[test]
fn m7_van_allen_equatorial_drift() {
    let (r_gc, a) = (10.0_f64, 0.1_f64);
    // Moment along −z: in the equatorial plane B = +|M|/r³ ẑ, 100 at R.
    let dipole = MagneticDipole {
        position: DVec3::ZERO,
        moment: DVec3::new(0.0, 0.0, -1e5),
        radius: 0.5,
    };
    let omega = 1e5 / r_gc.powi(3);
    let v = a * omega;
    let field = LevelField {
        dipoles: vec![dipole],
        ..LevelField::default()
    };
    // q = 1 moving +y at x: the force v B x̂ puts the gyration centre at x + a.
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c: f64::INFINITY,
        x0: DVec3::new(r_gc - a, 0.0, 0.0),
        p0: DVec3::new(0.0, v, 0.0),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 200.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    // Unwrapped longitude against time; least-squares slope.
    let mut pts = Vec::with_capacity(tr.samples.len());
    let mut last = 0.0_f64;
    let mut turns = 0.0;
    for s in &tr.samples {
        let phi = s.x.y.atan2(s.x.x);
        if phi - last > PI {
            turns -= TAU;
        } else if last - phi > PI {
            turns += TAU;
        }
        last = phi;
        pts.push((s.t, phi + turns));
    }
    #[allow(clippy::cast_precision_loss)]
    let n = pts.len() as f64;
    let (mt, mp) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let slope = pts.iter().map(|&(t, p)| (t - mt) * (p - mp)).sum::<f64>()
        / pts.iter().map(|&(t, _)| (t - mt).powi(2)).sum::<f64>();
    let jackson = 1.5 * (a / r_gc).powi(2) * omega;
    println!(
        "M7: |dφ/dt| = {:.6e}, Jackson {jackson:.6e}, ratio {:.4} over {:.2} rad ({} steps)",
        slope.abs(),
        slope.abs() / jackson,
        (slope * 200.0).abs(),
        tr.stats.n_accept
    );
    assert!((slope.abs() / jackson - 1.0).abs() < 0.02);
}

// --- M8–M10: ramped coils, vector potential and induction (Jackson §5.5, §5.15) ------

fn ring(radius: f64, rate: f64) -> CircularLoop {
    CircularLoop {
        center: DVec3::new(1.0, -0.5, 0.0),
        normal: DVec3::Z,
        radius,
        kappa: 1.0,
        wire_radius: 0.05,
        rate,
    }
}

fn square(rate: f64) -> PolygonCoil {
    PolygonCoil {
        vertices: vec![
            DVec3::new(-3.0, -2.0, 0.0),
            DVec3::new(4.0, -2.0, 0.0),
            DVec3::new(4.0, 3.0, 0.0),
            DVec3::new(-3.0, 3.0, 0.0),
        ],
        kappa: 1.0,
        wire_radius: 0.05,
        rate,
    }
}

/// M8: the vector potential's curl is the field (κ = 1), for a circular loop (in and off
/// its plane, near and far) and a polygon coil, by central differences (h = 1e-5).
#[test]
fn m8_curl_of_the_vector_potential_is_the_field() {
    let h = 1e-5;
    let curl = |a: &dyn Fn(DVec3) -> DVec3, x: DVec3| {
        let d = |axis: DVec3| (a(x + axis * h) - a(x - axis * h)) / (2.0 * h);
        let (dx, dy, dz) = (d(DVec3::X), d(DVec3::Y), d(DVec3::Z));
        DVec3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x)
    };
    let l = ring(3.0, 0.0);
    let p = square(0.0);
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(1.2, -0.4, 0.0),
        DVec3::new(2.5, 1.0, 0.0),
        DVec3::new(6.0, 2.0, 0.0),
        DVec3::new(0.3, 0.8, 1.1),
        DVec3::new(5.0, -3.0, -0.7),
    ] {
        let e1 =
            (curl(&|y| l.unit_vector_potential(y), x) - l.field(x)).length() / l.field(x).length();
        let e2 =
            (curl(&|y| p.unit_vector_potential(y), x) - p.field(x)).length() / p.field(x).length();
        println!("M8 at {x:?}: loop {e1:.1e}, polygon {e2:.1e}");
        worst = worst.max(e1).max(e2);
    }
    assert!(worst < 1e-8);
}

/// Gauss–Legendre nodes and weights on [0, 1] (n points, Newton on P_n).
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
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
            (0.5 * (u + 1.0), 0.5 * w)
        })
        .collect()
}

/// Flux of `b` (its z-component) through the disc of radius `r` around `c` in the plane,
/// by Gauss–Legendre in the radius and the trapezoidal rule (periodic) in the angle.
fn flux_through(b: &dyn Fn(DVec3) -> DVec3, c: DVec3, r: f64) -> f64 {
    let gl = gauss_legendre(40);
    let m = 256;
    let mut total = 0.0;
    for &(s, w) in &gl {
        let rho = s * r;
        for j in 0..m {
            let phi = TAU * f64::from(j) / f64::from(m);
            total += b(c + DVec3::new(rho * phi.cos(), rho * phi.sin(), 0.0)).z * rho * w * r * TAU
                / f64::from(m);
        }
    }
    total
}

/// Circulation of `a` around the circle of radius `r` around `c` (trapezoidal rule).
fn circulation(a: &dyn Fn(DVec3) -> DVec3, c: DVec3, r: f64) -> f64 {
    let m = 2048;
    (0..m)
        .map(|j| {
            let phi = TAU * f64::from(j) / f64::from(m);
            let t = DVec3::new(-phi.sin(), phi.cos(), 0.0);
            a(c + DVec3::new(r * phi.cos(), r * phi.sin(), 0.0)).dot(t) * r * TAU / f64::from(m)
        })
        .sum()
}

/// M9, Faraday's law (Jackson §5.15): for a ramped coil `E = −rate · A_unit`, so
/// `∮E·dl = −rate ∮A_unit·dl`, which must equal `−dΦ/dt = −rate Φ_unit`, the flux of the
/// field computed independently by quadrature. Circles inside a loop (concentric and
/// off-centre) and inside a polygon coil, clear of the wires.
#[test]
fn m9_faradays_law_for_ramped_coils() {
    let rate = 0.37;
    let l = ring(3.0, rate);
    let p = square(rate);
    let mut worst: f64 = 0.0;
    for (name, c, r) in [
        ("loop, concentric", l.center, 2.0),
        (
            "loop, off-centre",
            l.center + DVec3::new(0.6, 0.3, 0.0),
            1.5,
        ),
    ] {
        let emf = circulation(&|x| l.induced_e(x), c, r);
        let dphi = rate * flux_through(&|x| l.field(x), c, r);
        let err = (emf + dphi).abs() / dphi.abs();
        println!(
            "M9 {name}: ∮E·dl = {emf:.12}, −dΦ/dt = {:.12}, rel. error {err:.1e}",
            -dphi
        );
        worst = worst.max(err);
    }
    let c = DVec3::new(0.5, 0.4, 0.0);
    let emf = circulation(&|x| p.induced_e(x), c, 1.8);
    let dphi = rate * flux_through(&|x| p.field(x), c, 1.8);
    let err = (emf + dphi).abs() / dphi.abs();
    println!(
        "M9 polygon: ∮E·dl = {emf:.12}, −dΦ/dt = {:.12}, rel. error {err:.1e}",
        -dphi
    );
    worst = worst.max(err);
    assert!(worst < 1e-8);
}

/// M10: in the field of a ramped circular loop (axially symmetric, time-dependent) the
/// canonical angular momentum about the axis, `x p_y − y p_x + q (x A_y − y A_x)` with
/// `A = κ(t) A_unit`, is conserved exactly along any orbit in the loop's plane.
#[test]
fn m10_canonical_angular_momentum_in_a_ramped_loop() {
    let l = CircularLoop {
        center: DVec3::ZERO,
        ..ring(6.0, 0.02)
    };
    let (q, m) = (1.0, 1.0);
    let field = LevelField {
        loops: vec![l],
        ..LevelField::default()
    };
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: q,
            mass: m,
            radius: 0.0,
            moment: 0.0,
        },
        c: f64::INFINITY,
        x0: DVec3::new(2.0, 0.5, 0.0),
        p0: DVec3::new(0.1, 0.3, 0.0),
        detector: None,
        bounds: Some(cube(5.0)),
        t_max: 60.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    let canonical = |t: f64, x: DVec3, p: DVec3| {
        let a = l.unit_vector_potential(x) * l.kappa_at(t);
        x.x * p.y - x.y * p.x + q * (x.x * a.y - x.y * a.x)
    };
    let l0 = canonical(0.0, scn.x0, scn.p0);
    let scale = scn.x0.length() * scn.p0.length()
        + (q * l.unit_vector_potential(scn.x0).length() * scn.x0.length());
    let worst = tr
        .samples
        .iter()
        .map(|s| (canonical(s.t, s.x, s.p) - l0).abs() / scale)
        .fold(0.0, f64::max);
    let dp = (tr.end.p.length() - scn.p0.length()).abs() / scn.p0.length();
    println!(
        "M10: {} samples to t = {:.1}; max |ΔL_can| / scale = {worst:.1e}; |p| changed by {dp:.2} (the induced field works)",
        tr.samples.len(),
        tr.end.t
    );
    assert!(worst < 1e-10 && dp > 0.01);
}

/// Guiding centre `x + p×ẑ/(qB)` of a charge gyrating in a loop's plane.
fn guiding_centre(l: &CircularLoop, q: f64, x: DVec3, p: DVec3, t: f64) -> DVec3 {
    x + p.cross(DVec3::Z) / (q * l.field_at(x, t).z)
}

/// Time average of the magnetic moment `p²/B(guiding centre)` over `[t0, t1]`.
fn mean_moment(tr: &Trajectory, l: &CircularLoop, q: f64, t0: f64, t1: f64) -> f64 {
    let (mut sum, mut weight) = (0.0, 0.0);
    for w in tr.samples.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        if a.t >= t0 && b.t <= t1 {
            let gc = guiding_centre(l, q, a.x, a.p, a.t);
            sum += a.p.length_squared() / l.field_at(gc, a.t).z * (b.t - a.t);
            weight += b.t - a.t;
        }
    }
    sum / weight
}

/// M11, adiabatic invariance (Jackson §12.5). In a slowly rising field the magnetic
/// moment `p⊥²/B` (the flux through the gyro-orbit) is conserved while `p²` grows with
/// `B`, and the guiding centre drifts inward (E×B in the induced field) keeping the flux
/// through its own circle about the axis, `ρ² B`. A gyrating charge (gyroradius 0.5,
/// ω_c ≈ 1) with its guiding centre 2 from the axis of a large loop (radius 20, where
/// B is uniform to 0.8 %), while the current rises 4× over 150, 300 and 600 (at most
/// 0.13 of the field's rise per gyroperiod). The moment, averaged over three
/// gyroperiods at the start and at the end, changes by at most 1.2e-3 and the guiding
/// centre ends at ρ = 1.002–1.004 (flux conservation: 2/√4 = 1, up to the field's
/// non-uniformity).
#[test]
fn m11_adiabatic_invariance_of_the_magnetic_moment() {
    let (q, m) = (1.0, 1.0);
    let radius = 20.0;
    // B at the centre is 2πκ/R: κ₀ gives ω_c ≈ 1.
    let kappa0 = radius / TAU;
    for t_ramp in [150.0, 300.0, 600.0] {
        let l = CircularLoop {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius,
            kappa: kappa0,
            wire_radius: 0.05,
            rate: 3.0 * kappa0 / t_ramp,
        };
        let field = LevelField {
            loops: vec![l],
            ..LevelField::default()
        };
        let b0 = l.field_at(DVec3::new(2.0, 0.0, 0.0), 0.0).z;
        let p0 = 0.5 * q * b0;
        let scn = Scenario {
            field,
            obstacles: vec![],
            particle: Particle {
                charge: q,
                mass: m,
                radius: 0.0,
                moment: 0.0,
            },
            c: f64::INFINITY,
            // Guiding centre at (2, 0): start 0.5 above it, moving so that it circles it.
            x0: DVec3::new(2.0, 0.5, 0.0),
            p0: DVec3::new(p0, 0.0, 0.0),
            detector: None,
            bounds: Some(cube(15.0)),
            t_max: t_ramp,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        let period = TAU / (q * b0 / m);
        let mu0 = mean_moment(&tr, &l, q, 0.0, 3.0 * period);
        let mu1 = mean_moment(&tr, &l, q, t_ramp - 3.0 * period, t_ramp);
        let rho = guiding_centre(&l, q, tr.end.x, tr.end.p, tr.end.t).length();
        println!(
            "M11 ramp over {t_ramp}: B ×4, μ changed by {:.1e}, guiding centre at ρ = {rho:.4} (flux: 1)",
            mu1 / mu0 - 1.0
        );
        assert!((mu1 / mu0 - 1.0).abs() < 2e-3 && (rho - 1.0).abs() < 0.01);
    }
}
