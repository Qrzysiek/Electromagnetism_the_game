//! Validation tests A1–A4 for oscillating dipoles (PHYSICS.md §2.4, §9). Run with
//! `cargo test -p physics --test antenna -- --nocapture --test-threads=1`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use std::f64::consts::PI;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::antenna::OscillatingDipole;
use physics::field::{Coulomb, FieldSolver, FixedCharge, LevelField};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{RunSettings, Scenario, run};

/// Gauss–Legendre nodes and weights on [−1, 1] (Newton iteration on P_n).
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            // P_n(x) and P_n'(x) by the three-term recurrence.
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                #[allow(clippy::cast_precision_loss)]
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            dp = nf * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

// --- A1: radiated power -------------------------------------------------------------------

/// The time-averaged Poynting flux `S = (c²/4π) E × B` (units k = 1, μ₀ = 4π/c²) through
/// a sphere of any radius around the dipole equals the Larmor power
/// `⟨P⟩ = p₀² ω⁴ / (3c³)`: the near-field terms carry only reactive power, whose time
/// average vanishes exactly. Checked from inside the near zone (R = 0.2 λ) to deep in the
/// radiation zone (R = 20 λ). The quadratures are exact for the trigonometric
/// polynomials involved (Gauss–Legendre in cos θ, uniform in φ and in time).
#[test]
fn a1_radiated_power_equals_larmor_at_every_radius() {
    let c = 3.0;
    let dip = OscillatingDipole {
        position: DVec3::new(0.3, -0.2, 0.1),
        amplitude: DVec3::new(0.6, -1.1, 0.4),
        omega: 2.0,
        phase: 0.7,
        c,
        radius: 0.05,
    };
    let lambda = 2.0 * PI * c / dip.omega;
    let larmor = dip.larmor_power();
    let gl = gauss_legendre(24);
    let (n_phi, n_t) = (48, 16);
    let mut worst: f64 = 0.0;
    for r_over_lambda in [0.2, 1.0, 20.0] {
        let r = r_over_lambda * lambda;
        let mut power = 0.0;
        for k in 0..n_t {
            let t = 2.0 * PI / dip.omega * f64::from(k) / f64::from(n_t);
            for &(u, w) in &gl {
                let s = (1.0 - u * u).sqrt();
                for j in 0..n_phi {
                    let phi = 2.0 * PI * f64::from(j) / f64::from(n_phi);
                    let n = DVec3::new(s * phi.cos(), s * phi.sin(), u);
                    // Time is shifted by r/c so every radius sees the same emitted cycle
                    // (not needed for the average; it only makes the check uniform).
                    let f = dip.fields(dip.position + n * r, t + r / c);
                    let flux = c * c / (4.0 * PI) * f.e.cross(f.b).dot(n);
                    power += flux * r * r * w * (2.0 * PI / f64::from(n_phi));
                }
            }
        }
        power /= f64::from(n_t);
        let rel = (power - larmor).abs() / larmor;
        println!(
            "A1 R = {r_over_lambda} λ: ⟨P⟩ = {power:.15e}, Larmor {larmor:.15e}, rel. error {rel:.2e}"
        );
        worst = worst.max(rel);
    }
    assert!(worst < 1e-10, "A1 error {worst:.3e}");
}

// --- A2: static limit ---------------------------------------------------------------------

/// `ω = 0`: the field and the potential of a dipole `p` equal those of charges `±Q` at
/// `±d/2` along `p̂` with `Qd = |p|`, up to the quadrupole-free `O((d/r)²)` correction
/// (the pair has no quadrupole moment, so the next term is the octupole).
#[test]
fn a2_static_limit_is_the_charge_pair_field() {
    let p = DVec3::new(0.3, 0.8, -0.2);
    let dip = OscillatingDipole {
        position: DVec3::ZERO,
        amplitude: p,
        omega: 0.0,
        phase: 0.0,
        c: 4.0,
        radius: 0.05,
    };
    for d in [1e-2, 1e-3] {
        let q = p.length() / d;
        let pair = Coulomb::new(&[
            FixedCharge {
                position: p.normalize() * (d / 2.0),
                charge: q,
                radius: 0.0,
            },
            FixedCharge {
                position: -p.normalize() * (d / 2.0),
                charge: -q,
                radius: 0.0,
            },
        ]);
        let (mut worst, mut worst_phi): (f64, f64) = (0.0, 0.0);
        for x in [
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(-1.0, 1.5, 0.5),
            DVec3::new(0.2, 0.3, -3.0),
        ] {
            let (a, b) = (dip.fields(x, 1.3), pair.sample(x, 0.0));
            worst = worst.max((a.e - b.e).length() / b.e.length());
            // Potential relative to its size p/r² (it vanishes across the dipole).
            let scale = p.length() / x.length_squared();
            worst_phi = worst_phi.max((a.phi - b.phi).abs() / scale);
        }
        println!(
            "A2 d = {d}: max relative difference {worst:.2e} (E), {worst_phi:.2e} (potential)"
        );
        assert!(
            worst_phi < (d / 1.8_f64).powi(2),
            "A2 potential {worst_phi:.3e}"
        );
        // Octupole correction ~ (d/r)² with r ≥ 1.8.
        assert!(worst < (d / 1.8_f64).powi(2), "A2 {worst:.3e}");
    }
}

// --- A3: the 2D slice ---------------------------------------------------------------------

/// Antennas with in-plane moments (and fixed charges): at z = 0, E lies exactly in the
/// plane and B is exactly along z, so a particle starting in the plane stays there, bit
/// for bit.
#[test]
fn a3_plane_symmetry_with_antennas() {
    let c = 3.0;
    let antenna = |x: f64, y: f64, px: f64, py: f64| OscillatingDipole {
        position: DVec3::new(x, y, 0.0),
        amplitude: DVec3::new(px, py, 0.0),
        omega: 1.3,
        phase: 0.2,
        c,
        radius: 0.2,
    };
    let antennas = vec![antenna(3.0, 2.0, 2.0, 1.0), antenna(-2.0, -3.0, -1.0, 3.0)];
    let field = LevelField {
        coulomb: Coulomb::new(&[FixedCharge {
            position: DVec3::new(6.0, -1.0, 0.0),
            charge: 1.0,
            radius: 0.2,
        }]),
        antennas: antennas.clone(),
        time_offset: 0.4,
        ..LevelField::default()
    };
    assert!(!field.is_static());
    let scn = Scenario {
        field,
        obstacles: antennas
            .iter()
            .map(|a| {
                Shape::Sphere(Sphere {
                    center: a.position,
                    radius: a.radius,
                })
            })
            .collect(),
        particle: UNIT_PARTICLE,
        c,
        x0: DVec3::new(-8.0, 1.0, 0.0),
        p0: DVec3::new(1.0, -0.1, 0.0),
        detector: None,
        bounds: Some(cube(30.0)),
        t_max: 60.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
    println!("A3: {} steps, outcome {:?}", tr.stats.n_accept, tr.outcome);
    assert!(tr.stats.n_accept > 30);
    assert!(tr.samples.iter().all(|s| s.x.z == 0.0 && s.p.z == 0.0));
}

// --- A4: a static antenna conserves energy ------------------------------------------------

/// `ω = 0` (any c): the antenna is an electrostatic dipole with potential `φ = n·p/r²`, so
/// `T + qφ` is conserved along a flight past it. The kinetic energy swings by about half
/// the launch energy on this path (as in the audit's case: p = (0, 2) at (12, 11), passed
/// 3 cells below with T₀ = 0.5); the conserved energy must hold to the integration
/// accuracy.
#[test]
fn a4_static_antenna_conserves_energy() {
    for c in [f64::INFINITY, 5.0] {
        let antenna = OscillatingDipole {
            position: DVec3::new(12.0, 11.0, 0.0),
            amplitude: DVec3::new(0.0, 2.0, 0.0),
            omega: 0.0,
            phase: 0.0,
            c,
            radius: 0.3,
        };
        let field = LevelField {
            antennas: vec![antenna],
            ..LevelField::default()
        };
        assert!(field.is_static());
        let kin = physics::dynamics::Kinematics::new(1.0, c);
        let p0 = kin.momentum_from_kinetic_energy(0.5, DVec3::X);
        let scn = Scenario {
            field,
            obstacles: vec![Shape::Sphere(Sphere {
                center: antenna.position,
                radius: antenna.radius,
            })],
            particle: UNIT_PARTICLE,
            c,
            x0: DVec3::new(2.0, 8.0, 0.0),
            p0,
            detector: None,
            bounds: Some(cube(40.0)),
            t_max: 30.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
        let swing = tr
            .samples
            .iter()
            .map(|s| (kin.kinetic_energy(s.p) - 0.5).abs())
            .fold(0.0, f64::max);
        println!(
            "A4 c = {c}: kinetic energy swings by {:.3} T₀, energy error {:.2e} T₀",
            swing / 0.5,
            tr.energy_max_abs_error / 0.5
        );
        assert!(swing > 0.2 * 0.5, "the path must pass the dipole");
        assert!(
            tr.energy_max_abs_error < 1e-9 * 0.5,
            "{:.3e}",
            tr.energy_max_abs_error
        );
    }
}
