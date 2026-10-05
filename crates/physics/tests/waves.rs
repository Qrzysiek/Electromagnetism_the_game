//! Validation tests W1–W5 for external fields (PHYSICS.md §2.3, §9). Run with
//! `cargo test -p physics --test waves -- --nocapture --test-threads=1`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use std::f64::consts::PI;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::external::{External, PlaneWave};
use physics::field::{Coulomb, FieldSolver, FixedCharge, LevelField};
use physics::magnetic::MagneticDipole;
use physics::trajectory::{Outcome, RunSettings, Scenario, run};

const TOL: f64 = 1e-12;

fn wave_scenario(
    wave: PlaneWave,
    c: f64,
    x0: DVec3,
    p0: DVec3,
    t_max: f64,
) -> Scenario<LevelField> {
    Scenario {
        field: LevelField {
            external: vec![External::Wave(wave)],
            ..LevelField::default()
        },
        obstacles: vec![],
        particle: UNIT_PARTICLE,
        c,
        x0,
        p0,
        detector: None,
        bounds: Some(cube(1e6)),
        t_max,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    }
}

// --- W1: Newtonian particle in a uniform oscillating field -------------------------------

/// `c = ∞`: `E = E₀ ê cos(ωt + φ)` is uniform, and with `a = qE₀/m`
/// `v(t) = v₀ + (a/ω)[sin(ωt+φ) − sin φ] ê`,
/// `x(t) = x₀ + v₀t − (a/ω)[(cos(ωt+φ) − cos φ)/ω + t sin φ] ê`.
#[test]
fn w1_newtonian_oscillating_field() {
    let (e0, omega) = (0.8, 2.3);
    let mut worst: f64 = 0.0;
    for (i, phase) in [0.0, 0.9, 2.0, 4.4].into_iter().enumerate() {
        let wave = PlaneWave::in_plane(e0, 0.6, omega, phase, f64::INFINITY);
        let x0 = DVec3::new(0.5, -1.0, 0.0);
        let v0 = DVec3::new(0.7, 0.2, 0.0);
        let t_end = 11.0;
        let scn = wave_scenario(wave, f64::INFINITY, x0, v0, t_end);
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        assert_eq!(tr.outcome, Outcome::Timeout);
        assert!(
            tr.energy_max_abs_error.is_nan(),
            "no energy diagnostic in a wave"
        );
        let a = e0 / UNIT_PARTICLE.mass * UNIT_PARTICLE.charge;
        let t = tr.end.t;
        let s = omega * t + phase;
        let x = x0 + v0 * t
            - wave.polarization * (a / omega * ((s.cos() - phase.cos()) / omega + t * phase.sin()));
        let v = v0 + wave.polarization * (a / omega * (s.sin() - phase.sin()));
        // Scales: the quiver amplitude a/ω² and velocity a/ω.
        let ex = (tr.end.x - x).length() / (a / (omega * omega));
        let ev = (tr.end.p - v).length() / (a / omega);
        println!("W1 case {i}: t = {t:.3}, |Δx|/(a/ω²) = {ex:.2e}, |Δv|/(a/ω) = {ev:.2e}");
        worst = worst.max(ex).max(ev);
    }
    assert!(worst < 1e-9, "W1 error {worst:.3e}");
}

// --- W2: relativistic particle in a plane wave ---------------------------------------------

/// Exact invariants of motion in a plane wave (Landau & Lifshitz, Classical Theory of
/// Fields §47–48): with `A` the vector potential in the gauge φ = 0,
/// the light-front momentum `γmc − p·k̂` and the canonical transverse momentum
/// `p_⊥ + qA` are conserved. Together with the mass shell they determine `p` as a
/// function of the wave phase, so checking them is equivalent to checking the analytic
/// solution for the momentum.
#[test]
fn w2_relativistic_plane_wave_invariants() {
    let c = 2.0;
    let omega = 1.5;
    let mut worst: f64 = 0.0;
    // Normalized amplitudes a₀ = qE₀/(mcω) from weakly to strongly relativistic.
    for (a0, angle, phase, p0) in [
        (0.1, 0.0, 0.3, DVec3::ZERO),
        (1.0, 0.8, 1.7, DVec3::new(0.5, -0.3, 0.0)),
        (3.0, -2.2, 4.0, DVec3::new(-2.0, 1.0, 0.0)),
        (2.0, 2.9, 0.0, DVec3::new(0.0, 0.0, 1.5)),
    ] {
        let m = UNIT_PARTICLE.mass;
        let q = UNIT_PARTICLE.charge;
        let e0 = a0 * m * c * omega / q;
        let wave = PlaneWave::in_plane(e0, angle, omega, phase, c);
        let x0 = DVec3::new(1.0, 2.0, 0.0);
        let t_max = 40.0 * 2.0 * PI / omega;
        let scn = wave_scenario(wave, c, x0, p0, t_max);
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        let kin = Kinematics::new(m, c);
        let invariants = |t: f64, x: DVec3, p: DVec3| {
            let light_front = kin.gamma(p) * m * c - p.dot(wave.direction);
            let canonical =
                p - wave.direction * p.dot(wave.direction) + wave.vector_potential(x, t) * q;
            (light_front, canonical)
        };
        let (lf0, can0) = invariants(0.0, x0, p0);
        let mut dev: f64 = 0.0;
        for s in &tr.samples {
            let (lf, can) = invariants(s.t, s.x, s.p);
            dev = dev.max((lf - lf0).abs()).max((can - can0).length());
        }
        let p_max = tr.samples.iter().map(|s| s.p.length()).fold(0.0, f64::max);
        // Relative to the momentum scale of the motion.
        let rel = dev / (m * c).max(p_max);
        println!(
            "W2 a0 = {a0}: {} steps, max |p|/(mc) = {:.2}, invariant drift / p = {rel:.2e}",
            tr.stats.n_accept,
            p_max / (m * c)
        );
        assert!(tr.samples.iter().all(|s| kin.velocity(s.p).length() < c));
        worst = worst.max(rel);
    }
    assert!(worst < 1e-9, "W2 invariant drift {worst:.3e}");
}

// --- W3: the 2D slice is exact with external fields ----------------------------------------

/// Charges and magnets in the plane plus an in-plane wave and uniform stray fields
/// (E in the plane, B along z): a particle starting in the plane never leaves it, bit for bit.
#[test]
fn w3_plane_symmetry_with_external_fields() {
    let c = 3.0;
    let field = LevelField {
        coulomb: Coulomb::new(&[FixedCharge {
            position: DVec3::new(4.0, 1.0, 0.0),
            charge: -2.0,
            radius: 0.2,
        }]),
        dipoles: vec![MagneticDipole {
            position: DVec3::new(-3.0, 2.0, 0.0),
            moment: DVec3::new(0.0, 0.0, 4.0),
            radius: 0.2,
        }],
        external: vec![
            External::Wave(PlaneWave::in_plane(0.4, 1.1, 2.0, 0.5, c)),
            External::Uniform {
                e: DVec3::new(0.05, -0.02, 0.0),
                b: DVec3::new(0.0, 0.0, 0.1),
            },
        ],
        ..LevelField::default()
    };
    assert!(!field.is_static());
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(-6.0, -1.0, 0.0),
        p0: DVec3::new(1.2, 0.4, 0.0),
        detector: None,
        bounds: Some(cube(40.0)),
        t_max: 60.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    println!("W3: {} steps, outcome {:?}", tr.stats.n_accept, tr.outcome);
    assert!(tr.stats.n_accept > 50);
    assert!(tr.samples.iter().all(|s| s.x.z == 0.0 && s.p.z == 0.0));
}

// --- W4: energy conservation with static stray fields --------------------------------------

/// A uniform stray E (potential −E·x) and B, with a charge: `(γ−1)mc² + qφ` is
/// conserved, the same check as T1.
#[test]
fn w4_energy_conservation_with_static_stray_fields() {
    for c in [f64::INFINITY, 4.0] {
        let field = LevelField {
            coulomb: Coulomb::new(&[FixedCharge {
                position: DVec3::new(3.0, 0.5, 0.0),
                charge: 1.5,
                radius: 0.2,
            }]),
            external: vec![
                External::Uniform {
                    e: DVec3::new(0.02, 0.07, 0.0),
                    b: DVec3::new(0.0, 0.0, 0.15),
                },
                // A static (ω = 0) wave term is a uniform electric field.
                External::Wave(PlaneWave::in_plane(0.03, 0.4, 0.0, 0.8, c)),
            ],
            ..LevelField::default()
        };
        assert!(field.is_static());
        let scn = Scenario {
            field,
            obstacles: vec![],
            particle: UNIT_PARTICLE,
            c,
            x0: DVec3::new(-5.0, 0.0, 0.0),
            p0: DVec3::new(1.0, 0.1, 0.0),
            detector: None,
            bounds: Some(cube(60.0)),
            t_max: 80.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        let rel = tr.energy_max_abs_error / tr.kinetic_initial;
        println!(
            "W4 c = {c}: {} steps, outcome {:?}, max |ΔW|/T0 = {rel:.2e}",
            tr.stats.n_accept, tr.outcome
        );
        assert!(tr.stats.n_accept > 20);
        assert!(rel < 1e-10, "energy error {rel:.3e}");
    }
}

/// W5, the Lawson–Woodward theorem (Brau Ex. 4.6) and its violation by radiation reaction:
/// a charge (q = 0.1, m = 1) at rest in the path of a plane-wave pulse (c = 2, ω = 1,
/// a₀ = qE₀/(mωc) = 1) given by its vector potential `A_y = (E₀/ω) g(φ) sin φ`,
/// `φ = ω(t − x/c) − φ₀`, `g = exp(−φ²/2σ²)`, σ = 6π (a test-local field: E and B from A,
/// so A vanishes exactly before and after the pulse). Without radiation reaction the
/// canonical momentum `p⊥ + qA` and the light-front momentum `ρ = γ − p_x/mc` are conserved:
/// the charge ends at rest (a plane wave in vacuum gives no net energy), displaced by
/// `Δx = (q²/2m²cω) ∫A² dφ = (q²E₀²/2m²cω³)(σ√π/2)(1 − e^{−σ²})` and `Δy = −(q/mω)∫A dφ = 0`.
/// With it, Landau–Lifshitz's exact plane-wave solution (Di Piazza 2008, test R7) gives
/// `1/ρ = 1 + (2q⁴/3m³c⁵ω) ∫E² dφ`, `∫E² dφ = E₀²√π (σ/2 + 1/4σ)` up to e^{−σ²}: the pulse
/// leaves the charge moving forward (radiation pressure). Required: without, |p| below
/// 1e-9 of mca₀, Δx within 1e-8, Δy within 1e-8 of Δx; with, 1/ρ within 1e-9.
#[test]
fn w5_lawson_woodward_and_radiation_pressure() {
    use physics::field::FieldSample;
    struct Pulse {
        e0: f64,
        w: f64,
        sigma: f64,
        phi0: f64,
        c: f64,
    }
    impl FieldSolver for Pulse {
        fn sample(&self, x: DVec3, t: f64) -> FieldSample {
            let phi = self.w * (t - x.x / self.c) - self.phi0;
            let g = (-0.5 * (phi / self.sigma).powi(2)).exp();
            let dg = -phi / (self.sigma * self.sigma) * g;
            // E_y = −∂A_y/∂t = −ω dA/dφ, B = x̂ × E / c.
            let (s, co) = phi.sin_cos();
            let ey = -self.e0 * (dg * s + g * co);
            FieldSample {
                e: DVec3::new(0.0, ey, 0.0),
                b: DVec3::new(0.0, 0.0, ey / self.c),
                phi: 0.0,
            }
        }
    }
    let (q, m, c, w, a0) = (0.1_f64, 1.0_f64, 2.0_f64, 1.0_f64, 1.0_f64);
    let e0 = a0 * m * w * c / q;
    let sigma = 6.0 * PI;
    let pulse = Pulse {
        e0,
        w,
        sigma,
        phi0: 8.0 * sigma,
        c,
    };
    let sp = PI.sqrt();
    let dx_ref = q * q * e0 * e0 / (2.0 * m * m * c * w.powi(3))
        * (sigma * sp / 2.0)
        * (1.0 - (-sigma * sigma).exp());
    let inv_rho_ref = 1.0
        + 2.0 * q.powi(4) / (3.0 * m.powi(3) * c.powi(5) * w)
            * e0
            * e0
            * sp
            * (sigma / 2.0 + 1.0 / (4.0 * sigma));
    let kin = Kinematics::new(m, c);
    let mut ok = true;
    for reaction in [false, true] {
        let scn = Scenario {
            field: Pulse { ..pulse },
            obstacles: vec![],
            particle: Particle {
                charge: q,
                mass: m,
                radius: 0.0,
                moment: 0.0,
            },
            c,
            x0: DVec3::ZERO,
            p0: DVec3::ZERO,
            detector: None,
            bounds: None,
            t_max: 450.0,
            radiation_reaction: reaction,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        assert_eq!(tr.outcome, Outcome::Timeout);
        let p = tr.end.p;
        let rho = kin.gamma(p) - p.x / (m * c);
        if reaction {
            let rel = 1.0 / rho / inv_rho_ref - 1.0;
            let v = kin.velocity(p);
            println!(
                "W5 with radiation reaction: 1/ρ = {:.12} (Di Piazza {inv_rho_ref:.12}): {rel:.1e}; \
                 the charge leaves at v = ({:.6e}, {:.1e}) c ({} steps)",
                1.0 / rho,
                v.x / c,
                v.y / c,
                tr.stats.n_accept
            );
            ok &= rel.abs() < 1e-9;
        } else {
            let (ex, ey) = (tr.end.x.x / dx_ref - 1.0, tr.end.x.y / dx_ref);
            let ep = p.length() / (m * c * a0);
            println!(
                "W5 without: |p| = {:.1e} mca₀ after the pulse (ρ = 1 + {:.1e}); Δx = {:.12} \
                 (Lawson–Woodward drift {dx_ref:.12}): {ex:.1e}; Δy = {:.1e} Δx ({} steps)",
                ep,
                rho - 1.0,
                tr.end.x.x,
                ey,
                tr.stats.n_accept
            );
            ok &= ep < 1e-9 && ex.abs() < 1e-8 && ey.abs() < 1e-8;
        }
    }
    assert!(ok, "W5");
}
