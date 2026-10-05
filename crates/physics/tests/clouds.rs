//! Charge clouds: uniformly charged spheres particles can fly through (Thomson's atom,
//! PHYSICS.md §2.1). Run with `cargo test --release -p physics --test clouds -- --nocapture`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use std::f64::consts::PI;

use common::cube;
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{ChargeCloud, Coulomb, FieldSolver, FixedCharge};
use physics::trajectory::{Outcome, RunSettings, Scenario, run};

const TOL: f64 = 1e-12;

/// C1: the cloud's field is minus the gradient of its potential, obeys Gauss's law
/// inside (div E = 4πρ with ρ = 3Q/(4πR³), i.e. div E = 3Q/R³), and is continuous at the
/// surface; outside it equals a point charge's. A point charge beside it is unchanged.
#[test]
fn c1_cloud_field_potential_gauss_and_continuity() {
    let cloud = ChargeCloud {
        position: DVec3::new(1.0, -2.0, 0.0),
        charge: 2.5,
        radius: 3.0,
    };
    let point = FixedCharge {
        position: DVec3::new(9.0, 4.0, 0.0),
        charge: -1.2,
        radius: 0.3,
    };
    let f = Coulomb::with_clouds(&[point], &[cloud]);
    let h = 1e-5;
    let mut worst_grad: f64 = 0.0;
    let mut worst_div: f64 = 0.0;
    for x in [
        DVec3::new(1.5, -1.0, 0.3),
        DVec3::new(-0.7, -3.1, -1.0),
        DVec3::new(6.0, 1.0, 0.0),
    ] {
        let s = f.sample(x, 0.0);
        let phi = |y: DVec3| f.sample(y, 0.0).phi;
        let e_of = |y: DVec3| f.sample(y, 0.0).e;
        let grad = DVec3::new(
            phi(x + DVec3::X * h) - phi(x - DVec3::X * h),
            phi(x + DVec3::Y * h) - phi(x - DVec3::Y * h),
            phi(x + DVec3::Z * h) - phi(x - DVec3::Z * h),
        ) / (2.0 * h);
        worst_grad = worst_grad.max((s.e + grad).length() / s.e.length());
        let div = (e_of(x + DVec3::X * h).x - e_of(x - DVec3::X * h).x + e_of(x + DVec3::Y * h).y
            - e_of(x - DVec3::Y * h).y
            + e_of(x + DVec3::Z * h).z
            - e_of(x - DVec3::Z * h).z)
            / (2.0 * h);
        let inside = (x - cloud.position).length() < cloud.radius;
        let expected = if inside {
            3.0 * cloud.charge / cloud.radius.powi(3)
        } else {
            0.0
        };
        worst_div =
            worst_div.max((div - expected).abs() / (3.0 * cloud.charge / cloud.radius.powi(3)));
    }
    // Continuity across the surface (field and potential), and the point-charge value
    // outside.
    let dir = DVec3::new(0.6, 0.8, 0.0);
    let at = |r: f64| f.sample(cloud.position + dir * r, 0.0);
    let (a, b) = (
        at(cloud.radius * (1.0 - 1e-12)),
        at(cloud.radius * (1.0 + 1e-12)),
    );
    let jump = (a.e - b.e).length() / b.e.length() + (a.phi - b.phi).abs() / b.phi.abs();
    let only_cloud = Coulomb::with_clouds(&[], &[cloud]);
    let far = cloud.position + dir * 7.0;
    let pc = Coulomb::new(&[FixedCharge {
        position: cloud.position,
        charge: cloud.charge,
        radius: 0.0,
    }]);
    let outside = (only_cloud.sample(far, 0.0).e - pc.sample(far, 0.0).e).length();
    println!(
        "C1: |E + ∇φ|/|E| ≤ {worst_grad:.1e}; Gauss |div E − 4πρ| / (3Q/R³) ≤ {worst_div:.1e}; \
         surface jump {jump:.1e}; outside vs point charge {outside:.1e} (4π = {:.4})",
        4.0 * PI
    );
    assert!(worst_grad < 1e-8 && worst_div < 1e-5 && jump < 1e-10 && outside == 0.0);
}

/// C2, Jackson Problem 16.1 (and §16.7): a charge bound harmonically (here inside a
/// charge cloud, ω₀² = |qQ|/(mR³)) radiates, and its oscillation energy decays as
/// `e^{−Γt}` with `Γ = ω₀² τ`, `τ = 2q²/(3mc³)`. Landau–Lifshitz agrees with the
/// Abraham–Lorentz result to O(ω₀τ) (here 0.01): the fitted decay rate over five decay
/// times must agree within 2 %. And exactly: the end state after the five decay times
/// (~600 oscillations) against an independent integration of the full relativistic
/// Landau–Lifshitz equation written from Landau & Lifshitz §76
/// (`scripts/wolfram/c2_radiating_oscillator.wls`, Wolfram Engine 14.2, NDSolve at 30
/// digits; the same to 17 digits with a tighter goal): position and momentum within
/// 1e-6 of the end amplitudes (set before measuring).
#[test]
fn c2_radiating_oscillator_decays_at_gamma() {
    let (q, m, big_q, r, c) = (-1.0_f64, 1.0_f64, 1.0_f64, 4.0_f64, 2.0_f64);
    let cloud = ChargeCloud {
        position: DVec3::ZERO,
        charge: big_q,
        radius: r,
    };
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let gamma = omega0 * omega0 * tau;
    let field = Coulomb::with_clouds(&[], &[cloud]);
    // An in-plane ellipse of amplitude ~1 (well inside the cloud).
    let v0 = 0.6 * omega0;
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: q,
            mass: m,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(1.0, 0.0, 0.0),
        p0: Kinematics::new(m, c).momentum_from_kinetic_energy(0.5 * m * v0 * v0, DVec3::Y),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 5.0 / gamma,
        radiation_reaction: true,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    // End state (t = 5/Γ = 3840) from the Wolfram script: x, y, px, py.
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const END: [f64; 4] = [
        -0.059818079583771145,
        0.032384524846739647,
        -0.0069765044047088572,
        -0.0046503164591507392,
    ];
    let x_ref = DVec3::new(END[0], END[1], 0.0);
    let p_ref = DVec3::new(END[2], END[3], 0.0);
    let ex = (tr.end.x - x_ref).length() / x_ref.length();
    let ep = (tr.end.p - p_ref).length() / p_ref.length();
    println!(
        "C2 end state at t = {:.1}: position {ex:.1e}, momentum {ep:.1e} from the independent LL integration",
        tr.end.t
    );
    assert!((tr.end.t - 3840.0).abs() < 1e-9 && ex < 1e-6 && ep < 1e-6);
    let kin = Kinematics::new(m, c);
    let phi0 = scn.field.sample(DVec3::ZERO, 0.0).phi;
    // ln(oscillation energy) against t: least-squares slope.
    let pts: Vec<(f64, f64)> = tr
        .samples
        .iter()
        .map(|s| {
            let w = kin.kinetic_energy(s.p) + q * (scn.field.sample(s.x, 0.0).phi - phi0);
            (s.t, w.ln())
        })
        .collect();
    #[allow(clippy::cast_precision_loss)]
    let n = pts.len() as f64;
    let (mt, my) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let slope = pts.iter().map(|&(t, y)| (t - mt) * (y - my)).sum::<f64>()
        / pts.iter().map(|&(t, _)| (t - mt).powi(2)).sum::<f64>();
    println!(
        "C2: ω₀ = {omega0:.4}, ω₀τ = {:.4}: fitted decay rate {:.6e}, Jackson Γ = ω₀²τ = {gamma:.6e}, \
         ratio {:.4} ({} samples, max |F_RR|/|F_L| = {:.1e})",
        omega0 * tau,
        -slope,
        -slope / gamma,
        pts.len(),
        tr.reaction_ratio_max
    );
    assert!((-slope / gamma - 1.0).abs() < 0.02);
}

/// The bound charge of tests C3 and C4: an electron (q = −1, m = 1) at the centre of a cloud
/// (Q = 8, R = 2: ω₀ = 1), c = 8 (τ = 2q²/(3mc³): ω₀τ = 1.3e-3), radiation reaction on.
const BOUND: (f64, f64, f64, f64, f64) = (-1.0, 1.0, 8.0, 2.0, 8.0);

/// One flight of the bound charge (`BOUND`) driven by a plane wave of frequency `w` along x,
/// polarized along y, for 20 periods, started on the analytic steady state, at the amplitude
/// 1e-4 · min(1, c/ω) (v/c ≤ 1e-4). Returns the measured extinction cross section (the
/// wave's work per time, in the steady state minus the reaction's, over the intensity
/// `I = cE₀²/8π`), Jackson's (16.78), the flight's radiated (Liénard) energy per time over
/// I, and the number of steps.
fn bound_charge_in_wave(w: f64) -> (f64, f64, f64, u64) {
    use physics::external::{External, PlaneWave};
    use physics::field::LevelField;
    let (q, m, big_q, r, c) = BOUND;
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let gamma = omega0 * omega0 * tau;
    let sigma_t = 8.0 * PI / 3.0 * (q * q / (m * c * c)).powi(2);
    // Jackson (16.74) with Γ' = 0: X = χ E₀ for E = Re[E₀ e^{−iωt}] at the origin,
    // χ = (q/m)(1 − iωτ)/(ω₀² − ω² − iωΓ), Γ = ω₀²τ; the field amplitude that gives |X|.
    let (nr, ni) = (1.0, -w * tau);
    let (dr, di) = (omega0 * omega0 - w * w, -w * gamma);
    let d2 = dr * dr + di * di;
    let (chi_r, chi_i) = (
        (nr * dr + ni * di) / d2 * q / m,
        (ni * dr - nr * di) / d2 * q / m,
    );
    let e0 = 1e-4 * (c / w).min(1.0) / (chi_r * chi_r + chi_i * chi_i).sqrt();
    // x(0) = Re[X], v(0) = Re[−iωX] = ω Im[X], along the polarization (y).
    let (xr, xi) = (chi_r * e0, chi_i * e0);
    let v0 = w * xi;
    let field = LevelField {
        coulomb: Coulomb::with_clouds(
            &[],
            &[ChargeCloud {
                position: DVec3::ZERO,
                charge: big_q,
                radius: r,
            }],
        ),
        external: vec![External::Wave(PlaneWave::in_plane(e0, 0.0, w, 0.0, c))],
        ..LevelField::default()
    };
    let kin = Kinematics::new(m, c);
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: q,
            mass: m,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(0.0, xr, 0.0),
        p0: DVec3::new(
            0.0,
            kin.gamma_of_velocity(DVec3::new(0.0, v0, 0.0)) * m * v0,
            0.0,
        ),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 20.0 * 2.0 * PI / w,
        radiation_reaction: true,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    assert_eq!(tr.outcome, Outcome::Timeout, "flight at ω = {w}");
    let intensity = c * e0 * e0 / (8.0 * PI);
    (
        -tr.radiation_work / (scn.t_max * intensity),
        sigma_t * w.powi(4) / d2,
        tr.radiated_energy / (scn.t_max * intensity),
        tr.stats.n_accept,
    )
}

/// C3, Jackson §16.8: scattering of a plane wave by a bound charge (`bound_charge_in_wave`),
/// from Rayleigh's `σ_T (ω/ω₀)⁴` below the resonance through its peak `6πc²/ω₀²` (Jackson's
/// 6πƛ₀²) to Thomson's `σ_T` above it. Reference: Jackson (16.78) with no other damping
/// (Γ' = 0), `σ_t = σ_T ω⁴/((ω₀² − ω²)² + ω²Γ²)`, `Γ = ω₀²τ`, the exact steady state of his
/// (16.73), which is the Landau–Lifshitz reduction (his (16.10)). Required: within 1e-4 (the
/// neglected O(v/c)² is below 1e-8). Printed: the steady state of Abraham–Lorentz (`x⃛ →
/// iω³x`, the width τω² for Γ), which differs by up to ω₀τ in the wings of the line and by
/// (ωτ)² far from it, the order either reduction neglects; and the flight's radiated
/// (Liénard) energy, which is the extinction: nothing is absorbed (Γ' = 0).
#[test]
fn c3_bound_charge_scattering_cross_section() {
    let (q, m, big_q, r, c) = BOUND;
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let sigma_t = 8.0 * PI / 3.0 * (q * q / (m * c * c)).powi(2);
    let mut worst: f64 = 0.0;
    for w in [0.1, 0.3, 0.7, 0.999, 1.0, 1.001, 1.5, 3.0, 10.0] {
        let (measured, jackson, radiated, steps) = bound_charge_in_wave(w);
        let abraham_lorentz =
            sigma_t * w.powi(4) / ((omega0 * omega0 - w * w).powi(2) + (tau * w.powi(3)).powi(2));
        let rel = (measured / jackson - 1.0).abs();
        println!(
            "C3 ω/ω₀ = {w}: σ/σ_T = {:.6e} (Jackson (16.78) {:.6e}, off by {rel:.1e}; \
             Abraham–Lorentz {:.2e} relative; radiated (Liénard) {:.2e} relative; {steps} steps)",
            measured / sigma_t,
            jackson / sigma_t,
            abraham_lorentz / jackson - 1.0,
            radiated / measured - 1.0,
        );
        worst = worst.max(rel);
    }
    println!(
        "C3: σ_T = {sigma_t:.4e}; at the resonance 6πc²/ω₀² = {:.4e}; worst {worst:.1e}",
        6.0 * PI * c * c / (omega0 * omega0)
    );
    assert!(worst < 1e-4, "C3 {worst:.3e}");
}

/// C4, Jackson Problem 16.13: the dipole sum rule `∫₀^∞ σ_t dω = 2π²q²/(mc)`. Its premise, a
/// polarizability tending to the free charge's `−q²/(mω²)`, fails for the radiating
/// oscillator of §16.8 by the factor (1 − iωτ) of (16.74): its σ_t tends to Thomson's σ_T
/// (Fig. 16.2), and the integral diverges. What holds is the rule with the free charge's
/// Thomson scattering subtracted, exactly for this oscillator: `∫₀^∞ (σ_t − σ_T) dω =
/// (2π²q²/mc)(1 − ω₀²τ²)` (`scripts/wolfram/c4_dipole_sum_rule.wls`, Wolfram Engine 14.2:
/// symbolically for any Γ < 2ω₀; Abraham–Lorentz's steady state, whose plateau ends at
/// ω ~ 1/τ, integrates to 2.0000 times the rule, its runaway pole in the upper half plane
/// breaking the Kramers–Kronig premise). The measured cross section of C3 at the nodes of
/// Gauss–Legendre rules on three pieces, each smooth: from 0.1 ω₀ to the line in ω, on the
/// line (ω₀ ± 50Γ) in θ with `ω = ω₀ + (Γ/2) tan θ` (which makes it flat), above it to
/// 200 ω₀ in 1/ω. Outside, Jackson's σ_t stands in: below 0.1 ω₀ it is under 1e-4 σ_T (its
/// integral 2e-9 of the sum; the subtracted σ_T is exact), and a flight there would take
/// 20 periods of the drive at steps set by ω₀; beyond 200 ω₀ σ_t − σ_T is 8e-6 of the sum.
/// Required: within 1e-4 of the exact value (the measured σ_t is within 1e-6, C3). Printed:
/// the rules' own error on Jackson's σ_t, and how far the rule without the factor 1 − ω₀²τ²
/// is.
#[test]
fn c4_dipole_sum_rule() {
    let (q, m, big_q, r, c) = BOUND;
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let gamma = omega0 * omega0 * tau;
    let sigma_t = 8.0 * PI / 3.0 * (q * q / (m * c * c)).powi(2);
    let jackson =
        |w: f64| sigma_t * w.powi(4) / ((omega0 * omega0 - w * w).powi(2) + (w * gamma).powi(2));
    // From the Wolfram script (q = −1, m = 1, c = 8, ω₀ = 1).
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const EXACT: f64 = 2.4673969169886816951;
    let sum_rule = 2.0 * PI * PI * q * q / (m * c);
    let (half, w_min, w_max) = (50.0 * gamma, 0.1 * omega0, 200.0 * omega0);
    // Each piece: (Gauss–Legendre order, ends of the variable, ω and dω/d(variable)).
    type Map = fn(f64, f64, f64) -> (f64, f64);
    let pieces: [(usize, f64, f64, Map); 3] = [
        (32, w_min, omega0 - half, |s, _, _| (s, 1.0)),
        (
            96,
            -(2.0 * half / gamma).atan(),
            (2.0 * half / gamma).atan(),
            |s, w0, g| {
                let t = libm::tan(s);
                (w0 + 0.5 * g * t, 0.5 * g * (1.0 + t * t))
            },
        ),
        (32, 1.0 / w_max, 1.0 / (omega0 + half), |s, _, _| {
            (1.0 / s, 1.0 / (s * s))
        }),
    ];
    let (mut measured, mut analytic, mut flights, mut steps) = (0.0, 0.0, 0, 0);
    let mut line = 0.0;
    for (i, (n, a, b, map)) in pieces.into_iter().enumerate() {
        let (mut piece_m, mut piece_a) = (0.0, 0.0);
        for (x, wt) in common::gauss_legendre(n) {
            let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let (w, dw) = map(s, omega0, gamma);
            let (sigma, _, _, k) = bound_charge_in_wave(w);
            let weight = 0.5 * (b - a) * wt * dw;
            piece_m += weight * (sigma - sigma_t);
            piece_a += weight * (jackson(w) - sigma_t);
            flights += 1;
            steps += k;
        }
        println!(
            "C4 piece [{a:.4}, {b:.4}] ({n} flights): measured {piece_m:.10e}, Jackson {piece_a:.10e}"
        );
        measured += piece_m;
        analytic += piece_a;
        if i == 1 {
            line = piece_m;
        }
    }
    // Outside 0.1–200 ω₀: Jackson's σ_t, below in ω, beyond in 1/ω.
    let (mut below, mut beyond, mut below_sigma) = (0.0, 0.0, 0.0);
    for (x, wt) in common::gauss_legendre(32) {
        let w = 0.5 * w_min * (1.0 + x);
        below += 0.5 * w_min * wt * (jackson(w) - sigma_t);
        below_sigma += 0.5 * w_min * wt * jackson(w);
        let u = 0.5 / w_max * (1.0 + x);
        beyond += 0.5 / w_max * wt * (jackson(1.0 / u) - sigma_t) / (u * u);
    }
    measured += below + beyond;
    analytic += below + beyond;
    let rel = measured / EXACT - 1.0;
    println!(
        "C4: ∫(σ − σ_T) dω measured {measured:.12e} (exact {EXACT:.12e}, off by {rel:.2e}; the \
         rules on Jackson's σ_t {:.1e}; Jackson's σ_t below 0.1 ω₀ {:.1e} of it, σ_t − σ_T \
         beyond 200 ω₀ {:.1e}); over 2π²q²/(mc) {:.8} (exact 1 − ω₀²τ² = {:.8}: the rule \
         without the factor is off by {:.2e}); the line alone (ω₀ ± 50Γ) {:.6} of the rule; \
         {flights} flights, {steps} steps",
        analytic / EXACT - 1.0,
        below_sigma / EXACT,
        beyond / EXACT,
        measured / sum_rule,
        1.0 - (omega0 * tau).powi(2),
        measured / sum_rule - 1.0,
        line / sum_rule,
    );
    assert!(rel.abs() < 1e-4, "C4 {rel:.3e}");
}

/// The passing charge of tests C5 and C7: an electron (q = −1, m = 1) at rest at the centre
/// of a neutral cloud (Q = 1, R = 4: ω₀ = 1/8), and a charge z = 16 passing at v = 0.8c
/// (c = 100) at the impact parameter `b = ξγv/ω₀`, its field the engine's Liénard–Wiechert
/// field of uniform motion (an infinitely heavy projectile), from 300 b before the closest
/// approach to 300 b after it. Returns b, the energy left in the oscillator and the steps.
fn passing_charge_transfer(xi: f64) -> (f64, f64, u64) {
    use physics::field::FieldSample;
    use physics::lienard::{Worldline, fields};
    struct Line {
        x0: DVec3,
        v: DVec3,
    }
    impl Worldline for Line {
        fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
            (self.x0 + self.v * t, self.v, DVec3::ZERO)
        }
    }
    struct Passing {
        atom: Coulomb,
        z: f64,
        path: Line,
        c: f64,
    }
    impl FieldSolver for Passing {
        fn sample(&self, x: DVec3, t: f64) -> FieldSample {
            let a = self.atom.sample(x, t);
            let f = fields(&self.path, self.z, self.c, x, t);
            FieldSample {
                e: a.e + f.e(),
                b: a.b + f.b,
                phi: a.phi,
            }
        }
    }
    let (q, m, big_q, r, z, c, v) = PASSING;
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let gamma = 1.0 / (1.0 - (v / c).powi(2)).sqrt();
    let atom = || {
        Coulomb::with_clouds(
            &[],
            &[ChargeCloud {
                position: DVec3::ZERO,
                charge: big_q,
                radius: r,
            }],
        )
    };
    let phi0 = atom().sample(DVec3::ZERO, 0.0).phi;
    let kin = Kinematics::new(m, c);
    let b = xi * gamma * v / omega0;
    let span = 300.0 * b;
    let scn = Scenario {
        field: Passing {
            atom: atom(),
            z,
            path: Line {
                x0: DVec3::new(-span, b, 0.0),
                v: DVec3::new(v, 0.0, 0.0),
            },
            c,
        },
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
        t_max: 2.0 * span / v,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    assert_eq!(tr.outcome, Outcome::Timeout, "the pass at ξ = {xi}");
    let energy =
        kin.kinetic_energy(tr.end.p) + q * (scn.field.atom.sample(tr.end.x, 0.0).phi - phi0);
    (b, energy, tr.stats.n_accept)
}

/// q, m, Q, R, z, c, v of `passing_charge_transfer`.
const PASSING: (f64, f64, f64, f64, f64, f64, f64) = (-1.0, 1.0, 1.0, 4.0, 16.0, 100.0, 80.0);

/// C5, Jackson Problems 13.2–13.3 (Brau §5.2.3, excitation by a fast charged particle): a
/// charge z = 16 passes at v = 0.8c (c = 100) and impact parameter b an electron (q = −1,
/// m = 1) at rest at the centre of a neutral cloud (Q = 1, R = 4: ω₀ = 1/8)
/// (`passing_charge_transfer`). The energy left in the oscillator against Problem 13.3's
/// `ΔE = (2z²q²/(m b² v²)) [ξ² K₁(ξ)² + ξ² K₀(ξ)²/γ²]`, `ξ = ω₀b/(γv)` = 0.3, 1, 2, 4: from
/// nearly the impulse to a free charge (Pr. 13.1) to the adiabatic cut-off, where the
/// transfer falls as e^{−2ξ} (`scripts/wolfram/c5_bound_energy_transfer.wls`, Wolfram Engine
/// 14.2: the closed form and the Fourier integral of the field agree to 28 digits).
/// Required: within 1e-4. Neglected by the problem: the dipole approximation, (swing/b)² ~
/// 1e-9; the electron's speed, (v_e/c)² ~ 1e-10; the projectile's magnetic force, at first
/// order v_e v/c² ~ 1e-5 (the scale c = 100 keeps it there with swings of 1e-2); the field
/// left at the ends, ~1e-5.
#[test]
fn c5_energy_transfer_to_a_bound_charge() {
    let (q, m, _, _, z, _, v) = PASSING;
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const TRANSFER: [(f64, f64); 4] = [
        (0.3, 7.0433438139553840203e-7),
        (1.0, 2.9960611994406950320e-8),
        (2.0, 1.7038348612222998048e-9),
        (4.0, 1.4109718907197246090e-11),
    ];
    let mut worst: f64 = 0.0;
    for (xi, reference) in TRANSFER {
        let (b, energy, steps) = passing_charge_transfer(xi);
        let rel = energy / reference - 1.0;
        let impulse = 2.0 * (z * q).powi(2) / (m * (b * v).powi(2));
        println!(
            "C5 ξ = {xi} (b = {b:.2}): ΔE = {energy:.10e}, Jackson Pr. 13.3 {reference:.10e}: \
             {rel:.1e}; {:.4} of the impulse to a free charge; {steps} steps",
            energy / impulse
        );
        worst = worst.max(rel.abs());
    }
    assert!(worst < 1e-4, "C5 {worst:.3e}");
}

/// C7, Bohr's classical energy loss (Brau §7.3.1; Jackson §13.2, distant collisions): C5's
/// transfer summed over the impact parameters, `∫ 2π b ΔE(b) db`, the energy a passing
/// charge gives per unit length to a medium of one such electron per unit volume. Problem
/// 13.3's ΔE integrates in closed form: `(4πz²q²/(mv²)) [F(ξ₁) − F(ξ₂)]`, `F(ξ) = ξK₀K₁ −
/// (β²/2) ξ² (K₁² − K₀²)` (`scripts/wolfram/c5_bound_energy_transfer.wls`: equal to the
/// direct integral to 30 digits). Measured: the flights of C5 at the 32 nodes of a
/// Gauss–Legendre rule in ln ξ over ξ = 0.1 … 10 (the integrand `2π b² ΔE` in ln ξ is smooth
/// there; beyond 10 it is below e^{−20}). Required: within 1e-4 (as C5). Printed: Bohr's
/// approximation for small ξ₁, `F ≈ ln(1.123/ξ₁) − β²/2`: the adiabatic cut-off acts as a
/// largest impact parameter 1.123 γv/ω₀.
#[test]
fn c7_bohr_energy_loss() {
    let (_, _, big_q, r, z, c, v) = PASSING;
    let (q, m) = (PASSING.0, PASSING.1);
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let beta = v / c;
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const LOSS: f64 = 1.0554403419788837408;
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const F_LOW: f64 = 2.0997318468598346293;
    let (lo, hi) = (0.1_f64.ln(), 10.0_f64.ln());
    let (mut total, mut steps) = (0.0, 0);
    for (x, w) in common::gauss_legendre(32) {
        let xi = (0.5 * (lo + hi) + 0.5 * (hi - lo) * x).exp();
        let (b, energy, k) = passing_charge_transfer(xi);
        total += 0.5 * (hi - lo) * w * 2.0 * PI * b * b * energy;
        steps += k;
    }
    let rel = total / LOSS - 1.0;
    let prefactor = 4.0 * PI * (z * q).powi(2) / (m * v * v);
    println!(
        "C7: ∫ 2π b ΔE db over ξ = 0.1 … 10 = {total:.12e} (closed form {LOSS:.12e}): \
         {rel:.1e}; F(0.1) − F(10) = {:.6} (exact {F_LOW:.6}; Bohr's ln(1.123/ξ) − β²/2 = \
         {:.6}); the cut-off b = 1.123 γv/ω₀ = {:.1}; {steps} steps",
        total / prefactor,
        (1.123_f64 / 0.1).ln() - beta * beta / 2.0,
        1.123 * v / (omega0 * (1.0 - beta * beta).sqrt())
    );
    assert!(rel.abs() < 1e-4, "C7 {rel:.3e}");
}

/// C6, collective radiation damping (Brau §4.3.2 and §10.5; Landau & Lifshitz §75): N
/// electrons (q = −1, m = 1) at their planar equilibrium in a cloud (Q = 4, R = 4:
/// ω₀ = 1/4 for any N), at rest, set swinging along x by a smooth pulse of uniform field
/// (Gaussian, σ = 1, centred at t = 8: zero to 1e-14 at the launch, so the runner's uniform
/// past before the launch is exact; a kick at the launch leaves a kink in each world line,
/// which the reaction's field derivative turns into a worse defect at every other electron
/// and back: the square stalled at t ≈ 7 after a million steps). A uniform field moves only
/// the centre of mass, which oscillates at exactly ω₀ whatever the repulsion (Kohn's
/// theorem, exact for c = ∞): the cluster's dipole plasmon. At c = 3 (ω₀τ = 6.2e-3) each electron radiates (its own
/// Landau–Lifshitz reaction) and feels the others' radiation fields, the term
/// `(2/3c³) d⃛` of the total dipole moment (Landau & Lifshitz §75): the mode decays at N
/// times the single charge's rate Γ₁ = ω₀²τ as the cluster shrinks against the
/// wavelength. Reference, for in-phase dipoles at the cluster's size (separations 0.27–0.41
/// of λ/2π): the classical cooperative decay `Γ₁ Σⱼ F(k r_ij, θ_ij)`, `F(x, θ) =
/// (3/2)[sin²θ (sin x)/x + (1 − 3cos²θ)(cos x/x² − sin x/x³)]` (F = 1 for j = i; θ between
/// the motion and the separation), divided by `(1 + δ)²`: the moving charges' mutual
/// magnetic (Darwin) energy makes the mode heavier by `δ = Σ_{j≠i} q²(1 + cos²θ_ij)/(2mc²
/// r_ij)` (the k² part of the dipoles' coupling; c = 3 makes q²/mc² = 0.11), which also
/// lowers the frequency to `ω₀/√(1 + δ)`. Required, with the exact retarded interaction:
/// the fitted decay rate and the frequency (from the zero crossings) within 2 % of these
/// (Landau–Lifshitz is first order in ω₀τ, as for C2), fitted after the pulse (t > 14).
/// Printed: the default quasi-static
/// interaction, which continues each source's past as motion in a uniform field (no jerk,
/// so none of the others' `(2/3c³) d⃛`).
#[test]
#[allow(clippy::cast_precision_loss)] // small counts
fn c6_collective_radiation_damping() {
    use physics::beam::{BeamParticle, BeamScenario, Fates, run_beam_observed};
    use physics::field::FieldSample;
    struct Pulsed {
        atom: Coulomb,
        e0: f64,
    }
    impl FieldSolver for Pulsed {
        fn sample(&self, x: DVec3, t: f64) -> FieldSample {
            let a = self.atom.sample(x, t);
            let g = (-0.5 * (t - 8.0).powi(2)).exp();
            FieldSample {
                e: a.e + DVec3::new(self.e0 * g, 0.0, 0.0),
                b: a.b,
                phi: a.phi,
            }
        }
    }
    let (q, m, big_q, r, c) = (-1.0_f64, 1.0_f64, 4.0_f64, 4.0_f64, 3.0_f64);
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let gamma1 = omega0 * omega0 * tau;
    let k = omega0 / c;
    // Planar equilibria: a pair on the y axis (d³ = |q|R³/(4Q)), a square (±a, ±a)
    // (a³ = |q|R³ (1/4 + 1/(8√2))/Q).
    let d = (q.abs() * r.powi(3) / (4.0 * big_q)).cbrt();
    let a = (q.abs() * r.powi(3) * (0.25 + 1.0 / (8.0 * 2f64.sqrt())) / big_q).cbrt();
    let clusters: [(&str, Vec<DVec3>); 3] = [
        ("one", vec![DVec3::ZERO]),
        (
            "pair",
            vec![DVec3::new(0.0, d, 0.0), DVec3::new(0.0, -d, 0.0)],
        ),
        (
            "square",
            vec![
                DVec3::new(a, a, 0.0),
                DVec3::new(-a, a, 0.0),
                DVec3::new(-a, -a, 0.0),
                DVec3::new(a, -a, 0.0),
            ],
        ),
    ];
    let f = |x: f64, cos_t: f64| {
        let sin2 = 1.0 - cos_t * cos_t;
        1.5 * (sin2 * x.sin() / x
            + (1.0 - 3.0 * cos_t * cos_t) * (x.cos() / (x * x) - x.sin() / x.powi(3)))
    };
    // The pulse's impulse gives the centre of mass a swing of about 0.05 (σ = 1: the
    // factor e^{−ω₀²σ²/2} = 0.97 for a pulse not much shorter than the period).
    let e0 = 0.05 * omega0 / ((2.0 * PI).sqrt() * (-0.5 * omega0 * omega0).exp());
    let fit_from = 14.0;
    let mut worst: f64 = 0.0;
    for (name, sites) in clusters {
        let n = sites.len();
        // All electrons are equivalent: the sums for the first.
        let (mut sum, mut delta) = (1.0, 0.0);
        for &s in &sites[1..] {
            let sep = s - sites[0];
            let cos_t = sep.x / sep.length();
            sum += f(k * sep.length(), cos_t);
            delta += q * q * (1.0 + cos_t * cos_t) / (2.0 * m * c * c * sep.length());
        }
        let reference = gamma1 * sum / (1.0 + delta).powi(2);
        let frequency = omega0 / (1.0 + delta).sqrt();
        for retarded in [true, false] {
            let kin = Kinematics::new(m, c);
            let scn = BeamScenario {
                field: Pulsed {
                    atom: Coulomb::with_clouds(
                        &[],
                        &[ChargeCloud {
                            position: DVec3::ZERO,
                            charge: big_q,
                            radius: r,
                        }],
                    ),
                    e0,
                },
                obstacles: vec![],
                particles: sites
                    .iter()
                    .map(|&s| BeamParticle {
                        particle: Particle {
                            charge: q,
                            mass: m,
                            radius: 0.0,
                            moment: 0.0,
                        },
                        x0: s,
                        p0: DVec3::ZERO,
                        detector: None,
                        acceptance: None,
                    })
                    .collect(),
                c,
                bounds: None,
                t_max: fit_from + 3.0 / reference,
                interact: true,
                retarded,
                gates: Vec::new(),
                radiation_reaction: true,
                fates: Fates::default(),
            };
            let mut pts: Vec<(f64, f64)> = Vec::new();
            let mut crossings: Vec<f64> = Vec::new();
            let run = run_beam_observed(
                &scn,
                &RunSettings::with_tolerance(TOL),
                |dense, members, p_ref| {
                    if members.len() != n {
                        return;
                    }
                    let mean_x = |t: f64| {
                        (0..n).map(|k| dense.eval_component(6 * k, t)).sum::<f64>() / n as f64
                    };
                    let (t0, t1) = (dense.t_start(), dense.t_end());
                    // Zero crossings of the centre of mass, bisected on the dense output.
                    if t0 > fit_from && mean_x(t0).signum() != mean_x(t1).signum() {
                        let (mut lo, mut hi) = (t0, t1);
                        for _ in 0..60 {
                            let mid = 0.5 * (lo + hi);
                            if mean_x(mid).signum() == mean_x(lo).signum() {
                                lo = mid;
                            } else {
                                hi = mid;
                            }
                        }
                        crossings.push(0.5 * (lo + hi));
                    }
                    if t1 < fit_from {
                        return;
                    }
                    let (mut x, mut v) = (DVec3::ZERO, DVec3::ZERO);
                    for k in 0..n {
                        let comp = |i: usize| dense.eval_component(6 * k + i, t1);
                        x += DVec3::new(comp(0), comp(1), 0.0);
                        v += kin.velocity(DVec3::new(comp(3), comp(4), 0.0) * p_ref);
                    }
                    let (x, v) = (x / n as f64, v / n as f64);
                    // The mode's energy, with the mass the oscillation shows.
                    let energy = 0.5
                        * n as f64
                        * m
                        * ((1.0 + delta) * v.length_squared()
                            + omega0 * omega0 * x.length_squared());
                    pts.push((t1, energy.ln()));
                },
            );
            assert!(
                run.trajectories
                    .iter()
                    .all(|t| t.outcome == Outcome::Timeout),
                "C6 {name}: {:?}",
                run.trajectories
                    .iter()
                    .map(|t| (t.outcome, t.end.t))
                    .collect::<Vec<_>>()
            );
            let np = pts.len() as f64;
            let (mt, my) = (
                pts.iter().map(|p| p.0).sum::<f64>() / np,
                pts.iter().map(|p| p.1).sum::<f64>() / np,
            );
            let slope = pts.iter().map(|&(t, y)| (t - mt) * (y - my)).sum::<f64>()
                / pts.iter().map(|&(t, _)| (t - mt).powi(2)).sum::<f64>();
            let rel = -slope / reference - 1.0;
            let half_periods = (crossings.len() - 1) as f64;
            let measured_w = PI * half_periods / (crossings[crossings.len() - 1] - crossings[0]);
            let rel_w = measured_w / frequency - 1.0;
            println!(
                "C6 {name} (N = {n}), {}: decay rate {:.4} Γ₁ (reference {:.4} Γ₁: cooperative \
                 {sum:.4}, Darwin δ = {delta:.4}; point limit {n}): {rel:.2e}; frequency \
                 {:.6} ω₀ (reference {:.6}): {rel_w:.1e}; {} samples, {} crossings",
                if retarded {
                    "exact retarded"
                } else {
                    "quasi-static"
                },
                -slope / gamma1,
                reference / gamma1,
                measured_w / omega0,
                frequency / omega0,
                pts.len(),
                crossings.len()
            );
            if retarded {
                worst = worst.max(rel.abs()).max(rel_w.abs());
            }
        }
    }
    println!(
        "C6: ω₀ = {omega0}, ω₀τ = {:.2e}, kd = {:.3}, ka = {:.3}",
        omega0 * tau,
        k * d,
        k * a
    );
    assert!(worst < 0.02, "C6 {worst:.3e}");
}

/// C8, the polarization force of an atom (Brau §3.1.3 and Ex. 3.6): a slow charge z = 1/800
/// (m = 1) passes a neutral Thomson atom (cloud Q = 1, R = 4; an electron q = −1, m = 1 at
/// rest at its centre: ω₀ = 1/8, polarizability α = q²/(mω₀²) = R³ = 64) at impact parameter
/// b = 20, from 20 b before the closest approach to 20 b after it, both particles flown by
/// the beam runner (c = ∞). The charge displaces the electron and is attracted by the
/// dipole it induced: in the adiabatic limit by the potential `−α z²/(2r⁴)`, an impulse
/// `−3π α z²/(4 v b⁴)` toward the atom. Reference, at ξ = ω₀b/v = 10, 25, 50: the exact
/// linear response for the straight path, the oscillator driven from rest and its dipole's
/// force integrated over the same flight (`scripts/wolfram/c8_polarization_force.wls`,
/// Wolfram Engine 14.2, NDSolve at 30 digits), which approaches the adiabatic impulse as
/// ~1.9/ξ². Required: within 1e-4 (neglected: the induced charge's quadrupole, at first
/// order in swing/b = 1e-5; the charge's deflection, 6e-7 rad).
#[test]
fn c8_polarization_force_of_an_atom() {
    use physics::beam::{BeamParticle, BeamScenario, Fates, run_beam};
    let (q, m, big_q, r) = (-1.0_f64, 1.0_f64, 1.0_f64, 4.0_f64);
    let (z, mass, b) = (1.0 / 800.0, 1.0_f64, 20.0_f64);
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let alpha = q * q / (m * omega0 * omega0);
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const IMPULSE: [(f64, f64); 3] = [
        (10.0, -6.0094856779354391829e-9),
        (25.0, -1.4770866642470541742e-8),
        (50.0, -2.9474575648422653079e-8),
    ];
    let mut worst: f64 = 0.0;
    for (xi, reference) in IMPULSE {
        let v = omega0 * b / xi;
        let particle = |charge: f64, mass: f64, x0: DVec3, p0: DVec3| BeamParticle {
            particle: Particle {
                charge,
                mass,
                radius: 0.0,
                moment: 0.0,
            },
            x0,
            p0,
            detector: None,
            acceptance: None,
        };
        let scn = BeamScenario {
            field: Coulomb::with_clouds(
                &[],
                &[ChargeCloud {
                    position: DVec3::ZERO,
                    charge: big_q,
                    radius: r,
                }],
            ),
            obstacles: vec![],
            particles: vec![
                particle(
                    z,
                    mass,
                    DVec3::new(-20.0 * b, b, 0.0),
                    DVec3::new(mass * v, 0.0, 0.0),
                ),
                particle(q, m, DVec3::ZERO, DVec3::ZERO),
            ],
            c: f64::INFINITY,
            bounds: None,
            t_max: 40.0 * b / v,
            interact: true,
            retarded: false,
            gates: Vec::new(),
            radiation_reaction: false,
            fates: Fates::default(),
        };
        let run = run_beam(&scn, &RunSettings::with_tolerance(TOL));
        let got = run.trajectories[0].end.p.y;
        let rel = got / reference - 1.0;
        let adiabatic = -3.0 * PI * alpha * z * z / (4.0 * v * b.powi(4));
        println!(
            "C8 ξ = {xi} (v = {v}): impulse {got:.10e} (linear response {reference:.10e}): \
             {rel:.1e}; {:.4} of the adiabatic −3παz²/(4vb⁴); {} steps",
            got / adiabatic,
            run.stats.n_accept
        );
        worst = worst.max(rel.abs());
    }
    assert!(worst < 1e-4, "C8 {worst:.3e}");
}
