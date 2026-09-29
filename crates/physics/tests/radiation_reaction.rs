//! Validation tests R1–R3 for radiation reaction (PHYSICS.md §3.1, §9). Run with
//! `cargo test -p physics --test radiation_reaction -- --nocapture --test-threads=1`.

mod common;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle, ParticleOde};
use physics::field::{Coulomb, FixedCharge, UniformFields};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{RunSettings, Scenario, run, run_observed};

const TOL: f64 = 1e-12;

// --- R1: synchrotron damping in a uniform magnetic field ----------------------------------

/// Motion perpendicular to a uniform `B`. The Landau–Lifshitz force reduces exactly to
/// `−κ γ² v` with `κ = 2q⁴B²/(3m²c³)` (the derivative term vanishes, `E = 0`), so
/// `du/dt = −(κ/m) u sqrt(1 + u²)` for `u = |p|/(mc)`, with the exact solution
/// `u(t) = 1 / sinh(asinh(1/u₀) + κt/m)`: relativistic synchrotron damping.
#[test]
fn r1_synchrotron_damping_exact() {
    let (q, m, c, b): (f64, f64, f64, f64) = (1.0, 1.0, 2.0, 0.2);
    let kappa = 2.0 * q.powi(4) * b * b / (3.0 * m * m * c.powi(3));
    let mut worst: f64 = 0.0;
    for gamma0 in [1.05, 2.0, 10.0] {
        let u0 = (gamma0 * gamma0 - 1.0_f64).sqrt();
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
            x0: DVec3::ZERO,
            p0: DVec3::new(u0 * m * c, 0.0, 0.0),
            detector: None,
            bounds: Some(cube(1e6)),
            t_max: 600.0,
            radiation_reaction: true,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(TOL));
        let exact = |t: f64| 1.0 / ((1.0 / u0).asinh() + kappa / m * t).sinh();
        let mut dev: f64 = 0.0;
        for s in &tr.samples {
            let u = s.p.length() / (m * c);
            dev = dev.max((u - exact(s.t)).abs() / u0);
        }
        let u_end = tr.end.p.length() / (m * c);
        let balance = tr.energy_max_abs_error / tr.kinetic_initial;
        println!(
            "R1 γ₀ = {gamma0}: u {u0:.3} → {u_end:.4} over t = {:.0} ({} steps); max |Δu|/u₀ = {dev:.2e}; \
             energy balance |ΔW − work|/T₀ = {balance:.2e}; max |F_RR|/|F_L| = {:.2e}",
            tr.end.t, tr.stats.n_accept, tr.reaction_ratio_max
        );
        assert!(u_end < 0.5 * u0, "damping too weak to be a meaningful test");
        assert!(balance < 1e-10, "energy balance {balance:.3e}");
        worst = worst.max(dev);
    }
    assert!(worst < 1e-9, "R1 error {worst:.3e}");
}

// --- R2: radiated energy in Coulomb scattering --------------------------------------------

/// Scattering off a fixed charge. The work of the Landau–Lifshitz force equals minus the
/// Larmor–Liénard energy `∫P dt` plus the change of the Schott energy
/// `E_S = τ₀ m γ⁴ (v·a)` (τ₀ = 2q²/(3mc³)), which is not zero because the flight starts
/// and ends at a finite distance, up to corrections of relative order `τ₀ ω` (ω = v/b, the
/// inverse collision time) that are beyond the accuracy of the LL approximation itself.
/// The remaining difference must shrink with `τ₀`. (Without the Schott term the
/// difference is 1e-4 at every c and falls as 1/distance², which identifies it.)
#[test]
fn r2_radiated_energy_matches_larmor_in_scattering() {
    let mut rows = Vec::new();
    for c in [5.0, 10.0, 20.0] {
        let kappa = 1.0; // q Q
        let charge = FixedCharge {
            position: DVec3::ZERO,
            charge: kappa,
            radius: 0.05,
        };
        let (b, v) = (0.6, 0.8);
        let scn = Scenario {
            field: Coulomb::new(&[charge]),
            obstacles: vec![Shape::Sphere(Sphere {
                center: DVec3::ZERO,
                radius: 0.05,
            })],
            particle: UNIT_PARTICLE,
            c,
            x0: DVec3::new(-400.0, b, 0.0),
            p0: DVec3::new(
                Kinematics::new(1.0, c)
                    .momentum_from_kinetic_energy(0.5 * v * v, DVec3::X)
                    .x,
                0.0,
                0.0,
            ),
            detector: None,
            bounds: Some(cube(400.5)),
            t_max: 5000.0,
            radiation_reaction: true,
            acceptance: None,
            gates: Vec::new(),
        };
        // Accurate Larmor–Liénard integral: Simpson's rule on every step's dense output.
        let mut larmor = 0.0;
        let tr = run_observed(&scn, &RunSettings::with_tolerance(TOL), |step| {
            let ode: &ParticleOde<_> = step.ode;
            let power = |t: f64| {
                let (x, p) = step.state(t);
                let gamma = ode.kin.gamma(p);
                let vel = ode.kin.velocity(p);
                let f = ode.force(x, p, t);
                let a = (f - vel * (vel.dot(f) / (c * c))) / (gamma * ode.kin.mass);
                2.0 / 3.0
                    * gamma.powi(6)
                    * (a.length_squared() - vel.cross(a).length_squared() / (c * c))
                    / c.powi(3)
            };
            let (t0, t1) = (step.t_start(), step.t_end());
            let n = 16;
            let h = (t1 - t0) / f64::from(n);
            let mut sum = power(t0) + power(t1);
            for i in 1..n {
                sum += power(t0 + h * f64::from(i)) * if i % 2 == 1 { 4.0 } else { 2.0 };
            }
            larmor += sum * h / 3.0;
        });
        let tau0 = 2.0 / (3.0 * c.powi(3));
        let ode = ParticleOde::new(&scn.field, &scn.particle, c, 1.0);
        let schott = |t: f64, x: DVec3, p: DVec3| {
            let gamma = ode.kin.gamma(p);
            let vel = ode.kin.velocity(p);
            let f = ode.force(x, p, t);
            let a = (f - vel * (vel.dot(f) / (c * c))) / (gamma * ode.kin.mass);
            tau0 * ode.kin.mass * gamma.powi(4) * vel.dot(a)
        };
        let d_schott = schott(tr.end.t, tr.end.x, tr.end.p) - schott(0.0, scn.x0, scn.p0);
        let rel = (tr.radiation_work - (-larmor + d_schott)).abs() / larmor;
        println!(
            "R2 c = {c}: outcome {:?}, LL work {:.6e}, ∫P dt {larmor:.6e}, rel. diff {rel:.2e}, τ₀ω = {:.2e}, \
             |ΔW − work|/T₀ = {:.1e}",
            tr.outcome,
            tr.radiation_work,
            tau0 * v / b,
            tr.energy_max_abs_error / tr.kinetic_initial
        );
        assert!(tr.energy_max_abs_error / tr.kinetic_initial < 1e-10);
        rows.push((tau0 * v / b, rel));
    }
    for (tw, rel) in &rows {
        assert!(rel < tw, "difference {rel:.3e} not below τ₀ω = {tw:.3e}");
    }
}

// --- R3: no reaction in Newtonian mechanics -----------------------------------------------

/// `c = ∞`: no radiation, so switching radiation reaction on changes nothing, bit for bit.
#[test]
#[allow(clippy::float_cmp)] // exact zero by construction
fn r3_no_reaction_without_finite_c() {
    let make = |rr: bool| Scenario {
        field: UniformFields {
            e: DVec3::new(0.1, 0.0, 0.0),
            b: DVec3::new(0.0, 0.0, 0.3),
        },
        obstacles: vec![],
        particle: UNIT_PARTICLE,
        c: f64::INFINITY,
        x0: DVec3::ZERO,
        p0: DVec3::new(0.0, 1.0, 0.0),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 50.0,
        radiation_reaction: rr,
        acceptance: None,
        gates: Vec::new(),
    };
    let (a, b) = (
        run(&make(false), &RunSettings::with_tolerance(TOL)),
        run(&make(true), &RunSettings::with_tolerance(TOL)),
    );
    assert_eq!(a.samples, b.samples);
    assert_eq!(b.radiation_work, 0.0);
}

/// R4, Jackson Problem 16.2 (the classical atom): a charge −1 on a circular orbit around a
/// fixed charge +Z radiates and spirals in. Nonrelativistically the radius obeys
/// `r³ = r₀³ − 9 Z (cτ)² c t` with `τ = 2q²/(3mc³)`, i.e. `d(r³)/dt = −6 Z τ` (units with
/// Coulomb's constant 1, q = m = 1). At c = 8 (v/c = 0.056) the fitted slope of r³ over
/// ~12 orbits must agree within 1 %: the relativistic correction is O((v/c)²) ≈ 0.3 %, and
/// the slight eccentricity of the launch (Newtonian circular speed) scatters the samples.
#[test]
fn r4_classical_atom_orbit_decay() {
    let (z, c, r0) = (1.0_f64, 8.0_f64, 5.0_f64);
    let tau = 2.0 / (3.0 * c.powi(3));
    let nucleus = FixedCharge {
        position: DVec3::ZERO,
        charge: z,
        radius: 0.05,
    };
    let v = (z / r0).sqrt();
    let scn = Scenario {
        field: Coulomb::new(&[nucleus]),
        obstacles: vec![Shape::Sphere(Sphere {
            center: DVec3::ZERO,
            radius: 0.05,
        })],
        particle: Particle {
            charge: -1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(0.0, -r0, 0.0),
        p0: Kinematics::new(1.0, c).momentum_from_kinetic_energy(0.5 * v * v, DVec3::X),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 1500.0,
        radiation_reaction: true,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    // Least-squares slope of r³ against t over every accepted step.
    let pts: Vec<(f64, f64)> = tr
        .samples
        .iter()
        .map(|s| (s.t, s.x.length().powi(3)))
        .collect();
    #[allow(clippy::cast_precision_loss)]
    let n = pts.len() as f64;
    let (mt, my) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let slope = pts.iter().map(|&(t, y)| (t - mt) * (y - my)).sum::<f64>()
        / pts.iter().map(|&(t, _)| (t - mt).powi(2)).sum::<f64>();
    let jackson = -6.0 * z * tau;
    println!(
        "R4 c = {c}: d(r³)/dt = {slope:.6} (fit over {} samples), Jackson {jackson:.6}, \
         ratio {:.4}; max |F_RR|/|F_L| = {:.2e}",
        pts.len(),
        slope / jackson,
        tr.reaction_ratio_max
    );
    assert!((slope / jackson - 1.0).abs() < 0.01);
}

/// R5, Jackson Problem 14.5b: a charge q (mass m) fired head-on at a repulsive fixed
/// charge Q (v₀ = 0.8 at the launch, 400 away, c = 20: v/c = 0.04) radiates, nonrelativistically,
/// `(8/45)(q/Q) m v∞⁵/c³`. Exact reference (`scripts/wolfram/r5_head_on_collision.wls`,
/// Wolfram Engine 14.2): along a line `a = F/(γ³m)`, so Liénard's power is
/// `(2q²/3m²c³)(qQ/r²)²` at any speed, and its integral along the exact relativistic
/// (radiation-free) trajectory, in from the launch and out to the arena's edge, is
/// `W = 7.4254746851720199e-6` (Jackson's formula is 8.6e-5 lower: the relativistic
/// correction). Checks, targets set before measuring:
/// (1) without radiation reaction, the recorded Liénard energy equals W to 1e-6;
/// (2) with it, the Landau–Lifshitz work plus the change of the Schott term: in one
///     dimension LL reduces to `(2q³/3mc³) γ DE/Dt` (the other two terms cancel), whose
///     work is `−∫P dt + (2q³/3mc³)[γ v E]` exactly: equals the recorded Liénard energy to
///     1e-6 (the boundary term, from the weak Coulomb field at the ends, is 1.1e-4 of W);
/// (3) with it, W differs from the radiation-free W by the back-reaction, O(W/T₀) ~ 3e-5:
///     below 1e-4;
/// (4) Jackson's nonrelativistic formula within 1e-3 of W.
/// (The first version compared the recorded energy, then integrated by the trapezoidal
/// rule over steps, with Jackson's formula: 0.24 %, attributed to relativity; the exact
/// reference showed relativity accounts for 8.6e-5 and the quadrature for the rest.)
#[test]
fn r5_head_on_collision_radiates_jacksons_energy() {
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const W_EXACT: f64 = 7.4254746851720199e-6;
    let (q, big_q, m, v0, c) = (1.0_f64, 1.0_f64, 1.0_f64, 0.8_f64, 20.0_f64);
    let centre = FixedCharge {
        position: DVec3::ZERO,
        charge: big_q,
        radius: 0.05,
    };
    let scn = |radiation_reaction: bool| Scenario {
        field: Coulomb::new(&[centre]),
        obstacles: vec![Shape::Sphere(Sphere {
            center: DVec3::ZERO,
            radius: 0.05,
        })],
        particle: Particle {
            charge: q,
            mass: m,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(-400.0, 0.0, 0.0),
        p0: Kinematics::new(m, c).momentum_from_kinetic_energy(0.5 * m * v0 * v0, DVec3::X),
        detector: None,
        bounds: Some(cube(400.5)),
        t_max: 5000.0,
        radiation_reaction,
        acceptance: None,
        gates: Vec::new(),
    };
    let free = run(&scn(false), &RunSettings::with_tolerance(TOL));
    let with_rr = run(&scn(true), &RunSettings::with_tolerance(TOL));
    // The Schott boundary term (2q³/3mc³)[γ v E] between the launch and the end.
    let kin = Kinematics::new(m, c);
    let schott = |x: DVec3, p: DVec3| {
        let e_field = big_q * x / x.length().powi(3);
        2.0 * q.powi(3) / (3.0 * m * c.powi(3)) * kin.gamma(p) * kin.velocity(p).dot(e_field)
    };
    let boundary = schott(with_rr.end.x, with_rr.end.p) - schott(scn(true).x0, scn(true).p0);
    let ll = -with_rr.radiation_work + boundary;
    let e1 = free.radiated_energy / W_EXACT - 1.0;
    let e2 = ll / with_rr.radiated_energy - 1.0;
    let e3 = with_rr.radiated_energy / W_EXACT - 1.0;
    // Jackson's v₀ is the speed at infinity: the launch at distance 400 already has the
    // potential energy qQ/400.
    let v_inf = (v0 * v0 + 2.0 * q * big_q / (m * 400.0)).sqrt();
    let jackson = 8.0 / 45.0 * (q / big_q) * m * v_inf.powi(5) / c.powi(3);
    let e4 = jackson / W_EXACT - 1.0;
    println!(
        "R5: (1) without RR {:.9e} vs exact {W_EXACT:.9e}: {e1:.1e}; (2) LL work + Schott {ll:.9e} vs Liénard {:.9e}: {e2:.1e} (boundary {:.1e} of W); (3) back-reaction {e3:.1e}; (4) Jackson {e4:.1e}",
        free.radiated_energy,
        with_rr.radiated_energy,
        boundary / W_EXACT
    );
    assert!(e1.abs() < 1e-6 && e2.abs() < 1e-6 && e3.abs() < 1e-4 && e4.abs() < 1e-3);
}
