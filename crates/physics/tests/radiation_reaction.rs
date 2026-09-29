//! Validation tests R1–R6 for radiation reaction (PHYSICS.md §3.1, §9). Run with
//! `cargo test -p physics --test radiation_reaction -- --nocapture --test-threads=1`.

mod common;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle, ParticleOde};
use physics::field::{Coulomb, FieldSample, FieldSolver, FixedCharge, UniformFields};
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
/// fixed charge +Z radiates and spirals in; nonrelativistically `d(r³)/dt = −6 Z τ`,
/// `τ = 2q²/(3mc³)` (units with Coulomb's constant 1, q = m = 1). Exact reference at c = 8
/// (v/c = 0.056; `scripts/wolfram/r4_classical_atom.wls`, Wolfram Engine 14.2): the
/// adiabatic inspiral through relativistic circular orbits (γv² = Z/r), `dr/dt = −P/E′`
/// with `E = γc² − Z/r` and the relativistic Larmor power `(2/3c³)γ⁴(v²/r)²`, from r₀ = 5
/// to T = 1500: r(T) = 4.8377986383588659, a mean d(r³)/dt 0.48 % above Jackson's. The
/// particle is launched on that inspiral (the exact circular momentum and the inspiral's
/// radial velocity), so no eccentricity is excited. Required: r³(T) − r₀³ within 1e-4 of
/// the exact value (set before measuring: the neglected terms, the residual eccentricity
/// and Landau–Lifshitz's O((τω)²), are 1e-7 or less); Jackson's slope within 1 % of it.
/// (The first version launched with the Newtonian circular speed and fitted the slope of
/// r³ against Jackson's: 0.64 %, attributed to relativity plus the launch's eccentricity;
/// relativity is 0.48 %.)
#[test]
fn r4_classical_atom_orbit_decay() {
    // From the Wolfram script: tangential momentum γv, radial velocity dr/dt at r₀, r(T).
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const LAUNCH: (f64, f64, f64) = (
        0.44756311749347263,
        -0.00010465622167041062,
        4.8377986383588659,
    );
    let (z, c, r0, t_end) = (1.0_f64, 8.0_f64, 5.0_f64, 1500.0_f64);
    let tau = 2.0 / (3.0 * c.powi(3));
    let nucleus = FixedCharge {
        position: DVec3::ZERO,
        charge: z,
        radius: 0.05,
    };
    let (u, rdot, r_exact) = LAUNCH;
    // At (0, −r₀) moving +x: the radial direction is −y; γ is that of the tangential
    // motion (the radial velocity changes it by 1e-8).
    let gamma = (1.0 + u * u / (c * c)).sqrt();
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
        p0: DVec3::new(u, -gamma * rdot, 0.0),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: t_end,
        radiation_reaction: true,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    let d_measured = tr.end.x.length().powi(3) - r0.powi(3);
    let d_exact = r_exact.powi(3) - r0.powi(3);
    let err = d_measured / d_exact - 1.0;
    let jackson = -6.0 * z * tau * t_end;
    println!(
        "R4 c = {c}: r³(T) − r₀³ = {d_measured:.9} (t = {:.1}), exact {d_exact:.9}: {err:.1e}; Jackson {jackson:.6} ({:.2e} from exact); max |F_RR|/|F_L| = {:.2e}",
        tr.end.t,
        jackson / d_exact - 1.0,
        tr.reaction_ratio_max
    );
    assert!(err.abs() < 1e-4 && (jackson / d_exact - 1.0).abs() < 0.01);
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

// --- R6: a gap with radiation damping (Jackson Pr. 16.10–16.11) --------------------------

/// A gap of uniform field `E₀` between two charged grids at x = 0 and x = d, smeared over a
/// width w: `E_x = E₀ (tanh(x/w) − tanh((x − d)/w))/2`, with its potential.
struct Gap {
    e0: f64,
    d: f64,
    w: f64,
}

impl FieldSolver for Gap {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        let (u, v) = (x.x / self.w, (x.x - self.d) / self.w);
        let ln_cosh = |z: f64| z.abs() + (-2.0 * z.abs()).exp().ln_1p() - std::f64::consts::LN_2;
        FieldSample {
            e: DVec3::new(0.5 * self.e0 * (u.tanh() - v.tanh()), 0.0, 0.0),
            b: DVec3::ZERO,
            phi: -0.5 * self.e0 * self.w * (ln_cosh(u) - ln_cosh(v)),
        }
    }
}

/// R6, Jackson Pr. 16.10–16.11: a charge crosses a gap of uniform field (`Gap`) from
/// x = −20w to beyond d + 20w (where the field is below 5e-18 of the gap's). In one
/// dimension the Lorentz–Dirac equation, in the rapidity y (p = mc sinh y) and the proper
/// time s, is exactly the Abraham–Lorentz equation `mc (y′ − τ y″) = f`, `f = qE(x(s))`
/// (Pr. 16.8), whose physical (non-runaway) solution is the integro-differential form of
/// Pr. 16.10(a), `mc y′(s) = ∫₀^∞ e^{−u} f(s + τu) du`: the series `f + τ f′ + τ² f″ + …` of
/// Pr. 16.10(b), whose first two terms are exactly Landau–Lifshitz here. References
/// (`scripts/wolfram/r6_gap_damping.wls`, Wolfram Engine 14.2, 32 digits): the exact
/// solution (integrated backward in s from beyond the gap, where the runaway mode decays,
/// shooting on the exit rapidity; two methods agree), the Landau–Lifshitz flight and the
/// flight without radiation reaction, as the arrival time at x = d + 20w and the momentum
/// there. From them the effects of the damping: on the transit time T (from x = 0 to d,
/// the uniform motions outside the gap extrapolated to its edges) and on the exit velocity.
/// Cases: (A) nonrelativistic (c = 20, v 1/4 → 1/2 over d = 10) with sharp edges
/// (w = 0.01); (B) relativistic (c = 5, β 0.2 → 0.42 over d = 1) with steep edges
/// (w = 0.2, 0.1, 0.05), where Landau–Lifshitz's error is measurable. Required (set before
/// measuring):
/// (1) the flight without the reaction within 1e-10 of its reference;
/// (2) the effects of the reaction within 1e-4 (A) and 1e-7 (B) of the Landau–Lifshitz
///     reference's (the integrator's error, ~1e-11 of the flight, against effects of 1.5e-6
///     and 1e-3 of it);
/// (3) the Landau–Lifshitz effects within the reaction ratio shown in the flight details
///     (largest reaction force over the largest Lorentz force) of the exact ones: the model
///     note's claim that the force is valid while that ratio is small;
/// (4) in A, Jackson's first-order formulas of Pr. 16.11(b), `T′ = T − τ(1 − v₀/v₁)` and
///     `v₁′ = v₁ − (a²τ/v₁) T`, within 5e-3 (their corrections: the smooth edges, w/d = 1e-3,
///     and relativity, β² = 6e-4);
/// (5) in A, Pr. 16.11(c): the radiated (Liénard) energy plus the change of kinetic energy
///     equals the field's work, within 1e-5 of the radiated energy (Landau–Lifshitz's own
///     violation, of relative order τ² ∫Ḟ² dt / ∫F² dt, is 1e-8 here).
#[test]
fn r6_gap_with_radiation_damping() {
    // (c, E₀, d, w, v₀, [t, p] without reaction, [t, p] Landau–Lifshitz, [t, p] exact):
    // arrival time at x = d + 20w (t = 0 at x = −20w) and the momentum there.
    type Case = (f64, f64, f64, f64, f64, [f64; 2], [f64; 2], [f64; 2]);
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const REF: [Case; 4] = [
        (
            20.0,
            0.009375,
            10.0,
            0.01,
            0.25,
            [27.869597464329432, 0.50004638780987578],
            [27.869556084565951, 0.50004599749692859],
            [27.869556084498039, 0.50004599749695609],
        ),
        (
            5.0,
            2.0,
            1.0,
            0.2,
            1.0,
            [6.5340664296895947, 2.2987278158375551],
            [6.5344411256767111, 2.2936159769216889],
            [6.5344271975677316, 2.293619210186051],
        ),
        (
            5.0,
            2.0,
            1.0,
            0.1,
            1.0,
            [3.590821355683222, 2.2987278158375551],
            [3.589787367228894, 2.2929414230177335],
            [3.5897716776893105, 2.2929466760068178],
        ),
        (
            5.0,
            2.0,
            1.0,
            0.05,
            1.0,
            [2.116261248019058, 2.2987278158375551],
            [2.1143762296548531, 2.2925898407025205],
            [2.1143591587773803, 2.2925985112027798],
        ),
    ];
    let (q, m) = (1.0_f64, 1.0_f64);
    for (case, &(c, e0, d, w, v0, free_ref, ll_ref, exact_ref)) in REF.iter().enumerate() {
        let (l, tau) = (20.0 * w, 2.0 * q * q / (3.0 * m * c.powi(3)));
        let kin = Kinematics::new(m, c);
        let gap = Gap { e0, d, w };
        let fly = |radiation_reaction: bool| {
            let scn = Scenario {
                field: Gap { e0, d, w },
                obstacles: vec![],
                particle: Particle {
                    charge: q,
                    mass: m,
                    radius: 0.0,
                    moment: 0.0,
                },
                c,
                x0: DVec3::new(-l, 0.0, 0.0),
                p0: DVec3::new(kin.gamma_of_velocity(DVec3::X * v0) * m * v0, 0.0, 0.0),
                detector: None,
                bounds: None,
                t_max: free_ref[0] + 1.0,
                radiation_reaction,
                acceptance: None,
                gates: Vec::new(),
            };
            let tr = run(&scn, &RunSettings::with_tolerance(TOL));
            // The arrival at x = d + l, from the uniform motion beyond it.
            let v = kin.velocity(tr.end.p).x;
            ([tr.end.t - (tr.end.x.x - (d + l)) / v, tr.end.p.x], tr)
        };
        // Transit time (from x = 0 to x = d) and exit velocity of an arrival [t, p].
        let transit = |[t, p]: [f64; 2]| {
            let v1 = kin.velocity(DVec3::X * p).x;
            (t - l / v1 - l / v0, v1)
        };
        let effect = |with: [f64; 2], without: [f64; 2]| {
            let ((t1, v1), (t0, v0)) = (transit(with), transit(without));
            (t1 - t0, v1 - v0)
        };
        let (free, _) = fly(false);
        let (ll, tr) = fly(true);
        let e1 = ((free[0] / free_ref[0] - 1.0).abs()).max((free[1] / free_ref[1] - 1.0).abs());
        let ours = effect(ll, free);
        let reference = effect(ll_ref, free_ref);
        let exact = effect(exact_ref, free_ref);
        let rel = |a: (f64, f64), b: (f64, f64)| ((a.0 / b.0 - 1.0).abs(), (a.1 / b.1 - 1.0).abs());
        let e2 = rel(ours, reference);
        let e3 = rel(ours, exact);
        let ratio = tr.reaction_ratio_max;
        println!(
            "R6 c = {c}, w = {w}: without reaction {e1:.1e} from its reference; effects \
             ΔT = {:.6e}, Δv₁ = {:.6e}: {:.1e}, {:.1e} from the Landau–Lifshitz reference; \
             {:.1e}, {:.1e} from the exact (reaction ratio {ratio:.1e})",
            ours.0, ours.1, e2.0, e2.1, e3.0, e3.1
        );
        let e2_max = if case == 0 { 1e-4 } else { 1e-7 };
        assert!(
            e1 < 1e-10,
            "R6 case {case}: flight without reaction {e1:.3e}"
        );
        assert!(e2.0 < e2_max && e2.1 < e2_max, "R6 case {case}: {e2:?}");
        assert!(
            e3.0 < ratio && e3.1 < ratio,
            "R6 case {case}: {e3:?} against {ratio:.3e}"
        );
        if case == 0 {
            // (4) Jackson's formulas, with T, v₀, v₁ of the flight without reaction.
            let (t_free, v1) = transit(free);
            let a = q * e0 / m;
            let dt_jackson = -tau * (1.0 - v0 / v1);
            let dv_jackson = -(a * a * tau / v1) * t_free;
            let j = (ours.0 / dt_jackson - 1.0, ours.1 / dv_jackson - 1.0);
            // (5) Energy balance: kinetic energy gained plus radiated equals the work.
            let work =
                q * (gap.sample(DVec3::new(-l, 0.0, 0.0), 0.0).phi - gap.sample(tr.end.x, 0.0).phi);
            let gained = kin.kinetic_energy(tr.end.p) - tr.kinetic_initial;
            let balance = (gained + tr.radiated_energy - work) / tr.radiated_energy;
            println!(
                "R6 (A): T′ − T = {:.6e} (Jackson {dt_jackson:.6e}: {:.1e}), v₁′ − v₁ = {:.6e} \
                 (Jackson {dv_jackson:.6e}: {:.1e}); radiated {:.6e}, balance {balance:.1e}",
                ours.0, j.0, ours.1, j.1, tr.radiated_energy
            );
            assert!(j.0.abs() < 5e-3 && j.1.abs() < 5e-3, "R6 Jackson {j:?}");
            assert!(balance.abs() < 1e-5, "R6 energy balance {balance:.3e}");
        }
    }
}
