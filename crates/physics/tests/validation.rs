//! Physics validation tests T1–T10 (PHYSICS.md §9). Each compares against an analytic
//! result and prints the measured error; run with
//! `cargo test -p physics --test validation -- --nocapture --test-threads=1`.

#![allow(clippy::cast_precision_loss)] // small loop counters

mod common;

use std::f64::consts::PI;

use common::{
    CoulombOrbit, Rng, UNIT_PARTICLE, angle_between, central_charge, cube, periapsis_times,
};
use physics::DVec3;
use physics::dynamics::Kinematics;
use physics::field::{Coulomb, FixedCharge, UniformElectric};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{Outcome, RunSettings, Scenario, run, run_observed};

/// Tolerance used for the "verify" level in these tests.
const TOL: f64 = 1e-12;

// --- T1: energy conservation in a many-charge field ---------------------------------------

fn many_charge_scenario(c: f64, seed: u64) -> Scenario<Coulomb> {
    let mut rng = Rng::new(seed);
    let charges: Vec<FixedCharge> = (0..20)
        .map(|_| FixedCharge {
            position: DVec3::new(
                rng.range(-10.0, 10.0),
                rng.range(-10.0, 10.0),
                rng.range(-2.0, 2.0),
            ),
            charge: rng.range(-3.0, 3.0),
            radius: 0.2,
        })
        .collect();
    Scenario {
        field: Coulomb::new(&charges),
        obstacles: charges
            .iter()
            .map(|c| {
                Shape::Sphere(Sphere {
                    center: c.position,
                    radius: c.radius,
                })
            })
            .collect(),
        particle: UNIT_PARTICLE,
        c,
        x0: DVec3::new(-14.0, 0.3, 0.1),
        p0: DVec3::new(2.0, 0.0, 0.0),
        detector: None,
        bounds: Some(cube(30.0)),
        t_max: 200.0,
        radiation_reaction: false,
        acceptance: None,
    }
}

#[test]
fn t1_energy_conservation_many_charges() {
    for c in [f64::INFINITY, 5.0] {
        for seed in 1..=4 {
            let scn = many_charge_scenario(c, seed);
            let tr = run(&scn, &RunSettings::with_tolerance(TOL));
            let rel = tr.energy_max_abs_error / tr.kinetic_initial;
            println!(
                "T1 c = {c}, seed {seed}: outcome {:?} at t = {:.2}, {} steps, max |ΔW|/T0 = {rel:.2e}",
                tr.outcome, tr.end.t, tr.stats.n_accept
            );
            assert!(tr.stats.n_accept >= 20, "flight too short to be meaningful");
            assert!(rel < 1e-10, "energy error {rel:.3e}");
        }
    }
}

// --- T2/T3: Coulomb scattering ------------------------------------------------------------

/// Scatters a unit particle off charge `kappa` at the origin from 1e6 units away and
/// returns (measured deflection, analytic deflection).
fn scatter(kappa: f64, c: f64, b: f64, p: f64) -> (f64, f64) {
    let d = 1e6;
    let mut scn = central_charge(
        kappa,
        1e-3,
        c,
        DVec3::new(-d, b, 0.0),
        DVec3::new(p, 0.0, 0.0),
        1e9,
    );
    scn.bounds = Some(cube(2.0 * d));
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    assert_eq!(tr.outcome, Outcome::LeftBounds);
    let measured = angle_between(scn.p0, tr.end.p);
    let analytic = CoulombOrbit::from_state(kappa, 1.0, c, scn.x0, scn.p0).deflection();
    // The velocity direction at distance r differs from the asymptote by O(b κ/(m v² r²))
    // ≈ 1e-12 here, far below the tolerance.
    (measured, analytic)
}

#[test]
fn t2_rutherford_scattering() {
    for (b, p) in [(1.0, 1.0), (0.2, 1.0), (5.0, 0.7), (1.0, 3.0)] {
        let (measured, analytic) = scatter(1.0, f64::INFINITY, b, p);
        // Textbook form: tan(θ/2) = κ/(m v∞² b∞), with v∞ and b∞ from E and L.
        let e = p * p / 2.0 + 1.0 / (1e12 + b * b).sqrt();
        let v_inf = (2.0 * e).sqrt();
        let b_inf = b * p / v_inf;
        let rutherford = 2.0 * (1.0 / (v_inf * v_inf * b_inf)).atan();
        let rel = (measured - analytic).abs() / analytic;
        println!(
            "T2 b = {b}, p = {p}: θ = {measured:.12} rad, Rutherford {rutherford:.12}, relative error {rel:.2e}"
        );
        assert!((analytic - rutherford).abs() < 1e-14);
        assert!(rel < 1e-8, "relative error {rel:.3e}");
    }
}

#[test]
fn t3_relativistic_coulomb_scattering() {
    for (kappa, c, b, p) in [
        (1.0, 2.0, 1.0, 1.0),
        (1.0, 1.0, 1.0, 5.0),
        (-1.0, 2.0, 1.5, 1.0),
        (-0.5, 1.0, 1.0, 3.0),
    ] {
        let (measured, analytic) = scatter(kappa, c, b, p);
        let newtonian = CoulombOrbit::from_state(
            kappa,
            1.0,
            f64::INFINITY,
            DVec3::new(-1e6, b, 0.0),
            DVec3::new(p, 0.0, 0.0),
        )
        .deflection();
        let v_over_c = Kinematics::new(1.0, c)
            .velocity(DVec3::new(p, 0.0, 0.0))
            .length()
            / c;
        let rel = (measured - analytic).abs() / analytic;
        println!(
            "T3 κ = {kappa}, c = {c}, b = {b}, p = {p} (v/c = {v_over_c:.3}): θ = {measured:.12}, \
             analytic {analytic:.12} (Newtonian would be {newtonian:.6}), relative error {rel:.2e}"
        );
        assert!(rel < 1e-8, "relative error {rel:.3e}");
    }
}

// --- T4: Kepler orbit ---------------------------------------------------------------------

/// Unit particle at periapsis r = 0.5 of an attractive unit charge, with the Newtonian
/// speed of an e = 0.5, a = 1 orbit (period 2π when c = ∞).
fn kepler_scenario(c: f64, t_max: f64) -> Scenario<Coulomb> {
    central_charge(
        -1.0,
        1e-3,
        c,
        DVec3::new(0.5, 0.0, 0.0),
        DVec3::new(0.0, 3f64.sqrt(), 0.0),
        t_max,
    )
}

#[test]
fn t4_kepler_orbit() {
    let orbits = 10.0;
    let scn = kepler_scenario(f64::INFINITY, orbits * 2.0 * PI);
    let l0 = scn.x0.cross(scn.p0);
    let mut peri = Vec::new();
    let mut l_err: f64 = 0.0;
    let tr = run_observed(&scn, &RunSettings::with_tolerance(TOL), |step| {
        peri.extend(periapsis_times(step));
        let (x, p) = step.state(step.t_end());
        l_err = l_err.max((x.cross(p) - l0).length() / l0.length());
    });
    assert_eq!(tr.outcome, Outcome::Timeout);

    let period_exact = 2.0 * PI; // 2π sqrt(m a³ / |κ|), a = 1
    let period_err = peri
        .iter()
        .enumerate()
        .map(|(k, &t)| {
            (t - (k as f64 + 1.0) * period_exact).abs() / ((k as f64 + 1.0) * period_exact)
        })
        .fold(0.0, f64::max);
    let closure = (tr.end.x - scn.x0).length();
    let rel_e = tr.energy_max_abs_error / tr.energy_initial.abs();
    println!(
        "T4 Kepler, {orbits} orbits: {} periapsides, max relative period error {period_err:.2e}, \
         closure after {orbits} orbits {closure:.2e}, max |ΔL|/L {l_err:.2e}, max |ΔE|/|E| {rel_e:.2e}",
        peri.len()
    );
    // The 10th passage coincides with t_max = 10 T and may fall on either side of it.
    assert!((9..=10).contains(&peri.len()));
    assert!(period_err < 1e-9, "period error {period_err:.3e}");
    assert!(closure < 1e-8, "closure {closure:.3e}");
    assert!(l_err < 1e-10, "angular momentum error {l_err:.3e}");
}

// --- T5: relativistic perihelion precession -----------------------------------------------

#[test]
fn t5_relativistic_precession() {
    for c in [10.0, 3.0] {
        let revolutions = 20;
        let scn = kepler_scenario(c, (revolutions as f64 + 2.0) * 2.0 * PI * 1.5);
        let orbit = CoulombOrbit::from_state(-1.0, 1.0, c, scn.x0, scn.p0);
        let mut angles = Vec::new();
        let mut shape_err: f64 = 0.0;
        let mut phi_unwrapped = 0.0;
        let mut phi_prev = 0.0;
        run_observed(&scn, &RunSettings::with_tolerance(TOL), |step| {
            for t in periapsis_times(step) {
                let (x, _) = step.state(t);
                angles.push(x.y.atan2(x.x));
            }
            // Orbit shape r(φ) along the trajectory, with φ unwrapped from the start.
            let (x, _) = step.state(step.t_end());
            let phi = x.y.atan2(x.x);
            let mut d = phi - phi_prev;
            if d < -PI {
                d += 2.0 * PI;
            } else if d > PI {
                d -= 2.0 * PI;
            }
            phi_unwrapped += d;
            phi_prev = phi;
            shape_err = shape_err.max((x.length() - orbit.r(phi_unwrapped)).abs() / x.length());
        });
        assert!(angles.len() >= revolutions);
        // Successive periapsis angles differ by Δφ (mod 2π) with |Δφ| < π: unwrap the
        // differences and average over the revolutions.
        let mut total = 0.0;
        let mut prev = 0.0; // periapsis at t = 0 lies on the +x axis
        for &a in &angles[..revolutions] {
            let mut d = a - prev;
            while d <= -PI {
                d += 2.0 * PI;
            }
            while d > PI {
                d -= 2.0 * PI;
            }
            total += d;
            prev = a;
        }
        let measured = total / revolutions as f64;
        let analytic = 2.0 * PI * (1.0 / orbit.gamma - 1.0);
        let rel = (measured - analytic).abs() / analytic;
        println!(
            "T5 c = {c}: precession per revolution {measured:.12} rad, analytic {analytic:.12}, \
             relative error {rel:.2e}; max relative deviation from r(φ) {shape_err:.2e}"
        );
        assert!(rel < 1e-6, "relative error {rel:.3e}");
        assert!(shape_err < 1e-9, "shape error {shape_err:.3e}");
    }
}

// --- T6: hyperbolic motion in a uniform field ---------------------------------------------

#[test]
fn t6_hyperbolic_motion_uniform_field() {
    let (m, c, q_e) = (1.0, 1.0, 1.0);
    for p_perp in [0.0, 0.5, 3.0] {
        for t_end in [1.0, 10.0, 1000.0] {
            let scn = Scenario {
                field: UniformElectric {
                    e: DVec3::new(q_e, 0.0, 0.0),
                },
                obstacles: vec![],
                particle: UNIT_PARTICLE,
                c,
                x0: DVec3::ZERO,
                p0: DVec3::new(0.0, p_perp, 0.0),
                detector: None,
                bounds: None,
                t_max: t_end,
                radiation_reaction: false,
                acceptance: None,
            };
            let tr = run(&scn, &RunSettings::with_tolerance(TOL));
            let eps0 = (m * m * c.powi(4) + p_perp * p_perp * c * c).sqrt();
            let x_exact = ((eps0 * eps0 + (q_e * c * t_end).powi(2)).sqrt() - eps0) / q_e;
            let y_exact = p_perp * c / q_e * (q_e * c * t_end / eps0).asinh();
            let rel_x = (tr.end.x.x - x_exact).abs() / x_exact;
            let rel_y = if p_perp > 0.0 {
                (tr.end.x.y - y_exact).abs() / y_exact
            } else {
                tr.end.x.y.abs()
            };
            let rel_p = (tr.end.p - DVec3::new(q_e * t_end, p_perp, 0.0)).length() / (q_e * t_end);
            println!(
                "T6 p⊥ = {p_perp}, t = {t_end}: relative error x {rel_x:.2e}, y {rel_y:.2e}, p {rel_p:.2e}"
            );
            assert!(rel_x < 1e-11 && rel_y < 1e-11 && rel_p < 1e-12);
        }
    }
}

// --- T7 and T8 are property tests in `properties.rs`. -------------------------------------

// --- T9: grazing events -------------------------------------------------------------------

#[test]
fn t9_grazing_straight_line() {
    // No field: the particle moves on the line y = d; closest approach to the origin is d.
    let d = 1.0;
    for delta in [1e-6, 1e-9, 1e-12, 1e-14] {
        for (radius, should_hit) in [(d * (1.0 + delta), true), (d * (1.0 - delta), false)] {
            let scn = Scenario {
                field: Coulomb::new(&[]),
                obstacles: vec![Shape::Sphere(Sphere {
                    center: DVec3::ZERO,
                    radius,
                })],
                particle: UNIT_PARTICLE,
                c: f64::INFINITY,
                x0: DVec3::new(-10.0, d, 0.0),
                p0: DVec3::new(1.3, 0.0, 0.0),
                detector: None,
                bounds: Some(cube(20.0)),
                t_max: 100.0,
                radiation_reaction: false,
                acceptance: None,
            };
            let tr = run(&scn, &RunSettings::with_tolerance(1e-8));
            let hit = tr.outcome == Outcome::Collided(0);
            println!(
                "T9 straight line, δ = {delta:e}, hit expected {should_hit}: outcome {:?}",
                tr.outcome
            );
            assert_eq!(hit, should_hit, "δ = {delta:e}");
        }
    }
}

#[test]
fn t9_grazing_coulomb_orbit() {
    for (kappa, c) in [(1.0, f64::INFINITY), (1.0, 3.0), (-1.0, 3.0)] {
        let x0 = DVec3::new(-20.0, 1.0, 0.0);
        let p0 = DVec3::new(1.0, 0.0, 0.0);
        let r_min = CoulombOrbit::from_state(kappa, 1.0, c, x0, p0).r_min();

        // Measured closest approach with a negligible sphere.
        let mut scn = central_charge(kappa, 1e-6, c, x0, p0, 1e4);
        scn.bounds = Some(cube(40.0));
        let mut measured = f64::INFINITY;
        run_observed(&scn, &RunSettings::with_tolerance(TOL), |step| {
            for t in periapsis_times(step) {
                measured = measured.min(step.state(t).0.length());
            }
        });
        let rel = (measured - r_min).abs() / r_min;
        println!(
            "T9 Coulomb κ = {kappa}, c = {c}: r_min {r_min:.14}, measured {measured:.14}, relative error {rel:.2e}"
        );
        assert!(rel < 1e-10, "closest approach error {rel:.3e}");

        for delta in [1e-6, 1e-8] {
            for (radius, should_hit) in [
                (r_min * (1.0 + delta), true),
                (r_min * (1.0 - delta), false),
            ] {
                let mut scn = central_charge(kappa, radius, c, x0, p0, 1e4);
                scn.bounds = Some(cube(40.0));
                let tr = run(&scn, &RunSettings::with_tolerance(TOL));
                let hit = tr.outcome == Outcome::Collided(0);
                assert_eq!(
                    hit, should_hit,
                    "κ = {kappa}, c = {c}, δ = {delta:e}: {:?}",
                    tr.outcome
                );
            }
        }
    }
}

// --- T10: convergence with tolerance ------------------------------------------------------

/// Closure error after one Newtonian Kepler orbit (a = 1, period 2π) of eccentricity `e`.
fn kepler_closure_error(e: f64, tol: f64) -> f64 {
    let x0 = DVec3::new(1.0 - e, 0.0, 0.0);
    let p0 = DVec3::new(0.0, ((1.0 + e) / (1.0 - e)).sqrt(), 0.0);
    let scn = central_charge(-1.0, 1e-3, f64::INFINITY, x0, p0, 2.0 * PI);
    let tr = run(&scn, &RunSettings::with_tolerance(tol));
    (tr.end.x - x0).length()
}

#[test]
fn t10_error_scales_with_tolerance() {
    // Asymptotically DOP853's error estimate scales like the local error of the 8th-order
    // solution (h^9), so the global error should become proportional to tol: slope 1.
    // Measured in the asymptotic regime (≥ 20 steps per orbit).
    let points: Vec<(f64, f64)> = (9..=14)
        .map(|k| {
            let tol = 10f64.powi(-k);
            let err = kepler_closure_error(0.2, tol);
            println!(
                "T10 e = 0.2, tol {tol:.0e}: closure error {err:.3e} ({:.1} tol)",
                err / tol
            );
            (tol.log10(), err.log10())
        })
        .collect();
    let n = points.len() as f64;
    let (sx, sy) = points
        .iter()
        .fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
    let (mx, my) = (sx / n, sy / n);
    let num: f64 = points.iter().map(|(x, y)| (x - mx) * (y - my)).sum();
    let den: f64 = points.iter().map(|(x, _)| (x - mx).powi(2)).sum();
    let slope = num / den;
    println!("T10 slope d log(error) / d log(tol) = {slope:.3}");
    assert!((0.9..1.1).contains(&slope), "slope {slope:.3}");

    // Calibration guard: global error after one orbit stays within 100 tol, including the
    // pre-asymptotic regime of the more eccentric orbit.
    for e in [0.2, 0.5] {
        for k in 6..=13 {
            let tol = 10f64.powi(-k);
            let err = kepler_closure_error(e, tol);
            assert!(err < 100.0 * tol, "e = {e}, tol {tol:e}: error {err:.3e}");
        }
    }
}
