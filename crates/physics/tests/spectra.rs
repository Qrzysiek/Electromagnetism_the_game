//! Validation of the far-zone radiation measure (`spectrum.rs`, PHYSICS.md §3.4) against
//! Jackson's analytic results, on prescribed motions.

use std::f64::consts::PI;

use physics::DVec3;
use physics::spectrum::{Emission, band_energy, lienard_energy, spectrum};

/// Bessel function `J_m(x)` by its power series (moderate `x`).
fn bessel_j(m: i32, x: f64) -> f64 {
    if m < 0 {
        return if m % 2 == 0 { 1.0 } else { -1.0 } * bessel_j(-m, x);
    }
    let mut term = (x / 2.0).powi(m) / (1..=m).map(f64::from).product::<f64>();
    let mut sum = term;
    for k in 1..200 {
        term *= -(x * x / 4.0) / (f64::from(k) * f64::from(k + m));
        sum += term;
        if term.abs() < 1e-18 * sum.abs() {
            break;
        }
    }
    sum
}

/// Uniform circular motion (radius `rho`, angular frequency `w0`) sampled `per_turn`
/// times per turn over `turns` turns.
fn circle(rho: f64, w0: f64, turns: u32, per_turn: u32) -> Vec<Emission> {
    let n = turns * per_turn;
    (0..=n)
        .map(|i| {
            let t = f64::from(i) / f64::from(per_turn) * 2.0 * PI / w0;
            let (s, c) = (w0 * t).sin_cos();
            (
                t,
                DVec3::new(rho * c, rho * s, 0.0),
                DVec3::new(-rho * w0 * s, rho * w0 * c, 0.0),
                DVec3::new(-rho * w0 * w0 * c, -rho * w0 * w0 * s, 0.0),
            )
        })
        .collect()
}

/// S1, Jackson Pr. 14.15: a charge on a circle radiates into its orbital plane (θ = π/2)
/// at the harmonics `mω₀`, with power per steradian
/// `dP_m/dΩ = (q² m² ω₀² β²/2πc) J'_m(mβ)²`. Over `T` the energy in a band ±ω₀/2 around
/// each harmonic must be `T dP_m/dΩ`, up to the finite flight's leakage: each line is a
/// sinc² of width 2π/T whose tails outside the band carry a fraction ~2/(π² (ω₀/2) T)
/// (and they are shared with the neighbours), so the error falls as 1/T. β = 0.5, 40 and
/// 160 turns: at most 3.4e-3 and 9.1e-4.
#[test]
fn s1_harmonics_of_circular_motion() {
    let (q, c) = (1.0, 10.0);
    let (rho, w0) = (1.0, 5.0);
    let beta = rho * w0 / c;
    let n = DVec3::X;
    for (turns, bound) in [(40, 5e-3), (160, 1.2e-3)] {
        let samples = circle(rho, w0, turns, 400);
        let t_total = f64::from(turns) * 2.0 * PI / w0;
        let mut worst: f64 = 0.0;
        for m in 1..=4 {
            let mf = f64::from(m);
            let jp = 0.5 * (bessel_j(m - 1, mf * beta) - bessel_j(m + 1, mf * beta));
            let jackson =
                q * q * mf * mf * w0 * w0 * beta * beta / (2.0 * PI * c) * jp * jp * t_total;
            let e = band_energy(&samples, q, c, n, (mf - 0.5) * w0, (mf + 0.5) * w0, false);
            let err = e / jackson - 1.0;
            println!(
                "S1 {turns} turns, harmonic {m}: band energy {e:.6e}, Jackson {jackson:.6e}, rel. error {err:.1e}"
            );
            worst = worst.max(err.abs());
        }
        assert!(worst < bound);
    }
}

/// S2, Parseval (Jackson 14.60–14.65): the spectrum integrated over all frequencies equals
/// Liénard's energy per steradian, for a relativistic transient: a charge at β = 0.6
/// kicked sideways by a Gaussian pulse of acceleration, seen from several directions.
#[test]
fn s2_spectrum_integrates_to_the_lienard_energy() {
    let (q, c) = (1.0, 1.0);
    let v0 = 0.6;
    let (a0, sigma) = (0.05, 2.0);
    let dt = 0.002;
    // Velocity: v0 along x, v_y = ∫ a dt (the kick, ≪ c).
    let n_steps = 20_000;
    let mut samples: Vec<Emission> = Vec::with_capacity(n_steps + 1);
    let (mut x, mut vy) = (DVec3::new(-20.0 * v0, 0.0, 0.0), 0.0);
    for i in 0..=n_steps {
        #[allow(clippy::cast_precision_loss)]
        let t = -20.0 + i as f64 * dt;
        let ay = a0 * (-(t * t) / (2.0 * sigma * sigma)).exp();
        let v = DVec3::new(v0, vy, 0.0);
        samples.push((t, x, v, DVec3::new(0.0, ay, 0.0)));
        // Exact enough for the test: the same samples feed both sides.
        x += v * dt;
        vy += ay * dt;
    }
    let mut worst: f64 = 0.0;
    for deg in [0.0_f64, 20.0, 60.0, 135.0] {
        let n = DVec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
        let lienard = lienard_energy(&samples, q, c, n);
        // The pulse lasts ~σ(1 − n·β): its spectrum ends near a few /(σκ); 40/(σκ) is far out.
        let kappa = 1.0 - n.x * v0;
        let top = 40.0 / (sigma * kappa);
        let total = band_energy(&samples, q, c, n, 0.0, top, false);
        let err = total / lienard - 1.0;
        println!(
            "S2 at {deg}°: ∫ d²I/dωdΩ dω = {total:.8e}, Liénard {lienard:.8e}, rel. error {err:.1e}"
        );
        worst = worst.max(err.abs());
    }
    // A spot check that the spectrum is flat at low frequency (a net velocity change
    // radiates like a sudden kick for ω ≪ 1/σκ).
    let s = spectrum(&samples, q, c, DVec3::Y, 1e-3, 1e-3, 2, false);
    assert!((s[0] / s[1] - 1.0).abs() < 1e-3);
    assert!(worst < 1e-4);
}

/// S3, the measure along a computed flight (`trajectory::run` with a radiation goal): a
/// charge at β = 0.5 on a circle in a uniform field (radius 1, counter-clockwise) enters a
/// detector after about three quarters of a turn. Its measured energy per steradian into
/// an arc (all frequencies, a low band and a high band) must match the same measure on the
/// exact circular motion up to the detector entry (20 000 samples per turn and one at the
/// entry time).
#[test]
fn s3_measure_along_a_computed_flight() {
    use physics::dynamics::Particle;
    use physics::field::UniformFields;
    use physics::geometry::{Aabb, Region};
    use physics::spectrum::RadiationWindow;
    use physics::trajectory::{Acceptance, Outcome, RunSettings, Scenario, run};

    let (q, m, c): (f64, f64, f64) = (1.0, 1.0, 10.0);
    let v: f64 = 5.0;
    let gamma = 1.0 / (1.0 - (v / c) * (v / c)).sqrt();
    let w0 = v; // radius 1
    for band in [None, Some((2.0, 12.0)), Some((20.0, 60.0))] {
        let window = RadiationWindow {
            axis: DVec3::new(1.0, 1.0, 0.0).normalize(),
            half_angle: 15f64.to_radians(),
            band,
            energy: (0.0, 1e6),
            abrupt_stop: false,
        };
        let scn = Scenario {
            field: UniformFields {
                e: DVec3::ZERO,
                b: DVec3::new(0.0, 0.0, -gamma * m * v),
            },
            obstacles: vec![],
            particle: Particle {
                charge: q,
                mass: m,
                radius: 0.0,
                moment: 0.0,
            },
            c,
            x0: DVec3::new(1.0, 0.0, 0.0),
            p0: DVec3::new(0.0, gamma * m * v, 0.0),
            detector: Some(Region::Box(Aabb {
                min: DVec3::new(-0.1, -1.2, -1.0),
                max: DVec3::new(0.1, -0.8, 1.0),
            })),
            bounds: None,
            t_max: 10.0,
            radiation_reaction: false,
            acceptance: Some(Acceptance {
                radiation: Some(window),
                ..Acceptance::default()
            }),
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
        assert_eq!(tr.outcome, Outcome::Arrived);
        let measured = tr.radiation.expect("measured on arrival");
        // Exact motion up to the entry angle (x = −0.1 on the lower half).
        let theta_end = 2.0 * PI - (-0.1f64).acos();
        let t_end = theta_end / w0;
        // Ending exactly at the entry: the flight stops with the acceleration on, and that
        // edge feeds the high frequencies (a reference ending one sample early was off by
        // 2e-4 in the band 20–60).
        let mut exact: Vec<Emission> = circle(1.0, w0, 1, 20_000)
            .into_iter()
            .filter(|s| s.0 < t_end)
            .collect();
        let (s, co) = (w0 * t_end).sin_cos();
        exact.push((
            t_end,
            DVec3::new(co, s, 0.0),
            DVec3::new(-w0 * s, w0 * co, 0.0),
            DVec3::new(-w0 * w0 * co, -w0 * w0 * s, 0.0),
        ));
        let reference = window.measure(&exact, q, c);
        let err = measured / reference - 1.0;
        println!(
            "S3 band {band:?}: measured {measured:.8e}, exact motion {reference:.8e}, rel. error {err:.1e} (entry at t = {:.6}, exact {t_end:.6})",
            tr.end.t
        );
        // The band 20–60 (harmonics 4–12): pieces are set by the phase step there; its
        // error, 0.8–2.1e-4 for phase steps of 0.25–1 rad (not converging with the step: the
        // abrupt end at the entry, seen in the high frequencies), is held to the 1e-3 the
        // radiation goals need (their windows have factors of several to spare).
        let bound = if band.is_some_and(|(lo, _)| lo >= 20.0) {
            1e-3
        } else {
            1e-4
        };
        assert!(err.abs() < bound);
    }
}

/// S4, the sudden stop (Jackson §15.2; the `abrupt_stop` option). (a) A charge at β = 0.8
/// moving uniformly and stopped instantly at t = 0 radiates the flat spectrum
/// `(q²/4π²c) β² sin²θ / (1 − β cos θ)²`; (b) the same charge stopped smoothly within
/// τ = 0.002 (velocity `v₀ / (1 + e^{t/τ})`), integrated by the acceleration form alone,
/// must converge to it at frequencies far below the inverse stopping time: quadratically
/// in ωτ (the difference falls 16 ± 4 times when ω falls 4 times), and below 1e-5 at
/// ωτ = 1e-3. Measured: (a) ≤ 3e-16; (b) 7e-7 … 4.6e-6 at ωτ = 1e-3, 1.8e-4 … 1.2e-3 at
/// ωτ = 1.6e-2 (larger away from the direction of motion), ratios 16.0.
#[test]
fn s4_sudden_stop() {
    let (q, c): (f64, f64) = (1.0, 1.0);
    let beta: f64 = 0.8;
    let v0 = beta * c;
    let tau = 0.002;
    // (a) uniform motion to the stop.
    let uniform: Vec<Emission> = (0..=200)
        .map(|i| {
            let t = -2.0 + 0.01 * f64::from(i);
            (
                t,
                DVec3::new(v0 * t, 0.0, 0.0),
                DVec3::new(v0, 0.0, 0.0),
                DVec3::ZERO,
            )
        })
        .collect();
    // (b) the smooth stop: s = 1/(1 + e^{t/τ}), x = v₀ (t − τ ln(1 + e^{t/τ})).
    let n_steps = 60_000;
    let smooth: Vec<Emission> = (0..=n_steps)
        .map(|i| {
            let t = -2.0 + 2.1 * f64::from(i) / f64::from(n_steps);
            let e = (t / tau).exp();
            let s = 1.0 / (1.0 + e);
            let x = v0 * (t - tau * e.ln_1p());
            let a = -v0 * s * (1.0 - s) / tau;
            (
                t,
                DVec3::new(x, 0.0, 0.0),
                DVec3::new(v0 * s, 0.0, 0.0),
                DVec3::new(a, 0.0, 0.0),
            )
        })
        .collect();
    let mut worst_a: f64 = 0.0;
    for deg in [10.0_f64, 36.87, 60.0, 120.0] {
        let th = deg.to_radians();
        let n = DVec3::new(th.cos(), th.sin(), 0.0);
        let jackson = q * q / (4.0 * PI * PI * c) * beta * beta * th.sin().powi(2)
            / (1.0 - beta * th.cos()).powi(2);
        let mut diffs = Vec::new();
        for omega in [0.5, 2.0, 8.0] {
            let a = spectrum(&uniform, q, c, n, omega, 1.0, 1, true)[0];
            let b = spectrum(&smooth, q, c, n, omega, 1.0, 1, false)[0];
            let (ea, eb) = (a / jackson - 1.0, b / jackson - 1.0);
            println!(
                "S4 at {deg}°, ω = {omega}: abrupt {ea:.1e}, smooth stop (ωτ = {:.0e}) {eb:.1e}",
                omega * tau
            );
            worst_a = worst_a.max(ea.abs());
            diffs.push(eb.abs());
        }
        let (r1, r2) = (diffs[1] / diffs[0], diffs[2] / diffs[1]);
        println!("S4 at {deg}°: ratios {r1:.1}, {r2:.1} (quadratic: 16)");
        assert!(diffs[0] < 1e-5 && (12.0..20.0).contains(&r1) && (12.0..20.0).contains(&r2));
    }
    assert!(worst_a < 1e-12);
}
