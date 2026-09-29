//! Validation of the far-zone radiation measure (`spectrum.rs`, PHYSICS.md §3.4) against
//! Jackson's analytic results, on prescribed motions.

use std::f64::consts::PI;

use physics::DVec3;
use physics::spectrum::{
    Emission, band_energy, lienard_energy, spectrum, system_band_energy, system_lienard_energy,
    system_spectrum,
};

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
    circle_from(rho, w0, 0.0, turns, per_turn)
}

/// `circle`, with the charge ahead by the phase `theta` at t = 0.
fn circle_from(rho: f64, w0: f64, theta: f64, turns: u32, per_turn: u32) -> Vec<Emission> {
    let n = turns * per_turn;
    (0..=n)
        .map(|i| {
            let t = f64::from(i) / f64::from(per_turn) * 2.0 * PI / w0;
            let (s, c) = (w0 * t + theta).sin_cos();
            (
                t,
                DVec3::new(rho * c, rho * s, 0.0),
                DVec3::new(-rho * w0 * s, rho * w0 * c, 0.0),
                DVec3::new(-rho * w0 * w0 * c, -rho * w0 * w0 * s, 0.0),
            )
        })
        .collect()
}

/// S1, Jackson Pr. 14.15: a charge on a circle (β = 0.5) radiates into its orbital plane
/// at the harmonics `mω₀`, with power per steradian `(q² m² ω₀² β²/2πc) J'_m(mβ)²`. For a
/// flight of exactly N turns the energy in a band ±ω₀/2 around each harmonic is known
/// exactly: the integrand of Jackson 14.65 is periodic but for `e^{iωt}`, so the amplitude
/// is a sum over harmonics with finite-time factors (`scripts/wolfram/s1_circular_harmonics.wls`,
/// Wolfram Engine 14.2: Fourier coefficients by an exponentially accurate DFT, the band by
/// quadrature). Required: the measure to 1e-5 of the exact value (set before measuring:
/// 400 samples per turn), and the exact values within the finite flight's leakage of
/// Jackson's infinite-time power (a line is a sinc² of width 2π/T, whose tails outside
/// the band carry ~2/(π² (ω₀/2) T): 5e-3 at 40 turns, 1.2e-3 at 160).
#[test]
fn s1_harmonics_of_circular_motion() {
    // (turns, harmonic, exact band energy, infinite-time Jackson), from the Wolfram script.
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const EXACT: [(u32, i32, f64, f64); 8] = [
        (40, 1, 1.027697010635241, 1.0302753517029275),
        (40, 2, 0.8810862284884798, 0.8840475603754644),
        (40, 3, 0.5460370542835826, 0.5460828703215562),
        (40, 4, 0.29820585780936115, 0.2972098500096013),
        (160, 1, 4.118522708010444, 4.12110140681171),
        (160, 2, 3.5332287356623806, 3.5361902415018576),
        (160, 3, 2.184285704684461, 2.184331481286225),
        (160, 4, 1.1898354770073927, 1.1888394000384053),
    ];
    let (q, c) = (1.0, 10.0);
    let (rho, w0) = (1.0, 5.0);
    let beta = rho * w0 / c;
    let n = DVec3::X;
    let (mut worst, mut worst_jackson) = (0.0_f64, [0.0_f64; 2]);
    for (turns, m, exact, jackson_ref) in EXACT {
        let samples = circle(rho, w0, turns, 400);
        let mf = f64::from(m);
        // The infinite-time power, computed here too (checks the pasted column).
        let t_total = f64::from(turns) * 2.0 * PI / w0;
        let jp = 0.5 * (bessel_j(m - 1, mf * beta) - bessel_j(m + 1, mf * beta));
        let jackson = q * q * mf * mf * w0 * w0 * beta * beta / (2.0 * PI * c) * jp * jp * t_total;
        assert!((jackson / jackson_ref - 1.0).abs() < 1e-12);
        let e = band_energy(&samples, q, c, n, (mf - 0.5) * w0, (mf + 0.5) * w0, false);
        let (err, leak) = (e / exact - 1.0, exact / jackson - 1.0);
        println!(
            "S1 {turns} turns, harmonic {m}: measure {e:.12e}, exact {exact:.12e}, rel. error {err:.1e}; exact vs infinite-time {leak:.1e}"
        );
        worst = worst.max(err.abs());
        let slot = usize::from(turns == 160);
        worst_jackson[slot] = worst_jackson[slot].max(leak.abs());
    }
    println!(
        "S1: worst {worst:.1e}; leakage {:.1e} (40 turns), {:.1e} (160)",
        worst_jackson[0], worst_jackson[1]
    );
    assert!(worst < 1e-5 && worst_jackson[0] < 5e-3 && worst_jackson[1] < 1.2e-3);
}

/// S2, Parseval (Jackson 14.60–14.65): the spectrum integrated over all frequencies equals
/// Liénard's energy per steradian, for a relativistic transient: a charge at β = 0.6
/// kicked sideways by a pulse of acceleration `a₀ sech²(t/σ)`, seen from several
/// directions. The world line is exact (`v_y = a₀σ(1 + tanh(t/σ))`,
/// `y = a₀σ(t + σ ln cosh(t/σ))`): the first version stepped the position by Euler, which
/// is inconsistent with the velocities off the axis and limited the agreement to 6e-5.
/// Required: 1e-8 (set before measuring).
#[test]
fn s2_spectrum_integrates_to_the_lienard_energy() {
    let (q, c) = (1.0, 1.0);
    let v0 = 0.6;
    let (a0, sigma) = (0.05, 2.0);
    let n_steps = 20_000;
    let samples: Vec<Emission> = (0..=n_steps)
        .map(|i| {
            let t = -20.0 + 40.0 * f64::from(i) / f64::from(n_steps);
            let (th, ch) = ((t / sigma).tanh(), (t / sigma).cosh());
            let x = DVec3::new(v0 * t, a0 * sigma * (t + sigma * ch.ln()), 0.0);
            let v = DVec3::new(v0, a0 * sigma * (1.0 + th), 0.0);
            let a = DVec3::new(0.0, a0 * (1.0 - th * th), 0.0);
            (t, x, v, a)
        })
        .collect();
    let mut worst: f64 = 0.0;
    for deg in [0.0_f64, 20.0, 60.0, 135.0] {
        let n = DVec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
        let lienard = lienard_energy(&samples, q, c, n);
        // The pulse lasts ~σκ as received: its spectrum falls as e^{−πωσκ/2}; 40/(σκ) is
        // far out (κ ≥ 0.4 for these directions).
        let kappa = 1.0 - n.x * v0;
        let top = 40.0 / (sigma * kappa);
        let total = band_energy(&samples, q, c, n, 0.0, top, false);
        let err = total / lienard - 1.0;
        println!(
            "S2 at {deg}°: ∫ d²I/dωdΩ dω = {total:.12e}, Liénard {lienard:.12e}, rel. error {err:.1e}"
        );
        worst = worst.max(err.abs());
    }
    // A spot check that the spectrum is flat at low frequency (a net velocity change
    // radiates like a sudden kick for ω ≪ 1/σκ).
    let s = spectrum(&samples, q, c, DVec3::Y, 1e-3, 1e-3, 2, false);
    assert!((s[0] / s[1] - 1.0).abs() < 1e-3);
    assert!(worst < 1e-8);
}

/// S5, Jackson Pr. 14.23: charges `qⱼ` on one circle at fixed phases `θⱼ` radiate into
/// the harmonic `mω₀` as a single unit charge times the form factor `|Σⱼ qⱼ e^{imθⱼ}|²`:
/// N equal, equally spaced charges radiate only at multiples of `Nω₀`, with N² times one
/// charge's intensity. Checked at the harmonics of a flight of whole turns (where each
/// charge's amplitude is the first one's times `e^{−imθⱼ}` exactly, by periodicity) in two
/// in-plane directions. With the phases on the sampling grid (360 samples per turn) the
/// charges' samples are shifted copies and the identity holds to rounding: required 1e-10
/// of the single charge's value (set before measuring). Off the grid (θ = 0.3, 1.1 rad)
/// the sampling differs between the charges: required 1e-6 (S1's samples are accurate to
/// 1e-5 of the exact value at 400 per turn; the difference of two samplings is smaller).
#[test]
fn s5_ring_of_charges_form_factor() {
    let c = 10.0;
    let (rho, w0) = (1.0, 5.0);
    let (turns, per_turn) = (20, 360);
    let one = circle(rho, w0, turns, per_turn);
    let mut worst: [f64; 2] = [0.0; 2];
    // (name, charges and phases, 0 on the sampling grid / 1 off it)
    type Ring = Vec<(f64, f64)>;
    let cases: [(&str, Ring, usize); 5] = [
        ("N = 2", vec![(1.0, 0.0), (1.0, PI)], 0),
        (
            "N = 3",
            (0..3)
                .map(|j| (1.0, 2.0 * PI * f64::from(j) / 3.0))
                .collect(),
            0,
        ),
        (
            "N = 4",
            (0..4)
                .map(|j| (1.0, 2.0 * PI * f64::from(j) / 4.0))
                .collect(),
            0,
        ),
        (
            "unequal, on the grid",
            vec![
                (1.0, 0.0),
                (-0.5, 40f64.to_radians()),
                (2.0, 200f64.to_radians()),
            ],
            0,
        ),
        (
            "unequal, off the grid",
            vec![(1.0, 0.0), (-0.5, 0.3), (2.0, 1.1)],
            1,
        ),
    ];
    for (name, charges, grid) in &cases {
        let flights: Vec<Vec<Emission>> = charges
            .iter()
            .map(|&(_, th)| circle_from(rho, w0, th, turns, per_turn))
            .collect();
        let sources: Vec<(f64, &[Emission])> = charges
            .iter()
            .zip(&flights)
            .map(|(&(q, _), f)| (q, f.as_slice()))
            .collect();
        for deg in [0.0_f64, 35.0] {
            let n = DVec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
            let sys = system_spectrum(&sources, c, n, w0, w0, 6);
            let single = spectrum(&one, 1.0, c, n, w0, w0, 6, false);
            for m in 1..=6 {
                let (re, im) = charges.iter().fold((0.0, 0.0), |(re, im), &(q, th)| {
                    let ph = f64::from(m) * th;
                    (re + q * ph.cos(), im + q * ph.sin())
                });
                let factor = re * re + im * im;
                let k = usize::try_from(m - 1).unwrap();
                let err = (sys[k] - factor * single[k]).abs() / single[k];
                println!(
                    "S5 {name} at {deg}°, m = {m}: system {:.10e}, |F|² × single {:.10e} \
                     (|F|² = {factor:.6}), difference {err:.1e} of one charge",
                    sys[k],
                    factor * single[k]
                );
                worst[*grid] = worst[*grid].max(err);
            }
        }
    }
    assert!(worst[0] < 1e-10, "on the grid: {:.2e}", worst[0]);
    assert!(worst[1] < 1e-6, "off the grid: {:.2e}", worst[1]);
}

/// S6, Parseval for a system (Jackson 14.60–14.65 with the fields added): the system's
/// spectrum integrated over all frequencies equals its energy per steradian from the
/// time domain, `(1/4πc) ∫ |Σⱼ qⱼ gⱼ|² dτ`, for two charges kicked by pulses at different
/// times and places (their radiation overlapping in the phase time in some directions,
/// not in others), charges +1 and −0.7. Required 1e-8, as S2 (set before measuring). And
/// the time-domain system measure of one charge equals Liénard's (`lienard_energy`, a
/// different quadrature of the same integral): required 1e-8.
#[test]
fn s6_system_spectrum_integrates_to_the_system_energy() {
    let c = 1.0;
    let (a0, sigma) = (0.05, 2.0);
    // A charge at β = 0.6 along x from `x0`, kicked along y around `t_k`.
    let flight = |x0: DVec3, t_k: f64| -> Vec<Emission> {
        let v0 = 0.6;
        let n_steps = 20_000;
        (0..=n_steps)
            .map(|i| {
                let t = -20.0 + 40.0 * f64::from(i) / f64::from(n_steps);
                let u = t - t_k;
                let (th, ch) = ((u / sigma).tanh(), (u / sigma).cosh());
                let x = x0 + DVec3::new(v0 * t, a0 * sigma * (u + sigma * ch.ln()), 0.0);
                let v = DVec3::new(v0, a0 * sigma * (1.0 + th), 0.0);
                let a = DVec3::new(0.0, a0 * (1.0 - th * th), 0.0);
                (t, x, v, a)
            })
            .collect()
    };
    let first = flight(DVec3::ZERO, 0.0);
    let second = flight(DVec3::new(0.5, 1.5, 0.0), 1.5);
    let sources: [(f64, &[Emission]); 2] = [(1.0, &first), (-0.7, &second)];
    let (mut worst, mut worst_one) = (0.0_f64, 0.0_f64);
    for deg in [0.0_f64, 20.0, 60.0, 135.0] {
        let n = DVec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
        let time_domain = system_lienard_energy(&sources, c, n);
        let kappa = 1.0 - n.x * 0.6;
        let top = 40.0 / (sigma * kappa);
        let spectral = system_band_energy(&sources, c, n, 0.0, top);
        let err = spectral / time_domain - 1.0;
        let one =
            system_lienard_energy(&sources[..1], c, n) / lienard_energy(&first, 1.0, c, n) - 1.0;
        println!(
            "S6 at {deg}°: ∫ spectrum {spectral:.12e}, time domain {time_domain:.12e}, \
             rel. error {err:.1e}; one charge against Liénard {one:.1e}"
        );
        worst = worst.max(err.abs());
        worst_one = worst_one.max(one.abs());
    }
    assert!(worst < 1e-8, "{worst:.2e}");
    assert!(worst_one < 1e-8, "{worst_one:.2e}");
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
/// `(q²/4π²c) β² sin²θ / (1 − β cos θ)²`. (b) The same charge stopped smoothly within
/// τ = 0.002 (velocity `v₀ / (1 + e^{t/τ})`), integrated by the acceleration form alone:
/// its spectrum divided by (a) is `|R|²` with `R = ∫ W e^{iωψ} dt` (W = f′/(−f₀),
/// f = β sin θ/(1 − β cos θ), ψ = ∫₀ᵗ (1 − β cos θ) dt′), so `|R|² − 1 = −(ωτ)² V + O(ω⁴)`
/// with V the variance of ψ/τ under W. `scripts/wolfram/s4_sudden_stop.wls` computes V
/// and the exact `|R|² − 1` (Wolfram Engine 14.2, 40-digit quadrature), tabulated below.
/// Criteria: (a) 1e-12; (b) `|R|² − 1` to 1e-7 absolute (set before measuring: the smallest
/// value, 7e-7, to ~15 %, the largest to 1e-4 of itself).
/// Measured: (a) ≤ 3e-16; (b) ≤ 2.1e-10 (the first, linear Filon rule: 1.8e-8).
#[test]
fn s4_sudden_stop() {
    // (θ in degrees, V, exact |R|² − 1 at ωτ = 1e-3, 4e-3, 1.6e-2), from the Wolfram script.
    // Wolfram's 17 digits as printed (rounded to f64 by the compiler).
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const EXACT: [(f64, f64, [f64; 3]); 4] = [
        (
            10.0,
            0.69795801813479668,
            [
                -6.9795744079834772e-7,
                -0.000011167180494128837,
                -0.00017863942547117701,
            ],
        ),
        (
            36.87,
            1.1843553491307757,
            [
                -1.184354188323343e-6,
                -0.000018949388423757935,
                -0.00030311891370784225,
            ],
        ),
        (
            60.0,
            1.9739208802178717,
            [
                -1.9739183345964407e-6,
                -0.000031582082415207885,
                -0.00050515696251552457,
            ],
        ),
        (
            120.0,
            4.605815387175034,
            [
                -4.605802174204176e-6,
                -0.000073689663793984662,
                -0.0011782233346276761,
            ],
        ),
    ];
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
    let (mut worst_a, mut worst_b): (f64, f64) = (0.0, 0.0);
    for (deg, v, exact) in EXACT {
        let th = deg.to_radians();
        let n = DVec3::new(th.cos(), th.sin(), 0.0);
        let jackson = q * q / (4.0 * PI * PI * c) * beta * beta * th.sin().powi(2)
            / (1.0 - beta * th.cos()).powi(2);
        for (k, wt) in [1e-3, 4e-3, 1.6e-2].into_iter().enumerate() {
            let omega = wt / tau;
            let a = spectrum(&uniform, q, c, n, omega, 1.0, 1, true)[0];
            let b = spectrum(&smooth, q, c, n, omega, 1.0, 1, false)[0];
            let (ea, eb) = (a / jackson - 1.0, b / jackson - 1.0);
            println!(
                "S4 at {deg}°, ωτ = {wt:.0e}: abrupt {ea:.1e}; smooth |R|² − 1 = {eb:.9e}, exact {:.9e} (leading −(ωτ)²V {:.9e}), difference {:.1e}",
                exact[k],
                -wt * wt * v,
                eb - exact[k]
            );
            worst_a = worst_a.max(ea.abs());
            worst_b = worst_b.max((eb - exact[k]).abs());
        }
    }
    println!("S4: abrupt ≤ {worst_a:.1e}; smooth against Wolfram ≤ {worst_b:.1e}");
    assert!(worst_a < 1e-12 && worst_b < 1e-7);
}
