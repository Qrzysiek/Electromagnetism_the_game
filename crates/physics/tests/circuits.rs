//! Validation tests Z4–Z8 of the lumped circuits (`circuit.rs`, PHYSICS.md §2.9): modified
//! nodal analysis integrated by RADAU5, against closed forms. Common requirement, set
//! before measuring: at the tolerance 1e-12, every value within 1e-8 at evenly spaced
//! times, i.e. of the dense output, which inside a step is the collocation polynomial, of
//! lower order than the steps' ends (there RADAU5's global error is about its tolerance,
//! I2–I4).

#![allow(clippy::disallowed_methods)] // references

use physics::circuit::{Circuit, CircuitSolution, Component, Waveform};

const TOL: f64 = 1e-12;

/// The worst |got − want| over `count` evenly spaced times in `[t0, t1]` (the dense output)
/// and the times of the steps' ends.
fn worst(
    sol: &CircuitSolution,
    t0: f64,
    t1: f64,
    count: u32,
    got: impl Fn(&CircuitSolution, f64) -> f64,
    want: impl Fn(f64) -> f64,
) -> f64 {
    (0..=count)
        .map(|k| t0 + (t1 - t0) * f64::from(k) / f64::from(count))
        .map(|t| (got(sol, t) - want(t)).abs())
        .fold(0.0, f64::max)
}

/// Z4: a DC source charging a capacitor through a resistor (V₀ = 2, R = 0.5, C = 2, so
/// RC = 1), from an uncharged capacitor: `v_C = V₀ (1 − e^{−t/RC})`, the current
/// `(V₀/R) e^{−t/RC}`.
#[test]
fn z4_rc_charging() {
    let (v0, r, c) = (2.0, 0.5, 2.0);
    let circuit = Circuit {
        nodes: 2,
        components: vec![
            Component::Source {
                a: 1,
                b: 0,
                wave: Waveform::Dc(v0),
            },
            Component::Resistor { a: 1, b: 2, r },
            Component::Capacitor { a: 2, b: 0, c },
        ],
        ..Circuit::default()
    };
    let sol = circuit
        .solve(&vec![0.0; circuit.unknowns()], 10.0, TOL, v0, 0.05)
        .expect("solves");
    let vc = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| s.potential(2, t),
        |t| v0 * (1.0 - (-t / (r * c)).exp()),
    );
    // The source's current flows out of node 1 into the source: minus the charging current.
    let i = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| -s.current(0, t),
        |t| v0 / r * (-t / (r * c)).exp(),
    );
    println!(
        "Z4 RC: v_C within {vc:.1e}, the current within {i:.1e}; {} steps",
        sol.steps()
    );
    assert!(vc < 1e-8 && i < 1e-8, "Z4 {vc:.3e} {i:.3e}");
}

/// Z5: an LC tank (C = 0.5, L = 2: ω = 1) started with the capacitor at 1: `v = cos t`,
/// the inductor's current `C sin t`, and the energy `½ C v² + ½ L i² = 0.25` over eight
/// periods.
#[test]
fn z5_lc_oscillation() {
    let (c, l) = (0.5, 2.0);
    let circuit = Circuit {
        nodes: 1,
        components: vec![
            Component::Capacitor { a: 1, b: 0, c },
            Component::Inductor { a: 1, b: 0, l },
        ],
        ..Circuit::default()
    };
    let t_end = 16.0 * std::f64::consts::PI;
    let sol = circuit
        .solve(&[1.0, 0.0], t_end, TOL, 1.0, 0.05)
        .expect("solves");
    let v = worst(&sol, 0.0, t_end, 8001, |s, t| s.potential(1, t), f64::cos);
    let i = worst(
        &sol,
        0.0,
        t_end,
        8001,
        |s, t| s.current(1, t),
        |t| c * t.sin(),
    );
    let e = worst(
        &sol,
        0.0,
        t_end,
        8001,
        |s, t| 0.5 * c * s.potential(1, t).powi(2) + 0.5 * l * s.current(1, t).powi(2),
        |_| 0.25,
    );
    println!(
        "Z5 LC: v within {v:.1e}, i within {i:.1e}, energy within {e:.1e}; {} steps",
        sol.steps()
    );
    assert!(
        v < 1e-8 && i < 1e-8 && e < 1e-8,
        "Z5 {v:.3e} {i:.3e} {e:.3e}"
    );
}

/// Z6: charge sharing through a switch. C₁ = 1 charged to V₀ = 1, C₂ = 3 empty, a switch
/// (r = 0.2) between them closing at t = 1: both approach `C₁V₀/(C₁ + C₂)` with
/// `τ = r C₁C₂/(C₁ + C₂)`, and the heat in the switch, `∫ (v₁ − v₂)²/r dt`, is
/// `½ (C₁C₂/(C₁ + C₂)) V₀²` whatever r (here by Simpson's rule on the dense output,
/// 1e-9 relative required). Before the switch closes nothing moves.
#[test]
fn z6_switch_charge_sharing() {
    let (c1, c2, r, v0, ts) = (1.0, 3.0, 0.2, 1.0, 1.0);
    let circuit = Circuit {
        nodes: 2,
        components: vec![
            Component::Capacitor { a: 1, b: 0, c: c1 },
            Component::Capacitor { a: 2, b: 0, c: c2 },
            Component::Switch {
                a: 1,
                b: 2,
                r,
                closed: false,
                toggles: vec![ts],
            },
        ],
        ..Circuit::default()
    };
    let t_end = 5.0;
    let sol = circuit
        .solve(&[v0, 0.0], t_end, TOL, 1.0, 0.01)
        .expect("solves");
    let vf = c1 * v0 / (c1 + c2);
    let tau = r * c1 * c2 / (c1 + c2);
    let decay = |t: f64| {
        if t <= ts {
            1.0
        } else {
            (-(t - ts) / tau).exp()
        }
    };
    let w1 = worst(
        &sol,
        0.0,
        t_end,
        5001,
        |s, t| s.potential(1, t),
        |t| vf + (v0 - vf) * decay(t),
    );
    let w2 = worst(
        &sol,
        0.0,
        t_end,
        5001,
        |s, t| s.potential(2, t),
        |t| vf * (1.0 - decay(t)),
    );
    // The heat, from the switch's closing.
    let m = 200_000;
    let h = (t_end - ts) / f64::from(m);
    let p = |t: f64| (sol.potential(1, t) - sol.potential(2, t)).powi(2) / r;
    let mut heat = p(ts) + p(t_end);
    for k in 1..m {
        heat += if k % 2 == 1 { 4.0 } else { 2.0 } * p(ts + h * f64::from(k));
    }
    heat *= h / 3.0;
    let want = 0.5 * c1 * c2 / (c1 + c2) * v0 * v0;
    // The rest of the heat after t_end: (v₁ − v₂)² decays as e^{−2t/τ}.
    let tail = want * (-2.0 * (t_end - ts) / tau).exp();
    let rel = (heat + tail) / want - 1.0;
    println!(
        "Z6 switch: v₁ within {w1:.1e}, v₂ within {w2:.1e}; heat {heat:.12} + {tail:.1e} \
         (closed form {want:.12}): {rel:.1e}; {} steps",
        sol.steps()
    );
    assert!(w1 < 1e-8 && w2 < 1e-8, "Z6 {w1:.3e} {w2:.3e}");
    assert!(rel.abs() < 1e-9, "Z6 heat {rel:.3e}");
}

/// Z7: a trapezoid pulse (low 0, high 1, delay 1, rise 0.5, width 2, fall 0.5) through
/// R = 1 into C = 1, against superposed ramp responses: a ramp `k (t − t₀)` from t₀ gives
/// `k [(t − t₀) − RC (1 − e^{−(t − t₀)/RC})]`. The corners are breakpoints.
#[test]
fn z7_pulse_into_rc() {
    let wave = Waveform::Pulse {
        low: 0.0,
        high: 1.0,
        delay: 1.0,
        rise: 0.5,
        width: 2.0,
        fall: 0.5,
        period: 0.0,
    };
    let circuit = Circuit {
        nodes: 2,
        components: vec![
            Component::Source { a: 1, b: 0, wave },
            Component::Resistor { a: 1, b: 2, r: 1.0 },
            Component::Capacitor { a: 2, b: 0, c: 1.0 },
        ],
        ..Circuit::default()
    };
    let sol = circuit
        .solve(&vec![0.0; circuit.unknowns()], 8.0, TOL, 1.0, 0.05)
        .expect("solves");
    assert_eq!(sol.breakpoints, vec![1.0, 1.5, 3.5, 4.0]);
    let ramp = |t: f64, t0: f64, k: f64| {
        if t <= t0 {
            0.0
        } else {
            k * ((t - t0) - (1.0 - (-(t - t0)).exp()))
        }
    };
    let want =
        |t: f64| ramp(t, 1.0, 2.0) - ramp(t, 1.5, 2.0) - ramp(t, 3.5, 2.0) + ramp(t, 4.0, 2.0);
    let w = worst(&sol, 0.0, 8.0, 8001, |s, t| s.potential(2, t), want);
    println!("Z7 pulse into RC: within {w:.1e}; {} steps", sol.steps());
    assert!(w < 1e-8, "Z7 {w:.3e}");
}

/// Z8: two coupled inductors (L₁ = 1, L₂ = 2, mutual M = 0.8): a DC step V₀ = 1 through
/// R₁ = 1 into the primary, the secondary closed by R₂ = 2. `[L₁ M; M L₂] i' = [V₀ − R₁i₁;
/// −R₂i₂]`, solved in closed form by the eigenvectors of its 2×2 matrix: the primary
/// current rises to V₀/R₁ while the secondary's swings and dies. Also the currents' rates
/// of change (the derivative of the dense output, which drives a coil's induced field):
/// required 1e-7, a power of the step less accurate than the values.
#[test]
fn z8_mutual_inductance() {
    let (l1, l2, m, r1, r2, v0) = (1.0, 2.0, 0.8, 1.0, 2.0, 1.0);
    let circuit = Circuit {
        nodes: 3,
        components: vec![
            Component::Source {
                a: 1,
                b: 0,
                wave: Waveform::Dc(v0),
            },
            Component::Resistor { a: 1, b: 2, r: r1 },
            Component::Inductor { a: 2, b: 0, l: l1 },
            Component::Inductor { a: 3, b: 0, l: l2 },
            Component::Resistor { a: 3, b: 0, r: r2 },
        ],
        mutual: vec![(2, 3, m)],
        ..Circuit::default()
    };
    let sol = circuit
        .solve(&vec![0.0; circuit.unknowns()], 10.0, TOL, 1.0, 0.05)
        .expect("solves");
    // i' = K i + g with K = −L⁻¹ R, g = L⁻¹ (V₀, 0); the steady state i∞ = (V₀/R₁, 0).
    let det = l1 * l2 - m * m;
    let linv = [[l2 / det, -m / det], [-m / det, l1 / det]];
    let k = [
        [-linv[0][0] * r1, -linv[0][1] * r2],
        [-linv[1][0] * r1, -linv[1][1] * r2],
    ];
    let tr = k[0][0] + k[1][1];
    let dt = k[0][0] * k[1][1] - k[0][1] * k[1][0];
    let disc = (tr * tr / 4.0 - dt).sqrt();
    let (la, lb) = (tr / 2.0 + disc, tr / 2.0 - disc);
    // Eigenvectors (k01, λ − k00).
    let (va, vb) = ([k[0][1], la - k[0][0]], [k[0][1], lb - k[0][0]]);
    // i(0) − i∞ = −i∞ = ca va + cb vb.
    let d0 = [-v0 / r1, 0.0];
    let den = va[0] * vb[1] - va[1] * vb[0];
    let ca = (d0[0] * vb[1] - d0[1] * vb[0]) / den;
    let cb = (va[0] * d0[1] - va[1] * d0[0]) / den;
    let current = |t: f64, j: usize| {
        let inf = if j == 0 { v0 / r1 } else { 0.0 };
        inf + ca * va[j] * (la * t).exp() + cb * vb[j] * (lb * t).exp()
    };
    let w1 = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| s.current(2, t),
        |t| current(t, 0),
    );
    let w2 = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| s.current(3, t),
        |t| current(t, 1),
    );
    let rate =
        |t: f64, j: usize| ca * va[j] * la * (la * t).exp() + cb * vb[j] * lb * (lb * t).exp();
    let d1 = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| s.current_rate(2, t),
        |t| rate(t, 0),
    );
    let d2 = worst(
        &sol,
        0.0,
        10.0,
        4001,
        |s, t| s.current_rate(3, t),
        |t| rate(t, 1),
    );
    println!(
        "Z8 coupled inductors (λ = {la:.4}, {lb:.4}): the primary within {w1:.1e}, the \
         secondary within {w2:.1e}; their rates within {d1:.1e}, {d2:.1e}; {} steps",
        sol.steps()
    );
    assert!(w1 < 1e-8 && w2 < 1e-8, "Z8 {w1:.3e} {w2:.3e}");
    assert!(d1 < 1e-7 && d2 < 1e-7, "Z8 rates {d1:.3e} {d2:.3e}");
}

/// Z9: a source `V₀ sin ωt` (V₀ = 1, ω = 1.3) with internal resistance R₁ = 1 drives a
/// lossy tank: C = 1 in parallel with L = 1 in series with R₂ = 0.5, from rest. (a) After
/// the transient (its modes decay as e^{−0.75 t}: e^{−30} by t = 40) the tank's voltage
/// and the inductor's current against the phasors: `Z_tank = 1/(iωC + 1/(R₂ + iωL))`,
/// `V_tank = V₀ Z_tank/(R₁ + Z_tank)`, `I_L = V_tank/(R₂ + iωL)`, over t ∈ [40, 60];
/// required 1e-8 as above. (b) The energy balance over [0, 60]: the source's work
/// `∫ v₁ i dt` = the heat in R₁ and R₂ + `½ C v² + ½ L i_L²` at the end (Simpson's rule on
/// the dense output); required within 1e-8 of the work (the values' accuracy).
#[test]
fn z9_driven_tank_and_energy_balance() {
    let (v0, w, r1, c, l, r2) = (1.0, 1.3, 1.0, 1.0, 1.0, 0.5);
    let circuit = Circuit {
        nodes: 3,
        components: vec![
            Component::Source {
                a: 1,
                b: 0,
                wave: Waveform::Sine {
                    offset: 0.0,
                    amplitude: v0,
                    omega: w,
                    phase: 0.0,
                },
            },
            Component::Resistor { a: 1, b: 2, r: r1 },
            Component::Capacitor { a: 2, b: 0, c },
            Component::Inductor { a: 2, b: 3, l },
            Component::Resistor { a: 3, b: 0, r: r2 },
        ],
        ..Circuit::default()
    };
    let t_end = 60.0;
    let sol = circuit
        .solve(&vec![0.0; circuit.unknowns()], t_end, TOL, v0, 0.05)
        .expect("solves");
    // Complex arithmetic on (re, im) pairs.
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let inv = |a: (f64, f64)| {
        let d = a.0 * a.0 + a.1 * a.1;
        (a.0 / d, -a.1 / d)
    };
    let add = |a: (f64, f64), b: (f64, f64)| (a.0 + b.0, a.1 + b.1);
    let branch = (r2, w * l);
    let tank = inv(add((0.0, w * c), inv(branch)));
    let v_tank = mul((v0, 0.0), mul(tank, inv(add((r1, 0.0), tank))));
    let i_l = mul(v_tank, inv(branch));
    // v(t) = Im(V e^{iωt}).
    let at = |z: (f64, f64), t: f64| z.0 * (w * t).sin() + z.1 * (w * t).cos();
    let wv = worst(
        &sol,
        40.0,
        t_end,
        4001,
        |s, t| s.potential(2, t),
        |t| at(v_tank, t),
    );
    let wi = worst(
        &sol,
        40.0,
        t_end,
        4001,
        |s, t| s.current(3, t),
        |t| at(i_l, t),
    );
    // (b) The energy balance. The source's current flows from node 1 through it to ground:
    // it delivers −i.
    let m = 600_000;
    let h = t_end / f64::from(m);
    let simpson = |f: &dyn Fn(f64) -> f64| {
        let mut s = f(0.0) + f(t_end);
        for k in 1..m {
            s += if k % 2 == 1 { 4.0 } else { 2.0 } * f(h * f64::from(k));
        }
        s * h / 3.0
    };
    let work = simpson(&|t| -sol.potential(1, t) * sol.current(0, t));
    let heat = simpson(&|t| {
        (sol.potential(1, t) - sol.potential(2, t)).powi(2) / r1 + sol.current(3, t).powi(2) * r2
    });
    let stored =
        0.5 * c * sol.potential(2, t_end).powi(2) + 0.5 * l * sol.current(3, t_end).powi(2);
    let balance = (work - heat - stored) / work;
    println!(
        "Z9 driven tank: |V_tank| = {:.6}, the tank's voltage within {wv:.1e}, the inductor's \
         current within {wi:.1e}; work {work:.12}, heat {heat:.12}, stored {stored:.3e}: \
         balance {balance:.1e}; {} steps",
        v_tank.0.hypot(v_tank.1),
        sol.steps()
    );
    assert!(wv < 1e-8 && wi < 1e-8, "Z9 phasors {wv:.3e} {wi:.3e}");
    assert!(balance.abs() < 1e-8, "Z9 balance {balance:.3e}");
}
