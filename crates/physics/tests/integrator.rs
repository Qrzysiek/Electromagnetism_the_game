//! Validation of the DOP853 integrator against problems with known solutions and against
//! an independent implementation (`ode_solvers`). Measured values are recorded in
//! PHYSICS.md §5.

use std::f64::consts::PI;

use ode_solvers::{SVector, System};
use physics::integrator::{Dop853, OdeSystem, Settings};

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn integrate(sys: &impl OdeSystem, y0: &[f64], t_end: f64, s: Settings) -> Dop853 {
    let mut int = Dop853::new(sys, 0.0, y0, s);
    while !int
        .step(sys, t_end, f64::INFINITY)
        .expect("integration failed")
    {}
    int
}

/// Settings that make the controller take exactly `h` every step.
fn fixed_step(h: f64) -> Settings {
    Settings {
        rtol: 1e30,
        atol: 1e30,
        h_init: Some(h),
        h_max: h,
        max_steps: 10_000_000,
        ..Settings::default()
    }
}

// --- Arenstorf orbit (Hairer, Nørsett, Wanner, Solving ODEs I, §II.0): periodic. ---------

const MU: f64 = 0.012277471;
// Reference values are quoted with all published digits.
#[allow(clippy::excessive_precision)]
const ARENSTORF_T: f64 = 17.065_216_560_157_962_558_891_720_624_9;
#[allow(clippy::excessive_precision)]
const ARENSTORF_Y0: [f64; 4] = [0.994, 0.0, 0.0, -2.001_585_106_379_082_522_405_378_622_24];

fn arenstorf_rhs(y: &[f64], dy: &mut [f64]) {
    let mu1 = 1.0 - MU;
    let d1 = ((y[0] + MU).powi(2) + y[1].powi(2)).powf(1.5);
    let d2 = ((y[0] - mu1).powi(2) + y[1].powi(2)).powf(1.5);
    dy[0] = y[2];
    dy[1] = y[3];
    dy[2] = y[0] + 2.0 * y[3] - mu1 * (y[0] + MU) / d1 - MU * (y[0] - mu1) / d2;
    dy[3] = y[1] - 2.0 * y[2] - mu1 * y[1] / d1 - MU * y[1] / d2;
}

struct Arenstorf;

impl OdeSystem for Arenstorf {
    fn dim(&self) -> usize {
        4
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        arenstorf_rhs(y, dy);
    }
}

struct ArenstorfOde;

impl System<f64, SVector<f64, 4>> for ArenstorfOde {
    fn system(&self, _t: f64, y: &SVector<f64, 4>, dy: &mut SVector<f64, 4>) {
        arenstorf_rhs(y.as_slice(), dy.as_mut_slice());
    }
}

#[test]
fn arenstorf_orbit_closes() {
    let mut s = Settings::with_tolerance(1e-12, 1e-12);
    s.h_max = ARENSTORF_T;
    let int = integrate(&Arenstorf, &ARENSTORF_Y0, ARENSTORF_T, s);
    let err = max_abs_diff(int.y(), &ARENSTORF_Y0);
    println!("Arenstorf, tol 1e-12: error after one period = {err:.3e}");
    assert!(err < 5e-9, "error {err:.3e}");
}

#[test]
fn matches_independent_implementation_step_for_step() {
    for k in 6..=12 {
        let tol = 10f64.powi(-k);
        let mut s = Settings::with_tolerance(tol, tol);
        s.h_max = ARENSTORF_T;
        s.dense = false;
        let ours = integrate(&Arenstorf, &ARENSTORF_Y0, ARENSTORF_T, s);

        let mut theirs = ode_solvers::Dop853::from_param(
            ArenstorfOde,
            0.0,
            ARENSTORF_T,
            ARENSTORF_T,
            SVector::<f64, 4>::from_column_slice(&ARENSTORF_Y0),
            tol,
            tol,
            0.9,
            0.0,
            0.333,
            6.0,
            ARENSTORF_T,
            0.0,
            100_000,
            1000,
            ode_solvers::dop_shared::OutputType::Sparse,
        );
        let st = theirs.integrate().unwrap();
        let y_theirs = theirs.y_out().last().unwrap().as_slice().to_vec();

        let so = ours.stats();
        assert_eq!(so.n_accept, st.accepted_steps as u64, "tol {tol:e}");
        assert_eq!(so.n_reject, st.rejected_steps as u64, "tol {tol:e}");
        // The implementations round differently (ode_solvers forms y + Σ(h a) k instead of
        // y + h Σ a k); the Arenstorf orbit amplifies such differences to ~1e-10.
        let diff = max_abs_diff(ours.y(), &y_theirs);
        println!("tol {tol:e}: |ours - ode_solvers| = {diff:.2e}");
        assert!(diff < 1e-8, "tol {tol:e}: difference {diff:.3e}");
    }
}

// --- Convergence order ------------------------------------------------------------------

struct Oscillator;

impl OdeSystem for Oscillator {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1];
        dy[1] = -y[0];
    }
}

struct Kepler;

impl OdeSystem for Kepler {
    fn dim(&self) -> usize {
        4
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let r3 = (y[0] * y[0] + y[1] * y[1]).powf(1.5);
        dy[0] = y[2];
        dy[1] = y[3];
        dy[2] = -y[0] / r3;
        dy[3] = -y[1] / r3;
    }
}

#[test]
fn global_order_is_eight() {
    // Kepler orbit, e = 0.5, a = 1, GM = 1: period 2π; exact final state = initial state.
    let y0 = [0.5, 0.0, 0.0, 3f64.sqrt()];
    let err = |n: u32| {
        let int = integrate(&Kepler, &y0, 2.0 * PI, fixed_step(2.0 * PI / f64::from(n)));
        max_abs_diff(int.y(), &y0)
    };
    let (e1, e2) = (err(128), err(256));
    let order = (e1 / e2).log2();
    println!("Kepler fixed step: e(128) = {e1:.3e}, e(256) = {e2:.3e}, order = {order:.2}");
    assert!((7.5..8.5).contains(&order), "observed order {order:.2}");

    // Harmonic oscillator over 10 time units.
    let y0 = [1.0, 0.0];
    let err = |n: u32| {
        let int = integrate(&Oscillator, &y0, 10.0, fixed_step(10.0 / f64::from(n)));
        max_abs_diff(int.y(), &[10f64.cos(), -10f64.sin()])
    };
    let (e1, e2) = (err(20), err(40));
    let order = (e1 / e2).log2();
    println!("Oscillator fixed step: e(20) = {e1:.3e}, e(40) = {e2:.3e}, order = {order:.2}");
    assert!((7.5..8.5).contains(&order), "observed order {order:.2}");
}

#[test]
fn dense_output_local_order_is_eight() {
    let max_err = |h: f64| {
        let mut int = Dop853::new(&Oscillator, 0.0, &[1.0, 0.0], fixed_step(h));
        int.step(&Oscillator, 100.0, f64::INFINITY).unwrap();
        let d = int.dense();
        (1..200)
            .map(|i| {
                let t = h * f64::from(i) / 200.0;
                let e0 = (d.eval_component(0, t) - t.cos()).abs();
                let e1 = (d.eval_component(1, t) + t.sin()).abs();
                e0.max(e1)
            })
            .fold(0.0, f64::max)
    };
    for (h1, h2) in [(0.8, 0.4), (0.4, 0.2)] {
        let order = (max_err(h1) / max_err(h2)).log2();
        println!("dense output: h {h1} -> {h2}: local order {order:.2}");
        assert!((7.8..8.2).contains(&order), "observed order {order:.2}");
    }
}

#[test]
fn dense_output_matches_step_endpoints() {
    let mut int = Dop853::new(
        &Kepler,
        0.0,
        &[0.5, 0.0, 0.0, 3f64.sqrt()],
        Settings::default(),
    );
    for _ in 0..20 {
        let y_start = int.y().to_vec();
        int.step(&Kepler, 100.0, f64::INFINITY).unwrap();
        let d = int.dense();
        let mut a = vec![0.0; 4];
        let mut b = vec![0.0; 4];
        d.eval(d.t_start(), &mut a);
        d.eval(d.t_end(), &mut b);
        // Exact up to rounding: at θ = 1 the interpolant evaluates y + (y_new - y).
        let scale = y_start
            .iter()
            .chain(int.y())
            .fold(1.0f64, |m, v| m.max(v.abs()));
        let tol = 8.0 * f64::EPSILON * scale;
        // At θ = 0 the interpolant returns y exactly.
        assert_eq!(max_abs_diff(&a, &y_start).to_bits(), 0);
        let e = max_abs_diff(&b, int.y());
        assert!(e <= tol, "endpoint mismatch {e:.3e} > {tol:.3e}");
    }
}

// --- Relativistic hyperbolic motion (prototype of T6) -----------------------------------

/// Units m = c = qE = 1. State (x, p); exact solution p = t, x = sqrt(1 + t^2) - 1.
struct Hyperbolic;

impl OdeSystem for Hyperbolic {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = y[1] / (1.0 + y[1] * y[1]).sqrt();
        dy[1] = 1.0;
    }
}

#[test]
fn relativistic_hyperbolic_motion() {
    for t_end in [1.0, 10.0, 1000.0] {
        let int = integrate(
            &Hyperbolic,
            &[0.0, 0.0],
            t_end,
            Settings::with_tolerance(1e-12, 1e-12),
        );
        let x = (1.0 + t_end * t_end).sqrt() - 1.0;
        let rel = (int.y()[0] - x).abs() / x;
        println!("hyperbolic motion to t = {t_end}: relative error {rel:.2e}");
        assert!(rel < 1e-12, "t = {t_end}: relative error {rel:.3e}");
    }
}

// --- API behaviour ----------------------------------------------------------------------

/// y' = cos t, y(0) = 0: checks that stage times are used correctly.
struct Cosine;

impl OdeSystem for Cosine {
    fn dim(&self) -> usize {
        1
    }
    fn rhs(&self, t: f64, _y: &[f64], dy: &mut [f64]) {
        dy[0] = t.cos();
    }
}

#[test]
fn non_autonomous_right_hand_side() {
    let int = integrate(
        &Cosine,
        &[0.0],
        10.0,
        Settings::with_tolerance(1e-12, 1e-12),
    );
    let err = (int.y()[0] - 10f64.sin()).abs();
    assert!(err < 1e-11, "error {err:.3e}");
}

#[test]
fn step_respects_limit_and_cap() {
    let mut int = Dop853::new(
        &Kepler,
        0.0,
        &[0.5, 0.0, 0.0, 3f64.sqrt()],
        Settings::default(),
    );
    let cap = 0.01;
    loop {
        let t0 = int.t();
        let done = int.step(&Kepler, 1.0, cap).unwrap();
        // t0 + h is rounded, so allow one ulp of t.
        assert!(
            int.t() - t0 <= cap + 2.0 * f64::EPSILON,
            "step {}",
            int.t() - t0
        );
        if done {
            break;
        }
    }
    assert_eq!(int.t().to_bits(), 1.0f64.to_bits());
}
