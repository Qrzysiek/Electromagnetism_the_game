//! M1 integrator evaluation. Prints tables comparing our DOP853 port with the `ode_solvers`
//! crate on reference problems with known solutions. Run with
//! `cargo run --release -p integrator_eval`.

use std::f64::consts::PI;
use std::hint::black_box;
use std::time::Instant;

use ode_solvers::{SVector, System};
use physics::integrator::{Dop853, OdeSystem, Settings};

// ---------------------------------------------------------------------------------------
// Arenstorf orbit (restricted three-body problem), Hairer, Nørsett, Wanner, Solving ODEs I,
// §II.0, eq. (0.1). The orbit is periodic, so the exact state at t = T equals the initial
// state.
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

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

struct Run {
    y: Vec<f64>,
    steps: u64,
    rejects: u64,
    fcn: u64,
    micros: f64,
}

fn ours(sys: &impl OdeSystem, y0: &[f64], t_end: f64, mut s: Settings, reps: u32) -> Run {
    s.h_max = t_end; // same default as dop853.f / ode_solvers
    s.max_steps = 10_000_000;
    let mut last = None;
    let start = Instant::now();
    for _ in 0..reps {
        let mut int = Dop853::new(sys, 0.0, y0, s);
        while !int
            .step(sys, t_end, f64::INFINITY)
            .expect("integration failed")
        {}
        last = Some(int);
    }
    let micros = start.elapsed().as_secs_f64() * 1e6 / f64::from(reps);
    let int = last.unwrap();
    let st = int.stats();
    Run {
        y: black_box(int.y().to_vec()),
        steps: st.n_accept + st.n_reject,
        rejects: st.n_reject,
        fcn: st.n_fcn,
        micros,
    }
}

fn theirs(tol: f64, reps: u32) -> Run {
    let y0 = SVector::<f64, 4>::from_column_slice(&ARENSTORF_Y0);
    let mut out = None;
    let start = Instant::now();
    for _ in 0..reps {
        let mut s = ode_solvers::Dop853::from_param(
            ArenstorfOde,
            0.0,
            ARENSTORF_T,
            ARENSTORF_T,
            y0,
            tol,
            tol,
            0.9,
            0.0,
            0.333,
            6.0,
            ARENSTORF_T,
            0.0,
            10_000_000,
            1000,
            ode_solvers::dop_shared::OutputType::Sparse,
        );
        let stats = s.integrate().expect("integration failed");
        out = Some((s.y_out().last().unwrap().as_slice().to_vec(), stats));
    }
    let micros = start.elapsed().as_secs_f64() * 1e6 / f64::from(reps);
    let (y, st) = out.unwrap();
    Run {
        y: black_box(y),
        steps: (st.accepted_steps + st.rejected_steps) as u64,
        rejects: st.rejected_steps as u64,
        fcn: st.num_eval as u64,
        micros,
    }
}

fn arenstorf_table() {
    println!("## A. Arenstorf orbit, one period (rtol = atol = tol, dense output off)\n");
    println!(
        "| tol | error (ours) | error (ode_solvers) | |ours - theirs| | steps ours/theirs | \
         rejected ours/theirs | f-evals ours/theirs | time ours / theirs [us] |"
    );
    println!("|---|---|---|---|---|---|---|---|");
    for k in 6..=14 {
        let tol = 10f64.powi(-k);
        let mut s = Settings::with_tolerance(tol, tol);
        s.dense = false;
        let a = ours(&Arenstorf, &ARENSTORF_Y0, ARENSTORF_T, s, 50);
        let b = theirs(tol, 50);
        println!(
            "| 1e-{k} | {:.2e} | {:.2e} | {:.1e} | {}/{} | {}/{} | {}/{} | {:.0} / {:.0} |",
            max_abs_diff(&a.y, &ARENSTORF_Y0),
            max_abs_diff(&b.y, &ARENSTORF_Y0),
            max_abs_diff(&a.y, &b.y),
            a.steps,
            b.steps,
            a.rejects,
            b.rejects,
            a.fcn,
            b.fcn,
            a.micros,
            b.micros
        );
    }
    println!();
}

// ---------------------------------------------------------------------------------------
// Kepler problem, eccentricity 0.5, a = 1, GM = 1: period 2π. Fixed-step runs measure the
// convergence order of the method.
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

const KEPLER_E: f64 = 0.5;

fn kepler_y0() -> [f64; 4] {
    let e = KEPLER_E;
    [1.0 - e, 0.0, 0.0, ((1.0 + e) / (1.0 - e)).sqrt()]
}

/// Settings that make the controller take exactly `h` every step.
fn fixed_step(h: f64) -> Settings {
    Settings {
        rtol: 1e30,
        atol: 1e30,
        h_init: Some(h),
        h_max: h,
        max_steps: 10_000_000,
        dense: false,
        ..Settings::default()
    }
}

fn kepler_order_table() {
    println!("## B. Convergence order, fixed steps, Kepler e = 0.5, one period\n");
    println!("| steps | error | observed order log2(e(N/2)/e(N)) |");
    println!("|---|---|---|");
    let y0 = kepler_y0();
    let mut prev: Option<f64> = None;
    for n in [64u32, 128, 256, 512, 1024] {
        let h = 2.0 * PI / f64::from(n);
        let mut int = Dop853::new(&Kepler, 0.0, &y0, fixed_step(h));
        while !int.step(&Kepler, 2.0 * PI, f64::INFINITY).unwrap() {}
        let err = max_abs_diff(int.y(), &y0);
        let order = prev.map_or("-".to_string(), |p| format!("{:.2}", (p / err).log2()));
        println!("| {n} | {err:.3e} | {order} |");
        prev = Some(err);
    }
    println!();
}

// ---------------------------------------------------------------------------------------
// Dense output order: harmonic oscillator y = (cos t, -sin t). One step of size h from the
// exact initial value; maximum interpolation error over the step. The interpolant has
// order 7, so its local error is O(h^8).
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

fn dense_order_table() {
    println!("## C. Dense output order, harmonic oscillator, one step\n");
    println!("| h | max interpolation error | observed local order |");
    println!("|---|---|---|");
    let mut prev: Option<f64> = None;
    for h in [1.6, 0.8, 0.4, 0.2] {
        let mut s = fixed_step(h);
        s.dense = true;
        let mut int = Dop853::new(&Oscillator, 0.0, &[1.0, 0.0], s);
        int.step(&Oscillator, 100.0, f64::INFINITY).unwrap();
        let d = int.dense();
        let mut err: f64 = 0.0;
        for i in 1..200 {
            let t = h * f64::from(i) / 200.0;
            let e0 = (d.eval_component(0, t) - t.cos()).abs();
            let e1 = (d.eval_component(1, t) + t.sin()).abs();
            err = err.max(e0).max(e1);
        }
        let order = prev.map_or("-".to_string(), |p| format!("{:.2}", (p / err).log2()));
        println!("| {h} | {err:.3e} | {order} |");
        prev = Some(err);
    }
    println!();
}

// ---------------------------------------------------------------------------------------
// T6 prototype: relativistic motion from rest in a uniform field, units m = c = qE = 1.
// State (x, p); dx/dt = p / sqrt(1 + p^2), dp/dt = 1. Exact: p = t, x = sqrt(1 + t^2) - 1.
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

fn hyperbolic_table() {
    println!("## D. Relativistic hyperbolic motion (uniform E, from rest), rtol = atol = tol\n");
    println!("| tol | t_end | gamma(t_end) | relative error in x | steps |");
    println!("|---|---|---|---|---|");
    for tol in [1e-8, 1e-12] {
        for t_end in [1.0, 10.0, 1000.0] {
            let s = Settings::with_tolerance(tol, tol);
            let r = ours(&Hyperbolic, &[0.0, 0.0], t_end, s, 1);
            let x_exact = (1.0 + t_end * t_end).sqrt() - 1.0;
            let gamma = (1.0 + t_end * t_end).sqrt();
            println!(
                "| {tol:.0e} | {t_end} | {gamma:.1} | {:.2e} | {} |",
                (r.y[0] - x_exact).abs() / x_exact,
                r.steps
            );
        }
    }
    println!();
}

fn main() {
    arenstorf_table();
    kepler_order_table();
    dense_order_table();
    hyperbolic_table();
}
