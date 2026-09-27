//! Dormand–Prince 8(5,3) explicit Runge–Kutta method with 7th-order dense output.
//!
//! A step-at-a-time port of `DP86CO` from E. Hairer's `dop853.f` (Hairer, Nørsett, Wanner,
//! *Solving Ordinary Differential Equations I*, 2nd ed., §II.10). Stages, error estimator,
//! step-size controller, initial step guess and dense output follow the Fortran code
//! operation by operation. Differences from the Fortran driver:
//!
//! - The caller advances one accepted step at a time and gets the dense-output polynomial
//!   of that step, so events can be located on the interpolant.
//! - Each step accepts an extra upper bound on its size (`h_cap`), used by event safeguards.
//! - A step that reaches `t_limit` ends exactly at `t_limit`.
//! - Only forward integration (`t` increasing) is supported.
//! - The stiffness detector is omitted (it only issues a diagnostic in the original).

// Index loops mirror the Fortran operation order, which we reproduce deliberately.
#![allow(clippy::needless_range_loop)]

use super::dop853_coefficients::{A, B, BHH, C, D, E};

/// A first-order system `y' = f(t, y)`.
pub trait OdeSystem {
    fn dim(&self) -> usize;
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]);
}

/// Integrator parameters. Defaults are those of `dop853.f`.
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    pub rtol: f64,
    pub atol: f64,
    /// Maximum step size.
    pub h_max: f64,
    /// Initial step size; `None` uses Hairer's `HINIT` estimate.
    pub h_init: Option<f64>,
    pub safety: f64,
    /// Lower bound of `h_new / h_old`.
    pub fac_min: f64,
    /// Upper bound of `h_new / h_old`.
    pub fac_max: f64,
    /// Lund stabilization parameter (0 disables it).
    pub beta: f64,
    pub max_steps: u64,
    pub uround: f64,
    /// Compute the dense-output polynomial on every accepted step (3 extra evaluations).
    pub dense: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            rtol: 1e-6,
            atol: 1e-6,
            h_max: f64::INFINITY,
            h_init: None,
            safety: 0.9,
            fac_min: 0.333,
            fac_max: 6.0,
            beta: 0.0,
            max_steps: 100_000,
            uround: 2.3e-16,
            dense: true,
        }
    }
}

impl Settings {
    pub fn with_tolerance(rtol: f64, atol: f64) -> Self {
        Self {
            rtol,
            atol,
            ..Self::default()
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    StepSizeTooSmall,
    MaxStepsReached,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::StepSizeTooSmall => write!(f, "step size too small"),
            Error::MaxStepsReached => write!(f, "maximum number of steps reached"),
        }
    }
}

impl std::error::Error for Error {}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub n_fcn: u64,
    pub n_step: u64,
    pub n_accept: u64,
    pub n_reject: u64,
}

/// Dense-output polynomial of one accepted step, valid on `[t_start, t_end]`.
#[derive(Clone, Debug)]
pub struct Dense {
    t_old: f64,
    h: f64,
    n: usize,
    /// Eight coefficient blocks of length `n` (Hairer's `CONT`).
    cont: Vec<f64>,
}

impl Dense {
    fn new(n: usize) -> Self {
        Self {
            t_old: 0.0,
            h: 0.0,
            n,
            cont: vec![0.0; 8 * n],
        }
    }

    pub fn t_start(&self) -> f64 {
        self.t_old
    }

    pub fn t_end(&self) -> f64 {
        self.t_old + self.h
    }

    /// Component `i` of the interpolant at time `t` (Hairer's `CONTD8`).
    pub fn eval_component(&self, i: usize, t: f64) -> f64 {
        let n = self.n;
        let c = |r: usize| self.cont[r * n + i];
        let s = (t - self.t_old) / self.h;
        let s1 = 1.0 - s;
        let conpar = c(4) + s * (c(5) + s1 * (c(6) + s * c(7)));
        c(0) + s * (c(1) + s1 * (c(2) + s * (c(3) + s1 * conpar)))
    }

    /// Time derivative of component `i` of the interpolant at `t` (the derivative of the
    /// polynomial of `eval_component`, one order less accurate).
    pub fn eval_derivative_component(&self, i: usize, t: f64) -> f64 {
        let n = self.n;
        let c = |r: usize| self.cont[r * n + i];
        let s = (t - self.t_old) / self.h;
        let s1 = 1.0 - s;
        // Values and d/ds of the nested factors, innermost first.
        let q3 = c(6) + s * c(7);
        let dq3 = c(7);
        let q2 = c(5) + s1 * q3;
        let dq2 = -q3 + s1 * dq3;
        let q = c(4) + s * q2;
        let dq = q2 + s * dq2;
        let p4 = c(3) + s1 * q;
        let dp4 = -q + s1 * dq;
        let p3 = c(2) + s * p4;
        let dp3 = p4 + s * dp4;
        let p2 = c(1) + s1 * p3;
        let dp2 = -p3 + s1 * dp3;
        (p2 + s * dp2) / self.h
    }

    pub fn eval(&self, t: f64, out: &mut [f64]) {
        for (i, o) in out.iter_mut().enumerate() {
            *o = self.eval_component(i, t);
        }
    }
}

pub struct Dop853 {
    settings: Settings,
    n: usize,
    t: f64,
    y: Vec<f64>,
    /// Stage derivatives `k[1..=16]` (index 0 unused). `k[1]` is `f(t, y)`.
    k: Vec<Vec<f64>>,
    y1: Vec<f64>,
    y_new: Vec<f64>,
    incr: Vec<f64>,
    h: f64,
    facold: f64,
    reject: bool,
    dense: Dense,
    stats: Stats,
}

impl Dop853 {
    pub fn new(sys: &impl OdeSystem, t0: f64, y0: &[f64], settings: Settings) -> Self {
        let n = sys.dim();
        assert_eq!(y0.len(), n, "state dimension mismatch");
        let mut s = Self {
            n,
            t: t0,
            y: y0.to_vec(),
            k: vec![vec![0.0; n]; 17],
            y1: vec![0.0; n],
            y_new: vec![0.0; n],
            incr: vec![0.0; n],
            h: 0.0,
            facold: 1e-4,
            reject: false,
            dense: Dense::new(n),
            stats: Stats::default(),
            settings,
        };
        sys.rhs(s.t, &s.y, &mut s.k[1]);
        s.stats.n_fcn += 1;
        s.h = match s.settings.h_init {
            Some(h) => h,
            None => {
                s.stats.n_fcn += 1;
                s.hinit(sys)
            }
        };
        s
    }

    pub fn t(&self) -> f64 {
        self.t
    }

    pub fn y(&self) -> &[f64] {
        &self.y
    }

    /// Derivative `f(t, y)` at the current point.
    pub fn dy(&self) -> &[f64] {
        &self.k[1]
    }

    /// Step size proposed for the next step.
    pub fn h_next(&self) -> f64 {
        self.h
    }

    pub fn stats(&self) -> Stats {
        self.stats
    }

    /// Dense output of the last accepted step (only meaningful with `Settings::dense`).
    pub fn dense(&self) -> &Dense {
        &self.dense
    }

    /// Advances by one accepted step, never past `t_limit` and never longer than `h_cap`.
    /// Returns `true` when `t_limit` was reached.
    pub fn step(&mut self, sys: &impl OdeSystem, t_limit: f64, h_cap: f64) -> Result<bool, Error> {
        let n = self.n;
        let st = self.settings;
        let expo1 = 1.0 / 8.0 - st.beta * 0.2;
        let facc1 = 1.0 / st.fac_min;
        let facc2 = 1.0 / st.fac_max;
        #[allow(clippy::cast_precision_loss)] // n is a small state dimension
        let n_f = n as f64;

        loop {
            if self.stats.n_step > st.max_steps {
                return Err(Error::MaxStepsReached);
            }
            let mut h = self.h.min(h_cap);
            if 0.1 * h.abs() <= self.t.abs() * st.uround {
                return Err(Error::StepSizeTooSmall);
            }
            let mut last = false;
            if self.t + 1.01 * h - t_limit > 0.0 {
                // Stretch or shrink the step to land exactly on t_limit, but never beyond
                // the cap: `h_cap` is a hard bound used by event safeguards.
                let remaining = t_limit - self.t;
                if remaining <= h_cap {
                    h = remaining;
                    last = true;
                }
            }
            self.stats.n_step += 1;

            // The twelve stages.
            for s in 2..=12 {
                stage(&mut self.y1, &self.y, &self.k, &A[s], s, h);
                sys.rhs(self.t + C[s] * h, &self.y1, &mut self.k[s]);
            }
            self.stats.n_fcn += 11;

            // 8th-order increment and error estimate.
            let k = &self.k;
            let mut err = 0.0;
            let mut err2 = 0.0;
            for i in 0..n {
                let mut acc = 0.0;
                for j in [1, 6, 7, 8, 9, 10, 11, 12] {
                    acc += B[j] * k[j][i];
                }
                self.incr[i] = acc;
                self.y_new[i] = self.y[i] + h * acc;
                let sk = st.atol + st.rtol * self.y[i].abs().max(self.y_new[i].abs());
                let erri = acc - BHH[1] * k[1][i] - BHH[2] * k[9][i] - BHH[3] * k[12][i];
                err2 += (erri / sk) * (erri / sk);
                let mut erri = 0.0;
                for j in [1, 6, 7, 8, 9, 10, 11, 12] {
                    erri += E[j] * k[j][i];
                }
                err += (erri / sk) * (erri / sk);
            }
            let mut deno = err + 0.01 * err2;
            if deno <= 0.0 {
                deno = 1.0;
            }
            let err = h.abs() * err * (1.0 / (n_f * deno)).sqrt();

            // New step size.
            let fac11 = err.powf(expo1);
            let fac = fac11 / self.facold.powf(st.beta);
            let fac = facc2.max(facc1.min(fac / st.safety));
            let mut h_new = h / fac;

            if err <= 1.0 {
                // Accepted.
                self.facold = err.max(1e-4);
                self.stats.n_accept += 1;
                let t_new = if last { t_limit } else { self.t + h };
                sys.rhs(self.t + h, &self.y_new, &mut self.k[13]);
                self.stats.n_fcn += 1;

                if st.dense {
                    self.prepare_dense(sys, h);
                }

                self.k.swap(1, 13);
                std::mem::swap(&mut self.y, &mut self.y_new);
                self.t = t_new;

                if !last {
                    if h_new.abs() > st.h_max {
                        h_new = st.h_max;
                    }
                    if self.reject {
                        h_new = h_new.abs().min(h.abs());
                    }
                }
                self.reject = false;
                self.h = h_new;
                return Ok(last);
            }
            // Rejected.
            h_new = h / facc1.min(fac11 / st.safety);
            self.reject = true;
            if self.stats.n_accept >= 1 {
                self.stats.n_reject += 1;
            }
            self.h = h_new;
        }
    }

    /// Builds the dense-output coefficients of the step just accepted (before `y`, `k[1]`
    /// are overwritten). Requires `k[13] = f(t + h, y_new)`.
    fn prepare_dense(&mut self, sys: &impl OdeSystem, h: f64) {
        let n = self.n;
        // The next three function evaluations (stages 14, 15, 16).
        for s in 14..=16 {
            stage(&mut self.y1, &self.y, &self.k, &A[s], s, h);
            sys.rhs(self.t + C[s] * h, &self.y1, &mut self.k[s]);
        }
        self.stats.n_fcn += 3;

        let k = &self.k;
        let cont = &mut self.dense.cont;
        for i in 0..n {
            let ydiff = self.y_new[i] - self.y[i];
            let bspl = h * k[1][i] - ydiff;
            cont[i] = self.y[i];
            cont[n + i] = ydiff;
            cont[2 * n + i] = bspl;
            cont[3 * n + i] = ydiff - h * k[13][i] - bspl;
            for r in 4..=7 {
                let mut acc = 0.0;
                for j in [1, 6, 7, 8, 9, 10, 11, 12] {
                    acc += D[r][j] * k[j][i];
                }
                acc = acc + D[r][13] * k[13][i] + D[r][14] * k[14][i] + D[r][15] * k[15][i];
                acc += D[r][16] * k[16][i];
                cont[r * n + i] = h * acc;
            }
        }
        self.dense.t_old = self.t;
        self.dense.h = h;
    }

    /// Initial step size guess (Hairer's `HINIT` with `IORD = 8`).
    fn hinit(&mut self, sys: &impl OdeSystem) -> f64 {
        let st = &self.settings;
        let f0 = &self.k[1];
        let mut dnf = 0.0;
        let mut dny = 0.0;
        for i in 0..self.n {
            let sk = st.atol + st.rtol * self.y[i].abs();
            dnf += (f0[i] / sk) * (f0[i] / sk);
            dny += (self.y[i] / sk) * (self.y[i] / sk);
        }
        let mut h = if dnf <= 1e-10 || dny <= 1e-10 {
            1e-6
        } else {
            (dny / dnf).sqrt() * 0.01
        };
        h = h.min(st.h_max);
        for i in 0..self.n {
            self.y1[i] = self.y[i] + h * f0[i];
        }
        let mut f1 = vec![0.0; self.n];
        sys.rhs(self.t + h, &self.y1, &mut f1);
        let mut der2 = 0.0;
        for i in 0..self.n {
            let sk = st.atol + st.rtol * self.y[i].abs();
            der2 += ((f1[i] - f0[i]) / sk) * ((f1[i] - f0[i]) / sk);
        }
        der2 = der2.sqrt() / h;
        let der12 = der2.abs().max(dnf.sqrt());
        let h1 = if der12 <= 1e-15 {
            1e-6_f64.max(h.abs() * 1e-3)
        } else {
            (0.01 / der12).powf(1.0 / 8.0)
        };
        (100.0 * h.abs()).min(h1).min(st.h_max)
    }
}

/// `y1 = y + h * Σ_{j<s} a[j] k[j]`, summed in ascending `j`, skipping zero coefficients
/// (the same operation order as the Fortran code).
fn stage(y1: &mut [f64], y: &[f64], k: &[Vec<f64>], a: &[f64; 17], s: usize, h: f64) {
    for i in 0..y.len() {
        let mut acc = 0.0;
        for j in 1..s {
            if a[j] != 0.0 {
                acc += a[j] * k[j][i];
            }
        }
        y1[i] = y[i] + h * acc;
    }
}

#[cfg(test)]
mod tests {
    //! Checks of the generated coefficient table against Runge–Kutta order conditions
    //! (independent of how the table was produced).
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 2e-15 * (1.0 + b.abs())
    }

    #[test]
    fn row_sums_equal_nodes() {
        for s in (2..=12).chain(14..=16) {
            let sum: f64 = A[s].iter().sum();
            // Rows contain large coefficients that cancel; the rounding error of the
            // summation scales with sum |a|, not with the result.
            let scale: f64 = A[s].iter().map(|a| a.abs()).sum();
            let err = (sum - C[s]).abs();
            assert!(
                err <= 4.0 * f64::EPSILON * scale,
                "stage {s}: sum a = {sum}, c = {}",
                C[s]
            );
        }
    }

    #[test]
    fn weights_satisfy_quadrature_conditions_to_order_8() {
        for q in 1..=8 {
            let sum: f64 = (1..=12).map(|i| B[i] * C[i].powi(q - 1)).sum();
            assert!(close(sum, 1.0 / f64::from(q)), "q = {q}: {sum}");
        }
    }

    #[test]
    fn third_order_estimator_weights_are_exact_to_order_3() {
        let c = [C[1], C[9], C[12]];
        for q in 1..=3 {
            let sum: f64 = (0..3).map(|i| BHH[i + 1] * c[i].powi(q - 1)).sum();
            assert!(close(sum, 1.0 / f64::from(q)), "q = {q}: {sum}");
        }
    }

    /// The derivative of the dense output matches `f` at both ends of every step (the
    /// interpolant is built from them) and the exact derivative inside, for the harmonic
    /// oscillator `y = (cos t, −sin t)`.
    #[test]
    fn dense_derivative_matches_the_exact_derivative() {
        struct Osc;
        impl OdeSystem for Osc {
            fn dim(&self) -> usize {
                2
            }
            fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
                dy[0] = y[1];
                dy[1] = -y[0];
            }
        }
        let settings = Settings {
            rtol: 1e-12,
            atol: 1e-12,
            dense: true,
            ..Settings::default()
        };
        let mut int = Dop853::new(&Osc, 0.0, &[1.0, 0.0], settings);
        let mut worst: f64 = 0.0;
        let mut ends: f64 = 0.0;
        while !int.step(&Osc, 10.0, f64::INFINITY).unwrap() {
            let d = int.dense();
            for k in 0..=8 {
                let t = d.t_start() + (d.t_end() - d.t_start()) * f64::from(k) / 8.0;
                worst = worst
                    .max((d.eval_derivative_component(0, t) + t.sin()).abs())
                    .max((d.eval_derivative_component(1, t) + t.cos()).abs());
            }
            let t = d.t_end();
            ends = ends
                .max((d.eval_derivative_component(0, t) - int.dy()[0]).abs())
                .max((d.eval_derivative_component(1, t) - int.dy()[1]).abs());
        }
        assert!(worst < 1e-9, "inside: {worst:e}");
        assert!(ends < 1e-12, "at the ends: {ends:e}");
    }

    #[test]
    fn fifth_order_estimator_annihilates_polynomials_to_degree_4() {
        for q in 1..=5 {
            let sum: f64 = (1..=12).map(|i| E[i] * C[i].powi(q - 1)).sum();
            assert!(sum.abs() < 1e-14, "q = {q}: {sum}");
        }
    }
}
