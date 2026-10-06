//! RADAU5 (Hairer & Wanner, *Solving Ordinary Differential Equations II*, §IV.8): the
//! three-stage Radau IIA implicit Runge–Kutta method of order 5, for stiff systems and
//! index-1 differential-algebraic equations `M y' = f(t, y)`, with step-size control
//! (Gustafsson's predictive controller) and the collocation polynomial as dense output.
//!
//! A port of `radau5.f` (version of July 9, 1996, latest correction January 18, 2002) with
//! the dense linear algebra of `decsol.f` and `dc_decsol.f`, for full matrices only (no
//! banded Jacobians, no second-order structure, no Hessenberg form). The control flow of
//! RADCOR is kept, step for step: the simplified Newton iteration on the stages transformed
//! by `T⁻¹` into one real and one complex linear system (`(γ/h) M − J` and
//! `((α + iβ)/h) M − J`), its convergence rate θ deciding when the Jacobian and the
//! decompositions are recomputed, the embedded error estimate (ESTRAD, with its second
//! solve after a rejected or first step), and the tolerances transformed as RADAU5 does
//! (`rtol' = 0.1 rtol^{2/3}`, `atol' = rtol' atol/rtol`). Variables of index 2 and 3 are
//! supported as in RADAU5 (their error scaled by h and h²).
//!
//! Elementary functions come from `libm` (deterministic across platforms, CLAUDE.md).

/// A stiff system `M y' = f(t, y)`.
pub trait StiffSystem {
    fn dim(&self) -> usize;
    fn rhs(&self, t: f64, y: &[f64], f: &mut [f64]);
    /// The constant mass matrix `M`, column-major (`m[i + n j]` = row i, column j); None:
    /// the identity.
    fn mass(&self) -> Option<Vec<f64>> {
        None
    }
    /// The Jacobian `∂f/∂y` at `(t, y)`, column-major, if known analytically (returns
    /// whether it filled `jac`); otherwise it is taken by forward differences.
    fn jacobian(&self, _t: f64, _y: &[f64], _jac: &mut [f64]) -> bool {
        false
    }
}

/// Why the integration stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// The iteration matrix was singular five times in a row (RADAU5's IDID = −4).
    Singular,
    /// The step size fell below the rounding of t (IDID = −3).
    StepTooSmall,
    /// More than `n_max` steps (IDID = −2).
    MaxSteps,
}

/// Settings, with RADAU5's defaults (`Settings::new`).
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// The user's tolerances (transformed internally as RADAU5 does).
    pub rtol: f64,
    pub atol: f64,
    /// Initial step (RADAU5 starts from 1e-6 if 0).
    pub h0: f64,
    /// Largest step (0: the whole interval).
    pub h_max: f64,
    pub n_max: usize,
    /// Most Newton iterations per step.
    pub nit: usize,
    /// Safety factor of the step-size prediction.
    pub safe: f64,
    /// The Jacobian is recomputed when the Newton contraction θ exceeds this.
    pub thet: f64,
    /// Newton's stopping criterion (None: RADAU5's `max(10 u/rtol', min(0.03, √rtol'))`).
    pub fnewt: Option<f64>,
    /// The step is kept (and the decompositions reused) while `quot1 ≤ h_new/h ≤ quot2`.
    pub quot1: f64,
    pub quot2: f64,
    /// The step changes by at most these factors: `1/facr ≥ h_old/h_new ≥ 1/facl`.
    pub facl: f64,
    pub facr: f64,
    /// Gustafsson's predictive step-size controller (else the classical one).
    pub pred: bool,
    /// Start every Newton iteration from zero (else from the extrapolated collocation
    /// polynomial).
    pub startn: bool,
    /// Numbers of variables of index 1, 2 and 3, in this order (default: all index 1).
    pub index: Option<(usize, usize, usize)>,
}

impl Settings {
    pub fn new(rtol: f64, atol: f64) -> Self {
        Self {
            rtol,
            atol,
            h0: 0.0,
            h_max: 0.0,
            n_max: 100_000,
            nit: 7,
            safe: 0.9,
            thet: 0.001,
            fnewt: None,
            quot1: 1.0,
            quot2: 1.2,
            facl: 5.0,
            facr: 1.0 / 8.0,
            pred: true,
            startn: false,
            index: None,
        }
    }
}

/// Counters (RADAU5's IWORK(14..20); `n_fcn` without the evaluations for a numerical
/// Jacobian, as there).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub n_fcn: u64,
    pub n_jac: u64,
    pub n_step: u64,
    pub n_accept: u64,
    pub n_reject: u64,
    pub n_dec: u64,
    pub n_sol: u64,
}

/// RADAU5's rounding unit (its default WORK(1)).
const UROUND: f64 = 1e-16;

/// The collocation polynomial of the last accepted step (RADAU5's CONTR5).
#[derive(Clone, Debug)]
pub struct Dense {
    n: usize,
    t_old: f64,
    t: f64,
    c1m1: f64,
    c2m1: f64,
    cont: Vec<f64>,
}

impl Dense {
    pub fn t_start(&self) -> f64 {
        self.t_old
    }

    pub fn t_end(&self) -> f64 {
        self.t
    }

    /// Component `i` at time `x` (within the step).
    pub fn eval_component(&self, i: usize, x: f64) -> f64 {
        let n = self.n;
        let s = (x - self.t) / (self.t - self.t_old);
        let c = &self.cont;
        c[i] + s * (c[i + n] + (s - self.c2m1) * (c[i + 2 * n] + (s - self.c1m1) * c[i + 3 * n]))
    }

    /// The whole state at `x`.
    pub fn eval(&self, x: f64, y: &mut [f64]) {
        for (i, yi) in y.iter_mut().enumerate() {
            *yi = self.eval_component(i, x);
        }
    }
}

/// What the next call has to do first (RADCOR's labels 10, 20 and 30).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Next {
    Jacobian,
    Decompose,
    Step,
}

/// The RADAU5 integrator: `new`, then `step` until `t()` reaches the end.
pub struct Radau5 {
    n: usize,
    t: f64,
    t_end: f64,
    y: Vec<f64>,
    h: f64,
    hold: f64,
    hopt: f64,
    hacc: f64,
    erracc: f64,
    hmaxn: f64,
    posneg: f64,
    rtol: f64,
    atol: f64,
    fnewt: f64,
    cfac: f64,
    s: Settings,
    nind: (usize, usize, usize),
    // The method's constants.
    c1: f64,
    c2: f64,
    c1m1: f64,
    c2m1: f64,
    c1mc2: f64,
    dd: [f64; 3],
    u1: f64,
    alph: f64,
    beta: f64,
    tm: [[f64; 3]; 3],
    ti: [[f64; 3]; 3],
    // Work arrays.
    z1: Vec<f64>,
    z2: Vec<f64>,
    z3: Vec<f64>,
    y0: Vec<f64>,
    scal: Vec<f64>,
    f1: Vec<f64>,
    f2: Vec<f64>,
    f3: Vec<f64>,
    cont: Vec<f64>,
    fjac: Vec<f64>,
    fmas: Option<Vec<f64>>,
    e1: Vec<f64>,
    e2r: Vec<f64>,
    e2i: Vec<f64>,
    ip1: Vec<usize>,
    ip2: Vec<usize>,
    // Control.
    next: Next,
    first: bool,
    reject: bool,
    last: bool,
    caljac: bool,
    faccon: f64,
    theta: f64,
    hhfac: f64,
    nsing: usize,
    started: bool,
    stats: Stats,
}

impl Radau5 {
    /// An integrator for `sys` from `(t0, y0)` towards `t_end`.
    pub fn new(
        sys: &impl StiffSystem,
        t0: f64,
        y0: &[f64],
        t_end: f64,
        settings: Settings,
    ) -> Self {
        let n = sys.dim();
        assert_eq!(y0.len(), n);
        assert!(
            settings.atol > 0.0 && settings.rtol > 10.0 * UROUND,
            "tolerances are too small"
        );
        // The tolerances as RADAU5 transforms them.
        let expm = 2.0 / 3.0;
        let quot = settings.atol / settings.rtol;
        let rtol = 0.1 * libm::pow(settings.rtol, expm);
        let atol = rtol * quot;
        let fnewt = settings
            .fnewt
            .unwrap_or_else(|| (10.0 * UROUND / rtol).max(0.03_f64.min(rtol.sqrt())));
        let nind = settings.index.unwrap_or((n, 0, 0));
        assert_eq!(nind.0 + nind.1 + nind.2, n, "index counts must add up to n");
        let sq6 = 6.0_f64.sqrt();
        let c1 = (4.0 - sq6) / 10.0;
        let c2 = (4.0 + sq6) / 10.0;
        let cbrt81 = libm::cbrt(81.0);
        let cbrt9 = libm::cbrt(9.0);
        let u1 = (6.0 + cbrt81 - cbrt9) / 30.0;
        let alph = (12.0 - cbrt81 + cbrt9) / 60.0;
        let beta = (cbrt81 + cbrt9) * 3.0_f64.sqrt() / 60.0;
        let cno = alph * alph + beta * beta;
        #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
        let tm = [
            [
                9.1232394870892942792e-02,
                -0.14125529502095420843,
                -3.0029194105147424492e-02,
            ],
            [
                0.24171793270710701896,
                0.20412935229379993199,
                0.38294211275726193779,
            ],
            [0.96604818261509293619, 1.0, 0.0],
        ];
        #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
        let ti = [
            [
                4.3255798900631553510,
                0.33919925181580986954,
                0.54177053993587487119,
            ],
            [
                -4.1787185915519047273,
                -0.32768282076106238708,
                0.47662355450055045196,
            ],
            [
                -0.50287263494578687595,
                2.5719269498556054292,
                -0.59603920482822492497,
            ],
        ];
        let posneg = if t_end >= t0 { 1.0 } else { -1.0 };
        let span = (t_end - t0).abs();
        let hmaxn = if settings.h_max > 0.0 {
            settings.h_max.abs().min(span)
        } else {
            span
        };
        let mut h = settings.h0;
        if h.abs() <= 10.0 * UROUND {
            h = 1.0e-6;
        }
        h = h.abs().min(hmaxn) * posneg;
        let mut last = false;
        if (t0 + h * 1.0001 - t_end) * posneg >= 0.0 {
            h = t_end - t0;
            last = true;
        }
        let scal = y0.iter().map(|v| atol + rtol * v.abs()).collect();
        let mut y0v = vec![0.0; n];
        sys.rhs(t0, y0, &mut y0v);
        let mut cont = vec![0.0; 4 * n];
        cont[..n].copy_from_slice(y0);
        Self {
            n,
            t: t0,
            t_end,
            y: y0.to_vec(),
            h,
            hold: h,
            hopt: h,
            hacc: 0.0,
            erracc: 0.0,
            hmaxn,
            posneg,
            rtol,
            atol,
            fnewt,
            #[allow(clippy::cast_precision_loss)]
            cfac: settings.safe * (1 + 2 * settings.nit) as f64,
            s: settings,
            nind,
            c1,
            c2,
            c1m1: c1 - 1.0,
            c2m1: c2 - 1.0,
            c1mc2: c1 - c2,
            dd: [
                -(13.0 + 7.0 * sq6) / 3.0,
                (-13.0 + 7.0 * sq6) / 3.0,
                -1.0 / 3.0,
            ],
            u1: 1.0 / u1,
            alph: alph / cno,
            beta: beta / cno,
            tm,
            ti,
            z1: vec![0.0; n],
            z2: vec![0.0; n],
            z3: vec![0.0; n],
            y0: y0v,
            scal,
            f1: vec![0.0; n],
            f2: vec![0.0; n],
            f3: vec![0.0; n],
            cont,
            fjac: vec![0.0; n * n],
            fmas: sys.mass(),
            e1: vec![0.0; n * n],
            e2r: vec![0.0; n * n],
            e2i: vec![0.0; n * n],
            ip1: vec![0; n],
            ip2: vec![0; n],
            next: Next::Jacobian,
            first: true,
            reject: false,
            last,
            caljac: false,
            faccon: 1.0,
            theta: 0.0,
            hhfac: h,
            nsing: 0,
            started: false,
            stats: Stats {
                n_fcn: 1,
                ..Stats::default()
            },
        }
    }

    pub fn t(&self) -> f64 {
        self.t
    }

    pub fn y(&self) -> &[f64] {
        &self.y
    }

    /// The step size the next step will try.
    pub fn h(&self) -> f64 {
        self.h
    }

    pub fn stats(&self) -> Stats {
        self.stats
    }

    /// Whether the end has been reached.
    pub fn done(&self) -> bool {
        self.started && (self.t - self.t_end) * self.posneg >= 0.0
    }

    /// The collocation polynomial of the last accepted step.
    pub fn dense(&self) -> Dense {
        Dense {
            n: self.n,
            t_old: self.t - self.hold,
            t: self.t,
            c1m1: self.c1m1,
            c2m1: self.c2m1,
            cont: self.cont.clone(),
        }
    }

    /// Column-major index.
    fn ix(&self, i: usize, j: usize) -> usize {
        i + self.n * j
    }

    fn compute_jacobian(&mut self, sys: &impl StiffSystem) {
        self.stats.n_jac += 1;
        let n = self.n;
        if !sys.jacobian(self.t, &self.y, &mut self.fjac) {
            let mut f = vec![0.0; n];
            for i in 0..n {
                let ysafe = self.y[i];
                let delt = (UROUND * 1e-5_f64.max(ysafe.abs())).sqrt();
                self.y[i] = ysafe + delt;
                // Not counted in n_fcn (as RADAU5).
                sys.rhs(self.t, &self.y, &mut f);
                for (j, fj) in f.iter().enumerate() {
                    let k = self.ix(j, i);
                    self.fjac[k] = (fj - self.y0[j]) / delt;
                }
                self.y[i] = ysafe;
            }
        }
        self.caljac = true;
    }

    /// E1 = (γ/h) M − J and E2 = ((α + iβ)/h) M − J, decomposed (DECOMR, DECOMC). Err if
    /// one is singular.
    fn decompose(&mut self) -> Result<(), ()> {
        let n = self.n;
        let fac1 = self.u1 / self.h;
        let alphn = self.alph / self.h;
        let betan = self.beta / self.h;
        for j in 0..n {
            for i in 0..n {
                let k = self.ix(i, j);
                let m = match &self.fmas {
                    Some(fm) => fm[k],
                    None => f64::from(u8::from(i == j)),
                };
                self.e1[k] = m * fac1 - self.fjac[k];
                self.e2r[k] = m * alphn - self.fjac[k];
                self.e2i[k] = m * betan;
            }
        }
        dec(n, &mut self.e1, &mut self.ip1).map_err(|_| ())?;
        decc(n, &mut self.e2r, &mut self.e2i, &mut self.ip2).map_err(|_| ())?;
        self.stats.n_dec += 1;
        Ok(())
    }

    /// One accepted step (rejected attempts are retried inside).
    pub fn step(&mut self, sys: &impl StiffSystem) -> Result<(), Error> {
        let n = self.n;
        if self.done() {
            return Ok(());
        }
        self.started = true;
        loop {
            match self.next {
                Next::Jacobian => {
                    self.compute_jacobian(sys);
                    self.next = Next::Decompose;
                }
                Next::Decompose => {
                    if self.decompose().is_err() {
                        self.nsing += 1;
                        if self.nsing >= 5 {
                            return Err(Error::Singular);
                        }
                        self.unexpected_rejection();
                        continue;
                    }
                    self.next = Next::Step;
                }
                Next::Step => {
                    self.stats.n_step += 1;
                    if self.stats.n_step > self.s.n_max as u64 {
                        return Err(Error::MaxSteps);
                    }
                    if 0.1 * self.h.abs() <= self.t.abs() * UROUND {
                        return Err(Error::StepTooSmall);
                    }
                    if self.attempt(sys) {
                        return Ok(());
                    }
                }
            }
            let _ = n;
        }
    }

    /// RADCOR's label 78: Newton failed or the matrix was singular.
    fn unexpected_rejection(&mut self) {
        self.h *= 0.5;
        self.hhfac = 0.5;
        self.reject = true;
        self.last = false;
        self.next = if self.caljac {
            Next::Decompose
        } else {
            Next::Jacobian
        };
    }

    /// One attempt at a step of size h from the current point (RADCOR from label 30):
    /// true if accepted.
    #[allow(clippy::too_many_lines)]
    fn attempt(&mut self, sys: &impl StiffSystem) -> bool {
        let n = self.n;
        let (n1, n2, n3) = self.nind;
        if n2 > 0 {
            for v in &mut self.scal[n1..n1 + n2] {
                *v /= self.hhfac;
            }
        }
        if n3 > 0 {
            for v in &mut self.scal[n1 + n2..n1 + n2 + n3] {
                *v /= self.hhfac * self.hhfac;
            }
        }
        let h = self.h;
        let x = self.t;
        let xph = x + h;
        let (tm, ti) = (self.tm, self.ti);
        // Starting values for the Newton iteration.
        if self.first || self.s.startn {
            for v in [
                &mut self.z1,
                &mut self.z2,
                &mut self.z3,
                &mut self.f1,
                &mut self.f2,
                &mut self.f3,
            ] {
                v.fill(0.0);
            }
        } else {
            let c3q = h / self.hold;
            let c1q = self.c1 * c3q;
            let c2q = self.c2 * c3q;
            for i in 0..n {
                let ak1 = self.cont[i + n];
                let ak2 = self.cont[i + 2 * n];
                let ak3 = self.cont[i + 3 * n];
                let z1i = c1q * (ak1 + (c1q - self.c2m1) * (ak2 + (c1q - self.c1m1) * ak3));
                let z2i = c2q * (ak1 + (c2q - self.c2m1) * (ak2 + (c2q - self.c1m1) * ak3));
                let z3i = c3q * (ak1 + (c3q - self.c2m1) * (ak2 + (c3q - self.c1m1) * ak3));
                self.z1[i] = z1i;
                self.z2[i] = z2i;
                self.z3[i] = z3i;
                self.f1[i] = ti[0][0] * z1i + ti[0][1] * z2i + ti[0][2] * z3i;
                self.f2[i] = ti[1][0] * z1i + ti[1][1] * z2i + ti[1][2] * z3i;
                self.f3[i] = ti[2][0] * z1i + ti[2][1] * z2i + ti[2][2] * z3i;
            }
        }
        // The simplified Newton iteration.
        let nit = self.s.nit;
        let fac1 = self.u1 / h;
        let alphn = self.alph / h;
        let betan = self.beta / h;
        let mut newt = 0;
        self.faccon = libm::pow(self.faccon.max(UROUND), 0.8);
        self.theta = self.s.thet.abs();
        let mut dynold = 0.0;
        let mut thqold = 0.0;
        let mut stage = vec![0.0; n];
        loop {
            if newt >= nit {
                self.unexpected_rejection();
                return false;
            }
            // The right-hand sides at the stages.
            for (zk, ck) in [(0, self.c1), (1, self.c2), (2, 1.0)] {
                let z = match zk {
                    0 => &self.z1,
                    1 => &self.z2,
                    _ => &self.z3,
                };
                for i in 0..n {
                    stage[i] = self.y[i] + z[i];
                }
                let tk = if zk == 2 { xph } else { x + ck * h };
                let out = match zk {
                    0 => &mut self.z1,
                    1 => &mut self.z2,
                    _ => &mut self.z3,
                };
                sys.rhs(tk, &stage, out);
            }
            self.stats.n_fcn += 3;
            // Transformed, and the linear systems (SLVRAD).
            for i in 0..n {
                let (a1, a2, a3) = (self.z1[i], self.z2[i], self.z3[i]);
                self.z1[i] = ti[0][0] * a1 + ti[0][1] * a2 + ti[0][2] * a3;
                self.z2[i] = ti[1][0] * a1 + ti[1][1] * a2 + ti[1][2] * a3;
                self.z3[i] = ti[2][0] * a1 + ti[2][1] * a2 + ti[2][2] * a3;
            }
            match &self.fmas {
                None => {
                    for i in 0..n {
                        let s2 = -self.f2[i];
                        let s3 = -self.f3[i];
                        self.z1[i] -= self.f1[i] * fac1;
                        self.z2[i] += s2 * alphn - s3 * betan;
                        self.z3[i] += s3 * alphn + s2 * betan;
                    }
                }
                Some(fm) => {
                    for i in 0..n {
                        let (mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0);
                        for j in 0..n {
                            let bb = fm[i + n * j];
                            s1 -= bb * self.f1[j];
                            s2 -= bb * self.f2[j];
                            s3 -= bb * self.f3[j];
                        }
                        self.z1[i] += s1 * fac1;
                        self.z2[i] += s2 * alphn - s3 * betan;
                        self.z3[i] += s3 * alphn + s2 * betan;
                    }
                }
            }
            sol(n, &self.e1, &mut self.z1, &self.ip1);
            solc(
                n,
                &self.e2r,
                &self.e2i,
                &mut self.z2,
                &mut self.z3,
                &self.ip2,
            );
            self.stats.n_sol += 1;
            newt += 1;
            let mut dyno = 0.0;
            for i in 0..n {
                let d = self.scal[i];
                dyno +=
                    (self.z1[i] / d).powi(2) + (self.z2[i] / d).powi(2) + (self.z3[i] / d).powi(2);
            }
            #[allow(clippy::cast_precision_loss)]
            let dyno = (dyno / (3 * n) as f64).sqrt();
            // Bad convergence, or too many iterations to come.
            if newt > 1 && newt < nit {
                let thq = dyno / dynold;
                self.theta = if newt == 2 {
                    thq
                } else {
                    (thq * thqold).sqrt()
                };
                thqold = thq;
                if self.theta < 0.99 {
                    self.faccon = self.theta / (1.0 - self.theta);
                    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                    let dyth =
                        self.faccon * dyno * self.theta.powi((nit - 1 - newt) as i32) / self.fnewt;
                    if dyth >= 1.0 {
                        let qnewt = dyth.clamp(1.0e-4, 20.0);
                        #[allow(clippy::cast_precision_loss)]
                        let expo = -1.0 / (4.0 + (nit - 1 - newt) as f64);
                        self.hhfac = 0.8 * libm::pow(qnewt, expo);
                        self.h *= self.hhfac;
                        self.reject = true;
                        self.last = false;
                        self.next = if self.caljac {
                            Next::Decompose
                        } else {
                            Next::Jacobian
                        };
                        return false;
                    }
                } else {
                    self.unexpected_rejection();
                    return false;
                }
            }
            dynold = dyno.max(UROUND);
            for i in 0..n {
                let f1i = self.f1[i] + self.z1[i];
                let f2i = self.f2[i] + self.z2[i];
                let f3i = self.f3[i] + self.z3[i];
                self.f1[i] = f1i;
                self.f2[i] = f2i;
                self.f3[i] = f3i;
                self.z1[i] = tm[0][0] * f1i + tm[0][1] * f2i + tm[0][2] * f3i;
                self.z2[i] = tm[1][0] * f1i + tm[1][1] * f2i + tm[1][2] * f3i;
                self.z3[i] = tm[2][0] * f1i + f2i;
            }
            if self.faccon * dyno <= self.fnewt {
                break;
            }
        }
        // Error estimation (ESTRAD).
        let err = self.estimate(sys, h);
        // The new step size: 0.2 ≤ h_new/h ≤ 8 by default.
        #[allow(clippy::cast_precision_loss)]
        let fac = self.s.safe.min(self.cfac / (newt + 2 * nit) as f64);
        let mut quot = self.s.facr.max(self.s.facl.min(libm::pow(err, 0.25) / fac));
        let mut hnew = h / quot;
        if err < 1.0 {
            // Accepted.
            self.first = false;
            self.stats.n_accept += 1;
            if self.s.pred {
                // Gustafsson's predictive controller.
                if self.stats.n_accept > 1 {
                    let facgus =
                        (self.hacc / h) * libm::pow(err * err / self.erracc, 0.25) / self.s.safe;
                    let facgus = self.s.facr.max(self.s.facl.min(facgus));
                    quot = quot.max(facgus);
                    hnew = h / quot;
                }
                self.hacc = h;
                self.erracc = 1.0e-2_f64.max(err);
            }
            self.hold = h;
            self.t = xph;
            for i in 0..n {
                self.y[i] += self.z3[i];
                let (z1i, z2i) = (self.z1[i], self.z2[i]);
                self.cont[i + n] = (z2i - self.z3[i]) / self.c2m1;
                let ak = (z1i - z2i) / self.c1mc2;
                let acont3 = (ak - z1i / self.c1) / self.c2;
                self.cont[i + 2 * n] = (ak - self.cont[i + n]) / self.c1m1;
                self.cont[i + 3 * n] = self.cont[i + 2 * n] - acont3;
            }
            for i in 0..n {
                self.scal[i] = self.atol + self.rtol * self.y[i].abs();
                self.cont[i] = self.y[i];
            }
            self.caljac = false;
            if self.last {
                self.h = self.hopt;
                return true;
            }
            sys.rhs(self.t, &self.y, &mut self.y0);
            self.stats.n_fcn += 1;
            hnew = self.posneg * hnew.abs().min(self.hmaxn);
            self.hopt = h.min(hnew);
            if self.reject {
                hnew = self.posneg * hnew.abs().min(h.abs());
            }
            self.reject = false;
            if (self.t + hnew / self.s.quot1 - self.t_end) * self.posneg >= 0.0 {
                self.h = self.t_end - self.t;
                self.last = true;
            } else {
                let qt = hnew / h;
                self.hhfac = h;
                if self.theta <= self.s.thet && qt >= self.s.quot1 && qt <= self.s.quot2 {
                    self.next = Next::Step;
                    return true;
                }
                self.h = hnew;
            }
            self.hhfac = self.h;
            self.next = if self.theta <= self.s.thet {
                Next::Decompose
            } else {
                Next::Jacobian
            };
            true
        } else {
            // Rejected.
            self.reject = true;
            self.last = false;
            if self.first {
                self.h *= 0.1;
                self.hhfac = 0.1;
            } else {
                self.hhfac = hnew / h;
                self.h = hnew;
            }
            if self.stats.n_accept >= 1 {
                self.stats.n_reject += 1;
            }
            self.next = if self.caljac {
                Next::Decompose
            } else {
                Next::Jacobian
            };
            false
        }
    }

    /// The embedded error estimate (ESTRAD), scaled: below 1 accepts the step.
    fn estimate(&mut self, sys: &impl StiffSystem, h: f64) -> f64 {
        let n = self.n;
        let hee = [self.dd[0] / h, self.dd[1] / h, self.dd[2] / h];
        let mut cont = vec![0.0; n];
        match &self.fmas {
            None => {
                for (i, c) in cont.iter_mut().enumerate() {
                    self.f2[i] = hee[0] * self.z1[i] + hee[1] * self.z2[i] + hee[2] * self.z3[i];
                    *c = self.f2[i] + self.y0[i];
                }
            }
            Some(fm) => {
                for i in 0..n {
                    self.f1[i] = hee[0] * self.z1[i] + hee[1] * self.z2[i] + hee[2] * self.z3[i];
                }
                for i in 0..n {
                    let mut sum = 0.0;
                    for j in 0..n {
                        sum += fm[i + n * j] * self.f1[j];
                    }
                    self.f2[i] = sum;
                    cont[i] = sum + self.y0[i];
                }
            }
        }
        sol(n, &self.e1, &mut cont, &self.ip1);
        let rms = |c: &[f64], scal: &[f64]| {
            let mut e = 0.0;
            for i in 0..n {
                e += (c[i] / scal[i]).powi(2);
            }
            #[allow(clippy::cast_precision_loss)]
            let e = (e / n as f64).sqrt();
            e.max(1.0e-10)
        };
        let mut err = rms(&cont, &self.scal);
        if err >= 1.0 && (self.first || self.reject) {
            let mut ys = vec![0.0; n];
            for i in 0..n {
                ys[i] = self.y[i] + cont[i];
            }
            sys.rhs(self.t, &ys, &mut self.f1);
            self.stats.n_fcn += 1;
            for (i, c) in cont.iter_mut().enumerate() {
                *c = self.f1[i] + self.f2[i];
            }
            sol(n, &self.e1, &mut cont, &self.ip1);
            err = rms(&cont, &self.scal);
        }
        err
    }
}

/// LU decomposition with partial pivoting of the column-major n×n matrix `a` (DEC):
/// multipliers below the diagonal, negated, and the pivot rows in `ip`. Err(k) if the
/// matrix is singular at stage k.
fn dec(n: usize, a: &mut [f64], ip: &mut [usize]) -> Result<(), usize> {
    let ix = |i: usize, j: usize| i + n * j;
    for k in 0..n.saturating_sub(1) {
        let mut m = k;
        for i in k + 1..n {
            if a[ix(i, k)].abs() > a[ix(m, k)].abs() {
                m = i;
            }
        }
        ip[k] = m;
        let mut t = a[ix(m, k)];
        if m != k {
            a[ix(m, k)] = a[ix(k, k)];
            a[ix(k, k)] = t;
        }
        if t == 0.0 {
            return Err(k);
        }
        t = 1.0 / t;
        for i in k + 1..n {
            a[ix(i, k)] = -a[ix(i, k)] * t;
        }
        for j in k + 1..n {
            let t = a[ix(m, j)];
            a[ix(m, j)] = a[ix(k, j)];
            a[ix(k, j)] = t;
            if t != 0.0 {
                for i in k + 1..n {
                    a[ix(i, j)] += a[ix(i, k)] * t;
                }
            }
        }
    }
    if n > 0 && a[ix(n - 1, n - 1)] == 0.0 {
        return Err(n - 1);
    }
    Ok(())
}

/// Solves `A x = b` with the decomposition of `dec` (SOL); `b` becomes x.
fn sol(n: usize, a: &[f64], b: &mut [f64], ip: &[usize]) {
    let ix = |i: usize, j: usize| i + n * j;
    if n == 0 {
        return;
    }
    for k in 0..n - 1 {
        b.swap(ip[k], k);
        let t = b[k];
        for i in k + 1..n {
            b[i] += a[ix(i, k)] * t;
        }
    }
    for k in (1..n).rev() {
        b[k] /= a[ix(k, k)];
        let t = -b[k];
        for i in 0..k {
            b[i] += a[ix(i, k)] * t;
        }
    }
    b[0] /= a[ix(0, 0)];
}

/// `dec` for the complex matrix `ar + i ai` (DECC).
fn decc(n: usize, ar: &mut [f64], ai: &mut [f64], ip: &mut [usize]) -> Result<(), usize> {
    let ix = |i: usize, j: usize| i + n * j;
    for k in 0..n.saturating_sub(1) {
        let mut m = k;
        for i in k + 1..n {
            if ar[ix(i, k)].abs() + ai[ix(i, k)].abs() > ar[ix(m, k)].abs() + ai[ix(m, k)].abs() {
                m = i;
            }
        }
        ip[k] = m;
        let mut tr = ar[ix(m, k)];
        let mut ti = ai[ix(m, k)];
        if m != k {
            ar[ix(m, k)] = ar[ix(k, k)];
            ai[ix(m, k)] = ai[ix(k, k)];
            ar[ix(k, k)] = tr;
            ai[ix(k, k)] = ti;
        }
        if tr.abs() + ti.abs() == 0.0 {
            return Err(k);
        }
        let den = tr * tr + ti * ti;
        tr /= den;
        ti = -ti / den;
        for i in k + 1..n {
            let prodr = ar[ix(i, k)] * tr - ai[ix(i, k)] * ti;
            let prodi = ai[ix(i, k)] * tr + ar[ix(i, k)] * ti;
            ar[ix(i, k)] = -prodr;
            ai[ix(i, k)] = -prodi;
        }
        for j in k + 1..n {
            let tr = ar[ix(m, j)];
            let ti = ai[ix(m, j)];
            ar[ix(m, j)] = ar[ix(k, j)];
            ai[ix(m, j)] = ai[ix(k, j)];
            ar[ix(k, j)] = tr;
            ai[ix(k, j)] = ti;
            if tr.abs() + ti.abs() == 0.0 {
                continue;
            }
            if ti == 0.0 {
                for i in k + 1..n {
                    let prodr = ar[ix(i, k)] * tr;
                    let prodi = ai[ix(i, k)] * tr;
                    ar[ix(i, j)] += prodr;
                    ai[ix(i, j)] += prodi;
                }
                continue;
            }
            if tr == 0.0 {
                for i in k + 1..n {
                    let prodr = -ai[ix(i, k)] * ti;
                    let prodi = ar[ix(i, k)] * ti;
                    ar[ix(i, j)] += prodr;
                    ai[ix(i, j)] += prodi;
                }
                continue;
            }
            for i in k + 1..n {
                let prodr = ar[ix(i, k)] * tr - ai[ix(i, k)] * ti;
                let prodi = ai[ix(i, k)] * tr + ar[ix(i, k)] * ti;
                ar[ix(i, j)] += prodr;
                ai[ix(i, j)] += prodi;
            }
        }
    }
    if n > 0 && ar[ix(n - 1, n - 1)].abs() + ai[ix(n - 1, n - 1)].abs() == 0.0 {
        return Err(n - 1);
    }
    Ok(())
}

/// `sol` for the complex system (SOLC); `br + i bi` becomes the solution.
fn solc(n: usize, ar: &[f64], ai: &[f64], br: &mut [f64], bi: &mut [f64], ip: &[usize]) {
    let ix = |i: usize, j: usize| i + n * j;
    if n == 0 {
        return;
    }
    for k in 0..n - 1 {
        let m = ip[k];
        let tr = br[m];
        let ti = bi[m];
        br[m] = br[k];
        bi[m] = bi[k];
        br[k] = tr;
        bi[k] = ti;
        for i in k + 1..n {
            let prodr = ar[ix(i, k)] * tr - ai[ix(i, k)] * ti;
            let prodi = ai[ix(i, k)] * tr + ar[ix(i, k)] * ti;
            br[i] += prodr;
            bi[i] += prodi;
        }
    }
    let divide = |br: &mut [f64], bi: &mut [f64], k: usize| {
        let (a, b) = (ar[ix(k, k)], ai[ix(k, k)]);
        let den = a * a + b * b;
        let prodr = br[k] * a + bi[k] * b;
        let prodi = bi[k] * a - br[k] * b;
        br[k] = prodr / den;
        bi[k] = prodi / den;
    };
    for k in (1..n).rev() {
        divide(br, bi, k);
        let tr = -br[k];
        let ti = -bi[k];
        for i in 0..k {
            let prodr = ar[ix(i, k)] * tr - ai[ix(i, k)] * ti;
            let prodi = ai[ix(i, k)] * tr + ar[ix(i, k)] * ti;
            br[i] += prodr;
            bi[i] += prodi;
        }
    }
    divide(br, bi, 0);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The real and complex LU solves against products with known solutions.
    #[test]
    fn linear_solves() {
        let n = 4;
        // Column-major.
        let a0: Vec<f64> = vec![
            2.0, 1.0, -3.0, 0.5, //
            -1.0, 4.0, 2.0, 1.0, //
            0.0, 1.5, -2.0, 3.0, //
            5.0, -1.0, 1.0, -2.0,
        ];
        let x: Vec<f64> = vec![1.0, -2.0, 0.5, 3.0];
        let mut b = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                b[i] += a0[i + n * j] * x[j];
            }
        }
        let mut a = a0.clone();
        let mut ip = vec![0; n];
        dec(n, &mut a, &mut ip).expect("regular");
        sol(n, &a, &mut b, &ip);
        for i in 0..n {
            assert!((b[i] - x[i]).abs() < 1e-13, "{i}: {} vs {}", b[i], x[i]);
        }
        // Complex: (A + i B) z = c.
        let bm: Vec<f64> = (0..n * n)
            .map(|k| 0.3 * f64::from(u32::try_from(k % 5).unwrap()) - 0.4)
            .collect();
        let (zr, zi) = (vec![0.5, 1.0, -1.5, 2.0], vec![-1.0, 0.25, 2.0, 0.0]);
        let (mut cr, mut ci) = (vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            for j in 0..n {
                let (p, q) = (a0[i + n * j], bm[i + n * j]);
                cr[i] += p * zr[j] - q * zi[j];
                ci[i] += p * zi[j] + q * zr[j];
            }
        }
        let (mut ar, mut ai) = (a0.clone(), bm.clone());
        decc(n, &mut ar, &mut ai, &mut ip).expect("regular");
        solc(n, &ar, &ai, &mut cr, &mut ci, &ip);
        for i in 0..n {
            assert!(
                (cr[i] - zr[i]).abs() < 1e-13 && (ci[i] - zi[i]).abs() < 1e-13,
                "{i}"
            );
        }
    }

    /// Fidelity to `radau5.f`: the stiff relaxation `y₀' = −10⁶ (y₀ − cos y₁)`, `y₁' = 1`
    /// from 0 to t = 1, at the default settings, against the Fortran code (gfortran,
    /// `-ffp-contract=off`, with a driver of the same problem): the same counts of steps,
    /// accepted and rejected steps, decompositions and solves, function evaluations and
    /// Jacobians, and y₀(1) to 1e-14. (The error against the asymptotic solution
    /// `cos t + sin t/λ`, 2e-9 to 4e-9 from 1e-6 down to 1e-10 and 8e-12 at 1e-12, is
    /// RADAU5's own: it keeps the first Jacobian, whose −10⁶ sin y₁ coupling grows from 0,
    /// and its error estimate goes through it.)
    #[test]
    fn matches_radau5_f() {
        struct Relax;
        impl StiffSystem for Relax {
            fn dim(&self) -> usize {
                2
            }
            fn rhs(&self, _t: f64, y: &[f64], f: &mut [f64]) {
                f[0] = -1e6 * (y[0] - libm::cos(y[1]));
                f[1] = 1.0;
            }
        }
        // (tol, y₀(1), n_fcn, n_jac, n_step, n_accept, n_reject, n_dec, n_sol) from
        // radau5.f.
        #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
        let fortran = [
            (1e-6, 0.540303151717348773, 143, 1, 29, 28, 0, 19, 38),
            (1e-8, 0.540303149531565330, 250, 1, 49, 48, 0, 26, 67),
            (1e-10, 0.540303150641765706, 508, 1, 97, 95, 0, 39, 137),
            (1e-12, 0.540303147347222734, 1094, 1, 197, 195, 0, 53, 299),
        ];
        for (tol, y_f, fcn, jac, step, acc, rej, dec, sol) in fortran {
            let mut r = Radau5::new(&Relax, 0.0, &[0.0, 0.0], 1.0, Settings::new(tol, tol));
            while !r.done() {
                r.step(&Relax).expect("integrates");
            }
            let st = r.stats();
            assert_eq!(
                (
                    st.n_fcn,
                    st.n_jac,
                    st.n_step,
                    st.n_accept,
                    st.n_reject,
                    st.n_dec,
                    st.n_sol
                ),
                (fcn, jac, step, acc, rej, dec, sol),
                "tol {tol}"
            );
            assert!(
                (r.y()[0] - y_f).abs() < 1e-14,
                "tol {tol}: {} vs {y_f}",
                r.y()[0]
            );
        }
    }
}
