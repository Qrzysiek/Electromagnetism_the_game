//! Validation tests I1–I5 of the stiff integrator (RADAU5, `integrator/radau5.rs`;
//! PHYSICS.md §5.2): fidelity to Hairer's Fortran code, accuracy against independent
//! high-precision references, and a circuit as a differential-algebraic equation.

#![allow(clippy::disallowed_methods)] // references

use physics::integrator::radau5::{Radau5, Settings, StiffSystem};

/// Van der Pol's equation, ε = 1e-6 (Hairer & Wanner's dr1_radau5 driver).
struct VanDerPol;
impl StiffSystem for VanDerPol {
    fn dim(&self) -> usize {
        2
    }
    fn rhs(&self, _t: f64, y: &[f64], f: &mut [f64]) {
        f[0] = y[1];
        f[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / 1e-6;
    }
}

/// Robertson's chemical kinetics (ROBER, the Test Set for IVP Solvers).
struct Robertson;
impl StiffSystem for Robertson {
    fn dim(&self) -> usize {
        3
    }
    fn rhs(&self, _t: f64, y: &[f64], f: &mut [f64]) {
        f[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
        f[2] = 3.0e7 * y[1] * y[1];
        f[1] = -f[0] - f[2];
    }
}

/// HIRES (the Test Set): eight reactions of plant physiology.
struct Hires;
impl StiffSystem for Hires {
    fn dim(&self) -> usize {
        8
    }
    fn rhs(&self, _t: f64, y: &[f64], f: &mut [f64]) {
        f[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
        f[1] = 1.71 * y[0] - 8.75 * y[1];
        f[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
        f[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
        f[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
        f[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
        f[6] = 280.0 * y[5] * y[7] - 1.81 * y[6];
        f[7] = -f[6];
    }
}

/// Integrates to `t_end`; the final state and the counts (fcn, jac, step, accept, reject,
/// dec, sol).
fn solve(
    sys: &impl StiffSystem,
    y0: &[f64],
    t_end: f64,
    rtol: f64,
    atol: f64,
) -> (Vec<f64>, [u64; 7]) {
    let mut r = Radau5::new(sys, 0.0, y0, t_end, Settings::new(rtol, atol));
    while !r.done() {
        r.step(sys).expect("integrates");
    }
    let s = r.stats();
    (
        r.y().to_vec(),
        [
            s.n_fcn, s.n_jac, s.n_step, s.n_accept, s.n_reject, s.n_dec, s.n_sol,
        ],
    )
}

/// I1: fidelity to `radau5.f` on the three classical problems at rtol = 1e-6, 1e-8,
/// 1e-10 (atol the same; ROBER's 1e-6 times smaller), default settings and numerical
/// Jacobians: the Fortran code compiled with gfortran (`-O2 -ffp-contract=off`) and a
/// driver for each. Required: the same counts of function evaluations, Jacobians, steps,
/// accepted and rejected steps, decompositions and solves (every decision the same), and
/// the final state to 1e-9 relative. The two compilers round `f` and the finite-difference
/// Jacobian differently, and the simplified Newton iteration stops anywhere within
/// `fnewt·rtol'` of its limit (3e-8 at rtol = 1e-6), so the states differ by that much at
/// most. The first requirement, 1e-10, was set without this analysis; HIRES at rtol 1e-6
/// measured 1.7e-10.
#[test]
fn i1_matches_radau5_f() {
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let vdp = [
        (
            1e-6,
            vec![0.17061677081875102e1, -0.89280963355092913e0],
            [3894, 308, 487, 476, 8, 404, 1139],
        ),
        (
            1e-8,
            vec![0.17061674351089380e1, -0.89281001891180400e0],
            [8133, 482, 1027, 1019, 7, 834, 2371],
        ),
        (
            1e-10,
            vec![0.17061674374891436e1, -0.89281001659555193e0],
            [17084, 615, 2215, 2211, 2, 1705, 4957],
        ),
    ];
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let rober = [
        (
            1e-6,
            vec![
                0.20830769906625651e-07,
                0.83323081341356824e-13,
                0.99999997916914563e+00,
            ],
            [4294, 321, 471, 331, 1, 470, 1321],
        ),
        (
            1e-8,
            vec![
                0.20833376469098075e-07,
                0.83333507591683976e-13,
                0.99999997916653882e+00,
            ],
            [7092, 547, 772, 621, 2, 716, 2157],
        ),
        (
            1e-10,
            vec![
                0.20833401313843985e-07,
                0.83333606970672004e-13,
                0.99999997916651495e+00,
            ],
            [13478, 852, 1470, 1268, 2, 1096, 4070],
        ),
    ];
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let hires = [
        (
            1e-6,
            vec![
                0.73712636377581541e-03,
                0.14424761104315766e-03,
                0.58886368396438374e-04,
                0.11756425112130239e-02,
                0.23861995098027329e-02,
                0.62384451683828845e-02,
                0.28499028269742457e-02,
                0.28500971730257537e-02,
            ],
            [486, 27, 58, 51, 0, 50, 145],
        ),
        (
            1e-8,
            vec![
                0.73713142216297713e-03,
                0.14424860584469752e-03,
                0.58887324365183097e-04,
                0.11756517167289284e-02,
                0.23863591084184528e-02,
                0.62389658678273091e-02,
                0.28500097010432175e-02,
                0.28499902989567796e-02,
            ],
            [829, 36, 95, 94, 0, 60, 245],
        ),
        (
            1e-10,
            vec![
                0.73713126074725780e-03,
                0.14424857332528870e-03,
                0.58887297907294896e-04,
                0.11756513499812332e-02,
                0.23863562466577017e-02,
                0.62389679847477302e-02,
                0.28499986207993029e-02,
                0.28500013792006991e-02,
            ],
            [1647, 62, 196, 195, 1, 96, 484],
        ),
    ];
    let mut worst: f64 = 0.0;
    let mut check =
        |name: &str, tol: f64, got: (Vec<f64>, [u64; 7]), want: &[f64], counts: [u64; 7]| {
            let rel = got
                .0
                .iter()
                .zip(want)
                .map(|(a, b)| (a - b).abs() / b.abs())
                .fold(0.0, f64::max);
            println!(
                "I1 {name} rtol {tol:e}: counts {:?} (Fortran {counts:?}); final state {rel:.1e}",
                got.1
            );
            assert_eq!(got.1, counts, "{name} {tol}");
            worst = worst.max(rel);
        };
    for (tol, want, counts) in &vdp {
        check(
            "Van der Pol",
            *tol,
            solve(&VanDerPol, &[2.0, -0.66], 2.0, *tol, *tol),
            want,
            *counts,
        );
    }
    for (tol, want, counts) in &rober {
        check(
            "ROBER",
            *tol,
            solve(&Robertson, &[1.0, 0.0, 0.0], 1e11, *tol, *tol * 1e-6),
            want,
            *counts,
        );
    }
    for (tol, want, counts) in &hires {
        let mut y0 = [0.0; 8];
        y0[0] = 1.0;
        y0[7] = 0.0057;
        check(
            "HIRES",
            *tol,
            solve(&Hires, &y0, 321.8122, *tol, *tol),
            want,
            *counts,
        );
    }
    assert!(worst < 1e-9, "I1 {worst:.3e}");
}

/// I5: a circuit as a differential-algebraic equation `M y' = f` (modified nodal analysis,
/// as the circuits will be written): a source `V₀ sin ωt` drives R, L and C in series.
/// Unknowns: the node voltages v₁ (the source's), v₂ (between R and L), v₃ (the
/// capacitor's), the source current i_s and the inductor current i_L; the mass matrix
/// diag(0, 0, C, 0, L) makes v₁, v₂ and i_s algebraic (index 1). From rest, against the
/// closed form: `q = C v₃` with `q'' + (R/L) q' + q/(LC) = (V₀/L) sin ωt`, the driven
/// solution plus the damped free one (R = 0.5, L = C = 1, ω = 0.7, to t = 20, rtol = atol
/// = 1e-10). Required (set before measuring): every unknown at the step ends within 1e-8
/// (100× the tolerance), and the dense output (the collocation polynomial, of lower order
/// inside the step) at each step's midpoint within 1e-7.
#[test]
fn i5_series_rlc_dae() {
    struct Rlc {
        r: f64,
        l: f64,
        c: f64,
        v0: f64,
        w: f64,
    }
    impl StiffSystem for Rlc {
        fn dim(&self) -> usize {
            5
        }
        fn rhs(&self, t: f64, y: &[f64], f: &mut [f64]) {
            let (v1, v2, v3, is, il) = (y[0], y[1], y[2], y[3], y[4]);
            f[0] = is - (v1 - v2) / self.r; // KCL at node 1
            f[1] = (v1 - v2) / self.r - il; // KCL at node 2
            f[2] = il; // C v3' = i_L
            f[3] = v1 - self.v0 * (self.w * t).sin(); // the source
            f[4] = v2 - v3; // L i_L' = v2 − v3
        }
        fn mass(&self) -> Option<Vec<f64>> {
            let mut m = vec![0.0; 25];
            m[2 + 5 * 2] = self.c;
            m[4 + 5 * 4] = self.l;
            Some(m)
        }
    }
    let sys = Rlc {
        r: 0.5,
        l: 1.0,
        c: 1.0,
        v0: 1.0,
        w: 0.7,
    };
    // The closed form.
    let w0sq = 1.0 / (sys.l * sys.c);
    let g = sys.r / (2.0 * sys.l);
    let wd = (w0sq - g * g).sqrt();
    let d = (w0sq - sys.w * sys.w).powi(2) + (2.0 * g * sys.w).powi(2);
    let a_p = sys.v0 / sys.l * (w0sq - sys.w * sys.w) / d;
    let b_p = -sys.v0 / sys.l * 2.0 * g * sys.w / d;
    let a_h = -b_p;
    let b_h = (g * a_h - a_p * sys.w) / wd;
    let exact = |t: f64| {
        let e = (-g * t).exp();
        let (s, c) = ((sys.w * t).sin(), (sys.w * t).cos());
        let (sd, cd) = ((wd * t).sin(), (wd * t).cos());
        let q = a_p * s + b_p * c + e * (a_h * cd + b_h * sd);
        let dq = sys.w * (a_p * c - b_p * s)
            + e * (-g * (a_h * cd + b_h * sd) + wd * (-a_h * sd + b_h * cd));
        let v3 = q / sys.c;
        let v1 = sys.v0 * s;
        let v2 = v1 - sys.r * dq;
        [v1, v2, v3, dq, dq]
    };
    let mut r = Radau5::new(&sys, 0.0, &[0.0; 5], 20.0, Settings::new(1e-10, 1e-10));
    let (mut worst, mut worst_dense): (f64, f64) = (0.0, 0.0);
    let mut y = [0.0; 5];
    while !r.done() {
        r.step(&sys).expect("integrates");
        let want = exact(r.t());
        for (got, w) in r.y().iter().zip(want) {
            worst = worst.max((got - w).abs());
        }
        let dense = r.dense();
        let mid = 0.5 * (dense.t_start() + dense.t_end());
        dense.eval(mid, &mut y);
        let want = exact(mid);
        for (got, w) in y.iter().zip(want) {
            worst_dense = worst_dense.max((got - w).abs());
        }
    }
    println!(
        "I5: series RLC to t = 20: step ends within {worst:.1e}, the dense output at the \
         midpoints within {worst_dense:.1e}; {:?}",
        r.stats()
    );
    assert!(worst < 1e-8, "I5 {worst:.3e}");
    assert!(worst_dense < 1e-7, "I5 dense {worst_dense:.3e}");
}

/// I2–I4: accuracy on the three classical problems (Van der Pol at t = 2, ROBER at
/// t = 1e11, HIRES at t = 321.8122) against independent references: NDSolve
/// (`scripts/wolfram/i_stiff_references.wls`) at 34-digit working precision for Van der
/// Pol and HIRES (the 30-digit run agrees to 6e-13 and 1e-13), at 30 digits for ROBER (at
/// 34 NDSolve failed its own error test at t = 7e8); HIRES and ROBER agree with the Test
/// Set for IVP Solvers' published values to 14 digits. Criterion: RADAU5's design, its
/// tolerances transformed so that
/// the global error is about the tolerance: every component within 10 times the
/// tolerance's scale `atol + rtol |y|`, at rtol = 1e-6, 1e-8, 1e-10. (Not a threshold set
/// blind: the Fortran code's results, I1, had been seen.)
#[test]
fn i2_i4_accuracy_against_references() {
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let vdp_ref = [1.70616743754317598094, -0.89281001655112012742];
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let rober_ref = [
        2.08334014970137718425e-8,
        8.33336077033511521538e-14,
        0.99999997916651516938,
    ];
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let hires_ref = [
        7.3713125733256681729e-4,
        1.4424857263161847303e-4,
        5.8887297409675756881e-5,
        1.17565134328314921361e-3,
        2.38635619883133156652e-3,
        6.2389682527427992309e-3,
        2.84999839518576942912e-3,
        2.85000160481423057088e-3,
    ];
    let mut worst: f64 = 0.0;
    let mut check = |name: &str, tol: f64, atol: f64, got: &[f64], want: &[f64]| {
        let scaled = got
            .iter()
            .zip(want)
            .map(|(g, w)| (g - w).abs() / (atol + tol * w.abs()))
            .fold(0.0, f64::max);
        let rel = got
            .iter()
            .zip(want)
            .map(|(g, w)| (g - w).abs() / w.abs())
            .fold(0.0, f64::max);
        println!(
            "I {name} rtol {tol:e}: error {scaled:.2} of the tolerance's scale (relative {rel:.1e})"
        );
        worst = worst.max(scaled);
    };
    for tol in [1e-6, 1e-8, 1e-10] {
        let (y, _) = solve(&VanDerPol, &[2.0, -0.66], 2.0, tol, tol);
        check("2 Van der Pol", tol, tol, &y, &vdp_ref);
        let (y, _) = solve(&Robertson, &[1.0, 0.0, 0.0], 1e11, tol, tol * 1e-6);
        check("3 ROBER", tol, tol * 1e-6, &y, &rober_ref);
        let mut y0 = [0.0; 8];
        y0[0] = 1.0;
        y0[7] = 0.0057;
        let (y, _) = solve(&Hires, &y0, 321.8122, tol, tol);
        check("4 HIRES", tol, tol, &y, &hires_ref);
    }
    assert!(worst < 10.0, "I2-I4 {worst:.2}");
}
