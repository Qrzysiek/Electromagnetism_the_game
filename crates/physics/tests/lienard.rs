//! Validation tests L1–L4 for the Liénard–Wiechert fields (PHYSICS.md §2.5). Run with
//! `cargo test -p physics --test lienard -- --nocapture --test-threads=1`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

use physics::DVec3;
use physics::antenna::OscillatingDipole;
use physics::lienard::{Worldline, fields};

struct Uniform {
    x0: DVec3,
    v: DVec3,
}

impl Worldline for Uniform {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        (self.x0 + self.v * t, self.v, DVec3::ZERO)
    }
}

/// Circular motion of radius `r` and angular frequency `w` in the plane z = 0.
struct Circle {
    r: f64,
    w: f64,
}

impl Worldline for Circle {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let (s, c) = (self.w * t).sin_cos();
        (
            DVec3::new(c, s, 0.0) * self.r,
            DVec3::new(-s, c, 0.0) * (self.r * self.w),
            DVec3::new(c, s, 0.0) * (-self.r * self.w * self.w),
        )
    }
}

/// A charge circling at radius 1 (ω = 1.2) while oscillating across the plane
/// (`z = 0.3 sin 2.1t`): accelerated in three dimensions, with a varying jerk.
struct Bobbing;

impl Worldline for Bobbing {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let (s, c) = (1.2 * t).sin_cos();
        let (sz, cz) = (2.1 * t).sin_cos();
        (
            DVec3::new(c, s, 0.3 * sz),
            DVec3::new(-1.2 * s, 1.2 * c, 0.63 * cz),
            DVec3::new(-1.44 * c, -1.44 * s, -1.323 * sz),
        )
    }
}

/// Harmonic oscillation `z(t) = A cos(ωt)` along a unit direction.
struct Oscillator {
    amplitude: DVec3,
    w: f64,
}

impl Worldline for Oscillator {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let (s, c) = (self.w * t).sin_cos();
        (
            self.amplitude * c,
            self.amplitude * (-self.w * s),
            self.amplitude * (-self.w * self.w * c),
        )
    }
}

// --- L1: uniformly moving charge ---------------------------------------------------------

/// The field of a uniformly moving charge, from the retarded solution, equals the
/// Heaviside field directed from the *present* position:
/// `E = q (1 − β²) R_p / (R_p³ (1 − β² sin²ψ)^{3/2})`, `B = v × E / c²`
/// (ψ the angle between `R_p` and `v`).
#[test]
fn l1_uniform_motion_is_the_heaviside_field() {
    let (q, c) = (1.3, 2.0);
    let mut worst: f64 = 0.0;
    for (speed, dir) in [(0.1, 0.0), (1.2, 0.7), (1.9, 2.5)] {
        let w = Uniform {
            x0: DVec3::new(0.2, -0.1, 0.0),
            v: DVec3::new(f64::cos(dir), f64::sin(dir), 0.3).normalize() * speed,
        };
        for (x, t) in [
            (DVec3::new(3.0, 1.0, -0.5), 0.4),
            (DVec3::new(-2.0, 0.5, 1.0), 2.0),
            (DVec3::new(0.3, 7.0, 0.0), -1.0),
        ] {
            let f = fields(&w, q, c, x, t);
            let rp = x - w.state(t).0;
            let b2 = speed * speed / (c * c);
            let sin2 = rp.cross(w.v).length_squared() / (rp.length_squared() * speed * speed);
            let e = rp * (q * (1.0 - b2) / (rp.length().powi(3) * (1.0 - b2 * sin2).powf(1.5)));
            let b = w.v.cross(e) / (c * c);
            let err = ((f.e() - e).length() / e.length()).max((f.b - b).length() / e.length() * c);
            assert!(f.e_radiation.length() == 0.0);
            worst = worst.max(err);
        }
    }
    println!("L1: max relative error {worst:.2e}");
    assert!(worst < 1e-12, "L1 {worst:.3e}");
}

// --- L2: slow oscillation = oscillating dipole ---------------------------------------------

/// A charge oscillating with small amplitude `A` (`Aω ≪ c`) radiates like an oscillating
/// dipole `p₀ = qA` (an independent implementation, `antenna.rs`). In the radiation zone
/// the fields agree up to relative corrections `O(Aω/c)` and `O(A/r)`.
#[test]
fn l2_slow_oscillation_matches_the_dipole_antenna() {
    let (q, c, w) = (1.0, 3.0, 1.5);
    let lambda = 2.0 * std::f64::consts::PI * c / w;
    for amp in [1e-3, 1e-4] {
        let osc = Oscillator {
            amplitude: DVec3::new(0.0, amp, 0.0),
            w,
        };
        let dip = OscillatingDipole {
            position: DVec3::ZERO,
            amplitude: DVec3::new(0.0, q * amp, 0.0),
            omega: w,
            phase: 0.0,
            c,
            radius: 0.0,
        };
        let mut worst: f64 = 0.0;
        for x in [
            DVec3::new(20.0 * lambda, 0.0, 0.0),
            DVec3::new(-12.0 * lambda, 3.0, 9.0 * lambda),
        ] {
            for t in [0.0, 0.7, 2.1] {
                let lw = fields(&osc, q, c, x, t);
                let d = dip.fields(x, t);
                // The static Coulomb part q/r² of the point charge is not part of the
                // dipole field: compare the oscillating parts (subtract q n/r²).
                let e_lw = lw.e() - x.normalize() * (q / x.length_squared());
                let scale = d.e.length().max(1e-300);
                worst = worst.max((e_lw - d.e).length() / scale);
                worst = worst.max((lw.b - d.b).length() * c / scale);
            }
        }
        println!(
            "L2 A = {amp:e}: relative difference {worst:.2e} (Aω/c = {:.1e})",
            amp * w / c
        );
        assert!(worst < 20.0 * amp * w / c, "L2 {worst:.3e}");
    }
}

// --- L3: Maxwell equations for arbitrary motion ------------------------------------------

/// The Liénard–Wiechert field of a charge in relativistic circular motion (v = 0.8 c)
/// satisfies the vacuum Maxwell equations away from the charge (central differences).
#[test]
fn l3_circular_motion_satisfies_maxwell_equations() {
    let (q, c) = (1.0, 2.0);
    let orbit = Circle { r: 1.0, w: 1.6 };
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(3.0, 0.5, 0.2),
        DVec3::new(-6.0, 4.0, -2.0),
        DVec3::new(30.0, -25.0, 10.0),
    ] {
        let t = 1.1;
        let h = 1e-4;
        let e = |x: DVec3, t: f64| fields(&orbit, q, c, x, t).e();
        let b = |x: DVec3, t: f64| fields(&orbit, q, c, x, t).b;
        let d = |f: &dyn Fn(DVec3, f64) -> DVec3, axis: DVec3| {
            (f(x + axis * h, t) - f(x - axis * h, t)) / (2.0 * h)
        };
        let dt = |f: &dyn Fn(DVec3, f64) -> DVec3| (f(x, t + h) - f(x, t - h)) / (2.0 * h);
        let curl = |f: &dyn Fn(DVec3, f64) -> DVec3| {
            let (dx, dy, dz) = (d(f, DVec3::X), d(f, DVec3::Y), d(f, DVec3::Z));
            DVec3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x)
        };
        let div = |f: &dyn Fn(DVec3, f64) -> DVec3| {
            d(f, DVec3::X).x + d(f, DVec3::Y).y + d(f, DVec3::Z).z
        };
        let se = curl(&e).length().max(dt(&b).length());
        let sb = curl(&b).length().max(dt(&e).length() / (c * c));
        let r = x.length();
        let res = [
            div(&e).abs() * r / e(x, t).length(),
            div(&b).abs() * r / b(x, t).length(),
            (curl(&e) + dt(&b)).length() / se,
            (curl(&b) - dt(&e) / (c * c)).length() / sb,
        ];
        println!(
            "L3 at {x}: residuals {:.2e} {:.2e} {:.2e} {:.2e}",
            res[0], res[1], res[2], res[3]
        );
        worst = res.iter().copied().fold(worst, f64::max);
    }
    assert!(worst < 1e-6, "L3 {worst:.3e}");
}

// --- L4: the forms of Heaviside and Feynman (Jackson Pr. 6.2) ----------------------------

/// L4, Jackson Pr. 6.2: the fields of a point charge as observation-time derivatives of its
/// retarded direction and distance, Feynman's electric field
/// `E = q {R̂/R² + (R/c) d/dt(R̂/R²) + (1/c²) d²R̂/dt²}` and Heaviside's magnetic field
/// `B = (q/c²) {v × R̂/(κ²R²) + d/dt(v × R̂/κ)/(cR)}` (k = 1, μ₀/4π = 1/c²), against the
/// Liénard–Wiechert fields, for a charge circling while oscillating across the plane
/// (`Bobbing`: β up to 0.68 at c = 2), at six points 1.6 to 50 from its retarded
/// position. Reference:
/// `scripts/wolfram/l4_heaviside_feynman.wls` (Wolfram Engine 14.2), the derivatives
/// symbolic (`d/dt = κ⁻¹ d/dt′`), the retarded time to 50 digits; there the forms of
/// Pr. 6.2(b), from Jefimenko's equations, and the Liénard–Wiechert form agree with these to
/// 40 digits. Required: 1e-12 relative (set before measuring: the retarded time in f64,
/// amplified by up to κ⁻³ ≈ 30).
#[test]
fn l4_heaviside_feynman_forms() {
    // (x, t, E (Feynman), B (Heaviside)).
    type Point = ([f64; 3], f64, [f64; 3], [f64; 3]);
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const REF: [Point; 6] = [
        (
            [3.0, 0.5, 0.2],
            1.1,
            [
                0.14745743814834396,
                0.04357638188038207,
                0.044988831044971223,
            ],
            [
                0.0028772786907419333,
                -0.016461347517220598,
                0.0065138349740077087,
            ],
        ),
        (
            [-6.0, 4.0, -2.0],
            1.1,
            [
                -0.022131009254042458,
                -0.0060199289005101254,
                0.019981998676926141,
            ],
            [
                0.0051182479123875595,
                0.010823709762259083,
                0.0089295349282616275,
            ],
        ),
        (
            [30.0, -25.0, 10.0],
            1.1,
            [
                0.0023259568342506655,
                0.0029839932753845958,
                0.0013718562928862832,
            ],
            [
                -0.00079678737796664984,
                -0.00022149152407970998,
                0.0018327154809670343,
            ],
        ),
        (
            [0.5, -1.5, 1.0],
            2.3333333333333333,
            [
                0.084270609608741477,
                -0.019619516406187891,
                0.10697179985364153,
            ],
            [
                -0.048073825509557495,
                0.0084346771519695532,
                0.039418752179994262,
            ],
        ),
        (
            [0.0, 0.0, 50.0],
            2.3333333333333333,
            [
                -0.0013201835457572312,
                -0.0041641204682224627,
                0.0002198159438535444,
            ],
            [
                0.0020836376379434267,
                -0.00066089143744298233,
                -5.6749159048446928e-6,
            ],
        ),
        (
            [-2.0, -1.0, 0.0],
            3.25,
            [
                -1.9432697816406398,
                -0.069585123673196793,
                -0.82267404341025855,
            ],
            [
                0.31461063726615247,
                -0.42815288066893755,
                -0.70693888775047621,
            ],
        ),
    ];
    let (q, c) = (1.0, 2.0);
    let mut worst: f64 = 0.0;
    for (x, t, e, b) in REF {
        let f = fields(&Bobbing, q, c, DVec3::from_array(x), t);
        let (e, b) = (DVec3::from_array(e), DVec3::from_array(b));
        let de = (f.e() - e).length() / e.length();
        let db = (f.b - b).length() / b.length();
        println!("L4 at {x:?}, t = {t}: E {de:.1e}, B {db:.1e}");
        worst = worst.max(de).max(db);
    }
    assert!(worst < 1e-12, "L4 {worst:.3e}");
}
