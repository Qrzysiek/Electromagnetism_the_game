//! Validation tests L1–L3 for the Liénard–Wiechert fields (PHYSICS.md §2.5). Run with
//! `cargo test -p physics --test lienard -- --nocapture --test-threads=1`.

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
