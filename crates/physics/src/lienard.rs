//! Fields of a moving point charge: the Liénard–Wiechert fields (PHYSICS.md §2.5).
//!
//! With `R = x − r(t_r)`, `n = R/|R|`, `β = v/c`, `κ = 1 − n·β`, all evaluated at the
//! retarded time `t_r` defined by `|x − r(t_r)| = c (t − t_r)` (units `k = 1`,
//! force `q(E + v×B)`; Jackson §14.1):
//!
//! ```text
//! E = q (n − β)(1 − β²) / (κ³ R²)  +  (q/c) n × ((n − β) × β̇) / (κ³ R)
//! B = n × E / c
//! ```
//!
//! The first term is the velocity (generalized Coulomb) field, the second the radiation
//! field, proportional to the acceleration. Used to visualize the field a particle
//! carries and radiates; its effect back on the particle is radiation reaction
//! (`dynamics.rs`).

use glam::DVec3;

/// A world line: position, velocity and acceleration at any time. Speeds must stay
/// below `c`.
pub trait Worldline {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3);
}

/// Retarded time for the field point `x` at time `t`: the root of
/// `g(t_r) = c (t − t_r) − |x − r(t_r)|`, which is strictly decreasing for `|v| < c`.
/// Found by bracketing and bisection to full precision.
pub fn retarded_time(w: &impl Worldline, c: f64, x: DVec3, t: f64) -> f64 {
    let g = |tr: f64| c * (t - tr) - (x - w.state(tr).0).length();
    let mut hi = t;
    if g(hi) >= 0.0 {
        return hi;
    }
    // g(t − s) ≥ c s − |x − r(t)| − v_max s: expand the bracket until g > 0.
    let mut step = (x - w.state(t).0).length() / c;
    let mut lo = t - step;
    while g(lo) <= 0.0 {
        hi = lo;
        step *= 2.0;
        lo = t - step;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if g(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Liénard–Wiechert `(E, B)` of charge `q` on world line `w`, at `x` and time `t`.
/// Returns the velocity and radiation parts of `E` separately (`B = n × E / c` for
/// their sum).
pub fn fields(w: &impl Worldline, q: f64, c: f64, x: DVec3, t: f64) -> LwField {
    let tr = retarded_time(w, c, x, t);
    let (r, v, a) = w.state(tr);
    let d = x - r;
    let dist = d.length();
    let n = d / dist;
    let beta = v / c;
    let beta_dot = a / c;
    let kappa = 1.0 - n.dot(beta);
    let k3 = kappa * kappa * kappa;
    let e_vel = (n - beta) * (q * (1.0 - beta.length_squared()) / (k3 * dist * dist));
    let e_rad = n.cross((n - beta).cross(beta_dot)) * (q / (c * k3 * dist));
    let b = n.cross(e_vel + e_rad) / c;
    LwField {
        e_velocity: e_vel,
        e_radiation: e_rad,
        b,
        retarded_time: tr,
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LwField {
    pub e_velocity: DVec3,
    pub e_radiation: DVec3,
    pub b: DVec3,
    pub retarded_time: f64,
}

impl LwField {
    pub fn e(&self) -> DVec3 {
        self.e_velocity + self.e_radiation
    }
}

/// A world line known at samples (e.g. a computed flight), interpolated inside the
/// sampled interval (cubic Hermite for the position, linear for velocity and
/// acceleration) and continued with uniform motion before the first and after the last
/// sample.
#[derive(Clone, Debug, Default)]
pub struct SampledWorldline {
    t: Vec<f64>,
    x: Vec<DVec3>,
    v: Vec<DVec3>,
    a: Vec<DVec3>,
}

impl SampledWorldline {
    /// Samples `(t, x, v, a)` with strictly increasing `t` (at least one).
    pub fn new(samples: &[(f64, DVec3, DVec3, DVec3)]) -> Self {
        assert!(!samples.is_empty());
        Self {
            t: samples.iter().map(|s| s.0).collect(),
            x: samples.iter().map(|s| s.1).collect(),
            v: samples.iter().map(|s| s.2).collect(),
            a: samples.iter().map(|s| s.3).collect(),
        }
    }

    pub fn t_start(&self) -> f64 {
        self.t[0]
    }

    pub fn t_end(&self) -> f64 {
        *self.t.last().expect("non-empty")
    }
}

impl Worldline for SampledWorldline {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let n = self.t.len();
        if t <= self.t[0] {
            return (
                self.x[0] + self.v[0] * (t - self.t[0]),
                self.v[0],
                DVec3::ZERO,
            );
        }
        if t >= self.t[n - 1] {
            return (
                self.x[n - 1] + self.v[n - 1] * (t - self.t[n - 1]),
                self.v[n - 1],
                DVec3::ZERO,
            );
        }
        let i = self.t.partition_point(|&s| s <= t).clamp(1, n - 1);
        let (t0, t1) = (self.t[i - 1], self.t[i]);
        let h = t1 - t0;
        let s = (t - t0) / h;
        let (x0, x1, v0, v1) = (self.x[i - 1], self.x[i], self.v[i - 1], self.v[i]);
        // Cubic Hermite position; its derivative for consistency is not needed (the
        // velocity is interpolated linearly, like the acceleration).
        let h00 = (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s);
        let h10 = s * (1.0 - s) * (1.0 - s);
        let h01 = s * s * (3.0 - 2.0 * s);
        let h11 = s * s * (s - 1.0);
        let x = x0 * h00 + v0 * (h10 * h) + x1 * h01 + v1 * (h11 * h);
        (x, v0.lerp(v1, s), self.a[i - 1].lerp(self.a[i], s))
    }
}
