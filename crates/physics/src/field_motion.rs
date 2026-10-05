//! The motion of a charge in uniform fields, in the plane: `E` in the plane, `B` along z
//! (the fields of the plane `z = 0` always are), exact and relativistic (PHYSICS.md §3.3).
//! The quasi-static beam interaction continues each source's past along it: the motion
//! it would have had in the fields it feels now.
//!
//! In proper time `s`, with the 4-velocity `U = (γc, γv_x, γv_y)`,
//! `dU/ds = Λ U`, `Λ = (q/m) [[0, E_x/c, E_y/c], [E_x/c, 0, B], [E_y/c, −B, 0]]`.
//! Its characteristic polynomial is `λ (λ² + ω²)`, `ω² = (q/m)² (B² − E²/c²)`, so
//! `Λ³ = −ω² Λ` and
//!
//! ```text
//! U(s) = U₀ + S₁(s) Λ U₀ + C₁(s) Λ² U₀,        X(s) = X₀ + s U₀ + C₁(s) Λ U₀ + S₂(s) Λ² U₀
//! S₁ = sin(ωs)/ω,  C₁ = (1 − cos ωs)/ω²,  S₂ = (s − sin(ωs)/ω)/ω²
//! ```
//!
//! (hyperbolic functions for `E > cB`, powers of `s` for `E = cB`: all three are the
//! series `sᵏ Σₙ (−ω²s²)ⁿ/(2n + k)!`). The speed stays below c.
//!
//! The fields a source feels change along its path, which motion in uniform fields does
//! not follow: the jerk that misses (`JerkCorrection`) is added to the continued past as
//! a tapered cubic term in coordinate time, and the retarded point is searched on the
//! corrected world line (`FieldMotion::retarded_with`).

use glam::DVec3;

/// A continued past's jerk from the change of the fields along the source's path, which
/// the motion in uniform fields does not have (the quasi-static beam interaction,
/// PHYSICS.md §3.3). Added to the motion as `δa(τ) = Δȧ T g(τ/T)`, `g(u) = u (1 − u²)³`
/// for `|u| < 1` and 0 beyond, with the velocity and position changes `δv = Δȧ T² G₁(u)`,
/// `δx = Δȧ T³ G₂(u)` (`G₁' = g`, `G₂' = G₁`, both zero at `u = 0`): `Δȧ τ`, `Δȧ τ²/2`,
/// `Δȧ τ³/6` near the present (to relative `O(u²)`: `G₂ = u³/6 − 3u⁵/20 + …`), the jerk
/// fading within `T` (the acceleration continuous to its second derivative), the velocity
/// changed by at most `|Δȧ| T²/8`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JerkCorrection {
    /// `Δȧ`, the jerk added at the present.
    pub jerk: DVec3,
    /// `T`, the time over which it fades (positive and finite).
    pub scale: f64,
}

/// `G₂(1) = 1/6 − 3/20 + 1/14 − 1/72`.
const G2_ONE: f64 = 187.0 / 2520.0;

impl JerkCorrection {
    /// No correction.
    pub const NONE: Self = Self {
        jerk: DVec3::ZERO,
        scale: 1.0,
    };

    /// `(δx, δv, δa)` at the coordinate time `τ` from the present.
    pub fn at(&self, tau: f64) -> (DVec3, DVec3, DVec3) {
        let big_t = self.scale;
        let u = tau / big_t;
        let (g, g1, g2) = if u.abs() < 1.0 {
            let u2 = u * u;
            let w = 1.0 - u2;
            (
                u * w * w * w,
                // u²/2 − 3u⁴/4 + u⁶/2 − u⁸/8 = (1 − (1 − u²)⁴)/8.
                u2 * (0.5 - u2 * (0.75 - u2 * (0.5 - u2 / 8.0))),
                // u³/6 − 3u⁵/20 + u⁷/14 − u⁹/72.
                u * u2 * (1.0 / 6.0 - u2 * (0.15 - u2 * (1.0 / 14.0 - u2 / 72.0))),
            )
        } else {
            (0.0, 0.125, u.signum() * (G2_ONE + (u.abs() - 1.0) / 8.0))
        };
        (
            self.jerk * (big_t * big_t * big_t * g2),
            self.jerk * (big_t * big_t * g1),
            self.jerk * (big_t * g),
        )
    }
}

/// The motion of a charge in uniform in-plane fields (`FieldMotion::new`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldMotion {
    /// Coordinate time and position at `s = 0`.
    t0: f64,
    x0: DVec3,
    /// `U₀`, `Λ U₀`, `Λ² U₀` as `(ct, x, y)` components.
    u0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    /// `ω² = (q/m)² (B² − E²/c²)` (negative for `E > cB`).
    w2: f64,
    c: f64,
}

/// `(C₀, S₁, C₁, S₂)` at `s`: `C₀ = cos ωs` (`= 1 − ω² C₁`), and the functions above.
/// For `|ω²s²| < 0.01` five terms of their series (the next is below 3e-18 of the first);
/// otherwise the closed forms, `C₁` by the half angle (no cancellation).
fn functions(w2: f64, s: f64) -> (f64, f64, f64, f64) {
    let z = w2 * s * s;
    if z.abs() < 1e-2 {
        let s1 = s * (1.0 - z / 6.0 * (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0))));
        let c1 =
            s * s * 0.5 * (1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0))));
        let s2 = s * s * s / 6.0
            * (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0))));
        return (1.0 - w2 * c1, s1, c1, s2);
    }
    if w2 > 0.0 {
        let w = w2.sqrt();
        let (sn, cs) = libm::sincos(w * s);
        let half = libm::sin(0.5 * w * s);
        (cs, sn / w, 2.0 * half * half / w2, (s - sn / w) / w2)
    } else {
        let k = (-w2).sqrt();
        let (sh, ch) = (libm::sinh(k * s), libm::cosh(k * s));
        let half = libm::sinh(0.5 * k * s);
        (ch, sh / k, 2.0 * half * half / -w2, (sh / k - s) / -w2)
    }
}

impl FieldMotion {
    /// The motion of a charge `q`, mass `m` that is at `x0` with velocity `v0` at time
    /// `t0`, in the uniform fields `e`, `b`. `None` unless it stays in the plane z = 0
    /// (position, velocity and `E` in it, `B` along z), or for a neutral or massless one.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        q: f64,
        m: f64,
        c: f64,
        t0: f64,
        x0: DVec3,
        v0: DVec3,
        e: DVec3,
        b: DVec3,
    ) -> Option<Self> {
        let planar = x0.z == 0.0 && v0.z == 0.0 && e.z == 0.0 && b.x == 0.0 && b.y == 0.0;
        if !planar || q == 0.0 || m <= 0.0 || !c.is_finite() || v0.length() >= c {
            return None;
        }
        let gamma = 1.0 / (1.0 - v0.length_squared() / (c * c)).sqrt();
        let u0 = [gamma * c, gamma * v0.x, gamma * v0.y];
        let k = q / m;
        let lambda = |u: [f64; 3]| {
            [
                k * (e.x * u[1] + e.y * u[2]) / c,
                k * (e.x * u[0] / c + b.z * u[2]),
                k * (e.y * u[0] / c - b.z * u[1]),
            ]
        };
        let v1 = lambda(u0);
        let v2 = lambda(v1);
        Some(Self {
            t0,
            x0,
            u0,
            v1,
            v2,
            w2: k * k * (b.z * b.z - (e.x * e.x + e.y * e.y) / (c * c)),
            c,
        })
    }

    /// `U(s)`, `dU/ds` and `X(s) − X₀`.
    fn at(&self, s: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (c0, s1, c1, s2) = functions(self.w2, s);
        let mut u = [0.0; 3];
        let mut du = [0.0; 3];
        let mut dx = [0.0; 3];
        for i in 0..3 {
            u[i] = self.u0[i] + s1 * self.v1[i] + c1 * self.v2[i];
            du[i] = c0 * self.v1[i] + s1 * self.v2[i];
            dx[i] = s * self.u0[i] + c1 * self.v1[i] + s2 * self.v2[i];
        }
        (u, du, dx)
    }

    /// Coordinate time, position, velocity and acceleration at proper time `s`.
    pub fn state(&self, s: f64) -> (f64, DVec3, DVec3, DVec3) {
        let (u, du, dx) = self.at(s);
        self.state_from(u, du, dx)
    }

    /// `state` from an evaluation `at`.
    fn state_from(&self, u: [f64; 3], du: [f64; 3], dx: [f64; 3]) -> (f64, DVec3, DVec3, DVec3) {
        let c = self.c;
        let t = self.t0 + dx[0] / c;
        let x = self.x0 + DVec3::new(dx[1], dx[2], 0.0);
        let v = DVec3::new(u[1], u[2], 0.0) * (c / u[0]);
        // a = dv/dt = (c/U⁰) d/ds (c U/U⁰).
        let a = (DVec3::new(du[1], du[2], 0.0) * u[0] - DVec3::new(u[1], u[2], 0.0) * du[0])
            * (c * c / (u[0] * u[0] * u[0]));
        (t, x, v, a)
    }

    /// The jerk `da/dt` at `s = 0`.
    pub fn jerk(&self) -> DVec3 {
        let c = self.c;
        let (u0, u1, u2) = (self.u0, self.v1, self.v2);
        let sp = |w: [f64; 3]| DVec3::new(w[1], w[2], 0.0);
        // a = c² N/U⁰³ with N = U′ U⁰ − U U⁰′ (spatial), N′ = U″ U⁰ − U U⁰″.
        let n = sp(u1) * u0[0] - sp(u0) * u1[0];
        let dn = sp(u2) * u0[0] - sp(u0) * u2[0];
        let da_ds = (dn / u0[0].powi(3) - n * (3.0 * u1[0] / u0[0].powi(4))) * (c * c);
        da_ds * (c / u0[0])
    }

    /// The snap `d²a/dt²` at `s = 0`.
    pub fn snap(&self) -> DVec3 {
        let c = self.c;
        let sp = |w: [f64; 3]| DVec3::new(w[1], w[2], 0.0);
        // With W = U⁰ and u the spatial part: U' = Λ U₀, U'' = Λ² U₀, U''' = −ω² Λ U₀ at
        // s = 0; a = c² N W⁻³ with N = u'W − uW', and d/dt = (c/W) d/ds.
        let (w0, w1, w2) = (self.u0[0], self.v1[0], self.v2[0]);
        let w3 = -self.w2 * self.v1[0];
        let (u, u1, u2) = (sp(self.u0), sp(self.v1), sp(self.v2));
        let u3 = u1 * -self.w2;
        let n = u1 * w0 - u * w1;
        let dn = u2 * w0 - u * w2;
        let ddn = u3 * w0 + u2 * w1 - u1 * w2 - u * w3;
        (ddn / w0.powi(5) - dn * (7.0 * w1 / w0.powi(6)) - n * (3.0 * w2 / w0.powi(6))
            + n * (15.0 * w1 * w1 / w0.powi(7)))
            * c.powi(4)
    }
}

/// Liénard–Wiechert fields `(E, B)` at `x`, at the time `dt` after the motion's reference
/// time, of the charge `q` moving as `motion` (`FieldMotion::retarded`).
/// `None` where the motion has no retarded point (`FieldMotion::retarded`).
pub fn field_motion_fields(
    q: f64,
    x: DVec3,
    dt: f64,
    motion: &FieldMotion,
) -> Option<(DVec3, DVec3)> {
    let r = motion.retarded(x, dt)?;
    Some(motion.fields(q, x, &r))
}

/// The retarded point of a `FieldMotion` for an observer: proper time, coordinate time,
/// position, velocity, acceleration, and the excursion `|U(s) − U₀|/c` (how far the motion
/// has swung by then: about the change of velocity over c for slow particles, `2γβ
/// sin(θ/2)` after turning by θ).
#[derive(Clone, Copy, Debug)]
pub struct Retarded {
    pub s: f64,
    /// The emission time, in the same clock as the motion's reference time `t₀` (not
    /// relative to the observer).
    pub t: f64,
    pub x: DVec3,
    pub v: DVec3,
    pub a: DVec3,
    pub excursion: f64,
}

impl FieldMotion {
    /// `(U₀, Λ U₀, Λ² U₀, ω²)` (for the game's GPU field view, which evaluates the same
    /// closed form).
    pub fn parameters(&self) -> ([f64; 3], [f64; 3], [f64; 3], f64) {
        (self.u0, self.v1, self.v2, self.w2)
    }

    /// The motion's present velocity and acceleration (at `s = 0`).
    pub fn present(&self) -> (DVec3, DVec3) {
        let (_, _, v, a) = self.state(0.0);
        (v, a)
    }

    /// Liénard–Wiechert fields `(E, B)` at `x` of the charge `q` at its retarded point `r`.
    pub fn fields(&self, q: f64, x: DVec3, r: &Retarded) -> (DVec3, DVec3) {
        let f = crate::lienard::fields_from(q, self.c, x, r.t, r.x, r.v, r.a);
        (f.e(), f.b)
    }

    /// The retarded point for an observer at `x`, at the time `dt` after the reference time:
    /// the root of the strictly decreasing `G(s) = c (t₀ + dt − t(s)) − |x − x(s)|` (the
    /// world line is timelike). Newton's method from the uniform-motion guess, kept inside
    /// the bracket once both sides are known (bisection otherwise); if it has not converged
    /// in 12 steps, the bracket is searched by doubling first. `None` if there is no root:
    /// a charge accelerated forever by `E > cB` never reaches observers beyond its horizon
    /// (and far back its hyperbolic functions overflow): the search gives up after 64
    /// doublings or at a value that is not finite. Also `None` where `G` is not resolved:
    /// far back on such a world line `c (t₀ + dt − t(s))` and `|x − x(s)|` grow beyond
    /// 1e5 times the observer's distance and their difference is rounding noise (it once
    /// gave a "root" at an excursion of 1e73).
    pub fn retarded(&self, x: DVec3, dt: f64) -> Option<Retarded> {
        self.retarded_with(x, dt, &JerkCorrection::NONE)
    }

    /// `retarded` on the world line with the jerk correction `corr` added (in coordinate
    /// time from the reference time). The point's position, velocity and acceleration
    /// include the correction; its excursion is the motion's own.
    pub fn retarded_with(&self, x: DVec3, dt: f64, corr: &JerkCorrection) -> Option<Retarded> {
        let c = self.c;
        // `G`, its slope, its rounding error, and the evaluation.
        let eval = |s: f64| {
            let (u, du, dx) = self.at(s);
            let delta = corr.at(dx[0] / c);
            let d = x - (self.x0 + DVec3::new(dx[1], dx[2], 0.0) + delta.0);
            let dist = d.length();
            let value = c * dt - dx[0] - dist;
            // dG/ds = −U⁰ + n·dX/ds, dX/ds = u + δv U⁰/c (negative: the world line is
            // timelike).
            let slope = -u[0]
                + d.dot(DVec3::new(u[1], u[2], 0.0) + delta.1 * (u[0] / c)) / dist.max(1e-300);
            let noise = 4.0 * f64::EPSILON * ((c * dt).abs() + dx[0].abs() + dist);
            (value, slope, noise, (u, du, dx, delta))
        };
        #[allow(clippy::type_complexity)]
        let done_at = |s: f64, ev: ([f64; 3], [f64; 3], [f64; 3], (DVec3, DVec3, DVec3))| {
            let (u, du, dx, delta) = ev;
            let (t, xr, v, a) = self.state_from(u, du, dx);
            let d = [u[0] - self.u0[0], u[1] - self.u0[1], u[2] - self.u0[2]];
            Retarded {
                s,
                t,
                x: xr + delta.0,
                v: v + delta.1,
                a: a + delta.2,
                excursion: (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / c,
            }
        };
        // The uniform-motion guess, in proper time.
        let gamma0 = self.u0[0] / c;
        let v0 = DVec3::new(self.u0[1], self.u0[2], 0.0) / gamma0;
        let r = x - self.x0 - v0 * dt;
        let (rv, v2) = (r.dot(v0), v0.length_squared());
        let tau_u = dt - (rv + (rv * rv + (c * c - v2) * r.length_squared()).sqrt()) / (c * c - v2);
        let scale = (r.length() / c).max(dt.abs()).max(1e-300) / gamma0;
        // `G` is resolved where its rounding error is below 1e-10 of the observer's
        // distance (or light time).
        let unresolved = |noise: f64| noise > 1e-10 * r.length().max(c * dt.abs()).max(1e-300);
        let tol = |s: f64| 4.0 * f64::EPSILON * (s.abs() + scale);
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        let mut s = tau_u / gamma0;
        for iteration in 0..100 {
            if iteration == 12 && !(lo.is_finite() && hi.is_finite()) {
                // Slow: bracket by doubling from where it is.
                let mut step = scale;
                for _ in 0..64 {
                    if lo.is_finite() {
                        break;
                    }
                    let trial = hi.min(s) - step;
                    let (value, _, noise, _) = eval(trial);
                    if !value.is_finite() || unresolved(noise) {
                        return None;
                    }
                    if value >= 0.0 {
                        lo = trial;
                    } else {
                        hi = trial;
                        step *= 2.0;
                    }
                }
                for _ in 0..64 {
                    if hi.is_finite() {
                        break;
                    }
                    let trial = lo.max(s) + step;
                    let (value, _, noise, _) = eval(trial);
                    if !value.is_finite() || unresolved(noise) {
                        return None;
                    }
                    if value <= 0.0 {
                        hi = trial;
                    } else {
                        lo = trial;
                        step *= 2.0;
                    }
                }
                if !(lo.is_finite() && hi.is_finite()) {
                    return None;
                }
                s = 0.5 * (lo + hi);
            }
            let (value, slope, noise, ev) = eval(s);
            if !(value.is_finite() && slope.is_finite()) || unresolved(noise) {
                return None;
            }
            if value > 0.0 {
                lo = lo.max(s);
            } else {
                hi = hi.min(s);
            }
            let mut next = s - value / slope;
            if lo.is_finite() && hi.is_finite() && !(next > lo && next < hi) {
                next = 0.5 * (lo + hi);
            }
            if !next.is_finite() {
                next = if value > 0.0 { s + scale } else { s - scale };
            }
            if (next - s).abs() <= tol(s) || hi - lo <= tol(s) {
                return Some(done_at(s, ev));
            }
            s = next;
        }
        let (value, _, noise, ev) = eval(s);
        (value.is_finite() && !unresolved(noise)).then(|| done_at(s, ev))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The closed form against a fine Runge–Kutta integration of the Lorentz force in
    /// coordinate time, in magnetic-dominated, electric-dominated and nearly null fields,
    /// back and forth in proper time.
    #[test]
    fn closed_form_matches_integration() {
        let (q, m, c) = (0.7, 1.3, 5.0);
        let (x0, v0) = (DVec3::new(1.0, 2.0, 0.0), DVec3::new(3.0, -1.5, 0.0));
        for (e, b) in [
            (DVec3::new(0.8, -0.4, 0.0), 1.1),
            (DVec3::new(2.5, 1.0, 0.0), 0.2),
            (DVec3::new(5.5, 0.0, 0.0), 1.1),
        ] {
            let b = DVec3::new(0.0, 0.0, b);
            let motion = FieldMotion::new(q, m, c, 0.0, x0, v0, e, b).expect("planar");
            for s in [-1.3, -0.2, 0.7] {
                let (t, x, v, _) = motion.state(s);
                // RK4 in t for (x, p), 20 000 steps.
                let kin = |p: DVec3| p / (m * (1.0 + p.length_squared() / (m * c * m * c)).sqrt());
                let f = |(_, p): (DVec3, DVec3)| {
                    let vel = kin(p);
                    (vel, (e + vel.cross(b)) * q)
                };
                let n = 20_000;
                let h = t / f64::from(n);
                let gamma = 1.0 / (1.0 - v0.length_squared() / (c * c)).sqrt();
                let mut y = (x0, v0 * (gamma * m));
                for _ in 0..n {
                    let k1 = f(y);
                    let k2 = f((y.0 + k1.0 * (0.5 * h), y.1 + k1.1 * (0.5 * h)));
                    let k3 = f((y.0 + k2.0 * (0.5 * h), y.1 + k2.1 * (0.5 * h)));
                    let k4 = f((y.0 + k3.0 * h, y.1 + k3.1 * h));
                    y = (
                        y.0 + (k1.0 + (k2.0 + k3.0) * 2.0 + k4.0) * (h / 6.0),
                        y.1 + (k1.1 + (k2.1 + k3.1) * 2.0 + k4.1) * (h / 6.0),
                    );
                }
                assert!((x - y.0).length() < 1e-9, "{e} {b} {s}: {x} vs {}", y.0);
                assert!((v - kin(y.1)).length() < 1e-9, "{e} {b} {s}");
            }
        }
    }

    /// The acceleration at `s = 0` is the Lorentz force's, and the jerk is the derivative
    /// of the acceleration along the motion (central difference in coordinate time).
    #[test]
    fn acceleration_and_jerk() {
        let (q, m, c) = (-0.4, 0.9, 3.0);
        let (x0, v0) = (DVec3::new(0.5, -1.0, 0.0), DVec3::new(-1.2, 1.9, 0.0));
        let (e, b) = (DVec3::new(0.3, 0.9, 0.0), DVec3::new(0.0, 0.0, -0.7));
        let motion = FieldMotion::new(q, m, c, 2.0, x0, v0, e, b).expect("planar");
        let (_, _, _, a) = motion.state(0.0);
        let kin = crate::dynamics::Kinematics::new(m, c);
        let gamma = 1.0 / (1.0 - v0.length_squared() / (c * c)).sqrt();
        let lorentz = kin.acceleration(v0 * (gamma * m), (e + v0.cross(b)) * q);
        assert!(
            (a - lorentz).length() < 1e-13 * lorentz.length(),
            "{a} vs {lorentz}"
        );
        let h = 1e-4;
        let (tp, _, _, ap) = motion.state(h);
        let (tm, _, _, am) = motion.state(-h);
        let numeric = (ap - am) / (tp - tm);
        let jerk = motion.jerk();
        assert!(
            (jerk - numeric).length() < 1e-7 * jerk.length(),
            "{jerk} vs {numeric}"
        );
    }

    /// The snap is the derivative of the jerk along the motion: the jerks of the motions
    /// started from the states a little before and after (central difference in coordinate
    /// time), in magnetic- and electric-dominated fields.
    #[test]
    fn snap() {
        let (q, m, c) = (-0.4, 0.9, 3.0);
        let (x0, v0) = (DVec3::new(0.5, -1.0, 0.0), DVec3::new(-1.2, 1.9, 0.0));
        for (e, b) in [
            (DVec3::new(0.3, 0.9, 0.0), DVec3::new(0.0, 0.0, -0.7)),
            (DVec3::new(2.5, -1.0, 0.0), DVec3::new(0.0, 0.0, 0.1)),
        ] {
            let motion = FieldMotion::new(q, m, c, 2.0, x0, v0, e, b).expect("planar");
            let h = 1e-4;
            let jerk_at = |s: f64| {
                let (t, x, v, _) = motion.state(s);
                let m2 = FieldMotion::new(q, m, c, t, x, v, e, b).expect("planar");
                (t, m2.jerk())
            };
            let ((tp, jp), (tm, jm)) = (jerk_at(h), jerk_at(-h));
            let numeric = (jp - jm) / (tp - tm);
            let snap = motion.snap();
            assert!(
                (snap - numeric).length() < 1e-7 * snap.length(),
                "{snap} vs {numeric}"
            );
        }
    }

    /// The jerk correction's position, velocity and acceleration are each other's
    /// derivatives (central differences, inside the taper, across its end and beyond),
    /// start as `Δȧ τ³/6`, `Δȧ τ²/2`, `Δȧ τ`, and the velocity change stays below
    /// `|Δȧ| T²/8`.
    #[test]
    fn jerk_correction() {
        let corr = JerkCorrection {
            jerk: DVec3::new(0.3, -0.8, 0.0),
            scale: 2.5,
        };
        let h = 1e-5;
        for tau in [-7.0, -2.6, -2.5, -2.4, -1.0, -0.3, 0.0, 0.4, 2.45, 2.5, 3.0] {
            let (x, v, a) = corr.at(tau);
            let (xp, vp, _) = corr.at(tau + h);
            let (xm, vm, _) = corr.at(tau - h);
            assert!(((xp - xm) / (2.0 * h) - v).length() < 1e-8, "{tau}: v");
            assert!(((vp - vm) / (2.0 * h) - a).length() < 1e-8, "{tau}: a");
            assert!(
                v.length() <= corr.jerk.length() * corr.scale.powi(2) / 8.0 * (1.0 + 1e-15),
                "{tau}: |δv|"
            );
            assert!(x.is_finite() && v.is_finite() && a.is_finite());
        }
        let tau = 1e-3;
        let (x, v, a) = corr.at(tau);
        let rel = |got: DVec3, want: DVec3| (got - want).length() / want.length();
        assert!(rel(x, corr.jerk * (tau.powi(3) / 6.0)) < 1e-6);
        assert!(rel(v, corr.jerk * (tau * tau / 2.0)) < 1e-6);
        assert!(rel(a, corr.jerk * tau) < 1e-6);
        // The jerk of the corrected past at the present: Δȧ.
        let (_, _, ap) = corr.at(h);
        let (_, _, am) = corr.at(-h);
        assert!(rel((ap - am) / (2.0 * h), corr.jerk) < 1e-9);
    }

    /// The retarded point on the corrected world line: on the light cone of the observer,
    /// and with the correction's position, velocity and acceleration at its time.
    #[test]
    fn retarded_with_correction() {
        let (q, m, c) = (0.7, 1.3, 5.0);
        let (x0, v0) = (DVec3::new(1.0, 2.0, 0.0), DVec3::new(3.0, -1.5, 0.0));
        let (e, b) = (DVec3::new(0.8, -0.4, 0.0), DVec3::new(0.0, 0.0, 1.1));
        let motion = FieldMotion::new(q, m, c, 0.0, x0, v0, e, b).expect("planar");
        let corr = JerkCorrection {
            jerk: DVec3::new(-2.0, 1.0, 0.0),
            scale: 0.8,
        };
        for x in [DVec3::new(4.0, 1.0, 0.0), DVec3::new(-3.0, 6.0, 0.0)] {
            let plain = motion.retarded(x, 0.3).expect("a retarded point");
            let r = motion
                .retarded_with(x, 0.3, &corr)
                .expect("a retarded point");
            let (_, xs, vs, a_s) = motion.state(r.s);
            let (dx, dv, da) = corr.at(r.t);
            assert!(((r.x - (xs + dx)).length()) < 1e-13);
            assert!(((r.v - (vs + dv)).length()) < 1e-13);
            assert!(((r.a - (a_s + da)).length()) < 1e-12);
            let cone = c * (0.3 - r.t) - (x - r.x).length();
            assert!(cone.abs() < 1e-12, "{cone:e}");
            assert!(
                (r.t - plain.t).abs() > 1e-4,
                "the correction moves the point"
            );
        }
    }
}
