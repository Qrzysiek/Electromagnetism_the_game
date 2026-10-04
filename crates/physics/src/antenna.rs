//! Oscillating electric dipole (a small antenna) with its exact retarded fields
//! (PHYSICS.md §2.4).
//!
//! With `p(t) = p₀ cos(ωt + φ)`, unit vector `n` from the dipole to the field point at
//! distance `r`, and everything evaluated at the retarded time `t_r = t − r/c`
//! (internal units `k = 1`, `μ₀/4π = 1/c²`; Jackson, Classical Electrodynamics §9.2):
//!
//! ```text
//! E = [3n(n·p) − p] / r³ + [3n(n·ṗ) − ṗ] / (c r²) + [n × (n × p̈)] / (c² r)
//! B = (1/c²) [ṗ × n / r² + p̈ × n / (c r)]
//! ```
//!
//! This solves the vacuum Maxwell equations exactly everywhere outside the source point:
//! near (quasi-static), intermediate and radiation zones alike. For `c = ∞` it reduces to
//! the electrostatic dipole field of the instantaneous `p(t)` with `B = 0`.
//!
//! The potentials (Lorenz gauge, `E = −∇φ − ∂A/∂t`, `B = ∇ × A`) are
//! `φ = n·p / r² + n·ṗ / (c r)` and `A = ṗ / (c² r)`, retarded. The sample's `phi` is this
//! φ: for `ω = 0` the electrostatic dipole potential `n·p / r²`, which with the field makes
//! `T + qφ` the conserved energy of a static antenna's flight.

use glam::DVec3;

use crate::field::FieldSample;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OscillatingDipole {
    pub position: DVec3,
    /// Amplitude vector `p₀` (direction and size of the dipole moment).
    pub amplitude: DVec3,
    pub omega: f64,
    pub phase: f64,
    /// Speed of light of the world (may be infinite).
    pub c: f64,
    /// Radius of the solid antenna body (an obstacle).
    pub radius: f64,
}

impl OscillatingDipole {
    /// `(p, ṗ, p̈)` at time `t`.
    pub fn moment(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let (s, c) = libm::sincos(self.omega * t + self.phase);
        let w = self.omega;
        (
            self.amplitude * c,
            self.amplitude * (-w * s),
            self.amplitude * (-w * w * c),
        )
    }

    pub fn fields(&self, x: DVec3, t: f64) -> FieldSample {
        let d = x - self.position;
        let r = d.length();
        let n = d / r;
        if self.c.is_infinite() {
            let (p, _, _) = self.moment(t);
            return FieldSample {
                e: (n * (3.0 * n.dot(p)) - p) / (r * r * r),
                b: DVec3::ZERO,
                phi: n.dot(p) / (r * r),
            };
        }
        let c = self.c;
        let (p, pd, pdd) = self.moment(t - r / c);
        let e = (n * (3.0 * n.dot(p)) - p) / (r * r * r)
            + (n * (3.0 * n.dot(pd)) - pd) / (c * r * r)
            + n.cross(n.cross(pdd)) / (c * c * r);
        let b = (pd.cross(n) / (r * r) + pdd.cross(n) / (c * r)) / (c * c);
        let phi = n.dot(p) / (r * r) + n.dot(pd) / (c * r);
        FieldSample { e, b, phi }
    }

    /// The Lorenz-gauge vector potential `A = ṗ / (c² r)` (retarded; 0 for `c = ∞`).
    pub fn vector_potential(&self, x: DVec3, t: f64) -> DVec3 {
        if self.c.is_infinite() {
            return DVec3::ZERO;
        }
        let r = (x - self.position).length();
        let (_, pd, _) = self.moment(t - r / self.c);
        pd / (self.c * self.c * r)
    }

    /// Time-averaged radiated power (Larmor): `⟨P⟩ = p₀² ω⁴ / (3 c³)`.
    pub fn larmor_power(&self) -> f64 {
        self.amplitude.length_squared() * self.omega.powi(4) / (3.0 * self.c.powi(3))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dipole(c: f64) -> OscillatingDipole {
        OscillatingDipole {
            position: DVec3::new(0.5, -0.3, 0.1),
            amplitude: DVec3::new(0.8, 0.6, -0.3),
            omega: 1.7,
            phase: 0.4,
            c,
            radius: 0.1,
        }
    }

    /// All four vacuum Maxwell equations by central differences, in the near,
    /// intermediate and radiation zones (λ = 2πc/ω ≈ 11).
    #[test]
    fn oscillating_dipole_satisfies_maxwell_equations() {
        let c = 3.0;
        let a = dipole(c);
        for x in [
            DVec3::new(1.3, 0.2, -0.4),
            DVec3::new(-3.0, 4.0, 1.0),
            DVec3::new(20.0, -35.0, 12.0),
        ] {
            let t = 0.7;
            // Step small against both the distance and the wavelength. The residuals
            // shrink as h² (checked: 16× per 4× smaller h), i.e. they are pure
            // difference error; at this h they are ≤ 7e-8.
            let lambda = 2.0 * std::f64::consts::PI * c / a.omega;
            let h = 2.5e-5 * (x - a.position).length().min(lambda);
            let e = |x: DVec3, t: f64| a.fields(x, t).e;
            let b = |x: DVec3, t: f64| a.fields(x, t).b;
            let d = |f: &dyn Fn(DVec3, f64) -> DVec3, axis: DVec3| {
                (f(x + axis * h, t) - f(x - axis * h, t)) / (2.0 * h)
            };
            let dt = |f: &dyn Fn(DVec3, f64) -> DVec3| {
                let k = 1e-4;
                (f(x, t + k) - f(x, t - k)) / (2.0 * k)
            };
            let curl = |f: &dyn Fn(DVec3, f64) -> DVec3| {
                let (dx, dy, dz) = (d(f, DVec3::X), d(f, DVec3::Y), d(f, DVec3::Z));
                DVec3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x)
            };
            let div = |f: &dyn Fn(DVec3, f64) -> DVec3| {
                d(f, DVec3::X).x + d(f, DVec3::Y).y + d(f, DVec3::Z).z
            };
            // Scales: the size of the terms being compared.
            let r = (x - a.position).length();
            let se = curl(&e).length().max(dt(&b).length());
            let sb = curl(&b).length().max(dt(&e).length() / (c * c));
            assert!(div(&e).abs() < 1e-6 * e(x, t).length() / r);
            assert!(div(&b).abs() < 1e-6 * b(x, t).length() / r);
            let faraday = (curl(&e) + dt(&b)).length() / se;
            let ampere = (curl(&b) - dt(&e) / (c * c)).length() / sb;
            assert!(faraday < 1e-6, "Faraday {faraday:.2e} at {x}");
            assert!(ampere < 1e-6, "Ampère {ampere:.2e} at {x}");
        }
    }

    /// The potentials generate the fields: `E = −∇φ − ∂A/∂t`, `B = ∇ × A` (central
    /// differences), oscillating and static (ω = 0, where φ is the electrostatic dipole
    /// potential), at finite and infinite c.
    #[test]
    fn potentials_generate_the_fields() {
        for (c, omega) in [
            (3.0, 1.7),
            (3.0, 0.0),
            (f64::INFINITY, 1.7),
            (f64::INFINITY, 0.0),
        ] {
            let a = OscillatingDipole {
                omega,
                c,
                ..dipole(c)
            };
            for x in [DVec3::new(1.3, 0.2, -0.4), DVec3::new(-3.0, 4.0, 1.0)] {
                let t = 0.7;
                let h = 1e-5;
                let phi = |x: DVec3| a.fields(x, t).phi;
                let pot = |x: DVec3, t: f64| a.vector_potential(x, t);
                let grad = DVec3::new(
                    phi(x + DVec3::X * h) - phi(x - DVec3::X * h),
                    phi(x + DVec3::Y * h) - phi(x - DVec3::Y * h),
                    phi(x + DVec3::Z * h) - phi(x - DVec3::Z * h),
                ) / (2.0 * h);
                let da_dt = (pot(x, t + h) - pot(x, t - h)) / (2.0 * h);
                let d = |axis: DVec3| (pot(x + axis * h, t) - pot(x - axis * h, t)) / (2.0 * h);
                let (dx, dy, dz) = (d(DVec3::X), d(DVec3::Y), d(DVec3::Z));
                let curl = DVec3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x);
                let f = a.fields(x, t);
                let e_err = (-grad - da_dt - f.e).length() / f.e.length();
                assert!(e_err < 1e-8, "c {c}, ω {omega}, {x}: E {e_err:.2e}");
                let b_err = (curl - f.b).length();
                assert!(b_err <= 1e-8 * f.e.length() / c.min(1e300), "B {b_err:.2e}");
            }
        }
    }

    /// `c = ∞`: the static dipole field of the instantaneous moment.
    #[test]
    fn newtonian_limit_is_the_instantaneous_static_dipole() {
        let a = dipole(f64::INFINITY);
        let far = dipole(1e9);
        let x = DVec3::new(2.0, 1.0, -0.5);
        let (s, f) = (a.fields(x, 0.3), far.fields(x, 0.3));
        assert!((s.e - f.e).length() < 1e-8 * s.e.length());
        assert_eq!(s.b, DVec3::ZERO);
    }
}
