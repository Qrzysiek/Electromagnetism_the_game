//! External fields from sources outside the arena: uniform stray fields and plane waves
//! (PHYSICS.md §2.3). Both are exact solutions of the vacuum Maxwell equations, so adding
//! them introduces no approximation.

use glam::DVec3;

use crate::field::{FieldSample, FieldSolver};

/// A linearly polarized plane wave in vacuum:
/// `E = E₀ ê cos(ω (t − k̂·x / c) + φ)`, `B = k̂ × E / c`.
///
/// For `c = ∞` it becomes a spatially uniform field oscillating in time, with `B = 0`
/// (the long-wavelength limit, exact in Newtonian electrodynamics).
///
/// `ω = 0` is a static uniform electric field `E₀ ê cos φ` with `B = 0`: the wave vector
/// `k = ω k̂/c` vanishes, so Faraday's law (`k × E = ω B`) no longer ties a magnetic field
/// to it (`k̂` only names the direction ê is perpendicular to). Both `(E, 0)` and
/// `(E, k̂ × E/c)` solve the static equations; the term is the electric field the level
/// format and the interface call a static field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneWave {
    pub amplitude: f64,
    /// Unit propagation direction `k̂`.
    pub direction: DVec3,
    /// Unit polarization `ê`, perpendicular to `k̂`.
    pub polarization: DVec3,
    /// Angular frequency `ω` (0 gives a static uniform electric field `E₀ ê cos φ`).
    pub omega: f64,
    pub phase: f64,
    /// Speed of light of the world (may be infinite).
    pub c: f64,
}

impl PlaneWave {
    /// A wave travelling in the plane z = 0 at angle `angle` from +x, polarized in the
    /// plane along `ẑ × k̂`. Its `B` is along z, so it keeps a particle in the plane.
    pub fn in_plane(amplitude: f64, angle: f64, omega: f64, phase: f64, c: f64) -> Self {
        let (s, co) = libm::sincos(angle);
        Self {
            amplitude,
            direction: DVec3::new(co, s, 0.0),
            polarization: DVec3::new(-s, co, 0.0),
            omega,
            phase,
            c,
        }
    }

    /// Phase `ω (t − k̂·x / c) + φ`.
    pub fn phase_at(&self, x: DVec3, t: f64) -> f64 {
        let retard = if self.c.is_finite() {
            self.direction.dot(x) / self.c
        } else {
            0.0
        };
        self.omega * (t - retard) + self.phase
    }

    pub fn fields(&self, x: DVec3, t: f64) -> (DVec3, DVec3) {
        let e = self.polarization * (self.amplitude * libm::cos(self.phase_at(x, t)));
        let b = if self.c.is_finite() && self.omega != 0.0 {
            self.direction.cross(e) / self.c
        } else {
            DVec3::ZERO
        };
        (e, b)
    }

    /// Vector potential in the gauge `φ = 0`: `A = −(E₀/ω) ê sin(phase)`, so that
    /// `E = −∂A/∂t` and `B = ∇ × A`; for `ω = 0` (static E, no B) `A = −E₀ ê cos(φ) t`.
    /// Used by the tests (canonical momentum).
    pub fn vector_potential(&self, x: DVec3, t: f64) -> DVec3 {
        if self.omega == 0.0 {
            return self.polarization * (-self.amplitude * libm::cos(self.phase) * t);
        }
        self.polarization * (-self.amplitude / self.omega * libm::sin(self.phase_at(x, t)))
    }
}

/// One external field term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum External {
    /// Static uniform fields; the electric part has potential `−E·x`.
    Uniform {
        e: DVec3,
        b: DVec3,
    },
    Wave(PlaneWave),
}

impl External {
    pub fn is_static(&self) -> bool {
        match self {
            External::Uniform { .. } => true,
            External::Wave(w) => w.omega == 0.0,
        }
    }
}

impl FieldSolver for External {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        match self {
            External::Uniform { e, b } => FieldSample {
                e: *e,
                b: *b,
                phi: -e.dot(x),
            },
            External::Wave(w) => {
                let (e, b) = w.fields(x, t);
                // A static (ω = 0) wave term is the uniform electric field E₀ ê cos φ
                // (`PlaneWave`), with its potential −E·x.
                let phi = if w.omega == 0.0 { -e.dot(x) } else { 0.0 };
                FieldSample { e, b, phi }
            }
        }
    }

    fn is_static(&self) -> bool {
        External::is_static(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The wave satisfies the vacuum Maxwell equations (checked by central differences):
    /// ∇·E = 0, ∇·B = 0, ∇×E = −∂B/∂t, ∇×B = (1/c²) ∂E/∂t.
    #[test]
    #[allow(clippy::float_cmp)] // exact zeros by construction
    fn plane_wave_satisfies_maxwell_equations() {
        let c = 3.0;
        let w = PlaneWave::in_plane(2.0, 0.7, 1.3, 0.4, c);
        let x = DVec3::new(0.3, -1.2, 0.5);
        let t = 0.9;
        let h = 1e-4;
        let e = |x: DVec3, t: f64| w.fields(x, t).0;
        let b = |x: DVec3, t: f64| w.fields(x, t).1;
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
        let scale = 2.0 * 1.3 / c;
        assert!(div(&e).abs() < 1e-7 * scale);
        assert!(div(&b).abs() < 1e-7 * scale);
        assert!((curl(&e) + dt(&b)).length() < 1e-7 * scale);
        assert!((curl(&b) - dt(&e) / (c * c)).length() < 1e-7 * scale);
        // In-plane wave: E in the plane, B perpendicular to it.
        assert_eq!(e(x, t).z, 0.0);
        assert_eq!(b(x, t).x, 0.0);
        assert_eq!(b(x, t).y, 0.0);
    }

    #[test]
    fn vector_potential_generates_the_fields() {
        for omega in [2.0, 0.0] {
            let w = PlaneWave::in_plane(1.5, -0.4, omega, 1.1, 4.0);
            let (x, t, h) = (DVec3::new(1.0, 0.5, 0.0), 0.3, 1e-5);
            let a = |x: DVec3, t: f64| w.vector_potential(x, t);
            let e = -(a(x, t + h) - a(x, t - h)) / (2.0 * h);
            assert!((e - w.fields(x, t).0).length() < 1e-8);
            let d = |axis: DVec3| (a(x + axis * h, t) - a(x - axis * h, t)) / (2.0 * h);
            let (dx, dy) = (d(DVec3::X), d(DVec3::Y));
            let bz = dx.y - dy.x;
            assert!((bz - w.fields(x, t).1.z).abs() < 1e-8);
        }
    }

    /// `ω = 0`: a static uniform electric field (no magnetic field), with potential −E·x.
    #[test]
    #[allow(clippy::float_cmp)] // exact zeros by construction
    fn static_wave_term_is_a_uniform_electric_field() {
        let w = External::Wave(PlaneWave::in_plane(2.0, 0.6, 0.0, 0.3, 4.0));
        let (x, y) = (DVec3::new(1.0, -2.0, 0.0), DVec3::new(-3.0, 0.5, 0.0));
        let (a, b) = (w.sample(x, 0.4), w.sample(y, 7.0));
        assert_eq!(a.e, b.e);
        assert_eq!(a.b, DVec3::ZERO);
        assert!((a.e.length() - 2.0 * libm::cos(0.3)).abs() < 1e-15);
        assert!(((a.phi - b.phi) - (y - x).dot(a.e)).abs() < 1e-14);
    }
}
