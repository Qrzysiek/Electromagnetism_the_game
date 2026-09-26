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
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneWave {
    pub amplitude: f64,
    /// Unit propagation direction `k̂`.
    pub direction: DVec3,
    /// Unit polarization `ê`, perpendicular to `k̂`.
    pub polarization: DVec3,
    /// Angular frequency `ω` (0 gives a static uniform field `E₀ ê cos φ`).
    pub omega: f64,
    pub phase: f64,
    /// Speed of light of the world (may be infinite).
    pub c: f64,
}

impl PlaneWave {
    /// A wave travelling in the plane z = 0 at angle `angle` from +x, polarized in the
    /// plane along `ẑ × k̂`. Its `B` is along z, so it keeps a particle in the plane.
    pub fn in_plane(amplitude: f64, angle: f64, omega: f64, phase: f64, c: f64) -> Self {
        let (s, co) = angle.sin_cos();
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
        let e = self.polarization * (self.amplitude * self.phase_at(x, t).cos());
        let b = if self.c.is_finite() {
            self.direction.cross(e) / self.c
        } else {
            DVec3::ZERO
        };
        (e, b)
    }

    /// Vector potential in the gauge `φ = 0`: `A = −(E₀/ω) ê sin(phase)`, so that
    /// `E = −∂A/∂t` and `B = ∇ × A`. Used by the tests (canonical momentum).
    pub fn vector_potential(&self, x: DVec3, t: f64) -> DVec3 {
        self.polarization * (-self.amplitude / self.omega * self.phase_at(x, t).sin())
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
                // A static (ω = 0) wave term is a uniform field E₀ ê cos φ (with B = 0,
                // since k·x/c enters only multiplied by ω).
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
        let w = PlaneWave::in_plane(1.5, -0.4, 2.0, 1.1, 4.0);
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
