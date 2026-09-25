//! Relativistic equation of motion of a test particle (PHYSICS.md §3).
//!
//! State `(x, p)`, with `γ = sqrt(1 + |p|²/(m²c²))`, `v = p/(γm)`, `dp/dt = q(E + v×B)`.
//! `c = ∞` gives exact Newtonian mechanics (`γ = 1`), which the non-relativistic tests use.
//!
//! For error control the integrator works with the scaled state `(x, p / p_ref)`, so that
//! one tolerance is meaningful for both positions (in grid units) and momenta.

use glam::DVec3;

use crate::field::FieldSolver;
use crate::integrator::OdeSystem;

/// A test particle: rigid, spherically symmetric, non-polarizable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Particle {
    pub charge: f64,
    pub mass: f64,
    pub radius: f64,
}

/// Relativistic kinematics for one particle species in a world with speed of light `c`.
#[derive(Clone, Copy, Debug)]
pub struct Kinematics {
    pub mass: f64,
    pub c: f64,
    /// `1/(m c)²`; zero for `c = ∞`.
    inv_mc_sq: f64,
}

impl Kinematics {
    pub fn new(mass: f64, c: f64) -> Self {
        let mc = mass * c;
        Self {
            mass,
            c,
            inv_mc_sq: 1.0 / (mc * mc),
        }
    }

    pub fn gamma(&self, p: DVec3) -> f64 {
        (1.0 + p.length_squared() * self.inv_mc_sq).sqrt()
    }

    pub fn velocity(&self, p: DVec3) -> DVec3 {
        p / (self.gamma(p) * self.mass)
    }

    /// Kinetic energy `(γ − 1) m c²`, evaluated as `p² / (m (γ + 1))` to avoid
    /// cancellation at low speed (and to remain finite for `c = ∞`).
    pub fn kinetic_energy(&self, p: DVec3) -> f64 {
        p.length_squared() / (self.mass * (self.gamma(p) + 1.0))
    }

    /// Momentum of magnitude corresponding to kinetic energy `t`, along `direction`.
    pub fn momentum_from_kinetic_energy(&self, t: f64, direction: DVec3) -> DVec3 {
        // p² = T² / c² + 2 m T.
        let p2 = t * t * self.inv_mc_sq * self.mass * self.mass + 2.0 * self.mass * t;
        direction.normalize() * p2.sqrt()
    }
}

/// Right-hand side of the equation of motion in the scaled state `(x, p / p_ref)`.
pub struct ParticleOde<F> {
    pub field: F,
    pub charge: f64,
    pub kin: Kinematics,
    pub p_ref: f64,
}

impl<F: FieldSolver> ParticleOde<F> {
    pub fn new(field: F, particle: &Particle, c: f64, p_ref: f64) -> Self {
        Self {
            field,
            charge: particle.charge,
            kin: Kinematics::new(particle.mass, c),
            p_ref,
        }
    }

    pub fn pack(&self, x: DVec3, p: DVec3) -> [f64; 6] {
        let s = p / self.p_ref;
        [x.x, x.y, x.z, s.x, s.y, s.z]
    }

    pub fn position(y: &[f64]) -> DVec3 {
        DVec3::new(y[0], y[1], y[2])
    }

    pub fn momentum(&self, y: &[f64]) -> DVec3 {
        DVec3::new(y[3], y[4], y[5]) * self.p_ref
    }

    /// Force `q(E + v×B)` at state `(x, p)` and time `t`.
    pub fn force(&self, x: DVec3, p: DVec3, t: f64) -> DVec3 {
        let f = self.field.sample(x, t);
        let v = self.kin.velocity(p);
        (f.e + v.cross(f.b)) * self.charge
    }
}

impl<F: FieldSolver> OdeSystem for ParticleOde<F> {
    fn dim(&self) -> usize {
        6
    }

    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        let x = Self::position(y);
        let p = self.momentum(y);
        let v = self.kin.velocity(p);
        let f = self.field.sample(x, t);
        let force = (f.e + v.cross(f.b)) * self.charge;
        let dp = force / self.p_ref;
        dy[0] = v.x;
        dy[1] = v.y;
        dy[2] = v.z;
        dy[3] = dp.x;
        dy[4] = dp.y;
        dy[5] = dp.z;
    }
}
