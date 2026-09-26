//! Relativistic equation of motion of a test particle (PHYSICS.md §3).
//!
//! State `(x, p)`, with `γ = sqrt(1 + |p|²/(m²c²))`, `v = p/(γm)`, `dp/dt = q(E + v×B)`.
//! `c = ∞` gives exact Newtonian mechanics (`γ = 1`), which the non-relativistic tests use.
//!
//! For error control the integrator works with the scaled state `(x, p / p_ref)`, so that
//! one tolerance is meaningful for both positions (in grid units) and momenta.
//!
//! Optionally the Landau–Lifshitz radiation-reaction force is added (PHYSICS.md §3.1).
//! The state then has a 7th component, the work done by that force divided by `e_ref`, so
//! the energy balance can be checked to integration accuracy.

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
    /// Include the Landau–Lifshitz radiation-reaction force (no effect for `c = ∞`).
    pub radiation_reaction: bool,
    /// Energy unit of the radiation-work state component.
    pub e_ref: f64,
}

impl<F: FieldSolver> ParticleOde<F> {
    pub fn new(field: F, particle: &Particle, c: f64, p_ref: f64) -> Self {
        Self {
            field,
            charge: particle.charge,
            kin: Kinematics::new(particle.mass, c),
            p_ref,
            radiation_reaction: false,
            e_ref: 1.0,
        }
    }

    /// Enables radiation reaction (if `c` is finite), with energy unit `e_ref`.
    #[must_use]
    pub fn with_radiation_reaction(mut self, on: bool, e_ref: f64) -> Self {
        self.radiation_reaction = on && self.kin.c.is_finite();
        self.e_ref = e_ref;
        self
    }

    pub fn pack(&self, x: DVec3, p: DVec3) -> Vec<f64> {
        let s = p / self.p_ref;
        let mut y = vec![x.x, x.y, x.z, s.x, s.y, s.z];
        if self.radiation_reaction {
            y.push(0.0);
        }
        y
    }

    /// Work done so far by the radiation-reaction force (0 without it).
    pub fn radiation_work(&self, y: &[f64]) -> f64 {
        if self.radiation_reaction {
            y[6] * self.e_ref
        } else {
            0.0
        }
    }

    /// Landau–Lifshitz radiation-reaction force (Classical Theory of Fields §76, converted
    /// to force `q(E + v×B)` with `k = 1`, `μ₀/4π = 1/c²`):
    ///
    /// ```text
    /// f = (2q³/3mc³) γ [DE/Dt + v × DB/Dt]
    ///   + (2q⁴/3m²c⁴) [c E×B + c B×(B×v) + E (v·E)/c]
    ///   − (2q⁴/3m²c⁵) γ² v [(E + v×B)² − (E·v)²/c²]
    /// ```
    ///
    /// `D/Dt = ∂/∂t + v·∇` is the derivative along the world line, evaluated by a central
    /// difference over a displacement of 1e-5 cells (its relative error, ~1e-10, is far
    /// below the O(τ₀) accuracy of the Landau–Lifshitz approximation itself).
    pub fn radiation_reaction_force(&self, x: DVec3, p: DVec3, t: f64) -> DVec3 {
        let c = self.kin.c;
        if !c.is_finite() {
            return DVec3::ZERO;
        }
        let (q, m) = (self.charge, self.kin.mass);
        let gamma = self.kin.gamma(p);
        let v = self.kin.velocity(p);
        let f = self.field.sample(x, t);
        let (e, b) = (f.e, f.b);
        let h = 1e-5 / v.length().max(1.0);
        let fp = self.field.sample(x + v * h, t + h);
        let fm = self.field.sample(x - v * h, t - h);
        let de = (fp.e - fm.e) / (2.0 * h);
        let db = (fp.b - fm.b) / (2.0 * h);
        let a1 = 2.0 * q * q * q / (3.0 * m * c * c * c);
        let a2 = 2.0 * q.powi(4) / (3.0 * m * m * c.powi(4));
        let lorentz = e + v.cross(b);
        let ev = e.dot(v);
        (de + v.cross(db)) * (a1 * gamma)
            + (e.cross(b) * c + b.cross(b.cross(v)) * c + e * (ev / c)) * a2
            - v * (a2 / c * gamma * gamma * (lorentz.length_squared() - ev * ev / (c * c)))
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
        if self.radiation_reaction { 7 } else { 6 }
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
        if self.radiation_reaction {
            let f_rr = self.radiation_reaction_force(x, p, t);
            let dp_rr = f_rr / self.p_ref;
            dy[3] += dp_rr.x;
            dy[4] += dp_rr.y;
            dy[5] += dp_rr.z;
            dy[6] = v.dot(f_rr) / self.e_ref;
        }
    }
}
