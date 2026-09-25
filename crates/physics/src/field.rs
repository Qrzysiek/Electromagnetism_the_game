//! Field sources (PHYSICS.md §2).

use glam::DVec3;

/// Fields and potential at a point. Internal units: `k = 1/(4πε₀) = 1`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FieldSample {
    pub e: DVec3,
    pub b: DVec3,
    /// Electrostatic potential (only meaningful for static electric fields).
    pub phi: f64,
}

/// Source of external fields acting on test particles.
pub trait FieldSolver {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample;
}

impl<F: FieldSolver + ?Sized> FieldSolver for &F {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        (**self).sample(x, t)
    }
}

/// A fixed charge: a rigid sphere with a spherically symmetric charge distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FixedCharge {
    pub position: DVec3,
    pub charge: f64,
    pub radius: f64,
}

/// Exact Coulomb superposition of fixed charges, valid outside every sphere (the only
/// region particles can reach). Summed in the order the charges are stored.
#[derive(Clone, Debug, Default)]
pub struct Coulomb {
    positions: Vec<DVec3>,
    charges: Vec<f64>,
}

impl Coulomb {
    pub fn new(charges: &[FixedCharge]) -> Self {
        Self {
            positions: charges.iter().map(|c| c.position).collect(),
            charges: charges.iter().map(|c| c.charge).collect(),
        }
    }
}

impl FieldSolver for Coulomb {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        let mut e = DVec3::ZERO;
        let mut phi = 0.0;
        for (&r, &q) in self.positions.iter().zip(&self.charges) {
            let d = x - r;
            let inv_r = 1.0 / d.length();
            let q_inv_r = q * inv_r;
            phi += q_inv_r;
            e += d * (q_inv_r * inv_r * inv_r);
        }
        FieldSample {
            e,
            b: DVec3::ZERO,
            phi,
        }
    }
}

/// Uniform static electric field, for tests with analytic solutions.
#[derive(Clone, Copy, Debug)]
pub struct UniformElectric {
    pub e: DVec3,
}

impl FieldSolver for UniformElectric {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        FieldSample {
            e: self.e,
            b: DVec3::ZERO,
            phi: -self.e.dot(x),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coulomb_field_is_minus_gradient_of_potential() {
        let f = Coulomb::new(&[
            FixedCharge {
                position: DVec3::new(0.3, -1.0, 0.2),
                charge: 2.0,
                radius: 0.1,
            },
            FixedCharge {
                position: DVec3::new(-1.5, 0.5, -0.7),
                charge: -3.0,
                radius: 0.1,
            },
        ]);
        let x = DVec3::new(0.9, 0.4, -0.3);
        let h = 1e-5;
        let s = f.sample(x, 0.0);
        for axis in [DVec3::X, DVec3::Y, DVec3::Z] {
            let grad =
                (f.sample(x + axis * h, 0.0).phi - f.sample(x - axis * h, 0.0).phi) / (2.0 * h);
            assert!((s.e.dot(axis) + grad).abs() < 1e-8);
        }
    }
}
