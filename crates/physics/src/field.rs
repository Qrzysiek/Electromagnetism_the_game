//! Field sources (PHYSICS.md §2).

use glam::DVec3;

use crate::external::External;
use crate::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};

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

    /// True if the fields do not depend on time. Only then is `phi` a potential of `E`
    /// and the energy `(γ−1)mc² + qφ` conserved.
    fn is_static(&self) -> bool {
        true
    }
}

impl<F: FieldSolver + ?Sized> FieldSolver for &F {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        (**self).sample(x, t)
    }

    fn is_static(&self) -> bool {
        (**self).is_static()
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

impl Coulomb {
    /// Positions and charges, in summation order.
    pub fn charges(&self) -> impl Iterator<Item = (DVec3, f64)> + '_ {
        self.positions
            .iter()
            .copied()
            .zip(self.charges.iter().copied())
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

/// All sources of a level: fixed charges, magnets and coils, plus external fields
/// (uniform stray fields and plane waves, PHYSICS.md §2.3). Summed in this order.
#[derive(Clone, Debug, Default)]
pub struct LevelField {
    pub coulomb: Coulomb,
    pub dipoles: Vec<MagneticDipole>,
    pub loops: Vec<CircularLoop>,
    pub polygons: Vec<PolygonCoil>,
    pub external: Vec<External>,
}

impl LevelField {
    pub fn magnetic(&self, x: DVec3) -> DVec3 {
        let mut b = DVec3::ZERO;
        for d in &self.dipoles {
            b += d.field(x);
        }
        for l in &self.loops {
            b += l.field(x);
        }
        for p in &self.polygons {
            b += p.field(x);
        }
        b
    }
}

impl FieldSolver for LevelField {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        let mut s = self.coulomb.sample(x, t);
        s.b = self.magnetic(x);
        for ext in &self.external {
            let f = ext.sample(x, t);
            s.e += f.e;
            s.b += f.b;
            s.phi += f.phi;
        }
        s
    }

    fn is_static(&self) -> bool {
        self.external.iter().all(External::is_static)
    }
}

/// Uniform static electric and magnetic fields, for tests with analytic solutions.
#[derive(Clone, Copy, Debug)]
pub struct UniformFields {
    pub e: DVec3,
    pub b: DVec3,
}

impl FieldSolver for UniformFields {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        FieldSample {
            e: self.e,
            b: self.b,
            phi: -self.e.dot(x),
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
