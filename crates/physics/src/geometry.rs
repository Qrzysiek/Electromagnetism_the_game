//! Geometric primitives: spheres, boxes, and regions with signed distance functions.
//!
//! Every signed distance here is 1-Lipschitz in position, so along a trajectory its time
//! derivative is bounded by the particle speed. Event location relies on this.

use glam::DVec3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    pub center: DVec3,
    pub radius: f64,
}

impl Sphere {
    /// Signed distance from the surface: positive outside, negative inside.
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        (x - self.center).length() - self.radius
    }
}

/// Axis-aligned box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb {
    pub min: DVec3,
    pub max: DVec3,
}

impl Aabb {
    /// Signed distance from the surface: positive outside, negative inside (exact
    /// Euclidean distance outside; distance to the nearest face inside).
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        let center = (self.min + self.max) * 0.5;
        let half = (self.max - self.min) * 0.5;
        let q = (x - center).abs() - half;
        q.max(DVec3::ZERO).length() + q.max_element().min(0.0)
    }
}

/// A region of space (used for the detector B).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Region {
    Sphere(Sphere),
    Box(Aabb),
}

impl Region {
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        match self {
            Region::Sphere(s) => s.signed_distance(x),
            Region::Box(b) => b.signed_distance(x),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_signed_distance() {
        let b = Aabb {
            min: DVec3::new(-1.0, -2.0, -3.0),
            max: DVec3::new(1.0, 2.0, 3.0),
        };
        assert!((b.signed_distance(DVec3::ZERO) + 1.0).abs() < 1e-15);
        assert!((b.signed_distance(DVec3::new(4.0, 0.0, 0.0)) - 3.0).abs() < 1e-15);
        let corner = b.signed_distance(DVec3::new(4.0, 6.0, 3.0));
        assert!((corner - 5.0).abs() < 1e-15);
    }
}
