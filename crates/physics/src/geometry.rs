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

/// Solid torus: a circular wire of radius `minor` bent into a ring of radius `major`
/// around `center`, in the plane orthogonal to the unit vector `normal`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torus {
    pub center: DVec3,
    pub normal: DVec3,
    pub major: f64,
    pub minor: f64,
}

impl Torus {
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        let d = x - self.center;
        let axial = d.dot(self.normal);
        let radial = (d - self.normal * axial).length();
        libm::hypot(radial - self.major, axial) - self.minor
    }
}

/// Capsule: all points within `radius` of the segment from `a` to `b` (a straight wire).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Capsule {
    pub a: DVec3,
    pub b: DVec3,
    pub radius: f64,
}

impl Capsule {
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        let ab = self.b - self.a;
        let l2 = ab.length_squared();
        // A zero-length segment is a sphere (0/0 gave NaN, and the crossing search on a
        // NaN event value hung the flight: a polygon coil with a repeated vertex).
        let t = if l2 > 0.0 {
            ((x - self.a).dot(ab) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (x - (self.a + ab * t)).length() - self.radius
    }
}

/// A box with one axis along z, rotated in the plane (an electrode), inflated by
/// `margin` (the corners are rounded, as for a Minkowski sum with a sphere).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedBox {
    pub center: DVec3,
    /// Unit in-plane axis of the half-length.
    pub axis: DVec3,
    pub half: DVec3,
    pub margin: f64,
}

impl OrientedBox {
    /// Exact Euclidean signed distance (1-Lipschitz).
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        let d = x - self.center;
        let v = DVec3::new(-self.axis.y, self.axis.x, 0.0);
        let local = DVec3::new(d.dot(self.axis), d.dot(v), d.z);
        let q = local.abs() - self.half;
        let outside = q.max(DVec3::ZERO).length();
        let inside = q.x.max(q.y).max(q.z).min(0.0);
        outside + inside - self.margin
    }
}

/// A solid obstacle: touching it loses the particle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Shape {
    Sphere(Sphere),
    Torus(Torus),
    Capsule(Capsule),
    Box(OrientedBox),
}

impl Shape {
    /// Signed distance, 1-Lipschitz in position for every variant.
    pub fn signed_distance(&self, x: DVec3) -> f64 {
        match self {
            Shape::Sphere(s) => s.signed_distance(x),
            Shape::Torus(t) => t.signed_distance(x),
            Shape::Capsule(c) => c.signed_distance(x),
            Shape::Box(b) => b.signed_distance(x),
        }
    }
}

impl From<Sphere> for Shape {
    fn from(s: Sphere) -> Self {
        Shape::Sphere(s)
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

    #[test]
    fn torus_and_capsule_distances() {
        let t = Torus {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major: 5.0,
            minor: 0.2,
        };
        assert!((t.signed_distance(DVec3::new(0.0, 7.0, 0.0)) - 1.8).abs() < 1e-15);
        assert!((t.signed_distance(DVec3::new(5.0, 0.0, 1.0)) - 0.8).abs() < 1e-15);
        assert!((t.signed_distance(DVec3::ZERO) - 4.8).abs() < 1e-15);
        let c = Capsule {
            a: DVec3::ZERO,
            b: DVec3::new(4.0, 0.0, 0.0),
            radius: 0.1,
        };
        assert!((c.signed_distance(DVec3::new(2.0, 3.0, 0.0)) - 2.9).abs() < 1e-15);
        assert!((c.signed_distance(DVec3::new(7.0, 4.0, 0.0)) - 4.9).abs() < 1e-15);
    }
}

#[cfg(test)]
mod capsule_tests {
    use super::*;

    /// A zero-length capsule is a sphere (regression: NaN hung the crossing search).
    #[test]
    fn zero_length_capsule_is_a_sphere() {
        let c = Capsule {
            a: DVec3::new(1.0, 2.0, 0.0),
            b: DVec3::new(1.0, 2.0, 0.0),
            radius: 0.1,
        };
        let d = c.signed_distance(DVec3::new(4.0, 6.0, 0.0));
        assert!((d - 4.9).abs() < 1e-12, "{d}");
    }
}
