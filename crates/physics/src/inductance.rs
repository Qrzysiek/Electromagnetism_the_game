//! Inductances of the game's coils from their exact vector potentials (PHYSICS.md §2.8,
//! tests Z2–Z3), in units of `μ₀/4π` (the game's `1/c²`): the Neumann integral
//! `N = ∮∮ dl₁·dl₂/|x₁ − x₂|`, which a coil's vector potential per unit strength
//! (`κ = μ₀I/4π = 1`) gives as `∮ A·dl` around the other coil.

use std::f64::consts::TAU;

use glam::DVec3;

use crate::magnetic::{CircularLoop, PolygonCoil};

/// Relative agreement of two successive results (the point count doubled) at which the
/// trapezoidal rule stops; the rule converges exponentially for these smooth periodic
/// integrands, so the finer result is far more accurate than their difference.
const CONVERGED: f64 = 1e-14;

/// An orthonormal pair spanning the plane orthogonal to the unit vector `n`, with
/// `u × v = n` (for `n = ẑ`: x̂ and ŷ).
fn plane_basis(n: DVec3) -> (DVec3, DVec3) {
    let axis = if n.x.abs() <= n.y.abs() && n.x.abs() <= n.z.abs() {
        DVec3::X
    } else if n.y.abs() <= n.z.abs() {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let u = (axis - n * n.dot(axis)).normalize();
    (u, n.cross(u))
}

/// `∮ A·dl` around the circle of `radius` about `center` in the plane orthogonal to the
/// unit vector `normal`, counter-clockwise about it: the trapezoidal rule on `n` points.
pub fn around_circle(
    a: impl Fn(DVec3) -> DVec3,
    center: DVec3,
    normal: DVec3,
    radius: f64,
    n: usize,
) -> f64 {
    let (u, v) = plane_basis(normal);
    #[allow(clippy::cast_precision_loss)] // point counts are small
    let step = TAU / n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let t = step * i as f64;
        let (s, c) = (libm::sin(t), libm::cos(t));
        let x = center + (u * c + v * s) * radius;
        sum += a(x).dot((v * c - u * s) * radius);
    }
    sum * step
}

/// `f(n)` for n = 64, 128, … until two successive results agree to `CONVERGED` (at most
/// 2¹⁶ points); the last result.
fn doubled(f: impl Fn(usize) -> f64) -> f64 {
    let mut n = 64;
    let mut last = f(n);
    while n < 1 << 16 {
        n *= 2;
        let next = f(n);
        if (next - last).abs() <= CONVERGED * next.abs() {
            return next;
        }
        last = next;
    }
    last
}

/// Mutual inductance of two circular coils (`μ₀/4π` units): the flux of `a`'s field
/// linked by `b`, `∮_b A_a·dl` (Jackson §5.17; test Z2). Symmetric, and independent of
/// the coils' strengths.
pub fn mutual_circles(a: &CircularLoop, b: &CircularLoop) -> f64 {
    doubled(|n| {
        around_circle(
            |x| a.unit_vector_potential(x),
            b.center,
            b.normal,
            b.radius,
            n,
        )
    })
}

/// Self-inductance of a circular coil (`μ₀/4π` units) whose current flows uniformly over
/// its wire's surface, as on the perfect conductors of the game: the coil's own flux
/// through the circles on that surface, `(ρ, z) = (a + b cos ψ, b sin ψ)`, averaged over ψ
/// (Jackson §5.17 B and Pr. 5.32; test Z3:
/// `4πa[ln(8a/b) − 2 + (b/a)²((1/8) ln(8a/b) + 1/16) + …]`).
pub fn ring_self(l: &CircularLoop) -> f64 {
    let (a, b) = (l.radius, l.wire_radius);
    let (u, v) = plane_basis(l.normal);
    doubled(|n| {
        let mut sum = 0.0;
        for i in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let psi = TAU * i as f64 / n as f64;
            let rho = a + b * libm::cos(psi);
            let x = l.center + u * rho + l.normal * (b * libm::sin(psi));
            // By symmetry A is azimuthal and constant around the circle through x.
            sum += TAU * rho * l.unit_vector_potential(x).dot(v);
        }
        #[allow(clippy::cast_precision_loss)]
        let mean = sum / n as f64;
        mean
    })
}

/// Gauss–Legendre, 8 points on [−1, 1]: nodes and weights.
const GL8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

/// `∫ f(x)·dl` along the straight segment from `a` to `b`: composite 8-point Gauss–
/// Legendre on 1, 2, 4, … panels until two results agree to `CONVERGED` (at most 4096
/// panels; the integrands here are smooth on the segment).
fn along_segment(f: &impl Fn(DVec3) -> DVec3, a: DVec3, b: DVec3) -> f64 {
    let rule = |panels: usize| {
        let d = (b - a) / f64::from(u32::try_from(panels).expect("few panels"));
        let mut sum = 0.0;
        for k in 0..panels {
            let mid = a + d * (f64::from(u32::try_from(k).expect("few panels")) + 0.5);
            for &(x, w) in &GL8 {
                sum += w * (f(mid + d * (0.5 * x)).dot(d) + f(mid - d * (0.5 * x)).dot(d));
            }
        }
        0.5 * sum
    };
    let mut panels = 1;
    let mut last = rule(panels);
    while panels < 4096 {
        panels *= 2;
        let next = rule(panels);
        if (next - last).abs() <= CONVERGED * next.abs().max(1e-300) {
            return next;
        }
        last = next;
    }
    last
}

/// The segments of a closed polygon `(start, end)`, without those of no length (a
/// repeated vertex).
fn segments(p: &PolygonCoil) -> Vec<(DVec3, DVec3)> {
    let n = p.vertices.len();
    (0..n)
        .map(|i| (p.vertices[i], p.vertices[(i + 1) % n]))
        .filter(|(a, b)| a != b)
        .collect()
}

/// Self-inductance of a polygonal coil (`μ₀/4π` units) whose current flows uniformly over
/// its wire's surface, as `ring_self` for a ring: the coil's own flux (its centre line's
/// exact potential) through the closed curves on that surface, averaged over the angle ψ
/// around the wire. Such a curve is the polygon offset in the plane by `b cos ψ` (its
/// corners mitred: each side's offset line, meeting the next) at the height `b sin ψ`;
/// for a circle it is `ring_self`'s circle. (Test Z14: rectangles against Wolfram, and
/// fine regular polygons tend to the ring.)
pub fn polygon_self(p: &PolygonCoil) -> f64 {
    let s = segments(p);
    let n = s.len();
    let b = p.wire_radius;
    // In-plane unit normals of the sides (left of the current's direction).
    let normals: Vec<DVec3> = s
        .iter()
        .map(|&(a, c)| DVec3::Z.cross((c - a).normalize()))
        .collect();
    let potential = |x: DVec3| p.unit_vector_potential(x);
    doubled(|m| {
        let mut sum = 0.0;
        for i in 0..m {
            #[allow(clippy::cast_precision_loss)]
            let psi = TAU * i as f64 / m as f64;
            let (d, h) = (b * libm::cos(psi), b * libm::sin(psi));
            // Vertex k of the offset curve: where the offset lines of sides k − 1 and k
            // meet, `vertex + d (n₁ + n₂)/(1 + n₁·n₂)`.
            let q: Vec<DVec3> = (0..n)
                .map(|k| {
                    let (n1, n2) = (normals[(k + n - 1) % n], normals[k]);
                    s[k].0 + (n1 + n2) * (d / (1.0 + n1.dot(n2))) + DVec3::Z * h
                })
                .collect();
            sum += (0..n)
                .map(|k| along_segment(&potential, q[k], q[(k + 1) % n]))
                .sum::<f64>();
        }
        #[allow(clippy::cast_precision_loss)]
        let mean = sum / m as f64;
        mean
    })
}

/// Mutual inductance of a polygonal coil and another coil whose vector potential per unit
/// strength is `a_other` (`μ₀/4π` units): `∮_p A_other·dl` along the polygon's sides.
pub fn mutual_polygon(p: &PolygonCoil, a_other: impl Fn(DVec3) -> DVec3) -> f64 {
    segments(p)
        .iter()
        .map(|&(a, b)| along_segment(&a_other, a, b))
        .sum()
}

/// Mutual inductance of a polygonal coil and a circular one (`μ₀/4π` units).
pub fn mutual_polygon_circle(p: &PolygonCoil, c: &CircularLoop) -> f64 {
    mutual_polygon(p, |x| c.unit_vector_potential(x))
}

/// Mutual inductance of two polygonal coils (`μ₀/4π` units).
pub fn mutual_polygons(p: &PolygonCoil, q: &PolygonCoil) -> f64 {
    mutual_polygon(p, |x| q.unit_vector_potential(x))
}
