//! Inductances of the game's coils from their exact vector potentials (PHYSICS.md §2.8,
//! tests Z2–Z3), in units of `μ₀/4π` (the game's `1/c²`): the Neumann integral
//! `N = ∮∮ dl₁·dl₂/|x₁ − x₂|`, which a coil's vector potential per unit strength
//! (`κ = μ₀I/4π = 1`) gives as `∮ A·dl` around the other coil.

use std::f64::consts::TAU;

use glam::DVec3;

use crate::magnetic::CircularLoop;

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
