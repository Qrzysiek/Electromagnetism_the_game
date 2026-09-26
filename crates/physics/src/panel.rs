//! Exact potential and field of a uniformly charged flat triangle (PHYSICS.md §2.7).
//!
//! For a planar polygon with unit normal `n`, an observation point `r` at height
//! `h = n·(r − r₀)` above the plane, and for each edge `i` (unit tangent `l_i`, outward
//! in-plane normal `m_i = l_i × n`), with `ρ` the projection of `r` onto the plane:
//! - `P⁰_i = (ρ_i^± − ρ)·m_i` (signed distance of the projection from the edge line),
//! - `l_i^± = (ρ_i^± − ρ)·l_i` (positions of the edge end points along the edge),
//! - `R_i^± = |r − r_i^±|`, `(R⁰_i)² = (P⁰_i)² + h²`,
//! - `f_i = ln((R_i^+ + l_i^+)/(R_i^− + l_i^−))`,
//! - `β_i = atan(P⁰ l^+ / ((R⁰)² + |h| R^+)) − atan(P⁰ l^− / ((R⁰)² + |h| R^−))`,
//!
//! the integrals are (Wilton et al., IEEE Trans. Antennas Propag. 32, 276 (1984);
//! Graglia, ibid. 41, 1448 (1993)):
//!
//! ```text
//! ∫ 1/R dS  = Σ_i P⁰_i f_i − |h| Σ_i β_i
//! ∫ ∇(1/R) dS = −Σ_i m_i f_i − n sign(h) Σ_i β_i      (gradient in r)
//! ```
//!
//! so a triangle with surface charge density σ has `φ = σ ∫1/R`, `E = −σ ∫∇(1/R)`
//! (units k = 1). The formulas are exact for every point off the triangle's edges; on the
//! triangle itself (h = 0 inside) the potential is continuous and the in-plane field is
//! finite, while the normal field has the jump ±2πσ.

use glam::DVec3;

/// A flat triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Triangle {
    pub a: DVec3,
    pub b: DVec3,
    pub c: DVec3,
}

impl Triangle {
    pub fn normal(&self) -> DVec3 {
        (self.b - self.a).cross(self.c - self.a).normalize()
    }

    pub fn area(&self) -> f64 {
        0.5 * (self.b - self.a).cross(self.c - self.a).length()
    }

    pub fn centroid(&self) -> DVec3 {
        (self.a + self.b + self.c) / 3.0
    }

    /// Longest edge.
    pub fn size(&self) -> f64 {
        (self.b - self.a)
            .length()
            .max((self.c - self.b).length())
            .max((self.a - self.c).length())
    }

    /// `(∫ 1/R dS, ∫ ∇(1/R) dS)` at `r`, exactly.
    pub fn integrals(&self, r: DVec3) -> (f64, DVec3) {
        Panel::new(*self).integrals(r)
    }
}

/// A triangle with its geometry precomputed for repeated evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Panel {
    pub triangle: Triangle,
    n: DVec3,
    /// Per edge: start, end, unit tangent, outward in-plane normal.
    edges: [(DVec3, DVec3, DVec3, DVec3); 3],
}

impl Panel {
    pub fn new(t: Triangle) -> Self {
        let n = t.normal();
        let edge = |p: DVec3, q: DVec3| {
            let l = (q - p).normalize();
            (p, q, l, l.cross(n))
        };
        Self {
            triangle: t,
            n,
            edges: [edge(t.a, t.b), edge(t.b, t.c), edge(t.c, t.a)],
        }
    }

    /// `(∫ 1/R dS, ∫ ∇(1/R) dS)` at `r`, exactly.
    pub fn integrals(&self, r: DVec3) -> (f64, DVec3) {
        let n = self.n;
        let h = n.dot(r - self.triangle.a);
        let abs_h = h.abs();
        let rho = r - n * h;
        let mut pot = 0.0;
        let mut grad = DVec3::ZERO;
        let mut beta_sum = 0.0;
        for &(p, q, l, m) in &self.edges {
            let p0 = (p - rho).dot(m);
            let l_minus = (p - rho).dot(l);
            let l_plus = (q - rho).dot(l);
            let r_minus = (r - p).length();
            let r_plus = (r - q).length();
            let r0_sq = p0 * p0 + h * h;
            // f = ln((R+ + l+)/(R- + l-)); for points on the edge line beyond an end
            // the arguments vanish together: use the equivalent form with R − l.
            let f = if l_minus > 0.0 || r_minus + l_minus > 1e-12 * r_minus.max(1e-300) {
                ((r_plus + l_plus) / (r_minus + l_minus)).ln()
            } else {
                ((r_minus - l_minus) / (r_plus - l_plus)).ln()
            };
            let beta = if r0_sq > 0.0 {
                (p0 * l_plus / (r0_sq + abs_h * r_plus)).atan()
                    - (p0 * l_minus / (r0_sq + abs_h * r_minus)).atan()
            } else {
                0.0
            };
            pot += p0 * f;
            beta_sum += beta;
            grad -= m * f;
        }
        pot -= abs_h * beta_sum;
        grad -= n * (h.signum() * beta_sum);
        if h == 0.0 {
            // In the plane: the normal component is the average of the two sides (0).
            grad -= n * grad.dot(n);
        }
        (pot, grad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tri() -> Triangle {
        Triangle {
            a: DVec3::new(0.2, -0.1, 0.3),
            b: DVec3::new(1.4, 0.3, 0.1),
            c: DVec3::new(0.5, 1.2, -0.2),
        }
    }

    /// Reference by brute-force quadrature: subdivide into 4^k small triangles and use the
    /// 7-point degree-5 rule on each.
    fn quadrature(t: &Triangle, r: DVec3, level: u32) -> (f64, DVec3) {
        let mut tris = vec![*t];
        for _ in 0..level {
            let mut next = Vec::with_capacity(tris.len() * 4);
            for s in tris {
                let (ab, bc, ca) = ((s.a + s.b) / 2.0, (s.b + s.c) / 2.0, (s.c + s.a) / 2.0);
                next.push(Triangle {
                    a: s.a,
                    b: ab,
                    c: ca,
                });
                next.push(Triangle {
                    a: ab,
                    b: s.b,
                    c: bc,
                });
                next.push(Triangle {
                    a: ca,
                    b: bc,
                    c: s.c,
                });
                next.push(Triangle {
                    a: ab,
                    b: bc,
                    c: ca,
                });
            }
            tris = next;
        }
        // Dunavant degree-5 rule (barycentric weights).
        let a1 = 0.059_715_871_789_770;
        let b1 = 0.470_142_064_105_115;
        let a2 = 0.797_426_985_353_087;
        let b2 = 0.101_286_507_323_456;
        let pts = [
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225),
            (a1, b1, b1, 0.132_394_152_788_506),
            (b1, a1, b1, 0.132_394_152_788_506),
            (b1, b1, a1, 0.132_394_152_788_506),
            (a2, b2, b2, 0.125_939_180_544_827),
            (b2, a2, b2, 0.125_939_180_544_827),
            (b2, b2, a2, 0.125_939_180_544_827),
        ];
        let mut pot = 0.0;
        let mut grad = DVec3::ZERO;
        for s in &tris {
            let area = s.area();
            for &(u, v, w, wt) in &pts {
                let y = s.a * u + s.b * v + s.c * w;
                let d = r - y;
                let inv = 1.0 / d.length();
                pot += wt * area * inv;
                grad -= d * (wt * area * inv * inv * inv);
            }
        }
        (pot, grad)
    }

    #[test]
    fn matches_quadrature_off_the_triangle() {
        let t = tri();
        for r in [
            DVec3::new(0.7, 0.4, 0.8),
            DVec3::new(-1.0, 2.0, 0.5),
            DVec3::new(3.0, -2.0, -1.5),
            DVec3::new(0.6, 0.4, -0.05),
        ] {
            let (p, g) = t.integrals(r);
            let (pq, gq) = quadrature(&t, r, 6);
            assert!((p - pq).abs() < 1e-9 * pq.abs(), "{r}: {p} vs {pq}");
            assert!((g - gq).length() < 1e-8 * gq.length(), "{r}: {g} vs {gq}");
        }
    }

    #[test]
    fn gradient_is_the_derivative_of_the_potential() {
        let t = tri();
        let h = 1e-6;
        for r in [
            DVec3::new(0.7, 0.4, 0.3),
            DVec3::new(2.0, 0.0, 0.2),
            DVec3::new(0.5, 0.5, 0.2),
        ] {
            let g = t.integrals(r).1;
            for axis in [DVec3::X, DVec3::Y, DVec3::Z] {
                let d = (t.integrals(r + axis * h).0 - t.integrals(r - axis * h).0) / (2.0 * h);
                assert!(
                    (g.dot(axis) - d).abs() < 1e-7 * g.length().max(1.0),
                    "{r} {axis}"
                );
            }
        }
    }

    /// Crossing the triangle, the normal field jumps by 4π (for unit σ, ∫∇(1/R)
    /// changes by 4π n); the potential is continuous.
    #[test]
    fn jump_across_the_surface() {
        let t = tri();
        let n = t.normal();
        let x = t.centroid();
        let e = 1e-9;
        let (p1, g1) = t.integrals(x + n * e);
        let (p2, g2) = t.integrals(x - n * e);
        assert!((p1 - p2).abs() < 1e-7);
        let jump = (g2 - g1).dot(n);
        assert!((jump - 4.0 * std::f64::consts::PI).abs() < 1e-6, "{jump}");
    }
}
