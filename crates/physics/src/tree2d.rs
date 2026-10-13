//! A tree code for the z-invariant world's sums (PHYSICS.md §2.11): the field of many line
//! charges and charged segments at many points, and the line charges' potential, in
//! O(N log N) instead of O(N²).
//!
//! The sources are line charges (λ at z_k) and uniformly charged straight segments (density
//! σ from a to b). In the complex plane a group's potential is `Re Ω`,
//! `Ω(z) = −2 Σ q ln(z − z_k)` (a segment integrated along it), and outside a disc of radius
//! r about z_c it is the convergent multipole series
//! `Ω(z) = −2 [a₀ ln(z − z_c) − Σ_{n≥1} aₙ/(n (z − z_c)ⁿ)]`, `aₙ = Σ q (z_k − z_c)ⁿ` (a
//! segment's moments exact: `σ e^{inθ} 2 (L/2)^{n+1}/(n+1)` for even n about its middle,
//! zero for odd), with `E_x − i E_y = −Ω′(z) = 2 Σ aₙ/(z − z_c)^{n+1}`. A quadtree holds the
//! sources; a node whose radius is below `THETA` times its distance to the target is summed
//! by its series (`ORDER` terms: the omitted ones below `THETA^(ORDER+1)` ≈ 5e-7 of it),
//! every other source directly with its exact kernel: a line charge as a uniformly charged
//! rod of radius h (`2λ r⃗/max(r², h²)`), a segment by its exact integrals. A node is taken
//! as a series only when the target lies farther than its radius plus h, so the rods'
//! softening never enters a series: the tree reproduces the direct sums to that truncation.
//!
//! Deterministic: the tree's construction and each target's traversal follow a fixed order,
//! and the logarithms are `libm`'s.

use glam::DVec3;

use crate::zinv::{line_charge, segment_field};

/// Terms of the multipole series.
pub const ORDER: usize = 20;
/// A node is summed by its series when its radius is below this fraction of its distance.
pub const THETA: f64 = 0.5;
/// Most sources in a leaf.
const LEAF: usize = 32;
/// Deepest subdivision (coincident sources stop it).
const DEPTH: u32 = 40;

/// A complex number (re, im).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    const ZERO: C = C { re: 0.0, im: 0.0 };
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: C) -> C {
        C::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn scale(self, s: f64) -> C {
        C::new(self.re * s, self.im * s)
    }
    fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    fn inv(self) -> C {
        let d = self.abs2();
        C::new(self.re / d, -self.im / d)
    }
}

/// A source: a line charge at `a` (`b` = `a`) or a uniformly charged segment from `a` to `b`.
#[derive(Clone, Copy, Debug)]
pub struct Source {
    pub a: DVec3,
    pub b: DVec3,
}

impl Source {
    pub fn point(a: DVec3) -> Self {
        Self { a, b: a }
    }
    pub fn segment(a: DVec3, b: DVec3) -> Self {
        Self { a, b }
    }
    fn is_point(&self) -> bool {
        self.a == self.b
    }
    fn middle(&self) -> C {
        C::new(0.5 * (self.a.x + self.b.x), 0.5 * (self.a.y + self.b.y))
    }
    fn half_length(&self) -> f64 {
        0.5 * libm::sqrt((self.b.x - self.a.x).powi(2) + (self.b.y - self.a.y).powi(2))
    }
}

#[derive(Clone, Debug)]
struct Node {
    center: C,
    /// Every source within this distance of the centre (a segment wholly).
    radius: f64,
    /// Its sources: `order[start..end]`.
    start: usize,
    end: usize,
    /// Its children's indices (none: a leaf).
    children: Vec<usize>,
}

/// The tree of a set of sources (their geometry); `Tree::multipoles` sums their charges.
#[derive(Clone, Debug)]
pub struct Tree {
    sources: Vec<Source>,
    nodes: Vec<Node>,
    /// The sources in the tree's order.
    order: Vec<usize>,
    /// The line charges' radius (rods).
    rod: f64,
    binomial: Vec<Vec<f64>>,
}

/// A tree's multipole coefficients for one set of charges (per node, `ORDER + 1` each).
#[derive(Clone, Debug)]
pub struct Multipoles {
    coeffs: Vec<[C; ORDER + 1]>,
    charges: Vec<f64>,
}

impl Tree {
    /// The tree of `sources`; line charges are rods of radius `rod`.
    pub fn new(sources: Vec<Source>, rod: f64) -> Self {
        let n = sources.len();
        let mut order: Vec<usize> = (0..n).collect();
        let mut nodes = Vec::new();
        if n > 0 {
            let (mut lo, mut hi) = (C::new(f64::MAX, f64::MAX), C::new(f64::MIN, f64::MIN));
            for s in &sources {
                let m = s.middle();
                lo = C::new(lo.re.min(m.re), lo.im.min(m.im));
                hi = C::new(hi.re.max(m.re), hi.im.max(m.im));
            }
            let half = 0.5 * (hi.re - lo.re).max(hi.im - lo.im);
            let center = C::new(0.5 * (lo.re + hi.re), 0.5 * (lo.im + hi.im));
            build(&sources, &mut order, &mut nodes, 0, n, center, half, 0);
        }
        let binomial = (0..=ORDER)
            .map(|n| (0..=n).map(|k| binomial_coeff(n, k)).collect())
            .collect();
        Self {
            sources,
            nodes,
            order,
            rod,
            binomial,
        }
    }

    /// The coefficients of every node for `charges` (per source: a line charge's λ, a
    /// segment's density σ).
    ///
    /// # Panics
    /// If `charges` has another length than the sources.
    pub fn multipoles(&self, charges: &[f64]) -> Multipoles {
        assert_eq!(charges.len(), self.sources.len());
        let mut coeffs = vec![[C::ZERO; ORDER + 1]; self.nodes.len()];
        // Children are stored after their parents: sum from the last node back.
        for i in (0..self.nodes.len()).rev() {
            let node = &self.nodes[i];
            let mut a = [C::ZERO; ORDER + 1];
            if node.children.is_empty() {
                for &k in &self.order[node.start..node.end] {
                    let s = &self.sources[k];
                    let d = s.middle().sub(node.center);
                    if s.is_point() {
                        // q (z_k − z_c)^n.
                        let mut p = C::new(charges[k], 0.0);
                        for an in &mut a {
                            *an = an.add(p);
                            p = p.mul(d);
                        }
                    } else {
                        // The segment's moments about its middle, shifted to the centre.
                        let own = segment_moments(s, charges[k]);
                        self.shift_into(&own, d, &mut a);
                    }
                }
            } else {
                for &c in &node.children {
                    let d = self.nodes[c].center.sub(node.center);
                    let own = coeffs[c];
                    self.shift_into(&own, d, &mut a);
                }
            }
            coeffs[i] = a;
        }
        Multipoles {
            coeffs,
            charges: charges.to_vec(),
        }
    }

    /// Adds to `a` the coefficients `own` (about a point `d` from `a`'s centre):
    /// `aₙ += Σ_k C(n, k) ownₖ d^(n−k)`.
    fn shift_into(&self, own: &[C; ORDER + 1], d: C, a: &mut [C; ORDER + 1]) {
        let mut dp = [C::new(1.0, 0.0); ORDER + 1];
        for n in 1..=ORDER {
            dp[n] = dp[n - 1].mul(d);
        }
        for n in 0..=ORDER {
            let mut sum = C::ZERO;
            for k in 0..=n {
                sum = sum.add(own[k].mul(dp[n - k]).scale(self.binomial[n][k]));
            }
            a[n] = a[n].add(sum);
        }
    }

    /// The field at `x` of every source but `skip` (an index into the sources).
    pub fn field(&self, m: &Multipoles, x: DVec3, skip: Option<usize>) -> DVec3 {
        let mut e = DVec3::ZERO;
        if self.nodes.is_empty() {
            return e;
        }
        let z = C::new(x.x, x.y);
        let h2 = self.rod * self.rod;
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let node = &self.nodes[i];
            let d = z.sub(node.center);
            let dist2 = d.abs2();
            let far = node.radius < THETA * libm::sqrt(dist2)
                && libm::sqrt(dist2) > node.radius + self.rod;
            if far {
                // E_x − i E_y = 2 Σ aₙ w^(n+1), w = 1/(z − z_c).
                let w = d.inv();
                let mut wp = w;
                let mut s = C::ZERO;
                for an in &m.coeffs[i] {
                    s = s.add(an.mul(wp));
                    wp = wp.mul(w);
                }
                e += DVec3::new(2.0 * s.re, -2.0 * s.im, 0.0);
            } else if node.children.is_empty() {
                for &k in &self.order[node.start..node.end] {
                    if Some(k) == skip {
                        continue;
                    }
                    let s = &self.sources[k];
                    let q = m.charges[k];
                    if s.is_point() {
                        let r = DVec3::new(x.x - s.a.x, x.y - s.a.y, 0.0);
                        let r2 = (r.x * r.x + r.y * r.y).max(h2);
                        e += r * (2.0 * q / r2);
                    } else {
                        e += segment_field(s.a, s.b, x) * q;
                    }
                }
            } else {
                // Children in reverse, so that they are visited in order.
                stack.extend(node.children.iter().rev());
            }
        }
        e
    }

    /// The potential at `x` of the line charges as points (`−λ ln r²`, unsoftened: the
    /// induced charge's right-hand side) and of the segments.
    pub fn potential(&self, m: &Multipoles, x: DVec3) -> f64 {
        let mut phi = 0.0;
        if self.nodes.is_empty() {
            return phi;
        }
        let z = C::new(x.x, x.y);
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let node = &self.nodes[i];
            let d = z.sub(node.center);
            let dist = libm::sqrt(d.abs2());
            if node.radius < THETA * dist {
                // Re Ω = −2 [a₀ ln|z − z_c| − Re Σ_{n≥1} aₙ wⁿ / n].
                let w = d.inv();
                let a = &m.coeffs[i];
                let mut wp = w;
                let mut s = C::ZERO;
                for (n, an) in a.iter().enumerate().skip(1) {
                    #[allow(clippy::cast_precision_loss)]
                    let inv_n = 1.0 / n as f64;
                    s = s.add(an.mul(wp).scale(inv_n));
                    wp = wp.mul(w);
                }
                phi += -2.0 * (a[0].re * libm::log(dist) - a[0].im * libm::atan2(d.im, d.re))
                    + 2.0 * s.re;
            } else if node.children.is_empty() {
                for &k in &self.order[node.start..node.end] {
                    let s = &self.sources[k];
                    let q = m.charges[k];
                    if s.is_point() {
                        phi += line_charge(s.a, q, x).0;
                    } else {
                        phi += crate::zinv::segment_integrals(s.a, s.b, x).0 * q;
                    }
                }
            } else {
                stack.extend(node.children.iter().rev());
            }
        }
        phi
    }
}

/// A segment's moments about its middle: `σ e^{inθ} 2 (L/2)^{n+1}/(n+1)` for even n.
fn segment_moments(s: &Source, sigma: f64) -> [C; ORDER + 1] {
    let half = s.half_length();
    let (dx, dy) = (s.b.x - s.a.x, s.b.y - s.a.y);
    let e = C::new(dx / (2.0 * half), dy / (2.0 * half));
    let mut out = [C::ZERO; ORDER + 1];
    let mut en = C::new(1.0, 0.0);
    let mut hp = half;
    for (n, o) in out.iter_mut().enumerate() {
        if n % 2 == 0 {
            #[allow(clippy::cast_precision_loss)]
            let m = 2.0 * hp / (n as f64 + 1.0);
            *o = en.scale(sigma * m);
        }
        en = en.mul(e);
        hp *= half;
    }
    out
}

fn binomial_coeff(n: usize, k: usize) -> f64 {
    let mut c = 1.0;
    for i in 0..k {
        #[allow(clippy::cast_precision_loss)]
        {
            c = c * (n - i) as f64 / (i + 1) as f64;
        }
    }
    c
}

/// Subdivides `order[start..end]` (sources by their middles) into the node of `center` and
/// half-width `half`, appending it and its descendants to `nodes`; returns its index.
#[allow(clippy::too_many_arguments)]
fn build(
    sources: &[Source],
    order: &mut [usize],
    nodes: &mut Vec<Node>,
    start: usize,
    end: usize,
    center: C,
    half: f64,
    depth: u32,
) -> usize {
    let radius = order[start..end]
        .iter()
        .map(|&k| {
            let s = &sources[k];
            libm::sqrt(s.middle().sub(center).abs2()) + s.half_length()
        })
        .fold(0.0, f64::max);
    let index = nodes.len();
    nodes.push(Node {
        center,
        radius,
        start,
        end,
        children: Vec::new(),
    });
    if end - start <= LEAF || depth >= DEPTH || half <= 0.0 {
        return index;
    }
    // Quadrants: x below/above the centre, then y; a stable partition keeps the order fixed.
    let quadrant = |k: usize| {
        let m = sources[k].middle();
        usize::from(m.re >= center.re) + 2 * usize::from(m.im >= center.im)
    };
    order[start..end].sort_by_key(|&k| quadrant(k));
    let mut children = Vec::new();
    let mut lo = start;
    for q in 0..4 {
        let hi = lo
            + order[lo..end]
                .iter()
                .take_while(|&&k| quadrant(k) == q)
                .count();
        if hi > lo {
            let h = 0.5 * half;
            let c = C::new(
                center.re + if q % 2 == 1 { h } else { -h },
                center.im + if q / 2 == 1 { h } else { -h },
            );
            children.push(build(sources, order, nodes, lo, hi, c, h, depth + 1));
        }
        lo = hi;
    }
    nodes[index].children = children;
    index
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zinv::segment_field;

    /// A deterministic scatter in a 24 × 16 box.
    fn scatter(n: usize, seed: u64) -> Vec<DVec3> {
        let mut s = seed;
        let mut next = || {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            #[allow(clippy::cast_precision_loss)]
            let v = (s >> 11) as f64 / (1u64 << 53) as f64;
            v
        };
        (0..n)
            .map(|_| DVec3::new(24.0 * next(), 16.0 * next(), 0.0))
            .collect()
    }

    /// The tree against the direct sums: 2000 rods (radius 0.0625) and 300 segments of
    /// length 0.5 with charges of both signs; the field at every rod (its own excluded) and
    /// the potential at 200 points. Required: 1e-7 of the field's and the potential's
    /// scale, four orders below the mesh's error (1e-3) and below the direct path's segment
    /// series (`segment_field_fast`: 4.9e-5).
    #[test]
    fn tree_matches_direct_sums() {
        let rod = 0.0625;
        let pts = scatter(2000, 1);
        let ends = scatter(300, 2);
        let mut sources: Vec<Source> = pts.iter().map(|p| Source::point(*p)).collect();
        sources.extend(
            ends.iter()
                .map(|a| Source::segment(*a, *a + DVec3::new(0.3, 0.4, 0.0))),
        );
        #[allow(clippy::cast_precision_loss)]
        let charges: Vec<f64> = (0..sources.len())
            .map(|k| if k % 3 == 0 { -1.0 } else { 0.7 } * (1.0 + (k % 7) as f64 * 0.1))
            .collect();
        let tree = Tree::new(sources.clone(), rod);
        let m = tree.multipoles(&charges);
        let direct_field = |x: DVec3, skip: usize| {
            let mut e = DVec3::ZERO;
            for (k, (s, q)) in sources.iter().zip(&charges).enumerate() {
                if k == skip {
                    continue;
                }
                if s.is_point() {
                    let r = DVec3::new(x.x - s.a.x, x.y - s.a.y, 0.0);
                    e += r * (2.0 * q / (r.x * r.x + r.y * r.y).max(rod * rod));
                } else {
                    e += segment_field(s.a, s.b, x) * *q;
                }
            }
            e
        };
        let (mut worst, mut scale) = (0.0_f64, 0.0_f64);
        for (k, p) in pts.iter().enumerate() {
            let exact = direct_field(*p, k);
            let tree_e = tree.field(&m, *p, Some(k));
            worst = worst.max((tree_e - exact).length());
            scale = scale.max(exact.length());
        }
        let probes = scatter(200, 3);
        let (mut worst_phi, mut scale_phi) = (0.0_f64, 0.0_f64);
        for x in &probes {
            let mut exact = 0.0;
            for (s, q) in sources.iter().zip(&charges) {
                exact += if s.is_point() {
                    line_charge(s.a, *q, *x).0
                } else {
                    crate::zinv::segment_integrals(s.a, s.b, *x).0 * q
                };
            }
            worst_phi = worst_phi.max((tree.potential(&m, *x) - exact).abs());
            scale_phi = scale_phi.max(exact.abs());
        }
        println!(
            "tree: field {:.2e} of {scale:.3e}, potential {:.2e} of {scale_phi:.3e}",
            worst / scale,
            worst_phi / scale_phi
        );
        assert!(worst < 1e-7 * scale, "field {worst:e} of {scale:e}");
        assert!(
            worst_phi < 1e-7 * scale_phi,
            "potential {worst_phi:e} of {scale_phi:e}"
        );
    }
}
