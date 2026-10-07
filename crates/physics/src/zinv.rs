//! The z-invariant world (docs/TUBES.md, decision 1; PHYSICS.md §2.11): every source
//! is translation-invariant along z, so the plane z = 0 is a cross-section of the whole
//! and the fields are those of line charges. This is exact 3D electrostatics restricted
//! by a symmetry, the one the textbook diode laws assume (Child–Langmuir, Langmuir–
//! Blodgett), not a 2D model of the 3D slice the rest of the game uses.
//!
//! Units as everywhere (k = 1): a line charge λ (charge per unit length along z) has
//! `φ = −2λ ln r` and `E = 2λ r̂/r`. Electrodes are prisms along z; their cross-sections
//! are polygons cut into segments carrying a constant surface charge density σ (charge
//! per unit area), collocated at the segments' midpoints (the potential there equals the
//! electrode's), with exact segment integrals.
//!
//! **Reference.** A net line charge's potential grows like `−2Q ln r` without bound, so
//! potentials are defined relative to the circle of unit radius about each source (the
//! logarithm's zero). Results that do not depend on that choice (and the validation
//! uses only such) come from enclosed setups: a grounded or driven outer electrode, or
//! zero net charge.

use glam::DVec3;

use crate::bem::{Bias, Lu};

/// Potential and field at `x` of the segment from `a` to `b` (in the plane) carrying a
/// unit surface charge density, exactly: with `u` along the segment from `a`, `v` across
/// it (to the left), `φ = −[G(u) − G(u − L)]`, `G(w) = w ln(w² + v²) − 2w + 2v atan(w/v)`,
/// `E_u = ln(r_a² / r_b²)`, `E_v = 2 [atan((L − u)/v) + atan(u/v)]` (twice the angle the
/// segment subtends; on the segment itself its average, 0).
pub fn segment_integrals(a: DVec3, b: DVec3, x: DVec3) -> (f64, DVec3) {
    let d = b - a;
    let l = libm::sqrt(d.x * d.x + d.y * d.y);
    let t = DVec3::new(d.x / l, d.y / l, 0.0);
    let n = DVec3::new(-t.y, t.x, 0.0);
    let r = x - a;
    let (u, mut v) = (r.x * t.x + r.y * t.y, r.x * n.x + r.y * n.y);
    // On the segment's line within rounding (a collocation point computed as the
    // midpoint): exactly on it, where the normal field's average is 0.
    if v.abs() <= 1e-12 * l {
        v = 0.0;
    }
    let g = |w: f64| {
        let s = w * w + v * v;
        let log = if s > 0.0 { w * libm::log(s) } else { 0.0 };
        let arc = if v == 0.0 {
            0.0
        } else {
            2.0 * v * libm::atan(w / v)
        };
        log - 2.0 * w + arc
    };
    let phi = -(g(u) - g(u - l));
    let ra2 = u * u + v * v;
    let rb2 = (u - l) * (u - l) + v * v;
    let e_u = libm::log(ra2 / rb2);
    let e_v = if v == 0.0 {
        0.0
    } else {
        2.0 * (libm::atan((l - u) / v) + libm::atan(u / v))
    };
    (phi, t * e_u + n * e_v)
}

/// A line charge's potential and field at `x`.
pub fn line_charge(p: DVec3, lambda: f64, x: DVec3) -> (f64, DVec3) {
    let r = DVec3::new(x.x - p.x, x.y - p.y, 0.0);
    let r2 = r.x * r.x + r.y * r.y;
    (-lambda * libm::log(r2), r * (2.0 * lambda / r2))
}

/// The cross-section of an electrode prism.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Section {
    /// A rectangle: centre, in-plane angle, half-extents along and across it.
    Rect {
        center: DVec3,
        angle: f64,
        half_length: f64,
        half_thickness: f64,
    },
    /// A circle (a wire or a tube's wall seen end-on), as a regular polygon.
    Circle { center: DVec3, radius: f64 },
}

impl Section {
    /// The boundary as segments `(start, end)`, counter-clockwise, none longer than
    /// `size`; a rectangle's sides graded towards the corners (Chebyshev), where the
    /// density is singular.
    pub fn mesh(&self, size: f64) -> Vec<(DVec3, DVec3)> {
        match *self {
            Section::Rect {
                center,
                angle,
                half_length,
                half_thickness,
            } => {
                let (s, c) = libm::sincos(angle);
                let (u, v) = (DVec3::new(c, s, 0.0), DVec3::new(-s, c, 0.0));
                let corners = [
                    center - u * half_length - v * half_thickness,
                    center + u * half_length - v * half_thickness,
                    center + u * half_length + v * half_thickness,
                    center - u * half_length + v * half_thickness,
                ];
                let mut out = Vec::new();
                for k in 0..4 {
                    let (p, q) = (corners[k], corners[(k + 1) % 4]);
                    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                    let n = (((q - p).length() / size).ceil() as usize).max(2);
                    #[allow(clippy::cast_precision_loss)]
                    let at = |i: usize| {
                        let f = 0.5 * (1.0 - libm::cos(std::f64::consts::PI * i as f64 / n as f64));
                        p + (q - p) * f
                    };
                    for i in 0..n {
                        out.push((at(i), at(i + 1)));
                    }
                }
                out
            }
            Section::Circle { center, radius } => {
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let n = ((std::f64::consts::TAU * radius / size).ceil() as usize).max(8);
                #[allow(clippy::cast_precision_loss)]
                let at = |i: usize| {
                    let a = std::f64::consts::TAU * i as f64 / n as f64;
                    center + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * radius
                };
                (0..n).map(|i| (at(i), at(i + 1))).collect()
            }
        }
    }
}

/// An electrode prism along z: its cross-section and how it is held.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Electrode {
    pub section: Section,
    pub bias: Bias,
}

/// The electrodes of a z-invariant level with their surface charge for given line
/// charges: the source system (all electrodes at 0) plus unit systems combined through
/// the capacitance matrix per unit length, as for the 3D electrodes (`bem.rs`).
#[derive(Debug)]
pub struct Electrodes {
    pub electrodes: Vec<Electrode>,
    /// Segments `(start, end)` and the electrode each belongs to.
    pub segments: Vec<(DVec3, DVec3)>,
    owner: Vec<usize>,
    lu: Lu,
    unit: Vec<Vec<f64>>,
    /// Capacitance matrix per unit length: `capacitance[i][j]` the charge per unit
    /// length on `i` with `j` at unit potential and the others at 0.
    pub capacitance: Vec<Vec<f64>>,
    /// Surface charge density of each segment.
    pub sigma: Vec<f64>,
    /// Each electrode's potential.
    pub potentials: Vec<f64>,
}

impl Electrodes {
    /// Solves the electrodes for the line charges `sources` `(position, λ)`, with segments
    /// no longer than `size`.
    pub fn new(electrodes: Vec<Electrode>, sources: &[(DVec3, f64)], size: f64) -> Self {
        let mut segments = Vec::new();
        let mut owner = Vec::new();
        for (i, e) in electrodes.iter().enumerate() {
            let m = e.section.mesh(size);
            owner.extend(std::iter::repeat_n(i, m.len()));
            segments.extend(m);
        }
        let n = segments.len();
        let mid: Vec<DVec3> = segments.iter().map(|(a, b)| (*a + *b) * 0.5).collect();
        let mut a = vec![0.0; n * n];
        for (i, x) in mid.iter().enumerate() {
            crate::cancel::checkpoint();
            for (j, (p, q)) in segments.iter().enumerate() {
                a[i * n + j] = segment_integrals(*p, *q, *x).0;
            }
        }
        let lu = Lu::new(a, n);
        let m = electrodes.len();
        let unit: Vec<Vec<f64>> = (0..m)
            .map(|e| {
                let b: Vec<f64> = owner
                    .iter()
                    .map(|&o| if o == e { 1.0 } else { 0.0 })
                    .collect();
                lu.solve(&b)
            })
            .collect();
        let charge_on = |sigma: &[f64], e: usize| -> f64 {
            segments
                .iter()
                .zip(sigma)
                .zip(&owner)
                .filter(|(_, o)| **o == e)
                .map(|(((p, q), s), _)| s * (*q - *p).length())
                .sum()
        };
        let capacitance: Vec<Vec<f64>> = (0..m)
            .map(|i| (0..m).map(|j| charge_on(&unit[j], i)).collect())
            .collect();
        // The source system: the electrodes at 0 against the line charges.
        let b: Vec<f64> = mid
            .iter()
            .map(|x| {
                -sources
                    .iter()
                    .map(|&(p, l)| line_charge(p, l, *x).0)
                    .sum::<f64>()
            })
            .collect();
        let mut sigma = lu.solve(&b);
        let g: Vec<f64> = (0..m).map(|e| charge_on(&sigma, e)).collect();
        let mut alpha: Vec<f64> = electrodes
            .iter()
            .map(|e| match e.bias {
                Bias::Potential(v) => v,
                _ => 0.0,
            })
            .collect();
        let floating: Vec<usize> = (0..m)
            .filter(|&i| matches!(electrodes[i].bias, Bias::Charge(_)))
            .collect();
        if !floating.is_empty() {
            let k = floating.len();
            let mut mat = vec![0.0; k * k];
            let mut rhs = vec![0.0; k];
            for (r, &i) in floating.iter().enumerate() {
                for (c, &j) in floating.iter().enumerate() {
                    mat[r * k + c] = capacitance[i][j];
                }
                let Bias::Charge(q) = electrodes[i].bias else {
                    unreachable!()
                };
                let known: f64 = (0..m)
                    .filter(|j| !floating.contains(j))
                    .map(|j| capacitance[i][j] * alpha[j])
                    .sum();
                rhs[r] = q - g[i] - known;
            }
            let x = Lu::new(mat, k).solve(&rhs);
            for (r, &i) in floating.iter().enumerate() {
                alpha[i] = x[r];
            }
        }
        for (j, &v) in alpha.iter().enumerate() {
            if v != 0.0 {
                for (s, u) in sigma.iter_mut().zip(&unit[j]) {
                    *s += v * u;
                }
            }
        }
        Self {
            electrodes,
            segments,
            owner,
            lu,
            unit,
            capacitance,
            sigma,
            potentials: alpha,
        }
    }

    /// The electrodes' potential and field at `x` (exact segment integrals).
    pub fn sample(&self, x: DVec3) -> (f64, DVec3) {
        let mut phi = 0.0;
        let mut e = DVec3::ZERO;
        for ((p, q), s) in self.segments.iter().zip(&self.sigma) {
            let (f, g) = segment_integrals(*p, *q, x);
            phi += s * f;
            e += g * *s;
        }
        (phi, e)
    }

    /// Net charge per unit length on electrode `i`.
    pub fn charge(&self, i: usize) -> f64 {
        self.charges_of(&self.sigma)[i]
    }

    /// Each electrode's charge per unit length for the surface density `sigma` (e.g. the
    /// induced one; its rate of change is the electrode's Ramo current).
    pub fn charges_of(&self, sigma: &[f64]) -> Vec<f64> {
        let mut q = vec![0.0; self.electrodes.len()];
        for (((p, r), s), &o) in self.segments.iter().zip(sigma).zip(&self.owner) {
            q[o] += s * (*r - *p).length();
        }
        q
    }

    /// The surface charge that cancels the potential `phi_at(midpoint)` of other sources
    /// on the electrodes (all held at 0): the induced charge, e.g. of the particles.
    pub fn induced_by(&self, phi_at: impl Fn(DVec3) -> f64) -> Vec<f64> {
        let b: Vec<f64> = self
            .segments
            .iter()
            .map(|(a, c)| -phi_at((*a + *c) * 0.5))
            .collect();
        self.lu.solve(&b)
    }

    /// The electrode segment `j` belongs to.
    pub fn owner_of(&self, j: usize) -> usize {
        self.owner[j]
    }

    /// The unit system of electrode `j` (its density with `j` at 1, the others at 0).
    pub fn unit(&self, j: usize) -> &[f64] {
        &self.unit[j]
    }
}
