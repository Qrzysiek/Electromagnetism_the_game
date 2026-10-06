//! Box-shaped electrodes by the boundary element method (PHYSICS.md §2.7).
//!
//! An electrode is a rectangular box standing on the plane of the 2D slice: symmetric
//! about z = 0, of in-plane length and thickness and height `±hz` in z (a plate, a slab,
//! a wall). Its surface is cut into flat triangles carrying constant surface charge
//! density; the densities follow from collocation at the triangle centroids (the
//! surface potential equals the electrode's potential), with exact triangle integrals
//! (`panel.rs`).
//!
//! **Symmetry.** Only the upper half (z ≥ 0) is meshed; each panel stands for itself and
//! its mirror image. In the plane z = 0 the mirror's field is the panel's field with
//! E_z reversed, applied exactly, so E_z = 0 there bit for bit.
//!
//! **Mesh.** Each face is divided into a grid whose lines are graded towards the edges
//! (Chebyshev spacing), where the surface charge density is singular.
//!
//! **Bias.** As for spheres (`conductor.rs`): the source system (fixed charges, all
//! electrodes at 0) plus unit systems (one electrode at 1, the others at 0) combined
//! through the capacitance matrix; grounded, fixed potential or floating with a net
//! charge.
//!
//! **Accuracy.** Constant-density panels are a discretization: the error is measured by
//! comparing two mesh resolutions (preview and verification), which is how it enters
//! the verification of flights. The field of a given panel set is the exact field of a
//! real charge distribution, so energy is conserved to integration accuracy.
//!
//! **Evaluation.** Every panel uses the exact integrals, so the field is exactly the
//! gradient of the potential of the panel charges: energy is conserved to integration
//! accuracy. (A first version switched to a quadrature rule for distant panels; the
//! switch broke exact conservation, 1.5e-9 energy drift, and was removed.)

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use glam::DVec3;

use crate::field::{FieldSample, FieldSolver};
use crate::panel::{Panel, Triangle};

/// How an electrode is held (same meaning as for spheres).
pub use crate::conductor::Bias;

/// A box electrode: centre in the plane z = 0, in-plane direction angle, half-extents
/// along that direction (`half_length`), across it (`half_thickness`) and in z
/// (`half_height`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoxElectrode {
    pub center: DVec3,
    pub angle: f64,
    pub half_length: f64,
    pub half_thickness: f64,
    pub half_height: f64,
    pub bias: Bias,
}

/// Mesh resolution: largest panel size relative to the electrode's smallest dimension is
/// set by `cells_per_edge` along each face edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Resolution {
    Preview,
    Verify,
    Display,
}

impl Resolution {
    /// Target panel size in cells (upper bound; faces get at least 2 divisions).
    fn panel_size(self) -> f64 {
        match self {
            Resolution::Preview => 0.5,
            Resolution::Verify => 0.35,
            Resolution::Display => 1.0,
        }
    }
}

impl BoxElectrode {
    fn axes(&self) -> (DVec3, DVec3) {
        let (s, c) = libm::sincos(self.angle);
        (DVec3::new(c, s, 0.0), DVec3::new(-s, c, 0.0))
    }

    /// Upper-half surface (z ≥ 0) as graded triangles.
    pub fn mesh(&self, size: f64) -> Vec<Triangle> {
        let (u, v) = self.axes();
        let w = DVec3::Z;
        let (a, b, h) = (self.half_length, self.half_thickness, self.half_height);
        let c = self.center;
        let mut out = Vec::new();
        // Face given by origin corner and two edge vectors; graded both ways.
        let mut face = |o: DVec3, e1: DVec3, e2: DVec3, grade2_both: bool| {
            let n1 = divisions(e1.length(), size);
            let n2 = divisions(e2.length(), size);
            let g1 = graded(n1, true);
            let g2 = graded(n2, grade2_both);
            for i in 0..n1 {
                for j in 0..n2 {
                    let p = |s: f64, t: f64| o + e1 * s + e2 * t;
                    let (p00, p10, p01, p11) = (
                        p(g1[i], g2[j]),
                        p(g1[i + 1], g2[j]),
                        p(g1[i], g2[j + 1]),
                        p(g1[i + 1], g2[j + 1]),
                    );
                    out.push(Triangle {
                        a: p00,
                        b: p10,
                        c: p11,
                    });
                    out.push(Triangle {
                        a: p00,
                        b: p11,
                        c: p01,
                    });
                }
            }
        };
        // Top face z = h (normal +z): corners (−a,−b) → (+a,+b).
        face(
            c - u * a - v * b + w * h,
            u * (2.0 * a),
            v * (2.0 * b),
            true,
        );
        // Side faces from z = 0 to z = h; graded in z towards the top edge only (the
        // mid-plane z = 0 is not an edge of the full box).
        face(c + u * a - v * b, v * (2.0 * b), w * h, false);
        face(c - u * a + v * b, v * (-2.0 * b), w * h, false);
        face(c - u * a - v * b, u * (2.0 * a), w * h, false);
        face(c + u * a + v * b, u * (-2.0 * a), w * h, false);
        out
    }
}

fn divisions(length: f64, size: f64) -> usize {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = (length / size).ceil() as usize;
    n.max(2)
}

/// `n + 1` parameters in [0, 1]: Chebyshev-graded towards both ends, or towards 1 only.
fn graded(n: usize, both: bool) -> Vec<f64> {
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    (0..=n)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let kf = k as f64;
            if both {
                0.5 * (1.0 - libm::cos(std::f64::consts::PI * kf / nf))
            } else {
                // Grade towards 1 (the top edge): sin mapping.
                libm::sin(0.5 * std::f64::consts::PI * kf / nf)
            }
        })
        .collect()
}

fn mirror(t: &Triangle) -> Triangle {
    let m = |p: DVec3| DVec3::new(p.x, p.y, -p.z);
    // Reverse the orientation so the normal points outwards again.
    Triangle {
        a: m(t.a),
        b: m(t.c),
        c: m(t.b),
    }
}

/// Potential and gradient integrals of a panel and its mirror image at `x` (exact).
fn pair_integrals(t: &Panel, mirror_panel: &Panel, x: DVec3) -> (f64, DVec3) {
    let (p1, g1) = t.integrals(x);
    if x.z == 0.0 {
        // The mirror's contribution in the plane: same potential, E_z reversed.
        (2.0 * p1, DVec3::new(2.0 * g1.x, 2.0 * g1.y, 0.0))
    } else {
        let (p2, g2) = mirror_panel.integrals(x);
        (p1 + p2, g1 + g2)
    }
}

/// LU factorization with partial pivoting (row-major, in place).
#[derive(Debug)]
struct Lu {
    n: usize,
    a: Vec<f64>,
    perm: Vec<usize>,
}

impl Lu {
    #[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
    fn new(mut a: Vec<f64>, n: usize) -> Self {
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            crate::cancel::checkpoint();
            let p = (k..n)
                .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
                .expect("non-empty");
            if p != k {
                for j in 0..n {
                    a.swap(k * n + j, p * n + j);
                }
                perm.swap(k, p);
            }
            let pivot = a[k * n + k];
            for i in k + 1..n {
                let f = a[i * n + k] / pivot;
                a[i * n + k] = f;
                if f != 0.0 {
                    let (top, bottom) = a.split_at_mut(i * n);
                    let row_k = &top[k * n..k * n + n];
                    let row_i = &mut bottom[..n];
                    for j in k + 1..n {
                        row_i[j] -= f * row_k[j];
                    }
                }
            }
        }
        Self { n, a, perm }
    }

    #[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = self.perm.iter().map(|&i| b[i]).collect();
        for i in 0..n {
            let mut s = x[i];
            for j in 0..i {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = x[i];
            for j in i + 1..n {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s / self.a[i * n + i];
        }
        x
    }
}

/// Geometry-dependent part: panels (upper halves), their electrode, the LU of the
/// collocation matrix, the unit systems and the capacitance matrix.
#[derive(Debug)]
struct Geometry {
    panels: Vec<Triangle>,
    /// Precomputed panels and their mirror images.
    pre: Vec<(Panel, Panel)>,
    /// For pictures, per panel: centroid, squared radius of exact integration
    /// (`PICTURE_NEAR` sizes), the points of `PICTURE_RULE` and a third of the area.
    picture: Vec<PictureRule>,
    owner: Vec<usize>,
    lu: Lu,
    unit: Vec<Vec<f64>>,
    capacitance: Vec<Vec<f64>>,
    cost: crate::field::SetupCost,
}

fn geometry_key(electrodes: &[BoxElectrode], size: f64) -> Vec<u64> {
    let mut k: Vec<u64> = electrodes
        .iter()
        .flat_map(|e| {
            [
                e.center.x.to_bits(),
                e.center.y.to_bits(),
                e.angle.to_bits(),
                e.half_length.to_bits(),
                e.half_thickness.to_bits(),
                e.half_height.to_bits(),
            ]
        })
        .collect();
    k.push(size.to_bits());
    k
}

/// A panel prepared for the picture evaluation (`Electrodes::picture_sample`).
#[derive(Debug)]
struct PictureRule {
    centroid: DVec3,
    near2: f64,
    points: [DVec3; 3],
    third: f64,
}

impl PictureRule {
    fn new(t: &Triangle) -> Self {
        let near = PICTURE_NEAR * t.size();
        Self {
            centroid: t.centroid(),
            near2: near * near,
            points: PICTURE_RULE.map(|[u, v, k]| t.a * u + t.b * v + t.c * k),
            third: t.area() / 3.0,
        }
    }
}

fn cached_geometry(electrodes: &[BoxElectrode], size: f64) -> Arc<Geometry> {
    static CACHE: OnceLock<Mutex<HashMap<Vec<u64>, Arc<Geometry>>>> = OnceLock::new();
    let key = geometry_key(electrodes, size);
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(g) = cache.lock().expect("cache").get(&key) {
        return g.clone();
    }
    let start = std::time::Instant::now();
    let mut panels = Vec::new();
    let mut owner = Vec::new();
    for (i, e) in electrodes.iter().enumerate() {
        let m = e.mesh(size);
        owner.extend(std::iter::repeat_n(i, m.len()));
        panels.extend(m);
    }
    let n = panels.len();
    let pre: Vec<(Panel, Panel)> = panels
        .iter()
        .map(|t| (Panel::new(*t), Panel::new(mirror(t))))
        .collect();
    let mut a = vec![0.0; n * n];
    for (i, pi) in panels.iter().enumerate() {
        crate::cancel::checkpoint();
        let x = pi.centroid();
        for (j, (p, m)) in pre.iter().enumerate() {
            a[i * n + j] = pair_integrals(p, m, x).0;
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
        panels
            .iter()
            .zip(sigma)
            .zip(&owner)
            .filter(|(_, o)| **o == e)
            .map(|((t, s), _)| 2.0 * s * t.area())
            .sum()
    };
    let capacitance: Vec<Vec<f64>> = (0..m)
        .map(|i| (0..m).map(|j| charge_on(&unit[j], i)).collect())
        .collect();
    let cost = crate::field::SetupCost {
        seconds: start.elapsed().as_secs_f64(),
        bytes: 8 * lu.a.len(),
        unknowns: n,
    };
    let picture = panels.iter().map(PictureRule::new).collect();
    let g = Arc::new(Geometry {
        panels,
        pre,
        picture,
        owner,
        lu,
        unit,
        capacitance,
        cost,
    });
    let mut c = cache.lock().expect("cache");
    if c.len() > 32 {
        c.clear();
    }
    c.insert(key, g.clone());
    g
}

/// Box electrodes with their surface charge for given fixed sources.
#[derive(Clone, Debug, Default)]
pub struct Electrodes {
    pub electrodes: Vec<BoxElectrode>,
    geometry: Option<Arc<Geometry>>,
    /// Surface charge density of each (upper) panel.
    pub sigma: Vec<f64>,
    /// Potential of each electrode.
    pub potentials: Vec<f64>,
    /// Evaluate for pictures (display resolution: field lines, the field views): the
    /// panels nearer than `PICTURE_NEAR` of their sizes exactly, the others by a
    /// three-point rule (`picture_sample`). Physics always uses the exact integrals.
    pub picture: bool,
}

/// Panels nearer than this many of their sizes are integrated exactly in pictures; the
/// three-point rule for the others is then within 4.3e-6 of the plate potential of the
/// exact integrals (level 19, from 0.05 cells off the plates; with one point charge per
/// panel it was 4.7e-2, and 1.1e-3 with the near ones exact).
pub const PICTURE_NEAR: f64 = 2.0;

/// Barycentric points of the three-point rule on a triangle (degree 2, equal weights).
pub const PICTURE_RULE: [[f64; 3]; 3] = [
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];

impl Electrodes {
    /// Build cost of the geometry's factorization (zero without electrodes).
    pub fn setup_cost(&self) -> crate::field::SetupCost {
        self.geometry
            .as_ref()
            .map_or_else(Default::default, |g| g.cost)
    }

    pub fn new(
        electrodes: Vec<BoxElectrode>,
        sources: &[(DVec3, f64)],
        resolution: Resolution,
    ) -> Self {
        let mut e = Self::with_panel_size(electrodes, sources, resolution.panel_size());
        e.picture = resolution == Resolution::Display;
        e
    }

    /// As `new`, with an explicit largest panel size (cells).
    pub fn with_panel_size(
        electrodes: Vec<BoxElectrode>,
        sources: &[(DVec3, f64)],
        size: f64,
    ) -> Self {
        let m = electrodes.len();
        if m == 0 {
            return Self::default();
        }
        let geo = cached_geometry(&electrodes, size);
        // Source system: electrodes at 0, i.e. σ cancels the sources' potential.
        let b: Vec<f64> = geo
            .panels
            .iter()
            .map(|t| {
                let x = t.centroid();
                -sources
                    .iter()
                    .map(|&(p, q)| q / (x - p).length())
                    .sum::<f64>()
            })
            .collect();
        let mut sigma = geo.lu.solve(&b);
        let charge_on = |sigma: &[f64], e: usize| -> f64 {
            geo.panels
                .iter()
                .zip(sigma)
                .zip(&geo.owner)
                .filter(|(_, o)| **o == e)
                .map(|((t, s), _)| 2.0 * s * t.area())
                .sum()
        };
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
                    mat[r * k + c] = geo.capacitance[i][j];
                }
                let Bias::Charge(q) = electrodes[i].bias else {
                    unreachable!()
                };
                let known: f64 = (0..m)
                    .filter(|j| !floating.contains(j))
                    .map(|j| geo.capacitance[i][j] * alpha[j])
                    .sum();
                rhs[r] = q - g[i] - known;
            }
            let x = Lu::new(mat, k).solve(&rhs);
            for (r, &i) in floating.iter().enumerate() {
                alpha[i] = x[r];
            }
        }
        for (j, &a) in alpha.iter().enumerate() {
            if a != 0.0 {
                for (s, u) in sigma.iter_mut().zip(&geo.unit[j]) {
                    *s += a * u;
                }
            }
        }
        Self {
            electrodes,
            geometry: Some(geo),
            sigma,
            potentials: alpha,
            picture: false,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.electrodes.is_empty()
    }

    /// Geometry only (no charges): for obstacles and containment tests.
    pub fn shapes_only(electrodes: Vec<BoxElectrode>) -> Self {
        Self {
            electrodes,
            ..Self::default()
        }
    }

    /// Panels (upper halves) with their densities.
    pub fn panels(&self) -> impl Iterator<Item = (&Triangle, f64)> + '_ {
        self.geometry
            .iter()
            .flat_map(|g| g.panels.iter())
            .zip(self.sigma.iter().copied())
    }

    /// Net charge of electrode `e` (both halves).
    pub fn charge(&self, e: usize) -> f64 {
        let Some(g) = &self.geometry else { return 0.0 };
        g.panels
            .iter()
            .zip(&self.sigma)
            .zip(&g.owner)
            .filter(|(_, o)| **o == e)
            .map(|((t, s), _)| 2.0 * s * t.area())
            .sum()
    }

    /// Obstacles: each electrode's box plus a contact shell of `margin`.
    pub fn obstacles(&self, margin: f64) -> Vec<crate::geometry::Shape> {
        self.electrodes
            .iter()
            .map(|e| {
                crate::geometry::Shape::Box(crate::geometry::OrientedBox {
                    center: e.center,
                    axis: e.axes().0,
                    half: DVec3::new(e.half_length, e.half_thickness, e.half_height),
                    margin,
                })
            })
            .collect()
    }

    /// Whether `x` is inside an electrode (plus `margin`).
    pub fn contains(&self, x: DVec3, margin: f64) -> bool {
        self.electrodes.iter().any(|e| {
            let (u, v) = e.axes();
            let d = x - e.center;
            d.dot(u).abs() <= e.half_length + margin
                && d.dot(v).abs() <= e.half_thickness + margin
                && d.z.abs() <= e.half_height + margin
        })
    }
}

impl Electrodes {
    /// The picture evaluation (`picture`): each panel and its mirror image exactly when
    /// `x` is nearer than `PICTURE_NEAR` of its size to its centroid, otherwise by the
    /// three-point rule `PICTURE_RULE` (error of order (size/R)³ of the panel's share).
    fn picture_sample(
        &self,
        geo: &Geometry,
        x: DVec3,
        sigma: impl Fn(usize) -> f64,
    ) -> (f64, DVec3) {
        let mut phi = 0.0;
        let mut e = DVec3::ZERO;
        // In the plane the mirror image's points are as far as the panel's: doubled.
        let in_plane = x.z == 0.0;
        for (i, ((p, m), rule)) in geo.pre.iter().zip(&geo.picture).enumerate() {
            let s = sigma(i);
            if (x - rule.centroid).length_squared() < rule.near2 {
                let (pp, g) = pair_integrals(p, m, x);
                phi += s * pp;
                e -= g * s;
                continue;
            }
            let w = s * rule.third;
            for y in rule.points {
                let mut add = |y: DVec3, w: f64| {
                    let d = x - y;
                    let inv = 1.0 / d.length();
                    phi += w * inv;
                    e += d * (w * inv * inv * inv);
                };
                if in_plane {
                    add(y, 2.0 * w);
                } else {
                    add(y, w);
                    add(DVec3::new(y.x, y.y, -y.z), w);
                }
            }
        }
        if in_plane {
            // The mirror symmetry: no field across the plane.
            e.z = 0.0;
        }
        (phi, e)
    }

    /// Field and potential of the panels with the densities `sigma(i)`: exact integrals,
    /// or the picture evaluation.
    fn sample_densities(&self, x: DVec3, sigma: impl Fn(usize) -> f64) -> FieldSample {
        let mut phi = 0.0;
        let mut e = DVec3::ZERO;
        if let Some(geo) = &self.geometry
            && self.picture
        {
            (phi, e) = self.picture_sample(geo, x, sigma);
        } else if let Some(geo) = &self.geometry {
            for (i, (p, m)) in geo.pre.iter().enumerate() {
                let s = sigma(i);
                let (pp, g) = pair_integrals(p, m, x);
                phi += s * pp;
                e -= g * s;
            }
        }
        FieldSample {
            e,
            b: DVec3::ZERO,
            phi,
        }
    }

    /// As `sample`, with the potential of electrode `e` changed by `d` for each `(e, d)`
    /// of `shifts`: its unit system's surface charge times `d` added (the solution is
    /// linear in the potentials, §2.7). For electrodes driven by a circuit (`drive.rs`).
    pub fn sample_shifted(&self, x: DVec3, shifts: &[(usize, f64)]) -> FieldSample {
        let Some(geo) = &self.geometry else {
            return FieldSample::default();
        };
        self.sample_densities(x, |i| {
            let mut s = self.sigma[i];
            for &(e, d) in shifts {
                s += d * geo.unit[e][i];
            }
            s
        })
    }

    /// The coefficients of capacitance `C[i][j]`: the charge on electrode i with electrode
    /// j at potential 1 and the others grounded (Maxwell's matrix, §2.8). None without
    /// electrodes.
    pub fn capacitance(&self) -> Option<&[Vec<f64>]> {
        self.geometry.as_ref().map(|g| g.capacitance.as_slice())
    }
}

impl FieldSolver for Electrodes {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        self.sample_densities(x, |i| self.sigma[i])
    }
}
