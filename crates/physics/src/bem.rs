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
//! **Dielectrics.** Boxes of relative permittivity ε (`Dielectric`) are meshed the same
//! way and solved in the same system: their panels carry the bound (polarization) charge,
//! fixed by the continuity of the normal displacement across the surface. With
//! `λ = (ε − 1)/(ε + 1)` and the average normal field `Eₙ` at the centroid (the panel's
//! own charge contributing none there), `σ = (λ/2π) Eₙ` (k = 1; Jackson §4.4). The
//! electrodes' unit systems and capacitance matrix then include the dielectrics' response.
//! Each body's bound charge sums to zero exactly; the discretized equation alone does not
//! ensure it (the net-charge mode, a free conductor charge, is nearly singular as ε grows:
//! −0.6 % of a nearby charge at ε = 2, −48 % at ε = 10⁴). So each body adds that
//! constraint and one unknown γ, a uniform normal field on its surface (the field of a
//! monopole inside, which absorbs the inconsistency; tests N1–N3).
//! Particles do not enter dielectrics (obstacles, as electrodes).
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

/// The shape of a dielectric body: a box (as an electrode; its bias is ignored) or a
/// sphere centred in the plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BodyShape {
    Box(BoxElectrode),
    Sphere { center: DVec3, radius: f64 },
}

impl BodyShape {
    /// Upper-half surface (z ≥ 0) as triangles with outward normals.
    pub fn mesh(&self, size: f64) -> Vec<Triangle> {
        match self {
            BodyShape::Box(b) => b.mesh(size),
            BodyShape::Sphere { center, radius } => {
                sphere_mesh(*center, *radius, SPHERE_PANEL * size)
            }
        }
    }

    fn key(&self) -> Vec<u64> {
        match self {
            BodyShape::Box(e) => vec![
                0,
                e.center.x.to_bits(),
                e.center.y.to_bits(),
                e.angle.to_bits(),
                e.half_length.to_bits(),
                e.half_thickness.to_bits(),
                e.half_height.to_bits(),
            ],
            BodyShape::Sphere { center, radius } => {
                vec![1, center.x.to_bits(), center.y.to_bits(), radius.to_bits()]
            }
        }
    }
}

/// Panels on a dielectric sphere are this fraction of the resolution's size: its flat
/// facets meet at small kinks, where the normal-field condition converges only at order
/// about 1.3 (test N4: glass 2.5 % off at the verification's size), so it gets finer
/// panels than a box (about 1.3 % there, as the cube's 1.1 %).
pub const SPHERE_PANEL: f64 = 0.6;

/// A dielectric body of relative permittivity `permittivity` (≥ 0; 1 is vacuum).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dielectric {
    pub shape: BodyShape,
    pub permittivity: f64,
}

/// The upper half (z ≥ 0) of a sphere: the upper four faces of an octahedron, each cut
/// into `n²` triangles (`n` segments per edge, edges no longer than `size` on the
/// sphere), projected onto it. Nearly uniform panels, no slivers at the pole, and the
/// equator exactly at z = 0 (the mirror plane). Normals point outwards.
fn sphere_mesh(center: DVec3, radius: f64, size: f64) -> Vec<Triangle> {
    // An octahedron edge spans a quarter circle (π/2 · r) on the sphere.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = ((std::f64::consts::FRAC_PI_2 * radius / size).ceil() as usize).max(2);
    let top = DVec3::Z;
    let equator = [DVec3::X, DVec3::Y, -DVec3::X, -DVec3::Y];
    let on_sphere = |v: DVec3| center + v.normalize() * radius;
    let mut out = Vec::new();
    for k in 0..4 {
        let (a, b, c) = (top, equator[k], equator[(k + 1) % 4]);
        // Barycentric grid: point (i, j) = a + (b − a) i/n + (c − a) j/n, i + j ≤ n.
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        #[allow(clippy::cast_precision_loss)]
        let p = |i: usize, j: usize| {
            on_sphere(a + (b - a) * (i as f64 / nf) + (c - a) * (j as f64 / nf))
        };
        for i in 0..n {
            for j in 0..n - i {
                // Orientation (a, b, c) of the face (top, east, north): outward for a
                // counter-clockwise run seen from outside.
                let t1 = Triangle {
                    a: p(i, j),
                    b: p(i + 1, j),
                    c: p(i, j + 1),
                };
                out.push(orient_outward(t1, center));
                if i + j + 1 < n {
                    let t2 = Triangle {
                        a: p(i + 1, j),
                        b: p(i + 1, j + 1),
                        c: p(i, j + 1),
                    };
                    out.push(orient_outward(t2, center));
                }
            }
        }
    }
    // Flat triangles with their corners on the sphere enclose too little: their sag,
    // about h²/(8r), cost the polarizability (∝ r³) 15 % at panels of half the radius
    // (test N4). Scaled from the centre so that the polyhedron's volume is the sphere's
    // (the tetrahedra from the centre; the equator stays at z = 0).
    let half_volume: f64 = out
        .iter()
        .map(|t| (t.a - center).dot((t.b - center).cross(t.c - center)) / 6.0)
        .sum();
    let scale = libm::cbrt(std::f64::consts::FRAC_PI_3 * 2.0 * radius.powi(3) / half_volume);
    let grow = |v: DVec3| center + (v - center) * scale;
    out.iter()
        .map(|t| Triangle {
            a: grow(t.a),
            b: grow(t.b),
            c: grow(t.c),
        })
        .collect()
}

/// The triangle with its vertices ordered so that its normal points away from `center`.
fn orient_outward(t: Triangle, center: DVec3) -> Triangle {
    if t.normal().dot(t.centroid() - center) >= 0.0 {
        t
    } else {
        Triangle {
            a: t.a,
            b: t.c,
            c: t.b,
        }
    }
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
    pub fn panel_size(self) -> f64 {
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
    /// Per panel: `Some(λ)` for a dielectric's, None for an electrode's.
    lambda: Vec<Option<f64>>,
    /// For pictures, per panel: centroid, squared radius of exact integration
    /// (`PICTURE_NEAR` sizes), the points of `PICTURE_RULE` and a third of the area.
    picture: Vec<PictureRule>,
    /// The electrode each panel belongs to (`usize::MAX` for a dielectric's panel).
    owner: Vec<usize>,
    lu: Lu,
    unit: Vec<Vec<f64>>,
    capacitance: Vec<Vec<f64>>,
    cost: crate::field::SetupCost,
}

fn geometry_key(electrodes: &[BoxElectrode], dielectrics: &[Dielectric], size: f64) -> Vec<u64> {
    let shape = |e: &BoxElectrode| {
        [
            e.center.x.to_bits(),
            e.center.y.to_bits(),
            e.angle.to_bits(),
            e.half_length.to_bits(),
            e.half_thickness.to_bits(),
            e.half_height.to_bits(),
        ]
    };
    let mut k: Vec<u64> = electrodes.iter().flat_map(shape).collect();
    for d in dielectrics {
        // Marker, then the shape and the permittivity.
        k.push(u64::MAX);
        k.extend(d.shape.key());
        k.push(d.permittivity.to_bits());
    }
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

fn cached_geometry(
    electrodes: &[BoxElectrode],
    dielectrics: &[Dielectric],
    size: f64,
) -> Arc<Geometry> {
    static CACHE: OnceLock<Mutex<HashMap<Vec<u64>, Arc<Geometry>>>> = OnceLock::new();
    let key = geometry_key(electrodes, dielectrics, size);
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(g) = cache.lock().expect("cache").get(&key) {
        return g.clone();
    }
    let start = std::time::Instant::now();
    let mut panels = Vec::new();
    let mut owner = Vec::new();
    // Per panel: None for an electrode's (potential condition), Some(λ) for a
    // dielectric's (normal displacement condition).
    let mut lambda: Vec<Option<f64>> = Vec::new();
    for (i, e) in electrodes.iter().enumerate() {
        let m = e.mesh(size);
        owner.extend(std::iter::repeat_n(i, m.len()));
        lambda.extend(std::iter::repeat_n(None, m.len()));
        panels.extend(m);
    }
    // The dielectric body of each dielectric panel.
    let mut body: Vec<usize> = Vec::new();
    for (k, d) in dielectrics.iter().enumerate() {
        let m = d.shape.mesh(size);
        let l = (d.permittivity - 1.0) / (d.permittivity + 1.0);
        owner.extend(std::iter::repeat_n(usize::MAX, m.len()));
        lambda.extend(std::iter::repeat_n(Some(l), m.len()));
        body.extend(std::iter::repeat_n(k, m.len()));
        panels.extend(m);
    }
    let n = panels.len();
    let n_electrode_panels = n - body.len();
    // The bordered system: the panels' densities, then one γ per dielectric body.
    let big = n + dielectrics.len();
    let pre: Vec<(Panel, Panel)> = panels
        .iter()
        .map(|t| (Panel::new(*t), Panel::new(mirror(t))))
        .collect();
    let mut a = vec![0.0; big * big];
    for (i, pi) in panels.iter().enumerate() {
        crate::cancel::checkpoint();
        let x = pi.centroid();
        match lambda[i] {
            None => {
                for (j, (p, m)) in pre.iter().enumerate() {
                    a[i * big + j] = pair_integrals(p, m, x).0;
                }
            }
            Some(l) => {
                // σᵢ − (λ/2π) Σⱼ σⱼ Eₙ,ⱼ(xᵢ) = (λ/2π) E_src·nᵢ, with the field of a unit
                // density `−∫∇(1/R)`; the panel's own average normal field is 0 (on its
                // centroid h is a rounding error, so its normal part is removed).
                let normal = pi.normal();
                // The condition averaged over the panel (its flux), by `DIELECTRIC_RULE`:
                // centroid collocation alone converged slowly near the edges.
                let points =
                    DIELECTRIC_RULE.map(|([u, v, w], wt)| (pi.a * u + pi.b * v + pi.c * w, wt));
                for (j, (p, m)) in pre.iter().enumerate() {
                    let mut g = DVec3::ZERO;
                    for &(y, wt) in &points {
                        g += wt
                            * if j == i {
                                let (_, g1) = p.integrals(y);
                                let (_, g2) = m.integrals(y);
                                (g1 - normal * g1.dot(normal)) + g2
                            } else {
                                pair_integrals(p, m, y).1
                            };
                    }
                    a[i * big + j] = l / std::f64::consts::TAU * g.dot(normal);
                }
                a[i * big + i] += 1.0;
                // The body's γ: a uniform normal field on its surface.
                a[i * big + n + body[i - n_electrode_panels]] = 1.0;
            }
        }
    }
    // The constraints: each body's bound charge (both halves) sums to zero.
    for (j, (t, &k)) in panels[n_electrode_panels..].iter().zip(&body).enumerate() {
        a[(n + k) * big + n_electrode_panels + j] = 2.0 * t.area();
    }
    let lu = Lu::new(a, big);
    let m = electrodes.len();
    let unit: Vec<Vec<f64>> = (0..m)
        .map(|e| {
            let b: Vec<f64> = owner
                .iter()
                .map(|&o| if o == e { 1.0 } else { 0.0 })
                .collect();
            solve_bordered(&lu, &b, n)
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
        lambda,
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

/// The panels' densities from the bordered system's factorization: the right side padded
/// with the dielectric bodies' constraints (zero net charge), their γ dropped.
fn solve_bordered(lu: &Lu, b: &[f64], n: usize) -> Vec<f64> {
    let mut rhs = b.to_vec();
    rhs.resize(lu.n, 0.0);
    let mut x = lu.solve(&rhs);
    x.truncate(n);
    x
}

/// Box electrodes with their surface charge for given fixed sources.
#[derive(Clone, Debug, Default)]
pub struct Electrodes {
    pub electrodes: Vec<BoxElectrode>,
    /// Dielectric boxes, solved with the electrodes (their panels after the electrodes').
    pub dielectrics: Vec<Dielectric>,
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

/// The rule averaging a dielectric panel's condition: barycentric points and weights
/// (Dunavant's degree-4 rule, 6 points).
const DIELECTRIC_RULE: [([f64; 3], f64); 6] = {
    const A: f64 = 0.445_948_490_915_965;
    const WA: f64 = 0.223_381_589_678_011;
    const B: f64 = 0.091_576_213_509_771;
    const WB: f64 = 0.109_951_743_655_322;
    [
        ([A, A, 1.0 - 2.0 * A], WA),
        ([A, 1.0 - 2.0 * A, A], WA),
        ([1.0 - 2.0 * A, A, A], WA),
        ([B, B, 1.0 - 2.0 * B], WB),
        ([B, 1.0 - 2.0 * B, B], WB),
        ([1.0 - 2.0 * B, B, B], WB),
    ]
};

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

    /// Electrodes and dielectric boxes at a resolution.
    pub fn new_with_dielectrics(
        electrodes: Vec<BoxElectrode>,
        dielectrics: Vec<Dielectric>,
        sources: &[(DVec3, f64)],
        resolution: Resolution,
    ) -> Self {
        let mut e =
            Self::with_dielectrics(electrodes, dielectrics, sources, resolution.panel_size());
        e.picture = resolution == Resolution::Display;
        e
    }

    /// As `new`, with an explicit largest panel size (cells).
    pub fn with_panel_size(
        electrodes: Vec<BoxElectrode>,
        sources: &[(DVec3, f64)],
        size: f64,
    ) -> Self {
        Self::with_dielectrics(electrodes, Vec::new(), sources, size)
    }

    /// Electrodes and dielectric boxes, with an explicit largest panel size (cells).
    pub fn with_dielectrics(
        electrodes: Vec<BoxElectrode>,
        dielectrics: Vec<Dielectric>,
        sources: &[(DVec3, f64)],
        size: f64,
    ) -> Self {
        let m = electrodes.len();
        if m == 0 && dielectrics.is_empty() {
            return Self::default();
        }
        let geo = cached_geometry(&electrodes, &dielectrics, size);
        // Source system: electrodes at 0, i.e. σ cancels the sources' potential; on a
        // dielectric's panel `(λ/2π) E_src·n`.
        let b: Vec<f64> = geo
            .panels
            .iter()
            .zip(&geo.lambda)
            .map(|(t, l)| {
                let x = t.centroid();
                match l {
                    None => -sources
                        .iter()
                        .map(|&(p, q)| q / (x - p).length())
                        .sum::<f64>(),
                    Some(l) => {
                        // Averaged over the panel as the matrix's row.
                        let e: DVec3 = DIELECTRIC_RULE
                            .iter()
                            .map(|([u, v, w], wt)| (t.a * *u + t.b * *v + t.c * *w, *wt))
                            .flat_map(|(y, wt)| {
                                sources.iter().map(move |&(p, q)| {
                                    let r = y - p;
                                    r * (wt * q / (r.length() * r.length_squared()))
                                })
                            })
                            .sum::<DVec3>();
                        l / std::f64::consts::TAU * e.dot(t.normal())
                    }
                }
            })
            .collect();
        let mut sigma = solve_bordered(&geo.lu, &b, geo.panels.len());
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
            dielectrics,
            geometry: Some(geo),
            sigma,
            potentials: alpha,
            picture: false,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.electrodes.is_empty() && self.dielectrics.is_empty()
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

/// A ferrite body: magnetically soft (relative permeability `permeability`), electrically
/// an insulator, so it bends magnetic fields only (PHYSICS.md §2.7, "Ferrites").
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ferrite {
    pub shape: BodyShape,
    pub permeability: f64,
}

/// The geometry-dependent part of the ferrites: panels, their mirrors, the LU of the
/// collocation matrix and each panel's λ.
#[derive(Debug)]
struct FerriteGeometry {
    panels: Vec<Triangle>,
    pre: Vec<(Panel, Panel)>,
    lambda: Vec<f64>,
    lu: Lu,
    cost: crate::field::SetupCost,
}

/// Field integrals of a panel and its mirror image for a density odd under the mirror
/// (z → −z): the mirror's charge has the opposite sign. In the plane z = 0 the in-plane
/// components cancel exactly and B_z doubles.
fn odd_pair_gradient(t: &Panel, mirror_panel: &Panel, x: DVec3) -> DVec3 {
    let (_, g1) = t.integrals(x);
    if x.z == 0.0 {
        DVec3::new(0.0, 0.0, 2.0 * g1.z)
    } else {
        let (_, g2) = mirror_panel.integrals(x);
        g1 - g2
    }
}

fn ferrite_geometry(bodies: &[Ferrite], size: f64) -> Arc<FerriteGeometry> {
    static CACHE: OnceLock<Mutex<HashMap<Vec<u64>, Arc<FerriteGeometry>>>> = OnceLock::new();
    let mut key: Vec<u64> = bodies
        .iter()
        .flat_map(|b| {
            let mut k = b.shape.key();
            k.push(b.permeability.to_bits());
            k
        })
        .collect();
    key.push(size.to_bits());
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(g) = cache.lock().expect("cache").get(&key) {
        return g.clone();
    }
    let start = std::time::Instant::now();
    let mut panels = Vec::new();
    let mut lambda = Vec::new();
    for b in bodies {
        let m = b.shape.mesh(size);
        let l = (b.permeability - 1.0) / (b.permeability + 1.0);
        lambda.extend(std::iter::repeat_n(l, m.len()));
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
        let normal = pi.normal();
        let points = DIELECTRIC_RULE.map(|([u, v, w], wt)| (pi.a * u + pi.b * v + pi.c * w, wt));
        for (j, (p, m)) in pre.iter().enumerate() {
            let mut g = DVec3::ZERO;
            for &(y, wt) in &points {
                g += wt
                    * if j == i {
                        // The panel's own average normal field is 0 (its points lie on
                        // it up to rounding).
                        let (_, g1) = p.integrals(y);
                        let (_, g2) = m.integrals(y);
                        (g1 - normal * g1.dot(normal)) - g2
                    } else {
                        odd_pair_gradient(p, m, y)
                    };
            }
            a[i * n + j] = lambda[i] / std::f64::consts::TAU * g.dot(normal);
        }
        a[i * n + i] += 1.0;
    }
    let lu = Lu::new(a, n);
    let cost = crate::field::SetupCost {
        seconds: start.elapsed().as_secs_f64(),
        bytes: 8 * lu.a.len(),
        unknowns: n,
    };
    let g = Arc::new(FerriteGeometry {
        panels,
        pre,
        lambda,
        lu,
        cost,
    });
    let mut c = cache.lock().expect("cache");
    if c.len() > 32 {
        c.clear();
    }
    c.insert(key, g.clone());
    g
}

/// Ferrite bodies magnetized by static external fields: their bound magnetic surface
/// charge, from the continuity of the normal B (`σ = (λ/2π) H̄ₙ`, λ = (μ − 1)/(μ + 1),
/// the dielectric's equation with H for E), odd under the mirror z → −z because the
/// sources of the slice (in-plane coils, magnets along z) give B along z in the plane.
/// The total magnetic charge of a body is then zero by symmetry (no monopoles), and
/// the net-charge mode that the dielectrics needed bordering for does not exist.
#[derive(Clone, Debug, Default)]
pub struct Ferrites {
    pub bodies: Vec<Ferrite>,
    geometry: Option<Arc<FerriteGeometry>>,
    /// Magnetic surface charge density of each (upper) panel; the mirror carries −σ.
    pub sigma: Vec<f64>,
}

impl Ferrites {
    /// The bodies magnetized by the static field `b_ext`, with panels of at most `size`.
    pub fn new(bodies: Vec<Ferrite>, b_ext: impl Fn(DVec3) -> DVec3, size: f64) -> Self {
        if bodies.is_empty() {
            return Self::default();
        }
        let geo = ferrite_geometry(&bodies, size);
        let b: Vec<f64> = geo
            .panels
            .iter()
            .zip(&geo.lambda)
            .map(|(t, l)| {
                let h: DVec3 = DIELECTRIC_RULE
                    .iter()
                    .map(|([u, v, w], wt)| b_ext(t.a * *u + t.b * *v + t.c * *w) * *wt)
                    .sum();
                l / std::f64::consts::TAU * h.dot(t.normal())
            })
            .collect();
        let sigma = geo.lu.solve(&b);
        Self {
            bodies,
            geometry: Some(geo),
            sigma,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    /// Build cost of the factorization.
    pub fn setup_cost(&self) -> crate::field::SetupCost {
        self.geometry
            .as_ref()
            .map_or_else(Default::default, |g| g.cost)
    }

    /// Panels (upper halves) with their densities (the mirror carries the opposite).
    pub fn panels(&self) -> impl Iterator<Item = (&Triangle, f64)> + '_ {
        self.geometry
            .iter()
            .flat_map(|g| g.panels.iter())
            .zip(self.sigma.iter().copied())
    }

    /// The bodies' magnetic field at `x` (exact panel integrals).
    pub fn field(&self, x: DVec3) -> DVec3 {
        let Some(geo) = &self.geometry else {
            return DVec3::ZERO;
        };
        let mut b = DVec3::ZERO;
        for ((p, m), s) in geo.pre.iter().zip(&self.sigma) {
            b -= odd_pair_gradient(p, m, x) * *s;
        }
        b
    }
}
