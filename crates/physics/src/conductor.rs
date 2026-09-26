//! Conducting spheres (PHYSICS.md §2.6).
//!
//! **Building block: Dirichlet systems.** Every induced charge distribution is computed
//! for spheres held at prescribed potentials:
//! 1. Kelvin images to depth `D`: the image of a charge `q` at `p` in a grounded sphere
//!    (radius `a`, centre `c`) is `−q a/|p − c|` at `c + (p − c) a²/|p − c|²`; images are
//!    imaged again in the other spheres. This carries the singular near field exactly;
//!    what is left is `O(ρ^(D+1))`, `ρ` the largest image ratio between two spheres;
//! 2. the small smooth remainder is represented by the method of fundamental solutions
//!    (MFS): `K` equivalent charges on a shell of radius `s·a` inside each sphere, fitted
//!    by least squares to the boundary condition on `2K` surface points. The matrix
//!    depends only on the geometry: its QR factorization is computed once and cached.
//!
//! The MFS remainder is a model with a finite error (about 1e-9 of the source potential
//! on the surfaces). Preview and verification use two resolutions (`Resolution`), so
//! that this error enters the verification of every flight.
//!
//! **Bias.** The source system `G` (sources, all spheres grounded) and the unit systems
//! `U_j` (sphere `j` at potential 1, the others at 0; its seed is a charge `a_j` at the
//! centre, imaged like a source) are combined as `G + Σ α_j U_j`: sphere `i` then has
//! potential `α_i` and net charge `G_i + Σ_j C_ij α_j`, where `C` is the capacitance
//! matrix. Grounded: `α = 0`; fixed potential `V`: `α = V`; floating with net charge
//! `Q`: the `α` of the floating spheres solve the charge equations.
//!
//! **Accuracy.** The boundary residual (largest deviation of the surface potentials from
//! their values, relative to the sources' potential on the surfaces) is measured on an
//! independent, denser set of surface points. For the Dirichlet problem it bounds the
//! potential error everywhere outside the conductors (maximum principle).
//!
//! **The particle.** The charges induced by the moving particle are its Kelvin images to
//! depth `SELF_IMAGE_DEPTH` (grounded), corrected for floating spheres with the unit
//! systems and `C⁻¹`. Truncating by depth (not by size) and the symmetric `C` keep the
//! induced interaction symmetric, so the image force is conservative with energy
//! `½ q φ_self`; the neglected deeper images are bounded by `ρ^D` relative to the first.
//! Conductors respond instantaneously (electrostatics): exact for `c = ∞`, and for
//! finite `c` valid while the particle is slow compared with light over the size of the
//! setup.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use glam::DVec3;

use crate::field::{Coulomb, FieldSolver, FixedCharge};

/// How a conducting sphere is held.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Bias {
    /// Connected to ground: potential 0.
    Grounded,
    /// Isolated, with this net charge.
    Charge(f64),
    /// Held at this potential (relative to infinity) by a source.
    Potential(f64),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereConductor {
    pub center: DVec3,
    pub radius: f64,
    pub bias: Bias,
}

/// Point charges (position, charge).
pub type Charges = Vec<(DVec3, f64)>;

/// Resolution of the conductor model. Preview and verification use different ones, so
/// that the model's own error enters the verification of every flight (the difference
/// between the two is compared with the margins, like the integration error).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Resolution {
    Preview,
    Verify,
    /// Coarse, for pictures only (potential map, field lines): accurate to about 1e-4.
    Display,
}

impl Resolution {
    /// (equivalent charges per sphere, shell radius / sphere radius, image depth).
    fn parameters(self) -> (usize, f64, usize) {
        match self {
            Resolution::Preview => (440, 0.6, 7),
            Resolution::Verify => (600, 0.6, 8),
            Resolution::Display => (60, 0.6, 2),
        }
    }
}
/// Depth of the particle's image tree.
pub const SELF_IMAGE_DEPTH: usize = 6;

/// Kelvin images of `sources` in the spheres, to `depth` generations (a source is not
/// imaged in the sphere it lies in). Returns the images and the image charge in each
/// sphere.
pub fn kelvin_images(
    spheres: &[SphereConductor],
    sources: &[(DVec3, f64, Option<usize>)],
    depth: usize,
) -> (Charges, Vec<f64>) {
    let mut images = Vec::new();
    let mut totals = vec![0.0; spheres.len()];
    let mut generation: Vec<(DVec3, f64, Option<usize>)> = sources.to_vec();
    for _ in 0..depth {
        let mut next = Vec::new();
        for &(p, q, origin) in &generation {
            for (s, sphere) in spheres.iter().enumerate() {
                if Some(s) == origin {
                    continue;
                }
                let d = p - sphere.center;
                let r2 = d.length_squared();
                let a = sphere.radius;
                let qi = -q * a / r2.sqrt();
                let pi = sphere.center + d * (a * a / r2);
                totals[s] += qi;
                images.push((pi, qi));
                next.push((pi, qi, Some(s)));
            }
        }
        generation = next;
    }
    (images, totals)
}

/// Largest ratio |image| / |charge| for a charge in one sphere imaged in another:
/// `a_j / (|c_i − c_j| − a_i)`.
pub fn image_ratio(spheres: &[SphereConductor]) -> f64 {
    let mut rho: f64 = 0.0;
    for (i, a) in spheres.iter().enumerate() {
        for (j, b) in spheres.iter().enumerate() {
            if i != j {
                rho = rho.max(b.radius / ((a.center - b.center).length() - a.radius));
            }
        }
    }
    rho
}

/// `n` nearly uniform unit vectors (Fibonacci lattice), rotated by `twist`.
fn fibonacci(n: usize, twist: f64) -> Vec<DVec3> {
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let (fi, fn_) = (i as f64, n as f64);
            let z = 1.0 - 2.0 * (fi + 0.5) / fn_;
            let r = (1.0 - z * z).sqrt();
            let phi = golden * fi + twist;
            DVec3::new(r * phi.cos(), r * phi.sin(), z)
        })
        .collect()
}

/// MFS least-squares fit for Dirichlet data: equivalent charges, collocation points and
/// the Householder QR factorization of the matrix `1/|x_r − y_k|`.
#[derive(Debug)]
struct Fit {
    charges: Vec<(DVec3, usize)>,
    points: Vec<(DVec3, usize)>,
    depth: usize,
    qr: Vec<f64>,
    tau: Vec<f64>,
    rows: usize,
    cols: usize,
}

impl Fit {
    fn new(spheres: &[SphereConductor], resolution: Resolution) -> Self {
        let (k_charges, shell, depth) = resolution.parameters();
        let mut charges = Vec::new();
        let mut points = Vec::new();
        for (j, s) in spheres.iter().enumerate() {
            for u in fibonacci(k_charges, 0.0) {
                charges.push((s.center + u * (shell * s.radius), j));
            }
            for u in fibonacci(2 * k_charges, 0.37) {
                points.push((s.center + u * s.radius, j));
            }
        }
        let (rows, cols) = (points.len(), charges.len());
        let mut a = vec![0.0; rows * cols];
        for (k, &(y, _)) in charges.iter().enumerate() {
            for (r, &(x, _)) in points.iter().enumerate() {
                a[k * rows + r] = 1.0 / (x - y).length();
            }
        }
        let tau = householder_qr(&mut a, rows, cols);
        Self {
            charges,
            points,
            depth,
            qr: a,
            tau,
            rows,
            cols,
        }
    }

    /// Least-squares solution for right-hand side `b`.
    #[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
    fn solve(&self, mut b: Vec<f64>) -> Vec<f64> {
        let (m, n) = (self.rows, self.cols);
        for k in 0..n {
            let col = &self.qr[k * m..(k + 1) * m];
            let mut s = b[k];
            for i in k + 1..m {
                s += col[i] * b[i];
            }
            s *= self.tau[k];
            b[k] -= s;
            for i in k + 1..m {
                b[i] -= s * col[i];
            }
        }
        let mut x = vec![0.0; n];
        for k in (0..n).rev() {
            let mut s = b[k];
            for j in k + 1..n {
                s -= self.qr[j * m + k] * x[j];
            }
            x[k] = s / self.qr[k * m + k];
        }
        x
    }

    /// Dirichlet system: induced charges for `sources` (with the sphere each lies in, if
    /// any) such that sphere `i` is at potential `targets[i]`.
    fn system(
        &self,
        spheres: &[SphereConductor],
        sources: &[(DVec3, f64, Option<usize>)],
        targets: &[f64],
    ) -> Charges {
        let (mut induced, _) = kelvin_images(spheres, sources, self.depth);
        let mut known: Charges = sources.iter().map(|&(p, q, _)| (p, q)).collect();
        known.extend(induced.iter().copied());
        let b: Vec<f64> = self
            .points
            .iter()
            .map(|&(x, j)| {
                targets[j]
                    - known
                        .iter()
                        .map(|&(p, q)| q / (x - p).length())
                        .sum::<f64>()
            })
            .collect();
        let w = self.solve(b);
        induced.extend(self.charges.iter().zip(w).map(|(&(p, _), q)| (p, q)));
        induced
    }
}

/// In-place Householder QR of a column-major `m × n` matrix; returns the reflector
/// scalars `τ` (LAPACK convention, `v₀ = 1`).
fn householder_qr(a: &mut [f64], m: usize, n: usize) -> Vec<f64> {
    let mut tau = vec![0.0; n];
    for k in 0..n {
        let col = k * m;
        let norm = (k..m).map(|i| a[col + i] * a[col + i]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[col + k] > 0.0 { -norm } else { norm };
        let v0 = a[col + k] - alpha;
        for i in k + 1..m {
            a[col + i] /= v0;
        }
        tau[k] = -v0 / alpha;
        a[col + k] = alpha;
        for j in k + 1..n {
            let cj = j * m;
            let mut s = a[cj + k];
            for i in k + 1..m {
                s += a[col + i] * a[cj + i];
            }
            s *= tau[k];
            a[cj + k] -= s;
            for i in k + 1..m {
                a[cj + i] -= s * a[col + i];
            }
        }
    }
    tau
}

/// Solves `m x = b` for a small dense system (Gaussian elimination, partial pivoting).
#[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
fn solve_dense(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| m[i][k].abs().total_cmp(&m[j][k].abs()))
            .expect("non-empty");
        m.swap(k, p);
        b.swap(k, p);
        for i in k + 1..n {
            let f = m[i][k] / m[k][k];
            for j in k..n {
                m[i][j] -= f * m[k][j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| m[k][j] * x[j]).sum();
        x[k] = (b[k] - s) / m[k][k];
    }
    x
}

fn charge_in(spheres: &[SphereConductor], charges: &Charges, i: usize) -> f64 {
    let s = &spheres[i];
    charges
        .iter()
        .filter(|(p, _)| (*p - s.center).length() < s.radius)
        .map(|&(_, q)| q)
        .sum()
}

/// Geometry-dependent part: the MFS factorization, the unit systems and the capacitance
/// matrix. Cached by geometry.
#[derive(Debug)]
struct Geometry {
    fit: Fit,
    unit: Vec<Charges>,
    capacitance: Vec<Vec<f64>>,
    cost: crate::field::SetupCost,
}

fn geometry_key(spheres: &[SphereConductor], resolution: Resolution) -> Vec<u64> {
    let mut key: Vec<u64> = spheres
        .iter()
        .flat_map(|s| {
            [
                s.center.x.to_bits(),
                s.center.y.to_bits(),
                s.center.z.to_bits(),
                s.radius.to_bits(),
            ]
        })
        .collect();
    key.push(resolution as u64);
    key
}

fn cached_geometry(spheres: &[SphereConductor], resolution: Resolution) -> Arc<Geometry> {
    static CACHE: OnceLock<Mutex<HashMap<Vec<u64>, Arc<Geometry>>>> = OnceLock::new();
    let key = geometry_key(spheres, resolution);
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(g) = cache.lock().expect("cache").get(&key) {
        return g.clone();
    }
    let start = std::time::Instant::now();
    let fit = Fit::new(spheres, resolution);
    let n = spheres.len();
    let unit: Vec<Charges> = (0..n)
        .map(|j| {
            let mut targets = vec![0.0; n];
            targets[j] = 1.0;
            let seed = (spheres[j].center, spheres[j].radius, Some(j));
            let mut c = fit.system(spheres, &[seed], &targets);
            c.push((spheres[j].center, spheres[j].radius));
            c
        })
        .collect();
    let capacitance: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| charge_in(spheres, &unit[j], i)).collect())
        .collect();
    let cost = crate::field::SetupCost {
        seconds: start.elapsed().as_secs_f64(),
        bytes: 8 * (fit.qr.len() + fit.tau.len()),
        unknowns: fit.cols,
    };
    let g = Arc::new(Geometry {
        fit,
        unit,
        capacitance,
        cost,
    });
    let mut c = cache.lock().expect("cache");
    if c.len() > 64 {
        c.clear();
    }
    c.insert(key, g.clone());
    g
}

/// Conductors of a scene with their induced charges for the fixed sources, and what is
/// needed to compute the charges induced by a moving particle.
#[derive(Clone, Debug, Default)]
pub struct Conductors {
    pub spheres: Vec<SphereConductor>,
    /// Induced charges for the fixed sources and the bias.
    pub induced: Coulomb,
    geometry: Option<Arc<Geometry>>,
    /// Indices of floating spheres and the inverse of their capacitance block.
    floating: Vec<usize>,
    floating_inverse: Vec<Vec<f64>>,
    /// Potential of each sphere.
    pub alpha: Vec<f64>,
    /// Relative bound on the neglected part of the particle's image tree.
    pub self_truncation: f64,
}

impl Conductors {
    /// Induced charges of `spheres` for fixed point `sources` (outside all spheres).
    pub fn new(
        spheres: Vec<SphereConductor>,
        sources: &[(DVec3, f64)],
        resolution: Resolution,
    ) -> Self {
        let n = spheres.len();
        if n == 0 {
            return Self::default();
        }
        let geo = cached_geometry(&spheres, resolution);
        let src: Vec<(DVec3, f64, Option<usize>)> =
            sources.iter().map(|&(p, q)| (p, q, None)).collect();
        let grounded = geo.fit.system(&spheres, &src, &vec![0.0; n]);
        let g: Vec<f64> = (0..n).map(|i| charge_in(&spheres, &grounded, i)).collect();

        // α: the potentials of the spheres.
        let floating: Vec<usize> = (0..n)
            .filter(|&i| matches!(spheres[i].bias, Bias::Charge(_)))
            .collect();
        let mut alpha: Vec<f64> = spheres
            .iter()
            .map(|s| match s.bias {
                Bias::Potential(v) => v,
                _ => 0.0,
            })
            .collect();
        let block: Vec<Vec<f64>> = floating
            .iter()
            .map(|&i| floating.iter().map(|&j| geo.capacitance[i][j]).collect())
            .collect();
        let floating_inverse = invert(&block);
        if !floating.is_empty() {
            let rhs: Vec<f64> = floating
                .iter()
                .map(|&i| {
                    let Bias::Charge(q) = spheres[i].bias else {
                        unreachable!()
                    };
                    let known: f64 = (0..n)
                        .filter(|j| !floating.contains(j))
                        .map(|j| geo.capacitance[i][j] * alpha[j])
                        .sum();
                    q - g[i] - known
                })
                .collect();
            for (k, a) in mat_vec(&floating_inverse, &rhs).into_iter().enumerate() {
                alpha[floating[k]] = a;
            }
        }
        let mut induced = grounded;
        for (j, &a) in alpha.iter().enumerate() {
            if a != 0.0 {
                induced.extend(geo.unit[j].iter().map(|&(p, q)| (p, q * a)));
            }
        }

        // The source system and the unit systems share their equivalent-charge
        // positions: merge charges at identical positions (first-occurrence order, so
        // the summation order stays deterministic).
        let induced = merge(induced);
        let fixed: Vec<FixedCharge> = induced
            .iter()
            .map(|&(position, charge)| FixedCharge {
                position,
                charge,
                radius: 0.0,
            })
            .collect();
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let self_truncation = image_ratio(&spheres).powi(SELF_IMAGE_DEPTH as i32);
        Self {
            spheres,
            induced: Coulomb::new(&fixed),
            geometry: Some(geo),
            floating,
            floating_inverse,
            alpha,
            self_truncation,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.spheres.is_empty()
    }

    /// Build cost of the geometry's factorization (zero without spheres).
    pub fn setup_cost(&self) -> crate::field::SetupCost {
        self.geometry
            .as_ref()
            .map_or_else(Default::default, |g| g.cost)
    }

    /// Largest deviation of the surface potentials from their values, measured on 1000
    /// independent points per sphere, relative to the largest potential of the sources
    /// and bias on the surfaces (PHYSICS.md §2.6). Expensive; for checks, not per step.
    pub fn boundary_residual(&self, sources: &[(DVec3, f64)]) -> f64 {
        let mut all: Charges = sources.to_vec();
        all.extend(self.induced.charges());
        let potential = |x: DVec3| all.iter().map(|&(p, q)| q / (x - p).length()).sum::<f64>();
        let mut scale: f64 = self.alpha.iter().map(|a| a.abs()).fold(0.0, f64::max);
        let mut residual: f64 = 0.0;
        for (j, s) in self.spheres.iter().enumerate() {
            for u in fibonacci(1000, 2.3) {
                let x = s.center + u * s.radius;
                residual = residual.max((potential(x) - self.alpha[j]).abs());
                let v: f64 = sources
                    .iter()
                    .map(|&(p, q)| (q / (x - p).length()).abs())
                    .sum();
                scale = scale.max(v);
            }
        }
        residual / scale.max(1e-300)
    }

    /// Charges induced by a point charge `q` at `x`: its Kelvin images to depth
    /// `SELF_IMAGE_DEPTH` (grounded spheres), plus unit systems for the floating spheres.
    /// The charge a unit charge at `x` induces on grounded sphere `j` is, by Green's
    /// reciprocity, `−φ_Uj(x)` (the potential at `x` of unit system `j`); using it makes
    /// the floating correction `q φ_U(x)ᵀ C⁻¹ φ_U(x)` exactly symmetric, so the image
    /// force stays conservative.
    pub fn induced_by(&self, x: DVec3, q: f64) -> Charges {
        let (mut charges, _) = kelvin_images(&self.spheres, &[(x, q, None)], SELF_IMAGE_DEPTH);
        if let Some(geo) = &self.geometry
            && !self.floating.is_empty()
        {
            let rhs: Vec<f64> = self
                .floating
                .iter()
                .map(|&i| {
                    let phi: f64 = geo.unit[i].iter().map(|&(p, c)| c / (x - p).length()).sum();
                    q * phi
                })
                .collect();
            for (k, a) in mat_vec(&self.floating_inverse, &rhs)
                .into_iter()
                .enumerate()
            {
                let j = self.floating[k];
                charges.extend(geo.unit[j].iter().map(|&(p, c)| (p, c * a)));
            }
        }
        charges
    }

    /// Field and potential at `x` of the charges induced by a charge `q` at `x` itself
    /// (the image force is `q E`; the interaction energy is `½ q φ`).
    pub fn self_field(&self, x: DVec3, q: f64) -> (DVec3, f64) {
        if self.spheres.is_empty() || q == 0.0 {
            return (DVec3::ZERO, 0.0);
        }
        let mut e = DVec3::ZERO;
        let mut phi = 0.0;
        for (p, c) in self.induced_by(x, q) {
            let d = x - p;
            let inv_r = 1.0 / d.length();
            phi += c * inv_r;
            e += d * (c * inv_r * inv_r * inv_r);
        }
        (e, phi)
    }
}

/// Sums charges at bit-identical positions, keeping first-occurrence order.
fn merge(charges: Charges) -> Charges {
    let mut index: HashMap<[u64; 3], usize> = HashMap::new();
    let mut out: Charges = Vec::with_capacity(charges.len());
    for (p, q) in charges {
        let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        match index.get(&key) {
            Some(&i) => out[i].1 += q,
            None => {
                index.insert(key, out.len());
                out.push((p, q));
            }
        }
    }
    out
}

fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum())
        .collect()
}

fn invert(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = m.len();
    let columns: Vec<Vec<f64>> = (0..n)
        .map(|k| {
            let mut e = vec![0.0; n];
            e[k] = 1.0;
            solve_dense(m.to_vec(), e)
        })
        .collect();
    (0..n)
        .map(|i| (0..n).map(|k| columns[k][i]).collect())
        .collect()
}

impl FieldSolver for Conductors {
    fn sample(&self, x: DVec3, t: f64) -> crate::field::FieldSample {
        self.induced.sample(x, t)
    }
}
