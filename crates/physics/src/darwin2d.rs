//! The magnetic field of currents in the plane, in the Darwin (magnetoquasistatic) limit
//! (docs/DARWIN.md, PHYSICS.md §2.12). In the z-invariant world a current in the plane
//! makes a field along z; it and the induced field come from the transverse vector
//! potential `A = ∇ × χẑ` (Coulomb gauge), `B_z = −∇²χ`.
//!
//! Given the charges' canonical current `J_P = Σ q P/m` and their response `β = Σ q²/m` per
//! area (a carrier's kinetic current is `q (P − q A)/m`), `χ` minimises
//!
//! `W(χ) = ∫ B_z²/(2μ₀) + β |A|²/2 − J_P · A`,
//!
//! the magnetic energy plus the carriers' kinetic energy's dependence on `A`: so the
//! induced field is implicit (a metal larger than its skin depth `1/√(μ₀β)` screens it),
//! with no time derivative of the current. Its Euler–Lagrange equation,
//! `(1/μ₀) ∇⁴χ − ∇·(β ∇χ) = (∇ × J_P)_z`, is solved on a staggered grid by conjugate
//! gradients, preconditioned by the vacuum's operator `μ₀ L⁻²` (`L = −∇²`, exact through
//! sine transforms).
//!
//! The grid: `χ` and `B_z = Lχ` on the nodes, `A_x` on the vertical edges (between nodes
//! `(i, j)` and `(i, j+1)`), `A_y` on the horizontal ones, the currents and `β` deposited on
//! the same edges (so `∫ J·A` pairs exactly, and a loop of current is exactly the curl of a
//! step in the stream function: its field is `μ₀K` inside, nothing outside). The frame of
//! nodes around it holds `χ = 0` and `B_z = 0` (sine modes): flux escapes as from the
//! infinite plane, where an infinitely long solenoid's field vanishes outside it.
//!
//! Deterministic: a radix-2 FFT of our own (twiddles from `libm`), fixed orders, no fused
//! multiply-add.

/// A type-I discrete sine transform of `n` points, `n + 1` a power of two, through a
/// complex FFT of `2(n + 1)` points.
#[derive(Clone, Debug)]
struct Dst {
    n: usize,
    /// `e^{−2πik/N}`, `k < N/2`.
    twiddles: Vec<(f64, f64)>,
}

impl Dst {
    fn new(n: usize) -> Self {
        assert!(
            (n + 1).is_power_of_two(),
            "DST size {n}: n + 1 must be a power of two"
        );
        let big = 2 * (n + 1);
        #[allow(clippy::cast_precision_loss)]
        let twiddles = (0..big / 2)
            .map(|k| {
                let a = -std::f64::consts::TAU * k as f64 / big as f64;
                (libm::cos(a), libm::sin(a))
            })
            .collect();
        Self { n, twiddles }
    }

    /// `X_k = Σ_{j=1..n} x_j sin(π j k/(n+1))`, `k = 1..n` (indices from 0 in the slices).
    fn apply(&self, x: &[f64], out: &mut [f64], work: &mut [(f64, f64)]) {
        let n = self.n;
        let big = 2 * (n + 1);
        // The odd extension: y_0 = 0, y_j = x_j, y_{n+1} = 0, y_{N−j} = −x_j.
        work[..big].fill((0.0, 0.0));
        for j in 1..=n {
            work[j] = (x[j - 1], 0.0);
            work[big - j] = (-x[j - 1], 0.0);
        }
        fft(&mut work[..big], &self.twiddles);
        // FFT_k = −2i X_k.
        for k in 1..=n {
            out[k - 1] = -0.5 * work[k].1;
        }
    }
}

/// An in-place radix-2 FFT (decimation in time) with `e^{−2πik/N}` twiddles.
fn fft(a: &mut [(f64, f64)], twiddles: &[(f64, f64)]) {
    let n = a.len();
    // Bit reversal.
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let stride = n / len;
        for start in (0..n).step_by(len) {
            for k in 0..len / 2 {
                let (wr, wi) = twiddles[k * stride];
                let (xr, xi) = a[start + k + len / 2];
                let (tr, ti) = (wr * xr - wi * xi, wr * xi + wi * xr);
                let (ur, ui) = a[start + k];
                a[start + k] = (ur + tr, ui + ti);
                a[start + k + len / 2] = (ur - tr, ui - ti);
            }
        }
        len <<= 1;
    }
}

/// The Darwin field's grid and its last solution (the next solve's start).
#[derive(Clone, Debug)]
pub struct Darwin2d {
    /// Interior nodes along x and y (`nx + 1`, `ny + 1` powers of two).
    pub nx: usize,
    pub ny: usize,
    /// Node spacing.
    pub h: f64,
    /// The position of node `(0, 0)`, a frame node; the interior is `i = 1..=nx`,
    /// `j = 1..=ny`.
    pub origin: (f64, f64),
    /// `μ₀ = 4π/c²` (the game's Coulomb constant 1).
    pub mu0: f64,
    /// The stream function on the interior nodes, `(i − 1) + nx (j − 1)`.
    pub chi: Vec<f64>,
    dst_x: Dst,
    dst_y: Dst,
    /// The conjugate-gradient iterations of the last solve.
    pub iterations: usize,
}

/// Edge-centred fields: `x` on the vertical edges (`A_x`, `J_x` at `(i, j + ½)`, `i = 1..=nx`,
/// `j = 0..=ny`: index `(i − 1) + nx j`), `y` on the horizontal ones (`(i + ½, j)`,
/// `i = 0..=nx`, `j = 1..=ny`: index `i + (nx + 1)(j − 1)`).
#[derive(Clone, Debug, PartialEq)]
pub struct Edges {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Edges {
    #[must_use]
    pub fn zeros(grid: &Darwin2d) -> Self {
        Self {
            x: vec![0.0; grid.nx * (grid.ny + 1)],
            y: vec![0.0; (grid.nx + 1) * grid.ny],
        }
    }
}

impl Darwin2d {
    /// A grid of `nx × ny` interior nodes spaced `h`, its frame node `(0, 0)` at `origin`,
    /// for the speed of light `c`.
    #[must_use]
    pub fn new(nx: usize, ny: usize, h: f64, origin: (f64, f64), c: f64) -> Self {
        Self {
            nx,
            ny,
            h,
            origin,
            mu0: 4.0 * std::f64::consts::PI / (c * c),
            chi: vec![0.0; nx * ny],
            dst_x: Dst::new(nx),
            dst_y: Dst::new(ny),
            iterations: 0,
        }
    }

    /// Node `(i, j)`'s position (`0..=nx + 1`, `0..=ny + 1`).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn node(&self, i: usize, j: usize) -> (f64, f64) {
        (
            self.origin.0 + i as f64 * self.h,
            self.origin.1 + j as f64 * self.h,
        )
    }

    /// `χ` at node `(i, j)`, the frame's 0 included.
    fn at(&self, f: &[f64], i: usize, j: usize) -> f64 {
        if i == 0 || j == 0 || i > self.nx || j > self.ny {
            0.0
        } else {
            f[(i - 1) + self.nx * (j - 1)]
        }
    }

    /// `L f = −∇²f` on the interior (the frame at 0), the compact five-point stencil.
    fn laplacian(&self, f: &[f64], out: &mut [f64]) {
        let h2 = self.h * self.h;
        for j in 1..=self.ny {
            for i in 1..=self.nx {
                let c = self.at(f, i, j);
                let s = self.at(f, i - 1, j)
                    + self.at(f, i + 1, j)
                    + self.at(f, i, j - 1)
                    + self.at(f, i, j + 1);
                out[(i - 1) + self.nx * (j - 1)] = (4.0 * c - s) / h2;
            }
        }
    }

    /// `A = ∇ × χẑ`: `A_x = ∂_y χ`, `A_y = −∂_x χ` on the edges.
    fn curl(&self, f: &[f64], out: &mut Edges) {
        let h = self.h;
        for j in 0..=self.ny {
            for i in 1..=self.nx {
                out.x[(i - 1) + self.nx * j] = (self.at(f, i, j + 1) - self.at(f, i, j)) / h;
            }
        }
        for j in 1..=self.ny {
            for i in 0..=self.nx {
                out.y[i + (self.nx + 1) * (j - 1)] = -(self.at(f, i + 1, j) - self.at(f, i, j)) / h;
            }
        }
    }

    /// The adjoint of `curl` (`Σ_edges e · curl f = Σ_nodes f · curlᵀ e`): the discrete
    /// `(∇ × e)_z` on the interior nodes.
    fn curl_t(&self, e: &Edges, out: &mut [f64]) {
        let h = self.h;
        let nx = self.nx;
        for j in 1..=self.ny {
            for i in 1..=nx {
                let jx_below = e.x[(i - 1) + nx * (j - 1)];
                let jx_above = e.x[(i - 1) + nx * j];
                let jy_left = e.y[(i - 1) + (nx + 1) * (j - 1)];
                let jy_right = e.y[i + (nx + 1) * (j - 1)];
                out[(i - 1) + nx * (j - 1)] = (jx_below - jx_above) / h + (jy_right - jy_left) / h;
            }
        }
    }

    /// `μ₀ L⁻² f` through sine transforms (the vacuum's inverse: the preconditioner).
    fn vacuum_inverse(&self, f: &[f64], out: &mut [f64]) {
        let (nx, ny) = (self.nx, self.ny);
        let mut work = vec![(0.0, 0.0); 2 * (nx.max(ny) + 1)];
        let mut row = vec![0.0; nx.max(ny)];
        let mut col = vec![0.0; ny];
        let mut spec = vec![0.0; nx * ny];
        // Along x, row by row.
        for j in 0..ny {
            self.dst_x
                .apply(&f[nx * j..nx * (j + 1)], &mut row[..nx], &mut work);
            spec[nx * j..nx * (j + 1)].copy_from_slice(&row[..nx]);
        }
        // Along y, column by column.
        for i in 0..nx {
            for j in 0..ny {
                col[j] = spec[i + nx * j];
            }
            self.dst_y.apply(&col, &mut row[..ny], &mut work);
            for j in 0..ny {
                spec[i + nx * j] = row[j];
            }
        }
        // Divide by L's eigenvalues squared, times μ₀ and the transforms' normalisation.
        let h2 = self.h * self.h;
        #[allow(clippy::cast_precision_loss)]
        let (gx, gy) = ((nx + 1) as f64, (ny + 1) as f64);
        let norm = (2.0 / gx) * (2.0 / gy);
        for j in 0..ny {
            #[allow(clippy::cast_precision_loss)]
            let sy = libm::sin(std::f64::consts::PI * (j + 1) as f64 / (2.0 * gy));
            for i in 0..nx {
                #[allow(clippy::cast_precision_loss)]
                let sx = libm::sin(std::f64::consts::PI * (i + 1) as f64 / (2.0 * gx));
                let lambda = 4.0 * (sx * sx + sy * sy) / h2;
                spec[i + nx * j] *= norm * self.mu0 / (lambda * lambda);
            }
        }
        // Back: along y, then x.
        for i in 0..nx {
            for j in 0..ny {
                col[j] = spec[i + nx * j];
            }
            self.dst_y.apply(&col, &mut row[..ny], &mut work);
            for j in 0..ny {
                spec[i + nx * j] = row[j];
            }
        }
        for j in 0..ny {
            self.dst_x
                .apply(&spec[nx * j..nx * (j + 1)], &mut row[..nx], &mut work);
            out[nx * j..nx * (j + 1)].copy_from_slice(&row[..nx]);
        }
    }

    /// `K f = (1/μ₀) L L f + curlᵀ β curl f`.
    fn operator(&self, f: &[f64], beta: &Edges, out: &mut [f64], scratch: &mut Scratch) {
        self.laplacian(f, &mut scratch.lf);
        self.laplacian(&scratch.lf, out);
        let inv = 1.0 / self.mu0;
        for v in out.iter_mut() {
            *v *= inv;
        }
        self.curl(f, &mut scratch.edges);
        for (a, b) in scratch.edges.x.iter_mut().zip(&beta.x) {
            *a *= b;
        }
        for (a, b) in scratch.edges.y.iter_mut().zip(&beta.y) {
            *a *= b;
        }
        self.curl_t(&scratch.edges, &mut scratch.lf);
        for (o, v) in out.iter_mut().zip(&scratch.lf) {
            *o += v;
        }
    }

    /// Solves for `χ` given the canonical current `j` and the response `beta` on the edges
    /// (both per area), from the last solution, to a relative residual `tol`. Returns the
    /// iterations taken.
    pub fn solve(&mut self, j: &Edges, beta: &Edges, tol: f64) -> usize {
        let n = self.nx * self.ny;
        let mut rhs = vec![0.0; n];
        self.curl_t(j, &mut rhs);
        let mut scratch = Scratch {
            lf: vec![0.0; n],
            edges: Edges::zeros(self),
        };
        let mut x = std::mem::take(&mut self.chi);
        let mut ax = vec![0.0; n];
        self.operator(&x, beta, &mut ax, &mut scratch);
        let mut r: Vec<f64> = rhs.iter().zip(&ax).map(|(b, a)| b - a).collect();
        let norm_b = dot(&rhs, &rhs).sqrt();
        let mut z = vec![0.0; n];
        self.vacuum_inverse(&r, &mut z);
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let mut it = 0;
        while dot(&r, &r).sqrt() > tol * norm_b.max(1e-300) && it < 500 {
            self.operator(&p, beta, &mut ax, &mut scratch);
            let alpha = rz / dot(&p, &ax);
            for k in 0..n {
                x[k] += alpha * p[k];
                r[k] -= alpha * ax[k];
            }
            self.vacuum_inverse(&r, &mut z);
            let rz_new = dot(&r, &z);
            let beta_cg = rz_new / rz;
            rz = rz_new;
            for k in 0..n {
                p[k] = z[k] + beta_cg * p[k];
            }
            it += 1;
        }
        self.chi = x;
        self.iterations = it;
        it
    }

    /// `B_z = L χ` on the interior nodes.
    #[must_use]
    pub fn b_z(&self) -> Vec<f64> {
        let mut b = vec![0.0; self.nx * self.ny];
        self.laplacian(&self.chi, &mut b);
        b
    }

    /// `A_T = ∇ × χẑ` on the edges.
    #[must_use]
    pub fn a(&self) -> Edges {
        let mut a = Edges::zeros(self);
        self.curl(&self.chi, &mut a);
        a
    }

    /// The magnetic energy per unit length, `∫ B_z²/(2μ₀)`.
    #[must_use]
    pub fn magnetic_energy(&self) -> f64 {
        let b = self.b_z();
        let h2 = self.h * self.h;
        b.iter().map(|v| v * v).sum::<f64>() * h2 / (2.0 * self.mu0)
    }

    /// The four points of a staggered sub-grid around `(x, y)` with their bilinear weights
    /// and the weights' derivatives along x and y: the x-edges (`(i, j + ½)`) or the y-edges
    /// (`(i + ½, j)`). Points outside the grid (the frame's edges, where `A` is 0) are left
    /// out.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_possible_wrap
    )]
    fn stencil(&self, x: f64, y: f64, x_edges: bool) -> [(Option<usize>, f64, f64, f64); 4] {
        let (shift_x, shift_y) = if x_edges { (0.0, 0.5) } else { (0.5, 0.0) };
        let fi = (x - self.origin.0) / self.h - shift_x;
        let fj = (y - self.origin.1) / self.h - shift_y;
        let (i0, j0) = (libm::floor(fi), libm::floor(fj));
        let (a, b) = (fi - i0, fj - j0);
        let (i0, j0) = (i0 as i64, j0 as i64);
        let (nx, ny) = (self.nx as i64, self.ny as i64);
        let index = |i: i64, j: i64| -> Option<usize> {
            if x_edges {
                ((1..=nx).contains(&i) && (0..=ny).contains(&j))
                    .then(|| ((i - 1) + nx * j) as usize)
            } else {
                ((0..=nx).contains(&i) && (1..=ny).contains(&j))
                    .then(|| (i + (nx + 1) * (j - 1)) as usize)
            }
        };
        let inv = 1.0 / self.h;
        [
            (
                index(i0, j0),
                (1.0 - a) * (1.0 - b),
                -(1.0 - b) * inv,
                -(1.0 - a) * inv,
            ),
            (index(i0 + 1, j0), a * (1.0 - b), (1.0 - b) * inv, -a * inv),
            (index(i0, j0 + 1), (1.0 - a) * b, -b * inv, (1.0 - a) * inv),
            (index(i0 + 1, j0 + 1), a * b, b * inv, a * inv),
        ]
    }

    /// Adds a particle's vector `(vx, vy)` at `(x, y)` to `e` as a density (per area): the
    /// x-component on the x-edges, the y-component on the y-edges, with the weights
    /// `interpolate` reads back (so that `Σ q u · A(x)` pairs with `∫ J · A`).
    pub fn deposit_vector(&self, e: &mut Edges, x: f64, y: f64, vx: f64, vy: f64) {
        let area = self.h * self.h;
        for (k, w, _, _) in self.stencil(x, y, true) {
            if let Some(k) = k {
                e.x[k] += vx * w / area;
            }
        }
        for (k, w, _, _) in self.stencil(x, y, false) {
            if let Some(k) = k {
                e.y[k] += vy * w / area;
            }
        }
    }

    /// Adds a particle's scalar `s` (its `q²/m`) at `(x, y)` to both edge sets as a density.
    pub fn deposit_scalar(&self, e: &mut Edges, x: f64, y: f64, s: f64) {
        self.deposit_vector(e, x, y, s, s);
    }

    /// The edge field `a` at `(x, y)` (bilinear on each staggered sub-grid).
    #[must_use]
    pub fn interpolate(&self, a: &Edges, x: f64, y: f64) -> (f64, f64) {
        let mut out = (0.0, 0.0);
        for (k, w, _, _) in self.stencil(x, y, true) {
            if let Some(k) = k {
                out.0 += a.x[k] * w;
            }
        }
        for (k, w, _, _) in self.stencil(x, y, false) {
            if let Some(k) = k {
                out.1 += a.y[k] * w;
            }
        }
        out
    }

    /// The gradient of `v · a` at `(x, y)` for a fixed `v`: `(∂_x, ∂_y)` of
    /// `vx a_x + vy a_y`, from the bilinear interpolation's derivatives.
    #[must_use]
    pub fn grad_dot(&self, a: &Edges, x: f64, y: f64, vx: f64, vy: f64) -> (f64, f64) {
        let mut g = (0.0, 0.0);
        for (k, _, dx, dy) in self.stencil(x, y, true) {
            if let Some(k) = k {
                g.0 += vx * a.x[k] * dx;
                g.1 += vx * a.x[k] * dy;
            }
        }
        for (k, _, dx, dy) in self.stencil(x, y, false) {
            if let Some(k) = k {
                g.0 += vy * a.y[k] * dx;
                g.1 += vy * a.y[k] * dy;
            }
        }
        g
    }

    /// `B_z` at `(x, y)` (bilinear on the nodes; the frame's 0 beyond).
    #[must_use]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn b_at(&self, b: &[f64], x: f64, y: f64) -> f64 {
        let fi = (x - self.origin.0) / self.h;
        let fj = (y - self.origin.1) / self.h;
        if fi < 0.0 || fj < 0.0 {
            return 0.0;
        }
        let (i0, j0) = (libm::floor(fi), libm::floor(fj));
        let (a, c) = (fi - i0, fj - j0);
        let (i0, j0) = (i0 as usize, j0 as usize);
        let at = |i: usize, j: usize| self.at(b, i, j);
        (1.0 - a) * (1.0 - c) * at(i0, j0)
            + a * (1.0 - c) * at(i0 + 1, j0)
            + (1.0 - a) * c * at(i0, j0 + 1)
            + a * c * at(i0 + 1, j0 + 1)
    }

    /// The four edges of `stencil` as `(i, j)` with their weights, whether or not they lie in
    /// the grid (`edge_index` tells).
    #[allow(clippy::cast_possible_truncation)]
    fn stencil_cells(&self, x: f64, y: f64, x_edges: bool) -> [(i64, i64, f64); 4] {
        let (shift_x, shift_y) = if x_edges { (0.0, 0.5) } else { (0.5, 0.0) };
        let fi = (x - self.origin.0) / self.h - shift_x;
        let fj = (y - self.origin.1) / self.h - shift_y;
        let (i0, j0) = (libm::floor(fi), libm::floor(fj));
        let (a, b) = (fi - i0, fj - j0);
        let (i0, j0) = (i0 as i64, j0 as i64);
        [
            (i0, j0, (1.0 - a) * (1.0 - b)),
            (i0 + 1, j0, a * (1.0 - b)),
            (i0, j0 + 1, (1.0 - a) * b),
            (i0 + 1, j0 + 1, a * b),
        ]
    }

    /// An edge's index in `Edges::x` (`x_edges`) or `Edges::y`, None for the frame's.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss
    )]
    fn edge_index(&self, i: i64, j: i64, x_edges: bool) -> Option<usize> {
        let (nx, ny) = (self.nx as i64, self.ny as i64);
        if x_edges {
            ((1..=nx).contains(&i) && (0..=ny).contains(&j)).then(|| ((i - 1) + nx * j) as usize)
        } else {
            ((0..=nx).contains(&i) && (1..=ny).contains(&j))
                .then(|| (i + (nx + 1) * (j - 1)) as usize)
        }
    }

    /// The grid's self-inductance at a carrier in a uniform medium of response `beta` (per
    /// area; a metal's carriers, `Σ q²/m`): solved to `tol`, once per medium (`SelfResponse`).
    #[must_use]
    #[allow(clippy::cast_possible_wrap)]
    pub fn self_response(&self, beta: f64, tol: f64) -> SelfResponse {
        let (nx, ny) = (self.nx as i64, self.ny as i64);
        let (x_source, y_source) = ((nx / 2 + 1, ny / 2), (nx / 2, ny / 2 + 1));
        let mut medium = Edges::zeros(self);
        medium.x.fill(beta);
        medium.y.fill(beta);
        let unit = 1.0 / (self.h * self.h);
        let response = |x_edges: bool| {
            let (i, j) = if x_edges { x_source } else { y_source };
            let mut current = Edges::zeros(self);
            let k = self.edge_index(i, j, x_edges).expect("the grid's centre");
            if x_edges {
                current.x[k] = unit;
            } else {
                current.y[k] = unit;
            }
            let mut probe = self.clone();
            probe.chi.fill(0.0);
            probe.solve(&current, &medium, tol);
            probe.a()
        };
        SelfResponse {
            x_source,
            y_source,
            from_x: response(true),
            from_y: response(false),
        }
    }
}

/// A carrier's own vector potential at itself per unit of its canonical current `q u`, in a
/// uniform medium (the other carriers' response screening it): what the solver returns at a
/// carrier from the current it deposits there through the same stencil, its self-inductance on
/// the grid, `A_self = q g u`. A macroparticle's velocity jump (a collision, a wall) would
/// otherwise lose `½ q² Δv·gΔv` to it, a self-energy the many electrons it stands for, jumping
/// independently, do not have (PHYSICS.md §2.12). From two solves at the grid's centre, a unit
/// current on one x-edge and on one y-edge, read at any carrier by translation.
#[derive(Clone, Debug)]
pub struct SelfResponse {
    x_source: (i64, i64),
    y_source: (i64, i64),
    from_x: Edges,
    from_y: Edges,
}

impl SelfResponse {
    /// `g` at `(x, y)`: `g[a][b]` the component a of `A_self` per unit component b of `q u`.
    #[must_use]
    pub fn at(&self, grid: &Darwin2d, x: f64, y: f64) -> [[f64; 2]; 2] {
        let stencils = [
            grid.stencil_cells(x, y, true),
            grid.stencil_cells(x, y, false),
        ];
        let mut g = [[0.0; 2]; 2];
        for (b, (source, response)) in
            [(self.x_source, &self.from_x), (self.y_source, &self.from_y)]
                .into_iter()
                .enumerate()
        {
            for &(il, jl, wl) in &stencils[b] {
                if grid.edge_index(il, jl, b == 0).is_none() {
                    continue;
                }
                // The response to a unit current at (il, jl): the centre's, shifted.
                let (di, dj) = (source.0 - il, source.1 - jl);
                for (a, field) in [(0, &response.x), (1, &response.y)] {
                    for &(ik, jk, wk) in &stencils[a] {
                        if grid.edge_index(ik, jk, a == 0).is_none() {
                            continue;
                        }
                        if let Some(k) = grid.edge_index(ik + di, jk + dj, a == 0) {
                            g[a][b] += wk * wl * field[k];
                        }
                    }
                }
            }
        }
        g
    }
}

/// Work arrays of `Darwin2d::operator`.
struct Scratch {
    lf: Vec<f64>,
    edges: Edges,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;

    /// The sine transform against its definition, and its inverse `(2/(n+1)) DST`.
    #[test]
    fn sine_transform() {
        let n = 15;
        let d = Dst::new(n);
        let x: Vec<f64> = (0..n)
            .map(|k| libm::sin(0.7 * k as f64) + 0.1 * k as f64)
            .collect();
        let mut out = vec![0.0; n];
        let mut work = vec![(0.0, 0.0); 2 * (n + 1)];
        d.apply(&x, &mut out, &mut work);
        for k in 1..=n {
            let direct: f64 = (1..=n)
                .map(|j| {
                    x[j - 1] * libm::sin(std::f64::consts::PI * (j * k) as f64 / (n + 1) as f64)
                })
                .sum();
            assert!((out[k - 1] - direct).abs() < 1e-12, "k {k}");
        }
        let mut back = vec![0.0; n];
        d.apply(&out, &mut back, &mut work);
        for (a, b) in back.iter().zip(&x) {
            assert!((a * 2.0 / (n + 1) as f64 - b).abs() < 1e-12);
        }
    }

    /// A square loop of current: `J = curl ψ₀`, ψ₀ = K on the nodes inside (a sheet current
    /// K around the loop, one edge wide).
    fn loop_current(
        g: &Darwin2d,
        k: f64,
        (i0, j0): (usize, usize),
        (i1, j1): (usize, usize),
    ) -> Edges {
        let mut psi = vec![0.0; g.nx * g.ny];
        for j in j0..=j1 {
            for i in i0..=i1 {
                psi[(i - 1) + g.nx * (j - 1)] = k;
            }
        }
        let mut e = Edges::zeros(g);
        g.curl(&psi, &mut e);
        e
    }

    /// V-D1, V-D2: a loop's field is `μ₀K` inside and none outside, its energy `½ μ₀ A K²`
    /// (its inductance per unit length `μ₀` times its area), to the solver's tolerance:
    /// the staggered grid represents a loop of current exactly.
    #[test]
    fn loop_field_and_inductance() {
        let mut g = Darwin2d::new(63, 63, 0.25, (0.0, 0.0), 10.0);
        let k = 1.5;
        let (lo, hi) = ((20, 24), (44, 40));
        let j = loop_current(&g, k, lo, hi);
        let beta = Edges::zeros(&g);
        g.solve(&j, &beta, 1e-12);
        let b = g.b_z();
        let inside =
            |i: usize, jj: usize| (lo.0..=hi.0).contains(&i) && (lo.1..=hi.1).contains(&jj);
        let mut nodes_in = 0usize;
        for jj in 1..=g.ny {
            for i in 1..=g.nx {
                let v = b[(i - 1) + g.nx * (jj - 1)];
                let want = if inside(i, jj) { g.mu0 * k } else { 0.0 };
                assert!(
                    (v - want).abs() < 1e-9 * g.mu0 * k,
                    "({i}, {jj}): {v} vs {want}"
                );
                if inside(i, jj) {
                    nodes_in += 1;
                }
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let area = nodes_in as f64 * g.h * g.h;
        let energy = g.magnetic_energy();
        let want = 0.5 * g.mu0 * area * k * k;
        assert!((energy - want).abs() < 1e-9 * want, "{energy} vs {want}");
    }

    /// V-D3: nested loops: the inner one's flux from the outer one's current is
    /// `M K₁` with `M = μ₀` times the inner area, and the interaction energy `∫ J₂ · A₁`
    /// equals `M K₁ K₂`.
    #[test]
    fn nested_loops_mutual_inductance() {
        let mut g = Darwin2d::new(63, 63, 0.25, (0.0, 0.0), 7.0);
        let (k1, k2) = (2.0, 0.5);
        let outer = loop_current(&g, k1, (10, 10), (54, 54));
        let inner = loop_current(&g, k2, (25, 28), (36, 40));
        let beta = Edges::zeros(&g);
        g.solve(&outer, &beta, 1e-12);
        let a1 = g.a();
        let h2 = g.h * g.h;
        #[allow(clippy::cast_precision_loss)]
        let inner_area = ((36 - 25 + 1) * (40 - 28 + 1)) as f64 * h2;
        let m = g.mu0 * inner_area;
        let interaction: f64 = (a1.x.iter().zip(&inner.x).map(|(a, j)| a * j).sum::<f64>()
            + a1.y.iter().zip(&inner.y).map(|(a, j)| a * j).sum::<f64>())
            * h2;
        assert!(
            (interaction - m * k1 * k2).abs() < 1e-9 * m * k1 * k2,
            "{interaction} vs {}",
            m * k1 * k2
        );
    }

    /// A particle's deposit pairs with the interpolation: `Σ_edges e·a h²` equals
    /// `v · a(x)` for any edge field `a` (the energy `∫ J·A` and the force use the same
    /// weights), and `grad_dot` is the interpolation's gradient (finite differences).
    #[test]
    fn deposit_pairs_with_interpolation() {
        let g = Darwin2d::new(15, 15, 0.5, (-1.0, -2.0), 5.0);
        let mut a = Edges::zeros(&g);
        for (k, v) in a.x.iter_mut().enumerate() {
            *v = libm::sin(0.37 * k as f64) + 0.2;
        }
        for (k, v) in a.y.iter_mut().enumerate() {
            *v = libm::cos(0.53 * k as f64) - 0.1;
        }
        let h2 = g.h * g.h;
        for &(x, y) in &[(1.3, 0.7), (2.7, 3.1), (0.1, -1.6), (5.9, 4.2)] {
            let (vx, vy) = (0.8, -1.7);
            let mut e = Edges::zeros(&g);
            g.deposit_vector(&mut e, x, y, vx, vy);
            let pair: f64 = (e.x.iter().zip(&a.x).map(|(p, q)| p * q).sum::<f64>()
                + e.y.iter().zip(&a.y).map(|(p, q)| p * q).sum::<f64>())
                * h2;
            let (ax, ay) = g.interpolate(&a, x, y);
            assert!((pair - (vx * ax + vy * ay)).abs() < 1e-12, "({x}, {y})");
            let (gx, gy) = g.grad_dot(&a, x, y, vx, vy);
            let d = 1e-6;
            let f = |x: f64, y: f64| {
                let (ax, ay) = g.interpolate(&a, x, y);
                vx * ax + vy * ay
            };
            let (fx, fy) = (
                (f(x + d, y) - f(x - d, y)) / (2.0 * d),
                (f(x, y + d) - f(x, y - d)) / (2.0 * d),
            );
            assert!(
                (gx - fx).abs() < 1e-6 && (gy - fy).abs() < 1e-6,
                "({x}, {y})"
            );
        }
    }

    /// The modified Bessel function `I₀` by its series (to 1e-15 for `x` up to ~20).
    fn bessel_i0(x: f64) -> f64 {
        let q = 0.25 * x * x;
        let (mut term, mut sum) = (1.0, 1.0);
        for k in 1..60 {
            term *= q / (k * k) as f64;
            sum += term;
        }
        sum
    }

    /// V-D4: a disc of carriers (β, a skin depth `λ = 1/√(μ₀β)`, radius `a`) at the centre of
    /// a round loop (the field `B₀ = μ₀K` outside the disc) screens the field as London's
    /// `B₀ I₀(r/λ)/I₀(a/λ)` (a cylinder in an axial field). Round, so that the screening
    /// current `−βA` circles inside the disc (a slab that ends would need the electrostatic
    /// field to close its current: in the full model it does; here nothing would). The
    /// grid's error: `O(h²/λ²)` and the disc's staircase edge, 2 % allowed at
    /// `h/λ = 1/8`, `a/h = 24`, more than a cell from the edge.
    #[test]
    fn disc_screens_the_field() {
        let mut g = Darwin2d::new(127, 127, 0.25, (0.0, 0.0), 10.0);
        let k = 1.0;
        let (xc, yc) = g.node(64, 64);
        let r_loop = 12.0;
        // The loop: ψ₀ = K on the nodes within r_loop (a round sheet current).
        let mut psi = vec![0.0; g.nx * g.ny];
        for jj in 1..=g.ny {
            for i in 1..=g.nx {
                let (x, y) = g.node(i, jj);
                if libm::hypot(x - xc, y - yc) <= r_loop {
                    psi[(i - 1) + g.nx * (jj - 1)] = k;
                }
            }
        }
        let mut j = Edges::zeros(&g);
        g.curl(&psi, &mut j);
        let (lambda, a) = (2.0, 6.0);
        let beta_value = 1.0 / (g.mu0 * lambda * lambda);
        let mut beta = Edges::zeros(&g);
        for jj in 0..=g.ny {
            for i in 1..=g.nx {
                let (x, y) = g.node(i, jj);
                if libm::hypot(x - xc, y + 0.5 * g.h - yc) <= a {
                    beta.x[(i - 1) + g.nx * jj] = beta_value;
                }
            }
        }
        for jj in 1..=g.ny {
            for i in 0..=g.nx {
                let (x, y) = g.node(i, jj);
                if libm::hypot(x + 0.5 * g.h - xc, y - yc) <= a {
                    beta.y[i + (g.nx + 1) * (jj - 1)] = beta_value;
                }
            }
        }
        let it = g.solve(&j, &beta, 1e-11);
        let b = g.b_z();
        let b0 = g.mu0 * k;
        let mut worst: f64 = 0.0;
        for jj in 1..=g.ny {
            for i in 1..=g.nx {
                let (x, y) = g.node(i, jj);
                let r = libm::hypot(x - xc, y - yc);
                if r <= a - 1.0 {
                    let want = b0 * bessel_i0(r / lambda) / bessel_i0(a / lambda);
                    let got = b[(i - 1) + g.nx * (jj - 1)];
                    worst = worst.max(((got - want) / want).abs());
                }
            }
        }
        println!("V-D4: worst relative error {worst:.2e} ({it} iterations)");
        assert!(
            worst < 0.02,
            "worst relative error {worst} ({it} iterations)"
        );
        // Between the disc and the loop the field is B₀ (the screening current is the
        // disc's own, closed within it).
        let (x, y) = g.node(64, 64 + 36);
        assert!(libm::hypot(x - xc, y - yc) > a + 1.0);
        let between = b[63 + g.nx * (64 + 36 - 1)];
        assert!(
            (between - b0).abs() < 0.02 * b0,
            "between: {between} vs {b0}"
        );
    }

    /// The self-response (`SelfResponse`, the carriers' impulses' correction): a carrier's
    /// current deposited at an off-centre point of a uniform medium, solved directly, and its
    /// own vector potential read back through the same stencil, against the kernel's `g` read
    /// there by translation: within 1e-4. They differ by the frame alone: a current's
    /// transverse part (the Coulomb gauge's) has a long-range tail the medium does not screen,
    /// which the frame 8 away bends differently at the centre and at the point (measured up to
    /// 5.4e-6, unchanged with the solves converged to 1e-14); a wrong index or shift would be
    /// off by order 1, and `g`'s own accuracy is the uniform medium's, ~0.5 % in a metal. And
    /// `g` symmetric (the operator is) and positive on its diagonal.
    #[test]
    fn self_response_matches_a_direct_solve() {
        let grid = Darwin2d::new(63, 63, 0.25, (0.0, 0.0), 3.0);
        let beta = 40.0;
        let kernel = grid.self_response(beta, 1e-12);
        let mut medium = Edges::zeros(&grid);
        medium.x.fill(beta);
        medium.y.fill(beta);
        for &(x, y) in &[(8.13, 7.91), (7.6, 8.37), (8.0, 8.0)] {
            let g = kernel.at(&grid, x, y);
            for (b, (ux, uy)) in [(1.0, 0.0), (0.0, 1.0)].into_iter().enumerate() {
                let mut j = Edges::zeros(&grid);
                grid.deposit_vector(&mut j, x, y, ux, uy);
                let mut direct = grid.clone();
                direct.solve(&j, &medium, 1e-12);
                let (ax, ay) = direct.interpolate(&direct.a(), x, y);
                for (a, value) in [ax, ay].into_iter().enumerate() {
                    let scale = g[0][0].abs().max(g[1][1].abs());
                    assert!(
                        (value - g[a][b]).abs() < 1e-4 * scale,
                        "({x}, {y}) g[{a}][{b}]: direct {value}, kernel {}",
                        g[a][b]
                    );
                }
            }
            assert!(g[0][0] > 0.0 && g[1][1] > 0.0, "{g:?}");
            assert!((g[0][1] - g[1][0]).abs() < 1e-4 * g[0][0], "{g:?}");
        }
    }
}
