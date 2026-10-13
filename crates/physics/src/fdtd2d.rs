//! The z-invariant world's electromagnetic field (PHYSICS.md §2.11): Maxwell's equations
//! for the transverse-electric field (`E_x`, `E_y` in the plane, `B_z` along the axis) on a
//! Yee grid, in Gaussian units,
//! `∂B_z/∂t = −c (∂E_y/∂x − ∂E_x/∂y)`, `∂E/∂t = c ∇×(B_z ẑ) − 4π J`,
//! by the leapfrog of Yee (1966): `E_x` on the edges between nodes `(i, j)` and `(i+1, j)`,
//! `E_y` between `(i, j)` and `(i, j+1)`, `B_z` at the cell centres, half a step apart.
//!
//! Sources are line charges (per unit length) moving in the plane: each move is deposited
//! as a current by Esirkepov's scheme (2001) for linear (cloud-in-cell) shapes, which
//! satisfies the discrete continuity equation exactly, so that Gauss's law
//! `∇·E = 4πρ` holds at every step to rounding if it held at the start. Conductors are
//! perfect (a staircase of nodes): the edges between two metal nodes carry no field; charge
//! carried onto a conductor's surface stays there as its surface charge, redistributed by
//! the field. The grid's border is a perfect conductor, behind an optional absorbing layer
//! (the convolutional PML of Roden and Gedney, 2000) that lets outgoing waves leave.
//!
//! The start: the electrostatic field of given conductor potentials and charges, solved on
//! the same grid (conjugate gradients on the five-point Laplacian), so that it satisfies
//! the discrete Gauss's law and is a stationary solution.

use rayon::prelude::*;

/// The grading of the absorbing layer's conductivity (`σ ∝ depth^m`).
const PML_GRADING: i32 = 3;
/// The absorbing layer's design reflection at normal incidence.
const PML_REFLECTION: f64 = 1e-8;

/// A rectangular grid of `nx × ny` cells of side `h`, its first node at `origin`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub h: f64,
    pub origin: (f64, f64),
}

impl Grid {
    /// The node `(i, j)`'s position.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn node(&self, i: usize, j: usize) -> (f64, f64) {
        (
            self.origin.0 + i as f64 * self.h,
            self.origin.1 + j as f64 * self.h,
        )
    }

    fn node_index(&self, i: usize, j: usize) -> usize {
        i + (self.nx + 1) * j
    }

    /// Nodes in all.
    #[must_use]
    pub fn nodes(&self) -> usize {
        (self.nx + 1) * (self.ny + 1)
    }
}

/// The convolutional PML's coefficients along one axis at every half-cell position
/// (index `2k`: node `k`, `2k + 1`: the cell centre after it).
#[derive(Clone, Debug)]
struct Profile {
    b: Vec<f64>,
    a: Vec<f64>,
}

impl Profile {
    #[allow(clippy::cast_precision_loss)]
    fn new(n: usize, cells: usize, h: f64, c: f64, dt: f64) -> Self {
        let len = 2 * n + 1;
        let (mut b, mut a) = (vec![1.0; len], vec![0.0; len]);
        if cells == 0 {
            return Self { b, a };
        }
        let thick = cells as f64;
        let sigma_max =
            f64::from(PML_GRADING + 1) * c * -libm::log(PML_REFLECTION) / (2.0 * thick * h);
        for (k, (bk, ak)) in b.iter_mut().zip(a.iter_mut()).enumerate() {
            let u = k as f64 / 2.0;
            let depth = (thick - u).max(u - (n as f64 - thick)).max(0.0) / thick;
            if depth <= 0.0 {
                continue;
            }
            let sigma = sigma_max * libm::pow(depth, f64::from(PML_GRADING));
            *bk = libm::exp(-sigma * dt);
            *ak = *bk - 1.0;
        }
        Self { b, a }
    }
}

/// The field on a grid, its conductors and its absorbing layer.
#[derive(Clone, Debug)]
pub struct Fdtd {
    pub grid: Grid,
    /// Speed of light.
    pub c: f64,
    /// Time step.
    pub dt: f64,
    /// Time of `E` (and of the charges); `B_z` lags it by half a step.
    pub t: f64,
    /// `E_x` on the edges `(i, j)–(i+1, j)`: `nx × (ny + 1)`, index `i + nx j`.
    pub ex: Vec<f64>,
    /// `E_y` on the edges `(i, j)–(i, j+1)`: `(nx + 1) × ny`, index `i + (nx + 1) j`.
    pub ey: Vec<f64>,
    /// `B_z` at the cell centres: `nx × ny`, index `i + nx j`.
    pub bz: Vec<f64>,
    /// The current density over the step being taken (as `ex`, `ey`).
    jx: Vec<f64>,
    jy: Vec<f64>,
    /// Nodes inside a conductor.
    metal: Vec<bool>,
    /// The edges that carry no field (both ends in metal, or on the border).
    dead_x: Vec<bool>,
    dead_y: Vec<bool>,
    px: Profile,
    py: Profile,
    /// The PML's memory of `∂E_y/∂x`, `∂E_x/∂y` (for `B_z`), `∂B_z/∂y` (for `E_x`) and
    /// `∂B_z/∂x` (for `E_y`).
    psi_bx: Vec<f64>,
    psi_by: Vec<f64>,
    psi_ex: Vec<f64>,
    psi_ey: Vec<f64>,
}

impl Fdtd {
    /// An empty field on `grid` at time 0, its time step `courant` times the stability
    /// limit `h/(c√2)`, `pml` cells of absorbing layer inside each border (0: none, a
    /// closed box), `metal` the nodes inside conductors.
    ///
    /// # Panics
    /// If `metal` has another length than the nodes.
    #[must_use]
    pub fn new(grid: Grid, c: f64, courant: f64, pml: usize, metal: Vec<bool>) -> Self {
        assert_eq!(metal.len(), grid.nodes());
        let Grid { nx, ny, h, .. } = grid;
        let dt = courant * h / (c * std::f64::consts::SQRT_2);
        let dead_x = (0..nx * (ny + 1))
            .map(|k| {
                let (i, j) = (k % nx, k / nx);
                j == 0
                    || j == ny
                    || (metal[grid.node_index(i, j)] && metal[grid.node_index(i + 1, j)])
            })
            .collect();
        let dead_y = (0..(nx + 1) * ny)
            .map(|k| {
                let (i, j) = (k % (nx + 1), k / (nx + 1));
                i == 0
                    || i == nx
                    || (metal[grid.node_index(i, j)] && metal[grid.node_index(i, j + 1)])
            })
            .collect();
        Self {
            grid,
            c,
            dt,
            t: 0.0,
            ex: vec![0.0; nx * (ny + 1)],
            ey: vec![0.0; (nx + 1) * ny],
            bz: vec![0.0; nx * ny],
            jx: vec![0.0; nx * (ny + 1)],
            jy: vec![0.0; (nx + 1) * ny],
            metal,
            dead_x,
            dead_y,
            px: Profile::new(nx, pml, h, c, dt),
            py: Profile::new(ny, pml, h, c, dt),
            psi_bx: vec![0.0; nx * ny],
            psi_by: vec![0.0; nx * ny],
            psi_ex: vec![0.0; nx * (ny + 1)],
            psi_ey: vec![0.0; (nx + 1) * ny],
        }
    }

    /// Whether node `(i, j)` is inside a conductor.
    #[must_use]
    pub fn is_metal(&self, i: usize, j: usize) -> bool {
        self.metal[self.grid.node_index(i, j)]
    }

    /// The electrostatic potential at the nodes: `potential(i, j)` fixes a node (a
    /// conductor's, or the border's: `Some`), `rho` the charge density at the nodes; the
    /// rest solves `−∇²φ = 4πρ` (five-point Laplacian) by conjugate gradients to
    /// `tolerance` of the right-hand side's norm.
    #[must_use]
    #[allow(clippy::many_single_char_names)]
    pub fn solve_static(
        &self,
        potential: impl Fn(usize, usize) -> Option<f64>,
        rho: &[f64],
        tolerance: f64,
    ) -> Vec<f64> {
        let g = self.grid;
        let (nx, ny, h) = (g.nx, g.ny, g.h);
        let n = g.nodes();
        let mut phi = vec![0.0; n];
        let mut free = vec![false; n];
        for j in 0..=ny {
            for i in 0..=nx {
                let k = g.node_index(i, j);
                let border = i == 0 || j == 0 || i == nx || j == ny;
                match potential(i, j) {
                    Some(v) => phi[k] = v,
                    None if border => {}
                    None => free[k] = true,
                }
            }
        }
        // A x = b on the free nodes: (4φ − Σ neighbours)/h² = 4πρ + fixed neighbours/h².
        let apply = |x: &[f64], out: &mut [f64]| {
            out.par_iter_mut().enumerate().for_each(|(k, o)| {
                if !free[k] {
                    *o = 0.0;
                    return;
                }
                let s = x[k - 1] + x[k + 1] + x[k - (nx + 1)] + x[k + nx + 1];
                *o = (4.0 * x[k] - s) / (h * h);
            });
        };
        let mut b = vec![0.0; n];
        for k in 0..n {
            if !free[k] {
                continue;
            }
            let fixed: f64 = [k - 1, k + 1, k - (nx + 1), k + nx + 1]
                .iter()
                .filter(|&&m| !free[m])
                .map(|&m| phi[m])
                .sum();
            b[k] = 4.0 * std::f64::consts::PI * rho[k] + fixed / (h * h);
        }
        let mut x = vec![0.0; n];
        let mut r = b.clone();
        let mut p = r.clone();
        let mut ap = vec![0.0; n];
        let dot = |u: &[f64], v: &[f64]| -> f64 { u.iter().zip(v).map(|(a, b)| a * b).sum() };
        let norm_b = dot(&b, &b).sqrt().max(1e-300);
        let mut rr = dot(&r, &r);
        for _ in 0..20 * n {
            if rr.sqrt() <= tolerance * norm_b {
                break;
            }
            apply(&p, &mut ap);
            let alpha = rr / dot(&p, &ap);
            for k in 0..n {
                x[k] += alpha * p[k];
                r[k] -= alpha * ap[k];
            }
            let rr_new = dot(&r, &r);
            let beta = rr_new / rr;
            rr = rr_new;
            for k in 0..n {
                p[k] = r[k] + beta * p[k];
            }
        }
        for k in 0..n {
            if free[k] {
                phi[k] = x[k];
            }
        }
        phi
    }

    /// Sets `E = −∇φ` (differences along the edges) and `B_z = 0`: with `φ` from
    /// `solve_static`, a stationary field obeying the discrete Gauss's law.
    pub fn set_static(&mut self, phi: &[f64]) {
        let g = self.grid;
        let (nx, ny, h) = (g.nx, g.ny, g.h);
        for j in 0..=ny {
            for i in 0..nx {
                let k = i + nx * j;
                self.ex[k] = if self.dead_x[k] {
                    0.0
                } else {
                    -(phi[g.node_index(i + 1, j)] - phi[g.node_index(i, j)]) / h
                };
            }
        }
        for j in 0..ny {
            for i in 0..=nx {
                let k = i + (nx + 1) * j;
                self.ey[k] = if self.dead_y[k] {
                    0.0
                } else {
                    -(phi[g.node_index(i, j + 1)] - phi[g.node_index(i, j)]) / h
                };
            }
        }
        self.bz.iter_mut().for_each(|b| *b = 0.0);
    }

    /// The charge density at the nodes of a line charge `q` at `(x, y)` (its share by
    /// linear weights over the cell: `q S / h²`), added to `rho`.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn add_charge(&self, rho: &mut [f64], q: f64, x: f64, y: f64) {
        let g = self.grid;
        let (u, v) = ((x - g.origin.0) / g.h, (y - g.origin.1) / g.h);
        let (i0, j0) = (libm::floor(u), libm::floor(v));
        let (fu, fv) = (u - i0, v - j0);
        let (i0, j0) = (i0 as usize, j0 as usize);
        let w = q / (g.h * g.h);
        for (di, wx) in [(0, 1.0 - fu), (1, fu)] {
            for (dj, wy) in [(0, 1.0 - fv), (1, fv)] {
                rho[g.node_index(i0 + di, j0 + dj)] += w * wx * wy;
            }
        }
    }

    /// Deposits the move of a line charge `q` from `a` to `b` during the step being taken
    /// (Esirkepov's scheme, in pieces of at most half a cell): the current that carries its
    /// charge density from the one place to the other, exactly. Both ends must lie at
    /// least a cell inside the grid.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn deposit(&mut self, q: f64, a: (f64, f64), b: (f64, f64)) {
        let g = self.grid;
        let (du, dv) = ((b.0 - a.0) / g.h, (b.1 - a.1) / g.h);
        let pieces = (2.0 * du.abs().max(dv.abs())).ceil().max(1.0);
        let n = pieces as usize;
        let (u0, v0) = ((a.0 - g.origin.0) / g.h, (a.1 - g.origin.1) / g.h);
        for k in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let (s0, s1) = (k as f64 / pieces, (k + 1) as f64 / pieces);
            self.deposit_piece(
                q,
                (u0 + du * s0, v0 + dv * s0),
                (u0 + du * s1, v0 + dv * s1),
            );
        }
    }

    /// One piece (grid units, shorter than a cell): Esirkepov's `W` and the currents
    /// summed from the window's low side.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    fn deposit_piece(&mut self, q: f64, a: (f64, f64), b: (f64, f64)) {
        let g = self.grid;
        let (nx, h, dt) = (g.nx, g.h, self.dt);
        let bi = libm::floor(a.0.min(b.0)) as usize;
        let bj = libm::floor(a.1.min(b.1)) as usize;
        let shape = |u: f64, base: usize| -> [f64; 4] {
            let mut s = [0.0; 4];
            for (k, sk) in s.iter_mut().enumerate() {
                *sk = (1.0 - (u - (base + k) as f64).abs()).max(0.0);
            }
            s
        };
        let (sx0, sx1) = (shape(a.0, bi), shape(b.0, bi));
        let (sy0, sy1) = (shape(a.1, bj), shape(b.1, bj));
        let f = -q / (h * dt);
        for l in 0..4 {
            let mut acc = 0.0;
            for k in 0..3 {
                let wx = (sx1[k] - sx0[k]) * (sy0[l] + 0.5 * (sy1[l] - sy0[l]));
                acc += f * wx;
                let (i, j) = (bi + k, bj + l);
                if acc != 0.0 && i < nx && j <= g.ny {
                    self.jx[i + nx * j] += acc;
                }
            }
        }
        for k in 0..4 {
            let mut acc = 0.0;
            for l in 0..3 {
                let wy = (sy1[l] - sy0[l]) * (sx0[k] + 0.5 * (sx1[k] - sx0[k]));
                acc += f * wy;
                let (i, j) = (bi + k, bj + l);
                if acc != 0.0 && i <= nx && j < g.ny {
                    self.jy[i + (nx + 1) * j] += acc;
                }
            }
        }
    }

    /// Advances the field by one step: `B_z` by half a step past `E`'s time, then `E` by
    /// the curl of `B_z` and the currents deposited since the last step (cleared after).
    #[allow(clippy::many_single_char_names)]
    pub fn step(&mut self) {
        let Grid { nx, h, .. } = self.grid;
        let (cdt, dt) = (self.c * self.dt / h, self.dt);
        // B_z.
        {
            let (ex, ey) = (&self.ex, &self.ey);
            let (px, py) = (&self.px, &self.py);
            self.bz
                .par_chunks_mut(nx)
                .zip(self.psi_bx.par_chunks_mut(nx))
                .zip(self.psi_by.par_chunks_mut(nx))
                .enumerate()
                .for_each(|(j, ((bz, psx), psy))| {
                    for i in 0..nx {
                        let dx = ey[i + 1 + (nx + 1) * j] - ey[i + (nx + 1) * j];
                        let dy = ex[i + nx * (j + 1)] - ex[i + nx * j];
                        let (bxk, byk) = (2 * i + 1, 2 * j + 1);
                        psx[i] = px.b[bxk] * psx[i] + px.a[bxk] * dx;
                        psy[i] = py.b[byk] * psy[i] + py.a[byk] * dy;
                        bz[i] -= cdt * ((dx + psx[i]) - (dy + psy[i]));
                    }
                });
        }
        // E_x.
        {
            let bz = &self.bz;
            let py = &self.py;
            let four_pi = 4.0 * std::f64::consts::PI;
            self.ex
                .par_chunks_mut(nx)
                .zip(self.psi_ex.par_chunks_mut(nx))
                .zip(self.jx.par_chunks_mut(nx))
                .zip(self.dead_x.par_chunks(nx))
                .enumerate()
                .for_each(|(j, (((ex, ps), jx), dead))| {
                    for i in 0..nx {
                        if dead[i] {
                            ex[i] = 0.0;
                        } else {
                            let dy = bz[i + nx * j] - bz[i + nx * (j - 1)];
                            ps[i] = py.b[2 * j] * ps[i] + py.a[2 * j] * dy;
                            ex[i] += cdt * (dy + ps[i]) - four_pi * dt * jx[i];
                        }
                        jx[i] = 0.0;
                    }
                });
        }
        // E_y.
        {
            let bz = &self.bz;
            let px = &self.px;
            let four_pi = 4.0 * std::f64::consts::PI;
            self.ey
                .par_chunks_mut(nx + 1)
                .zip(self.psi_ey.par_chunks_mut(nx + 1))
                .zip(self.jy.par_chunks_mut(nx + 1))
                .zip(self.dead_y.par_chunks(nx + 1))
                .enumerate()
                .for_each(|(j, (((ey, ps), jy), dead))| {
                    for i in 0..=nx {
                        if dead[i] {
                            ey[i] = 0.0;
                        } else {
                            let dx = bz[i + nx * j] - bz[i - 1 + nx * j];
                            ps[i] = px.b[2 * i] * ps[i] + px.a[2 * i] * dx;
                            ey[i] -= cdt * (dx + ps[i]) + four_pi * dt * jy[i];
                        }
                        jy[i] = 0.0;
                    }
                });
        }
        self.t += self.dt;
    }

    /// `∇·E / 4π` at node `(i, j)` (an interior node): the charge density Gauss's law
    /// gives.
    #[must_use]
    pub fn gauss(&self, i: usize, j: usize) -> f64 {
        let Grid { nx, h, .. } = self.grid;
        let div = (self.ex[i + nx * j] - self.ex[i - 1 + nx * j]) / h
            + (self.ey[i + (nx + 1) * j] - self.ey[i + (nx + 1) * (j - 1)]) / h;
        div / (4.0 * std::f64::consts::PI)
    }

    /// The field energy per unit length, `Σ (E² + B²) h²/8π` (`B_z` half a step behind).
    #[must_use]
    pub fn energy(&self) -> f64 {
        let h2 = self.grid.h * self.grid.h;
        let s: f64 = self
            .ex
            .iter()
            .chain(&self.ey)
            .chain(&self.bz)
            .map(|v| v * v)
            .sum();
        s * h2 / (8.0 * std::f64::consts::PI)
    }

    /// The power per unit length flowing out of the rectangle `x0..x1 × y0..y1`:
    /// `∮ S·n ds`, `S = (c/4π) E × B_z ẑ = (c/4π)(E_y B_z, −E_x B_z)`, the fields sampled
    /// every `h` along its sides (`B_z` half a step behind `E`).
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub fn outflow(&self, x0: f64, y0: f64, x1: f64, y1: f64) -> f64 {
        let h = self.grid.h;
        let k = self.c / (4.0 * std::f64::consts::PI);
        let s = |x: f64, y: f64| {
            let b = self.bz_at(x, y);
            let (ex, ey) = self.e_at(x, y);
            (k * ey * b, -k * ex * b)
        };
        let nx = libm::round((x1 - x0) / h) as usize;
        let ny = libm::round((y1 - y0) / h) as usize;
        let mut p = 0.0;
        for i in 0..nx {
            let x = x0 + (i as f64 + 0.5) * h;
            p += (s(x, y1).1 - s(x, y0).1) * h;
        }
        for j in 0..ny {
            let y = y0 + (j as f64 + 0.5) * h;
            p += (s(x1, y).0 - s(x0, y).0) * h;
        }
        p
    }

    /// `B_z` at a point (bilinear between the cell centres).
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn bz_at(&self, x: f64, y: f64) -> f64 {
        let g = self.grid;
        let u = ((x - g.origin.0) / g.h - 0.5).clamp(0.0, (g.nx - 1) as f64 - 1e-9);
        let v = ((y - g.origin.1) / g.h - 0.5).clamp(0.0, (g.ny - 1) as f64 - 1e-9);
        let (i, j) = (libm::floor(u) as usize, libm::floor(v) as usize);
        let (fu, fv) = (u - i as f64, v - j as f64);
        let at = |i: usize, j: usize| self.bz[(i).min(g.nx - 1) + g.nx * (j).min(g.ny - 1)];
        (1.0 - fu) * (1.0 - fv) * at(i, j)
            + fu * (1.0 - fv) * at(i + 1, j)
            + (1.0 - fu) * fv * at(i, j + 1)
            + fu * fv * at(i + 1, j + 1)
    }

    /// `E` at a point: each component bilinear between its edges.
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn e_at(&self, x: f64, y: f64) -> (f64, f64) {
        let g = self.grid;
        let (u, v) = ((x - g.origin.0) / g.h, (y - g.origin.1) / g.h);
        let lerp = |f: &[f64], w: usize, u: f64, v: f64, ni: usize, nj: usize| {
            let u = u.clamp(0.0, ni as f64 - 1e-9);
            let v = v.clamp(0.0, nj as f64 - 1e-9);
            let (i, j) = (libm::floor(u) as usize, libm::floor(v) as usize);
            let (fu, fv) = (u - i as f64, v - j as f64);
            let at = |i: usize, j: usize| f[i.min(ni) + w * j.min(nj)];
            (1.0 - fu) * (1.0 - fv) * at(i, j)
                + fu * (1.0 - fv) * at(i + 1, j)
                + (1.0 - fu) * fv * at(i, j + 1)
                + fu * fv * at(i + 1, j + 1)
        };
        (
            lerp(&self.ex, g.nx, u - 0.5, v, g.nx - 1, g.ny),
            lerp(&self.ey, g.nx + 1, u, v - 0.5, g.nx, g.ny - 1),
        )
    }
}
