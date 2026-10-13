//! The electromagnetic field of a tube run (PHYSICS.md §2.11; for display):
//! `physics::fdtd2d` driven by the run's frames. The run's model is quasi-static (its
//! carriers feel the instantaneous field); this field is what its charges and currents
//! radiate at a finite speed of light, one-way coupled as the 3D levels' radiation views
//! are.
//!
//! Sources: every particle's motion between frames (matched by identity; a new particle
//! comes out of its emitter's surface, an absorbed one goes into the electrode it reached),
//! and each electrode's circuit current along a straight lead (the owner's convention,
//! 2026-10-12): from the innermost grounded ring enclosing it (an envelope), or, in an open
//! level (one with an output), from the nearest envelope outside it, or else from the
//! board's border; the shortest such segment that crosses no other electrode. A lead
//! carries the electrode's change of charge (from the frames' surface densities) less the
//! particles it took in and gave off, plus what the leads starting on it drew. The
//! electrodes are perfect conductors on the grid, so each one's charge follows the run's.
//! An open level's outermost envelope has no lead: its world conserves its total charge,
//! as a z-invariant world without currents from infinity must (the quasi-static model's
//! open 2D boundary lets the total drift; that drift does not radiate).

use physics::DVec3;
use physics::fdtd2d::{Fdtd, Grid};
use std::collections::HashMap;

use crate::tube::{KIND_NEGATIVE, KIND_POSITIVE, KIND_PROJECTILE, TubeRun};
use crate::{Level, Outline};

/// Absorbing-layer cells of an open tube (one with an output).
pub const EM_PML: usize = 16;
/// Grid cells per board cell.
pub const EM_PER_CELL: f64 = 4.0;
/// The display's speed of light of a tube without an output, in units of its carriers'
/// largest speed (`√(2 |q/m| ΔV)`): fast enough that the field is nearly quasi-static.
pub const EM_C_FACTOR: f64 = 10.0;
/// The time step's fraction of the stability limit.
const COURANT: f64 = 0.9;

/// A lead: from its parent's metal (or the border) to its electrode's.
pub type Lead = ((f64, f64), (f64, f64));

/// A charge moved from one point to another within a step.
type Move = (f64, (f64, f64), (f64, f64));

/// A path a charge takes over a frame interval (from `a` at its start to `b` at its end).
#[derive(Clone, Copy, Debug)]
struct Path {
    q: f64,
    a: (f64, f64),
    b: (f64, f64),
}

/// The sources over one frame interval.
#[derive(Clone, Debug)]
struct Interval {
    t0: f64,
    t1: f64,
    paths: Vec<Path>,
    /// Per electrode: the charge its lead carries over the interval.
    leads: Vec<f64>,
}

/// The field of a tube run, advanced along it.
#[derive(Clone, Debug)]
pub struct TubeEm {
    pub fdtd: Fdtd,
    outlines: Vec<Outline>,
    /// Per electrode: its segments (indices into the run's).
    segments: Vec<Vec<usize>>,
    /// Each grid node's electrode (None: vacuum or the border).
    node_owner: Vec<Option<usize>>,
    /// Per electrode: the electrode its lead starts on (None: the border, or no lead).
    parent: Vec<Option<usize>>,
    /// Per electrode: its lead (from the parent or the border, to it), if it has one.
    pub leads: Vec<Option<Lead>>,
    /// The electrodes children first (each lead's charge includes its children's).
    order: Vec<usize>,
    cathode: usize,
    ions: Option<usize>,
    negative_carriers: bool,
    /// Whether the grid is open (an absorbing layer inside its border).
    open: bool,
    /// The next frame interval to feed: frames `next`..`next + 1`.
    next: usize,
    interval: Option<Interval>,
    /// The energy per unit length that has flowed out through the contour (a rectangle
    /// three cells inside the absorbing layer, or the border), and its history `(t, W)`.
    pub radiated: f64,
    history: Vec<(f64, f64)>,
}

fn p2(v: DVec3) -> (f64, f64) {
    (v.x, v.y)
}

fn v3(p: (f64, f64)) -> DVec3 {
    DVec3::new(p.0, p.1, 0.0)
}

/// The nearest point of segment `a`–`b` to `p`.
fn nearest_on(a: DVec3, b: DVec3, p: DVec3) -> DVec3 {
    let d = b - a;
    let s = ((p - a).dot(d) / d.length_squared().max(1e-300)).clamp(0.0, 1.0);
    a + d * s
}

/// An outline's centre (a box's middle).
fn centre(o: &Outline) -> DVec3 {
    match *o {
        Outline::Rect(r) => (r[0] + r[1] + r[2] + r[3]) * 0.25,
        Outline::Circle { center, .. } | Outline::Ring { center, .. } => center,
    }
}

/// The farthest distance from `c` to an outline's material.
fn farthest(o: &Outline, c: DVec3) -> f64 {
    match *o {
        Outline::Rect(r) => r.iter().map(|p| (*p - c).length()).fold(0.0, f64::max),
        Outline::Circle { center, radius } => (center - c).length() + radius,
        Outline::Ring { center, outer, .. } => (center - c).length() + outer,
    }
}

/// The nearest grid node to `p` of electrode `owner` (None: of the border's two rings),
/// searched in growing squares (`p` itself if none is found).
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap
)]
fn nearest_node(
    grid: &Grid,
    node_owner: &[Option<usize>],
    owner: Option<usize>,
    p: (f64, f64),
    border: &dyn Fn(usize, usize) -> bool,
) -> (f64, f64) {
    let (u, v) = (
        (p.0 - grid.origin.0) / grid.h,
        (p.1 - grid.origin.1) / grid.h,
    );
    let (ci, cj) = (libm::round(u) as i64, libm::round(v) as i64);
    for r in 0..64_i64 {
        let mut best: Option<(f64, (f64, f64))> = None;
        for dj in -r..=r {
            for di in -r..=r {
                if di.abs() != r && dj.abs() != r {
                    continue;
                }
                let (i, j) = (ci + di, cj + dj);
                if i < 0 || j < 0 || i > grid.nx as i64 || j > grid.ny as i64 {
                    continue;
                }
                let (i, j) = (i as usize, j as usize);
                let k = i + (grid.nx + 1) * j;
                let ok = match owner {
                    Some(e) => node_owner[k] == Some(e),
                    None => border(i, j),
                };
                if ok {
                    let q = grid.node(i, j);
                    let d = (q.0 - p.0).powi(2) + (q.1 - p.1).powi(2);
                    if best.is_none_or(|(b, _)| d < b) {
                        best = Some((d, q));
                    }
                }
            }
        }
        if let Some((_, q)) = best {
            return q;
        }
    }
    p
}

impl TubeEm {
    /// The field of `level`'s tube run `run` (its first frame at least), at the display's
    /// speed of light (the output's, or `EM_C_FACTOR` times the carriers' largest speed),
    /// on a grid of `EM_PER_CELL` cells per board cell over the board, open (an absorbing
    /// layer inside the border) when the level has an output, otherwise a closed box.
    ///
    /// # Panics
    /// If the level is not a tube level or the run has no frame.
    #[must_use]
    #[allow(
        clippy::too_many_lines,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn new(level: &Level, run: &TubeRun) -> Self {
        let spec = level.tube.as_ref().expect("a tube level");
        let first = run.frames.first().expect("a frame");
        let (lo, hi) = run.potential_range;
        let speed = (2.0 * spec.charge_per_mass.abs() * (hi - lo).abs())
            .sqrt()
            .max(1.0);
        let c = spec.output.map_or(EM_C_FACTOR * speed, |o| o.c);
        let open = spec.output.is_some();
        let b = level.bounds();
        let h = 1.0 / EM_PER_CELL;
        let nx = ((b.max.x - b.min.x) / h).ceil() as usize;
        let ny = ((b.max.y - b.min.y) / h).ceil() as usize;
        let grid = Grid {
            nx,
            ny,
            h,
            origin: (b.min.x, b.min.y),
        };
        let outlines = level.electrode_outlines(&[]);
        let n_el = outlines.len();
        let mut segments = vec![Vec::new(); n_el];
        for (s, &o) in run.owners.iter().enumerate() {
            if o < n_el {
                segments[o].push(s);
            }
        }
        // The nodes: inside an electrode, or on the border (two rings: a lead ends there).
        let owner_of = |p: DVec3| outlines.iter().position(|o| o.distance_to(p) == 0.0);
        let node_owner: Vec<Option<usize>> = (0..grid.nodes())
            .map(|k| {
                let (i, j) = (k % (nx + 1), k / (nx + 1));
                owner_of(v3(grid.node(i, j)))
            })
            .collect();
        let border = |i: usize, j: usize| i <= 1 || j <= 1 || i + 1 >= nx || j + 1 >= ny;
        let metal: Vec<bool> = (0..grid.nodes())
            .map(|k| {
                let (i, j) = (k % (nx + 1), k / (nx + 1));
                node_owner[k].is_some() || border(i, j)
            })
            .collect();
        let pml = if open { EM_PML } else { 0 };
        let mut fdtd = Fdtd::new(grid, c, COURANT, pml, metal);
        // Envelopes: rings enclosing another electrode.
        let encloses = |p: usize, e: usize| match outlines[p] {
            Outline::Ring { center, inner, .. } => p != e && farthest(&outlines[e], center) < inner,
            _ => false,
        };
        let envelope: Vec<bool> = (0..n_el)
            .map(|p| (0..n_el).any(|e| encloses(p, e)))
            .collect();
        let parent: Vec<Option<usize>> = (0..n_el)
            .map(|e| {
                // The innermost enclosing ring.
                let inner = (0..n_el).filter(|&p| encloses(p, e)).min_by(|&p, &q| {
                    let r = |k: usize| match outlines[k] {
                        Outline::Ring { inner, .. } => inner,
                        _ => f64::INFINITY,
                    };
                    r(p).total_cmp(&r(q))
                });
                inner.or_else(|| {
                    if open && !envelope[e] {
                        (0..n_el)
                            .filter(|&p| envelope[p] && p != e)
                            .min_by(|&p, &q| {
                                outlines[p]
                                    .distance(&outlines[e])
                                    .total_cmp(&outlines[q].distance(&outlines[e]))
                            })
                    } else {
                        None
                    }
                })
            })
            .collect();
        // Leads: the shortest straight path from the electrode's surface to its parent's
        // material (or the border) crossing no other electrode, over 144 directions.
        let inside_grid = |p: DVec3| {
            let (x0, y0) = grid.origin;
            let (x1, y1) = (x0 + nx as f64 * h, y0 + ny as f64 * h);
            p.x > x0 + 1.5 * h && p.y > y0 + 1.5 * h && p.x < x1 - 1.5 * h && p.y < y1 - 1.5 * h
        };
        let leads: Vec<Option<Lead>> = (0..n_el)
            .map(|e| {
                if parent[e].is_none() && open && envelope[e] {
                    return None;
                }
                let c0 = centre(&outlines[e]);
                let reach = farthest(&outlines[e], c0) + 1.0;
                let ring_outer = match outlines[e] {
                    Outline::Ring { outer, .. } => Some(outer),
                    _ => None,
                };
                let mut best: Option<(f64, DVec3, DVec3)> = None;
                for k in 0..144 {
                    let a = f64::from(k) * std::f64::consts::TAU / 144.0;
                    let dir = DVec3::new(libm::cos(a), libm::sin(a), 0.0);
                    // The surface along this direction (a ring's outer one).
                    let s = if let Some(outer) = ring_outer {
                        c0 + dir * outer
                    } else {
                        if outlines[e].distance_to(c0) > 0.0 {
                            continue;
                        }
                        let (mut lo, mut hi) = (0.0, reach);
                        for _ in 0..50 {
                            let mid = 0.5 * (lo + hi);
                            if outlines[e].distance_to(c0 + dir * mid) > 0.0 {
                                hi = mid;
                            } else {
                                lo = mid;
                            }
                        }
                        c0 + dir * hi
                    };
                    // March to the parent or the border.
                    let step = 0.25 * h;
                    let mut d = 0.0;
                    let end = loop {
                        d += step;
                        let p = s + dir * d;
                        match parent[e] {
                            Some(q) if outlines[q].distance_to(p) == 0.0 => break Some(p),
                            None if !inside_grid(p) => break Some(p),
                            _ => {}
                        }
                        if !inside_grid(p)
                            || (0..n_el).any(|o| o != e && outlines[o].distance_to(p) == 0.0)
                        {
                            break None;
                        }
                    };
                    let Some(end) = end else { continue };
                    let len = (end - s).length();
                    if best.is_none_or(|(l, _, _)| len < l) {
                        best = Some((len, s, end));
                    }
                }
                // Each end on a node of its metal (snapped in `new`'s caller below).
                best.map(|(_, s, end)| (p2(end), p2(s)))
            })
            .collect();
        // Children before parents: by depth.
        let depth = |mut e: usize| {
            let mut d = 0;
            while let Some(p) = parent[e] {
                d += 1;
                e = p;
                if d > n_el {
                    break;
                }
            }
            d
        };
        let mut order: Vec<usize> = (0..n_el).collect();
        order.sort_by_key(|&e| std::cmp::Reverse(depth(e)));
        // Each lead's ends on nodes of their metal: all of the charge it carries leaves and
        // arrives there (a point near a thin vane would spread some onto vacuum nodes,
        // where it would stay as a spurious charge).
        let snap = |owner: Option<usize>, p: (f64, f64)| {
            nearest_node(&grid, &node_owner, owner, p, &border)
        };
        let leads: Vec<Option<Lead>> = leads
            .iter()
            .enumerate()
            .map(|(e, l)| l.map(|(from, to)| (snap(parent[e], from), snap(Some(e), to))))
            .collect();
        // The start: the electrostatic field of the first frame's potentials and charges.
        let mut rho = vec![0.0; grid.nodes()];
        for (x, q) in first.x.iter().zip(&first.q) {
            fdtd.add_charge(&mut rho, f64::from(*q), f64::from(x[0]), f64::from(x[1]));
        }
        let phi = fdtd.solve_static(
            |i, j| {
                let k = i + (nx + 1) * j;
                match node_owner[k] {
                    Some(e) => Some(f64::from(first.potentials.get(e).copied().unwrap_or(0.0))),
                    None if border(i, j) => Some(0.0),
                    None => None,
                }
            },
            &rho,
            1e-10,
        );
        fdtd.set_static(&phi);
        fdtd.t = first.t;
        Self {
            fdtd,
            outlines,
            segments,
            node_owner,
            parent,
            leads,
            order,
            cathode: spec.cathode,
            ions: spec.ions.map(|i| i.electrode),
            negative_carriers: spec.charge_per_mass < 0.0,
            open,
            next: 0,
            interval: None,
            radiated: 0.0,
            history: Vec::new(),
        }
    }

    /// The power per unit length flowing out through the contour now (`Fdtd::outflow`
    /// over a rectangle three cells inside the absorbing layer, or the border).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn outflow(&self) -> f64 {
        let g = self.fdtd.grid;
        let inset = (if self.open { EM_PML } else { 0 } + 3) as f64 * g.h;
        self.fdtd.outflow(
            g.origin.0 + inset,
            g.origin.1 + inset,
            g.origin.0 + g.nx as f64 * g.h - inset,
            g.origin.1 + g.ny as f64 * g.h - inset,
        )
    }

    /// The mean outflow over the last `span` of time (0 before there is any).
    #[must_use]
    pub fn mean_outflow(&self, span: f64) -> f64 {
        let Some(&(t1, w1)) = self.history.last() else {
            return 0.0;
        };
        let k = self.history.partition_point(|(t, _)| *t < t1 - span);
        let (t0, w0) = self.history[k.min(self.history.len() - 1)];
        if t1 > t0 { (w1 - w0) / (t1 - t0) } else { 0.0 }
    }

    /// Each electrode's charge per unit length in a frame (its segments' densities).
    fn charges(&self, run: &TubeRun, k: usize) -> Vec<f64> {
        let f = &run.frames[k];
        self.segments
            .iter()
            .map(|segs| {
                segs.iter()
                    .map(|&s| {
                        let (a, b) = run.segments[s];
                        f64::from(f.sigma.get(s).copied().unwrap_or(0.0)) * (b - a).length()
                    })
                    .sum()
            })
            .collect()
    }

    /// The node of electrode `e` nearest to its surface's nearest point to `p`: where a
    /// charge leaving or entering it at `p` comes from or goes to.
    fn surface_node(&self, run: &TubeRun, e: usize, p: DVec3) -> (f64, f64) {
        let s = self.segments[e]
            .iter()
            .map(|&s| {
                let (a, b) = run.segments[s];
                nearest_on(a, b, p)
            })
            .min_by(|a, b| (*a - p).length().total_cmp(&(*b - p).length()))
            .unwrap_or(p);
        let g = self.fdtd.grid;
        let (nx, ny) = (g.nx, g.ny);
        let border = |i: usize, j: usize| i <= 1 || j <= 1 || i + 1 >= nx || j + 1 >= ny;
        nearest_node(&g, &self.node_owner, Some(e), p2(s), &border)
    }

    /// The sources over frames `k`..`k + 1`.
    #[allow(clippy::too_many_lines)]
    fn interval(&self, run: &TubeRun, k: usize) -> Interval {
        let (f0, f1) = (&run.frames[k], &run.frames[k + 1]);
        let n_el = self.outlines.len();
        let next: HashMap<u32, usize> = f1.id.iter().enumerate().map(|(i, id)| (*id, i)).collect();
        let prev: HashMap<u32, usize> = f0.id.iter().enumerate().map(|(i, id)| (*id, i)).collect();
        let mut exchange = vec![0.0; n_el];
        let mut paths = Vec::new();
        let at = |x: [f32; 2]| DVec3::new(f64::from(x[0]), f64::from(x[1]), 0.0);
        for (i, id) in f0.id.iter().enumerate() {
            let (x, q) = (at(f0.x[i]), f64::from(f0.q[i]));
            match next.get(id) {
                Some(&j) => paths.push(Path {
                    q,
                    a: p2(x),
                    b: p2(at(f1.x[j])),
                }),
                None => {
                    // Absorbed by the nearest electrode (lost ones too: none leave a tube
                    // closed by an envelope; elsewhere their charge goes to the nearest).
                    let e = (0..n_el)
                        .min_by(|&a, &b| {
                            self.outlines[a]
                                .distance_to(x)
                                .total_cmp(&self.outlines[b].distance_to(x))
                        })
                        .unwrap_or(0);
                    exchange[e] += q;
                    paths.push(Path {
                        q,
                        a: p2(x),
                        b: self.surface_node(run, e, x),
                    });
                }
            }
        }
        for (j, id) in f1.id.iter().enumerate() {
            if prev.contains_key(id) {
                continue;
            }
            let (x, q, kind) = (at(f1.x[j]), f64::from(f1.q[j]), f1.kind[j]);
            let emitter = match kind {
                KIND_PROJECTILE => None,
                KIND_NEGATIVE if self.negative_carriers => Some(self.cathode),
                KIND_POSITIVE if !self.negative_carriers => Some(self.cathode),
                _ => self.ions,
            };
            let Some(e) = emitter else { continue };
            exchange[e] -= q;
            paths.push(Path {
                q,
                a: self.surface_node(run, e, x),
                b: p2(x),
            });
        }
        let (q0, q1) = (self.charges(run, k), self.charges(run, k + 1));
        let mut leads = vec![0.0; n_el];
        for &e in &self.order {
            let children: f64 = (0..n_el)
                .filter(|&c| self.parent[c] == Some(e) && self.leads[c].is_some())
                .map(|c| leads[c])
                .sum();
            leads[e] = q1[e] - q0[e] - exchange[e] + children;
        }
        Interval {
            t0: f0.t,
            t1: f1.t,
            paths,
            leads,
        }
    }

    /// The charge density at the grid's nodes the field's time implies: every path's
    /// charge where it is now (Gauss's law must give it in vacuum).
    #[must_use]
    pub fn charge_density(&self) -> Vec<f64> {
        let mut rho = vec![0.0; self.fdtd.grid.nodes()];
        if let Some(iv) = &self.interval {
            let g = ((self.fdtd.t - iv.t0) / (iv.t1 - iv.t0).max(1e-300)).clamp(0.0, 1.0);
            for p in &iv.paths {
                let x = (p.a.0 + (p.b.0 - p.a.0) * g, p.a.1 + (p.b.1 - p.a.1) * g);
                self.fdtd.add_charge(&mut rho, p.q, x.0, x.1);
            }
        }
        rho
    }

    /// Whether node `(i, j)` and its neighbours are all vacuum (Gauss's law applies there
    /// alone: metal nodes hold the conductors' charge).
    #[must_use]
    pub fn in_vacuum(&self, i: usize, j: usize) -> bool {
        (i.saturating_sub(1)..=i + 1)
            .all(|a| (j.saturating_sub(1)..=j + 1).all(|b| !self.fdtd.is_metal(a, b)))
    }

    /// Advances the field to time `t` (as far as the run's frames reach).
    pub fn advance_to(&mut self, run: &TubeRun, t: f64) {
        let dt = self.fdtd.dt;
        let end = t.min(run.frames.last().map_or(0.0, |f| f.t));
        while self.fdtd.t + dt <= end {
            let (mut s0, s1) = (self.fdtd.t, self.fdtd.t + dt);
            while s0 < s1 {
                if self.interval.as_ref().is_none_or(|iv| s0 >= iv.t1) {
                    if self.next + 1 >= run.frames.len() {
                        return;
                    }
                    self.interval = Some(self.interval(run, self.next));
                    self.next += 1;
                    continue;
                }
                let iv = self.interval.as_ref().expect("an interval");
                let e = s1.min(iv.t1);
                let span = (iv.t1 - iv.t0).max(1e-300);
                let (g0, g1) = ((s0 - iv.t0) / span, (e - iv.t0) / span);
                let lerp = |a: (f64, f64), b: (f64, f64), g: f64| {
                    (a.0 + (b.0 - a.0) * g, a.1 + (b.1 - a.1) * g)
                };
                let moves: Vec<Move> = iv
                    .paths
                    .iter()
                    .map(|p| (p.q, lerp(p.a, p.b, g0), lerp(p.a, p.b, g1)))
                    .chain(
                        iv.leads
                            .iter()
                            .zip(&self.leads)
                            .filter_map(|(q, l)| l.map(|(a, b)| (q * (g1 - g0), a, b))),
                    )
                    .collect();
                for (q, a, b) in moves {
                    if q != 0.0 {
                        self.fdtd.deposit(q, a, b);
                    }
                }
                s0 = e;
            }
            self.fdtd.step();
            self.radiated += self.outflow() * dt;
            self.history.push((self.fdtd.t, self.radiated));
        }
    }
}
