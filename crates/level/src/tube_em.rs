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
//! board's border; the shortest such segment that crosses no other electrode and no metal
//! (a conductor: Positive carriers' inner Hall contacts were wired through the other
//! strip). A lead
//! carries the electrode's change of charge (from the frames' surface densities) less the
//! particles it took in and gave off, plus what the leads starting on it drew. The
//! electrodes are perfect conductors on the grid, so each one's charge follows the run's.
//! An open level's outermost envelope has no lead: its world conserves its total charge,
//! as a z-invariant world without currents from infinity must (the quasi-static model's
//! open 2D boundary lets the total drift; that drift does not radiate).

use physics::DVec3;
use physics::fdtd2d::{Fdtd, Grid};
use std::collections::HashMap;

use crate::tube::{KIND_NEGATIVE, KIND_PROJECTILE, TubeRun};
use crate::{Level, Outline};

/// Absorbing-layer cells of an open tube (one with an output).
pub const EM_PML: usize = 16;

/// Free space between an open level's board and its absorbing layer (board cells).
pub const EM_MARGIN: f64 = 1.0;
/// Grid cells per board cell.
pub const EM_PER_CELL: f64 = 4.0;
/// The display's speed of light of a tube without an output, in units of its carriers'
/// largest speed (`√(2 |q/m| ΔV)`): fast enough that the field is nearly quasi-static.
pub const EM_C_FACTOR: f64 = 10.0;
/// A closed tube's speed of light also puts its box's lowest mode (`πc/L`, L the box's
/// longer side) at least this many times above the run's fastest frequency
/// (`fastest_frequency`): the box then responds to the run quasi-statically, not as a
/// resonator (Klystron's drive lay 9 % below the lowest mode and rang it).
pub const EM_MODE_FACTOR: f64 = 12.0;
/// A closed tube's field is shown smoothed in time by this many first-order low-passes
/// in cascade (`Smoothing`).
pub const EM_SMOOTH_STAGES: usize = 4;
/// Each stage's time constant in units of `1/ω₁`, ω₁ = πc/L the box's lowest mode.
pub const EM_SMOOTH_RATIO: f64 = 2.0;
/// Every level's field (open ones too) is smoothed at least over the run's sampling: each
/// stage's time constant at least this many frame intervals. The run's currents are known
/// at its frames only; between them they step from frame to frame, and the steps radiate
/// waves at the frames' rate (Magnetron's frames every 0.05: stripes along its antenna's
/// feed at ω ≈ 126), above the frames' Nyquist frequency π/T, where the run has no
/// information. 0.4 passes `2π/T` as 0.02, `π/T` as 0.15, and Magnetron's ω = 3 and 6
/// as 0.993 and 0.97.
pub const EM_FRAME_RATIO: f64 = 0.4;
/// The time step's fraction of the stability limit.
const COURANT: f64 = 0.9;

/// A lead: from its parent's metal (or the border) to its electrode's.
pub type Lead = ((f64, f64), (f64, f64));

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
    /// Per electrode: the charge its lead carries over the interval (from its start, the
    /// parent or the border, to the electrode).
    leads: Vec<f64>,
    /// The current of every path's and lead's whole move as if made in one step
    /// (`Fdtd::deposit_to`; `(edge, J)` of `E_x`'s and `E_y`'s edges): each step takes its
    /// share of the interval.
    jx: Vec<(usize, f64)>,
    jy: Vec<(usize, f64)>,
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
    /// The emitting electrodes and whether their carriers are negative (the cathode, the
    /// ion source, the metals' contacts): a newborn particle comes from the nearest of
    /// its sign.
    emitters: Vec<(usize, bool)>,
    /// The metals' backgrounds (the ions' charge) on the grid's nodes: fixed.
    background: Vec<f64>,
    /// Whether the grid is open (an absorbing layer inside its border).
    open: bool,
    /// The field smoothed in time: over the run's frames, and a closed box's over its
    /// lowest mode (`EM_FRAME_RATIO`, `EM_SMOOTH_RATIO`).
    smoothing: Option<Smoothing>,
    /// The next frame interval to feed: frames `next`..`next + 1`.
    next: usize,
    interval: Option<Interval>,
    /// The energy per unit length that has flowed out through the contour (a rectangle
    /// three cells inside the absorbing layer), and its history `(t, W)`: an open grid's
    /// (a closed box's stays 0).
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

/// The fastest frequency of a tube level's run in a field `b_z` along z, as far as its
/// level fixes it: its sources' (a sine's ω, a pulse's edges' `1/rise`) and its carriers'
/// cyclotron frequency `|q/m| B`.
#[must_use]
pub fn fastest_frequency(level: &Level, b_z: f64) -> f64 {
    let Some(spec) = level.tube.as_ref() else {
        return 0.0;
    };
    let sources = spec
        .circuits
        .iter()
        .map(|c| match c.source {
            crate::Source::Dc { .. } => 0.0,
            crate::Source::Sine { omega, .. } => omega.abs(),
            crate::Source::Pulse { rise, fall, .. } => [rise, fall]
                .into_iter()
                .filter(|&s| s > 0.0)
                .map(|s| 1.0 / s)
                .fold(0.0, f64::max),
        })
        .fold(0.0, f64::max);
    let ratio = spec
        .ions
        .iter()
        .map(|i| i.charge_per_mass.abs())
        .chain(spec.metals.iter().map(|m| m.charge_per_mass.abs()))
        .fold(spec.charge_per_mass.abs(), f64::max);
    sources.max(ratio * b_z.abs())
}

/// The field smoothed in time: `EM_SMOOTH_STAGES` first-order low-passes in cascade, each
/// step's field into the first (`y += α (x − y)`, `α = 1 − e^{−dt/τ}`), passing a
/// frequency ω as `(1 + (ωτ)²)^{−n/2}` and static fields unchanged. Every level's over the
/// run's frames (`EM_FRAME_RATIO`: its currents step from frame to frame, and the steps
/// radiate at the frames' rate). A closed box's more (`EM_SMOOTH_RATIO`): a lossless box
/// keeps every wave radiated in it, and the run's macroparticles radiate far more than its
/// real carriers would (their weight times: each one's impact on an electrode, its
/// graininess, its path's kinks between frames), so the raw field's fluctuations grow
/// without bound; they ring at the box's modes, `ω ≥ ω₁ = πc/L`. Being linear, it is the
/// field of the run's currents smoothed by the same kernel (a gamma distribution of mean
/// delay `nτ`): the field shown is taken that much later.
#[derive(Clone, Debug)]
struct Smoothing {
    alpha: f64,
    /// Per field (`E_x`, `E_y`, `B_z`): every stage's value, element by element
    /// (`EM_SMOOTH_STAGES` in a row: one pass takes a step through all of them).
    stages: [Vec<f64>; 3],
    /// Per field: the last stage's values (shown).
    out: [Vec<f64>; 3],
    /// The kernel's mean delay (each stage's `dt (1 − α)/α`).
    delay: f64,
    /// Per lead: its current through the same stages, the last one's (shown), the charge
    /// that has passed at that current since the start, and the largest shown current.
    lead_stages: Vec<[f64; EM_SMOOTH_STAGES]>,
    lead_out: Vec<f64>,
    lead_charge: Vec<f64>,
    lead_peak: f64,
}

/// Elements per parallel task of the smoothing's pass.
const SMOOTH_CHUNK: usize = 2048;

impl Smoothing {
    /// Stages of time constant `tau`, all starting at `fdtd`'s field (a static start
    /// stays as it is), with `leads` leads' currents starting at 0.
    fn new(fdtd: &Fdtd, tau: f64, leads: usize) -> Self {
        let alpha = 1.0 - libm::exp(-fdtd.dt / tau);
        #[allow(clippy::cast_precision_loss)]
        let delay = EM_SMOOTH_STAGES as f64 * fdtd.dt * (1.0 - alpha) / alpha;
        let fields = [&fdtd.ex, &fdtd.ey, &fdtd.bz];
        Self {
            alpha,
            stages: fields.map(|f| {
                f.iter()
                    .flat_map(|&v| std::iter::repeat_n(v, EM_SMOOTH_STAGES))
                    .collect()
            }),
            out: fields.map(Clone::clone),
            delay,
            lead_stages: vec![[0.0; EM_SMOOTH_STAGES]; leads],
            lead_out: vec![0.0; leads],
            lead_charge: vec![0.0; leads],
            lead_peak: 0.0,
        }
    }

    /// Takes in a step's lead currents (`dt` long), as the field.
    fn update_leads(&mut self, currents: &[f64], dt: f64) {
        let a = self.alpha;
        for (((ys, o), q), &i) in self
            .lead_stages
            .iter_mut()
            .zip(&mut self.lead_out)
            .zip(&mut self.lead_charge)
            .zip(currents)
        {
            let mut v = i;
            for y in ys {
                *y += a * (v - *y);
                v = *y;
            }
            *o = v;
            *q += v * dt;
            self.lead_peak = self.lead_peak.max(v.abs());
        }
    }

    /// Takes in a step's field: the three fields' chunks in one parallel pass (a pass per
    /// field paid the thread pool's start three times a step).
    fn update(&mut self, fdtd: &Fdtd) {
        use rayon::prelude::*;
        const N: usize = EM_SMOOTH_STAGES;
        let a = self.alpha;
        let inputs = [&fdtd.ex, &fdtd.ey, &fdtd.bz];
        let mut tasks = Vec::new();
        for ((stages, out), x) in self.stages.iter_mut().zip(&mut self.out).zip(inputs) {
            tasks.extend(
                stages
                    .chunks_mut(N * SMOOTH_CHUNK)
                    .zip(out.chunks_mut(SMOOTH_CHUNK))
                    .zip(x.chunks(SMOOTH_CHUNK)),
            );
        }
        tasks.into_par_iter().for_each(|((ys, out), x)| {
            for ((ys, o), &v) in ys.as_chunks_mut::<N>().0.iter_mut().zip(out).zip(x) {
                let mut v = v;
                for y in ys {
                    *y += a * (v - *y);
                    v = *y;
                }
                *o = v;
            }
        });
    }
}

/// The field a map shows (`TubeEm::shown`): its time and `E_x`, `E_y`, `B_z` on the
/// field's grid (`Fdtd`'s layout).
pub struct Shown<'a> {
    pub t: f64,
    pub ex: &'a [f64],
    pub ey: &'a [f64],
    pub bz: &'a [f64],
}

impl TubeEm {
    /// The field of `level`'s tube run `run` (its first frame at least), at the display's
    /// speed of light (the output's, or `EM_C_FACTOR` times the carriers' largest speed
    /// and enough for `EM_MODE_FACTOR`), on a grid of `EM_PER_CELL` cells per board cell
    /// over the board, open (an absorbing layer inside the border) when the level has an
    /// output, otherwise a closed box, shown smoothed in time (`Smoothing`).
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
        // The carriers' largest speed: what the potentials give them, or a metal's thermal
        // speed (4 √(kT/m): the fastest of its thousands of carriers), or a projectile's.
        let (lo, hi) = run.potential_range;
        let speed = (2.0 * spec.charge_per_mass.abs() * (hi - lo).abs())
            .sqrt()
            .max(1.0);
        let speed = spec
            .metals
            .iter()
            .map(|m| 4.0 * m.thermal_speed)
            .chain(spec.projectiles.iter().map(|p| p.speed))
            .fold(speed, f64::max);
        let open = spec.output.is_some();
        let b = level.bounds();
        let longer = (b.max.x - b.min.x).max(b.max.y - b.min.y);
        let c = spec.output.map_or_else(
            || {
                (EM_C_FACTOR * speed).max(
                    EM_MODE_FACTOR
                        * fastest_frequency(level, run.b_wave.map_or(run.b_z, |w| w.peak()))
                        * longer
                        / std::f64::consts::PI,
                )
            },
            |o| o.c,
        );
        let h = 1.0 / EM_PER_CELL;
        // An open level's grid reaches past the board: `EM_MARGIN` of free space, then the
        // absorbing layer, so that all the board shows is free field (with the layer on
        // the board its outer cells damped the waves: the fronts seemed to stick to the
        // edge).
        #[allow(clippy::cast_precision_loss)]
        let margin = if open {
            EM_MARGIN + EM_PML as f64 * h
        } else {
            0.0
        };
        let nx = ((b.max.x - b.min.x + 2.0 * margin) / h).ceil() as usize;
        let ny = ((b.max.y - b.min.y + 2.0 * margin) / h).ceil() as usize;
        let grid = Grid {
            nx,
            ny,
            h,
            origin: (b.min.x - margin, b.min.y - margin),
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
        // material (or the border) crossing no other electrode and no metal (its own
        // included), over 144 directions.
        let in_metal = |p: DVec3| {
            run.metals
                .iter()
                .any(|m| p.x >= m.lo.x && p.x <= m.hi.x && p.y >= m.lo.y && p.y <= m.hi.y)
        };
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
                            || in_metal(p)
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
        // The metals' backgrounds: each cell's share of a rectangle, at the share's centre.
        let mut background = vec![0.0; grid.nodes()];
        for m in &run.metals {
            let h = grid.h;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let (i0, i1, j0, j1) = (
                ((m.lo.x - grid.origin.0) / h).floor().max(0.0) as usize,
                ((m.hi.x - grid.origin.0) / h).ceil().max(0.0) as usize,
                ((m.lo.y - grid.origin.1) / h).floor().max(0.0) as usize,
                ((m.hi.y - grid.origin.1) / h).ceil().max(0.0) as usize,
            );
            for j in j0..j1.min(grid.ny) {
                for i in i0..i1.min(grid.nx) {
                    #[allow(clippy::cast_precision_loss)]
                    let (cx0, cy0) = grid.node(i, j);
                    let (ax, bx) = (cx0.max(m.lo.x), (cx0 + h).min(m.hi.x));
                    let (ay, by) = (cy0.max(m.lo.y), (cy0 + h).min(m.hi.y));
                    if bx > ax && by > ay {
                        fdtd.add_charge(
                            &mut background,
                            m.background * (bx - ax) * (by - ay),
                            0.5 * (ax + bx),
                            0.5 * (ay + by),
                        );
                    }
                }
            }
        }
        // The start: the electrostatic field of the first frame's potentials and charges.
        let mut rho = background.clone();
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
        // Smoothed over the run's frames at least, and a closed box over its lowest mode.
        let frame = run.frames.get(1).map_or(spec.step, |f| f.t - first.t);
        let tau = (EM_FRAME_RATIO * frame).max(if open {
            0.0
        } else {
            EM_SMOOTH_RATIO * longer / (std::f64::consts::PI * c)
        });
        let smoothing = Some(Smoothing::new(&fdtd, tau, n_el));
        Self {
            fdtd,
            outlines,
            segments,
            node_owner,
            parent,
            leads,
            order,
            emitters: std::iter::once((
                spec.cathode,
                spec.metal_of(spec.cathode)
                    .map_or(spec.charge_per_mass, |m| spec.metals[m].charge_per_mass)
                    < 0.0,
            ))
            .chain(spec.ions.map(|i| (i.electrode, false)))
            .chain(spec.metals.iter().flat_map(|m| {
                m.contacts
                    .iter()
                    .map(move |&e| (e, m.charge_per_mass < 0.0))
            }))
            .collect(),
            background,
            open,
            smoothing,
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
            if kind == KIND_PROJECTILE {
                continue;
            }
            // From the nearest emitter of its sign.
            let negative = kind == KIND_NEGATIVE;
            let Some(e) = self
                .emitters
                .iter()
                .filter(|(_, n)| *n == negative)
                .map(|(e, _)| *e)
                .min_by(|&a, &b| {
                    self.outlines[a]
                        .distance_to(x)
                        .total_cmp(&self.outlines[b].distance_to(x))
                })
            else {
                continue;
            };
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
        // The interval's current, deposited once (it was each step's share of every move,
        // a step's cost growing with the particles: Klystron's field ran slower than its
        // playback). A path's charge then passes from its start to its end in proportion
        // to the time (`charge_density`), not along the path: within a frame interval, a
        // few cells at most.
        let g = self.fdtd.grid;
        let mut dense_x = vec![0.0; g.nx * (g.ny + 1)];
        let mut dense_y = vec![0.0; (g.nx + 1) * g.ny];
        let moves = paths.iter().map(|p| (p.q, p.a, p.b)).chain(
            leads
                .iter()
                .zip(&self.leads)
                .filter_map(|(q, l)| l.map(|(a, b)| (*q, a, b))),
        );
        for (q, a, b) in moves {
            if q != 0.0 {
                self.fdtd.deposit_to(q, a, b, &mut dense_x, &mut dense_y);
            }
        }
        let sparse = |d: Vec<f64>| -> Vec<(usize, f64)> {
            d.into_iter()
                .enumerate()
                .filter(|(_, j)| *j != 0.0)
                .collect()
        };
        Interval {
            t0: f0.t,
            t1: f1.t,
            paths,
            leads,
            jx: sparse(dense_x),
            jy: sparse(dense_y),
        }
    }

    /// The charge density at the grid's nodes the field's time implies: every path's
    /// charge where it is now (Gauss's law must give it in vacuum).
    #[must_use]
    pub fn charge_density(&self) -> Vec<f64> {
        let mut rho = self.background.clone();
        if let Some(iv) = &self.interval {
            // Each path's charge passes from its start to its end in proportion to the time
            // (the interval's current is spread evenly over it).
            let g = ((self.fdtd.t - iv.t0) / (iv.t1 - iv.t0).max(1e-300)).clamp(0.0, 1.0);
            for p in &iv.paths {
                self.fdtd
                    .add_charge(&mut rho, p.q * (1.0 - g), p.a.0, p.a.1);
                self.fdtd.add_charge(&mut rho, p.q * g, p.b.0, p.b.1);
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

    /// Per electrode: the current its lead carries now (over the frame interval being fed;
    /// positive from the lead's start, its parent or the border, to the electrode).
    #[must_use]
    pub fn lead_currents(&self) -> Vec<f64> {
        self.interval.as_ref().map_or_else(
            || vec![0.0; self.outlines.len()],
            |iv| {
                let span = (iv.t1 - iv.t0).max(1e-300);
                iv.leads.iter().map(|q| q / span).collect()
            },
        )
    }

    /// Per electrode: its lead's current as the shown field sees it (smoothed the same
    /// way), the charge that has passed so, and the largest such current so far.
    #[must_use]
    pub fn shown_leads(&self) -> (&[f64], &[f64], f64) {
        self.smoothing
            .as_ref()
            .map_or((&[][..], &[][..], 0.0), |s| {
                (&s.lead_out[..], &s.lead_charge[..], s.lead_peak)
            })
    }

    /// How much later than the field's own time the shown field is (`Smoothing`; 0 when
    /// shown as it is).
    #[must_use]
    pub fn delay(&self) -> f64 {
        self.smoothing.as_ref().map_or(0.0, |s| s.delay)
    }

    /// The field the maps show: smoothed in time (`Smoothing`), at the field's time less
    /// the smoothing's delay.
    #[must_use]
    pub fn shown(&self) -> Shown<'_> {
        let f = &self.fdtd;
        match &self.smoothing {
            Some(s) => {
                let [ex, ey, bz] = &s.out;
                Shown {
                    t: f.t - s.delay,
                    ex,
                    ey,
                    bz,
                }
            }
            None => Shown {
                t: f.t,
                ex: &f.ex,
                ey: &f.ey,
                bz: &f.bz,
            },
        }
    }

    /// Advances the field until the shown field reaches time `t` (as far as the run's
    /// frames reach).
    pub fn show(&mut self, run: &TubeRun, t: f64) {
        self.advance_to(run, t + self.delay());
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
                self.fdtd.add_current(&iv.jx, &iv.jy, (e - s0) / span);
                s0 = e;
            }
            self.fdtd.step();
            if let Some(s) = &mut self.smoothing {
                s.update(&self.fdtd);
                let currents = self.interval.as_ref().map_or_else(Vec::new, |iv| {
                    let span = (iv.t1 - iv.t0).max(1e-300);
                    iv.leads.iter().map(|q| q / span).collect()
                });
                s.update_leads(&currents, dt);
            }
            // A closed box radiates nothing out (its border is a conductor).
            if self.open {
                self.radiated += self.outflow() * dt;
                self.history.push((self.fdtd.t, self.radiated));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Least-squares fit of `y ≈ a₀ + a cos ωt + b sin ωt`: `(a₀, a, −b)`, the static part
    /// and the complex amplitude `C = a − i b` of `Re[C e^{iωt}]`.
    #[allow(clippy::needless_range_loop)]
    fn fit(samples: &[(f64, f64)], omega: f64) -> (f64, f64, f64) {
        let mut m = [[0.0; 3]; 3];
        let mut r = [0.0; 3];
        for &(t, y) in samples {
            let f = [1.0, libm::cos(omega * t), libm::sin(omega * t)];
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] += f[i] * f[j];
                }
                r[i] += f[i] * y;
            }
        }
        // Gaussian elimination (the normal equations are well conditioned over whole
        // periods).
        for k in 0..3 {
            for i in k + 1..3 {
                let q = m[i][k] / m[k][k];
                for j in k..3 {
                    m[i][j] -= q * m[k][j];
                }
                r[i] -= q * r[k];
            }
        }
        let mut x = [0.0; 3];
        for i in (0..3).rev() {
            let s: f64 = (i + 1..3).map(|j| m[i][j] * x[j]).sum();
            x[i] = (r[i] - s) / m[i][i];
        }
        (x[0], x[1], -x[2])
    }

    /// F7: a closed box's smoothing (`Smoothing`) on F1's ringing TE₁₀ mode (a perfectly
    /// conducting square box, side 20, h = 0.5, c = 1) over a uniform `B_z` of 3 (a static
    /// solution), at `ωτ = EM_SMOOTH_RATIO`. After the transients (70 time constants) the
    /// shown field must be the static part unchanged and the mode times the cascade's exact
    /// discrete transfer `H(ω)ⁿ`, `H = α/(1 − (1 − α) e^{−iωΔt})` at the scheme's
    /// discrete ω: the complex amplitudes fitted over 10 periods, required 1e-6 (the
    /// transients' remainder and rounding), the static part 1e-9.
    #[test]
    fn f7_smoothing_of_a_cavity_mode() {
        let (n, h, c, l, b0) = (40, 0.5, 1.0, 20.0, 3.0);
        let grid = Grid {
            nx: n,
            ny: n,
            h,
            origin: (0.0, 0.0),
        };
        let mut f = Fdtd::new(grid, c, 0.9, 0, vec![false; grid.nodes()]);
        let k = PI / l;
        for j in 0..n {
            for i in 0..n {
                #[allow(clippy::cast_precision_loss)]
                let x = (i as f64 + 0.5) * h;
                f.bz[i + n * j] = b0 + libm::cos(k * x);
            }
        }
        let omega = 2.0 / f.dt * libm::asin(c * f.dt / h * libm::sin(k * h / 2.0));
        let tau = EM_SMOOTH_RATIO / omega;
        let mut s = Smoothing::new(&f, tau, 0);
        let period = 2.0 * PI / omega;
        let start = 70.0 * tau;
        let cell = 3 + n * 7;
        let (mut raw, mut shown) = (Vec::new(), Vec::new());
        while f.t < start + 10.0 * period {
            f.step();
            s.update(&f);
            if f.t >= start {
                raw.push((f.t, f.bz[cell]));
                shown.push((f.t, s.out[2][cell]));
            }
        }
        let (r0, rr, ri) = fit(&raw, omega);
        let (s0, sr, si) = fit(&shown, omega);
        // H(ω)ⁿ.
        let a = s.alpha;
        let (cr, ci) = (
            1.0 - (1.0 - a) * libm::cos(omega * f.dt),
            (1.0 - a) * libm::sin(omega * f.dt),
        );
        let d = cr * cr + ci * ci;
        let (hr, hi) = (a * cr / d, -a * ci / d);
        let (mut pr, mut pi) = (1.0, 0.0);
        for _ in 0..EM_SMOOTH_STAGES {
            (pr, pi) = (pr * hr - pi * hi, pr * hi + pi * hr);
        }
        // The measured transfer: shown / raw.
        let q = rr * rr + ri * ri;
        let (mr, mi) = ((sr * rr + si * ri) / q, (si * rr - sr * ri) / q);
        let error = libm::hypot(mr - pr, mi - pi) / libm::hypot(pr, pi);
        #[allow(clippy::cast_precision_loss)]
        let continuum = libm::pow(
            1.0 + EM_SMOOTH_RATIO * EM_SMOOTH_RATIO,
            -0.5 * EM_SMOOTH_STAGES as f64,
        );
        println!(
            "F7: |H^n| measured {:.8}, exact {:.8} (continuum {continuum:.8}), relative \
             error {error:.1e}; static part {s0} for {r0} (raw)",
            libm::hypot(mr, mi),
            libm::hypot(pr, pi),
        );
        assert!(error < 1e-6, "F7 transfer {error:e}");
        assert!((s0 / b0 - 1.0).abs() < 1e-9, "F7 static {s0}");
        assert!((r0 / b0 - 1.0).abs() < 1e-9, "F7 raw static {r0}");
    }
}
