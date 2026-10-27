//! Leads: the wires that bring each electrode of a tube level its circuit current, and the
//! straps that join a floating circuit's electrodes, as the level's geometry: the field
//! `tube_em` drives along them and the game draws them in every view. A lead comes from the
//! innermost grounded ring enclosing its electrode (an envelope), or, in an open level (one
//! with an output), from the nearest envelope outside it, or else from the board's border;
//! an open level's outermost envelope has none (the owner's convention, 2026-10-12).
//!
//! They are routed as circuit diagrams draw wires (the owner, 2026-10-27): along x and y with
//! right-angle bends, clear of every electrode and metal and of the wires routed before. On a
//! grid `ROUTE_SPACING` apart over the board, a shortest-path search (in whole grid steps)
//! charges a bend `BEND` steps, a node another wire took `CROSS` and a stub meeting its
//! surface on a slant `SLANT` (and on a box a step for every grid step off the face's middle),
//! so that a wire bends only to go round something, crosses another only where there is no
//! way round, meets a box's face (never its corner) in the middle and a round electrode on
//! its axis where it can (Magnetron's antenna feed: straight along the rod's axis); a node
//! closer
//! than `CLEARANCE` to an electrode or a metal is closed to every wire but those ending on it.
//! Only a wire's two ends leave the grid: a short stub onto each surface (perpendicular to a
//! flat face, radial to a round one). Straps are routed first, then the leads, the shortest
//! (by their ends' distance) first. Deterministic: integer costs, ties by the state's index.

use crate::geometry::segment_distance;
use crate::{Element, Level, Outline};
use physics::DVec3;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// The routing grid's spacing (board units: half a cell, so that a wire runs between two
/// electrodes a cell apart).
pub const ROUTE_SPACING: f64 = 0.5;
/// A wire keeps at least this far from every electrode and metal but those it ends on.
pub const CLEARANCE: f64 = 0.3;
/// A bend's cost in grid steps.
const BEND: u32 = 4;
/// The cost of a node another wire took (a crossing), in grid steps.
const CROSS: u32 = 60;
/// The cost of a stub that meets its surface on a slant (a round electrode off its axes, a
/// turned box), in grid steps: among routes alike the one meeting both ends square on. (A
/// box's face costs besides a step for every grid step off its middle: `end_cost`.)
const SLANT: u32 = 2;

/// A wire: its corners in order, the first and last on the surfaces it joins.
pub type Wire = Vec<(f64, f64)>;

/// A tube level's wiring for one placement.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LeadPlan {
    /// Per level electrode: the enclosing ring its lead comes from (None: the border, or no
    /// lead).
    pub parent: Vec<Option<usize>>,
    /// Per level electrode: whether it encloses another (an envelope).
    pub envelope: Vec<bool>,
    /// Per level electrode: its lead, from its parent's surface or the border to its own
    /// surface, if it has one (None also where no route exists).
    pub leads: Vec<Option<Wire>>,
    /// The floating circuits' straps: from a circuit's electrode to each one linked to it.
    pub straps: Vec<Wire>,
    /// With the carriers' own field (`TubeSpec::speed_of_light`; DEC-27): the ground
    /// terminal, on the grounded envelope (`tube::ENVELOPE_MARGIN` beyond the board) below
    /// the middle of the board's lower edge, where every circuit from the border returns.
    pub ground: Option<(f64, f64)>,
    /// With it, per level electrode with a lead from the border: the lead's return, along
    /// the envelope's inner face from the ground terminal, the shorter way round, then across
    /// the margin to where the lead starts. The envelope is then the circuits' ground wire as
    /// it is the electrostatics' ground, and each circuit a closed loop (an open lead's
    /// current, coming from nowhere, would have a field whose energy grows without bound with
    /// the distance: no inductance of its own).
    pub returns: Vec<Option<Wire>>,
}

/// The path along the border of the rectangle `lo`–`hi` from `a` (on it) to `b` (on it, or
/// near it: the path ends at its foot on the nearest side, then at `b`), the shorter way
/// round (counter-clockwise on a tie): its corners.
fn border_path(lo: DVec3, hi: DVec3, a: DVec3, b: DVec3) -> Vec<DVec3> {
    let (w, h) = (hi.x - lo.x, hi.y - lo.y);
    let perimeter = 2.0 * (w + h);
    // The foot of `p` on the nearest side, and its distance along the border
    // counter-clockwise from `lo`.
    let foot = |p: DVec3| {
        let sides = [
            (p.y - lo.y).abs(),
            (hi.x - p.x).abs(),
            (hi.y - p.y).abs(),
            (p.x - lo.x).abs(),
        ];
        let side = (0..4)
            .min_by(|&i, &j| sides[i].total_cmp(&sides[j]))
            .unwrap_or(0);
        let (x, y) = (p.x.clamp(lo.x, hi.x), p.y.clamp(lo.y, hi.y));
        match side {
            0 => (DVec3::new(x, lo.y, 0.0), x - lo.x),
            1 => (DVec3::new(hi.x, y, 0.0), w + (y - lo.y)),
            2 => (DVec3::new(x, hi.y, 0.0), w + h + (hi.x - x)),
            _ => (DVec3::new(lo.x, y, 0.0), 2.0 * w + h + (hi.y - y)),
        }
    };
    let (fa, sa) = foot(a);
    let (fb, sb) = foot(b);
    let corners = [
        (0.0, DVec3::new(lo.x, lo.y, 0.0)),
        (w, DVec3::new(hi.x, lo.y, 0.0)),
        (w + h, DVec3::new(hi.x, hi.y, 0.0)),
        (2.0 * w + h, DVec3::new(lo.x, hi.y, 0.0)),
    ];
    let forward = (sb - sa).rem_euclid(perimeter);
    let mut out = vec![fa];
    // The corners strictly between, in the order passed, over one turn either way.
    let mut passed: Vec<(f64, DVec3)> = Vec::new();
    for m in [-1.0, 0.0, 1.0] {
        for &(s, c) in &corners {
            let s = s + m * perimeter;
            let between = if forward <= 0.5 * perimeter {
                s > sa + 1e-9 && s < sa + forward - 1e-9
            } else {
                s < sa - 1e-9 && s > sa - (perimeter - forward) + 1e-9
            };
            if between {
                passed.push((s, c));
            }
        }
    }
    if forward <= 0.5 * perimeter {
        passed.sort_by(|x, y| x.0.total_cmp(&y.0));
    } else {
        passed.sort_by(|x, y| y.0.total_cmp(&x.0));
    }
    out.extend(passed.into_iter().map(|(_, c)| c));
    out.push(fb);
    out.push(b);
    out.dedup_by(|p, q| (*p - *q).length() < 1e-9);
    out
}

/// What a wire ends on.
#[derive(Clone, Copy, Debug, PartialEq)]
enum End {
    /// An electrode's outer surface.
    Electrode(usize),
    /// An enclosing ring's inner surface.
    Inside(usize),
    /// The board's border.
    Border,
}

/// An outline's centre (a box's middle).
pub(crate) fn centre(o: &Outline) -> DVec3 {
    match *o {
        Outline::Rect(r) => (r[0] + r[1] + r[2] + r[3]) * 0.25,
        Outline::Circle { center, .. } | Outline::Ring { center, .. } => center,
    }
}

/// The farthest distance from `c` to an outline's material.
pub(crate) fn farthest(o: &Outline, c: DVec3) -> f64 {
    match *o {
        Outline::Rect(r) => r.iter().map(|p| (*p - c).length()).fold(0.0, f64::max),
        Outline::Circle { center, radius } => (center - c).length() + radius,
        Outline::Ring { center, outer, .. } => (center - c).length() + outer,
    }
}

/// The point of an outline's material nearest to `p` (outside it): a box's face, a rod's or a
/// ring's circle (a ring's inner one from inside its hole).
fn surface_point(o: &Outline, p: DVec3) -> DVec3 {
    let radial = |center: DVec3, r: f64| {
        let d = p - center;
        let l = d.length();
        if l < 1e-12 {
            center + DVec3::new(r, 0.0, 0.0)
        } else {
            center + d * (r / l)
        }
    };
    match *o {
        Outline::Rect(r) => (0..4)
            .map(|i| {
                let (a, b) = (r[i], r[(i + 1) % 4]);
                let d = b - a;
                let s = ((p - a).dot(d) / d.length_squared().max(1e-300)).clamp(0.0, 1.0);
                a + d * s
            })
            .min_by(|a, b| (*a - p).length().total_cmp(&(*b - p).length()))
            .unwrap_or(p),
        Outline::Circle { center, radius } => radial(center, radius),
        Outline::Ring {
            center,
            inner,
            outer,
        } => {
            if (p - center).length() < 0.5 * (inner + outer) {
                radial(center, inner)
            } else {
                radial(center, outer)
            }
        }
    }
}

/// Each electrode's enclosing ring (its lead's start) and whether it is an envelope: the
/// innermost ring enclosing it, or in an open level the nearest envelope outside it (but
/// for an envelope), or none (the border).
pub(crate) fn topology(outlines: &[Outline], open: bool) -> (Vec<Option<usize>>, Vec<bool>) {
    let n = outlines.len();
    let encloses = |p: usize, e: usize| match outlines[p] {
        Outline::Ring { center, inner, .. } => p != e && farthest(&outlines[e], center) < inner,
        _ => false,
    };
    let envelope: Vec<bool> = (0..n).map(|p| (0..n).any(|e| encloses(p, e))).collect();
    let parent = (0..n)
        .map(|e| {
            let inner = (0..n).filter(|&p| encloses(p, e)).min_by(|&p, &q| {
                let r = |k: usize| match outlines[k] {
                    Outline::Ring { inner, .. } => inner,
                    _ => f64::INFINITY,
                };
                r(p).total_cmp(&r(q))
            });
            inner.or_else(|| {
                if open && !envelope[e] {
                    (0..n).filter(|&p| envelope[p] && p != e).min_by(|&p, &q| {
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
    (parent, envelope)
}

/// The routing grid over the board with its obstacles and the nodes wires have taken.
struct Router<'a> {
    origin: DVec3,
    nx: usize,
    ny: usize,
    outlines: &'a [Outline],
    /// Per node: inside some material (closed to every wire).
    inside: Vec<bool>,
    /// Per node: a metal within `CLEARANCE` (closed to every wire).
    metal_near: Vec<bool>,
    /// Per node: the electrodes within `CLEARANCE` (open only to a wire ending on them).
    near: Vec<Vec<usize>>,
    /// Per node: a wire took it.
    taken: Vec<bool>,
}

impl<'a> Router<'a> {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    fn new(lo: DVec3, hi: DVec3, outlines: &'a [Outline], metals: &[(DVec3, DVec3)]) -> Self {
        // In the plane of the outlines (the board's bounds reach z = ±1).
        let lo = DVec3::new(lo.x, lo.y, 0.0);
        let nx = libm::round((hi.x - lo.x) / ROUTE_SPACING).max(1.0) as usize;
        let ny = libm::round((hi.y - lo.y) / ROUTE_SPACING).max(1.0) as usize;
        let nodes = (nx + 1) * (ny + 1);
        let (mut inside, mut metal_near) = (vec![false; nodes], vec![false; nodes]);
        let mut near = vec![Vec::new(); nodes];
        for k in 0..nodes {
            let p = lo
                + DVec3::new(
                    (k % (nx + 1)) as f64 * ROUTE_SPACING,
                    (k / (nx + 1)) as f64 * ROUTE_SPACING,
                    0.0,
                );
            for (e, o) in outlines.iter().enumerate() {
                let d = o.distance_to(p);
                if d <= 0.0 {
                    inside[k] = true;
                } else if d < CLEARANCE {
                    near[k].push(e);
                }
            }
            for (a, b) in metals {
                let dx = (a.x - p.x).max(p.x - b.x).max(0.0);
                let dy = (a.y - p.y).max(p.y - b.y).max(0.0);
                let d = libm::hypot(dx, dy);
                if d <= 0.0 {
                    inside[k] = true;
                } else if d < CLEARANCE {
                    metal_near[k] = true;
                }
            }
        }
        Self {
            origin: lo,
            nx,
            ny,
            outlines,
            inside,
            metal_near,
            near,
            taken: vec![false; nodes],
        }
    }

    #[allow(clippy::cast_precision_loss)]
    fn point(&self, k: usize) -> DVec3 {
        self.origin
            + DVec3::new(
                (k % (self.nx + 1)) as f64 * ROUTE_SPACING,
                (k / (self.nx + 1)) as f64 * ROUTE_SPACING,
                0.0,
            )
    }

    /// Whether a wire between `a` and `b` may pass node `k`.
    fn open(&self, k: usize, a: End, b: End) -> bool {
        let own = |e: usize| {
            [a, b]
                .iter()
                .any(|x| matches!(*x, End::Electrode(f) | End::Inside(f) if f == e))
        };
        !self.inside[k] && !self.metal_near[k] && self.near[k].iter().all(|&e| own(e))
    }

    /// Whether node `k` may start or end a wire at `end`: next to its surface (on the right
    /// side of a ring), or on the border.
    fn touches(&self, k: usize, end: End) -> bool {
        let p = self.point(k);
        match end {
            End::Border => {
                let (i, j) = (k % (self.nx + 1), k / (self.nx + 1));
                i == 0 || j == 0 || i == self.nx || j == self.ny
            }
            End::Electrode(e) | End::Inside(e) => {
                let o = &self.outlines[e];
                if o.distance_to(p) > CLEARANCE + ROUTE_SPACING {
                    return false;
                }
                match (*o, end) {
                    (Outline::Ring { center, inner, .. }, End::Inside(_)) => {
                        (p - center).length() < inner
                    }
                    (Outline::Ring { center, outer, .. }, End::Electrode(_)) => {
                        (p - center).length() > outer
                    }
                    // A box square to the axes: only where the stub meets a face square on
                    // (a node off its corner would join it on a slant), and on the face, not
                    // at a corner (a wire there would run on along the box's edge, as if it
                    // were one with it).
                    (Outline::Rect(r), _)
                        if (0..4).all(|i| {
                            let d = r[(i + 1) % 4] - r[i];
                            d.x.abs() < 1e-9 || d.y.abs() < 1e-9
                        }) =>
                    {
                        let s = surface_point(o, p);
                        // Clear of the corners by the clearance, or by half a face shorter
                        // than twice that (a small box is met at its face's middle).
                        let shortest = (0..4)
                            .map(|i| (r[(i + 1) % 4] - r[i]).length())
                            .fold(f64::INFINITY, f64::min);
                        let keep = CLEARANCE.min(0.5 * shortest);
                        ((s.x - p.x).abs() < 1e-9 || (s.y - p.y).abs() < 1e-9)
                            && r.iter().all(|c| (s - *c).length() >= keep - 1e-9)
                    }
                    _ => true,
                }
            }
        }
    }

    /// The cost of the stub from node `k` onto `end`: `SLANT` on a slant (0 square on, along
    /// x or y, or on the border), and on a box a step for every grid step it meets the face
    /// off its middle (a wire meets a box's face as a circuit diagram's meets a terminal).
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn end_cost(&self, k: usize, end: End) -> u32 {
        let p = self.point(k);
        match end {
            End::Border => 0,
            End::Electrode(e) | End::Inside(e) => {
                let o = &self.outlines[e];
                let s = surface_point(o, p);
                let slant = if (s.x - p.x).abs() < 1e-9 || (s.y - p.y).abs() < 1e-9 {
                    0
                } else {
                    SLANT
                };
                let off = match *o {
                    Outline::Rect(r) => (0..4)
                        .map(|i| (r[i], r[(i + 1) % 4]))
                        .min_by(|(a, b), (c, d)| {
                            segment_distance(s, *a, *b).total_cmp(&segment_distance(s, *c, *d))
                        })
                        .map_or(0, |(a, b)| {
                            libm::round((s - (a + b) * 0.5).length() / ROUTE_SPACING) as u32
                        }),
                    _ => 0,
                };
                slant + off
            }
        }
    }

    /// The cheapest route from `a` to `b` (each end's stub included), its nodes then taken;
    /// None if there is none.
    #[allow(clippy::cast_possible_truncation)]
    fn route(&mut self, a: End, b: End) -> Option<Wire> {
        let (w, nodes) = (self.nx + 1, (self.nx + 1) * (self.ny + 1));
        // States: node × arrival direction (+x, −x, +y, −y, none yet).
        let mut cost = vec![u32::MAX; nodes * 5];
        let mut prev = vec![u32::MAX; nodes * 5];
        let mut heap = BinaryHeap::new();
        for k in 0..nodes {
            if self.open(k, a, b) && self.touches(k, a) {
                let c = self.end_cost(k, a);
                cost[k * 5 + 4] = c;
                heap.push(Reverse((c, (k * 5 + 4) as u32)));
            }
        }
        // The best finish so far: its cost with the end's stub, and its state.
        let mut found: Option<(u32, usize)> = None;
        while let Some(Reverse((c, s))) = heap.pop() {
            let s = s as usize;
            if found.is_some_and(|(f, _)| c >= f) {
                break;
            }
            if c > cost[s] {
                continue;
            }
            let (k, dir) = (s / 5, s % 5);
            if self.touches(k, b) {
                let f = c + self.end_cost(k, b);
                if found.is_none_or(|(g, _)| f < g) {
                    found = Some((f, s));
                }
            }
            let (i, j) = (k % w, k / w);
            let steps = [
                (i < self.nx).then(|| k + 1),
                (i > 0).then(|| k - 1),
                (j < self.ny).then(|| k + w),
                (j > 0).then(|| k - w),
            ];
            for (d, next) in steps.into_iter().enumerate() {
                let Some(n) = next else { continue };
                if !self.open(n, a, b) {
                    continue;
                }
                let mut nc = c + 1;
                if dir != 4 && dir != d {
                    nc += BEND;
                }
                if self.taken[n] {
                    nc += CROSS;
                }
                let ns = n * 5 + d;
                if nc < cost[ns] {
                    cost[ns] = nc;
                    prev[ns] = s as u32;
                    heap.push(Reverse((nc, ns as u32)));
                }
            }
        }
        let mut s = found?.1;
        let mut path = vec![s / 5];
        while prev[s] != u32::MAX {
            s = prev[s] as usize;
            path.push(s / 5);
        }
        path.reverse();
        for &k in &path {
            self.taken[k] = true;
        }
        // The corners: the ends and where the direction changes.
        let mut corners: Vec<DVec3> = Vec::new();
        for (m, &k) in path.iter().enumerate() {
            let turn = m > 0 && m + 1 < path.len() && {
                let (p, q) = (path[m - 1], path[m + 1]);
                (k as i64 - p as i64) != (q as i64 - k as i64)
            };
            if m == 0 || m + 1 == path.len() || turn {
                corners.push(self.point(k));
            }
        }
        let stub = |end: End, p: DVec3| match end {
            End::Border => None,
            End::Electrode(e) | End::Inside(e) => Some(surface_point(&self.outlines[e], p)),
        };
        let mut wire = Vec::with_capacity(corners.len() + 2);
        if let Some(s) = corners.first().and_then(|p| stub(a, *p)) {
            wire.push(s);
        }
        wire.extend(&corners);
        if let Some(s) = corners.last().and_then(|p| stub(b, *p)) {
            wire.push(s);
        }
        wire.dedup_by(|p, q| (*p - *q).length() < 1e-9);
        Some(wire.into_iter().map(|p| (p.x, p.y)).collect())
    }
}

impl Level {
    /// The tube level's wiring for `player`'s placement (`LeadPlan`): its electrodes' leads and
    /// its floating circuits' straps, routed among the electrodes (the player's plates
    /// included) and the metals.
    #[must_use]
    pub fn tube_leads(&self, player: &[Element]) -> LeadPlan {
        let Some(spec) = self.tube.as_ref() else {
            return LeadPlan::default();
        };
        let metals: Vec<(DVec3, DVec3)> = spec.metals.iter().map(|m| self.metal_rect(m)).collect();
        let outlines = self.electrode_outlines(player);
        let n = self.electrodes.len().min(outlines.len());
        let open = spec.output.is_some();
        let (parent, envelope) = topology(&outlines[..n], open);
        let b = self.bounds();
        let mut router = Router::new(b.min, b.max, &outlines, &metals);
        let straps = spec
            .circuits
            .iter()
            .filter(|c| c.resistance.is_none() && c.electrode < outlines.len())
            .flat_map(|c| {
                c.linked
                    .iter()
                    .filter(|&&e| e < outlines.len())
                    .map(move |&e| (c.electrode, e))
            })
            .collect::<Vec<_>>()
            .into_iter()
            .filter_map(|(a, e)| router.route(End::Electrode(a), End::Electrode(e)))
            .collect();
        // The leads, the shortest first: an electrode's distance to its parent, or to the
        // board's nearest edge.
        let to_border = |o: &Outline| {
            let c = centre(o);
            (c.x - b.min.x)
                .min(b.max.x - c.x)
                .min(c.y - b.min.y)
                .min(b.max.y - c.y)
        };
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&e| !(parent[e].is_none() && open && envelope[e]))
            .map(|e| {
                let d = parent[e].map_or_else(
                    || to_border(&outlines[e]),
                    |p| outlines[p].distance(&outlines[e]),
                );
                (d, e)
            })
            .collect();
        order.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let mut leads = vec![None; n];
        for (_, e) in order {
            // From inside the ring that encloses it, or from the outside of the envelope
            // nearest it (an open level's electrode outside every ring: Magnetron's antenna).
            let from = match parent[e] {
                Some(p) => match outlines[p] {
                    Outline::Ring { center, inner, .. }
                        if farthest(&outlines[e], center) < inner =>
                    {
                        End::Inside(p)
                    }
                    _ => End::Electrode(p),
                },
                None => End::Border,
            };
            leads[e] = router.route(from, End::Electrode(e));
        }
        // With the carriers' own field, the border's leads return to the ground terminal along
        // the grounded envelope's inner face, `ENVELOPE_MARGIN` beyond the board (DEC-27).
        let m = crate::tube::ENVELOPE_MARGIN;
        let (lo, hi) = (
            DVec3::new(b.min.x - m, b.min.y - m, 0.0),
            DVec3::new(b.max.x + m, b.max.y + m, 0.0),
        );
        let ground = spec
            .speed_of_light
            .map(|_| DVec3::new(0.5 * (lo.x + hi.x), lo.y, 0.0));
        let returns = (0..n)
            .map(|e| {
                let g = ground?;
                let start = leads[e].as_ref()?.first()?;
                parent[e].is_none().then(|| {
                    border_path(lo, hi, g, DVec3::new(start.0, start.1, 0.0))
                        .into_iter()
                        .map(|p| (p.x, p.y))
                        .collect()
                })
            })
            .collect();
        LeadPlan {
            parent,
            envelope,
            leads,
            straps,
            ground: ground.map(|g| (g.x, g.y)),
            returns,
        }
    }
}

impl LeadPlan {
    /// Per level electrode, with the carriers' own field: the whole path its circuit current
    /// comes along from the ground terminal, its return along the border then its lead (the
    /// field's sources, `physics::tube::Tube::leads`); None for an electrode without a lead
    /// from the border.
    #[must_use]
    pub fn circuit_paths(&self) -> Vec<Option<Vec<DVec3>>> {
        self.leads
            .iter()
            .zip(&self.returns)
            .map(|(lead, back)| {
                let (lead, back) = (lead.as_ref()?, back.as_ref()?);
                let mut path: Vec<DVec3> = back
                    .iter()
                    .chain(lead.iter().skip(1))
                    .map(|&(x, y)| DVec3::new(x, y, 0.0))
                    .collect();
                path.dedup_by(|p, q| (*p - *q).length() < 1e-9);
                Some(path)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests;
