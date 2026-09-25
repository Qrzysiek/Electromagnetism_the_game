//! Visual-only computations: field lines. Nothing here feeds back into gameplay.

use std::collections::{HashMap, VecDeque};

use physics::DVec3;
use physics::field::{Coulomb, FieldSolver};
use physics::trajectory::Scenario;

/// A traced field line with arrowheads along it.
pub struct FieldLine {
    pub points: Vec<DVec3>,
    /// Arrowhead positions and unit directions of **E** there.
    pub arrows: Vec<(DVec3, DVec3)>,
}

/// Largest direction change allowed per step (radians): keeps curves smooth.
const MAX_TURN: f64 = 0.02;

/// Spatial hash of the points of accepted lines, for "is any line within d?" queries.
struct Grid {
    cell: f64,
    map: HashMap<(i64, i64), Vec<(DVec3, usize)>>,
}

impl Grid {
    fn new(cell: f64) -> Self {
        Self {
            cell,
            map: HashMap::new(),
        }
    }

    #[allow(clippy::cast_possible_truncation)]
    fn key(&self, p: DVec3) -> (i64, i64) {
        (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
        )
    }

    fn insert(&mut self, p: DVec3, line: usize) {
        let k = self.key(p);
        self.map.entry(k).or_default().push((p, line));
    }

    /// Whether a point of a line other than `except` lies within `d` (`d ≤ cell`).
    fn near(&self, p: DVec3, d: f64, except: Option<usize>) -> bool {
        let (kx, ky) = self.key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(v) = self.map.get(&(kx + dx, ky + dy))
                    && v.iter()
                        .any(|&(q, l)| Some(l) != except && (q - p).length() < d)
                {
                    return true;
                }
            }
        }
        false
    }
}

/// Field lines in the plane z = 0 that start and end only on charges or at the arena
/// edge, as field lines of an electrostatic field do. Spacing is kept roughly even with
/// a variant of Jobard & Lefer's evenly spaced streamlines (1997): candidate lines are
/// seeded near charges, then offset by `spacing` perpendicular to accepted lines, then
/// from a grid sweep; each candidate is traced in full both ways, and rejected if, away
/// from the charges (where lines must converge), a noticeable part of it runs closer
/// than `spacing / 2` to an accepted line. Between close charges this leaves few lines.
///
/// In the 2D slice of a 3D field, line density does not represent field strength
/// (SPEC §4): lines show direction.
pub fn field_lines(scn: &Scenario<Coulomb>, spacing: f64) -> Vec<FieldLine> {
    if scn.obstacles.is_empty() {
        return Vec::new();
    }
    let b = scn.bounds.expect("bounds");
    let d_sep = spacing.max(0.2);
    let d_test = 0.5 * d_sep;
    // Near a charge lines converge; proximity there is not counted.
    let converge_zone = 1.5 * d_sep;
    let arrow_spacing = (2.5 * d_sep).max(2.0);

    let field_dir = |x: DVec3, sign: f64| -> Option<DVec3> {
        let e = scn.field.sample(x, 0.0).e;
        (e.length() > 0.0 && e.is_finite()).then(|| e.normalize() * sign)
    };
    let inside_bounds = |x: DVec3| x.x > b.min.x && x.y > b.min.y && x.x < b.max.x && x.y < b.max.y;
    let charge_distance = |x: DVec3| {
        scn.obstacles
            .iter()
            .map(|s| s.signed_distance(x))
            .fold(f64::INFINITY, f64::min)
    };

    // Traces from `seed` in direction `sign` (+1 along E) until the line reaches a charge
    // or the edge (or a point where E vanishes).
    let trace = |seed: DVec3, sign: f64| -> Vec<DVec3> {
        let mut pts = vec![seed];
        let mut x = seed;
        let mut ds: f64 = 0.02;
        let mut length = 0.0;
        'outer: for _ in 0..50_000 {
            let Some(d0) = field_dir(x, sign) else { break };
            let xn = loop {
                let k1 = d0;
                let Some(k2) = field_dir(x + k1 * (0.5 * ds), sign) else {
                    break 'outer;
                };
                let Some(k3) = field_dir(x + k2 * (0.5 * ds), sign) else {
                    break 'outer;
                };
                let Some(k4) = field_dir(x + k3 * ds, sign) else {
                    break 'outer;
                };
                let xn = x + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (ds / 6.0);
                let Some(dn) = field_dir(xn, sign) else {
                    break 'outer;
                };
                let turn = d0.angle_between(dn);
                if turn > MAX_TURN && ds > 1e-3 {
                    ds *= 0.5;
                    continue;
                }
                if turn < 0.3 * MAX_TURN {
                    ds = (ds * 1.5).min(0.25 * d_sep);
                }
                break xn;
            };
            length += (xn - x).length();
            x = xn;
            pts.push(x);
            if !inside_bounds(x) || charge_distance(x) < 0.0 || length > 2000.0 {
                break;
            }
        }
        pts
    };

    let mut grid = Grid::new(d_sep);
    let mut lines: Vec<FieldLine> = Vec::new();
    let mut queue: VecDeque<DVec3> = VecDeque::new();
    for s in &scn.obstacles {
        for i in 0..8 {
            let a = std::f64::consts::FRAC_PI_4 * f64::from(i) + 0.2;
            queue
                .push_back(s.center + DVec3::new(a.cos(), a.sin(), 0.0) * (s.radius + 0.6 * d_sep));
        }
    }
    let mut sweep: Vec<DVec3> = Vec::new();
    let mut y = b.min.y + 0.5 * d_sep;
    while y < b.max.y {
        let mut x = b.min.x + 0.5 * d_sep;
        while x < b.max.x {
            sweep.push(DVec3::new(x, y, 0.0));
            x += d_sep;
        }
        y += d_sep;
    }
    let mut sweep = sweep.into_iter();

    while let Some(seed) = queue.pop_front().or_else(|| sweep.next()) {
        if !inside_bounds(seed) || charge_distance(seed) < 0.0 {
            continue;
        }
        if charge_distance(seed) > converge_zone && grid.near(seed, d_sep, None) {
            continue;
        }
        let forward = trace(seed, 1.0);
        let backward = trace(seed, -1.0);
        let mut points: Vec<DVec3> = backward.into_iter().rev().collect();
        points.extend(forward.into_iter().skip(1));

        // Reject lines that crowd accepted ones away from the charges.
        let mut crowded = 0.0;
        let mut free = 0.0;
        for w in points.windows(2) {
            let seg = (w[1] - w[0]).length();
            if charge_distance(w[1]) > converge_zone {
                if grid.near(w[1], d_test, None) {
                    crowded += seg;
                } else {
                    free += seg;
                }
            }
        }
        if crowded > 0.1 * (crowded + free) || crowded > d_sep {
            continue;
        }

        let id = lines.len();
        let mut since_store = 0.0;
        let mut since_seed = 0.0;
        for w in points.windows(2) {
            let seg = (w[1] - w[0]).length();
            since_store += seg;
            since_seed += seg;
            if since_store >= 0.25 * d_sep {
                grid.insert(w[1], id);
                since_store = 0.0;
            }
            if since_seed >= d_sep {
                since_seed = 0.0;
                let t = (w[1] - w[0]).normalize_or_zero();
                let n = DVec3::new(-t.y, t.x, 0.0);
                queue.push_back(w[1] + n * d_sep);
                queue.push_back(w[1] - n * d_sep);
            }
        }
        let arrows = arrows_along(&points, arrow_spacing, |x| field_dir(x, 1.0));
        lines.push(FieldLine { points, arrows });
    }
    lines
}

/// Arrowheads every `spacing` of arc length (the first at half the spacing), each
/// pointing along E at its position.
fn arrows_along(
    points: &[DVec3],
    spacing: f64,
    e_dir: impl Fn(DVec3) -> Option<DVec3>,
) -> Vec<(DVec3, DVec3)> {
    let mut out = Vec::new();
    let mut next = 0.5 * spacing;
    let mut s = 0.0;
    for w in points.windows(2) {
        let seg = (w[1] - w[0]).length();
        while seg > 0.0 && s + seg >= next {
            let p = w[0].lerp(w[1], (next - s) / seg);
            if let Some(d) = e_dir(p) {
                out.push((p, d));
            }
            next += spacing;
        }
        s += seg;
    }
    out
}
