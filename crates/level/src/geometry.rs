//! The board's geometry: the electrodes' outlines and the rectangle tests the
//! placement checks share.

use super::*;

/// An electrode's cross-section for geometric checks: a box's corners, or (tube levels) a
/// rod's circle or a ring's annulus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outline {
    Rect([DVec3; 4]),
    Circle {
        center: DVec3,
        radius: f64,
    },
    Ring {
        center: DVec3,
        inner: f64,
        outer: f64,
    },
}

impl Outline {
    /// The nearest and farthest distances from `c` to this outline's material (the nearest
    /// 0 if `c` lies in it).
    fn radial_range(&self, c: DVec3) -> (f64, f64) {
        match *self {
            Outline::Rect(r) => {
                let near = if rect_contains(&r, c) {
                    0.0
                } else {
                    (0..4)
                        .map(|i| segment_distance(c, r[i], r[(i + 1) % 4]))
                        .fold(f64::INFINITY, f64::min)
                };
                let far = r.iter().map(|p| (*p - c).length()).fold(0.0, f64::max);
                (near, far)
            }
            Outline::Circle { center, radius } => {
                let d = (center - c).length();
                ((d - radius).max(0.0), d + radius)
            }
            Outline::Ring {
                center,
                inner,
                outer,
            } => {
                let d = (center - c).length();
                let near = if d < inner {
                    inner - d
                } else if d > outer {
                    d - outer
                } else {
                    0.0
                };
                (near, d + outer)
            }
        }
    }

    /// The distance between two outlines' material (0 where they touch or overlap). Exact
    /// for connected shapes: a shape reaching from a ring's hole to beyond it crosses the
    /// ring.
    pub fn distance(&self, other: &Outline) -> f64 {
        match (*self, *other) {
            (Outline::Rect(a), Outline::Rect(b)) => rect_distance(&a, &b),
            (Outline::Circle { center, radius }, o) | (o, Outline::Circle { center, radius }) => {
                (o.radial_range(center).0 - radius).max(0.0)
            }
            (
                Outline::Ring {
                    center,
                    inner,
                    outer,
                },
                o,
            )
            | (
                o,
                Outline::Ring {
                    center,
                    inner,
                    outer,
                },
            ) => {
                let (near, far) = o.radial_range(center);
                if far < inner {
                    inner - far
                } else if near > outer {
                    near - outer
                } else {
                    0.0
                }
            }
        }
    }

    /// The distance from point `p` to the material (0 inside it).
    pub fn distance_to(&self, p: DVec3) -> f64 {
        self.radial_range(p).0
    }

    /// Whether the outline lies within the axis-aligned rectangle `r` (corners as
    /// `Level::node_rect` gives them: lower left first, upper right third).
    pub(crate) fn within(&self, r: &[DVec3; 4]) -> bool {
        let (lo, hi) = (r[0], r[2]);
        let inside = |c: DVec3, e: f64| {
            c.x - e >= lo.x - 1e-9
                && c.y - e >= lo.y - 1e-9
                && c.x + e <= hi.x + 1e-9
                && c.y + e <= hi.y + 1e-9
        };
        match *self {
            Outline::Rect(corners) => corners.iter().all(|c| inside(*c, 0.0)),
            Outline::Circle { center, radius } => inside(center, radius),
            Outline::Ring { center, outer, .. } => inside(center, outer),
        }
    }
}

/// In-plane corners of an electrode box (z = 0), counter-clockwise.
pub(super) fn rect_corners(b: &physics::bem::BoxElectrode) -> [DVec3; 4] {
    let u = DVec3::new(libm::cos(b.angle), libm::sin(b.angle), 0.0) * b.half_length;
    let v = DVec3::new(-libm::sin(b.angle), libm::cos(b.angle), 0.0) * b.half_thickness;
    let c = DVec3::new(b.center.x, b.center.y, 0.0);
    [c - u - v, c + u - v, c + u + v, c - u + v]
}

/// Whether a point lies in a convex counter-clockwise polygon (boundary included).
pub(super) fn rect_contains(r: &[DVec3; 4], p: DVec3) -> bool {
    (0..4).all(|i| {
        let (a, b) = (r[i], r[(i + 1) % 4]);
        (b - a).x * (p - a).y - (b - a).y * (p - a).x >= -1e-12
    })
}

pub(super) fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

pub(super) fn segments_cross(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> bool {
    let cross = |o: DVec3, p: DVec3, q: DVec3| (p - o).x * (q - o).y - (p - o).y * (q - o).x;
    let (d1, d2) = (cross(c, d, a), cross(c, d, b));
    let (d3, d4) = (cross(a, b, c), cross(a, b, d));
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

/// Distance between two convex quadrilaterals in the plane (0 if they overlap).
pub(super) fn rect_distance(a: &[DVec3; 4], b: &[DVec3; 4]) -> f64 {
    let overlap = a.iter().any(|p| rect_contains(b, *p))
        || b.iter().any(|p| rect_contains(a, *p))
        || (0..4)
            .any(|i| (0..4).any(|j| segments_cross(a[i], a[(i + 1) % 4], b[j], b[(j + 1) % 4])));
    if overlap {
        return 0.0;
    }
    let one_way = |p: &[DVec3; 4], q: &[DVec3; 4]| {
        p.iter()
            .flat_map(|x| (0..4).map(move |j| segment_distance(*x, q[j], q[(j + 1) % 4])))
            .fold(f64::INFINITY, f64::min)
    };
    one_way(a, b).min(one_way(b, a))
}

impl Level {
    /// The outlines of every electrode of a placement: the level's (by their shapes), then
    /// the player's plates (boxes), in `all_box_electrodes`' order.
    pub fn electrode_outlines(&self, player: &[Element]) -> Vec<Outline> {
        self.electrodes
            .iter()
            .map(|e| {
                let center = self.grid.position(e.center);
                let center = DVec3::new(center.x, center.y, 0.0);
                match e.shape {
                    ElectrodeShape::Rod => Outline::Circle {
                        center,
                        radius: e.length / 2.0,
                    },
                    ElectrodeShape::Ring => Outline::Ring {
                        center,
                        inner: (e.length - e.thickness) / 2.0,
                        outer: (e.length + e.thickness) / 2.0,
                    },
                    ElectrodeShape::Box => {
                        Outline::Rect(rect_corners(&physics::bem::BoxElectrode {
                            center,
                            angle: e.angle_deg.to_radians(),
                            half_length: e.length / 2.0,
                            half_thickness: e.thickness / 2.0,
                            half_height: e.height / 2.0,
                            bias: Bias::Grounded,
                        }))
                    }
                }
            })
            .chain(
                player
                    .iter()
                    .filter(|e| e.kind == ElementKind::Plate)
                    .map(|e| Outline::Rect(rect_corners(&self.plate_box(e)))),
            )
            .collect()
    }
}
