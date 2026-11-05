//! A tube level's board: its electrodes as prisms, its metals, their contacts and
//! the player's cuts.

use super::*;

impl Level {
    /// The prisms of a placement: every box electrode (the level's at the player's
    /// supplies, then the player's plates), its cross-section.
    pub fn tube_prisms(&self, player: &[Element]) -> Vec<Prism> {
        let mut prisms: Vec<Prism> = self
            .all_box_electrodes(player)
            .into_iter()
            .enumerate()
            .map(|(i, b)| {
                let shape = self
                    .electrodes
                    .get(i)
                    .map_or(crate::ElectrodeShape::Box, |e| e.shape);
                let e = self.electrodes.get(i);
                let section = match (shape, e) {
                    (crate::ElectrodeShape::Rod, Some(e)) => Section::Circle {
                        center: b.center,
                        radius: e.length / 2.0,
                    },
                    (crate::ElectrodeShape::Ring, Some(e)) => Section::Ring {
                        center: b.center,
                        radius: e.length / 2.0,
                        thickness: e.thickness,
                    },
                    _ => Section::Rect {
                        center: b.center,
                        angle: b.angle,
                        half_length: b.half_length,
                        half_thickness: b.half_thickness,
                    },
                };
                Prism {
                    section,
                    bias: b.bias,
                }
            })
            .collect();
        // An electrode on a source (ideal or through a resistance) starts at the source's
        // value at t = 0: the supply was on before the cathode heated, so its capacitors
        // are charged. A floating one keeps its own bias.
        for c in self.tube_circuits(player) {
            if c.resistance.is_finite()
                && let Some(p) = prisms.get_mut(c.electrode)
            {
                p.bias = physics::conductor::Bias::Potential(c.source.value(0.0));
            }
        }
        // The envelope, grounded, last (ENVELOPE_MARGIN).
        let max = self.grid.position(self.grid.max_node());
        prisms.push(Prism {
            section: Section::Frame {
                center: DVec3::new(0.5 * max.x, 0.5 * max.y, 0.0),
                half_width: 0.5 * max.x + ENVELOPE_MARGIN,
                half_height: 0.5 * max.y + ENVELOPE_MARGIN,
                segment: ENVELOPE_SEGMENT,
            },
            bias: physics::conductor::Bias::Grounded,
        });
        prisms
    }

    /// A metal's rectangle (lower and upper corners).
    #[must_use]
    pub fn metal_rect(&self, m: &MetalSpec) -> (DVec3, DVec3) {
        let (a, b) = (
            self.grid.position(m.corners[0]),
            self.grid.position(m.corners[1]),
        );
        (
            DVec3::new(a.x.min(b.x), a.y.min(b.y), 0.0),
            DVec3::new(a.x.max(b.x), a.y.max(b.y), 0.0),
        )
    }

    /// A ring's hole (lower and upper corners; `MetalSpec::hole`), None for a solid metal.
    #[must_use]
    pub fn metal_hole(&self, m: &MetalSpec) -> Option<(DVec3, DVec3)> {
        let [a, b] = m.hole?;
        let (a, b) = (self.grid.position(a), self.grid.position(b));
        Some((
            DVec3::new(a.x.min(b.x), a.y.min(b.y), 0.0),
            DVec3::new(a.x.max(b.x), a.y.max(b.y), 0.0),
        ))
    }

    /// A metal's narrowest width: its rectangle's shorter side, a ring's thinnest wall.
    #[must_use]
    pub fn metal_width(&self, m: &MetalSpec) -> f64 {
        let (lo, hi) = self.metal_rect(m);
        let whole = (hi.x - lo.x).min(hi.y - lo.y);
        self.metal_hole(m).map_or(whole, |(a, b)| {
            (a.x - lo.x).min(hi.x - b.x).min(a.y - lo.y).min(hi.y - b.y)
        })
    }

    /// A metal's conductor's area: its rectangle's, less its hole's.
    #[must_use]
    pub fn metal_area(&self, m: &MetalSpec) -> f64 {
        let (lo, hi) = self.metal_rect(m);
        let whole = (hi.x - lo.x) * (hi.y - lo.y);
        self.metal_hole(m)
            .map_or(whole, |(a, b)| whole - (b.x - a.x) * (b.y - a.y))
    }

    /// Whether `x` lies in a metal's conductor: in its rectangle (its edges included), not
    /// inside its hole.
    #[must_use]
    pub fn metal_holds(&self, m: &MetalSpec, x: DVec3) -> bool {
        let (lo, hi) = self.metal_rect(m);
        let in_rect = x.x >= lo.x && x.x <= hi.x && x.y >= lo.y && x.y <= hi.y;
        in_rect
            && self
                .metal_hole(m)
                .is_none_or(|(a, b)| !(x.x > a.x && x.x < b.x && x.y > a.y && x.y < b.y))
    }

    /// An axis-aligned box electrode's corners (None for another shape or angle).
    pub(super) fn electrode_box(&self, e: usize) -> Option<(DVec3, DVec3)> {
        let el = self.electrodes.get(e)?;
        if el.shape != crate::ElectrodeShape::Box {
            return None;
        }
        let quarter = libm::round(el.angle_deg / 90.0);
        if (el.angle_deg - 90.0 * quarter).abs() > 1e-9 {
            return None;
        }
        #[allow(clippy::cast_possible_truncation)]
        let along_x = (quarter as i64).rem_euclid(2) == 0;
        let (hx, hy) = if along_x {
            (0.5 * el.length, 0.5 * el.thickness)
        } else {
            (0.5 * el.thickness, 0.5 * el.length)
        };
        let c = self.grid.position(el.center);
        Some((
            DVec3::new(c.x - hx, c.y - hy, 0.0),
            DVec3::new(c.x + hx, c.y + hy, 0.0),
        ))
    }

    /// How electrode `e` meets metal `m` as a contact, if it does: an axis-aligned box whose
    /// face lies along one side of the metal, within it, outside (`Toward` the metal: it
    /// gives off carriers that way only); or one in a ring's hole, flush with one of the
    /// hole's walls (`Toward` that wall) or filling the hole (`Around`: every face is
    /// flush with a wall and gives off carriers; a Corbino disc's inner contact). A box in
    /// the hole touching two or three walls without filling it is none (its other faces
    /// would give off carriers into the hole's vacuum).
    #[must_use]
    pub fn contact_side(&self, m: &MetalSpec, e: usize) -> Option<ContactSide> {
        let (b0, b1) = self.electrode_box(e)?;
        let (lo, hi) = self.metal_rect(m);
        let tol = 1e-9 * (1.0 + hi.x - lo.x + hi.y - lo.y);
        let within = |a0: f64, a1: f64, s0: f64, s1: f64| a0 >= s0 - tol && a1 <= s1 + tol;
        if (b1.x - lo.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            return Some(ContactSide::Toward(DVec3::X));
        } else if (b0.x - hi.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            return Some(ContactSide::Toward(-DVec3::X));
        } else if (b1.y - lo.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            return Some(ContactSide::Toward(DVec3::Y));
        } else if (b0.y - hi.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            return Some(ContactSide::Toward(-DVec3::Y));
        }
        let (h0, h1) = self.metal_hole(m)?;
        if !(within(b0.x, b1.x, h0.x, h1.x) && within(b0.y, b1.y, h0.y, h1.y)) {
            return None;
        }
        // The walls it is flush with, each with the direction into it.
        let flush: Vec<DVec3> = [
            ((b0.x - h0.x).abs() < tol, -DVec3::X),
            ((b1.x - h1.x).abs() < tol, DVec3::X),
            ((b0.y - h0.y).abs() < tol, -DVec3::Y),
            ((b1.y - h1.y).abs() < tol, DVec3::Y),
        ]
        .into_iter()
        .filter_map(|(f, d)| f.then_some(d))
        .collect();
        match flush.as_slice() {
            [d] => Some(ContactSide::Toward(*d)),
            [_, _, _, _] => Some(ContactSide::Around),
            _ => None,
        }
    }

    /// Contact `e`'s faces on metal `m` (`contact_side`), each a segment from one end to the
    /// other within the metal's side or the hole's wall: one for a contact along a side,
    /// four for one filling a ring's hole (the hole's walls); none if it is not a contact.
    #[must_use]
    pub fn contact_segments(&self, m: &MetalSpec, e: usize) -> Vec<[[f64; 2]; 2]> {
        let (Some(side), Some((b0, b1))) = (self.contact_side(m, e), self.electrode_box(e)) else {
            return Vec::new();
        };
        let (lo, hi) = self.metal_rect(m);
        let (y0, y1) = (b0.y.max(lo.y), b1.y.min(hi.y));
        let (x0, x1) = (b0.x.max(lo.x), b1.x.min(hi.x));
        // In a hole the box's own face is the wall's piece it covers; outside, the face it
        // presses against the metal's side.
        let inner = self.metal_hole(m).is_some_and(|(h0, h1)| {
            b0.x >= h0.x - 1e-9 && b1.x <= h1.x + 1e-9 && b0.y >= h0.y - 1e-9 && b1.y <= h1.y + 1e-9
        });
        let face = |d: DVec3| -> [[f64; 2]; 2] {
            if inner {
                if d == DVec3::X {
                    [[b1.x, b0.y], [b1.x, b1.y]]
                } else if d == -DVec3::X {
                    [[b0.x, b0.y], [b0.x, b1.y]]
                } else if d == DVec3::Y {
                    [[b0.x, b1.y], [b1.x, b1.y]]
                } else {
                    [[b0.x, b0.y], [b1.x, b0.y]]
                }
            } else if d == DVec3::X {
                [[lo.x, y0], [lo.x, y1]]
            } else if d == -DVec3::X {
                [[hi.x, y0], [hi.x, y1]]
            } else if d == DVec3::Y {
                [[x0, lo.y], [x1, lo.y]]
            } else {
                [[x0, hi.y], [x1, hi.y]]
            }
        };
        match side {
            ContactSide::Toward(d) => vec![face(d)],
            ContactSide::Around => [-DVec3::Y, DVec3::X, DVec3::Y, -DVec3::X]
                .into_iter()
                .map(face)
                .collect(),
        }
    }

    /// Metal `m`'s contacts' faces, all of them (`contact_segments`, in its contacts' order;
    /// `MetalView::contacts`).
    #[must_use]
    pub fn contact_faces(&self, m: &MetalSpec) -> Vec<[[f64; 2]; 2]> {
        m.contacts
            .iter()
            .flat_map(|&e| self.contact_segments(m, e))
            .collect()
    }

    /// Whether electrode `e` lies wholly in metal `m`'s hole (its box, or the square about a
    /// round one's centre: as `metal_issues` bounds it).
    #[must_use]
    pub fn in_metal_hole(&self, m: &MetalSpec, e: usize) -> bool {
        let Some((h0, h1)) = self.metal_hole(m) else {
            return false;
        };
        let (b0, b1) = self.electrode_bounds(e);
        b0.x >= h0.x - 1e-9 && b1.x <= h1.x + 1e-9 && b0.y >= h0.y - 1e-9 && b1.y <= h1.y + 1e-9
    }

    /// An electrode's bounding box in the plane: an axis-aligned box's corners, else the
    /// square about its centre with half the sum of its length and thickness as half its side
    /// (a turned box's or a round electrode's reach).
    pub(super) fn electrode_bounds(&self, e: usize) -> (DVec3, DVec3) {
        self.electrode_box(e).unwrap_or_else(|| {
            let c = self.grid.position(self.electrodes[e].center);
            let r = 0.5 * (self.electrodes[e].length + self.electrodes[e].thickness);
            (c - DVec3::new(r, r, 0.0), c + DVec3::new(r, r, 0.0))
        })
    }

    /// The depth of the layer a metal's contacts keep neutral: half a mean free path, at
    /// most a quarter of the metal's width (a ring's thinnest wall).
    #[must_use]
    pub fn contact_depth(&self, m: &MetalSpec) -> f64 {
        (0.5 * m.mean_free_path).min(0.25 * self.metal_width(m))
    }

    /// A cut's slit (`ElementKind::Cut`): from its node along `angle_deg` for `value` (its
    /// length); None if it has no length.
    #[must_use]
    pub fn cut_segment(&self, e: &Element) -> Option<(DVec3, DVec3)> {
        if e.kind != ElementKind::Cut || !(e.value.is_finite() && e.value > 0.0) {
            return None;
        }
        if !e.angle_deg.is_finite() {
            return None;
        }
        let a = self.grid.position(e.node);
        let angle = e.angle_deg.to_radians();
        let d = DVec3::new(libm::cos(angle), libm::sin(angle), 0.0);
        Some((
            DVec3::new(a.x, a.y, 0.0),
            DVec3::new(a.x, a.y, 0.0) + d * e.value,
        ))
    }

    /// The metal a cut lies in: both its ends in the metal's rectangle (on its edges too),
    /// its middle inside, in the conductor (a ring's: not in its hole, where another metal
    /// may lie).
    #[must_use]
    pub fn cut_metal(&self, e: &Element) -> Option<usize> {
        let spec = self.tube.as_ref()?;
        let (a, b) = self.cut_segment(e)?;
        let m = (a + b) * 0.5;
        spec.metals.iter().position(|metal| {
            let (lo, hi) = self.metal_rect(metal);
            let tol = 1e-9 * (1.0 + hi.x.abs() + hi.y.abs());
            let within = |p: DVec3| {
                p.x >= lo.x - tol && p.x <= hi.x + tol && p.y >= lo.y - tol && p.y <= hi.y + tol
            };
            within(a)
                && within(b)
                && m.x > lo.x
                && m.x < hi.x
                && m.y > lo.y
                && m.y < hi.y
                && self.metal_holds(metal, m)
        })
    }

    /// A cut from node `a` to node `b` (the player draws it from one to the other).
    #[must_use]
    pub fn cut_between(&self, a: Node, b: Node) -> Element {
        let (pa, pb) = (self.grid.position(a), self.grid.position(b));
        let d = pb - pa;
        Element {
            node: a,
            kind: ElementKind::Cut,
            value: libm::sqrt(d.x * d.x + d.y * d.y),
            angle_deg: libm::atan2(d.y, d.x).to_degrees(),
            omega: None,
            speed: None,
        }
    }

    /// The player's cuts: each with a length, lying in one of the metals (its ends in the
    /// metal, on its edges too).
    pub(crate) fn check_cuts(&self, player: &[Element]) -> Result<(), crate::PlacementError> {
        for e in player.iter().filter(|e| e.kind == ElementKind::Cut) {
            if self.cut_metal(e).is_none() {
                return Err(crate::PlacementError::CutNotInMetal(e.node));
            }
        }
        Ok(())
    }
}
