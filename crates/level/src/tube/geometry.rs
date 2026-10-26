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

    /// The direction from contact `e` into metal `m`, if `e` is an axis-aligned box whose
    /// face lies along one side of the metal, within it (it then gives off carriers into
    /// the metal only).
    #[must_use]
    pub fn contact_side(&self, m: &MetalSpec, e: usize) -> Option<DVec3> {
        let (b0, b1) = self.electrode_box(e)?;
        let (lo, hi) = self.metal_rect(m);
        let tol = 1e-9 * (1.0 + hi.x - lo.x + hi.y - lo.y);
        let within = |a0: f64, a1: f64, s0: f64, s1: f64| a0 >= s0 - tol && a1 <= s1 + tol;
        if (b1.x - lo.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            Some(DVec3::X)
        } else if (b0.x - hi.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            Some(-DVec3::X)
        } else if (b1.y - lo.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            Some(DVec3::Y)
        } else if (b0.y - hi.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            Some(-DVec3::Y)
        } else {
            None
        }
    }

    /// Metal `m`'s contacts' faces along its sides (`contact_side`), each a segment from one
    /// end to the other within the metal's side (`MetalView::contacts`).
    #[must_use]
    pub fn contact_faces(&self, m: &MetalSpec) -> Vec<[[f64; 2]; 2]> {
        let (lo, hi) = self.metal_rect(m);
        m.contacts
            .iter()
            .filter_map(|&e| {
                let side = self.contact_side(m, e)?;
                let (b0, b1) = self.electrode_box(e)?;
                let (y0, y1) = (b0.y.max(lo.y), b1.y.min(hi.y));
                let (x0, x1) = (b0.x.max(lo.x), b1.x.min(hi.x));
                Some(if side == DVec3::X {
                    [[lo.x, y0], [lo.x, y1]]
                } else if side == -DVec3::X {
                    [[hi.x, y0], [hi.x, y1]]
                } else if side == DVec3::Y {
                    [[x0, lo.y], [x1, lo.y]]
                } else {
                    [[x0, hi.y], [x1, hi.y]]
                })
            })
            .collect()
    }

    /// The depth of the layer a metal's contacts keep neutral: half a mean free path, at
    /// most a quarter of the metal's width.
    #[must_use]
    pub fn contact_depth(&self, m: &MetalSpec) -> f64 {
        let (lo, hi) = self.metal_rect(m);
        (0.5 * m.mean_free_path).min(0.25 * (hi.x - lo.x).min(hi.y - lo.y))
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
    /// its middle inside.
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
            within(a) && within(b) && m.x > lo.x && m.x < hi.x && m.y > lo.y && m.y < hi.y
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
