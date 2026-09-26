//! Editor state: player elements (charges, magnets, antennas) on the grid, cursor,
//! selected element.
//! Independent of the rendering engine so it can be unit tested.

use level::{ANTENNA_ANGLES, Element, ElementKind, Level, Node, PlacementError};

/// Kinds in palette order.
pub const KINDS: [ElementKind; 3] = [
    ElementKind::Charge,
    ElementKind::Magnet,
    ElementKind::Antenna,
];

pub struct Editor {
    /// The level as loaded (recommended grid).
    base: Level,
    /// The level at the current refinement.
    pub level: Level,
    pub placement: Vec<Element>,
    pub cursor: Node,
    /// Kind of new elements.
    pub kind: ElementKind,
    /// Index into the allowed values of `kind` for new elements.
    pub magnitude_index: usize,
    /// Sign of new elements (charges: sign of Q; magnets: moment along +z; antennas:
    /// phase).
    pub positive: bool,
    /// Orientation of new antennas, degrees (one of `ANTENNA_ANGLES`).
    pub angle_deg: f64,
    /// Last rejected action, for the status line.
    pub message: Option<String>,
    /// Incremented on every change of the physical setup.
    pub revision: u64,
}

/// Allowed magnitudes of a kind in a level.
pub fn magnitudes(level: &Level, kind: ElementKind) -> &[f64] {
    match kind {
        ElementKind::Charge => &level.limits.magnitudes,
        ElementKind::Magnet => &level.limits.magnet_strengths,
        ElementKind::Antenna => &level.limits.antenna_amplitudes,
    }
}

/// Maximum number of player elements of a kind.
pub fn max_of(level: &Level, kind: ElementKind) -> u32 {
    match kind {
        ElementKind::Charge => level.limits.max_charges,
        ElementKind::Magnet => level.limits.max_magnets,
        ElementKind::Antenna => level.limits.max_antennas,
    }
}

/// Whether both signs of a kind are available.
fn both_signs(level: &Level, kind: ElementKind) -> bool {
    kind != ElementKind::Charge || (level.limits.allow_positive && level.limits.allow_negative)
}

fn default_positive(level: &Level, kind: ElementKind) -> bool {
    kind != ElementKind::Charge || level.limits.allow_positive
}

impl Editor {
    pub fn new(level: Level) -> Self {
        let m = level.grid.max_node();
        let kind = KINDS
            .into_iter()
            .find(|&k| max_of(&level, k) > 0)
            .unwrap_or(ElementKind::Charge);
        let positive = default_positive(&level, kind);
        Self {
            base: level.clone(),
            level,
            placement: Vec::new(),
            cursor: [m[0] / 2, m[1] / 2, 0],
            kind,
            magnitude_index: 0,
            positive,
            angle_deg: ANTENNA_ANGLES[0],
            message: None,
            revision: 1,
        }
    }

    pub fn subdivision(&self) -> u32 {
        self.level.grid.subdivision
    }

    /// Value of a new element of the selected kind.
    pub fn selected_value(&self) -> f64 {
        let m = magnitudes(&self.level, self.kind)
            .get(self.magnitude_index)
            .copied()
            .unwrap_or(1.0);
        if self.positive { m } else { -m }
    }

    pub fn left(&self, kind: ElementKind) -> usize {
        let max = max_of(&self.level, kind) as usize;
        let used = self.placement.iter().filter(|e| e.kind == kind).count();
        max.saturating_sub(used)
    }

    /// Selects the kind of new elements (if the level allows it).
    pub fn set_kind(&mut self, kind: ElementKind) {
        if max_of(&self.level, kind) > 0 && kind != self.kind {
            self.kind = kind;
            self.magnitude_index = 0;
            self.positive = default_positive(&self.level, kind);
        }
    }

    /// Switches to the next kind the level allows.
    pub fn toggle_kind(&mut self) {
        let i = KINDS.iter().position(|&k| k == self.kind).unwrap_or(0);
        for step in 1..KINDS.len() {
            let k = KINDS[(i + step) % KINDS.len()];
            if max_of(&self.level, k) > 0 {
                self.set_kind(k);
                return;
            }
        }
    }

    /// Rotates the antenna under the cursor, or the orientation of new antennas, by one
    /// step of `ANTENNA_ANGLES` (45°).
    pub fn rotate(&mut self, step: isize) {
        let next = |a: f64| {
            let n = ANTENNA_ANGLES.len() as isize;
            let i = ANTENNA_ANGLES
                .iter()
                .position(|x| x.to_bits() == a.to_bits())
                .unwrap_or(0) as isize;
            ANTENNA_ANGLES[(i + step).rem_euclid(n) as usize]
        };
        if let Some(i) = self.player_index_at(self.cursor) {
            if self.placement[i].kind == ElementKind::Antenna {
                let mut trial = self.placement.clone();
                trial[i].angle_deg = next(trial[i].angle_deg);
                let _ = self.try_placement(trial);
            }
        } else {
            self.angle_deg = next(self.angle_deg);
        }
    }

    fn changed(&mut self) {
        self.revision += 1;
        self.message = None;
    }

    fn player_index_at(&self, node: Node) -> Option<usize> {
        self.placement.iter().position(|c| c.node == node)
    }

    pub fn move_cursor(&mut self, dx: i64, dy: i64) {
        let m = self.level.grid.max_node();
        self.cursor[0] = (self.cursor[0] + dx).clamp(0, m[0]);
        self.cursor[1] = (self.cursor[1] + dy).clamp(0, m[1]);
    }

    pub fn set_cursor(&mut self, node: Node) {
        if self.level.grid.contains(node) {
            self.cursor = node;
        }
    }

    /// Tries a new placement; keeps it if the level allows it.
    fn try_placement(&mut self, trial: Vec<Element>) -> Result<(), PlacementError> {
        match self.level.check_placement(&trial) {
            Ok(()) => {
                self.placement = trial;
                self.changed();
                Ok(())
            }
            Err(e) => {
                self.message = Some(describe(&e));
                Err(e)
            }
        }
    }

    /// Places the selected element at the cursor (replacing a player element there).
    pub fn place(&mut self) -> Result<(), PlacementError> {
        let mut trial = self.placement.clone();
        let e = Element {
            node: self.cursor,
            kind: self.kind,
            value: self.selected_value(),
            angle_deg: if self.kind == ElementKind::Antenna {
                self.angle_deg
            } else {
                0.0
            },
        };
        match self.player_index_at(self.cursor) {
            Some(i) => trial[i] = e,
            None => trial.push(e),
        }
        self.try_placement(trial)
    }

    pub fn remove(&mut self) {
        if let Some(i) = self.player_index_at(self.cursor) {
            self.placement.remove(i);
            self.changed();
        }
    }

    /// Flips the sign of the element under the cursor, or of the selection if none.
    pub fn flip_sign(&mut self) {
        if let Some(i) = self.player_index_at(self.cursor) {
            let mut trial = self.placement.clone();
            trial[i].value = -trial[i].value;
            let _ = self.try_placement(trial);
        } else {
            if both_signs(&self.level, self.kind) {
                self.positive = !self.positive;
            }
        }
    }

    /// Cycles the magnitude of the element under the cursor, or of the selection if none.
    pub fn cycle_magnitude(&mut self, step: isize) {
        let kind = self
            .player_index_at(self.cursor)
            .map_or(self.kind, |i| self.placement[i].kind);
        let list = magnitudes(&self.level, kind).to_vec();
        let n = list.len();
        if n == 0 {
            return;
        }
        let next = |i: usize| (i as isize + step).rem_euclid(n as isize) as usize;
        if let Some(i) = self.player_index_at(self.cursor) {
            let v = self.placement[i].value;
            // Magnitudes are exact values from the level's list.
            let current = list
                .iter()
                .position(|&m| m.to_bits() == v.abs().to_bits())
                .unwrap_or(0);
            let mut trial = self.placement.clone();
            trial[i].value = v.signum() * list[next(current)];
            let _ = self.try_placement(trial);
        } else {
            self.magnitude_index = next(self.magnitude_index);
        }
    }

    pub fn clear(&mut self) {
        if !self.placement.is_empty() {
            self.placement.clear();
            self.changed();
        }
    }

    /// Sets the grid refinement factor relative to the level's recommended grid. Elements
    /// keep their positions; returning to a coarser grid is only possible while every
    /// player element still lies on a node of it.
    pub fn set_refinement(&mut self, factor: u32) {
        let current = self.subdivision() / self.base.grid.subdivision;
        if factor == current || factor == 0 {
            return;
        }
        let from = i64::from(current);
        let to = i64::from(factor);
        let convert = |n: Node| -> Option<Node> {
            let mut out = [0; 3];
            for i in 0..3 {
                let v = n[i] * to;
                if v % from != 0 {
                    return None;
                }
                out[i] = v / from;
            }
            Some(out)
        };
        let Some(placement) = self
            .placement
            .iter()
            .map(|c| convert(c.node).map(|node| Element { node, ..*c }))
            .collect::<Option<Vec<_>>>()
        else {
            self.message = Some("An element is not on a node of the coarser grid.".into());
            return;
        };
        let mut level = self.base.clone();
        level.refine(factor, &mut []);
        self.cursor = convert(self.cursor).unwrap_or([
            self.cursor[0] * to / from,
            self.cursor[1] * to / from,
            0,
        ]);
        self.level = level;
        self.placement = placement;
        self.changed();
    }

    /// The level at its own (recommended) grid.
    pub fn base(&self) -> &Level {
        &self.base
    }

    /// Edits the level itself (sandbox). Works on the level's own grid: the refinement is
    /// reset to 1 first (player elements that do not fit it are removed).
    pub fn edit_level(&mut self, f: impl FnOnce(&mut Level)) {
        if self.refinement() != 1 {
            self.set_refinement(1);
            if self.refinement() != 1 {
                self.placement.clear();
                self.level = self.base.clone();
            }
        }
        f(&mut self.base);
        self.level = self.base.clone();
        let m = self.level.grid.max_node();
        self.cursor = [
            self.cursor[0].clamp(0, m[0]),
            self.cursor[1].clamp(0, m[1]),
            0,
        ];
        let n = magnitudes(&self.level, self.kind).len();
        self.magnitude_index = self.magnitude_index.min(n.saturating_sub(1));
        self.changed();
    }

    pub fn refinement(&self) -> u32 {
        self.subdivision() / self.base.grid.subdivision
    }
}

pub fn describe(e: &PlacementError) -> String {
    match e {
        PlacementError::TooManyCharges => "No charges left for this level.".into(),
        PlacementError::TooManyMagnets => "No magnets left for this level.".into(),
        PlacementError::TooManyAntennas => "No antennas left for this level.".into(),
        PlacementError::OutsideGrid(_) => "Outside the grid.".into(),
        PlacementError::OutsideRegion(_) => "Elements can only go in the marked region.".into(),
        PlacementError::NotInPlane(_) => "Must be in the plane.".into(),
        PlacementError::Occupied(_) => "That node is occupied.".into(),
        PlacementError::SignNotAllowed(_) => "That sign is not allowed here.".into(),
        PlacementError::MagnitudeNotAllowed(_) => "That value is not allowed here.".into(),
        PlacementError::AngleNotAllowed(_) => {
            "Antennas can point along 0°, 45°, 90° or 135°.".into()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn level() -> Level {
        Level::from_json(include_str!("../../../levels/01_first_bend.json")).unwrap()
    }

    #[test]
    fn place_edit_remove() {
        let mut e = Editor::new(level());
        e.set_cursor([5, 5, 0]);
        e.place().unwrap();
        assert_eq!(e.placement.len(), 1);
        let q0 = e.placement[0].value;
        e.flip_sign();
        assert_eq!(e.placement[0].value.to_bits(), (-q0).to_bits());
        e.cycle_magnitude(1);
        assert_ne!(e.placement[0].value.abs().to_bits(), q0.abs().to_bits());
        // The level allows one charge: a second one elsewhere is rejected.
        e.set_cursor([6, 5, 0]);
        assert!(e.place().is_err());
        assert!(e.message.is_some());
        e.set_cursor([5, 5, 0]);
        e.remove();
        assert!(e.placement.is_empty());
        // No magnets in this level: the kind stays Charge.
        e.toggle_kind();
        assert_eq!(e.kind, ElementKind::Charge);
    }

    #[test]
    fn refinement_round_trip_keeps_positions() {
        let mut e = Editor::new(level());
        e.set_cursor([4, 6, 0]);
        e.place().unwrap();
        let before = e.level.grid.position(e.placement[0].node);
        e.set_refinement(3);
        assert_eq!(e.level.grid.position(e.placement[0].node), before);
        e.set_refinement(1);
        assert_eq!(e.level.grid.position(e.placement[0].node), before);
        // An element on a fine-only node blocks coarsening.
        e.set_refinement(2);
        e.set_cursor([4 * 2, 6 * 2, 0]);
        e.remove();
        e.set_cursor([9, 13, 0]);
        e.place().unwrap();
        e.set_refinement(1);
        assert_eq!(e.refinement(), 2);
    }
}
