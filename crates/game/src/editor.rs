//! Editor state: player charges on the grid, cursor, selected charge value. Independent of
//! the rendering engine so it can be unit tested.

use level::{Charge, Level, Node, PlacementError};

pub struct Editor {
    /// The level as loaded (recommended grid).
    base: Level,
    /// The level at the current refinement.
    pub level: Level,
    pub placement: Vec<Charge>,
    pub cursor: Node,
    /// Index into `level.limits.magnitudes` for new charges.
    pub magnitude_index: usize,
    pub positive: bool,
    /// Last rejected action, for the status line.
    pub message: Option<String>,
    /// Incremented on every change of the physical setup.
    pub revision: u64,
}

impl Editor {
    pub fn new(level: Level) -> Self {
        let m = level.grid.max_node();
        let positive = level.limits.allow_positive;
        Self {
            base: level.clone(),
            level,
            placement: Vec::new(),
            cursor: [m[0] / 2, m[1] / 2, 0],
            magnitude_index: 0,
            positive,
            message: None,
            revision: 1,
        }
    }

    pub fn subdivision(&self) -> u32 {
        self.level.grid.subdivision
    }

    /// Value of a new charge.
    pub fn selected_charge(&self) -> f64 {
        let q = self
            .level
            .limits
            .magnitudes
            .get(self.magnitude_index)
            .copied()
            .unwrap_or(1.0);
        if self.positive { q } else { -q }
    }

    pub fn charges_left(&self) -> usize {
        (self.level.limits.max_charges as usize).saturating_sub(self.placement.len())
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
    fn try_placement(&mut self, trial: Vec<Charge>) -> Result<(), PlacementError> {
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

    /// Places the selected charge at the cursor (replacing a player charge there).
    pub fn place(&mut self) -> Result<(), PlacementError> {
        let mut trial = self.placement.clone();
        let c = Charge {
            node: self.cursor,
            charge: self.selected_charge(),
        };
        match self.player_index_at(self.cursor) {
            Some(i) => trial[i] = c,
            None => trial.push(c),
        }
        self.try_placement(trial)
    }

    pub fn remove(&mut self) {
        if let Some(i) = self.player_index_at(self.cursor) {
            self.placement.remove(i);
            self.changed();
        }
    }

    /// Flips the sign of the charge under the cursor, or of the selection if none.
    pub fn flip_sign(&mut self) {
        if let Some(i) = self.player_index_at(self.cursor) {
            let mut trial = self.placement.clone();
            trial[i].charge = -trial[i].charge;
            let _ = self.try_placement(trial);
        } else {
            let limits = &self.level.limits;
            if limits.allow_positive && limits.allow_negative {
                self.positive = !self.positive;
            }
        }
    }

    /// Cycles the magnitude of the charge under the cursor, or of the selection if none.
    pub fn cycle_magnitude(&mut self, step: isize) {
        let n = self.level.limits.magnitudes.len();
        let next = |i: usize| (i as isize + step).rem_euclid(n as isize) as usize;
        if let Some(i) = self.player_index_at(self.cursor) {
            let q = self.placement[i].charge;
            let current = self
                .level
                .limits
                .magnitudes
                .iter()
                // Magnitudes are exact values from the level's list.
                .position(|&m| m.to_bits() == q.abs().to_bits())
                .unwrap_or(0);
            let mut trial = self.placement.clone();
            trial[i].charge = q.signum() * self.level.limits.magnitudes[next(current)];
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

    /// Sets the grid refinement factor relative to the level's recommended grid. Charges
    /// keep their positions; returning to a coarser grid is only possible while every
    /// player charge still lies on a node of it.
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
            .map(|c| {
                convert(c.node).map(|node| Charge {
                    node,
                    charge: c.charge,
                })
            })
            .collect::<Option<Vec<_>>>()
        else {
            self.message = Some("A charge is not on a node of the coarser grid.".into());
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
    /// reset to 1 first (player charges that do not fit it are removed).
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
        let n = self.level.limits.magnitudes.len();
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
        PlacementError::OutsideGrid(_) => "Outside the grid.".into(),
        PlacementError::OutsideRegion(_) => "Charges can only go in the marked region.".into(),
        PlacementError::NotInPlane(_) => "Must be in the plane.".into(),
        PlacementError::Occupied(_) => "That node is occupied.".into(),
        PlacementError::SignNotAllowed(_) => "That sign is not allowed here.".into(),
        PlacementError::MagnitudeNotAllowed(_) => "That magnitude is not allowed here.".into(),
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
        let q0 = e.placement[0].charge;
        e.flip_sign();
        assert_eq!(e.placement[0].charge.to_bits(), (-q0).to_bits());
        e.cycle_magnitude(1);
        assert_ne!(e.placement[0].charge.abs().to_bits(), q0.abs().to_bits());
        // The level allows one charge: a second one elsewhere is rejected.
        e.set_cursor([6, 5, 0]);
        assert!(e.place().is_err());
        assert!(e.message.is_some());
        e.set_cursor([5, 5, 0]);
        e.remove();
        assert!(e.placement.is_empty());
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
        // A charge on a fine-only node blocks coarsening.
        e.set_refinement(2);
        e.set_cursor([9, 13, 0]);
        e.remove();
        e.set_cursor([4 * 2, 6 * 2, 0]);
        e.remove();
        e.set_cursor([9, 13, 0]);
        e.place().unwrap();
        e.set_refinement(1);
        assert_eq!(e.refinement(), 2);
    }
}
