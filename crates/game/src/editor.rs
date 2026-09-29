//! Editor state: player elements (charges, magnets, antennas, plates, power supplies)
//! on the grid, cursor, selected element.
//! Independent of the rendering engine so it can be unit tested.

use level::{
    ANTENNA_ANGLES, ConductorBias, Element, ElementKind, FREE_ANGLES, Level, Node, PLATE_ANGLES,
    PlacementError, value_range,
};

/// Kinds in palette order. Power supplies are not placed from the palette: they are
/// operated on their electrodes.
pub const KINDS: [ElementKind; 5] = [
    ElementKind::Charge,
    ElementKind::Magnet,
    ElementKind::Antenna,
    ElementKind::Plate,
    ElementKind::Free,
];

/// Kinds whose values are signed values from a list (plates, power supplies: potentials;
/// free charges: charges), not a magnitude with a separate sign.
pub fn is_signed(kind: ElementKind) -> bool {
    matches!(
        kind,
        ElementKind::Plate | ElementKind::Supply | ElementKind::Free
    )
}

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
    /// Orientation of new plates, degrees (one of `PLATE_ANGLES`).
    pub plate_angle_deg: f64,
    /// Index into the level's allowed antenna frequencies for new antennas.
    pub omega_index: usize,
    /// Launch direction of new free charges, degrees (one of `FREE_ANGLES`).
    pub free_angle_deg: f64,
    /// Index into the level's launch speeds for new free charges.
    pub speed_index: usize,
    /// Hardcore mode: magnitudes of new elements per kind (`KINDS` order), the
    /// frequency of new antennas and the speed of new free charges, anywhere in the
    /// level's ranges.
    pub continuous_magnitude: [f64; 5],
    pub continuous_omega: f64,
    pub continuous_speed: f64,
    /// Player element being moved (index) and its node when it was grabbed.
    pub grabbed: Option<(usize, Node)>,
    /// Element node minus cursor node when grabbed (a plate can be picked up anywhere).
    grab_offset: Node,
    /// Last rejected action, for the status line.
    pub message: Option<String>,
    /// Incremented on every change of the physical setup.
    pub revision: u64,
}

/// Allowed magnitudes of a kind in a level (signed potentials for plates and power
/// supplies, see `is_signed`).
pub fn magnitudes(level: &Level, kind: ElementKind) -> &[f64] {
    match kind {
        ElementKind::Charge => &level.limits.magnitudes,
        ElementKind::Magnet => &level.limits.magnet_strengths,
        ElementKind::Antenna => &level.limits.antenna_amplitudes,
        ElementKind::Plate => &level.limits.plate_voltages,
        ElementKind::Supply => &level.limits.supply_voltages,
        ElementKind::Free => &level.limits.free_charges,
    }
}

/// Maximum number of player elements of a kind.
pub fn max_of(level: &Level, kind: ElementKind) -> u32 {
    match kind {
        ElementKind::Charge => level.limits.max_charges,
        ElementKind::Magnet => level.limits.max_magnets,
        ElementKind::Antenna => level.limits.max_antennas,
        ElementKind::Plate => level.limits.max_plates,
        ElementKind::Supply => {
            u32::try_from(level.electrodes.iter().filter(|e| e.tunable).count()).unwrap_or(0)
        }
        ElementKind::Free => level.limits.max_free,
    }
}

/// Step of a signed potential per key press in hardcore mode: 1/20 of its range.
fn signed_step(list: &[f64]) -> f64 {
    value_range(list).map_or(0.0, |(lo, hi)| (hi - lo) / 20.0)
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
        let continuous_magnitude =
            KINDS.map(|k| magnitudes(&level, k).first().copied().unwrap_or(1.0));
        let continuous_omega = level.limits.antenna_omegas.first().copied().unwrap_or(1.0);
        let continuous_speed = level.limits.free_speeds.first().copied().unwrap_or(0.0);
        Self {
            base: level.clone(),
            level,
            placement: Vec::new(),
            cursor: [m[0] / 2, m[1] / 2, 0],
            kind,
            magnitude_index: 0,
            positive,
            angle_deg: ANTENNA_ANGLES[0],
            plate_angle_deg: PLATE_ANGLES[0],
            omega_index: 0,
            free_angle_deg: 0.0,
            speed_index: 0,
            continuous_magnitude,
            continuous_omega,
            continuous_speed,
            grabbed: None,
            grab_offset: [0; 3],
            message: None,
            revision: 1,
        }
    }

    pub fn subdivision(&self) -> u32 {
        self.level.grid.subdivision
    }

    /// Value of a new element of the selected kind.
    pub fn selected_value(&self) -> f64 {
        let m = if self.continuous() {
            self.continuous_magnitude[kind_index(self.kind)]
        } else {
            magnitudes(&self.level, self.kind)
                .get(self.magnitude_index)
                .copied()
                .unwrap_or(1.0)
        };
        if is_signed(self.kind) || self.positive {
            m
        } else {
            -m
        }
    }

    /// Hardcore mode: continuous values instead of the level's lists.
    pub fn continuous(&self) -> bool {
        self.level.limits.continuous
    }

    /// Switches hardcore mode. Leaving it snaps every player element to the nearest
    /// listed value, orientation and frequency.
    pub fn set_continuous(&mut self, on: bool) {
        if on == self.continuous() {
            return;
        }
        self.base.limits.continuous = on;
        self.level.limits.continuous = on;
        if on {
            // Start the sliders at the currently selected discrete values.
            for (i, k) in KINDS.into_iter().enumerate() {
                let list = magnitudes(&self.level, k);
                let idx = if k == self.kind {
                    self.magnitude_index
                } else {
                    0
                };
                if let Some(&m) = list.get(idx.min(list.len().saturating_sub(1))) {
                    self.continuous_magnitude[i] = m;
                }
            }
            if let Some(w) = self.selected_omega_discrete() {
                self.continuous_omega = w;
            }
            self.continuous_speed = self.selected_speed_discrete();
        } else {
            let level = self.level.clone();
            for e in &mut self.placement {
                snap(&level, e);
            }
            self.angle_deg = nearest_angle(self.angle_deg).0;
            self.plate_angle_deg = nearest_plate_angle(self.plate_angle_deg);
            self.free_angle_deg = nearest_free_angle(self.free_angle_deg);
        }
        self.changed();
    }

    /// Launch speed of new free charges: one of the level's speeds (hardcore: anything in
    /// their range).
    pub fn selected_speed(&self) -> f64 {
        if self.continuous() && !self.level.limits.free_speeds.is_empty() {
            return self.continuous_speed;
        }
        self.selected_speed_discrete()
    }

    fn selected_speed_discrete(&self) -> f64 {
        let list = &self.level.limits.free_speeds;
        list.get(self.speed_index.min(list.len().saturating_sub(1)))
            .copied()
            .unwrap_or(0.0)
    }

    /// Selects the speed of new free charges (index into the level's list).
    pub fn set_speed_index(&mut self, i: usize) {
        self.speed_index = i.min(self.level.limits.free_speeds.len().saturating_sub(1));
    }

    /// Sets the launch velocity of free charge `i` from a drag of its arrow: direction
    /// `angle_deg` and speed `speed`, continuous (the speed clamped to the range of the
    /// level's speeds, any direction). Kept if allowed.
    pub fn set_velocity(&mut self, i: usize, angle_deg: f64, speed: f64) {
        let Some(&e) = self.placement.get(i) else {
            return;
        };
        if e.kind != ElementKind::Free {
            return;
        }
        let speeds = &self.level.limits.free_speeds;
        let (lo, hi) = value_range(speeds).unwrap_or((0.0, 0.0));
        let (angle, v) = (angle_deg.rem_euclid(360.0), speed.clamp(lo, hi));
        let mut new = e;
        new.angle_deg = if v == 0.0 { 0.0 } else { angle };
        new.speed = Some(v);
        self.set_element(i, new);
    }

    /// Changes the player element at index `i` (hardcore sliders); kept if allowed.
    pub fn set_element(&mut self, i: usize, e: Element) {
        if self.placement.get(i) != Some(&e) {
            let mut trial = self.placement.clone();
            trial[i] = e;
            let _ = self.try_placement(trial);
        }
    }

    /// Index of the player element under the cursor.
    pub fn element_at_cursor(&self) -> Option<usize> {
        self.player_index_at(self.cursor)
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

    /// Rotates the antenna or plate under the cursor, or the orientation of new ones, by
    /// one step of `ANTENNA_ANGLES` (45°) or `PLATE_ANGLES` (90°).
    pub fn rotate(&mut self, step: isize) {
        let on = self.player_index_at(self.cursor);
        if on.map_or(self.kind, |i| self.placement[i].kind) == ElementKind::Free {
            let turn = |a: f64| {
                let n = FREE_ANGLES.len() as isize;
                let i = FREE_ANGLES
                    .iter()
                    .position(|x| x.to_bits() == nearest_free_angle(a).to_bits())
                    .unwrap_or(0) as isize;
                FREE_ANGLES[(i + step).rem_euclid(n) as usize]
            };
            match on {
                Some(i) => {
                    let mut e = self.placement[i];
                    e.angle_deg = turn(e.angle_deg);
                    self.set_element(i, e);
                }
                None => self.free_angle_deg = turn(self.free_angle_deg),
            }
            return;
        }
        let plate = on.map_or(self.kind, |i| self.placement[i].kind) == ElementKind::Plate;
        if plate {
            let turn = |a: f64| {
                if self.continuous() {
                    (a + 45.0 * f64::from(step_i32(step))).rem_euclid(180.0)
                } else {
                    let n = PLATE_ANGLES.len() as isize;
                    let i = PLATE_ANGLES
                        .iter()
                        .position(|x| x.to_bits() == a.to_bits())
                        .unwrap_or(0) as isize;
                    PLATE_ANGLES[(i + step).rem_euclid(n) as usize]
                }
            };
            match on {
                Some(i) => {
                    let mut e = self.placement[i];
                    e.angle_deg = turn(e.angle_deg);
                    self.set_element(i, e);
                }
                None => self.plate_angle_deg = turn(self.plate_angle_deg),
            }
            return;
        }
        if self.continuous() {
            let turn = |a: f64| (a + 45.0 * f64::from(step_i32(step))).rem_euclid(360.0);
            if let Some(i) = self.player_index_at(self.cursor) {
                if self.placement[i].kind == ElementKind::Antenna {
                    let mut e = self.placement[i];
                    e.angle_deg = turn(e.angle_deg);
                    self.set_element(i, e);
                }
            } else {
                self.angle_deg = turn(self.angle_deg);
            }
            return;
        }
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

    /// Frequency of new antennas: one of the level's allowed values (hardcore: anything
    /// in their range), or `None` (the level's RF generator) if it lists none.
    pub fn selected_omega(&self) -> Option<f64> {
        if self.continuous() && !self.level.limits.antenna_omegas.is_empty() {
            return Some(self.continuous_omega);
        }
        self.selected_omega_discrete()
    }

    fn selected_omega_discrete(&self) -> Option<f64> {
        let list = &self.level.limits.antenna_omegas;
        list.get(self.omega_index.min(list.len().saturating_sub(1)))
            .copied()
    }

    /// Changes the frequency of the antenna under the cursor, or of new antennas, to the
    /// next allowed value.
    pub fn cycle_omega(&mut self, step: isize) {
        let list = self.level.limits.antenna_omegas.clone();
        let n = list.len();
        if n == 0 {
            return;
        }
        if self.continuous() {
            let (lo, hi) = value_range(&list).expect("non-empty");
            let scale = |w: f64| (w * STEP.powi(step_i32(step))).clamp(lo, hi);
            if let Some(i) = self.player_index_at(self.cursor) {
                if self.placement[i].kind == ElementKind::Antenna {
                    let mut e = self.placement[i];
                    e.omega = e.omega.map(scale);
                    self.set_element(i, e);
                }
            } else {
                self.continuous_omega = scale(self.continuous_omega);
            }
            return;
        }
        let next = |i: usize| (i as isize + step).clamp(0, n as isize - 1) as usize;
        if let Some(i) = self.player_index_at(self.cursor) {
            if self.placement[i].kind == ElementKind::Antenna {
                let current = list
                    .iter()
                    .position(|w| Some(w.to_bits()) == self.placement[i].omega.map(f64::to_bits))
                    .unwrap_or(0);
                let mut trial = self.placement.clone();
                trial[i].omega = Some(list[next(current)]);
                let _ = self.try_placement(trial);
            }
        } else {
            self.omega_index = next(self.omega_index);
        }
    }

    fn changed(&mut self) {
        self.revision += 1;
        self.message = None;
    }

    /// Replaces the player's elements (nodes on the current, possibly refined, grid).
    pub fn set_placement(&mut self, placement: Vec<Element>) {
        self.placement = placement;
        self.changed();
    }

    /// The player element at a node: one on it, else a plate covering it, else the power
    /// supply of the tunable electrode under it.
    fn player_index_at(&self, node: Node) -> Option<usize> {
        if let Some(i) = self.placement.iter().position(|c| c.node == node) {
            return Some(i);
        }
        let p = self.level.grid.position(node);
        if let Some(i) = self.placement.iter().position(|e| {
            e.kind == ElementKind::Plate
                && physics::bem::Electrodes::shapes_only(vec![self.level.plate_box(e)])
                    .contains(p, PICK_MARGIN)
        }) {
            return Some(i);
        }
        let centre = self.tunable_at(node)?;
        self.placement
            .iter()
            .position(|e| e.kind == ElementKind::Supply && e.node == centre)
    }

    /// Centre node of the tunable level electrode under a node.
    pub fn tunable_at(&self, node: Node) -> Option<Node> {
        let p = self.level.grid.position(node);
        let boxes = self.level.box_electrodes();
        self.level
            .electrodes
            .iter()
            .zip(boxes)
            .find(|(e, b)| {
                e.tunable
                    && physics::bem::Electrodes::shapes_only(vec![*b]).contains(p, PICK_MARGIN)
            })
            .map(|(e, _)| e.center)
    }

    /// Switches on the power supply of the tunable electrode centred at `centre`, at the
    /// listed potential nearest to the electrode's own bias.
    fn add_supply(&mut self, centre: Node) -> Result<(), PlacementError> {
        let Some(el) = self.level.electrodes.iter().find(|e| e.center == centre) else {
            return Ok(());
        };
        let own = match el.bias {
            ConductorBias::Potential(v) => v,
            ConductorBias::Grounded | ConductorBias::Charge(_) => 0.0,
        };
        let Some(v) = nearest_linear(&self.level.limits.supply_voltages, own) else {
            return Ok(());
        };
        let mut trial = self.placement.clone();
        trial.push(Element::supply(centre, v));
        self.try_placement(trial)
    }

    pub fn move_cursor(&mut self, dx: i64, dy: i64) {
        let m = self.level.grid.max_node();
        self.cursor[0] = (self.cursor[0] + dx).clamp(0, m[0]);
        self.cursor[1] = (self.cursor[1] + dy).clamp(0, m[1]);
        self.follow_grabbed();
    }

    pub fn set_cursor(&mut self, node: Node) {
        if self.level.grid.contains(node) {
            self.cursor = node;
            self.follow_grabbed();
        }
    }

    /// Picks up the player element under the cursor to move it; returns whether one was
    /// there.
    pub fn grab(&mut self) -> bool {
        match self.player_index_at(self.cursor) {
            // Power supplies stay on their electrodes.
            Some(i) if self.placement[i].kind != ElementKind::Supply => {
                let node = self.placement[i].node;
                self.grabbed = Some((i, node));
                self.grab_offset = [0, 1, 2].map(|k| node[k] - self.cursor[k]);
                true
            }
            _ => false,
        }
    }

    /// Puts the grabbed element down where it is.
    pub fn drop_grabbed(&mut self) {
        self.grabbed = None;
    }

    /// Returns the grabbed element to where it was picked up.
    pub fn cancel_grab(&mut self) {
        if let Some((i, home)) = self.grabbed.take() {
            self.move_element(i, home);
        }
    }

    /// The grabbed element follows the cursor, to every node where it is allowed
    /// (over a blocked node it waits and jumps on when the cursor is past it).
    fn follow_grabbed(&mut self) {
        if let Some((i, _)) = self.grabbed {
            let target = [0, 1, 2].map(|k| self.cursor[k] + self.grab_offset[k]);
            self.move_element(i, target);
        }
    }

    fn move_element(&mut self, i: usize, node: Node) {
        if self.placement.get(i).is_none_or(|e| e.node == node) {
            return;
        }
        let mut trial = self.placement.clone();
        trial[i].node = node;
        if self.level.check_placement(&trial).is_ok() {
            self.placement = trial;
            self.changed();
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
    /// On a tunable electrode it operates its power supply instead: switches it on, or
    /// steps it to the next potential.
    pub fn place(&mut self) -> Result<(), PlacementError> {
        if let Some(centre) = self.tunable_at(self.cursor) {
            if self
                .placement
                .iter()
                .any(|e| e.kind == ElementKind::Supply && e.node == centre)
            {
                self.cycle_magnitude(1);
                return Ok(());
            }
            return self.add_supply(centre);
        }
        let mut trial = self.placement.clone();
        let e = Element {
            node: self.cursor,
            kind: self.kind,
            value: self.selected_value(),
            angle_deg: match self.kind {
                ElementKind::Antenna => self.angle_deg,
                ElementKind::Plate => self.plate_angle_deg,
                ElementKind::Free if self.selected_speed() != 0.0 => self.free_angle_deg,
                _ => 0.0,
            },
            omega: if self.kind == ElementKind::Antenna {
                self.selected_omega()
            } else {
                None
            },
            speed: (self.kind == ElementKind::Free).then(|| self.selected_speed()),
        };
        match self.player_index_at(self.cursor) {
            Some(i) => trial[i] = e,
            None => trial.push(e),
        }
        self.try_placement(trial)
    }

    pub fn remove(&mut self) {
        self.grabbed = None;
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
        } else if is_signed(self.kind) {
            // The opposite potential, where the level lists it.
            if self.continuous() {
                let k = kind_index(self.kind);
                let (lo, hi) =
                    value_range(magnitudes(&self.level, self.kind)).unwrap_or((0.0, 0.0));
                self.continuous_magnitude[k] = (-self.continuous_magnitude[k]).clamp(lo, hi);
            } else {
                let v = self.selected_value();
                if let Some(j) = magnitudes(&self.level, self.kind)
                    .iter()
                    .position(|m| m.to_bits() == (-v).to_bits())
                {
                    self.magnitude_index = j;
                }
            }
        } else if both_signs(&self.level, self.kind) {
            self.positive = !self.positive;
        }
    }

    /// Cycles the magnitude of the element under the cursor, or of the selection if none.
    pub fn cycle_magnitude(&mut self, step: isize) {
        let on = self.player_index_at(self.cursor);
        if on.is_none()
            && let Some(centre) = self.tunable_at(self.cursor)
        {
            let _ = self.add_supply(centre);
            return;
        }
        let kind = on.map_or(self.kind, |i| self.placement[i].kind);
        let list = magnitudes(&self.level, kind).to_vec();
        let n = list.len();
        if n == 0 {
            return;
        }
        if is_signed(kind) {
            self.cycle_signed(on, kind, &list, step);
            return;
        }
        if self.continuous() {
            let (lo, hi) = value_range(&list).expect("non-empty");
            let scale = |m: f64| (m * STEP.powi(step_i32(step))).clamp(lo, hi);
            if let Some(i) = self.player_index_at(self.cursor) {
                let mut e = self.placement[i];
                e.value = e.value.signum() * scale(e.value.abs());
                self.set_element(i, e);
            } else {
                let k = kind_index(kind);
                self.continuous_magnitude[k] = scale(self.continuous_magnitude[k]);
            }
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

    /// `cycle_magnitude` for signed potentials: the next listed value (hardcore: a step of
    /// 1/20 of the range).
    fn cycle_signed(&mut self, on: Option<usize>, kind: ElementKind, list: &[f64], step: isize) {
        let n = list.len();
        if self.continuous() {
            let (lo, hi) = value_range(list).expect("non-empty");
            let d = signed_step(list) * f64::from(step_i32(step));
            match on {
                Some(i) => {
                    let mut e = self.placement[i];
                    e.value = (e.value + d).clamp(lo, hi);
                    self.set_element(i, e);
                }
                None => {
                    let k = kind_index(kind);
                    self.continuous_magnitude[k] = (self.continuous_magnitude[k] + d).clamp(lo, hi);
                }
            }
            return;
        }
        let next = |i: usize| (i as isize + step).rem_euclid(n as isize) as usize;
        match on {
            Some(i) => {
                let v = self.placement[i].value;
                let current = list
                    .iter()
                    .position(|m| m.to_bits() == v.to_bits())
                    .unwrap_or(0);
                let mut trial = self.placement.clone();
                trial[i].value = list[next(current)];
                let _ = self.try_placement(trial);
            }
            None => self.magnitude_index = next(self.magnitude_index),
        }
    }

    /// Potential of the power supply of the tunable electrode centred at `centre`
    /// (`None`: switched off, the electrode keeps its own bias).
    pub fn supply(&self, centre: Node) -> Option<f64> {
        self.placement
            .iter()
            .find(|e| e.kind == ElementKind::Supply && e.node == centre)
            .map(|e| e.value)
    }

    /// Sets or switches off the power supply of the tunable electrode centred at
    /// `centre`; kept if allowed.
    pub fn set_supply(&mut self, centre: Node, value: Option<f64>) {
        if self.supply(centre).map(f64::to_bits) == value.map(f64::to_bits) {
            return;
        }
        let mut trial = self.placement.clone();
        let at = trial
            .iter()
            .position(|e| e.kind == ElementKind::Supply && e.node == centre);
        match (at, value) {
            (Some(i), Some(v)) => trial[i].value = v,
            (Some(i), None) => {
                trial.remove(i);
            }
            (None, Some(v)) => trial.push(Element::supply(centre, v)),
            (None, None) => {}
        }
        let _ = self.try_placement(trial);
    }

    pub fn clear(&mut self) {
        self.grabbed = None;
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

/// Factor per key press for continuous values (hardcore).
const STEP: f64 = 1.1;

/// How far outside a plate or electrode the cursor still picks it, in cells.
const PICK_MARGIN: f64 = 0.3;

/// Nearest listed value to `v` on a linear scale (signed potentials).
fn nearest_linear(list: &[f64], v: f64) -> Option<f64> {
    list.iter()
        .copied()
        .min_by(|a, b| (a - v).abs().total_cmp(&(b - v).abs()))
}

/// How a velocity arrow's length measures speed: in the speed itself, or in rapidity
/// `c·artanh(v/c)` (additive under boosts along a line; equal to v for slow particles, and
/// unbounded as v approaches c).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArrowMeasure {
    Speed,
    Rapidity,
}

/// The default measure: rapidity in relativistic levels (finite c and a listed speed
/// above 0.3 c), speed otherwise.
pub fn default_measure(level: &Level) -> ArrowMeasure {
    let c = level.c();
    let fast = level.limits.free_speeds.iter().any(|&v| v > 0.3 * c)
        || level
            .free_particles
            .iter()
            .any(|f| f.velocity[0].hypot(f.velocity[1]) > 0.3 * c);
    if c.is_finite() && fast {
        ArrowMeasure::Rapidity
    } else {
        ArrowMeasure::Speed
    }
}

/// The arrow length (in speed units, before `arrow_scale`) of a speed `v`.
pub fn arrow_length(level: &Level, v: f64, m: ArrowMeasure) -> f64 {
    let c = level.c();
    match m {
        ArrowMeasure::Rapidity if c.is_finite() => c * (v / c).min(1.0 - 1e-12).atanh(),
        _ => v,
    }
}

/// The speed of an arrow length (inverse of `arrow_length`).
pub fn arrow_speed(level: &Level, len: f64, m: ArrowMeasure) -> f64 {
    let c = level.c();
    match m {
        ArrowMeasure::Rapidity if c.is_finite() => c * (len / c).tanh(),
        _ => len,
    }
}

/// Cells of arrow per unit of arrow length: the longest listed launch is drawn 3 cells
/// long.
pub fn arrow_scale(level: &Level, m: ArrowMeasure) -> f64 {
    let longest = level
        .limits
        .free_speeds
        .iter()
        .copied()
        .chain(
            level
                .free_particles
                .iter()
                .map(|f| f.velocity[0].hypot(f.velocity[1])),
        )
        .map(|v| arrow_length(level, v, m))
        .fold(0.0, f64::max);
    if longest > 0.0 { 3.0 / longest } else { 1.0 }
}

/// Length of a free charge's velocity arrow in cells (`ARROW_MIN` at rest, where only
/// the handle is drawn).
fn arrow_cells(level: &Level, e: &Element, m: ArrowMeasure) -> f64 {
    let v = e.speed.unwrap_or(0.0);
    if v == 0.0 {
        ARROW_MIN
    } else {
        arrow_length(level, v, m) * arrow_scale(level, m)
    }
}

/// Velocity arrow of a free charge, in cell coordinates: its tip, ahead of the particle
/// in the direction of motion (drawn only when it moves).
pub fn arrow_tip(level: &Level, e: &Element, m: ArrowMeasure) -> [f64; 2] {
    let p = level.grid.position(e.node);
    let (sin, cos) = e.angle_deg.to_radians().sin_cos();
    let len = arrow_cells(level, e, m);
    [p.x + cos * len, p.y + sin * len]
}

/// The drag handle of a free charge: behind it, opposite to the arrow (a slingshot: pull
/// the handle back to launch the particle forward).
pub fn arrow_handle(level: &Level, e: &Element, m: ArrowMeasure) -> [f64; 2] {
    let p = level.grid.position(e.node);
    let (sin, cos) = e.angle_deg.to_radians().sin_cos();
    let len = arrow_cells(level, e, m);
    [p.x - cos * len, p.y - sin * len]
}

/// Where the handle of a free charge at rest sits, cells from it.
pub const ARROW_MIN: f64 = 0.8;

/// A drag ending closer than this to the particle (cells) sets it at rest.
pub const REST_ZONE: f64 = 0.4;

/// Index of the player's free charge whose drag handle is within `reach` cells of `x`.
pub fn arrow_handle_at(
    level: &Level,
    placement: &[Element],
    x: [f64; 2],
    reach: f64,
    m: ArrowMeasure,
) -> Option<usize> {
    placement
        .iter()
        .enumerate()
        .filter(|(_, e)| e.kind == ElementKind::Free)
        .map(|(i, e)| {
            let t = arrow_handle(level, e, m);
            (i, (t[0] - x[0]).hypot(t[1] - x[1]))
        })
        .filter(|&(_, d)| d <= reach)
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(i, _)| i)
}

/// The velocity (direction in degrees, speed) that pulling a free charge's handle from
/// the particle at `p` back to `x` asks for (before snapping): towards `p − x`, a slingshot.
pub fn dragged_velocity(level: &Level, p: [f64; 2], x: [f64; 2], m: ArrowMeasure) -> (f64, f64) {
    let (dx, dy) = (p[0] - x[0], p[1] - x[1]);
    let cells = dx.hypot(dy);
    // Close to the particle: at rest.
    let len = if cells < REST_ZONE {
        0.0
    } else {
        cells / arrow_scale(level, m)
    };
    (dy.atan2(dx).to_degrees(), arrow_speed(level, len, m))
}

/// Nearest allowed launch direction of a free charge (15° steps).
pub fn nearest_free_angle(a: f64) -> f64 {
    let a = a.rem_euclid(360.0);
    let dist = |x: f64| {
        let d = (a - x).abs();
        d.min(360.0 - d)
    };
    FREE_ANGLES
        .into_iter()
        .min_by(|x, y| dist(*x).total_cmp(&dist(*y)))
        .expect("angles")
}

/// Nearest allowed discrete plate orientation (a plate turned by 180° is the same).
fn nearest_plate_angle(a: f64) -> f64 {
    let a = a.rem_euclid(180.0);
    let dist = |x: f64| {
        let d = (a - x).abs();
        d.min(180.0 - d)
    };
    PLATE_ANGLES
        .into_iter()
        .min_by(|x, y| dist(*x).total_cmp(&dist(*y)))
        .expect("angles")
}

fn step_i32(step: isize) -> i32 {
    i32::try_from(step).unwrap_or(0)
}

/// Position of a kind in `KINDS`.
pub fn kind_index(kind: ElementKind) -> usize {
    KINDS
        .iter()
        .position(|&k| k == kind)
        .expect("all kinds listed")
}

/// Nearest allowed discrete antenna orientation, and whether the direction flips (an
/// orientation θ + 180° is the same antenna with the opposite sign).
fn nearest_angle(a: f64) -> (f64, bool) {
    let a = a.rem_euclid(360.0);
    let (base, flip) = if a >= 180.0 {
        (a - 180.0, true)
    } else {
        (a, false)
    };
    // Distance on the half circle (0° and 180° are the same line).
    let dist = |x: f64| {
        let d = (base - x).abs();
        d.min(180.0 - d)
    };
    let best = ANTENNA_ANGLES
        .into_iter()
        .min_by(|x, y| dist(*x).total_cmp(&dist(*y)))
        .expect("angles");
    // Wrapping from ~180° to 0° reverses the direction too.
    let wrapped = base - best > 90.0;
    (best, flip ^ wrapped)
}

/// Nearest listed value (on a log scale) to `v`.
fn nearest_listed(list: &[f64], v: f64) -> Option<f64> {
    list.iter()
        .copied()
        .min_by(|a, b| (a.ln() - v.ln()).abs().total_cmp(&(b.ln() - v.ln()).abs()))
}

/// Snaps an element to the level's discrete lists (leaving hardcore mode).
fn snap(level: &Level, e: &mut Element) {
    if is_signed(e.kind) {
        if let Some(v) = nearest_linear(magnitudes(level, e.kind), e.value) {
            e.value = v;
        }
        if e.kind == ElementKind::Plate {
            e.angle_deg = nearest_plate_angle(e.angle_deg);
        }
        // Free charges' velocities stay continuous in both modes.
        return;
    }
    if let Some(m) = nearest_listed(magnitudes(level, e.kind), e.value.abs()) {
        e.value = e.value.signum() * m;
    }
    if e.kind == ElementKind::Antenna {
        let (a, flip) = nearest_angle(e.angle_deg);
        e.angle_deg = a;
        if flip {
            e.value = -e.value;
        }
        e.omega = e
            .omega
            .and_then(|w| nearest_listed(&level.limits.antenna_omegas, w));
    }
}

pub fn describe(e: &PlacementError) -> String {
    match e {
        PlacementError::TooManyCharges => "No charges left for this level.".into(),
        PlacementError::TooManyMagnets => "No magnets left for this level.".into(),
        PlacementError::TooManyAntennas => "No antennas left for this level.".into(),
        PlacementError::TooManyPlates => "No plates left for this level.".into(),
        PlacementError::TooManyFree => "No free charges left for this level.".into(),
        PlacementError::SpeedNotAllowed(_) => "The level does not offer that speed.".into(),
        PlacementError::NoTunableElectrode(_) => {
            "Power supplies belong to the level's tunable electrodes.".into()
        }
        PlacementError::OutsideGrid(_) => "Outside the grid.".into(),
        PlacementError::OutsideRegion(_) => "Elements can only go in the marked region.".into(),
        PlacementError::NotInPlane(_) => "Must be in the plane.".into(),
        PlacementError::Occupied(_) => {
            "Occupied (plates keep 1 cell from other electrodes and stay clear of elements, \
             coils and detectors)."
                .into()
        }
        PlacementError::SignNotAllowed(_) => "That sign is not allowed here.".into(),
        PlacementError::MagnitudeNotAllowed(_) => "That value is not allowed here.".into(),
        PlacementError::FrequencyNotAllowed(_) => {
            "That antenna frequency is not available in this level.".into()
        }
        PlacementError::AngleNotAllowed(_) => {
            "Antennas point along 0°, 45°, 90° or 135°; plates lie along 0° or 90°; free \
             charges fly off in 15° steps."
                .into()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn level() -> Level {
        level::shipped("first_bend")
    }

    #[test]
    fn place_edit_remove() {
        // At most one charge, whatever the shipped level allows.
        let mut l = level();
        l.limits.max_charges = 1;
        let mut e = Editor::new(l);
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

    /// The deflection-plates level with its first plate tunable and room for one player
    /// plate.
    fn plate_level() -> Level {
        let mut l = level::shipped("deflection_plates");
        l.electrodes[0].tunable = true;
        l.limits.supply_voltages = vec![-6e4, -3e4, 0.0, 3e4, 6e4];
        l.limits.max_plates = 1;
        l.limits.plate_voltages = vec![-2e4, 0.0, 2e4];
        l
    }

    #[test]
    fn plates_place_rotate_tune_and_move() {
        let mut e = Editor::new(plate_level());
        e.set_kind(ElementKind::Plate);
        assert_eq!(e.kind, ElementKind::Plate);
        e.set_cursor([21, 10, 0]);
        e.place().unwrap();
        assert_eq!(e.placement, vec![Element::plate([21, 10, 0], -2e4, 0.0)]);
        e.rotate(1);
        assert_eq!(e.placement[0].angle_deg.to_bits(), 90f64.to_bits());
        e.cycle_magnitude(1);
        assert_eq!(e.placement[0].value.to_bits(), 0f64.to_bits());
        e.flip_sign(); // −0 is not listed: rejected, unchanged.
        e.cycle_magnitude(1);
        e.flip_sign();
        assert_eq!(e.placement[0].value.to_bits(), (-2e4f64).to_bits());
        // Picked up by its body (1 cell off centre along its length), it keeps the offset.
        e.set_cursor([21, 11, 0]);
        assert!(e.grab());
        e.move_cursor(-2, 0);
        e.drop_grabbed();
        assert_eq!(e.placement[0].node, [19, 10, 0]);
        // A second plate is over the limit.
        e.set_cursor([24, 4, 0]);
        assert!(matches!(e.place(), Err(PlacementError::TooManyPlates)));
        // Right click on its body removes it.
        e.set_cursor([19, 9, 0]);
        e.remove();
        assert!(e.placement.is_empty());
    }

    #[test]
    fn power_supplies_are_operated_on_their_electrode() {
        let mut e = Editor::new(plate_level());
        let centre = e.level.electrodes[0].center;
        // Anywhere on the tunable electrode: switches the supply on at the listed potential
        // nearest to the electrode's own bias (−30k), then steps it.
        e.set_cursor([centre[0] + 3, centre[1], 0]);
        e.place().unwrap();
        assert_eq!(e.supply(centre), Some(-3e4));
        e.place().unwrap();
        assert_eq!(e.supply(centre), Some(0.0));
        e.place().unwrap();
        e.flip_sign();
        assert_eq!(e.supply(centre), Some(-3e4));
        assert!(!e.grab(), "supplies stay on their electrode");
        e.set_supply(centre, Some(6e4));
        assert_eq!(e.supply(centre), Some(6e4));
        e.remove();
        assert_eq!(e.supply(centre), None);
        // The other electrode is not tunable: a click there places nothing.
        let other = e.level.electrodes[1].center;
        e.set_kind(ElementKind::Charge);
        e.set_cursor(other);
        assert!(e.place().is_err());
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
        // An element on a fine-only node blocks coarsening.
        e.set_refinement(2);
        e.set_cursor([4 * 2, 6 * 2, 0]);
        e.remove();
        e.set_cursor([9, 13, 0]);
        e.place().unwrap();
        e.set_refinement(1);
        assert_eq!(e.refinement(), 2);
    }

    #[test]
    fn grabbed_elements_follow_the_cursor_and_skip_blocked_nodes() {
        let mut e = Editor::new(level());
        e.set_cursor([10, 5, 0]);
        e.place().unwrap();
        let value = e.placement[0].value;
        assert!(e.grab());
        e.move_cursor(2, 1);
        assert_eq!(e.placement[0].node, [12, 6, 0]);
        assert_eq!(e.placement[0].value.to_bits(), value.to_bits());
        // The launch node is blocked: the element waits there and jumps on.
        let launch = e.level.shots[0].launch.node;
        e.set_cursor(launch);
        assert_eq!(e.placement[0].node, [12, 6, 0]);
        e.set_cursor([4, 4, 0]);
        assert_eq!(e.placement[0].node, [4, 4, 0]);
        e.cancel_grab();
        assert_eq!(e.placement[0].node, [10, 5, 0]);
        assert!(e.grabbed.is_none());
        assert_eq!(e.placement.len(), 1);
    }

    /// A level with free charges: speeds 0, 1 and 4 at c = 5.
    fn free_level() -> Level {
        let mut l = level::shipped("first_bend");
        l.limits.max_free = 2;
        l.limits.free_charges = vec![-1e-6, 1e-6];
        l.limits.free_speeds = vec![0.0, 1.0, 4.0];
        l
    }

    /// Dragging the arrow's handle to where the arrow of a speed ends gives back that
    /// speed (for both measures), and `set_velocity` snaps to the listed speeds and 15°
    /// directions.
    #[test]
    fn velocity_arrow_drag_round_trip_and_continuous_setting() {
        let l = free_level();
        for m in [ArrowMeasure::Speed, ArrowMeasure::Rapidity] {
            for v in [1.0, 4.0] {
                let e = Element::free([10, 10, 0], 1e-6, 30.0, v);
                let handle = arrow_handle(&l, &e, m);
                let p = l.grid.position(e.node);
                let (a, back) = dragged_velocity(&l, [p.x, p.y], handle, m);
                assert!((back - v).abs() < 1e-12, "{m:?}: {back} vs {v}");
                assert!((a - 30.0).abs() < 1e-9);
            }
        }
        // Rapidity: the fastest listed speed (0.8 c) gets the full 3 cells, and a speed
        // near c lies far beyond (unbounded).
        let r = |v: f64| arrow_length(&l, v, ArrowMeasure::Rapidity);
        assert!(r(4.999) > 2.5 * r(4.0));
        assert!((r(0.01) - 0.01).abs() < 1e-6, "slow: rapidity is the speed");

        let mut ed = Editor::new(l);
        ed.set_kind(ElementKind::Free);
        ed.set_cursor([10, 10, 0]);
        ed.place().unwrap();
        assert_eq!(ed.placement[0].kind, ElementKind::Free);
        // Continuous (the owner: snapping to the listed values was too clunky): any
        // direction, the speed clamped to the range of the level's speeds.
        ed.set_velocity(0, 37.0, 3.2);
        assert_eq!(ed.placement[0].angle_deg.to_bits(), 37.0_f64.to_bits());
        assert_eq!(ed.placement[0].speed, Some(3.2));
        ed.set_velocity(0, -20.0, 9.0);
        assert_eq!(ed.placement[0].angle_deg.to_bits(), 340.0_f64.to_bits());
        assert_eq!(ed.placement[0].speed, Some(4.0), "clamped to the fastest");
        ed.set_velocity(0, 100.0, 0.0);
        assert_eq!(ed.placement[0].speed, Some(0.0), "at rest");
        assert_eq!(ed.placement[0].angle_deg.to_bits(), 0.0_f64.to_bits());
        // The handle of the placed charge is found where it is drawn: behind it.
        let handle = arrow_handle(&ed.level, &ed.placement[0], ArrowMeasure::Speed);
        let p = ed.level.grid.position(ed.placement[0].node);
        let tip = arrow_tip(&ed.level, &ed.placement[0], ArrowMeasure::Speed);
        assert!(
            (handle[0] + tip[0] - 2.0 * p.x).abs() < 1e-12,
            "opposite the arrow"
        );
        assert_eq!(
            arrow_handle_at(&ed.level, &ed.placement, handle, 0.35, ArrowMeasure::Speed),
            Some(0)
        );
        // Pulling the handle back to the left launches it to the right.
        let (a, _) = dragged_velocity(&ed.level, [p.x, p.y], [p.x - 2.0, p.y], ArrowMeasure::Speed);
        assert!(a.abs() < 1e-9);
    }
}
