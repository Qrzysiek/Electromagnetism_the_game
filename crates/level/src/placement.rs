//! Where the player may put things: the placement check, the supplies' targets,
//! the plates' and the coils' room.

use super::*;

impl Level {
    /// Checks a player placement against the level's limits and occupied nodes.
    pub fn check_placement(&self, player: &[Element]) -> Result<(), PlacementError> {
        let count = |k: ElementKind| player.iter().filter(|e| e.kind == k).count();
        if count(ElementKind::Charge) > self.limits.max_charges as usize {
            return Err(PlacementError::TooManyCharges);
        }
        if count(ElementKind::Magnet) > self.limits.max_magnets as usize {
            return Err(PlacementError::TooManyMagnets);
        }
        if count(ElementKind::Antenna) > self.limits.max_antennas as usize {
            return Err(PlacementError::TooManyAntennas);
        }
        if count(ElementKind::Plate) > self.limits.max_plates as usize {
            return Err(PlacementError::TooManyPlates);
        }
        if count(ElementKind::Free) > self.limits.max_free as usize {
            return Err(PlacementError::TooManyFree);
        }
        if count(ElementKind::Cut) > self.limits.max_cuts as usize {
            return Err(PlacementError::TooManyCuts);
        }
        self.check_supplies(player)?;
        self.check_cuts(player)?;
        // Plates first: point elements are checked against every electrode.
        let boxes = self.check_plates(player)?;
        let mut occupied: Vec<Node> = self.elements.iter().map(|c| c.node).collect();
        occupied.extend(self.shots.iter().map(|s| s.launch.node));
        occupied.extend(self.free_particles.iter().map(|f| f.node));
        for e in player {
            if matches!(
                e.kind,
                ElementKind::Plate | ElementKind::Supply | ElementKind::Cut
            ) {
                continue;
            }
            if !self.grid.contains(e.node) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if self.limits.region.is_some_and(|r| !r.contains(e.node)) {
                return Err(PlacementError::OutsideRegion(e.node));
            }
            if self.grid.is_2d() && e.node[2] != 0 {
                return Err(PlacementError::NotInPlane(e.node));
            }
            if occupied.contains(&e.node) {
                return Err(PlacementError::Occupied(e.node));
            }
            // Not inside or touching a metal sphere.
            let reach = self
                .physics
                .charge_radius
                .max(self.physics.magnet_radius)
                .max(self.physics.antenna_radius);
            let p = self.grid.position(e.node);
            if physics::bem::Electrodes::shapes_only(boxes.clone())
                .contains(p, CONTACT_DISTANCE + reach)
            {
                return Err(PlacementError::Occupied(e.node));
            }
            if self.conductors.iter().any(|c| {
                (p - self.grid.position(c.center)).length() < c.radius + CONTACT_DISTANCE + reach
            }) {
                return Err(PlacementError::Occupied(e.node));
            }
            if self.dielectric_spheres.iter().any(|d| {
                (p - self.grid.position(d.center)).length() < d.radius + CONTACT_DISTANCE + reach
            }) {
                return Err(PlacementError::Occupied(e.node));
            }
            occupied.push(e.node);
            if e.value == 0.0 && e.kind != ElementKind::Free {
                return Err(PlacementError::SignNotAllowed(e.value));
            }
            // Magnitudes come from the level's list and are compared exactly (hardcore:
            // anything within its range).
            let magnitude = e.value.abs();
            match e.kind {
                ElementKind::Charge => {
                    let sign_ok = if e.value > 0.0 {
                        self.limits.allow_positive
                    } else {
                        self.limits.allow_negative
                    };
                    if !sign_ok {
                        return Err(PlacementError::SignNotAllowed(e.value));
                    }
                    if !self.limits.allows(&self.limits.magnitudes, magnitude) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                }
                ElementKind::Magnet => {
                    if !self.limits.allows(&self.limits.magnet_strengths, magnitude) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                }
                ElementKind::Antenna => {
                    if !self
                        .limits
                        .allows(&self.limits.antenna_amplitudes, magnitude)
                    {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                    let angle_ok = if self.limits.continuous {
                        e.angle_deg.is_finite()
                    } else {
                        ANTENNA_ANGLES
                            .iter()
                            .any(|a| a.to_bits() == e.angle_deg.to_bits())
                    };
                    if !angle_ok {
                        return Err(PlacementError::AngleNotAllowed(e.angle_deg));
                    }
                    let omegas = &self.limits.antenna_omegas;
                    let omega_ok = match e.omega {
                        None => omegas.is_empty(),
                        Some(w) => self.limits.allows(omegas, w),
                    };
                    if !omega_ok {
                        return Err(PlacementError::FrequencyNotAllowed(e.omega));
                    }
                }
                ElementKind::Free => {
                    if !self.limits.allows(&self.limits.free_charges, e.value) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                    // A free charge's launch velocity is always continuous: any speed in the
                    // range of the level's speeds, any direction (the owner: dragging the
                    // arrow between listed values was too clunky). The solver still searches
                    // the listed speeds and FREE_ANGLES only.
                    let speed_ok = e.speed.is_some_and(|v| {
                        v.is_finite()
                            && value_range(&self.limits.free_speeds)
                                .is_some_and(|(lo, hi)| (lo..=hi).contains(&v))
                    });
                    if !speed_ok {
                        return Err(PlacementError::SpeedNotAllowed(e.speed));
                    }
                    let angle_ok = e.angle_deg.is_finite();
                    if !angle_ok {
                        return Err(PlacementError::AngleNotAllowed(e.angle_deg));
                    }
                }
                ElementKind::Plate | ElementKind::Supply | ElementKind::Cut => {
                    unreachable!("checked above")
                }
            }
        }
        Ok(())
    }

    /// Power supplies: each on the centre of a different tunable electrode or coil, at a
    /// value on its slider.
    pub(super) fn check_supplies(&self, player: &[Element]) -> Result<(), PlacementError> {
        let targets = self.supply_targets();
        let mut supplied: Vec<Node> = Vec::new();
        for e in player.iter().filter(|e| e.kind == ElementKind::Supply) {
            let Some(&(_, target)) = targets.iter().find(|(c, _)| *c == e.node) else {
                return Err(PlacementError::NoTunableElectrode(e.node));
            };
            if supplied.contains(&e.node) {
                return Err(PlacementError::Occupied(e.node));
            }
            supplied.push(e.node);
            if !Limits::allows_potential(self.supply_list(target), e.value) {
                return Err(PlacementError::MagnitudeNotAllowed(e.value));
            }
        }
        Ok(())
    }

    /// What the player's power supplies operate: the tunable electrodes, each at its
    /// centre, and the tunable circular coils, each at its terminal (`Coil::supply_node`).
    pub fn supply_targets(&self) -> Vec<(Node, SupplyTarget)> {
        let electrodes = self
            .electrodes
            .iter()
            .enumerate()
            .filter(|(_, e)| e.tunable)
            .map(|(i, e)| (e.center, SupplyTarget::Electrode(i)));
        let coils = self
            .coils
            .iter()
            .enumerate()
            .filter_map(|(i, c)| c.supply_node().map(|n| (n, SupplyTarget::Coil(i))));
        let solenoid = self
            .tube
            .as_ref()
            .and_then(|t| t.solenoid)
            .map(|n| (n, SupplyTarget::Solenoid));
        electrodes.chain(coils).chain(solenoid).collect()
    }

    /// The electrodes a supply drives: its tunable electrode and those linked to it through
    /// its tube circuit (a grid's wires, a buncher's second plate, a gun's Wehnelt). None
    /// for a coil or a solenoid.
    pub fn supply_electrodes(&self, target: SupplyTarget) -> Vec<usize> {
        let SupplyTarget::Electrode(i) = target else {
            return Vec::new();
        };
        let mut out = vec![i];
        if let Some(t) = &self.tube {
            for c in t.circuits.iter().filter(|c| c.electrode == i) {
                out.extend(c.linked.iter().copied());
            }
        }
        out
    }

    /// The slider list of a supply's target: a voltage (an electrode's potential, or the
    /// source of a circuit driving a coil: a voltage source too) or a ramp rate (an
    /// undriven coil's).
    pub fn supply_list(&self, target: SupplyTarget) -> &[f64] {
        match target {
            SupplyTarget::Coil(i) if !self.coils[i].is_driven() => &self.limits.coil_rates,
            SupplyTarget::Electrode(_) | SupplyTarget::Coil(_) => &self.limits.supply_voltages,
            SupplyTarget::Solenoid => &self.limits.solenoid_fields,
        }
    }

    /// The player's power supply on the centre `node`, if placed.
    pub(super) fn supply_at<'a>(&self, node: Node, player: &'a [Element]) -> Option<&'a Element> {
        player
            .iter()
            .find(|s| s.kind == ElementKind::Supply && s.node == node)
    }

    /// Player plates: allowed orientation and potential, entirely inside the grid and the
    /// player region, at least `PLATE_CLEARANCE` from every other electrode, clear of
    /// metal spheres, coil wires, elements, launch points and detectors. Returns every
    /// electrode box of the placement (level electrodes first).
    pub(super) fn check_plates(
        &self,
        player: &[Element],
    ) -> Result<Vec<physics::bem::BoxElectrode>, PlacementError> {
        // Electrodes and dielectrics alike: plates keep clear of them, elements out.
        let mut boxes = self.box_electrodes();
        boxes.extend(self.body_boxes());
        let reach = self
            .physics
            .charge_radius
            .max(self.physics.magnet_radius)
            .max(self.physics.antenna_radius);
        let grid_rect = self.node_rect(Region2 {
            min: [0, 0, 0],
            max: self.grid.max_node(),
        });
        for e in player.iter().filter(|e| e.kind == ElementKind::Plate) {
            if !self.grid.contains(e.node) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if self.grid.is_2d() && e.node[2] != 0 {
                return Err(PlacementError::NotInPlane(e.node));
            }
            let angle_ok = if self.limits.continuous {
                e.angle_deg.is_finite()
            } else {
                PLATE_ANGLES
                    .iter()
                    .any(|a| a.to_bits() == e.angle_deg.to_bits())
            };
            if !angle_ok {
                return Err(PlacementError::AngleNotAllowed(e.angle_deg));
            }
            if !Limits::allows_potential(&self.limits.plate_voltages, e.value) {
                return Err(PlacementError::MagnitudeNotAllowed(e.value));
            }
            let b = self.plate_box(e);
            let corners = rect_corners(&b);
            if !corners.iter().all(|c| rect_contains(&grid_rect, *c)) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if let Some(r) = self.limits.region {
                let region = self.node_rect(r);
                if !corners.iter().all(|c| rect_contains(&region, *c)) {
                    return Err(PlacementError::OutsideRegion(e.node));
                }
            }
            let outline = Outline::Rect(corners);
            let blocked = self
                .electrode_outlines(&[])
                .iter()
                .any(|o| o.distance(&outline) < PLATE_CLEARANCE)
                || boxes[self.electrodes.len()..]
                    .iter()
                    .any(|o| rect_distance(&corners, &rect_corners(o)) < PLATE_CLEARANCE)
                || self.shots.iter().any(|s| {
                    rect_distance(
                        &corners,
                        &self.node_rect(Region2 {
                            min: s.detector.min,
                            max: s.detector.max,
                        }),
                    ) <= 0.0
                });
            let one = physics::bem::Electrodes::shapes_only(vec![b]);
            let points = self
                .elements
                .iter()
                .map(|x| (self.grid.position(x.node), reach))
                .chain(
                    self.shots
                        .iter()
                        .map(|s| (self.grid.position(s.launch.node), reach)),
                )
                .chain(self.coil_points().map(|p| (p, self.physics.wire_radius)))
                .chain(
                    self.conductors
                        .iter()
                        .map(|c| (self.grid.position(c.center), c.radius + PLATE_CLEARANCE)),
                )
                // Charge clouds: electrodes see them as point charges, so plates stay out.
                .chain(
                    self.clouds
                        .iter()
                        .map(|c| (self.grid.position(c.center), c.radius + CONTACT_DISTANCE)),
                );
            let touching = points
                .into_iter()
                .any(|(p, r)| one.contains(p, CONTACT_DISTANCE + r));
            if blocked || touching {
                return Err(PlacementError::Occupied(e.node));
            }
            boxes.push(b);
        }
        Ok(boxes)
    }

    /// Points along every coil wire, at most 0.1 cells apart (for clearance checks).
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub(super) fn coil_points(&self) -> impl Iterator<Item = DVec3> + '_ {
        self.coils.iter().flat_map(move |c| {
            let mut out = Vec::new();
            match c {
                Coil::Circle { center, radius, .. } => {
                    let c = self.grid.position(*center);
                    let n = ((std::f64::consts::TAU * radius / 0.1).ceil() as usize).max(8);
                    for i in 0..n {
                        let a = std::f64::consts::TAU * i as f64 / n as f64;
                        out.push(c + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * *radius);
                    }
                }
                Coil::Polygon { vertices, .. } => {
                    for i in 0..vertices.len() {
                        let a = self.grid.position(vertices[i]);
                        let b = self.grid.position(vertices[(i + 1) % vertices.len()]);
                        let n = (((b - a).length() / 0.1).ceil() as usize).max(1);
                        for k in 0..n {
                            out.push(a + (b - a) * (k as f64 / n as f64));
                        }
                    }
                }
            }
            out
        })
    }

    /// The rectangle of a box of nodes, in cells (in-plane corners).
    pub(super) fn node_rect(&self, r: Region2) -> [DVec3; 4] {
        let (a, b) = (self.grid.position(r.min), self.grid.position(r.max));
        [
            DVec3::new(a.x, a.y, 0.0),
            DVec3::new(b.x, a.y, 0.0),
            DVec3::new(b.x, b.y, 0.0),
            DVec3::new(a.x, b.y, 0.0),
        ]
    }
}
