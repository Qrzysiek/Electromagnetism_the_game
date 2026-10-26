//! What a level's setup and its model cannot compute (`setup_issues`,
//! `model_issues` and their parts), and the bounds on the neglected image force
//! and the moment's radiation.

use super::*;

/// In-plane distance from `p` to a coil's wire (its centre line).
pub(super) fn coil_wire_distance(level: &Level, coil: &Coil, p: DVec3) -> f64 {
    let p = DVec3::new(p.x, p.y, 0.0);
    match coil {
        Coil::Circle { center, radius, .. } => {
            let c = level.grid.position(*center);
            ((p - DVec3::new(c.x, c.y, 0.0)).length() - radius).abs()
        }
        Coil::Polygon { vertices, .. } => {
            let n = vertices.len();
            (0..n)
                .map(|i| {
                    let a = level.grid.position(vertices[i]);
                    let b = level.grid.position(vertices[(i + 1) % n]);
                    let (a, b) = (DVec3::new(a.x, a.y, 0.0), DVec3::new(b.x, b.y, 0.0));
                    let ab = b - a;
                    let len2 = ab.length_squared();
                    let s = if len2 > 0.0 {
                        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    (p - (a + ab * s)).length()
                })
                .fold(f64::INFINITY, f64::min)
        }
    }
}

/// Largest image force the electrodes would exert on the particle (neglected, PHYSICS.md
/// §2.7) at the points `(x, t)` of a flight, bounded by `q²/d²` (four times the force of
/// a flat grounded plane at distance `d`, which covers concave corners), relative to
/// `max(|qE|, F₀)` with `F₀ = T₀` per cell, the smallest force that matters for the
/// flight (relative to the Lorentz force alone the ratio is meaningless where that force
/// passes through zero). 0 without electrodes.
pub fn electrode_image_force_bound(
    scn: &Scenario<LevelField>,
    kinetic_initial: f64,
    points: impl IntoIterator<Item = (DVec3, f64)>,
) -> f64 {
    use physics::field::FieldSolver;
    if scn.field.electrodes.is_empty() {
        return 0.0;
    }
    let boxes = scn.field.electrodes.obstacles(0.0);
    let q = scn.particle.charge;
    let mut worst: f64 = 0.0;
    for (x, t) in points {
        let d = boxes
            .iter()
            .map(|b| b.signed_distance(x))
            .fold(f64::INFINITY, f64::min)
            .max(1e-3);
        let f = (scn.field.sample(x, t).e.length() * q.abs()).max(kinetic_initial);
        worst = worst.max(q * q / (d * d) / f);
    }
    worst
}

/// Estimated energy radiated by the particle's magnetic moment `m` along the states
/// `(x, p, t)` of a flight (neglected, PHYSICS.md §3.2). A moving moment carries the
/// electric dipole `v×m/c²`, which radiates `2|ȧ×m|²/(3c⁷)` (`ȧ`: the jerk); the magnetic
/// quadrupole term is of the same order, so the estimate is `m²|ȧ|²/c⁷`, with the jerk
/// from the force between consecutive states. 0 for `c = ∞` or without a moment.
pub fn moment_radiation_estimate(
    scn: &Scenario<LevelField>,
    states: impl IntoIterator<Item = (DVec3, DVec3, f64)>,
) -> f64 {
    let (m, c) = (scn.particle.moment, scn.c);
    if m == 0.0 || !c.is_finite() {
        return 0.0;
    }
    let ode = physics::dynamics::ParticleOde::new(&scn.field, &scn.particle, c, 1.0);
    let accel =
        |x: DVec3, p: DVec3, t: f64| ode.force(x, p, t) / (ode.kin.gamma(p) * scn.particle.mass);
    let mut w = 0.0;
    let mut prev: Option<(DVec3, f64)> = None;
    for (x, p, t) in states {
        let a = accel(x, p, t);
        if let Some((a0, t0)) = prev
            && t > t0
        {
            let jerk = (a - a0) / (t - t0);
            w += m * m * jerk.length_squared() / c.powi(7) * (t - t0);
        }
        prev = Some((a, t));
    }
    w
}

/// Largest neglected image force that is still below the numerical accuracy (relative to
/// the force that matters, see `electrode_image_force_bound`).
pub const IMAGE_FORCE_LIMIT: f64 = 1e-10;

impl Level {
    /// Geometry the physics cannot compute, which the level editor's free fields allow:
    /// dielectrics touching metal or each other (the boundary elements assume separate
    /// bodies), metal spheres overlapping or touching each other or an electrode (the image
    /// series and the boundary elements assume separate bodies; such setups failed or
    /// ran for minutes), and particles that start inside or touching metal or on a coil's
    /// wire (the field is singular there: NaN forces). Found by `scripts/pairs.py`. The
    /// game shows these instead of computing; the shipped levels have none.
    pub fn setup_issues(&self, player: &[Element]) -> Vec<String> {
        let mut out = Vec::new();
        if !self.is_tube() {
            for (i, e) in self.electrodes.iter().enumerate() {
                if e.shape != ElectrodeShape::Box {
                    out.push(format!(
                        "Electrode {}: rods and rings exist in tube levels only (the 3D \
                         levels' boundary elements mesh boxes).",
                        i + 1
                    ));
                }
            }
        }
        let pos = |n: Node| self.grid.position(n);
        let boxes = physics::bem::Electrodes::shapes_only(self.all_box_electrodes(player));
        for (i, a) in self.conductors.iter().enumerate() {
            for (j, b) in self.conductors.iter().enumerate().skip(i + 1) {
                if (pos(a.center) - pos(b.center)).length() < a.radius + b.radius + CONTACT_DISTANCE
                {
                    out.push(format!(
                        "Metal spheres {} and {} touch or overlap.",
                        i + 1,
                        j + 1
                    ));
                }
            }
            if boxes.contains(pos(a.center), a.radius + CONTACT_DISTANCE) {
                out.push(format!("Metal sphere {} touches an electrode.", i + 1));
            }
        }
        // Dielectrics: apart from metal and from each other.
        let dielectrics = self.body_boxes();
        let metal_boxes = self.all_box_electrodes(player);
        for (i, d) in dielectrics.iter().enumerate() {
            let corners = rect_corners(d);
            let touches = |o: &physics::bem::BoxElectrode| {
                rect_distance(&corners, &rect_corners(o)) < CONTACT_DISTANCE
            };
            if metal_boxes.iter().any(touches) {
                out.push(format!("Dielectric {} touches an electrode.", i + 1));
            }
            if dielectrics[i + 1..].iter().any(touches) {
                out.push(format!("Dielectric {} touches another dielectric.", i + 1));
            }
            let one = physics::bem::Electrodes::shapes_only(vec![*d]);
            if self
                .conductors
                .iter()
                .any(|c| one.contains(pos(c.center), c.radius + CONTACT_DISTANCE))
            {
                out.push(format!("Dielectric {} touches a metal sphere.", i + 1));
            }
        }
        // Dielectric spheres: apart from metal, from the boxes and from each other.
        let all_boxes = physics::bem::Electrodes::shapes_only(
            metal_boxes.iter().chain(&dielectrics).copied().collect(),
        );
        for (i, d) in self.dielectric_spheres.iter().enumerate() {
            let c = pos(d.center);
            let reach = d.radius + CONTACT_DISTANCE;
            let near_sphere = |q: DVec3, r: f64| (c - q).length() < reach + r;
            if all_boxes.contains(c, reach)
                || self
                    .conductors
                    .iter()
                    .any(|m| near_sphere(pos(m.center), m.radius))
                || self.dielectric_spheres[i + 1..]
                    .iter()
                    .any(|o| near_sphere(pos(o.center), o.radius))
            {
                out.push(format!(
                    "Dielectric sphere {} touches metal or another dielectric.",
                    i + 1
                ));
            }
        }
        let dielectric_shapes = physics::bem::Electrodes::shapes_only(dielectrics);
        let wire = self.physics.wire_radius + CONTACT_DISTANCE;
        let starts = self
            .shots
            .iter()
            .enumerate()
            .map(|(i, s)| (format!("Shot {}", i + 1), s.launch.node, 0.0))
            .chain(self.free_particles.iter().enumerate().map(|(i, f)| {
                (
                    format!("Free particle {}", i + 1),
                    f.node,
                    f.particle.radius,
                )
            }))
            .chain(
                player
                    .iter()
                    .filter(|e| e.kind == ElementKind::Free)
                    .enumerate()
                    .map(|(i, e)| {
                        (
                            format!("Your free charge {}", i + 1),
                            e.node,
                            self.limits.free_radius,
                        )
                    }),
            );
        for (who, node, radius) in starts {
            let p = pos(node);
            let in_sphere = self
                .conductors
                .iter()
                .any(|c| (p - pos(c.center)).length() < c.radius + radius + CONTACT_DISTANCE);
            if in_sphere || boxes.contains(p, radius + CONTACT_DISTANCE) {
                out.push(format!("{who} starts inside or touching metal."));
            }
            let in_dielectric_sphere = self
                .dielectric_spheres
                .iter()
                .any(|d| (p - pos(d.center)).length() < d.radius + radius + CONTACT_DISTANCE);
            if in_dielectric_sphere || dielectric_shapes.contains(p, radius + CONTACT_DISTANCE) {
                out.push(format!("{who} starts inside or touching a dielectric."));
            }
            if self
                .coils
                .iter()
                .any(|c| coil_wire_distance(self, c, p) < wire + radius)
            {
                out.push(format!("{who} starts on a coil's wire."));
            }
        }
        self.tube_issues(player, &mut out);
        out
    }

    /// Radiation goals (PHYSICS.md §3.4) need a finite `c` and sit on shot detectors only.
    /// Where several particles fly (free particles, free charges), a goal measures all of
    /// them together, their far fields added; a beam shot's particles arrive one by one, so
    /// a beam shot has none.
    pub(super) fn radiation_goal_issues(&self, out: &mut Vec<String>) {
        let has = |d: &Detector| d.acceptance.is_some_and(|a| a.radiation.is_some());
        if self.shots.iter().any(|s| has(&s.detector)) {
            if self.physics.c.is_none() {
                out.push("a radiation goal needs a finite speed of light".into());
            }
            if self
                .shots
                .iter()
                .any(|s| s.beam.is_some() && has(&s.detector))
            {
                out.push(
                    "a beam shot cannot have a radiation goal (its particles arrive one by one)"
                        .into(),
                );
            }
            if self.has_beams()
                && self.shots.iter().any(|s| {
                    s.detector
                        .acceptance
                        .and_then(|a| a.radiation)
                        .is_some_and(|r| r.abrupt_stop)
                })
            {
                out.push(
                    "an abrupt stop is measured on single flights only (not with other particles \
                     flying)"
                        .into(),
                );
            }
        }
        if self.shots.iter().any(|s| {
            s.detector
                .acceptance
                .and_then(|a| a.radiation)
                .is_some_and(|r| r.abrupt_stop && r.band.is_none())
        }) {
            out.push(
                "a radiation goal with an abrupt stop needs a band (the stop's spectrum is flat to infinite frequency)"
                    .into(),
            );
        }
        for (i, s) in self.shots.iter().enumerate() {
            let Some(acc) = s.detector.acceptance else {
                continue;
            };
            let Some(r) = acc.radiation else { continue };
            let Some([t1, t2]) = r.window else { continue };
            if !(t1.is_finite() && t2.is_finite() && 0.0 <= t1 && t1 < t2) {
                out.push(format!(
                    "shot {}: a steady receiver's window must run forward from t ≥ 0",
                    i + 1
                ));
            } else if t2 >= self.physics.t_max {
                out.push(format!(
                    "shot {}: a steady receiver's window must end before the time limit (the \
                     light time from the scatterer to the origin is part of it)",
                    i + 1
                ));
            }
            if r.abrupt_stop {
                out.push(format!(
                    "shot {}: a steady receiver measures no stop (the particle need not arrive)",
                    i + 1
                ));
            }
            if acc.direction.is_some() || acc.kinetic.is_some() {
                out.push(format!(
                    "shot {}: a steady receiver's detector has no conditions on the arrival \
                     (the particle need not arrive)",
                    i + 1
                ));
            }
        }
        if self.gates.iter().any(has) {
            out.push("gates cannot have radiation goals".into());
        }
        if self
            .free_particles
            .iter()
            .any(|f| f.detector.as_ref().is_some_and(has))
        {
            out.push("free particles' detectors cannot have radiation goals".into());
        }
    }

    /// Problems of the physical model of this level (combinations that are not
    /// supported), as messages. Empty if the level is consistent.
    /// Positions and radii of the level's bodies and launch points for the containment
    /// checks of `model_issues`: fixed elements with the radius of their kind (supplies
    /// and plates have rules of their own), the shots' and the free particles' launch
    /// points with the particle's radius (contact during the flight uses it too).
    pub(super) fn bodies(&self) -> impl Iterator<Item = (DVec3, f64)> + '_ {
        let elements = self.elements.iter().filter_map(|e| {
            let r = match e.kind {
                ElementKind::Charge => self.physics.charge_radius,
                ElementKind::Magnet => self.physics.magnet_radius,
                ElementKind::Antenna => self.physics.antenna_radius,
                ElementKind::Free => self.limits.free_radius,
                ElementKind::Plate | ElementKind::Supply | ElementKind::Cut => return None,
            };
            Some((self.grid.position(e.node), r))
        });
        let shots = self
            .shots
            .iter()
            .map(|s| (self.grid.position(s.launch.node), s.particle.radius));
        let free = self
            .free_particles
            .iter()
            .map(|f| (self.grid.position(f.node), f.particle.radius));
        elements.chain(shots).chain(free)
    }

    pub fn model_issues(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.drive_issues(&mut out);
        self.gate_issues(&mut out);
        self.radiation_goal_issues(&mut out);
        if self.has_beams() {
            self.beam_issues(&mut out);
        }
        // Magnetic moments: of the shots and of the free particles.
        let moments = self
            .shots
            .iter()
            .map(|s| s.particle.moment)
            .chain(self.free_particles.iter().map(|f| f.particle.moment))
            .any(|m| m != 0.0);
        if moments {
            // m ∇B_z is exact for E = 0, or for c = ∞ (PHYSICS.md §3.2). At finite c an
            // electric field adds velocity-dependent (Aharonov–Casher / hidden-momentum)
            // terms that are not modelled.
            let electric = self.elements.iter().any(|e| e.kind != ElementKind::Magnet)
                || self.limits.max_charges > 0
                || self.limits.max_antennas > 0
                || self.limits.max_plates > 0
                || self.limits.max_free > 0
                || !self.free_particles.is_empty()
                || !self.clouds.is_empty()
                || !self.conductors.is_empty()
                || !self.electrodes.is_empty()
                || self
                    .disturbances
                    .iter()
                    .any(|d| d.e != [0.0, 0.0] || !d.waves.is_empty());
            if self.physics.c.is_some() && electric {
                out.push(
                    "magnetic moments with electric fields at finite c need terms that are not \
                     modelled (use c = ∞, or magnetic fields only)"
                        .into(),
                );
            }
            if self.physics.radiation_reaction {
                out.push("radiation reaction is not modelled for magnetic moments".into());
            }
        }
        // A launch inside its own detector would arrive at t = 0.
        let in_detector = (0..self.shots.len()).any(|i| {
            self.detector_region(i)
                .signed_distance(self.grid.position(self.shots[i].launch.node))
                <= 0.0
        }) || self.free_particles.iter().any(|f| {
            f.detector.is_some_and(|d| {
                self.box_region(d)
                    .signed_distance(self.grid.position(f.node))
                    <= 0.0
            })
        });
        if in_detector {
            out.push("a launch point is inside its own detector".into());
        }
        if self.physics.radiation_reaction && self.physics.c.is_some() {
            // The metal responds electrostatically: its image follows the particle at
            // once and the radiation it would reflect is missing, which near its surface
            // changes the reaction by a factor of order one (PHYSICS.md §3.1).
            if self.has_metal(&[]) || self.limits.max_plates > 0 {
                out.push(
                    "radiation reaction near metal: the metal's electrostatic response does \
                     not reflect the radiation, which changes the reaction near it (not \
                     modelled)"
                        .into(),
                );
            }
            // A drained charge's field jumps where the light cone of its removal passes;
            // the reaction differentiates the field along the world line.
            if self.physics.instant_drain && self.interacts() {
                out.push(
                    "radiation reaction with an instant drain: the drained charge's field \
                     jumps where its light cone passes, which the reaction cannot follow"
                        .into(),
                );
            }
        }
        if self.limits.max_plates > 0 && self.limits.plate_voltages.is_empty() {
            out.push("player plates need at least one allowed potential".into());
        }
        let tunable: Vec<Node> = self
            .electrodes
            .iter()
            .filter(|e| e.tunable)
            .map(|e| e.center)
            .collect();
        if !tunable.is_empty() && self.limits.supply_voltages.is_empty() {
            out.push("tunable electrodes need at least one supply voltage".into());
        }
        let targets: Vec<Node> = self.supply_targets().iter().map(|t| t.0).collect();
        if (1..targets.len()).any(|i| targets[..i].contains(&targets[i])) {
            out.push("two tunable electrodes or coils share a centre node".into());
        }
        let tunable = |driven: bool| {
            self.coils
                .iter()
                .any(|c| matches!(c, Coil::Circle { tunable: true, .. }) && c.is_driven() == driven)
        };
        if tunable(false) && self.limits.coil_rates.is_empty() {
            out.push("tunable coils need at least one ramp rate".into());
        }
        if tunable(true) && self.limits.supply_voltages.is_empty() {
            out.push("tunable driven coils need at least one supply voltage".into());
        }
        if !self.ferrites.is_empty() {
            // Magnetized by the static fields only; a changing magnetization would also
            // induce an electric field, not modelled yet.
            if self.coils.iter().any(Coil::is_time_dependent) {
                out.push(
                    "ferrites respond to static magnets and coils only: no ramped or driven coils"
                        .into(),
                );
            }
            if !self.disturbances.is_empty() {
                out.push("ferrites do not respond to outside fields yet".into());
            }
            let moments = self.shots.iter().any(|s| s.particle.moment != 0.0)
                || self.free_particles.iter().any(|f| f.particle.moment != 0.0);
            if moments {
                out.push("ferrites and magnetic moments are not yet solved together".into());
            }
        }
        let dielectric = !self.dielectrics.is_empty() || !self.dielectric_spheres.is_empty();
        if dielectric {
            // The bound charge answers the fixed charges (as the electrodes' surface charge
            // does), not the metal spheres' images nor fields changing in time.
            if !self.conductors.is_empty() {
                out.push("metal spheres and dielectrics are not yet solved together".into());
            }
            let time_dependent = self.limits.max_antennas > 0
                || self.elements.iter().any(|e| e.kind == ElementKind::Antenna)
                || !self.disturbances.is_empty();
            if time_dependent {
                out.push(
                    "dielectrics respond to static charges only: no antennas or outside fields"
                        .into(),
                );
            }
        }
        if !self.electrodes.is_empty() || self.limits.max_plates > 0 {
            if !self.conductors.is_empty() {
                out.push("metal spheres and box electrodes are not yet solved together".into());
            }
            let time_dependent = self.limits.max_antennas > 0
                || self.elements.iter().any(|e| e.kind == ElementKind::Antenna)
                || !self.disturbances.is_empty();
            if time_dependent {
                out.push(
                    "electrodes respond electrostatically; antennas and disturbances would \
                     need a full-wave solution"
                        .into(),
                );
            }
            let boxes = physics::bem::Electrodes::shapes_only(self.box_electrodes());
            let inside = self
                .bodies()
                .any(|(x, r)| boxes.contains(x, CONTACT_DISTANCE + r));
            if inside {
                out.push("an element or launch point is inside or at an electrode".into());
            }
        }
        // Ramped and driven coils: their induced field is not conservative, which the metal
        // model (electrostatic) cannot screen; and time-dependent B is not combined with
        // magnetic moments (like antennas and waves).
        if self.coils.iter().any(Coil::is_time_dependent) {
            if !self.conductors.is_empty()
                || !self.electrodes.is_empty()
                || self.limits.max_plates > 0
            {
                out.push(
                    "ramped coils with metal: its electrostatic response cannot screen the \
                     induced field"
                        .into(),
                );
            }
            if moments {
                out.push("ramped coils with magnetic moments are not modelled".into());
            }
        }
        let c = self.c();
        let too_fast = self
            .free_particles
            .iter()
            .any(|f| libm::hypot(f.velocity[0], f.velocity[1]) >= c)
            || (self.limits.max_free > 0
                && self.limits.free_speeds.iter().any(|&v| v >= c || v < 0.0));
        if too_fast {
            out.push("free particles need speeds below c (and not negative)".into());
        }
        for cl in &self.clouds {
            let c = self.grid.position(cl.center);
            if cl.radius.is_nan() || cl.radius <= 0.0 || !cl.charge.is_finite() {
                out.push("a charge cloud needs a positive radius and a finite charge".into());
            }
            let touches_sphere = self.conductors.iter().any(|m| {
                (self.grid.position(m.center) - c).length()
                    <= m.radius + cl.radius + CONTACT_DISTANCE
            });
            let boxes = physics::bem::Electrodes::shapes_only(self.box_electrodes());
            if touches_sphere || boxes.contains(c, cl.radius + CONTACT_DISTANCE) {
                out.push(
                    "a charge cloud must stay clear of metal (the metal sees it as a point charge)"
                        .into(),
                );
            }
        }
        if !self.conductors.is_empty() {
            let time_dependent = self.limits.max_antennas > 0
                || self.elements.iter().any(|e| e.kind == ElementKind::Antenna)
                || !self.disturbances.is_empty();
            if time_dependent {
                out.push(
                    "metal spheres respond electrostatically; antennas and disturbances \
                     would need a full-wave solution"
                        .into(),
                );
            }
            for (i, a) in self.conductors.iter().enumerate() {
                for b in &self.conductors[i + 1..] {
                    let d = (self.grid.position(a.center) - self.grid.position(b.center)).length();
                    if d <= a.radius + b.radius + 2.0 * CONTACT_DISTANCE {
                        out.push("two metal spheres touch or overlap".into());
                    }
                }
                let c = self.grid.position(a.center);
                let inside = self
                    .bodies()
                    .any(|(x, r)| (x - c).length() < a.radius + CONTACT_DISTANCE + r);
                if inside {
                    out.push("an element or launch point is inside or at a metal sphere".into());
                }
            }
        }
        out
    }

    /// Circuits (PHYSICS.md §2.10): valid parameters, chains that fit their element, a
    /// circuit that can be solved, and the flights' integration (the beam runner does not
    /// stop at the circuit's breakpoints and knots).
    pub(super) fn drive_issues(&self, out: &mut Vec<String>) {
        if !self.has_drives() {
            return;
        }
        let before = out.len();
        if self.has_beams() {
            out.push("circuits are not combined with beams or free particles yet".into());
        }
        for e in &self.electrodes {
            let Some(d) = e.drive else { continue };
            if !d.is_valid() {
                out.push(
                    "a drive needs a positive resistance, valid source parameters and                      non-negative times, inductance and capacitance"
                        .into(),
                );
            }
            if matches!(e.bias, ConductorBias::Charge(_)) {
                out.push("a driven electrode starts at a potential, not floating".into());
            }
            if d.capacitance != 0.0 {
                out.push("an electrode's drive has no series capacitance (it is one)".into());
            }
            if d.switch.is_some() && d.inductance != 0.0 {
                out.push(
                    "a switch cannot open on an inductor's current (no inductance with a switch)"
                        .into(),
                );
            }
        }
        for c in &self.coils {
            let (Coil::Circle {
                drive: Some(d),
                rate,
                ..
            }
            | Coil::Polygon {
                drive: Some(d),
                rate,
                ..
            }) = c
            else {
                continue;
            };
            if !d.is_valid() {
                out.push(
                    "a drive needs a positive resistance, valid source parameters and                      non-negative times, inductance and capacitance"
                        .into(),
                );
            }
            if *rate != 0.0 {
                out.push("a driven coil has no ramp".into());
            }
            if d.inductance != 0.0 || d.switch.is_some() {
                out.push(
                    "a coil's drive has no series inductance (it is one) and no switch".into(),
                );
            }
        }
        if self.coils.iter().any(Coil::is_driven) && self.physics.c.is_none() {
            out.push("a driven coil needs a finite speed of light (its inductance is N/c²)".into());
        }
        if out.len() == before {
            let (field, _) = self.field_at(&self.reference_solution, Resolution::Preview);
            if let Err(e) = self.drives_for(&self.reference_solution, &field) {
                out.push(format!(
                    "the circuit cannot be solved ({e:?}: a loop of capacitors and sources, \
                     or a switch on an inductor's current)"
                ));
            }
        }
    }

    /// Gates: they must not overlap (each other) or contain a launch point.
    pub(super) fn gate_issues(&self, out: &mut Vec<String>) {
        if self.gates.is_empty() {
            return;
        }
        let rect = |d: &Detector| {
            let (a, b) = (self.grid.position(d.min), self.grid.position(d.max));
            (a.min(b), a.max(b))
        };
        let overlap = |(a0, a1): (DVec3, DVec3), (b0, b1): (DVec3, DVec3)| {
            a0.x <= b1.x && b0.x <= a1.x && a0.y <= b1.y && b0.y <= a1.y
        };
        for (i, g) in self.gates.iter().enumerate() {
            if self.gates[i + 1..]
                .iter()
                .any(|h| overlap(rect(g), rect(h)))
            {
                out.push("two gates overlap".into());
            }
            let (lo, hi) = rect(g);
            if self.shots.iter().any(|s| {
                let p = self.grid.position(s.launch.node);
                (lo.x..=hi.x).contains(&p.x) && (lo.y..=hi.y).contains(&p.y)
            }) {
                out.push("a launch point lies in a gate".into());
            }
        }
    }

    /// Beams (PHYSICS.md §3.3): what the beam runner models.
    pub(super) fn beam_issues(&self, out: &mut Vec<String>) {
        if self.interacts() {
            if self.has_metal(&[]) || self.limits.max_plates > 0 {
                out.push(
                    "interacting beams with metal: the charge one particle induces would act \
                     on the others (not modelled)"
                        .into(),
                );
            }
            // Interacting particles launched at the same point would have an infinite
            // interaction: a beam needs a position spread, and a single shot must not
            // share its launch point with other particles.
            let zero_spread = |s: &Shot| s.beam.is_none_or(|b| b.width == 0.0 && b.length == 0.0);
            let crowded = self.shots.iter().enumerate().any(|(i, a)| {
                zero_spread(a)
                    && (a.beam.is_some_and(|b| b.count > 1)
                        || self
                            .shots
                            .iter()
                            .enumerate()
                            .any(|(j, b)| j != i && b.launch.node == a.launch.node))
            });
            if crowded {
                out.push(
                    "interacting particles launched at the same point (give the beam a position spread)"
                        .into(),
                );
            }
            // Opposite point charges fall into each other (the classical point-charge
            // theory breaks down there); with a radius they collide as rigid spheres.
            let point_charges: Vec<f64> = self
                .shots
                .iter()
                .map(|s| s.particle)
                .chain(self.free_particles.iter().map(|f| f.particle))
                .filter(|p| p.radius == 0.0)
                .map(|p| p.charge)
                .chain(
                    if self.limits.max_free > 0 && self.limits.free_radius == 0.0 {
                        self.limits.free_charges.clone()
                    } else {
                        Vec::new()
                    },
                )
                .collect();
            let positive = point_charges.iter().any(|&q| q > 0.0);
            let negative = point_charges.iter().any(|&q| q < 0.0);
            if positive && negative {
                out.push(
                    "opposite point charges in one flight could fall into each other (give \
                     them a radius: they then collide as rigid spheres)"
                        .into(),
                );
            }
        }
        if self.shots.iter().any(|s| s.launch.time != 0.0) {
            out.push(
                "shots flying together (beams, or with dynamic particles) are launched at t = 0"
                    .into(),
            );
        }
        for s in &self.shots {
            if let Some(b) = s.beam {
                if b.count == 0 || !(b.transmission > 0.0 && b.transmission <= 1.0) {
                    out.push("a beam needs particles and a transmission in (0, 1]".into());
                }
                if b.energy_spread * b.reach() >= 1.0 || b.energy_spread < 0.0 {
                    out.push("a beam's energy spread would give non-positive energies".into());
                }
                if [b.angle_spread_deg, b.width, b.length]
                    .iter()
                    .any(|v| *v < 0.0)
                {
                    out.push("a beam's spreads must not be negative".into());
                }
            }
        }
    }
}
