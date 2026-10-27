//! What a tube level's setup cannot compute: its metals' and its tube's issues.

use super::*;

impl Level {
    /// The largest field along z the placement gives (the player's solenoid supply, else the
    /// level's own; a changing field's level is its amplitude).
    pub(super) fn field_peak(&self, spec: &TubeSpec, player: &[Element]) -> f64 {
        spec.solenoid
            .and_then(|n| {
                player
                    .iter()
                    .find(|e| e.kind == ElementKind::Supply && e.node == n)
            })
            .map_or(spec.b_z, |e| e.value)
            .abs()
    }

    /// Whether electrode `e` of a tube level is held by something outside the board: a bias
    /// other than a fixed charge, a supply, a drive or a tube circuit (a floating one keeps
    /// its charge).
    pub(super) fn electrode_held(&self, spec: &TubeSpec, e: usize) -> bool {
        let Some(el) = self.electrodes.get(e) else {
            return false;
        };
        !matches!(el.bias, crate::ConductorBias::Charge(_))
            || el.tunable
            || el.drive.is_some()
            || spec
                .circuits
                .iter()
                .any(|c| c.electrode == e || c.linked.contains(&e))
    }

    /// What a tube level's metals refuse (PHYSICS.md §2.11, the Drude metal): a
    /// rectangle off the board or crossing another metal or an electrode, contacts that
    /// do not lie along a side, a medium of cores too dense for the model, a plasma period
    /// or a mean free path the step does not resolve, too many carriers.
    #[allow(clippy::too_many_lines)]
    pub(super) fn metal_issues(&self, spec: &TubeSpec, player: &[Element], out: &mut Vec<String>) {
        let n = self.electrodes.len();
        let max = self.grid.position(self.grid.max_node());
        let mut carriers = 0.0;
        for (k, m) in spec.metals.iter().enumerate() {
            let (lo, hi) = self.metal_rect(m);
            let name = format!("Metal {}", k + 1);
            if !(hi.x > lo.x && hi.y > lo.y) {
                out.push(format!("{name}: its corners must span a rectangle."));
                continue;
            }
            if lo.x <= 0.0 || lo.y <= 0.0 || hi.x >= max.x || hi.y >= max.y {
                out.push(format!("{name} must lie inside the board."));
            }
            let positive = |v: f64| v.is_finite() && v > 0.0;
            if !(positive(m.density)
                && m.charge_per_mass.is_finite()
                && m.charge_per_mass != 0.0
                && positive(m.thermal_speed)
                && positive(m.mean_free_path)
                && positive(m.core_radius)
                && positive(m.carriers_per_area))
            {
                out.push(format!(
                    "{name}: its density, thermal speed, mean free path, core radius and \
                     carriers per area must be positive, its charge-to-mass ratio nonzero."
                ));
                continue;
            }
            let dense = m.core_radius / (2.0 * m.mean_free_path);
            if dense > MAX_CORE_DENSITY {
                out.push(format!(
                    "{name}: its cores are too dense for the model (a²n = a/2ℓ = {dense:.4}, \
                     at most {MAX_CORE_DENSITY}): make them smaller or the mean free path \
                     longer."
                ));
            }
            let width = (hi.x - lo.x).min(hi.y - lo.y);
            if m.mean_free_path > 0.25 * width {
                out.push(format!(
                    "{name}: its mean free path must be at most a quarter of its width (the \
                     carriers scatter off the cores, not off its walls)."
                ));
            }
            let wp = m.plasma_frequency();
            if wp * spec.step > MAX_PLASMA_STEP {
                out.push(format!(
                    "{name}: the time step does not resolve its plasma oscillation \
                     (ω_p dt = {:.2}, at most {MAX_PLASMA_STEP}).",
                    wp * spec.step
                ));
            }
            if 4.0 * m.thermal_speed * spec.step > m.mean_free_path {
                out.push(format!(
                    "{name}: the time step must be at most a quarter of the mean free time \
                     (ℓ over 4 thermal speeds)."
                ));
            }
            // The step must follow the carriers' gyration: at most a radian per step (the
            // Boris rotation's angle is then within 8 % of the true one; the verification's
            // half step within 2 %).
            let turn = m.charge_per_mass.abs() * self.field_peak(spec, player) * spec.step;
            if turn > 1.0 {
                out.push(format!(
                    "{name}: the field turns its carriers by {turn:.2} rad per step, at most 1 \
                     (the step cannot follow their gyration): a weaker field."
                ));
            }
            let area = (hi.x - lo.x) * (hi.y - lo.y);
            carriers += area * m.carriers_per_area;
            #[allow(clippy::cast_precision_loss)]
            let cores = area * m.core_density();
            if cores > MAX_CORES {
                out.push(format!("{name}: {cores:.0} cores, at most {MAX_CORES:.0}."));
            }
            // An isolated metal (eddy currents): every electrode floating. One held at a
            // potential in this open z-invariant world draws charge from the log gauge's
            // "infinity", and a floating conductor beside it runs away (docs/OPEN_PROBLEMS.md,
            // OP-1; PHYSICS.md 2.11).
            if m.contacts.is_empty() && (0..n).any(|e| self.electrode_held(spec, e)) {
                out.push(format!(
                    "{name} has no contact: an isolated metal needs every electrode floating \
                     (one held at a potential draws charge from this open world's infinity)."
                ));
            }
            for (i, &e) in m.contacts.iter().enumerate() {
                if e >= n || m.contacts[..i].contains(&e) {
                    out.push(format!("{name}: its contacts must be distinct electrodes."));
                } else if self.contact_side(m, e).is_none() {
                    out.push(format!(
                        "{name}: electrode {} is not a contact: a box along one of its sides, \
                         its face within the side.",
                        e + 1
                    ));
                }
                if spec.metals[..k].iter().any(|o| o.contacts.contains(&e)) {
                    out.push(format!("Electrode {} is a contact of two metals.", e + 1));
                }
            }
            if spec.metal_of(spec.cathode) == Some(k)
                && (spec.charge_per_mass - m.charge_per_mass).abs()
                    > 1e-12 * m.charge_per_mass.abs()
            {
                out.push(format!(
                    "{name}: the cathode is its contact: their carriers' charge-to-mass ratios \
                     must agree."
                ));
            }
            if spec.ions.is_some_and(|i| m.contacts.contains(&i.electrode)) {
                out.push(format!("{name}: the ion source cannot be its contact."));
            }
            // Nothing inside it: electrodes, other metals, projectiles, free charges.
            let overlaps =
                |a0: DVec3, a1: DVec3| a0.x < hi.x && a1.x > lo.x && a0.y < hi.y && a1.y > lo.y;
            for e in 0..n {
                let inside = match self.electrode_box(e) {
                    Some((b0, b1)) => overlaps(b0, b1),
                    None => {
                        let c = self.grid.position(self.electrodes[e].center);
                        let r = 0.5 * (self.electrodes[e].length + self.electrodes[e].thickness);
                        overlaps(c - DVec3::new(r, r, 0.0), c + DVec3::new(r, r, 0.0))
                    }
                };
                if inside {
                    out.push(format!("{name} overlaps electrode {}.", e + 1));
                }
            }
            for (j, o) in spec.metals.iter().enumerate().skip(k + 1) {
                let (a0, a1) = self.metal_rect(o);
                if overlaps(a0, a1) {
                    out.push(format!("{name} overlaps metal {}.", j + 1));
                }
            }
            let inside = |x: DVec3| x.x >= lo.x && x.x <= hi.x && x.y >= lo.y && x.y <= hi.y;
            if spec
                .projectiles
                .iter()
                .any(|p| inside(self.grid.position(p.node)))
                || player
                    .iter()
                    .filter(|e| e.kind == ElementKind::Free)
                    .any(|e| inside(self.grid.position(e.node)))
            {
                out.push(format!(
                    "{name}: no projectile or free charge may start inside it."
                ));
            }
        }
        if 2.0 * carriers * spec.particles > MAX_METAL_CARRIERS {
            out.push(format!(
                "Tube: the metals' carriers ({:.0} in the verification) exceed {MAX_METAL_CARRIERS:.0}.",
                2.0 * carriers * spec.particles
            ));
        }
    }

    /// What a tube level refuses: everything that breaks z-invariance, and an
    /// inconsistent spec.
    #[allow(clippy::too_many_lines)]
    pub(crate) fn tube_issues(&self, player: &[Element], out: &mut Vec<String>) {
        let Some(spec) = &self.tube else { return };
        let n = self.electrodes.len();
        if spec.cathode >= n {
            out.push("Tube: the cathode is not one of the level's electrodes.".into());
        }
        if let Some(c) = spec.speed_of_light {
            if !super::SPEED_OF_LIGHT_RANGE.contains(&c) {
                out.push(format!(
                    "Tube: the speed of light {c} lies outside {}–{}.",
                    super::SPEED_OF_LIGHT_RANGE.start(),
                    super::SPEED_OF_LIGHT_RANGE.end()
                ));
            }
            if spec.step > super::SELF_FIELD_MAX_STEP {
                out.push(format!(
                    "Tube: with the carriers' own field the step is at most {} (its splitting is \
                     first order, PHYSICS.md §2.12), not {}.",
                    super::SELF_FIELD_MAX_STEP,
                    spec.step
                ));
            }
            let max = self.grid.position(self.grid.max_node());
            let d = super::self_field_grid(max.x, max.y, c);
            if d.nx > super::SELF_FIELD_MAX_NODES || d.ny > super::SELF_FIELD_MAX_NODES {
                out.push(format!(
                    "Tube: the board is too large for the carriers' own field ({} × {} nodes, \
                     at most {} a side).",
                    d.nx,
                    d.ny,
                    super::SELF_FIELD_MAX_NODES
                ));
            }
        }
        for g in &spec.goals {
            // A carrier's own emitter takes back only those that return: not a goal.
            let emitter = match g.carriers {
                Carriers::Cathode => Some(spec.cathode),
                Carriers::Ions => spec.ions.map(|i| i.electrode),
                Carriers::Net => None,
            };
            if g.carriers == Carriers::Ions && spec.ions.is_none() {
                out.push("Tube: an ion current goal needs an ion source.".into());
            } else if g.electrode >= n || Some(g.electrode) == emitter {
                out.push(
                    "Tube: a current goal must be an electrode other than its carriers' \
                     emitter."
                        .into(),
                );
            }
            if !(g.min < g.max && 0.0 <= g.start && g.start < g.end && g.end.is_finite()) {
                out.push("Tube: a current goal needs min < max and 0 ≤ start < end.".into());
            }
        }
        if spec.goals.is_empty()
            && !spec.projectiles.iter().any(|p| p.detector.is_some())
            && spec.potential_goals.is_empty()
            && spec.heat_goals.is_empty()
        {
            out.push(
                "Tube: no goal (a current, a potential, or a projectile with a detector).".into(),
            );
        }
        if !(spec.charge_per_mass.is_finite() && spec.charge_per_mass != 0.0) {
            out.push("Tube: the carriers' charge-to-mass ratio must be finite and nonzero.".into());
        }
        self.metal_issues(spec, player, out);
        if let Some(ions) = spec.ions
            && (ions.electrode >= n
                || ions.electrode == spec.cathode
                || !(ions.charge_per_mass.is_finite() && ions.charge_per_mass > 0.0))
        {
            out.push(
                "Tube: the ion source must be another electrode, with a positive charge-to-mass \
                 ratio."
                    .into(),
            );
        }
        if spec
            .saturation
            .is_some_and(|j| !(j.is_finite() && j >= 0.0))
        {
            out.push(
                "Tube: the cathode's saturation current must not be negative (0: it emits \
                 nothing)."
                    .into(),
            );
        }
        if !(spec.particle_rate.is_finite() && spec.particle_rate >= 1.0) {
            out.push("Tube: the cathode's macroparticle rate must be at least 1.".into());
        }
        if !(spec.particles.is_finite() && (0.1..=100.0).contains(&spec.particles)) {
            out.push("Tube: the macroparticle factor must be between 0.1 and 100.".into());
        }
        if !(spec.emitter_mesh.is_finite() && (1.0..=16.0).contains(&spec.emitter_mesh)) {
            out.push("Tube: the emitters' mesh factor must be between 1 and 16.".into());
        }
        if let Some(o) = spec.output {
            let load = spec
                .circuits
                .get(o.circuit)
                .and_then(|c| c.resistance)
                .is_some_and(|r| r > 0.0 && r.is_finite());
            if !load || !(o.c.is_finite() && o.c > 0.0) {
                out.push(
                    "Tube: the output needs a circuit with a resistance (its load) and a \
                     positive wave speed."
                        .into(),
                );
            }
        }
        if spec.potential_goals.iter().any(|g| {
            g.reference
                .is_some_and(|r| r >= self.electrodes.len() || r == g.electrode)
        }) {
            out.push(
                "Tube: a potential goal's reference must be another of the level's electrodes."
                    .into(),
            );
        }
        // Every electrode in at most one circuit (as its electrode or linked to it).
        let mut taken: Vec<usize> = Vec::new();
        for (k, c) in spec.circuits.iter().enumerate() {
            let mine: Vec<usize> = std::iter::once(c.electrode)
                .chain(c.linked.iter().copied())
                .collect();
            let twice = mine
                .iter()
                .enumerate()
                .any(|(i, e)| taken.contains(e) || mine[..i].contains(e) || *e >= n)
                || c.return_to.is_some_and(|r| c.linked.contains(&r));
            // A coil needs a finite resistance (an open branch carries no current), and a
            // branch to another electrode a resistance or a coil (an ideal source has no
            // branch).
            let coil_ok = c.inductance.is_finite()
                && c.inductance >= 0.0
                && (c.inductance == 0.0 || c.resistance.is_some_and(f64::is_finite));
            let return_ok = c.return_to.is_none_or(|r| {
                r < n
                    && r != c.electrode
                    && (c.inductance > 0.0
                        || c.resistance.is_some_and(|x| x > 0.0 && x.is_finite()))
            });
            taken.extend(mine);
            let sane = c.electrode < n
                && !twice
                && c.source.is_valid()
                && c.resistance.is_none_or(|r| r.is_finite() && r >= 0.0)
                && c.capacitance.is_finite()
                && c.capacitance >= 0.0
                && (!c.supply_offset || matches!(c.source, crate::Source::Sine { .. }))
                && coil_ok
                && return_ok;
            if !sane {
                out.push(format!(
                    "Tube: circuit {} needs electrodes of its own, a valid source, a \
                     resistance, capacitance and inductance that are not negative (an offset \
                     supply needs a sine; a coil a finite resistance; a return to another \
                     electrode, not one linked to it, a resistance or a coil).",
                    k + 1
                ));
            }
        }
        for (k, g) in spec.potential_goals.iter().enumerate() {
            if !(g.electrode < n && g.min < g.max && 0.0 <= g.start && g.start < g.end) {
                out.push(format!(
                    "Tube: potential goal {} needs an electrode, min < max and 0 ≤ start < end.",
                    k + 1
                ));
            }
        }
        for (k, g) in spec.heat_goals.iter().enumerate() {
            if !(g.metal < spec.metals.len() && g.min < g.max && 0.0 <= g.start && g.start < g.end)
            {
                out.push(format!(
                    "Tube: heat goal {} needs a metal, min < max and 0 ≤ start < end.",
                    k + 1
                ));
            }
        }
        if spec.solenoid.is_some_and(|n| !self.grid.contains(n)) {
            out.push("Tube: the solenoid's terminal must be on the board.".into());
        }
        if !(spec.b_z.is_finite() && spec.step.is_finite() && spec.step > 0.0) {
            out.push("Tube: B_z must be finite and the step positive.".into());
        }
        for (j, p) in spec.projectiles.iter().enumerate() {
            let inside = self.grid.contains(p.node)
                && p.detector
                    .is_none_or(|[a, b]| self.grid.contains(a) && self.grid.contains(b));
            let sane = p.mass > 0.0
                && p.mass.is_finite()
                && p.charge.is_finite()
                && p.speed.is_finite()
                && p.speed >= 0.0
                && p.launch.is_finite()
                && p.launch >= 0.0;
            if !(inside && sane) {
                out.push(format!(
                    "Tube: projectile {} needs a node and detector on the board, a positive \
                     mass and a finite charge, speed and launch time.",
                    j + 1
                ));
            }
        }
        let breaks = [
            (!self.shots.is_empty(), "shots"),
            (
                !self.elements.is_empty(),
                "charges, magnets and other elements",
            ),
            (!self.coils.is_empty(), "coils"),
            (!self.conductors.is_empty(), "metal spheres"),
            (!self.clouds.is_empty(), "charge clouds"),
            (!self.free_particles.is_empty(), "free particles"),
            (
                !self.dielectrics.is_empty() || !self.dielectric_spheres.is_empty(),
                "dielectrics",
            ),
            (!self.ferrites.is_empty(), "ferrites"),
            (!self.disturbances.is_empty(), "disturbances"),
            (!self.gates.is_empty(), "gates"),
            (
                self.electrodes.iter().any(|e| e.drive.is_some()),
                "driven electrodes",
            ),
            (
                player.iter().any(|e| {
                    !matches!(
                        e.kind,
                        ElementKind::Plate
                            | ElementKind::Supply
                            | ElementKind::Free
                            | ElementKind::Cut
                    )
                }),
                "placed elements other than plates, supplies, free charges and cuts",
            ),
        ];
        // Newtonian: every carrier (falling through the largest potential difference)
        // and projectile (its launch speed plus that fall) below 0.1 c, where the
        // neglected v²/c² corrections (1 %) are below the macroparticles' error (5–8 %
        // between the two runs).
        let potential = |b: &physics::bem::BoxElectrode| match b.bias {
            physics::conductor::Bias::Potential(v) => v,
            _ => 0.0,
        };
        let boxes = self.all_box_electrodes(player);
        if let Some(c) = self.physics.c {
            let (s_lo, s_hi) = self.tube_source_extremes(player);
            let (lo, hi) = boxes
                .iter()
                .map(potential)
                .fold((s_lo, s_hi), |(a, b), v| (a.min(v), b.max(v)));
            let dv = hi - lo;
            let ions = spec.ions.map_or(0.0, |i| i.charge_per_mass);
            let free_mass = self.limits.free_mass.max(1e-300);
            let ratios = spec
                .projectiles
                .iter()
                .map(|p| (p.charge / p.mass).abs())
                .chain(
                    player
                        .iter()
                        .filter(|e| e.kind == ElementKind::Free)
                        .map(|e| (e.value / free_mass).abs()),
                )
                .chain([spec.charge_per_mass.abs(), ions]);
            let fall = ratios.fold(0.0_f64, f64::max);
            let launch = spec
                .projectiles
                .iter()
                .map(|p| p.speed)
                .chain(
                    player
                        .iter()
                        .filter(|e| e.kind == ElementKind::Free)
                        .map(|e| e.speed.unwrap_or(0.0)),
                )
                .fold(0.0, f64::max);
            let v = libm::sqrt(2.0 * fall * dv + launch * launch);
            if v > 0.1 * c {
                out.push(format!(
                    "Tube: the particles would reach {:.2} c; the tube model is Newtonian \
                     (up to 0.1 c): lower the voltages and speeds or leave c infinite.",
                    v / c
                ));
            }
        }
        // Every electrode inside the board: the envelope lies just outside it. The
        // electrodes' own outlines (a rod's circle, a ring's annulus, open inside).
        let outlines = self.electrode_outlines(player);
        let max = self.grid.position(self.grid.max_node());
        let board = [
            DVec3::ZERO,
            DVec3::new(max.x, 0.0, 0.0),
            DVec3::new(max.x, max.y, 0.0),
            DVec3::new(0.0, max.y, 0.0),
        ];
        for (i, o) in outlines.iter().enumerate() {
            if !o.within(&board) {
                out.push(format!(
                    "Tube: electrode {} reaches beyond the board (the tube's grounded \
                     envelope lies just outside it).",
                    i + 1
                ));
            }
        }
        // Prisms must not touch: the boundary elements of touching electrodes at different
        // potentials would hold an infinite field between them.
        for (i, a) in outlines.iter().enumerate() {
            for (j, b) in outlines.iter().enumerate().skip(i + 1) {
                if a.distance(b) < crate::CONTACT_DISTANCE {
                    out.push(format!(
                        "Electrodes {} and {} touch or overlap.",
                        i + 1,
                        j + 1
                    ));
                }
            }
        }
        // A projectile or a free charge starting inside metal.
        let in_metal = |x: DVec3| {
            outlines
                .iter()
                .any(|o| o.distance_to(x) < crate::CONTACT_DISTANCE)
        };
        for (j, p) in spec.projectiles.iter().enumerate() {
            if in_metal(self.grid.position(p.node)) {
                out.push(format!("Tube: projectile {} starts in metal.", j + 1));
            }
        }
        for (j, e) in player
            .iter()
            .filter(|e| e.kind == ElementKind::Free)
            .enumerate()
        {
            if in_metal(self.grid.position(e.node)) {
                out.push(format!("Tube: your free charge {} starts in metal.", j + 1));
            }
        }
        for (bad, what) in breaks {
            if bad {
                out.push(format!(
                    "Tube levels are z-invariant: {what} break the symmetry and are not supported there."
                ));
            }
        }
    }
}
