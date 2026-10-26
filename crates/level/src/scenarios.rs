//! The flights: each shot's scenario at a resolution (the display's, the
//! verification's) and the tolerances they are judged with.

use super::*;

impl Level {
    /// Metal spheres with the charges the fixed charges induce on them.
    pub(super) fn conductors_for(
        &self,
        charges: &[FixedCharge],
        resolution: Resolution,
    ) -> Conductors {
        if self.conductors.is_empty() {
            return Conductors::default();
        }
        let spheres = self
            .conductors
            .iter()
            .map(|c| SphereConductor {
                center: self.grid.position(c.center),
                radius: c.radius,
                bias: match c.bias {
                    ConductorBias::Grounded => Bias::Grounded,
                    ConductorBias::Charge(q) => Bias::Charge(q),
                    ConductorBias::Potential(v) => Bias::Potential(v),
                },
            })
            .collect();
        let sources: Vec<(DVec3, f64)> = charges.iter().map(|c| (c.position, c.charge)).collect();
        Conductors::new(spheres, &sources, resolution)
    }

    /// One flight's scenario (shot, disturbance) with metal at display resolution (for
    /// pictures: potential map, field lines, field views).
    pub fn display_scenario(
        &self,
        shot: usize,
        disturbance: usize,
        player: &[Element],
    ) -> Scenario<LevelField> {
        let (field, obstacles) = self.field_at(player, Resolution::Display);
        self.scenario_with(shot, disturbance, field, obstacles)
    }

    /// Scenarios of all flights with the field at `resolution`.
    pub fn scenarios_at(
        &self,
        player: &[Element],
        resolution: Resolution,
    ) -> Vec<Scenario<LevelField>> {
        let (field, obstacles) = self.field_at(player, resolution);
        (0..self.flight_count())
            .map(|i| {
                let (shot, d) = self.flight_of(i);
                self.scenario_with(shot, d, field.clone(), obstacles.clone())
            })
            .collect()
    }

    /// Verifies every flight: preview at preview resolution, the tighter run with the
    /// field at verification resolution (they differ only with metal: spheres, electrodes,
    /// player plates), so the field model's error enters the verdict too.
    pub fn verify_flights(&self, player: &[Element]) -> Vec<Verification> {
        let preview = self.scenarios(player);
        let fine = if self.has_metal(player) {
            self.scenarios_at(player, Resolution::Verify)
        } else {
            preview.clone()
        };
        preview
            .iter()
            .zip(&fine)
            .map(|(a, b)| verify_pair(a, b, self.tolerances()))
            .collect()
    }

    /// Flights per shot: one per disturbance, or one undisturbed flight.
    pub fn flights_per_shot(&self) -> usize {
        self.disturbances.len().max(1)
    }

    /// Total number of flights, `shots × flights_per_shot`.
    pub fn flight_count(&self) -> usize {
        self.shots.len() * self.flights_per_shot()
    }

    /// Shot and disturbance of flight `i` (flights are ordered shot-major).
    pub fn flight_of(&self, i: usize) -> (usize, usize) {
        (i / self.flights_per_shot(), i % self.flights_per_shot())
    }

    /// The physical scenario of one shot under its first disturbance (or undisturbed).
    pub fn scenario(&self, shot: usize, player: &[Element]) -> Scenario<LevelField> {
        let (field, obstacles) = self.field(player);
        self.scenario_with(shot, 0, field, obstacles)
    }

    /// Scenarios of all flights, shot-major (see `flight_of`), sharing one field.
    pub fn scenarios(&self, player: &[Element]) -> Vec<Scenario<LevelField>> {
        let (field, obstacles) = self.field(player);
        (0..self.flight_count())
            .map(|i| {
                let (shot, d) = self.flight_of(i);
                self.scenario_with(shot, d, field.clone(), obstacles.clone())
            })
            .collect()
    }

    pub(super) fn scenario_with(
        &self,
        shot: usize,
        disturbance: usize,
        mut field: LevelField,
        obstacles: Vec<Shape>,
    ) -> Scenario<LevelField> {
        if let Some(d) = self.disturbances.get(disturbance) {
            field.external = d.fields(self.c());
        }
        field.time_offset = self.shots[shot].launch.time;
        Scenario {
            field,
            obstacles,
            particle: self.particle(shot),
            c: self.c(),
            x0: self.grid.position(self.shots[shot].launch.node),
            p0: self.launch_momentum(shot),
            detector: Some(self.detector_region(shot)),
            bounds: Some(self.bounds()),
            t_max: self.physics.t_max,
            radiation_reaction: self.physics.radiation_reaction,
            acceptance: self.shots[shot]
                .detector
                .acceptance
                .map(DetectorAcceptance::to_physics),
            gates: (0..self.gates.len())
                .map(|i| physics::trajectory::Gate {
                    region: self.gate_region(i),
                    acceptance: self.gates[i].acceptance.map(DetectorAcceptance::to_physics),
                })
                .collect(),
        }
    }

    pub fn tolerances(&self) -> Tolerances {
        self.physics.tolerances.into()
    }
}
