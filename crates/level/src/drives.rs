//! The time-varying drives (PHYSICS.md §2.10): each element's circuit, their time
//! scales, and the bounds on their back-action and on the quasi-static
//! approximation.

use super::*;

impl Level {
    /// The player's power supply of a tunable electrode, if placed.
    pub(super) fn supply_of<'a>(
        &self,
        e: &Electrode,
        player: &'a [Element],
    ) -> Option<&'a Element> {
        player
            .iter()
            .find(|s| s.kind == ElementKind::Supply && e.tunable && s.node == e.center)
    }

    /// Whether a circuit drives an electrode or a coil (PHYSICS.md §2.10).
    pub fn has_drives(&self) -> bool {
        self.electrodes.iter().any(|e| e.drive.is_some()) || self.coils.iter().any(Coil::is_driven)
    }

    /// The circuits driving the level's electrodes and coils, solved over every flight's
    /// lab time (PHYSICS.md §2.10), for a field built by `field_at` (its electrodes, the
    /// level's then the player's, its loops and polygons, the coils of each shape in
    /// order). None without drives.
    pub fn drives_for(
        &self,
        player: &[Element],
        field: &LevelField,
    ) -> Result<Option<physics::drive::Drives>, physics::circuit::CircuitError> {
        let electrodes: Vec<(usize, physics::drive::Chain)> = self
            .electrodes
            .iter()
            .enumerate()
            .filter_map(|(i, e)| {
                let d = e.drive?;
                let source = match self.supply_of(e, player) {
                    Some(s) => d.source.with_level(s.value),
                    None => d.source,
                };
                Some((i, d.chain(source)))
            })
            .collect();
        let loops: Vec<(usize, physics::drive::Chain)> = self
            .coils
            .iter()
            .filter(|c| matches!(c, Coil::Circle { .. }))
            .enumerate()
            .filter_map(|(i, c)| match c {
                Coil::Circle { drive: Some(d), .. } => {
                    let source = match c.supply_node().and_then(|n| self.supply_at(n, player)) {
                        Some(s) => d.source.with_level(s.value),
                        None => d.source,
                    };
                    Some((i, d.chain(source)))
                }
                _ => None,
            })
            .collect();
        let polygons: Vec<(usize, physics::drive::Chain)> = self
            .coils
            .iter()
            .filter(|c| matches!(c, Coil::Polygon { .. }))
            .enumerate()
            .filter_map(|(i, c)| match c {
                Coil::Polygon { drive: Some(d), .. } => Some((i, d.chain(d.source))),
                _ => None,
            })
            .collect();
        if electrodes.is_empty() && loops.is_empty() && polygons.is_empty() {
            return Ok(None);
        }
        physics::drive::Drives::build(
            &field.electrodes,
            &electrodes,
            &field.loops,
            &loops,
            &field.polygons,
            &polygons,
            self.c(),
            self.drive_end(),
        )
        .map(Some)
    }

    /// The lab time the circuits are solved over: up to the last launch plus `t_max`.
    pub fn drive_end(&self) -> f64 {
        self.shots.iter().map(|s| s.launch.time).fold(0.0, f64::max) + self.physics.t_max
    }

    /// The largest magnitude each circuit quantity that drives the field reaches over the
    /// circuit's time, sampled at 2001 times: `(electrode potentials, coil strengths κ)`,
    /// in the order of `Drives::electrodes`, then `Drives::loops` followed by
    /// `Drives::polygons`.
    pub(super) fn drive_scales(&self, d: &physics::drive::Drives) -> (Vec<f64>, Vec<f64>) {
        let t_end = self.drive_end();
        let peak = |f: &dyn Fn(f64) -> f64| {
            (0..=2000)
                .map(|k| f(t_end * f64::from(k) / 2000.0).abs())
                .fold(0.0, f64::max)
        };
        (
            d.electrodes
                .iter()
                .map(|&(_, n)| peak(&|t| d.solution.potential(n, t)))
                .collect(),
            d.loops
                .iter()
                .chain(&d.polygons)
                .map(|&(_, k)| peak(&|t| d.solution.current(k, t)) * d.per_current)
                .collect(),
        )
    }

    /// A bound on the particles' neglected back-action on the circuits (one-way coupling,
    /// PHYSICS.md §2.10), relative to the drive. A charge q induces at most |q| on a plate
    /// of the circuit, changing its potential by at most |q|/C (C its coefficient of
    /// capacitance), against the largest potential the circuit's plates reach. A charge
    /// moving slower than c and kept a wire radius b from a driven coil of radius a links
    /// at most the flux |q|·2πa/(c b) (its vector potential is |q|v/(c² r)), which changes
    /// the coil's strength by at most that over its self-flux N (κ = I/c², L = N/c²),
    /// against the largest κ of the driven coils. A circuit that stays at 0 drives
    /// nothing: then only the charge the particle induces acts, which is its image force,
    /// bounded separately (§2.7). 0 without circuits.
    pub fn drive_back_action_bound(&self, player: &[Element]) -> f64 {
        let (field, _) = self.field(player);
        let Some(d) = &field.drives else {
            return 0.0;
        };
        let q = self
            .shots
            .iter()
            .map(|s| s.particle.charge.abs())
            .fold(0.0, f64::max);
        let (potentials, strengths) = self.drive_scales(d);
        let v = potentials.iter().fold(0.0_f64, |m, v| m.max(*v));
        let kappa = strengths.iter().fold(0.0_f64, |m, k| m.max(*k));
        let mut bound: f64 = 0.0;
        if let Some(cm) = field.electrodes.capacitance()
            && v > 0.0
        {
            for &(e, _) in &d.electrodes {
                bound = bound.max(q / (cm[e][e] * v));
            }
        }
        if kappa > 0.0 {
            for &(i, _) in &d.loops {
                let l = &field.loops[i];
                let n = physics::inductance::ring_self(l);
                bound = bound.max(
                    q * std::f64::consts::TAU * l.radius / (self.c() * l.wire_radius * n * kappa),
                );
            }
            // A polygon: the same with its perimeter for the ring's 2πa.
            for &(i, _) in &d.polygons {
                let p = &field.polygons[i];
                let n = physics::inductance::polygon_self(p);
                let k = p.vertices.len();
                let perimeter: f64 = (0..k)
                    .map(|j| (p.vertices[(j + 1) % k] - p.vertices[j]).length())
                    .sum();
                bound = bound.max(q * perimeter / (self.c() * p.wire_radius * n * kappa));
            }
        }
        bound
    }

    /// The circuits' quasi-static parameter (PHYSICS.md §2.10): the light time across the
    /// arena times the fastest relative rate of change of a driven potential or coil
    /// strength (finite differences on 2000 intervals). 0 without circuits or for c = ∞.
    pub fn drive_quasi_static_ratio(&self, player: &[Element]) -> f64 {
        let (field, _) = self.field(player);
        let Some(d) = &field.drives else {
            return 0.0;
        };
        let t_end = self.drive_end();
        let dt = t_end / 2000.0;
        let (potentials, strengths) = self.drive_scales(d);
        let rate = |f: &dyn Fn(f64) -> f64, scale: f64| {
            (0..2000)
                .map(|k| {
                    let t = dt * f64::from(k);
                    (f(t + dt) - f(t)).abs() / dt
                })
                .fold(0.0, f64::max)
                / scale
        };
        let mut fastest: f64 = 0.0;
        for (&(_, n), v) in d.electrodes.iter().zip(&potentials) {
            if *v > 0.0 {
                fastest = fastest.max(rate(&|t| d.solution.potential(n, t), *v));
            }
        }
        for (&(_, k), kappa) in d.loops.iter().chain(&d.polygons).zip(&strengths) {
            if *kappa > 0.0 {
                fastest = fastest.max(rate(&|t| d.solution.current(k, t) * d.per_current, *kappa));
            }
        }
        f64::from(self.grid.nx.max(self.grid.ny)) / self.c() * fastest
    }
}
