//! Computational cost of a level and its budgets (the sandbox's resource meters,
//! `generator check`). A level whose flights or linear systems are too expensive is slow
//! to play, to verify and to analyse, and leaves the range where the solvers were
//! tested; the budgets keep levels interactive and accurate.
//!
//! Times are measured wall-clock times on the machine at hand; memory, steps and the
//! boundary residual do not depend on the machine.

use std::time::Instant;

use physics::conductor::Resolution;
use physics::field::{LevelField, SetupCost};
use physics::trajectory::{RunSettings, Trajectory, run};

use crate::{Element, Level};

/// Measured cost of one setup (level and player elements).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cost {
    /// Flights per placement (shots × disturbances).
    pub flights: usize,
    /// Linear systems of metal spheres and electrodes at display, preview and
    /// verification resolution: build time of a new geometry, memory, unknowns.
    pub setup: SetupCost,
    /// Wall time of all preview flights so far: the delay after every change.
    pub preview_seconds: f64,
    pub previews: usize,
    /// Wall time of all verification flights so far.
    pub verify_seconds: f64,
    pub verifications: usize,
    /// Integrator steps (accepted and rejected) of the longest flight.
    pub max_steps: u64,
    /// Force evaluations of the preview flights (for the mean cost per evaluation).
    pub evaluations: u64,
    /// Metal boundary residual at verification resolution (PHYSICS.md §2.6); NaN without
    /// metal spheres.
    pub metal_residual: f64,
}

impl Cost {
    pub fn new(level: &Level) -> Self {
        Self {
            // A tube level has one run per resolution (its current), not flights.
            flights: if level.is_tube() {
                1
            } else {
                level.flight_count()
            },
            setup: SetupCost::default(),
            preview_seconds: 0.0,
            previews: 0,
            verify_seconds: 0.0,
            verifications: 0,
            max_steps: 0,
            evaluations: 0,
            metal_residual: f64::NAN,
        }
    }

    /// Adds the build cost of a field's linear systems (call once per resolution).
    pub fn add_field(&mut self, field: &LevelField) {
        self.setup = self.setup + field.setup_cost();
    }

    /// Adds the field at verification resolution: its build cost and the accuracy of the
    /// metal model.
    pub fn add_verification_field(&mut self, field: &LevelField) {
        self.add_field(field);
        if !field.conductors.is_empty() {
            let sources: Vec<_> = field.coulomb.charges().collect();
            self.metal_residual = field.conductors.boundary_residual(&sources);
        }
    }

    pub fn add_preview(&mut self, seconds: f64, traj: &Trajectory) {
        self.preview_seconds += seconds;
        self.previews += 1;
        self.evaluations += traj.stats.n_fcn;
        self.max_steps = self.max_steps.max(traj.stats.n_step);
    }

    /// A tube level's run (`Level::tube_run`): its wall time (which includes building its
    /// segments' linear system), its steps, and that system (a dense n × n
    /// factorization, 8n² bytes).
    pub fn add_tube_run(&mut self, seconds: f64, segments: usize, steps: u32, verification: bool) {
        let n = segments;
        self.setup = self.setup
            + SetupCost {
                seconds: 0.0,
                bytes: 8 * n * n,
                unknowns: n,
            };
        self.max_steps = self.max_steps.max(u64::from(steps));
        if verification {
            self.verify_seconds += seconds;
            self.verifications += 1;
        } else {
            self.preview_seconds += seconds;
            self.previews += 1;
            self.evaluations += u64::from(steps);
        }
    }

    pub fn add_verification(&mut self, seconds: f64, traj: &Trajectory) {
        self.verify_seconds += seconds;
        self.verifications += 1;
        self.max_steps = self.max_steps.max(traj.stats.n_step);
    }

    /// Mean wall time of one force evaluation in the preview flights.
    #[allow(clippy::cast_precision_loss)] // counts far below 2⁵²
    pub fn seconds_per_evaluation(&self) -> f64 {
        if self.evaluations == 0 {
            f64::NAN
        } else {
            self.preview_seconds / self.evaluations as f64
        }
    }

    /// The resource meters with their budgets. Preview and verification times are
    /// extrapolated to all flights while flights are still being computed.
    #[allow(clippy::cast_precision_loss)] // counts far below 2⁵²
    pub fn meters(&self) -> Vec<Meter> {
        let all = |seconds: f64, done: usize| {
            if done == 0 {
                f64::NAN
            } else {
                seconds * self.flights as f64 / done as f64
            }
        };
        let mut m = vec![
            Meter {
                label: "Preview",
                value: all(self.preview_seconds, self.previews),
                good: PREVIEW_GOOD,
                limit: PREVIEW_LIMIT,
                unit: Unit::Seconds,
                hint: "Wall time of all preview flights: the delay after every change. \
                       Grows with the number of shots and disturbances, the flight time \
                       and the cost of the field (sources, metal, electrodes).",
            },
            Meter {
                label: "Verification",
                value: all(self.verify_seconds, self.verifications),
                good: VERIFY_GOOD,
                limit: VERIFY_LIMIT,
                unit: Unit::Seconds,
                hint: "Wall time of all verification flights (tighter tolerance, finer \
                       metal and electrode model) until the verdict appears.",
            },
            Meter {
                label: "Steps",
                value: if self.previews == 0 {
                    f64::NAN
                } else {
                    self.max_steps as f64
                },
                good: STEPS_GOOD,
                limit: STEPS_LIMIT,
                unit: Unit::Count,
                hint: "Integrator steps of the longest flight. At 10⁶ steps a flight \
                       fails; the budget leaves room for other placements. Many steps \
                       come from long flight times, many turns in magnetic fields and \
                       close passes by charges.",
            },
        ];
        if self.setup.unknowns > 0 {
            m.push(Meter {
                label: "Metal setup",
                value: self.setup.seconds,
                good: SETUP_GOOD,
                limit: SETUP_LIMIT,
                unit: Unit::Seconds,
                hint: "Time to build the linear systems of the metal spheres and \
                       electrodes (all resolutions). Paid once per geometry: when the \
                       level loads and whenever metal is moved or resized.",
            });
            m.push(Meter {
                label: "Memory",
                value: self.setup.bytes as f64,
                good: MEMORY_GOOD,
                limit: MEMORY_LIMIT,
                unit: Unit::Bytes,
                hint: "Memory of the metal and electrode matrices (grows with the square \
                       of the surface area).",
            });
        }
        if !self.metal_residual.is_nan() {
            m.push(Meter {
                label: "Metal error",
                value: self.metal_residual,
                good: RESIDUAL_GOOD,
                limit: RESIDUAL_LIMIT,
                unit: Unit::Ratio,
                hint: "Largest deviation of the metal surfaces from their potential at \
                       verification resolution, relative to the potentials there \
                       (PHYSICS.md §2.6). Spheres very close to each other or to charges \
                       need more resolution than the solver has; shipped levels stay \
                       below 1e-10.",
            });
        }
        m
    }

    /// Meters over their limit.
    pub fn over_budget(&self) -> Vec<Meter> {
        self.meters()
            .into_iter()
            .filter(|m| m.load() == Load::Over)
            .collect()
    }
}

// Budgets. Preview: the game should answer a change within a fraction of a second.
// Verification: the verdict within seconds. Setup: moving metal must not freeze the
// editor for long. Steps: the integrator's hard limit is 10⁶ (RunSettings), the budget
// leaves room for other placements. Memory:
// well within a typical machine. Metal error: the threshold of the level tests.
pub const PREVIEW_GOOD: f64 = 0.2;
pub const PREVIEW_LIMIT: f64 = 1.0;
pub const VERIFY_GOOD: f64 = 3.0;
pub const VERIFY_LIMIT: f64 = 20.0;
pub const STEPS_GOOD: f64 = 1e4;
pub const STEPS_LIMIT: f64 = 2e5;
pub const SETUP_GOOD: f64 = 1.0;
pub const SETUP_LIMIT: f64 = 5.0;
pub const MEMORY_GOOD: f64 = 64e6;
pub const MEMORY_LIMIT: f64 = 512e6;
pub const RESIDUAL_GOOD: f64 = 1e-10;
pub const RESIDUAL_LIMIT: f64 = 1e-8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    Seconds,
    Bytes,
    Count,
    Ratio,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Load {
    Fine,
    High,
    Over,
}

/// One resource with its budget: fine up to `good`, high up to `limit`, over beyond.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Meter {
    pub label: &'static str,
    /// NaN while not yet measured.
    pub value: f64,
    pub good: f64,
    pub limit: f64,
    pub unit: Unit,
    pub hint: &'static str,
}

impl Meter {
    pub fn load(&self) -> Load {
        if self.value.is_nan() || self.value <= self.good {
            Load::Fine
        } else if self.value <= self.limit {
            Load::High
        } else {
            Load::Over
        }
    }

    /// Position on a logarithmic bar from `good/1000` (empty) to `3·limit` (full).
    pub fn fraction(&self) -> f64 {
        Self::position(self.value, self.good, self.limit)
    }

    /// Positions of the `good` and `limit` marks on the bar.
    pub fn marks(&self) -> (f64, f64) {
        (
            Self::position(self.good, self.good, self.limit),
            Self::position(self.limit, self.good, self.limit),
        )
    }

    fn position(v: f64, good: f64, limit: f64) -> f64 {
        let (lo, hi) = (good / 1000.0, 3.0 * limit);
        if v.is_nan() || v <= lo {
            return 0.0;
        }
        (libm::log(v / lo) / libm::log(hi / lo)).clamp(0.0, 1.0)
    }

    pub fn text(&self) -> String {
        format_value(self.value, self.unit)
    }
}

pub fn format_value(v: f64, unit: Unit) -> String {
    if v.is_nan() {
        return "…".into();
    }
    match unit {
        Unit::Seconds if v < 1e-5 => format!("{:.1} µs", v * 1e6),
        Unit::Seconds if v < 1e-3 => format!("{:.0} µs", v * 1e6),
        Unit::Seconds if v < 1.0 => format!("{:.0} ms", v * 1e3),
        Unit::Seconds => format!("{v:.1} s"),
        Unit::Bytes if v < 1e6 => format!("{:.0} kB", v / 1e3),
        Unit::Bytes if v < 1e9 => format!("{:.0} MB", v / 1e6),
        Unit::Bytes => format!("{:.1} GB", v / 1e9),
        Unit::Count if v < 1e4 => format!("{v:.0}"),
        Unit::Count => format!("{v:.1e}"),
        Unit::Ratio => format!("{v:.1e}"),
    }
}

impl Level {
    /// Measures the cost of a placement: builds the fields at every resolution and runs
    /// every flight at preview and verification tolerance (as the game does).
    pub fn measure_cost(&self, player: &[Element]) -> Cost {
        let mut cost = Cost::new(self);
        let tol = self.tolerances();
        let metal = self.has_metal(player);
        if self.is_tube() {
            for refine in [1, 2] {
                let t = Instant::now();
                let mut sim = self.tube_sim(player, refine);
                sim.advance_to(self.tube.map_or(0.0, |s| s.goal.end));
                cost.add_tube_run(
                    t.elapsed().as_secs_f64(),
                    sim.segment_count(),
                    sim.steps(),
                    refine == 2,
                );
            }
            return cost;
        }
        if self.has_beams() {
            return self.measure_beam_cost(player, cost);
        }
        let preview = self.scenarios(player);
        if metal {
            cost.add_field(&self.field_at(player, Resolution::Display).0);
            cost.add_field(&preview[0].field);
        }
        for scn in &preview {
            let t = Instant::now();
            let traj = run(scn, &RunSettings::with_tolerance(tol.preview));
            cost.add_preview(t.elapsed().as_secs_f64(), &traj);
        }
        let fine = if metal {
            let fine = self.scenarios_at(player, Resolution::Verify);
            cost.add_verification_field(&fine[0].field);
            fine
        } else {
            preview
        };
        for scn in &fine {
            let t = Instant::now();
            let traj = run(scn, &RunSettings::with_tolerance(tol.verify));
            cost.add_verification(t.elapsed().as_secs_f64(), &traj);
        }
        cost
    }

    /// `measure_cost` for beam levels: one flight per disturbance, all particles
    /// together (the step counts are those of the whole system); the verification is the
    /// verdict's (`verification_beam_scenarios`: exact retarded for interacting beams at
    /// finite c).
    fn measure_beam_cost(&self, player: &[Element], mut cost: Cost) -> Cost {
        let tol = self.tolerances();
        let passes = [
            (
                true,
                self.beam_scenarios(player, Resolution::Preview),
                tol.preview,
            ),
            (false, self.verification_beam_scenarios(player), tol.verify),
        ];
        for (preview, scenarios, rtol) in passes {
            for scn in scenarios {
                let t = Instant::now();
                let r = physics::beam::run_beam(&scn, &RunSettings::with_tolerance(rtol));
                let seconds = t.elapsed().as_secs_f64();
                if let Some(traj) = r.trajectories.first() {
                    if preview {
                        cost.add_preview(seconds, traj);
                    } else {
                        cost.add_verification(seconds, traj);
                    }
                }
            }
        }
        cost
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meter(value: f64) -> Meter {
        Meter {
            label: "x",
            value,
            good: 1.0,
            limit: 10.0,
            unit: Unit::Seconds,
            hint: "",
        }
    }

    #[test]
    fn load_and_bar_position() {
        assert_eq!(meter(f64::NAN).load(), Load::Fine);
        assert_eq!(meter(1.0).load(), Load::Fine);
        assert_eq!(meter(5.0).load(), Load::High);
        assert_eq!(meter(10.1).load(), Load::Over);
        // Logarithmic from good/1000 to 3·limit.
        assert!(meter(1e-4).fraction().abs() < 1e-15);
        assert!((meter(30.0).fraction() - 1.0).abs() < 1e-12);
        let (g, l) = meter(0.0).marks();
        assert!((g - libm::log(1000.0) / libm::log(30_000.0)).abs() < 1e-12);
        assert!(g < l && l < 1.0);
        assert!(meter(0.5).fraction() < g && meter(2.0).fraction() > g);
    }

    #[test]
    fn formatting() {
        assert_eq!(format_value(f64::NAN, Unit::Seconds), "…");
        assert_eq!(format_value(6e-7, Unit::Seconds), "0.6 µs");
        assert_eq!(format_value(2.5e-4, Unit::Seconds), "250 µs");
        assert_eq!(format_value(0.394, Unit::Seconds), "394 ms");
        assert_eq!(format_value(4.1, Unit::Seconds), "4.1 s");
        assert_eq!(format_value(80e6, Unit::Bytes), "80 MB");
        assert_eq!(format_value(338.0, Unit::Count), "338");
    }

    /// The measurement covers every flight; metal levels report their linear systems and
    /// the boundary residual, which matches the level test's (< 1e-10).
    #[test]
    fn measured_cost() {
        let l = crate::shipped("first_bend");
        let c = l.measure_cost(&l.reference_solution);
        assert_eq!(c.previews, l.flight_count());
        assert_eq!(c.verifications, l.flight_count());
        assert!(c.max_steps > 0 && c.evaluations > c.max_steps);
        assert_eq!(c.setup, SetupCost::default());
        assert!(c.metal_residual.is_nan());
        assert_eq!(c.meters().len(), 3);

        let l = crate::shipped("high_voltage_dome");
        let c = l.measure_cost(&l.reference_solution);
        assert!(c.setup.unknowns > 0 && c.setup.bytes >= 8 * c.setup.unknowns);
        assert!(c.metal_residual < 1e-10, "{}", c.metal_residual);
        assert_eq!(c.meters().len(), 6);
    }
}
