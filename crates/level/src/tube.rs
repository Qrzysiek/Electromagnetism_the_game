//! Tube levels (docs/TUBES.md; physics in PHYSICS.md §2.11): the level is z-invariant.
//! Its electrodes are prisms along z (their `height` is ignored), the cathode emits
//! electrons under the space-charge limit, and the goal is the current collected by an
//! electrode, averaged over a time window. Shots, charges, magnets, coils and the other
//! bodies break the symmetry and are refused (`Level::tube_issues`).
//!
//! **Resolutions.** The preview is a coarse run; the verification halves every scale of
//! it (segments, macroparticle weight, step, softening). The verdict takes the fine
//! run's current with an error of `ERROR_FACTOR` times the two runs' difference.
//! Richardson extrapolation was tried and rejected: it assumes first order, which V4
//! (the coaxial diode) showed but the game-scale planar diode does not yet (the runs'
//! differences shrink 7.9× from refine 1→2 to 2→4), where it overshot by 4 %.

use physics::DVec3;
use physics::tube::{Tube, TubeState};
use physics::zinv::{Electrode as Prism, Electrodes as Prisms, Section};
use serde::{Deserialize, Serialize};

use crate::{Element, ElementKind, Level};

/// A tube level's emitter and goal.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeSpec {
    /// The emitting electrode (an index into the level's `electrodes`).
    pub cathode: usize,
    /// Charge-to-mass ratio of the carriers (electrons: negative).
    #[serde(default = "electron_ratio")]
    pub charge_per_mass: f64,
    /// The cathode's coated face: only its faces whose outward normal is within 60° of
    /// this direction (degrees, in the plane) emit. None: every face.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub emit_toward_deg: Option<f64>,
    pub goal: CurrentGoal,
}

fn electron_ratio() -> f64 {
    -1.0
}

/// The current (charge per unit length per time, positive for carriers arriving)
/// collected by an electrode, averaged over `start..end`, must lie in `min..max`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurrentGoal {
    /// An index into the level's `electrodes`.
    pub electrode: usize,
    pub min: f64,
    pub max: f64,
    pub start: f64,
    pub end: f64,
}

/// Preview resolution: segment length (cells).
pub const PREVIEW_SEGMENT: f64 = 0.5;
/// Preview resolution: steps over the run (to the goal's `end`).
pub const PREVIEW_STEPS: u32 = 600;
/// Preview resolution: macroparticles per the cathode's charge before emission. In
/// steady state the space charge in flight is of the order of that charge (the planar
/// diode's is 4/3 of it), so this sets the number of particles in flight.
pub const PREVIEW_PARTICLES: f64 = 150.0;
/// Frames of the particles kept for display, over the run: every preview step, so that
/// slow playback stays smooth (150 looked jerky at the tube's playback rate; 600 frames
/// of ~200 particles are ~3 MB).
pub const FRAMES: u32 = PREVIEW_STEPS;
/// The fine run's error as a multiple of the two runs' difference. At first order the
/// error equals the difference; measured: V4 (order 0.95) 1.08×, the game-scale planar
/// diode at V = 40 0.15× (`tests/tubes.rs`, `planar_diode_convergence`, against
/// refine 4). 1.5 covers the first with a margin for orders down to about 0.6.
pub const ERROR_FACTOR: f64 = 1.5;

/// One run of a tube level.
#[derive(Clone, Debug, Default)]
pub struct TubeRun {
    /// Particle positions every few steps `(t, positions)`.
    pub frames: Vec<(f64, Vec<DVec3>)>,
    /// The goal electrode's collected carrier charge (positive for arrivals) over time.
    pub collected: Vec<(f64, f64)>,
    /// The goal's average current over its window.
    pub current: f64,
    /// Particles in flight at the end.
    pub in_flight: usize,
    /// Time steps taken.
    pub steps: u32,
    /// The final state, for the potential map: the segments `(start, end)` with their
    /// total surface density (the electrodes' and the induced), and every particle's
    /// line charge (the last frame's positions, each `particle_charge`).
    pub segments: Vec<(DVec3, DVec3)>,
    pub sigma: Vec<f64>,
    pub particle_charge: f64,
    /// Each prism's potential.
    pub potentials: Vec<f64>,
}

impl TubeRun {
    /// The final state's potential and field at `x`, as the maps show them: the
    /// segments' charge exactly, the particles as line charges softened over
    /// `PREVIEW_SEGMENT` (the preview's resolution of its space charge).
    pub fn field_at(&self, x: DVec3) -> (f64, DVec3) {
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        for ((a, b), s) in self.segments.iter().zip(&self.sigma) {
            let (f, g) = physics::zinv::segment_integrals(*a, *b, x);
            phi += s * f;
            e += g * *s;
        }
        let eps2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        let lambda = self.particle_charge;
        for p in self.frames.last().map_or(&[][..], |f| f.1.as_slice()) {
            let r = DVec3::new(x.x - p.x, x.y - p.y, 0.0);
            let r2 = r.x * r.x + r.y * r.y + eps2;
            phi -= lambda * libm::log(r2);
            e += r * (2.0 * lambda / r2);
        }
        (phi, e)
    }
}

/// A tube level's verdict.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TubeStatus {
    /// The current, within its error, lies in the goal's range.
    Met,
    /// It lies outside, within its error.
    Missed,
    /// The error straddles a bound of the range.
    Uncertain,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeVerdict {
    pub preview: f64,
    /// The verification's current: the verdict's value.
    pub fine: f64,
    pub error: f64,
    pub status: TubeStatus,
}

impl TubeVerdict {
    pub fn new(preview: f64, fine: f64, goal: &CurrentGoal) -> Self {
        let error = ERROR_FACTOR * (fine - preview).abs();
        let (lo, hi) = (fine - error, fine + error);
        let status = if lo >= goal.min && hi <= goal.max {
            TubeStatus::Met
        } else if hi < goal.min || lo > goal.max {
            TubeStatus::Missed
        } else {
            TubeStatus::Uncertain
        };
        Self {
            preview,
            fine,
            error,
            status,
        }
    }
}

impl Level {
    /// Whether this is a tube level.
    pub fn is_tube(&self) -> bool {
        self.tube.is_some()
    }

    /// The prisms of a placement: every box electrode (the level's at the player's
    /// supplies, then the player's plates), its cross-section.
    pub fn tube_prisms(&self, player: &[Element]) -> Vec<Prism> {
        self.all_box_electrodes(player)
            .into_iter()
            .map(|b| Prism {
                section: Section::Rect {
                    center: b.center,
                    angle: b.angle,
                    half_length: b.half_length,
                    half_thickness: b.half_thickness,
                },
                bias: b.bias,
            })
            .collect()
    }

    /// Runs the tube at `refine` = 1 (preview) or 2 (verification), calling `frame` with
    /// each displayed frame as it is made.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_run(
        &self,
        player: &[Element],
        refine: u32,
        mut frame: impl FnMut(&TubeRun),
    ) -> TubeRun {
        let spec = self.tube.expect("a tube level");
        let r = f64::from(refine);
        let size = PREVIEW_SEGMENT / r;
        let prisms = Prisms::new(self.tube_prisms(player), &[], size);
        let q0 = prisms.charge(spec.cathode).abs().max(1e-300);
        let max = self.grid.position(self.grid.max_node());
        let steps = PREVIEW_STEPS * refine;
        let tube = Tube {
            electrodes: prisms,
            cathode: spec.cathode,
            charge_per_mass: spec.charge_per_mass,
            weight: q0 / (PREVIEW_PARTICLES * r),
            softening: 0.5 * size,
            dt: spec.goal.end / f64::from(steps),
            arena: Some((DVec3::ZERO, DVec3::new(max.x, max.y, 0.0))),
            emit_toward: spec.emit_toward_deg.map(|a| {
                let (s, c) = libm::sincos(a.to_radians());
                DVec3::new(c, s, 0.0)
            }),
        };
        let sign = spec.charge_per_mass.signum();
        let mut s: TubeState = tube.state();
        let mut run = TubeRun::default();
        let every = (steps / FRAMES).max(1);
        let mut at_start = None;
        for k in 1..=steps {
            tube.step(&mut s);
            let collected = sign * s.collected[spec.goal.electrode];
            if at_start.is_none() && s.t >= spec.goal.start {
                at_start = Some((s.t, collected));
            }
            if k % every == 0 || k == steps {
                run.frames.push((s.t, s.x.clone()));
                run.collected.push((s.t, collected));
                frame(&run);
            }
        }
        let (t0, q0) = at_start.unwrap_or((0.0, 0.0));
        let (t1, q1) = *run.collected.last().unwrap_or(&(0.0, 0.0));
        run.current = if t1 > t0 { (q1 - q0) / (t1 - t0) } else { 0.0 };
        run.in_flight = s.x.len();
        run.steps = steps;
        run.sigma = tube.surface_density(&s);
        run.segments = tube.electrodes.segments.clone();
        run.particle_charge = sign * tube.weight;
        run.potentials = tube.electrodes.potentials.clone();
        run
    }

    /// The verdict of a placement: the preview and the verification run.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_verdict(&self, player: &[Element]) -> TubeVerdict {
        let goal = self.tube.expect("a tube level").goal;
        let coarse = self.tube_run(player, 1, |_| {}).current;
        let fine = self.tube_run(player, 2, |_| {}).current;
        TubeVerdict::new(coarse, fine, &goal)
    }

    /// The search objective (`solve::objective`) of a tube level, from the preview: the
    /// current's distance outside the goal's range in units of its width (0 inside), and
    /// `Arrived` inside, else `Rejected` (as a particle in its detector outside the
    /// detector's acceptance).
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_objective(&self, player: &[Element]) -> (f64, physics::trajectory::Outcome) {
        use physics::trajectory::Outcome;
        let g = self.tube.expect("a tube level").goal;
        let i = self.tube_run(player, 1, |_| {}).current;
        let outside = (g.min - i).max(i - g.max).max(0.0);
        if outside == 0.0 {
            (0.0, Outcome::Arrived)
        } else {
            (outside / (g.max - g.min).max(1e-300), Outcome::Rejected)
        }
    }

    /// What a tube level refuses: everything that breaks z-invariance, and an
    /// inconsistent spec.
    pub(crate) fn tube_issues(&self, player: &[Element], out: &mut Vec<String>) {
        let Some(spec) = self.tube else { return };
        let n = self.electrodes.len();
        if spec.cathode >= n {
            out.push("Tube: the cathode is not one of the level's electrodes.".into());
        }
        let g = spec.goal;
        if g.electrode >= n || g.electrode == spec.cathode {
            out.push("Tube: the goal must be an electrode other than the cathode.".into());
        }
        if !(g.min < g.max && 0.0 <= g.start && g.start < g.end && g.end.is_finite()) {
            out.push("Tube: the goal needs min < max and 0 ≤ start < end.".into());
        }
        if !(spec.charge_per_mass.is_finite() && spec.charge_per_mass != 0.0) {
            out.push("Tube: the carriers' charge-to-mass ratio must be finite and nonzero.".into());
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
                player
                    .iter()
                    .any(|e| !matches!(e.kind, ElementKind::Plate | ElementKind::Supply)),
                "placed elements other than plates and supplies",
            ),
        ];
        // Newtonian: the fastest carrier (falling through the largest potential
        // difference to the cathode) below 0.1 c, where the neglected v²/c² corrections
        // (1 %) are below the macroparticles' error (5–8 % between the two runs).
        let potential = |b: &physics::bem::BoxElectrode| match b.bias {
            physics::conductor::Bias::Potential(v) => v,
            _ => 0.0,
        };
        let boxes = self.all_box_electrodes(player);
        if let (Some(c), Some(cathode)) = (self.physics.c, boxes.get(spec.cathode)) {
            let vc = potential(cathode);
            let dv = boxes
                .iter()
                .map(|b| (potential(b) - vc).abs())
                .fold(0.0, f64::max);
            let v = libm::sqrt(2.0 * spec.charge_per_mass.abs() * dv);
            if v > 0.1 * c {
                out.push(format!(
                    "Tube: the electrons would reach {:.2} c; the tube model is Newtonian \
                     (up to 0.1 c): lower the voltages or leave c infinite.",
                    v / c
                ));
            }
        }
        // Prisms must not touch: the boundary elements of touching electrodes at different
        // potentials would hold an infinite field between them.
        for (i, a) in boxes.iter().enumerate() {
            for (j, b) in boxes.iter().enumerate().skip(i + 1) {
                if crate::rect_distance(&crate::rect_corners(a), &crate::rect_corners(b))
                    < crate::CONTACT_DISTANCE
                {
                    out.push(format!(
                        "Electrodes {} and {} touch or overlap.",
                        i + 1,
                        j + 1
                    ));
                }
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
