//! Tube levels (docs/TUBES.md; physics in PHYSICS.md §2.11): the level is z-invariant.
//! Its electrodes are prisms along z (their `height` is ignored), the cathode emits
//! electrons under the space-charge limit, optionally an electrode emits ions (the
//! bipolar diode), projectiles (single line charges) may be fired through, and a
//! uniform magnetic field along z may thread it all. The goals: the current collected
//! by an electrode, averaged over a time window, and each projectile with a detector
//! arriving there. Shots, charges, magnets, coils and the other bodies break the
//! symmetry and are refused (`Level::tube_issues`).
//!
//! **Resolutions.** The preview is a coarse run; the verification halves every scale of
//! it (segments, macroparticle weight, step, softening). The verdict takes the fine
//! run's current with an error of `ERROR_FACTOR` times the two runs' difference, and a
//! projectile arrives when it does in both runs. Richardson extrapolation was tried and
//! rejected: it assumes first order, which V4 (the coaxial diode) showed but the
//! game-scale planar diode does not yet (the runs' differences shrink 7.9× from
//! refine 1→2 to 2→4), where it overshot by 4 %.

use physics::DVec3;
use physics::tube::{Emitter, PROJECTILE_TAG, Projectile, Tube, TubeState};
use physics::zinv::{Electrode as Prism, Electrodes as Prisms, Section};
use serde::{Deserialize, Serialize};

use crate::{Element, ElementKind, Level, Node};

/// A tube level's emitters, projectiles, field and goal.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
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
    /// An electrode that emits positive ions under the space-charge limit (the bipolar
    /// diode).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ions: Option<IonSource>,
    /// Line charges fired through the tube.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub projectiles: Vec<ProjectileSpec>,
    /// A uniform magnetic field along z (a solenoid around the tube).
    #[serde(default, skip_serializing_if = "crate::is_zero")]
    pub b_z: f64,
    /// The preview's time step (the verification's is half).
    #[serde(default = "default_step", skip_serializing_if = "is_default_step")]
    pub step: f64,
    /// The current goal, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub goal: Option<CurrentGoal>,
}

fn electron_ratio() -> f64 {
    -1.0
}

/// The preview's default time step: 1/30 (the first tube levels' 600 steps over 20).
pub const DEFAULT_STEP: f64 = 1.0 / 30.0;

fn default_step() -> f64 {
    DEFAULT_STEP
}

#[allow(clippy::trivially_copy_pass_by_ref, clippy::float_cmp)]
fn is_default_step(v: &f64) -> bool {
    *v == DEFAULT_STEP
}

/// An ion emitter.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct IonSource {
    pub electrode: usize,
    /// Charge-to-mass ratio of the ions (positive).
    pub charge_per_mass: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub emit_toward_deg: Option<f64>,
}

/// A projectile: one line charge (charge and mass per unit length) launched from a node.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProjectileSpec {
    pub node: Node,
    pub direction_deg: f64,
    pub speed: f64,
    #[serde(default, skip_serializing_if = "crate::is_zero")]
    pub launch: f64,
    pub charge: f64,
    pub mass: f64,
    /// It must arrive in the box spanned by these nodes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detector: Option<[Node; 2]>,
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
/// Preview resolution: macroparticles per the emitter's charge before emission. In
/// steady state the space charge in flight is of the order of that charge (the planar
/// diode's is 4/3 of it), so this sets the number of particles in flight.
pub const PREVIEW_PARTICLES: f64 = 150.0;
/// The fine run's error as a multiple of the two runs' difference. At first order the
/// error equals the difference; measured: V4 (order 0.95) 1.08×, the game-scale planar
/// diode at V = 40 0.15× (`tests/tubes.rs`, `planar_diode_convergence`, against
/// refine 4). 1.5 covers the first with a margin for orders down to about 0.6.
pub const ERROR_FACTOR: f64 = 1.5;

/// What a displayed particle is.
pub const KIND_NEGATIVE: u8 = 0;
pub const KIND_POSITIVE: u8 = 1;
pub const KIND_PROJECTILE: u8 = 2;

/// One displayed frame of a tube run: the particles (positions, charges per unit length
/// and kinds, f32: display only) and the electrodes' total surface density (applied and
/// induced), so that the maps can show the field at any frame.
#[derive(Clone, Debug, Default)]
pub struct TubeFrame {
    pub t: f64,
    pub x: Vec<[f32; 2]>,
    pub q: Vec<f32>,
    pub kind: Vec<u8>,
    pub sigma: Vec<f32>,
}

/// One run of a tube level (or a chunk of one being streamed: `TubeSim::take_chunk`).
#[derive(Clone, Debug, Default)]
pub struct TubeRun {
    /// A frame every step.
    pub frames: Vec<TubeFrame>,
    /// The goal electrode's collected carrier charge (positive for arrivals) over time.
    pub collected: Vec<(f64, f64)>,
    /// The goal's average current over its window, once the run has passed its end.
    pub current: Option<f64>,
    /// Per projectile: its track `(t, x)` while in flight, and when it reached its
    /// detector.
    pub tracks: Vec<Vec<(f64, [f32; 2])>>,
    pub arrivals: Vec<Option<f64>>,
    /// Particles in flight at the last frame.
    pub in_flight: usize,
    /// Time steps taken.
    pub steps: u32,
    /// The segments `(start, end)` the frames' densities belong to.
    pub segments: Vec<(DVec3, DVec3)>,
    /// Each prism's potential.
    pub potentials: Vec<f64>,
}

impl TubeRun {
    /// Appends a streamed chunk (its frames and samples follow this run's).
    pub fn extend(&mut self, chunk: TubeRun) {
        self.frames.extend(chunk.frames);
        self.collected.extend(chunk.collected);
        self.current = self.current.or(chunk.current);
        if self.tracks.len() < chunk.tracks.len() {
            self.tracks.resize(chunk.tracks.len(), Vec::new());
        }
        for (t, c) in self.tracks.iter_mut().zip(chunk.tracks) {
            t.extend(c);
        }
        if self.arrivals.len() < chunk.arrivals.len() {
            self.arrivals.resize(chunk.arrivals.len(), None);
        }
        for (a, c) in self.arrivals.iter_mut().zip(chunk.arrivals) {
            *a = a.or(c);
        }
        self.in_flight = chunk.in_flight;
        self.steps = chunk.steps;
        if self.segments.is_empty() {
            self.segments = chunk.segments;
            self.potentials = chunk.potentials;
        }
    }

    /// The frame at time `t` (the last one at or before it; the first before it).
    pub fn frame_at(&self, t: f64) -> Option<&TubeFrame> {
        let k = self.frames.partition_point(|f| f.t <= t);
        self.frames.get(k.saturating_sub(1))
    }

    /// The potential and field of a frame at `x`, as the maps show them: the segments'
    /// charge exactly, the particles as line charges softened over `PREVIEW_SEGMENT` (the
    /// preview's resolution of its space charge).
    pub fn field_at(&self, frame: &TubeFrame, x: DVec3) -> (f64, DVec3) {
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        for ((a, b), s) in self.segments.iter().zip(&frame.sigma) {
            let (f, g) = physics::zinv::segment_integrals(*a, *b, x);
            phi += f64::from(*s) * f;
            e += g * f64::from(*s);
        }
        let eps2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        for (p, q) in frame.x.iter().zip(&frame.q) {
            let r = DVec3::new(x.x - f64::from(p[0]), x.y - f64::from(p[1]), 0.0);
            let r2 = r.x * r.x + r.y * r.y + eps2;
            let lambda = f64::from(*q);
            phi -= lambda * libm::log(r2);
            e += r * (2.0 * lambda / r2);
        }
        (phi, e)
    }
}

/// A tube run in progress: advanced in pieces (`advance_to`), its frames taken as they
/// are made (`take_chunk`), so that a long display run can be streamed.
pub struct TubeSim {
    tube: Tube,
    state: TubeState,
    spec: TubeSpec,
    /// Each projectile's detector box (world coordinates).
    detectors: Vec<Option<(DVec3, DVec3)>>,
    at_start: Option<(f64, f64)>,
    /// What has not been taken yet.
    pending: TubeRun,
    steps: u32,
    current: Option<f64>,
    arrivals: Vec<Option<f64>>,
    /// Per projectile with a detector: its closest approach to the detector so far.
    closest: Vec<f64>,
}

/// Distance from `x` to the box `(lo, hi)` (0 inside).
fn box_distance(x: DVec3, (lo, hi): (DVec3, DVec3)) -> f64 {
    let dx = (lo.x - x.x).max(x.x - hi.x).max(0.0);
    let dy = (lo.y - x.y).max(x.y - hi.y).max(0.0);
    libm::sqrt(dx * dx + dy * dy)
}

impl TubeSim {
    pub fn t(&self) -> f64 {
        self.state.t
    }

    /// Steps taken and the segments of its linear system (for the cost meters).
    pub fn steps(&self) -> u32 {
        self.steps
    }

    pub fn segment_count(&self) -> usize {
        self.tube.electrodes.segments.len()
    }

    /// The goal's current, once the run has passed the goal window's end.
    pub fn current(&self) -> Option<f64> {
        self.current
    }

    /// When each projectile reached its detector (None: not, or no detector).
    pub fn arrivals(&self) -> &[Option<f64>] {
        &self.arrivals
    }

    /// Each projectile's closest approach to its detector so far (0: arrived; infinite
    /// without a detector).
    pub fn closest(&self) -> &[f64] {
        &self.closest
    }

    /// Advances to time `t_end` (whole steps), recording a frame every step.
    pub fn advance_to(&mut self, t_end: f64) {
        let sign = self.spec.charge_per_mass.signum();
        let n_proj = self.spec.projectiles.len();
        while self.state.t < t_end - 0.5 * self.tube.dt {
            self.tube.step(&mut self.state);
            self.steps += 1;
            let s = &self.state;
            let collected = self
                .spec
                .goal
                .map_or(0.0, |g| sign * s.collected[g.electrode][0]);
            if let Some(g) = self.spec.goal {
                if self.at_start.is_none() && s.t >= g.start {
                    self.at_start = Some((s.t, collected));
                }
                if self.current.is_none() && s.t >= g.end - 0.5 * self.tube.dt {
                    let (t0, q0) = self.at_start.unwrap_or((0.0, 0.0));
                    self.current = Some(if s.t > t0 {
                        (collected - q0) / (s.t - t0)
                    } else {
                        0.0
                    });
                    self.pending.current = self.current;
                }
            }
            if self.pending.tracks.len() < n_proj {
                self.pending.tracks.resize(n_proj, Vec::new());
                self.pending.arrivals.resize(n_proj, None);
            }
            for j in 0..n_proj {
                let Some(x) = s.projectile(j) else { continue };
                #[allow(clippy::cast_possible_truncation)]
                self.pending.tracks[j].push((s.t, [x.x as f32, x.y as f32]));
                if let Some(d) = self.detectors[j] {
                    let dist = box_distance(x, d);
                    self.closest[j] = self.closest[j].min(dist);
                    if dist == 0.0 && self.arrivals[j].is_none() {
                        self.arrivals[j] = Some(s.t);
                        self.pending.arrivals[j] = Some(s.t);
                    }
                }
            }
            #[allow(clippy::cast_possible_truncation)]
            let frame = TubeFrame {
                t: s.t,
                x: s.x.iter().map(|p| [p.x as f32, p.y as f32]).collect(),
                q: s.charge.iter().map(|q| *q as f32).collect(),
                kind: s
                    .tag
                    .iter()
                    .zip(&s.charge)
                    .map(|(t, q)| {
                        if *t >= PROJECTILE_TAG {
                            KIND_PROJECTILE
                        } else if *q < 0.0 {
                            KIND_NEGATIVE
                        } else {
                            KIND_POSITIVE
                        }
                    })
                    .collect(),
                sigma: self
                    .tube
                    .surface_density(s)
                    .iter()
                    .map(|v| *v as f32)
                    .collect(),
            };
            self.pending.frames.push(frame);
            self.pending.collected.push((s.t, collected));
        }
        self.pending.in_flight = self.state.x.len();
        self.pending.steps = self.steps;
    }

    /// The frames and samples made since the last call (the first chunk also carries the
    /// segments and potentials).
    pub fn take_chunk(&mut self) -> TubeRun {
        std::mem::take(&mut self.pending)
    }
}

/// A tube level's verdict on one goal, or on all of them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TubeStatus {
    /// Met within the error (the current in range; the projectile arrives in both runs).
    Met,
    /// Missed within the error.
    Missed,
    /// The error straddles a bound of the range, or the runs disagree on an arrival.
    Uncertain,
}

/// The current goal's verdict.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurrentVerdict {
    pub preview: f64,
    /// The verification's current: the verdict's value.
    pub fine: f64,
    pub error: f64,
    pub status: TubeStatus,
}

impl CurrentVerdict {
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

/// The verdict on a placement: each goal's and the whole's.
#[derive(Clone, Debug, PartialEq)]
pub struct TubeVerdict {
    pub current: Option<CurrentVerdict>,
    /// Per projectile with a detector: (projectile, preview arrival, verification
    /// arrival, status).
    pub arrivals: Vec<(usize, Option<f64>, Option<f64>, TubeStatus)>,
    pub status: TubeStatus,
}

impl TubeVerdict {
    /// From the two runs' results `(current, arrivals)`.
    pub fn new(
        spec: &TubeSpec,
        preview: &(Option<f64>, Vec<Option<f64>>),
        fine: &(Option<f64>, Vec<Option<f64>>),
    ) -> Self {
        let current = spec
            .goal
            .map(|g| CurrentVerdict::new(preview.0.unwrap_or(0.0), fine.0.unwrap_or(0.0), &g));
        let arrivals: Vec<_> = spec
            .projectiles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.detector.is_some())
            .map(|(j, _)| {
                let (a, b) = (
                    preview.1.get(j).copied().flatten(),
                    fine.1.get(j).copied().flatten(),
                );
                let status = match (a.is_some(), b.is_some()) {
                    (true, true) => TubeStatus::Met,
                    (false, false) => TubeStatus::Missed,
                    _ => TubeStatus::Uncertain,
                };
                (j, a, b, status)
            })
            .collect();
        let all = current
            .iter()
            .map(|c| c.status)
            .chain(arrivals.iter().map(|a| a.3));
        let mut status = TubeStatus::Met;
        for s in all {
            status = match (status, s) {
                (TubeStatus::Missed, _) | (_, TubeStatus::Missed) => TubeStatus::Missed,
                (TubeStatus::Uncertain, _) | (_, TubeStatus::Uncertain) => TubeStatus::Uncertain,
                _ => TubeStatus::Met,
            };
        }
        Self {
            current,
            arrivals,
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

    /// A tube run at `refine` = 1 (preview) or 2 (verification), at t = 0.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_sim(&self, player: &[Element], refine: u32) -> TubeSim {
        let spec = self.tube.clone().expect("a tube level");
        let r = f64::from(refine);
        let size = PREVIEW_SEGMENT / r;
        let prisms = Prisms::new(self.tube_prisms(player), &[], size);
        let toward = |a: Option<f64>| {
            a.map(|a| {
                let (s, c) = libm::sincos(a.to_radians());
                DVec3::new(c, s, 0.0)
            })
        };
        let weight = |e: usize| prisms.charge(e).abs().max(1e-300) / (PREVIEW_PARTICLES * r);
        let mut emitters = vec![Emitter {
            electrode: spec.cathode,
            charge_per_mass: spec.charge_per_mass,
            weight: weight(spec.cathode),
            emit_toward: toward(spec.emit_toward_deg),
        }];
        if let Some(ions) = spec.ions {
            emitters.push(Emitter {
                electrode: ions.electrode,
                charge_per_mass: ions.charge_per_mass,
                weight: weight(ions.electrode),
                emit_toward: toward(ions.emit_toward_deg),
            });
        }
        let projectiles = spec
            .projectiles
            .iter()
            .map(|p| {
                let (s, c) = libm::sincos(p.direction_deg.to_radians());
                Projectile {
                    x0: self.grid.position(p.node),
                    v0: DVec3::new(c, s, 0.0) * p.speed,
                    launch: p.launch,
                    charge: p.charge,
                    mass: p.mass,
                }
            })
            .collect();
        let detectors: Vec<_> = spec
            .projectiles
            .iter()
            .map(|p| {
                p.detector
                    .map(|[a, b]| (self.grid.position(a), self.grid.position(b)))
            })
            .collect();
        let n_proj = detectors.len();
        let max = self.grid.position(self.grid.max_node());
        let tube = Tube {
            electrodes: prisms,
            emitters,
            projectiles,
            softening: 0.5 * size,
            dt: spec.step / r,
            arena: Some((DVec3::ZERO, DVec3::new(max.x, max.y, 0.0))),
            b_z: spec.b_z,
        };
        let state = tube.state();
        let pending = TubeRun {
            segments: tube.electrodes.segments.clone(),
            potentials: tube.electrodes.potentials.clone(),
            ..TubeRun::default()
        };
        TubeSim {
            tube,
            state,
            spec,
            detectors,
            at_start: None,
            pending,
            steps: 0,
            current: None,
            arrivals: vec![None; n_proj],
            closest: vec![f64::INFINITY; n_proj],
        }
    }

    /// The end of a tube level's verification runs: the current goal's window end, and
    /// with projectile detectors the level's `t_max` (they must arrive by then).
    pub fn tube_verify_end(&self) -> f64 {
        self.tube.as_ref().map_or(0.0, |t| {
            let goal = t.goal.map_or(0.0, |g| g.end);
            if t.projectiles.iter().any(|p| p.detector.is_some()) {
                goal.max(self.physics.t_max)
            } else {
                goal
            }
        })
    }

    /// The end of a tube level's display run: its `t_max`, at least the verification's.
    pub fn tube_display_end(&self) -> f64 {
        self.tube_verify_end().max(self.physics.t_max)
    }

    /// Runs the tube at `refine` to `t_end` (whole), as one run.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_run(&self, player: &[Element], refine: u32, t_end: f64) -> TubeRun {
        let mut sim = self.tube_sim(player, refine);
        sim.advance_to(t_end);
        sim.take_chunk()
    }

    /// The results a verdict needs from a run at `refine`: the goal's current and the
    /// projectiles' arrivals, run to `tube_verify_end`.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_results(&self, player: &[Element], refine: u32) -> (Option<f64>, Vec<Option<f64>>) {
        let mut sim = self.tube_sim(player, refine);
        sim.advance_to(self.tube_verify_end());
        (sim.current(), sim.arrivals().to_vec())
    }

    /// The verdict of a placement: the preview and the verification run.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_verdict(&self, player: &[Element]) -> TubeVerdict {
        let spec = self.tube.as_ref().expect("a tube level");
        TubeVerdict::new(
            spec,
            &self.tube_results(player, 1),
            &self.tube_results(player, 2),
        )
    }

    /// The search objective (`solve::objective`) of a tube level, from the preview: the
    /// current's distance outside the goal's range in units of its width (0 inside), plus
    /// each projectile's closest approach to its detector (cells), and `Arrived` when all
    /// are 0, else `Rejected` (as a particle in its detector outside its acceptance).
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_objective(&self, player: &[Element]) -> (f64, physics::trajectory::Outcome) {
        use physics::trajectory::Outcome;
        let spec = self.tube.as_ref().expect("a tube level");
        let mut sim = self.tube_sim(player, 1);
        sim.advance_to(self.tube_verify_end());
        let mut score = 0.0;
        if let Some(g) = spec.goal {
            let i = sim.current().unwrap_or(0.0);
            score += (g.min - i).max(i - g.max).max(0.0) / (g.max - g.min).max(1e-300);
        }
        score += sim.closest().iter().filter(|d| d.is_finite()).sum::<f64>();
        if score == 0.0 {
            (0.0, Outcome::Arrived)
        } else {
            (score, Outcome::Rejected)
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
        if let Some(g) = spec.goal {
            if g.electrode >= n || g.electrode == spec.cathode {
                out.push("Tube: the goal must be an electrode other than the cathode.".into());
            }
            if !(g.min < g.max && 0.0 <= g.start && g.start < g.end && g.end.is_finite()) {
                out.push("Tube: the goal needs min < max and 0 ≤ start < end.".into());
            }
        } else if !spec.projectiles.iter().any(|p| p.detector.is_some()) {
            out.push("Tube: no goal (a current goal or a projectile with a detector).".into());
        }
        if !(spec.charge_per_mass.is_finite() && spec.charge_per_mass != 0.0) {
            out.push("Tube: the carriers' charge-to-mass ratio must be finite and nonzero.".into());
        }
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
                player
                    .iter()
                    .any(|e| !matches!(e.kind, ElementKind::Plate | ElementKind::Supply)),
                "placed elements other than plates and supplies",
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
            let (lo, hi) = boxes
                .iter()
                .map(potential)
                .fold((0.0_f64, 0.0_f64), |(a, b), v| (a.min(v), b.max(v)));
            let dv = hi - lo;
            let ions = spec.ions.map_or(0.0, |i| i.charge_per_mass);
            let ratios = spec
                .projectiles
                .iter()
                .map(|p| (p.charge / p.mass).abs())
                .chain([spec.charge_per_mass.abs(), ions]);
            let fall = ratios.fold(0.0_f64, f64::max);
            let launch = spec.projectiles.iter().map(|p| p.speed).fold(0.0, f64::max);
            let v = libm::sqrt(2.0 * fall * dv + launch * launch);
            if v > 0.1 * c {
                out.push(format!(
                    "Tube: the particles would reach {:.2} c; the tube model is Newtonian \
                     (up to 0.1 c): lower the voltages and speeds or leave c infinite.",
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
        // A projectile starting inside metal.
        let shapes = physics::bem::Electrodes::shapes_only(boxes);
        for (j, p) in spec.projectiles.iter().enumerate() {
            if shapes.contains(self.grid.position(p.node), crate::CONTACT_DISTANCE) {
                out.push(format!("Tube: projectile {} starts in metal.", j + 1));
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
