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
    /// The cathode's temperature-limited emission: at most this current per unit length
    /// of its surface (`physics::tube::Emitter::saturation`). None: space-charge limited
    /// only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub saturation: Option<f64>,
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
    /// A solenoid the player powers: its terminal node, where a power supply (from
    /// `limits.solenoid_fields`) sets B along z in place of `b_z`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solenoid: Option<Node>,
    /// A temperature-limited cathode's macroparticles per unit time at its saturation, in
    /// the preview (the verification's twice as many; `TL_RATE` by default): a beam
    /// modulated at a high frequency needs many per period.
    #[serde(default = "default_rate", skip_serializing_if = "is_default_rate")]
    pub particle_rate: f64,
    /// The preview's time step (the verification's is half).
    #[serde(default = "default_step", skip_serializing_if = "is_default_step")]
    pub step: f64,
    /// The current goals (older level files: a single `goal`).
    #[serde(
        default,
        alias = "goal",
        deserialize_with = "one_or_many",
        skip_serializing_if = "Vec::is_empty"
    )]
    pub goals: Vec<CurrentGoal>,
    /// External circuits of electrodes (per unit length; PHYSICS.md 2.11).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub circuits: Vec<TubeCircuit>,
    /// Goals on electrodes' potentials over windows (their mean or their swing).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub potential_goals: Vec<PotentialGoal>,
}

/// An electrode's external circuit (per unit length along z): a source through a
/// resistance into the electrode, a capacitance from it to ground (`physics::tube::
/// Circuit`). Resistance 0: an ideal source the electrode follows; without a source
/// (`resistance` infinite, written as `null`): floating with its capacitance. A tunable
/// electrode's power supply sets the source's level (`Source::with_level`), or with
/// `supply_offset` a sine's offset (a grid's bias under its signal).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeCircuit {
    pub electrode: usize,
    pub source: crate::Source,
    /// None: no source (floating).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resistance: Option<f64>,
    #[serde(default, skip_serializing_if = "crate::is_zero")]
    pub capacitance: f64,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub supply_offset: bool,
    /// Further electrodes on the same ideal source (resistance 0 only): a grid's wires.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub linked: Vec<usize>,
    /// An inductance in series with the resistance (`physics::tube::Circuit`): a tank's
    /// coil. 0: none.
    #[serde(default, skip_serializing_if = "crate::is_zero")]
    pub inductance: f64,
    /// The branch's other end: this electrode instead of ground (a tank between two
    /// electrodes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub return_to: Option<usize>,
}

/// What a potential goal measures over its window.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Measure {
    /// The mean potential (a rectifier's output).
    Mean,
    /// Half its peak-to-peak swing.
    Swing,
    /// The amplitude of its oscillation at the goal's `omega` (an amplifier's output):
    /// a sine at ω fitted together with a linear trend over the window (a node still
    /// settling to its operating point drifts), exact for a sine on a line.
    Amplitude,
    /// The phase of its oscillation at the goal's `omega` relative to `sin(ω t)` (the
    /// phase of every source's sine without its own phase: a klystron's drive), in
    /// degrees: `V ≈ A sin(ω t + φ)`, from the same Fourier sums as the amplitude. Phases
    /// are compared modulo 360: a value is taken at its turn nearest the goal's middle.
    Phase,
}

impl PotentialGoal {
    /// A measured value as the goal compares it: a phase at its turn nearest the goal's
    /// middle (`Measure::Phase`), anything else as it is.
    pub fn compared(&self, v: f64) -> f64 {
        if self.measure == Measure::Phase {
            let middle = 0.5 * (self.min + self.max);
            v + 360.0 * libm::round((middle - v) / 360.0)
        } else {
            v
        }
    }
}

/// An electrode's potential, measured over `start..end`, must lie in `min..max`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PotentialGoal {
    pub electrode: usize,
    pub measure: Measure,
    /// The angular frequency of `Measure::Amplitude` and `Measure::Phase`.
    #[serde(default, skip_serializing_if = "crate::is_zero")]
    pub omega: f64,
    pub min: f64,
    pub max: f64,
    pub start: f64,
    pub end: f64,
}

fn electron_ratio() -> f64 {
    -1.0
}

/// The preview's default time step: 1/30 (the first tube levels' 600 steps over 20).
pub const DEFAULT_STEP: f64 = 1.0 / 30.0;

fn default_step() -> f64 {
    DEFAULT_STEP
}

fn default_rate() -> f64 {
    TL_RATE
}

#[allow(clippy::trivially_copy_pass_by_ref, clippy::float_cmp)]
fn is_default_rate(v: &f64) -> bool {
    *v == TL_RATE
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

/// The current (charge per unit length per time, positive for carriers arriving) of one
/// kind of carrier collected by an electrode, averaged over `start..end`, must lie in
/// `min..max`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurrentGoal {
    /// An index into the level's `electrodes`.
    pub electrode: usize,
    pub min: f64,
    pub max: f64,
    pub start: f64,
    pub end: f64,
    /// Whose current: the cathode's carriers (the default) or the ion source's ions.
    #[serde(default, skip_serializing_if = "Carriers::is_cathode")]
    pub carriers: Carriers,
}

/// The carriers a current goal counts (projectiles never count).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Carriers {
    /// The cathode's (electrons in every level so far).
    #[default]
    Cathode,
    /// The ion source's (`TubeSpec::ions`): e.g. the ion current into the cathode, its
    /// wear in a bipolar diode.
    Ions,
}

impl Carriers {
    #[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
    fn is_cathode(&self) -> bool {
        *self == Carriers::Cathode
    }

    /// The emitter index of these carriers in a run (`TubeState::collected`'s second
    /// index): the cathode is emitter 0, the ion source 1.
    fn emitter(self) -> usize {
        match self {
            Carriers::Cathode => 0,
            Carriers::Ions => 1,
        }
    }
}

/// The current goals: a list, or in older level files a single goal (`goal`).
fn one_or_many<'de, D>(d: D) -> Result<Vec<CurrentGoal>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum OneOrMany {
        One(CurrentGoal),
        Many(Vec<CurrentGoal>),
    }
    Ok(match Option::<OneOrMany>::deserialize(d)? {
        None => Vec::new(),
        Some(OneOrMany::One(g)) => vec![g],
        Some(OneOrMany::Many(v)) => v,
    })
}

/// A tube run's resolution: `refine` divides the segments, the step, the macroparticles'
/// weight and the softening (1: the preview, 2: the verification); `particles` and
/// `softening` further multiply the macroparticles' number and the softening alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeResolution {
    pub refine: u32,
    pub particles: f64,
    pub softening: f64,
}

impl TubeResolution {
    /// The game's resolutions: everything scaled by `refine` together.
    pub fn refined(refine: u32) -> Self {
        Self {
            refine,
            particles: 1.0,
            softening: 1.0,
        }
    }
}

/// Preview resolution: segment length (cells).
pub const PREVIEW_SEGMENT: f64 = 0.5;
/// Preview resolution: macroparticles per the emitter's charge before emission: the
/// largest charge of its carriers' sign its emitting face holds without particles, over a
/// period of the sources (an AC anode leaves it empty at t = 0). In steady state the space
/// charge in flight is of the order of that charge (the planar diode's is 4/3 of it), so
/// this sets the number of particles in flight.
pub const PREVIEW_PARTICLES: f64 = 150.0;
/// The macroparticle weight's floor, as a fraction of the emitter's charge with every
/// other electrode at the level's largest voltage (`Level::tube_reference_voltage`): an
/// emitter whose face holds (almost) no charge of its carriers' sign at any phase of the
/// sources still emits what passing charges induce on it, and a weight from its own charge
/// would be (almost) zero (a design scan once asked for ~1e300 macroparticles). 1/1000: at
/// the former 1/20 it overrode a grid-screened cathode's own charge, whose face holds far
/// less than the level's largest voltage implies (a triode with B+ at 100 on a ±120
/// slider: 27 macroparticles in flight instead of 152, and its current 28 % low).
pub const WEIGHT_FLOOR: f64 = 1.0 / 1000.0;

/// The tube's envelope (PHYSICS.md §2.11): a grounded shield around the board, this far
/// outside its edge (so that a projectile launched from the edge starts inside it). In the
/// z-invariant world an open set of electrodes carries a net charge that depends on the
/// logarithm's reference radius, and so did its far field: shifting every potential by
/// the same amount changed a triode's electrons lost off the board from 0.20 to 0.54 per
/// unit time. Inside a closed grounded contour only potential differences count; what
/// leaves the board ends on the envelope (counted as lost). A real tube's envelope
/// charges to about the cathode's potential, which is ground in every level.
pub const ENVELOPE_MARGIN: f64 = 0.5;

/// The envelope's segment length, at both resolutions (the preview's): it lies away from
/// the flow, and its own error at 0.6–1 cell from it is then 3.6e-4 of the potential
/// (V13), far below the macroparticles'. Segments of 1 gave 1.2e-2 (the graded sides'
/// middle segments are π/2 times longer): rejected.
pub const ENVELOPE_SEGMENT: f64 = 0.5;

/// The emitters' segments as a fraction of the other electrodes' (PHYSICS.md §2.11, V15):
/// macroparticles are released half of the emitter's segment out
/// (`physics::tube::EMISSION_OFFSET`), and the space-charge layer they skip is a
/// first-order error in that distance (V4). A quarter cut V4's error at the verification's
/// resolution from 7.8 % to 1.0 %, and the game-scale diode's verdict error from ±0.124 to
/// ±0.049, for 30–50 % more run time; narrow gaps (a triode's grid 1.25 from its cathode)
/// need it.
pub const EMITTER_MESH: f64 = 0.25;

/// The macroparticles' radius as a fraction of a segment (PHYSICS.md §2.11): each is a
/// uniformly charged rod, its field softened inside it and exactly Coulomb's outside. It
/// replaced the softening `r⃗/(r² + ε²)` at an eighth, which never quite became Coulomb's
/// (ε²/r²: 0.4 % one cell away), so that far groups of particles could not be summed as
/// multipoles exactly. Rods of a quarter (the same largest force) converged more slowly on
/// the calibration diode (ratio 0.76, against 0.51 for the eighth). The softening smooths the space
/// charge's field over that length, and in a narrow gap the smoothed space charge limits the
/// current too weakly: at half a segment (Plummer, the earliest value) the triode's output,
/// 1.25 between its cathode and grid, came out 3.89, 3.38, 3.17 at refine 1, 2, 3 (its
/// limit ~3.0); at an eighth 3.34, 3.12, 3.04.
pub const SOFTENING: f64 = 0.125;

/// A temperature-limited emitter's macroparticles per unit time at its saturation
/// (preview): its weight is at most `saturation × perimeter / TL_RATE`. The weight rule
/// of the space-charge limit (its charge before emission) left the Electron gun's beam
/// at 36 macroparticles in flight, sparse to the eye and grainy for the verdict; with
/// this rule ~150.
pub const TL_RATE: f64 = 50.0;
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
    /// Every electrode's potential.
    pub potentials: Vec<f32>,
}

/// One run of a tube level (or a chunk of one being streamed: `TubeSim::take_chunk`).
#[derive(Clone, Debug, Default)]
pub struct TubeRun {
    /// A frame every step.
    pub frames: Vec<TubeFrame>,
    /// Per current goal: its electrode's collected charge of its carriers (positive for
    /// arrivals) over time.
    pub collected: Vec<Vec<(f64, f64)>>,
    /// Per current goal: its average current over its window, once the run has passed
    /// its end.
    pub currents: Vec<Option<f64>>,
    /// Each potential goal's measure, once the run has passed its window.
    pub potential: Vec<Option<f64>>,
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
    /// The lowest and highest potential the electrodes reach: their own, and their circuits'
    /// sources over a period (the potential map's scale; at t = 0 alone an AC level's
    /// electrodes are all at 0).
    pub potential_range: (f64, f64),
}

impl TubeRun {
    /// Appends a streamed chunk (its frames and samples follow this run's).
    pub fn extend(&mut self, chunk: TubeRun) {
        self.frames.extend(chunk.frames);
        if self.collected.len() < chunk.collected.len() {
            self.collected.resize(chunk.collected.len(), Vec::new());
        }
        for (a, c) in self.collected.iter_mut().zip(chunk.collected) {
            a.extend(c);
        }
        if self.currents.len() < chunk.currents.len() {
            self.currents.resize(chunk.currents.len(), None);
        }
        for (a, c) in self.currents.iter_mut().zip(chunk.currents) {
            *a = a.or(c);
        }
        if self.potential.len() < chunk.potential.len() {
            self.potential.resize(chunk.potential.len(), None);
        }
        for (a, c) in self.potential.iter_mut().zip(chunk.potential) {
            *a = a.or(c);
        }
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
            self.potential_range = chunk.potential_range;
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
    /// Per current goal: the time and collected charge at its window's start.
    at_start: Vec<Option<(f64, f64)>>,
    /// What has not been taken yet.
    pending: TubeRun,
    steps: u32,
    currents: Vec<Option<f64>>,
    /// Each potential goal's window so far: sum, count, least, most, and the samples
    /// `(t, V)` (for `Measure::Amplitude`).
    windows: Vec<(f64, f64, f64, f64)>,
    series: Vec<Vec<(f64, f64)>>,
    potential: Vec<Option<f64>>,
    arrivals: Vec<Option<f64>>,
    /// Per projectile with a detector: its closest approach to the detector so far.
    closest: Vec<f64>,
}

/// The amplitude and phase (degrees, relative to `sin(ω t)`) of the oscillation at angular
/// frequency `omega` in the samples `(t, v)`: the least-squares fit of
/// `a + b (t − t̄) + c cos ωt + d sin ωt`, all four together (exact for a sine on a line,
/// over any window; fitting the line first and then the sine took part of the sine with
/// the line: 1e-3 of the amplitude over 10 periods, up to tens of per cent over one), and
/// `A = √(c² + d²)`, `φ = atan2(c, d)` for `v ≈ A sin(ω t + φ)`.
#[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
fn fourier_at(samples: &[(f64, f64)], omega: f64) -> (f64, f64) {
    #[allow(clippy::cast_precision_loss)]
    let n = samples.len() as f64;
    if samples.len() < 4 {
        return (0.0, 0.0);
    }
    let mean_t = samples.iter().map(|(t, _)| t).sum::<f64>() / n;
    // Normal equations of the four basis functions.
    let mut m = [[0.0_f64; 5]; 4];
    for (t, v) in samples {
        let (s, c) = libm::sincos(omega * t);
        let f = [1.0, t - mean_t, c, s];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += f[i] * f[j];
            }
            m[i][4] += f[i] * v;
        }
    }
    // Gaussian elimination with partial pivoting (a 4 × 4 system, fixed order).
    for k in 0..4 {
        let p = (k..4)
            .max_by(|&i, &j| m[i][k].abs().total_cmp(&m[j][k].abs()))
            .unwrap_or(k);
        m.swap(k, p);
        if m[k][k] == 0.0 {
            return (0.0, 0.0);
        }
        for i in k + 1..4 {
            let f = m[i][k] / m[k][k];
            for j in k..5 {
                m[i][j] -= f * m[k][j];
            }
        }
    }
    let mut x = [0.0_f64; 4];
    for i in (0..4).rev() {
        let rest: f64 = (i + 1..4).map(|j| m[i][j] * x[j]).sum();
        x[i] = (m[i][4] - rest) / m[i][i];
    }
    let (c, d) = (x[2], x[3]);
    (libm::sqrt(c * c + d * d), libm::atan2(c, d).to_degrees())
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

    /// Macroparticles in flight now.
    pub fn in_flight(&self) -> usize {
        self.state.x.len()
    }

    /// Each emitter's macroparticle weight.
    pub fn weights(&self) -> Vec<f64> {
        self.tube.emitters.iter().map(|e| e.weight).collect()
    }

    /// Each current goal's current, once the run has passed its window's end.
    pub fn currents(&self) -> &[Option<f64>] {
        &self.currents
    }

    /// The charge per unit length each electrode has collected so far, per emitter
    /// (signed: electrons negative), and what left the board.
    pub fn collected(&self) -> &[Vec<f64>] {
        &self.state.collected
    }

    pub fn lost(&self) -> f64 {
        self.state.lost
    }

    /// Emission paused at `physics::tube::MAX_PARTICLES` at least once.
    pub fn capped(&self) -> bool {
        self.state.capped
    }

    /// Where each projectile ended (None: still flying, or not launched yet).
    pub fn fates(&self) -> &[Option<physics::tube::Fate>] {
        &self.state.fates
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
        let n_proj = self.spec.projectiles.len();
        let n_goals = self.spec.goals.len();
        while self.state.t < t_end - 0.5 * self.tube.dt {
            self.tube.step(&mut self.state);
            self.steps += 1;
            let s = &self.state;
            if self.pending.collected.len() < n_goals {
                self.pending.collected.resize(n_goals, Vec::new());
                self.pending.currents.resize(n_goals, None);
            }
            for (k, g) in self.spec.goals.iter().enumerate() {
                // The carriers' charge collected, positive for arrivals.
                let sign = match g.carriers {
                    Carriers::Cathode => self.spec.charge_per_mass.signum(),
                    Carriers::Ions => 1.0,
                };
                let collected = sign
                    * s.collected[g.electrode]
                        .get(g.carriers.emitter())
                        .copied()
                        .unwrap_or(0.0);
                if self.at_start[k].is_none() && s.t >= g.start {
                    self.at_start[k] = Some((s.t, collected));
                }
                if self.currents[k].is_none() && s.t >= g.end - 0.5 * self.tube.dt {
                    let (t0, q0) = self.at_start[k].unwrap_or((0.0, 0.0));
                    self.currents[k] = Some(if s.t > t0 {
                        (collected - q0) / (s.t - t0)
                    } else {
                        0.0
                    });
                    self.pending.currents[k] = self.currents[k];
                }
                self.pending.collected[k].push((s.t, collected));
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
                potentials: s.potentials.iter().map(|v| *v as f32).collect(),
            };
            self.pending.frames.push(frame);
            // Each potential goal's measure over its window.
            if self.pending.potential.len() < self.potential.len() {
                self.pending.potential.resize(self.potential.len(), None);
            }
            for (k, g) in self.spec.potential_goals.iter().enumerate() {
                if s.t < g.start - 0.5 * self.tube.dt || self.potential[k].is_some() {
                    continue;
                }
                let v = s.potentials[g.electrode];
                let w = &mut self.windows[k];
                w.0 += v;
                w.1 += 1.0;
                w.2 = w.2.min(v);
                w.3 = w.3.max(v);
                if matches!(g.measure, Measure::Amplitude | Measure::Phase) {
                    self.series[k].push((s.t, v));
                }
                if s.t >= g.end - 0.5 * self.tube.dt {
                    let m = match g.measure {
                        Measure::Mean => w.0 / w.1.max(1.0),
                        Measure::Swing => 0.5 * (w.3 - w.2),
                        Measure::Amplitude => fourier_at(&self.series[k], g.omega).0,
                        Measure::Phase => fourier_at(&self.series[k], g.omega).1,
                    };
                    self.potential[k] = Some(m);
                    self.pending.potential[k] = Some(m);
                }
            }
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
    /// The verdict on a goal's range `min..max` from the two runs' values.
    pub fn new(preview: f64, fine: f64, min: f64, max: f64) -> Self {
        let error = ERROR_FACTOR * (fine - preview).abs();
        let (lo, hi) = (fine - error, fine + error);
        let status = if lo >= min && hi <= max {
            TubeStatus::Met
        } else if hi < min || lo > max {
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

/// What a verdict needs from one run.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TubeResults {
    /// Each current goal's current.
    pub currents: Vec<Option<f64>>,
    /// Each potential goal's measure.
    pub potential: Vec<Option<f64>>,
    /// When each projectile reached its detector.
    pub arrivals: Vec<Option<f64>>,
    /// Emission paused at the particle cap.
    pub capped: bool,
}

impl TubeSim {
    /// The run's results so far.
    pub fn results(&self) -> TubeResults {
        TubeResults {
            currents: self.currents.clone(),
            potential: self.potential.clone(),
            arrivals: self.arrivals.clone(),
            capped: self.state.capped,
        }
    }
}

/// The verdict on a placement: each goal's and the whole's.
#[derive(Clone, Debug, PartialEq)]
pub struct TubeVerdict {
    /// Each current goal's verdict.
    pub currents: Vec<CurrentVerdict>,
    /// Each potential goal's verdict (the same rule: the fine run's value, 1.5 times the
    /// runs' difference as its error).
    pub potential: Vec<CurrentVerdict>,
    /// Per projectile with a detector: (projectile, preview arrival, verification
    /// arrival, status).
    pub arrivals: Vec<(usize, Option<f64>, Option<f64>, TubeStatus)>,
    /// A run hit the particle cap (`physics::tube::MAX_PARTICLES`): uncertain.
    pub capped: bool,
    pub status: TubeStatus,
}

impl TubeVerdict {
    /// From the two runs' results.
    pub fn new(spec: &TubeSpec, preview: &TubeResults, fine: &TubeResults) -> Self {
        let currents: Vec<CurrentVerdict> = spec
            .goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    preview.currents.get(k).copied().flatten().unwrap_or(0.0),
                    fine.currents.get(k).copied().flatten().unwrap_or(0.0),
                    g.min,
                    g.max,
                )
            })
            .collect();
        let potential: Vec<CurrentVerdict> = spec
            .potential_goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    g.compared(preview.potential.get(k).copied().flatten().unwrap_or(0.0)),
                    g.compared(fine.potential.get(k).copied().flatten().unwrap_or(0.0)),
                    g.min,
                    g.max,
                )
            })
            .collect();
        let arrivals: Vec<_> = spec
            .projectiles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.detector.is_some())
            .map(|(j, _)| {
                let (a, b) = (
                    preview.arrivals.get(j).copied().flatten(),
                    fine.arrivals.get(j).copied().flatten(),
                );
                let status = match (a.is_some(), b.is_some()) {
                    (true, true) => TubeStatus::Met,
                    (false, false) => TubeStatus::Missed,
                    _ => TubeStatus::Uncertain,
                };
                (j, a, b, status)
            })
            .collect();
        let all = currents
            .iter()
            .chain(potential.iter())
            .map(|c| c.status)
            .chain(arrivals.iter().map(|a| a.3));
        let capped = preview.capped || fine.capped;
        let mut status = if capped {
            TubeStatus::Uncertain
        } else {
            TubeStatus::Met
        };
        for s in all {
            status = match (status, s) {
                (TubeStatus::Missed, _) | (_, TubeStatus::Missed) => TubeStatus::Missed,
                (TubeStatus::Uncertain, _) | (_, TubeStatus::Uncertain) => TubeStatus::Uncertain,
                _ => TubeStatus::Met,
            };
        }
        Self {
            currents,
            potential,
            arrivals,
            capped,
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
        let mut prisms: Vec<Prism> = self
            .all_box_electrodes(player)
            .into_iter()
            .enumerate()
            .map(|(i, b)| {
                let shape = self
                    .electrodes
                    .get(i)
                    .map_or(crate::ElectrodeShape::Box, |e| e.shape);
                let e = self.electrodes.get(i);
                let section = match (shape, e) {
                    (crate::ElectrodeShape::Rod, Some(e)) => Section::Circle {
                        center: b.center,
                        radius: e.length / 2.0,
                    },
                    (crate::ElectrodeShape::Ring, Some(e)) => Section::Ring {
                        center: b.center,
                        radius: e.length / 2.0,
                        thickness: e.thickness,
                    },
                    _ => Section::Rect {
                        center: b.center,
                        angle: b.angle,
                        half_length: b.half_length,
                        half_thickness: b.half_thickness,
                    },
                };
                Prism {
                    section,
                    bias: b.bias,
                }
            })
            .collect();
        // An electrode on a source (ideal or through a resistance) starts at the source's
        // value at t = 0: the supply was on before the cathode heated, so its capacitors
        // are charged. A floating one keeps its own bias.
        for c in self.tube_circuits(player) {
            if c.resistance.is_finite()
                && let Some(p) = prisms.get_mut(c.electrode)
            {
                p.bias = physics::conductor::Bias::Potential(c.source.value(0.0));
            }
        }
        // The envelope, grounded, last (ENVELOPE_MARGIN).
        let max = self.grid.position(self.grid.max_node());
        prisms.push(Prism {
            section: Section::Frame {
                center: DVec3::new(0.5 * max.x, 0.5 * max.y, 0.0),
                half_width: 0.5 * max.x + ENVELOPE_MARGIN,
                half_height: 0.5 * max.y + ENVELOPE_MARGIN,
                segment: ENVELOPE_SEGMENT,
            },
            bias: physics::conductor::Bias::Grounded,
        });
        prisms
    }

    /// A tube run at `refine` = 1 (preview) or 2 (verification), at t = 0.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_sim(&self, player: &[Element], refine: u32) -> TubeSim {
        self.tube_sim_at(player, TubeResolution::refined(refine))
    }

    /// A tube run at a resolution whose parts may be varied separately (convergence
    /// studies: `generator tube-scan --particles --softening`), at t = 0.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_sim_at(&self, player: &[Element], res: TubeResolution) -> TubeSim {
        let spec = self.tube.clone().expect("a tube level");
        let refine = res.refine;
        let r = f64::from(refine);
        let size = PREVIEW_SEGMENT / r;
        // The emitters meshed finer (EMITTER_MESH).
        let shapes = self.tube_prisms(player);
        let sizes: Vec<f64> = (0..shapes.len())
            .map(|e| {
                let emits = e == spec.cathode || spec.ions.is_some_and(|i| i.electrode == e);
                if emits { size * EMITTER_MESH } else { size }
            })
            .collect();
        let prisms = Prisms::with_sizes(shapes, &[], &sizes);
        let toward = |a: Option<f64>| {
            a.map(|a| {
                let (s, c) = libm::sincos(a.to_radians());
                DVec3::new(c, s, 0.0)
            })
        };
        // The macroparticles' weights are set once the tube is built (below); their floor.
        let v_ref = self.tube_reference_voltage();
        let floor = |e: usize| {
            let others: f64 = prisms.capacitance[e]
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != e)
                .map(|(_, c)| c.abs())
                .sum();
            WEIGHT_FLOOR * v_ref * others
        };
        // A temperature-limited cathode emits at most its saturation over its surface:
        // the weight also gives `TL_RATE` macroparticles per unit time at that rate (the
        // whole perimeter: an upper bound of the emitting faces, so at least that many).
        let perimeter: f64 = prisms
            .segments
            .iter()
            .enumerate()
            .filter(|(k, _)| prisms.owner_of(*k) == spec.cathode)
            .map(|(_, (a, b))| (*b - *a).length())
            .sum();
        let mut floors = vec![floor(spec.cathode)];
        let mut emitters = vec![Emitter {
            electrode: spec.cathode,
            charge_per_mass: spec.charge_per_mass,
            weight: 1.0,
            emit_toward: toward(spec.emit_toward_deg),
            saturation: spec.saturation,
        }];
        if let Some(ions) = spec.ions {
            floors.push(floor(ions.electrode));
            emitters.push(Emitter {
                electrode: ions.electrode,
                charge_per_mass: ions.charge_per_mass,
                weight: 1.0,
                emit_toward: toward(ions.emit_toward_deg),
                saturation: None,
            });
        }
        // The level's projectiles, then the player's free charges (thrown at t = 0).
        let mut projectiles: Vec<Projectile> = spec
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
        let mut detectors: Vec<_> = spec
            .projectiles
            .iter()
            .map(|p| {
                p.detector
                    .map(|[a, b]| (self.grid.position(a), self.grid.position(b)))
            })
            .collect();
        for e in player.iter().filter(|e| e.kind == ElementKind::Free) {
            let (s, c) = libm::sincos(e.angle_deg.to_radians());
            projectiles.push(Projectile {
                x0: self.grid.position(e.node),
                v0: DVec3::new(c, s, 0.0) * e.speed.unwrap_or(0.0),
                launch: 0.0,
                charge: e.value,
                mass: self.limits.free_mass,
            });
            detectors.push(None);
        }
        let n_proj = detectors.len();
        let max = self.grid.position(self.grid.max_node());
        // The solenoid's field: the player's supply on its terminal, else the level's.
        let b_z = spec
            .solenoid
            .and_then(|n| {
                player
                    .iter()
                    .find(|e| e.kind == ElementKind::Supply && e.node == n)
            })
            .map_or(spec.b_z, |e| e.value);
        let mut tube = Tube {
            electrodes: prisms,
            emitters,
            projectiles,
            softening: SOFTENING * size * res.softening,
            dt: spec.step / r,
            arena: Some((DVec3::ZERO, DVec3::new(max.x, max.y, 0.0))),
            b_z,
            circuits: self.tube_circuits(player),
        };
        // Each emitter's weight: 1/PREVIEW_PARTICLES of the largest charge its emitting face
        // holds (no particles) over 16 phases of the sources' longest period (t = 0 alone
        // without periodic sources), at least the floor; a temperature-limited cathode's at
        // most what gives TL_RATE macroparticles per unit time.
        let period = tube
            .circuits
            .iter()
            .filter_map(|c| match c.source {
                physics::circuit::Waveform::Sine { omega, .. } if omega > 0.0 => {
                    Some(std::f64::consts::TAU / omega)
                }
                physics::circuit::Waveform::Pulse { period, .. } if period > 0.0 => Some(period),
                _ => None,
            })
            .fold(0.0, f64::max);
        let phases: Vec<f64> = if period > 0.0 {
            (0..16).map(|k| period * f64::from(k) / 16.0).collect()
        } else {
            vec![0.0]
        };
        let samples: Vec<Vec<f64>> = phases
            .iter()
            .map(|&t| {
                let mut v = tube.electrodes.potentials.clone();
                for c in &tube.circuits {
                    if c.resistance.is_finite() {
                        v[c.electrode] = c.source.value(t);
                    }
                }
                v
            })
            .collect();
        for (i, fl) in floors.iter().enumerate() {
            let q = samples
                .iter()
                .map(|v| tube.emitting_charge(i, v))
                .fold(0.0, f64::max);
            let mut w = q.max(*fl).max(1e-300) / (PREVIEW_PARTICLES * r * res.particles);
            if i == 0
                && let Some(j) = spec.saturation
            {
                w = w.min(j * perimeter / (spec.particle_rate * r * res.particles));
            }
            tube.emitters[i].weight = w;
        }
        let state = tube.state();
        let potential_range = samples
            .iter()
            .flatten()
            .fold((0.0_f64, 0.0_f64), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
        let pending = TubeRun {
            segments: tube.electrodes.segments.clone(),
            potentials: tube.electrodes.potentials.clone(),
            potential_range,
            ..TubeRun::default()
        };
        let n_goals = spec.potential_goals.len();
        let n_currents = spec.goals.len();
        TubeSim {
            tube,
            state,
            spec,
            detectors,
            at_start: vec![None; n_currents],
            pending,
            steps: 0,
            currents: vec![None; n_currents],
            windows: vec![(0.0, 0.0, f64::INFINITY, f64::NEG_INFINITY); n_goals],
            series: vec![Vec::new(); n_goals],
            potential: vec![None; n_goals],
            arrivals: vec![None; n_proj],
            closest: vec![f64::INFINITY; n_proj],
        }
    }

    /// The physics circuits of a placement: each tube circuit with its source set by the
    /// player's supply on its electrode, if tunable and placed.
    pub fn tube_circuits(&self, player: &[Element]) -> Vec<physics::tube::Circuit> {
        let Some(spec) = &self.tube else {
            return Vec::new();
        };
        spec.circuits
            .iter()
            .map(|c| {
                let electrode = &self.electrodes[c.electrode];
                let supply = player.iter().find(|e| {
                    e.kind == ElementKind::Supply && electrode.tunable && e.node == electrode.center
                });
                let source = match (supply, c.source) {
                    (
                        Some(s),
                        crate::Source::Sine {
                            amplitude,
                            omega,
                            phase,
                            ..
                        },
                    ) if c.supply_offset => crate::Source::Sine {
                        offset: s.value,
                        amplitude,
                        omega,
                        phase,
                    },
                    (Some(s), source) => source.with_level(s.value),
                    (None, source) => source,
                };
                let primary = physics::tube::Circuit {
                    electrode: c.electrode,
                    source: source.waveform(),
                    resistance: c.resistance.unwrap_or(f64::INFINITY),
                    capacitance: c.capacitance,
                    inductance: c.inductance,
                    to: c.return_to,
                };
                // Linked electrodes follow the same ideal source.
                std::iter::once(primary).chain(c.linked.iter().map(move |&e| {
                    physics::tube::Circuit {
                        electrode: e,
                        ..primary
                    }
                }))
            })
            .flat_map(|it| it.collect::<Vec<_>>())
            .collect()
    }

    /// The extremes of every circuit source of a placement (over all time).
    fn tube_source_extremes(&self, player: &[Element]) -> (f64, f64) {
        let mut lo: f64 = 0.0;
        let mut hi: f64 = 0.0;
        for c in self.tube_circuits(player) {
            let (a, b) = match c.source {
                physics::circuit::Waveform::Dc(v) => (v, v),
                physics::circuit::Waveform::Sine {
                    offset, amplitude, ..
                } => (offset - amplitude.abs(), offset + amplitude.abs()),
                physics::circuit::Waveform::Pulse { low, high, .. } => {
                    (low.min(high), low.max(high))
                }
            };
            lo = lo.min(a);
            hi = hi.max(b);
        }
        (lo, hi)
    }

    /// The level's largest voltage: of its electrodes' biases and of the sliders' ranges
    /// (supplies, plates); at least 1. Sets the macroparticles' weight floor.
    pub fn tube_reference_voltage(&self) -> f64 {
        let biases = self.electrodes.iter().map(|e| match e.bias {
            crate::ConductorBias::Potential(v) => v.abs(),
            crate::ConductorBias::Grounded | crate::ConductorBias::Charge(_) => 0.0,
        });
        let lists = self
            .limits
            .supply_voltages
            .iter()
            .chain(&self.limits.plate_voltages)
            .map(|v| v.abs());
        let (lo, hi) = self.tube_source_extremes(&[]);
        biases
            .chain(lists)
            .chain([lo.abs(), hi.abs()])
            .fold(1.0, f64::max)
    }

    /// The end of a tube level's verification runs: the current goal's window end, and
    /// with projectile detectors the level's `t_max` (they must arrive by then).
    pub fn tube_verify_end(&self) -> f64 {
        self.tube.as_ref().map_or(0.0, |t| {
            let goal = t
                .goals
                .iter()
                .map(|g| g.end)
                .chain(t.potential_goals.iter().map(|g| g.end))
                .fold(0.0, f64::max);
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
    pub fn tube_results(&self, player: &[Element], refine: u32) -> TubeResults {
        let mut sim = self.tube_sim(player, refine);
        sim.advance_to(self.tube_verify_end());
        sim.results()
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
    /// current's and each potential goal's distance outside its range in units of its
    /// width (0 inside), plus each projectile's closest approach to its detector (cells),
    /// and `Arrived` when all are 0, else `Rejected` (as a particle in its detector outside
    /// its acceptance).
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_objective(&self, player: &[Element]) -> (f64, physics::trajectory::Outcome) {
        use physics::trajectory::Outcome;
        let spec = self.tube.as_ref().expect("a tube level");
        let mut sim = self.tube_sim(player, 1);
        sim.advance_to(self.tube_verify_end());
        let mut score = 0.0;
        for (g, i) in spec.goals.iter().zip(sim.currents()) {
            let i = i.unwrap_or(0.0);
            score += (g.min - i).max(i - g.max).max(0.0) / (g.max - g.min).max(1e-300);
        }
        let measured = sim.results().potential;
        for (k, g) in spec.potential_goals.iter().enumerate() {
            // Not measured (the run ended early): one width, finite for the annealing.
            score += measured.get(k).copied().flatten().map_or(1.0, |v| {
                let v = g.compared(v);
                (g.min - v).max(v - g.max).max(0.0) / (g.max - g.min).max(1e-300)
            });
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
        for g in &spec.goals {
            // A carrier's own emitter takes back only those that return: not a goal.
            let emitter = match g.carriers {
                Carriers::Cathode => Some(spec.cathode),
                Carriers::Ions => spec.ions.map(|i| i.electrode),
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
        {
            out.push(
                "Tube: no goal (a current, a potential, or a projectile with a detector).".into(),
            );
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
        if spec.saturation.is_some_and(|j| !(j.is_finite() && j > 0.0)) {
            out.push("Tube: the cathode's saturation current must be positive.".into());
        }
        if !(spec.particle_rate.is_finite() && spec.particle_rate >= 1.0) {
            out.push("Tube: the cathode's macroparticle rate must be at least 1.".into());
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
                || (!c.linked.is_empty() && (c.resistance != Some(0.0) || c.inductance != 0.0));
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
                     supply needs a sine; linked electrodes an ideal source, resistance 0 and no \
                     coil; a coil a finite resistance; a return to another electrode a \
                     resistance or a coil).",
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
                        ElementKind::Plate | ElementKind::Supply | ElementKind::Free
                    )
                }),
                "placed elements other than plates, supplies and free charges",
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The Fourier measures of a potential goal: `3 sin(ω t + φ)` plus a drift and an
    /// offset, sampled uniformly over whole periods, gives the amplitude 3 and the phase φ
    /// (degrees) to rounding, for phases all round the circle; a phase goal compares a value
    /// at its turn nearest the goal's middle.
    #[test]
    fn amplitude_and_phase() {
        let omega = 2.0 * std::f64::consts::PI / 1.5;
        for phi in [-170.0_f64, -90.0, -10.0, 0.0, 45.0, 135.0, 179.0] {
            let samples: Vec<(f64, f64)> = (0..600)
                .map(|k| {
                    let t = 10.0 + f64::from(k) * 0.025;
                    (
                        t,
                        3.0 * libm::sin(omega * t + phi.to_radians()) + 0.2 * t - 5.0,
                    )
                })
                .collect();
            let (a, p) = fourier_at(&samples, omega);
            assert!((a - 3.0).abs() < 1e-9, "amplitude {a} at {phi}");
            assert!((p - phi).abs() < 1e-7, "phase {p} for {phi}");
        }
        let goal = PotentialGoal {
            electrode: 0,
            measure: Measure::Phase,
            omega,
            min: 150.0,
            max: 210.0,
            start: 0.0,
            end: 1.0,
        };
        assert!((goal.compared(-170.0) - 190.0).abs() < 1e-12);
        assert!((goal.compared(170.0) - 170.0).abs() < 1e-12);
        let mean = PotentialGoal {
            measure: Measure::Mean,
            ..goal
        };
        assert!((mean.compared(-170.0) + 170.0).abs() < 1e-12);
    }
}
