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
use physics::lattice::Lattice;
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
    /// The solenoid's field over time (`physics::tube::SolenoidWave`, with the electric
    /// field it induces about the board's centre): this source's shape, its level set by
    /// the player's supply on the terminal, else `b_z` (a DC level, a sine's amplitude, a
    /// pulse's high: a ramp or a triangle). None: a steady field.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solenoid_drive: Option<crate::Source>,
    /// A temperature-limited cathode's macroparticles per unit time at its saturation, in
    /// the preview (the verification's twice as many; `TL_RATE` by default): a beam
    /// modulated at a high frequency needs many per period.
    #[serde(default = "default_rate", skip_serializing_if = "is_default_rate")]
    pub particle_rate: f64,
    /// This many times the macroparticles of every emitter (1 by default): a magnetized
    /// sheath needs many more than a diode's flow (PHYSICS.md §2.11; the tree code makes
    /// them affordable).
    #[serde(
        default = "default_particles",
        skip_serializing_if = "is_default_particles"
    )]
    pub particles: f64,
    /// The emitters' segments this many times finer than `EMITTER_MESH` (1 by default):
    /// their particles are released half a segment out, a first-order error a magnetized
    /// sheath magnifies (PHYSICS.md §2.11).
    #[serde(
        default = "default_particles",
        skip_serializing_if = "is_default_particles"
    )]
    pub emitter_mesh: f64,
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
    /// Goals on metals' heat over windows (the mean power the field gives their carriers).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub heat_goals: Vec<HeatGoal>,
    /// The tube's radiated output (a magnetron's microwaves): one circuit's load.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output: Option<TubeOutput>,
    /// Conductors' interiors (the Drude metal): resistors, Hall probes.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub metals: Vec<MetalSpec>,
}

/// The densest medium of cores a metal may have: `a² n_ion = a/(2ℓ)`. Beyond, the fixed
/// cores' own correction to the conductivity passes 3 % (PHYSICS.md §2.11, the Drude
/// metal: about `1.1 a²n ln(a²n)`).
pub const MAX_CORE_DENSITY: f64 = 0.005;

/// The plasma period resolved: at most this fraction of a radian of the plasma
/// oscillation per preview step (`ω_p dt`, `ω_p² = 4π ρ |q/m|`).
pub const MAX_PLASMA_STEP: f64 = 0.2;

/// The particle-frames a run with metals keeps for display, at most (about 70 MB):
/// `TubeSim` records a frame every so many steps to stay within it.
pub const FRAME_BUDGET: f64 = 4.0e6;

/// About this many of a metal run's carriers are recorded at every step between its frames
/// (`TubeFrame::paths`): the display traces a few dozen of them, and their paths bend within
/// a frame's span (Magnetoresistance at B = 7: a gyration period of 0.9 against frames 0.28
/// apart; drawn from the frames alone, a trail was two straight strokes).
pub const FINE_SAMPLE: usize = 160;

/// The most cores a metal may have (memory and the collision search's buckets).
pub const MAX_CORES: f64 = 2.0e6;

/// The most carriers the metals may have in the verification (a share of
/// `physics::tube::MAX_PARTICLES`, so that the contacts can still emit).
pub const MAX_METAL_CARRIERS: f64 = 15_000.0;

/// Offsets a metal's seed for its carriers' random numbers (its cores use the seed).
const METAL_SEED: u64 = 0x6d65_7461_6c73;

/// A conductor's interior (PHYSICS.md §2.11, the Drude metal; `physics::lattice`):
/// carriers over the ions' uniform background in an axis-aligned rectangle, scattering off
/// ion cores placed at random. Its walls reflect the carriers, except where a contact
/// touches it: a contact takes in the carriers that reach it and gives off carriers into
/// the metal under the space-charge limit (which keeps the metal next to it neutral).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MetalSpec {
    /// Two opposite corners of the rectangle.
    pub corners: [Node; 2],
    /// Its contacts (indices into the level's `electrodes`): boxes along its sides, each
    /// touching one side with its whole face.
    pub contacts: Vec<usize>,
    /// The carriers' charge density (its magnitude; the background's is opposite) and
    /// charge-to-mass ratio (electrons: negative).
    pub density: f64,
    #[serde(default = "electron_ratio")]
    pub charge_per_mass: f64,
    /// The carriers' thermal speed `√(kT/m)`.
    pub thermal_speed: f64,
    /// The mean free path `ℓ = 1/(2 a n_ion)` and the cores' radius `a`: the cores' number
    /// per area follows.
    pub mean_free_path: f64,
    pub core_radius: f64,
    /// Macroparticles per unit area in the preview (the verification's twice as many).
    pub carriers_per_area: f64,
    /// The cores' placement and the collisions' random numbers.
    #[serde(default)]
    pub seed: u64,
}

impl TubeSpec {
    /// The metal whose contact the electrode `e` is.
    #[must_use]
    pub fn metal_of(&self, e: usize) -> Option<usize> {
        self.metals.iter().position(|m| m.contacts.contains(&e))
    }

    /// A current goal's sign: positive for its carriers arriving (a lead's current: the
    /// carriers of the metal its electrode is a contact of).
    fn goal_sign(&self, g: &CurrentGoal) -> f64 {
        match g.carriers {
            Carriers::Cathode => self.charge_per_mass.signum(),
            Carriers::Ions => 1.0,
            Carriers::Net => self
                .metal_of(g.electrode)
                .map_or(self.charge_per_mass, |m| self.metals[m].charge_per_mass)
                .signum(),
        }
    }
}

impl MetalSpec {
    /// The cores' number per area.
    #[must_use]
    pub fn core_density(&self) -> f64 {
        1.0 / (2.0 * self.core_radius * self.mean_free_path)
    }

    /// The plasma frequency `√(4π ρ |q/m|)`.
    #[must_use]
    pub fn plasma_frequency(&self) -> f64 {
        libm::sqrt(4.0 * std::f64::consts::PI * self.density * self.charge_per_mass.abs())
    }
}

/// A tube's output (PHYSICS.md §2.11, the magnetron): the load of circuit `circuit` (its
/// resistance) is what an antenna radiates, the level's antenna electrode strapped to that
/// circuit's electrode. The game shows the run's electromagnetic field (`tube_em`) at the
/// speed of light `c`, a display parameter: the tube's own dynamics stay Newtonian (the
/// electrons at a fraction of `c`); an output's level is open (its board's edge absorbs
/// the waves).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct TubeOutput {
    /// The circuit (an index into `TubeSpec::circuits`) whose resistance is the load.
    pub circuit: usize,
    /// The speed of light of the radiated wave.
    pub c: f64,
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
    /// Further electrodes wired to this one: on an ideal source, a grid's wires; on any
    /// other circuit, straps (`physics::tube::Circuit::to`) making them one node with it
    /// (a magnetron's alternate vanes on one tank).
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

/// A metal's heat, the power the field gives its carriers (`TubeState::field_work`: the
/// Joule heat in a steady state), averaged over `start..end`, must lie in `min..max`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct HeatGoal {
    pub metal: usize,
    pub min: f64,
    pub max: f64,
    pub start: f64,
    pub end: f64,
}

/// An electrode's potential, measured over `start..end`, must lie in `min..max`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PotentialGoal {
    pub electrode: usize,
    /// Measured against this electrode's potential (a voltmeter across the two: a Hall
    /// probe's voltage); None: against ground.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<usize>,
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

fn default_particles() -> f64 {
    1.0
}

#[allow(clippy::trivially_copy_pass_by_ref, clippy::float_cmp)]
fn is_default_particles(v: &f64) -> bool {
    *v == 1.0
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
    /// Every carrier's: what the electrode takes in less what it gives off, the current in
    /// its lead (a metal's contact).
    Net,
}

impl Carriers {
    #[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
    fn is_cathode(&self) -> bool {
        *self == Carriers::Cathode
    }

    /// The emitter index of these carriers in a run (`TubeState::collected`'s second
    /// index): the cathode is emitter 0, the ion source 1; all of them for `Net`.
    fn emitter(self) -> Option<usize> {
        match self {
            Carriers::Cathode => Some(0),
            Carriers::Ions => Some(1),
            Carriers::Net => None,
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
/// `softening` further multiply the macroparticles' number and the softening alone, and
/// `emitter` divides the emitters' segments (their particles' release distance) alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeResolution {
    pub refine: u32,
    pub particles: f64,
    pub softening: f64,
    pub emitter: f64,
}

impl TubeResolution {
    /// The game's resolutions: everything scaled by `refine` together.
    pub fn refined(refine: u32) -> Self {
        Self {
            refine,
            particles: 1.0,
            softening: 1.0,
            emitter: 1.0,
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
    /// Each particle's identity (`physics::tube::TubeState::id`).
    pub id: Vec<u32>,
    pub sigma: Vec<f32>,
    /// Every electrode's potential.
    pub potentials: Vec<f32>,
    /// Per tube circuit (`TubeSpec::circuits`): the current through its branch (its
    /// coil's, if it has one; 0 for an ideal source).
    pub currents: Vec<f32>,
    /// The Joule heat the metals' lattices have taken so far (per unit length).
    pub heat: f64,
    /// The sampled carriers (`TubeRun::fine_every`) at every step since the previous frame,
    /// by identity (sorted): the last point is this frame's position, the one before it a
    /// step (`TubeRun::dt`) earlier, and so on. Empty where every step is a frame.
    pub paths: Vec<(u32, Vec<[f32; 2]>)>,
}

impl TubeFrame {
    /// Carrier `id`'s positions at every step since the previous frame (`paths`).
    #[must_use]
    pub fn path(&self, id: u32) -> Option<&[[f32; 2]]> {
        self.paths
            .binary_search_by_key(&id, |(i, _)| *i)
            .ok()
            .map(|k| self.paths[k].1.as_slice())
    }
}

/// A metal as the display shows it (`TubeSpec::metals`): its rectangle, its background's
/// charge per area (the ions'), its cores and its carriers' thermal speed `√(kT/m)`.
#[derive(Clone, Debug, Default)]
pub struct MetalView {
    pub lo: DVec3,
    pub hi: DVec3,
    pub background: f64,
    pub cores: Vec<[f32; 2]>,
    pub core_radius: f64,
    pub thermal_speed: f64,
}

/// Each metal's carriers' mean displacement (x, y) from frame `a` to frame `b`: the mean over
/// the carriers in it at `b` that were in flight at `a` (thousands: the thermal motion
/// averages to ~10 % of a frame's drift, and does not accumulate). Summed over frames, the
/// drift of the carriers as a whole (the game's chevrons; test `metal_drift_carries_the_current`).
#[allow(clippy::cast_precision_loss)]
#[must_use]
pub fn metal_steps(metals: &[MetalView], a: &TubeFrame, b: &TubeFrame) -> Vec<[f64; 2]> {
    let before: std::collections::HashMap<u32, [f32; 2]> =
        a.id.iter().zip(&a.x).map(|(id, x)| (*id, *x)).collect();
    let mut sums = vec![([0.0_f64; 2], 0_usize); metals.len()];
    for (x, id) in b.x.iter().zip(&b.id) {
        let inside = |m: &MetalView| {
            let (px, py) = (f64::from(x[0]), f64::from(x[1]));
            px >= m.lo.x && px <= m.hi.x && py >= m.lo.y && py <= m.hi.y
        };
        let (Some(m), Some(x0)) = (metals.iter().position(inside), before.get(id)) else {
            continue;
        };
        sums[m].0[0] += f64::from(x[0] - x0[0]);
        sums[m].0[1] += f64::from(x[1] - x0[1]);
        sums[m].1 += 1;
    }
    sums.into_iter()
        .map(|(sum, n)| {
            let n = n.max(1) as f64;
            [sum[0] / n, sum[1] / n]
        })
        .collect()
}

impl MetalView {
    /// The metal as a strip: its centre, long axis (x for a square), the normal across it,
    /// its length and width.
    #[must_use]
    pub fn strip(&self) -> (DVec3, DVec3, DVec3, f64, f64) {
        let (dx, dy) = (self.hi.x - self.lo.x, self.hi.y - self.lo.y);
        let centre = (self.lo + self.hi) * 0.5;
        if dx >= dy {
            (centre, DVec3::X, DVec3::Y, dx, dy)
        } else {
            (centre, DVec3::Y, DVec3::X, dy, dx)
        }
    }
}

/// The mean over `[a, b]` of a series `(t, v)` (times increasing) taken as linear between
/// its samples and constant beyond its ends; at `b` if the interval is empty, None if the
/// series is.
fn series_mean(series: &[(f64, f64)], a: f64, b: f64) -> Option<f64> {
    if series.is_empty() {
        return None;
    }
    let at = |t: f64| {
        let k = series.partition_point(|s| s.0 <= t);
        match (k.checked_sub(1).map(|i| series[i]), series.get(k)) {
            (Some(p), Some(q)) => p.1 + (q.1 - p.1) * (t - p.0) / (q.0 - p.0),
            (Some(p), None) => p.1,
            (None, Some(q)) => q.1,
            (None, None) => 0.0,
        }
    };
    if b <= a {
        return Some(at(b));
    }
    let k0 = series.partition_point(|s| s.0 <= a);
    let k1 = series.partition_point(|s| s.0 < b);
    let mut points = vec![(a, at(a))];
    points.extend_from_slice(&series[k0..k1.max(k0)]);
    points.push((b, at(b)));
    let integral: f64 = points
        .windows(2)
        .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
        .sum();
    Some(integral / (b - a))
}

/// The forces on a metal's average carrier in a magnetic field (`TubeRun::hall_forces`):
/// its middle's centre, long axis, the normal across it and its width, the carriers' mean
/// drift there, the magnetic push on that drift `q⟨v⟩ × B`, the electric push `q⟨E⟩` (the
/// mean field there) and the collisions' mean push, which in a steady drift balances the
/// other two (the carriers' mean momentum stays put): minus their sum.
#[derive(Clone, Copy, Debug)]
pub struct HallForces {
    pub centre: DVec3,
    pub axis: DVec3,
    pub normal: DVec3,
    pub width: f64,
    pub drift: DVec3,
    pub magnetic: DVec3,
    pub electric: DVec3,
    pub collisions: DVec3,
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
    /// Each heat goal's mean power, once the run has passed its window.
    pub heat: Vec<Option<f64>>,
    /// Per metal piece (`metals`): the work the field has done on its carriers so far,
    /// `(t, W)`, at each frame (`TubeState::field_work`; its rate is the heat).
    pub work: Vec<Vec<(f64, f64)>>,
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
    /// Each segment's electrode.
    pub owners: Vec<usize>,
    /// Each prism's potential.
    pub potentials: Vec<f64>,
    /// The lowest and highest potential the electrodes reach: their own, and their circuits'
    /// sources over a period (the potential map's scale; at t = 0 alone an AC level's
    /// electrodes are all at 0).
    pub potential_range: (f64, f64),
    /// The metals (their backgrounds are fixed charges of the maps).
    pub metals: Vec<MetalView>,
    /// The uniform magnetic field along z of the run (the level's, or its solenoid's
    /// supply).
    pub b_z: f64,
    /// A field that changes in time (`TubeSpec::solenoid_drive`), in place of `b_z`.
    pub b_wave: Option<physics::tube::SolenoidWave>,
    /// Per metal, in a field: the mean field in its middle at each frame, its x and its y
    /// component `(t, E)` (`middle_field`; the Hall diagram's electric push is its mean over
    /// the drift's time).
    pub middle: Vec<[Vec<(f64, f64)>; 2]>,
    /// The time step.
    pub dt: f64,
    /// The carriers recorded at every step (`TubeFrame::paths`): those whose identity is a
    /// multiple of this; 0 where every step is a frame.
    pub fine_every: u32,
}

impl TubeRun {
    /// B along z at `t`.
    #[must_use]
    pub fn b_at(&self, t: f64) -> f64 {
        self.b_wave.map_or(self.b_z, |w| w.b(t))
    }

    /// Whether the run has a magnetic field at all (steady or changing).
    #[must_use]
    pub fn has_field(&self) -> bool {
        self.b_z != 0.0 || self.b_wave.is_some()
    }

    /// Appends a streamed chunk (its frames and samples follow this run's).
    pub fn extend(&mut self, chunk: TubeRun) {
        self.frames.extend(chunk.frames);
        if self.collected.len() < chunk.collected.len() {
            self.collected.resize(chunk.collected.len(), Vec::new());
        }
        for (a, c) in self.collected.iter_mut().zip(chunk.collected) {
            a.extend(c);
        }
        if self.middle.len() < chunk.middle.len() {
            self.middle
                .resize(chunk.middle.len(), [Vec::new(), Vec::new()]);
        }
        for (a, [cx, cy]) in self.middle.iter_mut().zip(chunk.middle) {
            a[0].extend(cx);
            a[1].extend(cy);
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
        if self.heat.len() < chunk.heat.len() {
            self.heat.resize(chunk.heat.len(), None);
        }
        for (a, c) in self.heat.iter_mut().zip(chunk.heat) {
            *a = a.or(c);
        }
        if self.work.len() < chunk.work.len() {
            self.work.resize(chunk.work.len(), Vec::new());
        }
        for (a, c) in self.work.iter_mut().zip(chunk.work) {
            a.extend(c);
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

    /// The picture at time `t` between the stored frames, for the display (smooth at any
    /// playback speed): each particle in both neighbouring frames moved linearly between
    /// them, the surface densities and potentials interpolated; a particle absorbed before
    /// the next frame stays where it was, one emitted after it appears with that frame.
    /// Before the first frame the first, after the last the last.
    pub fn frame_between(&self, t: f64) -> Option<TubeFrame> {
        let k = self.frames.partition_point(|f| f.t <= t).saturating_sub(1);
        let a = self.frames.get(k)?;
        let Some(b) = self.frames.get(k + 1).filter(|b| b.t > a.t && t > a.t) else {
            return Some(a.clone());
        };
        #[allow(clippy::cast_possible_truncation)]
        let w = ((t - a.t) / (b.t - a.t)).clamp(0.0, 1.0) as f32;
        let next: std::collections::HashMap<u32, usize> =
            b.id.iter().enumerate().map(|(i, id)| (*id, i)).collect();
        let lerp = |u: f32, v: f32| u + (v - u) * w;
        let x =
            a.x.iter()
                .enumerate()
                .map(|(i, p)| match a.id.get(i).and_then(|id| next.get(id)) {
                    // A sampled carrier along its steps, the others straight on.
                    Some(&j) => b
                        .path(a.id[i])
                        .and_then(|path| self.along(*p, a.t, path, b.t, t))
                        .unwrap_or([lerp(p[0], b.x[j][0]), lerp(p[1], b.x[j][1])]),
                    None => *p,
                })
                .collect();
        Some(TubeFrame {
            t,
            x,
            q: a.q.clone(),
            kind: a.kind.clone(),
            id: a.id.clone(),
            sigma: a
                .sigma
                .iter()
                .zip(&b.sigma)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            potentials: a
                .potentials
                .iter()
                .zip(&b.potentials)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            currents: a
                .currents
                .iter()
                .zip(&b.currents)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            heat: a.heat + (b.heat - a.heat) * f64::from(w),
            paths: Vec::new(),
        })
    }

    /// The position at `t` of a carrier at `x0` at `t0` that then took the steps `path`
    /// (the last at `t1`, each a step `dt` after the one before): between its two steps
    /// around `t`. None if the path does not span the time from `t0` to `t1`.
    fn along(&self, x0: [f32; 2], t0: f64, path: &[[f32; 2]], t1: f64, t: f64) -> Option<[f32; 2]> {
        let n = path.len();
        if n == 0 || self.dt <= 0.0 {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        let start = t1 - self.dt * n as f64;
        if (start - t0).abs() > 0.25 * self.dt {
            return None;
        }
        // Steps after `start`: 0 at x0, k at path[k - 1].
        #[allow(clippy::cast_precision_loss)]
        let u = ((t - start) / self.dt).clamp(0.0, n as f64);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (libm::floor(u) as usize).min(n - 1);
        let p0 = if k == 0 { x0 } else { path[k - 1] };
        let p1 = path[k];
        #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
        let w = (u - k as f64) as f32;
        Some([p0[0] + (p1[0] - p0[0]) * w, p0[1] + (p1[1] - p0[1]) * w])
    }

    /// Carrier `id`'s path from `t0` to `t1`: its recorded steps where it is sampled
    /// (`TubeFrame::paths`), else its frames' positions, ending where it is at `t1`
    /// (`frame_between`); empty if it is in no frame of that time.
    #[must_use]
    pub fn trail(&self, id: u32, t0: f64, t1: f64) -> Vec<[f32; 2]> {
        let k0 = self.frames.partition_point(|f| f.t < t0);
        let k1 = self.frames.partition_point(|f| f.t <= t1);
        let mut out = Vec::new();
        // The frames up to t1, and the next one's steps before t1.
        for f in &self.frames[k0..(k1 + 1).min(self.frames.len())] {
            match f.path(id) {
                Some(path) =>
                {
                    #[allow(clippy::cast_precision_loss)]
                    for (i, x) in path.iter().enumerate() {
                        let ti = f.t - self.dt * (path.len() - 1 - i) as f64;
                        if ti >= t0 && ti <= t1 {
                            out.push(*x);
                        }
                    }
                }
                None if f.t <= t1 => {
                    if let Some(i) = f.id.iter().position(|j| *j == id) {
                        out.push(f.x[i]);
                    }
                }
                None => {}
            }
        }
        if let Some(now) = self.frame_between(t1)
            && let Some(i) = now.id.iter().position(|j| *j == id)
        {
            out.push(now.x[i]);
        }
        out
    }

    /// The frame at time `t` (the last one at or before it; the first before it).
    pub fn frame_at(&self, t: f64) -> Option<&TubeFrame> {
        let k = self.frames.partition_point(|f| f.t <= t);
        self.frames.get(k.saturating_sub(1))
    }

    /// The field of a frame at `x` as `field_at` gives it, faster: far segments by their
    /// multipoles (`physics::zinv::segment_field_fast`, within 5e-5): for the display's
    /// arrows and colour scale, many points every frame.
    pub fn e_at(&self, frame: &TubeFrame, x: DVec3) -> DVec3 {
        let mut e = DVec3::ZERO;
        for ((a, b), s) in self.segments.iter().zip(&frame.sigma) {
            e += physics::zinv::segment_field_fast(*a, *b, x) * f64::from(*s);
        }
        let h2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        for (p, q) in frame.x.iter().zip(&frame.q) {
            let r = DVec3::new(x.x - f64::from(p[0]), x.y - f64::from(p[1]), 0.0);
            e += r * (2.0 * f64::from(*q) / (r.x * r.x + r.y * r.y).max(h2));
        }
        for m in &self.metals {
            e += physics::zinv::rect_integrals(m.lo, m.hi, x).1 * m.background;
        }
        // A changing field's induced part (no potential: it circulates).
        if let Some(w) = &self.b_wave {
            e += w.induced(frame.t, x);
        }
        e
    }

    /// The mean field in metal `metal`'s middle on `frame` (`MetalView::strip`): the mean of
    /// `e_at` over the central three fifths along it (7 points) on lines half a unit apart
    /// from a unit inside each side (its sides hold the surface charge of its potential
    /// against its surroundings, within a Debye length; the Hall field is the bulk's), or on
    /// its centre line alone in a strip narrower than 3. Its carriers' graininess scatters
    /// the part across from frame to frame: on Hall probe's reference (100 macroparticles per
    /// area) by 0.11 against its 0.47, on the centre line alone by 0.19; on Positive
    /// carriers' lower strip (50 per area, a weaker current) by 0.20, on the centre line by
    /// 0.35, against 0.27.
    #[must_use]
    pub fn middle_field(&self, frame: &TubeFrame, metal: usize) -> DVec3 {
        let Some(m) = self.metals.get(metal) else {
            return DVec3::ZERO;
        };
        let (centre, axis, normal, length, width) = m.strip();
        #[allow(clippy::cast_possible_truncation)]
        let lines = libm::floor((0.5 * width - 1.0) / 0.5 + 1e-9).max(0.0) as i32;
        let (mut sum, mut n) = (DVec3::ZERO, 0.0);
        for l in -lines..=lines {
            for i in 0..7 {
                let x = centre
                    + axis * (0.3 * length * (f64::from(i) / 3.0 - 1.0))
                    + normal * (0.5 * f64::from(l));
                sum += self.e_at(frame, x);
                n += 1.0;
            }
        }
        sum / n
    }

    /// The Hall effect in metal `metal` of a run in a field along z: its carriers' mean
    /// drift from frame `before` to frame `now` in the strip's middle (the central three
    /// fifths along it, a quarter away from its sides), the field's magnetic push on that
    /// drift, `q ⟨v⟩ × B`, and the electric push across the strip, `q E⊥`: the field across
    /// its middle (`hall_field`) averaged over the same time (the run's series `hall`,
    /// linear between its frames: over a time unit its scatter falls to a few per cent),
    /// per unit charge of the carriers' sign. None without a field or with too few
    /// carriers.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn hall_forces(
        &self,
        now: &TubeFrame,
        before: &TubeFrame,
        metal: usize,
    ) -> Option<HallForces> {
        let m = self.metals.get(metal)?;
        let span = now.t - before.t;
        if !self.has_field() || span <= 0.0 {
            return None;
        }
        let (centre, axis, normal, length, width) = m.strip();
        let then: std::collections::HashMap<u32, [f32; 2]> = before
            .id
            .iter()
            .zip(&before.x)
            .map(|(id, x)| (*id, *x))
            .collect();
        let (mut sum, mut n, mut sign) = (DVec3::ZERO, 0.0, 0.0);
        for ((x, id), kind) in now.x.iter().zip(&now.id).zip(&now.kind) {
            let p = DVec3::new(f64::from(x[0]), f64::from(x[1]), 0.0);
            let d = p - centre;
            if d.dot(axis).abs() >= 0.3 * length || d.dot(normal).abs() >= 0.5 * width - 0.25 {
                continue;
            }
            let Some(x0) = then.get(id) else {
                continue;
            };
            sum += p - DVec3::new(f64::from(x0[0]), f64::from(x0[1]), 0.0);
            n += 1.0;
            sign += match *kind {
                KIND_NEGATIVE => -1.0,
                KIND_POSITIVE => 1.0,
                _ => 0.0,
            };
        }
        if n < 10.0 {
            return None;
        }
        let q = f64::signum(sign);
        let drift = sum / (n * span);
        // q v × B ẑ = q B (v_y, −v_x).
        // A changing field: B at the window's middle.
        let b = self.b_at(0.5 * (before.t + now.t));
        let magnetic = DVec3::new(drift.y, -drift.x, 0.0) * (q * b);
        let mean = |k: usize| {
            self.middle
                .get(metal)
                .and_then(|s| series_mean(&s[k], before.t, now.t))
        };
        let field = match (mean(0), mean(1)) {
            (Some(x), Some(y)) => DVec3::new(x, y, 0.0),
            _ => self.middle_field(now, metal),
        };
        let electric = field * q;
        Some(HallForces {
            centre,
            axis,
            normal,
            width,
            drift,
            magnetic,
            electric,
            collisions: -(electric + magnetic),
        })
    }

    /// The potential and field of a frame at `x`, as the maps show them: the segments'
    /// charge exactly, the particles as uniformly charged rods of radius `PREVIEW_SEGMENT`
    /// (the preview's resolution of its space charge; the physics' rods are thinner,
    /// `SOFTENING`), exactly Coulomb's beyond it. (Softened as `r⃗/(r² + ε²)`, which never
    /// becomes Coulomb's, a metal's carriers failed to cancel its background's exact field
    /// by `2λε²/r³` each: across Positive carriers' lower strip, from the upper's 128
    /// charges, a field that took the Hall arrows' electric push to 1.73 of the magnetic.)
    pub fn field_at(&self, frame: &TubeFrame, x: DVec3) -> (f64, DVec3) {
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        for ((a, b), s) in self.segments.iter().zip(&frame.sigma) {
            let (f, g) = physics::zinv::segment_integrals(*a, *b, x);
            phi += f64::from(*s) * f;
            e += g * f64::from(*s);
        }
        let h2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        let log_h2 = libm::log(h2);
        for (p, q) in frame.x.iter().zip(&frame.q) {
            let r = DVec3::new(x.x - f64::from(p[0]), x.y - f64::from(p[1]), 0.0);
            let r2 = r.x * r.x + r.y * r.y;
            let lambda = f64::from(*q);
            if r2 >= h2 {
                phi -= lambda * libm::log(r2);
                e += r * (2.0 * lambda / r2);
            } else {
                phi -= lambda * (log_h2 - 1.0 + r2 / h2);
                e += r * (2.0 * lambda / h2);
            }
        }
        for m in &self.metals {
            let (f, g) = physics::zinv::rect_integrals(m.lo, m.hi, x);
            phi += m.background * f;
            e += g * m.background;
        }
        if let Some(w) = &self.b_wave {
            e += w.induced(frame.t, x);
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
    /// The run's segments, metals and field without its frames (the first chunk's: the
    /// frames' Hall field, `TubeRun::hall_field`, is measured on them).
    layout: TubeRun,
    /// Each lattice's metal (`Level::metal_pieces`: a cut metal has several).
    lattice_metal: Vec<usize>,
    steps: u32,
    /// A frame every this many steps (`FRAME_BUDGET`).
    frame_stride: u32,
    /// The sampled carriers' steps since the last frame (`TubeFrame::paths`).
    fine: std::collections::BTreeMap<u32, Vec<[f32; 2]>>,
    currents: Vec<Option<f64>>,
    /// Each potential goal's window so far: sum, count, least, most, and the samples
    /// `(t, V)` (for `Measure::Amplitude`).
    windows: Vec<(f64, f64, f64, f64)>,
    series: Vec<Vec<(f64, f64)>>,
    potential: Vec<Option<f64>>,
    /// Per heat goal: the time and its metal's work at its window's start, and its result.
    heat_start: Vec<Option<(f64, f64)>>,
    heat: Vec<Option<f64>>,
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

    /// The field at `x` now, as the carriers feel it (the physics' softening).
    #[must_use]
    pub fn field_at(&self, x: DVec3) -> DVec3 {
        self.tube.field_at(&self.state, x)
    }

    /// The macroparticles' positions and charges now.
    #[must_use]
    pub fn particles(&self) -> (&[DVec3], &[f64]) {
        (&self.state.x, &self.state.charge)
    }

    /// Macroparticles in flight now.
    pub fn in_flight(&self) -> usize {
        self.state.x.len()
    }

    /// The electrodes' potentials now.
    pub fn potentials(&self) -> &[f64] {
        &self.state.potentials
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
                // The carriers' charge collected, positive for arrivals (a lead's: what
                // its electrode takes in less what it gives off).
                let sign = self.spec.goal_sign(g);
                let collected = sign
                    * match g.carriers.emitter() {
                        Some(i) => s.collected[g.electrode].get(i).copied().unwrap_or(0.0),
                        None => {
                            s.collected[g.electrode].iter().sum::<f64>()
                                - self
                                    .tube
                                    .emitters
                                    .iter()
                                    .zip(&s.released)
                                    .filter(|(em, _)| em.electrode == g.electrode)
                                    .map(|(_, q)| q)
                                    .sum::<f64>()
                        }
                    };
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
            // The sampled carriers every step (between frames).
            let every = self.layout.fine_every;
            if every > 0 {
                for ((x, id), tag) in s.x.iter().zip(&s.id).zip(&s.tag) {
                    if *tag < PROJECTILE_TAG && id.is_multiple_of(every) {
                        #[allow(clippy::cast_possible_truncation)]
                        self.fine
                            .entry(*id)
                            .or_default()
                            .push([x.x as f32, x.y as f32]);
                    }
                }
            }
            // A frame every `frame_stride` steps (a metal's thousands of carriers every
            // step would fill gigabytes over a long run).
            if self.steps.is_multiple_of(self.frame_stride) {
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
                    id: s.id.clone(),
                    currents: {
                        // Each tube circuit's first physics circuit (its linked electrodes'
                        // follow it).
                        let mut at = 0;
                        self.spec
                            .circuits
                            .iter()
                            .map(|c| {
                                let k = at;
                                at += 1 + c.linked.len();
                                let i = if self.tube.circuits[k].inductance > 0.0 {
                                    s.inductor_current[k]
                                } else {
                                    s.circuit_current[k]
                                };
                                i as f32
                            })
                            .collect()
                    },
                    heat: s.joule,
                    paths: std::mem::take(&mut self.fine).into_iter().collect(),
                };
                let metals = self.layout.metals.len();
                if metals > 0 {
                    self.pending.work.resize(metals, Vec::new());
                    for m in 0..metals {
                        let w = self.state.field_work.get(m).copied().unwrap_or(0.0);
                        self.pending.work[m].push((frame.t, w));
                    }
                }
                if self.layout.has_field() && metals > 0 {
                    self.pending.middle.resize(metals, [Vec::new(), Vec::new()]);
                    for m in 0..metals {
                        let e = self.layout.middle_field(&frame, m);
                        self.pending.middle[m][0].push((frame.t, e.x));
                        self.pending.middle[m][1].push((frame.t, e.y));
                    }
                }
                self.pending.frames.push(frame);
            }
            // Each heat goal's window.
            self.heat_goals();
            let s = &self.state;
            // Each potential goal's measure over its window.
            if self.pending.potential.len() < self.potential.len() {
                self.pending.potential.resize(self.potential.len(), None);
            }
            for (k, g) in self.spec.potential_goals.iter().enumerate() {
                if s.t < g.start - 0.5 * self.tube.dt || self.potential[k].is_some() {
                    continue;
                }
                let v = s.potentials[g.electrode] - g.reference.map_or(0.0, |r| s.potentials[r]);
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

    /// Metal `m`'s field work so far (the sum over its pieces).
    fn metal_work(&self, m: usize) -> f64 {
        self.lattice_metal
            .iter()
            .zip(&self.state.field_work)
            .filter(|(k, _)| **k == m)
            .map(|(_, w)| w)
            .sum()
    }

    /// The heat goals' windows after a step: each one's start, and its mean power at its end.
    fn heat_goals(&mut self) {
        let (t, dt) = (self.state.t, self.tube.dt);
        if self.pending.heat.len() < self.heat.len() {
            self.pending.heat.resize(self.heat.len(), None);
        }
        for k in 0..self.spec.heat_goals.len() {
            let g = self.spec.heat_goals[k];
            if self.heat[k].is_some() {
                continue;
            }
            let w = self.metal_work(g.metal);
            if self.heat_start[k].is_none() && t >= g.start - 0.5 * dt {
                self.heat_start[k] = Some((t, w));
            }
            if t >= g.end - 0.5 * dt {
                let (t0, w0) = self.heat_start[k].unwrap_or((0.0, 0.0));
                let p = if t > t0 { (w - w0) / (t - t0) } else { 0.0 };
                self.heat[k] = Some(p);
                self.pending.heat[k] = Some(p);
            }
        }
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
    /// Each heat goal's mean power.
    pub heat: Vec<Option<f64>>,
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
            heat: self.heat.clone(),
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
    /// Each heat goal's verdict (the same rule).
    pub heat: Vec<CurrentVerdict>,
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
        let heat: Vec<CurrentVerdict> = spec
            .heat_goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    preview.heat.get(k).copied().flatten().unwrap_or(0.0),
                    fine.heat.get(k).copied().flatten().unwrap_or(0.0),
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
            .chain(heat.iter())
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
            heat,
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

    /// A metal's rectangle (lower and upper corners).
    #[must_use]
    pub fn metal_rect(&self, m: &MetalSpec) -> (DVec3, DVec3) {
        let (a, b) = (
            self.grid.position(m.corners[0]),
            self.grid.position(m.corners[1]),
        );
        (
            DVec3::new(a.x.min(b.x), a.y.min(b.y), 0.0),
            DVec3::new(a.x.max(b.x), a.y.max(b.y), 0.0),
        )
    }

    /// An axis-aligned box electrode's corners (None for another shape or angle).
    fn electrode_box(&self, e: usize) -> Option<(DVec3, DVec3)> {
        let el = self.electrodes.get(e)?;
        if el.shape != crate::ElectrodeShape::Box {
            return None;
        }
        let quarter = libm::round(el.angle_deg / 90.0);
        if (el.angle_deg - 90.0 * quarter).abs() > 1e-9 {
            return None;
        }
        #[allow(clippy::cast_possible_truncation)]
        let along_x = (quarter as i64).rem_euclid(2) == 0;
        let (hx, hy) = if along_x {
            (0.5 * el.length, 0.5 * el.thickness)
        } else {
            (0.5 * el.thickness, 0.5 * el.length)
        };
        let c = self.grid.position(el.center);
        Some((
            DVec3::new(c.x - hx, c.y - hy, 0.0),
            DVec3::new(c.x + hx, c.y + hy, 0.0),
        ))
    }

    /// The direction from contact `e` into metal `m`, if `e` is an axis-aligned box whose
    /// face lies along one side of the metal, within it (it then gives off carriers into
    /// the metal only).
    #[must_use]
    pub fn contact_side(&self, m: &MetalSpec, e: usize) -> Option<DVec3> {
        let (b0, b1) = self.electrode_box(e)?;
        let (lo, hi) = self.metal_rect(m);
        let tol = 1e-9 * (1.0 + hi.x - lo.x + hi.y - lo.y);
        let within = |a0: f64, a1: f64, s0: f64, s1: f64| a0 >= s0 - tol && a1 <= s1 + tol;
        if (b1.x - lo.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            Some(DVec3::X)
        } else if (b0.x - hi.x).abs() < tol && within(b0.y, b1.y, lo.y, hi.y) {
            Some(-DVec3::X)
        } else if (b1.y - lo.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            Some(DVec3::Y)
        } else if (b0.y - hi.y).abs() < tol && within(b0.x, b1.x, lo.x, hi.x) {
            Some(-DVec3::Y)
        } else {
            None
        }
    }

    /// The depth of the layer a metal's contacts keep neutral: half a mean free path, at
    /// most a quarter of the metal's width.
    #[must_use]
    pub fn contact_depth(&self, m: &MetalSpec) -> f64 {
        let (lo, hi) = self.metal_rect(m);
        (0.5 * m.mean_free_path).min(0.25 * (hi.x - lo.x).min(hi.y - lo.y))
    }

    /// The metals cut by the player's cuts (`ElementKind::Cut`): per piece, its metal and
    /// rectangle, in the metals' order; an uncut metal is one piece. A cut along x (angle 0)
    /// is the line through its node parallel to x across the whole metal it lies in, along
    /// y (90) parallel to y; together they cut each metal into a grid of rectangles.
    #[must_use]
    pub fn metal_pieces(&self, player: &[Element]) -> Vec<(usize, DVec3, DVec3)> {
        let Some(spec) = &self.tube else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (k, m) in spec.metals.iter().enumerate() {
            let (lo, hi) = self.metal_rect(m);
            let (mut xs, mut ys) = (vec![lo.x, hi.x], vec![lo.y, hi.y]);
            for e in player.iter().filter(|e| e.kind == ElementKind::Cut) {
                let p = self.grid.position(e.node);
                if !(p.x > lo.x && p.x < hi.x && p.y > lo.y && p.y < hi.y) {
                    continue;
                }
                if cut_along_x(e.angle_deg) {
                    ys.push(p.y);
                } else {
                    xs.push(p.x);
                }
            }
            xs.sort_by(f64::total_cmp);
            xs.dedup();
            ys.sort_by(f64::total_cmp);
            ys.dedup();
            for j in 0..ys.len() - 1 {
                for i in 0..xs.len() - 1 {
                    out.push((
                        k,
                        DVec3::new(xs[i], ys[j], 0.0),
                        DVec3::new(xs[i + 1], ys[j + 1], 0.0),
                    ));
                }
            }
        }
        out
    }

    /// A cut's slit (`ElementKind::Cut`): its two ends, on the edges of the metal it lies in;
    /// None if it lies in none.
    #[must_use]
    pub fn cut_segment(&self, e: &Element) -> Option<(DVec3, DVec3)> {
        let spec = self.tube.as_ref()?;
        let p = self.grid.position(e.node);
        let (lo, hi) = spec
            .metals
            .iter()
            .map(|m| self.metal_rect(m))
            .find(|(lo, hi)| p.x > lo.x && p.x < hi.x && p.y > lo.y && p.y < hi.y)?;
        Some(if cut_along_x(e.angle_deg) {
            (DVec3::new(lo.x, p.y, 0.0), DVec3::new(hi.x, p.y, 0.0))
        } else {
            (DVec3::new(p.x, lo.y, 0.0), DVec3::new(p.x, hi.y, 0.0))
        })
    }

    /// The player's cuts: each strictly inside one of the metals, along x or y, and every
    /// piece at least four of its mean free paths wide (the metals' own rule: its carriers
    /// scatter off its cores, not its walls).
    pub(crate) fn check_cuts(&self, player: &[Element]) -> Result<(), crate::PlacementError> {
        let cuts: Vec<&Element> = player
            .iter()
            .filter(|e| e.kind == ElementKind::Cut)
            .collect();
        if cuts.is_empty() {
            return Ok(());
        }
        let Some(spec) = &self.tube else {
            return Err(crate::PlacementError::CutNotInMetal(cuts[0].node));
        };
        for e in &cuts {
            let p = self.grid.position(e.node);
            let quarter = libm::round(e.angle_deg / 90.0);
            let inside = spec.metals.iter().any(|m| {
                let (lo, hi) = self.metal_rect(m);
                p.x > lo.x && p.x < hi.x && p.y > lo.y && p.y < hi.y
            });
            if (e.angle_deg - 90.0 * quarter).abs() > 1e-9 || !inside {
                return Err(crate::PlacementError::CutNotInMetal(e.node));
            }
        }
        for (k, lo, hi) in self.metal_pieces(player) {
            let width = (hi.x - lo.x).min(hi.y - lo.y);
            if width < 4.0 * spec.metals[k].mean_free_path {
                let node = cuts
                    .iter()
                    .find(|e| {
                        let p = self.grid.position(e.node);
                        p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y
                    })
                    .map_or(cuts[0].node, |e| e.node);
                return Err(crate::PlacementError::PieceTooNarrow(node));
            }
        }
        Ok(())
    }

    /// A piece `lo`..`hi` of a metal with `many` times its preview's carriers; piece 0 of an
    /// uncut metal is the metal (its seeds unchanged), further pieces seeded apart.
    fn metal_lattice(&self, m: &MetalSpec, lo: DVec3, hi: DVec3, piece: u64, many: f64) -> Lattice {
        let seed = m.seed.wrapping_add(piece.wrapping_mul(0x9e37_79b9));
        let area = (hi.x - lo.x) * (hi.y - lo.y);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let cores = libm::round(area * m.core_density()) as usize;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let carriers = libm::round(area * m.carriers_per_area * many).max(1.0) as usize;
        Lattice::new(
            lo,
            hi,
            -m.charge_per_mass.signum() * m.density,
            Lattice::random_cores(lo, hi, cores, m.core_radius, seed),
            m.core_radius,
            m.charge_per_mass,
            m.thermal_speed,
            carriers,
            seed.wrapping_add(METAL_SEED),
        )
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
                let emits = e == spec.cathode
                    || spec.ions.is_some_and(|i| i.electrode == e)
                    || spec.metal_of(e).is_some();
                if emits {
                    size * EMITTER_MESH / (res.emitter * spec.emitter_mesh)
                } else {
                    size
                }
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
        // The metals, and their contacts as emitters of their carriers (at the carriers'
        // weight, toward the metal).
        let many = r * res.particles * spec.particles;
        // The metals cut into pieces (the player's cuts): a lattice each.
        let pieces = self.metal_pieces(player);
        let mut counted = vec![0_u64; spec.metals.len()];
        let lattices: Vec<Lattice> = pieces
            .iter()
            .map(|&(k, lo, hi)| {
                counted[k] += 1;
                self.metal_lattice(&spec.metals[k], lo, hi, counted[k] - 1, many)
            })
            .collect();
        let lattice_metal: Vec<usize> = pieces.iter().map(|p| p.0).collect();
        let contact = |e: usize| {
            let k = spec.metal_of(e)?;
            let m = &spec.metals[k];
            // Every piece of a metal has its carriers per area: the first one's weight.
            let first = lattice_metal.iter().position(|&j| j == k)?;
            Some(Emitter {
                electrode: e,
                charge_per_mass: m.charge_per_mass,
                weight: lattices[first].carrier_charge().abs(),
                emit_toward: self.contact_side(m, e),
                saturation: None,
                thermal_speed: m.thermal_speed,
                ohmic: Some(physics::tube::OhmicContact {
                    density: m.density,
                    depth: self.contact_depth(m),
                }),
            })
        };
        let mut floors = vec![floor(spec.cathode)];
        let mut emitters = vec![contact(spec.cathode).unwrap_or(Emitter {
            electrode: spec.cathode,
            charge_per_mass: spec.charge_per_mass,
            weight: 1.0,
            emit_toward: toward(spec.emit_toward_deg),
            saturation: spec.saturation,
            thermal_speed: 0.0,
            ohmic: None,
        })];
        if let Some(ions) = spec.ions {
            floors.push(floor(ions.electrode));
            emitters.push(Emitter {
                electrode: ions.electrode,
                charge_per_mass: ions.charge_per_mass,
                weight: 1.0,
                emit_toward: toward(ions.emit_toward_deg),
                saturation: None,
                thermal_speed: 0.0,
                ohmic: None,
            });
        }
        for m in &spec.metals {
            for &e in &m.contacts {
                if e != spec.cathode && spec.ions.is_none_or(|i| i.electrode != e) {
                    emitters.extend(contact(e));
                }
            }
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
            b_wave: spec.solenoid_drive.map(|s| physics::tube::SolenoidWave {
                waveform: s.with_level(b_z).waveform(),
                // The solenoid around the board: its axis at the board's centre.
                axis: DVec3::new(0.5 * max.x, 0.5 * max.y, 0.0),
            }),
            circuits: self.tube_circuits(player),
            lattices,
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
            if spec.metal_of(tube.emitters[i].electrode).is_some() {
                // A contact: its metal's carriers' weight.
                continue;
            }
            let q = samples
                .iter()
                .map(|v| tube.emitting_charge(i, v))
                .fold(0.0, f64::max);
            let mut w = q.max(*fl).max(1e-300) / (PREVIEW_PARTICLES * many);
            if i == 0
                && let Some(j) = spec.saturation
            {
                w = w.min(j * perimeter / (spec.particle_rate * many));
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
            owners: (0..tube.electrodes.segments.len())
                .map(|j| tube.electrodes.owner_of(j))
                .collect(),
            potentials: tube.electrodes.potentials.clone(),
            potential_range,
            #[allow(clippy::cast_possible_truncation)]
            metals: tube
                .lattices
                .iter()
                .map(|l| MetalView {
                    lo: l.lo,
                    hi: l.hi,
                    background: l.background,
                    cores: l.cores.iter().map(|c| [c.x as f32, c.y as f32]).collect(),
                    core_radius: l.core_radius,
                    thermal_speed: l.thermal_speed,
                })
                .collect(),
            b_z: tube.b_z,
            b_wave: tube.b_wave,
            ..TubeRun::default()
        };
        let n_goals = spec.potential_goals.len();
        let n_heat = spec.heat_goals.len();
        let n_currents = spec.goals.len();
        // A metal's run keeps a frame every few steps: its particles over the whole run
        // within `FRAME_BUDGET` (other tubes: every step).
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let frame_stride = if spec.metals.is_empty() {
            1
        } else {
            let steps = self.physics.t_max / tube.dt;
            (steps * state.x.len() as f64 / FRAME_BUDGET)
                .ceil()
                .max(1.0) as u32
        };
        // Between frames, a sample of the carriers at every step (identities are dealt out
        // in order: every so many of them).
        let mut pending = pending;
        pending.dt = tube.dt;
        pending.fine_every = if frame_stride > 1 {
            u32::try_from((state.x.len() / FINE_SAMPLE).max(1)).unwrap_or(u32::MAX)
        } else {
            0
        };
        TubeSim {
            tube,
            state,
            spec,
            detectors,
            at_start: vec![None; n_currents],
            layout: pending.clone(),
            lattice_metal,
            pending,
            steps: 0,
            frame_stride,
            fine: std::collections::BTreeMap::new(),
            currents: vec![None; n_currents],
            windows: vec![(0.0, 0.0, f64::INFINITY, f64::NEG_INFINITY); n_goals],
            series: vec![Vec::new(); n_goals],
            potential: vec![None; n_goals],
            heat_start: vec![None; n_heat],
            heat: vec![None; n_heat],
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
                // Linked electrodes follow the same ideal source, or are strapped to the
                // circuit's electrode.
                let ideal =
                    primary.resistance == 0.0 && primary.inductance == 0.0 && primary.to.is_none();
                std::iter::once(primary).chain(c.linked.iter().map(move |&e| {
                    if ideal {
                        physics::tube::Circuit {
                            electrode: e,
                            ..primary
                        }
                    } else {
                        physics::tube::Circuit {
                            electrode: e,
                            source: physics::circuit::Waveform::Dc(0.0),
                            resistance: 0.0,
                            capacitance: 0.0,
                            inductance: 0.0,
                            to: Some(primary.electrode),
                        }
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
                .chain(t.heat_goals.iter().map(|g| g.end))
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
        let heat = sim.results().heat;
        for (k, g) in spec.heat_goals.iter().enumerate() {
            score += heat.get(k).copied().flatten().map_or(1.0, |v| {
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

    /// What a tube level's metals refuse (PHYSICS.md §2.11, the Drude metal): a
    /// rectangle off the board or crossing another metal or an electrode, contacts that
    /// do not lie along a side, a medium of cores too dense for the model, a plasma period
    /// or a mean free path the step does not resolve, too many carriers.
    #[allow(clippy::too_many_lines)]
    fn metal_issues(&self, spec: &TubeSpec, player: &[Element], out: &mut Vec<String>) {
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
            let area = (hi.x - lo.x) * (hi.y - lo.y);
            carriers += area * m.carriers_per_area;
            #[allow(clippy::cast_precision_loss)]
            let cores = area * m.core_density();
            if cores > MAX_CORES {
                out.push(format!("{name}: {cores:.0} cores, at most {MAX_CORES:.0}."));
            }
            if m.contacts.is_empty() {
                out.push(format!("{name} needs a contact."));
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
        if spec.saturation.is_some_and(|j| !(j.is_finite() && j > 0.0)) {
            out.push("Tube: the cathode's saturation current must be positive.".into());
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

/// Whether a cut at `angle_deg` runs along x (0 or 180), else along y.
#[must_use]
pub fn cut_along_x(angle_deg: f64) -> bool {
    #[allow(clippy::cast_possible_truncation)]
    let quarter = libm::round(angle_deg / 90.0) as i64;
    quarter.rem_euclid(2) == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The picture between two frames: a particle in both moves linearly, one absorbed
    /// before the second stays where it was, the surface densities and potentials are
    /// interpolated; before the first frame the first, after the last the last.
    #[test]
    fn frames_interpolate() {
        let frame = |t: f64, x: Vec<[f32; 2]>, id: Vec<u32>, sigma: f32| TubeFrame {
            t,
            q: vec![-1.0; x.len()],
            kind: vec![KIND_NEGATIVE; x.len()],
            x,
            id,
            sigma: vec![sigma; 3],
            potentials: vec![sigma; 2],
            currents: Vec::new(),
            heat: 0.0,
            paths: Vec::new(),
        };
        let run = TubeRun {
            frames: vec![
                frame(1.0, vec![[0.0, 0.0], [5.0, 5.0]], vec![7, 8], 1.0),
                // Particle 8 absorbed, 9 emitted, 7 moved and now listed second.
                frame(2.0, vec![[9.0, 9.0], [1.0, 2.0]], vec![9, 7], 3.0),
            ],
            ..TubeRun::default()
        };
        let mid = run.frame_between(1.25).unwrap();
        assert!((mid.t - 1.25).abs() < 1e-12);
        assert_eq!(mid.x, vec![[0.25, 0.5], [5.0, 5.0]]);
        assert_eq!(mid.id, vec![7, 8]);
        assert_eq!(mid.sigma, vec![1.5; 3]);
        assert_eq!(mid.potentials, vec![1.5; 2]);
        assert_eq!(run.frame_between(0.5).unwrap().x, run.frames[0].x);
        assert_eq!(run.frame_between(3.0).unwrap().x, run.frames[1].x);
    }

    /// A carrier recorded at every step (`TubeFrame::paths`) moves along its steps between
    /// frames, not straight across: turning a corner between frames 0.3 apart (steps of
    /// 0.1: from (0, 0) to (1, 0), (1, 1) and (0, 1)), at t = 1.15 it is halfway up the
    /// corner's side, (1, 0.5) (straight across it would be at (0, 0.5)); its trail from
    /// t = 1.0 to 1.25 goes round the corner; an unrecorded carrier moves straight.
    #[test]
    fn sampled_carriers_follow_their_steps() {
        let frame = |t: f64, x: Vec<[f32; 2]>, paths: Vec<(u32, Vec<[f32; 2]>)>| TubeFrame {
            t,
            q: vec![-1.0; x.len()],
            kind: vec![KIND_NEGATIVE; x.len()],
            id: vec![7, 8],
            x,
            paths,
            ..TubeFrame::default()
        };
        let run = TubeRun {
            frames: vec![
                frame(1.0, vec![[0.0, 0.0], [0.0, 0.0]], Vec::new()),
                frame(
                    1.3,
                    vec![[0.0, 1.0], [0.0, 1.0]],
                    vec![(7, vec![[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])],
                ),
            ],
            dt: 0.1,
            fine_every: 7,
            ..TubeRun::default()
        };
        let mid = run.frame_between(1.15).unwrap();
        let close =
            |a: [f32; 2], b: [f32; 2]| (a[0] - b[0]).abs() < 1e-5 && (a[1] - b[1]).abs() < 1e-5;
        assert!(close(mid.x[0], [1.0, 0.5]), "sampled: {:?}", mid.x[0]);
        assert!(close(mid.x[1], [0.0, 0.5]), "straight: {:?}", mid.x[1]);
        let trail = run.trail(7, 1.0, 1.25);
        let expected = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.5, 1.0]];
        assert_eq!(trail.len(), expected.len(), "{trail:?}");
        for (a, b) in trail.iter().zip(expected) {
            assert!(close(*a, b), "{trail:?}");
        }
    }

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
            reference: None,
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
