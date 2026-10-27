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

mod build;
mod frames;
mod geometry;
mod issues;
mod run;
mod sim;
mod verdict;
mod verify;
use frames::lerp_cells;
pub use frames::metal_steps;

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
    /// of its surface (`physics::tube::Emitter::saturation`; 0: a cold electrode that emits
    /// nothing, the cathode a level with isolated metals only needs). None: space-charge limited
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
    /// The carriers' own magnetic field (the Darwin limit, PHYSICS.md §2.12) at this speed of
    /// light: their currents' field acting back on them (inductance, skin effect, transformers),
    /// its effects going as `1/c²`. None: none, the electrostatic model. Its step is first
    /// order: a preview step of at most `SELF_FIELD_MAX_STEP`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub speed_of_light: Option<f64>,
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
    /// Its function in the device (channel, load, Hall plate, …): shown with its number
    /// wherever it is named. Text only.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
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

/// The self-field's display samples apart (cells): two a cell, the frames' `self_b`.
pub const SELF_B_SPACING: f64 = 0.5;

/// A grid of display samples: `nx × ny` points `spacing` apart from `origin`, row by row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayGrid {
    pub origin: (f64, f64),
    pub spacing: f64,
    pub nx: usize,
    pub ny: usize,
}

impl DisplayGrid {
    /// The board `width × height` from the origin, `spacing` apart, its edges included.
    #[must_use]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub fn board(width: f64, height: f64, spacing: f64) -> Self {
        Self {
            origin: (0.0, 0.0),
            spacing,
            nx: libm::floor(width / spacing) as usize + 1,
            ny: libm::floor(height / spacing) as usize + 1,
        }
    }

    /// The points, row by row.
    #[allow(clippy::cast_precision_loss)]
    pub fn points(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        (0..self.ny).flat_map(move |j| {
            (0..self.nx).map(move |i| {
                (
                    self.origin.0 + i as f64 * self.spacing,
                    self.origin.1 + j as f64 * self.spacing,
                )
            })
        })
    }
}

/// The self-field's grid spacing (`physics::darwin2d`; V-D4–V-D6 were validated at it: a
/// metal's skin depth of ~2 is eight cells). The verification's refinement keeps it.
pub const SELF_FIELD_CELL: f64 = 0.25;
/// The self-field's grid beyond the board on every side: a loop's field stays inside it
/// (z-invariance), the Coulomb gauge's vector potential has a long tail.
pub const SELF_FIELD_MARGIN: f64 = 4.0;
/// The largest preview step with the self-field: its splitting is first order, and the
/// coupling's 5 % target holds up to dt ≈ 0.0125 (Y24: 4.0 % at 0.01, 7.8 % at 0.02;
/// PHYSICS.md §2.12). (0.015 until 2026-10-27, from Y24's 3.1 % at 0.01 before the
/// impulses' currents were fixed.)
pub const SELF_FIELD_MAX_STEP: f64 = 0.012;
/// The speeds of light a level may set (the editor's range, `check_editable`).
pub const SPEED_OF_LIGHT_RANGE: std::ops::RangeInclusive<f64> = 0.5..=1.0e4;
/// The self-field's largest grid (nodes per side): its solve each step is the cost.
pub const SELF_FIELD_MAX_NODES: usize = 511;

/// The self-field's grid for a board `width × height` from the origin at the speed of light
/// `c`: `SELF_FIELD_CELL` apart, the board and `SELF_FIELD_MARGIN` round it, its sides rounded
/// up to the sine transforms' powers of two, centred on the board.
#[must_use]
pub fn self_field_grid(width: f64, height: f64, c: f64) -> physics::darwin2d::Darwin2d {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let cells = |len: f64| {
        (libm::ceil((len + 2.0 * SELF_FIELD_MARGIN) / SELF_FIELD_CELL) as usize).next_power_of_two()
    };
    let (cx, cy) = (cells(width), cells(height));
    #[allow(clippy::cast_precision_loss)]
    let origin = (
        0.5 * width - 0.5 * cx as f64 * SELF_FIELD_CELL,
        0.5 * height - 0.5 * cy as f64 * SELF_FIELD_CELL,
    );
    physics::darwin2d::Darwin2d::new(cx - 1, cy - 1, SELF_FIELD_CELL, origin, c)
}

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
    /// With the carriers' own magnetic field (`speed_of_light`): its `B_z` on the run's
    /// display grid (`TubeRun::self_b_grid`), row by row; empty without it.
    pub self_b: Vec<f32>,
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
    /// Per metal (`TubeRun::metals`), its accounts by cell so far (`MetalCells`). Empty
    /// without metals.
    pub cells: Vec<MetalCells>,
}

/// A metal's accounts by cell at a frame (`physics::tube::CellAccounts`, `MetalView::cells`,
/// row by row, so far): the charge carried (x, y) and the heat goal's work (the field's on
/// the carriers; in an isolated metal the induced field's: `MetalView::isolated`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetalCells {
    pub flux: Vec<[f32; 2]>,
    pub work: Vec<f32>,
}

/// A metal's mean current density and power per area by cell over a span of time
/// (`TubeRun::metal_rates`, row by row).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetalRates {
    pub j: Vec<[f64; 2]>,
    pub power: Vec<f64>,
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
    /// Its accounts' cells (`physics::lattice::Lattice::account_size`), its cuts (the
    /// lattice's walls) and whether it is isolated (no contacts: its heat is the induced
    /// field's work, `TubeSim::metal_work`), for the maps inside it.
    pub cells: (usize, usize),
    pub walls: Vec<[[f64; 2]; 2]>,
    pub isolated: bool,
    /// Its contacts' faces along its sides (each from one end to the other): where the
    /// display takes carriers in and gives them off between frames (`TubeRun::frame_between`).
    pub contacts: Vec<[[f64; 2]; 2]>,
}

/// What the metal maps show at a time (`TubeRun::metal_maps`): the span averaged, whether
/// the current is the lock-in's (a swinging field), and per metal the values by cell with
/// their standard errors (`MetalRates` both; the errors' `j` per component).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetalMaps {
    pub from: f64,
    pub to: f64,
    pub lock_in: bool,
    pub values: Vec<MetalRates>,
    pub errors: Vec<MetalRates>,
}

/// How many sub-windows the metal maps' standard errors come from (the batch means').
const METAL_MAP_BATCHES: usize = 5;

/// The span of time the metal maps average over by default (`TubeRun::metal_rates`): a
/// cell's ~50 macroparticles then carry a thermal noise of about a fifth of an ordinary
/// drift (the thermal speed is five times it), a third of that after the 3 × 3 average.
pub const METAL_MAP_SPAN: f64 = 4.0;

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
    /// With the carriers' own magnetic field: the grid its frames' `self_b` are sampled on
    /// (the board, `SELF_B_SPACING` apart).
    pub self_b_grid: Option<DisplayGrid>,
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
    /// The contacts' potential steps so far, `(t, Σ)` over the electrodes, at each frame
    /// with metals, when booked (`TubeSim::book_potential_steps`; `TubeState::
    /// potential_steps`: the circuit's energy the field's work leaves out).
    pub potential_steps: Vec<(f64, f64)>,
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

/// A tube run in progress: advanced in pieces (`advance_to`), its frames taken as they
/// are made (`take_chunk`), so that a long display run can be streamed.
pub struct TubeSim {
    /// With the self-field: where its frames sample it (`TubeFrame::self_b`).
    self_b_grid: Option<DisplayGrid>,
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
    /// Each lattice's metal (one each: a cut is a wall inside it).
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

impl Level {
    /// Whether this is a tube level.
    pub fn is_tube(&self) -> bool {
        self.tube.is_some()
    }
}

#[cfg(test)]
mod tests;
