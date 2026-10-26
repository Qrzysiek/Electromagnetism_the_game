//! Vacuum tubes in the z-invariant world (docs/TUBES.md; PHYSICS.md §2.11): electrons
//! emitted from a cathode under the space-charge limit, flying to the other electrodes
//! in the field of the electrodes, of each other, and of the charge they induce on the
//! electrodes.
//!
//! - **Particles** are macroparticles: lines of charge along z, each `weight` charge per
//!   unit length, with the carriers' charge-to-mass ratio, so each follows an electron's
//!   trajectory exactly (only the graininess of the space charge depends on the weight:
//!   a stated approximation, measured by halving it). Newtonian (tube energies: v ~ 0.03
//!   c; the magnetic and retarded corrections, ~v²/c², are left out and stated).
//! - **Fields** at a particle: the electrodes' at their potentials (the static
//!   solution), the induced surface charge of all particles (solved each step from the
//!   factorization, all electrodes held at their potentials), and every other particle
//!   (2D Coulomb, `2λ r̂/r`): each a uniformly charged rod of radius `softening`, much
//!   below the spacing, so its field is exactly Coulomb outside the rod (`2λ r⃗/h²` inside).
//! - **Emission** (space-charge limited, by Gauss's law): each step, each segment of
//!   the cathode whose surface density has the carriers' sign emits that charge (its
//!   density times its length) as macroparticles, released at rest `EMISSION_OFFSET` out at
//!   points along it from a deterministic low-discrepancy sequence. Emitting the charge
//!   the surface holds drives the cathode's normal field to zero there, the limit.
//! - **Absorption**: a particle that crosses an electrode's boundary (any change of
//!   which side of it the particle is on) is collected by that electrode, which records
//!   the charge; the current is the collected charge per time. A carrier that comes back
//!   to its emitter is absorbed there on reaching the layer it was released from
//!   (`RETURN_LEAVE`, `RETURN_ABSORB`): released at rest a distance δ out, it lacks the
//!   energy of that first δ, so it could never come closer again and would bounce there.
//! - **Temperature-limited emission** (optional per emitter): at most a saturation
//!   current density, whatever the field (Richardson's limit of a hot cathode); the
//!   emission is the lesser of the two limits.
//! - **Species**: several emitters (an electron cathode, an ion-emitting anode: the
//!   bipolar diode), and projectiles, single line charges launched at given times that
//!   interact with everything and are absorbed where they hit.
//! - **Circuits** (per unit length along z, consistent with the symmetry): an electrode may
//!   be a circuit node, connected through a resistance to a source and through a
//!   capacitance to ground (resistance 0: an ideal source the electrode follows; infinite:
//!   floating). The node's charge (its conductor's, `Σ_j C_ij V_j` plus the particles'
//!   induced charge, and its capacitor's) changes by the resistor's current and by the
//!   charge of the particles it absorbs or emits; each step its potential is solved
//!   implicitly (backward Euler in the resistor) from that charge, the other electrodes'
//!   potentials and the induced charge, so the particles and the circuit act on each other
//!   exactly (the induced charge carries the Ramo currents). Floating electrodes (a
//!   fixed charge) are nodes without a circuit, so their charge is conserved as particles
//!   fly by.
//! - **Push**: Boris with a fixed step (a uniform B along z rotates the velocity; with
//!   B = 0 it is the kick-drift leapfrog), symplectic; launches and emission add
//!   particles at the start of a step.
//! - **A changing field** (`SolenoidWave`, optional): B_z follows a waveform, and by
//!   Faraday's law (`∇ × E = −∂B/∂t`; these units write the Lorentz force `q(E + v × B)`)
//!   it induces `E = −½ Ḃ ẑ × (x − axis)`, circulating about the solenoid's axis, which
//!   pushes every particle (a metal's eddy currents). The currents' own magnetic field is
//!   left out, as everywhere here: the eddy currents do not screen the change (valid where
//!   the skin depth is far larger than the metal).

use glam::DVec3;

use crate::bem::{Bias, Lu};
use crate::circuit::Waveform;
use rayon::prelude::*;

use crate::lattice::Lattice;
use crate::tree2d::{Source, Tree};
use crate::zinv::{
    Electrodes, PAR_CHUNK, Section, line_charge, segment_field_fast, segment_integrals,
};

mod circuits;
mod darwin;
mod emission;
mod fields;
mod metals;

/// From this many particles in flight the step sums the fields and the induced charge's
/// potentials by trees (`tree2d`) instead of directly: below it the direct sums are faster
/// (and the tree's segment field more exact than `segment_field_fast`'s two terms, by up to
/// 4.9e-5 of a segment's monopole, so results change slightly across it).
pub const TREE_MIN: usize = 2000;

/// Distance from the cathode at which macroparticles are released, in segment lengths:
/// at the mesh's resolution, where a macroparticle's own image (the field w/h of the
/// charge it induces on the cathode, an artifact of its finite weight: an infinitesimal
/// charge has none) is no stronger than its neighbours' fields. Released at 1e-3 of a
/// segment the self-image pulled them back into the cathode and the coaxial diode's
/// current came out 43 % (and 12.5 % refined) below Langmuir–Blodgett (test V4). Released
/// at rest half a segment out, the particles skip the first δ of the space-charge layer:
/// the current is high by an error of first order in the resolution (15 %, then 7.8 %
/// with everything halved; Richardson-extrapolated 0.5 %). Giving them the layer's
/// analytic speed at δ instead emptied the space charge near the cathode, which then
/// emitted more (+111 %, +67 %): rejected.
pub const EMISSION_OFFSET: f64 = 0.5;

/// A carrier coming back to its emitter (PHYSICS.md §2.11): once it has been farther
/// than `RETURN_LEAVE` release distances δ (`EMISSION_OFFSET` segments) from the emitter's
/// surface, it is absorbed by the emitter as soon as it is back within `RETURN_ABSORB` δ.
/// Released at rest at δ, a macroparticle lacks the energy a real carrier gains over that
/// first δ, so in a static field it returns exactly to δ with zero speed and turns again
/// (a real carrier would reach the surface): in a triode with a negative grid the
/// returning electrons bounced in front of the cathode forever, drifted along it and
/// left the board (more than the anode collected). Absorbing on the release layer's
/// return (between δ and 2δ, so not the exact turning point, where the outcome would
/// hinge on rounding) is the real carrier's fate, with an error of the release's order:
/// the last 0.5 δ of the way back (test V12, the orbits' time to 1e-3).
pub const RETURN_LEAVE: f64 = 2.0;
pub const RETURN_ABSORB: f64 = 1.5;

/// A space-charge-limited emitter: an electrode (or one coated face of it) releasing
/// carriers of one kind.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Emitter {
    pub electrode: usize,
    /// Charge-to-mass ratio of its carriers (electrons negative, ions positive).
    pub charge_per_mass: f64,
    /// Charge per unit length of a macroparticle (its magnitude).
    pub weight: f64,
    /// The emitting face (a coating): only segments whose outward normal is within 60°
    /// of this unit direction emit. None: every face.
    pub emit_toward: Option<DVec3>,
    /// Temperature-limited emission: at most this current per unit area (charge per
    /// unit length of the surface, per unit z-length, per time; the Richardson current
    /// of the cathode's temperature). Below it the space charge limits the emission,
    /// above it the temperature (the Langmuir knee). None: space-charge limited only.
    pub saturation: Option<f64>,
    /// Released at rest (0), or at this thermal speed `√(kT/m)` as the carriers of a
    /// gas at that temperature cross a wall into it (a metal's contact at the lattice's
    /// temperature): the normal speed Rayleigh-distributed, the tangential normal
    /// (`lattice::flux_velocity`).
    pub thermal_speed: f64,
    /// A metal's ohmic contact: in place of the space-charge limit, each emitting segment
    /// keeps the layer in front of it neutral (the boundary condition of device Monte
    /// Carlo codes), giving off what its carriers lack there.
    pub ohmic: Option<OhmicContact>,
}

/// An ohmic contact's rule (`Emitter::ohmic`): the layer `depth` deep in front of each
/// emitting segment holds `density` (the carriers' charge per area, its magnitude, the
/// metal's background's opposite); what it lacks is given off from the face as a gas at
/// the lattice's temperature (`thermal_speed`) crosses a wall (`lattice::flux_velocity`),
/// so that it gives back, on average, the energy the carriers it takes in bring. The
/// space-charge limit instead (the field at its face held at zero, as a cathode's) piled
/// carriers in front of an injecting contact: a layer with its own voltage; carriers
/// placed in the layer with a Maxwellian's velocities (kT each, where those reaching the
/// face bring 1.5 kT) made every contact a heat sink (PHYSICS.md §2.11).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OhmicContact {
    pub density: f64,
    pub depth: f64,
}

/// A projectile: one line charge (not a macroparticle) launched at a time, with its own
/// charge and mass per unit length; it interacts with everything and is absorbed by the
/// electrode it hits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Projectile {
    pub x0: DVec3,
    pub v0: DVec3,
    pub launch: f64,
    pub charge: f64,
    pub mass: f64,
}

/// An electrode's external circuit, per unit length along z: a branch of a source, a
/// resistance and an inductance in series, from ground (or from another electrode, `to`)
/// into the electrode, and a capacitance from the electrode to ground.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circuit {
    pub electrode: usize,
    pub source: Waveform,
    /// Resistance per unit length (potential × time / charge per unit length): 0 without
    /// an inductance is an ideal source (the electrode follows it), infinite no source
    /// (floating).
    pub resistance: f64,
    /// External capacitance per unit length to ground (charge per unit length / potential).
    pub capacitance: f64,
    /// Inductance in series with the resistance (potential × time² / charge per unit
    /// length): 0 none. A branch with one is integrated by the trapezoidal rule (an LC
    /// tank keeps its energy); one without, by backward Euler.
    pub inductance: f64,
    /// The branch's other end: this electrode instead of ground (a tank between two
    /// electrodes). Its source's voltage adds to that electrode's potential. With no
    /// resistance and no inductance the branch is a strap: a wire (its source unused), the
    /// electrode at that one's potential (a magnetron's strapped vanes: one circuit node).
    /// A strap's other end is not itself strapped.
    pub to: Option<usize>,
}

impl Circuit {
    /// An ideal source: the electrode follows it.
    fn is_ideal(&self) -> bool {
        self.resistance == 0.0 && self.inductance == 0.0 && self.to.is_none()
    }

    /// A strap: a wire to another electrode, whose potential this one follows.
    fn is_strap(&self) -> bool {
        self.resistance == 0.0 && self.inductance == 0.0 && self.to.is_some()
    }

    /// A branch with an inductance (and a finite resistance: an infinite one carries no
    /// current, whatever its inductance).
    fn is_inductive(&self) -> bool {
        self.inductance > 0.0 && self.resistance.is_finite()
    }
}

/// A circuit node: an electrode whose potential follows from its charge.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Node {
    electrode: usize,
    /// Step over resistance of its own resistive branch (0 without one, or with an
    /// inductive branch, which has its own unknown).
    conductance_dt: f64,
    capacitance: f64,
    source: Waveform,
}

/// Where a projectile ended.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Fate {
    /// Hit electrode `electrode` at time `t`, point `x`.
    Absorbed { electrode: usize, t: f64, x: DVec3 },
    /// Left the arena.
    Lost { t: f64, x: DVec3 },
}

/// Most macroparticles in flight: beyond it emission pauses and the state records it
/// (`TubeState::capped`; the level's verdict is then uncertain). The game's runs hold a
/// few hundred; a degenerate weight once asked for ~1e300 (a cathode with no charge of
/// its own emitting the charge a passing ion induced, at a weight of its own charge).
pub const MAX_PARTICLES: usize = 20_000;

/// The Darwin field's solve: its relative residual (`darwin2d::Darwin2d::solve`).
pub const DARWIN_TOLERANCE: f64 = 1e-10;

/// Tag of projectile `j` (`TubeState::tag`); emitted carriers carry their emitter's
/// index.
pub const PROJECTILE_TAG: u32 = 1 << 16;

/// A uniform field along z that changes in time (a long solenoid around the board on a
/// time-dependent supply): `B_z(t)` follows `waveform`, and induces the electric field
/// `E = −½ Ḃ ẑ × (x − axis)` (Faraday's law, `∇ × E = −Ḃ ẑ`), circulating about the
/// solenoid's axis. Inside an isolated conductor only its curl matters (the charges it
/// moves to the surface cancel any uniform part); a circuit through the field picks up the
/// EMF of the flux it encloses, which this gauge assigns about `axis`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolenoidWave {
    pub waveform: Waveform,
    pub axis: DVec3,
}

impl SolenoidWave {
    /// `B_z` at `t`.
    #[must_use]
    pub fn b(&self, t: f64) -> f64 {
        self.waveform.value(t)
    }

    /// The largest |B| it reaches.
    #[must_use]
    pub fn peak(&self) -> f64 {
        match self.waveform {
            Waveform::Dc(v) => v.abs(),
            Waveform::Sine {
                offset, amplitude, ..
            } => offset.abs() + amplitude.abs(),
            Waveform::Pulse { low, high, .. } => low.abs().max(high.abs()),
        }
    }

    /// The induced electric field at `x` at `t`.
    #[must_use]
    pub fn induced(&self, t: f64, x: DVec3) -> DVec3 {
        let half = 0.5 * self.waveform.slope(t);
        DVec3::new(half * (x.y - self.axis.y), -half * (x.x - self.axis.x), 0.0)
    }
}

/// The tube's setup.
#[derive(Debug)]
pub struct Tube {
    /// The electrodes at their potentials (the static solution).
    pub electrodes: Electrodes,
    pub emitters: Vec<Emitter>,
    pub projectiles: Vec<Projectile>,
    /// The macroparticles' radius h: each is a uniformly charged rod, its field `2λ r⃗/h²`
    /// inside and exactly `2λ r⃗/r²` outside (a compact softening: beyond h the mutual force
    /// is Coulomb's, so far groups of particles may be summed as multipoles exactly).
    pub softening: f64,
    /// Time step.
    pub dt: f64,
    /// The arena `(min, max)` corners: a particle that leaves it is lost (it flies off
    /// to infinity in the open tube; in the z-invariant world nothing brings it back
    /// unless an electrode encloses the arena, and then it never leaves). None: no bound.
    pub arena: Option<(DVec3, DVec3)>,
    /// A uniform magnetic field along z (a long solenoid around the tube: z-invariant).
    pub b_z: f64,
    /// A field that changes in time, with its induced electric field, in place of `b_z`.
    pub b_wave: Option<SolenoidWave>,
    /// External circuits of electrodes.
    pub circuits: Vec<Circuit>,
    /// Conductors' interiors (`lattice.rs`): their carriers fill them at the start (tagged
    /// `emitters.len() + k`), every carrier inside one scatters off its ion cores and
    /// walls, and their backgrounds add to the field.
    pub lattices: Vec<Lattice>,
    /// The charges' own magnetic field in the Darwin limit (`darwin2d`, docs/DARWIN.md,
    /// PHYSICS.md §2.12): its grid and speed of light. None: none (the electrostatic model
    /// alone, bit for bit as before).
    pub darwin: Option<crate::darwin2d::Darwin2d>,
}

/// A lattice's accounts by cell (`Lattice::account_cell`, row by row), so far: bookkeeping
/// for the display's maps of the current and the heat inside a metal. Nothing here acts on
/// the motion (the golden hashes stayed); the cells add up to `TubeState::field_work` and
/// `induced_work` (test `d9_cell_accounts_add_up`).
#[derive(Clone, Debug, Default)]
pub struct CellAccounts {
    /// The charge the carriers carry through the cell, `Σ q Δx` of the pieces of their steps'
    /// chords inside it (`Lattice::book_flux`; collisions and walls included in the steps):
    /// the current density's time integral over the cell.
    pub flux: Vec<[f64; 2]>,
    /// The field's work on the carriers in the cell (their share of `TubeState::field_work`),
    /// and the induced field's (of `induced_work`).
    pub work: Vec<f64>,
    pub induced: Vec<f64>,
}

/// The particles and the bookkeeping of a running tube.
#[derive(Clone, Debug, Default)]
pub struct TubeState {
    pub t: f64,
    pub x: Vec<DVec3>,
    pub v: Vec<DVec3>,
    pub charge: Vec<f64>,
    /// Charge-to-mass ratio of each particle.
    pub qm: Vec<f64>,
    /// Emitter index, or `PROJECTILE_TAG + j` for projectile j.
    pub tag: Vec<u32>,
    /// Each particle's identity, kept while it flies (the display interpolates between
    /// steps by it).
    pub id: Vec<u32>,
    /// The next particle's identity.
    next_id: u32,
    /// For each particle, per electrode, whether it is inside the electrode's section.
    side: Vec<Vec<bool>>,
    /// For each particle, whether it has left its emitter's release layer (beyond
    /// `RETURN_LEAVE` release distances; projectiles never).
    left: Vec<bool>,
    /// Per emitter, per segment of its electrode, charge owed but not yet emitted (less
    /// than a macroparticle).
    pending: Vec<Vec<f64>>,
    /// Counter of the emission sequence.
    emitted: u64,
    /// Charge per unit length collected by each electrode so far, per emitter.
    pub collected: Vec<Vec<f64>>,
    /// Charge per unit length each emitter has given off so far (signed): with
    /// `collected`, a contact's current is what it takes in less what it gives off.
    pub released: Vec<f64>,
    /// Charge per unit length that left the arena.
    pub lost: f64,
    /// Which projectiles have been launched, and where each ended.
    launched: Vec<bool>,
    pub fates: Vec<Option<Fate>>,
    /// Emission paused at `MAX_PARTICLES` at least once.
    pub capped: bool,
    /// The energy per unit length the lattices have taken in collisions (the Joule heat),
    /// and the collisions.
    pub joule: f64,
    pub collisions: u64,
    /// Per lattice: the work the field has done on the carriers inside it, `Σ (q dt/2)
    /// E·(v_old + v_new)` (the push's exact energy account: the magnetic rotation does no
    /// work). In a steady state the Joule heat, the lattice's take, with far less noise:
    /// each collision exchanges about kT at random, while the field's work fluctuates by
    /// `qEℓ/kT` times as much.
    pub field_work: Vec<f64>,
    /// Per lattice: the part of `field_work` the induced field (`SolenoidWave::induced`)
    /// does. In an isolated conductor it is the Joule heat: there the electrostatic field's
    /// work averages to nothing (it only stores and returns energy), while the leapfrog's
    /// finite step makes it heat the carriers through their own graininess (∝ dt²: on an
    /// isolated plate 12 × 12, 2.57 per time at the preview's resolution and 0.60 at the
    /// verification's, against an eddy-current heat of 1.24).
    pub induced_work: Vec<f64>,
    /// Per lattice: the current and those works by cell (`CellAccounts`).
    pub cells: Vec<CellAccounts>,
    /// Per electrode: the kinetic energy of the carriers it took in less that of those it
    /// gave off (an ohmic contact's), so far. With the lattices' heat, where the field's work
    /// on the carriers goes (test `d10_carriers_energy_account`).
    pub contact_energy: Vec<f64>,
    /// Per electrode, with lattices and `book_potential_steps`: the potential steps of the
    /// carriers it gave off and took in, `q (V − φ)` at each release and `q (φ − V)` at each
    /// absorption (φ the potential at the carrier of everything else, V its electrode's): the
    /// part of the circuit's energy that the field's work on the carriers leaves out (test
    /// `d10_carriers_energy_account`).
    pub potential_steps: Vec<f64>,
    /// Whether to book `potential_steps` (a diagnostic, off by default: each release and
    /// absorption then sums the potential over every particle, ~40 % more time in a metal
    /// level).
    pub book_potential_steps: bool,
    /// The density the lattices' backgrounds induce on the segments (constant).
    background: Vec<f64>,
    /// Every electrode's potential now.
    pub potentials: Vec<f64>,
    /// Per circuit node (`Tube::nodes`): its charge (its conductor's and its capacitor's).
    node_charge: Vec<f64>,
    /// Per circuit (in `Tube::circuits` order): the current from its source into the
    /// electrode over the last step (0 for an ideal source, whose current is not
    /// tracked).
    pub circuit_current: Vec<f64>,
    /// Per circuit: its inductance's current now (0 without one).
    pub inductor_current: Vec<f64>,
    /// The particles' induced surface density now (all electrodes grounded), kept from
    /// the end of the last step: the next step's emission and the display reuse it.
    induced: Vec<f64>,
    /// With the carriers' own magnetic field (`Tube::darwin`): each particle's canonical
    /// velocity `u = v + (q/m) A_T(x)` (its canonical momentum per mass), the field's solver
    /// with its last solution, and that solution's `A_T` on the grid's edges. Empty and None
    /// without it.
    pub u: Vec<DVec3>,
    pub darwin: Option<crate::darwin2d::Darwin2d>,
    pub self_a: Option<crate::darwin2d::Edges>,
}

impl TubeState {
    /// Where projectile `j` is (None before its launch and after its end).
    pub fn projectile(&self, j: usize) -> Option<DVec3> {
        let tag = PROJECTILE_TAG + u32::try_from(j).ok()?;
        self.tag.iter().position(|&t| t == tag).map(|k| self.x[k])
    }
}

/// Where the chord from `a` to `b` crosses the wall `w0`–`w1`, as the fraction of the chord
/// (None: it does not). A side counts as crossed from the wall itself to either side, so a
/// carrier on a wall cannot slip through it.
fn crossing(a: DVec3, b: DVec3, w0: DVec3, w1: DVec3) -> Option<f64> {
    let e = w1 - w0;
    let side = |p: DVec3| (p.x - w0.x) * e.y - (p.y - w0.y) * e.x;
    let (sa, sb) = (side(a), side(b));
    if (sa > 0.0 && sb > 0.0) || (sa < 0.0 && sb < 0.0) || (sa == 0.0 && sb == 0.0) {
        return None;
    }
    let t = sa / (sa - sb);
    let p = a + (b - a) * t;
    let len2 = e.x * e.x + e.y * e.y;
    let u = ((p.x - w0.x) * e.x + (p.y - w0.y) * e.y) / len2;
    (0.0..=1.0).contains(&u).then_some(t)
}

/// Whether `x` lies inside the section.
fn inside(section: &Section, x: DVec3) -> bool {
    match *section {
        Section::Rect {
            center,
            angle,
            half_length,
            half_thickness,
        } => {
            let (s, c) = libm::sincos(angle);
            let d = x - center;
            (d.x * c + d.y * s).abs() < half_length && (-d.x * s + d.y * c).abs() < half_thickness
        }
        Section::Circle { center, radius } => {
            let d = x - center;
            d.x * d.x + d.y * d.y < radius * radius
        }
        Section::Ring {
            center,
            radius,
            thickness,
        } => {
            let d = x - center;
            let r = libm::sqrt(d.x * d.x + d.y * d.y);
            (r - radius).abs() < 0.5 * thickness
        }
        Section::Frame {
            center,
            half_width,
            half_height,
            ..
        } => {
            let d = x - center;
            !(d.x.abs() < half_width && d.y.abs() < half_height)
        }
    }
}

/// The distance from `x` to the section (0 inside).
fn distance_outside(section: &Section, x: DVec3) -> f64 {
    match *section {
        Section::Rect {
            center,
            angle,
            half_length,
            half_thickness,
        } => {
            let (s, c) = libm::sincos(angle);
            let d = x - center;
            let u = ((d.x * c + d.y * s).abs() - half_length).max(0.0);
            let w = ((-d.x * s + d.y * c).abs() - half_thickness).max(0.0);
            libm::sqrt(u * u + w * w)
        }
        Section::Circle { center, radius } => {
            let d = x - center;
            (libm::sqrt(d.x * d.x + d.y * d.y) - radius).max(0.0)
        }
        Section::Ring {
            center,
            radius,
            thickness,
        } => {
            let d = x - center;
            ((libm::sqrt(d.x * d.x + d.y * d.y) - radius).abs() - 0.5 * thickness).max(0.0)
        }
        Section::Frame {
            center,
            half_width,
            half_height,
            ..
        } => {
            let d = x - center;
            (half_width - d.x.abs())
                .min(half_height - d.y.abs())
                .max(0.0)
        }
    }
}

impl Tube {
    /// A tube of `electrodes`, its particles softened to rods of radius `softening`, stepped
    /// by `dt`, with nothing else: no emitters, projectiles, arena, magnetic field, circuits,
    /// metals or self-field. Build the rest with the struct update syntax,
    /// `Tube { lattices, ..Tube::new(electrodes, softening, dt) }`, so that a field added to
    /// `Tube` is given its default here alone.
    #[must_use]
    pub fn new(electrodes: Electrodes, softening: f64, dt: f64) -> Self {
        Self {
            electrodes,
            emitters: Vec::new(),
            projectiles: Vec::new(),
            softening,
            dt,
            arena: None,
            b_z: 0.0,
            b_wave: None,
            circuits: Vec::new(),
            lattices: Vec::new(),
            darwin: None,
        }
    }

    pub fn state(&self) -> TubeState {
        let n_el = self.electrodes.electrodes.len();
        // The potentials at t = 0: the electrodes' own (their biases, floating ones as
        // solved), ideal sources at their value.
        let mut potentials = self.electrodes.potentials.clone();
        for c in &self.circuits {
            if c.is_ideal() {
                potentials[c.electrode] = c.source.value(0.0);
            }
        }
        self.follow_straps(&mut potentials);
        let cap = &self.electrodes.capacitance;
        let nodes = self.nodes();
        let node_charge = (0..nodes.len())
            .map(|r| {
                self.members(&nodes, r)
                    .map(|i| (0..n_el).map(|j| cap[i][j] * potentials[j]).sum::<f64>())
                    .sum::<f64>()
                    + nodes[r].capacitance * potentials[nodes[r].electrode]
            })
            .collect();
        let mut s = TubeState {
            // A quiet start: each segment's account begins at a different fraction of a
            // macroparticle (the golden-ratio sequence), so that segments emitting at the
            // same rate do not emit together in volleys.
            pending: self
                .emitters
                .iter()
                .map(|e| {
                    (0..self.segments_of(e.electrode).count())
                        .map(|k| {
                            #[allow(clippy::cast_precision_loss)]
                            let f = ((k + 1) as f64 * 0.618_033_988_749_894_9).fract();
                            f * e.weight
                        })
                        .collect()
                })
                .collect(),
            collected: vec![vec![0.0; self.emitters.len() + self.lattices.len()]; n_el],
            contact_energy: vec![0.0; n_el],
            potential_steps: vec![0.0; n_el],
            released: vec![0.0; self.emitters.len()],
            field_work: vec![0.0; self.lattices.len()],
            induced_work: vec![0.0; self.lattices.len()],
            cells: self
                .lattices
                .iter()
                .map(|l| {
                    let (nx, ny) = l.account_size();
                    CellAccounts {
                        flux: vec![[0.0; 2]; nx * ny],
                        work: vec![0.0; nx * ny],
                        induced: vec![0.0; nx * ny],
                    }
                })
                .collect(),
            launched: vec![false; self.projectiles.len()],
            fates: vec![None; self.projectiles.len()],
            potentials,
            node_charge,
            circuit_current: vec![0.0; self.circuits.len()],
            inductor_current: vec![0.0; self.circuits.len()],
            induced: vec![0.0; self.electrodes.segments.len()],
            darwin: self.darwin.clone(),
            self_a: self.darwin.as_ref().map(crate::darwin2d::Edges::zeros),
            ..TubeState::default()
        };
        // The conductors: their carriers, and the charge their backgrounds induce.
        if !self.lattices.is_empty() {
            let n_em = u32::try_from(self.emitters.len()).unwrap_or(0);
            for (k, l) in self.lattices.iter().enumerate() {
                let q = l.carrier_charge();
                let tag = n_em + u32::try_from(k).unwrap_or(0);
                for (x, v) in l.fill() {
                    self.add(&mut s, x, v, q, l.charge_per_mass, tag);
                }
            }
            s.background = self
                .electrodes
                .induced_by(|y| self.lattices.iter().map(|l| l.background_at(y).0).sum());
            s.induced = self.induced_density(&s);
            for (i, b) in s.induced.iter_mut().zip(&s.background) {
                *i += b;
            }
        } else {
            s.background = vec![0.0; self.electrodes.segments.len()];
        }
        s
    }

    /// One step: launches, emission, push (Boris: with B = 0 the kick-drift
    /// leapfrog), absorption.
    pub fn step(&self, s: &mut TubeState) {
        // The induced density is kept from the last step; the particles added during this
        // one (launched, emitted) add theirs.
        let mut induced = std::mem::take(&mut s.induced);
        let before = s.x.len();
        for (j, p) in self.projectiles.iter().enumerate() {
            if !s.launched[j] && s.t >= p.launch - 0.5 * self.dt {
                s.launched[j] = true;
                let tag = PROJECTILE_TAG + u32::try_from(j).unwrap_or(0);
                self.add(s, p.x0, p.v0, p.charge, p.charge / p.mass, tag);
            }
        }
        self.add_induced(&mut induced, s, before);
        let launched = s.x.len();
        let emitted = self.emit(s, &induced);
        self.add_induced(&mut induced, s, launched);
        let nodes = self.nodes();
        // The emitted charge left its electrode.
        for (e, q) in emitted.iter().enumerate() {
            if let Some(r) = self.node_of(&nodes, e) {
                s.node_charge[r] -= q;
            }
        }
        crate::cancel::checkpoint();
        // The potentials at the end of the step (sources, circuit nodes), then the field.
        let q_ind = self.electrodes.charges_of(&induced);
        self.update_potentials(s, s.t + self.dt, &q_ind);
        let sigma: Vec<f64> = self
            .static_density(&s.potentials)
            .iter()
            .zip(&induced)
            .map(|(a, b)| a + b)
            .collect();
        // In parallel, each particle's field whole by one thread (bit-identical).
        let fields = self.charges_fields(s, &sigma, s.x.len() >= TREE_MIN);
        let before = if self.lattices.is_empty() {
            Vec::new()
        } else {
            s.x.iter().copied().zip(s.v.iter().copied()).collect()
        };
        // A changing field: B at the positions' time, and its induced field below.
        let b_z = self.b_wave.map_or(self.b_z, |w| w.b(s.t));
        // With the charges' own magnetic field (Darwin) the positions move after the field's
        // solve (`darwin_step`).
        let moves = s.darwin.is_none();
        let v_start: Vec<DVec3> = if moves { Vec::new() } else { s.v.clone() };
        // The field's work on each lattice's carriers over the step.
        let mut work = vec![0.0; self.lattices.len()];
        let mut induced_work = vec![0.0; self.lattices.len()];
        for (((((x, v), qm), e), charge), tag) in
            s.x.iter_mut()
                .zip(&mut s.v)
                .zip(&s.qm)
                .zip(&fields)
                .zip(&s.charge)
                .zip(&s.tag)
        {
            let old = *v;
            let (e, induced) = self.with_external(*e, s.t, *x);
            let e = &e;
            let h = 0.5 * self.dt * qm;
            let minus = *v + *e * h;
            // Rotation about z by the magnetic field (Boris).
            let t = DVec3::new(0.0, 0.0, h * b_z);
            let prime = minus + minus.cross(t);
            let plus = minus + prime.cross(t * (2.0 / (1.0 + t.z * t.z)));
            *v = plus + *e * h;
            if !work.is_empty()
                && *tag < PROJECTILE_TAG
                && let Some(k) = self.lattices.iter().position(|l| l.contains(*x))
            {
                // Its cell's share too (the same sums, split by place).
                let c = self.lattices[k].account_cell(*x);
                let w = 0.5 * self.dt * charge * e.dot(old + *v);
                work[k] += w;
                s.cells[k].work[c] += w;
                if let Some(i) = induced {
                    let w = 0.5 * self.dt * charge * i.dot(old + *v);
                    induced_work[k] += w;
                    s.cells[k].induced[c] += w;
                }
            }
            if moves {
                *x += *v * self.dt;
            }
        }
        if !moves {
            self.darwin_step(s, &v_start);
        }
        for (total, w) in s.field_work.iter_mut().zip(&work) {
            *total += w;
        }
        for (total, w) in s.induced_work.iter_mut().zip(&induced_work) {
            *total += w;
        }
        s.t += self.dt;
        // The conductors: every carrier that was inside one scatters (lattice.rs).
        if !self.lattices.is_empty() {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let step = libm::round(s.t / self.dt) as u64;
            for (p, (a, v0)) in before.iter().enumerate() {
                if s.tag[p] >= PROJECTILE_TAG {
                    continue;
                }
                if let Some(k) = self.lattices.iter().position(|l| l.contains(*a)) {
                    let l = &self.lattices[k];
                    self.scatter(l, s, p, *a, *v0, step);
                    // The charge it carried over the step, along its chord.
                    l.book_flux(&mut s.cells[k].flux, *a, s.x[p], s.charge[p]);
                }
            }
        }
        // Absorption: a particle that crossed an electrode's boundary, came back to its
        // emitter (RETURN_ABSORB), or left.
        let n_el = self.electrodes.electrodes.len();
        let release: Vec<f64> = self
            .emitters
            .iter()
            .map(|em| {
                EMISSION_OFFSET
                    * self
                        .segments_of(em.electrode)
                        .map(|j| {
                            let (a, b) = self.electrodes.segments[j];
                            (b - a).length()
                        })
                        .fold(0.0, f64::max)
            })
            .collect();
        let mut k = 0;
        while k < s.x.len() {
            let crossed = (0..n_el)
                .find(|&e| inside(&self.electrodes.electrodes[e].section, s.x[k]) != s.side[k][e]);
            let mut returned = None;
            if let Some(em) = self.emitters.get(s.tag[k] as usize) {
                let delta = release[s.tag[k] as usize];
                let d = distance_outside(&self.electrodes.electrodes[em.electrode].section, s.x[k]);
                if d >= RETURN_LEAVE * delta {
                    s.left[k] = true;
                } else if s.left[k] && d < RETURN_ABSORB * delta {
                    returned = Some(em.electrode);
                }
            }
            let hit = crossed.or(returned);
            let out = self.arena.is_some_and(|(lo, hi)| {
                let x = s.x[k];
                x.x < lo.x || x.y < lo.y || x.x > hi.x || x.y > hi.y
            });
            if hit.is_some() || out {
                let tag = s.tag[k];
                if tag >= PROJECTILE_TAG {
                    let j = (tag - PROJECTILE_TAG) as usize;
                    s.fates[j] = Some(match hit {
                        Some(electrode) => Fate::Absorbed {
                            electrode,
                            t: s.t,
                            x: s.x[k],
                        },
                        None => Fate::Lost { t: s.t, x: s.x[k] },
                    });
                } else if let Some(e) = hit {
                    s.collected[e][tag as usize] += s.charge[k];
                    s.contact_energy[e] +=
                        0.5 * (s.charge[k] / s.qm[k]).abs() * s.v[k].length_squared();
                    if s.book_potential_steps && !self.lattices.is_empty() {
                        let phi = self.potential(s, &sigma, s.x[k], Some(k));
                        s.potential_steps[e] += s.charge[k] * (phi - s.potentials[e]);
                    }
                }
                // A circuit node keeps what it absorbs (projectiles' charge too).
                if let Some(r) = hit.and_then(|e| self.node_of(&nodes, e)) {
                    s.node_charge[r] += s.charge[k];
                }
                if hit.is_none() {
                    s.lost += s.charge[k];
                }
                s.x.swap_remove(k);
                s.v.swap_remove(k);
                s.charge.swap_remove(k);
                s.qm.swap_remove(k);
                s.tag.swap_remove(k);
                s.id.swap_remove(k);
                s.left.swap_remove(k);
                s.side.swap_remove(k);
                if !s.u.is_empty() {
                    s.u.swap_remove(k);
                }
            } else {
                k += 1;
            }
        }
        // The moved particles' induced density, for the next step and the display.
        s.induced = self.induced_density(s);
        for (i, b) in s.induced.iter_mut().zip(&s.background) {
            *i += b;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bem::Bias;
    use crate::lattice::{Lattice, Material};
    use crate::zinv::{Electrode, Electrodes, Section};

    /// The step's two ways of summing the charges' fields (`charges_fields`: the trees from
    /// `TREE_MIN` particles, the direct sums below) give the same fields for the same state, to
    /// the trees' accuracy: 1e-6 of the largest (`tree2d`'s series truncation is 5e-7 of a
    /// node's field; first set at 1e-9 from a comment here that claimed ~1e-11, which
    /// `tree2d`'s own bound contradicts; measured 1.2e-8). A metal of 2500 carriers over its
    /// background beside a floating electrode, after one step. Until 2026-10-26 the direct sums
    /// included the backgrounds and the push added them again: this test would have shown the
    /// backgrounds' whole field between the two.
    #[test]
    fn trees_and_direct_sums_agree() {
        let (lo, hi) = (DVec3::ZERO, DVec3::new(10.0, 4.0, 0.0));
        let cores = Lattice::random_cores(lo, hi, 2000, 0.01, 3);
        let lattice = Lattice::new(
            lo,
            hi,
            Material {
                background: 1.0,
                core_radius: 0.01,
                charge_per_mass: -1.0,
                thermal_speed: 1.0,
            },
            cores,
            2500,
            5,
        );
        let electrodes = Electrodes::new(
            vec![Electrode {
                section: Section::Rect {
                    center: DVec3::new(5.0, 6.0, 0.0),
                    angle: 0.0,
                    half_length: 2.0,
                    half_thickness: 0.25,
                },
                bias: Bias::Charge(0.0),
            }],
            &[],
            0.25,
        );
        let tube = Tube {
            arena: Some((DVec3::new(-2.0, -2.0, 0.0), DVec3::new(12.0, 9.0, 0.0))),
            lattices: vec![lattice],
            ..Tube::new(electrodes, 0.05, 0.02)
        };
        let mut s = tube.state();
        tube.step(&mut s);
        let sigma = tube.surface_density(&s);
        let trees = tube.charges_fields(&s, &sigma, true);
        let direct = tube.charges_fields(&s, &sigma, false);
        let top = direct.iter().map(|e| e.length()).fold(0.0, f64::max);
        let worst = trees
            .iter()
            .zip(&direct)
            .map(|(a, b)| (*a - *b).length())
            .fold(0.0, f64::max);
        println!("trees against direct sums: {worst:.2e} of the largest field {top:.3e}");
        assert!(worst < 1e-6 * top, "{worst} of {top}");
    }
}
