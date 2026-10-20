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
use crate::zinv::{Electrodes, PAR_CHUNK, Section, line_charge, segment_field_fast};

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

    /// The circuit nodes: electrodes with a circuit through a resistance or an inductance
    /// (not an ideal source), then floating electrodes (a fixed charge) without a circuit.
    fn nodes(&self) -> Vec<Node> {
        let mut out: Vec<Node> = self
            .circuits
            .iter()
            .filter(|c| !c.is_ideal() && !c.is_strap())
            .map(|c| Node {
                electrode: c.electrode,
                conductance_dt: if c.resistance.is_finite() && !c.is_inductive() {
                    self.dt / c.resistance
                } else {
                    0.0
                },
                capacitance: c.capacitance,
                source: c.source,
            })
            .collect();
        for (i, e) in self.electrodes.electrodes.iter().enumerate() {
            if matches!(e.bias, Bias::Charge(_)) && !self.circuits.iter().any(|c| c.electrode == i)
            {
                out.push(Node {
                    electrode: i,
                    conductance_dt: 0.0,
                    capacitance: 0.0,
                    source: Waveform::Dc(0.0),
                });
            }
        }
        out
    }

    /// The node index of an electrode, if it is a circuit node or strapped to one.
    fn node_of(&self, nodes: &[Node], electrode: usize) -> Option<usize> {
        let e = self.strap(electrode).unwrap_or(electrode);
        nodes.iter().position(|n| n.electrode == e)
    }

    /// The electrode a strapped one follows.
    fn strap(&self, electrode: usize) -> Option<usize> {
        self.circuits
            .iter()
            .find(|c| c.electrode == electrode && c.is_strap())
            .and_then(|c| c.to)
    }

    /// Node `r`'s electrodes: its own, then those strapped to it.
    fn members<'a>(&'a self, nodes: &'a [Node], r: usize) -> impl Iterator<Item = usize> + 'a {
        let own = nodes[r].electrode;
        std::iter::once(own).chain(
            self.circuits
                .iter()
                .filter(move |c| c.is_strap() && c.to == Some(own))
                .map(|c| c.electrode),
        )
    }

    /// Every strapped electrode at its other end's potential.
    fn follow_straps(&self, potentials: &mut [f64]) {
        for c in &self.circuits {
            if let (true, Some(to)) = (c.is_strap(), c.to) {
                potentials[c.electrode] = potentials[to];
            }
        }
    }

    /// The potentials at time `t`: ideal sources at their value, the nodes solved from
    /// their charges, with the particles' induced charge `q_ind` per electrode. A
    /// resistive branch delivers `dt (S + V_to − V)/R` over the step (backward Euler, at
    /// the step's end); an inductive one `dt (I₀ + I₁)/2`, its current from
    /// `L (I₁ − I₀) = dt/2 [(S + V_to − V − R I)₀ + (S + V_to − V − R I)₁]` (the trapezoidal
    /// rule: with the capacitances, the implicit midpoint rule, which keeps an LC tank's
    /// energy exactly). The unknowns are the nodes' potentials and the inductive branches'
    /// currents at `t`, from one linear system.
    #[allow(clippy::too_many_lines)]
    fn update_potentials(&self, s: &mut TubeState, t: f64, q_ind: &[f64]) {
        let old = s.potentials.clone();
        for c in &self.circuits {
            if c.is_ideal() {
                s.potentials[c.electrode] = c.source.value(t);
            }
        }
        self.follow_straps(&mut s.potentials);
        let nodes = self.nodes();
        if nodes.is_empty() {
            return;
        }
        let k = nodes.len();
        let inductive: Vec<usize> = (0..self.circuits.len())
            .filter(|&c| self.circuits[c].is_inductive())
            .collect();
        let size = k + inductive.len();
        let n_el = self.electrodes.electrodes.len();
        let cap = &self.electrodes.capacitance;
        let mut a = vec![0.0; size * size];
        let mut b = vec![0.0; size];
        // Charge balance of each node (its own resistive branch to ground as before), over
        // its electrodes (its own and those strapped to it).
        for (r, n) in nodes.iter().enumerate() {
            for c in 0..k {
                a[r * size + c] = self
                    .members(&nodes, r)
                    .flat_map(|i| self.members(&nodes, c).map(move |j| cap[i][j]))
                    .sum();
            }
            a[r * size + r] += n.capacitance + n.conductance_dt;
            let fixed: f64 = self
                .members(&nodes, r)
                .map(|i| {
                    (0..n_el)
                        .filter(|&j| self.node_of(&nodes, j).is_none())
                        .map(|j| cap[i][j] * s.potentials[j])
                        .sum::<f64>()
                })
                .sum();
            let induced: f64 = self.members(&nodes, r).map(|i| q_ind[i]).sum();
            b[r] = s.node_charge[r] + n.conductance_dt * n.source.value(t) - fixed - induced;
        }
        // A resistive branch returning to another electrode: its potential in the
        // branch's drop, and the charge out of it if it is a node.
        for c in &self.circuits {
            let Some(to) = c.to else { continue };
            if c.is_ideal() || c.is_strap() || c.is_inductive() || !c.resistance.is_finite() {
                continue;
            }
            let g = self.dt / c.resistance;
            let r = self
                .node_of(&nodes, c.electrode)
                .expect("a circuit's electrode is a node");
            match self.node_of(&nodes, to) {
                Some(q) => {
                    a[r * size + q] -= g;
                    a[q * size + q] += g;
                    a[q * size + r] -= g;
                    b[q] -= g * c.source.value(t);
                }
                None => b[r] += g * s.potentials[to],
            }
        }
        // The inductive branches: their charge into (and out of) the nodes, and their
        // own equations.
        let h = 0.5 * self.dt;
        for (j, &ci) in inductive.iter().enumerate() {
            let c = &self.circuits[ci];
            let u = k + j;
            let i0 = s.inductor_current[ci];
            let r = self
                .node_of(&nodes, c.electrode)
                .expect("a circuit's electrode is a node");
            a[r * size + u] -= h;
            b[r] += h * i0;
            let (q, v_to0, v_to1) = match c.to {
                Some(to) => (self.node_of(&nodes, to), old[to], s.potentials[to]),
                None => (None, 0.0, 0.0),
            };
            if let Some(q) = q {
                a[q * size + u] += h;
                b[q] -= h * i0;
                a[u * size + q] -= h;
            } else {
                b[u] += h * v_to1;
            }
            a[u * size + u] = c.inductance + h * c.resistance;
            a[u * size + r] += h;
            b[u] += c.inductance * i0
                + h * (c.source.value(t - self.dt) + v_to0 - old[c.electrode] - c.resistance * i0)
                + h * c.source.value(t);
        }
        let x = Lu::new(a, size).solve(&b);
        for (r, n) in nodes.iter().enumerate() {
            s.potentials[n.electrode] = x[r];
        }
        self.follow_straps(&mut s.potentials);
        // The charges the branches delivered.
        for (r, n) in nodes.iter().enumerate() {
            if n.conductance_dt == 0.0 {
                continue;
            }
            let Some(c) = self
                .circuits
                .iter()
                .position(|c| c.electrode == n.electrode && !c.is_ideal())
            else {
                continue;
            };
            let circuit = &self.circuits[c];
            let delivered = match circuit.to {
                None => n.conductance_dt * (n.source.value(t) - x[r]),
                Some(to) => n.conductance_dt * (n.source.value(t) + s.potentials[to] - x[r]),
            };
            s.node_charge[r] += delivered;
            if let Some(q) = circuit.to.and_then(|to| self.node_of(&nodes, to)) {
                s.node_charge[q] -= delivered;
            }
            s.circuit_current[c] = delivered / self.dt;
        }
        for (j, &ci) in inductive.iter().enumerate() {
            let c = &self.circuits[ci];
            let i1 = x[k + j];
            let delivered = h * (s.inductor_current[ci] + i1);
            let r = self
                .node_of(&nodes, c.electrode)
                .expect("a circuit's electrode is a node");
            s.node_charge[r] += delivered;
            if let Some(q) = c.to.and_then(|to| self.node_of(&nodes, to)) {
                s.node_charge[q] -= delivered;
            }
            s.inductor_current[ci] = i1;
            s.circuit_current[ci] = delivered / self.dt;
        }
    }

    /// The surface density of the electrodes at their potentials (no particles).
    fn static_density(&self, potentials: &[f64]) -> Vec<f64> {
        let mut sigma = vec![0.0; self.electrodes.segments.len()];
        for (j, v) in potentials.iter().enumerate() {
            if *v != 0.0 {
                for (s, u) in sigma.iter_mut().zip(self.electrodes.unit(j)) {
                    *s += v * u;
                }
            }
        }
        sigma
    }

    /// The particles' induced surface density (all electrodes grounded).
    fn induced_density(&self, s: &TubeState) -> Vec<f64> {
        if s.x.len() >= TREE_MIN {
            let tree = Tree::new(s.x.iter().map(|p| Source::point(*p)).collect(), 0.0);
            let m = tree.multipoles(&s.charge);
            return self.electrodes.induced_by(|y| tree.potential(&m, y));
        }
        self.electrodes.induced_by(|y| {
            s.x.iter()
                .zip(&s.charge)
                .map(|(p, l)| line_charge(*p, *l, y).0)
                .sum()
        })
    }

    fn segments_of(&self, electrode: usize) -> impl Iterator<Item = usize> + '_ {
        (0..self.electrodes.segments.len())
            .filter(move |&j| self.electrodes.owner_of(j) == electrode)
    }

    /// The total surface density: the electrodes' at their present potentials plus the
    /// particles' induced one (as of the state's last step).
    pub fn surface_density(&self, s: &TubeState) -> Vec<f64> {
        self.static_density(&s.potentials)
            .iter()
            .zip(&s.induced)
            .map(|(a, b)| a + b)
            .collect()
    }

    /// Adds to `induced` the surface density induced by the particles from index `from` on
    /// (those added during a step: launched or emitted).
    fn add_induced(&self, induced: &mut [f64], s: &TubeState, from: usize) {
        if s.x.len() <= from {
            return;
        }
        let extra = self.electrodes.induced_by(|y| {
            s.x[from..]
                .iter()
                .zip(&s.charge[from..])
                .map(|(p, l)| line_charge(*p, *l, y).0)
                .sum()
        });
        for (a, b) in induced.iter_mut().zip(&extra) {
            *a += b;
        }
    }

    /// The field at `x` of the surface charge `sigma` and of every particle but `skip`.
    fn field(&self, s: &TubeState, sigma: &[f64], x: DVec3, skip: Option<usize>) -> DVec3 {
        let mut e = DVec3::ZERO;
        for ((p, q), d) in self.electrodes.segments.iter().zip(sigma) {
            e += segment_field_fast(*p, *q, x) * *d;
        }
        let h2 = self.softening * self.softening;
        for (k, (p, l)) in s.x.iter().zip(&s.charge).enumerate() {
            if Some(k) == skip {
                continue;
            }
            let r = DVec3::new(x.x - p.x, x.y - p.y, 0.0);
            let r2 = (r.x * r.x + r.y * r.y).max(h2);
            e += r * (2.0 * l / r2);
        }
        for l in &self.lattices {
            e += l.background_at(x).1;
        }
        e
    }

    fn add(&self, s: &mut TubeState, x: DVec3, v: DVec3, charge: f64, qm: f64, tag: u32) {
        s.x.push(x);
        s.v.push(v);
        s.charge.push(charge);
        s.qm.push(qm);
        s.tag.push(tag);
        s.id.push(s.next_id);
        s.next_id = s.next_id.wrapping_add(1);
        s.left.push(false);
        s.side.push(
            self.electrodes
                .electrodes
                .iter()
                .map(|e| inside(&e.section, x))
                .collect(),
        );
    }

    /// The outward normal of `electrode`'s segment from `a` to `b` at its point `at`: the
    /// side away from the electrode's section.
    fn outward(&self, electrode: usize, a: DVec3, b: DVec3, at: DVec3) -> DVec3 {
        let t = (b - a).normalize();
        let n = DVec3::new(-t.y, t.x, 0.0);
        if inside(
            &self.electrodes.electrodes[electrode].section,
            at + n * 1e-6,
        ) {
            -n
        } else {
            n
        }
    }

    /// The charge of its carriers' sign (a magnitude) that an emitter's emitting face holds
    /// with the electrodes at `potentials` and no particles: what space-charge-limited
    /// emission releases first, the scale of the charge in flight (the macroparticles'
    /// weight, `level::tube`).
    pub fn emitting_charge(&self, emitter: usize, potentials: &[f64]) -> f64 {
        let em = &self.emitters[emitter];
        let sign = em.charge_per_mass.signum();
        let sigma = self.static_density(potentials);
        self.segments_of(em.electrode)
            .map(|j| {
                let (a, b) = self.electrodes.segments[j];
                let out = self.outward(em.electrode, a, b, (a + b) * 0.5);
                let q = sigma[j] * (b - a).length();
                if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) || q * sign <= 0.0 {
                    0.0
                } else {
                    q.abs()
                }
            })
            .sum()
    }

    /// Space-charge-limited emission of every emitter: each segment of its electrode
    /// whose surface density (the electrodes' at their potentials plus `induced`, the
    /// particles') has its carriers' sign releases that charge.
    fn emit(&self, s: &mut TubeState, induced: &[f64]) -> Vec<f64> {
        let mut emitted = vec![0.0; self.electrodes.electrodes.len()];
        let sigma: Vec<f64> = self
            .static_density(&s.potentials)
            .iter()
            .zip(induced)
            .map(|(a, b)| a + b)
            .collect();
        for (i, em) in self.emitters.iter().enumerate() {
            let sign = em.charge_per_mass.signum();
            let segs: Vec<usize> = self.segments_of(em.electrode).collect();
            if let Some(ohmic) = em.ohmic {
                self.emit_ohmic(s, i, em, ohmic, &segs, &mut emitted);
                continue;
            }
            for (k, &j) in segs.iter().enumerate() {
                let (a, b) = self.electrodes.segments[j];
                let charge = sigma[j] * (b - a).length();
                if charge * sign > 0.0 {
                    // The space-charge limit (the surface's charge), or the temperature
                    // limit over this step, whichever is less.
                    let cap = em
                        .saturation
                        .map_or(f64::INFINITY, |j_sat| j_sat * (b - a).length() * self.dt);
                    s.pending[i][k] += charge.abs().min(cap);
                }
                while s.pending[i][k] >= em.weight {
                    if s.x.len() >= MAX_PARTICLES {
                        s.capped = true;
                        break;
                    }
                    s.pending[i][k] -= em.weight;
                    // Low-discrepancy point along the segment (golden ratio).
                    s.emitted += 1;
                    #[allow(clippy::cast_precision_loss)]
                    let f = (s.emitted as f64 * 0.618_033_988_749_894_9).fract();
                    let mid = a + (b - a) * f;
                    let out = self.outward(em.electrode, a, b, mid);
                    if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) {
                        // Not the coated face: its charge stays.
                        s.pending[i][k] = 0.0;
                        break;
                    }
                    let x = mid + out * (EMISSION_OFFSET * (b - a).length());
                    let tag = u32::try_from(i).unwrap_or(u32::MAX);
                    let v = if em.thermal_speed > 0.0 {
                        crate::lattice::flux_velocity(out, em.thermal_speed, s.emitted)
                    } else {
                        DVec3::ZERO
                    };
                    self.add(s, x, v, sign * em.weight, em.charge_per_mass, tag);
                    emitted[em.electrode] += sign * em.weight;
                    s.released[i] += sign * em.weight;
                }
            }
        }
        emitted
    }

    /// An ohmic contact's emission (`OhmicContact`): per emitting segment of emitter `i`,
    /// the carriers its layer lacks.
    fn emit_ohmic(
        &self,
        s: &mut TubeState,
        i: usize,
        em: &Emitter,
        ohmic: OhmicContact,
        segs: &[usize],
        emitted: &mut [f64],
    ) {
        let sign = em.charge_per_mass.signum();
        let tag = u32::try_from(i).unwrap_or(u32::MAX);
        for &j in segs {
            let (a, b) = self.electrodes.segments[j];
            let len = (b - a).length();
            let t = (b - a) / len;
            let out = self.outward(em.electrode, a, b, (a + b) * 0.5);
            if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) {
                continue;
            }
            // The carriers of its sign in the layer now.
            let held: f64 =
                s.x.iter()
                    .zip(&s.charge)
                    .zip(&s.tag)
                    .filter(|((_, q), tag)| **tag < PROJECTILE_TAG && **q * sign > 0.0)
                    .filter(|((x, _), _)| {
                        let r = **x - a;
                        let (u, w) = (r.dot(t), r.dot(out));
                        (0.0..=len).contains(&u) && w > 0.0 && w <= ohmic.depth
                    })
                    .map(|((_, q), _)| q.abs())
                    .sum();
            let mut lack = ohmic.density * len * ohmic.depth - held;
            while lack >= em.weight {
                if s.x.len() >= MAX_PARTICLES {
                    s.capped = true;
                    break;
                }
                lack -= em.weight;
                s.emitted += 1;
                #[allow(clippy::cast_precision_loss)]
                let f = (s.emitted as f64 * 0.618_033_988_749_894_9).fract();
                // From the face, as a gas at the lattice's temperature crosses a wall: the
                // flux-weighted law, 1.5 kT each in the plane, what the carriers the contact
                // takes in carry on average.
                let v = crate::lattice::flux_velocity(out, em.thermal_speed, s.emitted);
                let x = a + t * (len * f) + out * (EMISSION_OFFSET * len);
                self.add(s, x, v, sign * em.weight, em.charge_per_mass, tag);
                emitted[em.electrode] += sign * em.weight;
                s.released[i] += sign * em.weight;
            }
        }
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
        let fields: Vec<DVec3> = if s.x.len() >= TREE_MIN {
            // Many particles: the segments' and the particles' trees (`tree2d`, the direct
            // sums to ~1e-11).
            let segments = &self.electrodes.tree;
            let segment_m = segments.multipoles(&sigma);
            let particles = Tree::new(
                s.x.iter().map(|p| Source::point(*p)).collect(),
                self.softening,
            );
            let particle_m = particles.multipoles(&s.charge);
            (0..s.x.len())
                .into_par_iter()
                .with_min_len(PAR_CHUNK)
                .map(|k| {
                    segments.field(&segment_m, s.x[k], None)
                        + particles.field(&particle_m, s.x[k], Some(k))
                })
                .collect()
        } else {
            (0..s.x.len())
                .into_par_iter()
                .with_min_len(PAR_CHUNK)
                .map(|k| self.field(s, &sigma, s.x[k], Some(k)))
                .collect()
        };
        let before = if self.lattices.is_empty() {
            Vec::new()
        } else {
            s.x.iter().copied().zip(s.v.iter().copied()).collect()
        };
        // A changing field: B at the positions' time, and its induced field below.
        let b_z = self.b_wave.map_or(self.b_z, |w| w.b(s.t));
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
            // The conductors' backgrounds (no sum at all without one: bit for bit as before).
            let e = if self.lattices.is_empty() {
                *e
            } else {
                *e + self
                    .lattices
                    .iter()
                    .map(|l| l.background_at(*x).1)
                    .sum::<DVec3>()
            };
            let induced = self.b_wave.map(|w| w.induced(s.t, *x));
            let e = match induced {
                Some(i) => e + i,
                None => e,
            };
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
            *x += *v * self.dt;
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

    /// The field at `x` now: the segments' charge, every particle's (softened) and the
    /// lattices' backgrounds, and a changing field's induced part.
    #[must_use]
    pub fn field_at(&self, s: &TubeState, x: DVec3) -> DVec3 {
        let sigma = self.surface_density(s);
        let e = self.field(s, &sigma, x, None);
        match &self.b_wave {
            Some(w) => e + w.induced(s.t, x),
            None => e,
        }
    }

    /// The carrier `p` of `lattice` after a move from `a` (its velocity `v0` before the
    /// step's kick): through the cores it entered (`Lattice::advance`: bounces, new
    /// thermal speeds, the lattice taking the energy), then the walls (reflecting, but
    /// where an electrode lies beyond: a contact, which absorbs the carrier below).
    fn scatter(
        &self,
        lattice: &Lattice,
        s: &mut TubeState,
        p: usize,
        a: DVec3,
        v0: DVec3,
        step: u64,
    ) {
        if let Some(m) = lattice.advance(a, v0, s.v[p], self.dt, s.id[p], step) {
            let mass = (s.charge[p] / s.qm[p]).abs();
            s.joule += mass * m.heat_per_mass;
            s.collisions += u64::from(m.collisions);
            s.x[p] = m.x;
            s.v[p] = m.v;
        }
        // Off its sides (out of the rectangle, not into a contact) and off its internal walls
        // (the step's chord from `from` crossing one), until it is inside.
        let mut from = a;
        for _ in 0..8 {
            let x = s.x[p];
            if lattice.contains(x) {
                let Some((t, w0, w1)) = lattice
                    .walls
                    .iter()
                    .filter_map(|&(w0, w1)| crossing(from, x, w0, w1).map(|t| (t, w0, w1)))
                    .min_by(|a, b| a.0.total_cmp(&b.0))
                else {
                    break;
                };
                // Specular: the position and the velocity mirrored in the wall's line.
                let e = w1 - w0;
                let n = DVec3::new(-e.y, e.x, 0.0) / libm::sqrt(e.x * e.x + e.y * e.y);
                let hit = from + (x - from) * t;
                let mirrored = x - n * (2.0 * (x - w0).dot(n));
                let v = s.v[p];
                s.v[p] = v - n * (2.0 * v.dot(n));
                s.x[p] = mirrored;
                // The next chord starts just off the wall, on the carrier's side.
                from = hit + (mirrored - hit) * 1e-6;
                continue;
            }
            if self
                .electrodes
                .electrodes
                .iter()
                .any(|e| inside(&e.section, x))
            {
                break;
            }
            let (lo, hi) = (lattice.lo, lattice.hi);
            let (mut x, mut v) = (x, s.v[p]);
            if x.x < lo.x {
                x.x = 2.0 * lo.x - x.x;
                v.x = -v.x;
            } else if x.x > hi.x {
                x.x = 2.0 * hi.x - x.x;
                v.x = -v.x;
            }
            if x.y < lo.y {
                x.y = 2.0 * lo.y - x.y;
                v.y = -v.y;
            } else if x.y > hi.y {
                x.y = 2.0 * hi.y - x.y;
                v.y = -v.y;
            }
            s.x[p] = x;
            s.v[p] = v;
        }
    }
}
