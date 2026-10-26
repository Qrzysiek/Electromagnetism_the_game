//! Level format: grid, shots (particle, launch, detector), fixed elements, coils, limits
//! and reference solutions (SPEC §3, §7.9).
//!
//! All positions are integer grid nodes. A node `n` lies at `n / subdivision` cell units.
//! Refining the grid by an integer factor multiplies every node and the subdivision, so
//! all existing nodes stay exactly where they were.
//!
//! Format history: version 1 had a single particle/launch/detector and `level_charges`
//! with `"charge"` values; it is migrated on load. Version 2 has `shots` and `elements`,
//! and optionally `disturbances`.
//!
//! **Flights:** every shot is flown once under each disturbance (once, undisturbed, if the
//! level has none). A setup solves the level when every flight arrives.

pub mod analysis;
pub mod beam;
pub mod cost;
pub mod model;
pub mod run_cache;
pub mod solve;
pub mod tube;
pub mod tube_em;

use physics::DVec3;
use physics::antenna::OscillatingDipole;
use physics::conductor::{Bias, Conductors, Resolution, SphereConductor};
use physics::dynamics::{Kinematics, Particle};
use physics::external::{External, PlaneWave};
use physics::field::{ChargeCloud, Coulomb, FixedCharge, LevelField};
use physics::geometry::{Aabb, Capsule, Region, Shape, Sphere, Torus};
use physics::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};
use physics::trajectory::Scenario;
use physics::verify::{Tolerances, Verification, verify_pair};
use serde::{Deserialize, Serialize};

mod drives;
mod fields;
mod geometry;
mod issues;
mod placement;
mod scenarios;
pub use geometry::Outline;
use geometry::{rect_contains, rect_corners, rect_distance};
pub use issues::{IMAGE_FORCE_LIMIT, electrode_image_force_bound, moment_radiation_estimate};

pub const FORMAT_VERSION: u32 = 2;
pub const ENGINE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Integer grid node `[x, y, z]`.
pub type Node = [i64; 3];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(from = "LevelFile")]
pub struct Level {
    pub format_version: u32,
    /// Physics engine version the level (and its reference solution) was made with.
    pub engine_version: String,
    pub name: String,
    pub description: String,
    /// What the level represents in the real world, and how its model differs from that
    /// (shown with the description; text only, it changes no run).
    #[serde(default, skip_serializing_if = "Reality::is_empty")]
    pub reality: Reality,
    pub grid: Grid,
    pub physics: WorldPhysics,
    /// Particles to deliver; one setup must bring every shot to its detector.
    pub shots: Vec<Shot>,
    /// Elements placed by the level (charges and magnets; obstacles).
    pub elements: Vec<Element>,
    /// Coils placed by the level.
    pub coils: Vec<Coil>,
    pub limits: Limits,
    /// A known solution (player elements), if any.
    pub reference_solution: Vec<Element>,
    /// External disturbances; one setup must work under each of them.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub disturbances: Vec<Disturbance>,
    /// Metal spheres placed by the level (PHYSICS.md §2.6).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub conductors: Vec<Conductor>,
    /// Charge clouds placed by the level: uniformly charged spheres that particles fly
    /// through (Thomson's atom, PHYSICS.md §2.1).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub clouds: Vec<Cloud>,
    /// Dynamic particles placed by the level (targets, partners): they fly with the
    /// shots as one interacting system; those with a detector must arrive.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub free_particles: Vec<FreeParticle>,
    /// Box electrodes (plates, slabs, walls) placed by the level (PHYSICS.md §2.7).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub electrodes: Vec<Electrode>,
    /// Dielectric boxes placed by the level (PHYSICS.md §2.7): their bound charge screens
    /// fields partially; particles do not enter them.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub dielectrics: Vec<Dielectric>,
    /// Dielectric spheres placed by the level (PHYSICS.md §2.7), as the boxes.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub dielectric_spheres: Vec<DielectricSphere>,
    /// Ferrite plates (high-μ insulators) placed by the level (PHYSICS.md §2.7): they
    /// bend magnetic fields only.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ferrites: Vec<FerritePlate>,
    /// Stages of the instrument: boxes every flight must pass, in order, before its
    /// detector counts, each with optional conditions on the entering particle
    /// (PHYSICS.md §6.2).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub gates: Vec<Detector>,
    /// A tube level (z-invariant; `tube/`, docs/TUBES.md): its cathode and current goal.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tube: Option<tube::TubeSpec>,
}

/// What a level represents and how far it is the real thing (the owner, 2026-10-21: leaving
/// it out gives a wrong impression of what a level shows). Written per level in
/// `scripts/reality.py`; every shipped level has it (test `every_level_states_its_reality`).
/// The numerical model's own approximations are `Level::model_notes`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Reality {
    /// The real thing the level shows, in a few words.
    pub represents: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub remarks: Vec<Remark>,
}

impl Reality {
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.represents.is_empty() && self.remarks.is_empty()
    }
}

/// One way the level is, or is not, the real thing: a short line and its explanation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Remark {
    pub kind: RemarkKind,
    pub line: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub detail: String,
}

/// How a remark relates the level to reality.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum RemarkKind {
    /// As in reality: what carries over.
    Kept,
    /// Idealized, scaled or approximated.
    #[default]
    Simplified,
    /// Absent from the model.
    Omitted,
}

impl RemarkKind {
    pub const ALL: [RemarkKind; 3] = [Self::Kept, Self::Simplified, Self::Omitted];

    /// The tag the panel shows.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Kept => "as real",
            Self::Simplified => "simplified",
            Self::Omitted => "left out",
        }
    }
}

/// A rectangular metal box standing on the plane (symmetric about it): a plate, slab or
/// wall electrode.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Electrode {
    /// Centre (grid node).
    pub center: Node,
    /// In-plane length along `angle_deg`, thickness across it, and height (z), in cells.
    pub length: f64,
    pub thickness: f64,
    pub height: f64,
    #[serde(default)]
    pub angle_deg: f64,
    /// The cross-section (tube levels only; the 3D levels' electrodes are boxes): a box,
    /// a rod (a circle of diameter `length`) or a ring (a tube's wall seen end-on: an
    /// annulus of mid-diameter `length` and wall `thickness`, open inside).
    #[serde(default, skip_serializing_if = "ElectrodeShape::is_box")]
    pub shape: ElectrodeShape,
    pub bias: ConductorBias,
    /// The player sets this electrode's potential with a power supply (an element of kind
    /// `Supply` on its centre, one of `limits.supply_voltages`); without one it keeps
    /// `bias`. A driven electrode's supply sets its source's level instead (`Source::
    /// with_level`); the bias is then its potential at t = 0.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub tunable: bool,
    /// A circuit driving its potential (PHYSICS.md §2.10).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub drive: Option<Drive>,
    /// Its function in the device (cathode, grid, source, gate, deflection plate, …):
    /// shown with its number wherever it is named. Text only.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
}

/// An electrode's cross-section (`Electrode::shape`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ElectrodeShape {
    #[default]
    Box,
    /// A circle of diameter `length` (a wire, a cylinder seen end-on).
    Rod,
    /// An annulus of mid-diameter `length` and wall `thickness` (a coaxial anode).
    Ring,
}

impl ElectrodeShape {
    #[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
    fn is_box(&self) -> bool {
        *self == ElectrodeShape::Box
    }
}

/// A dielectric box standing on the plane (symmetric about it), of relative permittivity
/// `permittivity` (≥ 1; PHYSICS.md §2.7).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Dielectric {
    /// Centre (grid node).
    pub center: Node,
    /// In-plane length along `angle_deg`, thickness across it, and height (z), in cells.
    pub length: f64,
    pub thickness: f64,
    pub height: f64,
    #[serde(default)]
    pub angle_deg: f64,
    pub permittivity: f64,
}

/// A ferrite box standing on the plane (symmetric about it): magnetically soft, of
/// relative permeability `permeability` (≥ 1), electrically an insulator.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct FerritePlate {
    pub center: Node,
    pub length: f64,
    pub thickness: f64,
    pub height: f64,
    #[serde(default)]
    pub angle_deg: f64,
    pub permeability: f64,
}

/// A dielectric sphere centred in the plane, of relative permittivity `permittivity`
/// (≥ 1; PHYSICS.md §2.7).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DielectricSphere {
    pub center: Node,
    /// Radius in cells.
    pub radius: f64,
    pub permittivity: f64,
}

/// A charge cloud: a sphere of uniform charge density that particles fly through. Inside,
/// a charge of the opposite sign is bound harmonically, `ω₀² = |qQ| / (m R³)`; outside, it
/// is a point charge.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cloud {
    pub center: Node,
    /// Radius in cells.
    pub radius: f64,
    /// Total charge.
    pub charge: f64,
}

/// A conducting (metal) sphere placed by the level.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Conductor {
    pub center: Node,
    /// Radius in cells.
    pub radius: f64,
    pub bias: ConductorBias,
}

/// How a metal sphere is held: grounded, isolated with a net charge, or at a potential.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "lowercase")]
pub enum ConductorBias {
    Grounded,
    Charge(f64),
    Potential(f64),
}

/// A particle touches a metal sphere at this distance from its surface (cells). A point
/// charge at the surface would feel an infinite image force (PHYSICS.md §2.6).
pub const CONTACT_DISTANCE: f64 = 0.02;

/// One realization of fields from outside the arena (PHYSICS.md §2.3): uniform stray
/// fields and plane waves travelling in the plane.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Disturbance {
    #[serde(default)]
    pub name: String,
    /// Uniform stray electric field `[E_x, E_y]` in the plane.
    #[serde(default)]
    pub e: [f64; 2],
    /// Uniform stray magnetic field `B_z` (perpendicular to the plane).
    #[serde(default)]
    pub bz: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub waves: Vec<Wave>,
}

/// A plane wave travelling in the plane, polarized in the plane (its B is along z):
/// `E = E₀ ê cos(ω (t − k̂·x/c) + φ)` with `ê = ẑ × k̂`. For `c = ∞` it is a uniform field
/// oscillating in time ("mains hum").
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Wave {
    /// Amplitude `E₀`.
    pub amplitude: f64,
    /// Propagation direction, degrees counter-clockwise from +x.
    pub direction_deg: f64,
    /// Angular frequency `ω` (0 gives a static uniform field).
    pub omega: f64,
    /// Phase `φ` at `t = 0`, `x = 0`, in degrees.
    pub phase_deg: f64,
}

impl Disturbance {
    /// The external field terms of this disturbance in a world with speed of light `c`.
    pub fn fields(&self, c: f64) -> Vec<External> {
        let mut out = Vec::new();
        if self.e != [0.0, 0.0] || self.bz != 0.0 {
            out.push(External::Uniform {
                e: DVec3::new(self.e[0], self.e[1], 0.0),
                b: DVec3::new(0.0, 0.0, self.bz),
            });
        }
        out.extend(self.waves.iter().map(|w| {
            External::Wave(PlaneWave::in_plane(
                w.amplitude,
                w.direction_deg.to_radians(),
                w.omega,
                w.phase_deg.to_radians(),
                c,
            ))
        }));
        out
    }
}

/// One particle to deliver: its species, launch and detector.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Shot {
    pub particle: ParticleSpec,
    pub launch: Launch,
    pub detector: Detector,
    /// Fired as a beam of many particles (`beam.rs`); otherwise a single particle.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beam: Option<beam::BeamSpec>,
}

/// The on-disk form, accepting both format versions.
#[derive(Deserialize)]
struct LevelFile {
    format_version: u32,
    engine_version: String,
    name: String,
    #[serde(default)]
    description: String,
    #[serde(default)]
    reality: Reality,
    grid: Grid,
    physics: WorldPhysics,
    #[serde(default)]
    shots: Vec<Shot>,
    // Version 1 single shot.
    particle: Option<ParticleSpec>,
    launch: Option<Launch>,
    detector: Option<Detector>,
    #[serde(default, alias = "level_charges")]
    elements: Vec<Element>,
    #[serde(default)]
    coils: Vec<Coil>,
    limits: Limits,
    #[serde(default)]
    reference_solution: Vec<Element>,
    #[serde(default)]
    disturbances: Vec<Disturbance>,
    #[serde(default)]
    conductors: Vec<Conductor>,
    #[serde(default)]
    clouds: Vec<Cloud>,
    #[serde(default)]
    free_particles: Vec<FreeParticle>,
    #[serde(default)]
    electrodes: Vec<Electrode>,
    #[serde(default)]
    dielectrics: Vec<Dielectric>,
    #[serde(default)]
    dielectric_spheres: Vec<DielectricSphere>,
    #[serde(default)]
    ferrites: Vec<FerritePlate>,
    #[serde(default)]
    gates: Vec<Detector>,
    #[serde(default)]
    tube: Option<tube::TubeSpec>,
}

impl From<LevelFile> for Level {
    fn from(f: LevelFile) -> Self {
        let mut shots = f.shots;
        if let (Some(particle), Some(launch), Some(detector)) = (f.particle, f.launch, f.detector) {
            shots.insert(
                0,
                Shot {
                    particle,
                    launch,
                    detector,
                    beam: None,
                },
            );
        }
        Level {
            format_version: FORMAT_VERSION.max(f.format_version),
            engine_version: f.engine_version,
            name: f.name,
            description: f.description,
            reality: f.reality,
            grid: f.grid,
            physics: f.physics,
            shots,
            elements: f.elements,
            coils: f.coils,
            limits: f.limits,
            reference_solution: f.reference_solution,
            disturbances: f.disturbances,
            conductors: f.conductors,
            clouds: f.clouds,
            free_particles: f.free_particles,
            electrodes: f.electrodes,
            dielectrics: f.dielectrics,
            dielectric_spheres: f.dielectric_spheres,
            ferrites: f.ferrites,
            gates: f.gates,
            tube: f.tube,
        }
    }
}

/// Grid of `nx × ny × nz` cells (`nz = 0` for a 2D level in the plane z = 0).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Grid {
    pub nx: u32,
    pub ny: u32,
    pub nz: u32,
    /// Nodes per cell edge (1 = the level's recommended grid).
    pub subdivision: u32,
}

impl Grid {
    pub fn is_2d(&self) -> bool {
        self.nz == 0
    }

    /// Largest node index along each axis.
    pub fn max_node(&self) -> Node {
        let s = i64::from(self.subdivision);
        [
            i64::from(self.nx) * s,
            i64::from(self.ny) * s,
            i64::from(self.nz) * s,
        ]
    }

    pub fn contains(&self, n: Node) -> bool {
        let m = self.max_node();
        (0..3).all(|i| (0..=m[i]).contains(&n[i]))
    }

    /// Position of a node in cell units.
    #[allow(clippy::cast_precision_loss)] // node indices are far below 2^53
    pub fn position(&self, n: Node) -> DVec3 {
        let s = f64::from(self.subdivision);
        DVec3::new(n[0] as f64 / s, n[1] as f64 / s, n[2] as f64 / s)
    }
}

fn default_radius() -> f64 {
    0.3
}

fn default_wire_radius() -> f64 {
    0.1
}

#[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
fn is_zero(v: &f64) -> bool {
    *v == 0.0
}

/// Orientations the player can give an antenna, degrees from +x (the opposite ones are
/// the same antenna with the sign of its amplitude flipped).
pub const ANTENNA_ANGLES: [f64; 4] = [0.0, 45.0, 90.0, 135.0];

/// Orientations of player plates (along x or along y); any angle in hardcore mode.
pub const PLATE_ANGLES: [f64; 2] = [0.0, 90.0];

/// Smallest gap between the surfaces of two electrodes when one of them is a player
/// plate, in cells: two panels of the preview mesh (0.5 cells, PHYSICS.md §2.7), so that
/// the surface charge between them is resolved.
pub const PLATE_CLEARANCE: f64 = 1.0;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WorldPhysics {
    /// Speed of light in internal units; `None` means Newtonian mechanics.
    pub c: Option<f64>,
    /// Radius of every fixed charge, in cell units.
    pub charge_radius: f64,
    /// Radius of every magnet (uniformly magnetized sphere), in cell units.
    #[serde(default = "default_radius")]
    pub magnet_radius: f64,
    /// Radius of coil wires, in cell units.
    #[serde(default = "default_wire_radius")]
    pub wire_radius: f64,
    /// Radius of every antenna body, in cell units.
    #[serde(default = "default_radius")]
    pub antenna_radius: f64,
    /// Angular frequency of the RF generator that drives every antenna, in phase
    /// (`p(t) = p₀ cos(ωt)`; 0 means static dipoles).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub rf_omega: f64,
    /// Include the particles' radiation reaction (Landau–Lifshitz, PHYSICS.md §3.1). Off:
    /// radiation is neglected and must be negligible.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub radiation_reaction: bool,
    /// The particles of beams interact (PHYSICS.md §3.3): Coulomb for c = ∞; for finite c
    /// quasi-static (fields of the present state continued with constant acceleration), or the exact
    /// retarded fields with `beam_retarded`.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub beam_interaction: bool,
    /// Exact retarded (Liénard–Wiechert) beam interaction at finite c (much slower).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub beam_retarded: bool,
    /// Beam particles absorbed by their detector vanish at once (the charge drained
    /// instantly) instead of fading as they fly into the screening cup behind the
    /// detector's mouth (the default, PHYSICS.md §3.3).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub instant_drain: bool,
    pub t_max: f64,
    pub tolerances: TolerancesSpec,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct TolerancesSpec {
    pub preview: f64,
    pub verify: f64,
}

impl From<TolerancesSpec> for Tolerances {
    fn from(t: TolerancesSpec) -> Self {
        Tolerances {
            preview: t.preview,
            verify: t.verify,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParticleSpec {
    pub charge: f64,
    pub mass: f64,
    pub radius: f64,
    /// Magnetic moment along z (a spin state perpendicular to the plane); its energy is
    /// `−moment · B_z` (PHYSICS.md §3.2).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub moment: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Launch {
    pub node: Node,
    pub direction: [f64; 3],
    pub kinetic_energy: f64,
    /// Lab time of the launch. It matters only with time-dependent fields (antennas,
    /// waves): the particle sees them at `time + t`.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub time: f64,
}

/// Detector: the box spanned by two nodes, optionally with conditions on the arriving
/// particle (as instruments and stages of instruments have).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Detector {
    pub min: Node,
    pub max: Node,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acceptance: Option<DetectorAcceptance>,
}

/// Conditions on the particle entering a detector (PHYSICS.md §6.1).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct DetectorAcceptance {
    /// Direction of motion: `[axis, half-angle]`, in degrees (axis from +x).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub direction: Option<[f64; 2]>,
    /// Kinetic energy `[min, max]`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kinetic: Option<[f64; 2]>,
    /// Radiation of the whole flight into an arc of directions (detectors only, not
    /// gates; needs a finite `c`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radiation: Option<RadiationGoal>,
}

/// A radiation goal (PHYSICS.md §3.4): the energy per steradian that the flight radiates
/// into an arc of in-plane directions (averaged over the arc), in all frequencies or in a
/// band, must lie in a window. The receiver is far away (the far zone).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct RadiationGoal {
    /// Arc of directions: `[axis, half-width]`, in degrees (axis from +x).
    pub direction: [f64; 2],
    /// Angular-frequency band `[ω_min, ω_max]`; absent: all frequencies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub band: Option<[f64; 2]>,
    /// Energy per steradian `[min, max]`.
    pub energy: [f64; 2],
    /// The detector is a target that stops the particle abruptly, and the stop's radiation
    /// counts (needs a band: its spectrum is flat to infinite frequency).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub abrupt_stop: bool,
    /// A steady receiver (a scatterer driven by a wave, PHYSICS.md §3.4): the window
    /// `[τ₁, τ₂]` of the receiver's time (arrival less the light time from the origin)
    /// over which it measures the mean power per steradian, which `energy` then bounds.
    /// The particle need not arrive: it is measured when its flight ends after the
    /// window (the flight times out), and must not end before.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window: Option<[f64; 2]>,
}

impl RadiationGoal {
    pub fn to_physics(self) -> physics::spectrum::RadiationWindow {
        let a = self.direction[0].to_radians();
        physics::spectrum::RadiationWindow {
            axis: DVec3::new(libm::cos(a), libm::sin(a), 0.0),
            half_angle: self.direction[1].to_radians(),
            band: self.band.map(|[lo, hi]| (lo, hi)),
            energy: (self.energy[0], self.energy[1]),
            abrupt_stop: self.abrupt_stop,
            steady: self.window.map(|[t1, t2]| (t1, t2)),
        }
    }
}

impl DetectorAcceptance {
    pub fn to_physics(self) -> physics::trajectory::Acceptance {
        physics::trajectory::Acceptance {
            direction: self.direction.map(|[axis, half]| {
                let a = axis.to_radians();
                (
                    DVec3::new(libm::cos(a), libm::sin(a), 0.0),
                    half.to_radians(),
                )
            }),
            kinetic: self.kinetic.map(|[lo, hi]| (lo, hi)),
            radiation: self.radiation.map(RadiationGoal::to_physics),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ElementKind {
    /// Fixed charge; `value` is its charge Q.
    #[default]
    Charge,
    /// Magnet: a uniformly magnetized sphere with moment along +z; `value` is
    /// `μ = μ₀ m_z / 4π` (field units, PHYSICS.md §2.2). In the plane its field is
    /// `B_z = −μ / r³`.
    Magnet,
    /// Antenna: a small oscillating electric dipole in the plane (PHYSICS.md §2.4),
    /// `p(t) = value · (cos α, sin α, 0) · cos(ω t)` with `α = angle_deg` and `ω` the
    /// element's own `omega`, or the level's RF generator `physics.rf_omega` if it has
    /// none. All antennas start in phase at t = 0; a negative value is the opposite
    /// phase.
    Antenna,
    /// Plate electrode placed by the player (PHYSICS.md §2.7): a box of the level's
    /// `limits.plate` size centred on the node, its length along `angle_deg`, held at the
    /// potential `value` (0: grounded).
    Plate,
    /// Power supply of a tunable level electrode: sets the potential `value` of the
    /// electrode centred on this node.
    Supply,
    /// Free charge placed by the player (PHYSICS.md §3.3): a particle of the level's
    /// `limits.free_mass` and `limits.free_radius` with charge `value`, launched from the
    /// node towards `angle_deg` with the element's `speed`. It moves under every force
    /// and interacts with the other particles (the level flies as one beam).
    Free,
    /// An insulating cut the player makes in a tube level's metal: a straight slit from the
    /// node along `angle_deg`, `value` long (any length and angle within one metal; drawn
    /// from one node to another, `tube::Level::cut_between`), which its carriers cannot
    /// cross (laminations: the eddy currents' loops must close round it).
    Cut,
}

/// An element on a grid node.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Element {
    pub node: Node,
    #[serde(default)]
    pub kind: ElementKind,
    #[serde(alias = "charge")]
    pub value: f64,
    /// Orientation in the plane, degrees from +x (antennas only).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub angle_deg: f64,
    /// Own angular frequency (antennas only); `None`: the level's RF generator.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub omega: Option<f64>,
    /// Launch speed, cells per time unit (free charges only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub speed: Option<f64>,
}

impl Element {
    pub fn charge(node: Node, q: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Charge,
            value: q,
            angle_deg: 0.0,
            omega: None,
            speed: None,
        }
    }

    pub fn magnet(node: Node, mu: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Magnet,
            value: mu,
            angle_deg: 0.0,
            omega: None,
            speed: None,
        }
    }

    pub fn plate(node: Node, potential: f64, angle_deg: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Plate,
            value: potential,
            angle_deg,
            omega: None,
            speed: None,
        }
    }

    pub fn supply(node: Node, potential: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Supply,
            value: potential,
            angle_deg: 0.0,
            omega: None,
            speed: None,
        }
    }

    pub fn antenna(node: Node, p0: f64, angle_deg: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Antenna,
            value: p0,
            angle_deg,
            omega: None,
            speed: None,
        }
    }

    /// A free charge `q` launched towards `angle_deg` with `speed`.
    pub fn free(node: Node, q: f64, angle_deg: f64, speed: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Free,
            value: q,
            angle_deg,
            omega: None,
            speed: Some(speed),
        }
    }

    /// The same antenna at its own frequency `omega`.
    #[must_use]
    pub fn with_omega(mut self, omega: f64) -> Self {
        self.omega = Some(omega);
        self
    }
}

impl Coil {
    /// Whether its current is ramped.
    pub fn is_ramped(&self) -> bool {
        match self {
            Coil::Circle { rate, .. } | Coil::Polygon { rate, .. } => *rate != 0.0,
        }
    }

    /// Whether a circuit drives its current.
    pub fn is_driven(&self) -> bool {
        matches!(
            self,
            Coil::Circle { drive: Some(_), .. } | Coil::Polygon { drive: Some(_), .. }
        )
    }

    /// Whether its current changes in time (ramped or driven): time-dependent B and an
    /// induced E.
    pub fn is_time_dependent(&self) -> bool {
        self.is_ramped() || self.is_driven()
    }

    /// The node a tunable circular coil's power supply sits on: its terminal, the node on
    /// its wire due east of the centre (`(cx + round(r), cy)`), so that concentric coils
    /// (a betatron's) have supplies of their own. None for other coils.
    pub fn supply_node(&self) -> Option<Node> {
        match self {
            Coil::Circle {
                center,
                radius,
                tunable: true,
                ..
            } => {
                #[allow(clippy::cast_possible_truncation)]
                let dx = libm::round(*radius) as i64;
                Some([center[0] + dx, center[1], center[2]])
            }
            _ => None,
        }
    }
}

/// A coil placed by the level, lying in the plane, with strength `kappa = μ₀ I / 4π`
/// (current counter-clockwise seen from +z for positive `kappa`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "lowercase")]
pub enum Coil {
    Circle {
        center: Node,
        /// Radius in cells.
        radius: f64,
        kappa: f64,
        /// Ramp rate dκ/dt: the strength is `kappa + rate·t` (lab time), with the induced
        /// field −∂A/∂t (quasi-static, PHYSICS.md §2.2). 0: a steady current.
        #[serde(default, skip_serializing_if = "is_zero")]
        rate: f64,
        /// A circuit driving its current from `kappa` at t = 0 (PHYSICS.md §2.10; not with
        /// a ramp).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        drive: Option<Drive>,
        /// The player sets its ramp rate with a power supply (an element of kind `Supply`
        /// on its terminal, `Coil::supply_node`, on the slider of `limits.coil_rates`), or,
        /// if a circuit drives it, its source's level (a voltage: `limits.supply_voltages`);
        /// without one it keeps `rate`.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        tunable: bool,
    },
    Polygon {
        vertices: Vec<Node>,
        kappa: f64,
        #[serde(default, skip_serializing_if = "is_zero")]
        rate: f64,
        /// A circuit driving its current, as for a circle (self-inductance by the wire's
        /// surface flux, `inductance::polygon_self`; not with a ramp).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        drive: Option<Drive>,
    },
}

/// A circuit driving a level electrode or a coil (PHYSICS.md §2.10, one way): a
/// source through a series resistance into the element. An electrode's chain may have an
/// inductance in series (an LC with the electrode's capacitance) or its resistance as a
/// timed switch; a coil's may have a capacitance in series. The element starts from its
/// static state (an electrode at its bias potential, a coil at its strength κ); the
/// chain's own inductor and capacitor start empty.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Drive {
    pub source: Source,
    /// Series resistance (> 0).
    pub resistance: f64,
    /// Series inductance (electrodes only; 0: none).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub inductance: f64,
    /// Series capacitance (coils only; 0: none).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub capacitance: f64,
    /// The resistance is a switch toggling once (electrodes without an inductance).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub switch: Option<Switch>,
}

/// A timed switch: closed (or open) at first, toggling at lab time `at`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Switch {
    #[serde(default)]
    pub closed: bool,
    pub at: f64,
}

/// A source's voltage over lab time (PHYSICS.md §2.9).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum Source {
    /// Constant.
    Dc { value: f64 },
    /// `offset + amplitude sin(omega t + phase)`.
    Sine {
        #[serde(default, skip_serializing_if = "is_zero")]
        offset: f64,
        amplitude: f64,
        omega: f64,
        #[serde(default, skip_serializing_if = "is_zero")]
        phase: f64,
    },
    /// A trapezoid: `low` until `delay`, rising over `rise` to `high`, held for `width`,
    /// falling over `fall`; repeated every `period` if that is positive.
    Pulse {
        #[serde(default, skip_serializing_if = "is_zero")]
        low: f64,
        high: f64,
        delay: f64,
        rise: f64,
        width: f64,
        fall: f64,
        #[serde(default, skip_serializing_if = "is_zero")]
        period: f64,
    },
}

impl Source {
    pub fn waveform(&self) -> physics::circuit::Waveform {
        use physics::circuit::Waveform;
        match *self {
            Source::Dc { value } => Waveform::Dc(value),
            Source::Sine {
                offset,
                amplitude,
                omega,
                phase,
            } => Waveform::Sine {
                offset,
                amplitude,
                omega,
                phase,
            },
            Source::Pulse {
                low,
                high,
                delay,
                rise,
                width,
                fall,
                period,
            } => Waveform::Pulse {
                low,
                high,
                delay,
                rise,
                width,
                fall,
                period,
            },
        }
    }

    /// The source with its level set by a power supply: a constant's value, a sine's
    /// amplitude, a pulse's top.
    #[must_use]
    pub fn with_level(mut self, v: f64) -> Self {
        match &mut self {
            Source::Dc { value } => *value = v,
            Source::Sine { amplitude, .. } => *amplitude = v,
            Source::Pulse { high, .. } => *high = v,
        }
        self
    }

    /// Whether every parameter is finite and the times and the frequency are not negative.
    pub fn is_valid(&self) -> bool {
        match *self {
            Source::Dc { value } => value.is_finite(),
            Source::Sine {
                offset,
                amplitude,
                omega,
                phase,
            } => {
                [offset, amplitude, phase].iter().all(|v| v.is_finite())
                    && omega.is_finite()
                    && omega >= 0.0
            }
            Source::Pulse {
                low,
                high,
                delay,
                rise,
                width,
                fall,
                period,
            } => {
                low.is_finite()
                    && high.is_finite()
                    && [delay, rise, width, fall, period]
                        .iter()
                        .all(|v| v.is_finite() && *v >= 0.0)
            }
        }
    }
}

impl Drive {
    /// The chain of the physics (`physics::drive::Chain`) with this source.
    pub fn chain(&self, source: Source) -> physics::drive::Chain {
        physics::drive::Chain {
            wave: source.waveform(),
            r: self.resistance,
            l: self.inductance,
            c: self.capacitance,
            switch: self.switch.map(|s| (s.closed, s.at)),
        }
    }

    /// Whether the parameters are valid: a positive resistance, inductance and capacitance
    /// not negative, a valid source, a switch time not negative.
    pub fn is_valid(&self) -> bool {
        self.resistance.is_finite()
            && self.resistance > 0.0
            && self.inductance.is_finite()
            && self.inductance >= 0.0
            && self.capacitance.is_finite()
            && self.capacitance >= 0.0
            && self.source.is_valid()
            && self.switch.is_none_or(|s| s.at.is_finite() && s.at >= 0.0)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Limits {
    pub max_charges: u32,
    /// Allowed charge magnitudes |Q|.
    pub magnitudes: Vec<f64>,
    pub allow_positive: bool,
    pub allow_negative: bool,
    /// Maximum number of magnets the player may place.
    #[serde(default)]
    pub max_magnets: u32,
    /// Allowed magnet strengths |μ| (either orientation).
    #[serde(default)]
    pub magnet_strengths: Vec<f64>,
    /// Maximum number of antennas the player may place.
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_antennas: u32,
    /// Allowed antenna amplitudes |p₀| (either phase; orientations `ANTENNA_ANGLES`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub antenna_amplitudes: Vec<f64>,
    /// Frequencies ω the player may give an antenna. Empty: player antennas run at the
    /// level's RF generator frequency.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub antenna_omegas: Vec<f64>,
    /// Hardcore: instead of the listed values, any value between the smallest and the
    /// largest listed one (magnitudes, strengths, amplitudes, frequencies), and any
    /// antenna orientation.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub continuous: bool,
    /// If set, player elements may only be placed inside this box of nodes (inclusive),
    /// like the electrode region of a real instrument. Player plates must lie inside it
    /// entirely.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub region: Option<Region2>,
    /// Maximum number of plates (electrodes) the player may place.
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_plates: u32,
    /// Potentials a player plate may be held at (signed; 0 is grounded): set with a
    /// slider, anything in the range of the list (`Limits::allows_potential`); the listed
    /// values are what the solver tries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub plate_voltages: Vec<f64>,
    /// Size of the player's plates.
    #[serde(default, skip_serializing_if = "PlateSize::is_default")]
    pub plate: PlateSize,
    /// Potentials the power supply of a tunable level electrode may be set to (signed): a
    /// slider over the range of the list, like `plate_voltages`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub supply_voltages: Vec<f64>,
    /// Ramp rates dκ/dt the power supply of a tunable undriven coil may be set to
    /// (signed): a slider over the range of the list, like `supply_voltages` (a driven
    /// coil's supply sets its source's voltage, from `supply_voltages`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub coil_rates: Vec<f64>,
    /// Fields B along z the power supply of a tube level's solenoid may be set to
    /// (signed: positive out of the plane): a slider over the range of the list
    /// (`tube::TubeSpec::solenoid`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub solenoid_fields: Vec<f64>,
    /// Maximum number of free charges the player may place (dynamic particles).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_free: u32,
    /// Maximum number of cuts the player may make through a tube level's metals.
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_cuts: u32,
    /// Charges a free charge may carry (signed).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub free_charges: Vec<f64>,
    /// Launch speeds of a free charge (0: released at rest); directions `FREE_ANGLES`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub free_speeds: Vec<f64>,
    /// Mass of the player's free charges.
    #[serde(
        default = "default_free_mass",
        skip_serializing_if = "is_default_free_mass"
    )]
    pub free_mass: f64,
    /// Radius of the player's free charges (they collide as rigid spheres).
    #[serde(default = "default_radius", skip_serializing_if = "is_default_radius")]
    pub free_radius: f64,
}

fn default_free_mass() -> f64 {
    1.0
}

#[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
fn is_default_free_mass(v: &f64) -> bool {
    v.to_bits() == default_free_mass().to_bits()
}

#[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
fn is_default_radius(v: &f64) -> bool {
    v.to_bits() == default_radius().to_bits()
}

/// Launch directions a free charge may be given, degrees from +x (15° steps).
pub const FREE_ANGLES: [f64; 24] = [
    0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0,
    210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0,
];

/// A dynamic particle placed by the level: it moves and interacts like every particle
/// of the flight, and may have its own detector (a goal) or none (part of the scene).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct FreeParticle {
    pub particle: ParticleSpec,
    pub node: Node,
    /// Launch velocity in the plane, cells per time unit.
    #[serde(default)]
    pub velocity: [f64; 2],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detector: Option<Detector>,
}

/// Size of a player plate, in cells: in-plane length (along its angle), thickness across
/// it, and height (z, symmetric about the plane).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlateSize {
    pub length: f64,
    pub thickness: f64,
    pub height: f64,
}

impl Default for PlateSize {
    fn default() -> Self {
        Self {
            length: 4.0,
            thickness: 0.4,
            height: 4.0,
        }
    }
}

impl PlateSize {
    #[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
    fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

#[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
fn is_zero_u32(v: &u32) -> bool {
    *v == 0
}

/// Smallest and largest value of a list (the range of hardcore mode).
/// Whether `path` in a levels directory is a level: a `.json` file other than the
/// directory's metadata (`golden_hashes.json`, `curriculum.json`).
/// The order of level files: by their number prefix (`NN_slug.json`, numerically: 100
/// after 99, which a plain name order put before 11), then by name; files without a
/// number last.
pub fn level_order_key(path: &std::path::Path) -> (u64, String) {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let number = name
        .split_once('_')
        .and_then(|(n, _)| n.parse::<u64>().ok())
        .unwrap_or(u64::MAX);
    (number, name)
}

pub fn is_level_file(path: &std::path::Path) -> bool {
    path.extension().is_some_and(|e| e == "json")
        && !path.ends_with("golden_hashes.json")
        && !path.ends_with("curriculum.json")
}

/// A shipped level of this workspace by its slug (`levels/NN_<slug>.json`, whatever its
/// number in the curriculum), for tests and tools.
///
/// # Panics
/// If no such level is shipped, or it does not parse.
pub fn shipped(slug: &str) -> Level {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let path = std::fs::read_dir(&dir)
        .expect("levels directory")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .find(|p| {
            p.file_stem().and_then(|s| s.to_str()).is_some_and(|s| {
                s.split_once('_')
                    .is_some_and(|(n, rest)| rest == slug && n.chars().all(|c| c.is_ascii_digit()))
            })
        })
        .unwrap_or_else(|| panic!("no shipped level '{slug}'"));
    Level::from_json(&std::fs::read_to_string(&path).expect("level file"))
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

pub fn value_range(list: &[f64]) -> Option<(f64, f64)> {
    let lo = list.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = list.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (lo <= hi).then_some((lo, hi))
}

impl Limits {
    /// Whether `v` is an allowed value from `list`: exactly one of its entries, or in
    /// hardcore mode anything in its range.
    pub fn allows(&self, list: &[f64], v: f64) -> bool {
        if self.continuous {
            value_range(list).is_some_and(|(lo, hi)| (lo..=hi).contains(&v))
        } else {
            list.iter().any(|m| m.to_bits() == v.to_bits())
        }
    }

    /// Whether `v` is a potential the player may set with a slider (power supplies and
    /// plates): anything in the range of `list`, in every mode. The owner (2026-10-06):
    /// voltages are tuned by sliders, over a range symmetric about 0 whose ends do not
    /// give the answer away (checked for the shipped levels by `slider_ranges_hide_the_
    /// answer`).
    pub fn allows_potential(list: &[f64], v: f64) -> bool {
        v.is_finite() && value_range(list).is_some_and(|(lo, hi)| (lo..=hi).contains(&v))
    }
}

/// What a power supply operates (`Level::supply_targets`): an electrode or a coil, by
/// its index in the level's list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SupplyTarget {
    Electrode(usize),
    Coil(usize),
    /// A tube level's solenoid (B along z).
    Solenoid,
}

/// A box of grid nodes, inclusive.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Region2 {
    pub min: Node,
    pub max: Node,
}

impl Region2 {
    pub fn contains(&self, n: Node) -> bool {
        (0..3).all(|i| (self.min[i]..=self.max[i]).contains(&n[i]))
    }
}

/// Why a player placement is not allowed.
#[derive(Clone, Debug, PartialEq)]
pub enum PlacementError {
    TooManyCharges,
    TooManyMagnets,
    TooManyAntennas,
    TooManyPlates,
    TooManyFree,
    TooManyCuts,
    /// A cut without a length, or not within one of a tube level's metals.
    CutNotInMetal(Node),
    /// A free charge's launch speed that the level does not offer.
    SpeedNotAllowed(Option<f64>),
    /// A power supply that is not on the centre of a tunable electrode.
    NoTunableElectrode(Node),
    OutsideGrid(Node),
    OutsideRegion(Node),
    NotInPlane(Node),
    Occupied(Node),
    SignNotAllowed(f64),
    MagnitudeNotAllowed(f64),
    AngleNotAllowed(f64),
    FrequencyNotAllowed(Option<f64>),
}

impl Level {
    pub fn from_json(s: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(s)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).expect("level serialization cannot fail")
    }

    /// The JSON of everything that can change a run: the level without its texts (name,
    /// description, reality, the electrodes' and metals' labels). The run cache and the
    /// test pass records key on it, so that rewording a level reruns nothing.
    #[must_use]
    pub fn physics_json(&self) -> String {
        let mut bare = Level {
            name: String::new(),
            description: String::new(),
            reality: Reality::default(),
            ..self.clone()
        };
        for e in &mut bare.electrodes {
            e.label.clear();
        }
        if let Some(t) = &mut bare.tube {
            for m in &mut t.metals {
                m.label.clear();
            }
        }
        bare.to_json()
    }

    /// Electrode `k` as the interface names it: its number, and its function if it has one
    /// ("Electrode 3 — gate").
    #[must_use]
    pub fn electrode_name(&self, k: usize) -> String {
        match self.electrodes.get(k).map(|e| e.label.as_str()) {
            Some(l) if !l.is_empty() => format!("Electrode {} — {l}", k + 1),
            _ => format!("Electrode {}", k + 1),
        }
    }

    /// Electrode `k` within a sentence: "electrode 2 (drain)", its function without the
    /// details in parentheses ("drain" of "drain (held at 10)").
    #[must_use]
    pub fn electrode_ref(&self, k: usize) -> String {
        let label = self.electrodes.get(k).map_or("", |e| e.label.as_str());
        let role = label.split(" (").next().unwrap_or("").trim();
        if role.is_empty() {
            format!("electrode {}", k + 1)
        } else {
            format!("electrode {} ({role})", k + 1)
        }
    }

    /// Metal `k` of a tube level as the interface names it, as `electrode_name`.
    #[must_use]
    pub fn metal_name(&self, k: usize) -> String {
        let label = self
            .tube
            .as_ref()
            .and_then(|t| t.metals.get(k))
            .map_or("", |m| m.label.as_str());
        if label.is_empty() {
            format!("Metal {}", k + 1)
        } else {
            format!("Metal {} — {label}", k + 1)
        }
    }

    pub fn c(&self) -> f64 {
        self.physics.c.unwrap_or(f64::INFINITY)
    }

    pub fn particle(&self, shot: usize) -> Particle {
        let p = self.shots[shot].particle;
        Particle {
            charge: p.charge,
            mass: p.mass,
            radius: p.radius,
            moment: p.moment,
        }
    }

    /// Initial momentum of a shot from its launch energy and direction.
    pub fn launch_momentum(&self, shot: usize) -> DVec3 {
        let s = &self.shots[shot];
        let d = s.launch.direction;
        Kinematics::new(s.particle.mass, self.c())
            .momentum_from_kinetic_energy(s.launch.kinetic_energy, DVec3::new(d[0], d[1], d[2]))
    }

    pub fn detector_region(&self, shot: usize) -> Region {
        self.box_region(self.shots[shot].detector)
    }

    /// The region of gate `i`.
    pub fn gate_region(&self, i: usize) -> Region {
        self.box_region(self.gates[i])
    }

    /// A box of nodes as a region (a slab through the plane in 2D).
    pub(crate) fn box_region(&self, d: Detector) -> Region {
        let a = self.grid.position(d.min);
        let b = self.grid.position(d.max);
        let mut min = a.min(b);
        let mut max = a.max(b);
        if self.grid.is_2d() {
            // The detector is a slab through the plane.
            min.z = -1.0;
            max.z = 1.0;
        }
        Region::Box(Aabb { min, max })
    }

    /// World bounds: the grid expanded by one cell on every side.
    pub fn bounds(&self) -> Aabb {
        let m = self.grid.position(self.grid.max_node());
        Aabb {
            min: DVec3::splat(-1.0),
            max: m + DVec3::ONE,
        }
    }

    /// Refines the grid by an integer factor, keeping every node in place.
    pub fn refine(&mut self, factor: u32, player: &mut [Element]) {
        assert!(factor >= 1);
        let f = i64::from(factor);
        let scale = |n: &mut Node| n.iter_mut().for_each(|v| *v *= f);
        self.grid.subdivision *= factor;
        for s in &mut self.shots {
            scale(&mut s.launch.node);
            scale(&mut s.detector.min);
            scale(&mut s.detector.max);
        }
        if let Some(r) = &mut self.limits.region {
            scale(&mut r.min);
            scale(&mut r.max);
        }
        for c in &mut self.coils {
            match c {
                Coil::Circle { center, .. } => scale(center),
                Coil::Polygon { vertices, .. } => vertices.iter_mut().for_each(scale),
            }
        }
        for c in &mut self.conductors {
            scale(&mut c.center);
        }
        for e in &mut self.electrodes {
            scale(&mut e.center);
        }
        for d in &mut self.dielectrics {
            scale(&mut d.center);
        }
        for d in &mut self.dielectric_spheres {
            scale(&mut d.center);
        }
        for f in &mut self.ferrites {
            scale(&mut f.center);
        }
        for e in self
            .elements
            .iter_mut()
            .chain(self.reference_solution.iter_mut())
            .chain(player.iter_mut())
        {
            scale(&mut e.node);
        }
    }

    /// Whether a shot's detector has a radiation goal (PHYSICS.md §3.4).
    pub fn has_radiation_goal(&self) -> bool {
        self.shots
            .iter()
            .any(|s| s.detector.acceptance.is_some_and(|a| a.radiation.is_some()))
    }
}

impl Electrode {
    /// A rod's radius, or a ring's outer radius (None for a box).
    pub fn outer_radius(&self) -> Option<f64> {
        match self.shape {
            ElectrodeShape::Box => None,
            ElectrodeShape::Rod => Some(self.length / 2.0),
            ElectrodeShape::Ring => Some((self.length + self.thickness) / 2.0),
        }
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod flight_tests;

#[cfg(test)]
mod model_tests;
