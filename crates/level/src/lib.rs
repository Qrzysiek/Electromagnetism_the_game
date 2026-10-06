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
pub mod solve;

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
    /// Stages of the instrument: boxes every flight must pass, in order, before its
    /// detector counts, each with optional conditions on the entering particle
    /// (PHYSICS.md §6.2).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub gates: Vec<Detector>,
}

/// A rectangular metal box standing on the plane (symmetric about it): a plate, slab or
/// wall electrode.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Electrode {
    /// Centre (grid node).
    pub center: Node,
    /// In-plane length along `angle_deg`, thickness across it, and height (z), in cells.
    pub length: f64,
    pub thickness: f64,
    pub height: f64,
    #[serde(default)]
    pub angle_deg: f64,
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
    gates: Vec<Detector>,
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
            gates: f.gates,
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
    /// Maximum number of free charges the player may place (dynamic particles).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_free: u32,
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

/// In-plane distance from `p` to a coil's wire (its centre line).
fn coil_wire_distance(level: &Level, coil: &Coil, p: DVec3) -> f64 {
    let p = DVec3::new(p.x, p.y, 0.0);
    match coil {
        Coil::Circle { center, radius, .. } => {
            let c = level.grid.position(*center);
            ((p - DVec3::new(c.x, c.y, 0.0)).length() - radius).abs()
        }
        Coil::Polygon { vertices, .. } => {
            let n = vertices.len();
            (0..n)
                .map(|i| {
                    let a = level.grid.position(vertices[i]);
                    let b = level.grid.position(vertices[(i + 1) % n]);
                    let (a, b) = (DVec3::new(a.x, a.y, 0.0), DVec3::new(b.x, b.y, 0.0));
                    let ab = b - a;
                    let len2 = ab.length_squared();
                    let s = if len2 > 0.0 {
                        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    (p - (a + ab * s)).length()
                })
                .fold(f64::INFINITY, f64::min)
        }
    }
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
        for e in self
            .elements
            .iter_mut()
            .chain(self.reference_solution.iter_mut())
            .chain(player.iter_mut())
        {
            scale(&mut e.node);
        }
    }

    /// Checks a player placement against the level's limits and occupied nodes.
    pub fn check_placement(&self, player: &[Element]) -> Result<(), PlacementError> {
        let count = |k: ElementKind| player.iter().filter(|e| e.kind == k).count();
        if count(ElementKind::Charge) > self.limits.max_charges as usize {
            return Err(PlacementError::TooManyCharges);
        }
        if count(ElementKind::Magnet) > self.limits.max_magnets as usize {
            return Err(PlacementError::TooManyMagnets);
        }
        if count(ElementKind::Antenna) > self.limits.max_antennas as usize {
            return Err(PlacementError::TooManyAntennas);
        }
        if count(ElementKind::Plate) > self.limits.max_plates as usize {
            return Err(PlacementError::TooManyPlates);
        }
        if count(ElementKind::Free) > self.limits.max_free as usize {
            return Err(PlacementError::TooManyFree);
        }
        self.check_supplies(player)?;
        // Plates first: point elements are checked against every electrode.
        let boxes = self.check_plates(player)?;
        let mut occupied: Vec<Node> = self.elements.iter().map(|c| c.node).collect();
        occupied.extend(self.shots.iter().map(|s| s.launch.node));
        occupied.extend(self.free_particles.iter().map(|f| f.node));
        for e in player {
            if matches!(e.kind, ElementKind::Plate | ElementKind::Supply) {
                continue;
            }
            if !self.grid.contains(e.node) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if self.limits.region.is_some_and(|r| !r.contains(e.node)) {
                return Err(PlacementError::OutsideRegion(e.node));
            }
            if self.grid.is_2d() && e.node[2] != 0 {
                return Err(PlacementError::NotInPlane(e.node));
            }
            if occupied.contains(&e.node) {
                return Err(PlacementError::Occupied(e.node));
            }
            // Not inside or touching a metal sphere.
            let reach = self
                .physics
                .charge_radius
                .max(self.physics.magnet_radius)
                .max(self.physics.antenna_radius);
            let p = self.grid.position(e.node);
            if physics::bem::Electrodes::shapes_only(boxes.clone())
                .contains(p, CONTACT_DISTANCE + reach)
            {
                return Err(PlacementError::Occupied(e.node));
            }
            if self.conductors.iter().any(|c| {
                (p - self.grid.position(c.center)).length() < c.radius + CONTACT_DISTANCE + reach
            }) {
                return Err(PlacementError::Occupied(e.node));
            }
            occupied.push(e.node);
            if e.value == 0.0 && e.kind != ElementKind::Free {
                return Err(PlacementError::SignNotAllowed(e.value));
            }
            // Magnitudes come from the level's list and are compared exactly (hardcore:
            // anything within its range).
            let magnitude = e.value.abs();
            match e.kind {
                ElementKind::Charge => {
                    let sign_ok = if e.value > 0.0 {
                        self.limits.allow_positive
                    } else {
                        self.limits.allow_negative
                    };
                    if !sign_ok {
                        return Err(PlacementError::SignNotAllowed(e.value));
                    }
                    if !self.limits.allows(&self.limits.magnitudes, magnitude) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                }
                ElementKind::Magnet => {
                    if !self.limits.allows(&self.limits.magnet_strengths, magnitude) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                }
                ElementKind::Antenna => {
                    if !self
                        .limits
                        .allows(&self.limits.antenna_amplitudes, magnitude)
                    {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                    let angle_ok = if self.limits.continuous {
                        e.angle_deg.is_finite()
                    } else {
                        ANTENNA_ANGLES
                            .iter()
                            .any(|a| a.to_bits() == e.angle_deg.to_bits())
                    };
                    if !angle_ok {
                        return Err(PlacementError::AngleNotAllowed(e.angle_deg));
                    }
                    let omegas = &self.limits.antenna_omegas;
                    let omega_ok = match e.omega {
                        None => omegas.is_empty(),
                        Some(w) => self.limits.allows(omegas, w),
                    };
                    if !omega_ok {
                        return Err(PlacementError::FrequencyNotAllowed(e.omega));
                    }
                }
                ElementKind::Free => {
                    if !self.limits.allows(&self.limits.free_charges, e.value) {
                        return Err(PlacementError::MagnitudeNotAllowed(e.value));
                    }
                    // A free charge's launch velocity is always continuous: any speed in the
                    // range of the level's speeds, any direction (the owner: dragging the
                    // arrow between listed values was too clunky). The solver still searches
                    // the listed speeds and FREE_ANGLES only.
                    let speed_ok = e.speed.is_some_and(|v| {
                        v.is_finite()
                            && value_range(&self.limits.free_speeds)
                                .is_some_and(|(lo, hi)| (lo..=hi).contains(&v))
                    });
                    if !speed_ok {
                        return Err(PlacementError::SpeedNotAllowed(e.speed));
                    }
                    let angle_ok = e.angle_deg.is_finite();
                    if !angle_ok {
                        return Err(PlacementError::AngleNotAllowed(e.angle_deg));
                    }
                }
                ElementKind::Plate | ElementKind::Supply => unreachable!("checked above"),
            }
        }
        Ok(())
    }

    /// Geometry the physics cannot compute, which the level editor's free fields allow:
    /// dielectrics touching metal or each other (the boundary elements assume separate
    /// bodies), metal spheres overlapping or touching each other or an electrode (the image
    /// series and the boundary elements assume separate bodies; such setups failed or
    /// ran for minutes), and particles that start inside or touching metal or on a coil's
    /// wire (the field is singular there: NaN forces). Found by `scripts/pairs.py`. The
    /// game shows these instead of computing; the shipped levels have none.
    pub fn setup_issues(&self, player: &[Element]) -> Vec<String> {
        let mut out = Vec::new();
        let pos = |n: Node| self.grid.position(n);
        let boxes = physics::bem::Electrodes::shapes_only(self.all_box_electrodes(player));
        for (i, a) in self.conductors.iter().enumerate() {
            for (j, b) in self.conductors.iter().enumerate().skip(i + 1) {
                if (pos(a.center) - pos(b.center)).length() < a.radius + b.radius + CONTACT_DISTANCE
                {
                    out.push(format!(
                        "Metal spheres {} and {} touch or overlap.",
                        i + 1,
                        j + 1
                    ));
                }
            }
            if boxes.contains(pos(a.center), a.radius + CONTACT_DISTANCE) {
                out.push(format!("Metal sphere {} touches an electrode.", i + 1));
            }
        }
        // Dielectrics: apart from metal and from each other.
        let dielectrics: Vec<physics::bem::BoxElectrode> =
            self.box_dielectrics().iter().map(|d| d.shape).collect();
        let metal_boxes = self.all_box_electrodes(player);
        for (i, d) in dielectrics.iter().enumerate() {
            let corners = rect_corners(d);
            let touches = |o: &physics::bem::BoxElectrode| {
                rect_distance(&corners, &rect_corners(o)) < CONTACT_DISTANCE
            };
            if metal_boxes.iter().any(touches) {
                out.push(format!("Dielectric {} touches an electrode.", i + 1));
            }
            if dielectrics[i + 1..].iter().any(touches) {
                out.push(format!("Dielectric {} touches another dielectric.", i + 1));
            }
            let one = physics::bem::Electrodes::shapes_only(vec![*d]);
            if self
                .conductors
                .iter()
                .any(|c| one.contains(pos(c.center), c.radius + CONTACT_DISTANCE))
            {
                out.push(format!("Dielectric {} touches a metal sphere.", i + 1));
            }
        }
        let dielectric_shapes = physics::bem::Electrodes::shapes_only(dielectrics);
        let wire = self.physics.wire_radius + CONTACT_DISTANCE;
        let starts = self
            .shots
            .iter()
            .enumerate()
            .map(|(i, s)| (format!("Shot {}", i + 1), s.launch.node, 0.0))
            .chain(self.free_particles.iter().enumerate().map(|(i, f)| {
                (
                    format!("Free particle {}", i + 1),
                    f.node,
                    f.particle.radius,
                )
            }))
            .chain(
                player
                    .iter()
                    .filter(|e| e.kind == ElementKind::Free)
                    .enumerate()
                    .map(|(i, e)| {
                        (
                            format!("Your free charge {}", i + 1),
                            e.node,
                            self.limits.free_radius,
                        )
                    }),
            );
        for (who, node, radius) in starts {
            let p = pos(node);
            let in_sphere = self
                .conductors
                .iter()
                .any(|c| (p - pos(c.center)).length() < c.radius + radius + CONTACT_DISTANCE);
            if in_sphere || boxes.contains(p, radius + CONTACT_DISTANCE) {
                out.push(format!("{who} starts inside or touching metal."));
            }
            if dielectric_shapes.contains(p, radius + CONTACT_DISTANCE) {
                out.push(format!("{who} starts inside or touching a dielectric."));
            }
            if self
                .coils
                .iter()
                .any(|c| coil_wire_distance(self, c, p) < wire + radius)
            {
                out.push(format!("{who} starts on a coil's wire."));
            }
        }
        out
    }

    /// Power supplies: each on the centre of a different tunable electrode or coil, at a
    /// value on its slider.
    fn check_supplies(&self, player: &[Element]) -> Result<(), PlacementError> {
        let targets = self.supply_targets();
        let mut supplied: Vec<Node> = Vec::new();
        for e in player.iter().filter(|e| e.kind == ElementKind::Supply) {
            let Some(&(_, target)) = targets.iter().find(|(c, _)| *c == e.node) else {
                return Err(PlacementError::NoTunableElectrode(e.node));
            };
            if supplied.contains(&e.node) {
                return Err(PlacementError::Occupied(e.node));
            }
            supplied.push(e.node);
            if !Limits::allows_potential(self.supply_list(target), e.value) {
                return Err(PlacementError::MagnitudeNotAllowed(e.value));
            }
        }
        Ok(())
    }

    /// What the player's power supplies operate: the tunable electrodes, each at its
    /// centre, and the tunable circular coils, each at its terminal (`Coil::supply_node`).
    pub fn supply_targets(&self) -> Vec<(Node, SupplyTarget)> {
        let electrodes = self
            .electrodes
            .iter()
            .enumerate()
            .filter(|(_, e)| e.tunable)
            .map(|(i, e)| (e.center, SupplyTarget::Electrode(i)));
        let coils = self
            .coils
            .iter()
            .enumerate()
            .filter_map(|(i, c)| c.supply_node().map(|n| (n, SupplyTarget::Coil(i))));
        electrodes.chain(coils).collect()
    }

    /// The slider list of a supply's target: a voltage (an electrode's potential, or the
    /// source of a circuit driving a coil: a voltage source too) or a ramp rate (an
    /// undriven coil's).
    pub fn supply_list(&self, target: SupplyTarget) -> &[f64] {
        match target {
            SupplyTarget::Coil(i) if !self.coils[i].is_driven() => &self.limits.coil_rates,
            SupplyTarget::Electrode(_) | SupplyTarget::Coil(_) => &self.limits.supply_voltages,
        }
    }

    /// The player's power supply on the centre `node`, if placed.
    fn supply_at<'a>(&self, node: Node, player: &'a [Element]) -> Option<&'a Element> {
        player
            .iter()
            .find(|s| s.kind == ElementKind::Supply && s.node == node)
    }

    /// Player plates: allowed orientation and potential, entirely inside the grid and the
    /// player region, at least `PLATE_CLEARANCE` from every other electrode, clear of
    /// metal spheres, coil wires, elements, launch points and detectors. Returns every
    /// electrode box of the placement (level electrodes first).
    fn check_plates(
        &self,
        player: &[Element],
    ) -> Result<Vec<physics::bem::BoxElectrode>, PlacementError> {
        // Electrodes and dielectrics alike: plates keep clear of them, elements out.
        let mut boxes = self.box_electrodes();
        boxes.extend(self.box_dielectrics().iter().map(|d| d.shape));
        let reach = self
            .physics
            .charge_radius
            .max(self.physics.magnet_radius)
            .max(self.physics.antenna_radius);
        let grid_rect = self.node_rect(Region2 {
            min: [0, 0, 0],
            max: self.grid.max_node(),
        });
        for e in player.iter().filter(|e| e.kind == ElementKind::Plate) {
            if !self.grid.contains(e.node) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if self.grid.is_2d() && e.node[2] != 0 {
                return Err(PlacementError::NotInPlane(e.node));
            }
            let angle_ok = if self.limits.continuous {
                e.angle_deg.is_finite()
            } else {
                PLATE_ANGLES
                    .iter()
                    .any(|a| a.to_bits() == e.angle_deg.to_bits())
            };
            if !angle_ok {
                return Err(PlacementError::AngleNotAllowed(e.angle_deg));
            }
            if !Limits::allows_potential(&self.limits.plate_voltages, e.value) {
                return Err(PlacementError::MagnitudeNotAllowed(e.value));
            }
            let b = self.plate_box(e);
            let corners = rect_corners(&b);
            if !corners.iter().all(|c| rect_contains(&grid_rect, *c)) {
                return Err(PlacementError::OutsideGrid(e.node));
            }
            if let Some(r) = self.limits.region {
                let region = self.node_rect(r);
                if !corners.iter().all(|c| rect_contains(&region, *c)) {
                    return Err(PlacementError::OutsideRegion(e.node));
                }
            }
            let blocked = boxes
                .iter()
                .any(|o| rect_distance(&corners, &rect_corners(o)) < PLATE_CLEARANCE)
                || self.shots.iter().any(|s| {
                    rect_distance(
                        &corners,
                        &self.node_rect(Region2 {
                            min: s.detector.min,
                            max: s.detector.max,
                        }),
                    ) <= 0.0
                });
            let one = physics::bem::Electrodes::shapes_only(vec![b]);
            let points = self
                .elements
                .iter()
                .map(|x| (self.grid.position(x.node), reach))
                .chain(
                    self.shots
                        .iter()
                        .map(|s| (self.grid.position(s.launch.node), reach)),
                )
                .chain(self.coil_points().map(|p| (p, self.physics.wire_radius)))
                .chain(
                    self.conductors
                        .iter()
                        .map(|c| (self.grid.position(c.center), c.radius + PLATE_CLEARANCE)),
                )
                // Charge clouds: electrodes see them as point charges, so plates stay out.
                .chain(
                    self.clouds
                        .iter()
                        .map(|c| (self.grid.position(c.center), c.radius + CONTACT_DISTANCE)),
                );
            let touching = points
                .into_iter()
                .any(|(p, r)| one.contains(p, CONTACT_DISTANCE + r));
            if blocked || touching {
                return Err(PlacementError::Occupied(e.node));
            }
            boxes.push(b);
        }
        Ok(boxes)
    }

    /// Points along every coil wire, at most 0.1 cells apart (for clearance checks).
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    fn coil_points(&self) -> impl Iterator<Item = DVec3> + '_ {
        self.coils.iter().flat_map(move |c| {
            let mut out = Vec::new();
            match c {
                Coil::Circle { center, radius, .. } => {
                    let c = self.grid.position(*center);
                    let n = ((std::f64::consts::TAU * radius / 0.1).ceil() as usize).max(8);
                    for i in 0..n {
                        let a = std::f64::consts::TAU * i as f64 / n as f64;
                        out.push(c + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * *radius);
                    }
                }
                Coil::Polygon { vertices, .. } => {
                    for i in 0..vertices.len() {
                        let a = self.grid.position(vertices[i]);
                        let b = self.grid.position(vertices[(i + 1) % vertices.len()]);
                        let n = (((b - a).length() / 0.1).ceil() as usize).max(1);
                        for k in 0..n {
                            out.push(a + (b - a) * (k as f64 / n as f64));
                        }
                    }
                }
            }
            out
        })
    }

    /// The rectangle of a box of nodes, in cells (in-plane corners).
    fn node_rect(&self, r: Region2) -> [DVec3; 4] {
        let (a, b) = (self.grid.position(r.min), self.grid.position(r.max));
        [
            DVec3::new(a.x, a.y, 0.0),
            DVec3::new(b.x, a.y, 0.0),
            DVec3::new(b.x, b.y, 0.0),
            DVec3::new(a.x, b.y, 0.0),
        ]
    }

    /// The box of a player plate.
    pub fn plate_box(&self, e: &Element) -> physics::bem::BoxElectrode {
        let s = self.limits.plate;
        physics::bem::BoxElectrode {
            center: self.grid.position(e.node),
            angle: e.angle_deg.to_radians(),
            half_length: s.length / 2.0,
            half_thickness: s.thickness / 2.0,
            half_height: s.height / 2.0,
            bias: if e.value == 0.0 {
                Bias::Grounded
            } else {
                Bias::Potential(e.value)
            },
        }
    }

    /// Every electrode of a placement: the level's (at the potentials of the player's
    /// power supplies where set), then the player's plates.
    pub fn all_box_electrodes(&self, player: &[Element]) -> Vec<physics::bem::BoxElectrode> {
        let mut boxes = self.box_electrodes();
        for (b, e) in boxes.iter_mut().zip(&self.electrodes) {
            if let Some(s) = self.supply_of(e, player)
                && e.drive.is_none()
            {
                b.bias = Bias::Potential(s.value);
            }
        }
        boxes.extend(
            player
                .iter()
                .filter(|e| e.kind == ElementKind::Plate)
                .map(|e| self.plate_box(e)),
        );
        boxes
    }

    /// The player's power supply of a tunable electrode, if placed.
    fn supply_of<'a>(&self, e: &Electrode, player: &'a [Element]) -> Option<&'a Element> {
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
    fn drive_scales(&self, d: &physics::drive::Drives) -> (Vec<f64>, Vec<f64>) {
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

    /// Whether a placement has metal (spheres, electrodes or player plates), whose model
    /// differs between preview and verification resolution.
    pub fn has_metal(&self, player: &[Element]) -> bool {
        !self.conductors.is_empty()
            || !self.electrodes.is_empty()
            || !self.dielectrics.is_empty()
            || player.iter().any(|e| e.kind == ElementKind::Plate)
    }

    /// The static field and the obstacles for a placement at preview resolution.
    pub fn field(&self, player: &[Element]) -> (LevelField, Vec<Shape>) {
        self.field_at(player, Resolution::Preview)
    }

    /// Whether a shot's detector has a radiation goal (PHYSICS.md §3.4).
    pub fn has_radiation_goal(&self) -> bool {
        self.shots
            .iter()
            .any(|s| s.detector.acceptance.is_some_and(|a| a.radiation.is_some()))
    }

    /// Radiation goals (PHYSICS.md §3.4) need a finite `c` and sit on shot detectors only.
    /// Where several particles fly (free particles, free charges), a goal measures all of
    /// them together, their far fields added; a beam shot's particles arrive one by one, so
    /// a beam shot has none.
    fn radiation_goal_issues(&self, out: &mut Vec<String>) {
        let has = |d: &Detector| d.acceptance.is_some_and(|a| a.radiation.is_some());
        if self.shots.iter().any(|s| has(&s.detector)) {
            if self.physics.c.is_none() {
                out.push("a radiation goal needs a finite speed of light".into());
            }
            if self
                .shots
                .iter()
                .any(|s| s.beam.is_some() && has(&s.detector))
            {
                out.push(
                    "a beam shot cannot have a radiation goal (its particles arrive one by one)"
                        .into(),
                );
            }
            if self.has_beams()
                && self.shots.iter().any(|s| {
                    s.detector
                        .acceptance
                        .and_then(|a| a.radiation)
                        .is_some_and(|r| r.abrupt_stop)
                })
            {
                out.push(
                    "an abrupt stop is measured on single flights only (not with other particles \
                     flying)"
                        .into(),
                );
            }
        }
        if self.shots.iter().any(|s| {
            s.detector
                .acceptance
                .and_then(|a| a.radiation)
                .is_some_and(|r| r.abrupt_stop && r.band.is_none())
        }) {
            out.push(
                "a radiation goal with an abrupt stop needs a band (the stop's spectrum is flat to infinite frequency)"
                    .into(),
            );
        }
        for (i, s) in self.shots.iter().enumerate() {
            let Some(acc) = s.detector.acceptance else {
                continue;
            };
            let Some(r) = acc.radiation else { continue };
            let Some([t1, t2]) = r.window else { continue };
            if !(t1.is_finite() && t2.is_finite() && 0.0 <= t1 && t1 < t2) {
                out.push(format!(
                    "shot {}: a steady receiver's window must run forward from t ≥ 0",
                    i + 1
                ));
            } else if t2 >= self.physics.t_max {
                out.push(format!(
                    "shot {}: a steady receiver's window must end before the time limit (the \
                     light time from the scatterer to the origin is part of it)",
                    i + 1
                ));
            }
            if r.abrupt_stop {
                out.push(format!(
                    "shot {}: a steady receiver measures no stop (the particle need not arrive)",
                    i + 1
                ));
            }
            if acc.direction.is_some() || acc.kinetic.is_some() {
                out.push(format!(
                    "shot {}: a steady receiver's detector has no conditions on the arrival \
                     (the particle need not arrive)",
                    i + 1
                ));
            }
        }
        if self.gates.iter().any(has) {
            out.push("gates cannot have radiation goals".into());
        }
        if self
            .free_particles
            .iter()
            .any(|f| f.detector.as_ref().is_some_and(has))
        {
            out.push("free particles' detectors cannot have radiation goals".into());
        }
    }

    /// Problems of the physical model of this level (combinations that are not
    /// supported), as messages. Empty if the level is consistent.
    /// Positions and radii of the level's bodies and launch points for the containment
    /// checks of `model_issues`: fixed elements with the radius of their kind (supplies
    /// and plates have rules of their own), the shots' and the free particles' launch
    /// points with the particle's radius (contact during the flight uses it too).
    fn bodies(&self) -> impl Iterator<Item = (DVec3, f64)> + '_ {
        let elements = self.elements.iter().filter_map(|e| {
            let r = match e.kind {
                ElementKind::Charge => self.physics.charge_radius,
                ElementKind::Magnet => self.physics.magnet_radius,
                ElementKind::Antenna => self.physics.antenna_radius,
                ElementKind::Free => self.limits.free_radius,
                ElementKind::Plate | ElementKind::Supply => return None,
            };
            Some((self.grid.position(e.node), r))
        });
        let shots = self
            .shots
            .iter()
            .map(|s| (self.grid.position(s.launch.node), s.particle.radius));
        let free = self
            .free_particles
            .iter()
            .map(|f| (self.grid.position(f.node), f.particle.radius));
        elements.chain(shots).chain(free)
    }

    pub fn model_issues(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.drive_issues(&mut out);
        self.gate_issues(&mut out);
        self.radiation_goal_issues(&mut out);
        if self.has_beams() {
            self.beam_issues(&mut out);
        }
        // Magnetic moments: of the shots and of the free particles.
        let moments = self
            .shots
            .iter()
            .map(|s| s.particle.moment)
            .chain(self.free_particles.iter().map(|f| f.particle.moment))
            .any(|m| m != 0.0);
        if moments {
            // m ∇B_z is exact for E = 0, or for c = ∞ (PHYSICS.md §3.2). At finite c an
            // electric field adds velocity-dependent (Aharonov–Casher / hidden-momentum)
            // terms that are not modelled.
            let electric = self.elements.iter().any(|e| e.kind != ElementKind::Magnet)
                || self.limits.max_charges > 0
                || self.limits.max_antennas > 0
                || self.limits.max_plates > 0
                || self.limits.max_free > 0
                || !self.free_particles.is_empty()
                || !self.clouds.is_empty()
                || !self.conductors.is_empty()
                || !self.electrodes.is_empty()
                || self
                    .disturbances
                    .iter()
                    .any(|d| d.e != [0.0, 0.0] || !d.waves.is_empty());
            if self.physics.c.is_some() && electric {
                out.push(
                    "magnetic moments with electric fields at finite c need terms that are not \
                     modelled (use c = ∞, or magnetic fields only)"
                        .into(),
                );
            }
            if self.physics.radiation_reaction {
                out.push("radiation reaction is not modelled for magnetic moments".into());
            }
        }
        // A launch inside its own detector would arrive at t = 0.
        let in_detector = (0..self.shots.len()).any(|i| {
            self.detector_region(i)
                .signed_distance(self.grid.position(self.shots[i].launch.node))
                <= 0.0
        }) || self.free_particles.iter().any(|f| {
            f.detector.is_some_and(|d| {
                self.box_region(d)
                    .signed_distance(self.grid.position(f.node))
                    <= 0.0
            })
        });
        if in_detector {
            out.push("a launch point is inside its own detector".into());
        }
        if self.physics.radiation_reaction && self.physics.c.is_some() {
            // The metal responds electrostatically: its image follows the particle at
            // once and the radiation it would reflect is missing, which near its surface
            // changes the reaction by a factor of order one (PHYSICS.md §3.1).
            if self.has_metal(&[]) || self.limits.max_plates > 0 {
                out.push(
                    "radiation reaction near metal: the metal's electrostatic response does \
                     not reflect the radiation, which changes the reaction near it (not \
                     modelled)"
                        .into(),
                );
            }
            // A drained charge's field jumps where the light cone of its removal passes;
            // the reaction differentiates the field along the world line.
            if self.physics.instant_drain && self.interacts() {
                out.push(
                    "radiation reaction with an instant drain: the drained charge's field \
                     jumps where its light cone passes, which the reaction cannot follow"
                        .into(),
                );
            }
        }
        if self.limits.max_plates > 0 && self.limits.plate_voltages.is_empty() {
            out.push("player plates need at least one allowed potential".into());
        }
        let tunable: Vec<Node> = self
            .electrodes
            .iter()
            .filter(|e| e.tunable)
            .map(|e| e.center)
            .collect();
        if !tunable.is_empty() && self.limits.supply_voltages.is_empty() {
            out.push("tunable electrodes need at least one supply voltage".into());
        }
        let targets: Vec<Node> = self.supply_targets().iter().map(|t| t.0).collect();
        if (1..targets.len()).any(|i| targets[..i].contains(&targets[i])) {
            out.push("two tunable electrodes or coils share a centre node".into());
        }
        let tunable = |driven: bool| {
            self.coils
                .iter()
                .any(|c| matches!(c, Coil::Circle { tunable: true, .. }) && c.is_driven() == driven)
        };
        if tunable(false) && self.limits.coil_rates.is_empty() {
            out.push("tunable coils need at least one ramp rate".into());
        }
        if tunable(true) && self.limits.supply_voltages.is_empty() {
            out.push("tunable driven coils need at least one supply voltage".into());
        }
        if !self.electrodes.is_empty() || self.limits.max_plates > 0 {
            if !self.conductors.is_empty() {
                out.push("metal spheres and box electrodes are not yet solved together".into());
            }
            let time_dependent = self.limits.max_antennas > 0
                || self.elements.iter().any(|e| e.kind == ElementKind::Antenna)
                || !self.disturbances.is_empty();
            if time_dependent {
                out.push(
                    "electrodes respond electrostatically; antennas and disturbances would \
                     need a full-wave solution"
                        .into(),
                );
            }
            let boxes = physics::bem::Electrodes::shapes_only(self.box_electrodes());
            let inside = self
                .bodies()
                .any(|(x, r)| boxes.contains(x, CONTACT_DISTANCE + r));
            if inside {
                out.push("an element or launch point is inside or at an electrode".into());
            }
        }
        // Ramped and driven coils: their induced field is not conservative, which the metal
        // model (electrostatic) cannot screen; and time-dependent B is not combined with
        // magnetic moments (like antennas and waves).
        if self.coils.iter().any(Coil::is_time_dependent) {
            if !self.conductors.is_empty()
                || !self.electrodes.is_empty()
                || self.limits.max_plates > 0
            {
                out.push(
                    "ramped coils with metal: its electrostatic response cannot screen the \
                     induced field"
                        .into(),
                );
            }
            if moments {
                out.push("ramped coils with magnetic moments are not modelled".into());
            }
        }
        let c = self.c();
        let too_fast = self
            .free_particles
            .iter()
            .any(|f| libm::hypot(f.velocity[0], f.velocity[1]) >= c)
            || (self.limits.max_free > 0
                && self.limits.free_speeds.iter().any(|&v| v >= c || v < 0.0));
        if too_fast {
            out.push("free particles need speeds below c (and not negative)".into());
        }
        for cl in &self.clouds {
            let c = self.grid.position(cl.center);
            if cl.radius.is_nan() || cl.radius <= 0.0 || !cl.charge.is_finite() {
                out.push("a charge cloud needs a positive radius and a finite charge".into());
            }
            let touches_sphere = self.conductors.iter().any(|m| {
                (self.grid.position(m.center) - c).length()
                    <= m.radius + cl.radius + CONTACT_DISTANCE
            });
            let boxes = physics::bem::Electrodes::shapes_only(self.box_electrodes());
            if touches_sphere || boxes.contains(c, cl.radius + CONTACT_DISTANCE) {
                out.push(
                    "a charge cloud must stay clear of metal (the metal sees it as a point charge)"
                        .into(),
                );
            }
        }
        if !self.conductors.is_empty() {
            let time_dependent = self.limits.max_antennas > 0
                || self.elements.iter().any(|e| e.kind == ElementKind::Antenna)
                || !self.disturbances.is_empty();
            if time_dependent {
                out.push(
                    "metal spheres respond electrostatically; antennas and disturbances \
                     would need a full-wave solution"
                        .into(),
                );
            }
            for (i, a) in self.conductors.iter().enumerate() {
                for b in &self.conductors[i + 1..] {
                    let d = (self.grid.position(a.center) - self.grid.position(b.center)).length();
                    if d <= a.radius + b.radius + 2.0 * CONTACT_DISTANCE {
                        out.push("two metal spheres touch or overlap".into());
                    }
                }
                let c = self.grid.position(a.center);
                let inside = self
                    .bodies()
                    .any(|(x, r)| (x - c).length() < a.radius + CONTACT_DISTANCE + r);
                if inside {
                    out.push("an element or launch point is inside or at a metal sphere".into());
                }
            }
        }
        out
    }

    /// The field and the obstacles for a placement: level elements first, then player
    /// elements, in the given order (the summation order of the field). Metal spheres
    /// are computed at the given resolution (PHYSICS.md §2.6).
    pub fn field_at(&self, player: &[Element], resolution: Resolution) -> (LevelField, Vec<Shape>) {
        let all: Vec<&Element> = self.elements.iter().chain(player).collect();
        let charges: Vec<FixedCharge> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Charge)
            .map(|e| FixedCharge {
                position: self.grid.position(e.node),
                charge: e.value,
                radius: self.physics.charge_radius,
            })
            .collect();
        let dipoles: Vec<MagneticDipole> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Magnet)
            .map(|e| MagneticDipole {
                position: self.grid.position(e.node),
                moment: DVec3::new(0.0, 0.0, e.value),
                radius: self.physics.magnet_radius,
            })
            .collect();
        let antennas: Vec<OscillatingDipole> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Antenna)
            .map(|e| {
                let a = e.angle_deg.to_radians();
                OscillatingDipole {
                    position: self.grid.position(e.node),
                    amplitude: DVec3::new(libm::cos(a), libm::sin(a), 0.0) * e.value,
                    omega: e.omega.unwrap_or(self.physics.rf_omega),
                    phase: 0.0,
                    c: self.c(),
                    radius: self.physics.antenna_radius,
                }
            })
            .collect();
        let wire = self.physics.wire_radius;
        let mut loops = Vec::new();
        let mut polygons = Vec::new();
        let mut obstacles: Vec<Shape> = charges
            .iter()
            .map(|c| {
                Shape::Sphere(Sphere {
                    center: c.position,
                    radius: c.radius,
                })
            })
            .collect();
        obstacles.extend(dipoles.iter().map(|d| {
            Shape::Sphere(Sphere {
                center: d.position,
                radius: d.radius,
            })
        }));
        obstacles.extend(antennas.iter().map(|a| {
            Shape::Sphere(Sphere {
                center: a.position,
                radius: a.radius,
            })
        }));
        for coil in &self.coils {
            match coil {
                Coil::Circle {
                    center,
                    radius,
                    kappa,
                    rate,
                    drive,
                    tunable: _,
                } => {
                    // A tunable coil's supply sets its ramp rate (a driven one's: its
                    // source's level, `drives_for`).
                    let supplied = coil.supply_node().and_then(|n| self.supply_at(n, player));
                    let rate = match supplied {
                        Some(s) if drive.is_none() => s.value,
                        _ => *rate,
                    };
                    let l = CircularLoop {
                        center: self.grid.position(*center),
                        normal: DVec3::Z,
                        radius: *radius,
                        kappa: *kappa,
                        wire_radius: wire,
                        rate,
                    };
                    obstacles.push(Shape::Torus(Torus {
                        center: l.center,
                        normal: l.normal,
                        major: l.radius,
                        minor: wire,
                    }));
                    loops.push(l);
                }
                Coil::Polygon {
                    vertices,
                    kappa,
                    rate,
                    drive: _,
                } => {
                    let v: Vec<DVec3> = vertices.iter().map(|n| self.grid.position(*n)).collect();
                    for i in 0..v.len() {
                        obstacles.push(Shape::Capsule(Capsule {
                            a: v[i],
                            b: v[(i + 1) % v.len()],
                            radius: wire,
                        }));
                    }
                    polygons.push(PolygonCoil {
                        vertices: v,
                        kappa: *kappa,
                        wire_radius: wire,
                        rate: *rate,
                    });
                }
            }
        }
        // Electrode boxes (the level's, then the player's plates) with the contact shell.
        let boxes = self.all_box_electrodes(player);
        obstacles.extend(
            physics::bem::Electrodes::shapes_only(boxes.clone()).obstacles(CONTACT_DISTANCE),
        );
        let dielectrics = self.box_dielectrics();
        obstacles.extend(
            physics::bem::Electrodes::shapes_only(dielectrics.iter().map(|d| d.shape).collect())
                .obstacles(CONTACT_DISTANCE),
        );
        for c in &self.conductors {
            obstacles.push(Shape::Sphere(Sphere {
                center: self.grid.position(c.center),
                radius: c.radius + CONTACT_DISTANCE,
            }));
        }
        let clouds: Vec<ChargeCloud> = self
            .clouds
            .iter()
            .map(|c| ChargeCloud {
                position: self.grid.position(c.center),
                charge: c.charge,
                radius: c.radius,
            })
            .collect();
        // Metal and electrodes see the clouds as point charges: the validation keeps
        // every cloud clear of them, where its field is a point charge's.
        let sources: Vec<FixedCharge> = charges
            .iter()
            .copied()
            .chain(clouds.iter().map(|c| FixedCharge {
                position: c.position,
                charge: c.charge,
                radius: 0.0,
            }))
            .collect();
        let mut field = LevelField {
            coulomb: Coulomb::with_clouds(&charges, &clouds),
            dipoles,
            loops,
            polygons,
            antennas,
            external: Vec::new(),
            conductors: self.conductors_for(&sources, resolution),
            electrodes: electrodes_for(boxes, dielectrics, &sources, resolution),
            time_offset: 0.0,
            drives: None,
        };
        // A circuit that cannot be solved is a model issue (`drive_issues`): the field is
        // then without it.
        field.drives = self
            .drives_for(player, &field)
            .ok()
            .flatten()
            .map(std::sync::Arc::new);
        (field, obstacles)
    }

    /// Circuits (PHYSICS.md §2.10): valid parameters, chains that fit their element, a
    /// circuit that can be solved, and the flights' integration (the beam runner does not
    /// stop at the circuit's breakpoints and knots).
    fn drive_issues(&self, out: &mut Vec<String>) {
        if !self.has_drives() {
            return;
        }
        let before = out.len();
        if self.has_beams() {
            out.push("circuits are not combined with beams or free particles yet".into());
        }
        for e in &self.electrodes {
            let Some(d) = e.drive else { continue };
            if !d.is_valid() {
                out.push(
                    "a drive needs a positive resistance, valid source parameters and                      non-negative times, inductance and capacitance"
                        .into(),
                );
            }
            if matches!(e.bias, ConductorBias::Charge(_)) {
                out.push("a driven electrode starts at a potential, not floating".into());
            }
            if d.capacitance != 0.0 {
                out.push("an electrode's drive has no series capacitance (it is one)".into());
            }
            if d.switch.is_some() && d.inductance != 0.0 {
                out.push(
                    "a switch cannot open on an inductor's current (no inductance with a switch)"
                        .into(),
                );
            }
        }
        for c in &self.coils {
            let (Coil::Circle {
                drive: Some(d),
                rate,
                ..
            }
            | Coil::Polygon {
                drive: Some(d),
                rate,
                ..
            }) = c
            else {
                continue;
            };
            if !d.is_valid() {
                out.push(
                    "a drive needs a positive resistance, valid source parameters and                      non-negative times, inductance and capacitance"
                        .into(),
                );
            }
            if *rate != 0.0 {
                out.push("a driven coil has no ramp".into());
            }
            if d.inductance != 0.0 || d.switch.is_some() {
                out.push(
                    "a coil's drive has no series inductance (it is one) and no switch".into(),
                );
            }
        }
        if self.coils.iter().any(Coil::is_driven) && self.physics.c.is_none() {
            out.push("a driven coil needs a finite speed of light (its inductance is N/c²)".into());
        }
        if out.len() == before {
            let (field, _) = self.field_at(&self.reference_solution, Resolution::Preview);
            if let Err(e) = self.drives_for(&self.reference_solution, &field) {
                out.push(format!(
                    "the circuit cannot be solved ({e:?}: a loop of capacitors and sources, \
                     or a switch on an inductor's current)"
                ));
            }
        }
    }

    /// Gates: they must not overlap (each other) or contain a launch point.
    fn gate_issues(&self, out: &mut Vec<String>) {
        if self.gates.is_empty() {
            return;
        }
        let rect = |d: &Detector| {
            let (a, b) = (self.grid.position(d.min), self.grid.position(d.max));
            (a.min(b), a.max(b))
        };
        let overlap = |(a0, a1): (DVec3, DVec3), (b0, b1): (DVec3, DVec3)| {
            a0.x <= b1.x && b0.x <= a1.x && a0.y <= b1.y && b0.y <= a1.y
        };
        for (i, g) in self.gates.iter().enumerate() {
            if self.gates[i + 1..]
                .iter()
                .any(|h| overlap(rect(g), rect(h)))
            {
                out.push("two gates overlap".into());
            }
            let (lo, hi) = rect(g);
            if self.shots.iter().any(|s| {
                let p = self.grid.position(s.launch.node);
                (lo.x..=hi.x).contains(&p.x) && (lo.y..=hi.y).contains(&p.y)
            }) {
                out.push("a launch point lies in a gate".into());
            }
        }
    }

    /// Beams (PHYSICS.md §3.3): what the beam runner models.
    fn beam_issues(&self, out: &mut Vec<String>) {
        if self.interacts() {
            if self.has_metal(&[]) || self.limits.max_plates > 0 {
                out.push(
                    "interacting beams with metal: the charge one particle induces would act \
                     on the others (not modelled)"
                        .into(),
                );
            }
            // Interacting particles launched at the same point would have an infinite
            // interaction: a beam needs a position spread, and a single shot must not
            // share its launch point with other particles.
            let zero_spread = |s: &Shot| s.beam.is_none_or(|b| b.width == 0.0 && b.length == 0.0);
            let crowded = self.shots.iter().enumerate().any(|(i, a)| {
                zero_spread(a)
                    && (a.beam.is_some_and(|b| b.count > 1)
                        || self
                            .shots
                            .iter()
                            .enumerate()
                            .any(|(j, b)| j != i && b.launch.node == a.launch.node))
            });
            if crowded {
                out.push(
                    "interacting particles launched at the same point (give the beam a position spread)"
                        .into(),
                );
            }
            // Opposite point charges fall into each other (the classical point-charge
            // theory breaks down there); with a radius they collide as rigid spheres.
            let point_charges: Vec<f64> = self
                .shots
                .iter()
                .map(|s| s.particle)
                .chain(self.free_particles.iter().map(|f| f.particle))
                .filter(|p| p.radius == 0.0)
                .map(|p| p.charge)
                .chain(
                    if self.limits.max_free > 0 && self.limits.free_radius == 0.0 {
                        self.limits.free_charges.clone()
                    } else {
                        Vec::new()
                    },
                )
                .collect();
            let positive = point_charges.iter().any(|&q| q > 0.0);
            let negative = point_charges.iter().any(|&q| q < 0.0);
            if positive && negative {
                out.push(
                    "opposite point charges in one flight could fall into each other (give \
                     them a radius: they then collide as rigid spheres)"
                        .into(),
                );
            }
        }
        if self.shots.iter().any(|s| s.launch.time != 0.0) {
            out.push(
                "shots flying together (beams, or with dynamic particles) are launched at t = 0"
                    .into(),
            );
        }
        for s in &self.shots {
            if let Some(b) = s.beam {
                if b.count == 0 || !(b.transmission > 0.0 && b.transmission <= 1.0) {
                    out.push("a beam needs particles and a transmission in (0, 1]".into());
                }
                if b.energy_spread * b.reach() >= 1.0 || b.energy_spread < 0.0 {
                    out.push("a beam's energy spread would give non-positive energies".into());
                }
                if [b.angle_spread_deg, b.width, b.length]
                    .iter()
                    .any(|v| *v < 0.0)
                {
                    out.push("a beam's spreads must not be negative".into());
                }
            }
        }
    }

    /// The level's dielectric boxes for the physics.
    pub fn box_dielectrics(&self) -> Vec<physics::bem::BoxDielectric> {
        self.dielectrics
            .iter()
            .map(|d| physics::bem::BoxDielectric {
                shape: physics::bem::BoxElectrode {
                    center: self.grid.position(d.center),
                    angle: d.angle_deg.to_radians(),
                    half_length: d.length / 2.0,
                    half_thickness: d.thickness / 2.0,
                    half_height: d.height / 2.0,
                    bias: Bias::Charge(0.0),
                },
                permittivity: d.permittivity,
            })
            .collect()
    }

    /// The level's electrodes as physics boxes, at their own bias.
    pub fn box_electrodes(&self) -> Vec<physics::bem::BoxElectrode> {
        self.electrodes
            .iter()
            .map(|e| physics::bem::BoxElectrode {
                center: self.grid.position(e.center),
                angle: e.angle_deg.to_radians(),
                half_length: e.length / 2.0,
                half_thickness: e.thickness / 2.0,
                half_height: e.height / 2.0,
                bias: match e.bias {
                    ConductorBias::Grounded => Bias::Grounded,
                    ConductorBias::Charge(q) => Bias::Charge(q),
                    ConductorBias::Potential(v) => Bias::Potential(v),
                },
            })
            .collect()
    }
}

/// Largest image force the electrodes would exert on the particle (neglected, PHYSICS.md
/// §2.7) at the points `(x, t)` of a flight, bounded by `q²/d²` (four times the force of
/// a flat grounded plane at distance `d`, which covers concave corners), relative to
/// `max(|qE|, F₀)` with `F₀ = T₀` per cell, the smallest force that matters for the
/// flight (relative to the Lorentz force alone the ratio is meaningless where that force
/// passes through zero). 0 without electrodes.
pub fn electrode_image_force_bound(
    scn: &Scenario<LevelField>,
    kinetic_initial: f64,
    points: impl IntoIterator<Item = (DVec3, f64)>,
) -> f64 {
    use physics::field::FieldSolver;
    if scn.field.electrodes.is_empty() {
        return 0.0;
    }
    let boxes = scn.field.electrodes.obstacles(0.0);
    let q = scn.particle.charge;
    let mut worst: f64 = 0.0;
    for (x, t) in points {
        let d = boxes
            .iter()
            .map(|b| b.signed_distance(x))
            .fold(f64::INFINITY, f64::min)
            .max(1e-3);
        let f = (scn.field.sample(x, t).e.length() * q.abs()).max(kinetic_initial);
        worst = worst.max(q * q / (d * d) / f);
    }
    worst
}

/// Estimated energy radiated by the particle's magnetic moment `m` along the states
/// `(x, p, t)` of a flight (neglected, PHYSICS.md §3.2). A moving moment carries the
/// electric dipole `v×m/c²`, which radiates `2|ȧ×m|²/(3c⁷)` (`ȧ`: the jerk); the magnetic
/// quadrupole term is of the same order, so the estimate is `m²|ȧ|²/c⁷`, with the jerk
/// from the force between consecutive states. 0 for `c = ∞` or without a moment.
pub fn moment_radiation_estimate(
    scn: &Scenario<LevelField>,
    states: impl IntoIterator<Item = (DVec3, DVec3, f64)>,
) -> f64 {
    let (m, c) = (scn.particle.moment, scn.c);
    if m == 0.0 || !c.is_finite() {
        return 0.0;
    }
    let ode = physics::dynamics::ParticleOde::new(&scn.field, &scn.particle, c, 1.0);
    let accel =
        |x: DVec3, p: DVec3, t: f64| ode.force(x, p, t) / (ode.kin.gamma(p) * scn.particle.mass);
    let mut w = 0.0;
    let mut prev: Option<(DVec3, f64)> = None;
    for (x, p, t) in states {
        let a = accel(x, p, t);
        if let Some((a0, t0)) = prev
            && t > t0
        {
            let jerk = (a - a0) / (t - t0);
            w += m * m * jerk.length_squared() / c.powi(7) * (t - t0);
        }
        prev = Some((a, t));
    }
    w
}

/// Largest neglected image force that is still below the numerical accuracy (relative to
/// the force that matters, see `electrode_image_force_bound`).
pub const IMAGE_FORCE_LIMIT: f64 = 1e-10;

/// Box electrodes with the surface charge the fixed charges induce on them.
fn electrodes_for(
    boxes: Vec<physics::bem::BoxElectrode>,
    dielectrics: Vec<physics::bem::BoxDielectric>,
    charges: &[FixedCharge],
    resolution: Resolution,
) -> physics::bem::Electrodes {
    if boxes.is_empty() && dielectrics.is_empty() {
        return physics::bem::Electrodes::default();
    }
    let sources: Vec<(DVec3, f64)> = charges.iter().map(|c| (c.position, c.charge)).collect();
    let res = match resolution {
        Resolution::Preview => physics::bem::Resolution::Preview,
        Resolution::Verify => physics::bem::Resolution::Verify,
        Resolution::Display => physics::bem::Resolution::Display,
    };
    physics::bem::Electrodes::new_with_dielectrics(boxes, dielectrics, &sources, res)
}

/// In-plane corners of an electrode box (z = 0), counter-clockwise.
fn rect_corners(b: &physics::bem::BoxElectrode) -> [DVec3; 4] {
    let u = DVec3::new(libm::cos(b.angle), libm::sin(b.angle), 0.0) * b.half_length;
    let v = DVec3::new(-libm::sin(b.angle), libm::cos(b.angle), 0.0) * b.half_thickness;
    let c = DVec3::new(b.center.x, b.center.y, 0.0);
    [c - u - v, c + u - v, c + u + v, c - u + v]
}

/// Whether a point lies in a convex counter-clockwise polygon (boundary included).
fn rect_contains(r: &[DVec3; 4], p: DVec3) -> bool {
    (0..4).all(|i| {
        let (a, b) = (r[i], r[(i + 1) % 4]);
        (b - a).x * (p - a).y - (b - a).y * (p - a).x >= -1e-12
    })
}

fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

fn segments_cross(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> bool {
    let cross = |o: DVec3, p: DVec3, q: DVec3| (p - o).x * (q - o).y - (p - o).y * (q - o).x;
    let (d1, d2) = (cross(c, d, a), cross(c, d, b));
    let (d3, d4) = (cross(a, b, c), cross(a, b, d));
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

/// Distance between two convex quadrilaterals in the plane (0 if they overlap).
fn rect_distance(a: &[DVec3; 4], b: &[DVec3; 4]) -> f64 {
    let overlap = a.iter().any(|p| rect_contains(b, *p))
        || b.iter().any(|p| rect_contains(a, *p))
        || (0..4)
            .any(|i| (0..4).any(|j| segments_cross(a[i], a[(i + 1) % 4], b[j], b[(j + 1) % 4])));
    if overlap {
        return 0.0;
    }
    let one_way = |p: &[DVec3; 4], q: &[DVec3; 4]| {
        p.iter()
            .flat_map(|x| (0..4).map(move |j| segment_distance(*x, q[j], q[(j + 1) % 4])))
            .fold(f64::INFINITY, f64::min)
    };
    one_way(a, b).min(one_way(b, a))
}

impl Level {
    /// Metal spheres with the charges the fixed charges induce on them.
    fn conductors_for(&self, charges: &[FixedCharge], resolution: Resolution) -> Conductors {
        if self.conductors.is_empty() {
            return Conductors::default();
        }
        let spheres = self
            .conductors
            .iter()
            .map(|c| SphereConductor {
                center: self.grid.position(c.center),
                radius: c.radius,
                bias: match c.bias {
                    ConductorBias::Grounded => Bias::Grounded,
                    ConductorBias::Charge(q) => Bias::Charge(q),
                    ConductorBias::Potential(v) => Bias::Potential(v),
                },
            })
            .collect();
        let sources: Vec<(DVec3, f64)> = charges.iter().map(|c| (c.position, c.charge)).collect();
        Conductors::new(spheres, &sources, resolution)
    }

    /// One flight's scenario (shot, disturbance) with metal at display resolution (for
    /// pictures: potential map, field lines, field views).
    pub fn display_scenario(
        &self,
        shot: usize,
        disturbance: usize,
        player: &[Element],
    ) -> Scenario<LevelField> {
        let (field, obstacles) = self.field_at(player, Resolution::Display);
        self.scenario_with(shot, disturbance, field, obstacles)
    }

    /// Scenarios of all flights with the field at `resolution`.
    pub fn scenarios_at(
        &self,
        player: &[Element],
        resolution: Resolution,
    ) -> Vec<Scenario<LevelField>> {
        let (field, obstacles) = self.field_at(player, resolution);
        (0..self.flight_count())
            .map(|i| {
                let (shot, d) = self.flight_of(i);
                self.scenario_with(shot, d, field.clone(), obstacles.clone())
            })
            .collect()
    }

    /// Verifies every flight: preview at preview resolution, the tighter run with the
    /// field at verification resolution (they differ only with metal: spheres, electrodes,
    /// player plates), so the field model's error enters the verdict too.
    pub fn verify_flights(&self, player: &[Element]) -> Vec<Verification> {
        let preview = self.scenarios(player);
        let fine = if self.has_metal(player) {
            self.scenarios_at(player, Resolution::Verify)
        } else {
            preview.clone()
        };
        preview
            .iter()
            .zip(&fine)
            .map(|(a, b)| verify_pair(a, b, self.tolerances()))
            .collect()
    }

    /// Flights per shot: one per disturbance, or one undisturbed flight.
    pub fn flights_per_shot(&self) -> usize {
        self.disturbances.len().max(1)
    }

    /// Total number of flights, `shots × flights_per_shot`.
    pub fn flight_count(&self) -> usize {
        self.shots.len() * self.flights_per_shot()
    }

    /// Shot and disturbance of flight `i` (flights are ordered shot-major).
    pub fn flight_of(&self, i: usize) -> (usize, usize) {
        (i / self.flights_per_shot(), i % self.flights_per_shot())
    }

    /// The physical scenario of one shot under its first disturbance (or undisturbed).
    pub fn scenario(&self, shot: usize, player: &[Element]) -> Scenario<LevelField> {
        let (field, obstacles) = self.field(player);
        self.scenario_with(shot, 0, field, obstacles)
    }

    /// Scenarios of all flights, shot-major (see `flight_of`), sharing one field.
    pub fn scenarios(&self, player: &[Element]) -> Vec<Scenario<LevelField>> {
        let (field, obstacles) = self.field(player);
        (0..self.flight_count())
            .map(|i| {
                let (shot, d) = self.flight_of(i);
                self.scenario_with(shot, d, field.clone(), obstacles.clone())
            })
            .collect()
    }

    fn scenario_with(
        &self,
        shot: usize,
        disturbance: usize,
        mut field: LevelField,
        obstacles: Vec<Shape>,
    ) -> Scenario<LevelField> {
        if let Some(d) = self.disturbances.get(disturbance) {
            field.external = d.fields(self.c());
        }
        field.time_offset = self.shots[shot].launch.time;
        Scenario {
            field,
            obstacles,
            particle: self.particle(shot),
            c: self.c(),
            x0: self.grid.position(self.shots[shot].launch.node),
            p0: self.launch_momentum(shot),
            detector: Some(self.detector_region(shot)),
            bounds: Some(self.bounds()),
            t_max: self.physics.t_max,
            radiation_reaction: self.physics.radiation_reaction,
            acceptance: self.shots[shot]
                .detector
                .acceptance
                .map(DetectorAcceptance::to_physics),
            gates: (0..self.gates.len())
                .map(|i| physics::trajectory::Gate {
                    region: self.gate_region(i),
                    acceptance: self.gates[i].acceptance.map(DetectorAcceptance::to_physics),
                })
                .collect(),
        }
    }

    pub fn tolerances(&self) -> Tolerances {
        self.physics.tolerances.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_level() -> Level {
        Level {
            format_version: FORMAT_VERSION,
            engine_version: ENGINE_VERSION.to_string(),
            name: "test".into(),
            description: String::new(),
            grid: Grid {
                nx: 20,
                ny: 10,
                nz: 0,
                subdivision: 1,
            },
            physics: WorldPhysics {
                c: Some(5.0),
                charge_radius: 0.25,
                magnet_radius: 0.3,
                wire_radius: 0.1,
                antenna_radius: 0.3,
                rf_omega: 1.5,
                radiation_reaction: true,
                beam_interaction: false,
                beam_retarded: false,
                instant_drain: false,
                t_max: 100.0,
                tolerances: TolerancesSpec {
                    preview: 1e-10,
                    verify: 1e-12,
                },
            },
            shots: vec![Shot {
                particle: ParticleSpec {
                    charge: 1.0,
                    mass: 1.0,
                    radius: 0.0,
                    moment: 0.0,
                },
                launch: Launch {
                    node: [0, 5, 0],
                    direction: [1.0, 0.1, 0.0],
                    kinetic_energy: 0.3,
                    time: 0.25,
                },
                detector: Detector {
                    min: [19, 3, 0],
                    max: [20, 7, 0],
                    acceptance: Some(DetectorAcceptance {
                        direction: Some([0.0, 20.0]),
                        kinetic: Some([0.1, 1.0]),
                        radiation: None,
                    }),
                },
                beam: None,
            }],
            // 0.1 + 0.2 is not exactly representable in decimal.
            elements: vec![Element::charge([10, 5, 0], 0.1 + 0.2)],
            coils: vec![Coil::Circle {
                center: [10, 5, 0],
                radius: 8.0,
                kappa: 0.5,
                rate: 0.25,
                drive: None,
                tunable: false,
            }],
            limits: Limits {
                max_charges: 2,
                magnitudes: vec![1.0, 2.0],
                allow_positive: true,
                allow_negative: true,
                max_magnets: 1,
                magnet_strengths: vec![5.0],
                max_antennas: 1,
                antenna_amplitudes: vec![3.0],
                antenna_omegas: vec![0.5, 2.0],
                continuous: false,
                region: None,
                max_plates: 1,
                plate_voltages: vec![0.0, 5e4],
                plate: PlateSize {
                    length: 3.0,
                    thickness: 0.5,
                    height: 2.0,
                },
                supply_voltages: vec![-1e4, 1e4],
                coil_rates: vec![],
                max_free: 2,
                free_charges: vec![-1e-6, 2e-6],
                free_speeds: vec![0.0, 0.7],
                free_mass: 3.0,
                free_radius: 0.4,
            },
            reference_solution: vec![],
            disturbances: vec![Disturbance {
                name: "hum".into(),
                e: [0.01, -0.02],
                bz: 0.003,
                waves: vec![Wave {
                    amplitude: 0.1,
                    direction_deg: 30.0,
                    omega: 2.0,
                    phase_deg: 90.0,
                }],
            }],
            conductors: vec![],
            clouds: vec![Cloud {
                center: [18, 8, 0],
                radius: 1.0,
                charge: 1.5,
            }],
            free_particles: vec![FreeParticle {
                particle: ParticleSpec {
                    charge: 1e-6,
                    mass: 4.0,
                    radius: 0.3,
                    moment: 0.0,
                },
                node: [18, 1, 0],
                velocity: [0.1, -0.2],
                detector: Some(Detector {
                    min: [1, 1, 0],
                    max: [2, 3, 0],
                    acceptance: None,
                }),
            }],
            electrodes: vec![],
            dielectrics: vec![],
            gates: vec![],
        }
    }

    #[test]
    fn json_round_trip_is_exact() {
        let l = sample_level();
        let back = Level::from_json(&l.to_json()).unwrap();
        assert_eq!(l, back);
        assert_eq!(
            l.elements[0].value.to_bits(),
            back.elements[0].value.to_bits()
        );
    }

    #[test]
    fn version_1_files_are_migrated() {
        let v1 = r#"{
            "format_version": 1, "engine_version": "0.1.0", "name": "old",
            "grid": {"nx": 20, "ny": 10, "nz": 0, "subdivision": 1},
            "physics": {"c": 5.0, "charge_radius": 0.3, "t_max": 100.0,
                        "tolerances": {"preview": 1e-10, "verify": 1e-12}},
            "particle": {"charge": 1e-6, "mass": 1.0, "radius": 0.0},
            "launch": {"node": [0, 5, 0], "direction": [1.0, 0.0, 0.0], "kinetic_energy": 0.5},
            "detector": {"min": [18, 3, 0], "max": [20, 7, 0]},
            "level_charges": [{"node": [10, 5, 0], "charge": 2e6}],
            "limits": {"max_charges": 1, "magnitudes": [1e6], "allow_positive": true,
                       "allow_negative": true},
            "reference_solution": [{"node": [4, 4, 0], "charge": -1e6}]
        }"#;
        let l = Level::from_json(v1).unwrap();
        assert_eq!(l.format_version, FORMAT_VERSION);
        assert_eq!(l.shots.len(), 1);
        assert_eq!(l.shots[0].launch.node, [0, 5, 0]);
        assert_eq!(l.elements, vec![Element::charge([10, 5, 0], 2e6)]);
        assert_eq!(l.reference_solution, vec![Element::charge([4, 4, 0], -1e6)]);
        assert_eq!(l.limits.max_magnets, 0);
        assert!(l.disturbances.is_empty());
        assert_eq!(l.flight_count(), 1);
        // Saving writes version 2, which loads back identically.
        assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
    }

    /// A static level with one tunable electrode (at (5, 8), along x), room for two player
    /// plates and no coils.
    fn plate_level() -> Level {
        let mut l = sample_level();
        l.coils.clear();
        l.disturbances.clear();
        // Single flights: no dynamic particles (with them every shot flies at t = 0).
        l.free_particles.clear();
        l.limits.max_free = 0;
        l.limits.max_antennas = 0;
        // Metal does not reflect radiation: no radiation reaction next to it.
        l.physics.radiation_reaction = false;
        l.limits.max_plates = 2;
        l.limits.plate_voltages = vec![0.0, 1e4];
        l.limits.supply_voltages = vec![-2e4, 2e4];
        l.limits.plate = PlateSize::default();
        l.electrodes = vec![Electrode {
            center: [5, 8, 0],
            length: 4.0,
            thickness: 0.4,
            height: 4.0,
            angle_deg: 0.0,
            bias: ConductorBias::Grounded,
            tunable: true,
            drive: None,
        }];
        l
    }

    #[test]
    fn plate_and_supply_placement() {
        let l = plate_level();
        let plate = |x, y, v, a| Element::plate([x, y, 0], v, a);
        let err = |p: &[Element]| l.check_placement(p);
        assert_eq!(err(&[plate(14, 2, 1e4, 0.0)]), Ok(()));
        assert_eq!(err(&[plate(14, 2, 0.0, 90.0)]), Ok(()), "grounded plate");
        assert!(matches!(
            err(&[plate(14, 2, 1e4, 45.0)]),
            Err(PlacementError::AngleNotAllowed(_))
        ));
        // Potentials are sliders: anything in the list's range, nothing beyond it.
        assert_eq!(err(&[plate(14, 2, 5e3, 0.0)]), Ok(()), "unlisted, in range");
        assert!(matches!(
            err(&[plate(14, 2, 2e4, 0.0)]),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        // Sticking out of the grid (x from −1 to 3).
        assert!(matches!(
            err(&[plate(1, 2, 1e4, 0.0)]),
            Err(PlacementError::OutsideGrid(_))
        ));
        // On the launch point, on the level charge, across the detector.
        for p in [
            plate(2, 5, 1e4, 0.0),
            plate(10, 5, 1e4, 0.0),
            plate(18, 5, 1e4, 0.0),
        ] {
            assert!(
                matches!(err(&[p]), Err(PlacementError::Occupied(_))),
                "{p:?}"
            );
        }
        // Clearance: gap 0.6 < 1 between two plates or to the level electrode; 1.6 is fine.
        assert!(matches!(
            err(&[plate(14, 2, 1e4, 0.0), plate(14, 3, 1e4, 0.0)]),
            Err(PlacementError::Occupied(_))
        ));
        assert_eq!(
            err(&[plate(14, 2, 1e4, 0.0), plate(14, 4, 1e4, 0.0)]),
            Ok(())
        );
        assert!(matches!(
            err(&[plate(5, 7, 1e4, 0.0)]),
            Err(PlacementError::Occupied(_))
        ));
        assert_eq!(err(&[plate(5, 6, 1e4, 0.0)]), Ok(()));
        // Crossing plates overlap without any corner inside the other.
        assert!(matches!(
            err(&[plate(14, 4, 1e4, 0.0), plate(14, 4, 0.0, 90.0)]),
            Err(PlacementError::Occupied(_))
        ));
        // A charge on or in a player plate.
        for c in [[14, 2, 0], [15, 2, 0]] {
            assert!(matches!(
                err(&[plate(14, 2, 1e4, 0.0), Element::charge(c, 1.0)]),
                Err(PlacementError::Occupied(_))
            ));
        }
        assert!(matches!(
            err(&[
                plate(14, 2, 1e4, 0.0),
                plate(14, 4, 1e4, 0.0),
                plate(14, 6, 1e4, 0.0)
            ]),
            Err(PlacementError::TooManyPlates)
        ));
        // Power supplies: on the tunable electrode's centre, once, at a potential on the
        // slider (the list's range).
        assert_eq!(err(&[Element::supply([5, 8, 0], 2e4)]), Ok(()));
        assert_eq!(err(&[Element::supply([5, 8, 0], 1e4)]), Ok(()), "unlisted");
        assert!(matches!(
            err(&[Element::supply([6, 8, 0], 2e4)]),
            Err(PlacementError::NoTunableElectrode(_))
        ));
        assert!(matches!(
            err(&[Element::supply([5, 8, 0], 3e4)]),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        assert!(matches!(
            err(&[
                Element::supply([5, 8, 0], 2e4),
                Element::supply([5, 8, 0], -2e4)
            ]),
            Err(PlacementError::Occupied(_))
        ));
        assert_eq!(l.model_issues(), Vec::<String>::new());
    }

    /// A player plate is the same electrode as a level electrode at the same place and
    /// potential, and a power supply sets the potential of its electrode: identical
    /// fields and obstacles.
    #[test]
    fn plates_and_supplies_are_electrodes() {
        use physics::field::FieldSolver;
        let l = plate_level();
        let player = [
            Element::plate([14, 2, 0], 1e4, 90.0),
            Element::supply([5, 8, 0], -2e4),
        ];
        let mut fixed = l.clone();
        fixed.electrodes[0].tunable = false;
        fixed.electrodes[0].bias = ConductorBias::Potential(-2e4);
        fixed.electrodes.push(Electrode {
            center: [14, 2, 0],
            length: 4.0,
            thickness: 0.4,
            height: 4.0,
            angle_deg: 90.0,
            bias: ConductorBias::Potential(1e4),
            tunable: false,
            drive: None,
        });
        let (a, oa) = l.field(&player);
        let (b, ob) = fixed.field(&[]);
        assert_eq!(oa, ob);
        for x in [DVec3::new(3.0, 3.0, 0.0), DVec3::new(12.0, 6.0, 0.5)] {
            let (fa, fb) = (a.sample(x, 0.0), b.sample(x, 0.0));
            assert_eq!(fa.phi.to_bits(), fb.phi.to_bits());
            assert_eq!(fa.e, fb.e);
        }
        // The supply matters: without it the electrode is grounded.
        let (c, _) = l.field(&player[..1]);
        let x = DVec3::new(5.0, 6.0, 0.0);
        assert!((c.sample(x, 0.0).phi - a.sample(x, 0.0).phi).abs() > 1e3);
        assert!(l.has_metal(&[]) && !sample_level().has_metal(&[]));
        assert!(sample_level().has_metal(&player[..1]));
    }

    /// A driven electrode (PHYSICS.md §2.10): the JSON round trip keeps the drive; the
    /// circuit charges it from its bias towards the source's level, which a power supply
    /// on a tunable electrode sets; and a drive's issues are listed.
    #[test]
    fn driven_electrodes() {
        let mut l = plate_level();
        l.electrodes[0].drive = Some(Drive {
            source: Source::Dc { value: 1e4 },
            resistance: 0.05,
            inductance: 0.0,
            capacitance: 0.0,
            switch: None,
        });
        assert!(l.has_drives());
        assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
        let end = l.drive_end();
        let potential = |player: &[Element]| {
            let (f, _) = l.field(player);
            let d = f.drives.clone().expect("driven");
            let n = d
                .electrodes
                .iter()
                .find(|&&(e, _)| e == 0)
                .expect("in it")
                .1;
            (d.solution.potential(n, 0.0), d.solution.potential(n, end))
        };
        // From the grounded bias to the source's level (RC ≪ t_max).
        let (v0, v1) = potential(&[]);
        assert!(v0.abs() < 1e-6 && (v1 - 1e4).abs() < 1e-2, "{v0} {v1}");
        // A supply sets the level; the bias stays the start.
        let (v0, v1) = potential(&[Element::supply([5, 8, 0], -2e4)]);
        assert!(v0.abs() < 1e-6 && (v1 + 2e4).abs() < 1e-2, "{v0} {v1}");
        let drive_issue = |l: &Level, text: &str| l.model_issues().iter().any(|i| i.contains(text));
        assert!(!drive_issue(&l, "drive") && !drive_issue(&l, "circuit"));
        let bound = l.drive_back_action_bound(&[]);
        assert!(bound > 0.0 && bound.is_finite(), "{bound}");
        let mut bad = l.clone();
        bad.electrodes[0].drive.as_mut().expect("driven").resistance = 0.0;
        assert!(drive_issue(&bad, "positive resistance"));
        let mut switched = l.clone();
        {
            let d = switched.electrodes[0].drive.as_mut().expect("driven");
            d.inductance = 1.0;
            d.switch = Some(Switch {
                closed: false,
                at: 1.0,
            });
        }
        assert!(drive_issue(&switched, "switch cannot open"));
        let mut floating = l.clone();
        floating.electrodes[0].bias = ConductorBias::Charge(0.0);
        assert!(drive_issue(&floating, "not floating"));
    }

    /// A driven coil (DC through R): from its strength κ at t = 0 towards the source's
    /// current over c², as the RL circuit's closed form; it needs a finite c and no ramp.
    #[test]
    fn driven_coils() {
        let mut l = sample_level();
        l.free_particles.clear();
        l.limits.max_free = 0;
        let Coil::Circle { rate, drive, .. } = &mut l.coils[0] else {
            unreachable!("a circle")
        };
        *rate = 0.0;
        *drive = Some(Drive {
            source: Source::Dc { value: 1.0 },
            resistance: 1.0,
            inductance: 0.0,
            capacitance: 0.0,
            switch: None,
        });
        assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
        let c = l.physics.c.expect("finite c");
        let (f, _) = l.field(&[]);
        let d = f.drives.clone().expect("driven");
        let (_, k) = d.loops[0];
        let kappa = |t: f64| d.solution.current(k, t) * d.per_current;
        // κ∞ + (κ₀ − κ∞) e^{−t/τ}: κ₀ = 0.5, κ∞ = (V/R)/c², τ = L/R with L its
        // self-inductance over c².
        let tau = physics::inductance::ring_self(&f.loops[0]) / (c * c);
        let infinity = 1.0 / (c * c);
        let end = l.drive_end();
        for t in [0.0, 0.1 * end, 0.5 * end, end] {
            let want = infinity + (0.5 - infinity) * libm::exp(-t / tau);
            assert!((kappa(t) - want).abs() < 1e-9, "{t}: {} {want}", kappa(t));
        }
        let drive_issue = |l: &Level, text: &str| l.model_issues().iter().any(|i| i.contains(text));
        assert!(!drive_issue(&l, "drive") && !drive_issue(&l, "circuit"));
        let mut newtonian = l.clone();
        newtonian.physics.c = None;
        assert!(drive_issue(&newtonian, "finite speed of light"));
        let mut ramped = l.clone();
        if let Coil::Circle { rate, .. } = &mut ramped.coils[0] {
            *rate = 0.1;
        }
        assert!(drive_issue(&ramped, "no ramp"));
    }

    /// Coil supplies: a power supply on a tunable circular coil's terminal sets its ramp
    /// rate (a value on `coil_rates`' slider), or, if a circuit drives it, its source's
    /// level; not on an untunable coil, not beyond the slider.
    #[test]
    fn coil_supplies_set_the_ramp_or_the_source() {
        let mut l = sample_level();
        l.free_particles.clear();
        l.limits.max_free = 0;
        l.limits.coil_rates = vec![-0.5, 0.5];
        // The coil (centre (10, 5), radius 8): its terminal at (18, 5).
        let centre = [18, 5, 0];
        let supply = |v: f64| [Element::supply(centre, v)];
        // Not tunable yet: refused.
        assert!(matches!(
            l.check_placement(&supply(0.3)),
            Err(PlacementError::NoTunableElectrode(_))
        ));
        let Coil::Circle { tunable, .. } = &mut l.coils[0] else {
            unreachable!("a circle")
        };
        *tunable = true;
        assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
        assert_eq!(l.check_placement(&supply(0.3)), Ok(()));
        assert!(matches!(
            l.check_placement(&supply(0.7)),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        // The ramp: the supply's rate, else the coil's own (0.25).
        assert_eq!(
            l.field(&supply(0.3)).0.loops[0].rate.to_bits(),
            0.3f64.to_bits()
        );
        assert_eq!(l.field(&[]).0.loops[0].rate.to_bits(), 0.25f64.to_bits());
        // Driven: the supply sets the DC source's value, κ∞ = (V/R)/c².
        let Coil::Circle { rate, drive, .. } = &mut l.coils[0] else {
            unreachable!("a circle")
        };
        *rate = 0.0;
        *drive = Some(Drive {
            source: Source::Dc { value: 1.0 },
            resistance: 1.0,
            inductance: 0.0,
            capacitance: 0.0,
            switch: None,
        });
        // Its supply is now a voltage (`supply_voltages`, here ±1e4).
        assert_eq!(l.supply_targets()[0].0, centre);
        assert_eq!(
            l.supply_list(l.supply_targets()[0].1),
            &l.limits.supply_voltages[..]
        );
        let c = l.physics.c.expect("finite c");
        let end = l.drive_end();
        for (placement, v) in [(&supply(0.4)[..], 0.4), (&[][..], 1.0)] {
            let (f, _) = l.field(placement);
            let d = f.drives.clone().expect("driven");
            let (_, k) = d.loops[0];
            let tau = physics::inductance::ring_self(&f.loops[0]) / (c * c);
            let infinity = v / (c * c);
            let want = infinity + (0.5 - infinity) * libm::exp(-end / tau);
            let got = d.solution.current(k, end) * d.per_current;
            assert!((got - want).abs() < 1e-9, "{v}: {got} {want}");
        }
    }

    /// Refining keeps metal in place too (metal sphere and electrode centres are nodes),
    /// and power supplies stay on their electrodes.
    #[test]
    fn refinement_keeps_metal_in_place() {
        let mut l = plate_level();
        l.conductors.push(Conductor {
            center: [15, 7, 0],
            radius: 1.0,
            bias: ConductorBias::Grounded,
        });
        let mut player = vec![Element::supply([5, 8, 0], 2e4)];
        let before = l.scenario(0, &player).obstacles;
        l.refine(2, &mut player);
        assert_eq!(before, l.scenario(0, &player).obstacles);
        assert_eq!(l.check_placement(&player), Ok(()));
    }

    #[test]
    fn refinement_keeps_positions() {
        let mut l = sample_level();
        let before = l.scenario(0, &[]);
        let mut player = vec![Element::charge([3, 4, 0], 1.0)];
        let p_before = l.grid.position(player[0].node);
        l.refine(4, &mut player);
        let after = l.scenario(0, &[]);
        assert_eq!(before.x0, after.x0);
        assert_eq!(before.obstacles, after.obstacles);
        assert_eq!(p_before, l.grid.position(player[0].node));
    }

    #[test]
    fn placement_limits() {
        let l = sample_level();
        let ok = Element::charge([5, 5, 0], -2.0);
        assert_eq!(l.check_placement(&[ok]), Ok(()));
        assert!(matches!(
            l.check_placement(&[Element::charge([5, 5, 0], 3.0)]),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        assert!(matches!(
            l.check_placement(&[Element::charge([10, 5, 0], 1.0)]),
            Err(PlacementError::Occupied(_))
        ));
        assert!(matches!(
            l.check_placement(&[ok, ok, ok]),
            Err(PlacementError::TooManyCharges)
        ));
        let magnet = Element::magnet([6, 6, 0], -5.0);
        assert_eq!(l.check_placement(&[ok, magnet]), Ok(()));
        assert!(matches!(
            l.check_placement(&[magnet, Element::magnet([7, 6, 0], 5.0)]),
            Err(PlacementError::TooManyMagnets)
        ));
        let antenna = Element::antenna([8, 2, 0], -3.0, 135.0).with_omega(2.0);
        assert_eq!(l.check_placement(&[ok, magnet, antenna]), Ok(()));
        assert!(matches!(
            l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 30.0).with_omega(2.0)]),
            Err(PlacementError::AngleNotAllowed(_))
        ));
        // The level lists frequencies, so the generator frequency (None) or an
        // unlisted one is not allowed.
        assert!(matches!(
            l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 0.0)]),
            Err(PlacementError::FrequencyNotAllowed(None))
        ));
        assert!(matches!(
            l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 0.0).with_omega(1.0)]),
            Err(PlacementError::FrequencyNotAllowed(Some(_)))
        ));
        assert!(matches!(
            l.check_placement(&[
                antenna,
                Element::antenna([9, 2, 0], 3.0, 0.0).with_omega(0.5)
            ]),
            Err(PlacementError::TooManyAntennas)
        ));
        // Hardcore: anything within the ranges, any orientation.
        let mut hard = l.clone();
        hard.limits.continuous = true;
        let odd = Element::antenna([8, 2, 0], -3.0, 17.0).with_omega(1.3);
        assert!(l.check_placement(&[odd]).is_err());
        assert_eq!(
            hard.check_placement(&[odd, Element::charge([5, 5, 0], 1.7)]),
            Ok(())
        );
        assert!(matches!(
            hard.check_placement(&[Element::charge([5, 5, 0], 2.1)]),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        assert!(matches!(
            hard.check_placement(&[odd.with_omega(2.5)]),
            Err(PlacementError::FrequencyNotAllowed(_))
        ));
        // Each antenna oscillates at its own frequency.
        let (field, _) = l.field(&[antenna]);
        assert_eq!(field.antennas[0].omega.to_bits(), 2.0_f64.to_bits());
    }
}

#[cfg(test)]
mod flight_tests {
    use super::*;

    #[test]
    fn every_shot_flies_under_every_disturbance() {
        let mut l = shipped("twin_beams");
        assert_eq!(l.flight_count(), 2);
        assert!(l.scenarios(&[]).iter().all(|s| s.field.external.is_empty()));
        l.disturbances = vec![
            Disturbance::default(),
            Disturbance {
                e: [0.5, 0.0],
                ..Disturbance::default()
            },
            Disturbance {
                waves: vec![Wave {
                    amplitude: 1.0,
                    direction_deg: 0.0,
                    omega: 1.0,
                    phase_deg: 0.0,
                }],
                ..Disturbance::default()
            },
        ];
        let scn = l.scenarios(&[]);
        assert_eq!(scn.len(), 6);
        assert_eq!(l.flight_of(4), (1, 1));
        assert_eq!(scn[4].x0, l.scenario(1, &[]).x0);
        assert!(scn[3].field.external.is_empty());
        assert_eq!(scn[4].field.external.len(), 1);
        assert!(matches!(scn[5].field.external[0], External::Wave(_)));
    }
}

#[cfg(test)]
mod model_tests {
    use super::*;

    fn has(level: &Level, text: &str) -> bool {
        level.model_issues().iter().any(|i| i.contains(text))
    }

    /// A magnetic moment at finite c next to an electric source is rejected, whether the
    /// moment is a shot's or a free particle's, and whether the source is a fixed charge
    /// or a charge cloud (the Aharonov–Casher term is not modelled).
    #[test]
    fn moments_with_electric_sources_are_rejected() {
        let base = shipped("stern_gerlach");
        assert!(base.physics.c.is_some());
        assert!(base.model_issues().is_empty());
        let moment = "magnetic moments with electric fields";
        let mut cloud = base.clone();
        cloud.clouds.push(Cloud {
            center: [2, 2, 0],
            radius: 1.0,
            charge: 1.0,
        });
        assert!(has(&cloud, moment), "a cloud is an electric source");
        let mut free = base.clone();
        let spec = free.shots[0].particle;
        for s in &mut free.shots {
            s.particle.moment = 0.0;
        }
        free.elements.push(Element::charge([2, 2, 0], 1.0));
        assert!(!has(&free, moment));
        free.free_particles.push(FreeParticle {
            particle: spec,
            node: [3, 3, 0],
            velocity: [0.0, 0.0],
            detector: None,
        });
        assert!(has(&free, moment), "a free particle's moment counts");
    }

    /// Containment checks use each body's own radius (they once added 0.3 to every
    /// element), and a launch point inside its own detector is listed.
    #[test]
    fn containment_uses_the_bodies_radii() {
        let mut l = shipped("high_voltage_dome");
        assert!(l.model_issues().is_empty());
        let sphere = l.conductors[0];
        let c = l.grid.position(sphere.center);
        // A magnet whose surface is 0.1 cells from the sphere's.
        let gap = 0.1;
        l.physics.magnet_radius = 0.6;
        let d = sphere.radius + l.physics.magnet_radius + gap;
        let s = f64::from(l.grid.subdivision);
        #[allow(clippy::cast_possible_truncation)]
        let node = [((c.x + d) * s).round() as i64, (c.y * s).round() as i64, 0];
        let at = l.grid.position(node);
        assert!(((at - c).length() - d).abs() < 0.5 / s, "node rounding");
        l.elements.push(Element::magnet(node, 1.0));
        // Inside with its own radius 0.6, outside with the 0.3 once assumed for all.
        assert!((at - c).length() < sphere.radius + CONTACT_DISTANCE + 0.6);
        assert!((at - c).length() > sphere.radius + CONTACT_DISTANCE + 0.3);
        assert!(has(&l, "inside or at a metal sphere"));
        l.physics.magnet_radius = 0.3;
        assert!(!has(&l, "inside or at a metal sphere"));
        // Now a launch inside its detector.
        let mut m = shipped("first_bend");
        assert!(m.model_issues().is_empty());
        m.shots[0].launch.node = m.shots[0].detector.min;
        assert!(has(&m, "inside its own detector"));
    }

    /// Radiation reaction near metal is rejected (the metal's response is electrostatic).
    #[test]
    fn radiation_reaction_near_metal_is_rejected() {
        let mut l = shipped("image_charge");
        assert!(l.model_issues().is_empty());
        l.physics.c = Some(50.0);
        l.physics.radiation_reaction = true;
        assert!(has(&l, "radiation reaction near metal"));
    }
}
