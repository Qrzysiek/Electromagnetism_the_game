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
pub mod solve;

use physics::DVec3;
use physics::antenna::OscillatingDipole;
use physics::conductor::{Bias, Conductors, Resolution, SphereConductor};
use physics::dynamics::{Kinematics, Particle};
use physics::external::{External, PlaneWave};
use physics::field::{Coulomb, FixedCharge, LevelField};
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
    /// Box electrodes (plates, slabs, walls) placed by the level (PHYSICS.md §2.7).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub electrodes: Vec<Electrode>,
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
    electrodes: Vec<Electrode>,
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
            electrodes: f.electrodes,
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
}

impl DetectorAcceptance {
    pub fn to_physics(self) -> physics::trajectory::Acceptance {
        physics::trajectory::Acceptance {
            direction: self.direction.map(|[axis, half]| {
                let a = axis.to_radians();
                (DVec3::new(a.cos(), a.sin(), 0.0), half.to_radians())
            }),
            kinetic: self.kinetic.map(|[lo, hi]| (lo, hi)),
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
}

impl Element {
    pub fn charge(node: Node, q: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Charge,
            value: q,
            angle_deg: 0.0,
            omega: None,
        }
    }

    pub fn magnet(node: Node, mu: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Magnet,
            value: mu,
            angle_deg: 0.0,
            omega: None,
        }
    }

    pub fn antenna(node: Node, p0: f64, angle_deg: f64) -> Self {
        Self {
            node,
            kind: ElementKind::Antenna,
            value: p0,
            angle_deg,
            omega: None,
        }
    }

    /// The same antenna at its own frequency `omega`.
    #[must_use]
    pub fn with_omega(mut self, omega: f64) -> Self {
        self.omega = Some(omega);
        self
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
    },
    Polygon {
        vertices: Vec<Node>,
        kappa: f64,
    },
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
    /// like the electrode region of a real instrument.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub region: Option<Region2>,
}

#[allow(clippy::trivially_copy_pass_by_ref)] // serde's skip_serializing_if signature
fn is_zero_u32(v: &u32) -> bool {
    *v == 0
}

/// Smallest and largest value of a list (the range of hardcore mode).
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
        let d = self.shots[shot].detector;
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
        let mut occupied: Vec<Node> = self.elements.iter().map(|c| c.node).collect();
        occupied.extend(self.shots.iter().map(|s| s.launch.node));
        for e in player {
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
            if !self.electrodes.is_empty()
                && physics::bem::Electrodes::shapes_only(self.box_electrodes())
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
            if e.value == 0.0 {
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
            }
        }
        Ok(())
    }

    /// The static field and the obstacles for a placement at preview resolution.
    pub fn field(&self, player: &[Element]) -> (LevelField, Vec<Shape>) {
        self.field_at(player, Resolution::Preview)
    }

    /// Problems of the physical model of this level (combinations that are not
    /// supported), as messages. Empty if the level is consistent.
    pub fn model_issues(&self) -> Vec<String> {
        let mut out = Vec::new();
        if !self.electrodes.is_empty() {
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
            let inside =
                self.elements
                    .iter()
                    .any(|e| boxes.contains(self.grid.position(e.node), CONTACT_DISTANCE + 0.3))
                    || self.shots.iter().any(|s| {
                        boxes.contains(self.grid.position(s.launch.node), CONTACT_DISTANCE)
                    });
            if inside {
                out.push("an element or launch point is inside or at an electrode".into());
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
                let inside = self.elements.iter().any(|e| {
                    (self.grid.position(e.node) - c).length() < a.radius + CONTACT_DISTANCE + 0.3
                }) || self.shots.iter().any(|s| {
                    (self.grid.position(s.launch.node) - c).length() < a.radius + CONTACT_DISTANCE
                });
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
                    amplitude: DVec3::new(a.cos(), a.sin(), 0.0) * e.value,
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
                } => {
                    let l = CircularLoop {
                        center: self.grid.position(*center),
                        normal: DVec3::Z,
                        radius: *radius,
                        kappa: *kappa,
                        wire_radius: wire,
                    };
                    obstacles.push(Shape::Torus(Torus {
                        center: l.center,
                        normal: l.normal,
                        major: l.radius,
                        minor: wire,
                    }));
                    loops.push(l);
                }
                Coil::Polygon { vertices, kappa } => {
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
                    });
                }
            }
        }
        // Electrode boxes with the contact shell.
        obstacles.extend(
            physics::bem::Electrodes::shapes_only(self.box_electrodes())
                .obstacles(CONTACT_DISTANCE),
        );
        for c in &self.conductors {
            obstacles.push(Shape::Sphere(Sphere {
                center: self.grid.position(c.center),
                radius: c.radius + CONTACT_DISTANCE,
            }));
        }
        let field = LevelField {
            coulomb: Coulomb::new(&charges),
            dipoles,
            loops,
            polygons,
            antennas,
            external: Vec::new(),
            conductors: self.conductors_for(&charges, resolution),
            electrodes: self.electrodes_for(&charges, resolution),
            time_offset: 0.0,
        };
        (field, obstacles)
    }

    /// Box electrodes with the surface charge the fixed charges induce on them.
    fn electrodes_for(
        &self,
        charges: &[FixedCharge],
        resolution: Resolution,
    ) -> physics::bem::Electrodes {
        if self.electrodes.is_empty() {
            return physics::bem::Electrodes::default();
        }
        let sources: Vec<(DVec3, f64)> = charges.iter().map(|c| (c.position, c.charge)).collect();
        let res = match resolution {
            Resolution::Preview => physics::bem::Resolution::Preview,
            Resolution::Verify => physics::bem::Resolution::Verify,
            Resolution::Display => physics::bem::Resolution::Display,
        };
        physics::bem::Electrodes::new(self.box_electrodes(), &sources, res)
    }

    /// The electrodes as physics boxes.
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

    /// One shot's scenario with metal spheres at display resolution (for pictures:
    /// potential map, field lines, field views).
    pub fn display_scenario(&self, shot: usize, player: &[Element]) -> Scenario<LevelField> {
        let (field, obstacles) = self.field_at(player, Resolution::Display);
        self.scenario_with(shot, 0, field, obstacles)
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
    /// field at verification resolution (they differ only with metal spheres), so the
    /// field model's error enters the verdict too.
    pub fn verify_flights(&self, player: &[Element]) -> Vec<Verification> {
        let preview = self.scenarios(player);
        let fine = if self.conductors.is_empty() && self.electrodes.is_empty() {
            preview.clone()
        } else {
            self.scenarios_at(player, Resolution::Verify)
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
                    }),
                },
            }],
            // 0.1 + 0.2 is not exactly representable in decimal.
            elements: vec![Element::charge([10, 5, 0], 0.1 + 0.2)],
            coils: vec![Coil::Circle {
                center: [10, 5, 0],
                radius: 8.0,
                kappa: 0.5,
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
            electrodes: vec![],
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
        let mut l = Level::from_json(include_str!("../../../levels/06_twin_beams.json")).unwrap();
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
