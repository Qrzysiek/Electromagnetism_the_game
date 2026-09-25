//! Level format: grid, charges, launch parameters, limits and reference solutions
//! (SPEC §3, §7.9).
//!
//! All positions are integer grid nodes. A node `n` lies at `n / subdivision` cell units.
//! Refining the grid by an integer factor multiplies every node and the subdivision, so
//! all existing nodes stay exactly where they were.

pub mod solve;

use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Aabb, Region, Sphere};
use physics::trajectory::Scenario;
use physics::verify::Tolerances;
use serde::{Deserialize, Serialize};

pub const FORMAT_VERSION: u32 = 1;
pub const ENGINE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Integer grid node `[x, y, z]`.
pub type Node = [i64; 3];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Level {
    pub format_version: u32,
    /// Physics engine version the level (and its reference solution) was made with.
    pub engine_version: String,
    pub name: String,
    #[serde(default)]
    pub description: String,
    pub grid: Grid,
    pub physics: WorldPhysics,
    pub particle: ParticleSpec,
    pub launch: Launch,
    pub detector: Detector,
    /// Charges placed by the level (obstacles).
    pub level_charges: Vec<Charge>,
    pub limits: Limits,
    /// A known solution (player charges), if any.
    #[serde(default)]
    pub reference_solution: Vec<Charge>,
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

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WorldPhysics {
    /// Speed of light in internal units; `None` means Newtonian mechanics.
    pub c: Option<f64>,
    /// Radius of every fixed charge, in cell units.
    pub charge_radius: f64,
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
}

/// Detector B: the box spanned by two nodes.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Detector {
    pub min: Node,
    pub max: Node,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Charge {
    pub node: Node,
    pub charge: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Limits {
    pub max_charges: u32,
    /// Allowed magnitudes |Q|.
    pub magnitudes: Vec<f64>,
    pub allow_positive: bool,
    pub allow_negative: bool,
    /// If set, player charges may only be placed inside this box of nodes (inclusive),
    /// like the electrode region of a real instrument.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub region: Option<Region2>,
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
    OutsideGrid(Node),
    OutsideRegion(Node),
    NotInPlane(Node),
    Occupied(Node),
    SignNotAllowed(f64),
    MagnitudeNotAllowed(f64),
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

    pub fn particle(&self) -> Particle {
        Particle {
            charge: self.particle.charge,
            mass: self.particle.mass,
            radius: self.particle.radius,
        }
    }

    /// Initial momentum from the launch energy and direction.
    pub fn launch_momentum(&self) -> DVec3 {
        let d = self.launch.direction;
        Kinematics::new(self.particle.mass, self.c())
            .momentum_from_kinetic_energy(self.launch.kinetic_energy, DVec3::new(d[0], d[1], d[2]))
    }

    pub fn detector_region(&self) -> Region {
        let a = self.grid.position(self.detector.min);
        let b = self.grid.position(self.detector.max);
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
    pub fn refine(&mut self, factor: u32, player: &mut [Charge]) {
        assert!(factor >= 1);
        let f = i64::from(factor);
        let scale = |n: &mut Node| n.iter_mut().for_each(|v| *v *= f);
        self.grid.subdivision *= factor;
        scale(&mut self.launch.node);
        scale(&mut self.detector.min);
        scale(&mut self.detector.max);
        if let Some(r) = &mut self.limits.region {
            scale(&mut r.min);
            scale(&mut r.max);
        }
        for c in self
            .level_charges
            .iter_mut()
            .chain(self.reference_solution.iter_mut())
            .chain(player.iter_mut())
        {
            scale(&mut c.node);
        }
    }

    /// Checks a player placement against the level's limits and occupied nodes.
    pub fn check_placement(&self, player: &[Charge]) -> Result<(), PlacementError> {
        if player.len() > self.limits.max_charges as usize {
            return Err(PlacementError::TooManyCharges);
        }
        let mut occupied: Vec<Node> = self.level_charges.iter().map(|c| c.node).collect();
        occupied.push(self.launch.node);
        for c in player {
            if !self.grid.contains(c.node) {
                return Err(PlacementError::OutsideGrid(c.node));
            }
            if self.limits.region.is_some_and(|r| !r.contains(c.node)) {
                return Err(PlacementError::OutsideRegion(c.node));
            }
            if self.grid.is_2d() && c.node[2] != 0 {
                return Err(PlacementError::NotInPlane(c.node));
            }
            if occupied.contains(&c.node) {
                return Err(PlacementError::Occupied(c.node));
            }
            occupied.push(c.node);
            let sign_ok = if c.charge > 0.0 {
                self.limits.allow_positive
            } else {
                self.limits.allow_negative
            };
            if !sign_ok || c.charge == 0.0 {
                return Err(PlacementError::SignNotAllowed(c.charge));
            }
            if !self.limits.magnitudes.contains(&c.charge.abs()) {
                return Err(PlacementError::MagnitudeNotAllowed(c.charge));
            }
        }
        Ok(())
    }

    /// The physical scenario for a given player placement: level charges first, then
    /// player charges, in the given order (the summation order of the field).
    pub fn scenario(&self, player: &[Charge]) -> Scenario<Coulomb> {
        let charges: Vec<FixedCharge> = self
            .level_charges
            .iter()
            .chain(player)
            .map(|c| FixedCharge {
                position: self.grid.position(c.node),
                charge: c.charge,
                radius: self.physics.charge_radius,
            })
            .collect();
        Scenario {
            obstacles: charges
                .iter()
                .map(|c| Sphere {
                    center: c.position,
                    radius: c.radius,
                })
                .collect(),
            field: Coulomb::new(&charges),
            particle: self.particle(),
            c: self.c(),
            x0: self.grid.position(self.launch.node),
            p0: self.launch_momentum(),
            detector: Some(self.detector_region()),
            bounds: Some(self.bounds()),
            t_max: self.physics.t_max,
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
                t_max: 100.0,
                tolerances: TolerancesSpec {
                    preview: 1e-10,
                    verify: 1e-12,
                },
            },
            particle: ParticleSpec {
                charge: 1.0,
                mass: 1.0,
                radius: 0.0,
            },
            launch: Launch {
                node: [0, 5, 0],
                direction: [1.0, 0.1, 0.0],
                kinetic_energy: 0.3,
            },
            detector: Detector {
                min: [19, 3, 0],
                max: [20, 7, 0],
            },
            level_charges: vec![Charge {
                node: [10, 5, 0],
                charge: 0.1 + 0.2, // not exactly representable in decimal
            }],
            limits: Limits {
                max_charges: 2,
                magnitudes: vec![1.0, 2.0],
                allow_positive: true,
                allow_negative: true,
                region: None,
            },
            reference_solution: vec![],
        }
    }

    #[test]
    fn json_round_trip_is_exact() {
        let l = sample_level();
        let back = Level::from_json(&l.to_json()).unwrap();
        assert_eq!(l, back);
        assert_eq!(
            l.level_charges[0].charge.to_bits(),
            back.level_charges[0].charge.to_bits()
        );
    }

    #[test]
    fn refinement_keeps_positions() {
        let mut l = sample_level();
        let before = l.scenario(&[]);
        let mut player = vec![Charge {
            node: [3, 4, 0],
            charge: 1.0,
        }];
        let p_before = l.grid.position(player[0].node);
        l.refine(4, &mut player);
        let after = l.scenario(&[]);
        assert_eq!(before.x0, after.x0);
        assert_eq!(before.obstacles, after.obstacles);
        assert_eq!(p_before, l.grid.position(player[0].node));
    }

    #[test]
    fn placement_limits() {
        let l = sample_level();
        let ok = Charge {
            node: [5, 5, 0],
            charge: -2.0,
        };
        assert_eq!(l.check_placement(&[ok]), Ok(()));
        let bad_mag = Charge {
            node: [5, 5, 0],
            charge: 3.0,
        };
        assert!(matches!(
            l.check_placement(&[bad_mag]),
            Err(PlacementError::MagnitudeNotAllowed(_))
        ));
        let occupied = Charge {
            node: [10, 5, 0],
            charge: 1.0,
        };
        assert!(matches!(
            l.check_placement(&[occupied]),
            Err(PlacementError::Occupied(_))
        ));
        assert!(matches!(
            l.check_placement(&[ok, ok, ok]),
            Err(PlacementError::TooManyCharges)
        ));
    }
}
