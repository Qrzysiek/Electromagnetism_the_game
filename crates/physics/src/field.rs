//! Field sources (PHYSICS.md §2).

use std::sync::Arc;

use glam::DVec3;

use crate::antenna::OscillatingDipole;
use crate::bem::Electrodes;
use crate::conductor::Conductors;
use crate::drive::Drives;
use crate::external::External;
use crate::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};

/// Fields and potential at a point. Internal units: `k = 1/(4πε₀) = 1`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FieldSample {
    pub e: DVec3,
    pub b: DVec3,
    /// Electrostatic potential of the static sources (only meaningful for static electric
    /// fields; a `LevelField` adds nothing for its time-dependent sources).
    pub phi: f64,
}

/// Source of external fields acting on test particles.
pub trait FieldSolver: Sync {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample;

    /// True if the fields do not depend on time. Only then is `phi` a potential of `E`
    /// and the energy `(γ−1)mc² + qφ` conserved.
    fn is_static(&self) -> bool {
        true
    }

    /// Field and potential at `x` of the charges that a charge `q` at `x` induces in the
    /// scene's conductors (PHYSICS.md §2.6). The force on the particle is `q E`, its
    /// interaction energy `½ q φ`. Zero without conductors.
    fn self_field(&self, _x: DVec3, _q: f64) -> (DVec3, f64) {
        (DVec3::ZERO, 0.0)
    }

    /// Gradient of `B_z` at a point of the plane z = 0, for the force `m ∇B_z` on a
    /// magnetic moment `m ẑ` (PHYSICS.md §3.2). Covers the static magnetic sources
    /// (magnets, coils); zero for fields without them.
    fn grad_bz(&self, _x: DVec3, _t: f64) -> DVec3 {
        DVec3::ZERO
    }

    /// The times (of the flight, ascending) at which the fields change abruptly: a
    /// circuit's switch toggles and pulse corners (PHYSICS.md §2.10). The integration
    /// stops just before each and restarts on it, so that no step straddles one.
    fn breakpoints(&self) -> Vec<f64> {
        Vec::new()
    }

    /// The times (of the flight, ascending) at which the fields are continuous but their
    /// time derivatives are not: a circuit's integration steps (PHYSICS.md §2.10). The
    /// integration's steps end on them, so that each integrates smooth fields.
    fn knots(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl<F: FieldSolver + ?Sized> FieldSolver for &F {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        (**self).sample(x, t)
    }

    fn is_static(&self) -> bool {
        (**self).is_static()
    }

    fn self_field(&self, x: DVec3, q: f64) -> (DVec3, f64) {
        (**self).self_field(x, q)
    }

    fn grad_bz(&self, x: DVec3, t: f64) -> DVec3 {
        (**self).grad_bz(x, t)
    }

    fn breakpoints(&self) -> Vec<f64> {
        (**self).breakpoints()
    }

    fn knots(&self) -> Vec<f64> {
        (**self).knots()
    }
}

/// A fixed charge: a rigid sphere with a spherically symmetric charge distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FixedCharge {
    pub position: DVec3,
    pub charge: f64,
    pub radius: f64,
}

/// A charge cloud: a sphere of uniform charge density that particles can fly through
/// (J. J. Thomson's atom). Inside, `E = Q d / R³` and `φ = Q (3R² − |d|²) / (2R³)`, a
/// linear restoring field: a charge q of the opposite sign oscillates harmonically with
/// `ω₀² = |qQ| / (m R³)` (Jackson's bound charge, §16.7, Pr. 16.1). Outside it is a point
/// charge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChargeCloud {
    pub position: DVec3,
    pub charge: f64,
    pub radius: f64,
}

/// Exact Coulomb superposition of fixed charges, valid outside every rigid sphere (the
/// only region particles can reach), and of charge clouds, exact inside and outside.
/// Summed in the order stored: the fixed charges, then the clouds.
#[derive(Clone, Debug, Default)]
pub struct Coulomb {
    positions: Vec<DVec3>,
    charges: Vec<f64>,
    /// Radius of each source's uniform charge sphere if it is a cloud; 0 for a point
    /// charge (a rigid sphere, whose inside is never reached).
    radii: Vec<f64>,
}

impl Coulomb {
    pub fn new(charges: &[FixedCharge]) -> Self {
        Self::with_clouds(charges, &[])
    }

    /// Fixed charges and charge clouds.
    pub fn with_clouds(charges: &[FixedCharge], clouds: &[ChargeCloud]) -> Self {
        Self {
            positions: charges
                .iter()
                .map(|c| c.position)
                .chain(clouds.iter().map(|c| c.position))
                .collect(),
            charges: charges
                .iter()
                .map(|c| c.charge)
                .chain(clouds.iter().map(|c| c.charge))
                .collect(),
            radii: std::iter::repeat_n(0.0, charges.len())
                .chain(clouds.iter().map(|c| c.radius))
                .collect(),
        }
    }
}

impl Coulomb {
    /// Positions and charges, in summation order (clouds as point charges: their field
    /// outside themselves, which is all that metal outside them sees).
    pub fn charges(&self) -> impl Iterator<Item = (DVec3, f64)> + '_ {
        self.positions
            .iter()
            .copied()
            .zip(self.charges.iter().copied())
    }

    /// Positions, charges and radii of the clouds' uniform spheres (0 for a point
    /// charge), in summation order: what pictures of the potential need.
    pub fn sources(&self) -> impl Iterator<Item = (DVec3, f64, f64)> + '_ {
        self.charges()
            .zip(self.radii.iter().copied())
            .map(|((p, q), r)| (p, q, r))
    }
}

impl FieldSolver for Coulomb {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        let mut e = DVec3::ZERO;
        let mut phi = 0.0;
        for ((&r, &q), &radius) in self.positions.iter().zip(&self.charges).zip(&self.radii) {
            let d = x - r;
            let dist2 = d.length_squared();
            if dist2 < radius * radius {
                // Inside a cloud: the uniform sphere's field.
                let inv_r3 = 1.0 / (radius * radius * radius);
                phi += q * (3.0 * radius * radius - dist2) * 0.5 * inv_r3;
                e += d * (q * inv_r3);
                continue;
            }
            let inv_r = 1.0 / dist2.sqrt();
            let q_inv_r = q * inv_r;
            phi += q_inv_r;
            e += d * (q_inv_r * inv_r * inv_r);
        }
        FieldSample {
            e,
            b: DVec3::ZERO,
            phi,
        }
    }
}

/// One-off cost of a field's linear systems (metal spheres, electrodes): the time it took
/// to build them (they are cached per geometry afterwards, so this is the cost of a new
/// geometry), their matrix memory and their number of unknowns.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SetupCost {
    pub seconds: f64,
    pub bytes: usize,
    pub unknowns: usize,
}

impl std::ops::Add for SetupCost {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            seconds: self.seconds + o.seconds,
            bytes: self.bytes + o.bytes,
            unknowns: self.unknowns + o.unknowns,
        }
    }
}

/// All sources of a level: fixed charges, magnets, coils and antennas, plus external
/// fields (uniform stray fields and plane waves, PHYSICS.md §2.3–2.4). Summed in this
/// order.
#[derive(Clone, Debug, Default)]
pub struct LevelField {
    pub coulomb: Coulomb,
    pub dipoles: Vec<MagneticDipole>,
    pub loops: Vec<CircularLoop>,
    pub polygons: Vec<PolygonCoil>,
    pub antennas: Vec<OscillatingDipole>,
    pub external: Vec<External>,
    /// Conducting spheres and the charges the fixed sources induce on them.
    pub conductors: Conductors,
    /// Box electrodes (BEM) with their surface charge.
    pub electrodes: Electrodes,
    /// Added to the flight time before evaluating time-dependent sources: a particle
    /// launched at lab time `t₀` sees the fields at `t₀ + t`.
    pub time_offset: f64,
    /// The circuit driving electrodes and coils (`drive.rs`), solved in lab time; the
    /// driven coils' own `kappa` and `rate` are then not used.
    pub drives: Option<Arc<Drives>>,
}

impl LevelField {
    /// Build cost of the field's linear systems (metal spheres, electrodes, circuits).
    pub fn setup_cost(&self) -> SetupCost {
        self.conductors.setup_cost()
            + self.electrodes.setup_cost()
            + self
                .drives
                .as_ref()
                .map_or_else(SetupCost::default, |d| d.cost)
    }

    /// Magnetic field of the magnets and coils at t = 0.
    pub fn magnetic(&self, x: DVec3) -> DVec3 {
        self.magnetic_at(x, 0.0)
    }

    /// Magnetic field of the magnets and coils at lab time `t` (ramped and driven coils
    /// change).
    pub fn magnetic_at(&self, x: DVec3, t: f64) -> DVec3 {
        let mut b = DVec3::ZERO;
        for d in &self.dipoles {
            b += d.field(x);
        }
        for (i, l) in self.loops.iter().enumerate() {
            b += match self.driven_loop(i, t) {
                Some((kappa, _)) => CircularLoop { kappa, ..*l }.field(x),
                None => l.field_at(x, t),
            };
        }
        for (i, p) in self.polygons.iter().enumerate() {
            b += match self.driven_polygon(i, t) {
                Some((kappa, _)) => PolygonCoil {
                    kappa,
                    rate: 0.0,
                    ..p.clone()
                }
                .field(x),
                None => p.field_at(x, t),
            };
        }
        b
    }

    /// Whether a coil's current is ramped or driven (time-dependent B and an induced E).
    pub fn has_ramps(&self) -> bool {
        self.loops.iter().any(|l| l.rate != 0.0)
            || self.polygons.iter().any(|p| p.rate != 0.0)
            || self
                .drives
                .as_ref()
                .is_some_and(|d| !d.loops.is_empty() || !d.polygons.is_empty())
    }

    /// The strength and rate of loop `i` at lab time `t` if a circuit drives it.
    fn driven_loop(&self, i: usize, t: f64) -> Option<(f64, f64)> {
        self.drives.as_ref()?.loop_strength(i, t)
    }

    /// The strength and rate of polygonal coil `i` at lab time `t` if a circuit drives it.
    fn driven_polygon(&self, i: usize, t: f64) -> Option<(f64, f64)> {
        self.drives.as_ref()?.polygon_strength(i, t)
    }
}

impl FieldSolver for LevelField {
    fn sample(&self, x: DVec3, t: f64) -> FieldSample {
        let mut s = self.coulomb.sample(x, t);
        if !self.conductors.is_empty() {
            let c = self.conductors.induced.sample(x, t);
            s.e += c.e;
            s.phi += c.phi;
        }
        let t_lab = t + self.time_offset;
        if !self.electrodes.is_empty() {
            let c = match &self.drives {
                Some(d) if !d.electrodes.is_empty() => self
                    .electrodes
                    .sample_shifted(x, &d.electrode_shifts(t_lab)),
                _ => self.electrodes.sample(x, t),
            };
            s.e += c.e;
            s.phi += c.phi;
        }
        s.b = self.magnetic_at(x, t_lab);
        // Induced field of ramped and driven coils, −∂A/∂t (quasi-static, PHYSICS.md
        // §2.2).
        for (i, l) in self.loops.iter().enumerate() {
            s.e += match self.driven_loop(i, t_lab) {
                Some((_, rate)) => l.unit_vector_potential(x) * (-rate),
                None => l.induced_e(x),
            };
        }
        for (i, p) in self.polygons.iter().enumerate() {
            s.e += match self.driven_polygon(i, t_lab) {
                Some((_, rate)) => p.unit_vector_potential(x) * (-rate),
                None => p.induced_e(x),
            };
        }
        for a in &self.antennas {
            let f = a.fields(x, t_lab);
            s.e += f.e;
            s.b += f.b;
            // `phi` is the potential of the static sources (as for the waves): a static
            // antenna's electrostatic dipole potential; an oscillating one has none.
            if a.omega == 0.0 {
                s.phi += f.phi;
            }
        }
        for ext in &self.external {
            let f = ext.sample(x, t_lab);
            s.e += f.e;
            s.b += f.b;
            s.phi += f.phi;
        }
        s
    }

    fn is_static(&self) -> bool {
        !self.has_ramps()
            && self.drives.is_none()
            && self.external.iter().all(External::is_static)
            && self.antennas.iter().all(|a| a.omega == 0.0)
    }

    fn breakpoints(&self) -> Vec<f64> {
        self.drives.as_ref().map_or_else(Vec::new, |d| {
            d.breakpoints()
                .iter()
                .map(|&b| b - self.time_offset)
                .filter(|&b| b > 0.0)
                .collect()
        })
    }

    fn knots(&self) -> Vec<f64> {
        self.drives.as_ref().map_or_else(Vec::new, |d| {
            d.solution
                .step_ends()
                .into_iter()
                .map(|k| k - self.time_offset)
                .filter(|&k| k > 0.0)
                .collect()
        })
    }

    fn self_field(&self, x: DVec3, q: f64) -> (DVec3, f64) {
        self.conductors.self_field(x, q)
    }

    /// Magnets and coils (uniform stray B has no gradient). Time-dependent magnetic
    /// fields (antennas, waves) are not included: levels do not combine them with
    /// magnetic moments (`Level::model_issues`).
    fn grad_bz(&self, x: DVec3, t: f64) -> DVec3 {
        let mut g = DVec3::ZERO;
        for d in &self.dipoles {
            g += d.grad_bz(x);
        }
        let t_lab = t + self.time_offset;
        for (i, l) in self.loops.iter().enumerate() {
            g += if let Some((kappa, _)) = self.driven_loop(i, t_lab) {
                CircularLoop { kappa, ..*l }.grad_bz_in_plane(x)
            } else if l.rate == 0.0 {
                l.grad_bz_in_plane(x)
            } else {
                CircularLoop {
                    kappa: l.kappa_at(t_lab),
                    ..*l
                }
                .grad_bz_in_plane(x)
            };
        }
        for (i, p) in self.polygons.iter().enumerate() {
            g += if let Some((kappa, _)) = self.driven_polygon(i, t_lab) {
                PolygonCoil { kappa, ..p.clone() }.grad_bz_in_plane(x)
            } else if p.rate == 0.0 {
                p.grad_bz_in_plane(x)
            } else {
                PolygonCoil {
                    kappa: p.kappa_at(t_lab),
                    ..p.clone()
                }
                .grad_bz_in_plane(x)
            };
        }
        g
    }
}

/// Uniform static electric and magnetic fields, for tests with analytic solutions.
#[derive(Clone, Copy, Debug)]
pub struct UniformFields {
    pub e: DVec3,
    pub b: DVec3,
}

impl FieldSolver for UniformFields {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        FieldSample {
            e: self.e,
            b: self.b,
            phi: -self.e.dot(x),
        }
    }
}

/// Uniform static electric field, for tests with analytic solutions.
#[derive(Clone, Copy, Debug)]
pub struct UniformElectric {
    pub e: DVec3,
}

impl FieldSolver for UniformElectric {
    fn sample(&self, x: DVec3, _t: f64) -> FieldSample {
        FieldSample {
            e: self.e,
            b: DVec3::ZERO,
            phi: -self.e.dot(x),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coulomb_field_is_minus_gradient_of_potential() {
        let f = Coulomb::new(&[
            FixedCharge {
                position: DVec3::new(0.3, -1.0, 0.2),
                charge: 2.0,
                radius: 0.1,
            },
            FixedCharge {
                position: DVec3::new(-1.5, 0.5, -0.7),
                charge: -3.0,
                radius: 0.1,
            },
        ]);
        let x = DVec3::new(0.9, 0.4, -0.3);
        let h = 1e-5;
        let s = f.sample(x, 0.0);
        for axis in [DVec3::X, DVec3::Y, DVec3::Z] {
            let grad =
                (f.sample(x + axis * h, 0.0).phi - f.sample(x - axis * h, 0.0).phi) / (2.0 * h);
            assert!((s.e.dot(axis) + grad).abs() < 1e-8);
        }
    }
}
