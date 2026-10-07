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
//!   (2D Coulomb, `2λ r̂/r`, softened within `softening`, much below the spacing).
//! - **Emission** (space-charge limited, by Gauss's law): each step, each segment of
//!   the cathode whose surface density has the carriers' sign emits that charge (its
//!   density times its length) as macroparticles, released at rest `EMISSION_OFFSET` out at
//!   points along it from a deterministic low-discrepancy sequence. Emitting the charge
//!   the surface holds drives the cathode's normal field to zero there, the limit.
//! - **Absorption**: a particle that crosses an electrode's boundary (any change of
//!   which side of it the particle is on) is collected by that electrode, which records
//!   the charge; the current is the collected charge per time.
//! - **Push**: leapfrog (kick-drift) with a fixed step, symplectic; the emission adds
//!   particles at the start of a step.

use glam::DVec3;

use crate::zinv::{Electrodes, Section, line_charge, segment_integrals};

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

/// The tube's setup.
#[derive(Debug)]
pub struct Tube {
    /// The electrodes at their potentials (the static solution).
    pub electrodes: Electrodes,
    /// The emitting electrode (the cathode).
    pub cathode: usize,
    /// Charge-to-mass ratio of the carriers (negative for electrons).
    pub charge_per_mass: f64,
    /// Charge per unit length of a macroparticle (its magnitude).
    pub weight: f64,
    /// Softening length of the particles' mutual force.
    pub softening: f64,
    /// Time step.
    pub dt: f64,
}

/// The particles and the bookkeeping of a running tube.
#[derive(Clone, Debug, Default)]
pub struct TubeState {
    pub t: f64,
    pub x: Vec<DVec3>,
    pub v: Vec<DVec3>,
    pub charge: Vec<f64>,
    /// For each particle, per electrode, whether it is inside the electrode's section.
    side: Vec<Vec<bool>>,
    /// Per cathode segment, charge owed but not yet emitted (less than a macroparticle).
    pending: Vec<f64>,
    /// Counter of the emission sequence.
    emitted: u64,
    /// Charge per unit length collected by each electrode so far.
    pub collected: Vec<f64>,
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
    }
}

impl Tube {
    pub fn state(&self) -> TubeState {
        let cathode_segments = self.cathode_segments().count();
        TubeState {
            pending: vec![0.0; cathode_segments],
            collected: vec![0.0; self.electrodes.electrodes.len()],
            ..TubeState::default()
        }
    }

    fn cathode_segments(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.electrodes.segments.len()).filter(|&j| self.owner(j) == self.cathode)
    }

    fn owner(&self, j: usize) -> usize {
        self.electrodes.owner_of(j)
    }

    /// The total surface density: the static solution plus the particles' induced one.
    pub fn surface_density(&self, s: &TubeState) -> Vec<f64> {
        let induced = self.electrodes.induced_by(|y| {
            s.x.iter()
                .zip(&s.charge)
                .map(|(p, l)| line_charge(*p, *l, y).0)
                .sum()
        });
        self.electrodes
            .sigma
            .iter()
            .zip(&induced)
            .map(|(a, b)| a + b)
            .collect()
    }

    /// The field at `x` of the surface charge `sigma` and of every particle but `skip`.
    fn field(&self, s: &TubeState, sigma: &[f64], x: DVec3, skip: Option<usize>) -> DVec3 {
        let mut e = DVec3::ZERO;
        for ((p, q), d) in self.electrodes.segments.iter().zip(sigma) {
            e += segment_integrals(*p, *q, x).1 * *d;
        }
        let eps2 = self.softening * self.softening;
        for (k, (p, l)) in s.x.iter().zip(&s.charge).enumerate() {
            if Some(k) == skip {
                continue;
            }
            let r = DVec3::new(x.x - p.x, x.y - p.y, 0.0);
            let r2 = r.x * r.x + r.y * r.y + eps2;
            e += r * (2.0 * l / r2);
        }
        e
    }

    /// One step: emission, push, absorption.
    pub fn step(&self, s: &mut TubeState) {
        let sign = self.charge_per_mass.signum();
        // Emission from the cathode under the space-charge limit.
        let sigma = self.surface_density(s);
        let segs: Vec<usize> = self.cathode_segments().collect();
        let n_el = self.electrodes.electrodes.len();
        for (k, &j) in segs.iter().enumerate() {
            let (a, b) = self.electrodes.segments[j];
            let charge = sigma[j] * (b - a).length();
            if charge * sign > 0.0 {
                s.pending[k] += charge.abs();
            }
            while s.pending[k] >= self.weight {
                s.pending[k] -= self.weight;
                // Low-discrepancy point along the segment (golden ratio), just outside.
                s.emitted += 1;
                #[allow(clippy::cast_precision_loss)]
                let f = (s.emitted as f64 * 0.618_033_988_749_894_9).fract();
                let t = (b - a).normalize();
                let n = DVec3::new(-t.y, t.x, 0.0);
                let mid = a + (b - a) * f;
                // Outward: the side away from the cathode's section.
                let out = if inside(
                    &self.electrodes.electrodes[self.cathode].section,
                    mid + n * 1e-6,
                ) {
                    -n
                } else {
                    n
                };
                let delta = EMISSION_OFFSET * (b - a).length();
                let x = mid + out * delta;
                s.x.push(x);
                s.v.push(DVec3::ZERO);
                s.charge.push(sign * self.weight);
                s.side.push(
                    self.electrodes
                        .electrodes
                        .iter()
                        .map(|e| inside(&e.section, x))
                        .collect(),
                );
            }
        }
        // Kick and drift.
        let sigma = self.surface_density(s);
        let accel: Vec<DVec3> = (0..s.x.len())
            .map(|k| self.field(s, &sigma, s.x[k], Some(k)) * self.charge_per_mass)
            .collect();
        for ((x, v), a) in s.x.iter_mut().zip(&mut s.v).zip(&accel) {
            *v += *a * self.dt;
            *x += *v * self.dt;
        }
        s.t += self.dt;
        // Absorption: a particle that crossed an electrode's boundary.
        let mut k = 0;
        while k < s.x.len() {
            let hit = (0..n_el)
                .find(|&e| inside(&self.electrodes.electrodes[e].section, s.x[k]) != s.side[k][e]);
            if let Some(e) = hit {
                s.collected[e] += s.charge[k];
                s.x.swap_remove(k);
                s.v.swap_remove(k);
                s.charge.swap_remove(k);
                s.side.swap_remove(k);
            } else {
                k += 1;
            }
        }
    }
}
