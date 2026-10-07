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
//! - **Species**: several emitters (an electron cathode, an ion-emitting anode: the
//!   bipolar diode), and projectiles, single line charges launched at given times that
//!   interact with everything and are absorbed where they hit.
//! - **Push**: Boris with a fixed step (a uniform B along z rotates the velocity; with
//!   B = 0 it is the kick-drift leapfrog), symplectic; launches and emission add
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

/// Where a projectile ended.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Fate {
    /// Hit electrode `electrode` at time `t`, point `x`.
    Absorbed { electrode: usize, t: f64, x: DVec3 },
    /// Left the arena.
    Lost { t: f64, x: DVec3 },
}

/// Tag of projectile `j` (`TubeState::tag`); emitted carriers carry their emitter's
/// index.
pub const PROJECTILE_TAG: u32 = 1 << 16;

/// The tube's setup.
#[derive(Debug)]
pub struct Tube {
    /// The electrodes at their potentials (the static solution).
    pub electrodes: Electrodes,
    pub emitters: Vec<Emitter>,
    pub projectiles: Vec<Projectile>,
    /// Softening length of the particles' mutual force.
    pub softening: f64,
    /// Time step.
    pub dt: f64,
    /// The arena `(min, max)` corners: a particle that leaves it is lost (it flies off
    /// to infinity in the open tube; in the z-invariant world nothing brings it back
    /// unless an electrode encloses the arena, and then it never leaves). None: no bound.
    pub arena: Option<(DVec3, DVec3)>,
    /// A uniform magnetic field along z (a long solenoid around the tube: z-invariant).
    pub b_z: f64,
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
    /// For each particle, per electrode, whether it is inside the electrode's section.
    side: Vec<Vec<bool>>,
    /// Per emitter, per segment of its electrode, charge owed but not yet emitted (less
    /// than a macroparticle).
    pending: Vec<Vec<f64>>,
    /// Counter of the emission sequence.
    emitted: u64,
    /// Charge per unit length collected by each electrode so far, per emitter.
    pub collected: Vec<Vec<f64>>,
    /// Charge per unit length that left the arena.
    pub lost: f64,
    /// Which projectiles have been launched, and where each ended.
    launched: Vec<bool>,
    pub fates: Vec<Option<Fate>>,
}

impl TubeState {
    /// Where projectile `j` is (None before its launch and after its end).
    pub fn projectile(&self, j: usize) -> Option<DVec3> {
        let tag = PROJECTILE_TAG + u32::try_from(j).ok()?;
        self.tag.iter().position(|&t| t == tag).map(|k| self.x[k])
    }
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
    }
}

impl Tube {
    pub fn state(&self) -> TubeState {
        let n_el = self.electrodes.electrodes.len();
        TubeState {
            pending: self
                .emitters
                .iter()
                .map(|e| vec![0.0; self.segments_of(e.electrode).count()])
                .collect(),
            collected: vec![vec![0.0; self.emitters.len()]; n_el],
            launched: vec![false; self.projectiles.len()],
            fates: vec![None; self.projectiles.len()],
            ..TubeState::default()
        }
    }

    fn segments_of(&self, electrode: usize) -> impl Iterator<Item = usize> + '_ {
        (0..self.electrodes.segments.len())
            .filter(move |&j| self.electrodes.owner_of(j) == electrode)
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

    fn add(&self, s: &mut TubeState, x: DVec3, v: DVec3, charge: f64, qm: f64, tag: u32) {
        s.x.push(x);
        s.v.push(v);
        s.charge.push(charge);
        s.qm.push(qm);
        s.tag.push(tag);
        s.side.push(
            self.electrodes
                .electrodes
                .iter()
                .map(|e| inside(&e.section, x))
                .collect(),
        );
    }

    /// Space-charge-limited emission of every emitter: each segment of its electrode
    /// whose surface density has its carriers' sign releases that charge.
    fn emit(&self, s: &mut TubeState) {
        let sigma = self.surface_density(s);
        for (i, em) in self.emitters.iter().enumerate() {
            let sign = em.charge_per_mass.signum();
            let segs: Vec<usize> = self.segments_of(em.electrode).collect();
            for (k, &j) in segs.iter().enumerate() {
                let (a, b) = self.electrodes.segments[j];
                let charge = sigma[j] * (b - a).length();
                if charge * sign > 0.0 {
                    s.pending[i][k] += charge.abs();
                }
                while s.pending[i][k] >= em.weight {
                    s.pending[i][k] -= em.weight;
                    // Low-discrepancy point along the segment (golden ratio).
                    s.emitted += 1;
                    #[allow(clippy::cast_precision_loss)]
                    let f = (s.emitted as f64 * 0.618_033_988_749_894_9).fract();
                    let t = (b - a).normalize();
                    let n = DVec3::new(-t.y, t.x, 0.0);
                    let mid = a + (b - a) * f;
                    // Outward: the side away from the electrode's section.
                    let out = if inside(
                        &self.electrodes.electrodes[em.electrode].section,
                        mid + n * 1e-6,
                    ) {
                        -n
                    } else {
                        n
                    };
                    if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) {
                        // Not the coated face: its charge stays.
                        s.pending[i][k] = 0.0;
                        break;
                    }
                    let x = mid + out * (EMISSION_OFFSET * (b - a).length());
                    let tag = u32::try_from(i).unwrap_or(u32::MAX);
                    self.add(s, x, DVec3::ZERO, sign * em.weight, em.charge_per_mass, tag);
                }
            }
        }
    }

    /// One step: launches, emission, push (Boris: with B = 0 the kick-drift
    /// leapfrog), absorption.
    pub fn step(&self, s: &mut TubeState) {
        for (j, p) in self.projectiles.iter().enumerate() {
            if !s.launched[j] && s.t >= p.launch - 0.5 * self.dt {
                s.launched[j] = true;
                let tag = PROJECTILE_TAG + u32::try_from(j).unwrap_or(0);
                self.add(s, p.x0, p.v0, p.charge, p.charge / p.mass, tag);
            }
        }
        self.emit(s);
        crate::cancel::checkpoint();
        let sigma = self.surface_density(s);
        let fields: Vec<DVec3> = (0..s.x.len())
            .map(|k| self.field(s, &sigma, s.x[k], Some(k)))
            .collect();
        for (((x, v), qm), e) in s.x.iter_mut().zip(&mut s.v).zip(&s.qm).zip(&fields) {
            let h = 0.5 * self.dt * qm;
            let minus = *v + *e * h;
            // Rotation about z by the magnetic field (Boris).
            let t = DVec3::new(0.0, 0.0, h * self.b_z);
            let prime = minus + minus.cross(t);
            let plus = minus + prime.cross(t * (2.0 / (1.0 + t.z * t.z)));
            *v = plus + *e * h;
            *x += *v * self.dt;
        }
        s.t += self.dt;
        // Absorption: a particle that crossed an electrode's boundary, or left.
        let n_el = self.electrodes.electrodes.len();
        let mut k = 0;
        while k < s.x.len() {
            let hit = (0..n_el)
                .find(|&e| inside(&self.electrodes.electrodes[e].section, s.x[k]) != s.side[k][e]);
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
                if hit.is_none() {
                    s.lost += s.charge[k];
                }
                s.x.swap_remove(k);
                s.v.swap_remove(k);
                s.charge.swap_remove(k);
                s.qm.swap_remove(k);
                s.tag.swap_remove(k);
                s.side.swap_remove(k);
            } else {
                k += 1;
            }
        }
    }
}
