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

use glam::DVec3;

use crate::bem::{Bias, Lu};
use crate::circuit::Waveform;
use rayon::prelude::*;

use crate::zinv::{Electrodes, PAR_CHUNK, Section, line_charge, segment_field_fast};

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
    /// electrodes). Its source's voltage adds to that electrode's potential.
    pub to: Option<usize>,
}

impl Circuit {
    /// An ideal source: the electrode follows it.
    fn is_ideal(&self) -> bool {
        self.resistance == 0.0 && self.inductance == 0.0
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
    /// External circuits of electrodes.
    pub circuits: Vec<Circuit>,
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
    /// Charge per unit length that left the arena.
    pub lost: f64,
    /// Which projectiles have been launched, and where each ended.
    launched: Vec<bool>,
    pub fates: Vec<Option<Fate>>,
    /// Emission paused at `MAX_PARTICLES` at least once.
    pub capped: bool,
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
        let cap = &self.electrodes.capacitance;
        let node_charge = self
            .nodes()
            .iter()
            .map(|n| {
                let i = n.electrode;
                (0..n_el).map(|j| cap[i][j] * potentials[j]).sum::<f64>()
                    + n.capacitance * potentials[i]
            })
            .collect();
        TubeState {
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
            collected: vec![vec![0.0; self.emitters.len()]; n_el],
            launched: vec![false; self.projectiles.len()],
            fates: vec![None; self.projectiles.len()],
            potentials,
            node_charge,
            circuit_current: vec![0.0; self.circuits.len()],
            inductor_current: vec![0.0; self.circuits.len()],
            induced: vec![0.0; self.electrodes.segments.len()],
            ..TubeState::default()
        }
    }

    /// The circuit nodes: electrodes with a circuit through a resistance or an inductance
    /// (not an ideal source), then floating electrodes (a fixed charge) without a circuit.
    fn nodes(&self) -> Vec<Node> {
        let mut out: Vec<Node> = self
            .circuits
            .iter()
            .filter(|c| !c.is_ideal())
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

    /// The node index of an electrode, if it is a circuit node.
    fn node_of(nodes: &[Node], electrode: usize) -> Option<usize> {
        nodes.iter().position(|n| n.electrode == electrode)
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
        // Charge balance of each node (its own resistive branch to ground as before).
        for (r, n) in nodes.iter().enumerate() {
            let i = n.electrode;
            for (c, m) in nodes.iter().enumerate() {
                a[r * size + c] = cap[i][m.electrode];
            }
            a[r * size + r] += n.capacitance + n.conductance_dt;
            let fixed: f64 = (0..n_el)
                .filter(|&j| Self::node_of(&nodes, j).is_none())
                .map(|j| cap[i][j] * s.potentials[j])
                .sum();
            b[r] = s.node_charge[r] + n.conductance_dt * n.source.value(t) - fixed - q_ind[i];
        }
        // A resistive branch returning to another electrode: its potential in the
        // branch's drop, and the charge out of it if it is a node.
        for c in &self.circuits {
            let Some(to) = c.to else { continue };
            if c.is_ideal() || c.is_inductive() || !c.resistance.is_finite() {
                continue;
            }
            let g = self.dt / c.resistance;
            let r = Self::node_of(&nodes, c.electrode).expect("a circuit's electrode is a node");
            match Self::node_of(&nodes, to) {
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
            let r = Self::node_of(&nodes, c.electrode).expect("a circuit's electrode is a node");
            a[r * size + u] -= h;
            b[r] += h * i0;
            let (q, v_to0, v_to1) = match c.to {
                Some(to) => (Self::node_of(&nodes, to), old[to], s.potentials[to]),
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
            if let Some(q) = circuit.to.and_then(|to| Self::node_of(&nodes, to)) {
                s.node_charge[q] -= delivered;
            }
            s.circuit_current[c] = delivered / self.dt;
        }
        for (j, &ci) in inductive.iter().enumerate() {
            let c = &self.circuits[ci];
            let i1 = x[k + j];
            let delivered = h * (s.inductor_current[ci] + i1);
            let r = Self::node_of(&nodes, c.electrode).expect("a circuit's electrode is a node");
            s.node_charge[r] += delivered;
            if let Some(q) = c.to.and_then(|to| Self::node_of(&nodes, to)) {
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
                    self.add(s, x, DVec3::ZERO, sign * em.weight, em.charge_per_mass, tag);
                    emitted[em.electrode] += sign * em.weight;
                }
            }
        }
        emitted
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
            if let Some(r) = Self::node_of(&nodes, e) {
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
        let fields: Vec<DVec3> = (0..s.x.len())
            .into_par_iter()
            .with_min_len(PAR_CHUNK)
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
                if let Some(r) = hit.and_then(|e| Self::node_of(&nodes, e)) {
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
                s.left.swap_remove(k);
                s.side.swap_remove(k);
            } else {
                k += 1;
            }
        }
        // The moved particles' induced density, for the next step and the display.
        s.induced = self.induced_density(s);
    }
}
