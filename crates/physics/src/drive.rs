//! Circuits driving the electrodes and coils of a level (PHYSICS.md §2.10). One-way
//! coupling: the circuit is solved first, over all the time the flights can take, and the
//! particles fly in the fields it sets: the electrodes' potentials, the coils' currents
//! and their rates. Their back-action on the circuit (the charge they deposit, the
//! currents they induce) is left out, a marked approximation.

use crate::bem::{Bias, Electrodes};
use crate::circuit::{Circuit, CircuitError, CircuitSolution, Component, Waveform};
use crate::inductance;
use crate::magnetic::CircularLoop;

/// RADAU5's relative tolerance for driving circuits: their potentials and currents are
/// then within about 1e-9 of their scale (tests Z4–Z9).
pub const TOLERANCE: f64 = 1e-12;

/// A series chain driving an electrode or a coil: a source `V(t)` (its terminal `−` on
/// ground) through a resistance `r` (> 0) into the element, with an inductance `l`
/// (electrodes; 0: none) or a capacitance `c` (coils; 0: none) in series. The element
/// closes the chain: an electrode through its capacitances to ground and to the other
/// electrodes, a coil through its inductance to ground. An electrode's chain without an
/// inductance may have its resistance as a switch: `(closed at first, the time it
/// toggles)`; open, it leaves the electrode floating with its charge. (Opening a switch
/// on an inductor's current needs a resistance the ideal components leave out: such
/// circuits are reported as singular.)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Chain {
    pub wave: Waveform,
    pub r: f64,
    pub l: f64,
    pub c: f64,
    pub switch: Option<(bool, f64)>,
}

/// The solved circuit of a level and how it drives the field's sources.
#[derive(Clone, Debug)]
pub struct Drives {
    /// The circuit's solution, in lab time.
    pub solution: CircuitSolution,
    /// The electrodes whose potentials follow the circuit, with their nodes: the driven
    /// ones and, with them, the floating ones, whose charges the circuit conserves.
    pub electrodes: Vec<(usize, usize)>,
    /// Those electrodes' potentials in the static solution (`Electrodes::potentials`,
    /// the state at t = 0), from which their surface charge is shifted.
    pub initial: Vec<f64>,
    /// The driven circular coils: (index in `LevelField::loops`, the coil's inductor in
    /// the circuit).
    pub loops: Vec<(usize, usize)>,
    /// A coil's strength per unit current: `μ₀/4π = 1/c²`.
    pub per_current: f64,
    /// What building and solving the circuit took (for the sandbox's meters).
    pub cost: crate::field::SetupCost,
}

impl Drives {
    /// Builds and solves the circuit of the driven electrodes (`(index, chain)`, with the
    /// static solution `electrodes`) and coils (`(index in loops, chain)`) over `[0, t_end]`
    /// (lab time). The electrodes' capacitances are the boundary-element Maxwell matrix;
    /// a coil's self-inductance is that of a uniform surface current on its wire and the
    /// coils' mutual inductances their linked fluxes (`inductance.rs`), divided by `c²`
    /// (a coil in a circuit needs a finite c). At t = 0 the electrodes are at their
    /// static potentials, the coils carry their static strengths, and the chains' own
    /// inductors and capacitors are empty.
    pub fn build(
        electrodes: &Electrodes,
        driven_electrodes: &[(usize, Chain)],
        loops: &[CircularLoop],
        driven_loops: &[(usize, Chain)],
        c_light: f64,
        t_end: f64,
    ) -> Result<Self, CircuitError> {
        let started = std::time::Instant::now();
        let mut circuit = Circuit::default();
        // Initial node potentials (node k at k − 1) and the initial currents of the
        // components that carry one (inductors, sources), by component.
        let mut v0: Vec<f64> = Vec::new();
        let mut i0: Vec<(usize, f64)> = Vec::new();
        let node = |circuit: &mut Circuit, v0: &mut Vec<f64>, v: f64| {
            circuit.nodes += 1;
            v0.push(v);
            circuit.nodes
        };
        let mut scale_v: f64 = 0.0;
        let mut scale_i: f64 = 0.0;

        // The electrodes' nodes and their capacitance block.
        let mut members: Vec<(usize, usize)> = Vec::new();
        if !driven_electrodes.is_empty() {
            let capacitance = electrodes.capacitance().ok_or(CircuitError::Singular)?;
            for (e, el) in electrodes.electrodes.iter().enumerate() {
                let driven = driven_electrodes.iter().any(|&(d, _)| d == e);
                if driven || matches!(el.bias, Bias::Charge(_)) {
                    let v = electrodes.potentials[e];
                    scale_v = scale_v.max(v.abs());
                    members.push((e, node(&mut circuit, &mut v0, v)));
                }
            }
            let m = members.len();
            let mut block = vec![0.0; m * m];
            for (r, &(er, _)) in members.iter().enumerate() {
                for (s, &(es, _)) in members.iter().enumerate() {
                    block[r + m * s] = capacitance[er][es];
                }
            }
            circuit
                .blocks
                .push((members.iter().map(|&(_, n)| n).collect(), block));
        }
        // A chain from a new source node into `target`; returns nothing.
        let chain = |circuit: &mut Circuit,
                     v0: &mut Vec<f64>,
                     i0: &mut Vec<(usize, f64)>,
                     ch: &Chain,
                     target: usize,
                     series_l: f64,
                     series_c: f64| {
            let s = node(circuit, v0, ch.wave.value(0.0));
            i0.push((circuit.components.len(), 0.0));
            circuit.components.push(Component::Source {
                a: s,
                b: 0,
                wave: ch.wave,
            });
            if series_l > 0.0 {
                let mid = node(circuit, v0, 0.0);
                circuit.components.push(Component::Resistor {
                    a: s,
                    b: mid,
                    r: ch.r,
                });
                i0.push((circuit.components.len(), 0.0));
                circuit.components.push(Component::Inductor {
                    a: mid,
                    b: target,
                    l: series_l,
                });
            } else if series_c > 0.0 {
                let mid = node(circuit, v0, v0[target - 1]);
                circuit.components.push(Component::Resistor {
                    a: s,
                    b: mid,
                    r: ch.r,
                });
                circuit.components.push(Component::Capacitor {
                    a: mid,
                    b: target,
                    c: series_c,
                });
            } else if let Some((closed, at)) = ch.switch {
                circuit.components.push(Component::Switch {
                    a: s,
                    b: target,
                    r: ch.r,
                    closed,
                    toggles: vec![at],
                });
            } else {
                circuit.components.push(Component::Resistor {
                    a: s,
                    b: target,
                    r: ch.r,
                });
            }
        };
        for (e, ch) in driven_electrodes {
            if ch.switch.is_some() && ch.l > 0.0 {
                return Err(CircuitError::Singular);
            }
            let target = members
                .iter()
                .find(|&&(m, _)| m == *e)
                .expect("a driven electrode is a member")
                .1;
            scale_v = scale_v.max(wave_scale(&ch.wave));
            chain(&mut circuit, &mut v0, &mut i0, ch, target, ch.l, 0.0);
        }
        // The coils: each an inductor from its node to ground.
        let per_current = 1.0 / (c_light * c_light);
        let mut coil_loops: Vec<(usize, usize)> = Vec::new();
        for (i, ch) in driven_loops {
            if ch.switch.is_some() {
                return Err(CircuitError::Singular);
            }
            let l = &loops[*i];
            let top = node(&mut circuit, &mut v0, 0.0);
            let current = l.kappa / per_current;
            scale_i = scale_i.max(current.abs()).max(wave_scale(&ch.wave) / ch.r);
            scale_v = scale_v.max(wave_scale(&ch.wave));
            coil_loops.push((*i, circuit.components.len()));
            i0.push((circuit.components.len(), current));
            circuit.components.push(Component::Inductor {
                a: top,
                b: 0,
                l: inductance::ring_self(l) * per_current,
            });
            chain(&mut circuit, &mut v0, &mut i0, ch, top, 0.0, ch.c);
        }
        for (a, &(ia, ka)) in coil_loops.iter().enumerate() {
            for &(ib, kb) in &coil_loops[a + 1..] {
                let m = inductance::mutual_circles(&loops[ia], &loops[ib]) * per_current;
                circuit.mutual.push((ka, kb, m));
            }
        }
        // The initial state: potentials, then the currents in the components' order.
        let mut x0 = v0;
        i0.sort_by_key(|&(k, _)| k);
        x0.extend(i0.iter().map(|&(_, i)| i));
        let scale = [scale_v, scale_i]
            .into_iter()
            .filter(|&s| s > 0.0)
            .fold(f64::INFINITY, f64::min);
        let scale = if scale.is_finite() { scale } else { 1.0 };
        let solution = circuit.solve(&x0, t_end, TOLERANCE, scale, t_end / 64.0)?;
        let n = x0.len();
        let cost = crate::field::SetupCost {
            seconds: started.elapsed().as_secs_f64(),
            // RADAU5's matrices (the Jacobian, the mass matrix, the real and the complex
            // decomposition) and the dense outputs (four coefficients per unknown).
            bytes: 8 * (5 * n * n + 4 * n * solution.steps()),
            unknowns: n,
        };
        Ok(Self {
            solution,
            initial: members
                .iter()
                .map(|&(e, _)| electrodes.potentials[e])
                .collect(),
            electrodes: members,
            loops: coil_loops,
            per_current,
            cost,
        })
    }

    /// The changes of the circuit's electrodes' potentials from the static solution at
    /// lab time `t`: `(electrode, change)`.
    pub fn electrode_shifts(&self, t: f64) -> Vec<(usize, f64)> {
        self.electrodes
            .iter()
            .zip(&self.initial)
            .map(|(&(e, n), &v)| (e, self.solution.potential(n, t) - v))
            .collect()
    }

    /// The strength `κ` of loop `i` and its rate `κ̇` at lab time `t`, if the loop is
    /// driven.
    pub fn loop_strength(&self, i: usize, t: f64) -> Option<(f64, f64)> {
        let &(_, k) = self.loops.iter().find(|&&(l, _)| l == i)?;
        Some((
            self.solution.current(k, t) * self.per_current,
            self.solution.current_rate(k, t) * self.per_current,
        ))
    }

    /// The times at which the drives change abruptly (lab time): switch toggles and pulse
    /// corners, where the integration of the flights restarts.
    pub fn breakpoints(&self) -> &[f64] {
        &self.solution.breakpoints
    }
}

/// The largest magnitude a waveform reaches.
fn wave_scale(w: &Waveform) -> f64 {
    match *w {
        Waveform::Dc(v) => v.abs(),
        Waveform::Sine {
            offset, amplitude, ..
        } => offset.abs() + amplitude.abs(),
        Waveform::Pulse { low, high, .. } => low.abs().max(high.abs()),
    }
}
