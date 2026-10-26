//! The circuits: their nodes (electrodes strapped together), the supplies' and floating
//! electrodes' potentials, the static surface charge.

use super::*;

impl Tube {
    /// The circuit nodes: electrodes with a circuit through a resistance or an inductance
    /// (not an ideal source), then floating electrodes (a fixed charge) without a circuit.
    pub(super) fn nodes(&self) -> Vec<Node> {
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
    pub(super) fn node_of(&self, nodes: &[Node], electrode: usize) -> Option<usize> {
        let e = self.strap(electrode).unwrap_or(electrode);
        nodes.iter().position(|n| n.electrode == e)
    }

    /// The electrode a strapped one follows.
    pub(super) fn strap(&self, electrode: usize) -> Option<usize> {
        self.circuits
            .iter()
            .find(|c| c.electrode == electrode && c.is_strap())
            .and_then(|c| c.to)
    }

    /// Node `r`'s electrodes: its own, then those strapped to it.
    pub(super) fn members<'a>(
        &'a self,
        nodes: &'a [Node],
        r: usize,
    ) -> impl Iterator<Item = usize> + 'a {
        let own = nodes[r].electrode;
        std::iter::once(own).chain(
            self.circuits
                .iter()
                .filter(move |c| c.is_strap() && c.to == Some(own))
                .map(|c| c.electrode),
        )
    }

    /// Every strapped electrode at its other end's potential.
    pub(super) fn follow_straps(&self, potentials: &mut [f64]) {
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
    pub(super) fn update_potentials(&self, s: &mut TubeState, t: f64, q_ind: &[f64]) {
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
    pub(super) fn static_density(&self, potentials: &[f64]) -> Vec<f64> {
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
}
