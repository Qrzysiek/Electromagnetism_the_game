//! Lumped linear circuits (docs/CIRCUITS.md, PHYSICS.md §2.9): modified nodal analysis,
//! integrated by RADAU5.
//!
//! Units as everywhere (`k = 1`, `μ₀/4π = 1/c²`): `Q = C V`, `V = R I`, `V = L dI/dt`.
//!
//! **Equations.** Node 0 is ground; the unknowns are the potentials of nodes 1..=n, the
//! inductors' currents (each flowing from its node `a` to its node `b`) and the voltage
//! sources' currents (each flowing out of its node `a` through the source to `b`). Each
//! node's Kirchhoff current law, the inductors' `L dI/dt = v_a − v_b` (with the mutual
//! inductances), and the sources' `v_a − v_b = V(t)` make `M x' = A x + b(t)`:
//!
//! ```text
//! M = [C 0 0; 0 L 0; 0 0 0],   A = [−G −E_L −E_V; E_Lᵀ 0 0; E_Vᵀ 0 0],   b = [0; 0; −V(t)]
//! ```
//!
//! with C the nodes' capacitance matrix (capacitors, and blocks such as the electrodes'
//! Maxwell matrix), G the conductances, E the incidence matrices. With a resistor in series
//! with every source and capacitor loop it is an index-1 DAE.
//!
//! **Breakpoints.** Switches change the topology at set times, pulses have corners: the
//! integration restarts at each breakpoint. The charges and fluxes `M x` are continuous
//! across it (no impulsive currents in index-1 circuits), and the state after it is the
//! consistent one with those charges and fluxes: `[M_d; A_a] x = [M_d x⁻; −b_a(t)]`, the
//! nonzero rows of M and the algebraic rows of the new system.

use crate::integrator::radau5::{self, Dense, Radau5, Settings, StiffSystem};

/// A source's voltage as a function of time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Waveform {
    Dc(f64),
    /// `offset + amplitude sin(ω t + φ)`.
    Sine {
        offset: f64,
        amplitude: f64,
        omega: f64,
        phase: f64,
    },
    /// A trapezoid: `low` until `delay`, rising linearly over `rise` to `high`, held for
    /// `width`, falling over `fall`; repeated every `period` if that is positive.
    Pulse {
        low: f64,
        high: f64,
        delay: f64,
        rise: f64,
        width: f64,
        fall: f64,
        period: f64,
    },
}

impl Waveform {
    pub fn value(&self, t: f64) -> f64 {
        match *self {
            Waveform::Dc(v) => v,
            Waveform::Sine {
                offset,
                amplitude,
                omega,
                phase,
            } => offset + amplitude * libm::sin(omega * t + phase),
            Waveform::Pulse {
                low,
                high,
                delay,
                rise,
                width,
                fall,
                period,
            } => {
                let mut s = t - delay;
                if s < 0.0 {
                    return low;
                }
                if period > 0.0 {
                    s -= period * libm::floor(s / period);
                }
                if s < rise {
                    low + (high - low) * s / rise
                } else if s < rise + width {
                    high
                } else if s < rise + width + fall {
                    high + (low - high) * (s - rise - width) / fall
                } else {
                    low
                }
            }
        }
    }

    /// The times in `(0, t_end)` where the waveform has a corner.
    pub fn corners(&self, t_end: f64) -> Vec<f64> {
        let Waveform::Pulse {
            delay,
            rise,
            width,
            fall,
            period,
            ..
        } = *self
        else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut start = delay;
        loop {
            for c in [
                start,
                start + rise,
                start + rise + width,
                start + rise + width + fall,
            ] {
                if c > 0.0 && c < t_end {
                    out.push(c);
                }
            }
            if period <= 0.0 || start + period >= t_end {
                break;
            }
            start += period;
        }
        out
    }
}

/// A two-terminal component between nodes `a` and `b` (0 is ground).
#[derive(Clone, Debug, PartialEq)]
pub enum Component {
    Resistor {
        a: usize,
        b: usize,
        r: f64,
    },
    Capacitor {
        a: usize,
        b: usize,
        c: f64,
    },
    Inductor {
        a: usize,
        b: usize,
        l: f64,
    },
    /// `v_a − v_b = V(t)`.
    Source {
        a: usize,
        b: usize,
        wave: Waveform,
    },
    /// A resistor `r` while closed, nothing while open; `toggles` are the times at which it
    /// changes state, starting from `closed`.
    Switch {
        a: usize,
        b: usize,
        r: f64,
        closed: bool,
        toggles: Vec<f64>,
    },
}

/// A circuit: `nodes` nodes besides ground, its components, capacitance blocks (a
/// symmetric matrix over some nodes, e.g. the electrodes' Maxwell matrix: their charges
/// `Q = C V`), and mutual inductances between pairs of its inductors (indices into
/// `components`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Circuit {
    pub nodes: usize,
    pub components: Vec<Component>,
    pub blocks: Vec<(Vec<usize>, Vec<f64>)>,
    pub mutual: Vec<(usize, usize, f64)>,
}

/// Where the unknowns are: node `k` (1..=n) at `k − 1`, then a current for each inductor
/// and each source, in the order of the components.
#[derive(Clone, Debug, PartialEq)]
struct Layout {
    n: usize,
    /// The current's index of each component (inductors and sources).
    current: Vec<Option<usize>>,
}

impl Circuit {
    fn layout(&self) -> Layout {
        let mut n = self.nodes;
        let current = self
            .components
            .iter()
            .map(|c| match c {
                Component::Inductor { .. } | Component::Source { .. } => {
                    n += 1;
                    Some(n - 1)
                }
                _ => None,
            })
            .collect();
        Layout { n, current }
    }

    /// The number of unknowns: the nodes' potentials, then the inductors' and the sources'
    /// currents.
    pub fn unknowns(&self) -> usize {
        self.layout().n
    }

    /// The times in `(0, t_end)` at which switches toggle or pulses have corners, sorted.
    pub fn breakpoints(&self, t_end: f64) -> Vec<f64> {
        let mut out: Vec<f64> = Vec::new();
        for c in &self.components {
            match c {
                Component::Source { wave, .. } => out.extend(wave.corners(t_end)),
                Component::Switch { toggles, .. } => {
                    out.extend(toggles.iter().copied().filter(|&t| t > 0.0 && t < t_end));
                }
                _ => {}
            }
        }
        out.sort_by(f64::total_cmp);
        out.dedup();
        out
    }

    /// Whether each component, if a switch, is closed just after time `t`.
    fn switch_states(&self, t: f64) -> Vec<bool> {
        self.components
            .iter()
            .map(|c| match c {
                Component::Switch {
                    closed, toggles, ..
                } => {
                    let flips = toggles.iter().filter(|&&s| s <= t).count();
                    *closed ^ (flips % 2 == 1)
                }
                _ => false,
            })
            .collect()
    }

    /// The equations with the switches as `closed` says.
    fn mna(&self, layout: &Layout, closed: &[bool]) -> Mna {
        let n = layout.n;
        let mut mass = vec![0.0; n * n];
        let mut a = vec![0.0; n * n];
        let ix = |i: usize, j: usize| i + n * j;
        // An admittance `g` between nodes `na` and `nb` (0: ground) into the matrix `m`.
        let stamp = |m: &mut [f64], na: usize, nb: usize, g: f64| {
            if na > 0 {
                m[ix(na - 1, na - 1)] += g;
            }
            if nb > 0 {
                m[ix(nb - 1, nb - 1)] += g;
            }
            if na > 0 && nb > 0 {
                m[ix(na - 1, nb - 1)] -= g;
                m[ix(nb - 1, na - 1)] -= g;
            }
        };
        // A branch current `i` from node `na` to `nb`: it leaves `na`'s current law and
        // enters `nb`'s, and its own row reads `v_a − v_b`.
        let branch = |m: &mut [f64], i: usize, na: usize, nb: usize| {
            if na > 0 {
                m[ix(na - 1, i)] -= 1.0;
                m[ix(i, na - 1)] += 1.0;
            }
            if nb > 0 {
                m[ix(nb - 1, i)] += 1.0;
                m[ix(i, nb - 1)] -= 1.0;
            }
        };
        let mut sources = Vec::new();
        for (k, c) in self.components.iter().enumerate() {
            match *c {
                Component::Resistor { a: na, b: nb, r } => stamp(&mut a, na, nb, -1.0 / r),
                Component::Switch {
                    a: na, b: nb, r, ..
                } => {
                    if closed[k] {
                        stamp(&mut a, na, nb, -1.0 / r);
                    }
                }
                Component::Capacitor { a: na, b: nb, c } => stamp(&mut mass, na, nb, c),
                Component::Inductor { a: na, b: nb, l } => {
                    let i = layout.current[k].expect("an inductor has a current");
                    branch(&mut a, i, na, nb);
                    mass[ix(i, i)] += l;
                }
                Component::Source { a: na, b: nb, wave } => {
                    let i = layout.current[k].expect("a source has a current");
                    branch(&mut a, i, na, nb);
                    sources.push((i, wave));
                }
            }
        }
        for (nodes, matrix) in &self.blocks {
            let m = nodes.len();
            for (r, &nr) in nodes.iter().enumerate() {
                for (s, &ns) in nodes.iter().enumerate() {
                    if nr > 0 && ns > 0 {
                        mass[ix(nr - 1, ns - 1)] += matrix[r + m * s];
                    }
                }
            }
        }
        for &(k1, k2, m) in &self.mutual {
            let i1 = layout.current[k1].expect("a mutual inductance joins inductors");
            let i2 = layout.current[k2].expect("a mutual inductance joins inductors");
            mass[ix(i1, i2)] += m;
            mass[ix(i2, i1)] += m;
        }
        Mna {
            n,
            mass,
            a,
            sources,
        }
    }
}

/// The circuit's equations `M x' = A x + b(t)` between two breakpoints.
struct Mna {
    n: usize,
    mass: Vec<f64>,
    a: Vec<f64>,
    /// The source rows and their waveforms (`b = −V(t)` there).
    sources: Vec<(usize, Waveform)>,
}

impl StiffSystem for Mna {
    fn dim(&self) -> usize {
        self.n
    }

    fn rhs(&self, t: f64, y: &[f64], f: &mut [f64]) {
        let n = self.n;
        for (i, fi) in f.iter_mut().enumerate() {
            let mut s = 0.0;
            for (j, yj) in y.iter().enumerate() {
                s += self.a[i + n * j] * yj;
            }
            *fi = s;
        }
        for &(i, w) in &self.sources {
            f[i] -= w.value(t);
        }
    }

    fn mass(&self) -> Option<Vec<f64>> {
        Some(self.mass.clone())
    }

    fn jacobian(&self, _t: f64, _y: &[f64], jac: &mut [f64]) -> bool {
        jac.copy_from_slice(&self.a);
        true
    }
}

impl Mna {
    /// The consistent state at time `t` whose charges and fluxes `M x` are those of `x0`:
    /// on the row space of M (found by elimination, `T M` in row echelon form) `M x = M x0`,
    /// on its left null space the algebraic equations `A x + b(t) = 0`. None if that system
    /// is singular (a circuit of index above 1).
    fn consistent(&self, x0: &[f64], t: f64) -> Option<Vec<f64>> {
        let n = self.n;
        let ix = |i: usize, j: usize| i + n * j;
        let mut r = self.mass.clone();
        let mut tm: Vec<f64> = (0..n * n)
            .map(|k| f64::from(u8::from(k % n == k / n)))
            .collect();
        let scale = r.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1e-300);
        let mut row = 0;
        let mut differential = vec![false; n];
        for col in 0..n {
            if row == n {
                break;
            }
            let (mut best, mut best_val) = (row, 0.0);
            for i in row..n {
                if r[ix(i, col)].abs() > best_val {
                    best = i;
                    best_val = r[ix(i, col)].abs();
                }
            }
            if best_val <= 1e-13 * scale {
                continue;
            }
            for j in 0..n {
                r.swap(ix(best, j), ix(row, j));
                tm.swap(ix(best, j), ix(row, j));
            }
            for i in 0..n {
                if i != row {
                    let f = r[ix(i, col)] / r[ix(row, col)];
                    if f != 0.0 {
                        for j in 0..n {
                            r[ix(i, j)] -= f * r[ix(row, j)];
                            tm[ix(i, j)] -= f * tm[ix(row, j)];
                        }
                    }
                }
            }
            differential[row] = true;
            row += 1;
        }
        let mut b = vec![0.0; n];
        for &(i, w) in &self.sources {
            b[i] = -w.value(t);
        }
        let mx0: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| self.mass[ix(i, j)] * x0[j]).sum())
            .collect();
        let mut sys = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for (row, &diff) in differential.iter().enumerate() {
            let m = if diff { &self.mass } else { &self.a };
            for j in 0..n {
                sys[ix(row, j)] = (0..n).map(|k| tm[ix(row, k)] * m[ix(k, j)]).sum();
            }
            rhs[row] = (0..n)
                .map(|k| tm[ix(row, k)] * if diff { mx0[k] } else { -b[k] })
                .sum();
        }
        solve_dense(n, &mut sys, &mut rhs).then_some(rhs)
    }
}

/// Solves the column-major n×n system `a x = b` by Gaussian elimination with partial
/// pivoting; `b` becomes x. False if singular.
fn solve_dense(n: usize, a: &mut [f64], b: &mut [f64]) -> bool {
    let ix = |i: usize, j: usize| i + n * j;
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1e-300);
    for k in 0..n {
        let mut p = k;
        for i in k + 1..n {
            if a[ix(i, k)].abs() > a[ix(p, k)].abs() {
                p = i;
            }
        }
        if a[ix(p, k)].abs() <= 1e-14 * scale {
            return false;
        }
        if p != k {
            for j in 0..n {
                a.swap(ix(p, j), ix(k, j));
            }
            b.swap(p, k);
        }
        for i in k + 1..n {
            let f = a[ix(i, k)] / a[ix(k, k)];
            if f != 0.0 {
                for j in k..n {
                    a[ix(i, j)] -= f * a[ix(k, j)];
                }
                b[i] -= f * b[k];
            }
        }
    }
    for k in (0..n).rev() {
        let mut s = b[k];
        for j in k + 1..n {
            s -= a[ix(k, j)] * b[j];
        }
        b[k] = s / a[ix(k, k)];
    }
    true
}

/// The circuit's solution over `[0, t_end]`: the collocation polynomials of RADAU5's steps.
#[derive(Clone, Debug)]
pub struct CircuitSolution {
    layout: Layout,
    segments: Vec<Dense>,
    /// The consistent state at t = 0.
    start: Vec<f64>,
    /// The breakpoints (the integration restarted there).
    pub breakpoints: Vec<f64>,
    t_end: f64,
    /// The inductors (component, node a, node b) and the inverse of their inductance
    /// matrix (column-major), for the currents' rates.
    inductors: Vec<(usize, usize, usize)>,
    inverse_inductance: Vec<f64>,
}

impl CircuitSolution {
    /// The step whose polynomial gives the state at `t` (≥ 0): right-continuous, so that
    /// at a breakpoint it is the step after it (the integration of the flights restarts
    /// there). The steps are told apart by their ends only: a step's recorded start is its
    /// end minus its length, which may differ from the previous end by an ulp.
    fn segment(&self, t: f64) -> Option<&Dense> {
        let k = self.segments.partition_point(|d| d.t_end() <= t);
        self.segments
            .get(k.min(self.segments.len().saturating_sub(1)))
    }

    fn value(&self, i: usize, t: f64) -> f64 {
        match self.segment(t) {
            Some(d) if t >= 0.0 => d.eval_component(i, t.min(d.t_end())),
            _ => self.start[i],
        }
    }

    /// The potential of node `node` at time `t` (0 for ground; before the start the
    /// initial state, after the end the last).
    pub fn potential(&self, node: usize, t: f64) -> f64 {
        if node == 0 {
            0.0
        } else {
            self.value(node - 1, t)
        }
    }

    /// The current of component `k` (an inductor or a source) at time `t`. An inductor's
    /// is its value at the start of RADAU5's step plus the integral of its rate
    /// (`current_rate`) over the step, so that the rate is exactly its derivative (a coil's
    /// induced field then keeps Faraday's law with its field exactly); at the steps' ends
    /// it is RADAU5's value, by Radau quadrature's orthogonality.
    pub fn current(&self, k: usize, t: f64) -> f64 {
        let i = self.layout.current[k].expect("a component with a current");
        let Some(p) = self.inductors.iter().position(|&(c, ..)| c == k) else {
            return self.value(i, t);
        };
        let Some(d) = self.segment(t) else {
            return self.start[i];
        };
        if t < 0.0 {
            return self.start[i];
        }
        let x = t.min(d.t_end());
        let integral = |node: usize| {
            if node == 0 {
                0.0
            } else {
                d.integral_component(node - 1, x)
            }
        };
        let m = self.inductors.len();
        let mut current = d.eval_component(i, d.t_start());
        for (j, &(_, a, b)) in self.inductors.iter().enumerate() {
            current += self.inverse_inductance[p + m * j] * (integral(a) - integral(b));
        }
        current
    }

    /// The rate of change of inductor `k`'s current at time `t`, from the inductors' law
    /// `L dI/dt = v_a − v_b` (with the mutual inductances) at the dense output's potentials,
    /// so as accurate as they are; 0 outside `[0, t_end]`, where the state is held. (The
    /// derivative of the dense output is a power of the step less accurate. Within a step
    /// the two differ by a cubic vanishing at the three collocation points, which
    /// integrates to zero over the step, Radau quadrature being exact to degree 4: the
    /// current's change over each step is the integral of this rate, and `current` is
    /// that integral.)
    pub fn current_rate(&self, k: usize, t: f64) -> f64 {
        if !(0.0..=self.t_end).contains(&t) {
            return 0.0;
        }
        let m = self.inductors.len();
        let p = self
            .inductors
            .iter()
            .position(|&(c, ..)| c == k)
            .expect("an inductor");
        self.inductors
            .iter()
            .enumerate()
            .map(|(j, &(_, a, b))| {
                self.inverse_inductance[p + m * j] * (self.potential(a, t) - self.potential(b, t))
            })
            .sum()
    }

    /// The times at which RADAU5's steps ended: the polynomials join there, continuous
    /// but with their derivatives changing.
    pub fn step_ends(&self) -> Vec<f64> {
        self.segments.iter().map(Dense::t_end).collect()
    }

    /// The whole state at time `t`.
    pub fn state(&self, t: f64) -> Vec<f64> {
        (0..self.layout.n).map(|i| self.value(i, t)).collect()
    }

    /// The number of RADAU5 steps taken.
    pub fn steps(&self) -> usize {
        self.segments.len()
    }
}

/// Why a circuit could not be solved.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CircuitError {
    /// Its equations are of index above 1 (a loop of capacitors and voltage sources, a cut
    /// set of inductors), or singular.
    Singular,
    Integration(radau5::Error),
}

impl Circuit {
    /// Solves the circuit from `x0` (its charges and fluxes `M x0` kept, the rest made
    /// consistent) to `t_end`: RADAU5 at the relative tolerance `tol` (absolute `tol ·
    /// scale`), steps at most `h_max`, restarting at every breakpoint with the charges and
    /// fluxes carried over.
    pub fn solve(
        &self,
        x0: &[f64],
        t_end: f64,
        tol: f64,
        scale: f64,
        h_max: f64,
    ) -> Result<CircuitSolution, CircuitError> {
        let layout = self.layout();
        assert_eq!(x0.len(), layout.n, "one initial value per unknown");
        let breakpoints = self.breakpoints(t_end);
        let mut cuts = vec![0.0];
        cuts.extend(&breakpoints);
        cuts.push(t_end);
        let mut x = x0.to_vec();
        let mut segments = Vec::new();
        let mut start = Vec::new();
        for w in cuts.windows(2) {
            let (ta, tb) = (w[0], w[1]);
            let mna = self.mna(&layout, &self.switch_states(ta));
            x = mna.consistent(&x, ta).ok_or(CircuitError::Singular)?;
            if start.is_empty() {
                start.clone_from(&x);
            }
            let mut settings = Settings::new(tol, tol * scale);
            settings.h_max = h_max;
            let mut r = Radau5::new(&mna, ta, &x, tb, settings);
            while !r.done() {
                r.step(&mna).map_err(CircuitError::Integration)?;
                segments.push(r.dense());
            }
            x = r.y().to_vec();
        }
        let (inductors, inverse_inductance) = self.inverse_inductance()?;
        Ok(CircuitSolution {
            layout,
            segments,
            start,
            breakpoints,
            t_end,
            inductors,
            inverse_inductance,
        })
    }

    /// The inductors (component, node a, node b) and the inverse of their inductance
    /// matrix, column-major.
    #[allow(clippy::type_complexity)]
    fn inverse_inductance(&self) -> Result<(Vec<(usize, usize, usize)>, Vec<f64>), CircuitError> {
        let inductors: Vec<(usize, usize, usize)> = self
            .components
            .iter()
            .enumerate()
            .filter_map(|(k, c)| match *c {
                Component::Inductor { a, b, .. } => Some((k, a, b)),
                _ => None,
            })
            .collect();
        let m = inductors.len();
        let position = |k: usize| {
            inductors
                .iter()
                .position(|&(c, ..)| c == k)
                .expect("a mutual inductance joins inductors")
        };
        let mut matrix = vec![0.0; m * m];
        for (p, &(k, ..)) in inductors.iter().enumerate() {
            if let Component::Inductor { l, .. } = self.components[k] {
                matrix[p + m * p] += l;
            }
        }
        for &(k1, k2, mutual) in &self.mutual {
            let (p1, p2) = (position(k1), position(k2));
            matrix[p1 + m * p2] += mutual;
            matrix[p2 + m * p1] += mutual;
        }
        let mut inverse = vec![0.0; m * m];
        for j in 0..m {
            let mut a = matrix.clone();
            let mut column = vec![0.0; m];
            column[j] = 1.0;
            if !solve_dense(m, &mut a, &mut column) {
                return Err(CircuitError::Singular);
            }
            inverse[m * j..m * (j + 1)].copy_from_slice(&column);
        }
        Ok((inductors, inverse))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The inductors' rates (from their law at the dense output's potentials) integrate
    /// over each RADAU5 step to the step's change of the currents: Radau quadrature's
    /// orthogonality (see `current_rate`), exact for a linear circuit up to rounding.
    /// Coupled inductors driven through resistors (test Z8's circuit, with a pulse);
    /// required 1e-12 of the currents' scale.
    #[test]
    fn rates_integrate_to_the_step_changes() {
        let circuit = Circuit {
            nodes: 3,
            components: vec![
                Component::Source {
                    a: 1,
                    b: 0,
                    wave: Waveform::Pulse {
                        low: 0.0,
                        high: 1.0,
                        delay: 0.5,
                        rise: 0.25,
                        width: 2.0,
                        fall: 0.25,
                        period: 0.0,
                    },
                },
                Component::Resistor { a: 1, b: 2, r: 1.0 },
                Component::Inductor { a: 2, b: 0, l: 1.0 },
                Component::Inductor { a: 3, b: 0, l: 2.0 },
                Component::Resistor { a: 3, b: 0, r: 2.0 },
            ],
            mutual: vec![(2, 3, 0.8)],
            ..Circuit::default()
        };
        let sol = circuit
            .solve(&vec![0.0; circuit.unknowns()], 6.0, 1e-12, 1.0, 0.05)
            .expect("solves");
        // Two-point Gauss–Legendre: exact for the cubic rate within a step.
        let g = 0.5 / 3.0_f64.sqrt();
        let mut worst: f64 = 0.0;
        let mut t0 = 0.0;
        for t1 in sol.step_ends() {
            let (mid, half) = (0.5 * (t0 + t1), 0.5 * (t1 - t0));
            for k in [2, 3] {
                let integral = half
                    * (sol.current_rate(k, mid - 2.0 * g * half)
                        + sol.current_rate(k, mid + 2.0 * g * half));
                let change = sol.current(k, t1) - sol.current(k, t0);
                worst = worst.max((integral - change).abs());
            }
            t0 = t1;
        }
        println!(
            "rates against the step changes: {worst:.1e} over {} steps",
            sol.steps()
        );
        assert!(worst < 1e-12, "{worst:.3e}");
    }
}
