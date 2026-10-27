//! A tube level's simulation: its steps, its measurements and its results.

use super::*;

/// The amplitude and phase (degrees, relative to `sin(ω t)`) of the oscillation at angular
/// frequency `omega` in the samples `(t, v)`: the least-squares fit of
/// `a + b (t − t̄) + c cos ωt + d sin ωt`, all four together (exact for a sine on a line,
/// over any window; fitting the line first and then the sine took part of the sine with
/// the line: 1e-3 of the amplitude over 10 periods, up to tens of per cent over one), and
/// `A = √(c² + d²)`, `φ = atan2(c, d)` for `v ≈ A sin(ω t + φ)`.
#[allow(clippy::needless_range_loop)] // index form mirrors the algorithm
pub(super) fn fourier_at(samples: &[(f64, f64)], omega: f64) -> (f64, f64) {
    #[allow(clippy::cast_precision_loss)]
    let n = samples.len() as f64;
    if samples.len() < 4 {
        return (0.0, 0.0);
    }
    let mean_t = samples.iter().map(|(t, _)| t).sum::<f64>() / n;
    // Normal equations of the four basis functions.
    let mut m = [[0.0_f64; 5]; 4];
    for (t, v) in samples {
        let (s, c) = libm::sincos(omega * t);
        let f = [1.0, t - mean_t, c, s];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += f[i] * f[j];
            }
            m[i][4] += f[i] * v;
        }
    }
    // Gaussian elimination with partial pivoting (a 4 × 4 system, fixed order).
    for k in 0..4 {
        let p = (k..4)
            .max_by(|&i, &j| m[i][k].abs().total_cmp(&m[j][k].abs()))
            .unwrap_or(k);
        m.swap(k, p);
        if m[k][k] == 0.0 {
            return (0.0, 0.0);
        }
        for i in k + 1..4 {
            let f = m[i][k] / m[k][k];
            for j in k..5 {
                m[i][j] -= f * m[k][j];
            }
        }
    }
    let mut x = [0.0_f64; 4];
    for i in (0..4).rev() {
        let rest: f64 = (i + 1..4).map(|j| m[i][j] * x[j]).sum();
        x[i] = (m[i][4] - rest) / m[i][i];
    }
    let (c, d) = (x[2], x[3]);
    (libm::sqrt(c * c + d * d), libm::atan2(c, d).to_degrees())
}

/// Distance from `x` to the box `(lo, hi)` (0 inside).
pub(super) fn box_distance(x: DVec3, (lo, hi): (DVec3, DVec3)) -> f64 {
    let dx = (lo.x - x.x).max(x.x - hi.x).max(0.0);
    let dy = (lo.y - x.y).max(x.y - hi.y).max(0.0);
    libm::sqrt(dx * dx + dy * dy)
}

impl TubeSim {
    pub fn t(&self) -> f64 {
        self.state.t
    }

    /// Books the contacts' potential steps from now on (`TubeRun::potential_steps`; a
    /// diagnostic: ~40 % more time).
    pub fn book_potential_steps(&mut self) {
        self.state.book_potential_steps = true;
    }

    /// Steps taken and the segments of its linear system (for the cost meters).
    pub fn steps(&self) -> u32 {
        self.steps
    }

    pub fn segment_count(&self) -> usize {
        self.tube.electrodes.segments.len()
    }

    /// The field at `x` now, as the carriers feel it (the physics' softening).
    #[must_use]
    pub fn field_at(&self, x: DVec3) -> DVec3 {
        self.tube.field_at(&self.state, x)
    }

    /// The macroparticles' positions and charges now.
    #[must_use]
    pub fn particles(&self) -> (&[DVec3], &[f64]) {
        (&self.state.x, &self.state.charge)
    }

    /// Macroparticles in flight now.
    pub fn in_flight(&self) -> usize {
        self.state.x.len()
    }

    /// The electrodes' potentials now.
    pub fn potentials(&self) -> &[f64] {
        &self.state.potentials
    }

    /// Each emitter's macroparticle weight.
    pub fn weights(&self) -> Vec<f64> {
        self.tube.emitters.iter().map(|e| e.weight).collect()
    }

    /// Each current goal's current, once the run has passed its window's end.
    pub fn currents(&self) -> &[Option<f64>] {
        &self.currents
    }

    /// The charge per unit length each electrode has collected so far, per emitter
    /// (signed: electrons negative), and what left the board.
    pub fn collected(&self) -> &[Vec<f64>] {
        &self.state.collected
    }

    pub fn lost(&self) -> f64 {
        self.state.lost
    }

    /// Emission paused at `physics::tube::MAX_PARTICLES` at least once.
    pub fn capped(&self) -> bool {
        self.state.capped
    }

    /// Where each projectile ended (None: still flying, or not launched yet).
    pub fn fates(&self) -> &[Option<physics::tube::Fate>] {
        &self.state.fates
    }

    /// When each projectile reached its detector (None: not, or no detector).
    pub fn arrivals(&self) -> &[Option<f64>] {
        &self.arrivals
    }

    /// Each projectile's closest approach to its detector so far (0: arrived; infinite
    /// without a detector).
    pub fn closest(&self) -> &[f64] {
        &self.closest
    }

    /// Advances to time `t_end` (whole steps), recording a frame every step.
    pub fn advance_to(&mut self, t_end: f64) {
        let n_proj = self.spec.projectiles.len();
        let n_goals = self.spec.goals.len();
        while self.state.t < t_end - 0.5 * self.tube.dt {
            self.tube.step(&mut self.state);
            self.steps += 1;
            let s = &self.state;
            if self.pending.collected.len() < n_goals {
                self.pending.collected.resize(n_goals, Vec::new());
                self.pending.currents.resize(n_goals, None);
            }
            for (k, g) in self.spec.goals.iter().enumerate() {
                // The carriers' charge collected, positive for arrivals (a lead's: what
                // its electrode takes in less what it gives off).
                let sign = self.spec.goal_sign(g);
                let collected = sign
                    * match g.carriers.emitter() {
                        Some(i) => s.collected[g.electrode].get(i).copied().unwrap_or(0.0),
                        None => {
                            s.collected[g.electrode].iter().sum::<f64>()
                                - self
                                    .tube
                                    .emitters
                                    .iter()
                                    .zip(&s.released)
                                    .filter(|(em, _)| em.electrode == g.electrode)
                                    .map(|(_, q)| q)
                                    .sum::<f64>()
                        }
                    };
                if self.at_start[k].is_none() && s.t >= g.start {
                    self.at_start[k] = Some((s.t, collected));
                }
                if self.currents[k].is_none() && s.t >= g.end - 0.5 * self.tube.dt {
                    let (t0, q0) = self.at_start[k].unwrap_or((0.0, 0.0));
                    self.currents[k] = Some(if s.t > t0 {
                        (collected - q0) / (s.t - t0)
                    } else {
                        0.0
                    });
                    self.pending.currents[k] = self.currents[k];
                }
                self.pending.collected[k].push((s.t, collected));
            }
            if self.pending.tracks.len() < n_proj {
                self.pending.tracks.resize(n_proj, Vec::new());
                self.pending.arrivals.resize(n_proj, None);
            }
            for j in 0..n_proj {
                let Some(x) = s.projectile(j) else { continue };
                #[allow(clippy::cast_possible_truncation)]
                self.pending.tracks[j].push((s.t, [x.x as f32, x.y as f32]));
                if let Some(d) = self.detectors[j] {
                    let dist = box_distance(x, d);
                    self.closest[j] = self.closest[j].min(dist);
                    if dist == 0.0 && self.arrivals[j].is_none() {
                        self.arrivals[j] = Some(s.t);
                        self.pending.arrivals[j] = Some(s.t);
                    }
                }
            }
            // The sampled carriers every step (between frames).
            let every = self.layout.fine_every;
            if every > 0 {
                for ((x, id), tag) in s.x.iter().zip(&s.id).zip(&s.tag) {
                    if *tag < PROJECTILE_TAG && id.is_multiple_of(every) {
                        #[allow(clippy::cast_possible_truncation)]
                        self.fine
                            .entry(*id)
                            .or_default()
                            .push([x.x as f32, x.y as f32]);
                    }
                }
            }
            // A frame every `frame_stride` steps (a metal's thousands of carriers every
            // step would fill gigabytes over a long run).
            if self.steps.is_multiple_of(self.frame_stride) {
                #[allow(clippy::cast_possible_truncation)]
                let frame = TubeFrame {
                    self_b: match (&s.darwin, &self.self_b_grid) {
                        (Some(d), Some(g)) => {
                            let b = d.b_z();
                            g.points().map(|(x, y)| d.b_at(&b, x, y) as f32).collect()
                        }
                        _ => Vec::new(),
                    },
                    t: s.t,
                    x: s.x.iter().map(|p| [p.x as f32, p.y as f32]).collect(),
                    q: s.charge.iter().map(|q| *q as f32).collect(),
                    kind: s
                        .tag
                        .iter()
                        .zip(&s.charge)
                        .map(|(t, q)| {
                            if *t >= PROJECTILE_TAG {
                                KIND_PROJECTILE
                            } else if *q < 0.0 {
                                KIND_NEGATIVE
                            } else {
                                KIND_POSITIVE
                            }
                        })
                        .collect(),
                    sigma: self
                        .tube
                        .surface_density(s)
                        .iter()
                        .map(|v| *v as f32)
                        .collect(),
                    potentials: s.potentials.iter().map(|v| *v as f32).collect(),
                    id: s.id.clone(),
                    currents: {
                        // Each tube circuit's first physics circuit (its linked electrodes'
                        // follow it).
                        let mut at = 0;
                        self.spec
                            .circuits
                            .iter()
                            .map(|c| {
                                let k = at;
                                at += 1 + c.linked.len();
                                let i = if self.tube.circuits[k].inductance > 0.0 {
                                    s.inductor_current[k]
                                } else {
                                    s.circuit_current[k]
                                };
                                i as f32
                            })
                            .collect()
                    },
                    heat: s.joule,
                    paths: std::mem::take(&mut self.fine).into_iter().collect(),
                    cells: s
                        .cells
                        .iter()
                        .zip(&self.layout.metals)
                        .map(|(c, m)| MetalCells {
                            flux: c.flux.iter().map(|f| [f[0] as f32, f[1] as f32]).collect(),
                            work: (if m.isolated { &c.induced } else { &c.work })
                                .iter()
                                .map(|w| *w as f32)
                                .collect(),
                        })
                        .collect(),
                };
                let metals = self.layout.metals.len();
                if metals > 0 {
                    self.pending.work.resize(metals, Vec::new());
                    for m in 0..metals {
                        let w = self.state.field_work.get(m).copied().unwrap_or(0.0);
                        self.pending.work[m].push((frame.t, w));
                    }
                    if self.state.book_potential_steps {
                        let steps = self.state.potential_steps.iter().sum();
                        self.pending.potential_steps.push((frame.t, steps));
                    }
                }
                if self.layout.has_field() && metals > 0 {
                    self.pending.middle.resize(metals, [Vec::new(), Vec::new()]);
                    for m in 0..metals {
                        let e = self.layout.middle_field(&frame, m);
                        self.pending.middle[m][0].push((frame.t, e.x));
                        self.pending.middle[m][1].push((frame.t, e.y));
                    }
                }
                self.pending.frames.push(frame);
            }
            // Each heat goal's window.
            self.heat_goals();
            let s = &self.state;
            // Each potential goal's measure over its window.
            if self.pending.potential.len() < self.potential.len() {
                self.pending.potential.resize(self.potential.len(), None);
            }
            for (k, g) in self.spec.potential_goals.iter().enumerate() {
                if s.t < g.start - 0.5 * self.tube.dt || self.potential[k].is_some() {
                    continue;
                }
                let v = s.potentials[g.electrode] - g.reference.map_or(0.0, |r| s.potentials[r]);
                let w = &mut self.windows[k];
                w.0 += v;
                w.1 += 1.0;
                w.2 = w.2.min(v);
                w.3 = w.3.max(v);
                if matches!(g.measure, Measure::Amplitude | Measure::Phase) {
                    self.series[k].push((s.t, v));
                }
                if s.t >= g.end - 0.5 * self.tube.dt {
                    let m = match g.measure {
                        Measure::Mean => w.0 / w.1.max(1.0),
                        Measure::Swing => 0.5 * (w.3 - w.2),
                        Measure::Amplitude => fourier_at(&self.series[k], g.omega).0,
                        Measure::Phase => fourier_at(&self.series[k], g.omega).1,
                    };
                    self.potential[k] = Some(m);
                    self.pending.potential[k] = Some(m);
                }
            }
        }
        self.pending.in_flight = self.state.x.len();
        self.pending.steps = self.steps;
    }

    /// Metal `m`'s heat so far, the sum over its pieces: the field's work on its carriers
    /// (`TubeState::field_work`), or, in a metal without contacts, the induced field's alone
    /// (`TubeState::induced_work`: there the electrostatic field's work is only the time
    /// step's numerical heating).
    fn metal_work(&self, m: usize) -> f64 {
        let isolated = self
            .spec
            .metals
            .get(m)
            .is_some_and(|x| x.contacts.is_empty());
        let work = if isolated {
            &self.state.induced_work
        } else {
            &self.state.field_work
        };
        self.lattice_metal
            .iter()
            .zip(work)
            .filter(|(k, _)| **k == m)
            .map(|(_, w)| w)
            .sum()
    }

    /// The heat goals' windows after a step: each one's start, and its mean power at its end.
    fn heat_goals(&mut self) {
        let (t, dt) = (self.state.t, self.tube.dt);
        if self.pending.heat.len() < self.heat.len() {
            self.pending.heat.resize(self.heat.len(), None);
        }
        for k in 0..self.spec.heat_goals.len() {
            let g = self.spec.heat_goals[k];
            if self.heat[k].is_some() {
                continue;
            }
            let w = self.metal_work(g.metal);
            if self.heat_start[k].is_none() && t >= g.start - 0.5 * dt {
                self.heat_start[k] = Some((t, w));
            }
            if t >= g.end - 0.5 * dt {
                let (t0, w0) = self.heat_start[k].unwrap_or((0.0, 0.0));
                let p = if t > t0 { (w - w0) / (t - t0) } else { 0.0 };
                self.heat[k] = Some(p);
                self.pending.heat[k] = Some(p);
            }
        }
    }

    /// The frames and samples made since the last call (the first chunk also carries the
    /// segments and potentials).
    pub fn take_chunk(&mut self) -> TubeRun {
        std::mem::take(&mut self.pending)
    }
}

impl TubeSim {
    /// The run's results so far.
    pub fn results(&self) -> TubeResults {
        TubeResults {
            currents: self.currents.clone(),
            potential: self.potential.clone(),
            heat: self.heat.clone(),
            arrivals: self.arrivals.clone(),
            capped: self.state.capped,
        }
    }
}
