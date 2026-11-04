//! The carriers' own magnetic field's part of a step (docs/DARWIN.md).

use super::*;
use crate::darwin2d::{Darwin2d, Edges};

impl Tube {
    /// The self-field's part of a step (docs/DARWIN.md §2.1): the canonical velocities take
    /// the step's other impulses (`v − v_start`); the field is solved from them, the induced
    /// field implicitly through the carriers' response `β = q²/m`; the velocities follow as
    /// `v = u − (q/m) A(x)`, then the canonical force `(q/m) ∇(v·A)` from the field just
    /// solved at the same positions kicks both; the particles move. A particle whose velocity
    /// was set since the last solve (`VelocitySet`) enters it with its current fixed, not
    /// responding, its `u` following from the solution as `v + (q/m) A(x)`: one given off at
    /// the velocity it was given; one a collision or a wall changed at the velocity it would
    /// have had without the jump, plus the jump, from a first solve with every carrier
    /// responding and its `u` the one before the jump (the other carriers' induction, an
    /// inductor's back-EMF, is then felt as by the rest). (Until 2026-10-26 the kick came
    /// first, from the last step's field, solved where the particles had been: each one's own
    /// vector potential then peaked a move behind it and dragged it, a loss growing as the
    /// step and the grid's sharpness, PHYSICS.md §2.12.)
    pub(super) fn darwin_step(&self, s: &mut TubeState, v_start: &[DVec3]) {
        let Some(mut d) = s.darwin.take() else {
            return;
        };
        for (((u, v), v0), set) in s.u.iter_mut().zip(&s.v).zip(v_start).zip(&s.set) {
            if *set != VelocitySet::Born {
                *u += *v - *v0;
            }
        }
        if s.set.iter().any(|v| matches!(v, VelocitySet::Jumped(_))) {
            let (j, beta) = deposit(&d, s);
            d.solve(&j, &beta, DARWIN_TOLERANCE);
            let a = d.a();
            for k in 0..s.x.len() {
                if let VelocitySet::Jumped(dv) = s.set[k] {
                    let (ax, ay) = d.interpolate(&a, s.x[k].x, s.x[k].y);
                    s.v[k] = s.u[k] - DVec3::new(ax, ay, 0.0) * s.qm[k] + dv;
                    s.set[k] = VelocitySet::Born;
                }
            }
        }
        let (j, beta) = deposit(&d, s);
        d.solve(&j, &beta, DARWIN_TOLERANCE);
        let a = d.a();
        if let Some(leads) = s.leads.as_mut() {
            leads.update_rest(&d, &a);
        }
        // Each carrier's kick and move, in parallel (its own arithmetic).
        let dt = self.dt;
        s.x.par_iter_mut()
            .zip(&mut s.v)
            .zip(&mut s.u)
            .zip(&mut s.set)
            .zip(&s.qm)
            .for_each(|((((x, v), u), set), &qm)| {
                let (ax, ay) = d.interpolate(&a, x.x, x.y);
                let vel = if *set == VelocitySet::Free {
                    *u - DVec3::new(ax, ay, 0.0) * qm
                } else {
                    *set = VelocitySet::Free;
                    *u = *v + DVec3::new(ax, ay, 0.0) * qm;
                    *v
                };
                let (gx, gy) = d.grad_dot(&a, x.x, x.y, vel.x, vel.y);
                let kick = DVec3::new(gx, gy, 0.0) * (qm * dt);
                *u += kick;
                *v = vel + kick;
                *x += *v * dt;
            });
        s.self_a = Some(a);
        s.darwin = Some(d);
    }
}

/// A metal's carrier's share of the field's work over a step with the self-field, booked
/// once the step's velocities are known (`Tube::book_work`).
pub(super) struct PendingWork {
    pub particle: usize,
    pub lattice: usize,
    pub cell: usize,
    pub charge: f64,
    /// Its mass, `q/(q/m)`.
    pub mass: f64,
    /// The field at the step's start (the push's), and its induced part.
    pub field: DVec3,
    pub induced: Option<DVec3>,
    /// The velocity at the step's start, and after the push (before the self-field's part),
    /// and the position at the step's start.
    pub old: DVec3,
    pub pushed: DVec3,
    pub at: DVec3,
}

/// The averaging window of the macroscopic currents of `TubeState::mutual_work`: short
/// against a drive's period (its central differences' error `(ωT)²/6`, 0.4 % at ω = 3), long
/// against a step (the macroparticles' grid-crossing noise averaged out).
pub const MUTUAL_WINDOW: f64 = 0.05;

/// The macroscopic currents' accounts (`Tube::book_mutual`).
#[derive(Clone, Debug, Default)]
pub(super) struct MutualAccount {
    /// Per lattice, this window's sum of its carriers' currents (per area, on the edges),
    /// and the window's steps.
    sums: Vec<Edges>,
    steps: usize,
    /// The last closed window's mean currents and fields, and the fields of the one before
    /// (per lattice).
    current: Vec<Edges>,
    field: Vec<Edges>,
    before: Vec<Edges>,
    /// Per lattice, each x-edge's and y-edge's account cell (`Lattice::account_cell` of its
    /// middle).
    cells: Vec<(Vec<usize>, Vec<usize>)>,
}

impl Tube {
    /// The field's work on the metals' carriers over a step with the self-field:
    /// `(q dt/2) E·(v_old + v_new)` with `v_new` the velocity after the self-field's part
    /// (`darwin_step`), the kinetic current the field drives. Not the push's own `v_new`:
    /// the induction then takes back part of the push's impulse (all of it in a screened
    /// metal), and the push's velocity would book `(q²/2m) E·(∂A/∂t) dt` too much per carrier
    /// and time, ~`dt/2τ` of a screened metal's unscreened Joule heat (the self-field's part
    /// of the step books it back: the energy's account holds either way; until 2026-11-02,
    /// 8 % of the unscreened heat in excess in a ring 20 × 20 at dt 0.012, docs/HANDOFF.md).
    ///
    /// And the induced fields' work, the drive's and the carriers' own (`inductive_work`):
    /// the drive's at the push, `(q dt/2) E_ind·(v_old + v_push)`, the push's exact energy
    /// account of it, and the self-field's part of the step's change of the kinetic energy,
    /// `½ m (v_new² − v_push²)`. Together the carriers' energy gained from the induced fields
    /// over the step, to round-off (the electrostatic field's part of the push, its numerical
    /// heating, left out): the induction that takes back the push's impulse is the
    /// self-field's negative work, booked with it.
    pub(super) fn book_work(
        &self,
        s: &mut TubeState,
        pending: &[PendingWork],
        work: &mut [f64],
        induced_work: &mut [f64],
        inductive_work: &mut [f64],
    ) {
        for w in pending {
            let v = s.v[w.particle];
            let mean = w.old + v;
            let done = 0.5 * self.dt * w.charge * w.field.dot(mean);
            work[w.lattice] += done;
            s.cells[w.lattice].work[w.cell] += done;
            let mut inductive = 0.5 * w.mass * (v.length_squared() - w.pushed.length_squared());
            if let Some(i) = w.induced {
                let done = 0.5 * self.dt * w.charge * i.dot(mean);
                induced_work[w.lattice] += done;
                s.cells[w.lattice].induced[w.cell] += done;
                inductive += 0.5 * self.dt * w.charge * i.dot(w.old + w.pushed);
            }
            inductive_work[w.lattice] += inductive;
            s.cells[w.lattice].inductive[w.cell] += inductive;
        }
    }
}

impl Tube {
    /// The power the lattices' macroscopic currents pass to each other
    /// (`TubeState::mutual_work`), with two lattices or more: each lattice's carriers' `q v`
    /// at the middle of their step deposited on the self-field's grid and averaged over a
    /// window of `MUTUAL_WINDOW` (`J̄`), its field alone in the vacuum (`Ā`,
    /// `Darwin2d::vacuum_field`); when window m + 1 closes, lattice k takes from each other
    /// lattice j `−∫ J̄_kᵐ·(Ā_jᵐ⁺¹ − Ā_jᵐ⁻¹)/2 dA` for window m, edge by edge into the cell
    /// of the edge's middle.
    pub(super) fn book_mutual(&self, s: &mut TubeState, pending: &[PendingWork]) {
        let n = self.lattices.len();
        if n < 2 {
            return;
        }
        let Some(d) = s.darwin.take() else {
            return;
        };
        let mut acc = s.mutual.take().unwrap_or_default();
        if acc.sums.is_empty() {
            let (xs, ys) = d.edge_positions();
            acc.sums = (0..n).map(|_| Edges::zeros(&d)).collect();
            acc.cells = self
                .lattices
                .iter()
                .map(|l| {
                    let cell = |p: &(f64, f64)| l.account_cell(DVec3::new(p.0, p.1, 0.0));
                    (xs.iter().map(cell).collect(), ys.iter().map(cell).collect())
                })
                .collect();
        }
        for w in pending {
            let x = (w.at + s.x[w.particle]) * 0.5;
            let v = s.v[w.particle] * w.charge;
            d.deposit_vector(&mut acc.sums[w.lattice], x.x, x.y, v.x, v.y);
        }
        acc.steps += 1;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let window = libm::round(MUTUAL_WINDOW / self.dt).max(1.0) as usize;
        if acc.steps >= window {
            #[allow(clippy::cast_precision_loss)]
            let per = 1.0 / acc.steps as f64;
            let mean: Vec<Edges> = acc
                .sums
                .iter()
                .map(|e| Edges {
                    x: e.x.iter().map(|v| v * per).collect(),
                    y: e.y.iter().map(|v| v * per).collect(),
                })
                .collect();
            let field: Vec<Edges> = mean.iter().map(|j| d.vacuum_field(j)).collect();
            if !acc.before.is_empty() {
                let area = d.h * d.h;
                for k in 0..n {
                    let (cx, cy) = &acc.cells[k];
                    let jk = &acc.current[k];
                    for j in (0..n).filter(|j| *j != k) {
                        let (now, then) = (&field[j], &acc.before[j]);
                        let mut total = 0.0;
                        for (e, c) in cx.iter().enumerate() {
                            let w = -0.5 * area * jk.x[e] * (now.x[e] - then.x[e]);
                            s.cells[k].mutual[*c] += w;
                            total += w;
                        }
                        for (e, c) in cy.iter().enumerate() {
                            let w = -0.5 * area * jk.y[e] * (now.y[e] - then.y[e]);
                            s.cells[k].mutual[*c] += w;
                            total += w;
                        }
                        s.mutual_work[k] += total;
                    }
                }
            }
            if !acc.current.is_empty() {
                acc.before = std::mem::take(&mut acc.field);
            }
            acc.field = field;
            acc.current = mean;
            for e in &mut acc.sums {
                e.x.fill(0.0);
                e.y.fill(0.0);
            }
            acc.steps = 0;
        }
        s.mutual = Some(acc);
        s.darwin = Some(d);
    }
}

/// The current and the response on the grid's edges: a free particle's canonical current
/// `q u` with its response `q²/m` (the solve's implicit induction), a jumped one's likewise (its
/// `u` the one before the jump), a given-off one's current `q v`, fixed; the leads' currents,
/// fixed (their circuits took their induction, `leads.rs`).
fn deposit(d: &Darwin2d, s: &TubeState) -> (Edges, Edges) {
    let mut j = Edges::zeros(d);
    let mut beta = Edges::zeros(d);
    for k in 0..s.x.len() {
        let (x, q, qm) = (s.x[k], s.charge[k], s.qm[k]);
        if s.set[k] == VelocitySet::Born {
            let v = s.v[k];
            d.deposit_vector(&mut j, x.x, x.y, q * v.x, q * v.y);
        } else {
            let u = s.u[k];
            d.deposit_vector(&mut j, x.x, x.y, q * u.x, q * u.y);
            d.deposit_scalar(&mut beta, x.x, x.y, q * qm);
        }
    }
    if let Some(leads) = &s.leads {
        leads.deposit(d, &mut j);
    }
    (j, beta)
}
