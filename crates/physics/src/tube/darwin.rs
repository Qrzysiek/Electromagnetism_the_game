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
        for k in 0..s.x.len() {
            let (x, qm) = (s.x[k], s.qm[k]);
            let (ax, ay) = d.interpolate(&a, x.x, x.y);
            let v = if s.set[k] == VelocitySet::Free {
                s.u[k] - DVec3::new(ax, ay, 0.0) * qm
            } else {
                s.set[k] = VelocitySet::Free;
                s.u[k] = s.v[k] + DVec3::new(ax, ay, 0.0) * qm;
                s.v[k]
            };
            let (gx, gy) = d.grad_dot(&a, x.x, x.y, v.x, v.y);
            let kick = DVec3::new(gx, gy, 0.0) * (qm * self.dt);
            s.u[k] += kick;
            s.v[k] = v + kick;
            s.x[k] = x + s.v[k] * self.dt;
        }
        s.self_a = Some(a);
        s.darwin = Some(d);
    }
}

/// The current and the response on the grid's edges: a free particle's canonical current
/// `q u` with its response `q²/m` (the solve's implicit induction), a jumped one's likewise (its
/// `u` the one before the jump), a given-off one's current `q v`, fixed.
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
    (j, beta)
}
