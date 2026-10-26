//! The carriers' own magnetic field's part of a step (docs/DARWIN.md).

use super::*;

impl Tube {
    /// The self-field's part of a step (docs/DARWIN.md §2.1): the canonical velocities take
    /// the step's other impulses (`v − v_start`); the field is solved from them, the induced
    /// field implicitly through the carriers' response `β = q²/m`; the velocities follow as
    /// `v = u − (q/m) A(x)`, then the canonical force `(q/m) ∇(v·A)` from the field just
    /// solved at the same positions kicks both; the particles move. (Until 2026-10-26 the
    /// kick came first, from the last step's field, solved where the particles had been: each
    /// one's own vector potential then peaked a move behind it and dragged it, a loss growing
    /// as the step and the grid's sharpness, PHYSICS.md §2.12.)
    pub(super) fn darwin_step(&self, s: &mut TubeState, v_start: &[DVec3]) {
        let Some(mut d) = s.darwin.take() else {
            return;
        };
        for ((u, v), v0) in s.u.iter_mut().zip(&s.v).zip(v_start) {
            *u += *v - *v0;
        }
        let mut j = crate::darwin2d::Edges::zeros(&d);
        let mut beta = crate::darwin2d::Edges::zeros(&d);
        for k in 0..s.x.len() {
            let (x, u, q, qm) = (s.x[k], s.u[k], s.charge[k], s.qm[k]);
            d.deposit_vector(&mut j, x.x, x.y, q * u.x, q * u.y);
            d.deposit_scalar(&mut beta, x.x, x.y, q * qm);
        }
        d.solve(&j, &beta, DARWIN_TOLERANCE);
        let a = d.a();
        for k in 0..s.x.len() {
            let (x, qm) = (s.x[k], s.qm[k]);
            let (ax, ay) = d.interpolate(&a, x.x, x.y);
            let v = s.u[k] - DVec3::new(ax, ay, 0.0) * qm;
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
