//! A metal's carriers after a move: through its cores, off its walls and cuts, into its
//! contacts.

use super::*;

impl Tube {
    /// The metal a carrier at `x` is in: the one whose conductor holds it (`Lattice::holds`;
    /// a metal may lie in another's hole), else the first whose rectangle does (one a rounding
    /// left just past its walls stays with its metal, as before). (Until 2026-10-27 the first
    /// rectangle alone: a ring's carriers in another's hole took that one, no cores and none of
    /// their own walls; test `nested_rings_keep_their_carriers`.)
    pub(super) fn metal_of(&self, x: DVec3) -> Option<usize> {
        self.lattices
            .iter()
            .position(|l| l.holds(x))
            .or_else(|| self.lattices.iter().position(|l| l.contains(x)))
    }

    /// The carrier `p` of `lattice` after a move from `a` (its velocity `v0` before the
    /// step's kick): through the cores it entered (`Lattice::advance`: bounces, new
    /// thermal speeds, the lattice taking the energy), then the walls (reflecting, but
    /// where an electrode lies beyond: a contact, which absorbs the carrier below).
    /// Returns where it ends and its velocity, its heat and collisions (when it met cores)
    /// and its velocity's mark for the self-field (when the jump set one): the state is the
    /// caller's to change (`Tube::step`, in the carriers' order).
    pub(super) fn scattered(
        &self,
        k: usize,
        s: &TubeState,
        p: usize,
        a: DVec3,
        v0: DVec3,
        step: u64,
    ) -> Scattered {
        let lattice = &self.lattices[k];
        let v_entry = s.v[p];
        let (mut at, mut vel) = (s.x[p], s.v[p]);
        let (mut heat, mut collisions) = (None, 0);
        if let Some(m) = lattice.advance(a, v0, vel, self.dt, s.id[p], step) {
            let mass = (s.charge[p] / s.qm[p]).abs();
            heat = Some(mass * m.heat_per_mass);
            collisions = u64::from(m.collisions);
            at = m.x;
            vel = m.v;
        }
        // Off its sides (out of the rectangle, not into a contact) and off its internal walls
        // (the step's chord from `from` crossing one), until it is inside.
        let mut from = a;
        for _ in 0..8 {
            let x = at;
            if lattice.contains(x) {
                let Some((t, w0, w1)) = lattice
                    .walls
                    .iter()
                    .filter_map(|&(w0, w1)| crossing(from, x, w0, w1).map(|t| (t, w0, w1)))
                    .min_by(|a, b| a.0.total_cmp(&b.0))
                else {
                    break;
                };
                // Specular: the position and the velocity mirrored in the wall's line.
                let e = w1 - w0;
                let n = DVec3::new(-e.y, e.x, 0.0) / libm::sqrt(e.x * e.x + e.y * e.y);
                let hit = from + (x - from) * t;
                let mirrored = x - n * (2.0 * (x - w0).dot(n));
                let v = vel;
                vel = v - n * (2.0 * v.dot(n));
                at = mirrored;
                // The next chord starts just off the wall, on the carrier's side.
                from = hit + (mirrored - hit) * 1e-6;
                continue;
            }
            if self
                .electrodes
                .electrodes
                .iter()
                .any(|e| inside(&e.section, x))
            {
                break;
            }
            let (lo, hi) = (lattice.lo, lattice.hi);
            let (mut x, mut v) = (x, vel);
            if x.x < lo.x {
                x.x = 2.0 * lo.x - x.x;
                v.x = -v.x;
            } else if x.x > hi.x {
                x.x = 2.0 * hi.x - x.x;
                v.x = -v.x;
            }
            if x.y < lo.y {
                x.y = 2.0 * lo.y - x.y;
                v.y = -v.y;
            } else if x.y > hi.y {
                x.y = 2.0 * hi.y - x.y;
                v.y = -v.y;
            }
            at = x;
            vel = v;
        }
        // With the self-field, the jump the collisions or the walls made: the next solve
        // fixes the carrier's current at the velocity it would have had without it, plus the
        // jump, the other carriers responding (`darwin_step`); `u` keeps the velocity before
        // it until then. (Until 2026-10-27 `u` took the jump `Δv` with a uniform medium's
        // self-inductance g added, `Δv + q(q/m) g Δv`: the metal's edges, unscreened, make
        // the true one larger, 13 % on average in Y25's disc at c = 3, and each collision
        // lost the difference; before that `Δv` alone, losing `½ q² Δv·gΔv`.)
        let set = (!s.u.is_empty() && vel != v_entry).then(|| {
            let dv = vel - v_entry;
            match s.set[p] {
                VelocitySet::Born => VelocitySet::Born,
                VelocitySet::Jumped(before) => VelocitySet::Jumped(before + dv),
                VelocitySet::Free => VelocitySet::Jumped(dv),
            }
        });
        Scattered {
            x: at,
            v: vel,
            heat,
            collisions,
            set,
        }
    }
}

/// A carrier's scattering over a step (`Tube::scattered`).
pub(super) struct Scattered {
    pub x: DVec3,
    pub v: DVec3,
    /// The heat it gave the lattice (`Some` when it met cores, even 0) and its collisions.
    pub heat: Option<f64>,
    pub collisions: u64,
    /// Its velocity's new mark for the self-field, if the step's jumps changed it.
    pub set: Option<VelocitySet>,
}
