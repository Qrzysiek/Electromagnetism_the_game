//! The fields: the surface charge the particles induce, the charges' fields by trees or
//! direct sums, the backgrounds and induced field added once (`with_external`), the
//! potential, the field at a point.

use super::*;

impl Tube {
    /// The particles' induced surface density (all electrodes grounded).
    pub(super) fn induced_density(&self, s: &TubeState) -> Vec<f64> {
        if s.x.len() >= TREE_MIN {
            let tree = Tree::new(s.x.iter().map(|p| Source::point(*p)).collect(), 0.0);
            let m = tree.multipoles(&s.charge);
            return self.electrodes.induced_by(|y| tree.potential(&m, y));
        }
        self.electrodes.induced_by(|y| {
            s.x.iter()
                .zip(&s.charge)
                .map(|(p, l)| line_charge(*p, *l, y).0)
                .sum()
        })
    }

    pub(super) fn segments_of(&self, electrode: usize) -> impl Iterator<Item = usize> + '_ {
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
    pub(super) fn add_induced(&self, induced: &mut [f64], s: &TubeState, from: usize) {
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

    /// Every particle's field from the surface charge `sigma` and every other particle, by
    /// the trees (`tree2d`: the direct sums to its series' truncation, 5e-7 of a node's field)
    /// or by the direct sums
    /// (`charges_field`); the step takes the trees from `TREE_MIN` particles on. Neither
    /// includes the backgrounds or a changing field's induced part: `with_external` adds them
    /// to both alike. (Until 2026-10-26 the direct sums included the backgrounds and the push
    /// added them again: below `TREE_MIN` they pulled twice, PHYSICS.md §2.11; test
    /// `trees_and_direct_sums_agree`.)
    pub(super) fn charges_fields(&self, s: &TubeState, sigma: &[f64], trees: bool) -> Vec<DVec3> {
        if trees {
            let segments = &self.electrodes.tree;
            let segment_m = segments.multipoles(sigma);
            let particles = Tree::new(
                s.x.iter().map(|p| Source::point(*p)).collect(),
                self.softening,
            );
            let particle_m = particles.multipoles(&s.charge);
            (0..s.x.len())
                .into_par_iter()
                .with_min_len(PAR_CHUNK)
                .map(|k| {
                    segments.field(&segment_m, s.x[k], None)
                        + particles.field(&particle_m, s.x[k], Some(k))
                })
                .collect()
        } else {
            (0..s.x.len())
                .into_par_iter()
                .with_min_len(PAR_CHUNK)
                .map(|k| self.charges_field(s, sigma, s.x[k], Some(k)))
                .collect()
        }
    }

    /// The whole field at `x` at time `t` from the charges' part `charges`: the lattices'
    /// backgrounds and a changing field's induced part added (no sum at all without a
    /// lattice: bit for bit as before them), and that induced part apart (its work is booked
    /// alone). The one place the push and the display (`field_at`) both take it from.
    pub(super) fn with_external(&self, charges: DVec3, t: f64, x: DVec3) -> (DVec3, Option<DVec3>) {
        let e = if self.lattices.is_empty() {
            charges
        } else {
            charges
                + self
                    .lattices
                    .iter()
                    .map(|l| l.background_at(x).1)
                    .sum::<DVec3>()
        };
        let induced = self.b_wave.map(|w| w.induced(t, x));
        let e = match induced {
            Some(i) => e + i,
            None => e,
        };
        (e, induced)
    }

    /// The field at `x` of the surface charge `sigma` and of every particle but `skip`, by
    /// direct sums (`charges_fields`' below `TREE_MIN`).
    pub(super) fn charges_field(
        &self,
        s: &TubeState,
        sigma: &[f64],
        x: DVec3,
        skip: Option<usize>,
    ) -> DVec3 {
        let mut e = DVec3::ZERO;
        for ((p, q), d) in self.electrodes.segments.iter().zip(sigma) {
            e += segment_field_fast(*p, *q, x) * *d;
        }
        let h2 = self.softening * self.softening;
        for (k, (p, l)) in s.x.iter().zip(&s.charge).enumerate() {
            if Some(k) == skip {
                continue;
            }
            let r = DVec3::new(x.x - p.x, x.y - p.y, 0.0);
            let r2 = (r.x * r.x + r.y * r.y).max(h2);
            e += r * (2.0 * l / r2);
        }
        e
    }

    /// The potential at `x` of the segments' charge, the particles and the backgrounds: the
    /// segments' exact integrals, each particle's rod (`−λ ln r²` outside its radius, `−λ (r²/h² − 1 + ln h²)` inside, whose
    /// field `charges_field` uses) and the lattices' backgrounds.
    pub(super) fn potential(
        &self,
        s: &TubeState,
        sigma: &[f64],
        x: DVec3,
        skip: Option<usize>,
    ) -> f64 {
        let mut phi = 0.0;
        for ((a, b), d) in self.electrodes.segments.iter().zip(sigma) {
            phi += segment_integrals(*a, *b, x).0 * *d;
        }
        let h2 = self.softening * self.softening;
        for (k, (p, l)) in s.x.iter().zip(&s.charge).enumerate() {
            if Some(k) == skip {
                continue;
            }
            let (dx, dy) = (x.x - p.x, x.y - p.y);
            let r2 = dx * dx + dy * dy;
            phi -= if r2 >= h2 {
                l * libm::log(r2)
            } else {
                l * (r2 / h2 - 1.0 + libm::log(h2))
            };
        }
        for l in &self.lattices {
            phi += l.background_at(x).0;
        }
        phi
    }

    /// The field at `x` now: the segments' charge, every particle's (softened) and the
    /// lattices' backgrounds, and a changing field's induced part.
    #[must_use]
    pub fn field_at(&self, s: &TubeState, x: DVec3) -> DVec3 {
        let sigma = self.surface_density(s);
        self.with_external(self.charges_field(s, &sigma, x, None), s.t, x)
            .0
    }
}
