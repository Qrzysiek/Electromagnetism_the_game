//! Particles added and given off: the emitters' space-charge-limited emission, the ohmic
//! contacts' carriers.

use super::*;

impl Tube {
    pub(super) fn add(
        &self,
        s: &mut TubeState,
        x: DVec3,
        v: DVec3,
        charge: f64,
        qm: f64,
        tag: u32,
    ) {
        // With the self-field: born at v, its current fixed at `q v` in the next solve, which
        // sets its canonical velocity (`darwin_step`); until then `v + (q/m) A` from the last.
        // (Until 2026-10-27 `u = v + (q/m) A`, then with a uniform medium's self-inductance
        // added: a contact's carriers were born slow, and the lattice reheated them.)
        if let (Some(d), Some(a)) = (s.darwin.as_ref(), s.self_a.as_ref()) {
            let (ax, ay) = d.interpolate(a, x.x, x.y);
            s.u.push(v + DVec3::new(ax, ay, 0.0) * qm);
            s.set.push(VelocitySet::Born);
        }
        s.x.push(x);
        s.v.push(v);
        s.charge.push(charge);
        s.qm.push(qm);
        s.tag.push(tag);
        s.id.push(s.next_id);
        s.next_id = s.next_id.wrapping_add(1);
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
    pub(super) fn outward(&self, electrode: usize, a: DVec3, b: DVec3, at: DVec3) -> DVec3 {
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
    pub(super) fn emit(&self, s: &mut TubeState, induced: &[f64]) -> Vec<f64> {
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
            if let Some(ohmic) = em.ohmic {
                self.emit_ohmic(s, &sigma, i, em, ohmic, &segs, &mut emitted);
                continue;
            }
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
                    let v = if em.thermal_speed > 0.0 {
                        crate::lattice::flux_velocity(out, em.thermal_speed, s.emitted)
                    } else {
                        DVec3::ZERO
                    };
                    self.add(s, x, v, sign * em.weight, em.charge_per_mass, tag);
                    emitted[em.electrode] += sign * em.weight;
                    s.released[i] += sign * em.weight;
                }
            }
        }
        emitted
    }

    /// An ohmic contact's emission (`OhmicContact`): per emitting segment of emitter `i`,
    /// the carriers its layer lacks (`sigma`: the segments' density now, for the potential
    /// step of each release).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn emit_ohmic(
        &self,
        s: &mut TubeState,
        sigma: &[f64],
        i: usize,
        em: &Emitter,
        ohmic: OhmicContact,
        segs: &[usize],
        emitted: &mut [f64],
    ) {
        let sign = em.charge_per_mass.signum();
        let tag = u32::try_from(i).unwrap_or(u32::MAX);
        // The other ohmic contacts' layers (each emitting segment's: its start, direction,
        // outward normal, length and depth), to share where they overlap this one's.
        let others = self.ohmic_layers(i);
        for &j in segs {
            let (a, b) = self.electrodes.segments[j];
            let len = (b - a).length();
            let t = (b - a) / len;
            let out = self.outward(em.electrode, a, b, (a + b) * 0.5);
            if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) {
                continue;
            }
            let layer = Layer {
                a,
                t,
                out,
                len,
                depth: ohmic.depth,
            };
            // The other contacts' layers overlapping this one (two contacts meeting at a
            // metal's corner): a carrier there is held by the face it is nearest. (Until
            // 2026-11-04 by both: the contact that emitted first each step filled the shared
            // corner and the other never did, a current between two grounded contacts; a
            // Corbino disc's four outer contacts took 17 % apart, D12.)
            let near: Vec<&Layer> = others.iter().filter(|o| o.overlaps(&layer)).collect();
            // Its share of its layer: the whole, or where it holds the carriers (a fixed
            // grid of 16 × 16 points; a corner's square split along its diagonal: half).
            let share = if near.is_empty() {
                1.0
            } else {
                layer.share(&near)
            };
            // The carriers of its sign in the layer now.
            let held: f64 =
                s.x.iter()
                    .zip(&s.charge)
                    .zip(&s.tag)
                    .filter(|((_, q), tag)| **tag < PROJECTILE_TAG && **q * sign > 0.0)
                    .filter(|((x, _), _)| {
                        layer.depth_of(**x).is_some_and(|w| {
                            near.iter().all(|o| o.depth_of(**x).is_none_or(|v| w <= v))
                        })
                    })
                    .map(|((_, q), _)| q.abs())
                    .sum();
            let mut lack = if near.is_empty() {
                ohmic.density * len * ohmic.depth - held
            } else {
                ohmic.density * len * ohmic.depth * share - held
            };
            while lack >= em.weight {
                if s.x.len() >= MAX_PARTICLES {
                    s.capped = true;
                    break;
                }
                lack -= em.weight;
                s.emitted += 1;
                #[allow(clippy::cast_precision_loss)]
                let f = (s.emitted as f64 * 0.618_033_988_749_894_9).fract();
                // From the face, as a gas at the lattice's temperature crosses a wall: the
                // flux-weighted law, 1.5 kT each in the plane, what the carriers the contact
                // takes in carry on average.
                let v = crate::lattice::flux_velocity(out, em.thermal_speed, s.emitted);
                let face = a + t * (len * f);
                let x = face + out * (EMISSION_OFFSET * len);
                // In the currents' accounts it comes from the face, as a carrier the contact
                // takes back short of it goes on to the face (`Tube::step`).
                if let Some(m) = self.metal_of(x) {
                    self.lattices[m].book_flux(&mut s.cells[m].flux, face, x, sign * em.weight);
                }
                if s.book_potential_steps {
                    let phi = self.potential(s, sigma, x, None);
                    s.potential_steps[em.electrode] +=
                        sign * em.weight * (s.potentials[em.electrode] - phi);
                }
                self.add(s, x, v, sign * em.weight, em.charge_per_mass, tag);
                emitted[em.electrode] += sign * em.weight;
                s.released[i] += sign * em.weight;
                s.contact_energy[em.electrode] -=
                    0.5 * (em.weight / em.charge_per_mass).abs() * v.length_squared();
            }
        }
    }

    /// The layers of every ohmic contact's emitting segments but emitter `skip`'s
    /// (`emit_ohmic`).
    fn ohmic_layers(&self, skip: usize) -> Vec<Layer> {
        let mut layers = Vec::new();
        for (k, em) in self.emitters.iter().enumerate() {
            let Some(ohmic) = em.ohmic.filter(|_| k != skip) else {
                continue;
            };
            for j in self.segments_of(em.electrode) {
                let (a, b) = self.electrodes.segments[j];
                let len = (b - a).length();
                let out = self.outward(em.electrode, a, b, (a + b) * 0.5);
                if em.emit_toward.is_some_and(|d| out.dot(d) < 0.5) {
                    continue;
                }
                layers.push(Layer {
                    a,
                    t: (b - a) / len,
                    out,
                    len,
                    depth: ohmic.depth,
                });
            }
        }
        layers
    }
}

/// An ohmic contact's layer in front of one emitting segment (`OhmicContact`): the points
/// `a + u t + w out` with `0 ≤ u ≤ len`, `0 < w ≤ depth`.
struct Layer {
    a: DVec3,
    t: DVec3,
    out: DVec3,
    len: f64,
    depth: f64,
}

impl Layer {
    /// How far in front of the face `x` lies, if it is in the layer.
    fn depth_of(&self, x: DVec3) -> Option<f64> {
        let r = x - self.a;
        let (u, w) = (r.dot(self.t), r.dot(self.out));
        ((0.0..=self.len).contains(&u) && w > 0.0 && w <= self.depth).then_some(w)
    }

    /// The part of it nearer its own face than any of `near`'s layers that hold the point:
    /// over a fixed grid of 16 × 16 points (the cells' middles), the fraction it holds.
    fn share(&self, near: &[&Self]) -> f64 {
        const N: u32 = 16;
        let mut held = 0u32;
        for i in 0..N {
            for k in 0..N {
                let u = self.len * (f64::from(i) + 0.5) / f64::from(N);
                let w = self.depth * (f64::from(k) + 0.5) / f64::from(N);
                let x = self.a + self.t * u + self.out * w;
                if near.iter().all(|o| o.depth_of(x).is_none_or(|v| w <= v)) {
                    held += 1;
                }
            }
        }
        f64::from(held) / f64::from(N * N)
    }

    /// Its corners' bounding box.
    fn bounds(&self) -> (DVec3, DVec3) {
        let corners = [
            self.a,
            self.a + self.t * self.len,
            self.a + self.out * self.depth,
            self.a + self.t * self.len + self.out * self.depth,
        ];
        corners.iter().fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(lo, hi), c| (lo.min(*c), hi.max(*c)),
        )
    }

    /// Whether its bounding box and `other`'s overlap with some area (layers of contacts
    /// that meet at a corner; a shared edge alone is not an overlap).
    fn overlaps(&self, other: &Self) -> bool {
        let ((a0, a1), (b0, b1)) = (self.bounds(), other.bounds());
        a0.x < b1.x && b0.x < a1.x && a0.y < b1.y && b0.y < a1.y
    }
}
