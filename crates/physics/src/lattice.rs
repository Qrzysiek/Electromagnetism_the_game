//! A conductor's interior in the z-invariant world (PHYSICS.md §2.11, the Drude metals):
//! carriers in an axis-aligned rectangle over a uniform background of the opposite charge
//! (the ions' charge, jellium: its field in closed form, `zinv::rect_integrals`), bouncing
//! off fixed ion cores (hard discs placed at random: an amorphous metal, whose free paths
//! are exponential as the Boltzmann equation of the model assumes). A collision is a
//! specular bounce off the core with a new speed drawn from the lattice's temperature: the
//! lattice takes the energy the field gave (the Joule heat). The speed is drawn from
//! `∝ v² e^{−v²/2σ²}`, `σ² = kT/m`: the collision rate grows as v, so this is the speed
//! distribution of the carriers that collide in a Maxwellian gas, which stays Maxwellian
//! at rest. Random numbers come from a hash of (seed, carrier, step, collision):
//! reproducible, whatever order the carriers are handled in.

use glam::DVec3;

use crate::zinv::rect_integrals;

/// A conductor: its rectangle, background, cores and carriers.
#[derive(Clone, Debug, PartialEq)]
pub struct Lattice {
    pub lo: DVec3,
    pub hi: DVec3,
    /// The background's charge per area (the ions'; the carriers' is opposite).
    pub background: f64,
    /// The ion cores' centres and radius.
    pub cores: Vec<DVec3>,
    pub core_radius: f64,
    /// The carriers' charge-to-mass ratio and thermal speed `σ = √(kT/m)`.
    pub charge_per_mass: f64,
    pub thermal_speed: f64,
    /// Macroparticles at the start (together the background's charge, opposite).
    pub carriers: usize,
    pub seed: u64,
    /// The cores by cell (`cell` wide, `nx × ny` from `lo`).
    cell: f64,
    nx: usize,
    ny: usize,
    buckets: Vec<Vec<u32>>,
}

/// splitmix64: a well-mixed 64-bit hash.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// A uniform number in (0, 1] from a hash.
#[allow(clippy::cast_precision_loss)]
fn unit(h: u64) -> f64 {
    ((h >> 11) as f64 + 1.0) / (1u64 << 53) as f64
}

/// Two standard normals from two uniforms (Box–Muller).
fn normals(u1: f64, u2: f64) -> (f64, f64) {
    let r = libm::sqrt(-2.0 * libm::log(u1));
    let (s, c) = libm::sincos(std::f64::consts::TAU * u2);
    (r * c, r * s)
}

/// The most collisions a carrier has in one step (the rest of the step is then a free
/// move: at the time steps used, a third collision within one step is rare).
const MAX_COLLISIONS: u32 = 4;

/// A carrier's step through the cores (`Lattice::advance`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Advance {
    /// The position at the step's end and the velocity to store (the half step's, as the
    /// leapfrog push keeps it).
    pub x: DVec3,
    pub v: DVec3,
    /// The kinetic energy per mass the lattice took (`(v_in² − v_out²)/2` summed).
    pub heat_per_mass: f64,
    pub collisions: u32,
}

impl Lattice {
    /// A conductor (`lo`..`hi`), its background, cores and carriers; `cores` must lie in it.
    #[must_use]
    #[allow(
        clippy::too_many_arguments,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn new(
        lo: DVec3,
        hi: DVec3,
        background: f64,
        cores: Vec<DVec3>,
        core_radius: f64,
        charge_per_mass: f64,
        thermal_speed: f64,
        carriers: usize,
        seed: u64,
    ) -> Self {
        let area = (hi.x - lo.x) * (hi.y - lo.y);
        let spacing = libm::sqrt(area / cores.len().max(1) as f64);
        let cell = (2.0 * spacing).max(4.0 * core_radius).max(1e-9);
        let nx = ((hi.x - lo.x) / cell).ceil().max(1.0) as usize;
        let ny = ((hi.y - lo.y) / cell).ceil().max(1.0) as usize;
        let mut buckets = vec![Vec::new(); nx * ny];
        for (k, c) in cores.iter().enumerate() {
            let i = (((c.x - lo.x) / cell) as usize).min(nx - 1);
            let j = (((c.y - lo.y) / cell) as usize).min(ny - 1);
            buckets[i + nx * j].push(u32::try_from(k).unwrap_or(u32::MAX));
        }
        Self {
            lo,
            hi,
            background,
            cores,
            core_radius,
            charge_per_mass,
            thermal_speed,
            carriers,
            seed,
            cell,
            nx,
            ny,
            buckets,
        }
    }

    /// `count` cores placed independently and uniformly over `lo`..`hi`, at least `radius`
    /// from the walls (they may overlap: the Boolean model, whose free paths are exactly
    /// exponential with the rate `2 a n v`; a lattice, even jittered, is more regular, and
    /// one jittered less than a cell keeps open corridors, free paths without end).
    #[must_use]
    pub fn random_cores(lo: DVec3, hi: DVec3, count: usize, radius: f64, seed: u64) -> Vec<DVec3> {
        (0..count as u64)
            .map(|k| {
                let h = mix(seed ^ mix(k ^ 0x000c_04e5));
                DVec3::new(
                    lo.x + radius + unit(h) * (hi.x - lo.x - 2.0 * radius),
                    lo.y + radius + unit(mix(h)) * (hi.y - lo.y - 2.0 * radius),
                    0.0,
                )
            })
            .collect()
    }

    /// Whether `x` lies in the conductor (its edges included).
    #[must_use]
    pub fn contains(&self, x: DVec3) -> bool {
        x.x >= self.lo.x && x.x <= self.hi.x && x.y >= self.lo.y && x.y <= self.hi.y
    }

    /// The background's potential and field at `x`.
    #[must_use]
    pub fn background_at(&self, x: DVec3) -> (f64, DVec3) {
        let (p, e) = rect_integrals(self.lo, self.hi, x);
        (self.background * p, e * self.background)
    }

    /// Each carrier's charge at the start: the background's over them, opposite.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn carrier_charge(&self) -> f64 {
        -self.background * (self.hi.x - self.lo.x) * (self.hi.y - self.lo.y)
            / self.carriers.max(1) as f64
    }

    /// The first core the straight move from `a` to `b` enters: the fraction of the move
    /// and the core's centre.
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_possible_wrap,
        clippy::cast_precision_loss
    )]
    pub fn first_hit(&self, a: DVec3, b: DVec3) -> Option<(f64, DVec3)> {
        let d = b - a;
        let r = self.core_radius;
        let cell = |v: f64, lo: f64, n: usize| {
            ((v - lo) / self.cell).floor().clamp(0.0, (n - 1) as f64) as usize
        };
        let (i0, i1) = (
            cell(a.x.min(b.x) - r, self.lo.x, self.nx),
            cell(a.x.max(b.x) + r, self.lo.x, self.nx),
        );
        let (j0, j1) = (
            cell(a.y.min(b.y) - r, self.lo.y, self.ny),
            cell(a.y.max(b.y) + r, self.lo.y, self.ny),
        );
        let dd = d.length_squared();
        let mut best: Option<(f64, DVec3)> = None;
        for j in j0..=j1 {
            for i in i0..=i1 {
                for &k in &self.buckets[i + self.nx * j] {
                    let c = self.cores[k as usize];
                    let m = a - c;
                    // |m + s d|² = r²: s² dd + 2 s (m·d) + |m|² − r² = 0, entering (m·d < 0).
                    let (bh, cc) = (m.dot(d), m.length_squared() - r * r);
                    if cc < 0.0 || bh >= 0.0 || dd == 0.0 {
                        continue;
                    }
                    let disc = bh * bh - dd * cc;
                    if disc < 0.0 {
                        continue;
                    }
                    let s = (-bh - libm::sqrt(disc)) / dd;
                    if (0.0..=1.0).contains(&s) && best.is_none_or(|(t, _)| s < t) {
                        best = Some((s, c));
                    }
                }
            }
        }
        best
    }

    /// A new speed after a collision (carrier `id`'s `k`th at step `step`):
    /// `∝ v² e^{−v²/2σ²}`, the length of three normal components.
    #[must_use]
    pub fn draw_speed(&self, id: u32, step: u64, k: u32) -> f64 {
        let h = mix(self.seed ^ mix(u64::from(id) ^ mix(step ^ (u64::from(k) << 56))));
        let (h2, h3, h4) = (mix(h), mix(mix(h)), mix(mix(mix(h))));
        let (g1, g2) = normals(unit(h), unit(h2));
        let (g3, _) = normals(unit(h3), unit(h4));
        self.thermal_speed * libm::sqrt(g1 * g1 + g2 * g2 + g3 * g3)
    }

    /// A carrier's step (carrier `id` at step `step`, of `dt`) from `a`, its stored
    /// velocity `v0` before the step's kick and `v1` after it: the leapfrog push moves it
    /// to `a + v1 dt`; here the move is followed through the cores it enters (at most
    /// `MAX_COLLISIONS`). The kick `g = v1 − v0` is spread over the step (the velocity at
    /// the fraction φ of it is `v1 + g (φ − ½)`, the path a parabola), so a collision
    /// reflects the velocity of its moment, the part of the kick after it is kept, and the
    /// field's work up to it equals the kinetic energy gained. Each collision: a specular
    /// bounce off the core (found on the chord of the path, with the core's normal there;
    /// the carrier placed on the path: they differ by `g dt/8`), a new thermal speed.
    /// `None` without a collision (the push's move stands).
    #[must_use]
    pub fn advance(
        &self,
        a: DVec3,
        v0: DVec3,
        v1: DVec3,
        dt: f64,
        id: u32,
        step: u64,
    ) -> Option<Advance> {
        let g = v1 - v0;
        let (mut x, mut u, mut phi) = (a, v1 - g * 0.5, 0.0);
        let (mut heat, mut k) = (0.0, 0);
        loop {
            let rest = 1.0 - phi;
            let end = x + (u + g * (0.5 * rest)) * (rest * dt);
            let hit = if k < MAX_COLLISIONS {
                self.first_hit(x, end)
            } else {
                None
            };
            let Some((f, c)) = hit else {
                if k == 0 {
                    return None;
                }
                return Some(Advance {
                    x: end,
                    v: u + g * (rest - 0.5),
                    heat_per_mass: heat,
                    collisions: k,
                });
            };
            // The hit found on the chord, its normal there (on the surface: a normal taken
            // off it would lean against the field, a bias of `g dt/8a`), the carrier placed
            // on the path (a parabola: the field's work up to it is then exactly the
            // kinetic energy gained).
            let span = f * rest;
            let n = (x + (end - x) * f - c).normalize_or_zero();
            let at = x + (u + g * (0.5 * span)) * (span * dt);
            let u_hit = u + g * span;
            let dir = (u_hit - n * (2.0 * u_hit.dot(n))).normalize_or_zero();
            let speed = self.draw_speed(id, step, k);
            heat += 0.5 * (u_hit.length_squared() - speed * speed);
            u = dir * speed;
            x = at + n * (1e-9 * self.core_radius);
            phi += span;
            k += 1;
        }
    }

    /// The carriers at the start: on a jittered grid over the conductor (a quiet start),
    /// out of the cores, with Maxwellian velocities (each component normal, σ).
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn fill(&self) -> Vec<(DVec3, DVec3)> {
        let (w, h) = (self.hi.x - self.lo.x, self.hi.y - self.lo.y);
        let n = self.carriers.max(1);
        let nx = libm::round(libm::sqrt(n as f64 * w / h)).max(1.0) as usize;
        let ny = n.div_ceil(nx);
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            let (i, j) = (k % nx, k / nx);
            let hk = mix(self.seed ^ mix(0x5eed ^ k as u64));
            let mut x = DVec3::new(
                self.lo.x + (i as f64 + unit(hk)) / nx as f64 * w,
                self.lo.y + (j as f64 + unit(mix(hk))) / ny as f64 * h,
                0.0,
            );
            // Out of a core it fell in: to its rim, outward.
            if let Some(c) = self
                .cores
                .iter()
                .find(|c| (x - **c).length() < self.core_radius)
            {
                let dir = (x - *c).normalize_or(DVec3::X);
                x = *c + dir * (self.core_radius * 1.000_001);
            }
            let (g1, g2) = normals(unit(mix(mix(hk))), unit(mix(mix(mix(hk)))));
            out.push((x, DVec3::new(g1, g2, 0.0) * self.thermal_speed));
        }
        out
    }
}
