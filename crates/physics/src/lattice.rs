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
    /// Insulating walls inside it (straight segments: cuts), which its carriers cannot cross
    /// (`physics::tube`: reflected off them as off its sides). Empty: none.
    pub walls: Vec<(DVec3, DVec3)>,
    /// Holes (`with_hole`): rectangles inside it with no cores, no carriers and no
    /// background, their edges walls (a ring: an inductor's turn). Empty: none.
    pub holes: Vec<(DVec3, DVec3)>,
    /// A round conductor (`with_disc`): the disc's centre and radius, inside the rectangle;
    /// its outline walls. None: the rectangle.
    pub disc: Option<(DVec3, f64)>,
    /// The cores by cell (`cell` wide, `nx × ny` from `lo`).
    cell: f64,
    nx: usize,
    ny: usize,
    buckets: Vec<Vec<u32>>,
}

/// A round conductor's outline: the sides of its polygon (`Lattice::with_disc`).
pub const DISC_SIDES: usize = 96;

/// Potential and field at `x` of the disc of `radius` about `center` carrying a unit charge
/// per area, in the convention of `zinv::rect_integrals` (`φ = −∫∫ ln r² dA`): outside,
/// that of its charge `πR²` on its axis, `φ = −πR² ln r²`, `E = 2πR² r̂/r`; inside, by Gauss,
/// `E = 2π r` (r from the centre) and `φ = −πR² ln R² + π(R² − r²)` (test
/// `disc_background_is_the_disc`).
#[must_use]
pub fn disc_integrals(center: DVec3, radius: f64, x: DVec3) -> (f64, DVec3) {
    let pi = std::f64::consts::PI;
    let d = DVec3::new(x.x - center.x, x.y - center.y, 0.0);
    let (r2, big) = (d.x * d.x + d.y * d.y, radius * radius);
    if r2 >= big {
        (-pi * big * libm::log(r2), d * (2.0 * pi * big / r2))
    } else {
        (-pi * big * libm::log(big) + pi * (big - r2), d * (2.0 * pi))
    }
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

/// The velocity of a carrier of a gas at thermal speed `sigma` crossing a wall along the
/// unit normal `n` (the flux-weighted half-Maxwellian: the normal speed
/// Rayleigh-distributed, `∝ v e^{−v²/2σ²}`, the tangential speed normal), from a hash of
/// `key`: what a contact at the lattice's temperature gives off.
#[must_use]
pub fn flux_velocity(n: DVec3, sigma: f64, key: u64) -> DVec3 {
    let h = mix(key ^ 0x00f1_0a5e);
    let normal = sigma * libm::sqrt(-2.0 * libm::log(unit(h)));
    let (g, _) = normals(unit(mix(h)), unit(mix(mix(h))));
    let t = DVec3::new(-n.y, n.x, 0.0);
    n * normal + t * (sigma * g)
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
            walls: Vec::new(),
            holes: Vec::new(),
            disc: None,
        }
    }

    /// The cores by cell again, after some went.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn rebucket(&mut self) {
        self.buckets = vec![Vec::new(); self.nx * self.ny];
        for (k, c) in self.cores.iter().enumerate() {
            let i = (((c.x - self.lo.x) / self.cell) as usize).min(self.nx - 1);
            let j = (((c.y - self.lo.y) / self.cell) as usize).min(self.ny - 1);
            self.buckets[i + self.nx * j].push(u32::try_from(k).unwrap_or(u32::MAX));
        }
    }

    /// The lattice with a rectangular hole `lo`..`hi` inside it, which makes it a ring: the
    /// hole's edges become walls (its carriers reflect off them as off a cut), the cores
    /// centred in it go (those across its edge stay: they scatter on the metal's side), and
    /// neither carriers (`fill`) nor background (`background_at`, `carrier_charge`) are put
    /// in it.
    #[must_use]
    pub fn with_hole(mut self, lo: DVec3, hi: DVec3) -> Self {
        let corners = [
            lo,
            DVec3::new(hi.x, lo.y, 0.0),
            hi,
            DVec3::new(lo.x, hi.y, 0.0),
        ];
        for k in 0..4 {
            self.walls.push((corners[k], corners[(k + 1) % 4]));
        }
        self.holes.push((lo, hi));
        let inside = |c: &DVec3| c.x > lo.x && c.x < hi.x && c.y > lo.y && c.y < hi.y;
        self.cores.retain(|c| !inside(c));
        self.rebucket();
        self
    }

    /// The lattice cut to the disc of `radius` about `center` (inside its rectangle): a round
    /// conductor. Its outline is a regular polygon of `DISC_SIDES` walls with the disc's area
    /// (within 0.04 % of the radius of the circle at 96 sides; its carriers reflect off them as
    /// off a cut), the cores centred outside it go, no carriers are put outside it (`fill`),
    /// and its background is the disc's, in closed form (`background_at`).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn with_disc(mut self, center: DVec3, radius: f64) -> Self {
        let n = DISC_SIDES as f64;
        let step = 2.0 * std::f64::consts::PI / n;
        // The vertices' radius that gives the polygon the disc's area: N/2 R_v² sin(2π/N) = πR².
        let rv = radius * libm::sqrt(step / libm::sin(step));
        let vertex = |k: usize| {
            let a = step * k as f64;
            center + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * rv
        };
        for k in 0..DISC_SIDES {
            self.walls.push((vertex(k), vertex((k + 1) % DISC_SIDES)));
        }
        self.disc = Some((center, radius));
        let cores = std::mem::take(&mut self.cores);
        self.cores = cores.into_iter().filter(|c| self.in_outline(*c)).collect();
        self.rebucket();
        self
    }

    /// Whether `x` lies inside its disc's outline (true without a disc): on the inner side of
    /// every side of the polygon.
    #[allow(clippy::cast_precision_loss)]
    fn in_outline(&self, x: DVec3) -> bool {
        let Some((center, radius)) = self.disc else {
            return true;
        };
        let n = DISC_SIDES as f64;
        let step = 2.0 * std::f64::consts::PI / n;
        let apothem = radius * libm::sqrt(step / libm::sin(step)) * libm::cos(0.5 * step);
        let d = x - center;
        (0..DISC_SIDES).all(|k| {
            let a = step * (k as f64 + 0.5);
            d.x * libm::cos(a) + d.y * libm::sin(a) <= apothem
        })
    }

    /// Whether `x` lies outside the conductor though in its rectangle: in a hole, or out of
    /// its disc's outline.
    #[must_use]
    pub fn outside(&self, x: DVec3) -> bool {
        self.in_hole(x) || !self.in_outline(x)
    }

    /// Whether `x` lies in one of its holes.
    #[must_use]
    pub fn in_hole(&self, x: DVec3) -> bool {
        self.holes
            .iter()
            .any(|(lo, hi)| x.x > lo.x && x.x < hi.x && x.y > lo.y && x.y < hi.y)
    }

    /// Its conductor's area: the rectangle's (or its disc's) less its holes'.
    #[must_use]
    pub fn area(&self) -> f64 {
        let whole = match self.disc {
            Some((_, r)) => std::f64::consts::PI * r * r,
            None => (self.hi.x - self.lo.x) * (self.hi.y - self.lo.y),
        };
        whole
            - self
                .holes
                .iter()
                .map(|(lo, hi)| (hi.x - lo.x) * (hi.y - lo.y))
                .sum::<f64>()
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

    /// The cells of the display's accounts (`physics::tube::CellAccounts`): `nx × ny` of
    /// about a unit each, from `lo`.
    #[must_use]
    pub fn account_size(&self) -> (usize, usize) {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = |w: f64| libm::ceil(w - 1e-9).max(1.0) as usize;
        (n(self.hi.x - self.lo.x), n(self.hi.y - self.lo.y))
    }

    /// Books the charge `q` carried along the chord from `a` to `b` (a step) into the
    /// account cells it crosses (`flux`, row by row): each the part of `q (b − a)` along the
    /// piece of the chord inside it. (Booked whole to the cell where the step began, the
    /// steps across a cell's edge biased the cells along a wall: each edge between two rows
    /// gave its thermal crossings' overshoot to the row below and took it from the row
    /// above, which cancels inside but left the rows along the walls a spurious flux into
    /// them, ±0.14 per area in the levels' metal, ~70 % of a strip's current density.)
    #[allow(clippy::cast_precision_loss)]
    pub fn book_flux(&self, flux: &mut [[f64; 2]], a: DVec3, b: DVec3, q: f64) {
        let d = b - a;
        let (ca, cb) = (self.account_cell(a), self.account_cell(b));
        if ca == cb {
            flux[ca][0] += q * d.x;
            flux[ca][1] += q * d.y;
            return;
        }
        // Where the chord crosses the cells' inner edges (fractions of it), in order.
        let (nx, ny) = self.account_size();
        let mut cuts = [0.0_f64; 10];
        let mut n = 0;
        let mut lines = |u0: f64, u1: f64, lo: f64, w: f64, count: usize| {
            for k in 1..count {
                let line = lo + w * k as f64;
                if (line - u0) * (line - u1) < 0.0 && n < cuts.len() {
                    cuts[n] = (line - u0) / (u1 - u0);
                    n += 1;
                }
            }
        };
        lines(a.x, b.x, self.lo.x, (self.hi.x - self.lo.x) / nx as f64, nx);
        lines(a.y, b.y, self.lo.y, (self.hi.y - self.lo.y) / ny as f64, ny);
        let cuts = &mut cuts[..n];
        cuts.sort_by(f64::total_cmp);
        let mut from = 0.0;
        for to in cuts.iter().copied().chain(std::iter::once(1.0)) {
            let c = self.account_cell(a + d * (0.5 * (from + to)));
            flux[c][0] += q * d.x * (to - from);
            flux[c][1] += q * d.y * (to - from);
            from = to;
        }
    }

    /// The account cell (row by row from `lo`) that `x` lies in; the nearest one outside.
    #[must_use]
    pub fn account_cell(&self, x: DVec3) -> usize {
        let (nx, ny) = self.account_size();
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let at = |u: f64, lo: f64, hi: f64, n: usize| {
            (libm::floor((u - lo) / (hi - lo) * n as f64).max(0.0) as usize).min(n - 1)
        };
        at(x.x, self.lo.x, self.hi.x, nx) + nx * at(x.y, self.lo.y, self.hi.y, ny)
    }

    /// The background's potential and field at `x`.
    #[must_use]
    pub fn background_at(&self, x: DVec3) -> (f64, DVec3) {
        let (mut p, mut e) = match self.disc {
            Some((center, radius)) => disc_integrals(center, radius, x),
            None => rect_integrals(self.lo, self.hi, x),
        };
        for (lo, hi) in &self.holes {
            let (ph, eh) = rect_integrals(*lo, *hi, x);
            p -= ph;
            e -= eh;
        }
        (self.background * p, e * self.background)
    }

    /// Each carrier's charge at the start: the background's over them, opposite.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn carrier_charge(&self) -> f64 {
        // A plain rectangle's product in its original order, `(−ρ w) h`: `−ρ (w h)` rounds
        // differently at a density that is not a power of two (Eddy currents in silicon's 2.7
        // changed its run, and its golden hash, by an ulp of every carrier's charge).
        let whole = if self.holes.is_empty() && self.disc.is_none() {
            -self.background * (self.hi.x - self.lo.x) * (self.hi.y - self.lo.y)
        } else {
            -self.background * self.area()
        };
        whole / self.carriers.max(1) as f64
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
            // Out of a hole (or off its disc) it fell in: drawn again uniformly over the
            // rectangle until in the conductor (none without either: bit for bit as before).
            let mut redraw = 0_u64;
            while self.outside(x) && redraw < 1000 {
                redraw += 1;
                let hr = mix(hk ^ mix(0x0401e ^ redraw));
                x = DVec3::new(self.lo.x + unit(hr) * w, self.lo.y + unit(mix(hr)) * h, 0.0);
            }
            // Out of a core it fell in: to its rim, outward (inward instead where outward is
            // out of the conductor: a core across a hole's edge or the outline).
            if let Some(c) = self
                .cores
                .iter()
                .find(|c| (x - **c).length() < self.core_radius)
            {
                let dir = (x - *c).normalize_or(DVec3::X);
                x = *c + dir * (self.core_radius * 1.000_001);
                if self.outside(x) {
                    x = *c - dir * (self.core_radius * 1.000_001);
                }
            }
            let (g1, g2) = normals(unit(mix(mix(hk))), unit(mix(mix(mix(hk)))));
            out.push((x, DVec3::new(g1, g2, 0.0) * self.thermal_speed));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The disc's closed form against the disc summed from 2000 thin rectangles across it
    /// (each exact, `rect_integrals`; the staircase's error O(h²) ~ 1e-6 away from the rim):
    /// relative 1e-5, inside and outside.
    #[test]
    fn disc_background_is_the_disc() {
        let (center, radius) = (DVec3::new(1.0, -2.0, 0.0), 3.0);
        let n = 2000;
        let h = 2.0 * radius / f64::from(n);
        let stairs = |x: DVec3| {
            let (mut p, mut e) = (0.0, DVec3::ZERO);
            for i in 0..n {
                let u = -radius + h * (f64::from(i) + 0.5);
                let half = libm::sqrt(radius * radius - u * u);
                let lo = center + DVec3::new(u - 0.5 * h, -half, 0.0);
                let hi = center + DVec3::new(u + 0.5 * h, half, 0.0);
                let (pi, ei) = rect_integrals(lo, hi, x);
                p += pi;
                e += ei;
            }
            (p, e)
        };
        let mut worst: f64 = 0.0;
        for (r, a) in [(0.0, 0.0), (0.9, 0.3), (2.1, 2.0), (4.5, -1.0), (9.0, 4.0)] {
            let x = center + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * r;
            let (p, e) = disc_integrals(center, radius, x);
            let (ps, es) = stairs(x);
            let scale = std::f64::consts::PI * radius * radius;
            let dp = (p - ps).abs() / (scale * (1.0 + libm::log(r * r + 1.0).abs()));
            let de = (e - es).length() / (2.0 * std::f64::consts::PI * radius);
            println!(
                "r {r}: φ {p:.6} vs {ps:.6}, E {:.6} vs {:.6}",
                e.length(),
                es.length()
            );
            worst = worst.max(dp).max(de);
        }
        println!("worst relative {worst:.2e}");
        assert!(worst < 1e-5, "worst {worst}");
    }

    /// A round lattice: its carriers and cores inside its outline, its area the disc's, its
    /// carriers' charge the background's, opposite.
    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn disc_lattice_stays_inside() {
        let (center, radius) = (DVec3::new(2.0, 3.0, 0.0), 4.0);
        let half = DVec3::new(radius * 1.001, radius * 1.001, 0.0);
        let (lo, hi) = (center - half, center + half);
        let cores = Lattice::random_cores(lo, hi, 400, 0.01, 5);
        let lattice =
            Lattice::new(lo, hi, 1.0, cores, 0.01, -1.0, 1.0, 3000, 7).with_disc(center, radius);
        let pi = std::f64::consts::PI;
        assert!((lattice.area() - pi * radius * radius).abs() < 1e-12);
        let total = lattice.carrier_charge() * lattice.carriers as f64;
        assert!(
            (total + pi * radius * radius).abs() < 1e-9,
            "charge {total}"
        );
        // The outline's vertices lie 0.036 % beyond the radius at 96 sides.
        let rv = radius * 1.000_357;
        let fill = lattice.fill();
        assert_eq!(fill.len(), 3000);
        for (x, _) in &fill {
            assert!(!lattice.outside(*x), "carrier at {x}");
            assert!((*x - center).length() <= rv, "carrier at {x}");
        }
        assert!(lattice.cores.iter().all(|c| (*c - center).length() <= rv));
        // Uniform: as many in the inner half of the area as in the outer.
        let inner = fill
            .iter()
            .filter(|(x, _)| (*x - center).length() < radius / std::f64::consts::SQRT_2)
            .count();
        assert!((inner as f64 / 3000.0 - 0.5).abs() < 0.03, "inner {inner}");
        assert_eq!(lattice.walls.len(), DISC_SIDES);
    }
}
