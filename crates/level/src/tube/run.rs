//! A recorded run: the frames between its steps, the particles' trails, what the
//! contacts took and gave off, and the fields and forces it shows.

use super::*;

/// The point a fraction `s` of the way from `a` to `b`.
#[allow(clippy::cast_possible_truncation)]
pub(super) fn towards(a: [f32; 2], b: [f32; 2], s: f64) -> [f32; 2] {
    let s = s as f32;
    [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s]
}

/// The distance from `a` to `b`.
pub(super) fn distance(a: [f32; 2], b: [f32; 2]) -> f64 {
    libm::hypot(f64::from(b[0] - a[0]), f64::from(b[1] - a[1]))
}

/// The mean over `[a, b]` of a series `(t, v)` (times increasing) taken as linear between
/// its samples and constant beyond its ends; at `b` if the interval is empty, None if the
/// series is.
pub(super) fn series_mean(series: &[(f64, f64)], a: f64, b: f64) -> Option<f64> {
    if series.is_empty() {
        return None;
    }
    let at = |t: f64| {
        let k = series.partition_point(|s| s.0 <= t);
        match (k.checked_sub(1).map(|i| series[i]), series.get(k)) {
            (Some(p), Some(q)) => p.1 + (q.1 - p.1) * (t - p.0) / (q.0 - p.0),
            (Some(p), None) => p.1,
            (None, Some(q)) => q.1,
            (None, None) => 0.0,
        }
    };
    if b <= a {
        return Some(at(b));
    }
    let k0 = series.partition_point(|s| s.0 <= a);
    let k1 = series.partition_point(|s| s.0 < b);
    let mut points = vec![(a, at(a))];
    points.extend_from_slice(&series[k0..k1.max(k0)]);
    points.push((b, at(b)));
    let integral: f64 = points
        .windows(2)
        .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
        .sum();
    Some(integral / (b - a))
}

impl TubeRun {
    /// B along z at `t`.
    #[must_use]
    pub fn b_at(&self, t: f64) -> f64 {
        self.b_wave.map_or(self.b_z, |w| w.b(t))
    }

    /// Whether the run has a magnetic field at all (steady or changing).
    #[must_use]
    pub fn has_field(&self) -> bool {
        self.b_z != 0.0 || self.b_wave.is_some()
    }

    /// Appends a streamed chunk (its frames and samples follow this run's).
    pub fn extend(&mut self, chunk: TubeRun) {
        if self.self_b_grid.is_none() {
            self.self_b_grid = chunk.self_b_grid;
        }
        self.frames.extend(chunk.frames);
        if self.collected.len() < chunk.collected.len() {
            self.collected.resize(chunk.collected.len(), Vec::new());
        }
        for (a, c) in self.collected.iter_mut().zip(chunk.collected) {
            a.extend(c);
        }
        if self.middle.len() < chunk.middle.len() {
            self.middle
                .resize(chunk.middle.len(), [Vec::new(), Vec::new()]);
        }
        for (a, [cx, cy]) in self.middle.iter_mut().zip(chunk.middle) {
            a[0].extend(cx);
            a[1].extend(cy);
        }
        if self.currents.len() < chunk.currents.len() {
            self.currents.resize(chunk.currents.len(), None);
        }
        for (a, c) in self.currents.iter_mut().zip(chunk.currents) {
            *a = a.or(c);
        }
        if self.potential.len() < chunk.potential.len() {
            self.potential.resize(chunk.potential.len(), None);
        }
        for (a, c) in self.potential.iter_mut().zip(chunk.potential) {
            *a = a.or(c);
        }
        if self.heat.len() < chunk.heat.len() {
            self.heat.resize(chunk.heat.len(), None);
        }
        for (a, c) in self.heat.iter_mut().zip(chunk.heat) {
            *a = a.or(c);
        }
        if self.work.len() < chunk.work.len() {
            self.work.resize(chunk.work.len(), Vec::new());
        }
        for (a, c) in self.work.iter_mut().zip(chunk.work) {
            a.extend(c);
        }
        self.potential_steps.extend(chunk.potential_steps);
        if self.tracks.len() < chunk.tracks.len() {
            self.tracks.resize(chunk.tracks.len(), Vec::new());
        }
        for (t, c) in self.tracks.iter_mut().zip(chunk.tracks) {
            t.extend(c);
        }
        if self.arrivals.len() < chunk.arrivals.len() {
            self.arrivals.resize(chunk.arrivals.len(), None);
        }
        for (a, c) in self.arrivals.iter_mut().zip(chunk.arrivals) {
            *a = a.or(c);
        }
        self.in_flight = chunk.in_flight;
        self.steps = chunk.steps;
        if self.segments.is_empty() {
            self.segments = chunk.segments;
            self.potentials = chunk.potentials;
            self.potential_range = chunk.potential_range;
        }
    }

    /// The picture at time `t` between the stored frames, for the display (smooth at any
    /// playback speed): each particle in both neighbouring frames moved linearly between
    /// them, the surface densities and potentials interpolated. A metal's carrier that a
    /// contact takes before the next frame goes into it and one it gives off after the
    /// frame comes out of it (`taken`, `given_off`): frames a metal's run keeps steps apart,
    /// where they stood still and then vanished, and the new ones appeared all at once (the
    /// owner, 2026-10-26: the carriers by the contacts "locked", moving "in bursts"). Any
    /// other particle absorbed before the next frame stays where it was, one emitted after
    /// it appears with that frame (their frames are steps). Before the first frame the
    /// first, after the last the last.
    pub fn frame_between(&self, t: f64) -> Option<TubeFrame> {
        let k = self.frames.partition_point(|f| f.t <= t).saturating_sub(1);
        let a = self.frames.get(k)?;
        let Some(b) = self.frames.get(k + 1).filter(|b| b.t > a.t && t > a.t) else {
            return Some(a.clone());
        };
        #[allow(clippy::cast_possible_truncation)]
        let w = ((t - a.t) / (b.t - a.t)).clamp(0.0, 1.0) as f32;
        let next: std::collections::HashMap<u32, usize> =
            b.id.iter().enumerate().map(|(i, id)| (*id, i)).collect();
        let lerp = |u: f32, v: f32| u + (v - u) * w;
        let n = a.x.len();
        let (mut x, mut q, mut kind, mut id) = (
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        );
        for (i, p) in a.x.iter().enumerate() {
            let now = match a.id.get(i).map(|me| (*me, next.get(me))) {
                // A sampled carrier along its steps, the others straight on.
                Some((me, Some(&j))) => Some(
                    b.path(me)
                        .and_then(|path| self.along(*p, a.t, path, b.t, t))
                        .unwrap_or([lerp(p[0], b.x[j][0]), lerp(p[1], b.x[j][1])]),
                ),
                Some((me, None)) => self.taken(*p, me, a.t, b, t),
                None => Some(*p),
            };
            if let Some(now) = now {
                x.push(now);
                q.push(a.q.get(i).copied().unwrap_or_default());
                kind.push(a.kind.get(i).copied().unwrap_or_default());
                id.push(a.id.get(i).copied().unwrap_or_default());
            }
        }
        let before: std::collections::HashSet<u32> = a.id.iter().copied().collect();
        for (j, me) in b.id.iter().enumerate() {
            if before.contains(me) {
                continue;
            }
            if let Some(now) = self.given_off(b.x[j], *me, a.t, b, t) {
                x.push(now);
                q.push(b.q[j]);
                kind.push(b.kind[j]);
                id.push(*me);
            }
        }
        Some(TubeFrame {
            self_b: a
                .self_b
                .iter()
                .zip(&b.self_b)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            t,
            x,
            q,
            kind,
            id,
            sigma: a
                .sigma
                .iter()
                .zip(&b.sigma)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            potentials: a
                .potentials
                .iter()
                .zip(&b.potentials)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            currents: a
                .currents
                .iter()
                .zip(&b.currents)
                .map(|(u, v)| lerp(*u, *v))
                .collect(),
            heat: a.heat + (b.heat - a.heat) * f64::from(w),
            paths: Vec::new(),
            cells: lerp_cells(&a.cells, &b.cells, w),
        })
    }

    /// The position at `t` of a carrier at `x0` at `t0` that then took the steps `path`
    /// (the last at `t1`, each a step `dt` after the one before): between its two steps
    /// around `t`. None if the path does not span the time from `t0` to `t1`.
    pub(super) fn along(
        &self,
        x0: [f32; 2],
        t0: f64,
        path: &[[f32; 2]],
        t1: f64,
        t: f64,
    ) -> Option<[f32; 2]> {
        let n = path.len();
        if n == 0 || self.dt <= 0.0 {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        let start = t1 - self.dt * n as f64;
        if (start - t0).abs() > 0.25 * self.dt {
            return None;
        }
        // Steps after `start`: 0 at x0, k at path[k - 1].
        #[allow(clippy::cast_precision_loss)]
        let u = ((t - start) / self.dt).clamp(0.0, n as f64);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (libm::floor(u) as usize).min(n - 1);
        let p0 = if k == 0 { x0 } else { path[k - 1] };
        let p1 = path[k];
        #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
        let w = (u - k as f64) as f32;
        Some([p0[0] + (p1[0] - p0[0]) * w, p0[1] + (p1[1] - p0[1]) * w])
    }

    /// Where a metal's carrier at `x0` at `t0` that a contact takes before frame `b` is at
    /// `t`: along its recorded steps if it is sampled (they end at its last step), else
    /// straight to the nearest point of its metal's contact faces at the carriers' mean
    /// speed, there by `b` at the latest; None once taken. Outside a metal with contacts,
    /// where it was.
    pub(super) fn taken(
        &self,
        x0: [f32; 2],
        id: u32,
        t0: f64,
        b: &TubeFrame,
        t: f64,
    ) -> Option<[f32; 2]> {
        if let Some(path) = b.path(id).filter(|p| !p.is_empty()) {
            #[allow(clippy::cast_precision_loss)]
            let last = t0 + self.dt * path.len() as f64;
            return if t <= last {
                self.along(x0, t0, path, last, t)
            } else {
                None
            };
        }
        let Some((face, speed)) = self.nearest_contact(x0) else {
            return Some(x0);
        };
        let span = (distance(x0, face) / speed).min(b.t - t0);
        let s = if span > 0.0 {
            (t - t0) / span
        } else {
            f64::INFINITY
        };
        (s <= 1.0).then(|| towards(x0, face, s))
    }

    /// Where a metal's carrier at `x1` in frame `b`, given off by a contact since `t0`, is
    /// at `t`: from its first recorded step if it is sampled, else from the nearest point of
    /// its metal's contact faces at the carriers' mean speed, leaving no earlier than `t0`;
    /// None before. Outside a metal with contacts, None (it appears with `b`).
    pub(super) fn given_off(
        &self,
        x1: [f32; 2],
        id: u32,
        t0: f64,
        b: &TubeFrame,
        t: f64,
    ) -> Option<[f32; 2]> {
        if let Some(path) = b.path(id).filter(|p| !p.is_empty()) {
            #[allow(clippy::cast_precision_loss)]
            let first = b.t - self.dt * (path.len() - 1) as f64;
            if t < first || self.dt <= 0.0 {
                return None;
            }
            let u = (t - first) / self.dt;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let k = (libm::floor(u) as usize).min(path.len() - 1);
            let next = path.get(k + 1).copied().unwrap_or(path[k]);
            #[allow(clippy::cast_precision_loss)]
            return Some(towards(path[k], next, u - k as f64));
        }
        let (face, speed) = self.nearest_contact(x1)?;
        let span = (distance(face, x1) / speed).min(b.t - t0);
        let born = b.t - span;
        if t < born {
            return None;
        }
        let s = if span > 0.0 { (t - born) / span } else { 1.0 };
        Some(towards(face, x1, s.min(1.0)))
    }

    /// The nearest point to `x` of the contact faces of the metal it is in, and that metal's
    /// carriers' mean speed (`√(π/2)` times the thermal speed `√(kT/m)`, a 2D Maxwellian's).
    pub(super) fn nearest_contact(&self, x: [f32; 2]) -> Option<([f32; 2], f64)> {
        let p = [f64::from(x[0]), f64::from(x[1])];
        let m = self
            .metals
            .iter()
            .find(|m| p[0] >= m.lo.x && p[0] <= m.hi.x && p[1] >= m.lo.y && p[1] <= m.hi.y)?;
        let d2 = |u: &[f64; 2]| (u[0] - p[0]) * (u[0] - p[0]) + (u[1] - p[1]) * (u[1] - p[1]);
        let nearest = m
            .contacts
            .iter()
            .map(|[a, b]| {
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                let len2 = dx * dx + dy * dy;
                let s = if len2 > 0.0 {
                    (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                [a[0] + s * dx, a[1] + s * dy]
            })
            .min_by(|u, v| d2(u).total_cmp(&d2(v)))?;
        let speed = m.thermal_speed * libm::sqrt(std::f64::consts::FRAC_PI_2);
        #[allow(clippy::cast_possible_truncation)]
        (speed > 0.0).then_some(([nearest[0] as f32, nearest[1] as f32], speed))
    }

    /// Carrier `id`'s path from `t0` to `t1`: its recorded steps where it is sampled
    /// (`TubeFrame::paths`), else its frames' positions, ending where it is at `t1`
    /// (`frame_between`); empty if it is in no frame of that time.
    #[must_use]
    pub fn trail(&self, id: u32, t0: f64, t1: f64) -> Vec<[f32; 2]> {
        let k0 = self.frames.partition_point(|f| f.t < t0);
        let k1 = self.frames.partition_point(|f| f.t <= t1);
        let mut out = Vec::new();
        // The frames up to t1, and the next one's steps before t1.
        for f in &self.frames[k0..(k1 + 1).min(self.frames.len())] {
            match f.path(id) {
                Some(path) =>
                {
                    #[allow(clippy::cast_precision_loss)]
                    for (i, x) in path.iter().enumerate() {
                        let ti = f.t - self.dt * (path.len() - 1 - i) as f64;
                        if ti >= t0 && ti <= t1 {
                            out.push(*x);
                        }
                    }
                }
                None if f.t <= t1 => {
                    if let Some(i) = f.id.iter().position(|j| *j == id) {
                        out.push(f.x[i]);
                    }
                }
                None => {}
            }
        }
        if let Some(now) = self.frame_between(t1)
            && let Some(i) = now.id.iter().position(|j| *j == id)
        {
            out.push(now.x[i]);
        }
        out
    }

    /// The frame at time `t` (the last one at or before it; the first before it).
    pub fn frame_at(&self, t: f64) -> Option<&TubeFrame> {
        let k = self.frames.partition_point(|f| f.t <= t);
        self.frames.get(k.saturating_sub(1))
    }

    /// The field of a frame at `x` as `field_at` gives it, faster: far segments by their
    /// multipoles (`physics::zinv::segment_field_fast`, within 5e-5): for the display's
    /// arrows and colour scale, many points every frame.
    pub fn e_at(&self, frame: &TubeFrame, x: DVec3) -> DVec3 {
        let mut e = DVec3::ZERO;
        for ((a, b), s) in self.segments.iter().zip(&frame.sigma) {
            e += physics::zinv::segment_field_fast(*a, *b, x) * f64::from(*s);
        }
        let h2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        for (p, q) in frame.x.iter().zip(&frame.q) {
            let r = DVec3::new(x.x - f64::from(p[0]), x.y - f64::from(p[1]), 0.0);
            e += r * (2.0 * f64::from(*q) / (r.x * r.x + r.y * r.y).max(h2));
        }
        for m in &self.metals {
            e += physics::zinv::rect_integrals(m.lo, m.hi, x).1 * m.background;
        }
        // A changing field's induced part (no potential: it circulates).
        if let Some(w) = &self.b_wave {
            e += w.induced(frame.t, x);
        }
        e
    }

    /// The mean field in metal `metal`'s middle on `frame` (`MetalView::strip`): the mean of
    /// `e_at` over the central three fifths along it (7 points) on lines half a unit apart
    /// from a unit inside each side (its sides hold the surface charge of its potential
    /// against its surroundings, within a Debye length; the Hall field is the bulk's), or on
    /// its centre line alone in a strip narrower than 3. Its carriers' graininess scatters
    /// the part across from frame to frame: on Hall probe's reference (100 macroparticles per
    /// area) by 0.11 against its 0.47, on the centre line alone by 0.19; on Positive
    /// carriers' lower strip (50 per area, a weaker current) by 0.20, on the centre line by
    /// 0.35, against 0.27.
    #[must_use]
    pub fn middle_field(&self, frame: &TubeFrame, metal: usize) -> DVec3 {
        let Some(m) = self.metals.get(metal) else {
            return DVec3::ZERO;
        };
        let (centre, axis, normal, length, width) = m.strip();
        #[allow(clippy::cast_possible_truncation)]
        let lines = libm::floor((0.5 * width - 1.0) / 0.5 + 1e-9).max(0.0) as i32;
        let (mut sum, mut n) = (DVec3::ZERO, 0.0);
        for l in -lines..=lines {
            for i in 0..7 {
                let x = centre
                    + axis * (0.3 * length * (f64::from(i) / 3.0 - 1.0))
                    + normal * (0.5 * f64::from(l));
                sum += self.e_at(frame, x);
                n += 1.0;
            }
        }
        sum / n
    }

    /// The Hall effect in metal `metal` of a run in a field along z: its carriers' mean
    /// drift from frame `before` to frame `now` in the strip's middle (the central three
    /// fifths along it, a quarter away from its sides), the field's magnetic push on that
    /// drift, `q ⟨v⟩ × B`, and the electric push across the strip, `q E⊥`: the field across
    /// its middle (`hall_field`) averaged over the same time (the run's series `hall`,
    /// linear between its frames: over a time unit its scatter falls to a few per cent),
    /// per unit charge of the carriers' sign. None without a field or with too few
    /// carriers.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn hall_forces(
        &self,
        now: &TubeFrame,
        before: &TubeFrame,
        metal: usize,
    ) -> Option<HallForces> {
        let m = self.metals.get(metal)?;
        let span = now.t - before.t;
        if !self.has_field() || span <= 0.0 {
            return None;
        }
        let (centre, axis, normal, length, width) = m.strip();
        let then: std::collections::HashMap<u32, [f32; 2]> = before
            .id
            .iter()
            .zip(&before.x)
            .map(|(id, x)| (*id, *x))
            .collect();
        let (mut sum, mut n, mut sign) = (DVec3::ZERO, 0.0, 0.0);
        for ((x, id), kind) in now.x.iter().zip(&now.id).zip(&now.kind) {
            let p = DVec3::new(f64::from(x[0]), f64::from(x[1]), 0.0);
            let d = p - centre;
            if d.dot(axis).abs() >= 0.3 * length || d.dot(normal).abs() >= 0.5 * width - 0.25 {
                continue;
            }
            let Some(x0) = then.get(id) else {
                continue;
            };
            sum += p - DVec3::new(f64::from(x0[0]), f64::from(x0[1]), 0.0);
            n += 1.0;
            sign += match *kind {
                KIND_NEGATIVE => -1.0,
                KIND_POSITIVE => 1.0,
                _ => 0.0,
            };
        }
        if n < 10.0 {
            return None;
        }
        let q = f64::signum(sign);
        let drift = sum / (n * span);
        // q v × B ẑ = q B (v_y, −v_x).
        // A changing field: B at the window's middle.
        let b = self.b_at(0.5 * (before.t + now.t));
        let magnetic = DVec3::new(drift.y, -drift.x, 0.0) * (q * b);
        let mean = |k: usize| {
            self.middle
                .get(metal)
                .and_then(|s| series_mean(&s[k], before.t, now.t))
        };
        let field = match (mean(0), mean(1)) {
            (Some(x), Some(y)) => DVec3::new(x, y, 0.0),
            _ => self.middle_field(now, metal),
        };
        let electric = field * q;
        Some(HallForces {
            centre,
            axis,
            normal,
            width,
            drift,
            magnetic,
            electric,
            collisions: -(electric + magnetic),
        })
    }

    /// The potential and field of a frame at `x`, as the maps show them: the segments'
    /// charge exactly, the particles as uniformly charged rods of radius `PREVIEW_SEGMENT`
    /// (the preview's resolution of its space charge; the physics' rods are thinner,
    /// `SOFTENING`), exactly Coulomb's beyond it. (Softened as `r⃗/(r² + ε²)`, which never
    /// becomes Coulomb's, a metal's carriers failed to cancel its background's exact field
    /// by `2λε²/r³` each: across Positive carriers' lower strip, from the upper's 128
    /// charges, a field that took the Hall arrows' electric push to 1.73 of the magnetic.)
    pub fn field_at(&self, frame: &TubeFrame, x: DVec3) -> (f64, DVec3) {
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        for ((a, b), s) in self.segments.iter().zip(&frame.sigma) {
            let (f, g) = physics::zinv::segment_integrals(*a, *b, x);
            phi += f64::from(*s) * f;
            e += g * f64::from(*s);
        }
        let h2 = PREVIEW_SEGMENT * PREVIEW_SEGMENT;
        let log_h2 = libm::log(h2);
        for (p, q) in frame.x.iter().zip(&frame.q) {
            let r = DVec3::new(x.x - f64::from(p[0]), x.y - f64::from(p[1]), 0.0);
            let r2 = r.x * r.x + r.y * r.y;
            let lambda = f64::from(*q);
            if r2 >= h2 {
                phi -= lambda * libm::log(r2);
                e += r * (2.0 * lambda / r2);
            } else {
                phi -= lambda * (log_h2 - 1.0 + r2 / h2);
                e += r * (2.0 * lambda / h2);
            }
        }
        for m in &self.metals {
            let (f, g) = physics::zinv::rect_integrals(m.lo, m.hi, x);
            phi += m.background * f;
            e += g * m.background;
        }
        if let Some(w) = &self.b_wave {
            e += w.induced(frame.t, x);
        }
        (phi, e)
    }
}
