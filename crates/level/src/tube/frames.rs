//! What a run records for the display: its frames, the metals' cells and rates,
//! their interpolation and smoothing, the maps.

use super::*;

impl TubeFrame {
    /// Carrier `id`'s positions at every step since the previous frame (`paths`).
    #[must_use]
    pub fn path(&self, id: u32) -> Option<&[[f32; 2]]> {
        self.paths
            .binary_search_by_key(&id, |(i, _)| *i)
            .ok()
            .map(|k| self.paths[k].1.as_slice())
    }
}

/// Each metal's carriers' mean displacement (x, y) from frame `a` to frame `b`: the mean over
/// the carriers in it at `b` that were in flight at `a` (thousands: the thermal motion
/// averages to ~10 % of a frame's drift, and does not accumulate). Summed over frames, the
/// drift of the carriers as a whole (the game's chevrons; test `metal_drift_carries_the_current`).
#[allow(clippy::cast_precision_loss)]
#[must_use]
pub fn metal_steps(metals: &[MetalView], a: &TubeFrame, b: &TubeFrame) -> Vec<[f64; 2]> {
    let before: std::collections::HashMap<u32, [f32; 2]> =
        a.id.iter().zip(&a.x).map(|(id, x)| (*id, *x)).collect();
    let mut sums = vec![([0.0_f64; 2], 0_usize); metals.len()];
    for (x, id) in b.x.iter().zip(&b.id) {
        let inside = |m: &MetalView| {
            let (px, py) = (f64::from(x[0]), f64::from(x[1]));
            px >= m.lo.x && px <= m.hi.x && py >= m.lo.y && py <= m.hi.y
        };
        let (Some(m), Some(x0)) = (metals.iter().position(inside), before.get(id)) else {
            continue;
        };
        sums[m].0[0] += f64::from(x[0] - x0[0]);
        sums[m].0[1] += f64::from(x[1] - x0[1]);
        sums[m].1 += 1;
    }
    sums.into_iter()
        .map(|(sum, n)| {
            let n = n.max(1) as f64;
            [sum[0] / n, sum[1] / n]
        })
        .collect()
}

impl MetalView {
    /// The metal as a strip: its centre, long axis (x for a square), the normal across it,
    /// its length and width.
    #[must_use]
    pub fn strip(&self) -> (DVec3, DVec3, DVec3, f64, f64) {
        let (dx, dy) = (self.hi.x - self.lo.x, self.hi.y - self.lo.y);
        let centre = (self.lo + self.hi) * 0.5;
        if dx >= dy {
            (centre, DVec3::X, DVec3::Y, dx, dy)
        } else {
            (centre, DVec3::Y, DVec3::X, dy, dx)
        }
    }
}

/// The metals' accounts a fraction `w` of the way from `a` to `b` (`TubeFrame::cells`).
pub(super) fn lerp_cells(a: &[MetalCells], b: &[MetalCells], w: f32) -> Vec<MetalCells> {
    let lerp = |u: f32, v: f32| u + (v - u) * w;
    a.iter()
        .zip(b)
        .map(|(u, v)| MetalCells {
            flux: u
                .flux
                .iter()
                .zip(&v.flux)
                .map(|(p, q)| [lerp(p[0], q[0]), lerp(p[1], q[1])])
                .collect(),
            work: u
                .work
                .iter()
                .zip(&v.work)
                .map(|(p, q)| lerp(*p, *q))
                .collect(),
        })
        .collect()
}

impl MetalView {
    /// Cell `c`'s centre (`cells`, row by row).
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn cell_centre(&self, c: usize) -> [f64; 2] {
        let (nx, ny) = (self.cells.0.max(1), self.cells.1.max(1));
        let (w, h) = (
            (self.hi.x - self.lo.x) / nx as f64,
            (self.hi.y - self.lo.y) / ny as f64,
        );
        [
            self.lo.x + ((c % nx) as f64 + 0.5) * w,
            self.lo.y + ((c / nx) as f64 + 0.5) * h,
        ]
    }

    /// The cell `x` lies in, if it is inside the metal.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    #[must_use]
    pub fn cell_at(&self, x: [f64; 2]) -> Option<usize> {
        let (nx, ny) = self.cells;
        let (u, v) = (
            (x[0] - self.lo.x) / (self.hi.x - self.lo.x),
            (x[1] - self.lo.y) / (self.hi.y - self.lo.y),
        );
        if nx == 0 || ny == 0 || !(0.0..=1.0).contains(&u) || !(0.0..=1.0).contains(&v) {
            return None;
        }
        let i = ((u * nx as f64) as usize).min(nx - 1);
        let j = ((v * ny as f64) as usize).min(ny - 1);
        Some(i + nx * j)
    }

    /// Whether a cut separates `a` from `b` (the segment between them meets one).
    #[must_use]
    pub fn separated(&self, a: [f64; 2], b: [f64; 2]) -> bool {
        let side = |o: [f64; 2], p: [f64; 2], q: [f64; 2]| {
            (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
        };
        self.walls.iter().any(|[w0, w1]| {
            let (d1, d2) = (side(*w0, *w1, a), side(*w0, *w1, b));
            let (d3, d4) = (side(a, b, *w0), side(a, b, *w1));
            d1 * d2 <= 0.0 && d3 * d4 <= 0.0 && (d1 != 0.0 || d2 != 0.0)
        })
    }
}

/// Each cell's mean with those of its eight neighbours that no cut separates from it
/// (`TubeRun::metal_rates`).
#[must_use]
pub(super) fn smoothed(m: &MetalView, rates: &MetalRates) -> MetalRates {
    let (nx, ny) = m.cells;
    let n = rates.j.len().min(rates.power.len());
    if n < nx * ny || nx * ny == 0 {
        return rates.clone();
    }
    let mut out = MetalRates {
        j: vec![[0.0; 2]; n],
        power: vec![0.0; n],
    };
    for c in 0..n {
        let (i, k) = (c % nx, c / nx);
        let centre = m.cell_centre(c);
        let mut count = 0.0;
        for kk in k.saturating_sub(1)..=(k + 1).min(ny - 1) {
            for ii in i.saturating_sub(1)..=(i + 1).min(nx - 1) {
                let d = ii + nx * kk;
                if d != c && m.separated(centre, m.cell_centre(d)) {
                    continue;
                }
                out.j[c][0] += rates.j[d][0];
                out.j[c][1] += rates.j[d][1];
                out.power[c] += rates.power[d];
                count += 1.0;
            }
        }
        out.j[c][0] /= count;
        out.j[c][1] /= count;
        out.power[c] /= count;
    }
    out
}

impl TubeRun {
    /// Each metal's accounts at `t` (`TubeFrame::cells`, between the frames around it).
    pub(super) fn cells_at(&self, t: f64) -> Option<Vec<MetalCells>> {
        let k = self.frames.partition_point(|f| f.t <= t);
        let a = &self.frames[k.checked_sub(1)?];
        let Some(b) = self.frames.get(k) else {
            return Some(a.cells.clone());
        };
        #[allow(clippy::cast_possible_truncation)]
        let w = if b.t > a.t {
            ((t - a.t) / (b.t - a.t)) as f32
        } else {
            0.0
        };
        Some(lerp_cells(&a.cells, &b.cells, w))
    }

    /// The swinging solenoid's period, if its field swings (a sine, a repeating pulse).
    #[must_use]
    pub fn drive_period(&self) -> Option<f64> {
        match self.b_wave?.waveform {
            physics::circuit::Waveform::Sine {
                omega, amplitude, ..
            } if omega > 0.0 && amplitude != 0.0 => Some(std::f64::consts::TAU / omega),
            physics::circuit::Waveform::Pulse { period, .. } if period > 0.0 => Some(period),
            _ => None,
        }
    }

    /// The metal maps at `t` (`MetalMaps`). A steady drive: the current density and the
    /// power per area over the `span` before `t` (`metal_rates`, smoothed). A swinging field
    /// (`drive_period`): over the whole half periods before `t` (up to four: each a whole
    /// cycle of the heat, which goes as `Ḃ²`), the power as before and the current by
    /// lock-in (`lock_in_rates`: in a plain mean it would cancel). Either way the span is
    /// cut into sub-windows (`METAL_MAP_BATCHES`, or the half periods) whose spread gives
    /// each value's standard error. None before two sub-windows have passed.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    #[must_use]
    pub fn metal_maps(&self, t: f64, span: f64) -> Option<MetalMaps> {
        let first = self.frames.first()?.t;
        let (n, sub, lock_in) = match self.drive_period() {
            Some(p) if t - first >= p => {
                let half = 0.5 * p;
                (((t - first) / half).floor().min(4.0) as usize, half, true)
            }
            _ => (METAL_MAP_BATCHES, span / METAL_MAP_BATCHES as f64, false),
        };
        let from = t - sub * n as f64;
        let mut batches = Vec::with_capacity(n);
        for i in 0..n {
            let a = from + sub * i as f64;
            let rates = if lock_in {
                self.lock_in_rates(a, a + sub)?
            } else {
                self.metal_rates(a, a + sub, false)?
            };
            batches.push(
                self.metals
                    .iter()
                    .zip(&rates)
                    .map(|(m, r)| smoothed(m, r))
                    .collect::<Vec<_>>(),
            );
        }
        if batches.len() < 2 {
            return None;
        }
        let k = batches.len() as f64;
        let mut values = Vec::new();
        let mut errors = Vec::new();
        for m in 0..self.metals.len() {
            let cells = batches[0].get(m).map_or(0, |r| r.power.len());
            let mut value = MetalRates {
                j: vec![[0.0; 2]; cells],
                power: vec![0.0; cells],
            };
            let mut error = value.clone();
            for c in 0..cells {
                let get = |f: &dyn Fn(&MetalRates) -> f64| -> (f64, f64) {
                    let v: Vec<f64> = batches.iter().map(|b| f(&b[m])).collect();
                    let mean = v.iter().sum::<f64>() / k;
                    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (k - 1.0);
                    (mean, libm::sqrt(var / k))
                };
                let (jx, ex) = get(&|r| r.j[c][0]);
                let (jy, ey) = get(&|r| r.j[c][1]);
                let (p, ep) = get(&|r| r.power[c]);
                value.j[c] = [jx, jy];
                value.power[c] = p;
                error.j[c] = [ex, ey];
                error.power[c] = ep;
            }
            values.push(value);
            errors.push(error);
        }
        Some(MetalMaps {
            from,
            to: t,
            lock_in,
            values,
            errors,
        })
    }

    /// Under a swinging field, each metal's current density by lock-in over `[t0, t1]`,
    /// unsmoothed: `J = ⟨j Ḃ⟩/⟨Ḃ²⟩` per cell (from the frames' accounts between frames and
    /// the field's change over each), the current per unit rate of change of the field,
    /// times the rms of `Ḃ`: the rms of the part of the current in step with `Ḃ`, the way it
    /// flows while the field rises (Ohmic eddy currents are that part alone). The power as
    /// `metal_rates`'.
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn lock_in_rates(&self, t0: f64, t1: f64) -> Option<Vec<MetalRates>> {
        let mut rates = self.metal_rates(t0, t1, false)?;
        // The times: the span's ends and the frames between.
        let mut times = vec![t0];
        times.extend(
            self.frames
                .iter()
                .map(|f| f.t)
                .filter(|t| *t > t0 && *t < t1),
        );
        times.push(t1);
        let cells: Vec<Vec<MetalCells>> = times
            .iter()
            .map(|t| self.cells_at(*t))
            .collect::<Option<_>>()?;
        let (mut den, mut num) = (0.0, vec![Vec::<[f64; 2]>::new(); self.metals.len()]);
        for (w, (a, b)) in times
            .windows(2)
            .zip(cells.windows(2).map(|c| (&c[0], &c[1])))
        {
            let dt = w[1] - w[0];
            if dt <= 0.0 {
                continue;
            }
            let db = self.b_at(w[1]) - self.b_at(w[0]);
            den += db * db / dt;
            for (m, (u, v)) in a.iter().zip(b).enumerate() {
                num[m].resize(u.flux.len(), [0.0; 2]);
                for (c, (p, q)) in u.flux.iter().zip(&v.flux).enumerate() {
                    num[m][c][0] += (f64::from(q[0]) - f64::from(p[0])) * db / dt;
                    num[m][c][1] += (f64::from(q[1]) - f64::from(p[1])) * db / dt;
                }
            }
        }
        if den <= 0.0 {
            return None;
        }
        let rms = libm::sqrt(den / (t1 - t0));
        for ((r, m), num) in rates.iter_mut().zip(&self.metals).zip(&num) {
            let (nx, ny) = m.cells;
            let area = (m.hi.x - m.lo.x) * (m.hi.y - m.lo.y) / (nx * ny).max(1) as f64;
            for (j, n) in r.j.iter_mut().zip(num) {
                *j = [n[0] / (area * den) * rms, n[1] / (area * den) * rms];
            }
        }
        Some(rates)
    }

    /// Each metal's mean current density and power per area by cell (`MetalView::cells`)
    /// over `[t0, t1]` (within the frames), from the frames' accounts (`TubeFrame::cells`):
    /// the power the heat goal measures (the field's work on the carriers there; an
    /// isolated metal's, the induced field's). With `smooth`, each cell's mean with those of
    /// its eight neighbours that no cut separates from it (`METAL_MAP_SPAN`). None without
    /// two frames in the span.
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn metal_rates(&self, t0: f64, t1: f64, smooth: bool) -> Option<Vec<MetalRates>> {
        let (first, last) = (self.frames.first()?.t, self.frames.last()?.t);
        let (t0, t1) = (t0.max(first), t1.min(last));
        if t1 <= t0 {
            return None;
        }
        let (a, b) = (self.cells_at(t0)?, self.cells_at(t1)?);
        Some(
            self.metals
                .iter()
                .zip(a.iter().zip(&b))
                .map(|(m, (a, b))| {
                    let (nx, ny) = m.cells;
                    let area = (m.hi.x - m.lo.x) * (m.hi.y - m.lo.y) / (nx * ny).max(1) as f64;
                    let per = 1.0 / ((t1 - t0) * area);
                    let j: Vec<[f64; 2]> = a
                        .flux
                        .iter()
                        .zip(&b.flux)
                        .map(|(p, q)| {
                            [
                                (f64::from(q[0]) - f64::from(p[0])) * per,
                                (f64::from(q[1]) - f64::from(p[1])) * per,
                            ]
                        })
                        .collect();
                    let power: Vec<f64> = a
                        .work
                        .iter()
                        .zip(&b.work)
                        .map(|(p, q)| (f64::from(*q) - f64::from(*p)) * per)
                        .collect();
                    let rates = MetalRates { j, power };
                    if smooth { smoothed(m, &rates) } else { rates }
                })
                .collect(),
        )
    }
}
