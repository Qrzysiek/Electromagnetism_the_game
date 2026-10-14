//! A meter's time average of a tube level's maps, the player's option (the owner: the
//! fields are volatile, a real instrument integrates). A first-order low-pass, as an RC
//! meter or a power meter's detector: `x̄ ← x̄ + (x − x̄)(1 − e^{−Δt/τ})` over the run's
//! time. The potential and the electric field are linear in their sources, so the maps
//! average the sources (the particles' charge binned on the map's cells, the electrodes'
//! surface densities and potentials); the electromagnetic maps average B_z and the
//! Poynting vector S per sample (the mean S is what a power meter reads).

use std::collections::BTreeMap;

use level::tube::TubeFrame;

use crate::em::EmSnapshot;

/// The cells the particles' charge is averaged on (as `potential::map_frame` bins them).
const BIN: f64 = 0.25;

#[derive(Default)]
pub struct MapAverage {
    revision: Option<u64>,
    t: f64,
    /// Per cell and kind: charge, and charge times x and y (its centre).
    cells: BTreeMap<(i64, i64, u8), [f64; 3]>,
    sigma: Vec<f64>,
    potentials: Vec<f64>,
    em_revision: Option<u64>,
    em_t: f64,
    em: Vec<[f64; 3]>,
}

impl MapAverage {
    /// The frame averaged with time constant `tau` up to `frame` (run `revision`).
    #[allow(clippy::cast_possible_truncation)]
    pub fn frame(&mut self, revision: u64, frame: &TubeFrame, tau: f64) -> TubeFrame {
        let mut now: BTreeMap<(i64, i64, u8), [f64; 3]> = BTreeMap::new();
        for ((x, q), kind) in frame.x.iter().zip(&frame.q).zip(&frame.kind) {
            let (px, py, q) = (f64::from(x[0]), f64::from(x[1]), f64::from(*q));
            let key = ((px / BIN).floor() as i64, (py / BIN).floor() as i64, *kind);
            let c = now.entry(key).or_insert([0.0; 3]);
            c[0] += q;
            c[1] += q * px;
            c[2] += q * py;
        }
        let sigma: Vec<f64> = frame.sigma.iter().map(|s| f64::from(*s)).collect();
        let potentials: Vec<f64> = frame.potentials.iter().map(|v| f64::from(*v)).collect();
        let restart =
            self.revision != Some(revision) || frame.t < self.t || self.sigma.len() != sigma.len();
        if restart {
            self.revision = Some(revision);
            self.cells = now;
            self.sigma = sigma;
            self.potentials = potentials;
            self.t = frame.t;
        } else if frame.t > self.t {
            let w = 1.0 - libm_exp(-(frame.t - self.t) / tau.max(1e-9));
            for c in self.cells.values_mut() {
                for v in c.iter_mut() {
                    *v *= 1.0 - w;
                }
            }
            for (key, c) in now {
                let a = self.cells.entry(key).or_insert([0.0; 3]);
                for (v, n) in a.iter_mut().zip(c) {
                    *v += w * n;
                }
            }
            // Cells whose charge has faded away.
            let total: f64 = self.cells.values().map(|c| c[0].abs()).sum();
            self.cells.retain(|_, c| c[0].abs() > 1e-9 * total);
            for (a, n) in self.sigma.iter_mut().zip(&sigma) {
                *a += w * (n - *a);
            }
            for (a, n) in self.potentials.iter_mut().zip(&potentials) {
                *a += w * (n - *a);
            }
            self.t = frame.t;
        }
        let mut out = TubeFrame {
            t: frame.t,
            sigma: self.sigma.iter().map(|s| *s as f32).collect(),
            potentials: self.potentials.iter().map(|v| *v as f32).collect(),
            currents: frame.currents.clone(),
            ..TubeFrame::default()
        };
        for (k, ((_, _, kind), c)) in self.cells.iter().enumerate() {
            if c[0] == 0.0 {
                continue;
            }
            out.x.push([(c[1] / c[0]) as f32, (c[2] / c[0]) as f32]);
            out.q.push(c[0] as f32);
            out.kind.push(*kind);
            out.id.push(u32::try_from(k).unwrap_or(u32::MAX));
        }
        out
    }

    /// The electromagnetic field's snapshot averaged with time constant `tau`.
    #[allow(clippy::cast_possible_truncation)]
    pub fn snapshot(&mut self, snap: &EmSnapshot, tau: f64) -> EmSnapshot {
        let restart = self.em_revision != Some(snap.revision)
            || snap.t < self.em_t
            || self.em.len() != snap.samples.len();
        if restart {
            self.em_revision = Some(snap.revision);
            self.em = snap
                .samples
                .iter()
                .map(|s| [f64::from(s[0]), f64::from(s[1]), f64::from(s[2])])
                .collect();
            self.em_t = snap.t;
        } else if snap.t > self.em_t {
            let w = 1.0 - libm_exp(-(snap.t - self.em_t) / tau.max(1e-9));
            for (a, s) in self.em.iter_mut().zip(&snap.samples) {
                for (v, n) in a.iter_mut().zip(s) {
                    *v += w * (f64::from(*n) - *v);
                }
            }
            self.em_t = snap.t;
        }
        EmSnapshot {
            samples: self
                .em
                .iter()
                .map(|a| [a[0] as f32, a[1] as f32, a[2] as f32])
                .collect(),
            ..snap.clone()
        }
    }
}

/// `e^x` for the display (the game's own arithmetic; not the authoritative path).
fn libm_exp(x: f64) -> f64 {
    x.exp()
}
