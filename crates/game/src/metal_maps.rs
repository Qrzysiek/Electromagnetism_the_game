//! A metal's maps on the board: their grid, span and colour scale.

use super::*;

/// The metal maps' values as drawn: each cell faded by how far it stands above its noise
/// (`level::tube::MetalMaps`' standard errors): within one standard error of nothing it
/// vanishes, beyond three it is drawn in full. What the noise cannot tell from nothing is
/// not drawn as something (a laminated core's heat, whose fall is the point).
pub(super) fn faded(maps: &level::tube::MetalMaps) -> Vec<level::tube::MetalRates> {
    let weight = |v: f64, e: f64| {
        if e > 0.0 {
            ((v / e - 1.0) / 2.0).clamp(0.0, 1.0)
        } else {
            1.0
        }
    };
    maps.values
        .iter()
        .zip(&maps.errors)
        .map(|(v, e)| level::tube::MetalRates {
            j: v.j
                .iter()
                .zip(&e.j)
                .map(|(j, ej)| {
                    let w = weight(
                        j[0].hypot(j[1]),
                        (0.5 * (ej[0] * ej[0] + ej[1] * ej[1])).sqrt(),
                    );
                    [j[0] * w, j[1] * w]
                })
                .collect(),
            power: v
                .power
                .iter()
                .zip(&e.power)
                .map(|(p, ep)| p * weight(p.abs(), *ep))
                .collect(),
        })
        .collect()
}

/// The span the metal maps average over: twice the maps' time constant τ when one is set
/// (a boxcar of 2τ is about as quiet as a first-order average over τ), else
/// `level::tube::METAL_MAP_SPAN`.
pub fn metal_map_span(tau: f64) -> f64 {
    if tau > 0.0 {
        2.0 * tau
    } else {
        level::tube::METAL_MAP_SPAN
    }
}

/// The metal maps' grid over the board `min`..`max` (setup `revision`, time `t`), a quarter
/// cell apart (the shader interpolates between samples: the cells stay sharp): each sample
/// in a metal its cell's power (`heat`, as the magnetic map's signed B_z) or current density
/// (as the energy flow's vector), zero outside; and the current's arrows, one per cell.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(super) fn metal_map_grid(
    run: &level::tube::TubeRun,
    rates: &[level::tube::MetalRates],
    (min, max): (physics::DVec3, physics::DVec3),
    heat: bool,
    (revision, t): (u64, f64),
) -> (em::EmSnapshot, TubeArrows) {
    let spacing = 0.25;
    let nx = (((max.x - min.x) / spacing).ceil() as usize).max(2);
    let ny = (((max.y - min.y) / spacing).ceil() as usize).max(2);
    let mut samples = vec![[0.0_f32; 3]; nx * ny];
    for (k, s) in samples.iter_mut().enumerate() {
        let p = [
            min.x + ((k % nx) as f64 + 0.5) * spacing,
            min.y + ((k / nx) as f64 + 0.5) * spacing,
        ];
        for (m, r) in run.metals.iter().zip(rates) {
            if let Some(c) = m.cell_at(p).filter(|c| *c < r.power.len()) {
                *s = if heat {
                    [r.power[c] as f32, 0.0, 0.0]
                } else {
                    [0.0, r.j[c][0] as f32, r.j[c][1] as f32]
                };
                break;
            }
        }
    }
    let arrows = if heat {
        Vec::new()
    } else {
        run.metals
            .iter()
            .zip(rates)
            .flat_map(|(m, r)| {
                (0..r.j.len()).map(move |c| {
                    let x = m.cell_centre(c);
                    (
                        Vec2::new(x[0] as f32, x[1] as f32),
                        Vec2::new(r.j[c][0] as f32, r.j[c][1] as f32),
                    )
                })
            })
            .collect()
    };
    let snap = em::EmSnapshot {
        revision,
        t,
        origin: (min.x, min.y),
        spacing,
        nx,
        ny,
        samples,
        c: 0.0,
        radiated: 0.0,
        at_end: false,
        delay: 0.0,
        leads: Vec::new(),
        lead_currents: Vec::new(),
        lead_charges: Vec::new(),
        lead_peak: 0.0,
    };
    (snap, arrows)
}

/// The 95th percentile over every metal's cells of `f(rates, cell)` (a metal map's colour
/// unit: the strongest few cells saturate).
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(super) fn metal_percentile(
    rates: &[level::tube::MetalRates],
    f: &dyn Fn(&level::tube::MetalRates, usize) -> f64,
) -> f64 {
    let mut v: Vec<f64> = rates
        .iter()
        .flat_map(|r| (0..r.power.len()).map(move |c| f(r, c)))
        .collect();
    v.sort_by(f64::total_cmp);
    v.get((0.95 * v.len() as f64) as usize)
        .copied()
        .unwrap_or(0.0)
}
