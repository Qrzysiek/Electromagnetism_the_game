//! What a run's results mean for the level's goals: each current's verdict and
//! the run's.

use super::*;

impl CurrentVerdict {
    /// The verdict on a goal's range `min..max` from the two runs' values.
    pub fn new(preview: f64, fine: f64, min: f64, max: f64) -> Self {
        let error = ERROR_FACTOR * (fine - preview).abs();
        let (lo, hi) = (fine - error, fine + error);
        let status = if lo >= min && hi <= max {
            TubeStatus::Met
        } else if hi < min || lo > max {
            TubeStatus::Missed
        } else {
            TubeStatus::Uncertain
        };
        Self {
            preview,
            fine,
            error,
            status,
        }
    }
}

impl TubeVerdict {
    /// From the two runs' results.
    pub fn new(spec: &TubeSpec, preview: &TubeResults, fine: &TubeResults) -> Self {
        let currents: Vec<CurrentVerdict> = spec
            .goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    preview.currents.get(k).copied().flatten().unwrap_or(0.0),
                    fine.currents.get(k).copied().flatten().unwrap_or(0.0),
                    g.min,
                    g.max,
                )
            })
            .collect();
        let potential: Vec<CurrentVerdict> = spec
            .potential_goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    g.compared(preview.potential.get(k).copied().flatten().unwrap_or(0.0)),
                    g.compared(fine.potential.get(k).copied().flatten().unwrap_or(0.0)),
                    g.min,
                    g.max,
                )
            })
            .collect();
        let heat: Vec<CurrentVerdict> = spec
            .heat_goals
            .iter()
            .enumerate()
            .map(|(k, g)| {
                CurrentVerdict::new(
                    preview.heat.get(k).copied().flatten().unwrap_or(0.0),
                    fine.heat.get(k).copied().flatten().unwrap_or(0.0),
                    g.min,
                    g.max,
                )
            })
            .collect();
        let arrivals: Vec<_> = spec
            .projectiles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.detector.is_some())
            .map(|(j, _)| {
                let (a, b) = (
                    preview.arrivals.get(j).copied().flatten(),
                    fine.arrivals.get(j).copied().flatten(),
                );
                let status = match (a.is_some(), b.is_some()) {
                    (true, true) => TubeStatus::Met,
                    (false, false) => TubeStatus::Missed,
                    _ => TubeStatus::Uncertain,
                };
                (j, a, b, status)
            })
            .collect();
        let all = currents
            .iter()
            .chain(potential.iter())
            .chain(heat.iter())
            .map(|c| c.status)
            .chain(arrivals.iter().map(|a| a.3));
        let capped = preview.capped || fine.capped;
        let mut status = if capped {
            TubeStatus::Uncertain
        } else {
            TubeStatus::Met
        };
        for s in all {
            status = match (status, s) {
                (TubeStatus::Missed, _) | (_, TubeStatus::Missed) => TubeStatus::Missed,
                (TubeStatus::Uncertain, _) | (_, TubeStatus::Uncertain) => TubeStatus::Uncertain,
                _ => TubeStatus::Met,
            };
        }
        Self {
            currents,
            potential,
            heat,
            arrivals,
            capped,
            status,
        }
    }
}
