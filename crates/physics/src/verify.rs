//! Outcome verification: "no numerical luck" (SPEC §2.3, PHYSICS.md §7).
//!
//! A flight is computed at two tolerances (preview and verify, 100× apart). The outcome is
//! **verified** when both runs agree and every margin — the closest approach to each
//! boundary that was not crossed, the penetration depth into the one that was, and the
//! time left before `t_max` — exceeds `SAFETY` times its numerical error estimate plus a
//! floor. The error estimate of a margin is its difference between the two runs, which
//! bounds the error of the looser run and therefore overestimates that of the tighter one
//! by roughly the tolerance ratio.

use crate::field::FieldSolver;
use crate::trajectory::{MARGIN_SAFE, Outcome, RunSettings, Scenario, Trajectory, run};

/// Required ratio of margin to its error estimate.
pub const SAFETY: f64 = 10.0;
/// Absolute margin floor in grid units (about the rounding floor measured in M1/M2).
pub const FLOOR: f64 = 1e-9;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerances {
    pub preview: f64,
    pub verify: f64,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self {
            preview: 1e-10,
            verify: 1e-12,
        }
    }
}

/// Which boundary a margin refers to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    Obstacle(usize),
    Bounds,
    Detector,
    /// Time remaining before `t_max` at the event.
    TimeLimit,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Status {
    Verified,
    /// The two runs disagree on the outcome.
    OutcomeMismatch {
        preview: Outcome,
        verify: Outcome,
    },
    /// A margin is too small relative to its numerical error.
    SmallMargin {
        boundary: Boundary,
        margin: f64,
        error_estimate: f64,
    },
    /// The integration failed in at least one run.
    Failed,
}

impl Status {
    pub fn is_verified(&self) -> bool {
        matches!(self, Status::Verified)
    }
}

#[derive(Clone, Debug)]
pub struct Verification {
    pub status: Status,
    pub preview: Trajectory,
    pub verified: Trajectory,
}

impl Verification {
    /// The outcome of the tighter run.
    pub fn outcome(&self) -> Outcome {
        self.verified.outcome
    }
}

/// Runs a scenario at both tolerances and classifies the result.
pub fn verify<F: FieldSolver>(scn: &Scenario<F>, tol: Tolerances) -> Verification {
    verify_pair(scn, scn, tol)
}

/// Verification where the tighter run also uses a finer field model (`fine`, e.g.
/// conductors at verification resolution, PHYSICS.md §2.6), so that the model's own
/// error enters the comparison like the integration error.
pub fn verify_pair<F: FieldSolver>(
    scn: &Scenario<F>,
    fine: &Scenario<F>,
    tol: Tolerances,
) -> Verification {
    let preview = run(scn, &RunSettings::with_tolerance(tol.preview));
    let verified = run(fine, &RunSettings::with_tolerance(tol.verify));
    let status = classify(&preview, &verified, scn.t_max);
    Verification {
        status,
        preview,
        verified,
    }
}

/// Compares two runs of the same scenario (`b` the tighter one).
pub fn classify(a: &Trajectory, b: &Trajectory, t_max: f64) -> Status {
    if matches!(a.outcome, Outcome::Failed(_)) || matches!(b.outcome, Outcome::Failed(_)) {
        return Status::Failed;
    }
    if a.outcome != b.outcome {
        return Status::OutcomeMismatch {
            preview: a.outcome,
            verify: b.outcome,
        };
    }
    let (ma, mb) = match (&a.margins, &b.margins) {
        (Some(ma), Some(mb)) => (ma, mb),
        _ => panic!("verification requires margins"),
    };

    let mut boundaries: Vec<(Boundary, f64, f64)> = Vec::new();
    for (i, (&x, &y)) in ma.obstacles.iter().zip(&mb.obstacles).enumerate() {
        boundaries.push((Boundary::Obstacle(i), x, y));
    }
    if let (Some(x), Some(y)) = (ma.bounds, mb.bounds) {
        boundaries.push((Boundary::Bounds, x, y));
    }
    if let (Some(x), Some(y)) = (ma.detector, mb.detector) {
        boundaries.push((Boundary::Detector, x, y));
    }

    let mut worst: Option<Status> = None;
    let mut worst_ratio = f64::INFINITY;
    let mut consider = |boundary: Boundary, margin: f64, error_estimate: f64, floor: f64| {
        let ratio = margin.abs() / (SAFETY * error_estimate + floor);
        if ratio <= 1.0 && ratio < worst_ratio {
            worst_ratio = ratio;
            worst = Some(Status::SmallMargin {
                boundary,
                margin,
                error_estimate,
            });
        }
    };

    for (boundary, x, y) in boundaries {
        if x.abs() >= MARGIN_SAFE && y.abs() >= MARGIN_SAFE {
            continue;
        }
        consider(boundary, y, (x - y).abs(), FLOOR);
    }
    if b.outcome != Outcome::Timeout {
        consider(
            Boundary::TimeLimit,
            t_max - b.end.t,
            (a.end.t - b.end.t).abs(),
            FLOOR * t_max,
        );
    }
    worst.unwrap_or(Status::Verified)
}
