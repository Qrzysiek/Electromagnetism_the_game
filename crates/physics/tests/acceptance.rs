//! Tests D1–D3 for detector acceptance (direction and energy windows, PHYSICS.md §6.1).
//! Run with `cargo test -p physics --test acceptance -- --nocapture`.

mod common;

use common::{UNIT_PARTICLE, cube};
use physics::DVec3;
use physics::field::UniformElectric;
use physics::geometry::{Aabb, Region};
use physics::trajectory::{Acceptance, Outcome, RunSettings, Scenario, run};
use physics::verify::{Boundary, Status, Tolerances, verify};

fn scenario(e: f64, acceptance: Acceptance) -> Scenario<UniformElectric> {
    Scenario {
        field: UniformElectric {
            e: DVec3::new(e, 0.0, 0.0),
        },
        obstacles: vec![],
        particle: UNIT_PARTICLE,
        c: f64::INFINITY,
        x0: DVec3::ZERO,
        p0: DVec3::new(0.6, 0.8, 0.0),
        detector: Some(Region::Box(Aabb {
            min: DVec3::new(10.0, -100.0, -1.0),
            max: DVec3::new(12.0, 100.0, 1.0),
        })),
        bounds: Some(cube(200.0)),
        t_max: 100.0,
        radiation_reaction: false,
        acceptance: Some(acceptance),
        gates: Vec::new(),
    }
}

/// Field-free straight flight at angle atan2(0.8, 0.6) = 53.13° from +x: accepted by a
/// window around +x exactly when its half-angle exceeds that, with margin equal to the
/// difference.
#[test]
fn d1_direction_window() {
    let angle = 0.8f64.atan2(0.6);
    for (half, arrives) in [(angle + 0.01, true), (angle - 0.01, false)] {
        let acc = Acceptance {
            direction: Some((DVec3::X, half)),
            kinetic: None,
        };
        let tr = run(&scenario(0.0, acc), &RunSettings::with_tolerance(1e-12));
        let m = tr.margins.unwrap().acceptance.unwrap();
        println!("D1 half-angle {half:.4}: {:?}, margin {m:.3e}", tr.outcome);
        assert_eq!(tr.outcome == Outcome::Arrived, arrives);
        assert_eq!(tr.outcome == Outcome::Rejected, !arrives);
        assert!((m - (half - angle)).abs() < 1e-12);
    }
}

/// Uniform field along x, detector face at x = 10: the kinetic energy on entry is
/// T0 + qE·10 exactly (energy conservation). The relative energy margin of a window
/// [lo, hi] is min(T − lo, hi − T)/hi.
#[test]
fn d2_energy_window() {
    let (e, t0) = (0.03, 0.5);
    let t = t0 + e * 10.0;
    for (lo, hi, arrives) in [(0.7, 0.9, true), (0.9, 1.1, false), (0.5, 0.79, false)] {
        let acc = Acceptance {
            direction: None,
            kinetic: Some((lo, hi)),
        };
        let tr = run(&scenario(e, acc), &RunSettings::with_tolerance(1e-12));
        let m = tr.margins.unwrap().acceptance.unwrap();
        let exact = (t - lo).min(hi - t) / hi;
        println!(
            "D2 window [{lo}, {hi}]: {:?}, margin {m:.6} (exact {exact:.6})",
            tr.outcome
        );
        assert_eq!(tr.outcome == Outcome::Arrived, arrives);
        assert!((m - exact).abs() < 1e-10);
    }
}

/// A flight exactly at the edge of the window is not declared verified: the margin (or
/// the outcome) is not reliable at the level of the numerical error.
#[test]
fn d3_edge_of_the_window_is_not_verified() {
    let angle = 0.8f64.atan2(0.6);
    let acc = Acceptance {
        direction: Some((DVec3::X, angle + 1e-14)),
        kinetic: None,
    };
    let v = verify(&scenario(0.0, acc), Tolerances::default());
    println!("D3: {:?}", v.status);
    assert!(
        matches!(
            v.status,
            Status::SmallMargin {
                boundary: Boundary::Acceptance,
                ..
            } | Status::OutcomeMismatch { .. }
        ),
        "{:?}",
        v.status
    );
}
