//! Outcome verification (PHYSICS.md §7): clear results are verified, results decided by
//! less than the numerical error are marginal.

mod common;

use common::{CoulombOrbit, central_charge, cube};
use physics::DVec3;
use physics::trajectory::Outcome;
use physics::verify::{Boundary, Status, Tolerances, verify};

fn grazing(delta: f64) -> (Outcome, Status, f64) {
    let (kappa, c) = (1.0, 3.0);
    let x0 = DVec3::new(-20.0, 1.0, 0.0);
    let p0 = DVec3::new(1.0, 0.0, 0.0);
    let r_min = CoulombOrbit::from_state(kappa, 1.0, c, x0, p0).r_min();
    let mut scn = central_charge(kappa, r_min * (1.0 + delta), c, x0, p0, 1e4);
    scn.bounds = Some(cube(40.0));
    let v = verify(&scn, Tolerances::default());
    let m = v.verified.margins.as_ref().unwrap().obstacles[0];
    (v.outcome(), v.status, m)
}

#[test]
fn clear_hit_and_clear_miss_are_verified() {
    for delta in [1e-3, -1e-3, 1e-6, -1e-6] {
        let (outcome, status, margin) = grazing(delta);
        println!("δ = {delta:e}: {outcome:?}, margin {margin:.3e}, {status:?}");
        assert_eq!(outcome == Outcome::Collided(0), delta > 0.0);
        assert!(status.is_verified(), "δ = {delta:e}: {status:?}");
    }
}

#[test]
fn margins_match_analytic_closest_approach() {
    // Radius r_min (1 − δ): the obstacle margin is the gap r_min δ.
    for delta in [1e-2, 1e-4, 1e-6] {
        let (_, _, margin) = grazing(-delta);
        let x0 = DVec3::new(-20.0, 1.0, 0.0);
        let r_min = CoulombOrbit::from_state(1.0, 1.0, 3.0, x0, DVec3::new(1.0, 0.0, 0.0)).r_min();
        let expected = r_min * delta;
        let rel = (margin - expected).abs() / expected;
        println!(
            "δ = {delta:e}: margin {margin:.6e}, expected {expected:.6e}, relative error {rel:.1e}"
        );
        assert!(rel < 1e-5, "relative error {rel:.2e}");
    }
}

#[test]
fn razor_thin_outcomes_are_marginal() {
    for delta in [1e-11, -1e-11, 1e-13, -1e-13] {
        let (outcome, status, margin) = grazing(delta);
        println!("δ = {delta:e}: {outcome:?}, margin {margin:.3e}, {status:?}");
        let marginal = matches!(
            status,
            Status::SmallMargin {
                boundary: Boundary::Obstacle(0),
                ..
            } | Status::OutcomeMismatch { .. }
        );
        assert!(marginal, "δ = {delta:e}: {status:?}");
    }
}
