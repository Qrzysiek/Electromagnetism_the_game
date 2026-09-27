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

/// Regression: a clear exit through the map edge was flagged marginal because the
/// penetration-depth continuation sometimes stopped after one step (it compared the
/// integrator state with the dense output at the same time, which differ by rounding).
#[test]
fn clear_exit_through_bounds_is_verified() {
    use physics::dynamics::{Kinematics, Particle};
    use physics::field::{Coulomb, FixedCharge};
    use physics::geometry::{Aabb, Region, Shape, Sphere};
    use physics::trajectory::Scenario;

    let charge = FixedCharge {
        position: DVec3::new(15.0, 12.0, 0.0),
        charge: 1e6,
        radius: 0.3,
    };
    let kin = Kinematics::new(1.0, 5.0);
    let scn = Scenario {
        field: Coulomb::new(&[charge]),
        obstacles: vec![Shape::Sphere(Sphere {
            center: charge.position,
            radius: charge.radius,
        })],
        particle: Particle {
            charge: 1e-6,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c: 5.0,
        x0: DVec3::new(0.0, 10.0, 0.0),
        p0: kin.momentum_from_kinetic_energy(0.5, DVec3::X),
        detector: Some(Region::Box(Aabb {
            min: DVec3::new(27.0, 8.0, -1.0),
            max: DVec3::new(30.0, 12.0, 1.0),
        })),
        bounds: Some(Aabb {
            min: DVec3::splat(-1.0),
            max: DVec3::new(31.0, 21.0, 1.0),
        }),
        t_max: 200.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let v = verify(&scn, Tolerances::default());
    println!(
        "{:?} {:?} margins {:?}",
        v.outcome(),
        v.status,
        v.verified.margins
    );
    assert_eq!(v.outcome(), Outcome::LeftBounds);
    assert!(v.status.is_verified(), "{:?}", v.status);
}
