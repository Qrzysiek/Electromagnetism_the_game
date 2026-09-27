//! Gates: regions a flight must pass, in order, before the detector counts (PHYSICS.md
//! §6.2). Straight field-free flights have exact references.
//!
//! `cargo test --release -p physics --test gates -- --nocapture`

use physics::DVec3;
use physics::dynamics::Particle;
use physics::field::Coulomb;
use physics::geometry::{Aabb, Region};
use physics::trajectory::{Acceptance, Gate, Outcome, RunSettings, Scenario, run};
use physics::verify::{Boundary, Status, verify};

fn cell(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::Box(Aabb {
        min: DVec3::new(x0, y0, -1.0),
        max: DVec3::new(x1, y1, 1.0),
    })
}

fn gate(x0: f64, y0: f64, x1: f64, y1: f64) -> Gate {
    Gate {
        region: cell(x0, y0, x1, y1),
        acceptance: None,
    }
}

/// A particle flying along +x at height `y` with speed 1 (Newtonian), detector at
/// x ∈ [20, 22].
fn flight(y: f64, angle: f64, gates: Vec<Gate>) -> Scenario<Coulomb> {
    Scenario {
        field: Coulomb::new(&[]),
        obstacles: vec![],
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c: f64::INFINITY,
        x0: DVec3::new(0.0, y, 0.0),
        p0: DVec3::new(angle.cos(), angle.sin(), 0.0),
        detector: Some(cell(20.0, -5.0, 22.0, 5.0)),
        bounds: None,
        t_max: 100.0,
        radiation_reaction: false,
        acceptance: None,
        gates,
    }
}

const TOL: f64 = 1e-12;

/// Two gates passed in order: arrival. The gate margins are the depths reached inside
/// (exact: the distance from the line y = 0.3 to the nearest face of each box).
#[test]
fn gates_passed_in_order() {
    let scn = flight(
        0.3,
        0.0,
        vec![gate(5.0, -1.0, 6.0, 2.0), gate(10.0, -2.0, 12.0, 1.0)],
    );
    let t = run(&scn, &RunSettings::with_tolerance(TOL));
    assert_eq!(t.outcome, Outcome::Arrived);
    let m = t.margins.unwrap();
    // Gate 1 (x 5..6, y −1..2): deepest at x = 5.5, depth min(0.5, 1.3, 1.7) = 0.5.
    // Gate 2 (x 10..12, y −2..1): deepest point depth min(1, 0.7, 2.3) = 0.7.
    println!("gate margins {:?}", m.gates);
    assert!((m.gates[0] + 0.5).abs() < 1e-12 && (m.gates[1] + 0.7).abs() < 1e-12);
}

/// Missing a gate, or passing it out of order, does not count.
#[test]
fn skipped_and_out_of_order_gates() {
    // Gate 2 lies off the path.
    let t = run(
        &flight(
            0.3,
            0.0,
            vec![gate(5.0, -1.0, 6.0, 2.0), gate(10.0, 3.0, 12.0, 4.0)],
        ),
        &RunSettings::with_tolerance(TOL),
    );
    assert_eq!(t.outcome, Outcome::SkippedGate(1));
    // Closest approach to the missed gate: 3 − 0.3.
    assert!((t.margins.unwrap().gates[1] - 2.7).abs() < 1e-12);
    // The second gate lies first on the path: passing it then does not count.
    let t = run(
        &flight(
            0.3,
            0.0,
            vec![gate(10.0, -1.0, 11.0, 2.0), gate(5.0, -1.0, 6.0, 2.0)],
        ),
        &RunSettings::with_tolerance(TOL),
    );
    assert_eq!(t.outcome, Outcome::SkippedGate(1));
}

/// A gate with a direction window: entering at 10° against a ±5° window along +x is not
/// a pass; at 2° it is. The acceptance margin is exact (5° − the angle).
#[test]
fn gate_direction_window() {
    let window = Gate {
        region: cell(8.0, -10.0, 9.0, 10.0),
        acceptance: Some(Acceptance {
            direction: Some((DVec3::X, 5f64.to_radians())),
            kinetic: None,
        }),
    };
    for (deg, expected) in [(10.0, Outcome::SkippedGate(0)), (2.0, Outcome::Arrived)] {
        let t = run(
            &flight(0.0, f64::to_radians(deg), vec![window]),
            &RunSettings::with_tolerance(TOL),
        );
        let m = t.margins.unwrap().gate_acceptance[0].unwrap();
        println!(
            "entering at {deg}°: {:?}, acceptance margin {m:.6}",
            t.outcome
        );
        assert_eq!(t.outcome, expected);
        assert!((m - (5f64 - deg).to_radians()).abs() < 1e-12);
    }
}

/// Verification: a flight grazing a gate's edge by 1e-11 is too close to call; one with
/// a clear pass is verified.
#[test]
fn gate_margins_enter_verification() {
    let tol = physics::verify::Tolerances::default();
    let clear = verify(&flight(0.3, 0.0, vec![gate(5.0, -1.0, 6.0, 2.0)]), tol);
    assert_eq!(clear.status, Status::Verified);
    let graze = verify(
        &flight(-1.0 + 1e-11, 0.0, vec![gate(5.0, -1.0, 6.0, 2.0)]),
        tol,
    );
    println!("grazing: {:?}", graze.status);
    assert!(matches!(
        graze.status,
        Status::SmallMargin {
            boundary: Boundary::Gate(0),
            ..
        } | Status::OutcomeMismatch { .. }
    ));
}
