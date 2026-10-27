use super::sim::fourier_at;
use super::*;

/// The picture between two frames: a particle in both moves linearly, one absorbed
/// before the second stays where it was, the surface densities and potentials are
/// interpolated; before the first frame the first, after the last the last.
#[test]
fn frames_interpolate() {
    let frame = |t: f64, x: Vec<[f32; 2]>, id: Vec<u32>, sigma: f32| TubeFrame {
        t,
        q: vec![-1.0; x.len()],
        kind: vec![KIND_NEGATIVE; x.len()],
        x,
        id,
        sigma: vec![sigma; 3],
        potentials: vec![sigma; 2],
        currents: Vec::new(),
        heat: 0.0,
        paths: Vec::new(),
        cells: Vec::new(),
    };
    let run = TubeRun {
        frames: vec![
            frame(1.0, vec![[0.0, 0.0], [5.0, 5.0]], vec![7, 8], 1.0),
            // Particle 8 absorbed, 9 emitted, 7 moved and now listed second.
            frame(2.0, vec![[9.0, 9.0], [1.0, 2.0]], vec![9, 7], 3.0),
        ],
        ..TubeRun::default()
    };
    let mid = run.frame_between(1.25).unwrap();
    assert!((mid.t - 1.25).abs() < 1e-12);
    assert_eq!(mid.x, vec![[0.25, 0.5], [5.0, 5.0]]);
    assert_eq!(mid.id, vec![7, 8]);
    assert_eq!(mid.sigma, vec![1.5; 3]);
    assert_eq!(mid.potentials, vec![1.5; 2]);
    assert_eq!(run.frame_between(0.5).unwrap().x, run.frames[0].x);
    assert_eq!(run.frame_between(3.0).unwrap().x, run.frames[1].x);
}

/// A carrier recorded at every step (`TubeFrame::paths`) moves along its steps between
/// frames, not straight across: turning a corner between frames 0.3 apart (steps of
/// 0.1: from (0, 0) to (1, 0), (1, 1) and (0, 1)), at t = 1.15 it is halfway up the
/// corner's side, (1, 0.5) (straight across it would be at (0, 0.5)); its trail from
/// t = 1.0 to 1.25 goes round the corner; an unrecorded carrier moves straight.
#[test]
fn sampled_carriers_follow_their_steps() {
    let frame = |t: f64, x: Vec<[f32; 2]>, paths: Vec<(u32, Vec<[f32; 2]>)>| TubeFrame {
        t,
        q: vec![-1.0; x.len()],
        kind: vec![KIND_NEGATIVE; x.len()],
        id: vec![7, 8],
        x,
        paths,
        ..TubeFrame::default()
    };
    let run = TubeRun {
        frames: vec![
            frame(1.0, vec![[0.0, 0.0], [0.0, 0.0]], Vec::new()),
            frame(
                1.3,
                vec![[0.0, 1.0], [0.0, 1.0]],
                vec![(7, vec![[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])],
            ),
        ],
        dt: 0.1,
        fine_every: 7,
        ..TubeRun::default()
    };
    let mid = run.frame_between(1.15).unwrap();
    let close = |a: [f32; 2], b: [f32; 2]| (a[0] - b[0]).abs() < 1e-5 && (a[1] - b[1]).abs() < 1e-5;
    assert!(close(mid.x[0], [1.0, 0.5]), "sampled: {:?}", mid.x[0]);
    assert!(close(mid.x[1], [0.0, 0.5]), "straight: {:?}", mid.x[1]);
    let trail = run.trail(7, 1.0, 1.25);
    let expected = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.5, 1.0]];
    assert_eq!(trail.len(), expected.len(), "{trail:?}");
    for (a, b) in trail.iter().zip(expected) {
        assert!(close(*a, b), "{trail:?}");
    }
}

/// A metal's carriers that its contacts take or give off between frames move between
/// them (the owner, 2026-10-26: they stood still, then vanished, and the new ones
/// appeared all at once). A metal 10 × 2 with contacts on both short sides, thermal
/// speed 2 (mean speed 2√(π/2) ≈ 2.507), frames 1 apart, steps 0.1. Taken: an unsampled
/// carrier at (1, 1) goes straight to the left face, (0, 1), in 1/2.507 = 0.399: at
/// t = 1.2 halfway, at 1.5 gone; a sampled one follows its two recorded steps (at 1.15
/// halfway between them, 9.8) and is gone after them. Given off: an unsampled carrier
/// at (1.5, 1) in the next frame leaves the face at 2 − 1.5/2.507 = 1.402 (absent at
/// 1.3, halfway at 1.701); a sampled one from its first recorded step, 1.7 (absent at
/// 1.65, between its first two steps at 1.75).
#[test]
fn contacts_take_and_give_off_carriers_between_frames() {
    let frame = |t: f64, x: Vec<[f32; 2]>, id: Vec<u32>, paths| TubeFrame {
        t,
        q: vec![-1.0; x.len()],
        kind: vec![KIND_NEGATIVE; x.len()],
        x,
        id,
        paths,
        ..TubeFrame::default()
    };
    let run = TubeRun {
        frames: vec![
            frame(1.0, vec![[1.0, 1.0], [9.5, 1.0]], vec![1, 14], Vec::new()),
            frame(
                2.0,
                vec![[1.5, 1.0], [8.0, 1.0]],
                vec![3, 21],
                vec![
                    (14, vec![[9.7, 1.0], [9.9, 1.0]]),
                    (21, vec![[9.6, 1.0], [9.0, 1.0], [8.5, 1.0], [8.0, 1.0]]),
                ],
            ),
        ],
        metals: vec![MetalView {
            lo: DVec3::ZERO,
            hi: DVec3::new(10.0, 2.0, 0.0),
            thermal_speed: 2.0,
            contacts: vec![[[0.0, 0.0], [0.0, 2.0]], [[10.0, 0.0], [10.0, 2.0]]],
            ..MetalView::default()
        }],
        dt: 0.1,
        fine_every: 7,
        ..TubeRun::default()
    };
    let at = |t: f64, id: u32| {
        let f = run.frame_between(t).unwrap();
        f.id.iter().position(|i| *i == id).map(|k| f.x[k])
    };
    let near = |x: Option<[f32; 2]>, want: f32| {
        x.is_some_and(|x| (x[0] - want).abs() < 2e-3 && (x[1] - 1.0).abs() < 1e-6)
    };
    assert!(
        near(at(1.2, 1), 0.499),
        "taken, unsampled: {:?}",
        at(1.2, 1)
    );
    assert_eq!(at(1.5, 1), None, "taken, unsampled, gone");
    assert!(
        near(at(1.15, 14), 9.8),
        "taken, sampled: {:?}",
        at(1.15, 14)
    );
    assert_eq!(at(1.25, 14), None, "taken, sampled, gone");
    assert_eq!(at(1.3, 3), None, "given off, unsampled, not yet");
    assert!(
        near(at(1.701, 3), 0.75),
        "given off, unsampled: {:?}",
        at(1.701, 3)
    );
    assert_eq!(at(1.65, 21), None, "given off, sampled, not yet");
    assert!(
        near(at(1.75, 21), 9.3),
        "given off, sampled: {:?}",
        at(1.75, 21)
    );
}

/// The Fourier measures of a potential goal: `3 sin(ω t + φ)` plus a drift and an
/// offset, sampled uniformly over whole periods, gives the amplitude 3 and the phase φ
/// (degrees) to rounding, for phases all round the circle; a phase goal compares a value
/// at its turn nearest the goal's middle.
#[test]
fn amplitude_and_phase() {
    let omega = 2.0 * std::f64::consts::PI / 1.5;
    for phi in [-170.0_f64, -90.0, -10.0, 0.0, 45.0, 135.0, 179.0] {
        let samples: Vec<(f64, f64)> = (0..600)
            .map(|k| {
                let t = 10.0 + f64::from(k) * 0.025;
                (
                    t,
                    3.0 * libm::sin(omega * t + phi.to_radians()) + 0.2 * t - 5.0,
                )
            })
            .collect();
        let (a, p) = fourier_at(&samples, omega);
        assert!((a - 3.0).abs() < 1e-9, "amplitude {a} at {phi}");
        assert!((p - phi).abs() < 1e-7, "phase {p} for {phi}");
    }
    let goal = PotentialGoal {
        electrode: 0,
        reference: None,
        measure: Measure::Phase,
        omega,
        min: 150.0,
        max: 210.0,
        start: 0.0,
        end: 1.0,
    };
    assert!((goal.compared(-170.0) - 190.0).abs() < 1e-12);
    assert!((goal.compared(170.0) - 170.0).abs() < 1e-12);
    let mean = PotentialGoal {
        measure: Measure::Mean,
        ..goal
    };
    assert!((mean.compared(-170.0) + 170.0).abs() < 1e-12);
}

/// The self-field's grid (`self_field_grid`): it covers the board and `SELF_FIELD_MARGIN` round
/// it, `SELF_FIELD_CELL` apart, its sides rounded up to powers of two, centred on the board.
#[test]
fn self_field_grid_covers_the_board() {
    let d = self_field_grid(40.0, 20.0, 3.0);
    assert_eq!((d.nx + 1, d.ny + 1), (256, 128));
    assert!((d.h - SELF_FIELD_CELL).abs() < 1e-15);
    #[allow(clippy::cast_precision_loss)]
    let far = (
        d.origin.0 + (d.nx + 1) as f64 * d.h,
        d.origin.1 + (d.ny + 1) as f64 * d.h,
    );
    assert!(d.origin.0 <= -SELF_FIELD_MARGIN && d.origin.1 <= -SELF_FIELD_MARGIN);
    assert!(far.0 >= 40.0 + SELF_FIELD_MARGIN && far.1 >= 20.0 + SELF_FIELD_MARGIN);
    assert!(((d.origin.0 + far.0) / 2.0 - 20.0).abs() < 1e-12);
    assert!(((d.origin.1 + far.1) / 2.0 - 10.0).abs() < 1e-12);
}

/// A tube level's speed of light: its issues (the step at most `SELF_FIELD_MAX_STEP`, the
/// speed within `SPEED_OF_LIGHT_RANGE`), and the engine's tube built with the self-field
/// exactly when the level has one.
#[test]
fn speed_of_light_turns_the_self_field_on() {
    let mut level = crate::shipped("ohms_law");
    let issues = |l: &Level| {
        let mut out = Vec::new();
        l.tube_issues(&l.reference_solution, &mut out);
        out
    };
    assert!(issues(&level).is_empty(), "{:?}", issues(&level));
    assert!(
        level
            .tube_sim(&level.reference_solution, 1)
            .tube
            .darwin
            .is_none()
    );
    let spec = level.tube.as_mut().expect("a tube level");
    spec.speed_of_light = Some(3.0);
    assert!(issues(&level).iter().any(|i| i.contains("step is at most")));
    level.tube.as_mut().expect("a tube level").step = 0.01;
    assert!(issues(&level).is_empty(), "{:?}", issues(&level));
    let sim = level.tube_sim(&level.reference_solution, 1);
    let d = sim.tube.darwin.as_ref().expect("the self-field");
    assert!((d.mu0 - 4.0 * std::f64::consts::PI / 9.0).abs() < 1e-12);
    level.tube.as_mut().expect("a tube level").speed_of_light = Some(0.1);
    assert!(issues(&level).iter().any(|i| i.contains("speed of light")));
}
