//! Beams of interacting particles (PHYSICS.md §3.3).
//!
//! `cargo test --release -p physics --test beams -- --nocapture --test-threads=1`

mod common;

use common::Rng;
use physics::DVec3;
use physics::beam::{BeamParticle, BeamScenario, run_beam};
use physics::dynamics::Particle;
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Aabb, Region, Shape, Sphere};
use physics::trajectory::{MARGIN_SAFE, Outcome, RunSettings, Scenario, run};
use physics::verify::classify;

const TOL: f64 = 1e-12;

fn particle(charge: f64, mass: f64) -> Particle {
    Particle {
        charge,
        mass,
        radius: 0.0,
        moment: 0.0,
    }
}

fn beam(
    field: Coulomb,
    obstacles: Vec<Shape>,
    particles: Vec<BeamParticle>,
    interact: bool,
    t_max: f64,
) -> BeamScenario<Coulomb> {
    BeamScenario {
        field,
        obstacles,
        particles,
        c: f64::INFINITY,
        detector: None,
        acceptance: None,
        bounds: None,
        t_max,
        interact,
    }
}

/// B1: two interacting particles. The relative coordinate r = x_A − x_B obeys the
/// one-body problem of reduced mass μ = m_A m_B / (m_A + m_B) in the field of a fixed
/// charge q_B, solved independently by the (validated) single-particle runner; the
/// centre of mass moves uniformly.
#[test]
fn b1_two_body_reduces_to_one_body() {
    let (qa, ma, qb, mb) = (0.8, 1.0, 1.3, 3.0);
    let (xa, va) = (DVec3::new(-6.0, 0.7, 0.0), DVec3::new(1.1, 0.05, 0.0));
    let (xb, vb) = (DVec3::new(2.0, -0.2, 0.0), DVec3::new(-0.3, 0.1, 0.0));
    let t_max = 20.0;
    let scn = beam(
        Coulomb::new(&[]),
        vec![],
        vec![
            BeamParticle {
                particle: particle(qa, ma),
                x0: xa,
                p0: va * ma,
            },
            BeamParticle {
                particle: particle(qb, mb),
                x0: xb,
                p0: vb * mb,
            },
        ],
        true,
        t_max,
    );
    let run2 = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    let (a, b) = (&run2.trajectories[0].end, &run2.trajectories[1].end);

    let mu = ma * mb / (ma + mb);
    let one = Scenario {
        field: Coulomb::new(&[FixedCharge {
            position: DVec3::ZERO,
            charge: qb,
            radius: 0.0,
        }]),
        obstacles: vec![],
        particle: particle(qa, mu),
        c: f64::INFINITY,
        x0: xa - xb,
        p0: (va - vb) * mu,
        detector: None,
        bounds: None,
        t_max,
        radiation_reaction: false,
        acceptance: None,
    };
    let r1 = run(&one, &RunSettings::with_tolerance(TOL)).end;
    let rel = ((a.x - b.x) - r1.x).length() / r1.x.length();
    let cm0 = (xa * ma + xb * mb) / (ma + mb);
    let vcm = (va * ma + vb * mb) / (ma + mb);
    let cm = (a.x * ma + b.x * mb) / (ma + mb);
    let cm_err = (cm - (cm0 + vcm * t_max)).length();
    let closest = (xa - xb).length();
    println!(
        "B1: relative position differs from the one-body solution by {rel:.1e}; centre of \
         mass off uniform motion by {cm_err:.1e}; energy drift {:.1e} (start distance {closest:.1})",
        run2.energy_max_rel_error
    );
    assert_eq!(a.t.to_bits(), t_max.to_bits(), "flies to the end");
    assert!(rel < 1e-9 && cm_err < 1e-9 && run2.energy_max_rel_error < 1e-10);
}

/// B2: a Coulomb explosion of 12 like charges: energy (kinetic + pair interaction),
/// momentum and angular momentum are conserved.
#[test]
fn b2_coulomb_explosion_conserves_energy_momentum_angular_momentum() {
    let mut rng = Rng::new(7);
    let particles: Vec<BeamParticle> = (0..12)
        .map(|_| {
            let x = DVec3::new(rng.next_f64() * 4.0 - 2.0, rng.next_f64() * 4.0 - 2.0, 0.0);
            let v = DVec3::new(rng.next_f64() - 0.5, rng.next_f64() - 0.5, 0.0) * 0.2;
            let m = 0.5 + rng.next_f64();
            BeamParticle {
                particle: particle(0.3, m),
                x0: x,
                p0: v * m,
            }
        })
        .collect();
    let p_tot0: DVec3 = particles.iter().map(|b| b.p0).sum();
    let l0: f64 = particles.iter().map(|b| b.x0.cross(b.p0).z).sum();
    let scn = beam(Coulomb::new(&[]), vec![], particles, true, 15.0);
    let r = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    let p_tot: DVec3 = r.trajectories.iter().map(|t| t.end.p).sum();
    let l: f64 = r
        .trajectories
        .iter()
        .map(|t| t.end.x.cross(t.end.p).z)
        .sum();
    let p_scale: f64 = r.trajectories.iter().map(|t| t.end.p.length()).sum();
    let l_scale: f64 = r
        .trajectories
        .iter()
        .map(|t| t.end.x.cross(t.end.p).length())
        .sum();
    let dp = (p_tot - p_tot0).length() / p_scale;
    let dl = (l - l0).abs() / l_scale;
    println!(
        "B2: energy drift {:.1e}, momentum {dp:.1e}, angular momentum {dl:.1e} ({} steps)",
        r.energy_max_rel_error, r.stats.n_step
    );
    assert!(r.energy_max_rel_error < 1e-10 && dp < 1e-12 && dl < 1e-12);
}

/// A scene with a fixed charge (an obstacle), bounds and a detector, and three
/// particles with different fates.
fn scene(interact: bool) -> BeamScenario<Coulomb> {
    let charge = FixedCharge {
        position: DVec3::new(10.0, 0.0, 0.0),
        charge: -2.0,
        radius: 0.3,
    };
    let launch = |y: f64, vy: f64, v: f64| BeamParticle {
        particle: particle(0.5, 1.0),
        x0: DVec3::new(0.0, y, 0.0),
        p0: DVec3::new(v, vy, 0.0),
    };
    BeamScenario {
        field: Coulomb::new(&[charge]),
        obstacles: vec![Shape::Sphere(Sphere {
            center: charge.position,
            radius: charge.radius,
        })],
        particles: vec![
            launch(0.0, 0.0, 1.0),   // straight into the charge
            launch(3.0, 0.05, 1.2),  // bent towards it, on to the detector
            launch(-3.0, -0.4, 1.0), // out of bounds
        ],
        c: f64::INFINITY,
        detector: Some(Region::Box(Aabb {
            min: DVec3::new(20.0, -2.0, -1.0),
            max: DVec3::new(22.0, 4.0, 1.0),
        })),
        acceptance: None,
        bounds: Some(Aabb {
            min: DVec3::new(-1.0, -6.0, -1.0),
            max: DVec3::new(23.0, 8.0, 1.0),
        }),
        t_max: 60.0,
        interact,
    }
}

/// B3: without interaction a beam is its particles flown one by one: outcomes, event
/// times and all margins (closest approaches, and the penetration depths that the
/// ghosts follow past removal) agree with the single-particle runner.
#[test]
fn b3_non_interacting_beam_matches_single_flights() {
    let mut scn = scene(false);
    // Add a particle that grazes the fixed charge (closest approach between 0.02 and
    // 0.2 cells, where margins are refined), found by scanning the launch height.
    let single = |b: &BeamParticle| Scenario {
        field: scn.field.clone(),
        obstacles: scn.obstacles.clone(),
        particle: b.particle,
        c: scn.c,
        x0: b.x0,
        p0: b.p0,
        detector: scn.detector,
        bounds: scn.bounds,
        t_max: scn.t_max,
        radiation_reaction: false,
        acceptance: None,
    };
    let grazing = (0..200)
        .map(|k| BeamParticle {
            particle: particle(0.5, 1.0),
            x0: DVec3::new(0.0, 1.0 + 0.01 * f64::from(k), 0.0),
            p0: DVec3::new(1.2, -0.05, 0.0),
        })
        .find(|b| {
            let t = run(&single(b), &RunSettings::with_tolerance(TOL));
            let m = t.margins.unwrap().obstacles[0];
            t.outcome != Outcome::Collided(0) && (0.02..0.2).contains(&m)
        })
        .expect("a grazing launch");
    scn.particles.push(grazing);
    let r = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    let mut worst: f64 = 0.0;
    for (i, (b, t)) in scn.particles.iter().zip(&r.trajectories).enumerate() {
        let one = Scenario {
            field: scn.field.clone(),
            obstacles: scn.obstacles.clone(),
            particle: b.particle,
            c: scn.c,
            x0: b.x0,
            p0: b.p0,
            detector: scn.detector,
            bounds: scn.bounds,
            t_max: scn.t_max,
            radiation_reaction: false,
            acceptance: None,
        };
        let s = run(&one, &RunSettings::with_tolerance(TOL));
        assert_eq!(t.outcome, s.outcome, "particle {i}");
        let (mt, ms) = (t.margins.as_ref().unwrap(), s.margins.as_ref().unwrap());
        // Margins beyond MARGIN_SAFE are lower bounds, deliberately not refined (and
        // ignored by the verification): compare them clamped.
        let clamp = |m: f64| m.clamp(-MARGIN_SAFE, MARGIN_SAFE);
        let dm = mt
            .all()
            .zip(ms.all())
            .map(|(a, b)| (clamp(a) - clamp(b)).abs())
            .fold(0.0, f64::max);
        let dt = (t.end.t - s.end.t).abs();
        println!(
            "B3 particle {i}: {:?} at t = {:.4}; event time differs by {dt:.1e}, margins by {dm:.1e}",
            t.outcome, t.end.t
        );
        worst = worst.max(dt).max(dm);
    }
    let fates: Vec<Outcome> = r.trajectories.iter().map(|t| t.outcome).collect();
    assert_eq!(
        &fates[..3],
        &[Outcome::Collided(0), Outcome::Arrived, Outcome::LeftBounds]
    );
    assert!(worst < 1e-8, "{worst:.3e}");
}

/// B4: the interacting version is deterministic (bit-identical reruns), each particle's
/// outcome is verified by the usual preview/verify comparison, and the interaction
/// changes the flights measurably.
#[test]
fn b4_interacting_beam_is_deterministic_and_verified() {
    let scn = scene(true);
    let a = run_beam(&scn, &RunSettings::with_tolerance(1e-10));
    let b = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    let b2 = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    for (x, y) in b.trajectories.iter().zip(&b2.trajectories) {
        assert_eq!(x.end.x.x.to_bits(), y.end.x.x.to_bits());
        assert_eq!(x.end.p.y.to_bits(), y.end.p.y.to_bits());
    }
    let free = run_beam(&scene(false), &RunSettings::with_tolerance(TOL));
    for (i, (pa, pb)) in a.trajectories.iter().zip(&b.trajectories).enumerate() {
        let status = classify(pa, pb, scn.t_max);
        let shift = (pb.end.x - free.trajectories[i].end.x).length();
        println!(
            "B4 particle {i}: {:?}, {status:?}; interaction moved its end point by {shift:.3}",
            pb.outcome
        );
        assert!(status.is_verified(), "particle {i}: {status:?}");
    }
    println!(
        "B4: energy drift {:.1e}, {} restarts",
        b.energy_max_rel_error, b.restarts
    );
    assert!(b.energy_max_rel_error < 1e-10);
}
