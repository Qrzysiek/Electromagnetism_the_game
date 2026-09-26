//! Checks on the shipped levels.

use std::path::PathBuf;

use level::Level;
use physics::trajectory::{Outcome, RunSettings, run};

fn shipped_levels() -> Vec<(String, Level)> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut out: Vec<(String, Level)> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| {
            p.extension().is_some_and(|e| e == "json") && !p.ends_with("golden_hashes.json")
        })
        .map(|p| {
            let name = p.file_stem().unwrap().to_string_lossy().into_owned();
            (
                name,
                Level::from_json(&std::fs::read_to_string(&p).unwrap()).unwrap(),
            )
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[test]
fn reference_solutions_are_verified() {
    for (name, level) in shipped_levels() {
        assert_eq!(
            level.check_placement(&level.reference_solution),
            Ok(()),
            "{name}"
        );
        for (i, v) in level
            .verify_flights(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let (shot, d) = level.flight_of(i);
            let flight = format!("{name} shot {} disturbance {}", shot + 1, d + 1);
            assert_eq!(v.outcome(), Outcome::Arrived, "{flight}");
            assert!(v.status.is_verified(), "{flight}: {:?}", v.status);
            // Resource budget (level::cost), machine-independent part: steps per flight.
            #[allow(clippy::cast_precision_loss)]
            let steps = v.preview.stats.n_step.max(v.verified.stats.n_step) as f64;
            assert!(steps <= level::cost::STEPS_LIMIT, "{flight}: {steps} steps");
        }
    }
}

/// The model neglects radiation (PHYSICS.md §8). A level is only physically consistent if
/// the energy radiated along its reference flight, by the Liénard formula, is below the
/// numerical accuracy. The radiated fraction per close pass is about r_cl / r with
/// r_cl = q²/(mc²), which is why levels use weakly charged particles and strongly charged
/// electrodes.
#[test]
fn neglected_radiation_is_below_numerical_accuracy() {
    for (name, level) in shipped_levels() {
        for (i, scn) in level
            .scenarios(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let tr = run(
                scn,
                &RunSettings::with_tolerance(level.physics.tolerances.verify),
            );
            let fraction = tr.radiated_energy / tr.kinetic_initial;
            if level.physics.radiation_reaction {
                // Radiation is part of the model; the Landau–Lifshitz treatment must be
                // valid instead (radiation reaction small against the Lorentz force).
                assert!(
                    tr.reaction_ratio_max < 0.05,
                    "{name}: |F_RR|/|F_L| = {:.3e}",
                    tr.reaction_ratio_max
                );
                continue;
            }
            let (shot, d) = level.flight_of(i);
            println!(
                "{name} shot {} disturbance {}: radiated / T0 = {fraction:.2e}",
                shot + 1,
                d + 1
            );
            assert!(fraction < 1e-10, "{name}: radiated fraction {fraction:.3e}");
        }
    }
}

/// Levels that include radiation reaction must need it: with radiation switched off, the
/// reference solution must not solve them (otherwise the radiation is decoration).
#[test]
fn radiation_levels_need_radiation() {
    for (name, mut level) in shipped_levels() {
        if !level.physics.radiation_reaction {
            continue;
        }
        level.physics.radiation_reaction = false;
        let all_arrive = level
            .verify_flights(&level.reference_solution)
            .iter()
            .all(|v| v.outcome() == Outcome::Arrived);
        assert!(!all_arrive, "{name} is solved without radiation reaction");
    }
}

/// Levels with metal spheres: a consistent model (no unsupported combinations), the
/// conductor model accurate at verification resolution for the reference placement
/// (boundary residual < 1e-10, PHYSICS.md §2.6), and its matrices within the memory
/// budget (level::cost).
#[test]
fn metal_levels_are_accurate_and_consistent() {
    for (name, level) in shipped_levels() {
        assert_eq!(level.model_issues(), Vec::<String>::new(), "{name}");
        if level.conductors.is_empty() {
            continue;
        }
        let (field, _) = level.field_at(
            &level.reference_solution,
            physics::conductor::Resolution::Verify,
        );
        let sources: Vec<_> = field.coulomb.charges().collect();
        let residual = field.conductors.boundary_residual(&sources);
        #[allow(clippy::cast_precision_loss)]
        let bytes = field.setup_cost().bytes as f64;
        assert!(bytes <= level::cost::MEMORY_LIMIT, "{name}: {bytes} bytes");
        println!("{name}: conductor boundary residual {residual:.1e}");
        assert!(residual < 1e-10, "{name}: residual {residual:.3e}");
    }
}

/// Levels with box electrodes: the particle's image force is not computed for
/// electrodes (PHYSICS.md §2.7), so it must be negligible. Upper bound q²/d² (four times
/// the force of a flat grounded plane at distance d, which covers concave corners),
/// relative to max(|F_Lorentz|, F₀) along every reference flight: < 1e-10. F₀ = T₀ per
/// cell is the smallest force that matters for the flight (it changes the kinetic energy
/// by T₀ over one cell); relative to the Lorentz force alone the ratio is meaningless
/// where that force passes through zero.
#[test]
fn electrode_image_force_is_negligible() {
    use physics::field::FieldSolver;
    for (name, level) in shipped_levels() {
        if level.electrodes.is_empty() {
            continue;
        }
        let (field, _) = level.field_at(
            &level.reference_solution,
            physics::conductor::Resolution::Verify,
        );
        // Resource budget (level::cost): the BEM matrices.
        #[allow(clippy::cast_precision_loss)]
        let bytes = field.setup_cost().bytes as f64;
        assert!(bytes <= level::cost::MEMORY_LIMIT, "{name}: {bytes} bytes");
        let boxes = physics::bem::Electrodes::shapes_only(level.box_electrodes()).obstacles(0.0);
        let mut worst: f64 = 0.0;
        for (i, scn) in level
            .scenarios(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let tr = physics::trajectory::run(
                scn,
                &physics::trajectory::RunSettings::with_tolerance(level.physics.tolerances.preview),
            );
            let shot = &level.shots[level.flight_of(i).0];
            let q = shot.particle.charge;
            let f0 = shot.launch.kinetic_energy;
            for s in &tr.samples {
                let d = boxes
                    .iter()
                    .map(|b| b.signed_distance(s.x))
                    .fold(f64::INFINITY, f64::min)
                    .max(1e-3);
                let f = (scn.field.sample(s.x, s.t).e.length() * q.abs()).max(f0);
                worst = worst.max(q * q / (d * d) / f);
            }
        }
        println!("{name}: image force bound / max(|F|, F0) <= {worst:.1e}");
        assert!(worst < 1e-10, "{name}: {worst:.3e}");
    }
}
