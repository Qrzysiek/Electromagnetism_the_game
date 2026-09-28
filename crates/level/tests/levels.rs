//! Checks on the shipped levels.

use std::path::PathBuf;

use level::Level;
use physics::trajectory::{Outcome, RunSettings, run};

fn shipped_levels() -> Vec<(String, Level)> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut out: Vec<(String, Level)> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| level::is_level_file(p))
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
        if level.has_beams() {
            // Every beam shot reaches its verified transmission in every flight.
            let v = level.verify_beams(&level.reference_solution);
            for (d, f) in v.iter().enumerate() {
                for s in 0..level.shots.len() {
                    let (ok, n) = f.transmitted(s);
                    println!("{name} flight {d} shot {s}: {ok}/{n} arrive, verified");
                }
            }
            assert!(level.beams_solved(&v), "{name}");
            // The model neglects each particle's radiation (PHYSICS.md §3.3): the energy it
            // radiates (Liénard, with the other particles' fields) must be below the
            // numerical accuracy, as for single flights.
            // Quasi-static interaction at finite c (PHYSICS.md §3.3): the exact retarded
            // interaction must give the same outcome for every particle, and the difference
            // between the models must match the indicator. Both are flown without radiation
            // reaction (it would multiply the cost of the exact run eightfold), so that the
            // comparison measures the interaction models only: comparing the level's own
            // flights (with radiation reaction) against an exact run without it mixed in
            // the radiation reaction's effect (found on level 51, where it dominated).
            if level.physics.c.is_some()
                && level.physics.beam_interaction
                && !level.physics.beam_retarded
            {
                let rs = RunSettings::with_tolerance(level.physics.tolerances.preview);
                let flights = |retarded: bool, interact: bool| {
                    let mut l = level.clone();
                    l.physics.beam_retarded = retarded;
                    l.physics.beam_interaction = interact;
                    l.physics.radiation_reaction = false;
                    l.beam_scenarios(
                        &level.reference_solution,
                        physics::conductor::Resolution::Preview,
                    )
                    .iter()
                    .map(|scn| physics::beam::run_beam(scn, &rs))
                    .collect::<Vec<_>>()
                };
                let (quasi, exact, alone) = (
                    flights(false, true),
                    flights(true, true),
                    flights(false, false),
                );
                for (((f, q), r), lone) in v.iter().zip(&quasi).zip(&exact).zip(&alone) {
                    let effect = q
                        .trajectories
                        .iter()
                        .zip(&lone.trajectories)
                        .map(|(a, b)| (a.end.x - b.end.x).length())
                        .fold(0.0f64, f64::max);
                    let mut worst: f64 = 0.0;
                    for (i, ((a, b), c)) in f
                        .verified
                        .trajectories
                        .iter()
                        .zip(&r.trajectories)
                        .zip(&q.trajectories)
                        .enumerate()
                    {
                        assert_eq!(a.outcome, b.outcome, "{name}: particle {i}");
                        worst = worst.max((c.end.x - b.end.x).length());
                    }
                    let indicator = q
                        .neglected_retardation
                        .iter()
                        .fold(0.0f64, |m, &x| m.max(x));
                    println!(
                        "{name}: quasi-static and retarded end points differ by up to {worst:.1e} \
                         cells, {:.1e} of the interaction's effect {effect:.2}; indicator up to \
                         {indicator:.1e}",
                        worst / effect
                    );
                    // The indicator is an order-of-magnitude estimate of the relative error.
                    assert!(
                        worst / effect < 3.0 * indicator,
                        "{name}: indicator too small"
                    );
                }
            }
            // With radiation reaction the Landau–Lifshitz treatment must be valid instead.
            if level.physics.c.is_some() {
                for f in &v {
                    for t in &f.verified.trajectories {
                        if level.physics.radiation_reaction {
                            let r = t.reaction_ratio_max;
                            assert!(r < 0.05, "{name}: |F_RR|/|F_L| = {r:.3e}");
                        } else {
                            let fraction = t.radiated_energy / t.kinetic_initial;
                            assert!(fraction < 1e-10, "{name}: radiated fraction {fraction:.3e}");
                        }
                    }
                }
            }
            continue;
        }
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
            // Charge (Liénard) plus the estimate for a magnetic moment (PHYSICS.md §3.2).
            let moment =
                level::moment_radiation_estimate(scn, tr.samples.iter().map(|s| (s.x, s.p, s.t)));
            let fraction = (tr.radiated_energy + moment) / tr.kinetic_initial;
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

/// Levels that include radiation reaction must need it (otherwise the radiation is
/// decoration): with radiation switched off, either the reference solution no longer
/// solves them, or the radiation that would be neglected is above the 1e-10 of the
/// launch energy allowed for neglected radiation, so the model needs it to be consistent
/// (beams; and levels whose goal is the radiation itself, where a particle that radiates
/// enough to be measured must also feel it).
#[test]
fn radiation_levels_need_radiation() {
    for (name, mut level) in shipped_levels() {
        if !level.physics.radiation_reaction {
            continue;
        }
        level.physics.radiation_reaction = false;
        if level.has_beams() {
            let v = level.verify_beams(&level.reference_solution);
            let solved = level.beams_solved(&v);
            let radiated = v
                .iter()
                .flat_map(|f| &f.verified.trajectories)
                .map(|t| t.radiated_energy / t.kinetic_initial)
                .fold(0.0f64, f64::max);
            println!(
                "{name}: without radiation reaction solved {solved}, radiates up to {radiated:.1e}"
            );
            assert!(
                !solved || radiated >= 1e-10,
                "{name} needs no radiation reaction"
            );
            continue;
        }
        let flights = level.verify_flights(&level.reference_solution);
        let all_arrive = flights.iter().all(|v| v.outcome() == Outcome::Arrived);
        let radiated = flights
            .iter()
            .map(|v| v.verified.radiated_energy / v.verified.kinetic_initial)
            .fold(0.0f64, f64::max);
        println!(
            "{name}: without radiation reaction solved {all_arrive}, radiates up to {radiated:.1e}"
        );
        assert!(
            !all_arrive || radiated >= 1e-10,
            "{name} is solved without radiation reaction, which it does not need"
        );
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
/// electrodes (PHYSICS.md §2.7), so it must be negligible along every reference flight:
/// `level::electrode_image_force_bound` < `IMAGE_FORCE_LIMIT` (1e-10).
#[test]
fn electrode_image_force_is_negligible() {
    for (name, level) in shipped_levels() {
        if !level.has_metal(&level.reference_solution) {
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
            let t0 = level.shots[level.flight_of(i).0].launch.kinetic_energy;
            let points = tr.samples.iter().map(|s| (s.x, s.t));
            worst = worst.max(level::electrode_image_force_bound(scn, t0, points));
        }
        println!("{name}: image force bound / max(|F|, F0) <= {worst:.1e}");
        assert!(worst < level::IMAGE_FORCE_LIMIT, "{name}: {worst:.3e}");
    }
}

/// Realistic iterations (chapter 15, SPEC §3) build an earlier level's idealised solution
/// in and add a real effect that breaks it: on its own, without the player's elements,
/// the built-in design must not solve the level (otherwise the effect is decoration).
#[test]
fn realistic_levels_break_the_idealised_design() {
    const REALISTIC: &[&str] = &[
        "chromatic_aberration",
        "real_analyzer",
        "crt_earth_field",
        "beam_pipe",
        "calutron_space_charge",
        "soft_landing_current",
    ];
    let mut found = 0;
    for (name, level) in shipped_levels() {
        if !REALISTIC.iter().any(|r| name.ends_with(r)) {
            continue;
        }
        found += 1;
        let solved = if level.has_beams() {
            level.beams_solved(&level.verify_beams(&[]))
        } else {
            level
                .verify_flights(&[])
                .iter()
                .all(|v| v.outcome() == Outcome::Arrived && v.status.is_verified())
        };
        println!("{name}: idealised design alone solves it: {solved}");
        assert!(
            !solved,
            "{name}: the idealised design still works on its own"
        );
    }
    assert_eq!(found, REALISTIC.len(), "all realistic levels shipped");
}
