//! Checks on the shipped levels.

use std::path::PathBuf;

use level::Level;
use physics::trajectory::Outcome;

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
        // Tube levels include the particles' induced charge (PHYSICS.md §2.11).
        if !level.has_metal(&level.reference_solution) || level.is_tube() {
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

/// Circuits drive plates and coils one way (PHYSICS.md §2.10): in every shipped level
/// the bound on the particles' neglected back-action, `Level::drive_back_action_bound`,
/// is below `IMAGE_FORCE_LIMIT` (1e-10), as for the neglected image force.
#[test]
fn drive_back_action_is_negligible() {
    for (name, level) in shipped_levels() {
        if !level.has_drives() {
            continue;
        }
        let bound = level.drive_back_action_bound(&level.reference_solution);
        println!("{name}: circuit back-action bound {bound:.1e}");
        assert!(bound < level::IMAGE_FORCE_LIMIT, "{name}: {bound:.3e}");
    }
}

/// No shipped level has geometry the physics cannot compute (`Level::setup_issues`):
/// metal apart, no particle starting in metal or on a wire.
#[test]
fn shipped_setups_are_computable() {
    for (name, level) in shipped_levels() {
        let issues = level.setup_issues(&level.reference_solution);
        assert!(issues.is_empty(), "{name}: {issues:?}");
    }
}

/// Voltages are set with sliders (the owner, 2026-10-06): each power supply's and
/// plate's range is symmetric about 0, and no reference potential lies where a player
/// would try first: within 5 % of the range's half from 0, from either end or from the
/// middle of either half (`scripts/levels.py`, `guessable`; `generator window` shows
/// each level's whole solving interval, which the level designs keep clear of them).
#[test]
fn slider_ranges_hide_the_answer() {
    use level::ElementKind;
    for (name, level) in shipped_levels() {
        let l = &level.limits;
        for (what, list) in [
            ("supply voltages", &l.supply_voltages),
            ("coil rates", &l.coil_rates),
            ("plate voltages", &l.plate_voltages),
        ] {
            if let Some((lo, hi)) = level::value_range(list) {
                // Exactly symmetric (a grounded-only list [0] has lo = 0 = -hi).
                #[allow(clippy::float_cmp)]
                let symmetric = lo == -hi;
                assert!(
                    symmetric,
                    "{name}: {what} slider [{lo}, {hi}] not symmetric"
                );
            }
        }
        // Each reference value against its own slider: a plate's, or the list of what its
        // power supply operates (an electrode's potentials, a coil's ramp rates).
        let targets = level.supply_targets();
        for e in &level.reference_solution {
            let list: &[f64] = match e.kind {
                ElementKind::Plate => &l.plate_voltages,
                ElementKind::Supply => targets
                    .iter()
                    .find(|(c, _)| *c == e.node)
                    .map_or(&[][..], |&(_, t)| level.supply_list(t)),
                _ => continue,
            };
            let Some((_, hi)) = level::value_range(list) else {
                continue;
            };
            if hi == 0.0 {
                continue; // nothing to tune (grounded plates only)
            }
            let v = e.value.abs();
            for p in [0.0, 0.5 * hi, hi] {
                assert!(
                    (v - p).abs() >= 0.05 * hi,
                    "{name}: {:?} {} near {p} on the slider ±{hi}",
                    e.kind,
                    e.value
                );
            }
        }
    }
}

/// A player's supply switched off leaves its electrode's circuit driving nothing: a DC
/// source at 0, a sine with no swing (or, where the supply sets its offset, none), a pulse
/// whose top is its rest. The owner (2026-10-26): Ringing tank's step, the Inverter's pulse
/// and the Rectifier's and the Klystron's sources had their answers built into the circuit
/// a switched-off supply leaves in place, so they solved (or half solved) with nothing
/// done. Switched on, a supply starts at 0, which the slider rule keeps out of every
/// solving window (`slider_ranges_hide_the_answer`).
#[test]
fn switched_off_supplies_drive_nothing() {
    use level::Source;
    let idle = |s: &Source, offset_is_level: bool| match *s {
        Source::Dc { value } => value == 0.0,
        Source::Sine {
            offset, amplitude, ..
        } => {
            if offset_is_level {
                offset == 0.0
            } else {
                amplitude == 0.0
            }
        }
        // Exactly: a pulse whose top is its rest, as written in the level.
        #[allow(clippy::float_cmp)]
        Source::Pulse { low, high, .. } => high == low,
    };
    let mut checked = 0;
    for (name, level) in shipped_levels() {
        for (k, e) in level.electrodes.iter().enumerate() {
            if !e.tunable {
                continue;
            }
            if let Some(d) = &e.drive {
                checked += 1;
                assert!(
                    idle(&d.source, false),
                    "{name}: electrode {}'s drive {:?} with its supply off",
                    k + 1,
                    d.source
                );
            }
            for c in level.tube.iter().flat_map(|t| &t.circuits) {
                if c.electrode == k {
                    checked += 1;
                    assert!(
                        idle(&c.source, c.supply_offset),
                        "{name}: electrode {}'s circuit {:?} with its supply off",
                        k + 1,
                        c.source
                    );
                }
            }
        }
    }
    assert!(checked > 10, "only {checked} tunable circuits found");
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
