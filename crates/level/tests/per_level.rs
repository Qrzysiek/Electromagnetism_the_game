//! The heavy checks of the shipped levels, one test per check and level (DEC-16), named
//! `<check>::<level file's stem>` and ignored by default: `cargo test` runs the fast tier
//! alone, `scripts/affected.py` the ones a change can affect, and
//! `cargo test --release -p level --test per_level -- --ignored [<filter>]` any by name
//! (a filter is a part of the name: `golden::`, `reference::109_`, `122_inverter`).
//!
//! - `golden::` (T11, PHYSICS.md §9): the reference solution flown at both tolerances, every
//!   bit of the result (all samples, margins, outcome, step statistics) hashed and compared
//!   with `levels/golden_hashes.json`, committed from one platform and checked on every CI
//!   platform. `UPDATE_GOLDEN=1` (optionally with `GOLDEN_ONLY=<part of a name>`) instead
//!   computes the matching levels' hashes and rewrites their entries, keeping the rest: only
//!   after an intentional change to the physics or the levels.
//! - `repeated::`: the same hash computed twice in one process.
//! - `reference::`: every shot's reference solution verified, and the level's resource
//!   budget kept.
//! - `radiation::`: the radiation the model neglects below 1e-10 of the launch energy.
//! - `needs_radiation::` (levels with radiation reaction): the reaction is needed.
//! - `image_force::` (levels with box electrodes, not tubes): the neglected image force on
//!   the particle negligible.
//!
//! A check that passed is skipped while the test program and the level's physics are
//! unchanged (`level::run_cache::TestPasses`; `EM_TEST_FULL=1` runs it anyway): a full run
//! interrupted resumes where it stopped.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use level::Level;
use libtest_mimic::{Arguments, Trial};
use physics::trajectory::{Outcome, RunSettings, Trajectory, run};

/// The checks' names: `<check>::<level>`.
const GOLDEN: &str = "golden";
const REPEATED: &str = "repeated";
const REFERENCE: &str = "reference";
const RADIATION: &str = "radiation";
const NEEDS_RADIATION: &str = "needs_radiation";

/// Trials run at once by default (`--test-threads` sets it): they share the process's
/// rayon threads, and four keep the memory under four levels' budgets (`level::cost`).
const DEFAULT_THREADS: usize = 4;

fn main() -> ExitCode {
    let mut args = Arguments::from_args();
    if std::env::var_os("UPDATE_GOLDEN").is_some() && !args.list {
        update_golden();
        return ExitCode::SUCCESS;
    }
    if std::env::var_os("EM_TEST_INPUTS").is_some() {
        print_inputs();
        return ExitCode::SUCCESS;
    }
    args.test_threads = args.test_threads.or(Some(DEFAULT_THREADS));
    let mut trials = Vec::new();
    for path in level_paths(None) {
        let name = stem(&path);
        let needs = Level::from_json(&std::fs::read_to_string(&path).unwrap())
            .unwrap()
            .physics
            .radiation_reaction;
        let mut checks = vec![GOLDEN, REPEATED, REFERENCE, RADIATION];
        if needs {
            checks.push(NEEDS_RADIATION);
        }
        for check in checks {
            let (path, name) = (path.clone(), name.clone());
            trials.push(
                Trial::test(format!("{check}::{name}"), move || {
                    run_check(check, &name, &path);
                    Ok(())
                })
                .with_ignored_flag(true),
            );
        }
    }
    // Returned, not `exit()`: the process must end normally for an instrumented build to
    // write its coverage (`scripts/affected.py`).
    libtest_mimic::run(&args, trials).exit_code()
}

/// One check of one level, skipped if it passed under this program for the level's
/// physics (`TestPasses`). The golden check takes a recorded hash for the computed one but
/// compares it with the golden file all the same: the record says what the program
/// computes, not that the file still holds it.
fn run_check(check: &str, name: &str, path: &Path) {
    let bytes = std::fs::read(path).unwrap();
    let passes = level::run_cache::TestPasses::open(pass_record(check));
    let recorded = passes.as_ref().and_then(|t| t.passed(&bytes));
    if let Some(value) = &recorded
        && check != GOLDEN
    {
        println!("{check}::{name}: passed before ({value})");
        return;
    }
    let level = Level::from_json(std::str::from_utf8(&bytes).unwrap()).unwrap();
    let value = match check {
        GOLDEN => {
            let hash = recorded.unwrap_or_else(|| level_hash(&level));
            let golden =
                std::fs::read_to_string(golden_path()).expect("golden_hashes.json missing");
            println!("{name}: {hash}");
            assert!(
                golden.contains(&format!("\"{name}\": \"{hash}\"")),
                "{name}: hash {hash} not in golden file (platform-dependent result or physics change)"
            );
            hash
        }
        REPEATED => {
            let (a, b) = (level_hash(&level), level_hash(&level));
            assert_eq!(a, b, "{name}");
            a
        }
        REFERENCE => {
            check_reference(name, &level);
            "verified".to_string()
        }
        RADIATION => {
            check_radiation(name, &level);
            "neglected".to_string()
        }
        NEEDS_RADIATION => {
            check_needs_radiation(name, &level);
            "needed".to_string()
        }
        _ => unreachable!("{check}"),
    };
    if let Some(t) = &passes {
        t.record(&bytes, &value);
    }
}

/// The pass records' names, those of the tests these checks were before (their records
/// stay valid until the program changes anyway).
fn pass_record(check: &str) -> &'static str {
    match check {
        GOLDEN => "t11-golden",
        REPEATED => "t11-repeated",
        REFERENCE => "reference",
        RADIATION => "radiation",
        NEEDS_RADIATION => "needs-radiation",
        _ => unreachable!("{check}"),
    }
}

fn levels_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels")
}

fn golden_path() -> PathBuf {
    levels_dir().join("golden_hashes.json")
}

/// The shipped levels' files, sorted; with `only`, those whose name contains it.
fn level_paths(only: Option<&str>) -> Vec<PathBuf> {
    let mut paths: Vec<_> = std::fs::read_dir(levels_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| level::is_level_file(p))
        .filter(|p| only.is_none_or(|o| stem(p).contains(o)))
        .collect();
    paths.sort();
    paths
}

/// A level's file name (its stem: the golden file's key and the trials' name).
fn stem(p: &Path) -> String {
    p.file_stem().unwrap().to_string_lossy().into_owned()
}

/// `EM_TEST_INPUTS=1`: each level's physics (`Level::physics_json`, hashed as the pass
/// records key it), one `<level> <hash>` line per level: the input other than code of its
/// checks, for `scripts/affected.py` (DEC-16). Rewording a level changes none.
fn print_inputs() {
    for path in level_paths(None) {
        let level = Level::from_json(&std::fs::read_to_string(&path).unwrap()).unwrap();
        let physics = blake3::hash(level.physics_json().as_bytes());
        println!("{} {}", stem(&path), physics.to_hex());
    }
}

/// `UPDATE_GOLDEN=1`: the hashes of every level (`GOLDEN_ONLY=<part of a name>`: of the
/// matching levels, the file keeping the others' entries) computed and written, one level
/// after another (each printed when done).
fn update_golden() {
    let only = std::env::var("GOLDEN_ONLY").ok();
    let passes = level::run_cache::TestPasses::open(pass_record(GOLDEN));
    let mut hashes = BTreeMap::new();
    let mut files = Vec::new();
    for path in level_paths(only.as_deref()) {
        let (name, bytes) = (stem(&path), std::fs::read(&path).unwrap());
        let hash = passes
            .as_ref()
            .and_then(|t| t.passed(&bytes))
            .unwrap_or_else(|| {
                level_hash(&Level::from_json(std::str::from_utf8(&bytes).unwrap()).unwrap())
            });
        println!("{name}: {hash}");
        hashes.insert(name, hash.clone());
        files.push((bytes, hash));
    }
    if only.is_some() {
        let old = std::fs::read_to_string(golden_path()).unwrap_or_default();
        for line in old.lines() {
            let parts: Vec<&str> = line.split('"').collect();
            if parts.len() >= 5 {
                hashes
                    .entry(parts[1].to_string())
                    .or_insert_with(|| parts[3].to_string());
            }
        }
    }
    let body: Vec<String> = hashes
        .iter()
        .map(|(k, v)| format!("  \"{k}\": \"{v}\""))
        .collect();
    std::fs::write(golden_path(), format!("{{\n{}\n}}\n", body.join(",\n"))).unwrap();
    if let Some(t) = &passes {
        for (bytes, hash) in &files {
            t.record(bytes, hash);
        }
    }
}

// The checks' bodies, moved unchanged from levels.rs (reference, radiation, the need for
// radiation reaction as one level's check) and determinism.rs (the hash).

/// End-point shifts (cells) that no view shows and no margin of the shipped levels comes
/// near: below it the interaction's effect is too small to compare the beam models by.
const NEGLIGIBLE_EFFECT: f64 = 1e-2;

/// One level's reference solution, verified (`reference::<level>`).
#[allow(clippy::too_many_lines)]
fn check_reference(name: &str, level: &Level) {
    assert_eq!(
        level.check_placement(&level.reference_solution),
        Ok(()),
        "{name}"
    );
    if level.is_tube() {
        // Tube levels: the verified current meets the goal (`level::tube`).
        let v = level.tube_verdict(&level.reference_solution);
        println!("{name}: tube verdict {v:?}");
        assert_eq!(v.status, level::tube::TubeStatus::Met, "{name}");
        return;
    }
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
        // Both at the verification tolerance: at the preview tolerance their integration
        // errors (5e-6 cells on Soft landing at full current) masked the difference of
        // the models.
        if level.physics.c.is_some() && level.interacts() && !level.physics.beam_retarded {
            let rs = RunSettings::with_tolerance(level.physics.tolerances.verify);
            // The interaction is switched on the flights themselves: a level with free
            // particles or free charges always interacts (`Level::interacts`), whatever
            // its `beam_interaction` says (the first version switched that flag and
            // compared such a level with itself).
            let flights = |retarded: bool, interact: bool| {
                let mut l = level.clone();
                l.physics.beam_retarded = retarded;
                l.physics.radiation_reaction = false;
                l.beam_scenarios(
                    &level.reference_solution,
                    physics::conductor::Resolution::Preview,
                )
                .into_iter()
                .map(|mut scn| {
                    scn.interact = interact;
                    physics::beam::run_beam(&scn, &rs)
                })
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
                     cells, {:.1e} of the interaction's effect {effect:.1e}; indicator up to \
                     {indicator:.1e}",
                    worst / effect
                );
                if effect < NEGLIGIBLE_EFFECT {
                    // The interaction hardly moves the end points (a pair repelling
                    // symmetrically, as in a quiet ring): the ratio says nothing about
                    // the indicator, and the difference itself must be negligible.
                    assert!(worst < NEGLIGIBLE_EFFECT, "{name}: {worst:.2e} cells");
                } else {
                    // The indicator is an order-of-magnitude estimate of the relative
                    // error.
                    assert!(
                        worst / effect < 3.0 * indicator,
                        "{name}: indicator too small"
                    );
                }
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
        return;
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

/// One level's neglected radiation (`radiation::<level>`).
fn check_radiation(name: &str, level: &Level) {
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

/// Levels that include radiation reaction must need it (otherwise the radiation is
/// decoration): with radiation switched off, either the reference solution no longer
/// solves them, or the radiation that would be neglected is above the 1e-10 of the
/// launch energy allowed for neglected radiation, so the model needs it to be consistent
/// (beams; and levels whose goal is the radiation itself, where a particle that radiates
/// enough to be measured must also feel it).
fn check_needs_radiation(name: &str, level: &Level) {
    let mut level = level.clone();
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
        return;
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

/// FNV-1a, 64-bit.
struct Hasher(u64);

impl Hasher {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn u64(&mut self, v: u64) {
        for b in v.to_le_bytes() {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    fn f64(&mut self, v: f64) {
        self.u64(v.to_bits());
    }

    fn trajectory(&mut self, t: &Trajectory) {
        self.u64(format!("{:?}", t.outcome).len() as u64);
        for b in format!("{:?}", t.outcome).bytes() {
            self.u64(u64::from(b));
        }
        for s in &t.samples {
            self.f64(s.t);
            for v in [s.x, s.p] {
                self.f64(v.x);
                self.f64(v.y);
                self.f64(v.z);
            }
        }
        self.u64(t.stats.n_fcn);
        self.u64(t.stats.n_accept);
        self.u64(t.stats.n_reject);
        self.f64(t.energy_max_abs_error);
        self.f64(t.radiated_energy);
        if let Some(m) = &t.margins {
            for v in m.all() {
                self.f64(v);
            }
        }
    }
}

/// A level's hash: both runs of its reference, every bit.
fn level_hash(level: &Level) -> String {
    {
        let mut h = Hasher::new();
        if level.is_tube() {
            // Tube levels: both runs, every frame's particles and the currents.
            for refine in [1, 2] {
                let end = level.tube_verify_end();
                let run = level.tube_run(&level.reference_solution, refine, end);
                for f in &run.frames {
                    h.f64(f.t);
                    for (x, q) in f.x.iter().zip(&f.q) {
                        h.u64(u64::from(x[0].to_bits()));
                        h.u64(u64::from(x[1].to_bits()));
                        h.u64(u64::from(q.to_bits()));
                    }
                    for s in &f.sigma {
                        h.u64(u64::from(s.to_bits()));
                    }
                }
                // Each current goal's collected charge and current. Without one, zeros
                // and NaN, as the single-goal format recorded them, so that the hashes
                // carried over to the list of goals.
                if run.collected.is_empty() {
                    for f in &run.frames {
                        h.f64(f.t);
                        h.f64(0.0);
                    }
                    h.f64(f64::NAN);
                }
                for (series, current) in run.collected.iter().zip(&run.currents) {
                    for (t, q) in series {
                        h.f64(*t);
                        h.f64(*q);
                    }
                    h.f64(current.unwrap_or(f64::NAN));
                }
                for track in &run.tracks {
                    for (t, x) in track {
                        h.f64(*t);
                        h.u64(u64::from(x[0].to_bits()));
                        h.u64(u64::from(x[1].to_bits()));
                    }
                }
            }
        } else if level.has_beams() {
            // Beam levels: every particle of every flight.
            for v in level.verify_beams(&level.reference_solution) {
                for t in v
                    .preview
                    .trajectories
                    .iter()
                    .chain(&v.verified.trajectories)
                {
                    h.trajectory(t);
                }
            }
        } else {
            for v in level.verify_flights(&level.reference_solution) {
                h.trajectory(&v.preview);
                h.trajectory(&v.verified);
            }
        }
        format!("{:016x}", h.0)
    }
}
