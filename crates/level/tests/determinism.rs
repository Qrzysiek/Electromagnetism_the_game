//! T11 (PHYSICS.md §9): trajectories are bit-identical across runs and platforms.
//!
//! For every level in `levels/`, the reference solution is flown at both tolerances and
//! every bit of the result (all samples, margins, outcome, step statistics) is hashed.
//! The hashes are compared with `levels/golden_hashes.json`, committed from one platform
//! and checked on every CI platform. Regenerate with `UPDATE_GOLDEN=1 cargo test -p level`
//! (only after an intentional change to the physics, which must bump the engine version).

use std::collections::BTreeMap;
use std::path::PathBuf;

use level::Level;
use physics::trajectory::Trajectory;

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

fn levels_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels")
}

/// Every level's file (`GOLDEN_ONLY=<part of a name>`: only the matching levels, a new
/// level's hash without the hour the whole set takes), sorted.
fn level_paths() -> Vec<PathBuf> {
    let only = std::env::var("GOLDEN_ONLY").ok();
    let mut paths: Vec<_> = std::fs::read_dir(levels_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| level::is_level_file(p))
        .filter(|p| {
            only.as_ref().is_none_or(|o| {
                p.file_stem()
                    .is_some_and(|s| s.to_string_lossy().contains(o.as_str()))
            })
        })
        .collect();
    paths.sort();
    paths
}

/// A level's file name (its stem: the golden file's key).
fn stem(p: &std::path::Path) -> String {
    p.file_stem().unwrap().to_string_lossy().into_owned()
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

/// The levels run twice must hash alike, each level checked again only where its file or
/// this test program changed since it last passed (`level::run_cache::TestPasses`).
#[test]
fn t11_repeated_runs_are_identical() {
    let passes = level::run_cache::TestPasses::open("t11-repeated");
    for p in level_paths() {
        let bytes = std::fs::read(&p).unwrap();
        if passes.as_ref().and_then(|t| t.passed(&bytes)).is_some() {
            continue;
        }
        let level = Level::from_json(std::str::from_utf8(&bytes).unwrap()).unwrap();
        let (a, b) = (level_hash(&level), level_hash(&level));
        assert_eq!(a, b, "{}", stem(&p));
        if let Some(t) = &passes {
            t.record(&bytes, &a);
        }
    }
}

#[test]
fn t11_matches_golden_hashes() {
    // A level's hash computed under this test program for its file, kept: the runs are
    // deterministic (`level::run_cache::TestPasses`).
    let passes = level::run_cache::TestPasses::open("t11-golden");
    let files: Vec<(String, Vec<u8>)> = level_paths()
        .iter()
        .map(|p| (stem(p), std::fs::read(p).unwrap()))
        .collect();
    let mut hashes = BTreeMap::new();
    for (name, bytes) in &files {
        let hash = passes
            .as_ref()
            .and_then(|t| t.passed(bytes))
            .unwrap_or_else(|| {
                let level = Level::from_json(std::str::from_utf8(bytes).unwrap()).unwrap();
                level_hash(&level)
            });
        hashes.insert(name.clone(), hash);
    }
    let record = |hashes: &BTreeMap<String, String>| {
        if let Some(t) = &passes {
            for (name, bytes) in &files {
                t.record(bytes, &hashes[name]);
            }
        }
    };
    let golden_path = levels_dir().join("golden_hashes.json");
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        // With GOLDEN_ONLY, the other levels' hashes are kept as they are.
        if std::env::var_os("GOLDEN_ONLY").is_some() {
            let old = std::fs::read_to_string(&golden_path).unwrap_or_default();
            for line in old.lines() {
                let parts: Vec<&str> = line.split('"').collect();
                if parts.len() >= 5 {
                    hashes
                        .entry(parts[1].to_string())
                        .or_insert_with(|| parts[3].to_string());
                }
            }
        }
        let mut s = String::from("{\n");
        let body: Vec<String> = hashes
            .iter()
            .map(|(k, v)| format!("  \"{k}\": \"{v}\""))
            .collect();
        s += &body.join(",\n");
        s += "\n}\n";
        std::fs::write(&golden_path, s).unwrap();
        record(&hashes);
        return;
    }
    let golden = std::fs::read_to_string(&golden_path).expect("golden_hashes.json missing");
    for (name, hash) in &hashes {
        let expected = format!("\"{name}\": \"{hash}\"");
        println!("{name}: {hash}");
        assert!(
            golden.contains(&expected),
            "{name}: hash {hash} not in golden file (platform-dependent result or physics change)"
        );
    }
    record(&hashes);
}
