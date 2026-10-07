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

fn level_hashes() -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let mut paths: Vec<_> = std::fs::read_dir(levels_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| level::is_level_file(p))
        .collect();
    paths.sort();
    for p in paths {
        let level = Level::from_json(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let mut h = Hasher::new();
        if level.is_tube() {
            // Tube levels: both runs, every frame's particles and the currents.
            for refine in [1, 2] {
                let run = level.tube_run(&level.reference_solution, refine, |_| {});
                for (t, xs) in &run.frames {
                    h.f64(*t);
                    for x in xs {
                        h.f64(x.x);
                        h.f64(x.y);
                    }
                }
                for (t, q) in &run.collected {
                    h.f64(*t);
                    h.f64(*q);
                }
                h.f64(run.current);
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
        let name = p.file_stem().unwrap().to_string_lossy().into_owned();
        out.insert(name, format!("{:016x}", h.0));
    }
    out
}

#[test]
fn t11_repeated_runs_are_identical() {
    assert_eq!(level_hashes(), level_hashes());
}

#[test]
fn t11_matches_golden_hashes() {
    let hashes = level_hashes();
    let golden_path = levels_dir().join("golden_hashes.json");
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        let mut s = String::from("{\n");
        let body: Vec<String> = hashes
            .iter()
            .map(|(k, v)| format!("  \"{k}\": \"{v}\""))
            .collect();
        s += &body.join(",\n");
        s += "\n}\n";
        std::fs::write(&golden_path, s).unwrap();
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
}
