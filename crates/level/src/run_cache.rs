//! A hash-locked cache of tube runs for the offline tools (the owner, 2026-10-20: "a must").
//! Each run's record is stored under the hash of everything that could change it: the level
//! (all of it, as JSON), the placement, the resolution and the run's end, and the engine, a
//! fingerprint the tool sets (`generator`: its own executable's hash, so that any change of
//! the physics, the level code, a library or the compiler gives another key). The runs are
//! deterministic (the golden hashes), so a record is bit for bit the run's; its floats are
//! stored as their bits (infinities and NaN too). With `verify` every hit is run again and
//! compared bit for bit, and a difference stops the tool. Off unless a tool enables it: the
//! game and the tests run as before.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::{Value, json};

use crate::Element;
use crate::tube::{TubeResolution, TubeResults};

/// The records' format: part of every key, so that a change of what is stored is a change
/// of key.
const SCHEMA: &str = "em-run-cache-1";

/// One tube run to the verification's end: what the verdict, the search objective and the
/// cost meters take from it.
#[derive(Clone, Debug, PartialEq)]
pub struct TubeRecord {
    pub results: TubeResults,
    /// Each projectile's closest approach to its detector.
    pub closest: Vec<f64>,
    pub segments: usize,
    pub steps: u32,
    /// How long the run took when it was computed (the cost meters'; not compared).
    pub seconds: f64,
}

fn bits(v: f64) -> Value {
    json!(v.to_bits())
}

fn opt(v: Option<f64>) -> Value {
    v.map_or(Value::Null, bits)
}

fn float(v: &Value) -> Option<f64> {
    v.as_u64().map(f64::from_bits)
}

fn floats(v: &Value) -> Option<Vec<f64>> {
    v.as_array()?.iter().map(float).collect()
}

fn options(v: &Value) -> Option<Vec<Option<f64>>> {
    v.as_array()?
        .iter()
        .map(|x| {
            if x.is_null() {
                Some(None)
            } else {
                float(x).map(Some)
            }
        })
        .collect()
}

impl TubeRecord {
    fn to_json(&self) -> Value {
        let r = &self.results;
        json!({
            "currents": r.currents.iter().copied().map(opt).collect::<Vec<_>>(),
            "potential": r.potential.iter().copied().map(opt).collect::<Vec<_>>(),
            "heat": r.heat.iter().copied().map(opt).collect::<Vec<_>>(),
            "arrivals": r.arrivals.iter().copied().map(opt).collect::<Vec<_>>(),
            "capped": r.capped,
            "closest": self.closest.iter().copied().map(bits).collect::<Vec<_>>(),
            "segments": self.segments,
            "steps": self.steps,
            "seconds": bits(self.seconds),
        })
    }

    fn from_json(v: &Value) -> Option<Self> {
        Some(Self {
            results: TubeResults {
                currents: options(&v["currents"])?,
                potential: options(&v["potential"])?,
                heat: options(&v["heat"])?,
                arrivals: options(&v["arrivals"])?,
                capped: v["capped"].as_bool()?,
            },
            closest: floats(&v["closest"])?,
            segments: usize::try_from(v["segments"].as_u64()?).ok()?,
            steps: u32::try_from(v["steps"].as_u64()?).ok()?,
            seconds: float(&v["seconds"])?,
        })
    }

    /// Whether two records are the same run bit for bit (all but the time it took).
    fn same_run(&self, other: &Self) -> bool {
        let mut a = self.to_json();
        let mut b = other.to_json();
        a["seconds"] = Value::Null;
        b["seconds"] = Value::Null;
        a == b
    }
}

/// The cache: its directory, the engine's fingerprint and whether hits are run again.
pub struct RunCache {
    dir: PathBuf,
    engine: String,
    verify: bool,
    hits: AtomicU64,
    misses: AtomicU64,
}

impl RunCache {
    #[must_use]
    pub fn new(dir: PathBuf, engine: String, verify: bool) -> Self {
        Self {
            dir,
            engine,
            verify,
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
        }
    }

    /// The key of a run: the hash of the schema, the engine, the kind of run, the level's
    /// JSON, the placement and the run's parameters (as bits).
    #[must_use]
    pub fn key(
        &self,
        kind: &str,
        level_json: &str,
        placement: &[Element],
        parts: &[u64],
    ) -> String {
        let mut h = blake3::Hasher::new();
        for s in [SCHEMA, &self.engine, kind, level_json] {
            h.update(&(s.len() as u64).to_le_bytes());
            h.update(s.as_bytes());
        }
        let placement = serde_json::to_string(placement).expect("a placement serializes");
        h.update(&(placement.len() as u64).to_le_bytes());
        h.update(placement.as_bytes());
        for p in parts {
            h.update(&p.to_le_bytes());
        }
        h.finalize().to_hex().to_string()
    }

    /// The key of a tube run of `level` with `player` at `res` to `end`.
    #[must_use]
    pub fn tube_key(
        &self,
        kind: &str,
        level: &crate::Level,
        player: &[Element],
        res: TubeResolution,
        end: f64,
        extra: &[f64],
    ) -> String {
        let mut parts = vec![
            u64::from(res.refine),
            res.particles.to_bits(),
            res.softening.to_bits(),
            res.emitter.to_bits(),
            end.to_bits(),
        ];
        parts.extend(extra.iter().map(|x| x.to_bits()));
        self.key(kind, &level.physics_json(), player, &parts)
    }

    fn path(&self, key: &str) -> PathBuf {
        self.dir.join(format!("{key}.json"))
    }

    fn read(&self, key: &str) -> Option<Value> {
        let text = std::fs::read_to_string(self.path(key)).ok()?;
        serde_json::from_str(&text).ok()
    }

    /// Writes atomically: a temporary file renamed over the record (parallel tools share the
    /// directory; a record lost to a race is computed again, never read half written).
    fn write(&self, key: &str, value: &Value) {
        if std::fs::create_dir_all(&self.dir).is_err() {
            return;
        }
        let tmp = self.dir.join(format!("{key}.{}.tmp", std::process::id()));
        if std::fs::write(&tmp, value.to_string()).is_ok()
            && std::fs::rename(&tmp, self.path(key)).is_err()
        {
            let _ = std::fs::remove_file(&tmp);
        }
    }

    /// The record of the run `key`, from the cache or computed (and stored).
    ///
    /// # Panics
    /// With `verify`, if a stored record differs from the run computed again.
    pub fn tube_record(&self, key: &str, compute: impl FnOnce() -> TubeRecord) -> TubeRecord {
        if let Some(stored) = self.read(key).as_ref().and_then(TubeRecord::from_json) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            if self.verify {
                let fresh = compute();
                assert!(
                    stored.same_run(&fresh),
                    "run cache: the record {key} differs from the run computed again:\n\
                     stored {stored:?}\nfresh  {fresh:?}"
                );
            }
            return stored;
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        let fresh = compute();
        self.write(key, &fresh.to_json());
        fresh
    }

    /// A text result (a tool's printed line) of the run `key`, from the cache or computed.
    ///
    /// # Panics
    /// With `verify`, if a stored text differs from the one computed again.
    pub fn text(&self, key: &str, compute: impl FnOnce() -> String) -> String {
        if let Some(stored) = self
            .read(key)
            .and_then(|v| v["text"].as_str().map(String::from))
        {
            self.hits.fetch_add(1, Ordering::Relaxed);
            if self.verify {
                let fresh = compute();
                assert!(
                    stored == fresh,
                    "run cache: the text {key} differs from the run computed again:\n\
                     stored {stored}\nfresh  {fresh}"
                );
            }
            return stored;
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        let fresh = compute();
        self.write(key, &json!({ "text": fresh }));
        fresh
    }

    /// Hits and misses so far.
    #[must_use]
    pub fn stats(&self) -> (u64, u64) {
        (
            self.hits.load(Ordering::Relaxed),
            self.misses.load(Ordering::Relaxed),
        )
    }

    #[must_use]
    pub fn dir(&self) -> &Path {
        &self.dir
    }
}

/// The running executable's fingerprint (the hash of its bytes): the engine part of a
/// tool's keys, so that a rebuilt tool computes again (a rebuild that changes nothing
/// physical costs reruns, never a stale record).
///
/// # Errors
/// If the executable cannot be read.
pub fn exe_fingerprint() -> std::io::Result<String> {
    let bytes = std::fs::read(std::env::current_exe()?)?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}

/// The heavy per-level tests' record of passes (the owner, 2026-10-20: run only the tests a
/// change can affect). A level's check that passed is not run again while the test program
/// and the level's file are byte-identical: the runs are deterministic, so it could only pass
/// again. A change to the physics or the level code rebuilds the test program, and every
/// level is checked again; a new or changed level file is checked alone. `EM_TEST_FULL=1`
/// checks everything. Kept in `target/test-passes` (the workspace's, whatever the build
/// directory).
pub struct TestPasses {
    dir: PathBuf,
    program: String,
    test: String,
}

impl TestPasses {
    /// The record of the test `test` in this test program; None with `EM_TEST_FULL=1` or if
    /// the program cannot be read (every level is then checked).
    #[must_use]
    pub fn open(test: &str) -> Option<Self> {
        if std::env::var("EM_TEST_FULL").is_ok_and(|v| v == "1") {
            return None;
        }
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/test-passes");
        Self::at(dir, exe_fingerprint().ok()?, test)
    }

    /// The record kept in `dir` for `program` (a fingerprint) and `test`.
    #[must_use]
    pub fn at(dir: PathBuf, program: String, test: &str) -> Option<Self> {
        Some(Self {
            dir,
            program,
            test: test.to_owned(),
        })
    }

    fn path(&self, level_file: &[u8]) -> PathBuf {
        // The level's physics (`Level::physics_json`): rewording it checks nothing again.
        let physics = std::str::from_utf8(level_file)
            .ok()
            .and_then(|s| crate::Level::from_json(s).ok())
            .map(|l| l.physics_json());
        let level_file = physics.as_ref().map_or(level_file, |p| p.as_bytes());
        let mut h = blake3::Hasher::new();
        for part in [self.program.as_bytes(), self.test.as_bytes(), level_file] {
            h.update(&(part.len() as u64).to_le_bytes());
            h.update(part);
        }
        self.dir.join(format!("{}.pass", h.finalize().to_hex()))
    }

    /// What the check of the level whose file is `level_file` recorded when it passed under
    /// this program (e.g. its hash), or None if it has not.
    #[must_use]
    pub fn passed(&self, level_file: &[u8]) -> Option<String> {
        std::fs::read_to_string(self.path(level_file)).ok()
    }

    /// Records that the check of the level `level_file` passed, with `value`.
    pub fn record(&self, level_file: &[u8], value: &str) {
        if std::fs::create_dir_all(&self.dir).is_err() {
            return;
        }
        let path = self.path(level_file);
        let tmp = path.with_extension(format!("{}.tmp", std::process::id()));
        if std::fs::write(&tmp, value).is_ok() && std::fs::rename(&tmp, &path).is_err() {
            let _ = std::fs::remove_file(&tmp);
        }
    }
}

static CACHE: OnceLock<RunCache> = OnceLock::new();

/// Turns the cache on for this process (once); false if it already was.
pub fn enable(cache: RunCache) -> bool {
    CACHE.set(cache).is_ok()
}

/// The process's cache, if a tool enabled it.
#[must_use]
pub fn get() -> Option<&'static RunCache> {
    CACHE.get()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!("em-run-cache-test-{}-{name}", std::process::id()))
    }

    fn record() -> TubeRecord {
        TubeRecord {
            results: TubeResults {
                currents: vec![Some(2.514_914_123_456_789), None],
                potential: vec![Some(-0.0), Some(f64::NAN)],
                heat: vec![Some(f64::MIN_POSITIVE)],
                arrivals: vec![None, Some(12.5)],
                capped: true,
            },
            closest: vec![f64::INFINITY, 0.1 + 0.2],
            segments: 812,
            steps: 4571,
            seconds: 78.6,
        }
    }

    /// Every float, -0, NaN and the infinities included, comes back bit for bit.
    #[test]
    fn records_round_trip_bit_for_bit() {
        let r = record();
        let back = TubeRecord::from_json(&r.to_json()).expect("parses");
        assert!(r.same_run(&back));
        assert_eq!(
            back.results.potential[0].map(f64::to_bits),
            Some((-0.0_f64).to_bits())
        );
        assert!(back.results.potential[1].is_some_and(f64::is_nan));
        assert_eq!(back.closest[1].to_bits(), (0.1_f64 + 0.2).to_bits());
        assert_eq!(back.seconds.to_bits(), 78.6_f64.to_bits());
    }

    /// The key changes with each of its inputs and only with them.
    #[test]
    fn keys_follow_every_input() {
        let a = RunCache::new(scratch("keys"), "engine-a".into(), false);
        let b = RunCache::new(scratch("keys"), "engine-b".into(), false);
        let e = |v: f64| -> Element {
            serde_json::from_value(json!({ "node": [21, 15, 0], "kind": "supply", "value": v }))
                .expect("an element")
        };
        let base = a.key("tube", "{\"level\": 1}", &[e(30.0)], &[1, 2]);
        assert_eq!(base, a.key("tube", "{\"level\": 1}", &[e(30.0)], &[1, 2]));
        for other in [
            b.key("tube", "{\"level\": 1}", &[e(30.0)], &[1, 2]),
            a.key("scan", "{\"level\": 1}", &[e(30.0)], &[1, 2]),
            a.key("tube", "{\"level\": 2}", &[e(30.0)], &[1, 2]),
            a.key(
                "tube",
                "{\"level\": 1}",
                &[e(30.000_000_000_000_004)],
                &[1, 2],
            ),
            a.key("tube", "{\"level\": 1}", &[e(30.0)], &[1, 3]),
            a.key("tube", "{\"level\": 1}", &[e(30.0)], &[2, 1]),
        ] {
            assert_ne!(base, other);
        }
    }

    /// A pass is found again only for the same program, test and level file.
    #[test]
    fn passes_follow_program_test_and_level() {
        let dir = scratch("passes");
        let _ = std::fs::remove_dir_all(&dir);
        let a = TestPasses::at(dir.clone(), "program-a".into(), "golden").expect("a record");
        assert_eq!(a.passed(b"level one"), None);
        a.record(b"level one", "0123abcd");
        assert_eq!(a.passed(b"level one").as_deref(), Some("0123abcd"));
        assert_eq!(a.passed(b"level one, changed"), None);
        let b = TestPasses::at(dir.clone(), "program-b".into(), "golden").expect("a record");
        assert_eq!(b.passed(b"level one"), None);
        let c = TestPasses::at(dir.clone(), "program-a".into(), "repeated").expect("a record");
        assert_eq!(c.passed(b"level one"), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A level's texts (name, description, reality) key neither a pass nor a run; its
    /// physics does.
    #[test]
    fn rewording_a_level_keeps_its_keys() {
        let dir = scratch("reworded");
        let _ = std::fs::remove_dir_all(&dir);
        let passes = TestPasses::at(dir.clone(), "program".into(), "golden").expect("a record");
        let cache = RunCache::new(dir.clone(), "engine".into(), false);
        let level = crate::shipped("gate");
        let mut reworded = level.clone();
        reworded.name = "Another name".into();
        reworded.description = "Another description.".into();
        reworded.reality.represents = "something else".into();
        reworded.reality.remarks.clear();
        reworded.electrodes[0].label = "another function".into();
        if let Some(t) = &mut reworded.tube {
            t.metals[0].label = "another function".into();
        }
        let mut changed = level.clone();
        changed.physics.t_max += 1.0;
        let res = TubeResolution::refined(1);
        let key = |l: &crate::Level| cache.tube_key("tube", l, &[], res, 10.0, &[]);
        passes.record(level.to_json().as_bytes(), "pass");
        assert_eq!(
            passes.passed(reworded.to_json().as_bytes()).as_deref(),
            Some("pass")
        );
        assert_eq!(key(&reworded), key(&level));
        assert_eq!(passes.passed(changed.to_json().as_bytes()), None);
        assert_ne!(key(&changed), key(&level));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A stored record is returned without running; a run again with `verify` that differs
    /// from it stops the tool.
    #[test]
    fn hits_skip_the_run_and_verify_catches_a_difference() {
        let dir = scratch("hits");
        let _ = std::fs::remove_dir_all(&dir);
        let cache = RunCache::new(dir.clone(), "engine".into(), false);
        let first = cache.tube_record("k", record);
        let second = cache.tube_record("k", || panic!("a hit must not run"));
        assert!(first.same_run(&second));
        assert_eq!(cache.stats(), (1, 1));
        let checking = RunCache::new(dir.clone(), "engine".into(), true);
        let same = checking.tube_record("k", record);
        assert!(same.same_run(&first));
        let differs = std::panic::catch_unwind(|| {
            checking.tube_record("k", || {
                let mut r = record();
                r.steps += 1;
                r
            })
        });
        assert!(differs.is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
