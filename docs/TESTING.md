# Testing

Two tiers (DEC-16): a fast tier that always runs whole, and a heavy tier of one test per
check and level, of which a change runs only the tests it can affect.

## The fast tier

`cargo test --workspace` (or `cargo nextest run --workspace`, which shows each test as it
finishes): every unit and integration test but the per-level checks. 346 tests, 3.6 min
under nextest on the owner's 16-thread machine (2026-10-27); the slowest are the tube
engine's tests of `level/tests/tubes.rs`, 2 to 3.5 min each, run side by side. Before every
commit, with `cargo fmt --all` and `cargo clippy --workspace --all-targets -- -D warnings`;
CI runs it on Windows and Linux.

## The heavy tier: one test per check and level

`crates/level/tests/per_level.rs` makes a test of each check of each shipped level, named
`<check>::<level file's stem>` (521 tests), ignored by default:

| check | what it asserts |
|---|---|
| `golden::` | T11 (PHYSICS.md §9): the reference flown at both tolerances, every bit hashed, equal to `levels/golden_hashes.json` |
| `repeated::` | the same hash twice in one process |
| `reference::` | every shot's reference verified, the level's step and memory budgets kept (tube levels: the goal met) |
| `radiation::` | the radiation the model neglects below 1e-10 of the launch energy |
| `needs_radiation::` | (levels with radiation reaction) the reaction is needed |

Run any by name: a filter is a part of the name.

    cargo test --release -p level --test per_level -- --ignored reference::109_
    cargo test --release -p level --test per_level -- --ignored 122_inverter
    cargo nextest run --release -p level --test per_level --run-ignored only -E 'test(/^golden::/)'

Four run at a time (they share the process's rayon threads; nextest gives each four
threads, `.config/nextest.toml`). A check that passed is skipped while its level's physics
(`Level::physics_json`) and the test program are unchanged (`level::run_cache::TestPasses`,
`target/test-passes`), so an interrupted run resumes; `EM_TEST_FULL=1` runs it anyway. The
golden check compares even a recorded hash with the golden file.

Golden hashes change only when the physics or a level changes intentionally:
`UPDATE_GOLDEN=1 cargo test --release -p level --test per_level` writes every level's
(`GOLDEN_ONLY=<part of a name>`: only the matching levels', keeping the rest).

## The heavy tests a change affects

    python scripts/affected.py status     # what the working tree selects, and why
    python scripts/affected.py run        # their verdicts (plain program), then their coverage
    python scripts/affected.py run --no-coverage   # verdicts only; coverage left for later
    python scripts/affected.py refresh    # the coverage left for later
    python scripts/affected.py note       # the change note for the commit message

`run` gives the verdicts first, with the plain test program at full speed (four tests at a
time on four threads each, as the program runs them itself), then takes the same tests'
coverage with the instrumented one. A test whose coverage is older than its verdict (`run
--no-coverage`, or a coverage run cut short) counts as depending on every engine file until
`refresh` brings it up to date: a change of code may have sent it into files it did not
execute before.

Each heavy test's record holds what its last run depended on:

- the workspace files it executed: LLVM source-based coverage of an instrumented build
  (`-C instrument-coverage`, in `target/affected/build`), one process per test;
- the state of every engine file at the time (`crates/physics/src`, `crates/level/src`
  without their unit tests, and `per_level.rs`), each as two hashes, both without comments
  or layout: the file's code, and its declarations (the code without its functions' bodies);
- the build's global inputs: `Cargo.lock`, the `Cargo.toml` files of the workspace, `physics`
  and `level`, `rust-toolchain.toml`, the compiler's version;
- its level's physics (`Level::physics_json`, hashed, as the program computes it); for
  `golden::`, the golden entry.

A test is selected again when:

| change | selects |
|---|---|
| a function's body | the tests that executed its file |
| declarations: constants, statics, types, signatures, macros, attributes, `const fn` bodies; a new or removed engine file | every test (they can act without being executed, where coverage cannot see) |
| a global input | every test |
| a level's physics (not its texts) | that level's tests |
| a golden entry | that level's `golden::` |
| a test that failed, or was never recorded | that test |
| any function's body, for a test whose coverage is older than its verdict | that test |

Comments and formatting select nothing, nor do changes outside the engine (the game, the
generator, the documents); `scripts/levels.py` reaches the heavy tests through the level
files it writes. Measured on two levels (a charges level and a vacuum tube, 2026-10-27):
a body in the tube engine selected the tube's three checks that ran it, not its radiation
check (no flights) nor the charges level's; a body in the flight integrator the charges
level's four; a constant all eight; a comment, a reworded description and nothing at all,
none; a level's launch energy its four; a golden entry its one.

Why it is sound: the runs are deterministic (PHYSICS.md §9), so a test whose inputs are bit
for bit those of a passing run passes again. A test's inputs are the code it runs, the code
that acts without running (declarations, compile-time evaluation), the build, and the files
it reads; the record covers each. What it relies on: the heavy checks read no file but their
level's and the golden file (a check that starts reading another must add it to the record),
and the engine reads no environment variable that changes a result.

Costs: only the workspace's crates are instrumented (`scripts/coverage_rustc.py`, cargo's
`RUSTC_WORKSPACE_WRAPPER`; the dependencies are `Cargo.lock`'s). The instrumented build runs
a tube level about 2.5 times slower than the plain one on one thread (14 s against 5.8 s;
with every crate instrumented, glam's and libm's inner loops included, 38 s), and slower
still on several, the counters being memory the threads share (23 s at 4 threads; the plain
build 1.6 s at 16). The coverage runs therefore give each test one thread, as many tests
at once as the machine has threads, the longest first (by their last coverage run, else ten
times their plain run, else by their level's kind): about 2.5 times the plain tier's CPU time,
spread over every core. Grouped by check instead (a first plan), each group ended with a few
metal tube levels on a few cores while the rest waited. A fresh checkout records everything
with `python scripts/affected.py run --coverage-only` (verdict and coverage in one
instrumented run); each test's record is saved as it finishes, so an interruption loses only
the tests in flight.

## Full runs

A full run (`python scripts/affected.py run --all`, or every heavy test by
`cargo test --release -p level --test per_level -- --ignored`) is due when a global input
changes (every record then differs anyway), on request, and when a bug escapes the selected
tests: then log what the records missed (DEC-16's history) before fixing it. There is no
fixed cadence.

## CI

`.github/workflows/ci.yml`: format, clippy and the fast tier on Windows and Linux; the
golden hashes of one level per feature on Linux (T11 across platforms: the hashes are
committed from Windows); the dependencies' licences and advisories, unused dependencies and
the documentation's links (DEC-17).
