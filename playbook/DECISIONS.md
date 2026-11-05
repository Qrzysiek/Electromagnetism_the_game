# Process decisions, for reuse

The decisions behind the rules in [`CLAUDE.template.md`](CLAUDE.template.md), restated for
any repository: each choice, why it was made (with the evidence from this project, where a
rule usually came from a failure), the alternatives turned down, and what to adapt. The
`DEC-n` references point to this repository's own log, `docs/DECISIONS.md`, which has the
full history; the `P-n` numbers are this file's and are cited by the template.

## P-0: The decision log's format

Start the new repository's `docs/DECISIONS.md` with this:

```markdown
# Design decisions

Every design choice: the model, algorithm or design taken, its purpose, what else was
considered and why not, and where it lives. When a choice is changed, its entry is not
rewritten: the change is appended to its history with the reasons.

## DEC-<n>: <short title>

- **Status:** in use | replaced by DEC-<m> | deferred | dropped (since <date>)
- **Area:** <physics | numerics | display | interface | tools | process | …>
- **Choice:** <the model, algorithm or design, precisely>
- **Purpose:** <what it is for; whose request, with the date>
- **Alternatives:** <what else was considered, and why not (measured where possible)>
- **Where:** <files, functions, tests, docs>
- **History:**
  - <date>: chosen. <the reasons, the measurements that decided it>
  - <date>: changed from <a> to <b>. <why: what was found inappropriate, the evidence>
```

Entries are numbered in order and never renumbered.

## P-1: An append-only design-decision log (DEC-8)

- **Choice:** every model, algorithm or design chosen for a feature, a tool or a piece of
  content gets an entry in the same commit; a change appends a dated history line and
  updates the status, and never rewrites the entry.
- **Why:** the owner wanted to see each choice, why it was made, and why it changed when it
  did. Rewritten entries lose the evidence that a choice was tried and failed, and the next
  session tries it again. Here the history lines carry the measurements that overturned a
  choice: the first test-selection rule missed compile-time constants (DEC-16), and the
  first runner skipped hung items, which the owner revised (DEC-23).
- **Alternatives:** commit messages alone (not findable by topic, no status); a wiki
  (drifts from the code; not in the commit that makes the change).
- **Adapt:** start empty with an entry for the log itself. Older choices are added when
  they are next touched, not reconstructed.

## P-2: The owner's rulings cited where they are applied

- **Choice:** a rule that comes from the owner carries the date and, where it matters, the
  owner's words: "(the owner, 2026-11-02: longer is too slow for sensible feedback)". A
  later ruling cites what it supersedes.
- **Why:** an agent cannot tell a deliberate rule from an accident unless the rule says
  where it came from. A dated reason also shows when a rule has been overtaken. For
  example, "every few levels, a full run" was superseded by triggered full runs (DEC-16),
  and the 30 s preview rule was clarified the next day into separate bounds.
- **Adapt:** as is.

## P-3: A handoff document instead of transcripts

- **Choice:** `docs/HANDOFF.md` holds the state (what is on main, what is verified and how),
  the open questions to the owner, and the next steps. Sessions read it first and update it
  last.
- **Why:** sessions end and contexts get compacted. A transcript is too long to reread, and
  a summary of one loses the reasons. The handoff is the one page that must be true.
- **Adapt:** as is. Keep it short; move finished material into the permanent documents.

## P-4: Open problems logged, not lost (the `open-problems` skill)

- **Choice:** a question the literature does not answer, and that cannot be derived in
  reasonable time with the tools at hand, goes into `docs/OPEN_PROBLEMS.md`. Each entry
  gives the precise question, the evidence, what was searched and tried, and the
  workaround in use.
- **Why:** without the log, the workaround looks like the answer, and the search gets
  repeated.
- **Adapt:** the skill is user-level (`~/.claude/skills/open-problems/`), so it already
  works in every repository.

## P-5: Audits archived with their resolution

- **Choice:** each external review goes into `audits/<date>_<fix commit>/` as delivered,
  with a `RESOLUTION.md` giving a verdict and a fix for every remark. `audits/README.md`
  records the audited commit and the fixing commit, so the next audit reviews only
  `git diff <fixed in>..HEAD`.
- **Why:** the remarks would otherwise be fixed once and forgotten, and the next audit would
  start from nothing.
- **Adapt:** as is.

## P-6: Numbers in documents written by the tests that measure them (DEC-22)

- **Choice:** a test records each number a document quotes (`measured(id, value, digits)`
  writes `target/measured/<id>`). The document marks each quote with
  `<!--m:<id>-->0.027<!--/m-->`, which is invisible when rendered. `scripts/measured.py
  update` writes the recorded values in, and `check` lists the quotes that differ.
- **Why:** numbers copied by hand from test output drifted: two table rows were rewritten
  twice as the method changed.
- **Alternatives:** generating whole tables (their setups are prose that no test knows);
  tests asserting the documented numbers (a test would then fail on a documentation lag,
  not on the code).
- **Adapt:** `scripts/measured.py` is generic apart from its document list. The helper is
  about ten lines in each test language.

## P-7: Two test tiers; the heavy tier selected by what each test executed (DEC-16)

- **Choice:**
  - The fast tier (seconds to minutes) runs whole before every commit and in CI.
  - The heavy tier has one test per check and per content item, ignored by default and
    runnable by name.
  - `scripts/affected.py` keeps, for each heavy test, the source files its last run
    executed (LLVM source-based coverage, one process per test) with their hashes, plus the
    non-code inputs (data files, golden values, the lock file, the toolchain). It runs every
    test whose recorded inputs differ, and writes a change note (what changed, what it
    affects, what ran) for the commit message.
  - Full runs are triggered, not scheduled: when a global input changes, on request, or
    when a bug escapes the selection (log what the records missed).
- **Why:** the heavy checks had grown to about 1.5 h and reran in full after any engine
  change. Selection is sound only because the runs are deterministic (P-12): a test whose
  inputs are bit for bit unchanged can only give its recorded result.
- **Lessons learned the hard way:**
  - Constants are folded in at compile time, so "the files it executed" misses them. A
    change outside function bodies (a declaration, a `const fn`) therefore selects every
    test; a change inside a body selects only the tests that executed its file.
  - Coverage is per file, so a 4000-line file makes every test depend on it. The largest
    files were split by concern (P-11).
  - A verdict-only run leaves every record's coverage stale, and a stale record is selected
    by any change. For a day every change ran all 525 tests. Always run with coverage, or
    refresh right after.
  - A test that computes exactly what another does (`repeated::` = `golden::` twice) takes
    the other's coverage instead of an instrumented run of its own (a metal tube level's
    took 34 minutes).
  - The records are shared by every worktree, in the main checkout. Per-worktree records
    made each new worktree select everything again.
  - A hang limit based on an estimate needs an estimate of the same work. A plain run that
    took a recorded result (0.2 s) set a 600 s limit on an instrumented run that computes
    for 30 minutes, and killed it (2026-11-03).
  - The instrumentation's slowdown has to be measured, not assumed, and it is a property of
    the workload: 4–32× here, alike within an item, not by its kind. An assumed 10× left
    the longest tests at two-thirds of their limit. Use the item's own measured ratio, and
    the highest seen for a new item.
- **Alternatives:** static selection from the module graph (trait impls, macros and generics
  make it unsafe to under-approximate in Rust); a hand-written map from features to modules
  (goes stale silently); selection per crate (too coarse); learned selection (needs a
  history); Bazel (a migration out of proportion).
- **Adapt:** `affected.py` is tied to Rust, LLVM coverage (`scripts/coverage_rustc.py` as
  `RUSTC_WORKSPACE_WRAPPER`) and a data-driven heavy test program (`libtest-mimic`, one
  trial per check and item). For another language keep the idea (Ekstazi, Gligoric et al.,
  ISSTA 2015; Microsoft's Test Impact Analysis), the declaration/body rule and the lessons
  above.

## P-8: Long computations estimated, tracked, chunked and resumable (DEC-23, `docs/PROGRESS.md`)

- **Choice:**
  - Estimate per item and in total before launching (within 50 % is fine), and say it.
  - Run through `scripts/progress.py run` with a regex counting finished items. `/progress`
    shows done/total, time per item against the estimate, and the ETA.
  - Anything of an hour or more is split into chunks of at most an hour, each saving its
    results, chained to run by themselves (`progress.py chain`).
  - Runs of many independent items go through `scripts/runner.py`. Each item runs in its
    own process, N at a time. A ledger keyed by each item's command and inputs records
    every outcome as it comes. `ctl <job> pause|resume|stop|skip ITEM|jobs N` works while
    it runs.
  - An item past four times its estimate is hung, and the job stops there. After the fix
    the same command continues: the stopping item runs first, and only finished items whose
    inputs changed run again.
  - Items are timed by the machine's awake time, not the wall clock: a night's sleep once
    counted as 6.7 h of a test's run, killed it as hung and interrupted its 15 neighbours.
- **Why:** before this, a run with one hung or failing item had to be killed and redone
  whole, twice in one session. "Hangs on one item for several hours" also went unseen,
  because progress was counted per test binary.
- **Alternatives:** nextest alone (per-test timeouts, but no ledger across runs and tests
  only); GNU parallel (not on Windows, no fingerprints); skipping a hung item and going on
  (the owner: the next items may hang too).
- **Adapt:** `runner.py` is generic. `progress.py` is generic except its `estimate`
  subcommand (this game's cost model): drop it or write your own. Write `docs/PROGRESS.md`'s
  table of item patterns and timings anew. Keep an unattended run's machine awake
  (`SetThreadExecutionState` on Windows) and release it at exit.

## P-9: Scripts type-checked and tested like the code (DEC-21)

- **Choice:** `python -m mypy` (mypy's defaults, which check the bodies of annotated
  functions) and `python -m pytest` on every script, locally and in CI. The tests cover what
  could go wrong silently, above all the tools the process trusts (the selection's parser).
- **Why:** the scripts had grown to 9300 lines without checks, and test selection depends on
  one of them.
- **Alternatives:** strict mypy (thousands of annotations on data-like code); pyright (an
  equivalent tool, one more to install).
- **Adapt:** `mypy.ini` and `pytest.ini` as is.

## P-10: A pinned toolchain and CI hygiene (DEC-17)

- **Choice:** `rust-toolchain.toml` pins the compiler and its components; an upgrade is a
  deliberate bump followed by every test. CI installs that toolchain, and its hygiene job
  runs cargo-deny (licences, advisories, sources), cargo-machete (unused dependencies) and
  the documentation with warnings as errors.
- **Why:** golden hashes, exact caches and test selection all assume the compiler does not
  change under them.
- **Adapt:** for other ecosystems, a lock file plus a pinned interpreter or compiler, and
  their audit tools.

## P-11: Structure for safety, and refactors as verified pure moves (DEC-18, DEC-19, DEC-20)

- **Choice:**
  - A quantity is computed in one place and used by every consumer. A doubled background
    came from fields assembled in two places.
  - A constructor leaves every optional part empty, and literals name only what they set
    (`..Tube::new(...)`). One new field had needed 18 literals edited.
  - A named struct replaces several alike positional numbers: a constructor with nine
    positional arguments, four of them alike floats, had invited swaps.
  - Test helpers are shared, not copied.
  - The largest files are split into modules by concern as pure moves: visibility kept
    exactly, paths re-exported, made by a mover that lexes the language, and verified bit
    for bit (golden outputs, every test).
- **Why:** each rule answers a bug class this project actually hit. Splitting by concern
  also makes coverage-based selection finer (P-7).
- **Alternatives:** full builder types (more code for the same safety); new crates for the
  splits (a crate boundary changes privacy, which a pure move must not).
- **Adapt:** as principles.

## P-12: A deterministic authoritative path

- **Choice:** the computation whose results count (here the trajectories) is f64 on the CPU
  and bit-reproducible:
  - no fused multiply-add;
  - a fixed summation order, in parallel reductions too;
  - elementary functions (`sin`, `exp`, `atan2`, …) from a portable library (the `libm`
    crate), never the platform's C library. `clippy.toml` disallows the std methods.

  The GPU is for visuals only.
- **Why:** the golden hashes differed between Windows and Linux CI until the elementary
  functions moved to `libm`: `atan2` in an acceptance margin. Determinism is what makes
  golden hashes, the hash-keyed run cache (keyed by every input including the executable)
  and test selection exact rather than approximate.
- **Adapt:** for a numerical project, take it early; retrofitting changed every hash once.

## P-13: Validation discipline (the owner, 2026-10-20)

- **Choice:**
  - Validation tests compare against analytic or independent results, such as a CAS at
    high precision or an independent implementation.
  - Thresholds are set from the target accuracy before measuring. A threshold change is
    documented with its reason, and never loosened to make a test pass.
  - Symmetry and conservation tests are mandatory for every algorithm that implements the
    model:
    - the system's symmetries;
    - those of the equations of motion and their conservation laws;
    - the background's (translations, rotations, boosts);
    - the discrete ones (P, T, C, their products) and the model's scalings.

    Each is tested against the exact transformed solution, and the model document records
    how well and over what time each holds.
  - Every number in an explanatory document is computed by its own script, independently of
    the program, and compared with the program's output. A failed comparison is fixed, not
    widened.
- **Why:** a numerical method can silently break a symmetry, and a threshold fitted to the
  measurement proves nothing.
- **Adapt:** name the model document (here `PHYSICS.md`) and its symmetry section.

## P-14: Disk hygiene (DEC-37)

- **Choice:**
  - Builds go only to the checkout's standard build directories; temporary ones go under
    `target/temp/`.
  - A session's files go to `target/scratch/`; whatever a tracked file cites from there
    moves to `prototypes/` in the same commit.
  - `scripts/cleanup.py` (`status` read-only; `run`) deletes:
    - temporary directories and anything under `target/` its `KNOWN` does not name;
    - builds over their caps (about twice their fresh size; cargo never deletes stale
      builds);
    - cache records older than the program that wrote them;
    - uncited scratch untouched for 14 days;
    - removed worktrees' leftovers.
  - It runs after each merge and whenever a build directory passes 20 GB.
- **Why:** the project reached 84 GB. 66 GB were extra build directories that sessions had
  made and never removed. The rest was stale builds inside the profiles: 26 versions of one
  library. 29 citations in tracked files pointed into one worktree's unversioned scratch.
- **Alternatives:** cargo's unstable gc; cargo-sweep (access times, unreliable on
  Windows); a check at every session start (slow, and it could delete what a running session
  uses).
- **Adapt:** `cleanup.py`'s `KNOWN`, `CAPS`, the rebuild after a cap and the cache rule are
  this project's; the structure is generic.

## P-15: Worktrees, fast-forward merges, and finishing steps

- **Choice:**
  - Work happens on a branch in a worktree.
  - Finished work is merged to main fast-forward only.
  - What the owner runs from the main checkout is rebuilt if code changed.
  - Main and the branch are pushed, then the cleanup runs.
  - The stash is shared by every worktree and session: stash only with a unique tag and
    apply by SHA, or prefer a WIP commit.
  - Commit messages end with the change note (P-7) and the co-author line.
- **Why:** the owner plays from the main checkout. A fast-forward keeps main exactly the
  tested commit.
- **Adapt:** as is.

## P-16: Visual checks without touching the desktop (DEC-29)

- **Choice:** the program has a capture mode. Environment variables open a given state and
  view, it saves a screenshot of its own window at a given frame (or when the results are
  in), and it exits. The agent never sends clicks or keys to the desktop. README
  screenshots are retaken by one script (WebP, a quarter of PNG's size).
- **Why:** the owner uses the machine while the agent works, and a screenshot taken by hand
  is not reproducible. The capture mode also carries regression checks: switching views
  at given frames, minimizing and restoring, abandoning stale computations, and a log of
  the frame rate and the time shown each second (the animation must play on while results
  are computed: the game's had stood still through a verification, DEC-39).
- **Adapt:** design the capture mode in from the start; add a variable for each state a
  check needs.

## P-17: Interactive content: complete, explained, honest, non-trivial

For an interactive simulation or educational game (this project's levels):

- **Reproducible in the sandbox:** the editor edits every field of the content format. Each
  type is destructured exhaustively, with no `..`, so a new field fails to compile until it
  gets a control. The editor's ranges are shared constants that a test validates on every
  shipped file. A new feature needs its control in the same change.
- **Everything explained, without clutter:** every mark drawn that is not hoverable has a
  legend row. Panel text is one short line with the detail on hover. Controls carry hover
  text.
- **Says what it represents** (the owner, 2026-10-21): the real thing shown, what carries
  over, what is simplified or left out. Also its scales: the real-world magnitudes, and the
  dimensionless ratios that decide the behaviour. Each item has a test.
- **No trivial solutions** (DEC-12, DEC-30):
  - nothing solves with the controls untouched or at their first touch;
  - slider ranges hide the answer: no solving interval at 0, an end, or a half's middle,
    and on a logarithmic slider none within 5 % of those;
  - one slider scaling a response in proportion is an introduction at most.
- **Showcase everything feasible:** every control spans all the behaviours the model can
  produce within reason, and the views make them visible. The answer is picked from among
  them, with a comfortable solving interval (about 5 % of a slider).
- **Curriculum:** every new element or concept gets an easy introduction first. Within a
  chapter, difficulty and the number of elements needed rise on average.
- **Timing** (the owner, 2026-11-02/03):
  - the preliminary answer within 30 s;
  - the animation streams without stutter at speed 1, however long it is;
  - the verification within a few minutes;
  - each bound may be exceeded by at most 100 % where that makes an item a lot better.

  What does not fit accurately becomes another format (here an "interlude").
- **One source of truth for content:** a generator script builds the content files, which
  are never edited by hand. Every shipped item has a verified reference solution, checked
  in the heavy tier.
- **Adapt:** the rules are about the player's understanding, not about physics. The
  thresholds (30 s, 5 %) are this owner's.
