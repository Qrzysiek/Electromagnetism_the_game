# Handoff (2026-10-28): where the work stands and what comes next

Read this instead of the previous session's transcript. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; decisions: `docs/DECISIONS.md` (DEC-26,
DEC-27 are the latest).

## State

`main` = `worktree-tubes-next` (worktree `.claude/worktrees/tubes-next`). Landed on 2026-10-27/28:

- `1e3fdf7` nested metals keep their carriers (`Tube::metal_of`, `Lattice::holds`): V-D7, the
  transformer test, passes (11.6 %); its old 40° failure was this bug.
- `08f57fd` self-field circuits (DEC-27, approved by the owner): with a level's speed of light,
  each supplied electrode with a border lead is a circuit node; the lead's current is a field
  source and its EMF is in the branch (`crates/physics/src/tube/leads.rs`); every circuit
  closes along the grounded envelope to a ground terminal (`LeadPlan::returns`,
  `circuit_paths`). Refused with the self-field until modelled: leads from rings, circuits
  between electrodes, straps. Tests: leads.rs unit tests, Y27 (`darwin.rs`, `--ignored`).
- `76e2bc8` lead beads in every view (`level::tube_em::LeadFlow`, test
  `lead_flow_keeps_kirchhoffs_law`).
- `b2185e7`, `4597a78` DEC-27's approval; model notes for self-field levels (`level::model`).

## 1. The heavy-test verification (running since 2026-10-28 02:50)

`scripts/refresh_coverage.py` (job `coverage-refresh`, 874 items: 521 verdicts then 353
coverage runs, ~11.4 h, expected end ~14:15) on `4597a78`, in its own window, building in the
worktree's `target/`. Do not build in that worktree while it runs.

- Progress: `python scripts/progress.py status` (or `/progress`). If it shows `lost` or
  stopped early, run `scripts\refresh_coverage.cmd` again: finished tests are kept.
- Results: `python scripts/affected.py status` (failures, stale records). A failing test runs
  alone: `cargo test --release -p level --test per_level -- --ignored <name>` (one name per
  invocation: the harness takes one filter).
- Expectation: every verdict passes. The session's engine changes are bit for bit unchanged for
  levels without a speed of light (none ships one), checked on golden 106, 107, 108, 109, 113.
  A golden-hash mismatch is then a real regression: suspect `1e3fdf7` (which lattice a carrier
  belongs to, `crates/physics/src/tube/metals.rs`) or `08f57fd` (`circuits.rs`: nodes and the
  circuit system with no lead branches must be as before). Never update golden hashes to make
  them pass.
- After it: note the outcome in the next commit's `affected.py note`.

## 2. Ring-shaped metals in the level format (task #34; needed by the Transformer level)

Physics already has rings: `Lattice::with_hole` (walls, no cores, carriers or background inside)
and nested metals since `1e3fdf7`. Steps:

1. Format: `MetalSpec::hole: Option<[Node; 2]>` (`crates/level/src/tube/mod.rs`, serde default,
   skipped when None): the hole's grid corners strictly inside `corners`, walls ≥ 1 cell.
2. Builder (`crates/level/src/tube/build.rs`, the lattices ~line 103, `metal_lattice` ~line 9):
   apply `.with_hole(...)`; the carriers' count from the net area (area − hole), not the
   rectangle's.
3. Editor (`crates/game/src/level_editor/tube.rs` ~line 312, destructured exhaustively, no `..`):
   a hole checkbox and two `node()` pickers; `check.rs` (~line 174): corners on the grid, inside
   the metal by ≥ 1 cell (sandbox completeness: `all_levels_are_reproducible_in_the_sandbox`).
4. Issues (`tube/issues.rs`): the hole inside its metal, clear of its contacts' faces; a metal
   may lie wholly inside another's hole (relax the metals' overlap check for that case only).
5. Geometry and display: everything using `metal_rect` must leave the hole open: the leads'
   router obstacles (`crates/level/src/leads.rs`), `TubeEm::new`'s background deposit
   (`tube_em.rs`), the metal maps and cell accounts, the board's drawing of a metal.
6. Tests: a ring's lattice (no carriers in the hole, the count from the net area); a ring in a
   swinging field circulates (a short run); sandbox reproducibility; add a ring metal to
   `scripts/pairs.py`'s `ELEMENTS`; a capture (`EM_CAPTURE`, `EM_LEVEL_FILE`).
7. Docs: PHYSICS.md (metals), a DECISIONS.md entry, LEVELS.md if a level uses it.

## 3. The first self-field level: Transformer (task #61, after #34)

From V-D7 (`crates/physics/tests/darwin.rs`, `transformer`): two nested closed rings (primary
14 × 14 in a secondary 20 × 20, walls 1, the tests' metal) in the solenoid's swinging field
(`solenoid_drive`), speed of light 3, step ≤ 0.012 (`SELF_FIELD_MAX_STEP`).

- Controls must interact (no one-slider scaling; the owner's rule): e.g. the drive's amplitude
  and a cut in the secondary or the primary's size, the goal (the secondary's current, read as
  `B_z` in the gap, or its heat) needing both; sliders symmetric about 0, a comfortable solving
  interval (≥ ~5 % of the range) avoiding 0, the ends and the halves' middles
  (`generator window`).
- Show: the self-field magnetic map (the primary screening its hole, the secondary's field in
  the gap). The model notes already state the self-field (`4597a78`).
- Run times: V-D7's runs take ~4 min per 9 periods at ω 0.8 with ~7000 carriers; keep the
  level's runs to a few periods, carriers per area modest; check the preview's time.
- Checklist (CLAUDE.md): `scripts/designs/` (its arc) and `curriculum.py`; `scripts/reality.py`
  (headline: a transformer's mutual inductance; simplified: z-invariant one-turn coils, slow
  light); scales; `python scripts/levels.py --only <slug>`; golden hash with
  `GOLDEN_ONLY=<slug>`; reference and radiation tests; the level document
  (`scripts/level_docs.py new`); `docs/LEVELS.md`; a capture.

## Working notes

- Finished work: `git fetch -q origin`, `git merge-base --is-ancestor origin/main HEAD`, then
  `git push -q origin worktree-tubes-next:main` and `git push -q origin worktree-tubes-next`.
  The owner pulls and builds `cargo build --release -p game -p generator`.
- Before committing: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`, the affected heavy tests; `python -m mypy` / `pytest` for scripts;
  `python scripts/measured.py update`, `python scripts/physics_pdf.py` after PHYSICS.md edits.
- The Bash guard refuses compound commands with computed values: write small Python scripts in
  `target/scratch` and run them. Windows file locks (Errno 22): retry.
- Never send clicks or keys to the desktop; captures through the game's own `EM_CAPTURE`.
