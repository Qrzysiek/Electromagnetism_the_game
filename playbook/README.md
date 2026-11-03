# Playbook: the working rules and process decisions, for another repository

This project's way of working with an agent (Claude Code) on a long, accuracy-critical
project, extracted so another repository can start with it:

- [`CLAUDE.template.md`](CLAUDE.template.md): the rules the agent follows, in layers. A is
  for any repository, B for numerical or scientific code, C for a Rust workspace, D for an
  interactive simulation or educational game, and E for Claude Code on Windows.
- [`DECISIONS.md`](DECISIONS.md): the process and tooling decisions behind those rules
  (P-0 to P-17). Each gives the choice, why it was made, what was turned down, and what to
  adapt. Most rules came from a failure here, and the entry says which.

The game's own decisions (its physics, levels and display) stay in `docs/DECISIONS.md`.
Only what carries over to other projects is restated here. When a rule or process decision
that would carry over is added or changed in `CLAUDE.md` or `docs/DECISIONS.md`, mirror it
here in the same commit (DEC-38).

## Applying it to another repository

1. Copy `CLAUDE.template.md` to the new repository as `CLAUDE.md`. Keep the layers that
   apply, fill the `<…>` placeholders, and delete the rest along with the opening comment.
2. Start `docs/DECISIONS.md` from P-0's format. Copy the P-entries you keep as its first
   entries, each dated with "adopted from <this repository>, P-n". The `P-n` citations in
   the rules then become `DEC-n`.
3. Create the documents the rules name: `docs/HANDOFF.md`, `docs/OPEN_PROBLEMS.md`,
   `docs/PROGRESS.md` (its table of item patterns and timings, written anew),
   `audits/README.md`, and `prototypes/README.md`.
4. Copy the tools the kept rules use, adapting them as the table says.
5. Copy `.claude/skills/progress/` (the `/progress` report, read-only). The `open-problems`
   and `wolfram-engine` skills are user-level (`~/.claude/skills/`) and already work in every
   repository.
6. Put machine facts (interpreter versions, installed packages) in the agent's memory, not
   in the rules.

## The tools

| File | What it does | What to adapt |
|---|---|---|
| `scripts/runner.py` (+ `scripts/tests/test_runner.py`) | runs many independent items, N at a time, with a ledger keyed by their inputs; pause, resume, skip and stop while running; a hang stops the job (P-8) | nothing |
| `scripts/progress.py` | wraps a long run with a progress count, ETA and chunk chains; `status` for `/progress` (P-8) | drop or rewrite `estimate` (this game's cost model) |
| `scripts/measured.py` (+ `test_measured.py`) and `measured()` in `crates/physics/tests/common/mod.rs` | writes test-measured numbers into marked spots in the documents (P-6) | its list of documents; the helper in the tests' language |
| `scripts/cleanup.py` (+ `test_cleanup.py`) | reports disk use; deletes leftovers, over-cap builds, stale caches and old uncited scratch (P-14) | `KNOWN`, `CAPS`, the rebuild after a cap, the run-cache rule |
| `scripts/affected.py`, `scripts/coverage_rustc.py` (+ `test_affected.py`) | selects heavy tests by the files each executed (LLVM coverage), writes the commit's change note (P-7) | Rust- and project-specific: the engine's crates, the heavy test program, the non-code inputs, `DERIVED`; in another language, reuse the design |
| `mypy.ini`, `pytest.ini` | type checks and tests for the scripts (P-9) | nothing |
| `rust-toolchain.toml`, `deny.toml`, `clippy.toml`, `.config/nextest.toml` | pinned toolchain, dependency hygiene, deterministic math (P-10, P-12) | versions; the disallowed list fits any deterministic f64 code |
| `.claude/skills/progress/SKILL.md` | `/progress`, a read-only report of the long runs | nothing |

Copy the tools rather than linking to them: each repository's copy then evolves with its
own rules.
