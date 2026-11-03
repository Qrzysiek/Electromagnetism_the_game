<!--
A CLAUDE.md template from the playbook of "Electromagnetism – the game" (playbook/README.md).
Layers: A any repository; B numerical or scientific code; C a Rust workspace; D an
interactive simulation or educational game; E Claude Code on Windows. Keep the layers that
apply, fill the <…>, delete the rest and this comment. P-n: the reason, in
playbook/DECISIONS.md (copy the ones you keep into the new repository's docs).
-->
# Project rules

<What the project is, in a sentence.> `<SPEC.md>` describes what it does, `<MODEL.md>` how
it works and how it is validated, `README.md` how to run it.

## A. Any repository

### Decisions and documents
- Every model, algorithm or design chosen for a feature, a tool or a piece of content gets
  an entry in `docs/DECISIONS.md` in the same commit (P-0, P-1). The entry gives the choice,
  its purpose, the alternatives considered and why not, and where it lives. When a choice is
  changed or dropped, its entry is not rewritten: append a dated history line with what
  changed and why (the evidence), and update its status.
- A rule that comes from the owner cites the date and the reason ("the owner, <date>: …");
  a later ruling cites what it supersedes (P-2).
- `docs/HANDOFF.md` holds where the work stands, what is verified and how, the open
  questions to the owner, and the next steps. Read it at a session's start instead of old
  transcripts; update it at the end (P-3).
- Open problems (questions the literature does not answer and that could not be derived in
  reasonable time) go into `docs/OPEN_PROBLEMS.md`, with what was searched and tried and the
  workaround in use (the `open-problems` skill; P-4).
- Audits are archived in `audits/` with a verdict per remark in each folder's
  `RESOLUTION.md`; the next audit reviews `git diff <last fix>..HEAD` (P-5).
- Numbers that documents quote from tests are marked (`<!--m:<id>-->…<!--/m-->`) and written
  by `python scripts/measured.py update` from what the tests record (P-6).

### How to work
- Find the cause first: explain a failure, a slowdown or a surprising number before fixing
  or working around it, and name the cause in the commit.
- Measure before adopting. A tool or optimisation proposed for speed is first measured on
  this machine against the status quo, and the numbers go into its decision entry, adopted
  or not.
- Use the process tools as intended. Never bypass one for speed without restoring what it
  needs before the next change. If a tool is wrong, fix the tool. Don't overdo the process:
  full runs happen on a trigger, not a cadence.
- Ask the owner only about decisions that are theirs (trade-offs, scope, taste). Give the
  options with their measured costs and a recommendation, and record the answer where the
  rule lives.
- Never send clicks or keys to the desktop. Visual checks go through the program's capture
  mode (P-16).

### Tests
- The fast tier (`<fast test command>`, minutes) runs before every commit and in CI. The
  heavy tier has one test per check and item, ignored by default and runnable by name
  (`<heavy test command> <part of a name>`), and is selected by what each test executed
  (`python scripts/affected.py run`, with coverage; P-7). Run everything when a global input
  changes (toolchain, lock file, build profile), on request, or when a bug escapes the
  selection (log what the records missed).
- Golden outputs change only when the behaviour changes intentionally. Updating them is a
  mode of the test program (`<UPDATE_GOLDEN=1 …>`).
- Scripts are type-checked and tested: `python -m mypy`, `python -m pytest` (P-9).

### Before committing
- `<formatter>`, `<linter with warnings as errors>`, `<fast tier>`, then the selected heavy
  tests (`affected.py run`). For script changes, also mypy and pytest. Paste
  `python scripts/affected.py note` (what changed, what it affects, what ran) into the
  commit message.
- Commit messages give what changed and why, with the evidence (measured numbers), and end
  with the co-author line used in the history.
- Work on a branch in a worktree. Finished work merges to main fast-forward only. Then
  rebuild what the owner runs from the main checkout, push main and the branch, and run
  `python scripts/cleanup.py run` (P-15).

### Long computations (`docs/PROGRESS.md`; `/progress` shows where they are; P-8)
- Estimate per item and in total before launching (within 50 % is fine), and say it.
- Run it through `scripts/progress.py run <job> --total N --each S --match "<regex>"`.
- Split anything of an hour or more into chunks of at most an hour that save their results
  and run one after another by themselves (`progress.py chain`). An interruption must not
  lose more than an hour of work.
- A run of many independent items goes through `scripts/runner.py`:
  - each outcome is saved as it comes, keyed by the item's inputs;
  - `runner.py ctl <job> pause|resume|stop|skip ITEM` works while it runs;
  - a hang stops the job with a record of where;
  - after the fix, the same command continues and redoes only what the fix changed.
- A hang limit is a multiple of an estimate, so make sure the estimate measures the same
  work (a cached or recorded result is not a computation's time).
- Keep the machine awake during an unattended run, and let it sleep after only when asked.

### Disk (P-14)
- Builds go only to `<the build directories>`. A temporary build directory goes under
  `target/temp/` and is removed when done.
- A session's scripts and outputs go to `target/scratch/` (not versioned). A scratch file
  that a tracked file cites moves to `prototypes/` in the same commit.
- `python scripts/cleanup.py` shows what takes space. Run `cleanup.py run` after each merge
  to main and whenever a checkout's `target/` passes <20 GB>. It deletes temporary and
  unregistered entries, builds over their caps, stale cache records, uncited scratch older
  than 14 days, and removed worktrees' leftovers. A tool that writes a new folder under
  `target/` registers it in the script's `KNOWN`.

## B. Numerical or scientific code
- Any change to the model or the numerics updates `<MODEL.md>` in the same commit: the
  model, the method, and the measured test results.
- Validation tests compare against analytic or independent results (a CAS at high
  precision, an independent implementation; the `wolfram-engine` skill). Set thresholds from
  the target accuracy before measuring. If a threshold must change, document why; never
  loosen one just to make a test pass (P-13).
- Symmetry and conservation tests are mandatory for every algorithm that implements the
  model:
  - the system's symmetries (a symmetric setup gives a symmetric solution; any breaking is
    explained);
  - those of the equations of motion and the conservation laws they imply (energy,
    momentum, angular momentum, charge);
  - the background's: translations, rotations, boosts;
  - the discrete ones (parity, time reversal, charge conjugation, their products) and the
    model's scalings.

  Each is tested against the exact transformed solution with a threshold from the target
  accuracy. `<MODEL.md>` (§ Symmetries) records how well and over what time each holds
  (P-13).
- The authoritative path is deterministic (P-12): f64 on the CPU, no fused multiply-add, a
  fixed summation order, and elementary functions from a portable library, never the
  platform's C library. The GPU is for visuals only. Caches of runs are keyed by every input,
  the program's executable included.
- Every number in an explanatory document is computed by its own script, independently of
  the program, and compared with the program's output. A failed comparison is fixed, never
  widened.
- A solver or feature with a cost reports it against the budgets (`<cost module>`).

## C. A Rust workspace
- The toolchain is pinned (`rust-toolchain.toml`): upgrade it deliberately, by bumping the
  version and running every test. CI checks the dependencies (`deny.toml`: licences,
  advisories), unused dependencies (cargo-machete) and the documentation's links (P-10).
- `clippy.toml` disallows the standard library's elementary functions on the authoritative
  path (use the `libm` crate; P-12). Run `cargo clippy --workspace --all-targets -- -D
  warnings`.
- Make the compiler enforce completeness. Destructure the types that an editor, a serializer
  or a check must handle exhaustively, with no `..`, so a new field fails to compile until it
  is handled or explicitly marked as not handled, with a reason.
- Structure for safety (P-11): compute a quantity in one place for every consumer; give
  constructors every optional part empty, with struct-update literals; use a named struct
  instead of several alike positional numbers; share test helpers in `tests/common`. Split a
  large file by concern as a pure move (visibility kept exactly, paths re-exported), verified
  bit for bit.

## D. An interactive simulation or educational game (P-17)
- Every shipped item has a verified reference solution, checked in the heavy tier.
  `<generator script>` is the single source of truth for the content files; never edit
  those by hand.
- Everything shipped can be reproduced by hand in the sandbox. The editor edits every field
  of the format (exhaustive destructuring, C), and its ranges are shared constants that a
  test checks on every shipped file. A new element or feature needs its sandbox control in
  the same change.
- Every element of the interface is explained, without clutter. A mark drawn on the view
  that is not a hoverable object gets a legend row. Panel text is a short line with the
  detail on hover, and controls carry hover text.
- Every item says what it represents (the real thing, what carries over, what is simplified
  or left out) and shows its scales (typical real magnitudes; the dimensionless ratios that
  decide its behaviour). A test checks each item has them.
- No trivial solutions:
  - nothing solves with the controls off or at their first touch;
  - slider ranges hide the answer: no solving interval contains 0, an end, or a half's
    middle, and on a logarithmic slider none lies within 5 % of its length from those;
  - one slider scaling a response in proportion is an introduction at most; a puzzle needs
    interacting controls or goals that pin them separately.
- Showcase everything feasible. Each control spans all the behaviours the model can produce
  within reason, and the views make them visible. The answer is picked from among them, with
  a comfortable solving interval (at least about 5 % of a slider).
- Curriculum: every new element or concept gets an easy introduction first. Within a
  chapter, difficulty and the number of elements needed rise on average.
- Timing (on the owner's machine):
  - the preliminary answer within <30 s>;
  - the animation streams without stutter at speed 1, however long it is;
  - the verification within <a few minutes>;
  - a bound may be exceeded by at most 100 % where that makes the item a lot better.

  Simplify what does not fit; what cannot be simplified accurately is simulation material
  (`<the other format>`).
- Visual checks run in the capture mode: `<ENV_CAPTURE=out.png ENV_ITEM=n … program>`
  opens a state, saves a screenshot of its own window at a frame, and exits. README
  screenshots are retaken by `<script>` (P-16).

## E. Claude Code on Windows
- Write scripts with the Write tool, not shell heredocs: the Bash tool halved backslashes
  in heredocs. Python writes text files with `newline="\n"`.
- The git stash is shared by every worktree and session. Stash only with a unique tag
  (`git stash push -u -m "<tag>"`) and apply by SHA, or prefer a WIP commit.
- Bulk deletions with `rm -rf` are denied in auto mode. Use the build tool's clean
  (`cargo clean --target-dir <dir>`, `--profile dev`) or the cleanup script.
- An unattended run's launcher keeps Windows awake:
  `ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)`. To sleep the PC after a
  run, only when asked:
  `[System.Windows.Forms.Application]::SetSuspendState('Suspend', $false, $false)`.
- Machine facts (interpreter versions, which `python` has which packages) go into the
  agent's memory, not into these rules.
