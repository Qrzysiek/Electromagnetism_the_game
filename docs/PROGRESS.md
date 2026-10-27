# Long computations

How this repository's long computations are run so that they can be followed (`/progress`)
and lose at most an hour of work when interrupted (the owner, 2026-10-20). The rules are in
`CLAUDE.md`; this is the reference. The tool: `scripts/progress.py` (`run`, `chain`,
`status`, `watch`, `estimate`; its docstring).

## Launching a long computation (more than a few minutes)

1. **Estimate before starting**, per item and in total; within 50 % is good enough. Tube levels
   with metals: `python scripts/progress.py estimate <level.json>` (carriers × steps ×
   1.5e-6 s); otherwise the timings below, or time one item first.
2. **Run it through the wrapper**, counting this computation's items with a pattern that
   matches one line per finished item (each computation has its own: see the table):
   `python scripts/progress.py run <job> --total N --each S --match "<regex>" [--idle] -- <command>`
   Its output goes to `target/progress/<job>.log`.
3. **Say it in the message**: what runs, N items × ~S s ≈ total, and that `/progress` shows
   where it is.
4. **When woken or asked, check**: measured per item against the estimate. Over 1.5 times
   slower: find out why before waiting longer.

## Chunks: never more than an hour of work at risk

A computation that would take an hour or more is split into consecutive chunks of at most
an hour each (less is better), each saving its results before the next starts, so that an
interruption (the owner needs the machine, a crash, a restart) loses at most one chunk. The
chunks run one after another by themselves, without the owner's input:

    python scripts/progress.py chain <job> --plan target/progress/<job>.plan.json [--idle]

The plan: `{"chunks": [{"name": "...", "cmd": [...], "total": N, "each": S, "match": "..."}, ...]}`.
Each chunk runs as `run` would (its log `target/progress/<job>-<k>.log`, kept: often it is the
chunk's result); each one that finishes is recorded in `target/progress/<job>.chain.json`; a
failed chunk stops the chain (look at its log, fix, and run the same command again: it
resumes at the first chunk not finished, as after any interruption). Launch the chain itself
in the background. A chunk that writes its own results should write them as they come, so
that even a chunk cut short keeps most of its work. Ways to split this repository's long runs:
- `generator window`: one slider at a time (`--element K`) and its range in parts
  (`--from`/`--to`), each part's intervals appended to a results file.
- `generator analyze`: one level per run (and fewer `--samples`/`--runs` if a level alone
  exceeds an hour).
- Golden hashes for every level (`UPDATE_GOLDEN=1 cargo test --release -p level --test
  per_level`, about an hour): in groups with `GOLDEN_ONLY=<part of the names>` (the file
  keeps the other entries).
- The per-level checks (`-- --ignored`): by check or by a range of levels, each a filter
  (`golden::`, `reference::1`); a check that passed is skipped when rerun, so a chain of
  them resumes.
- `generator tube-scan`: a placements file per chunk (or one placement per process, run in
  parallel when the machine is free).
- Design scripts: write each result as it comes (append to a JSON lines file) and skip the
  ones already written when restarted.

## Runs of independent items: `scripts/runner.py` (DEC-23)

Anything that runs many independent items (tests, a scan's points, design runs) goes through
`scripts/runner.py`: from a plan (`runner.py run <job> --plan plan.json`: the items' commands
and their inputs, the job's own inputs such as the program's binary) or as a library
(`affected.py` does). Each item's outcome is saved as it finishes, keyed by the fingerprint of
its command and inputs. While it runs, from another terminal:

    python scripts/runner.py ctl <job> pause | resume | stop | skip ITEM | jobs N
    python scripts/runner.py state <job>             # outcomes, and where it stopped

An item past four times its expected time is hung: the job stops there (the next items may
hang too; `--on-hang continue` goes on), and the incident log names it with its output's end.
Fix it and run the same command again: the item that stopped the job runs first, finished
items whose inputs are unchanged are not rerun, and those the fix changed are. `state --save
FILE` keeps a ledger, `--load FILE` restores one. The heavy tests' jobs are
`affected-verdicts` and `affected-coverage`.

## Item patterns and timings (owner's 16-thread machine, full speed, 2026-10-20)

| computation | item | `--match` | per item |
|---|---|---|---|
| `cargo test --workspace` | a test binary | `^test result:` | the fast tier: the per-level checks are ignored (DEC-16) |
| `cargo test --release -p level --test per_level -- --ignored` | one check of one level (521) | `^test \S+ \.\.\. ` | seconds for most, minutes for metal tubes and beams |
| `python scripts/affected.py run [--filter F] [--part K/N]` | one heavy test, instrumented | `^\[\d+/\d+\] ` | ~6.5x the plain build on one thread; a test per hardware thread |
| `UPDATE_GOLDEN=1 cargo test --release -p level --test per_level` | a level | `^\S+: [0-9a-f]{16}$` | its golden hash: both runs of its reference |
| `generator window` (tube level) | one slider value | `^\[window\]` | one preview run (`estimate`) |
| `generator tube-scan` | one placement | `^\d+: ` | one run at `--refine` (`estimate`: preview or verification) |
| `generator check` (tube level) | the level | `^  cost:` | preview + verification (Laminations 36 + 124 s) |
| `generator analyze` | a tenth of a phase | `^\[analyze\] .*: \d+/\d+ after` | tube levels: ~1 preview per new placement |
| `python scripts/levels.py --only …` | a level | `reference (verified|NOT VERIFIED)` | its `check`, plus a search if it has no reference |
| `python scripts/pairs.py --only X` | the whole run | `^physics worker:` | a game build (~1.5 min) and a second per case |

Metal tube runs: 1.4–1.6e-6 s per carrier-step (`progress.py estimate`); the verification has
twice the carriers and twice the steps: ~4 times the preview. Several runs at once share the
cores: their times add up roughly as if run one after another.
