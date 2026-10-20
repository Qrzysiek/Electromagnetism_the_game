---
name: progress
description: Show where this repository's long computations are (test suites, level scans, generator window/analyze/check, design runs, chains of chunks): done/total, time per item against its estimate, ETA. Read-only. Use when the user types /progress.
---

# /progress

Read-only: report, change nothing (no jobs started, stopped, registered or edited, no files
written).

Run `python scripts/progress.py status` (`--all` also lists jobs finished over an hour ago;
the command only reads `target/progress/`) and report each job in a line: done/total,
elapsed, measured time per item against the estimate, ETA; for a chain, its chunks done and
the ETA of all that remain. Point out what looks wrong:
- `SLOW`: more than 1.5 times the estimate per item; say why if it is evident (other jobs
  running, idle priority, a configuration far from a solution).
- `lost`: the process ended without the wrapper recording it (killed, crashed); a lost chain
  resumes when its command is run again.
- a log quiet for much longer than one item should take: possibly stuck.

If a long computation is running that the wrapper does not track, say so (with its output
file and how far it seems to be from its last lines) instead of registering it.
