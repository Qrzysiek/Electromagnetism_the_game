"""Brings every heavy test's verdict and coverage up to date, starting at once: `affected.py
run` (the selected tests' verdicts with the plain program at full speed, then their coverage
with the instrumented one, longest first), tracked by `progress.py` (`/progress`; the job
`coverage-refresh`). It takes the whole machine for hours (2026-10-27, all 521 tests stale:
about 1.5 h of verdicts and 10 h of coverage); on Windows it keeps the machine from sleeping
while it runs (not the display), unless `--allow-sleep`. Interrupted, start it again: the
runner keeps what is done and redoes only what changed (DEC-23). The records it writes are the
main checkout's, shared by every worktree (`affected.py`). Any arguments after `--` go to
`affected.py run` (e.g. `-- --jobs 2`).

    python scripts/refresh_coverage.py          # or double-click scripts/refresh_coverage.cmd
"""

import argparse
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The items `affected.py run` counts: every heavy test's verdict, then the coverage of those
# whose coverage is not another's (`repeated::` takes `golden::`'s): an estimate for the
# progress bar, not a limit.
ITEMS = 521 + 353
SECONDS_EACH = 47


def keep_awake():
    """While this process runs, Windows does not sleep (`SetThreadExecutionState`; released at
    exit). Elsewhere nothing."""
    if sys.platform != "win32":
        return
    import ctypes

    es_continuous, es_system_required = 0x80000000, 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(es_continuous | es_system_required)


def main():
    parser = argparse.ArgumentParser(
        description="Bring every heavy test's verdict and coverage up to date, starting at once "
        "(hours; the whole machine)."
    )
    parser.add_argument("--allow-sleep", action="store_true",
                        help="let the machine sleep as usual while it runs")
    parser.add_argument("rest", nargs=argparse.REMAINDER,
                        help="after --: arguments for affected.py run")
    args = parser.parse_args()
    rest = [a for a in args.rest if a != "--"]
    if not args.allow_sleep:
        keep_awake()
    command = [
        sys.executable, os.path.join("scripts", "progress.py"), "run", "coverage-refresh",
        "--total", str(ITEMS), "--each", str(SECONDS_EACH), "--match", r"^\[\d+/\d+\] ",
        "--", sys.executable, os.path.join("scripts", "affected.py"), "run", *rest,
    ]
    print("[refresh] " + " ".join(command), flush=True)
    return subprocess.call(command, cwd=ROOT)


if __name__ == "__main__":
    sys.exit(main())
