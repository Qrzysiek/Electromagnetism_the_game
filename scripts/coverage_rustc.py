"""rustc for the workspace's own crates, instrumented for coverage (`scripts/affected.py`,
DEC-16): cargo runs it as `RUSTC_WORKSPACE_WRAPPER`, for the workspace members only, so the
dependencies (glam's vector algebra, libm's functions, rayon), which only `Cargo.lock`
decides, run without counters in the engine's inner loops.

    coverage_rustc.py <rustc> <args...>
"""
import subprocess
import sys

args = sys.argv[1:]
# Cargo also asks the wrapper for rustc's version and target information (`-vV`, `-` with
# `--print`): pass those through untouched.
if "-" not in args and not any(a.startswith("--print") for a in args) and "-vV" not in args:
    args = args + ["-C", "instrument-coverage"]
sys.exit(subprocess.run(args).returncode)
