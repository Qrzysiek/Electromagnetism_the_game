"""Audit every tube level's field display rate (the magnetic and energy-flow maps).

The game's playback waits for the field (`tube_views.rs`), so a field slower than one time
unit a second stutters the animation even when the preview is fast. Runs `generator
tube-field <level> --every <end>` (the reference's preview run to the display end, then the
field over it in one go), parses its two timings, and appends one JSON line per level to
target/scratch/field_audit.jsonl (levels already there are skipped: a restart resumes).
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "target" / "release" / "generator.exe"
OUT = ROOT / "target" / "scratch" / "field_audit.jsonl"


def display_end(d: dict) -> float:
    """`Level::tube_display_end`: the verification's end (the goals' latest) or t_max."""
    t = d["physics"]["t_max"]
    tube = d["tube"]
    for key in ("heat_goals", "goals", "potential_goals", "oscillation_goals"):
        for g in tube.get(key) or []:
            if isinstance(g, dict) and "end" in g:
                t = max(t, g["end"])
    return t


def main() -> None:
    done = set()
    if OUT.exists():
        for line in OUT.read_text(encoding="utf-8").splitlines():
            done.add(json.loads(line)["file"])
    files = sorted(
        (p for p in (ROOT / "levels").glob("*.json") if re.match(r"\d+_", p.name)),
        key=lambda p: int(p.name.split("_")[0]),
    )
    for p in files:
        d = json.loads(p.read_text(encoding="utf-8"))
        if not d.get("tube") or p.name in done:
            continue
        end = display_end(d)
        r = subprocess.run(
            [str(GEN), "tube-field", str(p), "--every", str(end), "--end", str(end)],
            cwd=ROOT,
            env=dict(os.environ),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        run = re.search(r"^run ([\d.]+) s", r.stdout, re.M)
        field = re.search(r"^field ([\d.]+) s", r.stdout, re.M)
        grid = re.search(r"field grid (\d+) x (\d+)", r.stdout)
        run_s = float(run.group(1)) if run else float("nan")
        field_s = float(field.group(1)) if field else float("nan")
        rec = {
            "file": p.name,
            "end": end,
            "run_s": run_s,
            "field_s": field_s,
            "grid": [int(grid.group(1)), int(grid.group(2))] if grid else None,
            "run_rate": end / run_s if run_s > 0 else None,
            "field_rate": end / field_s if field_s > 0 else None,
            "exit": r.returncode,
        }
        with OUT.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        fr = f"{rec['field_rate']:.2f}" if rec["field_rate"] else "?"
        print(
            f"[field] {p.name}: run {run_s:.1f} s, field {field_s:.1f} s to t = {end}, "
            f"{fr} units/s, exit {r.returncode}",
            flush=True,
        )
        if r.returncode != 0:
            print(r.stderr[-2000:], file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
