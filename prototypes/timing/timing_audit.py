"""Audit every tube level's preview time and streaming rate (HANDOFF.md: the timing rules).

Runs `generator check` on each tube level with a fresh run cache (each run computed once,
its time measured alone), parses the cost meters, and appends one JSON line per level to
target/scratch/timing_audit.jsonl (levels already there are skipped: a restart resumes).
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "target" / "release" / "generator.exe"
OUT = ROOT / "target" / "scratch" / "timing_audit.jsonl"
CACHE = ROOT / "target" / "scratch" / "timing_cache"


def seconds(text: str) -> float:
    m = re.match(r"([\d.]+)\s*(ms|s|min|µs)", text.strip())
    if not m:
        return float("nan")
    v, u = float(m.group(1)), m.group(2)
    return v * {"ms": 1e-3, "s": 1.0, "min": 60.0, "µs": 1e-6}[u]


def main() -> None:
    done = set()
    if OUT.exists():
        for line in OUT.read_text(encoding="utf-8").splitlines():
            done.add(json.loads(line)["file"])
    files = sorted(
        (p for p in (ROOT / "levels").glob("*.json") if re.match(r"\d+_", p.name)),
        key=lambda p: int(p.name.split("_")[0]),
    )
    tubes = []
    for p in files:
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("tube"):
            tubes.append((p, d))
    env = dict(os.environ, EM_RUN_CACHE=str(CACHE))
    for p, d in tubes:
        if p.name in done:
            continue
        r = subprocess.run(
            [str(GEN), "check", str(p)],
            cwd=ROOT,
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        cost = next((l for l in r.stdout.splitlines() if l.startswith("  cost:")), "")
        meters = dict(
            (m.group(1), m.group(2))
            for m in re.finditer(r"(Preview|Verification|Steps) ([^(]+)\(", cost)
        )
        t_max = d["physics"]["t_max"]
        preview = seconds(meters.get("Preview", ""))
        verification = seconds(meters.get("Verification", ""))
        rec = {
            "file": p.name,
            "name": d["name"],
            "t_max": t_max,
            "metal": bool(d["tube"].get("metals")),
            "c": d["tube"].get("speed_of_light"),
            "preview_s": preview,
            "verification_s": verification,
            "steps": meters.get("Steps", "").strip(),
            "rate": t_max / preview if preview > 0 else None,
            "verdict": [l.strip() for l in r.stdout.splitlines() if l.startswith("  shot")],
            "exit": r.returncode,
        }
        with OUT.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        rate = f"{rec['rate']:.2f}" if rec["rate"] else "?"
        print(
            f"[level] {p.name}: preview {preview:.1f} s, verification {verification:.1f} s, "
            f"t_max {t_max}, {rate} units/s, exit {r.returncode}",
            flush=True,
        )
        if r.returncode != 0:
            print(r.stderr[-2000:], file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
