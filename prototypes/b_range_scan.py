"""The three metal levels' goals across the solenoid's full range (the reference's supplies),
three wrapped tube-scan jobs in parallel (`/progress` shows them)."""
import json
import os
import subprocess
import sys

ROOT = r"C:\Users\qrzys\Documents\Projects\Electromagnetism_the_game\.claude\worktrees\tubes-next"
S = os.path.join(ROOT, "target", "scratch")
GEN = os.path.join(ROOT, r"target\claude-gen\release\generator.exe")


def supply(x, y, v):
    return {"node": [x, y, 0], "kind": "supply", "value": v}


jobs = {
    "brange-hall": ("112_hall_probe.json", [supply(21, 8, 30.0)], (23, 15),
                    [0.5, 1, 2, 3, 4, 6, 8, 12, 16]),
    "brange-positive": ("113_positive_carriers.json", [supply(21, 14, 30.0), supply(21, 5, -12.0)],
                        (23, 19), [-0.5, -1, -2, -3, -4, -6, -8, -12, -16]),
    "brange-mr": ("114_magnetoresistance.json", [supply(21, 15, 7.0)], (23, 19),
                  [1, 2, 3, 4, 5, 6, 8, 10, 13, 16]),
}
procs = []
for job, (level, base, node, fields) in jobs.items():
    placements = [base + [supply(node[0], node[1], b)] for b in fields]
    path = os.path.join(S, f"{job}.json")
    json.dump(placements, open(path, "w", encoding="utf-8"))
    # Three in parallel share the cores: ~3 times one preview (~20 s alone) per item.
    procs.append(subprocess.Popen(
        [sys.executable, os.path.join(ROOT, "scripts", "progress.py"), "run", job,
         "--total", str(len(fields)), "--each", "60", "--match", r"^\d+: ", "--",
         GEN, "tube-scan", os.path.join(ROOT, "levels", level), path], cwd=ROOT))
for p in procs:
    p.wait()
print("done")
