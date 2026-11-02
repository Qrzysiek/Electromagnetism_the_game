"""The self-inductance level's solving interval in omega: `generator check` on copies of the
level with the reference supply's omega replaced, all at once. Usage:
python target/scratch/inductance/omega_window.py LEVEL.json w1 w2 ..."""

import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
GEN = os.path.join(ROOT, "target", "release", "generator.exe")
OUT = os.path.join(ROOT, "target", "scratch", "inductance", "window")
os.makedirs(OUT, exist_ok=True)

level_path, omegas = sys.argv[1], [float(w) for w in sys.argv[2:]]
base = json.load(open(level_path, encoding="utf-8"))


def run(w):
    lvl = json.loads(json.dumps(base))
    for e in lvl["reference_solution"]:
        if e["kind"] == "supply":
            e["omega"] = w
    path = os.path.join(OUT, f"w{w:g}.json")
    json.dump(lvl, open(path, "w", encoding="utf-8"), indent=1)
    r = subprocess.run([GEN, "check", path], capture_output=True, text=True, encoding="utf-8",
                       errors="replace", cwd=ROOT)
    line = next((l for l in r.stdout.splitlines() if "shot" in l), r.stdout[-300:] + r.stderr[-300:])
    print(f"w={w:g}: {line.strip()}", flush=True)
    return w, line


with ThreadPoolExecutor(len(omegas)) as pool:
    results = list(pool.map(run, omegas))
json.dump(results, open(os.path.join(OUT, "results.json"), "w", encoding="utf-8"), indent=1)
