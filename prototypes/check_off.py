"""Checks levels with every player control off (an empty placement): each must fail.

    python target/scratch/check_off.py 105_rectifier 107_klystron ...
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GEN = os.path.join(ROOT, "target", "claude", "release", "generator.exe")
for stem in sys.argv[1:]:
    d = json.load(open(os.path.join(ROOT, "levels", f"{stem}.json"), encoding="utf-8"))
    d["reference_solution"] = []
    path = os.path.join(ROOT, "target", "scratch", f"off_{stem}.json")
    json.dump(d, open(path, "w", encoding="utf-8"))
    out = subprocess.run([GEN, "check", path], capture_output=True, text=True).stdout
    lines = [l.strip()[:150] for l in out.splitlines() if "shot " in l]
    solved = bool(lines) and all("Met" in l and "Verified" in l for l in lines)
    print(f"== {stem}: {'SOLVES with everything off' if solved else 'fails, as it must'}")
    for l in lines:
        print("   ", l)
