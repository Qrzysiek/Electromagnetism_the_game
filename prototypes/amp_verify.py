"""Runs a prototype level's placements at the preview (one process, all placements) and at
--refine 2 (one process per placement) in parallel; prints each run's potential measures.

    python target/scratch/amp_verify.py tag [refine]
"""
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S = os.path.join(ROOT, "target", "scratch")
GEN = os.path.join(ROOT, "target", "claude", "release", "generator.exe")
tag = sys.argv[1]
refine = sys.argv[2] if len(sys.argv) > 2 else "2"
lvl = os.path.join(S, f"amp_{tag}.json")
places = json.load(open(os.path.join(S, f"amp_{tag}_ref.json")))
t0 = time.time()
runs = [("preview", os.path.join(S, f"amp_{tag}_ref.json"), [])]
for p in places:
    f = os.path.join(S, f"amp_{tag}_ref{p[0]['value']:g}.json")
    json.dump([p], open(f, "w"))
    runs.append((f"refine {refine} at {p[0]['value']:g}", f, ["--refine", refine]))
procs = []
for name, f, extra in runs:
    out = os.path.join(S, f"amp_{tag}_{len(procs)}.out")
    procs.append((name, out, subprocess.Popen([GEN, "tube-scan", lvl, f, *extra],
                                              stdout=open(out, "w"),
                                              stderr=subprocess.DEVNULL)))
for name, out, p in procs:
    p.wait()
    print(f"{name} ({time.time() - t0:.0f} s):")
    for line in open(out):
        m = re.search(r"^(\d+):.*?(potential \[[^\]]*\])", line)
        if m:
            print(f"  {m.group(1)}: {m.group(2)}")
