"""Runs bridge_proto.py with the fields split over parallel processes; prints per B the
bridge's output V0 - V2 and the two left ends' potentials.

    python target/scratch/bridge_run.py tag V R0 R2 procs B1 [B2 ...]
"""
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S = os.path.join(ROOT, "target", "scratch")
GEN = os.path.join(ROOT, "target", "claude", "release", "generator.exe")
tag, v, r0, r2, procs, *fields = sys.argv[1:]
procs = int(procs)
groups = [fields[i::procs] for i in range(procs)]
t0 = time.time()
runs = []
for i, g in enumerate(groups):
    if not g:
        continue
    t = f"{tag}{i}"
    subprocess.run([sys.executable, os.path.join(S, "bridge_proto.py"), t, v, r0, r2, *g],
                   check=True, stdout=subprocess.DEVNULL)
    out = os.path.join(S, f"br_{t}.out")
    runs.append((g, out, subprocess.Popen(
        [GEN, "tube-scan", os.path.join(S, f"br_{t}.json"), os.path.join(S, f"br_{t}_ref.json")],
        stdout=open(out, "w"), stderr=subprocess.DEVNULL)))
num = r"Some\(([-0-9.e]+)\)"
rows = []
for g, out, p in runs:
    p.wait()
    lines = [x for x in open(out) if "potential" in x]
    for b, line in zip(g, lines):
        m = re.search(rf"potential \[{num}, {num}, {num}\]", line)
        rows.append((float(b), m.groups() if m else line[:100]))
print(f"{time.time() - t0:.0f} s; V {v}, R0 {r0}, R2 {r2}")
for b, m in sorted(rows):
    if isinstance(m, tuple):
        d, a, c = (float(x) for x in m)
        print(f"  B {b:6.2f}: output {d:8.4f}  (V0 {a:7.4f}, V2 {c:7.4f})")
    else:
        print(f"  B {b:6.2f}: {m}")
