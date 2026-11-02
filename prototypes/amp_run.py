"""Runs several amp_proto.py variants in parallel; prints per bias the drain's mean, its
amplitude, the gain (amplitude / A) and the phase.

    python target/scratch/amp_run.py "biases" variant [variant ...]
    variant: tag:L:R_L:V_DD:A:T:periods:settle   (biases comma-separated)
"""
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S = os.path.join(ROOT, "target", "scratch")
GEN = os.path.join(ROOT, "target", "claude", "release", "generator.exe")
biases = sys.argv[1].split(",")
procs = []
t0 = time.time()
for v in sys.argv[2:]:
    tag, *params = v.split(":")
    subprocess.run([sys.executable, os.path.join(S, "amp_proto.py"), tag, *params, *biases],
                   check=True, stdout=subprocess.DEVNULL)
    out = open(os.path.join(S, f"amp_{tag}.out"), "w")
    procs.append((tag, params, subprocess.Popen(
        [GEN, "tube-scan", os.path.join(S, f"amp_{tag}.json"),
         os.path.join(S, f"amp_{tag}_ref.json")], stdout=out, stderr=subprocess.DEVNULL)))
for _, _, p in procs:
    p.wait()
print(f"{time.time() - t0:.0f} s")
num = r"Some\(([-0-9.e]+)\)"
for tag, params, p in procs:
    amp = float(params[3])
    print(f"{tag} (L R_L V_DD A T periods settle = {' '.join(params)})")
    text = open(os.path.join(S, f"amp_{tag}.out")).read()
    for b, line in zip(biases, [x for x in text.splitlines() if "potential" in x]):
        m = re.search(rf"potential \[{num}, {num}, {num}\]", line)
        if m:
            mean, a, ph = (float(g) for g in m.groups())
            print(f"  bias {b:>6}: mean {mean:7.2f}  amplitude {a:6.3f}  gain {a / amp:5.2f}  "
                  f"phase {ph:7.1f}")
        else:
            print(f"  bias {b:>6}: {line[:120]}")
