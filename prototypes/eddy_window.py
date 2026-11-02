"""Eddy currents at omega 0.125: the plate's heat at three drives, preview and verification,
in parallel (one tube-scan per drive and resolution)."""
import json
import os
import subprocess

ROOT = r"C:\Users\qrzys\Documents\Projects\Electromagnetism_the_game\.claude\worktrees\tubes-next"
GEN = os.path.join(ROOT, r"target\claude-gen\release\generator.exe")
S = os.path.join(ROOT, "target", "scratch")
level = os.path.join(S, "draft_eddy_currents.json")
procs = []
for a in (0.55, 0.68, 0.8):
    path = os.path.join(S, f"ew_{a}.json")
    json.dump([[{"node": [23, 19, 0], "kind": "supply", "value": a}]],
              open(path, "w", encoding="utf-8"))
    for refine in ("1", "2"):
        out = open(os.path.join(S, f"ew_{a}_{refine}.txt"), "w", encoding="utf-8")
        procs.append((a, refine, out, subprocess.Popen(
            [GEN, "tube-scan", level, path, "--refine", refine], cwd=ROOT, stdout=out,
            stderr=subprocess.STDOUT)))
for a, refine, out, p in procs:
    p.wait()
    out.close()
    text = open(os.path.join(S, f"ew_{a}_{refine}.txt"), encoding="utf-8").read()
    heat = [l for l in text.splitlines() if "heat [" in l]
    theory = 0.47 * 728.75 * (a * 0.125) ** 2 / 2
    print(f"A {a} refine {refine}: theory {theory:.4f} | {heat[-1][-30:] if heat else text[-300:]}")
