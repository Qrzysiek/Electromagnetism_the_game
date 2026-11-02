"""Laminations 20 x 12: each layout's heat, preview and verification, against Ohm's law."""
import json
import os
import re

ROOT = r"C:\Users\qrzys\Documents\Projects\Electromagnetism_the_game\.claude\worktrees\tubes-next"
P = os.path.join(ROOT, "target", "progress")
S = os.path.join(ROOT, "target", "scratch")
names = list(json.load(open(os.path.join(S, "lam20_layouts.json"), encoding="utf-8")))
ns = {}
exec(open(os.path.join(S, "eddy_theory.py"), encoding="utf-8").read().split("sigma = 0.47")[0], ns)
integral = ns["integral"]


def I(w, b):
    return integral(min(w, b), max(w, b))


theory = {"uncut": I(20, 12), "along1": 2 * I(20, 6), "along2": 3 * I(20, 4),
          "along3": 4 * I(20, 3), "across1": 2 * I(10, 12), "along_across": 4 * I(10, 6)}
b2 = (0.68 * 0.125) ** 2 / 2


def heats(k):
    out = {}
    for line in open(os.path.join(P, f"lam20-{k}.log"), encoding="utf-8"):
        m = re.match(r"^(\d+): .*heat \[Some\(([-\d.e]+)\)\]", line)
        if m:
            out[int(m.group(1))] = float(m.group(2))
    return out


pre, fine = heats(1), heats(2)
u = fine.get(0)
for i, n in enumerate(names):
    t = theory.get(n)
    t_s = f"{0.47 * b2 * t:.3f}" if t else "  -  "
    print(f"{n:13s} preview {pre.get(i, float('nan')):.3f} verification {fine.get(i, float('nan')):.3f}"
          f" (vs uncut {fine.get(i, float('nan')) / u:.3f}) Ohm {t_s}")
