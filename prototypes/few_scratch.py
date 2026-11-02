"""A scratch build of Few carriers or many from the pending levels.py (design only: levels/ stays
untouched while the suite reads it), and placements to scan.

    python target/scratch/few_scratch.py [B values] [--supplies U1,U2]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = open(os.path.join(ROOT, "target", "scratch", "levels_few_or_many.py.pending"),
           encoding="utf-8").read()
g = {"__file__": os.path.join(ROOT, "scripts", "levels.py"), "__name__": "levels_pending"}
exec(compile(src, "levels_pending", "exec"), g)
lvl = g["resolve_auto"](g["few_or_many"]())
out = os.path.join(ROOT, "target", "scratch", "few.json")
json.dump(lvl, open(out, "w", encoding="utf-8"), indent=1)
ref = lvl["reference_solution"]
print("reference", json.dumps(ref))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
u = [30.0, 30.0]
if "--supplies" in sys.argv:
    u = [float(x) for x in sys.argv[sys.argv.index("--supplies") + 1].split(",")]
    args = [a for a in args if a != sys.argv[sys.argv.index("--supplies") + 1]]
if args:
    places = []
    for b in [float(x) for x in args]:
        p = json.loads(json.dumps(ref))
        p[0]["value"], p[1]["value"], p[2]["value"] = u[0], b, u[1]
        places.append(p)
    json.dump(places, open(os.path.join(ROOT, "target", "scratch", "few_b.json"), "w"))
    print(len(places), "placements to target/scratch/few_b.json")
