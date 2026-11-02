"""A scratch build of a levels.py function (design only), and placements to scan.

    python target/scratch/level_scratch.py <function> [--set k=v,k=v] [--vary k v1 v2 ...]

Writes target/scratch/<function>.json and, with --vary, target/scratch/<function>_vary.json:
placements of the reference with element k (1-based) at each value; --set fixes other
elements first (k=v, 1-based).
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = open(os.path.join(ROOT, "scripts", "levels.py"), encoding="utf-8").read()
g = {"__file__": os.path.join(ROOT, "scripts", "levels.py"), "__name__": "levels_scratch"}
exec(compile(src, "levels_scratch", "exec"), g)
name = sys.argv[1]
lvl = g["resolve_auto"](g[name]())
out = os.path.join(ROOT, "target", "scratch", f"{name}.json")
json.dump(lvl, open(out, "w", encoding="utf-8"), indent=1)
ref = lvl["reference_solution"]
print("reference", json.dumps(ref))
base = json.loads(json.dumps(ref))
if "--set" in sys.argv:
    for kv in sys.argv[sys.argv.index("--set") + 1].split(","):
        k, v = kv.split("=")
        base[int(k) - 1]["value"] = float(v)
if "--vary" in sys.argv:
    i = sys.argv.index("--vary")
    k = int(sys.argv[i + 1])
    values = [float(v) for v in sys.argv[i + 2:] if not v.startswith("--")]
    places = []
    for v in values:
        p = json.loads(json.dumps(base))
        p[k - 1]["value"] = v
        places.append(p)
    path = os.path.join(ROOT, "target", "scratch", f"{name}_vary.json")
    json.dump(places, open(path, "w"))
    print(len(places), "placements to", path)
