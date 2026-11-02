"""Wheatstone bridge prototype (metals arc): Magnetoresistance's two strips (a long 16 x 4, a
short wide 4 x 10) on one fixed supply V at their right ends, their left ends grounded through
R0 (the long strip's) and R2 (the short one's); the output is the left ends' difference
V0 - V2 over t = 10..80, the player's solenoid B. Writes target/scratch/br_<tag>.json and
br_<tag>_ref.json (one placement per B).

    python target/scratch/bridge_proto.py tag V R0 R2 B1 [B2 ...]   (EM_END: the window's end, 80)
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.argv.append("--text")  # no solver at import
import levels as L  # noqa: E402

sys.argv.pop()
tag = sys.argv[1]
v, r0, r2 = (float(a) for a in sys.argv[2:5])
fields = [float(a) for a in sys.argv[5:]]
end = float(os.environ.get("EM_END", "80"))
lvl = L.level(
    "Bridge (prototype)", "scratch",
    grid=(24, 20), c=None, t_max=end + 20.0,
    electrodes=[L.contact(3, 15, 4), L.contact(21, 15, 4),
                L.contact(9, 6, 10), L.contact(15, 6, 10)],
    solenoid_fields=L.slider(16.0, fields[0]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP, solenoid=(23, 19),
                circuits=[L.tube_circuit(1, L.dc(v), resistance=0.0, linked=[3]),
                          L.tube_circuit(0, L.dc(0.0), resistance=r0),
                          L.tube_circuit(2, L.dc(0.0), resistance=r2)],
                potentials=[L.potential_goal(0, "mean", -1e3, 1e3, 10.0, end, reference=2),
                            L.potential_goal(0, "mean", -1e3, 1e3, 10.0, end),
                            L.potential_goal(2, "mean", -1e3, 1e3, 10.0, end)],
                metals=[L.metal((4, 13), (20, 17), [0, 1], seed=11, **L.METAL),
                        L.metal((10, 1), (14, 11), [2, 3], seed=12, **L.METAL)]),
    reference=[L.supply(23, 19, fields[0])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"br_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(23, 19, b)] for b in fields],
          open(os.path.join(out, f"br_{tag}_ref.json"), "w"))
print(f"br_{tag}.json: V {v}, R0 {r0}, R2 {r2}, B {fields}")
