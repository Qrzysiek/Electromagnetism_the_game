"""Scratch: the Magnetoresistance introduction (a long strip and a short wide one on one
supply in the player's solenoid), for tube-scan.

usage: design_mr_intro.py <placements "V,B;V,B;..."> [end] [width_b] [seed]
"""
import importlib.util
import json
import sys

spec = importlib.util.spec_from_file_location("levels", "scripts/levels.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)

places = [tuple(float(x) for x in p.split(",")) for p in sys.argv[1].split(";")]
end = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
wb = int(sys.argv[3]) if len(sys.argv) > 3 else 10
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 11
yb = 1 + (10 - wb) // 2

lvl = L.level(
    "Magnetoresistance", "design", grid=(24, 20), c=None, t_max=end,
    # 0, 1: the long strip's contacts; 2, 3: the short wide strip's.
    electrodes=[L.contact(3, 15, 4), L.contact(21, 15, 4, tunable=True),
                L.contact(9, yb + wb // 2, wb), L.contact(15, yb + wb // 2, wb)],
    supplies=L.slider(16.0, 6.0, 4.0),
    solenoid_fields=L.slider(7.0, 2.0, 3.0),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP, solenoid=(23, 19),
                circuits=[L.tube_circuit(1, L.dc(0.0), resistance=0.0, linked=[3])],
                goals=[L.current_goal(0, -100.0, 100.0, 10.0, end, carriers="net"),
                       L.current_goal(2, -100.0, 100.0, 10.0, end, carriers="net")],
                metals=[L.metal((4, 13), (20, 17), [0, 1], seed=seed, **L.METAL),
                        L.metal((10, yb), (14, yb + wb), [2, 3], seed=seed + 1, **L.METAL)]),
    reference=[L.supply(21, 15, places[0][0]), L.supply(23, 19, places[0][1])])
lvl = L.resolve_auto(lvl)
json.dump(lvl, open("target/scratch/mr_intro_design.json", "w", encoding="utf-8", newline="\n"),
          indent=1)
json.dump([[{"node": [21, 15, 0], "kind": "supply", "value": v},
            {"node": [23, 19, 0], "kind": "supply", "value": b}] for v, b in places],
          open("target/scratch/mr_intro_places.json", "w", encoding="utf-8", newline="\n"))
print("written", len(places))
