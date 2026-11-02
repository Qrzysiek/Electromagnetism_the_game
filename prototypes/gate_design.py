"""Gate (transistors arc intro) design scans: a channel 8 x 2 (density 1, mean free path 0.5)
between a grounded source and a drain at a fixed 10, gates 0.5 above and below on the
player's supply. The channel's charging time R C ~ L^2: the prototype's 16-long strip had
not settled by t = 30. Placements: the gate over +20 .. -100; the current over t = 20-60.

    python target/scratch/gate_design.py
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = open(os.path.join(ROOT, "scripts", "levels.py"), encoding="utf-8").read()
g = {"__file__": os.path.join(ROOT, "scripts", "levels.py"), "__name__": "levels_scratch"}
exec(compile(src, "levels_scratch", "exec"), g)
L = g
m = dict(density=1.0, thermal_speed=2.0, mean_free_path=0.5, core_radius=0.005,
         carriers_per_area=100.0)
lvl = L["level"](
    "Gate (design)", "scratch",
    grid=(20, 16), c=None, t_max=80.0,
    # 0: the source (grounded); 1: the drain (a fixed source); 2, 3: the gates (one supply).
    electrodes=[L["contact"](5, 8, 2), L["contact"](15, 8, 2),
                L["plate"](10, 10, 6, thickness=1.0, tunable=True),
                L["plate"](10, 6, 6, thickness=1.0)],
    supplies=L["slider"](100.0, 20.0, -40.0),
    tube=L["tube"](0, None, 0, 0, 0, 0, step=0.04,
                   circuits=[L["tube_circuit"](1, L["dc"](10.0), resistance=0.0),
                             L["tube_circuit"](2, L["dc"](0.0), resistance=0.0, linked=[3])],
                   goals=[L["current_goal"](1, 0.1, 10.0, 20.0, 60.0, carriers="net")],
                   metals=[L["metal"]((6, 7), (14, 9), [0, 1], seed=3, **m)]),
    reference=[L["supply"](10, 10, 0.0)])
lvl = L["resolve_auto"](lvl)
json.dump(lvl, open(os.path.join(ROOT, "target", "scratch", "gate.json"), "w"), indent=1)
places = []
for vg in [20.0, 0.0, -20.0, -40.0, -60.0, -80.0, -100.0]:
    p = json.loads(json.dumps(lvl["reference_solution"]))
    p[0]["value"] = vg
    places.append(p)
json.dump(places, open(os.path.join(ROOT, "target", "scratch", "gate_vary.json"), "w"))
print("ok")
