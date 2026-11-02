"""A scratch prototype (not a level): the field effect in the Drude metal. A strip 16 x 2
(density 1, mean free path 0.5: Debye length 0.56) between a grounded source contact and a
drain on a supply, with gates above and below it (1 away, on one supply). A negative gate
pushes the carriers out of the strip (depletion): the current falls. Placements: the drain
at 10, the gate at 0, -10, -20, -30, -40.

    python target/scratch/fet_proto.py
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
    "Field effect (prototype)", "scratch",
    grid=(24, 16), c=None, t_max=60.0,
    # 0: the source (grounded); 1: the drain (a supply); 2, 3: the gates (one supply).
    electrodes=[L["contact"](3, 8, 2), L["contact"](21, 8, 2, tunable=True),
                L["plate"](12, 11, 14, thickness=2.0, tunable=True),
                L["plate"](12, 5, 14, thickness=2.0)],
    supplies=L["slider"](40.0, 10.0, 20.0),
    tube=L["tube"](0, None, 0, 0, 0, 0, step=0.04,
                   circuits=[L["tube_circuit"](2, L["dc"](0.0), resistance=0.0, linked=[3])],
                   goals=[L["current_goal"](1, 0.1, 10.0, 10.0, 30.0, carriers="net")],
                   metals=[L["metal"]((4, 7), (20, 9), [0, 1], seed=3, **m)]),
    reference=[L["supply"](21, 8, 10.0), L["supply"](12, 11, 0.0)])
lvl = L["resolve_auto"](lvl)
json.dump(lvl, open(os.path.join(ROOT, "target", "scratch", "fet.json"), "w"), indent=1)
places = []
for vg in [0.0, -10.0, -20.0, -30.0, -40.0]:
    p = json.loads(json.dumps(lvl["reference_solution"]))
    p[1]["value"] = vg
    places.append(p)
json.dump(places, open(os.path.join(ROOT, "target", "scratch", "fet_vary.json"), "w"))
print("ok", json.dumps(lvl["reference_solution"]))
