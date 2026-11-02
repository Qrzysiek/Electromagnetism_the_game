"""scripts/levels.py: Laminations' core 20 x 12 (30 carriers per area), three free cuts."""
p = "scripts/levels.py"
s = open(p, encoding="utf-8").read()


def sub(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)


sub('''        max_cuts=2,
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, b_z=0.68,
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 0.0, 0.22, 25.0, 125.0)],
                  metals=[metal((4, 6), (20, 14), [], seed=122, **METAL)]),
        reference=[cut(4, 10, 20, 10), cut(12, 6, 12, 14)])''',
    '''        max_cuts=3,
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, b_z=0.68,
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 0.0, 0.6, 25.0, 125.0)],
                  metals=[metal((2, 4), (22, 16), [], seed=122,
                                **dict(METAL, carriers_per_area=30.0))]),
        reference=[cut(2, 8, 22, 8), cut(2, 12, 22, 12)])''')
sub('''the narrower the pieces, the more. You have two cuts: make the core take at most "
        "0.22 of heat per unit time (averaged over $t = 25$ to $125$).",''',
    '''the narrower the pieces, the more. Drag across the core to cut it, any length and "
        "angle; you have three cuts. Make the core take at most 0.6 of heat per unit time "
        "(averaged over $t = 25$ to $125$).",''')
open(p, "w", encoding="utf-8", newline="\n").write(s)

import json
import os
S = os.path.join("target", "scratch")


def cut(x0, y0, x1, y1):
    import math
    dx, dy = x1 - x0, y1 - y0
    e = {"node": [x0, y0, 0], "kind": "cut", "value": math.hypot(dx, dy)}
    a = math.degrees(math.atan2(dy, dx))
    if a:
        e["angle_deg"] = a
    return e


layouts = {
    "uncut": [],
    "along1": [cut(2, 10, 22, 10)],
    "along2": [cut(2, 8, 22, 8), cut(2, 12, 22, 12)],
    "along3": [cut(2, 7, 22, 7), cut(2, 10, 22, 10), cut(2, 13, 22, 13)],
    "across1": [cut(12, 4, 12, 16)],
    "along_across": [cut(2, 10, 22, 10), cut(12, 4, 12, 16)],
    "staggered": [cut(2, 8, 18, 8), cut(6, 12, 22, 12)],
    "half_along2": [cut(2, 8, 12, 8), cut(2, 12, 12, 12)],
    "diagonal": [cut(2, 4, 22, 16)],
}
json.dump(layouts, open(os.path.join(S, "lam20_layouts.json"), "w", encoding="utf-8"), indent=1)
json.dump(list(layouts.values()), open(os.path.join(S, "lam20_placements.json"), "w",
                                       encoding="utf-8"))
print("ok", list(layouts))
