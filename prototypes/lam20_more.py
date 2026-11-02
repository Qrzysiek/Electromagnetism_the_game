"""Laminations 20 x 12: more three-cut layouts (is the solving set wide?), verification."""
import json
import math
import os

S = os.path.dirname(os.path.abspath(__file__))


def cut(x0, y0, x1, y1):
    dx, dy = x1 - x0, y1 - y0
    e = {"node": [x0, y0, 0], "kind": "cut", "value": math.hypot(dx, dy)}
    a = math.degrees(math.atan2(dy, dx))
    if a:
        e["angle_deg"] = a
    return e


layouts = {
    "along2_across": [cut(2, 8, 22, 8), cut(2, 12, 22, 12), cut(12, 4, 12, 16)],
    "along3_uneven": [cut(2, 6, 22, 6), cut(2, 10, 22, 10), cut(2, 14, 22, 14)],
    "along3_offset": [cut(2, 7, 22, 7), cut(2, 11, 22, 11), cut(2, 14, 22, 14)],
}
json.dump(layouts, open(os.path.join(S, "lam20_more_layouts.json"), "w", encoding="utf-8"), indent=1)
json.dump(list(layouts.values()), open(os.path.join(S, "lam20_more.json"), "w", encoding="utf-8"))
print("ok")
