"""Prototype levels for the thermal-noise design: a strip, a grounded contact at its left end,
a floating contact (a lumped capacitance to ground) at its right end, optionally a grounded
plate beside the floating contact. Written from the Hall probe's file.

    python proto.py NAME L W C END [--plate GAP LENGTH] [--seed S] [--cpa N]
"""
import argparse
import copy
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


def build(length, width, cap, end, plate=None, seed=1, cpa=100.0, mfp=1.0, kick=None):
    base = json.load(open(os.path.join(ROOT, "levels", "112_hall_probe.json"), encoding="utf-8"))
    lvl = copy.deepcopy(base)
    x0, y0 = 6, 7
    length, width = int(length), int(width)
    ym = y0 + width // 2
    x1 = x0 + length
    box = {"thickness": 2.0, "height": 4.0, "angle_deg": 90.0, "bias": {"kind": "grounded"}}
    electrodes = [
        dict(box, center=[x0 - 1, ym, 0], length=width, label="grounded end"),
        dict(box, center=[x1 + 1, ym, 0], length=width, label="floating contact"),
    ]
    if plate:
        gap, plen = plate
        electrodes.append(dict(box, center=[int(x1 + 2 + gap + 0.5), ym, 0], length=plen, thickness=1.0,
                               label="grounded plate"))
    if kick:
        volts, at = kick
        electrodes.append(dict(box, center=[x1 + 1, y0 + width + 1, 0], length=2.0, thickness=1.0,
                               angle_deg=0.0, label="kicker"))
    lvl["electrodes"] = electrodes
    lvl["name"] = "noise prototype"
    lvl["reference_solution"] = []
    lvl["limits"]["solenoid_fields"] = []
    t = lvl["tube"]
    t["solenoid"] = None
    t["metals"] = [dict(t["metals"][0], corners=[[x0, y0, 0], [x1, y0 + width, 0]], contacts=[0, 1],
                        seed=seed, carriers_per_area=cpa, mean_free_path=mfp, label="strip")]
    t["circuits"] = [{"electrode": 1, "source": {"kind": "dc", "value": 0.0}, "capacitance": cap}]
    if kick:
        t["circuits"].append({"electrode": len(electrodes) - 1, "resistance": 0.0,
                              "source": {"kind": "pulse", "low": 0.0, "high": kick[0], "delay": kick[1],
                                         "rise": 0.2, "width": 100000.0, "fall": 0.0}})
    t["goals"] = [{"electrode": 0, "min": -100.0, "max": 100.0, "start": 10.0, "end": end,
                   "carriers": "net"}]
    t["potential_goals"] = [{"electrode": 1, "measure": "mean", "min": -10.0, "max": 10.0,
                             "start": 10.0, "end": end}]
    return lvl


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("length", type=float)
    ap.add_argument("width", type=float)
    ap.add_argument("cap", type=float)
    ap.add_argument("end", type=float)
    ap.add_argument("--plate", type=float, nargs=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpa", type=float, default=100.0)
    ap.add_argument("--mfp", type=float, default=1.0)
    ap.add_argument("--kick", type=float, nargs=2, help="VOLTS AT: a step on an electrode beside the contact")
    a = ap.parse_args()
    lvl = build(a.length, a.width, a.cap, a.end, a.plate, a.seed, a.cpa, a.mfp, a.kick)
    path = os.path.join(HERE, a.name + ".json")
    json.dump(lvl, open(path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(path)
