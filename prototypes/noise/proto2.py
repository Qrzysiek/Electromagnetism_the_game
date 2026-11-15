"""The thermal-noise level's prototype (DEC-42): two strips grounded at their left ends, each
with a floating contact at its right end. A is the reference; B, the player's, has a kicker
above its contact that steps at `kick_at`. The player may place a grounded plate and make a
cut. Goals: B's noise's spread over A's over `10..kick_at`, B's decay over `kick_at..end`.

    python proto2.py NAME [--length L] [--width W] [--cap C] [--kick V] [--at T] [--end T]
                          [--cpa N] [--sub S] [--seed S]
"""
import argparse
import copy
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


def build(length=4, width=4, cap=0.0, kick=40.0, at=400.0, end=420.0, cpa=25.0, sub=1, seed=1,
          ratio=(0.5, 0.7), decay=(2.0, 3.0), square=None):
    base = json.load(open(os.path.join(ROOT, "levels", "112_hall_probe.json"), encoding="utf-8"))
    lvl = copy.deepcopy(base)
    lvl["name"] = "noise prototype 2"
    lvl["grid"] = {"nx": 20, "ny": 16, "nz": 0, "subdivision": sub}
    box = {"thickness": 2.0, "height": 4.0, "angle_deg": 90.0, "bias": {"kind": "grounded"}}
    x0 = 4
    x1 = x0 + length
    ya, yb = 2, 9  # the strips' lower edges
    electrodes, metals = [], []
    for k, y0 in enumerate((ya, yb)):
        ym = y0 + width // 2
        name = "AB"[k]
        electrodes += [
            dict(box, center=[x0 - 1, ym, 0], length=width, label=f"strip {name}'s grounded end"),
            dict(box, center=[x1 + 1, ym, 0], length=width, label=f"contact {name} (floating)"),
        ]
        metal = dict(base["tube"]["metals"][0], corners=[[x0, y0, 0], [x1, y0 + width, 0]],
                     contacts=[2 * k, 2 * k + 1], seed=seed + k, carriers_per_area=cpa,
                     label=f"strip {name}")
        metals.append(metal)
    electrodes.append(dict(box, center=[x1 + 1, yb + width + 1, 0], length=2.0, thickness=1.0,
                           angle_deg=0.0, label="kicker"))
    lvl["electrodes"] = electrodes
    lvl["reference_solution"] = []
    lvl["limits"] = {"max_charges": 0, "magnitudes": [], "allow_positive": True,
                     "allow_negative": True, "max_magnets": 0, "magnet_strengths": [],
                     "max_plates": 1, "plate_voltages": [0.0],
                     "plate": {"length": 4.0, "thickness": 0.5, "height": 4.0}, "max_cuts": 1}
    t = lvl["tube"]
    t["solenoid"] = None
    t["metals"] = metals
    t["circuits"] = [
        {"electrode": 1, "source": {"kind": "dc", "value": 0.0}, "capacitance": cap},
        {"electrode": 3, "source": {"kind": "dc", "value": 0.0}, "capacitance": cap},
        {"electrode": 4, "resistance": 0.0,
         "source": ({"kind": "pulse", "low": 0.0, "high": kick, "delay": at, "rise": 0.2,
                     "width": 100000.0, "fall": 0.0} if square is None else
                    {"kind": "pulse", "low": -kick, "high": kick, "delay": 0.0, "rise": 0.1,
                     "width": square / 2 - 0.1, "fall": 0.1, "period": square})},
    ]
    t["goals"] = []
    t["potential_goals"] = [
        {"electrode": 3, "reference": 1, "measure": "spread_ratio", "min": ratio[0],
         "max": ratio[1], "start": 10.0, "end": at},
        {"electrode": 3, "measure": "decay", "min": decay[0], "max": decay[1], "start": at,
         "end": end},
    ]
    return lvl


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--length", type=int, default=4)
    ap.add_argument("--width", type=int, default=4)
    ap.add_argument("--cap", type=float, default=0.0)
    ap.add_argument("--kick", type=float, default=40.0)
    ap.add_argument("--at", type=float, default=400.0)
    ap.add_argument("--end", type=float, default=420.0)
    ap.add_argument("--cpa", type=float, default=25.0)
    ap.add_argument("--sub", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--square", type=float, help="a square wave of this period, +-kick")
    a = ap.parse_args()
    lvl = build(a.length, a.width, a.cap, a.kick, a.at, a.end, a.cpa, a.sub, a.seed, square=a.square)
    path = os.path.join(HERE, a.name + ".json")
    json.dump(lvl, open(path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(path)
