"""Field-tuned damping (AC arc prototype): Ringing tank's parallel tank with Magnetoresistance's
short wide strip (4 along the current, 10 across) as its resistor, in the player's solenoid:
the left contact fed by a step V at t = 5 through a coil L (no series resistance) with a
capacitor C to ground, the right contact grounded. Measures the left contact's mean and its
amplitude at omega (EM_OMEGA, default 1/sqrt(LC)) over t = start..end. Writes
target/scratch/ftd_<tag>.json and ftd_<tag>_ref.json (one placement per B).

    python target/scratch/ftd_proto.py tag V L C start end t_max B1 [B2 ...]
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.argv.append("--text")  # no solver at import
import levels as L  # noqa: E402

sys.argv.pop()
tag = sys.argv[1]
v, ind, cap, start, end, t_max = (float(a) for a in sys.argv[2:8])
fields = [float(a) for a in sys.argv[8:]]
omega = float(os.environ.get("EM_OMEGA", 1.0 / math.sqrt(ind * cap)))
lvl = L.level(
    "Field-tuned damping (prototype)", "scratch",
    grid=(24, 20), c=None, t_max=t_max,
    electrodes=[L.contact(9, 10, 10), L.contact(15, 10, 10)],
    solenoid_fields=L.slider(16.0, fields[0]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP, solenoid=(23, 19),
                circuits=[L.tube_circuit(0, L.pulse(v, 1.0, delay=5.0), resistance=0.0,
                                         inductance=ind, capacitance=cap)],
                potentials=[L.potential_goal(0, "mean", -1e3, 1e3, start, end),
                            L.potential_goal(0, "amplitude", 0.0, 1e3, start, end,
                                             omega=omega)],
                metals=[L.metal((10, 5), (14, 15), [0, 1], seed=12,
                                **dict(L.METAL, carriers_per_area=float(os.environ.get("EM_CPA", "50"))))]),
    reference=[L.supply(23, 19, fields[0])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"ftd_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(23, 19, b)] for b in fields],
          open(os.path.join(out, f"ftd_{tag}_ref.json"), "w"))
print(f"ftd_{tag}.json: V {v}, L {ind}, C {cap}, omega {omega:.4f} (period "
      f"{2 * math.pi / omega:.1f}), window {start}-{end}, t_max {t_max}, B {fields}")
