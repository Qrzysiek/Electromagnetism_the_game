"""RLC with a strip (AC arc prototype): a strip's left contact fed by a step (0 to V at t = 5,
rising over 1) through a small resistance R_s and an inductance L, with a capacitance C to
ground; the strip's right contact grounded. The left end rings at ~1/sqrt(LC), damped by the
strip (a parallel RLC: Q = R sqrt(C/L)). Goals measured: the left end's mean and its
oscillation amplitude at omega0 over t = start..end. Writes target/scratch/rlc_<tag>.json and
rlc_<tag>_ref.json (one placement per V).

    python target/scratch/rlc_proto.py tag strip_len L C R_s start end t_max V1 [V2 ...]
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
n = int(sys.argv[2])
ind, cap, rs, start, end, t_max = (float(a) for a in sys.argv[3:9])
volts = [float(a) for a in sys.argv[9:]]
omega = float(os.environ.get("EM_OMEGA", 1.0 / math.sqrt(ind * cap)))
x0, x1 = 12 - n // 2, 12 + n // 2
lvl = L.level(
    "RLC (prototype)", "scratch",
    grid=(24, 20), c=None, t_max=t_max,
    electrodes=[L.contact(x0 - 1, 10, 4, tunable=True), L.contact(x1 + 1, 10, 4)],
    supplies=L.slider(40.0, volts[0]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP,
                circuits=[L.tube_circuit(0, L.pulse(volts[0], 1.0, delay=5.0), resistance=rs,
                                         inductance=ind, capacitance=cap)],
                potentials=[L.potential_goal(0, "mean", -1e3, 1e3, start, end),
                            L.potential_goal(0, "amplitude", 0.0, 1e3, start, end,
                                             omega=omega)],
                metals=[L.metal((x0, 8), (x1, 12), [0, 1], seed=11, **L.METAL)]),
    reference=[L.supply(x0 - 1, 10, volts[0])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"rlc_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(x0 - 1, 10, v)] for v in volts],
          open(os.path.join(out, f"rlc_{tag}_ref.json"), "w"))
print(f"rlc_{tag}.json: strip {n} x 4, L {ind}, C {cap}, R_s {rs}, omega0 {omega:.4f} "
      f"(period {2 * math.pi / omega:.1f}), window {start}-{end}, t_max {t_max}, V {volts}")
