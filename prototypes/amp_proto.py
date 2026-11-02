"""Amplifier prototype (transistors arc): Gate's channel (length L), a sine of amplitude A and
period T on the gates around the bias (the player's supply), the drain fed from V_DD through
R_L. Writes target/scratch/amp_<tag>.json and amp_<tag>_ref.json (one placement per bias);
the potential goals measure the drain's mean, its amplitude at the signal's frequency and
its phase (degrees against sin wt) over t = settle .. settle + periods T.

    python target/scratch/amp_proto.py tag L R_L V_DD A T periods settle bias1 [bias2 ...]
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
length, rl, vdd, amp, period = (float(a) for a in sys.argv[2:7])
periods, settle = int(sys.argv[7]), float(sys.argv[8])
biases = [float(a) for a in sys.argv[9:]]
omega = 2.0 * math.pi / period
n = int(length)
x0, x1 = 10 - n // 2, 10 + n // 2
start, end = settle, settle + periods * period
lvl = L.level(
    "Amplifier (prototype)", "scratch",
    grid=(20, 16), c=None, t_max=end + 10.0,
    electrodes=[L.contact(x0 - 1, 8, 2), L.contact(x1 + 1, 8, 2),
                L.plate(10, 11, max(n - 2, 2), thickness=1.0, tunable=True),
                L.plate(10, 5, max(n - 2, 2), thickness=1.0)],
    supplies=L.slider(60.0, biases[0]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.CHANNEL_STEP,
                circuits=[L.tube_circuit(1, L.dc(vdd), resistance=rl),
                          L.tube_circuit(2, L.sine(amp, omega), resistance=0.0,
                                         supply_offset=True, linked=[3])],
                potentials=[L.potential_goal(1, "mean", -1e3, 1e3, start, end),
                            L.potential_goal(1, "amplitude", 0.0, 1e3, start, end,
                                             omega=omega),
                            L.potential_goal(1, "phase", -180.0, 180.0, start, end,
                                             omega=omega)],
                metals=[L.metal((x0, 7), (x1, 9), [0, 1], seed=3, **L.CHANNEL)]),
    reference=[L.supply(10, 11, biases[0])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"amp_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(10, 11, b)] for b in biases],
          open(os.path.join(out, f"amp_{tag}_ref.json"), "w"))
print(f"amp_{tag}.json: L {n}, R_L {rl}, V_DD {vdd}, A {amp}, T {period}, window "
      f"{start}-{end}, biases {biases}")
