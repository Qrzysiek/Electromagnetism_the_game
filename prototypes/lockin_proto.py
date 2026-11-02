"""Lock-in Hall (AC arc prototype): Hall probe's strip (16 x 4 at the board's centre, 100
carriers per area) with its right end on the player's supply and floating side contacts at
mid-length, in the player's swinging solenoid B0 sin(omega t) (EDDY_OMEGA). Measures, over two
whole periods from t = 10: the top contact against the bottom one, its amplitude and phase at
omega (degrees against sin omega t), and the current into the right end. Writes
target/scratch/li_<tag>.json and li_<tag>_ref.json (one placement per (V, B0)).

    python target/scratch/lockin_proto.py tag V,B0 [V,B0 ...]
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
pairs = [tuple(float(x) for x in p.split(",")) for p in sys.argv[2:]]
omega = L.EDDY_OMEGA
start = 10.0
end = start + 2 * 2 * math.pi / omega
hall = dict(L.METAL, carriers_per_area=100.0)
lvl = L.level(
    "Lock-in Hall (prototype)", "scratch",
    grid=(24, 16), c=None, t_max=end + 5.0,
    electrodes=[L.contact(3, 8, 4), L.contact(21, 8, 4, tunable=True),
                L.contact(12, 11, 2, vertical=False), L.contact(12, 5, 2, vertical=False)],
    supplies=L.slider(40.0, pairs[0][0]),
    solenoid_fields=L.slider(16.0, pairs[0][1]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP, solenoid=(23, 15),
                solenoid_drive=L.sine(1.0, omega),
                circuits=[L.tube_circuit(2, L.dc(0.0), capacitance=1.0),
                          L.tube_circuit(3, L.dc(0.0), capacitance=1.0)],
                goals=[L.current_goal(1, -1e3, 1e3, start, end, carriers="net")],
                potentials=[L.potential_goal(2, "amplitude", 0.0, 1e3, start, end, omega=omega,
                                             reference=3),
                            L.potential_goal(2, "phase", -180.0, 180.0, start, end,
                                             omega=omega, reference=3),
                            L.potential_goal(2, "mean", -1e3, 1e3, start, end, reference=3)],
                metals=[L.metal((4, 6), (20, 10), [0, 1, 2, 3], seed=11, **hall)]),
    reference=[L.supply(21, 8, pairs[0][0]), L.supply(23, 15, pairs[0][1])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"li_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(21, 8, v), L.supply(23, 15, b)] for v, b in pairs],
          open(os.path.join(out, f"li_{tag}_ref.json"), "w"))
print(f"li_{tag}.json: omega {omega} (period {2 * math.pi / omega:.2f}), window {start}-{end:.2f}, "
      f"(V, B0) {pairs}")
