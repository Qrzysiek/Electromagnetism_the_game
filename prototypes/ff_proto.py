"""Flip-flop prototype (transistors arc): two of Gate's transistors, A (channel y 3-5) and B
(y 13-15), each drain fed from V_DD through R_L, each pair of gates wired to the other's drain
through a level-shifting source -S (the player's supplies: S_A on A's gates, S_B on B's) and
a resistance r. Reset: drain A's supply drops to 0 for t = 0..44 (B's gates go low, B off, A
on); set: drain B's supply drops for t = 150..184 (A off, B on). Goals (measures): each drain's
mean over t = 100..150 (reset held) and 250..300 (set held).
Writes target/scratch/ff_<tag>.json and ff_<tag>_ref.json (one placement per (S_A, S_B)).

    python target/scratch/ff_proto.py tag V_DD R_L r t_max SA,SB [SA,SB ...]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.argv.append("--text")  # no solver at import
import levels as L  # noqa: E402

sys.argv.pop()
tag = sys.argv[1]
vdd, rl, r, t_max = (float(a) for a in sys.argv[2:6])
pairs = [tuple(float(v) for v in p.split(",")) for p in sys.argv[6:]]
reset = L.pulse(0.0, 2.0, delay=0.0, width=40.0, fall=2.0, low=vdd)
setp = L.pulse(0.0, 2.0, delay=150.0, width=30.0, fall=2.0, low=vdd)
lvl = L.level(
    "Flip-flop (prototype)", "scratch",
    grid=(20, 20), c=None, t_max=t_max,
    # 0, 1: A's source and drain; 2, 3: A's gates; 4, 5: B's source and drain; 6, 7: B's
    # gates.
    electrodes=[L.contact(5, 4, 2), L.contact(15, 4, 2),
                L.plate(10, 7, 6, thickness=1.0, tunable=True), L.plate(10, 1, 6, thickness=1.0),
                L.contact(5, 14, 2), L.contact(15, 14, 2),
                L.plate(10, 17, 6, thickness=1.0, tunable=True),
                L.plate(10, 11, 6, thickness=1.0)],
    supplies=L.slider(240.0, -pairs[0][0], -pairs[0][1]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.CHANNEL_STEP,
                circuits=[L.tube_circuit(1, reset, resistance=rl),
                          L.tube_circuit(5, setp, resistance=rl),
                          L.tube_circuit(2, L.dc(0.0), resistance=r, linked=[3], return_to=5),
                          L.tube_circuit(6, L.dc(0.0), resistance=r, linked=[7], return_to=1)],
                potentials=[L.potential_goal(1, "mean", -1e3, 1e3, 100.0, 150.0),
                            L.potential_goal(5, "mean", -1e3, 1e3, 100.0, 150.0),
                            L.potential_goal(1, "mean", -1e3, 1e3, 250.0, 300.0),
                            L.potential_goal(5, "mean", -1e3, 1e3, 250.0, 300.0)],
                metals=[L.metal((6, 3), (14, 5), [0, 1], seed=3, **L.CHANNEL),
                        L.metal((6, 13), (14, 15), [4, 5], seed=5, **L.CHANNEL)]),
    reference=[L.supply(10, 7, -pairs[0][0]), L.supply(10, 17, -pairs[0][1])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"ff_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(10, 7, -a), L.supply(10, 17, -b)] for a, b in pairs],
          open(os.path.join(out, f"ff_{tag}_ref.json"), "w"))
print(f"ff_{tag}.json: V_DD {vdd}, R_L {rl}, r {r}, t_max {t_max}, (S_A, S_B) {pairs}")
