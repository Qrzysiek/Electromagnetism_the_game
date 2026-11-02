"""Inverter prototype (transistors arc): the Gate channel with its drain fed from V_DD through a
load resistor; the gate held at -40 (off) until a pulse brings it to the player's level at
t = 60 for 70. The drain's potential is the output: high while the channel is off, low while
it conducts.

    python target/scratch/fet_inverter.py [V_DD] [R_L] [gate_high]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = open(os.path.join(ROOT, "scripts", "levels.py"), encoding="utf-8").read()
g = {"__file__": os.path.join(ROOT, "scripts", "levels.py"), "__name__": "levels_scratch"}
exec(compile(src, "levels_scratch", "exec"), g)
L = g
vdd, rl, high = (float(a) for a in (sys.argv[1:4] if len(sys.argv) >= 4 else ("20", "30", "0")))
lvl = L["level"](
    "Inverter (prototype)", "scratch",
    grid=(20, 16), c=None, t_max=200.0,
    # 0: the source (grounded); 1: the drain, fed from V_DD through R_L; 2, 3: the gates.
    electrodes=[L["contact"](5, 8, 2), L["contact"](15, 8, 2, tunable=True),
                L["plate"](10, 11, 6, thickness=1.0, tunable=True),
                L["plate"](10, 5, 6, thickness=1.0)],
    supplies=L["slider"](60.0, vdd, high),
    tube=L["tube"](0, None, 0, 0, 0, 0, step=L["CHANNEL_STEP"],
                   circuits=[L["tube_circuit"](1, L["dc"](vdd), resistance=rl),
                             L["tube_circuit"](2, L["pulse"](high, 2.0, delay=60.0, width=70.0,
                                                             fall=2.0, low=-40.0),
                                               resistance=0.0, linked=[3])],
                   goals=[L["current_goal"](1, 0.0, 10.0, 90.0, 130.0, carriers="net")],
                   metals=[L["metal"]((6, 7), (14, 9), [0, 1], seed=3, **L["CHANNEL"])]),
    reference=[L["supply"](15, 8, vdd), L["supply"](10, 11, high)])
lvl = L["resolve_auto"](lvl)
json.dump(lvl, open(os.path.join(ROOT, "target", "scratch", "inverter.json"), "w"), indent=1)
json.dump([lvl["reference_solution"]],
          open(os.path.join(ROOT, "target", "scratch", "inverter_ref.json"), "w"))
print("ok")
