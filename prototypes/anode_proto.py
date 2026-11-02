"""Resistive anode (combination prototype): a temperature-limited cathode (a small plate at the
bottom, at -V_K) shoots electrons up past a pair of deflection plates (the left one on the
player's supply, the right one grounded) onto a metal strip 16 x 4 (y 13-17) whose two end
contacts are grounded; the beam's current divides between the ends as the strip's resistances
on either side of the landing point (charge division). Measures the net currents into both
end contacts over t = start..end. Writes target/scratch/an_<tag>.json and an_<tag>_ref.json
(one placement per supply).

    python target/scratch/anode_proto.py tag V_K saturation start end t_max V1 [V2 ...]
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
vk, sat, start, end, t_max = (float(a) for a in sys.argv[2:7])
volts = [float(a) for a in sys.argv[7:]]
lvl = L.level(
    "Resistive anode (prototype)", "scratch",
    grid=(24, 20), c=None, t_max=t_max,
    # 0: the cathode; 1, 2: the strip's left and right end contacts; 3, 4: the deflection
    # plates (left: your supply; right: grounded).
    electrodes=[L.plate(12, 3, 1, thickness=0.5, kind="potential", value=-vk),
                L.contact(3, 15, 4), L.contact(21, 15, 4),
                L.plate(8, 8, 4, thickness=0.5, angle_deg=90.0, tunable=True),
                L.plate(16, 8, 4, thickness=0.5, angle_deg=90.0)],
    supplies=L.slider(40.0, volts[0]),
    tube=L.tube(0, None, 0, 0, 0, 0, step=L.METAL_STEP, saturation=sat,
                emit_toward_deg=90.0,
                goals=[L.current_goal(1, -1e3, 1e3, start, end, carriers="net"),
                       L.current_goal(2, -1e3, 1e3, start, end, carriers="net")],
                metals=[L.metal((4, 13), (20, 17), [1, 2], seed=11, **L.METAL)]),
    reference=[L.supply(8, 8, volts[0])])
out = os.path.join(ROOT, "target", "scratch")
json.dump(lvl, open(os.path.join(out, f"an_{tag}.json"), "w"), indent=1)
json.dump([[L.supply(8, 8, v)] for v in volts],
          open(os.path.join(out, f"an_{tag}_ref.json"), "w"))
print(f"an_{tag}.json: V_K {vk}, saturation {sat}, window {start}-{end}, t_max {t_max}, "
      f"supplies {volts}")
