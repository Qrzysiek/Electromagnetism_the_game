"""A 20 x 12 core: local Ohm's law for full-length cut layouts, relative to the uncut core."""
import os

S = os.path.dirname(os.path.abspath(__file__))
ns = {}
exec(open(os.path.join(S, "eddy_theory.py"), encoding="utf-8").read().split("sigma = 0.47")[0], ns)
integral = ns["integral"]


def I(w, b):
    return integral(min(w, b), max(w, b))


uncut = I(20, 12)
for name, v in [
    ("uncut", uncut),
    ("1 along (y = 10)", 2 * I(20, 6)),
    ("2 along (y = 8, 12)", 3 * I(20, 4)),
    ("3 along (y = 7, 10, 13)", 4 * I(20, 3)),
    ("1 across (x = 12)", 2 * I(10, 12)),
    ("2 across (x = 9, 16)", 2 * I(7, 12) + I(6, 12)),
    ("1 along + 1 across", 4 * I(10, 6)),
    ("2 along + 1 across", 6 * I(10, 4)),
]:
    print(f"{name:26s} {v / uncut:.3f}")
b2 = (0.68 * 0.125) ** 2 / 2
print("uncut heat at B0 = 0.68, omega 0.125:", 0.47 * b2 * uncut)
