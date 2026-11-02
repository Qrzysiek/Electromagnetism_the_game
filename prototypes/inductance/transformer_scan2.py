import sys
sys.path.insert(0, "target/scratch/inductance")
from transformer2d import scan
omegas = [0.03 * (100 ** (k / 10)) for k in range(11)]
scan("v1: outer 24 wall 4, inner 14 wall 2, density 1, speed 1, mfp 1, c 3.2",
     [(12.0, 4.0), (7.0, 2.0)], 1.0, 1.0, 1.0, 3.2, omegas)
scan("v2: outer 24 wall 3, inner 16 wall 2, density 1, speed 1, mfp 0.75, c 3.2",
     [(12.0, 3.0), (8.0, 2.0)], 1.0, 1.0, 0.75, 3.2, omegas)
