import sys
import numpy as np
sys.path.insert(0, "target/scratch/inductance")
from transformer2d import heats
from eddy2d import sigma
density, speed, c = 2.0, 1.0, 4.5
for rings, mfps, label in (([(5.0, 1.0), (3.0, 1.0)], [0.25, 0.25], "outer 10 w1, inner 6 w1"),
                           ([(6.0, 1.0), (4.0, 1.0)], [0.25, 0.25], "outer 12 w1, inner 8 w1")):
    print(f"== {label}: omega_p {np.sqrt(4*np.pi*density):.2f}, c/omega_p {c/np.sqrt(4*np.pi*density):.2f}", flush=True)
    for w in [0.1 * (100 ** (k / 12)) for k in range(13)]:
        sig = [sigma(density, speed, m, w) for m in mfps]
        p, bc = heats(rings, sig, c, w, 0.05)
        print(f"w {w:6.3f} c/w {c / w:5.1f}: outer {p[0]:.4g} inner {p[1]:.4g} ratio {p[1] / p[0]:.4f}  |B centre| {abs(bc):.3f}", flush=True)
