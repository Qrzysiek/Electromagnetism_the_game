"""The small Transformer candidate: exact heats of the two rings against omega."""
import sys
import numpy as np
sys.path.insert(0, "target/scratch/inductance")
from transformer2d import heats
from eddy2d import sigma

density, speed, c = 2.0, 1.0, 4.5
rings = [(5.0, 2.0), (2.0, 1.0)]  # outer 10 wide wall 2, inner 4 wide wall 1 (half sides)
mfps = [0.5, 0.25]
print(f"omega_p {np.sqrt(4 * np.pi * density):.2f}, c/omega_p {c / np.sqrt(4 * np.pi * density):.2f}")
for w in [0.1 * (100 ** (k / 12)) for k in range(13)]:
    sig = [sigma(density, speed, m, w) for m in mfps]
    p, bc = heats(rings, sig, c, w, 0.05)
    a0 = heats([rings[0]], [sig[0]], c, w, 0.05)[0][0]
    a1 = heats([rings[1]], [sig[1]], c, w, 0.05)[0][0]
    print(f"w {w:6.3f} c/w {c / w:5.1f}: outer {p[0]:.4g} inner {p[1]:.4g} ratio {p[1] / p[0]:.4f}"
          f"  alone {a0:.4g} {a1:.4g}  |B centre| {abs(bc):.3f}", flush=True)
