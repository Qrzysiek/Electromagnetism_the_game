import sys, math
sys.path.insert(0, "target/scratch/inductance")
from eddy2d import ring_power, sigma
side, wall, speed, mfp = 20.0, 2.0, 2.0, 0.5
for density, c in ((60.0, 24.0), (41.0, 27.4)):
    wp = math.sqrt(4 * math.pi * density)
    print(f"== density {density} c {c}: omega_p {wp:.2f} c/omega_p {c / wp:.2f} dt(0.2) {0.2 / wp:.4f}", flush=True)
    for w in [0.1 * (100 ** (k / 16)) for k in range(17)]:
        s = sigma(density, speed, mfp, w)
        p, bc = ring_power(side, wall, s, c, w, 0.1)
        print(f"  w {w:6.3f} c/w {c / w:5.1f}: P {p:9.4g} (B0 = 1)  |B hole| {abs(bc):.3f}", flush=True)
