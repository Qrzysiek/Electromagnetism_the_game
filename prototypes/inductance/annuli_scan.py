import math, sys
sys.path.insert(0, "target/scratch/inductance")
from annuli import mobility, report
c, b0 = 3.0, 0.2
mu0 = 4 * math.pi / (c * c)
for omega in (0.5, 1.0):
    sigma = mobility(1 / 0.5, omega)
    for prim, sec in [((2.0, 3.5), (4.5, 6.0)), ((2.5, 4.5), (5.0, 6.0)), ((3.0, 4.5), (5.0, 6.0)),
                      ((2.0, 4.0), (5.0, 6.0)), ((3.0, 4.0), (4.5, 6.0)), ((3.5, 4.5), (5.0, 6.0))]:
        both, ext = report([(prim[0], prim[1], sigma), (sec[0], sec[1], sigma)], omega, mu0, b0)
        print(f"w {omega} inner {prim} outer {sec}: heats {both[0]:.4f} {both[1]:.4f}; drive's work "
              f"{ext[0]:.4f} {ext[1]:.4f}; ratio inner {both[0]/ext[0]:.2f} outer {both[1]/ext[1]:.3f}", flush=True)
