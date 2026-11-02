"""The inductance of two antiparallel straight leads (current sheets along z in the
z-invariant world), length l, d apart, in the infinite plane, in the Darwin (Coulomb-gauge)
field: L = (μ₀/2π) ∫∫ [½ ln(1 + d²/Δ²) + d²/(Δ² + d²)] dy dy', Δ = y − y'
   = (μ₀/2π) · 2 ∫₀ˡ (l − u) f(u) du,
from W = (μ₀/2)∫ψ², ψ the stream function, through the biharmonic Green's function
r² ln r/(8π) of (−∇²)⁻² acting on the leads' curl (a dipole layer along each). Long leads:
μ₀ d l. Evaluated with mpmath to 30 digits; also checked against ∫ψ² summed on a fine grid
from the angle each lead subtends (ψ = (I/2π)(Ω₁ − Ω₂))."""
import math
import sys

import mpmath as mp

mp.mp.dps = 30


def loop(d, l, mu0):
    f = lambda u: 0.5 * mp.log(1 + d * d / (u * u)) + d * d / (u * u + d * d)
    integral = 2 * mp.quad(lambda u: (l - u) * f(u), [0, d, l])
    return mu0 / (2 * mp.pi) * integral


def grid_check(d, l, h=0.01, extent=60.0):
    """∫ψ² dA (ψ per unit current) by the midpoint rule over a box, the leads at x = 0, d
    from y = 0 to l, the box `extent` beyond them; the far tail ~1/r⁴ neglected."""
    import numpy as np

    xs = np.arange(-extent, d + extent, h) + 0.5 * h
    ys = np.arange(-extent, l + extent, h) + 0.5 * h
    total = 0.0
    for y in ys:
        def omega(x0):
            ax, ay = x0 - xs, 0.0 - y
            bx, by = x0 - xs, l - y
            return np.arctan2(ax * by - ay * bx, ax * bx + ay * by)

        psi = (omega(0.0) - omega(d)) / (2 * math.pi)
        total += float(np.sum(psi * psi)) * h * h
    return total


if __name__ == "__main__":
    c = 3.0
    mu0 = 4 * math.pi / c**2
    for d, l in [(4.0, 9.5), (2.0, 24.0), (4.0, 24.0)]:
        L = loop(d, l, mu0)
        print(f"d {d}, l {l}: L = {mp.nstr(L, 12)} (μ₀ d l = {mu0 * d * l:.6f}; ratio {float(L) / (mu0 * d * l):.6f})")
    if "--grid" in sys.argv:
        d, l = 4.0, 9.5
        s = grid_check(d, l, h=0.02, extent=40.0)
        print(f"grid check d {d}, l {l}: μ₀∫ψ² = {mu0 * s:.6f}")
