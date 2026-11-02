"""Two concentric round conductors (a disc or annulus inside an annulus) in a uniform swinging
field B0 e^{iwt} along z (sigma(w) for e^{iwt}: Im < 0), z-invariant, Darwin limit: the exact heat of each.

Phi' = 2 pi r B, B' = i w mu0 sigma(r) Phi / (2 pi r) (Faraday round a circle of radius r,
Ampere J = -B'/mu0, Ohm J = sigma E, E = -i w Phi/(2 pi r)); from r = 0 (Phi = 0, B = 1),
scaled so B = B0 outside. Heat per conductor 1/2 Re sigma int |E|^2 2 pi r dr.
Checked against the cylinder's Bessel formula (V-D8's) for one disc.
Usage: python annuli.py"""

import cmath
import math

COS_THETA = -1.0 / 3.0


def mobility(kappa, omega):
    n, top = 28000, 14.0
    h = top / n
    p = q = 0j
    for k in range(n + 1):
        v = h * k
        w = 1.0 if k in (0, n) else (4.0 if k % 2 else 2.0)
        d = kappa * kappa * v * v + omega * omega
        if d == 0.0:
            continue
        f = complex(kappa * kappa * v / d, -kappa * omega / d)
        e = math.exp(-0.5 * v * v)
        p += w * v ** 3 * e / (2 * math.pi) * f
        q += w * v * math.sqrt(2 / math.pi) * v * v * e * f
    p, q = p * h / 3, q * h / 3
    return p / (1 - COS_THETA * q) * (math.pi / kappa)


def heats(regions, omega, mu0, b0, steps_per_unit=4000):
    """regions: [(r_in, r_out, sigma)] sorted, non-overlapping; vacuum between. The heat of
    each region (per unit length along z)."""
    r_end = regions[-1][1]

    def sigma_at(r):
        for a, b, s in regions:
            if a <= r <= b:
                return s
        return 0j

    def deriv(r, phi, bz, s):
        dphi = 2 * math.pi * r * bz
        db = 0j if r == 0 else 1j * omega * mu0 * s * phi / (2 * math.pi * r)
        return dphi, db

    # RK4 over each piece (vacuum and conductors), recording |E|^2 r for the heat.
    edges = sorted({0.0} | {a for a, _, _ in regions} | {b for _, b, _ in regions})
    phi, bz = 0j, 1 + 0j
    out = []
    for a, b in zip(edges, edges[1:]):
        s = sigma_at(0.5 * (a + b))
        n = max(2, int((b - a) * steps_per_unit))
        n += n % 2
        h = (b - a) / n
        r = a
        integrand = []  # |E|^2 2 pi r at each node
        for i in range(n + 1):
            e = 1j * omega * phi / (2 * math.pi * r) if r > 0 else 0j
            integrand.append(abs(e) ** 2 * 2 * math.pi * r)
            if i == n:
                break
            k1 = deriv(r, phi, bz, s)
            k2 = deriv(r + h / 2, phi + h / 2 * k1[0], bz + h / 2 * k1[1], s)
            k3 = deriv(r + h / 2, phi + h / 2 * k2[0], bz + h / 2 * k2[1], s)
            k4 = deriv(r + h, phi + h * k3[0], bz + h * k3[1], s)
            phi += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            bz += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            r = a + h * (i + 1)
        simpson = sum(w * f for w, f in zip(
            [1] + [4 if i % 2 else 2 for i in range(1, n)] + [1], integrand)) * h / 3
        if s != 0:
            out.append(0.5 * s.real * simpson)
    scale = abs(b0 / bz) ** 2
    return [p * scale for p in out], b0 / bz


def bessel_i(nu, z, terms=80):
    return sum((z / 2) ** (2 * m + nu) / (math.factorial(m) * math.factorial(m + nu))
               for m in range(terms))


if __name__ == "__main__":
    # V-D8's cylinder: c = 3, omega 1, B0 0.2, radius 6, mean free path 0.5.
    c, omega, b0, radius = 3.0, 1.0, 0.2, 6.0
    mu0 = 4 * math.pi / (c * c)
    sigma = mobility(1 / 0.5, omega)
    k = cmath.sqrt(1j * omega * mu0 * sigma)
    rim = bessel_i(0, k * radius)
    n = 2000
    h = radius / n
    integral = 0.0
    for i in range(n + 1):
        r = h * i
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        f = (bessel_i(1, k * r) / (k * r) if r > 0 else 0.5) / rim
        integral += w * r ** 3 * abs(f) ** 2
    integral *= h / 3
    exact = math.pi * sigma.real * omega ** 2 * b0 ** 2 * integral
    shot, _ = heats([(0.0, radius, sigma)], omega, mu0, b0)
    print(f"sigma {sigma:.6f}; cylinder Bessel {exact:.6f}, shooting {shot[0]:.6f}, "
          f"ratio {shot[0] / exact:.8f}")


def profile(regions, omega, mu0, b0, n=20000):
    """(r, Phi, B, sigma) on a fine grid by RK4 (as `heats`), scaled to B0 outside."""
    r_end = regions[-1][1]

    def sigma_at(r):
        for a, b, s in regions:
            if a <= r <= b:
                return s
        return 0j
    h = r_end / n
    phi, bz, r = 0j, 1 + 0j, 0.0
    rows = []
    for i in range(n + 1):
        rows.append((r, phi, bz, sigma_at(r)))
        if i == n:
            break
        s = sigma_at(r + h / 2)
        f = lambda rr, p, b: (2 * math.pi * rr * b,
                              0j if rr == 0 else 1j * omega * mu0 * s * p / (2 * math.pi * rr))
        k1 = f(r, phi, bz)
        k2 = f(r + h / 2, phi + h / 2 * k1[0], bz + h / 2 * k1[1])
        k3 = f(r + h / 2, phi + h / 2 * k2[0], bz + h / 2 * k2[1])
        k4 = f(r + h, phi + h * k3[0], bz + h * k3[1])
        phi += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        bz += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        r += h
    scale = b0 / bz
    return [(r, p * scale, b * scale, s) for r, p, b, s in rows]


def report(regions, omega, mu0, b0):
    rows = profile(regions, omega, mu0, b0)
    heat = [0.0] * len(regions)
    ext = [0.0] * len(regions)
    for (r0, p0, _, s0), (r1, p1, _, _) in zip(rows, rows[1:]):
        r = 0.5 * (r0 + r1)
        for k, (a, b, s) in enumerate(regions):
            if a <= r <= b:
                phi = 0.5 * (p0 + p1)
                e = -1j * omega * phi / (2 * math.pi * r)
                j = s * e
                e_ext = -1j * omega * b0 * r / 2
                dA = 2 * math.pi * r * (r1 - r0)
                heat[k] += 0.5 * (e * j.conjugate()).real * dA
                ext[k] += 0.5 * (e_ext * j.conjugate()).real * dA
    return heat, ext


def scan():
    c, omega, b0 = 3.0, 1.0, 0.2
    mu0 = 4 * math.pi / (c * c)
    sigma = mobility(1 / 0.5, omega)
    for prim, sec in [((0.0, 3.0), (4.0, 6.0)), ((0.0, 3.0), (4.5, 6.0)), ((2.0, 3.5), (4.5, 6.0)),
                      ((0.0, 2.5), (3.5, 6.0))]:
        both, ext = report([(prim[0], prim[1], sigma), (sec[0], sec[1], sigma)], omega, mu0, b0)
        alone, _ = report([(prim[0], prim[1], sigma)], omega, mu0, b0)
        sec_alone, _ = report([(sec[0], sec[1], sigma)], omega, mu0, b0)
        print(f"primary {prim} secondary {sec}: heats {both[0]:.4f} {both[1]:.4f}; drive's work "
              f"{ext[0]:.4f} {ext[1]:.4f}; primary alone {alone[0]:.4f}, secondary alone {sec_alone[0]:.4f}")


if __name__ == "__main__":
    scan()
