"""The exact magnetoquasistatic answer for a z-invariant square ring in a uniform axial AC
field (an independent check of the game's Darwin self-field): B_z obeys
div(D grad B) = i w B with the magnetic diffusivity D = 1/(mu0 sigma(w)) in the wall; in the
hole (sigma = 0) B is uniform, which a huge D enforces (the flux through the hole's edge is
then Faraday's law round it); B = B_ext = 1 on the ring's outer edge (z-invariance: the
ring's own field vanishes outside). Time-averaged power per unit length
P = (1/2) Re(1/sigma) int |J|^2, J = (1/mu0) curl(B z).

Usage: py -3.14 eddy2d.py  (prints the ratio P(c)/P(c_ref) at each omega)."""

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

COS_THETA = -1.0 / 3.0


def mobility(kappa, omega):
    """The Boltzmann mobility of the game's metal (thermal speed 1; PHYSICS.md 2.11),
    crates/physics/tests/common/mod.rs `mobility`."""
    n, top = 28000, 14.0
    v = np.linspace(0.0, top, n + 1)
    w = np.ones(n + 1)
    w[1:-1:2], w[2:-1:2] = 4.0, 2.0
    d = kappa * kappa * v * v + omega * omega
    k = kappa / (kappa * v + 1j * omega) if omega != 0 else np.where(v > 0, 1 / np.maximum(v, 1e-300), 0)
    e = np.exp(-0.5 * v * v)
    fp = v ** 3 * e / (2 * np.pi)
    fq = v * np.sqrt(2 / np.pi) * v * v * e
    h = top / n
    p = np.sum(w * fp * k) * h / 3
    q = np.sum(w * fq * k) * h / 3
    return p / (1 - COS_THETA * q) * (np.pi / kappa)


def sigma(density, speed, mfp, omega, ac=True):
    """The metal's conductivity at the drive's frequency: the thermal-speed-1 mobility
    rescaled (times scale as 1/speed), |q/m| = 1."""
    kappa = 1.0 / mfp
    if not ac:
        omega = 0.0
    return density * mobility(kappa, omega / speed) / speed


def ring_power(side, wall, sig, c, omega, h=0.05):
    mu0 = 4 * np.pi / c ** 2
    n = int(round(side / h))  # cells per side; nodes 0..n, boundary nodes fixed
    xs = (np.arange(n + 1) * h) - side / 2
    # Cell (edge midpoint) diffusivity: wall or hole.
    inner = side / 2 - wall

    def in_hole(x, y):
        return (abs(x) < inner) & (abs(y) < inner)

    d_wall = 1.0 / (mu0 * sig)
    d_hole = 1e6 * abs(d_wall)
    m = n - 1  # interior nodes per side
    idx = lambda i, j: (i - 1) + m * (j - 1)
    rows, cols, vals = [], [], []
    rhs = np.zeros(m * m, dtype=complex)
    for j in range(1, n):
        for i in range(1, n):
            k = idx(i, j)
            diag = -1j * omega * h * h
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xm = xs[i] + 0.5 * di * h
                ym = xs[j] + 0.5 * dj * h
                dd = d_hole if in_hole(xm, ym) else d_wall
                diag -= dd
                ii, jj = i + di, j + dj
                if 1 <= ii <= n - 1 and 1 <= jj <= n - 1:
                    rows.append(k)
                    cols.append(idx(ii, jj))
                    vals.append(dd)
                else:
                    rhs[k] -= dd * 1.0  # B = 1 on the outer edge
            rows.append(k)
            cols.append(k)
            vals.append(diag)
    a = sp.csr_matrix((vals, (rows, cols)), shape=(m * m, m * m))
    b = spla.spsolve(a.tocsc(), rhs)
    full = np.ones((n + 1, n + 1), dtype=complex)
    full[1:n, 1:n] = b.reshape(m, m).T  # full[i, j]
    # Power: on each wall edge, J component from the difference across it.
    p = 0.0
    rho = (1 / sig).real
    for i in range(n):
        for j in range(n + 1):
            xm, ym = xs[i] + 0.5 * h, xs[j]
            if not in_hole(xm, ym) and abs(ym) <= side / 2:
                jy = -(full[i + 1, j] - full[i, j]) / h / mu0
                wgt = 0.5 if j in (0, n) else 1.0
                p += wgt * 0.5 * rho * abs(jy) ** 2 * h * h
    for i in range(n + 1):
        for j in range(n):
            xm, ym = xs[i], xs[j] + 0.5 * h
            if not in_hole(xm, ym):
                jx = (full[i, j + 1] - full[i, j]) / h / mu0
                wgt = 0.5 if i in (0, n) else 1.0
                p += wgt * 0.5 * rho * abs(jx) ** 2 * h * h
    centre = full[n // 2, n // 2]
    return p, centre


if __name__ == "__main__":
    side, wall, density, speed, mfp, c = 20.0, 2.0, 20.0, 2.0, 0.5, 7.3
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
    for omega in (0.1, 0.2, 0.4, 0.8, 1.6, 3.2):
        for ac in (False, True):
            s = sigma(density, speed, mfp, omega, ac)
            p, bc = ring_power(side, wall, s, c, omega, h)
            p0, _ = ring_power(side, wall, s, 100.0, omega, h)
            print(f"w={omega:<4} ac={ac!s:5} sigma={s:.4f} P={p:.4g} P(c=100)={p0:.4g} "
                  f"ratio={p / p0:.4f} |B_hole|={abs(bc):.3f}", flush=True)
