"""The exact magnetoquasistatic heats of nested square rings in a uniform axial AC field
(eddy2d.py generalised): div(D grad B) = i w B, D = 1/(mu0 sigma(w)) in the walls, huge in the
vacuum (B uniform there), B = 1 (the drive's amplitude) outside the outermost ring.
Each ring's time-averaged heat per unit length (1/2) Re(1/sigma) int |J|^2, for B0 = 1.

Rings are (half side, wall) centred on the origin, nested, apart.
Usage: py -3.14 transformer2d.py"""

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from eddy2d import sigma  # noqa: E402


def heats(rings, sig, c, omega, h=0.1):
    """Each ring's heat (B0 = 1) and B at the centre; `sig` one conductivity or one per ring."""
    sigs = list(sig) if isinstance(sig, (list, tuple)) else [sig] * len(rings)
    mu0 = 4 * np.pi / c ** 2
    half = max(r[0] for r in rings)
    n = int(round(2 * half / h))
    xs = np.arange(n + 1) * h - half

    def which(x, y):
        """The ring a point is in the wall of, else -1."""
        ax, ay = np.abs(x), np.abs(y)
        for k, (hs, w) in enumerate(rings):
            outer = (ax <= hs) & (ay <= hs)
            inner = (ax < hs - w) & (ay < hs - w)
            if outer and not inner:
                return k
        return -1

    d_walls = [1.0 / (mu0 * s) for s in sigs]
    d_vac = 1e6 * max(abs(d) for d in d_walls)
    m = n - 1
    rows, cols, vals = [], [], []
    rhs = np.zeros(m * m, dtype=complex)
    idx = lambda i, j: (i - 1) + m * (j - 1)
    for j in range(1, n):
        for i in range(1, n):
            k = idx(i, j)
            diag = -1j * omega * h * h
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xm, ym = xs[i] + 0.5 * di * h, xs[j] + 0.5 * dj * h
                kk = which(xm, ym)
                dd = d_walls[kk] if kk >= 0 else d_vac
                diag -= dd
                ii, jj = i + di, j + dj
                if 1 <= ii <= n - 1 and 1 <= jj <= n - 1:
                    rows.append(k)
                    cols.append(idx(ii, jj))
                    vals.append(dd)
                else:
                    rhs[k] -= dd
            rows.append(k)
            cols.append(k)
            vals.append(diag)
    a = sp.csr_matrix((vals, (rows, cols)), shape=(m * m, m * m))
    b = spla.spsolve(a.tocsc(), rhs)
    full = np.ones((n + 1, n + 1), dtype=complex)
    full[1:n, 1:n] = b.reshape(m, m).T
    rhos = [(1 / s).real for s in sigs]
    p = [0.0] * len(rings)
    for i in range(n):
        for j in range(n + 1):
            k = which(xs[i] + 0.5 * h, xs[j])
            if k >= 0:
                jy = -(full[i + 1, j] - full[i, j]) / h / mu0
                wgt = 0.5 if j in (0, n) else 1.0
                p[k] += wgt * 0.5 * rhos[k] * abs(jy) ** 2 * h * h
    for i in range(n + 1):
        for j in range(n):
            k = which(xs[i], xs[j] + 0.5 * h)
            if k >= 0:
                jx = (full[i, j + 1] - full[i, j]) / h / mu0
                wgt = 0.5 if i in (0, n) else 1.0
                p[k] += wgt * 0.5 * rhos[k] * abs(jx) ** 2 * h * h
    return p, full[n // 2, n // 2]


def scan(label, rings, density, speed, mfp, c, omegas, h=0.1):
    print(f"== {label}: rings {rings}, density {density}, speed {speed}, mfp {mfp}, c {c}; "
          f"omega_p {np.sqrt(4 * np.pi * density):.2f}, c/omega_p {c / np.sqrt(4 * np.pi * density):.2f}",
          flush=True)
    for w in omegas:
        s = sigma(density, speed, mfp, w)
        p, bc = heats(rings, s, c, w, h)
        alone = [heats([r], s, c, w, h)[0][0] for r in rings]
        print(f"  w {w:<6.3f} c/w {c / w:6.1f}: heats " + " ".join(f"{x:.4g}" for x in p)
              + "  ratio inner/outer " + f"{p[-1] / p[0]:.4f}"
              + "  alone " + " ".join(f"{x:.4g}" for x in alone)
              + f"  |B centre| {abs(bc):.3f}", flush=True)


if __name__ == "__main__":
    omegas = [0.03 * (100 ** (k / 10)) for k in range(11)]  # 0.03 .. 3
    # Outer ring 24 x 24 wall 4 (half 12), inner 10 x 10 wall 2 (half 5).
    scan("thick outer, density 1, c 3.2", [(12.0, 4.0), (5.0, 2.0)], 1.0, 2.0, 1.0, 3.2, omegas)
