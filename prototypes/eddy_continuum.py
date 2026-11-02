"""OP-1: the continuum model of an isolated metal rectangle in a ramping field, to compare
with the particles (crates/physics/tests/drude.rs D8 and the strip explorations).

Charge conservation on a finite-volume grid (no current through the walls), Ohm's law with
the model's Hall tensor (the Boltzmann mobility of drude.rs, D1/D4), the carriers' thermal
diffusion (the electrochemical field -grad rho_net, kT/q = -1 for q/m = -1, rho_c0 = -1:
the Debye length 1/sqrt(4 pi) = 0.28), the potential of the net charge with the
z-invariant kernel -2 Q ln r, and the induced field E = -1/2 Bdot z x (x - axis).

usage: eddy_continuum.py W B PERIOD B0 END [H] [DT]
"""
import math
import sys

import numpy as np

W, B, PERIOD, B0, END = (float(a) for a in sys.argv[1:6])
H = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
DT = float(sys.argv[7]) if len(sys.argv) > 7 else 0.02
# A grounded electrode, a square 0.5 wide at (EX, EY) (outside the metal): EX = none for none.
import os
EX = os.environ.get("EM_EX")
ELECTRODE = None if EX is None else (float(EX), float(os.environ.get("EM_EY", "-1.5")))
KAPPA, COS_THETA = 2.0, -1.0 / 3.0


def mobility(kappa, omega):
    """drude.rs `mobility`: (v_x + i v_y)/a of carriers (thermal speed 1) colliding at kappa v,
    rotating at -omega (dv/dt = a - i omega v), Simpson on [0, 14]."""
    n, top = 28000, 14.0
    h = top / n
    v = np.linspace(0.0, top, n + 1)
    w = np.where((np.arange(n + 1) % 2) == 1, 4.0, 2.0)
    w[0] = w[-1] = 1.0
    d = kappa * kappa * v * v + omega * omega
    d[d == 0] = np.inf
    re, im = kappa * kappa * v / d, -kappa * omega / d
    e = np.exp(-0.5 * v * v)
    fp = v ** 3 * e / (2 * math.pi)
    fq = v * math.sqrt(2 / math.pi) * v * v * e
    p = complex(np.sum(w * fp * re) * h / 3, np.sum(w * fp * im) * h / 3)
    q = complex(np.sum(w * fq * re) * h / 3, np.sum(w * fq * im) * h / 3)
    return math.pi / kappa * p / (1 - COS_THETA * q)


nx, ny = int(round(W / H)), int(round(B / H))
xc = (np.arange(nx) + 0.5) * H
yc = (np.arange(ny) + 0.5) * H
X, Y = np.meshgrid(xc, yc, indexing="ij")
ax, ay = 0.5 * W, 0.5 * B
# Green's matrix: phi_i = sum_j -2 Q_j ln |r_i - r_j|, the self term the square cell's mean ln r.
px, py = X.ravel(), Y.ravel()
dx = px[:, None] - px[None, :]
dy = py[:, None] - py[None, :]
r2 = dx * dx + dy * dy
np.fill_diagonal(r2, 1.0)
K = -np.log(r2)  # -2 ln r
sub = (np.arange(40) + 0.5) / 40 - 0.5
SX, SY = np.meshgrid(sub, sub)
mean_ln = np.mean(0.5 * np.log(SX * SX + SY * SY))  # unit square, about -1.0466
np.fill_diagonal(K, -2.0 * (math.log(H) + mean_ln))
K *= H * H  # times the cell's area: phi = K @ rho

# The electrode: its charge Q_e (on a 4 x 4 grid of points over its square) so that its centre's
# potential is 0: phi_e = sum_strip -2 Q_j ln|r_e - r_j| + Q_e * self_e = 0.
if ELECTRODE is not None:
    exx, eyy = ELECTRODE
    ed = np.sqrt((px - exx) ** 2 + (py - eyy) ** 2)
    k_to_e = -2.0 * np.log(ed) * H * H  # potential at the electrode per unit density in a cell
    k_from_e = -2.0 * np.log(ed)  # potential at the cells per unit charge on the electrode
    self_e = -2.0 * (math.log(0.5) + mean_ln)  # a uniformly charged square 0.5 wide, at its centre
cache = {}


def tensor(bz):
    """sigma_xx, sigma_H at B (q/m = -1: omega_c = -B), rounded to cache."""
    key = round(bz, 4)
    if key not in cache:
        m = mobility(KAPPA, -key)
        cache[key] = (m.real, m.imag)
    return cache[key]


def wave(t):
    """The triangle between -B0 and B0 (rising first), and its slope."""
    s = t - PERIOD * math.floor(t / PERIOD)
    half = 0.5 * PERIOD
    if s < half:
        return -B0 + 2 * B0 * s / half, 2 * B0 / half
    return B0 - 2 * B0 * (s - half) / half, -2 * B0 / half


rho = np.zeros((nx, ny))
# Face positions.
fx_x = (np.arange(nx - 1) + 1.0) * H  # x-faces between i and i+1
bins = 20
left = np.zeros(bins)
power = np.zeros(bins)
count = np.zeros(bins)
charge_e = np.zeros(bins)
q_e = 0.0
t = 0.0
while t < END:
    bz, bdot = wave(t)
    sxx, sh = tensor(bz)
    phi_flat = K @ rho.ravel()
    if ELECTRODE is not None:
        q_e = -np.dot(k_to_e, rho.ravel()) / self_e
        phi_flat = phi_flat + q_e * k_from_e
    phi = phi_flat.reshape(nx, ny)
    u = phi + rho  # potential plus the diffusion's: G = E_ind - grad(u)
    # The cells' gradients of u (central inside, one-sided at the walls).
    gux = np.gradient(u, H, axis=0)
    guy = np.gradient(u, H, axis=1)
    eix = 0.5 * bdot * (Y - ay)
    eiy = -0.5 * bdot * (X - ax)
    gx_c, gy_c = eix - gux, eiy - guy
    # x-faces (between i and i+1): normal G_x, tangential G_y averaged.
    yx = Y[:-1, :]
    gx = 0.5 * bdot * (yx - ay) - (u[1:, :] - u[:-1, :]) / H
    gyf = 0.5 * (gy_c[1:, :] + gy_c[:-1, :])
    jx = sxx * gx - sh * gyf
    # y-faces (between j and j+1).
    xy = X[:, :-1]
    gy = -0.5 * bdot * (xy - ax) - (u[:, 1:] - u[:, :-1]) / H
    gxf = 0.5 * (gx_c[:, 1:] + gx_c[:, :-1])
    jy = sxx * gy + sh * gxf
    div = np.zeros((nx, ny))
    div[:-1, :] += jx / H
    div[1:, :] -= jx / H
    div[:, :-1] += jy / H
    div[:, 1:] -= jy / H
    # The measures: the left half's current along y (faces within it) and the power.
    s = t - PERIOD * math.floor(t / PERIOD)
    rising = s < 0.5 * PERIOD
    since = s if rising else s - 0.5 * PERIOD
    k = min(int(since / (0.5 * PERIOD) * bins), bins - 1)
    sign = 1.0 if rising else -1.0
    inner_left = xy < 0.5 * W
    left[k] += sign * np.sum(jy[inner_left]) * H * H
    power[k] += (np.sum(jx * 0.5 * bdot * (yx - ay)) + np.sum(jy * (-0.5 * bdot * (xy - ax)))) * H * H
    count[k] += 1
    charge_e[k] += q_e
    rho -= DT * div
    t += DT


def torsion_area(w, b):
    s = sum(math.tanh(n * math.pi * b / (2 * w)) / n ** 5 for n in range(1, 400, 2))
    return w ** 3 * b / 12 * (1 - 192 * w / (math.pi ** 5 * b) * s)


def torsion_midline(w, b):
    s = sum((1 if (n // 2) % 2 == 0 else -1) * math.tanh(n * math.pi * b / (2 * w)) / n ** 4
            for n in range(1, 400, 2))
    return w * w * b / 8 - 8 * w ** 3 / math.pi ** 4 * s


sigma = mobility(KAPPA, 0.0).real
bdot = 4 * B0 / PERIOD
row = lambda v: " ".join(f"{v[k] / count[k]:.2f}" for k in range(bins))
print(f"electrode {ELECTRODE}: continuum {W} x {B}, h {H}, ramps of {0.5 * PERIOD}, b0 {B0}: left-half current "
      f"(theory {sigma * bdot * torsion_midline(W, B):.3f}): {row(left)}")
print(f"  electrode's charge: {row(charge_e)}")
print(f"  power (theory {sigma * bdot * bdot * torsion_area(W, B):.3f}): {row(power)}")
