"""OP-1: the continuum model on a grid fine enough for the Debye layers (h ~ 0.1 against
lambda_D = 0.28), the free-space log kernel by zero-padded FFT convolution; as
eddy_continuum.py otherwise (Ohm's law with the model's Hall tensor, diffusion, charge
conservation, no current through the walls), with an optional grounded electrode (EM_EX,
EM_EY: a square 0.5 wide held at 0). Also the current along y in the wall layers (0.5).

usage: eddy_continuum_fft.py W B PERIOD B0 END [H] [DT]
"""
import math
import os
import sys

import numpy as np

W, B, PERIOD, B0, END = (float(a) for a in sys.argv[1:6])
H = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
DT = float(sys.argv[7]) if len(sys.argv) > 7 else 0.01
KAPPA, COS_THETA = 2.0, -1.0 / 3.0
EX = os.environ.get("EM_EX")
ELECTRODE = None if EX is None else (float(EX), float(os.environ.get("EM_EY", "-1.5")))


def mobility(kappa, omega):
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
# The kernel on the padded grid: -2 ln r times the cell's area, its centre the cell's mean.
sub = (np.arange(40) + 0.5) / 40 - 0.5
SX, SY = np.meshgrid(sub, sub)
mean_ln = np.mean(0.5 * np.log(SX * SX + SY * SY))
ix = np.arange(2 * nx)
iy = np.arange(2 * ny)
dxs = np.where(ix < nx, ix, ix - 2 * nx) * H
dys = np.where(iy < ny, iy, iy - 2 * ny) * H
DX, DY = np.meshgrid(dxs, dys, indexing="ij")
r2 = DX * DX + DY * DY
r2[0, 0] = 1.0
kern = -np.log(r2)
kern[0, 0] = -2.0 * (math.log(H) + mean_ln)
kern *= H * H
KF = np.fft.rfft2(kern)


def potential(rho):
    pad = np.zeros((2 * nx, 2 * ny))
    pad[:nx, :ny] = rho
    return np.fft.irfft2(np.fft.rfft2(pad) * KF, s=(2 * nx, 2 * ny))[:nx, :ny]


if ELECTRODE is not None:
    exx, eyy = ELECTRODE
    ed = np.sqrt((X - exx) ** 2 + (Y - eyy) ** 2)
    k_to_e = -2.0 * np.log(ed) * H * H
    k_from_e = -2.0 * np.log(ed)
    self_e = -2.0 * (math.log(0.5) + mean_ln)

cache = {}


def tensor(bz):
    key = round(bz, 4)
    if key not in cache:
        m = mobility(KAPPA, -key)
        cache[key] = (m.real, m.imag)
    return cache[key]


def wave(t):
    s = t - PERIOD * math.floor(t / PERIOD)
    half = 0.5 * PERIOD
    if s < half:
        return -B0 + 2 * B0 * s / half, 2 * B0 / half
    return B0 - 2 * B0 * (s - half) / half, -2 * B0 / half


rho = np.zeros((nx, ny))
bins = 10
cols = 5  # left layer, left inside, right inside, right layer, power
acc = np.zeros((bins, cols))
qe_acc = np.zeros(bins)
count = np.zeros(bins)
layer = int(round(0.5 / H))
t = 0.0
q_e = 0.0
while t < END:
    bz, bdot = wave(t)
    sxx, sh = tensor(bz)
    phi = potential(rho)
    if ELECTRODE is not None:
        q_e = -np.sum(k_to_e * rho) / self_e
        phi = phi + q_e * k_from_e
    u = phi + rho
    gux = np.gradient(u, H, axis=0)
    guy = np.gradient(u, H, axis=1)
    gx_c = 0.5 * bdot * (Y - ay) - gux
    gy_c = -0.5 * bdot * (X - ax) - guy
    yx = Y[:-1, :]
    gx = 0.5 * bdot * (yx - ay) - (u[1:, :] - u[:-1, :]) / H
    gyf = 0.5 * (gy_c[1:, :] + gy_c[:-1, :])
    jx = sxx * gx - sh * gyf
    xy = X[:, :-1]
    gy = -0.5 * bdot * (xy - ax) - (u[:, 1:] - u[:, :-1]) / H
    gxf = 0.5 * (gx_c[:, 1:] + gx_c[:, :-1])
    jy = sxx * gy + sh * gxf
    div = np.zeros((nx, ny))
    div[:-1, :] += jx / H
    div[1:, :] -= jx / H
    div[:, :-1] += jy / H
    div[:, 1:] -= jy / H
    s = t - PERIOD * math.floor(t / PERIOD)
    rising = s < 0.5 * PERIOD
    since = s if rising else s - 0.5 * PERIOD
    k = min(int(since / (0.5 * PERIOD) * bins), bins - 1)
    sign = 1.0 if rising else -1.0
    a = H * H
    acc[k, 0] += sign * np.sum(jy[:layer, :]) * a
    acc[k, 1] += sign * np.sum(jy[layer:nx // 2, :]) * a
    acc[k, 2] += sign * np.sum(jy[nx // 2:nx - layer, :]) * a
    acc[k, 3] += sign * np.sum(jy[nx - layer:, :]) * a
    acc[k, 4] += (np.sum(jx * 0.5 * bdot * (yx - ay)) + np.sum(jy * (-0.5 * bdot * (xy - ax)))) * a
    qe_acc[k] += q_e
    count[k] += 1
    rho -= DT * div
    t += DT

print(f"continuum (fft) {W} x {B}, h {H}, dt {DT}, ramps of {0.5 * PERIOD}, b0 {B0}, electrode {ELECTRODE}:")
for k in range(bins):
    c = count[k]
    print(f"  bin {k}: along y: left layer {acc[k, 0] / c:.2f}, left inside {acc[k, 1] / c:.2f}, "
          f"right inside {acc[k, 2] / c:.2f}, right layer {acc[k, 3] / c:.2f}; power {acc[k, 4] / c:.3f}; "
          f"electrode's charge {qe_acc[k] / c:.3f}")
