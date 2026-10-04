"""Independent checks of the areas the formula audit left open.

Does not import the game. Prints one residual or a short verdict per check.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PHYS = ROOT / "crates" / "physics" / "src"
LEVELS = ROOT / "levels"


def show(name: str, value: float) -> None:
    print(f"{name}: {value:.6e}")


# ---------- DOP853 order conditions (Hairer), from the generated table ----------

def parse_rust_array(text: str, name: str) -> np.ndarray:
    # Non-greedy up to `=`. A greedy `.*` binds every later array to D.
    m = re.search(rf"pub const {name}:.*?=\s*\[(.*?)\];", text, re.S)
    if not m:
        raise SystemExit(f"missing {name}")
    nums = re.findall(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][+-]?\d+)?", m.group(1))
    return np.array([float(x) for x in nums], dtype=np.float64)


def dop853() -> None:
    text = (PHYS / "integrator" / "dop853_coefficients.rs").read_text(encoding="utf-8")
    c = parse_rust_array(text, "C")
    b = parse_rust_array(text, "B")
    e = parse_rust_array(text, "E")
    bhh = parse_rust_array(text, "BHH")
    # A is 17 rows; take the first bracket group after `pub const A`.
    am = re.search(r"pub const A:.*?=\s*\[(.*?)\];", text, re.S)
    rows = re.findall(r"\[(.*?)\]", am.group(1), re.S)
    a = np.array(
        [[float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][+-]?\d+)?", row)] for row in rows]
    )
    row_err = 0.0
    for s in list(range(2, 13)) + [14, 15, 16]:
        scale = np.abs(a[s]).sum()
        row_err = max(row_err, abs(a[s].sum() - c[s]) / max(scale, 1.0))
    quad = 0.0
    for q in range(1, 9):
        s = sum(b[i] * c[i] ** (q - 1) for i in range(1, 13))
        quad = max(quad, abs(s - 1.0 / q))
    err_ann = 0.0
    for q in range(1, 6):
        s = sum(e[i] * c[i] ** (q - 1) for i in range(1, 13))
        err_ann = max(err_ann, abs(s))
    bhh_err = 0.0
    nodes = [c[1], c[9], c[12]]
    for q in range(1, 4):
        s = sum(bhh[i + 1] * nodes[i] ** (q - 1) for i in range(3))
        bhh_err = max(bhh_err, abs(s - 1.0 / q))
    show("dop853 row-sum relative", row_err)
    show("dop853 quadrature to order 8", quad)
    show("dop853 5th estimator annihilates deg 4", err_ann)
    show("dop853 3rd estimator to order 3", bhh_err)


# ---------- flat triangle, Wilton / Graglia, against barycentric quadrature ----------

def triangle_integrals(a, b, c, r):
    """Published edge formula. Returns (int 1/R, int grad_r(1/R))."""
    a, b, c, r = map(np.asarray, (a, b, c, r))
    n = np.cross(b - a, c - a)
    n = n / np.linalg.norm(n)
    h = np.dot(n, r - a)
    rho = r - n * h
    edges = [(a, b), (b, c), (c, a)]
    pot = 0.0
    grad = np.zeros(3)
    beta_sum = 0.0
    for p, q in edges:
        l = q - p
        l = l / np.linalg.norm(l)
        m = np.cross(l, n)
        p0 = np.dot(p - rho, m)
        lm = np.dot(p - rho, l)
        lp = np.dot(q - rho, l)
        rm = np.linalg.norm(r - p)
        rp = np.linalg.norm(r - q)
        r0_sq = p0 * p0 + h * h
        if lm > 0.0 or rm + lm > 1e-12 * max(rm, 1e-300):
            f = math.log((rp + lp) / (rm + lm))
        else:
            f = math.log((rm - lm) / (rp - lp))
        if r0_sq > 0.0:
            beta = math.atan(p0 * lp / (r0_sq + abs(h) * rp)) - math.atan(
                p0 * lm / (r0_sq + abs(h) * rm)
            )
        else:
            beta = 0.0
        pot += p0 * f
        beta_sum += beta
        grad -= m * f
    pot -= abs(h) * beta_sum
    grad -= n * (math.copysign(1.0, h) * beta_sum if h != 0.0 else 0.0)
    if h == 0.0:
        grad -= n * np.dot(grad, n)
    return pot, grad


def quad_triangle(a, b, c, r, n=48):
    """Gauss–Legendre on the barycentric square, independent of the edge formula."""
    xs, ws = np.polynomial.legendre.leggauss(n)
    a, b, c, r = map(np.asarray, (a, b, c, r))
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    pot = 0.0
    grad = np.zeros(3)
    for u, wu in zip(xs, ws):
        for v, wv in zip(xs, ws):
            # map [-1,1]^2 onto u>=0, v>=0, u+v<=1
            uu = 0.5 * (u + 1.0)
            vv = 0.5 * (v + 1.0) * (1.0 - uu)
            jac = 0.25 * (1.0 - uu)
            x = a + uu * (b - a) + vv * (c - a)
            d = r - x
            dist = np.linalg.norm(d)
            w = wu * wv * jac * 2.0 * area
            pot += w / dist
            # ∇_r (1/|r−x|) = −(r−x)/|r−x|³. The previous sign made this the
            # negative of the edge formula, a relative difference of 2.
            grad -= w * d / dist**3
    return pot, grad


def panels() -> None:
    a = np.array([0.2, -0.1, 0.3])
    b = np.array([1.4, 0.3, 0.1])
    c = np.array([0.5, 1.2, -0.2])
    r = np.array([0.4, 0.6, 1.5])
    p1, g1 = triangle_integrals(a, b, c, r)
    p2, g2 = quad_triangle(a, b, c, r, n=64)
    show("panel potential vs quadrature", abs(p1 - p2) / abs(p2))
    show("panel gradient vs quadrature", np.linalg.norm(g1 - g2) / np.linalg.norm(g2))
    # The analytic gradient must be the derivative of the analytic potential.
    h = 1e-6
    num = np.zeros(3)
    for i in range(3):
        step = np.zeros(3)
        step[i] = h
        pp, _ = triangle_integrals(a, b, c, r + step)
        pm, _ = triangle_integrals(a, b, c, r - step)
        num[i] = (pp - pm) / (2.0 * h)
    show("panel grad vs d(potential)", np.linalg.norm(g1 - num) / np.linalg.norm(num))
    # Jump of the normal integral. The crate's test is (g_minus - g_plus)·n = 4π.
    mid = (a + b + c) / 3.0
    n = np.cross(b - a, c - a)
    n = n / np.linalg.norm(n)
    eps = 1e-6
    _, gp = triangle_integrals(a, b, c, mid + n * eps)
    _, gm = triangle_integrals(a, b, c, mid - n * eps)
    jump = np.dot(gm - gp, n)
    show("panel normal jump minus 4 pi", abs(jump - 4.0 * math.pi))
    # Display substitution: σ·area at the centroid versus the integral, one cell away.
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    centroid = (a + b + c) / 3.0
    obs = centroid + np.array([1.0, 0.0, 0.0])
    exact, _ = triangle_integrals(a, b, c, obs)
    mono = area / np.linalg.norm(obs - centroid)
    show("display monopole vs panel, relative", abs(mono - exact) / abs(exact))


# ---------- geometry, from the definitions ----------

def sphere_sd(center, radius, x):
    return np.linalg.norm(x - center) - radius


def torus_sd(center, normal, major, minor, x):
    d = x - center
    axial = np.dot(d, normal)
    radial = np.linalg.norm(d - normal * axial)
    return math.hypot(radial - major, axial) - minor


def capsule_sd(a, b, radius, x):
    ab = b - a
    l2 = np.dot(ab, ab)
    t = 0.0 if l2 == 0.0 else float(np.clip(np.dot(x - a, ab) / l2, 0.0, 1.0))
    return np.linalg.norm(x - (a + ab * t)) - radius


def aabb_sd(lo, hi, x):
    center = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    q = np.abs(x - center) - half
    return np.linalg.norm(np.maximum(q, 0.0)) + min(float(q.max()), 0.0)


def oriented_box_sd(center, axis, half, margin, x):
    axis = axis / np.linalg.norm(axis)
    v = np.array([-axis[1], axis[0], 0.0])
    d = x - center
    local = np.array([np.dot(d, axis), np.dot(d, v), d[2]])
    q = np.abs(local) - half
    outside = np.linalg.norm(np.maximum(q, 0.0))
    inside = min(float(q.max()), 0.0)
    return outside + inside - margin


def geometry() -> None:
    # Each residual is the transcribed SDF against a distance built another way.
    center = np.array([1.0, -2.0, 0.5])
    x = np.array([4.0, -2.0, 0.5])
    show("sphere sdf", abs(sphere_sd(center, 1.5, x) - math.hypot(3.0, 0.0) + 1.5))
    show(
        "torus exterior",
        abs(torus_sd(np.zeros(3), np.array([0.0, 0.0, 1.0]), 5.0, 0.2, np.array([0.0, 7.0, 0.0])) - (2.0 - 0.2)),
    )
    show(
        "torus above ring",
        abs(torus_sd(np.zeros(3), np.array([0.0, 0.0, 1.0]), 5.0, 0.2, np.array([5.0, 0.0, 1.0])) - (1.0 - 0.2)),
    )
    a = np.zeros(3)
    b = np.array([4.0, 0.0, 0.0])
    show("capsule side", abs(capsule_sd(a, b, 0.1, np.array([2.0, 3.0, 0.0])) - (3.0 - 0.1)))
    show("capsule end", abs(capsule_sd(a, b, 0.1, np.array([7.0, 4.0, 0.0])) - (5.0 - 0.1)))
    show(
        "zero capsule",
        abs(capsule_sd(np.array([1.0, 2.0, 0.0]), np.array([1.0, 2.0, 0.0]), 0.1, np.array([4.0, 6.0, 0.0])) - (5.0 - 0.1)),
    )
    lo = np.array([-1.0, -2.0, -3.0])
    hi = np.array([1.0, 2.0, 3.0])
    show("aabb inside", abs(aabb_sd(lo, hi, np.zeros(3)) - (-1.0)))
    show("aabb corner", abs(aabb_sd(lo, hi, np.array([4.0, 6.0, 3.0])) - 5.0))
    half = np.array([2.0, 1.0, 0.5])
    along_x = oriented_box_sd(np.zeros(3), np.array([1.0, 0.0, 0.0]), half, 0.0, np.array([5.0, 0.0, 0.0]))
    along_y = oriented_box_sd(np.zeros(3), np.array([0.0, 1.0, 0.0]), half, 0.0, np.array([0.0, 5.0, 0.0]))
    # The rotated query must match the axis-aligned distance of the unrotated point.
    show("oriented box rotation", abs(along_x - along_y) + abs(along_x - aabb_sd(-half, half, np.array([5.0, 0.0, 0.0]))))


# ---------- certified crossing ----------

def first_crossing(g, ta, tb, ga, gb, v_max, budget=10_000):
    if math.isnan(ga) or math.isnan(gb):
        return None
    if gb > 0.0 and ga + gb - v_max * (tb - ta) > 0.0:
        return None
    tm = 0.5 * (ta + tb)
    if tm <= ta or tm >= tb or budget <= 0:
        return tb if gb <= 0.0 else None
    gm = g(tm)
    if gm <= 0.0:
        earlier = first_crossing(g, ta, tm, ga, gm, v_max, budget - 1)
        return tm if earlier is None else earlier
    left = first_crossing(g, ta, tm, ga, gm, v_max, budget - 1)
    if left is not None:
        return left
    return first_crossing(g, tm, tb, gm, gb, v_max, budget - 1)


def events() -> None:
    root = first_crossing(lambda t: 1.0 - t, 0.0, 2.0, 1.0, -1.0, 1.0)
    show("crossing at 1", abs(root - 1.0))
    dip = lambda t: (t - 0.5) ** 2 - 1e-6
    t_dip = first_crossing(dip, 0.0, 1.0, dip(0.0), dip(1.0), 2.0)
    show("dip found", abs(t_dip - (0.5 - 1e-3)))
    near = lambda t: (t - 0.5) ** 2 + 1e-6
    missed = first_crossing(near, 0.0, 1.0, near(0.0), near(1.0), 2.0)
    show("near miss certified clear", 0.0 if missed is None else 1.0)
    # A bound smaller than the true slope certifies a real dip as clear.
    hidden = first_crossing(dip, 0.0, 1.0, dip(0.0), dip(1.0), 0.1)
    show("underestimated v_max misses the dip", 0.0 if hidden is None else 1.0)
    nan = first_crossing(lambda t: t, 0.0, 1.0, float("nan"), 1.0, 1.0)
    show("nan endpoint dropped", 0.0 if nan is None else 1.0)


# ---------- elastic bisection against a scalar root ----------

def kinetic(m, c, p):
    return math.sqrt(m * m * c * c * c * c + p * p * c * c) - m * c * c


def elastic_root(m1, m2, c, p1, p2):
    """1D, n = +1. Positive root of the energy balance, by the code's bracket."""
    v1 = p1 * c * c / math.sqrt(m1 * m1 * c**4 + p1 * p1 * c * c)
    v2 = p2 * c * c / math.sqrt(m2 * m2 * c**4 + p2 * p2 * c * c)
    closing = v1 - v2
    if closing <= 0.0:
        return None
    before = kinetic(m1, c, p1) + kinetic(m2, c, p2)

    def f(j):
        return kinetic(m1, c, p1 - j) + kinetic(m2, c, p2 + j) - before

    newtonian = 2.0 * m1 * m2 / (m1 + m2) * closing
    lo, hi = 0.5 * newtonian, 2.0 * newtonian
    while f(lo) > 0.0:
        lo *= 0.5
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), f


def elastic() -> None:
    # Unequal masses, target moving. Compare the bisection to a bisection started
    # from a much wider bracket (the same convex function, an independent bracket).
    root, f = elastic_root(3.0, 5.0, 4.0, 2.0, -1.0)
    lo, hi = 0.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    show("elastic bisection vs wide bracket", abs(root - 0.5 * (lo + hi)) / root)
    show("elastic residual of the balance", abs(f(root)) / kinetic(3.0, 4.0, 2.0))
    show("elastic not approaching", 0.0 if elastic_root(3.0, 5.0, 4.0, 1.0, 2.0) is None else 1.0)


# ---------- Liénard radiation term: power of a charge at rest ----------

def larmor_numeric() -> None:
    q, c, a = 1.3, 5.0, 0.4
    acc = np.array([0.0, a, 0.0])
    # Fibonacci directions. E_rad = q/(c^2 R) ((n·a)n - a), B = n×E/c, S = (c^2/4π) E×B.
    n_dir = 2000
    golden = math.pi * (3.0 - math.sqrt(5.0))
    power = 0.0
    r = 2.0
    for i in range(n_dir):
        z = 1.0 - 2.0 * (i + 0.5) / n_dir
        rho = math.sqrt(max(0.0, 1.0 - z * z))
        phi = golden * i
        n = np.array([rho * math.cos(phi), rho * math.sin(phi), z])
        e = (q / (c * c * r)) * ((np.dot(n, acc) * n) - acc)
        b = np.cross(n, e) / c
        s = (c * c / (4.0 * math.pi)) * np.cross(e, b)
        # equal-area directions: each represents 4π/n
        power += np.dot(s, n) * r * r * (4.0 * math.pi / n_dir)
    claim = (2.0 / 3.0) * q * q * a * a / c**3
    show("radiation term vs Larmor", abs(power - claim) / claim)


# ---------- uniform-field motion, hyperbolic branch, against RK4 ----------

def field_motion() -> None:
    q, m, c = 0.7, 1.3, 5.0
    e = np.array([2.5, 1.0])  # E > c B
    b = 0.2
    v0 = np.array([1.0, -0.4])
    gamma = 1.0 / math.sqrt(1.0 - np.dot(v0, v0) / c**2)
    u0 = np.array([gamma * c, gamma * v0[0], gamma * v0[1]])
    # Λ rows on (U0, Ux, Uy)
    scale = q / m
    lam = scale * np.array(
        [
            [0.0, e[0] / c, e[1] / c],
            [e[0] / c, 0.0, b],
            [e[1] / c, -b, 0.0],
        ]
    )
    w2 = (q / m) ** 2 * (b * b - np.dot(e, e) / c**2)
    s = 0.4
    k = math.sqrt(-w2)
    sh, ch = math.sinh(k * s), math.cosh(k * s)
    half = math.sinh(0.5 * k * s)
    c0 = ch
    s1 = sh / k
    c1 = 2.0 * half * half / (-w2)
    s2 = (sh / k - s) / (-w2)
    lam_u = lam @ u0
    lam2_u = lam @ lam_u
    u_closed = u0 + s1 * lam_u + c1 * lam2_u
    # RK4 of dU/ds = Λ U
    n = 4000
    h = s / n
    u = u0.copy()
    for _ in range(n):
        k1 = lam @ u
        k2 = lam @ (u + 0.5 * h * k1)
        k3 = lam @ (u + 0.5 * h * k2)
        k4 = lam @ (u + h * k3)
        u = u + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    show("hyperbolic field motion vs RK4", np.linalg.norm(u_closed - u) / np.linalg.norm(u))
    # C0 should be 1 - w2 C1
    show("hyperbolic C0 identity", abs(c0 - (1.0 - w2 * c1)))


# ---------- retarded time of uniform motion, Newton's method ----------

def retarded_uniform() -> None:
    c = 5.0
    x0 = np.array([0.0, 0.0, 0.0])
    v0 = np.array([1.2, -0.3, 0.0])
    x = np.array([4.0, 1.0, 0.0])
    dt = 3.0
    # Closed form: |x - x0 - v (dt - tau)| = c tau, tau = t - t_r > 0, the past root.
    r = x - x0 - v0 * dt
    rv = np.dot(r, v0)
    v2 = np.dot(v0, v0)
    disc = rv * rv + (c * c - v2) * np.dot(r, r)
    tau = -(rv) / (c * c - v2)  # placeholder, use the physical root
    # Solve quadratic (c^2 - v^2) tau^2 + 2 (r·v) tau - |r|^2 = 0? 
    # |r + v tau| = c tau if the emission is dt - tau and position is x0+v0*(dt-tau),
    # so |x - x0 - v0*dt + v0*tau| = c tau, |r + v tau| = c tau.
    # (r + v tau)^2 = c^2 tau^2
    # |r|^2 + 2 tau r·v + tau^2 v^2 = c^2 tau^2
    # (v^2 - c^2) tau^2 + 2 (r·v) tau + |r|^2 = 0
    aa = v2 - c * c
    bb = 2.0 * rv
    cc = np.dot(r, r)
    disc = bb * bb - 4 * aa * cc
    roots = [(-bb + math.sqrt(disc)) / (2 * aa), (-bb - math.sqrt(disc)) / (2 * aa)]
    tau = min(t for t in roots if t > 0.0)
    # Newton on G(s) = c*dt - c*gamma*s wait: proper time s, x(s) = x0 + v0 * (gamma? )
    # For uniform motion, coordinate time from the reference is tau_coord = gamma * s? 
    # U0 = gamma c, dx0/ds = U0 so t(s) = t0 + (gamma c / c) *? dx[0] = U^0 * s = gamma c * s,
    # and the code's value is c*dt - dx[0] - |x - x(s)|, with x(s) = x0 + v * (gamma s)
    # because spatial U = gamma v, so Δx = gamma v s = v * Δt, Δt = gamma s.
    gamma = 1.0 / math.sqrt(1.0 - v2 / c**2)
    s_true = (dt - tau) * 0.0  # emission proper time relative to t0, observer at t0+dt
    # t_r = t0 + dt - tau, s = (t_r - t0) / gamma = (dt - tau) / gamma
    s_true = (dt - tau) / gamma

    def g(s):
        dt_s = gamma * s
        pos = x0 + v0 * dt_s
        return c * dt - (gamma * c) * s - np.linalg.norm(x - pos)

    s = 0.0
    for _ in range(20):
        # slope ≈ -gamma c + n·(gamma v) which is negative
        eps = 1e-8
        slope = (g(s + eps) - g(s - eps)) / (2 * eps)
        s = s - g(s) / slope
    show("uniform retarded Newton", abs(s - s_true))
    show("uniform retarded residual", abs(g(s_true)))


# ---------- Filon quartic ----------

def filon_moments(d: float):
    e = complex(math.cos(d), math.sin(d))
    out = [0j] * 5
    if abs(d) < 0.5:
        term = 1 + 0j
        for j in range(16):
            if j > 0 and abs(term) < 1e-17:
                break
            for k in range(5):
                out[k] += term / (j + k + 1)
            term *= 1j * d / (j + 1)
        return out
    out[0] = (e - 1) / (1j * d)
    for k in range(1, 5):
        out[k] = (e - k * out[k - 1]) / (1j * d)
    return out


def newton_powers(u, g):
    n = len(u)
    dd = list(g)
    for level in range(1, n):
        for m in range(n - 1, level - 1, -1):
            dd[m] = (dd[m] - dd[m - 1]) / (u[m] - u[m - level])
    coeffs = [0.0] * n
    coeffs[0] = dd[-1]
    degree = 0
    for m in range(n - 2, -1, -1):
        for k in range(degree + 1, 0, -1):
            coeffs[k] = coeffs[k - 1] - coeffs[k] * u[m]
        coeffs[0] = dd[m] - coeffs[0] * u[m]
        degree += 1
    return coeffs


def filon_quartic() -> None:
    # A known quartic on uneven nodes. Integrate c(u) e^{i D u} from 0 to 1.
    rng = np.random.default_rng(0)
    truth = np.array([0.3, -1.1, 0.4, 0.2, -0.5])
    u = np.array([0.0, 0.17, 0.41, 0.73, 1.0])
    g = np.array([sum(truth[k] * ui**k for k in range(5)) for ui in u])
    coeffs = np.array(newton_powers(u.tolist(), g.tolist()))
    show("quartic interpolation", np.max(np.abs(coeffs - truth)))
    d = 3.5
    phi = filon_moments(d)
    got = sum(truth[k] * phi[k] for k in range(5))
    # analytic: integrate term by term with the same recurrence's definition
    claim = sum(truth[k] * phi[k] for k in range(5))
    # independent quadrature
    xs, ws = np.polynomial.legendre.leggauss(64)
    quad = 0j
    for x, w in zip(xs, ws):
        uu = 0.5 * (x + 1.0)
        val = sum(truth[k] * uu**k for k in range(5))
        quad += 0.5 * w * val * np.exp(1j * d * uu)
    show("filon quartic vs quadrature", abs(got - quad) / abs(quad))


# ---------- pipe mode ----------

def pipe() -> None:
    # Rectangular mouth of width w. φ = sin(π x / w) e^{-π z / w} is the lowest
    # Dirichlet mode: zero on x = 0 and x = w, and Laplace's equation.
    w = 2.0
    k = math.pi / w
    xs = np.linspace(0.05, w - 0.05, 7)
    zs = np.linspace(0.0, 3.0, 5)
    worst = 0.0
    h = 1e-5
    for x in xs:
        for z in zs:
            def phi(xx, zz, k=k, w=w):
                return math.sin(k * xx) * math.exp(-k * zz)

            lap = (
                phi(x + h, z) + phi(x - h, z) + phi(x, z + h) + phi(x, z - h) - 4 * phi(x, z)
            ) / h**2
            worst = max(worst, abs(lap) / max(abs(phi(x, z)), 1e-30))
    edge = abs(math.sin(k * 0.0)) + abs(math.sin(k * w))
    # Next integer mode decays as 2k, strictly faster.
    show("pipe mode laplace relative", worst)
    show("pipe mode wall value", edge)
    show("next mode decays faster", 0.0 if 2 * k > k else 1.0)


# ---------- MFS for one grounded sphere, against Kelvin ----------

def mfs() -> None:
    a = 1.0
    src = np.array([3.0, 0.4, -0.2])
    q = 1.7
    n = 80
    golden = math.pi * (3.0 - math.sqrt(5.0))

    def fib(n, twist):
        out = []
        for i in range(n):
            z = 1.0 - 2.0 * (i + 0.5) / n
            rho = math.sqrt(max(0.0, 1.0 - z * z))
            phi = golden * i + twist
            out.append([rho * math.cos(phi), rho * math.sin(phi), z])
        return np.array(out)

    y = 0.6 * a * fib(n, 0.0)
    x = a * fib(2 * n, 0.7)
    # Cancel the external potential on the sphere.
    rhs = -q / np.linalg.norm(x - src, axis=1)
    mat = 1.0 / np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)
    coef, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    # Boundary residual on a third lattice.
    xt = a * fib(300, 1.3)
    pred = mat_eval = coef @ (1.0 / np.linalg.norm(xt[:, None, :] - y[None, :, :], axis=2)).T
    exact_ext = q / np.linalg.norm(xt - src, axis=1)
    resid = np.max(np.abs(pred + exact_ext))
    scale = np.max(np.abs(exact_ext))
    show("mfs boundary residual relative", resid / scale)
    # Force on the source: field of the equivalent charges, versus the Kelvin image.
    r = src  # centre at 0
    r2 = np.dot(r, r)
    qi = -q * a / math.sqrt(r2)
    pi = r * (a * a / r2)
    d = src - pi
    e_kelvin = qi * d / np.linalg.norm(d) ** 3
    e_mfs = np.zeros(3)
    for ck, yk in zip(coef, y):
        d = src - yk
        e_mfs += ck * d / np.linalg.norm(d) ** 3
    show("mfs force vs Kelvin", np.linalg.norm(e_mfs - e_kelvin) / np.linalg.norm(e_kelvin))


# ---------- verification constant ----------
# The value printed here is pub const MARGIN_SAFE in trajectory.rs.

def rust_f64(path: Path, name: str) -> float:
    text = path.read_text(encoding="utf-8")
    found = re.search(rf"pub const {name}:\s*f64\s*=\s*([0-9.eE+-]+)\s*;", text)
    if not found:
        raise SystemExit(f"missing {name} in {path.name}")
    return float(found.group(1))


def verify_rule() -> None:
    margin_safe = rust_f64(PHYS / "trajectory.rs", "MARGIN_SAFE")
    show("verify MARGIN_SAFE", margin_safe)
    # Distance from the preview bound 0.2 up to the constant just read.
    show("verify preview short of MARGIN_SAFE", margin_safe - 0.2)


# ---------- launch accelerations of the fast single shots ----------

def dipole_b(moment_z, origin, x):
    r = x - origin
    r2 = np.dot(r, r)
    inv = 1.0 / math.sqrt(r2)
    # m along z: B = (3 (m·rhat) rhat - m) / r^3, m = (0,0,μ)
    m = np.array([0.0, 0.0, moment_z])
    rhat = r * inv
    return (3.0 * np.dot(m, rhat) * rhat - m) * (inv / r2)


def coulomb_e(q, origin, x):
    d = x - origin
    r2 = np.dot(d, d)
    return q * d / (r2 * math.sqrt(r2))


def beta_of(t, m, c):
    # ε = T/(m c^2), β = sqrt(ε(ε+2))/(1+ε)
    eps = t / (m * c * c)
    return math.sqrt(eps * (eps + 2.0)) / (1.0 + eps)


def fast_shots() -> None:
    # Real curriculum filenames. Level 30 shot 0, then the six γ = 3 shots.
    wanted = {
        "30_beta_spectrometer.json": [0],
        "80_jackson_beaming.json": [0],
        "81_jackson_critical_frequency.json": [0],
        "83_jackson_quiet_turn.json": [0],
        "84_jackson_thomson.json": [0],
        "85_jackson_braking.json": [0],
        "88_jackson_undulator.json": [0],
    }
    for name, shots in wanted.items():
        path = next(LEVELS.glob(name))
        level = json.loads(path.read_text(encoding="utf-8"))
        c = level["physics"]["c"]
        grid = level["grid"]
        sub = grid.get("subdivision", 1)

        def pos(node):
            return np.array([node[0] / sub, node[1] / sub, 0.0], dtype=float)

        refs = level.get("reference_solution") or []
        for si in shots:
            shot = level["shots"][si]
            p = shot["particle"]
            m = p["mass"]
            q = p["charge"]
            tkin = shot["launch"]["kinetic_energy"]
            beta = beta_of(tkin, m, c)
            direction = np.array(shot["launch"]["direction"], dtype=float)
            direction = direction / np.linalg.norm(direction)
            x = pos(shot["launch"]["node"])
            e = np.zeros(3)
            b = np.zeros(3)
            supported = True
            for el in refs:
                kind = el["kind"]
                if kind == "magnet":
                    b += dipole_b(el["value"], pos(el["node"]), x)
                elif kind == "charge":
                    e += coulomb_e(el["value"], pos(el["node"]), x)
                else:
                    supported = False
            if level.get("coils") or level.get("elements"):
                # Reference is the placed set; static elements would add more.
                if level.get("coils") or any(
                    e.get("kind") not in ("magnet", "charge") for e in level.get("elements", [])
                ):
                    supported = supported and not level.get("coils")
            gamma = 1.0 / math.sqrt(1.0 - beta * beta)
            v = direction * beta * c
            # p_vec magnitude from T: p^2 = T^2/c^2 + 2 m T
            pmag = math.sqrt(tkin * tkin / (c * c) + 2.0 * m * tkin)
            # force
            force = q * (e + np.cross(v, b))
            inv_c2 = 1.0 / (c * c)
            acc = (force - v * (np.dot(v, force) * inv_c2)) / (gamma * m)
            if not supported or np.linalg.norm(acc) == 0.0:
                print(f"fast {name} shot {si}: beta {beta:.6f} acceleration not a pure magnet/charge field")
                continue
            ahat = acc / np.linalg.norm(acc)
            # Speed after the taper's full 0.1c along â, and along −â (the past of a brake).
            def speed(sign):
                vv = v + ahat * (sign * 0.1 * c)
                return np.linalg.norm(vv) / c

            parallel = np.dot(ahat, direction)
            print(
                f"fast {name} shot {si}: beta {beta:.6f} a_parallel {parallel:.6f} "
                f"future {speed(1):.6f}c past {speed(-1):.6f}c"
            )


def newton_control(g, gp, s0, scale, niter=100):
    """The retarded-time loop's bracket, copied onto a scalar monotone function."""
    lo, hi = -math.inf, math.inf
    s = s0
    bracket_at = None
    for iteration in range(niter):
        if iteration == 12 and not (math.isfinite(lo) and math.isfinite(hi)):
            bracket_at = 12
            step = scale
            for _ in range(64):
                if math.isfinite(lo):
                    break
                trial = min(hi, s) - step
                value = g(trial)
                if value >= 0.0:
                    lo = trial
                else:
                    hi = trial
                    step *= 2.0
            for _ in range(64):
                if math.isfinite(hi):
                    break
                trial = max(lo, s) + step
                value = g(trial)
                if value <= 0.0:
                    hi = trial
                else:
                    lo = trial
                    step *= 2.0
            if not (math.isfinite(lo) and math.isfinite(hi)):
                return None, bracket_at
            s = 0.5 * (lo + hi)
        value = g(s)
        slope = gp(s)
        if value > 0.0:
            lo = max(lo, s)
        else:
            hi = min(hi, s)
        nxt = s - value / slope
        if math.isfinite(lo) and math.isfinite(hi) and not (lo < nxt < hi):
            nxt = 0.5 * (lo + hi)
        if abs(nxt - s) <= 1e-12 or (math.isfinite(lo) and math.isfinite(hi) and hi - lo <= 1e-12):
            return s, bracket_at
        s = nxt
    return s, bracket_at


def newton_bracket() -> None:
    # Correct slope: one Newton step lands on the root, so the iteration-12
    # bracket is never opened.
    s, bracket_at = newton_control(lambda t: 100.0 - t, lambda t: -1.0, 0.0, 1.0)
    show("newton direct residual", abs(100.0 - s))
    show("newton direct skips bracket", 0.0 if bracket_at is None else 1.0)
    # A slope 20× too shallow creeps. The first twelve samples stay positive
    # and the doubling search has to open. The loop then returns its last
    # sample after 100 iterations, the same fallback as FieldMotion::retarded.
    s2, bracket_at2 = newton_control(lambda t: 100.0 - t, lambda t: -20.0, 0.0, 1.0)
    show("newton bracket engaged", 0.0 if bracket_at2 == 12 else 1.0)
    show("newton shallow-slope last sample", abs(100.0 - s2))
    # Always positive and finite: 64 doublings never see a sign change.
    missing, _ = newton_control(lambda t: 1.0, lambda t: -1e-6, 0.0, 1.0)
    show("newton no-root returns none", 0.0 if missing is None else 1.0)


def j0(x: float) -> float:
    # Power series. Enough digits to locate the first zero.
    term = 1.0
    total = 1.0
    y = -0.25 * x * x
    for k in range(1, 40):
        term *= y / (k * k)
        total += term
    return total


def pipe_bessel() -> None:
    # First positive root of J0, by bisection on (2, 3).
    lo, hi = 2.0, 3.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if j0(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    j01 = 0.5 * (lo + hi)
    # A spherical detector uses width 2R, so the game's k is π/(2R).
    # The circular pipe's lowest Dirichlet decay is j01/R.
    game_over_circular = (math.pi / 2.0) / j01
    show("bessel J0 at the root", abs(j0(j01)))
    show("sphere-cup k versus circular pipe", abs(game_over_circular - 1.0))


def agm_k(m, niter, dtype):
    one = dtype(1.0)
    half = dtype(0.5)
    a = one
    b = np.sqrt(one - dtype(m))
    for _ in range(niter):
        an = half * (a + b)
        bn = np.sqrt(a * b)
        c = half * (a - b)
        a, b = an, bn
        if abs(float(c)) <= 1e-7 * abs(float(a)) and dtype is float:
            break
    return float(np.pi / (dtype(2.0) * a))


def shader_agm() -> None:
    # Shader: 12 iterations, f32. CPU: 40 iterations, f64, early exit at f64 epsilon.
    worst = 0.0
    for m in (1e-6, 0.3, 0.9, 0.99, 1.0 - 1e-4):
        ref = agm_k(m, 40, float)
        shader = agm_k(np.float32(m).item(), 12, np.float32)
        worst = max(worst, abs(shader - ref) / ref)
    show("shader AGM vs 40-iter", worst)


def graded(n, both):
    out = []
    for k in range(n + 1):
        if both:
            out.append(0.5 * (1.0 - math.cos(math.pi * k / n)))
        else:
            out.append(math.sin(0.5 * math.pi * k / n))
    return out


def electrode_mesh(center, angle, half_l, half_t, half_h, size):
    u = np.array([math.cos(angle), math.sin(angle), 0.0])
    v = np.array([-math.sin(angle), math.cos(angle), 0.0])
    w = np.array([0.0, 0.0, 1.0])
    tris = []

    def face(origin, e1, e2, grade2_both):
        n1 = max(2, math.ceil(np.linalg.norm(e1) / size))
        n2 = max(2, math.ceil(np.linalg.norm(e2) / size))
        g1 = graded(n1, True)
        g2 = graded(n2, grade2_both)
        for i in range(n1):
            for j in range(n2):
                def p(s, t, i=i, j=j):
                    return origin + e1 * s + e2 * t

                p00, p10 = p(g1[i], g2[j]), p(g1[i + 1], g2[j])
                p01, p11 = p(g1[i], g2[j + 1]), p(g1[i + 1], g2[j + 1])
                tris.append((p00, p10, p11))
                tris.append((p00, p11, p01))

    a, b, h = half_l, half_t, half_h
    face(center - u * a - v * b + w * h, u * (2 * a), v * (2 * b), True)
    face(center + u * a - v * b, v * (2 * b), w * h, False)
    face(center - u * a + v * b, v * (-2 * b), w * h, False)
    face(center - u * a - v * b, u * (2 * a), w * h, False)
    face(center + u * a + v * b, u * (-2 * a), w * h, False)
    return tris


def display_electrode() -> None:
    # Level 19, the upper grounded plate. Display panel size is 1 cell.
    level = json.loads((LEVELS / "19_power_supply.json").read_text(encoding="utf-8"))
    e = level["electrodes"][0]
    center = np.array(e["center"], dtype=float)
    tris = electrode_mesh(
        center,
        math.radians(e["angle_deg"]),
        e["length"] / 2.0,
        e["thickness"] / 2.0,
        e["height"] / 2.0,
        1.0,
    )
    # One cell outside the +y face, in the midplane. The mirror doubles the potential.
    obs = center + np.array([0.0, e["thickness"] / 2.0 + 1.0, 0.0])
    exact = 0.0
    mono = 0.0
    for a, b, c in tris:
        pot, _ = triangle_integrals(a, b, c, obs)
        area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
        centroid = (a + b + c) / 3.0
        exact += 2.0 * pot
        mono += 2.0 * area / np.linalg.norm(obs - centroid)
    show("level19 display panels", float(len(tris)))
    show("level19 display monopole vs panels", abs(mono - exact) / abs(exact))


def curriculum_order() -> None:
    files = sorted(LEVELS.glob("[0-9]*.json"), key=lambda p: p.name)
    slugs = []
    data = json.loads((LEVELS / "curriculum.json").read_text(encoding="utf-8"))
    for arc in data["arcs"]:
        for tier in arc["tiers"]:
            slugs.extend(tier["levels"])
    file_slugs = [p.name.split("_", 1)[1].removesuffix(".json") for p in files]
    order_mismatch = sum(a != b for a, b in zip(file_slugs, slugs)) + abs(len(file_slugs) - len(slugs))
    text = (ROOT / "docs" / "LEVELS.md").read_text(encoding="utf-8")
    rows = [ln for ln in text.splitlines() if ln.startswith("| ") and ln[2:4].strip(" |").isdigit()]
    names = []
    for path in files:
        names.append(json.loads(path.read_text(encoding="utf-8"))["name"])
    name_mismatch = 0
    for row, name in zip(rows, names):
        cells = [c.strip() for c in row.strip("|").split("|")]
        if cells[1] != name:
            name_mismatch += 1
    show("curriculum slug mismatches", float(order_mismatch))
    show("LEVELS.md name mismatches", float(name_mismatch + abs(len(rows) - len(names))))


def main() -> None:
    dop853()
    panels()
    geometry()
    events()
    elastic()
    larmor_numeric()
    field_motion()
    retarded_uniform()
    filon_quartic()
    pipe()
    mfs()
    verify_rule()
    fast_shots()
    newton_bracket()
    pipe_bessel()
    shader_agm()
    display_electrode()
    curriculum_order()
    print("remainder checks finished")


if __name__ == "__main__":
    main()
