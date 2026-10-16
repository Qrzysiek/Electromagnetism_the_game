"""Companion of docs/levels/first_bend/doc.tex (Level 1, First bend).

Every number of the document, from the level file, by the document's formulas (their labels
in the comments): the exact orbit at 40 significant digits (mpmath), and the first event of
every single-charge placement along its exact orbit (numpy). Compared with the game: the
reference flight, three more placements in detail, and the verdicts of all single-charge
placements (placements.json's family).
"""
import math

import mpmath as mp
import numpy as np

DIGITS = 40
# Comparisons with the verified flights: the game's verified positions are accurate to about
# 1e-10 cells on exact orbits (Foundations, section 6); ten times that, in cells and time
# units. The document's numbers are exact to far more digits than that.
TOL = 1e-9
# The game integrates the radiated power with three-point Gauss-Legendre on each of its
# steps (PHYSICS.md 8: 3e-9 relative on a test with an exact reference); a factor of ten
# thousand above that still tells a wrong formula (a factor 2/3, a missing gamma^6) apart.
TOL_RAD = 1e-4


class Board:
    """The level's board: launch, detector, world, charge radius, time limit (Foundations 5)."""

    def __init__(self, lvl, mpf):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = mpf(shot["particle"]["charge"])
        self.m = mpf(shot["particle"]["mass"])
        self.a = mpf(shot["particle"]["radius"])
        self.c = mpf(lvl["physics"]["c"])
        self.T0 = mpf(shot["launch"]["kinetic_energy"])
        self.R = mpf(lvl["physics"]["charge_radius"])
        self.tmax = mpf(lvl["physics"]["t_max"])
        self.x0 = [mpf(v) / s for v in shot["launch"]["node"][:2]]
        d = [mpf(v) for v in shot["launch"]["direction"][:2]]
        n = (d[0] ** 2 + d[1] ** 2) ** 0.5
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [mpf(v) / s for v in shot["detector"]["min"][:2]]
        self.dmax = [mpf(v) / s for v in shot["detector"]["max"][:2]]
        self.wmin = [mpf(-1), mpf(-1)]
        self.wmax = [mpf(lvl["grid"]["nx"]) / s + 1, mpf(lvl["grid"]["ny"]) / s + 1]
        # The launch (Foundations eq:launch).
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = (self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0) ** 0.5


class Orbit:
    """The exact orbit about one charge Q at rc (section 4): in the angle theta swept in the
    direction of motion, u = 1/r = B + (u0 - B) C(theta) + u0' S(theta) (eq:general)."""

    def __init__(self, b, Q, rc, lib):
        self.b, self.lib = b, lib
        self.kappa = b.q * Q
        self.rc = rc
        rho = [b.x0[0] - rc[0], b.x0[1] - rc[1]]
        self.r0 = (rho[0] ** 2 + rho[1] ** 2) ** 0.5
        self.phi0 = mp.atan2(rho[1], rho[0]) if lib is mp else math.atan2(rho[1], rho[0])
        p = [b.p0 * b.d[0], b.p0 * b.d[1]]
        Lz = rho[0] * p[1] - rho[1] * p[0]
        self.sigma = 1 if Lz > 0 else -1
        self.L = abs(Lz)
        self.pr0 = (p[0] * rho[0] + p[1] * rho[1]) / self.r0
        self.W = b.gamma0 * b.m * b.c ** 2 + self.kappa / self.r0          # eq:energy
        self.u0 = 1 / self.r0
        if self.L < 1e-12:  # radial motion (eq:radial-motion): no orbit in theta
            return
        self.G2 = 1 - self.kappa ** 2 / (self.L ** 2 * b.c ** 2)           # eq:orbit
        self.B = -self.kappa * self.W / (b.c ** 2 * self.L ** 2 * self.G2)  # eq:general
        self.du0 = -self.pr0 / self.L                                       # eq:radial

    def u(self, th):
        lib, g2 = self.lib, self.G2
        if g2 > 0:
            g = g2 ** 0.5
            return self.B + (self.u0 - self.B) * lib.cos(g * th) + self.du0 * lib.sin(g * th) / g
        g = (-g2) ** 0.5
        return self.B + (self.u0 - self.B) * lib.cosh(g * th) + self.du0 * lib.sinh(g * th) / g

    def du(self, th):
        lib, g2 = self.lib, self.G2
        if g2 > 0:
            g = g2 ** 0.5
            return -(self.u0 - self.B) * g * lib.sin(g * th) + self.du0 * lib.cos(g * th)
        g = (-g2) ** 0.5
        return (self.u0 - self.B) * g * lib.sinh(g * th) + self.du0 * lib.cosh(g * th)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.rc[0] + r * self.lib.cos(phi), self.rc[1] + r * self.lib.sin(phi))

    def dtdth(self, th):                                                    # eq:dtdth
        u = self.u(th)
        return (self.W - self.kappa * u) / (self.b.c ** 2 * self.L * u * u)

    def momentum(self, th):                                                 # eq:momentum
        lib = self.lib
        phi = self.phi0 + self.sigma * th
        pr = -self.L * self.du(th)
        pphi = self.sigma * self.L * self.u(th)
        return (pr * lib.cos(phi) - pphi * lib.sin(phi), pr * lib.sin(phi) + pphi * lib.cos(phi)), pr, pphi


def segment_box(x1, y1, x2, y2, lo, hi):
    """Where the chords (x1,y1)->(x2,y2) enter the box [lo, hi] (fraction in [0, 1]), or inf."""
    t0 = np.zeros_like(x1)
    t1 = np.ones_like(x1)
    for p1, p2, l, h in ((x1, x2, lo[0], hi[0]), (y1, y2, lo[1], hi[1])):
        dp = p2 - p1
        with np.errstate(divide="ignore", invalid="ignore"):
            ta = np.where(dp != 0, (l - p1) / dp, np.where((p1 >= l) & (p1 <= h), -np.inf, np.inf))
            tb = np.where(dp != 0, (h - p1) / dp, np.where((p1 >= l) & (p1 <= h), np.inf, -np.inf))
        t0 = np.maximum(t0, np.minimum(ta, tb))
        t1 = np.minimum(t1, np.maximum(ta, tb))
    return np.where(t0 <= t1, t0, np.inf)


def segment_circle(x1, y1, x2, y2, cx, cy, rad):
    """Where the chords first come within rad of (cx, cy) (fraction in [0, 1]), or inf."""
    dx, dy = x2 - x1, y2 - y1
    fx, fy = x1 - cx, y1 - cy
    A = dx * dx + dy * dy
    Bq = 2 * (fx * dx + fy * dy)
    C = fx * fx + fy * fy - rad * rad
    disc = Bq * Bq - 4 * A * C
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (-Bq - np.sqrt(np.maximum(disc, 0))) / (2 * A)
    inside0 = C <= 0
    return np.where(inside0, 0.0, np.where((disc >= 0) & (t >= 0) & (t <= 1), t, np.inf))


def first_event(orb, b, dth=4e-4):
    """The first event along the exact orbit, from chords of `dth` (a few 1e-7 cells off the
    arc): (code, theta interval, face) with code A, L, C or T (Foundations 5)."""
    if orb.L < 1e-12:  # on the launch line: radial motion (eq:radial-motion)
        if orb.pr0 < 0:
            if orb.kappa <= 0:
                return "C", None, None
            r_turn = orb.kappa / (orb.W - b.m * b.c ** 2)
            return ("C" if r_turn <= b.R + b.a else "L"), None, None
        return "L", None, None
    th0, t0 = 0.0, 0.0
    while True:
        th = th0 + dth * np.arange(0, 20001)
        u = orb.u(th)
        escaped = u <= 0
        if escaped.any():
            k = int(np.argmax(escaped))
            th, u = th[:max(k, 2)], u[:max(k, 2)]
        x, y = orb.pos(th)
        dtd = orb.dtdth(th)
        t = t0 + np.concatenate([[0.0], np.cumsum((dtd[1:] + dtd[:-1]) * 0.5 * dth)])
        x1, y1, x2, y2 = x[:-1], y[:-1], x[1:], y[1:]
        f_det = segment_box(x1, y1, x2, y2, b.dmin, b.dmax)
        f_col = segment_circle(x1, y1, x2, y2, orb.rc[0], orb.rc[1], b.R + b.a)
        out = (x2 < b.wmin[0]) | (x2 > b.wmax[0]) | (y2 < b.wmin[1]) | (y2 > b.wmax[1])
        f_out = np.where(out, 1.0, np.inf)
        f_late = np.where(t[1:] > b.tmax, 1.0, np.inf)
        first = np.minimum.reduce([f_det, f_col, f_out, f_late])
        hit = np.isfinite(first)
        if hit.any():
            k = int(np.argmax(hit))
            # Within one chord: the earliest; ties ranked obstacle, world, detector.
            cands = [(f_col[k], 0, "C"), (f_out[k], 1, "L"), (f_det[k], 2, "A"), (f_late[k], 3, "T")]
            f, _, code = min(cands)
            return code, (float(th[k]), float(th[k + 1])), (x1[k], y1[k], x2[k], y2[k])
        if escaped.any():
            return "L", None, None
        th0, t0 = float(th[-1]), float(t[-1])


def entry_exact(orb, b, bracket):
    """The exact entry into the detector within a chord's theta interval: the face crossed
    and the theta where the orbit meets it (Newton's method, eq:newton of the Foundations)."""
    lo, hi = mp.mpf(bracket[0]), mp.mpf(bracket[1])
    faces = [(0, b.dmin[0]), (0, b.dmax[0]), (1, b.dmin[1]), (1, b.dmax[1])]
    best = None
    for axis, value in faces:
        f = lambda th: orb.pos(th)[axis] - value
        if f(lo) * f(hi) > 0:
            continue
        th = mp.findroot(f, (lo, hi), solver="anderson")
        p = orb.pos(th)
        other = 1 - axis
        if b.dmin[other] - mp.mpf("1e-20") <= p[other] <= b.dmax[other] + mp.mpf("1e-20"):
            if best is None or th < best[0]:
                best = (th, axis, value)
    if best is None:
        raise ValueError("no face crossed in the bracket")
    return best


def compute(doc):
    mp.mp.dps = DIGITS
    lvl = doc.level
    b = Board(lvl, mp.mpf)
    ref = lvl["reference_solution"][0]
    s = lvl["grid"]["subdivision"]
    Q = mp.mpf(ref["value"])
    rc = [mp.mpf(v) / s for v in ref["node"][:2]]

    # ---- Setup -------------------------------------------------------------------------
    doc.value("q", b.q, 3, sci=True)
    doc.value("m", b.m, 3)
    doc.value("c", b.c, 3)
    doc.value("Tzero", b.T0, 3)
    doc.value("R", b.R, 3)
    doc.value("Q", Q, 3, sci=True)
    doc.value("tmax", b.tmax, 4)
    orb = Orbit(b, Q, rc, mp)
    doc.value("kappa", orb.kappa, 6)

    # ---- Step 1: the launch (Foundations eq:launch) -------------------------------------------
    doc.value("gammaZero", b.gamma0, 10)
    doc.value("pZero", b.p0, 15)
    v0 = b.p0 / (b.gamma0 * b.m)
    doc.value("vZero", v0, 10)
    doc.value("betaZero", v0 / b.c, 6)
    doc.value("vNewton", (2 * b.T0 / b.m) ** 0.5, 6)

    # ---- Step 2: the conserved quantities (eq:energy, eq:angmom) -------------------------
    doc.value("rZero", orb.r0, 6)
    doc.value("W", orb.W, 15)
    doc.value("L", orb.L, 15)
    doc.check(orb.sigma < 0, "the particle turns clockwise about the charge (L < 0)")
    doc.check(abs(orb.pr0) < mp.mpf("1e-35"), "the launch is perpendicular to the radius")

    # ---- Step 3: the orbit's constants (eq:orbit, eq:general, eq:amplitude) --------------
    eta = orb.kappa / (b.c * orb.L)
    doc.value("eta", eta, 10)
    G = orb.G2 ** 0.5
    doc.value("Gamma", G, 15)
    doc.value("B", orb.B, 15)
    A = orb.u0 - orb.B
    A_energy = (orb.B ** 2 + (orb.W ** 2 - b.m ** 2 * b.c ** 4) / (b.c ** 2 * orb.L ** 2 * orb.G2)) ** 0.5
    doc.value("A", A, 15)
    doc.value("Adiff", abs(A - A_energy), 2)
    doc.check(abs(A - A_energy) < mp.mpf("1e-35"), "A from the launch equals A from the energy (eq:amplitude)")
    doc.check(orb.B < 0 < A + orb.B, "B < 0 < A + B: a repulsive, unbound orbit")
    psi_inf = mp.acos(-orb.B / A) / G
    doc.value("psiInf", psi_inf, 10)
    doc.value("psiInfDeg", psi_inf * 180 / mp.pi, 6)

    # ---- Step 4: crossing x = x_D (Newton's method on eq:crossing) -------------------------
    xD = b.dmin[0]
    doc.value("xD", xD, 4)
    f = lambda psi: mp.sin(psi) - xD * (A * mp.cos(G * psi) + orb.B)
    fp = lambda psi: mp.cos(psi) + xD * A * G * mp.sin(G * psi)
    psi = mp.mpf(1)
    rows = []
    for k in range(8):
        val = f(psi)
        rows.append([str(k), psi, val])
        if abs(val) < mp.mpf("1e-30"):
            break
        psi = psi - val / fp(psi)
    doc.table("newton", [r"$k$", r"$\theta_k$", r"$f(\theta_k)$"], rows, align="rll",
              digits=[None, 22, 3], sci=[False, False, True])
    psi_e = mp.findroot(f, psi)
    doc.value("psiE", psi_e, 20)
    doc.value("psiEdeg", psi_e * 180 / mp.pi, 8)
    xe, ye = orb.pos(psi_e)
    ue = orb.u(psi_e)
    doc.value("uE", ue, 15)
    doc.value("rE", 1 / ue, 15)
    doc.value("yE", ye, 15)
    doc.check(abs(xe - xD) < mp.mpf("1e-35"), "the crossing is at x = x_D")
    doc.check(b.dmin[1] <= ye <= b.dmax[1], "y_e lies on the detector's left face")
    doc.check(ye < b.wmax[1], "the flight stays below the world's top edge")

    # ---- Step 5: the time (eq:time, eq:i1, eq:i2) --------------------------------------------
    sE = G * psi_e
    k2 = A ** 2 - orb.B ** 2
    I1 = 2 / mp.sqrt(k2) * mp.atanh(mp.sqrt((A - orb.B) / (A + orb.B)) * mp.tan(sE / 2))
    I2 = A * mp.sin(sE) / (k2 * (A * mp.cos(sE) + orb.B)) - orb.B / k2 * I1
    t_e = (orb.W * I2 - orb.kappa * I1) / (b.c ** 2 * orb.L * G)
    t_quad = mp.quad(orb.dtdth, [0, psi_e])
    doc.value("sE", sE, 15)
    doc.value("IOne", I1, 15)
    doc.value("ITwo", I2, 15)
    doc.value("tE", t_e, 15)
    doc.value("tEquadDiff", abs(t_e - t_quad), 2)
    doc.check(abs(t_e - t_quad) < mp.mpf("1e-30"), "the closed form of the time equals its quadrature")

    # ---- Step 6: the momentum at the entry (eq:momentum) ----------------------------------------
    gamma_e = (orb.W - orb.kappa * ue) / (b.m * b.c ** 2)
    pe_abs = b.m * b.c * mp.sqrt(gamma_e ** 2 - 1)
    (pxe, pye), pre, pphie = orb.momentum(psi_e)
    doc.value("gammaE", gamma_e, 15)
    doc.value("pE", pe_abs, 15)
    doc.value("prE", pre, 15)
    doc.value("pphiE", pphie, 15)
    doc.value("pxE", pxe, 15)
    doc.value("pyE", pye, 15)
    doc.value("TE", (gamma_e - 1) * b.m * b.c ** 2, 10)
    doc.value("angleE", mp.atan2(pye, pxe) * 180 / mp.pi, 8)
    doc.check(abs(mp.sqrt(pxe ** 2 + pye ** 2) - pe_abs) < mp.mpf("1e-35"),
              "the momentum's components agree with its magnitude from the energy")

    # ---- Step 7: the radiated energy (Foundations eq:larmor) --------------------------------------
    def power_dt(psi):
        x, y = orb.pos(psi)
        rx, ry = x - rc[0], y - rc[1]
        r3 = (rx * rx + ry * ry) ** mp.mpf(1.5)
        F = (orb.kappa * rx / r3, orb.kappa * ry / r3)
        (px, py), _, _ = orb.momentum(psi)
        gam = (orb.W - orb.kappa * orb.u(psi)) / (b.m * b.c ** 2)
        v = (px / (gam * b.m), py / (gam * b.m))
        vF = v[0] * F[0] + v[1] * F[1]
        acc = ((F[0] - v[0] * vF / b.c ** 2) / (gam * b.m), (F[1] - v[1] * vF / b.c ** 2) / (gam * b.m))
        a2 = acc[0] ** 2 + acc[1] ** 2
        vxa = v[0] * acc[1] - v[1] * acc[0]
        P = mp.mpf(2) / 3 * b.q ** 2 / b.c ** 3 * gam ** 6 * (a2 - vxa ** 2 / b.c ** 2)
        return P * orb.dtdth(psi)

    E_rad = mp.quad(power_dt, [0, psi_e / 2, psi_e])
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / b.T0, 3, sci=True)

    # ---- The comparison with the game -------------------------------------------------------------
    g = doc.flight(0)
    doc.compare("Outcome", "outcome", "Arrived", g["outcome"])
    doc.compare("T", r"time of arrival $t_\mathrm{e}$", t_e, g["end"]["t"], TOL)
    doc.compare("X", r"$x_\mathrm{e}$", xe, g["end"]["x"][0], TOL)
    doc.compare("Y", r"$y_\mathrm{e}$", ye, g["end"]["x"][1], TOL)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", pxe, g["end"]["p"][0], TOL)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", pye, g["end"]["p"][1], TOL)
    doc.compare("Rad", r"radiated energy (rel.)", E_rad, g["radiated_energy"], TOL_RAD, rel=True, digits=8)

    # ---- The reference orbit, for the figure ------------------------------------------------------
    doc.data("orbit", [[float(v) for v in orb.pos(psi_e * k / 200)] for k in range(201)])

    # ---- Other solutions: every single-charge placement (section 7) -------------------------------
    bf = Board(lvl, float)
    codes = []
    solving = {}
    for placement in doc.family():
        e = placement[0]
        node, Qv = e["node"], e["value"]
        rcf = [node[0] / s, node[1] / s]
        dist = math.hypot(bf.x0[0] - rcf[0], bf.x0[1] - rcf[1])
        if dist <= bf.R + bf.a:
            codes.append("X")   # the launch inside the charge: not a legal placement
            continue
        code, _, _ = first_event(Orbit(bf, Qv, rcf, np), bf)
        codes.append(code)
        if code == "A":
            solving.setdefault(Qv, []).append(rcf)
    ours = "".join(codes)
    doc.compare_verdicts("Map", "every single-charge placement: first event", ours)
    values = doc.family()
    distinct = sorted({p[0]["value"] for p in values}, key=lambda v: (v < 0, abs(v)))
    names = {1e6: "PlusOne", 2e6: "PlusTwo", 4e6: "PlusFour",
             -1e6: "MinusOne", -2e6: "MinusTwo", -4e6: "MinusFour"}
    for Qv in distinct:
        pts = solving.get(Qv, [])
        doc.value(f"n{names[Qv]}", len(pts))
        doc.data(f"map_{names[Qv].lower()}", pts if pts else [[-5, -5]])
    doc.value("nSolving", sum(len(v) for v in solving.values()))
    doc.value("nFamily", len(values))
    doc.value("nNodes", len(values) // len(distinct))

    # ---- Three more placements in detail (placements.json) ----------------------------------------
    for k in range(1, doc.n_results()):
        e = doc.result(k)["placement"][0]
        rck = [mp.mpf(v) / s for v in e["node"][:2]]
        ok = Orbit(b, mp.mpf(e["value"]), rck, mp)
        code, bracket, _ = first_event(Orbit(bf, e["value"], [float(v) for v in rck], np), bf)
        gk = doc.flight(k)
        tag = ["", "One", "Two", "Three", "Four", "Five"][k]
        doc.text(f"det{tag}Node", f"({e['node'][0]}, {e['node'][1]})")
        doc.value(f"det{tag}Q", mp.mpf(e["value"]), 3, sci=True)
        if code != "A":
            doc.compare(f"Det{tag}", f"placement {k}: outcome", {"L": "LeftBounds", "C": "Collided",
                        "T": "Timeout"}[code], gk["outcome"].split("(")[0])
            continue
        th, axis, value = entry_exact(ok, b, bracket)
        p = ok.pos(th)
        t_k = mp.quad(ok.dtdth, [0, th])
        if k == 1:
            # How thin the corner is: the gap to the right face, and the angle the orbit
            # spends inside the detector before it would leave through that face.
            doc.value("cornerGap", b.dmax[0] - p[0], 2)
            th_out = mp.findroot(lambda t: ok.pos(t)[0] - b.dmax[0], th + mp.mpf("1e-5"))
            doc.value("sliverAngle", th_out - th, 2, sci=True)
            doc.value("sliverDepth", max(ok.pos(th + (th_out - th) * j / 50)[1] for j in range(51)) - b.dmin[1], 2, sci=True)
        doc.value(f"det{tag}T", t_k, 12)
        doc.value(f"det{tag}X", p[0], 12)
        doc.value(f"det{tag}Y", p[1], 12)
        doc.text(f"det{tag}Face", ["left", "right", "bottom", "top"][2 * axis + (value == b.dmax[axis])])
        doc.compare(f"Det{tag}T", rf"placement {k}: $t_\mathrm{{e}}$", t_k, gk["end"]["t"], TOL)
        doc.compare(f"Det{tag}X", rf"placement {k}: $x_\mathrm{{e}}$", p[0], gk["end"]["x"][0], TOL)
        doc.compare(f"Det{tag}Y", rf"placement {k}: $y_\mathrm{{e}}$", p[1], gk["end"]["x"][1], TOL)
