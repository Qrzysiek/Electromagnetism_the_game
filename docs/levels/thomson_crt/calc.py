"""Companion of docs/levels/thomson_crt/doc.tex (Level 7, Thomson's cathode-ray tube).

One repulsive charge and two launch energies. Each flight is Level 1's orbit
(first_bend eq:general). The launch is not an apsis, so the time and the radiated
energy are Simpson's rule (Foundations eq:simpson) on d t/d theta (first_bend eq:dtdth).
The face is polished by Newton's method (Foundations eq:newton).

The Newtonian impulse from infinity, theta = kappa/(T b), is an estimate of the
deflection. It is computed here so the document can quote how far it misses the
arrival angle. The comparison uses the orbit at the face.

The step-doubling error of the time quadrature is required to be below 1e-12 before
any comparison, so it does not set the tolerance. A verified flight is accurate to
about 1e-10 (Foundations, section 6). Ten times that.
"""
import math

import mpmath as mp
import numpy as np


TOL = 1e-9
# Relative tolerance on the radiated energy. The quadrature's own relative error is
# checked below, and is far smaller; 1e-4 is the room used for this integral before.
TOL_RAD = 1e-4
N_SIMP = 20000


class Shot:
    def __init__(self, lvl, shot, mpf):
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
        nrm = (d[0] ** 2 + d[1] ** 2) ** 0.5
        self.d = [d[0] / nrm, d[1] / nrm]
        self.dmin = [mpf(v) / s for v in shot["detector"]["min"][:2]]
        self.dmax = [mpf(v) / s for v in shot["detector"]["max"][:2]]
        self.wmin = [mpf(-1), mpf(-1)]
        self.wmax = [mpf(lvl["grid"]["nx"]) / s + 1, mpf(lvl["grid"]["ny"]) / s + 1]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = (self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0) ** 0.5


class Orbit:
    """Level 1's orbit about one charge (first_bend eq:general)."""

    def __init__(self, shot, Q, rc, lib):
        self.s, self.lib, self.rc = shot, lib, rc
        self.kappa = shot.q * Q
        rho = [shot.x0[0] - rc[0], shot.x0[1] - rc[1]]
        self.r0 = (rho[0] ** 2 + rho[1] ** 2) ** 0.5
        self.phi0 = lib.atan2(rho[1], rho[0])
        p = [shot.p0 * shot.d[0], shot.p0 * shot.d[1]]
        Lz = rho[0] * p[1] - rho[1] * p[0]
        self.sigma = 1 if Lz > 0 else -1
        self.L = abs(Lz)
        self.pr0 = (p[0] * rho[0] + p[1] * rho[1]) / self.r0
        self.W = shot.gamma0 * shot.m * shot.c ** 2 + self.kappa / self.r0
        self.u0 = 1 / self.r0
        self.G2 = 1 - self.kappa ** 2 / (self.L ** 2 * shot.c ** 2)
        self.B = -self.kappa * self.W / (shot.c ** 2 * self.L ** 2 * self.G2)
        self.du0 = -self.pr0 / self.L

    def u(self, th):
        g = self.G2 ** 0.5
        return self.B + (self.u0 - self.B) * self.lib.cos(g * th) + self.du0 * self.lib.sin(g * th) / g

    def du(self, th):
        g = self.G2 ** 0.5
        return -(self.u0 - self.B) * g * self.lib.sin(g * th) + self.du0 * self.lib.cos(g * th)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.rc[0] + r * self.lib.cos(phi), self.rc[1] + r * self.lib.sin(phi))

    def dtdth(self, th):
        u = self.u(th)
        return (self.W - self.kappa * u) / (self.s.c ** 2 * self.L * u * u)

    def momentum(self, th):
        phi = self.phi0 + self.sigma * th
        pr = -self.L * self.du(th)
        pp = self.sigma * self.L * self.u(th)
        return (pr * self.lib.cos(phi) - pp * self.lib.sin(phi),
                pr * self.lib.sin(phi) + pp * self.lib.cos(phi))

    def dx(self, th):
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        return -self.lib.sin(phi) * self.sigma / u - self.lib.cos(phi) * du / (u * u)

    def dy(self, th):
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        return self.lib.cos(phi) * self.sigma / u - self.lib.sin(phi) * du / (u * u)


def newton_face(orb, axis, value, th0, n=12):
    deriv = orb.dx if axis == 0 else orb.dy
    th, rows = th0, []
    for k in range(n):
        f = orb.pos(th)[axis] - value
        rows.append((k, th, f))
        if abs(f) < mp.mpf("1e-28"):
            break
        th = th - f / deriv(th)
    return th, rows


def simpson(fn, a, b, n):
    h = (b - a) / n
    xs = [a + h * k for k in range(n + 1)]
    vals = [float(fn(x)) for x in xs]
    acc = vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2])
    return h / 3 * acc


def scan(orb, shot, dth=1e-3):
    th = 0.0
    while th < 4:
        ths = th + dth * np.arange(0, 4001)
        u = np.array([float(orb.u(t)) for t in ths])
        pos = [orb.pos(t) for t in ths]
        xs = np.array([float(p[0]) for p in pos])
        ys = np.array([float(p[1]) for p in pos])
        rs = 1 / np.maximum(u, 1e-15)
        for k in range(len(ths) - 1):
            if rs[k + 1] <= float(shot.R + shot.a):
                return "Collided", mp.mpf(ths[k]), None
            out = (xs[k + 1] < float(shot.wmin[0]) or xs[k + 1] > float(shot.wmax[0])
                   or ys[k + 1] < float(shot.wmin[1]) or ys[k + 1] > float(shot.wmax[1]))
            inside = (float(shot.dmin[0]) <= xs[k + 1] <= float(shot.dmax[0])
                      and float(shot.dmin[1]) <= ys[k + 1] <= float(shot.dmax[1]))
            if out:
                return "LeftBounds", mp.mpf(ths[k]), _crossed_face(
                    xs[k], ys[k], xs[k + 1], ys[k + 1], shot.wmin, shot.wmax)
            if inside:
                return "Arrived", mp.mpf(ths[k]), _crossed_face(
                    xs[k], ys[k], xs[k + 1], ys[k + 1], shot.dmin, shot.dmax)
        th = float(ths[-1])
    return "Timeout", None, None


def _crossed_face(x1, y1, x2, y2, lo, hi):
    tests = [(0, float(lo[0]), x1, x2), (0, float(hi[0]), x1, x2),
             (1, float(lo[1]), y1, y2), (1, float(hi[1]), y1, y2)]
    best = None
    for axis, value, p1, p2 in tests:
        if (p1 - value) * (p2 - value) <= 0 and p1 != p2:
            frac = (value - p1) / (p2 - p1)
            if 0 <= frac <= 1 and (best is None or frac < best[0]):
                best = (frac, axis, value)
    return (best[1], mp.mpf(best[2]))


def radiated(orb, shot, rc, th):
    def power_dt(psi):
        xx, yy = orb.pos(psi)
        rx, ry = xx - rc[0], yy - rc[1]
        r3 = (rx * rx + ry * ry) ** mp.mpf("1.5")
        F = (orb.kappa * rx / r3, orb.kappa * ry / r3)
        pxx, pyy = orb.momentum(psi)
        gam = (orb.W - orb.kappa * orb.u(psi)) / (shot.m * shot.c ** 2)
        v = (pxx / (gam * shot.m), pyy / (gam * shot.m))
        vF = v[0] * F[0] + v[1] * F[1]
        ax = (F[0] - v[0] * vF / shot.c ** 2) / (gam * shot.m)
        ay = (F[1] - v[1] * vF / shot.c ** 2) / (gam * shot.m)
        vxa = v[0] * ay - v[1] * ax
        P = (mp.mpf(2) / 3 * shot.q ** 2 / shot.c ** 3 * gam ** 6
             * (ax * ax + ay * ay - vxa ** 2 / shot.c ** 2))
        return P * orb.dtdth(psi)

    e1 = simpson(power_dt, 0, th, 4000)
    e2 = simpson(power_dt, 0, th, 8000)
    return e2, abs(e2 - e1) / abs(e2) / 15


def closest(orb, th_end):
    g = mp.sqrt(orb.G2)
    th = mp.atan(orb.du0 / (g * (orb.u0 - orb.B))) / g
    cands = [mp.mpf(0), th_end]
    for base in (th, th + mp.pi / g, th - mp.pi / g):
        if mp.mpf(0) <= base <= th_end:
            cands.append(base)
    best = max(cands, key=lambda t: orb.u(t))
    return best, 1 / orb.u(best)


def solve(shot, shot_f, Q, rc):
    orb = Orbit(shot, Q, rc, mp)
    code, th_bracket, face = scan(Orbit(shot_f, float(Q), [float(c) for c in rc], math), shot_f)
    th, rows = newton_face(orb, face[0], face[1], th_bracket + mp.mpf("5e-4"))
    residual = abs(orb.pos(th)[face[0]] - face[1])
    x, y = orb.pos(th)
    if face[0] == 0:
        x = face[1]
    else:
        y = face[1]
    px, py = orb.momentum(th)
    t1 = simpson(orb.dtdth, 0, th, N_SIMP)
    t2 = simpson(orb.dtdth, 0, th, 2 * N_SIMP)
    rad, raderr = radiated(orb, shot, rc, th)
    th_min, r_min = closest(orb, th)
    samples = []
    for k in range(81):
        if k == 0:
            samples.append((float(shot.x0[0]), float(shot.x0[1])))
        else:
            p = orb.pos(th * k / 80)
            samples.append((float(p[0]), float(p[1])))
    return {
        "orb": orb, "code": code, "theta": th, "newton": rows, "face": face,
        "residual": residual, "guess": th_bracket + mp.mpf("5e-4"),
        "x": x, "y": y, "px": px, "py": py,
        "t": t2, "terr": abs(t2 - t1) / 15,
        "rad": rad, "raderr": raderr, "path": samples,
        "rmin": r_min, "thmin": th_min,
        "ang": mp.atan2(py, px),
        "T": (mp.sqrt(1 + (px ** 2 + py ** 2) / (shot.m ** 2 * shot.c ** 2)) - 1) * shot.m * shot.c ** 2,
    }


def compare_flight(doc, tag, what, sol, k, flight):
    g = doc.flight(k, flight=flight)
    doc.compare(f"{tag}Out", f"{what}: outcome", sol["code"], g["outcome"])
    doc.compare(f"{tag}T", rf"{what}: $t$", sol["t"], g["end"]["t"], TOL)
    doc.compare(f"{tag}X", rf"{what}: $x$", sol["x"], g["end"]["x"][0], TOL)
    doc.compare(f"{tag}Y", rf"{what}: $y$", sol["y"], g["end"]["x"][1], TOL)
    doc.compare(f"{tag}Px", rf"{what}: $p_x$", sol["px"], g["end"]["p"][0], TOL)
    doc.compare(f"{tag}Py", rf"{what}: $p_y$", sol["py"], g["end"]["p"][1], TOL)
    doc.compare(f"{tag}Rad", f"{what}: radiated energy (rel.)",
                sol["rad"], g["radiated_energy"], TOL_RAD, rel=True)


def compute(doc):
    mp.mp.dps = 40
    lvl = doc.level
    ref = lvl["reference_solution"][0]
    sdiv = lvl["grid"]["subdivision"]
    Q = mp.mpf(ref["value"])
    rc = [mp.mpf(v) / sdiv for v in ref["node"][:2]]
    shots = [Shot(lvl, sh, mp.mpf) for sh in lvl["shots"]]
    shots_f = [Shot(lvl, sh, float) for sh in lvl["shots"]]
    slow, fast = shots
    doc.check(slow.T0 < fast.T0, "the first shot is the slower one")
    doc.check(slow.x0 == fast.x0 and slow.d == fast.d, "both shots share the launch ray")

    doc.value("q", slow.q, 3, sci=True)
    doc.value("m", slow.m, 3)
    doc.value("c", slow.c, 3)
    doc.value("R", slow.R, 3)
    doc.value("tmax", slow.tmax, 4)
    doc.value("Q", Q, 3, sci=True)
    doc.value("kappa", slow.q * Q, 4)
    doc.value("Tslow", slow.T0, 3)
    doc.value("Tfast", fast.T0, 3)
    doc.value("xLaunch", slow.x0[0], 2)
    doc.value("yLaunch", slow.x0[1], 2)
    doc.value("xC", rc[0], 2)
    doc.value("yC", rc[1], 2)
    doc.value("b", slow.x0[1] - rc[1], 4)
    doc.value("xWorld", slow.wmax[0], 2)
    doc.value("yWorld", slow.wmax[1], 2)
    doc.value("nSimp", 2 * N_SIMP)
    doc.value("nRad", 8000)
    reg = lvl["limits"]["region"]
    doc.value("xRlo", mp.mpf(reg["min"][0]) / sdiv, 2)
    doc.value("yRlo", mp.mpf(reg["min"][1]) / sdiv, 2)
    doc.value("xRhi", mp.mpf(reg["max"][0]) / sdiv, 2)
    doc.value("yRhi", mp.mpf(reg["max"][1]) / sdiv, 2)
    for shot, name, face in ((slow, "Slow", "xS"), (fast, "Fast", "xF")):
        doc.value(f"gamma{name}", shot.gamma0, 10)
        doc.value(f"p{name}", shot.p0, 12)
        doc.value(f"{face}", shot.dmin[0], 2)
        doc.value(f"y{name}Lo", shot.dmin[1], 2)
        doc.value(f"y{name}Hi", shot.dmax[1], 2)
        doc.value(f"x{name}Hi", shot.dmax[0], 2)

    solved = []
    for shot, shot_f, name in zip(shots, shots_f, ("Slow", "Fast")):
        sol = solve(shot, shot_f, Q, rc)
        solved.append(sol)
        orb = sol["orb"]
        doc.value(f"r{name}", orb.r0, 8)
        doc.value(f"phi{name}", orb.phi0, 8)
        doc.value(f"L{name}", orb.L, 10)
        doc.value(f"Lz{name}", orb.sigma * orb.L, 10)
        doc.value(f"W{name}", orb.W, 10)
        doc.value(f"Gamma{name}", mp.sqrt(orb.G2), 10)
        doc.value(f"B{name}", orb.B, 8)
        doc.value(f"u{name}", orb.u0, 8)
        doc.value(f"du{name}", orb.du0, 8)
        doc.value(f"pr{name}", orb.pr0, 8)
        doc.value(f"theta{name}", sol["theta"], 10)
        doc.value(f"t{name}", sol["t"], 12)
        doc.value(f"x{name}", sol["x"], 12)
        doc.value(f"y{name}", sol["y"], 12)
        doc.value(f"px{name}", sol["px"], 12)
        doc.value(f"py{name}", sol["py"], 12)
        doc.value(f"ang{name}", mp.degrees(sol["ang"]), 6)
        doc.value(f"Tarr{name}", sol["T"], 8)
        doc.value(f"terr{name}", sol["terr"], 2, sci=True)
        doc.value(f"rmin{name}", sol["rmin"], 6)
        doc.value(f"thMin{name}", sol["thmin"], 6)
        doc.value(f"Erad{name}", sol["rad"], 6, sci=True)
        doc.value(f"EradFrac{name}", sol["rad"] / shot.T0, 3, sci=True)
        doc.value(f"radErr{name}", sol["raderr"], 2, sci=True)
        doc.value(f"newtonRes{name}", sol["residual"], 3, sci=True)
        doc.value(f"dy{name}", sol["y"] - shot.x0[1], 6)
        # Newtonian impulse from infinity: theta = kappa / (T b), in radians.
        # b is the vertical offset of a ray along +x. Derived in the document.
        imp = orb.kappa / (shot.T0 * (shot.x0[1] - rc[1]))
        doc.value(f"imp{name}", mp.degrees(imp), 6)
        doc.value(f"impGap{name}", abs(mp.degrees(imp) - mp.degrees(sol["ang"])), 4)
        doc.check(sol["code"] == "Arrived", f"the {name} shot arrives")
        doc.check(sol["face"][0] == 0 and sol["face"][1] == shot.dmin[0],
                  f"the {name} shot meets the screen face")
        doc.check(shot.dmin[1] <= sol["y"] <= shot.dmax[1], f"the {name} ordinate is inside its spot")
        doc.check(sol["terr"] < 1e-12, f"the {name} time quadrature is finer than the tolerance")
        doc.check(sol["residual"] < mp.mpf("1e-18"), f"Newton lands on the {name} face")
        doc.check(sol["rmin"] > shot.R, f"the {name} shot misses the charge")
        doc.check(0 < sol["thmin"] < sol["theta"], f"the {name} closest approach is on the flight")
        doc.check(abs(orb.du0) > mp.mpf("1e-3"), f"the {name} launch is not an apsis")
        doc.data(f"path{name.lower()}", sol["path"])

    # The screen rises in proportion to 1/T when the turn is at one place.
    ratio_dy = (solved[0]["y"] - slow.x0[1]) / (solved[1]["y"] - fast.x0[1])
    ratio_T = fast.T0 / slow.T0
    doc.value("ratioDy", ratio_dy, 6)
    doc.value("ratioT", ratio_T, 6)
    doc.value("ratioGap", abs(ratio_dy - ratio_T), 4)
    doc.value("gGap", 1 - mp.sqrt(solved[0]["orb"].G2), 3, sci=True)

    rows = []
    for k, th, f in solved[0]["newton"]:
        rows.append([str(k), th, f])
    doc.table("newton", [r"$k$", r"$\theta_k$", r"$x(\theta_k) - x_S$"], rows,
              align="rll", digits=[None, 12, 3], sci=[False, False, True])
    doc.value("newtonSteps", len(rows) - 1)
    doc.value("thGuess", solved[0]["guess"], 6)
    samples = []
    for sol, launch in zip(solved, ("slow", "fast")):
        orb = sol["orb"]
        for frac in (0, 0.5, 1):
            th = sol["theta"] * frac
            x, y = shots[0].x0 if frac == 0 else orb.pos(th)
            # Both launches share the point; the fast shot's launch is the same point.
            if frac == 0:
                x, y = shots[solved.index(sol)].x0
            samples.append([launch, th, x, y])
    doc.table("samples", [r"shot", r"$\theta$", r"$x$", r"$y$"], samples,
              align="lrrr", digits=[None, 4, 4, 4])

    length = mp.mpf("2.2")
    doc.value("arrowX", slow.x0[0] + length, 4)
    doc.value("arrowY", slow.x0[1], 4)

    # ---- The neighbouring magnitudes, same node -------------------------------------------
    tags = ("Weak", "Hard")
    for k, tag in enumerate(tags, start=1):
        e = doc.result(k)["placement"][0]
        Qk = mp.mpf(e["value"])
        rck = [mp.mpf(v) / sdiv for v in e["node"][:2]]
        doc.text(f"{tag}Node", f"({float(rck[0]):g}, {float(rck[1]):g})")
        doc.value(f"{tag}Q", Qk, 3, sci=True)
        doc.value(f"{tag}K", slow.q * Qk, 4)
        doc.check(rck[0] == rc[0] and rck[1] == rc[1], f"the {tag} charge uses the reference node")
        for i, name in enumerate(("Slow", "Fast")):
            sol = solve(shots[i], shots_f[i], Qk, rck)
            doc.value(f"{tag}{name}T", sol["t"], 8)
            doc.value(f"{tag}{name}X", sol["x"], 8)
            doc.value(f"{tag}{name}Y", sol["y"], 8)
            doc.value(f"{tag}{name}Ang", mp.degrees(sol["ang"]), 4)
            doc.value(f"{tag}{name}Err", sol["terr"], 2, sci=True)
            doc.value(f"{tag}{name}Erad", sol["rad"], 3, sci=True)
            doc.value(f"{tag}{name}RadErr", sol["raderr"], 2, sci=True)
            doc.check(sol["code"] == "LeftBounds", f"the {tag} {name} shot leaves the world")
            doc.check(sol["terr"] < 1e-12, f"the {tag} {name} quadrature is finer than the tolerance")
            doc.check(sol["residual"] < mp.mpf("1e-18"), f"Newton lands on the {tag} {name} edge")
            doc.data(f"path{tag.lower()}{name.lower()}", sol["path"])
            compare_flight(doc, f"{tag}{name}", f"{tag.lower()} {name.lower()}", sol, k, i)

    # Reference comparisons last, after the quadrature checks above.
    compare_flight(doc, "Slow", "slow", solved[0], 0, 0)
    compare_flight(doc, "Fast", "fast", solved[1], 0, 1)
