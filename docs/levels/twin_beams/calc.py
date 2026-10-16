"""Companion of docs/levels/twin_beams/doc.tex (Level 4, Twin beams).

One repulsive charge and two launches. Each flight is Level 1's orbit (first_bend eq:general),
evaluated independently of the game. The detector face is found by Newton's method
(Foundations eq:newton); the time and the radiated energy are Simpson's rule
(Foundations eq:simpson) on d t/d theta (first_bend eq:dtdth). The step-doubling error of
that quadrature is computed before any comparison and is far below 1e-10, so the tolerance
is the exact-orbit one, 1e-9.
"""
import math

import mpmath as mp
import numpy as np


# Exact orbit against a verified flight: the game is accurate to about 1e-10
# (Foundations, section 6). Ten times that. The quadrature error is checked to be smaller.
TOL = 1e-9
TOL_RAD = 1e-4
# Simpson panels on [0, theta_e]. The error is |S_n - S_{2n}|/15 (Foundations eq:simpson).
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
        n = (d[0] ** 2 + d[1] ** 2) ** 0.5
        self.d = [d[0] / n, d[1] / n]
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
        """d x/d theta, for Newton's method on the face x = const."""
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        # x = cx + cos(phi)/u, d phi/d theta = sigma.
        return -self.lib.sin(phi) * self.sigma / u - self.lib.cos(phi) * du / (u * u)

    def dy(self, th):
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        return self.lib.cos(phi) * self.sigma / u - self.lib.sin(phi) * du / (u * u)


def newton_face(orb, axis, value, th0, n=8):
    """Newton's method on position[axis](theta) = value (Foundations eq:newton)."""
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
    """Simpson's rule on n panels (n even) of fn over [a, b]."""
    h = (b - a) / n
    xs = [a + h * k for k in range(n + 1)]
    # Float samples: the orbit is smooth and float64 is far inside the tolerance.
    vals = [float(fn(x)) for x in xs]
    acc = vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2])
    return h / 3 * acc


def scan(orb, shot, dth=1e-3):
    """First event along the orbit, from chords of dth. Returns (code, theta, face).

    code is Arrived, LeftBounds or Collided. `face` is (axis, value) of the boundary.
    """
    th = 0.0
    while th < 3:
        ths = th + dth * np.arange(0, 4001)
        u = np.array([float(orb.u(t)) for t in ths])
        pos = [orb.pos(t) for t in ths]
        xs = np.array([float(p[0]) for p in pos])
        ys = np.array([float(p[1]) for p in pos])
        rs = 1 / np.maximum(u, 1e-15)
        for k in range(len(ths) - 1):
            hit_c = rs[k + 1] <= float(shot.R + shot.a)
            out = (xs[k + 1] < float(shot.wmin[0]) or xs[k + 1] > float(shot.wmax[0])
                   or ys[k + 1] < float(shot.wmin[1]) or ys[k + 1] > float(shot.wmax[1]))
            inside = (float(shot.dmin[0]) <= xs[k + 1] <= float(shot.dmax[0])
                      and float(shot.dmin[1]) <= ys[k + 1] <= float(shot.dmax[1]))
            if hit_c:
                return "Collided", mp.mpf(ths[k]), None
            if out:
                return "LeftBounds", mp.mpf(ths[k]), _crossed_face(
                    xs[k], ys[k], xs[k + 1], ys[k + 1], shot.wmin, shot.wmax)
            if inside:
                return "Arrived", mp.mpf(ths[k]), _crossed_face(
                    xs[k], ys[k], xs[k + 1], ys[k + 1], shot.dmin, shot.dmax)
        th = float(ths[-1])
    return "Timeout", None, None


def _crossed_face(x1, y1, x2, y2, lo, hi):
    """Which side of the box the chord crosses, as (axis, value)."""
    tests = [(0, float(lo[0]), x1, x2), (0, float(hi[0]), x1, x2),
             (1, float(lo[1]), y1, y2), (1, float(hi[1]), y1, y2)]
    best = None
    for axis, value, p1, p2 in tests:
        if (p1 - value) * (p2 - value) <= 0 and p1 != p2:
            frac = (value - p1) / (p2 - p1)
            if 0 <= frac <= 1 and (best is None or frac < best[0]):
                best = (frac, axis, value)
    return (best[1], mp.mpf(best[2]))


def rutherford(kappa, T, b):
    """Newtonian scattering angle from infinity, in radians (geiger_marsden eq:chi-tan)."""
    return 2 * math.atan(float(kappa) / (2 * float(T) * float(b)))


def asymptote_chi(orb):
    """Relativistic scattering angle of this orbit, asymptote to asymptote, in radians."""
    def root(th0, step):
        th, up = th0, orb.u(th0)
        for _ in range(20000):
            th2 = th + step
            u2 = orb.u(th2)
            if up > 0 >= u2:
                inside, outside = th, th2
                for _ in range(60):
                    mid = (inside + outside) / 2
                    if orb.u(mid) > 0:
                        inside = mid
                    else:
                        outside = mid
                return (inside + outside) / 2
            th, up = th2, u2
            if abs(th) > 4:
                return None
        return None
    th_out = root(mp.mpf(0), mp.mpf("1e-3"))
    th_in = root(mp.mpf(0), mp.mpf("-1e-3"))
    return math.pi - float(th_out - th_in)


def solve(shot, shot_f, Q, rc):
    """The flight of one shot about one charge: event, state, time, radiation.

    `shot` is the high-precision launch; `shot_f` is the same launch in binary floats,
    used only to bracket the event.
    """
    orb = Orbit(shot, Q, rc, mp)
    code, th_bracket, face = scan(
        Orbit(shot_f, float(Q), [float(c) for c in rc], math), shot_f)
    th, rows = newton_face(orb, face[0], face[1], th_bracket + mp.mpf("5e-4"))
    # The bracket start can lie just before the face; Newton from slightly inside the
    # bracket is started once, and the polished root is checked to land on the face.
    x, y = orb.pos(th)
    px, py = orb.momentum(th)
    n = N_SIMP
    t1 = simpson(orb.dtdth, 0, th, n)
    t2 = simpson(orb.dtdth, 0, th, 2 * n)
    terr = abs(t2 - t1) / 15

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
        P = mp.mpf(2) / 3 * shot.q ** 2 / shot.c ** 3 * gam ** 6 * (ax * ax + ay * ay - vxa ** 2 / shot.c ** 2)
        return P * orb.dtdth(psi)

    # Radiation on a coarser panel count: the integral is smooth and tiny.
    e1 = simpson(power_dt, 0, th, 4000)
    e2 = simpson(power_dt, 0, th, 8000)
    samples = [orb.pos(th * k / 80) for k in range(81)]
    # Closest approach: the unique maximum of u, where u' = 0 (doc eq:uprime).
    g = mp.sqrt(orb.G2)
    th_min = mp.atan(orb.du0 / (g * (orb.u0 - orb.B))) / g
    r_min = 1 / orb.u(th_min)
    rs = [1 / orb.u(th * k / 400) for k in range(401)]
    return {
        "orb": orb, "code": code, "theta": th, "newton": rows, "face": face,
        "x": x, "y": y, "px": px, "py": py, "t": t2, "terr": terr,
        "rad": e2, "raderr": abs(e2 - e1) / abs(e2) / 15,
        "path": [(float(p[0]), float(p[1])) for p in samples],
        "rmin": r_min, "thmin": th_min,
        "rminGap": abs(float(r_min) - min(float(r) for r in rs)),
        "chi": asymptote_chi(orb),
    }


def angle_deg(px, py):
    return math.degrees(math.atan2(float(py), float(px)))


def kinetic(shot, px, py):
    g = mp.sqrt(1 + (px ** 2 + py ** 2) / (shot.m ** 2 * shot.c ** 2))
    return (g - 1) * shot.m * shot.c ** 2


def compute(doc):
    mp.mp.dps = 40
    lvl = doc.level
    ref = lvl["reference_solution"][0]
    sdiv = lvl["grid"]["subdivision"]
    Q = mp.mpf(ref["value"])
    rc = [mp.mpf(v) / sdiv for v in ref["node"][:2]]
    shots = [Shot(lvl, sh, mp.mpf) for sh in lvl["shots"]]
    shots_f = [Shot(lvl, sh, float) for sh in lvl["shots"]]
    up, down = shots

    doc.value("q", up.q, 3, sci=True)
    doc.value("m", up.m, 3)
    doc.value("c", up.c, 3)
    doc.value("Tzero", up.T0, 3)
    doc.value("R", up.R, 3)
    doc.value("tmax", up.tmax, 4)
    doc.value("Q", Q, 3, sci=True)
    doc.value("kappa", up.q * Q, 4)
    doc.text("rC", f"({float(rc[0]):g}, {float(rc[1]):g})")
    doc.value("gammaZero", up.gamma0, 10)
    doc.value("pZero", up.p0, 12)
    doc.value("nSimp", 2 * N_SIMP)
    doc.value("xFace", up.dmin[0], 2)

    names = ("Up", "Down")
    solved = []
    for shot, name in zip(shots, names):
        sol = solve(shot, shots_f[names.index(name)], Q, rc)
        solved.append(sol)
        orb = sol["orb"]
        doc.value(f"r{name}", orb.r0, 8)
        doc.value(f"L{name}", orb.L, 10)
        doc.value(f"W{name}", orb.W, 10)
        doc.value(f"Gamma{name}", mp.sqrt(orb.G2), 10)
        doc.value(f"B{name}", orb.B, 8)
        doc.value(f"u{name}", orb.u0, 8)
        doc.value(f"du{name}", orb.du0, 8)
        doc.value(f"b{name}", abs(float(rc[1] - shot.x0[1])), 4)
        chi_f = rutherford(orb.kappa, shot.T0, abs(float(rc[1] - shot.x0[1])))
        doc.value(f"chiForm{name}", math.degrees(chi_f), 6)
        doc.value(f"chiRel{name}", math.degrees(sol["chi"]), 6)
        doc.value(f"chiGap{name}", abs(math.degrees(chi_f) - math.degrees(sol["chi"])), 4)
        doc.value(f"angLeft{name}", abs(math.degrees(sol["chi"]) - abs(angle_deg(sol["px"], sol["py"]))), 4)
        # Launch potential over T_0, and how far Gamma falls short of 1: the two reasons
        # the infinity formula (geiger_marsden eq:chi-tan) misses chiRel. Computed before
        # any comparison.
        doc.value(f"pot{name}", abs(orb.kappa) / (orb.r0 * shot.T0), 4)
        doc.value(f"gGap{name}", 1 - mp.sqrt(orb.G2), 3, sci=True)
        doc.value(f"theta{name}", sol["theta"], 10)
        doc.value(f"t{name}", sol["t"], 12)
        doc.value(f"x{name}", sol["x"], 12)
        doc.value(f"y{name}", sol["y"], 12)
        doc.value(f"px{name}", sol["px"], 12)
        doc.value(f"py{name}", sol["py"], 12)
        doc.value(f"ang{name}", angle_deg(sol["px"], sol["py"]), 6)
        doc.value(f"T{name}", kinetic(shot, sol["px"], sol["py"]), 8)
        doc.value(f"terr{name}", sol["terr"], 2, sci=True)
        doc.value(f"rmin{name}", sol["rmin"], 6)
        doc.value(f"thMin{name}", sol["thmin"], 6)
        doc.value(f"rminGap{name}", sol["rminGap"], 2, sci=True)
        doc.check(0 < float(sol["thmin"]) < float(sol["theta"]),
                  f"the {name} closest approach lies on the flight")
        doc.value(f"Erad{name}", sol["rad"], 6, sci=True)
        doc.value(f"EradFrac{name}", sol["rad"] / float(shot.T0), 3, sci=True)
        doc.value(f"radErr{name}", sol["raderr"], 2, sci=True)
        doc.check(sol["code"] == "Arrived", f"the {name} shot arrives")
        doc.check(sol["terr"] < 1e-12, f"the {name} time quadrature is finer than 1e-9")
        doc.check(sol["rmin"] > float(shot.R), f"the {name} shot misses the charge")
        doc.check(abs(float(sol["x"]) - float(shot.dmin[0])) < 1e-9, f"the {name} shot meets x = 27")
        doc.check(abs(orb.L - abs(rc[1] - shot.x0[1]) * shot.p0) < mp.mpf("1e-9"),
                  f"the {name} angular momentum is the offset times the launch momentum")
        doc.data(f"path{name.lower()}", sol["path"])

    # Newton's method on the upper shot, for the hand procedure.
    rows = []
    for k, th, f in solved[0]["newton"]:
        rows.append([str(k), th, f])
    doc.table("newton", [r"$k$", r"$\theta_k$", r"$x(\theta_k) - 27$"], rows,
              align="rll", digits=[None, 12, 3], sci=[False, False, True])
    doc.value("newtonSteps", len(rows) - 1)

    consts = []
    for sol, shot in zip(solved, shots):
        orb = sol["orb"]
        launch = f"({float(shot.x0[0]):g}, {float(shot.x0[1]):g})"
        consts.append([launch, abs(float(rc[1] - shot.x0[1])),
                       orb.L, mp.sqrt(orb.G2), sol["theta"], sol["t"], sol["y"]])
    doc.table("orbits", [r"launch", r"$b$", r"$L$", r"$\Gamma$", r"$\theta_\mathrm{e}$",
                         r"$t_\mathrm{e}$", r"$y_\mathrm{e}$"],
              consts, align="lrrrrrr", digits=[None, 2, 5, 6, 6, 5, 5])
    samples = []
    for sol, shot, launch in zip(solved, shots, ("upper", "lower")):
        orb = sol["orb"]
        for frac in (0, 0.25, 0.5, 0.75, 1):
            th = sol["theta"] * frac
            # At theta = 0 the position is the launch. The sampled formula leaves a
            # roundoff of about 10^{-40} in x, which a table should not print.
            x, y = shot.x0 if frac == 0 else orb.pos(th)
            samples.append([launch, th, x, y])
    doc.table("samples", [r"shot", r"$\theta$", r"$x$", r"$y$"], samples,
              align="lrrr", digits=[None, 4, 4, 4])

    # ---- Comparison: both reference shots -------------------------------------------------
    for i, name, what in ((0, "Up", "upper"), (1, "Down", "lower")):
        g = doc.flight(0, flight=i)
        sol = solved[i]
        doc.compare(f"{name}Out", f"{what} shot: outcome", sol["code"], g["outcome"])
        doc.compare(f"{name}T", rf"{what} shot: $t_\mathrm{{e}}$", sol["t"], g["end"]["t"], TOL)
        doc.compare(f"{name}X", rf"{what} shot: $x_\mathrm{{e}}$", sol["x"], g["end"]["x"][0], TOL)
        doc.compare(f"{name}Y", rf"{what} shot: $y_\mathrm{{e}}$", sol["y"], g["end"]["x"][1], TOL)
        doc.compare(f"{name}Px", rf"{what} shot: $p_x$", sol["px"], g["end"]["p"][0], TOL)
        doc.compare(f"{name}Py", rf"{what} shot: $p_y$", sol["py"], g["end"]["p"][1], TOL)
        doc.compare(f"{name}Rad", f"{what} shot: radiated energy (rel.)",
                    sol["rad"], g["radiated_energy"], TOL_RAD, rel=True)

    # ---- The symmetric charge and the charge that deflects too much -----------------------
    tags = ("Sym", "Hard")
    for k, tag in enumerate(tags, start=1):
        e = doc.result(k)["placement"][0]
        Qk = mp.mpf(e["value"])
        rck = [mp.mpf(v) / sdiv for v in e["node"][:2]]
        doc.text(f"{tag}Node", f"({e['node'][0]}, {e['node'][1]})")
        doc.value(f"{tag}Q", Qk, 3, sci=True)
        doc.value(f"{tag}K", up.q * Qk, 4)
        mirror = []
        for i, side in enumerate(("Up", "Down")):
            sol = solve(shots[i], shots_f[i], Qk, rck)
            g = doc.flight(k, flight=i)
            doc.compare(f"{tag}{side}Out", f"{tag.lower()} {side.lower()}: outcome",
                        sol["code"], g["outcome"])
            doc.value(f"{tag}{side}T", sol["t"], 8)
            doc.value(f"{tag}{side}X", sol["x"], 8)
            doc.value(f"{tag}{side}Y", sol["y"], 8)
            doc.value(f"{tag}{side}Ang", angle_deg(sol["px"], sol["py"]), 6)
            doc.value(f"{tag}{side}Err", sol["terr"], 2, sci=True)
            if sol["code"] == "Arrived":
                doc.compare(f"{tag}{side}Y", rf"{tag.lower()} {side.lower()}: $y_\mathrm{{e}}$",
                            sol["y"], g["end"]["x"][1], TOL)
                doc.compare(f"{tag}{side}T", rf"{tag.lower()} {side.lower()}: $t_\mathrm{{e}}$",
                            sol["t"], g["end"]["t"], TOL)
            else:
                doc.compare(f"{tag}{side}X", rf"{tag.lower()} {side.lower()}: $x$",
                            sol["x"], g["end"]["x"][0], TOL)
                doc.compare(f"{tag}{side}Y", rf"{tag.lower()} {side.lower()}: $y$",
                            sol["y"], g["end"]["x"][1], TOL)
            mirror.append(sol)
        doc.check(abs(float(mirror[0]["t"] - mirror[1]["t"])) < 1e-12,
                  f"the {tag.lower()} flights take the same time")
        doc.check(abs(float(mirror[0]["y"] + mirror[1]["y"]) - 20) < 1e-9,
                  f"the {tag.lower()} ordinates are mirrors through y = 10")
        doc.check(abs(angle_deg(mirror[0]["px"], mirror[0]["py"])
                      + angle_deg(mirror[1]["px"], mirror[1]["py"])) < 1e-9,
                  f"the {tag.lower()} flights are mirror images")


