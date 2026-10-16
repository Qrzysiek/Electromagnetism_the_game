"""Companion of docs/levels/injection/doc.tex (Level 6, Injection).

One attractive charge. The flight is Level 1's orbit (first_bend eq:general) about that
charge. The launch is not an apsis, so the time and the radiated energy are Simpson's
rule (Foundations eq:simpson) on d t/d theta (first_bend eq:dtdth). The detector face is
the root of the orbit, polished by Newton's method (Foundations eq:newton). The direction
cone is the angle between the arrival velocity and the detector's axis.

The step-doubling error of the time quadrature is required to be below 1e-12 before any
comparison, so it does not set the tolerance. A verified flight is accurate to about
1e-10 (Foundations, section 6). Ten times that.
"""
import math

import mpmath as mp
import numpy as np


TOL = 1e-9
# Relative tolerance on the radiated energy. The quadrature's own relative error is
# checked below, and is far smaller than this; 1e-4 is the room used for the same
# integral in the earlier one-charge documents.
TOL_RAD = 1e-4
# Simpson panels on [0, theta_e]. The reported value uses 2*N_SIMP; the error is
# |S_n - S_{2n}|/15 (Foundations eq:simpson).
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
        self.dir_norm = (d[0] ** 2 + d[1] ** 2) ** 0.5
        self.d = [d[0] / self.dir_norm, d[1] / self.dir_norm]
        self.dmin = [mpf(v) / s for v in shot["detector"]["min"][:2]]
        self.dmax = [mpf(v) / s for v in shot["detector"]["max"][:2]]
        self.wmin = [mpf(-1), mpf(-1)]
        self.wmax = [mpf(lvl["grid"]["nx"]) / s + 1, mpf(lvl["grid"]["ny"]) / s + 1]
        acc = shot["detector"]["acceptance"]["direction"]
        self.axis_deg = mpf(acc[0])
        self.half_deg = mpf(acc[1])
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = (self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0) ** 0.5
        self.v0 = self.p0 / (self.gamma0 * self.m)


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
    vals = [float(fn(x)) for x in xs]
    acc = vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2])
    return h / 3 * acc


def scan(orb, shot, dth=1e-3):
    """First geometric event along the orbit. Returns (code, theta, face).

    code is Arrived, LeftBounds or Collided, before the direction cone. `face` is
    (axis, value) of the boundary the chord reaches first.
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


def deviation(px, py, axis):
    """Angle in radians from `axis` (a unit vector) to the momentum, in [0, pi].

    v is parallel to p, so |v × ê| / (v · ê) is |p × ê| / (p · ê).
    """
    cross = abs(px * axis[1] - py * axis[0])
    dot = px * axis[0] + py * axis[1]
    return mp.atan2(cross, dot)


def kinetic(shot, px, py):
    g = mp.sqrt(1 + (px ** 2 + py ** 2) / (shot.m ** 2 * shot.c ** 2))
    return (g - 1) * shot.m * shot.c ** 2


def radiated(orb, shot, rc, th):
    """Simpson integral of the Larmor power along the orbit, and its relative error."""

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

    e1 = simpson(power_dt, 0, th, 4000)
    e2 = simpson(power_dt, 0, th, 8000)
    return e2, abs(e2 - e1) / abs(e2) / 15


def closest(orb, th_end):
    """The minimum of r on [0, th_end]: a critical point of u, or an endpoint."""
    g = mp.sqrt(orb.G2)
    th = mp.atan(orb.du0 / (g * (orb.u0 - orb.B))) / g
    cands = [mp.mpf(0), th_end]
    for base in (th, th + mp.pi / g, th - mp.pi / g):
        if mp.mpf(0) <= base <= th_end:
            cands.append(base)
    best = max(cands, key=lambda t: orb.u(t))
    return best, 1 / orb.u(best)


def with_cone(geometric, dev, half):
    """The detector counts an entry only when the deviation lies inside the cone."""
    if geometric != "Arrived":
        return geometric
    return "Arrived" if dev <= half else "Rejected"


def solve(shot, shot_f, Q, rc, axis, half):
    """One charge: the face, the state, the time, the radiation, and the cone."""
    orb = Orbit(shot, Q, rc, mp)
    geometric, th_bracket, face = scan(
        Orbit(shot_f, float(Q), [float(c) for c in rc], math), shot_f)
    th, rows = newton_face(orb, face[0], face[1], th_bracket + mp.mpf("5e-4"))
    residual = abs(orb.pos(th)[face[0]] - face[1])
    x, y = orb.pos(th)
    if face[0] == 0:
        x = face[1]
    else:
        y = face[1]
    px, py = orb.momentum(th)
    dev = deviation(px, py, axis)
    n = N_SIMP
    t1 = simpson(orb.dtdth, 0, th, n)
    t2 = simpson(orb.dtdth, 0, th, 2 * n)
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
        "orb": orb, "code": with_cone(geometric, dev, half), "geometric": geometric,
        "theta": th, "newton": rows, "face": face, "residual": residual,
        "guess": th_bracket + mp.mpf("5e-4"),
        "x": x, "y": y, "px": px, "py": py, "dev": dev,
        "t": t2, "terr": abs(t2 - t1) / 15,
        "rad": rad, "raderr": raderr,
        "path": samples, "rmin": r_min, "thmin": th_min,
        "T": kinetic(shot, px, py),
    }


def straight(shot, axis, half):
    """No charge: uniform motion. The first box entered, then the cone."""
    vx, vy = shot.v0 * shot.d[0], shot.v0 * shot.d[1]
    det = _box_hit(shot.x0, (vx, vy), shot.dmin, shot.dmax, entering=True)
    world = _box_hit(shot.x0, (vx, vy), shot.wmin, shot.wmax, entering=False)
    if det is not None and (world is None or det[0] < world[0]):
        t, face, x, y = det
        geometric = "Arrived"
    else:
        t, face, x, y = world
        geometric = "LeftBounds"
    px, py = shot.p0 * shot.d[0], shot.p0 * shot.d[1]
    dev = deviation(px, py, axis)
    path = []
    for k in range(81):
        tk = t * k / 80
        path.append((float(shot.x0[0] + vx * tk), float(shot.x0[1] + vy * tk)))
    return {
        "code": with_cone(geometric, dev, half), "geometric": geometric,
        "t": t, "x": x, "y": y, "px": px, "py": py, "dev": dev,
        "face": face, "path": path, "T": kinetic(shot, px, py),
        "world_t": world[0],
    }


def _box_hit(x0, v, lo, hi, entering):
    """First t>0 at which the ray crosses into `entering` or out of (not entering) the box.

    Returns (t, (axis, value), x, y) with the face coordinate set to the face value.
    """
    best = None
    for axis in (0, 1):
        for side, face in ((True, lo[axis]), (False, hi[axis])):
            inward = (v[axis] > 0) if side else (v[axis] < 0)
            if v[axis] == 0 or inward != entering:
                continue
            t = (face - x0[axis]) / v[axis]
            if t <= 0:
                continue
            pos = [x0[0] + v[0] * t, x0[1] + v[1] * t]
            other = 1 - axis
            if lo[other] <= pos[other] <= hi[other] and (best is None or t < best[0]):
                pos[axis] = face
                best = (t, (axis, face), pos[0], pos[1])
    return best


def axis_vector(shot):
    """Unit axis of the cone. The file stores the axis angle in degrees."""
    a = shot.axis_deg * mp.pi / 180
    return (mp.cos(a), mp.sin(a)), shot.half_deg * mp.pi / 180


def compare_flight(doc, tag, what, sol, k, axis_f, radiated_rel):
    g = doc.flight(k)
    doc.compare(f"{tag}Out", f"{what}: outcome", sol["code"], g["outcome"])
    doc.compare(f"{tag}T", rf"{what}: $t_\mathrm{{e}}$", sol["t"], g["end"]["t"], TOL)
    doc.compare(f"{tag}X", rf"{what}: $x_\mathrm{{e}}$", sol["x"], g["end"]["x"][0], TOL)
    doc.compare(f"{tag}Y", rf"{what}: $y_\mathrm{{e}}$", sol["y"], g["end"]["x"][1], TOL)
    doc.compare(f"{tag}Px", rf"{what}: $p_x$", sol["px"], g["end"]["p"][0], TOL)
    doc.compare(f"{tag}Py", rf"{what}: $p_y$", sol["py"], g["end"]["p"][1], TOL)
    # The game stores momentum, not the angle. The deviation is this document's test,
    # applied to both momenta. Only the game momentum is read here.
    game_dev = math.atan2(
        abs(g["end"]["p"][0] * axis_f[1] - g["end"]["p"][1] * axis_f[0]),
        g["end"]["p"][0] * axis_f[0] + g["end"]["p"][1] * axis_f[1])
    doc.compare(f"{tag}Dev", rf"{what}: deviation $\beta$", sol["dev"], game_dev, TOL)
    if radiated_rel:
        doc.compare(f"{tag}Rad", f"{what}: radiated energy (rel.)",
                    sol["rad"], g["radiated_energy"], TOL_RAD, rel=True)
    else:
        # No force, so the power is identically zero. An absolute comparison: the game
        # value is zero and a relative difference is undefined.
        doc.compare(f"{tag}Rad", f"{what}: radiated energy", 0, g["radiated_energy"], 0)


def compute(doc):
    mp.mp.dps = 40
    lvl = doc.level
    ref = lvl["reference_solution"][0]
    sdiv = lvl["grid"]["subdivision"]
    Q = mp.mpf(ref["value"])
    rc = [mp.mpf(v) / sdiv for v in ref["node"][:2]]
    shot = Shot(lvl, lvl["shots"][0], mp.mpf)
    shot_f = Shot(lvl, lvl["shots"][0], float)
    axis, half = axis_vector(shot)
    axis_f = (float(axis[0]), float(axis[1]))

    doc.value("q", shot.q, 3, sci=True)
    doc.value("m", shot.m, 3)
    doc.value("c", shot.c, 3)
    doc.value("Tzero", shot.T0, 3)
    doc.value("R", shot.R, 3)
    doc.value("tmax", shot.tmax, 4)
    doc.value("Q", Q, 3, sci=True)
    doc.value("kappa", shot.q * Q, 4)
    doc.value("gammaZero", shot.gamma0, 10)
    doc.value("pZero", shot.p0, 12)
    doc.value("vZero", shot.v0, 10)
    doc.value("xLaunch", shot.x0[0], 2)
    doc.value("yLaunch", shot.x0[1], 2)
    doc.value("xC", rc[0], 2)
    doc.value("yC", rc[1], 2)
    doc.value("xD", shot.dmin[0], 2)
    doc.value("xDhi", shot.dmax[0], 2)
    doc.value("yDlo", shot.dmin[1], 2)
    doc.value("yDhi", shot.dmax[1], 2)
    doc.value("xWorld", shot.wmax[0], 2)
    doc.value("yWorld", shot.wmax[1], 2)
    doc.value("axisDeg", shot.axis_deg, 2)
    doc.value("halfDeg", shot.half_deg, 2)
    doc.value("halfRad", half, 10)
    doc.value("nSimp", 2 * N_SIMP)
    doc.value("nRad", 8000)
    launch = mp.atan2(shot.d[1], shot.d[0])
    doc.value("launchAng", mp.degrees(launch), 12)
    doc.value("dirResidual", abs(launch - mp.mpf(20) * mp.pi / 180), 3, sci=True)
    doc.value("dirNormGap", abs(shot.dir_norm - 1), 3, sci=True)

    length = mp.mpf("2.4")
    doc.value("arrowX", shot.x0[0] + length * shot.d[0], 6)
    doc.value("arrowY", shot.x0[1] + length * shot.d[1], 6)
    # The cone drawn on the figure is the accepted wedge, from the middle of the entrance.
    cone_o = (shot.dmin[0], (shot.dmin[1] + shot.dmax[1]) / 2)
    reach = mp.mpf(4)
    doc.value("coneOx", cone_o[0], 4)
    doc.value("coneOy", cone_o[1], 4)
    doc.value("coneXp", cone_o[0] + reach * mp.cos(half), 6)
    doc.value("coneYp", cone_o[1] + reach * mp.sin(half), 6)
    doc.value("coneXm", cone_o[0] + reach * mp.cos(-half), 6)
    doc.value("coneYm", cone_o[1] + reach * mp.sin(-half), 6)

    sol = solve(shot, shot_f, Q, rc, axis, half)
    orb = sol["orb"]
    doc.value("rZero", orb.r0, 10)
    doc.value("phiZero", orb.phi0, 10)
    doc.value("Lz", orb.sigma * orb.L, 12)
    doc.value("L", orb.L, 12)
    doc.value("sigma", orb.sigma)
    doc.value("prZero", orb.pr0, 10)
    doc.value("W", orb.W, 12)
    doc.value("Gamma", mp.sqrt(orb.G2), 12)
    doc.value("B", orb.B, 10)
    doc.value("uZero", orb.u0, 10)
    doc.value("duZero", orb.du0, 10)
    doc.value("gGap", 1 - mp.sqrt(orb.G2), 3, sci=True)
    doc.value("thetaE", sol["theta"], 12)
    doc.value("newtonRes", sol["residual"], 3, sci=True)
    doc.value("newtonSteps", len(sol["newton"]) - 1)
    doc.value("thGuess", sol["guess"], 6)
    doc.value("tE", sol["t"], 12)
    doc.value("xE", sol["x"], 12)
    doc.value("yE", sol["y"], 12)
    doc.value("pxE", sol["px"], 12)
    doc.value("pyE", sol["py"], 12)
    doc.value("angE", mp.degrees(sol["dev"]), 8)
    doc.value("devE", sol["dev"], 10)
    doc.value("marginE", mp.degrees(half - sol["dev"]), 6)
    doc.value("TE", sol["T"], 8)
    doc.value("terr", sol["terr"], 2, sci=True)
    doc.value("rmin", sol["rmin"], 6)
    doc.value("thMin", sol["thmin"], 6)
    doc.value("rOverR", sol["rmin"] / shot.R, 4)
    doc.value("Erad", sol["rad"], 6, sci=True)
    doc.value("EradFrac", sol["rad"] / shot.T0, 3, sci=True)
    doc.value("radErr", sol["raderr"], 2, sci=True)
    doc.value("potFrac", abs(orb.kappa) / (orb.r0 * shot.T0), 4)
    doc.data("path", sol["path"])

    doc.check(orb.G2 > 0, "Gamma^2 is positive")
    doc.check(abs(orb.du0) > mp.mpf("1e-3"), "the launch is not an apsis")
    doc.check(abs(orb.pos(0)[0] - shot.x0[0]) < mp.mpf("1e-25"), "the orbit starts at the launch")
    doc.check(sol["geometric"] == "Arrived", "the reference enters the detector")
    # The left face is crossed while the ray is still below the box, so it is not the entry.
    th_left, _ = newton_face(orb, 0, shot.dmin[0], sol["theta"] - mp.mpf("0.05"))
    y_left = orb.pos(th_left)[1]
    doc.value("thetaLeft", th_left, 8)
    doc.value("yAtLeft", y_left, 8)
    doc.value("leftRes", abs(orb.pos(th_left)[0] - shot.dmin[0]), 3, sci=True)
    doc.check(sol["face"][0] == 1 and sol["face"][1] == shot.dmin[1],
              "the reference meets the floor of the detector")
    doc.check(abs(orb.pos(th_left)[0] - shot.dmin[0]) < mp.mpf("1e-18"),
              "the left-face root is on x = 28")
    doc.check(th_left < sol["theta"], "the left face is crossed before the floor")
    doc.check(y_left < shot.dmin[1], "at the left face the reference is still below the detector")
    doc.check(shot.dmin[0] <= sol["x"] <= shot.dmax[0], "the floor crossing lies on the detector")
    doc.check(sol["code"] == "Arrived", "the reference arrival lies inside the cone")
    doc.check(sol["dev"] < half, "the reference deviation is below the half-angle")
    doc.check(sol["residual"] < mp.mpf("1e-18"), "Newton lands on the face")
    doc.check(sol["terr"] < 1e-12, "the time quadrature is finer than the tolerance")
    doc.check(sol["rmin"] > 10 * shot.R, "the closest approach is far outside the sphere")
    doc.check(0 < sol["thmin"] < sol["theta"], "the closest approach lies on the flight")
    doc.check(abs(orb.du(sol["thmin"])) < mp.mpf("1e-18"), "the closest approach is a critical point")
    doc.check(sol["t"] < shot.tmax, "the reference arrives before the time limit")
    doc.check(abs(shot.axis_deg) == 0, "the cone axis is +x")
    doc.check("kinetic" not in lvl["shots"][0]["detector"]["acceptance"], "there is no kinetic window")

    rows = []
    for k, th, f in sol["newton"]:
        rows.append([str(k), th, f])
    doc.table("newton", [r"$k$", r"$\theta_k$", r"$y(\theta_k) - y_D$"], rows,
              align="rll", digits=[None, 12, 3], sci=[False, False, True])
    samples = []
    for frac in (0, 0.25, 0.5, 0.75, 1):
        th = sol["theta"] * frac
        x, y = shot.x0 if frac == 0 else orb.pos(th)
        samples.append([th, x, y])
    doc.table("samples", [r"$\theta$", r"$x$", r"$y$"], samples,
              align="rrr", digits=[6, 5, 5])

    # ---- No charge: a straight ray, then the same cone ------------------------------------
    free = straight(shot, axis, half)
    doc.value("freeT", free["t"], 12)
    doc.value("freeX", free["x"], 12)
    doc.value("freeY", free["y"], 12)
    doc.value("freePx", free["px"], 12)
    doc.value("freePy", free["py"], 12)
    doc.value("freeAng", mp.degrees(free["dev"]), 8)
    doc.value("freeDev", free["dev"], 10)
    doc.value("freeMargin", mp.degrees(half - free["dev"]), 6)
    doc.value("freeKin", free["T"], 8)
    doc.data("pathfree", free["path"])
    doc.check(free["geometric"] == "Arrived", "the straight ray enters the detector")
    doc.check(free["face"][0] == 0 and free["face"][1] == shot.dmin[0],
              "the straight ray meets the left face")
    doc.check(shot.dmin[1] <= free["y"] <= shot.dmax[1], "the left-face crossing lies on the detector")
    doc.check(free["code"] == "Rejected", "the straight ray is outside the cone")
    doc.check(free["t"] < free["world_t"], "the detector comes before the edge of the world")
    doc.check(abs(free["dev"] - launch) < mp.mpf("1e-30"),
              "with no charge the arrival direction is the launch direction")

    # ---- Half the reference charge, read from the recorded placement ----------------------
    placed = doc.result(2)["placement"][0]
    Qh = mp.mpf(placed["value"])
    rch = [mp.mpf(v) / sdiv for v in placed["node"][:2]]
    doc.check(rch[0] == rc[0] and rch[1] == rc[1], "the weaker charge sits at the reference node")
    weak = solve(shot, shot_f, Qh, rch, axis, half)
    doc.text("halfNode", f"({float(rch[0]):g}, {float(rch[1]):g})")
    doc.value("halfQ", Qh, 3, sci=True)
    doc.value("halfK", shot.q * Qh, 4)
    doc.value("halfTheta", weak["theta"], 8)
    doc.value("halfT", weak["t"], 12)
    doc.value("halfX", weak["x"], 12)
    doc.value("halfY", weak["y"], 12)
    doc.value("halfPx", weak["px"], 12)
    doc.value("halfPy", weak["py"], 12)
    doc.value("halfAng", mp.degrees(weak["dev"]), 6)
    doc.value("halfDev", weak["dev"], 10)
    doc.value("halfMargin", mp.degrees(half - weak["dev"]), 6)
    doc.value("halfKin", weak["T"], 8)
    doc.value("halfTerr", weak["terr"], 2, sci=True)
    doc.value("halfRmin", weak["rmin"], 6)
    doc.value("halfErad", weak["rad"], 6, sci=True)
    doc.value("halfEradFrac", weak["rad"] / shot.T0, 3, sci=True)
    doc.value("halfRadErr", weak["raderr"], 2, sci=True)
    doc.value("halfRes", weak["residual"], 3, sci=True)
    doc.data("pathhalf", weak["path"])
    doc.check(weak["geometric"] == "Arrived", "the weaker charge enters the detector")
    doc.check(weak["face"][0] == 0 and weak["face"][1] == shot.dmin[0],
              "the weaker charge meets the left face")
    doc.check(shot.dmin[1] <= weak["y"] <= shot.dmax[1],
              "the weaker charge's left-face crossing lies on the detector")
    doc.check(weak["code"] == "Rejected", "the weaker charge is still outside the cone")
    doc.check(weak["terr"] < 1e-12, "the weaker charge's time quadrature is finer than the tolerance")
    doc.check(weak["rmin"] > 10 * shot.R, "the weaker charge stays far from the sphere")
    doc.check(weak["residual"] < mp.mpf("1e-18"), "Newton lands on the weaker charge's face")

    # Comparisons last, after every error estimate above.
    compare_flight(doc, "Ref", "reference", sol, 0, axis_f, True)
    compare_flight(doc, "Free", "no charge", free, 1, axis_f, False)
    compare_flight(doc, "Half", "weaker charge", weak, 2, axis_f, True)
