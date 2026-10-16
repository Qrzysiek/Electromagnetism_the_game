"""Companion of docs/levels/two_stages/doc.tex (Level 5, Two stages).

The reference charge sits at the origin, and the launch is perpendicular to the radius, so
the flight is Level 1's orbit (first_bend eq:reference) with Level 1's closed-form time
(first_bend eq:time). The new rule is the gate: the trajectory must enter that box before
the detector counts (doc eq:gate). A positive charge at the origin keeps x increasing and y
from falling (first_bend sec:first-event), which locates the gate and the detector by one
root each. The empty placement is a straight line.
"""
import math

import mpmath as mp


# Exact orbit against a verified flight (Foundations, section 6). Set before comparing.
TOL = 1e-9
TOL_RAD = 1e-4


class Shot:
    def __init__(self, lvl, shot):
        s = lvl["grid"]["subdivision"]
        self.q = mp.mpf(shot["particle"]["charge"])
        self.m = mp.mpf(shot["particle"]["mass"])
        self.a = mp.mpf(shot["particle"]["radius"])
        self.c = mp.mpf(lvl["physics"]["c"])
        self.T0 = mp.mpf(shot["launch"]["kinetic_energy"])
        self.R = mp.mpf(lvl["physics"]["charge_radius"])
        self.tmax = mp.mpf(lvl["physics"]["t_max"])
        self.x0 = [mp.mpf(v) / s for v in shot["launch"]["node"][:2]]
        d = [mp.mpf(v) for v in shot["launch"]["direction"][:2]]
        n = (d[0] ** 2 + d[1] ** 2) ** 0.5
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [mp.mpf(v) / s for v in shot["detector"]["min"][:2]]
        self.dmax = [mp.mpf(v) / s for v in shot["detector"]["max"][:2]]
        self.wmin = [mp.mpf(-1), mp.mpf(-1)]
        self.wmax = [mp.mpf(lvl["grid"]["nx"]) / s + 1, mp.mpf(lvl["grid"]["ny"]) / s + 1]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = (self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0) ** 0.5
        g = lvl["gates"][0]
        self.gmin = [mp.mpf(v) / s for v in g["min"][:2]]
        self.gmax = [mp.mpf(v) / s for v in g["max"][:2]]


class Orbit:
    """Level 1's orbit about one charge (first_bend eq:general)."""

    def __init__(self, shot, Q, rc):
        self.s, self.rc = shot, rc
        self.kappa = shot.q * Q
        rho = [shot.x0[0] - rc[0], shot.x0[1] - rc[1]]
        self.r0 = (rho[0] ** 2 + rho[1] ** 2) ** 0.5
        self.phi0 = mp.atan2(rho[1], rho[0])
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
        self.G = self.G2 ** mp.mpf("0.5")
        self.A = self.u0 - self.B

    def u(self, th):
        g = self.G
        return self.B + (self.u0 - self.B) * mp.cos(g * th) + self.du0 * mp.sin(g * th) / g

    def du(self, th):
        g = self.G
        return -(self.u0 - self.B) * g * mp.sin(g * th) + self.du0 * mp.cos(g * th)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.rc[0] + r * mp.cos(phi), self.rc[1] + r * mp.sin(phi))

    def dtdth(self, th):
        u = self.u(th)
        return (self.W - self.kappa * u) / (self.s.c ** 2 * self.L * u * u)

    def dx(self, th):
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        return -mp.sin(phi) * self.sigma / u - mp.cos(phi) * du / (u * u)

    def dy(self, th):
        phi = self.phi0 + self.sigma * th
        u, du = self.u(th), self.du(th)
        return mp.cos(phi) * self.sigma / u - mp.sin(phi) * du / (u * u)

    def momentum(self, th):
        phi = self.phi0 + self.sigma * th
        pr = -self.L * self.du(th)
        pp = self.sigma * self.L * self.u(th)
        return (pr * mp.cos(phi) - pp * mp.sin(phi),
                pr * mp.sin(phi) + pp * mp.cos(phi))


def newton(f, fp, th0, n=8):
    """Newton's method (Foundations eq:newton). Returns (root, rows)."""
    th, rows = th0, []
    for k in range(n):
        val = f(th)
        rows.append((k, th, val))
        if abs(val) < mp.mpf("1e-28"):
            break
        th = th - val / fp(th)
    return th, rows


def root_axis(orb, axis, value, th0):
    """The root of position[axis] = value, started at th0."""
    deriv = orb.dx if axis == 0 else orb.dy
    f = lambda th: orb.pos(th)[axis] - value
    return newton(f, deriv, th0)


def closed_time(orb, th):
    """Level 1's elementary time, for an apsis (first_bend eq:time, eq:i1, eq:i2)."""
    s = orb.G * th
    k2 = orb.A ** 2 - orb.B ** 2
    I1 = 2 / mp.sqrt(k2) * mp.atanh(
        mp.sqrt((orb.A - orb.B) / (orb.A + orb.B)) * mp.tan(s / 2))
    I2 = orb.A * mp.sin(s) / (k2 * (orb.A * mp.cos(s) + orb.B)) - orb.B / k2 * I1
    t = (orb.W * I2 - orb.kappa * I1) / (orb.s.c ** 2 * orb.L * orb.G)
    return t, I1, I2


def radiated(orb, th):
    """Simpson on the Liénard power, 8000 panels against 4000 (Foundations eq:simpson)."""
    shot, rc = orb.s, orb.rc

    def power_dt(psi):
        xx, yy = orb.pos(psi)
        rx, ry = xx - rc[0], yy - rc[1]
        r3 = (rx * rx + ry * ry) ** mp.mpf("1.5")
        F = (orb.kappa * rx / r3, orb.kappa * ry / r3)
        px, py = orb.momentum(psi)
        gam = (orb.W - orb.kappa * orb.u(psi)) / (shot.m * shot.c ** 2)
        v = (px / (gam * shot.m), py / (gam * shot.m))
        vF = v[0] * F[0] + v[1] * F[1]
        ax = (F[0] - v[0] * vF / shot.c ** 2) / (gam * shot.m)
        ay = (F[1] - v[1] * vF / shot.c ** 2) / (gam * shot.m)
        vxa = v[0] * ay - v[1] * ax
        P = mp.mpf(2) / 3 * shot.q ** 2 / shot.c ** 3 * gam ** 6 * (
            ax * ax + ay * ay - vxa ** 2 / shot.c ** 2)
        return P * orb.dtdth(psi)

    def simp(n):
        h = th / n
        acc = power_dt(0) + power_dt(th)
        for k in range(1, n):
            acc += (4 if k % 2 else 2) * power_dt(h * k)
        return h / 3 * acc

    e1, e2 = simp(4000), simp(8000)
    return e2, abs(e2 - e1) / abs(e2) / 15


def gate_entry(orb, shot):
    """First entry of a path whose x and y both increase, into the one gate.

    The path starts below and to the left of the box. It meets the box at the later of
    the left face and the bottom face, and that point is an entry when it is still
    inside the other two sides (doc eq:gate).
    """
    th_left, left_rows = root_axis(orb, 0, shot.gmin[0], mp.mpf("0.4"))
    y_left = orb.pos(th_left)[1]
    th_bottom, bot_rows = root_axis(orb, 1, shot.gmin[1], mp.mpf("0.4"))
    x_bottom = orb.pos(th_bottom)[0]
    # Judge the face by the free coordinate. The root's own coordinate is the face
    # value: a Newton residual must not push a point that lies on the face outside it.
    if shot.gmin[1] <= y_left <= shot.gmax[1]:
        th, face, x, y, nrow = th_left, "left", shot.gmin[0], y_left, left_rows
        entered = True
    elif shot.gmin[0] <= x_bottom <= shot.gmax[0]:
        th, face, x, y, nrow = th_bottom, "bottom", x_bottom, shot.gmin[1], bot_rows
        entered = True
    else:
        th, face, x, y, nrow = th_bottom, "missed", x_bottom, shot.gmin[1], bot_rows
        entered = False
    return {
        "entered": entered, "theta": th, "face": face, "x": x, "y": y,
        "th_left": th_left, "y_left": y_left, "th_bottom": th_bottom, "x_bottom": x_bottom,
        "newton": nrow,
    }


def ending(orb, shot):
    """The first ending event, using x increasing and y non-decreasing.

    The detector's left face counts when the ordinate there lies on it. Its bottom face
    counts when the abscissa there lies on it. The world's right edge is the other
    candidate. The earlier one wins. A gate is not an ending.
    """
    th_left, det_rows = root_axis(orb, 0, shot.dmin[0], mp.mpf(1))
    y_left = orb.pos(th_left)[1]
    th_world, _ = root_axis(orb, 0, shot.wmax[0], mp.mpf(1))
    # x and y both increase, so the detector is entered through its left face or,
    # if the ray is still below that face, through its bottom, and only while that
    # crossing lies on the face.
    hit = None
    if shot.dmin[1] <= y_left <= shot.dmax[1]:
        hit = th_left
    elif y_left < shot.dmin[1]:
        th_bot, _ = root_axis(orb, 1, shot.dmin[1], th_left)
        x_bot = orb.pos(th_bot)[0]
        if shot.dmin[0] <= x_bot <= shot.dmax[0] and th_bot < th_world:
            hit = th_bot
    if hit is not None:
        return "Arrived", hit, det_rows
    return "LeftBounds", th_world, det_rows


def kinetic(shot, px, py):
    g = mp.sqrt(1 + (px ** 2 + py ** 2) / (shot.m ** 2 * shot.c ** 2))
    return (g - 1) * shot.m * shot.c ** 2


def flight_of(shot, Q, rc):
    orb = Orbit(shot, Q, rc)
    code, th, rows = ending(orb, shot)
    gate = gate_entry(orb, shot)
    # A detector entry with the gate still ahead is a skipped gate (doc eq:gate).
    if code == "Arrived" and (not gate["entered"] or gate["theta"] > th):
        code = "SkippedGate(0)"
    x, y = orb.pos(th)
    px, py = orb.momentum(th)
    t, I1, I2 = closed_time(orb, th)
    t_gate = closed_time(orb, gate["theta"])[0] if gate["entered"] else None
    rad, raderr = radiated(orb, th)
    path = [orb.pos(th * k / 80) for k in range(81)]
    return {
        "orb": orb, "code": code, "theta": th, "newton": rows, "gate": gate,
        "x": x, "y": y, "px": px, "py": py, "t": t, "I1": I1, "I2": I2,
        "t_gate": t_gate, "rad": rad, "raderr": raderr,
        "path": [(float(p[0]), float(p[1])) for p in path],
    }


def compute(doc):
    mp.mp.dps = 40
    lvl = doc.level
    shot = Shot(lvl, lvl["shots"][0])
    sdiv = lvl["grid"]["subdivision"]
    ref = lvl["reference_solution"][0]
    Q = mp.mpf(ref["value"])
    rc = [mp.mpf(v) / sdiv for v in ref["node"][:2]]

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
    doc.value("xD", shot.dmin[0], 2)
    doc.value("yDlo", shot.dmin[1], 2)
    doc.value("yDhi", shot.dmax[1], 2)
    doc.value("xGlo", shot.gmin[0], 2)
    doc.value("xGhi", shot.gmax[0], 2)
    doc.value("yGlo", shot.gmin[1], 2)
    doc.value("yGhi", shot.gmax[1], 2)
    doc.value("xWorld", shot.wmax[0], 2)

    sol = flight_of(shot, Q, rc)
    orb = sol["orb"]
    doc.value("rZero", orb.r0, 6)
    doc.value("W", orb.W, 12)
    doc.value("L", orb.L, 12)
    doc.value("Gamma", orb.G, 12)
    doc.value("B", orb.B, 12)
    doc.value("A", orb.A, 12)
    A_energy = (orb.B ** 2 + (orb.W ** 2 - shot.m ** 2 * shot.c ** 4)
                / (shot.c ** 2 * orb.L ** 2 * orb.G2)) ** mp.mpf("0.5")
    doc.value("Adiff", abs(orb.A - A_energy), 2, sci=True)
    doc.check(abs(orb.pr0) < mp.mpf("1e-30"), "the launch is an apsis")
    doc.check(abs(orb.A - A_energy) < mp.mpf("1e-30"), "A from the energy agrees")
    doc.check(sol["code"] == "Arrived", "the reference arrives")
    doc.check(sol["gate"]["entered"] and sol["gate"]["theta"] < sol["theta"],
              "the reference enters the gate before the detector")

    gate = sol["gate"]
    doc.text("gateFace", gate["face"])
    doc.value("thetaG", gate["theta"], 10)
    doc.value("xG", gate["x"], 10)
    doc.value("yG", gate["y"], 10)
    doc.value("tG", sol["t_gate"], 10)
    doc.value("yAtLeft", gate["y_left"], 8)
    doc.value("xAtBottom", gate["x_bottom"], 8)
    doc.value("thetaE", sol["theta"], 12)
    doc.value("tE", sol["t"], 12)
    doc.value("yE", sol["y"], 12)
    doc.value("pxE", sol["px"], 12)
    doc.value("pyE", sol["py"], 12)
    doc.value("angE", float(mp.atan2(sol["py"], sol["px"]) * 180 / mp.pi), 6)
    doc.value("TE", kinetic(shot, sol["px"], sol["py"]), 8)
    doc.value("IOne", sol["I1"], 10)
    doc.value("ITwo", sol["I2"], 10)
    doc.value("Erad", sol["rad"], 6, sci=True)
    doc.value("EradFrac", sol["rad"] / shot.T0, 3, sci=True)
    doc.value("radErr", sol["raderr"], 2, sci=True)
    doc.data("path", sol["path"])

    rows = [[str(k), th, f] for k, th, f in sol["newton"]]
    doc.table("newton", [r"$k$", r"$\theta_k$", r"$x(\theta_k) - x_D$"], rows,
              align="rll", digits=[None, 12, 3], sci=[False, False, True])
    doc.value("newtonSteps", len(rows) - 1)
    if gate["newton"]:
        grows = [[str(k), th, f] for k, th, f in gate["newton"]]
        residual = r"$x(\theta_k) - x_G$" if gate["face"] == "left" else r"$y(\theta_k) - y_G$"
        doc.table("gateNewton", [r"$k$", r"$\theta_k$", residual], grows,
                  align="rll", digits=[None, 10, 3], sci=[False, False, True])

    g = doc.flight(0)
    doc.compare("Out", "outcome", sol["code"], g["outcome"])
    doc.compare("T", r"$t_\mathrm{e}$", sol["t"], g["end"]["t"], TOL)
    doc.compare("X", r"$x_\mathrm{e}$", sol["x"], g["end"]["x"][0], TOL)
    doc.compare("Y", r"$y_\mathrm{e}$", sol["y"], g["end"]["x"][1], TOL)
    doc.compare("Px", r"$p_x$", sol["px"], g["end"]["p"][0], TOL)
    doc.compare("Py", r"$p_y$", sol["py"], g["end"]["p"][1], TOL)
    doc.compare("Rad", "radiated energy (rel.)", sol["rad"], g["radiated_energy"], TOL_RAD, rel=True)

    # Straight flight: no charge. Radiation is identically zero.
    v0 = shot.p0 / (shot.gamma0 * shot.m)
    t_free = shot.wmax[0] / v0
    doc.value("vZero", v0, 10)
    doc.value("freeT", t_free, 12)
    doc.value("freeX", shot.wmax[0], 2)
    doc.value("freeY", shot.x0[1], 2)
    doc.check(shot.x0[1] < shot.gmin[1], "the straight ray passes under the gate")
    doc.check(shot.x0[1] < shot.dmin[1], "the straight ray passes under the detector")
    g1 = doc.flight(1)
    doc.compare("FreeOut", "no charge: outcome", "LeftBounds", g1["outcome"])
    doc.compare("FreeT", r"no charge: $t$", t_free, g1["end"]["t"], TOL)
    doc.compare("FreeY", r"no charge: $y$", shot.x0[1], g1["end"]["x"][1], TOL)
    doc.data("pathfree", [(0.0, float(shot.x0[1])), (float(shot.wmax[0]), float(shot.x0[1]))])

    # The weaker charge at the same node.
    e = doc.result(2)["placement"][0]
    Qh = mp.mpf(e["value"])
    rch = [mp.mpf(v) / sdiv for v in e["node"][:2]]
    half = flight_of(shot, Qh, rch)
    doc.value("halfQ", Qh, 3, sci=True)
    doc.value("halfK", shot.q * Qh, 4)
    doc.value("halfT", half["t"], 10)
    doc.value("halfX", half["x"], 10)
    doc.value("halfY", half["y"], 10)
    doc.value("halfAng", float(mp.atan2(half["py"], half["px"]) * 180 / mp.pi), 6)
    doc.value("halfErad", half["rad"], 6, sci=True)
    doc.value("halfRadErr", half["raderr"], 2, sci=True)
    doc.text("halfGate", "entered" if half["gate"]["entered"] else "missed")
    doc.text("halfFace", half["gate"]["face"])
    doc.value("halfXg", half["gate"]["x"], 8)
    doc.value("halfYg", half["gate"]["y"], 8)
    doc.value("halfYLeft", half["gate"]["y_left"], 6)
    # Ordinate as the ray crosses the detector's abscissa, whether or not it is on the face.
    th27, _ = root_axis(half["orb"], 0, shot.dmin[0], mp.mpf(1))
    doc.value("halfY27", half["orb"].pos(th27)[1], 8)
    th14, _ = root_axis(half["orb"], 1, shot.dmin[1], th27)
    doc.value("halfX14", half["orb"].pos(th14)[0], 8)
    doc.check(half["orb"].pos(th14)[0] > shot.dmax[0],
              "the weaker ray reaches the detector's floor only past its right side")
    doc.check(half["code"] == "LeftBounds", "the weaker charge leaves the world")
    doc.check(half["gate"]["entered"], "the weaker charge does pass the gate")
    doc.check(half["orb"].pos(th27)[1] < shot.dmin[1],
              "at the detector's abscissa the weaker ray is still below it")
    doc.data("pathhalf", half["path"])
    g2 = doc.flight(2)
    doc.compare("HalfOut", "weaker charge: outcome", half["code"], g2["outcome"])
    doc.compare("HalfT", r"weaker charge: $t$", half["t"], g2["end"]["t"], TOL)
    doc.compare("HalfX", r"weaker charge: $x$", half["x"], g2["end"]["x"][0], TOL)
    doc.compare("HalfY", r"weaker charge: $y$", half["y"], g2["end"]["x"][1], TOL)
    doc.compare("HalfRad", "weaker charge: radiated energy (rel.)",
                half["rad"], g2["radiated_energy"], TOL_RAD, rel=True)
