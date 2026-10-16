"""Companion of docs/levels/hemispherical_analyzer/doc.tex (Level 12, Hemispherical analyzer).

The exact orbit of one charge (first_bend eq:general). Every launch is an apsis, so
u = A cos(Gamma theta) + B and the time is the elementary integral (first_bend eq:time),
with the arctan form of I1 because B > |A|, split at s = pi where Gamma theta exceeds pi.
The first event is the earliest face a chord meets; Newton's method then puts the orbit
on that face. Radiation is Simpson's rule on the power.
"""
import math

import mpmath as mp

DIGITS = 40
# The orbit, the time and the face root are closed forms or Newton's method taken to a
# residual below 1e-28. The game's verified positions are about 1e-10
# (Foundations, section 6). The state tolerance is 1e-9, ten times that, and it is the
# floor for an exact result. The radiated energy is the game's own quadrature
# (PHYSICS.md 8: a few 1e-9 relative on an exact test); the radiation tolerance is the
# larger of 1e-4 and ten times this file's Simpson uncertainty. Both are fixed before
# any comparison.
TOL = 1e-9
TOL_RAD = 1e-4
DTH = 1e-4


class Board:
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
    """first_bend eq:general. theta increases in the direction of motion."""

    def __init__(self, b, Q, rc, lib):
        self.b, self.lib = b, lib
        self.kappa = b.q * Q
        self.rc = rc
        rho = [b.x0[0] - rc[0], b.x0[1] - rc[1]]
        self.r0 = (rho[0] ** 2 + rho[1] ** 2) ** 0.5
        self.phi0 = lib.atan2(rho[1], rho[0])
        p = [b.p0 * b.d[0], b.p0 * b.d[1]]
        Lz = rho[0] * p[1] - rho[1] * p[0]
        self.Lz = Lz
        self.sigma = 1 if Lz > 0 else -1
        self.L = abs(Lz)
        self.pr0 = (p[0] * rho[0] + p[1] * rho[1]) / self.r0
        self.W = b.gamma0 * b.m * b.c ** 2 + self.kappa / self.r0
        self.u0 = 1 / self.r0
        self.G2 = 1 - self.kappa ** 2 / (self.L ** 2 * b.c ** 2)
        self.B = -self.kappa * self.W / (b.c ** 2 * self.L ** 2 * self.G2)
        self.du0 = -self.pr0 / self.L
        self.A = self.u0 - self.B

    def u(self, th):
        g = self.G2 ** 0.5
        return self.B + self.A * self.lib.cos(g * th) + self.du0 * self.lib.sin(g * th) / g

    def du(self, th):
        g = self.G2 ** 0.5
        return -self.A * g * self.lib.sin(g * th) + self.du0 * self.lib.cos(g * th)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.rc[0] + r * self.lib.cos(phi), self.rc[1] + r * self.lib.sin(phi))

    def dpos(self, th):
        phi = self.phi0 + self.sigma * th
        u = self.u(th)
        du = self.du(th)
        r = 1 / u
        dr = -du / (u * u)
        c, s = self.lib.cos(phi), self.lib.sin(phi)
        # d(phi)/d(theta) = sigma
        return (dr * c - r * s * self.sigma, dr * s + r * c * self.sigma)

    def dtdth(self, th):
        u = self.u(th)
        return (self.W - self.kappa * u) / (self.b.c ** 2 * self.L * u * u)

    def momentum(self, th):
        lib = self.lib
        phi = self.phi0 + self.sigma * th
        pr = -self.L * self.du(th)
        pphi = self.sigma * self.L * self.u(th)
        px = pr * lib.cos(phi) - pphi * lib.sin(phi)
        py = pr * lib.sin(phi) + pphi * lib.cos(phi)
        return (px, py), pr, pphi


def time_closed(orb, th):
    """first_bend eq:time. I2 is eq:i2. I1 is the arctan integral, split at s = pi."""
    A, B = orb.A, orb.B
    G = orb.G2 ** 0.5
    s = G * th
    Delta = mp.sqrt(B ** 2 - A ** 2)

    def piece(bcoef, s_local):
        return 2 / Delta * mp.atan(mp.tan(s_local / 2) * mp.sqrt((B - bcoef) / (B + bcoef)))

    if s < mp.pi:
        I1 = piece(A, s)
    else:
        I1 = mp.pi / Delta + piece(-A, s - mp.pi)
    disc = A ** 2 - B ** 2
    I2 = A * mp.sin(s) / (disc * (A * mp.cos(s) + B)) - B / disc * I1
    t = (orb.W * I2 - orb.kappa * I1) / (orb.b.c ** 2 * orb.L * G)
    return t, I1, I2, s


def box_hit(x1, y1, x2, y2, lo, hi, entering):
    """The earliest face of the box that this chord crosses in the asked direction."""
    best = None
    faces = (("left", 0, lo[0]), ("right", 0, hi[0]), ("bottom", 1, lo[1]), ("top", 1, hi[1]))
    for name, axis, value in faces:
        p1 = x1 if axis == 0 else y1
        p2 = x2 if axis == 0 else y2
        d = p2 - p1
        if d == 0.0:
            continue
        t = (value - p1) / d
        if t < 0.0 or t > 1.0:
            continue
        o1 = y1 if axis == 0 else x1
        o2 = y2 if axis == 0 else x2
        other = o1 + t * (o2 - o1)
        olo = lo[1] if axis == 0 else lo[0]
        ohi = hi[1] if axis == 0 else hi[0]
        if other < olo - 1e-8 or other > ohi + 1e-8:
            continue
        inward = (d > 0.0 and name in ("left", "bottom")) or (d < 0.0 and name in ("right", "top"))
        if inward != entering:
            continue
        if best is None or t < best[0]:
            best = (t, name, axis, value)
    return best


def circle_hit(x1, y1, x2, y2, cx, cy, rad):
    dx, dy = x2 - x1, y2 - y1
    fx, fy = x1 - cx, y1 - cy
    A = dx * dx + dy * dy
    Bq = 2 * (fx * dx + fy * dy)
    C = fx * fx + fy * fy - rad * rad
    if C <= 0.0:
        return 0.0
    disc = Bq * Bq - 4 * A * C
    if disc < 0.0 or A == 0.0:
        return None
    t = (-Bq - math.sqrt(disc)) / (2 * A)
    if 0.0 <= t <= 1.0:
        return t
    return None


def first_event(orb, b, dth):
    """Chords of dth. Returns the winner and the largest chord-to-midpoint gap."""
    th = 0.0
    sag = 0.0
    # A bound orbit with B > |A| cannot escape; stop well past one turn.
    limit = 2.0 * math.pi / math.sqrt(orb.G2) + 1.0
    while th < limit:
        th2 = th + dth
        x1, y1 = orb.pos(th)
        x2, y2 = orb.pos(th2)
        xm, ym = orb.pos(0.5 * (th + th2))
        sag = max(sag, math.hypot(xm - 0.5 * (x1 + x2), ym - 0.5 * (y1 + y2)))
        t_col = circle_hit(x1, y1, x2, y2, orb.rc[0], orb.rc[1], b.R + b.a)
        world = box_hit(x1, y1, x2, y2, b.wmin, b.wmax, entering=False)
        det = box_hit(x1, y1, x2, y2, b.dmin, b.dmax, entering=True)
        cands = []
        if t_col is not None:
            cands.append((t_col, 0, "C", "sphere", None, None))
        if world is not None:
            cands.append((world[0], 1, "L", world[1], world[2], world[3]))
        if det is not None:
            cands.append((det[0], 2, "A", det[1], det[2], det[3]))
        if cands:
            cands.sort()
            _t, _rank, code, face, axis, value = cands[0]
            return {"code": code, "face": face, "axis": axis, "value": value,
                    "bracket": (th, th2), "sag": sag}
        th = th2
    raise RuntimeError("no event before one full turn")


def newton_face(orb, bracket, axis, value):
    lo, hi = mp.mpf(bracket[0]), mp.mpf(bracket[1])
    f = lambda th: orb.pos(th)[axis] - value
    if f(lo) * f(hi) > 0:
        raise RuntimeError("Newton bracket has no sign change")
    th = (lo + hi) / 2
    rows = []
    residual = f(th)
    for _ in range(12):
        df = orb.dpos(th)[axis]
        step = residual / df
        rows.append((th, residual, step))
        if abs(residual) < mp.mpf("1e-28"):
            break
        th = th - step
        residual = f(th)
    return th, residual, rows


def simpson(func, a, b, n):
    h = (b - a) / n
    total = func(a) + func(b)
    for i in range(1, n):
        total += (4 if i % 2 else 2) * func(a + i * h)
    return total * h / 3


def radiated(orbf, th_e, n):
    """Simpson on the Larmor power, evaluated in binary64 on the same orbit."""
    def power_dt(th):
        x, y = orbf.pos(th)
        rx, ry = x - orbf.rc[0], y - orbf.rc[1]
        r3 = (rx * rx + ry * ry) ** 1.5
        F = (orbf.kappa * rx / r3, orbf.kappa * ry / r3)
        (px, py), _, _ = orbf.momentum(th)
        gam = (orbf.W - orbf.kappa * orbf.u(th)) / (orbf.b.m * orbf.b.c ** 2)
        v = (px / (gam * orbf.b.m), py / (gam * orbf.b.m))
        vF = v[0] * F[0] + v[1] * F[1]
        acc = ((F[0] - v[0] * vF / orbf.b.c ** 2) / (gam * orbf.b.m),
               (F[1] - v[1] * vF / orbf.b.c ** 2) / (gam * orbf.b.m))
        a2 = acc[0] ** 2 + acc[1] ** 2
        vxa = v[0] * acc[1] - v[1] * acc[0]
        P = (2.0 / 3.0) * orbf.b.q ** 2 / orbf.b.c ** 3 * gam ** 6 * (a2 - vxa * vxa / orbf.b.c ** 2)
        return P * orbf.dtdth(th)

    fine = simpson(power_dt, 0.0, th_e, n)
    coarse = simpson(power_dt, 0.0, th_e, n // 2)
    rel = abs(fine - coarse) / abs(fine) / 15.0
    return fine, rel


def solve(b, Q, rc, dth):
    orb = Orbit(b, Q, rc, mp)
    # The chord search uses the same formulas at float precision.
    bf = FloatBoard(b)
    orbf = Orbit(bf, float(Q), [float(v) for v in rc], math)
    ev = first_event(orbf, bf, dth)
    if ev["code"] == "C":
        raise RuntimeError("a collision is not expected on these nine flights")
    th, residual, rows = newton_face(orb, ev["bracket"], ev["axis"], mp.mpf(ev["value"]))
    t, I1, I2, s = time_closed(orb, th)
    t_quad = mp.quad(lambda z: orb.dtdth(z), [0, th])
    pos = orb.pos(th)
    # The face coordinate is the face. The residual is the orbit's miss before that.
    snapped = [pos[0], pos[1]]
    snapped[ev["axis"]] = mp.mpf(ev["value"])
    (px, py), pr, pphi = orb.momentum(th)
    # Where x = x_Dhi, if the orbit reaches it before this event.
    cross = None
    x_face = float(b.dmax[0])
    th_f = 0.0
    while th_f < float(th):
        x1 = orbf.pos(th_f)[0]
        x2 = orbf.pos(th_f + dth)[0]
        if x1 >= x_face >= x2 and th_f > 1.0:
            cross = newton_face(orb, (th_f, th_f + dth), 0, b.dmax[0])
            break
        th_f += dth
    return {"orb": orb, "ev": ev, "th": th, "residual": residual, "rows": rows,
            "t": t, "I1": I1, "I2": I2, "s": s, "tQuad": t_quad,
            "x": snapped[0], "y": snapped[1], "px": px, "py": py, "pr": pr, "pphi": pphi,
            "cross": cross, "sag": ev["sag"], "orbf": orbf}


class FloatBoard:
    """The same board at binary64, for the chord search and the radiation sum."""

    def __init__(self, b):
        self.q = float(b.q)
        self.m = float(b.m)
        self.a = float(b.a)
        self.c = float(b.c)
        self.T0 = float(b.T0)
        self.R = float(b.R)
        self.tmax = float(b.tmax)
        self.gamma0 = float(b.gamma0)
        self.p0 = float(b.p0)
        self.x0 = [float(v) for v in b.x0]
        self.d = [float(v) for v in b.d]
        self.dmin = [float(v) for v in b.dmin]
        self.dmax = [float(v) for v in b.dmax]
        self.wmin = [float(v) for v in b.wmin]
        self.wmax = [float(v) for v in b.wmax]


def compute(doc):
    mp.mp.dps = DIGITS
    lvl = doc.level
    s = lvl["grid"]["subdivision"]
    places = [lvl["reference_solution"],
              doc.result(1)["placement"],
              doc.result(2)["placement"]]
    place_tags = ("Ref", "Weak", "Strong")
    shot_tags = ("Low", "Mid", "High")
    runs = {}
    sag = 0.0
    rad_err = 0.0
    for pi, place in enumerate(places):
        Q = mp.mpf(place[0]["value"])
        rc = [mp.mpf(v) / s for v in place[0]["node"][:2]]
        for si, shot in enumerate(lvl["shots"]):
            b = Board(lvl, shot, mp.mpf)
            sol = solve(b, Q, rc, DTH)
            n_rad = 4096
            e_rad, e_err = radiated(sol["orbf"], float(sol["th"]), n_rad)
            sol["rad"] = e_rad
            sol["radErr"] = e_err
            sol["board"] = b
            sag = max(sag, sol["sag"])
            rad_err = max(rad_err, float(e_err))
            doc.check(abs(sol["t"] - sol["tQuad"]) < mp.mpf("1e-30"),
                      f"{place_tags[pi]} {shot_tags[si]} time matches its quadrature")
            doc.check(abs(sol["residual"]) < mp.mpf("1e-24"),
                      f"{place_tags[pi]} {shot_tags[si]} sits on its face")
            doc.check(abs(sol["orb"].B) > abs(sol["orb"].A),
                      f"{place_tags[pi]} {shot_tags[si]} has B > |A|")
            doc.check(sol["s"] < 2 * mp.pi, f"{place_tags[pi]} {shot_tags[si]} stays below s = 2 pi")
            runs[(pi, si)] = sol
    # State tolerance: the exact-result floor, fixed above. Radiation tolerance: the
    # larger of 1e-4 and ten times the largest Simpson uncertainty, fixed here, before
    # any comparison.
    tol = TOL
    tol_rad = max(TOL_RAD, 10.0 * rad_err)

    refb = runs[(0, 1)]["board"]
    ref = runs[(0, 1)]["orb"]
    doc.value("q", refb.q, 3, sci=True)
    doc.value("m", refb.m, 3)
    doc.value("c", refb.c, 3)
    doc.value("R", refb.R, 3)
    doc.value("aRad", refb.a, 3)
    doc.value("tmax", refb.tmax, 3)
    doc.value("xLaunch", refb.x0[0], 3)
    doc.value("yLaunch", refb.x0[1], 3)
    doc.value("xC", ref.rc[0], 3)
    doc.value("yC", ref.rc[1], 3)
    doc.value("rZero", ref.r0, 3)
    doc.value("xWorld", refb.wmax[0], 3)
    doc.value("yWorld", refb.wmax[1], 3)
    doc.value("xDlo", refb.dmin[0], 3)
    doc.value("xDhi", refb.dmax[0], 3)
    doc.value("Q", ref.kappa / refb.q, 3, sci=True)
    doc.value("kappa", ref.kappa, 4)
    doc.value("kAbs", abs(ref.kappa), 3)
    doc.value("Qweak", runs[(1, 0)]["orb"].kappa / refb.q, 3, sci=True)
    doc.value("Qstrong", runs[(2, 0)]["orb"].kappa / refb.q, 3, sci=True)
    doc.value("kWeak", runs[(1, 0)]["orb"].kappa, 3)
    doc.value("kStrong", runs[(2, 0)]["orb"].kappa, 3)
    doc.value("kAbsWeak", abs(runs[(1, 0)]["orb"].kappa), 3)
    doc.value("kAbsStrong", abs(runs[(2, 0)]["orb"].kappa), 3)
    doc.value("phiZero", ref.phi0, 8)
    doc.value("sigma", ref.sigma)
    doc.value("duZero", ref.du0, 2, sci=True)
    doc.value("prZero", ref.pr0, 2, sci=True)
    doc.value("dth", DTH, 2, sci=True)
    doc.value("sagMax", sag, 2, sci=True)
    doc.value("nChord", int(round(2 * math.pi / DTH)))
    doc.value("nSimp", 4096)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("tolRad", tol_rad, 2, sci=True)
    doc.value("radErrMax", rad_err, 2, sci=True)
    doc.value("maxCharges", lvl["limits"]["max_charges"])
    region = lvl["limits"]["region"]
    doc.value("xRegLo", region["min"][0], 3)
    doc.value("xRegHi", region["max"][0], 3)
    doc.value("yRegLo", region["min"][1], 3)
    doc.value("yRegHi", region["max"][1], 3)

    # The three reference orbits, one per energy.
    for si, tag in enumerate(shot_tags):
        sol = runs[(0, si)]
        o, b = sol["orb"], sol["board"]
        v0 = b.p0 / (b.gamma0 * b.m)
        circ = b.p0 * v0 * o.r0
        circ2 = b.gamma0 * b.m * v0 ** 2 * o.r0
        doc.value("T" + tag, b.T0, 3)
        doc.value("gamma" + tag, b.gamma0, 8)
        doc.value("p" + tag, b.p0, 10)
        doc.value("v" + tag, v0, 8)
        doc.value("circ" + tag, circ, 8)
        doc.value("circAlt" + tag, circ2, 8)
        doc.value("circGap" + tag, circ - abs(o.kappa), 6)
        doc.value("VZero" + tag, o.kappa / o.r0, 6)
        doc.value("Wbind" + tag, b.T0 + o.kappa / o.r0, 6)
        doc.value("W" + tag, o.W, 10)
        doc.value("L" + tag, o.L, 10)
        doc.value("G" + tag, mp.sqrt(o.G2), 8)
        doc.value("A" + tag, o.A, 8)
        doc.value("B" + tag, o.B, 8)
        doc.value("oneG" + tag, 1 - mp.sqrt(o.G2), 6)
        A_energy = mp.sqrt(o.B ** 2 + (o.W ** 2 - b.m ** 2 * b.c ** 4) / (b.c ** 2 * o.L ** 2 * o.G2))
        doc.value("Adiff" + tag, abs(abs(o.A) - A_energy), 2, sci=True)
        doc.value("yDlo" + tag, b.dmin[1], 3)
        doc.value("yDhi" + tag, b.dmax[1], 3)
        doc.value("th" + tag, sol["th"], 10)
        doc.value("past" + tag, (sol["th"] - mp.pi) * 180 / mp.pi, 6)
        doc.value("s" + tag, sol["s"], 8)
        doc.value("faceRes" + tag, sol["residual"], 2, sci=True)
        doc.value("ang" + tag, mp.atan2(sol["py"], sol["px"]) * 180 / mp.pi, 6)
        ue = o.u(sol["th"])
        doc.value("uE" + tag, ue, 8)
        doc.value("rE" + tag, 1 / ue, 8)
        # Minimum radius on the way to the face. A < 0: u grows, the end is the closest
        # if Gamma theta < pi; A > 0: the launch is the closest.
        us = [o.u(mp.mpf(0)), ue]
        if sol["s"] >= mp.pi:
            us.append(o.B - o.A)
        rmin = 1 / max(us)
        doc.value("rMin" + tag, rmin, 6)
        doc.value("clear" + tag, rmin - b.R, 6)
        Tarr = (o.W - b.m * b.c ** 2) - o.kappa * ue
        doc.value("Tarr" + tag, Tarr, 6)
        doc.data("path" + tag.lower(),
                 [[float(v) for v in o.pos(sol["th"] * k / 200)] for k in range(201)])

    # The middle energy, in the detail the derivation evaluates by hand.
    mid = runs[(0, 1)]
    doc.value("IOne", mid["I1"], 10)
    doc.value("ITwo", mid["I2"], 10)
    doc.value("tQuadDiff", abs(mid["t"] - mid["tQuad"]), 2, sci=True)
    doc.value("prE", mid["pr"], 8)
    doc.value("pphiE", mid["pphi"], 8)
    doc.table("newton",
              [r"$\theta$", r"$f(\theta)$", r"$\Delta\theta$"],
              [(float(th), float(res), float(step)) for th, res, step in mid["rows"]],
              digits=[12, 3, 3],
              sci=[False, True, True])
    doc.value("nNewton", len(mid["rows"]))

    codes = {"A": "Arrived", "L": "LeftBounds", "C": "Collided"}
    rad_frac = 0.0
    for pi, ptag in enumerate(place_tags):
        for si, stag in enumerate(shot_tags):
            sol = runs[(pi, si)]
            name = ptag + stag
            doc.text(name + "Outcome", codes[sol["ev"]["code"]])
            doc.text(name + "Face", sol["ev"]["face"])
            doc.value(name + "T", sol["t"], 8)
            doc.value(name + "X", sol["x"], 8)
            doc.value(name + "Y", sol["y"], 8)
            doc.value(name + "Px", sol["px"], 8)
            doc.value(name + "Py", sol["py"], 8)
            doc.value(name + "Ang", mp.atan2(sol["py"], sol["px"]) * 180 / mp.pi, 6)
            doc.value(name + "Rad", sol["rad"], 4, sci=True)
            doc.value(name + "RadErr", sol["radErr"], 2, sci=True)
            rad_frac = max(rad_frac, float(sol["rad"] / sol["board"].T0))
            if sol["cross"] is None:
                doc.text(name + "Cross", "no")
                doc.value(name + "CrossY", 0, 3)
            else:
                cy = sol["orb"].pos(sol["cross"][0])[1]
                doc.text(name + "Cross", "yes")
                doc.value(name + "CrossY", cy, 6)
    doc.value("radFracMax", rad_frac, 2, sci=True)
    doc.value("arrowX", refb.x0[0] + mp.mpf("1.6"), 3)
    doc.value("arrowY", refb.x0[1], 3)

    # The two misses drawn on the figure: the weak high energy, which meets the floor,
    # and the strong low energy, which meets the left edge of the world.
    for key, dest in (( (1, 2), "pathfloor"), ( (2, 0), "pathleft")):
        sol = runs[key]
        doc.data(dest, [[float(v) for v in sol["orb"].pos(sol["th"] * k / 160)] for k in range(161)])

    def put(pi, si):
        sol = runs[(pi, si)]
        g = doc.flight(pi, flight=si)
        ge = g["end"]
        title = {"Ref": "ref.", "Weak": "weak", "Strong": "strong"}[place_tags[pi]]
        title += " " + {"Low": "low", "Mid": "mid", "High": "high"}[shot_tags[si]]
        tag = place_tags[pi] + shot_tags[si]
        doc.compare(tag + "Outcome", title + " outcome", codes[sol["ev"]["code"]], g["outcome"])
        doc.compare(tag + "T", title + " $t$", sol["t"], ge["t"], tol)
        doc.compare(tag + "X", title + " $x$", sol["x"], ge["x"][0], tol)
        doc.compare(tag + "Y", title + " $y$", sol["y"], ge["x"][1], tol)
        if pi == 0:
            doc.compare(tag + "Px", title + " $p_x$", sol["px"], ge["p"][0], tol)
            doc.compare(tag + "Py", title + " $p_y$", sol["py"], ge["p"][1], tol)
            doc.compare(tag + "Rad", title + " $E$ (rel.)", sol["rad"],
                        g["radiated_energy"], tol_rad, rel=True)

    for pi in range(3):
        for si in range(3):
            put(pi, si)
