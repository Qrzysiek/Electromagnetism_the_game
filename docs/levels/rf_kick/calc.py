"""Companion of docs/levels/rf_kick/doc.tex (Level 39, RF kick).

Every number of the document, from the level file. The antenna is a point dipole
p(t) = p0 cos(omega t) with its exact retarded fields (eq:fields, derived in the document from
the potentials eq:potentials); the field formula is checked against Maxwell's equations by
central differences (sec:maxwell-check). The flight: the relativistic equations of motion
(Foundations eq:eom) marched by the classical Runge-Kutta method at a fixed step (Foundations
eq:rk4) with two more states, the field's work on the ion (eq:work) and Lienard's radiated
energy (Foundations eq:larmor); the entry located on the cubic Hermite interpolant by
bisection; the error against the step 2h. The kick to first order in the field (eq:first-order)
by the same method along the undeflected line. Compared with the game: the reference, the empty
board and three more placements (placements.json): the antenna at (4, 8) in either phase and
at (9, 8).
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-6      # the radiated energy: the game's three-point Gauss-Legendre per step (sixth
                    # order, steps ~0.65, the integrand varying at 2 omega) is ~1e-7 relative; ten times
H = 0.004
TABLE_EVERY = 4.0   # the march's table: a row every 4 time units


class Antenna:
    """eq:fields: the retarded fields of p(t) = p0 cos(omega t) (3-vectors, any point)."""

    def __init__(self, e, s, omega, c):
        a = math.radians(e["angle_deg"])
        self.pos = (e["node"][0] / s, e["node"][1] / s, 0.0)
        self.amp = (math.cos(a) * e["value"], math.sin(a) * e["value"], 0.0)
        self.omega, self.c = omega, c

    def moment(self, t):
        co, si, w = math.cos(self.omega * t), math.sin(self.omega * t), self.omega
        return ([a * co for a in self.amp], [-w * a * si for a in self.amp],
                [-w * w * a * co for a in self.amp])

    def fields(self, x, t):
        d = [x[k] - self.pos[k] for k in range(3)]
        r = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        n = [v / r for v in d]
        c = self.c
        p, pd, pdd = self.moment(t - r / c)                   # the retarded time t - r/c
        dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        cross = lambda a, b: (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        npd, npp, nq = dot(n, p), dot(n, pd), dot(n, pdd)
        e = [(3 * n[k] * npd - p[k]) / r ** 3 + (3 * n[k] * npp - pd[k]) / (c * r * r)
             + (n[k] * nq - pdd[k]) / (c * c * r) for k in range(3)]   # n x (n x a) = n(n.a) - a
        b1, b2 = cross(pd, n), cross(pdd, n)
        b = [(b1[k] / (r * r) + b2[k] / (c * r)) / (c * c) for k in range(3)]
        return e, b

    def potentials(self, x, t):
        """eq:potentials: phi = n.[p]/r^2 + n.[p']/(c r), A = [p']/(c^2 r)."""
        d = [x[k] - self.pos[k] for k in range(3)]
        r = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        n = [v / r for v in d]
        p, pd, _ = self.moment(t - r / self.c)
        phi = sum(n[k] * p[k] for k in range(3)) / (r * r) + sum(n[k] * pd[k] for k in range(3)) / (self.c * r)
        return phi, [pd[k] / (self.c ** 2 * r) for k in range(3)]


class Setup:
    def __init__(self, lvl, placement):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.c = ph["c"]
        self.omega = ph["rf_omega"]
        self.R = ph["antenna_radius"]
        self.tmax = ph["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        els = [e for e in lvl["elements"] if e["kind"] == "antenna"] + list(placement)
        self.antennas = [Antenna(e, s, e.get("omega") or self.omega, self.c) for e in els]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)                 # Foundations eq:launch
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)

    def field(self, x, y, t):
        ex = ey = bz = 0.0
        for a in self.antennas:
            e, b = a.fields((x, y, 0.0), t)
            ex, ey, bz = ex + e[0], ey + e[1], bz + b[2]
        return ex, ey, bz

    def rhs(self, s, t):
        """State (x, y, px, py, work, radiated): Foundations eq:eom, eq:work, eq:larmor."""
        x, y, px, py = s[:4]
        m, c, q = self.m, self.c, self.q
        g = math.sqrt(1 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        ex, ey, bz = self.field(x, y, t)
        fx, fy = q * (ex + vy * bz), q * (ey - vx * bz)
        # dv/dt = (F - v (v.F)/c^2) / (gamma m), then Lienard's power
        vf = vx * fx + vy * fy
        ax, ay = (fx - vx * vf / c ** 2) / (g * m), (fy - vy * vf / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        power = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        return (vx, vy, fx, fy, q * (ex * vx + ey * vy), power)

    def kinetic(self, px, py):
        p2 = px * px + py * py
        g = math.sqrt(1 + p2 / (self.m * self.c) ** 2)
        return p2 / (self.m * (g + 1))                      # (gamma - 1) m c^2


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


def rk4_step(st, s, f, t, h):
    k1 = f
    k2 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)), t + h / 2)
    k3 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)), t + h / 2)
    k4 = st.rhs(tuple(a + h * b for a, b in zip(s, k3)), t + h)
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def fly(st, h, table=False):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0, 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    wdev = 0.0
    path = [(s[0], s[1])]
    rows = [(0.0,) + s[:4]] if table else []
    preds = [(f"Collided({i})", 0, lambda p, a=a: math.hypot(p[0] - a.pos[0], p[1] - a.pos[1]) <= st.R)
             for i, a in enumerate(st.antennas)]
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
    every = round(TABLE_EVERY / h)
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, t, hs)
        f2 = st.rhs(s2, t + hs)
        pts = [hermite(s, f, s2, f2, hs, u) for u in us]
        events = []
        for name, rank, pred in preds:
            for k, p in enumerate(pts):
                if pred(p):
                    events.append((us[k], rank, name, pred))
                    break
        if events:
            u_hit, _r, outcome, pred = min(events, key=lambda e: (e[0], e[1]))
            lo, hi = u_hit - 1 / 64, u_hit
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if pred(hermite(s, f, s2, f2, hs, mid)):
                    hi = mid
                else:
                    lo = mid
            te = t + hi * hs
            se = rk4_step(st, s, f, t, te - t)
            pe = hermite(s, f, s2, f2, hs, hi)
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else (
                face_of(pe, st.wmin, st.wmax) if outcome == "LeftBounds" else "")
            path.append(pe)
            wdev = max(wdev, abs(st.kinetic(se[2], se[3]) - st.T0 - se[4]))
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "steps": n + 1,
                    "wdev": wdev, "work": se[4], "rad": se[5], "path": path, "rows": rows}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.kinetic(s[2], s[3]) - st.T0 - s[4]))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if table and n % every == 0:
            rows.append((t,) + s[:4])
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "",
                    "steps": n, "wdev": wdev, "work": s[4], "rad": s[5], "path": path, "rows": rows}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def first_order(st, h):
    """eq:first-order: along the undeflected line x = x0 + v0 t, y = y0, the transverse force
    q(E_y - v0 B_z) integrated twice, by RK4 on (p_y1, y1), to the time the line reaches the
    detector's near face. Returns that time, p_y1 and y1 there."""
    v0 = st.p0 / (st.gamma0 * st.m)
    t_star = (st.dmin[0] - st.x0[0]) / v0

    def g(t, s):
        ex, ey, bz = st.field(st.x0[0] + v0 * t, st.x0[1], t)
        return (st.q * (ey - v0 * bz), s[0] / (st.gamma0 * st.m))

    n = round(t_star / h)
    hh = t_star / n
    s, t = (0.0, st.x0[1]), 0.0
    for _ in range(n):
        k1 = g(t, s)
        k2 = g(t + hh / 2, tuple(a + hh / 2 * b for a, b in zip(s, k1)))
        k3 = g(t + hh / 2, tuple(a + hh / 2 * b for a, b in zip(s, k2)))
        k4 = g(t + hh, tuple(a + hh * b for a, b in zip(s, k3)))
        s = tuple(a + hh / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))
        t += hh
    return t_star, s[0], s[1]


def maxwell_residuals(ant, x, t, d=1e-4):
    """sec:maxwell-check: div E, div B, curl E + dB/dt, curl B - dE/dt / c^2 by central
    differences, each relative to the size of its largest term."""
    def ef(xx, tt):
        return ant.fields(xx, tt)

    def shifted(k, sgn):
        y = list(x)
        y[k] += sgn * d
        return ef(y, t)

    der = [[[(shifted(k, 1)[w][i] - shifted(k, -1)[w][i]) / (2 * d) for k in range(3)]
            for i in range(3)] for w in range(2)]                # der[w][component][direction]
    e_p, e_m = ef(x, t + d), ef(x, t - d)
    dt = [[(e_p[w][i] - e_m[w][i]) / (2 * d) for i in range(3)] for w in range(2)]
    c2 = ant.c ** 2

    def curl(w):
        g = der[w]
        return (g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1])

    # Each residual relative to the largest single derivative entering it: the near field's
    # curl is a cancellation of terms ~100 times its result, and the differences' error
    # (about (d/r)^2 times a few) is relative to those terms.
    scale_e = max(abs(der[0][i][k]) for i in range(3) for k in range(3))
    scale_b = max(abs(der[1][i][k]) for i in range(3) for k in range(3))
    ce, cb = curl(0), curl(1)
    div_e = abs(sum(der[0][i][i] for i in range(3))) / scale_e
    div_b = abs(sum(der[1][i][i] for i in range(3))) / scale_b
    far = max(abs(ce[i] + dt[1][i]) for i in range(3)) / max(scale_e, max(abs(v) for v in dt[1]))
    amp = max(abs(cb[i] - dt[0][i] / c2) for i in range(3)) / max(scale_b, max(abs(v) / c2 for v in dt[0]))
    return div_e, div_b, far, amp


def potential_residual(ant, x, t, d=1e-4):
    """eq:fields from eq:potentials: E = -grad phi - dA/dt, B = curl A by central differences,
    relative to the largest component of each field."""
    def shifted(k, sgn):
        y = list(x)
        y[k] += sgn * d
        return ant.potentials(y, t)

    grad = [(shifted(k, 1)[0] - shifted(k, -1)[0]) / (2 * d) for k in range(3)]
    da = [[(shifted(k, 1)[1][i] - shifted(k, -1)[1][i]) / (2 * d) for k in range(3)] for i in range(3)]
    a_p, a_m = ant.potentials(x, t + d)[1], ant.potentials(x, t - d)[1]
    e_pot = [-grad[i] - (a_p[i] - a_m[i]) / (2 * d) for i in range(3)]
    b_pot = (da[2][1] - da[1][2], da[0][2] - da[2][0], da[1][0] - da[0][1])
    e, b = ant.fields(x, t)
    re = max(abs(e_pot[i] - e[i]) for i in range(3)) / max(map(abs, e))
    rb = max(abs(b_pot[i] - b[i]) for i in range(3)) / max(map(abs, b))
    return re, rb


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(lvl["physics"]["c"] is not None and not lvl.get("disturbances") and len(ref) == 1
              and not lvl["physics"].get("radiation_reaction", False),
              "relativistic, no disturbance, one placed antenna, no radiation reaction")
    st = Setup(lvl, ref)
    (ant,) = st.antennas
    e = ref[0]
    for k_, v in (("q", st.q), ("m", st.m), ("Tzero", st.T0), ("c", st.c), ("omega", st.omega),
                  ("Ra", st.R), ("xZero", st.x0[0]), ("yZero", st.x0[1]), ("xDetLo", st.dmin[0]),
                  ("yDetLo", st.dmin[1]), ("xDetHi", st.dmax[0]), ("yDetHi", st.dmax[1]),
                  ("xWorld", st.wmax[0]), ("yWorld", st.wmax[1]), ("xA", ant.pos[0]), ("yA", ant.pos[1]),
                  ("pA", e["value"]), ("angleA", e["angle_deg"]), ("H", H), ("tMax", st.tmax)):
        doc.value(k_, v, 4, sci=k_ in ("q", "pA"))
    doc.value("period", 2 * math.pi / st.omega, 4)
    doc.value("wavelength", 2 * math.pi * st.c / st.omega, 3)
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p0, 10)
    v0 = st.p0 / (st.gamma0 * st.m)
    doc.value("vZero", v0, 10)
    doc.value("kr", st.omega / st.c * math.hypot(st.x0[0] - ant.pos[0], st.x0[1] - ant.pos[1]), 3)
    doc.value("rLaunch", math.hypot(st.x0[0] - ant.pos[0], st.x0[1] - ant.pos[1]), 4)
    doc.value("ampCross", math.cos(math.radians(e["angle_deg"])) * e["value"], 2, sci=True)

    # The field at the launch, t = 0, term by term (the hand procedure's first step).
    d = [st.x0[0] - ant.pos[0], st.x0[1] - ant.pos[1]]
    r0 = math.hypot(*d)
    tr = -r0 / st.c
    doc.value("tRet", tr, 6)
    doc.value("phaseRet", st.omega * tr, 6)
    p, pd, pdd = ant.moment(tr)
    ny = d[1] / r0
    near = (3 * ny * ny * p[1] - p[1]) / r0 ** 3
    ind = (3 * ny * ny * pd[1] - pd[1]) / (st.c * r0 * r0)
    rad = (ny * ny * pdd[1] - pdd[1]) / (st.c * st.c * r0)
    ex, ey, bz = st.field(st.x0[0], st.x0[1], 0.0)
    doc.check(abs(near + ind + rad - ey) < 1e-9 * abs(ey), "the launch field term by term")
    doc.value("pRet", p[1], 6, sci=True)
    doc.value("pdRet", pd[1], 6, sci=True)
    doc.value("eNear", near, 6, sci=True)
    doc.value("eInd", ind, 6, sci=True)
    doc.value("eRad", rad, 6, sci=True)
    doc.value("eYLaunch", ey, 6, sci=True)
    doc.value("bZLaunch", bz, 4, sci=True)
    doc.value("fLaunch", st.q * (ey - v0 * bz), 6)
    doc.value("accLaunch", st.q * (ey - v0 * bz) / (st.gamma0 * st.m), 6)

    # The field formula against Maxwell's equations, at three points and times.
    worst = 0.0
    for xx, tt in (((7.3, 13.1, 0.7), 3.1), ((20.0, 2.0, -1.5), 11.7), ((1.0, 6.2, 0.2), 0.4)):
        worst = max(worst, *maxwell_residuals(ant, xx, tt))
    doc.check(worst < 1e-6, "the field satisfies Maxwell's equations")
    doc.value("maxwellRes", worst, 1, sci=True)
    worst = 0.0
    for xx, tt in (((7.3, 13.1, 0.7), 3.1), ((20.0, 2.0, -1.5), 11.7), ((1.0, 6.2, 0.2), 0.4)):
        worst = max(worst, *potential_residual(ant, xx, tt))
    doc.check(worst < 1e-6, "the fields are the potentials' derivatives")
    doc.value("potentialRes", worst, 1, sci=True)

    # The reference: the march, its error, the work and radiation checks.
    run, coarse = fly(st, H, table=True), fly(st, 2 * H)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    err = err_of(run, coarse)
    doc.value("tRef", run["t"], 10)
    doc.value("xRef", run["x"][0], 10)
    doc.value("yRef", run["x"][1], 10)
    doc.value("pxRef", run["p"][0], 10)
    doc.value("pyRef", run["p"][1], 10)
    doc.text("faceRef", run["face"])
    doc.value("errRef", err, 2, sci=True)
    doc.value("wDevRef", run["wdev"], 2, sci=True)
    doc.value("workRef", run["work"], 6)
    doc.value("radRef", run["rad"], 3, sci=True)
    doc.value("radFrac", run["rad"] / st.T0, 2, sci=True)
    doc.check(run["rad"] < 1e-10 * st.T0, "the neglected radiation is below 1e-10 T0")
    doc.check(run["wdev"] < 1e-11, "the work equals the change of the kinetic energy")
    doc.data("pathref", run["path"])
    doc.table("march", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$"],
              [r for r in run["rows"]] + [(run["t"], run["x"][0], run["x"][1], run["p"][0], run["p"][1])],
              digits=[4, 8, 8, 8, 8])

    # The kick to first order in the field.
    t_star, py1, y1 = first_order(st, H)
    t_star2, py1b, y1b = first_order(st, 2 * H)
    doc.check(abs(y1 - y1b) < 1e-9, "the first-order integrals converge")
    doc.value("tStar", t_star, 6)
    doc.value("pyFirst", py1, 6)
    doc.value("yFirst", y1, 6)
    doc.value("firstDiff", abs(y1 - run["x"][1]), 2)
    doc.value("firstRel", abs(y1 - run["x"][1]) / abs(run["x"][1] - st.x0[1]), 2)
    doc.value("liftRef", run["x"][1] - st.x0[1], 4)

    g = doc.flight(0)
    gp = doc.flight(0, run="preview")
    gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][1] - g["end"]["x"][1]),
              abs(gp["end"]["p"][0] - g["end"]["p"][0]), abs(gp["end"]["p"][1] - g["end"]["p"][1]))
    tol = max(TOL, 10 * err, gap)
    doc.compare("Out", "outcome", run["outcome"], g["outcome"])
    doc.compare("T", r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
    doc.compare("X", r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
    doc.compare("Y", r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
    doc.compare("Rad", "radiated energy (relative)", run["rad"], g["radiated_energy"], RAD_REL, rel=True)
    largest = max(abs(run["t"] - g["end"]["t"]), abs(run["x"][1] - g["end"]["x"][1]),
                  abs(run["p"][0] - g["end"]["p"][0]), abs(run["p"][1] - g["end"]["p"][1]))

    # The empty board: a straight line along y = y0 to the world's right edge.
    rb = fly(Setup(lvl, []), H)
    gb = doc.flight(1)
    doc.compare("BareOut", "empty board: outcome", rb["outcome"], gb["outcome"])
    doc.compare("BareT", r"empty board: $t_\mathrm{e}$", rb["t"], gb["end"]["t"], TOL)
    doc.value("bareT", rb["t"], 8)
    doc.check(abs(rb["t"] - (st.wmax[0] - st.x0[0]) / v0) < 1e-9, "the empty board: t = (x_max - x0)/v0")

    # Three more placements: the antenna at (4, 8) in either phase, and at (9, 8).
    for k, tag in ((2, "Anti"), (3, "Late"), (4, "Same")):
        placement = doc.result(k)["placement"]
        (pe,) = placement
        sk = Setup(lvl, placement)
        doc.value(f"x{tag}", pe["node"][0], 4)
        doc.value(f"y{tag}", pe["node"][1], 4)
        doc.value(f"p{tag}", pe["value"], 4, sci=True)
        # when the undeflected line passes the antenna, and the drive's phase then (retarded)
        t_pass = (pe["node"][0] - st.x0[0]) / v0
        r_pass = abs(st.x0[1] - pe["node"][1])
        doc.value(f"tPass{tag}", t_pass, 3)
        doc.value(f"cosPass{tag}", math.cos(st.omega * (t_pass - r_pass / st.c)), 2)
        rk, ck = fly(sk, H), fly(sk, 2 * H)
        ek = err_of(rk, ck)
        gk = doc.flight(k)
        gkp = doc.flight(k, run="preview")
        gapk = max(abs(gkp["end"]["t"] - gk["end"]["t"]), abs(gkp["end"]["x"][0] - gk["end"]["x"][0]),
                   abs(gkp["end"]["x"][1] - gk["end"]["x"][1]))
        tolk = max(TOL, 10 * ek, gapk)
        if tag == "Same":
            doc.game_value("gapSame", gapk, 2, sci=True)
        title = f"$p_0 = {latex_p(pe['value'])}$ at ({pe['node'][0]}, {pe['node'][1]}): "
        doc.compare(f"{tag}Out", title + "outcome", rk["outcome"], gk["outcome"])
        doc.compare(f"{tag}X", title + r"$x_\mathrm{e}$", rk["x"][0], gk["end"]["x"][0], tolk)
        doc.compare(f"{tag}Y", title + r"$y_\mathrm{e}$", rk["x"][1], gk["end"]["x"][1], tolk)
        doc.text(f"out{tag}", rk["outcome"])
        doc.value(f"xEnd{tag}", rk["x"][0], 5)
        doc.value(f"yEnd{tag}", rk["x"][1], 5)
        doc.data("path" + tag.lower(), rk["path"])
        largest = max(largest, abs(rk["x"][0] - gk["end"]["x"][0]), abs(rk["x"][1] - gk["end"]["x"][1]))
    doc.game_value("largest", largest, 2, sci=True)


def latex_p(v):
    """A placed amplitude for the table's labels in LaTeX: -10^{6}, 2 x 10^{6}."""
    mant, exp = f"{v:.1e}".split("e")
    mant = mant.rstrip("0").rstrip(".")
    lead = {"1": "", "-1": "-"}.get(mant, mant + r" \times ")
    return f"{lead}10^{{{int(exp)}}}"
