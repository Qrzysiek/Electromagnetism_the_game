"""Companion of docs/levels/rf_separator/doc.tex (Level 45, RF separator).

Every number of the document, from the level file. The antenna's exact retarded fields (rf_kick
eq:fields); each ion's flight from its launch time t0, the fields taken at t0 + t (eq:offset);
the relativistic equations of motion (Foundations eq:eom) marched by the classical Runge-Kutta
method at a fixed step (Foundations eq:rk4) with Lienard's radiated energy as one more state;
the entry by cubic Hermite interpolation and bisection; the error against the step 2h. The exact
symmetry (eq:conjugate): the second ion, half a period late, is the first with its charge
reversed, checked by marching both. The kicks to first order in the field (rf_kick
eq:first-order), equal and opposite. Compared with the game: the reference and the empty board,
both ions.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-6      # the radiated energy, as in rf_kick: the game's Gauss-Legendre per step, ~1e-7; ten times
H = 0.004


class Antenna:
    """rf_kick eq:fields: the retarded fields of p(t) = p0 cos(omega t)."""

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
        p, pd, pdd = self.moment(t - r / c)
        dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        cross = lambda a, b: (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        npd, npp, nq = dot(n, p), dot(n, pd), dot(n, pdd)
        e = [(3 * n[k] * npd - p[k]) / r ** 3 + (3 * n[k] * npp - pd[k]) / (c * r * r)
             + (n[k] * nq - pdd[k]) / (c * c * r) for k in range(3)]
        b1, b2 = cross(pd, n), cross(pdd, n)
        b = [(b1[k] / (r * r) + b2[k] / (c * r)) / (c * c) for k in range(3)]
        return e, b


class Setup:
    def __init__(self, lvl, placement, shot_i, t0=None, charge_sign=1.0):
        shot = lvl["shots"][shot_i]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = charge_sign * shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.t0 = shot["launch"].get("time", 0.0) if t0 is None else t0
        self.c = ph["c"]
        self.omega = ph["rf_omega"]
        self.R = ph["antenna_radius"]
        self.Rq = ph["charge_radius"]
        self.tmax = ph["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        els = list(lvl["elements"]) + list(placement)
        assert all(e["kind"] == "antenna" for e in els)
        self.antennas = [Antenna(e, s, e.get("omega") or self.omega, self.c) for e in els]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)

    def field(self, x, y, t):
        """eq:offset: the fields at the lab time t0 + t."""
        ex = ey = bz = 0.0
        for a in self.antennas:
            e, b = a.fields((x, y, 0.0), self.t0 + t)
            ex, ey, bz = ex + e[0], ey + e[1], bz + b[2]
        return ex, ey, bz

    def rhs(self, s, t):
        x, y, px, py = s[:4]
        m, c, q = self.m, self.c, self.q
        g = math.sqrt(1 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        ex, ey, bz = self.field(x, y, t)
        fx, fy = q * (ex + vy * bz), q * (ey - vx * bz)
        vf = vx * fx + vy * fy
        ax, ay = (fx - vx * vf / c ** 2) / (g * m), (fy - vy * vf / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        power = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        return (vx, vy, fx, fy, power)


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


def fly(st, h):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    path = [(s[0], s[1])]
    preds = [(f"Collided({i})", 0, lambda p, a=a: math.hypot(p[0] - a.pos[0], p[1] - a.pos[1]) <= st.R)
             for i, a in enumerate(st.antennas)]
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
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
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else ""
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "path": path, "rad": se[4]}
        s, f, t, n = s2, f2, t + hs, n + 1
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "path": path, "rad": s[4]}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def first_order(st, h):
    """rf_kick eq:first-order with the launch time: along x = x0 + v0 t, y = y0, the transverse
    force q(E_y - v0 B_z) at t0 + t integrated twice to the detector's near face."""
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


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(lvl["physics"]["c"] is not None and len(lvl["shots"]) == 2 and len(ref) == 1
              and not lvl.get("disturbances") and not lvl["physics"].get("radiation_reaction"),
              "relativistic, two shots, one antenna, no disturbance, no radiation reaction")
    sts = [Setup(lvl, ref, i) for i in range(2)]
    a, b = sts
    (ant,) = a.antennas
    e = ref[0]
    period = 2 * math.pi / a.omega
    for k_, v, dg in (("Tzero", a.T0, 4), ("m", a.m, 4), ("c", a.c, 4), ("omega", a.omega, 4), ("Ra", a.R, 4),
                      ("xZero", a.x0[0], 4), ("yZero", a.x0[1], 4), ("xWorld", a.wmax[0], 4),
                      ("yWorld", a.wmax[1], 4), ("xA", ant.pos[0], 4), ("yA", ant.pos[1], 4),
                      ("angleA", e["angle_deg"], 4), ("H", H, 4), ("tMax", a.tmax, 4),
                      ("tLaunch", b.t0, 10), ("period", period, 6), ("vZero", a.p0 / (a.gamma0 * a.m), 10),
                      ("gammaZero", a.gamma0, 10), ("pZero", a.p0, 10)):
        doc.value(k_, v, dg)
    doc.value("q", a.q, 4, sci=True)
    doc.value("pA", e["value"], 4, sci=True)
    for i, (st, nm) in enumerate(((a, "One"), (b, "Two"))):
        for k_, v in (("xLo", st.dmin[0]), ("yLo", st.dmin[1]), ("xHi", st.dmax[0]), ("yHi", st.dmax[1])):
            doc.value(f"{k_}{nm}", v, 4)
    doc.value("halfRatio", b.t0 / period, 12)
    doc.check(abs(a.omega * b.t0 - math.pi) < 1e-14, "the second launch half a period later")
    lim = lvl["limits"]
    for k_, v in (("xRegLo", lim["region"]["min"][0]), ("yRegLo", lim["region"]["min"][1]),
                  ("xRegHi", lim["region"]["max"][0]), ("yRegHi", lim["region"]["max"][1])):
        doc.value(k_, v, 4)

    largest = 0.0
    runs = []
    for i, (st, nm) in enumerate(((a, "One"), (b, "Two"))):
        run, coarse = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", f"ion {i + 1} arrives")
        err = err_of(run, coarse)
        runs.append(run)
        doc.value(f"t{nm}", run["t"], 10)
        doc.value(f"x{nm}", run["x"][0], 10)
        doc.value(f"y{nm}", run["x"][1], 10)
        doc.value(f"px{nm}", run["p"][0], 8)
        doc.value(f"py{nm}", run["p"][1], 8)
        doc.text(f"face{nm}", run["face"])
        doc.value(f"err{nm}", err, 2, sci=True)
        doc.value(f"rad{nm}", run["rad"], 2, sci=True)
        doc.data("path" + nm.lower(), run["path"])
        g, gp = doc.flight(0, flight=i), doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]), abs(gp["end"]["p"][0] - g["end"]["p"][0]),
                  abs(gp["end"]["p"][1] - g["end"]["p"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"ion {i + 1}: "
        doc.compare(f"Out{nm}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{nm}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{nm}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{nm}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Px{nm}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
        doc.compare(f"Py{nm}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        doc.compare(f"Rad{nm}", title + "radiated energy (relative)", run["rad"], g["radiated_energy"], RAD_REL, rel=True)
        largest = max(largest, abs(run["t"] - g["end"]["t"]), abs(run["x"][1] - g["end"]["x"][1]),
                      abs(run["p"][0] - g["end"]["p"][0]), abs(run["p"][1] - g["end"]["p"][1]))
        # the empty board: straight to the world's right edge
        gb = doc.flight(1, flight=i)
        t_free = (st.wmax[0] - st.x0[0]) / (st.p0 / (st.gamma0 * st.m))
        doc.compare(f"BareOut{nm}", f"empty board, ion {i + 1}: outcome", "LeftBounds", gb["outcome"])
        doc.compare(f"BareT{nm}", f"empty board, ion {i + 1}: " + r"$t_\mathrm{e}$", t_free, gb["end"]["t"], TOL)
    doc.game_value("largest", largest, 2, sci=True)
    doc.value("bareT", (a.wmax[0] - a.x0[0]) / (a.p0 / (a.gamma0 * a.m)), 8)

    # The symmetry: ion 2 is ion 1 with its charge reversed, launched at t = 0 (built from the
    # second shot, which differs from the first only in its launch time and its detector).
    s0, s1 = lvl["shots"]
    doc.check(s0["particle"] == s1["particle"] and all(s0["launch"][k] == s1["launch"][k]
                                                       for k in ("node", "direction", "kinetic_energy")),
              "the two shots differ only in launch time and detector")
    conj = fly(Setup(lvl, ref, 1, t0=0.0, charge_sign=-1.0), H)
    dconj = max(abs(conj["t"] - runs[1]["t"]), abs(conj["x"][0] - runs[1]["x"][0]),
                abs(conj["x"][1] - runs[1]["x"][1]), abs(conj["p"][0] - runs[1]["p"][0]),
                abs(conj["p"][1] - runs[1]["p"][1]))
    doc.check(dconj < 1e-10, "ion 2 equals ion 1 with its charge reversed")
    doc.value("conjDiff", dconj, 1, sci=True)

    # First order: equal and opposite kicks; the measured asymmetry is of second order.
    t1, py1, y1 = first_order(a, H)
    t2, py2, y2 = first_order(b, H)
    doc.check(abs((y1 - a.x0[1]) + (y2 - b.x0[1])) < 1e-12, "the first-order lifts are opposite")
    doc.value("tStar", t1, 6)
    doc.value("liftFirst", y1 - a.x0[1], 6)
    doc.value("pyFirst", py1, 6)
    up, down = runs[0]["x"][1] - a.x0[1], b.x0[1] - runs[1]["x"][1]
    doc.value("liftUp", up, 6)
    doc.value("liftDown", down, 6)
    doc.value("asym", up - down, 3)
    doc.value("asymRel", (up - down) / (0.5 * (up + down)), 2)
