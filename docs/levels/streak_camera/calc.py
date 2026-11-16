"""Companion of docs/levels/streak_camera/doc.tex (Level 47, Streak camera).

Every number of the document, from the level file. The antenna's exact retarded fields (rf_kick
eq:fields) at the level's frequency; each ion's flight from its launch time, the fields at
t0 + t (rf_separator eq:offset); the relativistic equations of motion (Foundations eq:eom)
marched by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) with Lienard's
radiated energy as one more state; the entry by cubic Hermite interpolation and bisection; the
error against the step 2h. To first order in the field (rf_kick eq:first-order) each ion's lift
is the real part of one complex response times exp(i omega t0) (tune_the_rf eq:phasor); launched
a third of a period apart, the three lifts sum to zero (eq:sum), so equal spacing puts the middle
one on the beam line (eq:spacing). Compared with the game: the reference, the empty board and
the reference antenna at half its amplitude, every ion.
"""
import cmath
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-6      # the radiated energy, as in rf_kick: the game's Gauss-Legendre per step, ~1e-7; ten times
H = 0.004
NAMES = ("One", "Two", "Three")


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
    def __init__(self, lvl, placement, shot_i, t0=None):
        shot = lvl["shots"][shot_i]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.t0 = shot["launch"].get("time", 0.0) if t0 is None else t0
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
        els = list(lvl["elements"]) + list(placement)
        assert all(e["kind"] == "antenna" for e in els)
        self.antennas = [Antenna(e, s, e.get("omega") or self.omega, self.c) for e in els]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)
        self.v0 = self.p0 / (self.gamma0 * self.m)

    def field(self, x, y, t):
        """rf_separator eq:offset: the fields at the lab time t0 + t."""
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


def fly(st, h, x_line=None):
    """The march to the first event (Foundations sec:board); with x_line, also the height at the
    crossing of the line x = x_line (the first-order comparison), on the same interpolant."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    path = [(s[0], s[1])]
    y_line = None
    preds = [(f"Collided({i})", 0, lambda p, a=a: math.hypot(p[0] - a.pos[0], p[1] - a.pos[1]) <= st.R)
             for i, a in enumerate(st.antennas)]
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, t, hs)
        f2 = st.rhs(s2, t + hs)
        if x_line is not None and y_line is None and s2[0] >= x_line:
            lo, hi = 0.0, 1.0
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if hermite(s, f, s2, f2, hs, mid)[0] >= x_line:
                    hi = mid
                else:
                    lo = mid
            y_line = hermite(s, f, s2, f2, hs, hi)[1]
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
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "path": path,
                    "rad": se[4], "y_line": y_line}
        s, f, t, n = s2, f2, t + hs, n + 1
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "path": path,
                    "rad": s[4], "y_line": y_line}


def err_of(a, b):
    """Foundations sec:numerics: the h run's error from the 2h run, |y_2h - y_h| / 15."""
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def first_order(st, h):
    """rf_kick eq:first-order with the launch time: along x = x0 + v0 t, y = y0, the transverse
    force q(E_y - v0 B_z) at t0 + t integrated twice to the detectors' near face (RK4, the step
    nearest h that divides the flight time)."""
    v0 = st.v0
    t_star = (st.dmin[0] - st.x0[0]) / v0

    def g(t, s):
        ex, ey, bz = st.field(st.x0[0] + v0 * t, st.x0[1], t)
        return (st.q * (ey - v0 * bz), s[0] / (st.gamma0 * st.m))

    n = round(t_star / h)
    hh = t_star / n
    s, t = (0.0, 0.0), 0.0
    for _ in range(n):
        k1 = g(t, s)
        k2 = g(t + hh / 2, tuple(a + hh / 2 * b for a, b in zip(s, k1)))
        k3 = g(t + hh / 2, tuple(a + hh / 2 * b for a, b in zip(s, k2)))
        k4 = g(t + hh, tuple(a + hh * b for a, b in zip(s, k3)))
        s = tuple(a + hh / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))
        t += hh
    return t_star, s[1]


def response(lvl, placement, h):
    """tune_the_rf eq:phasor: K = L(0) + i L(-pi/2 omega), from two first-order integrals."""
    w = lvl["physics"]["rf_omega"]
    _, yc = first_order(Setup(lvl, placement, 0, t0=0.0), h)
    _, ys = first_order(Setup(lvl, placement, 0, t0=-math.pi / (2 * w)), h)
    return complex(yc, ys)


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    lim = lvl["limits"]
    doc.check(lvl["physics"]["c"] is not None and len(lvl["shots"]) == 3 and len(ref) == 1
              and not lvl["elements"] and not lvl.get("disturbances")
              and not lvl["physics"].get("radiation_reaction")
              and ref[0]["kind"] == "antenna" and "omega" not in ref[0],
              "relativistic, three shots, one antenna at the level's frequency")
    shots = lvl["shots"]
    doc.check(all(sh["particle"] == shots[0]["particle"] and all(sh["launch"][k] == shots[0]["launch"][k]
                                                                for k in ("node", "direction", "kinetic_energy"))
                  for sh in shots),
              "the three shots differ only in launch time and detector")
    sts = [Setup(lvl, ref, i) for i in range(3)]
    a = sts[0]
    (ant,) = a.antennas
    e = ref[0]
    w = a.omega
    period = 2 * math.pi / w
    for k, st in enumerate(sts):
        doc.check(abs(st.t0 - k * period / 3) < 1e-12, f"ion {k + 1} leaves {k}/3 of a period after the first")
    dt = sts[1].t0 - sts[0].t0
    for k_, v, dg in (("Tzero", a.T0, 4), ("m", a.m, 4), ("c", a.c, 4), ("Ra", a.R, 4),
                      ("xZero", a.x0[0], 4), ("yZero", a.x0[1], 4), ("xWorld", a.wmax[0], 4),
                      ("yWorld", a.wmax[1], 4), ("xA", ant.pos[0], 4), ("yA", ant.pos[1], 4),
                      ("angleA", e["angle_deg"], 4), ("omega", w, 4), ("H", H, 4), ("tMax", a.tmax, 4),
                      ("dt", dt, 6), ("dtTwo", sts[2].t0, 6), ("period", period, 6),
                      ("vZero", a.v0, 10), ("gammaZero", a.gamma0, 10), ("pZero", a.p0, 10)):
        doc.value(k_, v, dg)
    doc.value("q", a.q, 4, sci=True)
    doc.value("pA", e["value"], 4, sci=True)
    doc.text("amplitudes", ", ".join(rf"\num{{{x:g}}}" for x in lim["antenna_amplitudes"]))
    for st, nm in zip(sts, NAMES):
        for k_, v in (("xLo", st.dmin[0]), ("yLo", st.dmin[1]), ("xHi", st.dmax[0]), ("yHi", st.dmax[1])):
            doc.value(f"{k_}{nm}", v, 4)
        doc.value(f"aim{nm}", 0.5 * (st.dmin[1] + st.dmax[1]) - st.x0[1], 4)
    for k_, v in (("xRegLo", lim["region"]["min"][0]), ("yRegLo", lim["region"]["min"][1]),
                  ("xRegHi", lim["region"]["max"][0]), ("yRegHi", lim["region"]["max"][1])):
        doc.value(k_, v, 4)
    aims = [0.5 * (st.dmin[1] + st.dmax[1]) - st.x0[1] for st in sts]
    doc.check(abs(aims[1]) < 1e-12 and abs(aims[0] + aims[2]) < 1e-12,
              "the detectors are centred on the beam line and symmetric about it (eq:spacing)")
    half_width = 0.5 * (sts[0].dmax[1] - sts[0].dmin[1])
    doc.value("halfWidth", half_width, 4)

    # The response (tune_the_rf eq:phasor) and the three first-order lifts at the detectors' face.
    t_star, _ = first_order(a, H)
    doc.value("tStar", t_star, 6)
    K = response(lvl, ref, H)
    lifts = [(K * cmath.exp(1j * w * st.t0)).real for st in sts]
    for st, nm, L in zip(sts, NAMES, lifts):
        _, direct = first_order(st, H)
        doc.check(abs(direct - L) < 1e-10, f"the phasor gives ion {nm}'s first-order lift")
        doc.value(f"L{nm}", L, 4)
    # eq:sum: 1 + exp(2 pi i/3) + exp(4 pi i/3) = 0.
    doc.check(abs(sum(lifts)) < 1e-10, "the three first-order lifts sum to zero (eq:sum)")
    doc.value("liftSum", abs(sum(lifts)), 1, sci=True)
    alpha = cmath.phase(K)
    doc.value("Kabs", abs(K), 4)
    doc.value("Karg", math.degrees(alpha), 4)
    doc.check(alpha < 0, "the reference's response lags (the text writes exp(-i |alpha|))")
    doc.value("KargNeg", -math.degrees(alpha), 4)
    doc.value("KargIdeal", -30.0, 4)
    doc.value("Kneed", aims[0] / math.cos(math.pi / 6), 4)
    doc.value("KneedLo", (aims[0] - half_width) / math.cos(math.pi / 6), 4)
    doc.value("KneedHi", (aims[0] + half_width) / math.cos(math.pi / 6), 4)
    # The listed amplitudes at the reference's place, orientation and sign: K scales with p0.
    rows = []
    for amp in lim["antenna_amplitudes"]:
        f = amp / abs(e["value"])
        rows.append((rf"\num{{{math.copysign(amp, e['value']):g}}}", abs(K) * f, lifts[0] * f,
                     lifts[1] * f, lifts[2] * f))
    doc.table("amps", [r"$p_0$", r"$\abs{K}$", r"$L_1$", r"$L_2$", r"$L_3$"], rows, digits=[4, 4, 4, 4, 4])
    half = dict(e, value=e["value"] / 2)
    doc.check(abs(half["value"]) in lim["antenna_amplitudes"], "half the reference's amplitude is on the list")
    doc.value("pHalf", half["value"], 4, sci=True)

    # eq:orders: the marched lift at the face, the amplitude scaled by u, over u; it tends to the
    # first-order lift as u -> 0, linearly in u (the second order), which a Richardson step from
    # the two smallest u removes.
    us = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
    per_u = []
    for u in us:
        row = []
        for i in range(3):
            su = Setup(lvl, [dict(e, value=e["value"] * u)], i)
            r = fly(su, H, su.dmin[0])
            doc.check(r["y_line"] is not None, "the scaled flight crosses the detectors' face")
            row.append((r["y_line"] - su.x0[1]) / u)
        per_u.append(row)
    rows = [(rf"\num{{{u:g}}}", *row) for u, row in zip(us, per_u)]
    extra = [2 * per_u[-1][i] - per_u[-2][i] for i in range(3)]
    rows.append((r"$u \to 0$ (extrapolated)", *extra))
    rows.append((r"first order, \eqref{eq:lifts}", *lifts))
    doc.table("orders", [r"$u$", r"$y_1/u$", r"$y_2/u$", r"$y_3/u$"], rows, digits=[4, 5, 5, 5])
    gap0 = max(abs(extra[i] - lifts[i]) for i in range(3))
    doc.check(gap0 < 0.01 * abs(K), "the marched lifts tend to the first-order ones as the amplitude vanishes")
    doc.value("extraGap", gap0, 1, sci=True)
    for i, nm in enumerate(NAMES):
        doc.value(f"sTwo{nm}", (per_u[-1][i] - lifts[i]) / us[-1], 2)

    largest = 0.0
    for i, (st, nm) in enumerate(zip(sts, NAMES)):
        run, coarse = fly(st, H, st.dmin[0]), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", f"ion {i + 1} arrives")
        err = err_of(run, coarse)
        doc.value(f"t{nm}", run["t"], 10)
        doc.value(f"x{nm}", run["x"][0], 10)
        doc.value(f"y{nm}", run["x"][1], 10)
        doc.value(f"px{nm}", run["p"][0], 8)
        doc.value(f"py{nm}", run["p"][1], 8)
        doc.text(f"face{nm}", run["face"])
        doc.value(f"err{nm}", err, 2, sci=True)
        doc.value(f"rad{nm}", run["rad"], 2, sci=True)
        doc.check(run["y_line"] is not None, f"ion {i + 1} crosses the detectors' face")
        doc.value(f"yLine{nm}", run["y_line"] - st.x0[1], 4)
        doc.value(f"resid{nm}", run["y_line"] - st.x0[1] - lifts[i], 3)
        doc.data("path" + nm.lower(), run["path"])
        if nm == "Two":
            doc.value("marginTwo", run["x"][1] - st.dmin[1], 3)
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
        largest = max(largest, abs(run["t"] - g["end"]["t"]), abs(run["x"][0] - g["end"]["x"][0]),
                      abs(run["x"][1] - g["end"]["x"][1]), abs(run["p"][0] - g["end"]["p"][0]),
                      abs(run["p"][1] - g["end"]["p"][1]))
        # The empty board: a straight line along y = y0, into the detector's near face if the
        # line crosses it, else to the world's right edge.
        gb = doc.flight(1, flight=i)
        on_line = st.dmin[1] <= st.x0[1] <= st.dmax[1]
        bare_out = "Arrived" if on_line else "LeftBounds"
        t_free = ((st.dmin[0] if on_line else st.wmax[0]) - st.x0[0]) / st.v0
        doc.compare(f"BareOut{nm}", f"empty board, ion {i + 1}: outcome", bare_out, gb["outcome"])
        doc.compare(f"BareT{nm}", f"empty board, ion {i + 1}: " + r"$t_\mathrm{e}$", t_free, gb["end"]["t"], TOL)
        doc.check(on_line == (nm == "Two"), "only the middle detector lies on the beam line")
        # Half the amplitude: the outer ions fall short of their detectors.
        sh = Setup(lvl, [half], i)
        hr, hc = fly(sh, H), fly(sh, 2 * H)
        herr = err_of(hr, hc)
        gh, ghp = doc.flight(2, flight=i), doc.flight(2, flight=i, run="preview")
        hgap = max(abs(ghp["end"]["t"] - gh["end"]["t"]), abs(ghp["end"]["x"][1] - gh["end"]["x"][1]))
        htol = max(TOL, 10 * herr, hgap)
        doc.text(f"halfOut{nm}", {"Arrived": "arrives", "LeftBounds": "leaves the world"}.get(hr["outcome"], hr["outcome"]))
        doc.value(f"halfY{nm}", hr["x"][1], 6)
        doc.value(f"halfX{nm}", hr["x"][0], 6)
        title = f"half amplitude, ion {i + 1}: "
        doc.compare(f"HalfOut{nm}", title + "outcome", hr["outcome"], gh["outcome"])
        doc.compare(f"HalfT{nm}", title + r"$t_\mathrm{e}$", hr["t"], gh["end"]["t"], htol)
        doc.compare(f"HalfY{nm}", title + r"$y_\mathrm{e}$", hr["x"][1], gh["end"]["x"][1], htol)
    doc.game_value("largest", largest, 2, sci=True)
    doc.value("bareT", (a.wmax[0] - a.x0[0]) / a.v0, 8)
