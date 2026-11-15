"""Companion of docs/levels/ramp_the_coil/doc.tex (Level 40, Ramp the coil).

Every number of the document, from the level file. Newtonian mechanics (c = infinity). The coil
is a circular loop whose current rises as kappa(t) = rate t: its field B_z is First coil's closed
form (first_coil eq:loop-field, K and E by the arithmetic-geometric mean, first_coil eq:agm)
scaled by kappa(t), and it induces E = -rate A with A the loop's vector potential per unit kappa
(eq:vector-potential; checked against its defining integral by the trapezoidal rule and against
B_z = (1/rho) d(rho A)/d rho by central differences). The flight: Newton's equation marched by the
classical Runge-Kutta method at a fixed step (Foundations eq:rk4) with the field's work as one
more state; the entry by cubic Hermite interpolation and bisection; the error against the step
2h; the canonical angular momentum about the coil's axis (eq:canonical) conserved. The kick to
first order in the rate along the undeflected line (eq:first-order). The solving interval of the
rate: its edges by bisection on the rate. Compared with the game: the reference, the supply at 0
(no supply placed), and the rate just inside and just outside each edge (placements.json).
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.004
H_EDGE = 0.02       # the edges' bisection: step error ~1e-9 in y, far below the edges' sensitivity
TABLE_EVERY = 4.0


def agm_ke(m):
    """K(m), E(m), m = k^2, by the arithmetic-geometric mean (first_coil eq:agm)."""
    a, b = 1.0, math.sqrt(1.0 - m)
    csum, p = 0.5 * m, 1.0
    for _ in range(60):
        c = 0.5 * (a - b)
        a, b = 0.5 * (a + b), math.sqrt(a * b)
        p *= 2.0
        csum += 0.5 * p * c * c
        if abs(c) < 1e-17:
            break
    K = math.pi / (2.0 * a)
    return K, K * (1.0 - csum)


def h_series(m):
    """h(m) = (1 - m/2) K - E = (pi/2) sum_{n>=2} m^n [2n c_n/(2n-1) - c_{n-1}/2],
    c_n = ((2n)!/(4^n n!^2))^2 (eq:h-series), for small m where the direct form cancels."""
    c_prev, total, n = 0.25, 0.0, 2          # c_1 = 1/4
    while True:
        c_n = c_prev * ((2 * n - 1) / (2 * n)) ** 2
        term = m ** n * (2 * n * c_n / (2 * n - 1) - c_prev / 2)
        total += term
        if abs(term) < 1e-18 * abs(total):
            return 0.5 * math.pi * total
        c_prev, n = c_n, n + 1


class Loop:
    """A circular loop of radius a in the plane, per unit kappa (fields at z = 0)."""

    def __init__(self, cx, cy, a):
        self.cx, self.cy, self.a = cx, cy, a

    def parts(self, rho):
        a = self.a
        alpha2, beta = (a - rho) ** 2, a + rho
        m = 4 * a * rho / (beta * beta)
        return alpha2, beta, m

    def bz(self, rho):
        """first_coil eq:loop-field with kappa = 1."""
        alpha2, beta, m = self.parts(rho)
        K, E = agm_ke(m)
        return 2.0 / (alpha2 * beta) * ((self.a ** 2 - rho ** 2) * E + alpha2 * K)

    def a_phi(self, rho):
        """eq:vector-potential with kappa = 1: A_phi = 2 beta h(m) / rho."""
        _, beta, m = self.parts(rho)
        if m < 0.05:
            h = h_series(m)
        else:
            K, E = agm_ke(m)
            h = (1 - m / 2) * K - E
        return 2.0 * beta * h / rho

    def a_phi_integral(self, rho, n):
        """eq:a-integral: A_phi = a int_0^2pi cos t dt / sqrt(rho^2 + a^2 - 2 a rho cos t),
        by the trapezoidal rule with n nodes (spectrally accurate for a periodic integrand)."""
        a = self.a
        s = sum(math.cos(2 * math.pi * k / n) / math.sqrt(rho * rho + a * a - 2 * a * rho * math.cos(2 * math.pi * k / n))
                for k in range(n))
        return a * s * 2 * math.pi / n


class Setup:
    def __init__(self, lvl, placement):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.wire = ph["wire_radius"]
        self.tmax = ph["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        (coil,) = lvl["coils"]
        self.loop = Loop(coil["center"][0] / s, coil["center"][1] / s, coil["radius"] / s)
        self.kappa0 = coil["kappa"]
        sup = [e for e in placement if e["kind"] == "supply"]
        self.rate = sup[0]["value"] if sup else coil.get("rate", 0.0)
        self.p0 = math.sqrt(2 * self.m * self.T0)

    def fields(self, x, y, t):
        """E (induced, -rate A) and B_z at (x, y, t)."""
        dx, dy = x - self.loop.cx, y - self.loop.cy
        rho = math.hypot(dx, dy)
        bz = (self.kappa0 + self.rate * t) * self.loop.bz(rho)
        ephi = -self.rate * self.loop.a_phi(rho)
        return -ephi * dy / rho, ephi * dx / rho, bz

    def rhs(self, s, t):
        """State (x, y, px, py, work): Newton's equation with q(E + v x B), and q E.v."""
        x, y, px, py = s[:4]
        vx, vy = px / self.m, py / self.m
        ex, ey, bz = self.fields(x, y, t)
        return (vx, vy, self.q * (ex + vy * bz), self.q * (ey - vx * bz), self.q * (ex * vx + ey * vy))

    def canonical(self, s, t):
        """eq:canonical: (r x p)_z + q rho A_phi(rho, t) about the coil's axis."""
        dx, dy = s[0] - self.loop.cx, s[1] - self.loop.cy
        rho = math.hypot(dx, dy)
        return dx * s[3] - dy * s[2] + self.q * rho * (self.kappa0 + self.rate * t) * self.loop.a_phi(rho)


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


def fly(st, h, table=False, checks=True):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    j0 = st.canonical(s, 0.0)
    jdev = wdev = 0.0
    path = [(s[0], s[1])]
    rows = [(0.0,) + s[:4]] if table else []
    lp = st.loop
    preds = [("Collided(coil)", 0, lambda p: abs(math.hypot(p[0] - lp.cx, p[1] - lp.cy) - lp.a) <= st.wire),
             ("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])),
             ("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0)]
    us = [k / 64 for k in range(1, 65)]
    every = round(TABLE_EVERY / h)
    kin = lambda s: (s[2] ** 2 + s[3] ** 2) / (2 * st.m)
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
            if checks:
                jdev = max(jdev, abs(st.canonical(se, te) - j0))
                wdev = max(wdev, abs(kin(se) - st.T0 - se[4]))
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "steps": n + 1,
                    "jdev": jdev, "wdev": wdev, "j0": j0, "work": se[4], "path": path, "rows": rows}
        s, f, t, n = s2, f2, t + hs, n + 1
        if checks:
            jdev = max(jdev, abs(st.canonical(s, t) - j0))
            wdev = max(wdev, abs(kin(s) - st.T0 - s[4]))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if table and n % every == 0:
            rows.append((t,) + s[:4])
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "steps": n,
                    "jdev": jdev, "wdev": wdev, "j0": j0, "work": s[4], "path": path, "rows": rows}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def first_order(st, h):
    """eq:first-order per unit rate: along x = x0 + v0 t, y = y0, the transverse force per unit
    rate, q (E_y - v0 B_z)/rate, integrated twice (RK4 on (p_y1, y1)) to the time the line
    reaches the detector's near face. Returns that time and dy/d(rate) there."""
    v0 = st.p0 / st.m
    t_star = (st.dmin[0] - st.x0[0]) / v0
    lp = st.loop

    def g(t, s):
        x, y = st.x0[0] + v0 * t, st.x0[1]
        dx, dy = x - lp.cx, y - lp.cy
        rho = math.hypot(dx, dy)
        ey = -lp.a_phi(rho) * dx / rho                  # E_y per unit rate
        return (st.q * (ey - v0 * t * lp.bz(rho)), s[0] / st.m)

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


def with_rate(lvl, rate):
    (coil,) = lvl["coils"]
    s = lvl["grid"]["subdivision"]
    node = [round(coil["center"][0] + coil["radius"] * s), coil["center"][1], 0]
    return [{"node": node, "kind": "supply", "value": rate}]


def edge(lvl, inside, outside):
    """The rate where the outcome changes, by bisection between a solving and a failing rate."""
    for _ in range(60):
        mid = 0.5 * (inside + outside)
        if fly(Setup(lvl, with_rate(lvl, mid)), H_EDGE, checks=False)["outcome"] == "Arrived":
            inside = mid
        else:
            outside = mid
        if abs(inside - outside) < 1e-6:
            break
    return inside, outside


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    (coil,) = lvl["coils"]
    doc.check(lvl["physics"]["c"] is None and len(ref) == 1 and ref[0]["kind"] == "supply"
              and coil["shape"] == "circle" and coil["tunable"] and not coil.get("rate"),
              "Newtonian, one tunable circular coil, its supply the only placement")
    st = Setup(lvl, ref)
    lp = st.loop
    rates = lvl["limits"]["coil_rates"]
    for k_, v in (("q", st.q), ("m", st.m), ("Tzero", st.T0), ("xZero", st.x0[0]), ("yZero", st.x0[1]),
                  ("xDetLo", st.dmin[0]), ("yDetLo", st.dmin[1]), ("xDetHi", st.dmax[0]),
                  ("yDetHi", st.dmax[1]), ("xWorld", st.wmax[0]), ("yWorld", st.wmax[1]),
                  ("xC", lp.cx), ("yC", lp.cy), ("a", lp.a), ("wire", st.wire), ("rate", st.rate),
                  ("rateMin", min(rates)), ("rateMax", max(rates)), ("H", H), ("tMax", st.tmax),
                  ("xS", ref[0]["node"][0]), ("yS", ref[0]["node"][1])):
        doc.value(k_, v, 6)
    doc.value("pZero", st.p0, 10)
    v0 = st.p0 / st.m
    doc.value("vZero", v0, 10)
    doc.value("rhoZero", math.hypot(st.x0[0] - lp.cx, st.x0[1] - lp.cy), 6)
    doc.value("rhoLine", abs(st.x0[1] - lp.cy), 6)
    doc.check(abs(ref[0]["node"][0] - (coil["center"][0] + coil["radius"])) < 1e-12
              and ref[0]["node"][1] == coil["center"][1], "the supply on the coil's terminal, due east")

    # The vector potential: closed form against its integral, and against B_z.
    worst_int = worst_curl = 0.0
    rows = []
    for rho in (0.5, 2.0, 7.0, 10.0, 16.0):
        closed = lp.a_phi(rho)
        integral = lp.a_phi_integral(rho, 400)
        worst_int = max(worst_int, abs(closed - integral) / abs(closed))
        d = 1e-5
        curl = ((rho + d) * lp.a_phi(rho + d) - (rho - d) * lp.a_phi(rho - d)) / (2 * d) / rho
        worst_curl = max(worst_curl, abs(curl - lp.bz(rho)) / abs(lp.bz(rho)))
        rows.append((rho, closed, lp.bz(rho)))
    doc.check(worst_int < 1e-12, "A_phi: the closed form equals its integral")
    doc.check(worst_curl < 1e-8, "B_z = (1/rho) d(rho A_phi)/d rho")
    doc.value("aIntegral", worst_int, 1, sci=True)
    doc.value("aCurl", worst_curl, 1, sci=True)
    doc.table("atab", [r"$\rho$", r"$A_\phi/\kappa$", r"$B_z/\kappa$"], rows, digits=[3, 10, 10])
    m_sw = 0.05
    rho_sw = next(r for r in (x / 1e4 for x in range(1, 40000)) if lp.parts(r)[2] >= m_sw)
    _, beta_sw, msw = lp.parts(rho_sw)
    Ksw, Esw = agm_ke(msw)
    sdiff = abs(h_series(msw) - ((1 - msw / 2) * Ksw - Esw)) / h_series(msw)
    doc.check(sdiff < 1e-10, "the series and the direct h(m) agree where the series takes over")
    doc.value("seriesDiff", sdiff, 1, sci=True)

    # The march.
    run, coarse = fly(st, H, table=True), fly(st, 2 * H)
    doc.check(run["outcome"] == "Arrived", "the reference arrives")
    err = err_of(run, coarse)
    for k_, v, dg in (("tRef", run["t"], 10), ("xRef", run["x"][0], 10), ("yRef", run["x"][1], 10),
                      ("pxRef", run["p"][0], 10), ("pyRef", run["p"][1], 10), ("jZero", run["j0"], 10),
                      ("workRef", run["work"], 6)):
        doc.value(k_, v, dg)
    doc.text("faceRef", run["face"])
    doc.value("errRef", err, 2, sci=True)
    doc.value("jDev", run["jdev"], 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    doc.check(run["jdev"] < 1e-10 and run["wdev"] < 1e-10, "J conserved and the work balanced")
    doc.value("kappaEnd", st.rate * run["t"], 6)
    rho_e = math.hypot(run["x"][0] - lp.cx, run["x"][1] - lp.cy)
    doc.value("rhoEnd", rho_e, 6)
    doc.value("speedEnd", math.hypot(*run["p"]) / st.m, 8)
    # the closest approach to the wire, on the march's recorded points (every 20 steps)
    doc.value("wireGap", min(abs(math.hypot(x - lp.cx, y - lp.cy) - lp.a) for x, y in run["path"]), 2)
    doc.data("pathref", run["path"])
    doc.table("march", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$"],
              run["rows"] + [(run["t"], run["x"][0], run["x"][1], run["p"][0], run["p"][1])],
              digits=[4, 8, 8, 8, 8])
    # The canonical momentum gives |p| from the position alone? No: it fixes one combination;
    # report J at the entry from its parts.
    dx, dy = run["x"][0] - lp.cx, run["x"][1] - lp.cy
    doc.value("jMech", dx * run["p"][1] - dy * run["p"][0], 8)
    doc.value("jField", st.q * rho_e * st.rate * run["t"] * lp.a_phi(rho_e), 8)

    # First order in the rate.
    t_star, dydr = first_order(st, H)
    _, dydr2 = first_order(st, 2 * H)
    doc.check(abs(dydr - dydr2) < 1e-14, "the first-order integral converges")
    doc.value("tStar", t_star, 6)
    doc.value("dydr", dydr, 6, sci=True)
    target = 0.5 * (st.dmin[1] + st.dmax[1])
    doc.value("yMid", target, 4)
    doc.value("rateFirst", (target - st.x0[1]) / dydr, 4)
    y_first = st.x0[1] + st.rate * dydr
    doc.value("yFirst", y_first, 6)
    doc.value("firstDiff", y_first - run["x"][1], 3)

    # The solving interval: from the reference outward to a failing rate, then bisection.
    def bracket(step):
        r = st.rate
        while fly(Setup(lvl, with_rate(lvl, r)), H_EDGE, checks=False)["outcome"] == "Arrived":
            r += step
        return r
    lo_in, lo_out = edge(lvl, st.rate, bracket(-250.0))
    hi_in, hi_out = edge(lvl, st.rate, bracket(250.0))
    doc.value("edgeLo", lo_in, 7)
    doc.value("edgeHi", hi_in, 7)
    doc.value("width", (hi_in - lo_in), 4)
    doc.value("widthFrac", 100 * (hi_in - lo_in) / (max(rates) - min(rates)), 3)
    for tag, r in (("Lo", lo_in), ("Hi", hi_in)):
        e = fly(Setup(lvl, with_rate(lvl, r)), H, checks=False)
        doc.text(f"face{tag}", e["face"])
        doc.value(f"xEdge{tag}", e["x"][0], 6)
        doc.value(f"yEdge{tag}", e["x"][1], 6)
        doc.data("path" + tag.lower(), e["path"])

    # Compared with the game: the reference, then placements.json.
    g, gp = doc.flight(0), doc.flight(0, run="preview")
    gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][1] - g["end"]["x"][1]),
              abs(gp["end"]["p"][0] - g["end"]["p"][0]), abs(gp["end"]["p"][1] - g["end"]["p"][1]))
    tol = max(TOL, 10 * err, gap)
    doc.compare("Out", "outcome", run["outcome"], g["outcome"])
    doc.compare("T", r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
    doc.compare("X", r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
    doc.compare("Y", r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
    largest = max(abs(run["t"] - g["end"]["t"]), abs(run["x"][1] - g["end"]["x"][1]),
                  abs(run["p"][0] - g["end"]["p"][0]), abs(run["p"][1] - g["end"]["p"][1]))

    names = {"Off": "no supply"}
    for k in range(1, doc.n_results()):
        placement = doc.result(k)["placement"]
        sk = Setup(lvl, placement)
        tag = "Off" if not placement else f"R{k}"
        rk, ck = fly(sk, H, checks=False), fly(sk, 2 * H, checks=False)
        ek = err_of(rk, ck)
        gk, gkp = doc.flight(k), doc.flight(k, run="preview")
        gapk = max(abs(gkp["end"]["t"] - gk["end"]["t"]), abs(gkp["end"]["x"][0] - gk["end"]["x"][0]),
                   abs(gkp["end"]["x"][1] - gk["end"]["x"][1]))
        tolk = max(TOL, 10 * ek, gapk)
        title = (names.get(tag) or f"rate ${sk.rate:.6g}$") + ": "
        doc.compare(f"{tag}Out", title + "outcome", rk["outcome"], gk["outcome"])
        doc.compare(f"{tag}X", title + r"$x_\mathrm{e}$", rk["x"][0], gk["end"]["x"][0], tolk)
        doc.compare(f"{tag}Y", title + r"$y_\mathrm{e}$", rk["x"][1], gk["end"]["x"][1], tolk)
        largest = max(largest, abs(rk["x"][0] - gk["end"]["x"][0]), abs(rk["x"][1] - gk["end"]["x"][1]))
        if tag == "Off":
            doc.value("offT", rk["t"], 8)
            doc.text("offOut", rk["outcome"])
    doc.game_value("largest", largest, 2, sci=True)
