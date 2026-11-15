"""Companion of docs/levels/stern_gerlach/doc.tex (Level 31, Stern-Gerlach).

Every number of the document, from the level file, by the document's two methods:
- the two spin states' flights by the classical Runge-Kutta method at a fixed step
  (Foundations eq:rk4) on eq:eom-moment, the entry into the detector by cubic Hermite
  interpolation and bisection, the error by doubling the step;
- each flight reduced to quadratures in the central potential U = m mu / r^3 of the one
  magnet (eq:potential) by its two conserved quantities, the energy (eq:energy-moment) and
  the angular momentum about the magnet (eq:angular): the turning point by bisection
  (eq:turning), the angle and the time as integrals over r (eq:phi-of-r, eq:t-of-r) by
  mpmath's quadrature, the entry by bisection on r.
Compared with the game: both shots of the reference and of the empty board.
"""
import math

import mpmath


TOL = 1e-9          # ten times the larger of our estimate (~1e-12) and the game's (1e-10)
H = 0.005


class Setup:
    def __init__(self, lvl, placement, shot_index):
        shot = lvl["shots"][shot_index]
        s = lvl["grid"]["subdivision"]
        part = shot["particle"]
        self.q = part["charge"]
        self.mom = part.get("moment", 0.0)
        self.m = part["mass"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.Rm = lvl["physics"]["magnet_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        mags = [e for e in lvl["elements"] if e["kind"] == "magnet"] + \
               [e for e in placement if e["kind"] == "magnet"]
        self.magnets = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in mags]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def U(self, x, y):
        """eq:potential: -m B_z with B_z = -mu/r^3."""
        return sum(self.mom * mu / math.hypot(x - mx, y - my) ** 3 for mx, my, mu in self.magnets)

    def rhs(self, s):
        """eq:eom-moment: F = m grad B_z = 3 m mu r / r^5 for each magnet."""
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        fx = fy = 0.0
        for mx, my, mu in self.magnets:
            rx, ry = x - mx, y - my
            r2 = rx * rx + ry * ry
            k = 3 * self.mom * mu / (r2 * r2 * math.sqrt(r2))
            fx += k * rx
            fy += k * ry
        return (px / (g * self.m), py / (g * self.m), fx, fy)

    def energy(self, s):
        x, y, px, py = s
        return math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2) + self.U(x, y)

    def angular(self, s):
        mx, my, _ = self.magnets[0]
        return (s[0] - mx) * s[3] - (s[1] - my) * s[2]


def rk4_step(st, s, f, h):
    k1 = f
    k2 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)))
    k3 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)))
    k4 = st.rhs(tuple(a + h * b for a, b in zip(s, k3)))
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


def fly(st, h):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    L0 = st.angular(s) if st.magnets else 0.0
    path, wdev, Ldev = [(s[0], s[1])], 0.0, 0.0
    preds = []
    for idx, (mx, my, _) in enumerate(st.magnets):
        preds.append((f"Collided({idx})", 0, lambda p, mx=mx, my=my: math.hypot(p[0] - mx, p[1] - my) <= st.Rm))
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                                  and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, hs)
        f2 = st.rhs(s2)
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
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, hs, hi)
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else (
                face_of(pe, st.wmin, st.wmax) if outcome == "LeftBounds" else "")
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "steps": n + 1, "face": face,
                    "path": path, "wdev": wdev, "Ldev": Ldev, "w0": w0}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if st.magnets:
            Ldev = max(Ldev, abs(st.angular(s) - L0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "steps": n,
                    "face": "", "path": path, "wdev": wdev, "Ldev": Ldev, "w0": w0}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1])) / 15


def quadrature_flight(st):
    """The flight in the central potential of the one magnet, from W and L (section 4.3)."""
    mpmath.mp.dps = 30
    mx, my, mu = (mpmath.mpf(v) for v in st.magnets[0])
    mom, m, c = mpmath.mpf(st.mom), mpmath.mpf(st.m), mpmath.mpf(st.c)
    x0, y0 = mpmath.mpf(st.x0[0]), mpmath.mpf(st.x0[1])
    p0 = mpmath.sqrt(mpmath.mpf(st.T0) ** 2 / c ** 2 + 2 * m * mpmath.mpf(st.T0))
    px0, py0 = p0 * mpmath.mpf(st.d[0]), p0 * mpmath.mpf(st.d[1])
    k = mom * mu
    rx0, ry0 = x0 - mx, y0 - my
    r0 = mpmath.hypot(rx0, ry0)
    W = mpmath.sqrt(m ** 2 * c ** 4 + p0 ** 2 * c ** 2) + k / r0 ** 3       # eq:energy-moment
    L = rx0 * py0 - ry0 * px0                                                 # eq:angular
    p2 = lambda r: ((W - k / r ** 3) / c) ** 2 - (m * c) ** 2
    pr2 = lambda r: p2(r) - (L / r) ** 2
    gam = lambda r: (W - k / r ** 3) / (m * c * c)
    pr0 = (rx0 * px0 + ry0 * py0) / r0
    # eq:turning: the largest r below r0 where p_r = 0, by bisection (p_r^2 > 0 at r0).
    lo, hi = mpmath.mpf("0.01"), r0
    rs = [lo + (hi - lo) * j / 4000 for j in range(4001)]
    vals = [pr2(r) for r in rs]
    j = max(i for i in range(4000) if vals[i] <= 0 < vals[i + 1])
    a, b = rs[j], rs[j + 1]
    for _ in range(200):
        mid = (a + b) / 2
        if pr2(mid) > 0:
            b = mid
        else:
            a = mid
    rt = b

    def integral(f, R):
        return mpmath.quad(lambda u: f(rt + u * u) * 2 * u, [0, mpmath.sqrt(R - rt)])

    # |p_r^2|: at the turning point p_r^2 rounds to +-1e-31; its sign there is rounding.
    dphi = lambda r: L / (r * r * mpmath.sqrt(abs(pr2(r))))       # eq:phi-of-r
    dt = lambda r: gam(r) * m / mpmath.sqrt(abs(pr2(r)))          # eq:t-of-r
    I0, J0 = integral(dphi, r0), integral(dt, r0)
    phi0 = mpmath.atan2(ry0, rx0)
    inward = pr0 < 0

    def point(r):
        phi = phi0 + (I0 if inward else -I0) + integral(dphi, r)
        t = (J0 if inward else -J0) + integral(dt, r)
        prr, pphi = mpmath.sqrt(pr2(r)), L / r
        return (mx + r * mpmath.cos(phi), my + r * mpmath.sin(phi), t,
                prr * mpmath.cos(phi) - pphi * mpmath.sin(phi),
                prr * mpmath.sin(phi) + pphi * mpmath.cos(phi))

    lx, ly, hx, hy = (mpmath.mpf(v) for v in (st.dmin[0], st.dmin[1], st.dmax[0], st.dmax[1]))

    def inside(r):
        x, y, *_ = point(r)
        return lx <= x <= hx and ly <= y <= hy

    r_a, r = rt, rt + mpmath.mpf("0.25")
    r_b = None
    while r < 60:
        if inside(r):
            r_b = r
            break
        r_a, r = r, r + mpmath.mpf("0.25")
    for _ in range(100):
        mid = (r_a + r_b) / 2
        if inside(mid):
            r_b = mid
        else:
            r_a = mid
    x, y, t, pxe, pye = point(r_b)
    return {"W": W, "L": L, "k": k, "rt": rt, "r0": r0, "pr0": pr0, "x": (x, y), "t": t,
            "p": (pxe, pye), "rE": r_b}


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(len(lvl["shots"]) == 2 and len(ref) == 1, "two spin states and one magnet")
    tags = ["Up", "Down"]
    for i, tag in enumerate(tags):
        st = Setup(lvl, ref, i)
        if i == 0:
            mx, my, mu = st.magnets[0]
            doc.value("m", st.m, 3)
            doc.value("H", H, 4)
            doc.value("c", st.c, 3)
            doc.value("Tzero", st.T0, 3)
            doc.value("gammaZero", st.gamma0, 10)
            doc.value("pZero", st.p0, 12)
            doc.value("vZero", st.p0 / (st.gamma0 * st.m), 10)
            doc.value("Rm", st.Rm, 3)
            doc.value("xZero", st.x0[0], 3)
            doc.value("yZero", st.x0[1], 3)
            doc.value("xMag", mx, 3)
            doc.value("yMag", my, 3)
            doc.value("muMag", mu, 3, sci=True)
            doc.value("xWorld", st.wmax[0], 3)
            doc.value("yWorld", st.wmax[1], 3)
            reg = lvl["limits"]["region"]
            for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                          ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
                doc.value(k_, v, 3)
            r0 = math.hypot(st.x0[0] - mx, st.x0[1] - my)
            doc.value("rZero", r0, 6)
            doc.value("bZero", -mu / r0 ** 3, 6)
        doc.value(f"mom{tag}", st.mom, 3, sci=True)
        doc.value(f"k{tag}", st.mom * st.magnets[0][2], 4)
        for k_, v in (("xDetLo", st.dmin[0]), ("yDetLo", st.dmin[1]), ("xDetHi", st.dmax[0]),
                      ("yDetHi", st.dmax[1])):
            doc.value(k_ + tag, v, 3)
        doc.value(f"UZero{tag}", st.U(*st.x0), 8)
        run, half = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", tag + ": the reference flight arrives")
        err = err_of(run, half)
        doc.check(err < 1e-9, tag + ": the step-halving estimate is below 1e-9")
        doc.value(f"t{tag}", run["t"], 12)
        doc.value(f"x{tag}", run["x"][0], 12)
        doc.value(f"y{tag}", run["x"][1], 12)
        doc.value(f"px{tag}", run["p"][0], 10)
        doc.value(f"py{tag}", run["p"][1], 10)
        doc.text(f"face{tag}", run["face"])
        doc.value(f"ang{tag}", math.degrees(math.atan2(run["p"][1], run["p"][0])), 5)
        doc.value(f"steps{tag}", run["steps"])
        doc.value(f"err{tag}", err, 2, sci=True)
        doc.value(f"wDev{tag}", run["wdev"], 2, sci=True)
        doc.value(f"LDev{tag}", run["Ldev"], 2, sci=True)
        doc.data("path" + tag.lower(), run["path"])
        qf = quadrature_flight(st)
        doc.value(f"W{tag}", float(qf["W"]), 12)
        doc.value(f"L{tag}", float(qf["L"]), 12)
        doc.value(f"prZero{tag}", float(qf["pr0"]), 10)
        doc.value(f"rT{tag}", float(qf["rt"]), 10)
        doc.value(f"rE{tag}", float(qf["rE"]), 8)
        dq = max(abs(float(qf["t"]) - run["t"]), abs(float(qf["x"][0]) - run["x"][0]),
                 abs(float(qf["x"][1]) - run["x"][1]), abs(float(qf["p"][0]) - run["p"][0]),
                 abs(float(qf["p"][1]) - run["p"][1]))
        doc.value(f"quadDiff{tag}", dq, 2, sci=True)
        doc.check(dq < 1e-9, tag + ": the quadratures agree with the march")
        doc.value(f"tQ{tag}", float(qf["t"]), 12)
        doc.value(f"yQ{tag}", float(qf["x"][1]), 12)

        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        title = "spin up: " if i == 0 else "spin down: "
        doc.compare(f"Out{tag}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{tag}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{tag}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{tag}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Px{tag}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
        doc.compare(f"Py{tag}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)

        bare = Setup(lvl, [], i)
        rb = fly(bare, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{tag}", ("empty board, " + title) + "outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BareT{tag}", ("empty board, " + title) + r"$t_\mathrm{e}$", rb["t"], gb["end"]["t"], TOL)
        if i == 0:
            doc.value("bareT", rb["t"], 10)
            doc.value("bareX", rb["x"][0], 4)
            doc.value("bareY", rb["x"][1], 4)
