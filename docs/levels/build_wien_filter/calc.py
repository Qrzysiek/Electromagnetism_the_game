"""Companion of docs/levels/build_wien_filter/doc.tex (Level 36, Build a Wien filter).

Every number of the document, from the level file: the three ions' flights by the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4) in the placed charges' Coulomb field and
the placed magnets' in-plane dipole fields (first_magnet eq:dipole-plane), the entries by cubic
Hermite interpolation and bisection, the middle ion's arrival direction against the detector's
cone, the errors by doubling the step; the balance E_y/B_z along the axis (wien_filter
eq:balance). Compared with the game: the three reference flights and the empty board's.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.004


class Setup:
    def __init__(self, lvl, placement, i):
        shot = lvl["shots"][i]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.Rm = lvl["physics"]["magnet_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        det = shot["detector"]
        self.dmin = [v / s for v in det["min"][:2]]
        self.dmax = [v / s for v in det["max"][:2]]
        acc = det.get("acceptance") or {}
        self.cone = acc.get("direction")            # [axis in degrees, half-angle in degrees]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        els = list(lvl["elements"]) + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in els if e["kind"] == "charge"]
        self.magnets = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in els if e["kind"] == "magnet"]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def bz(self, x, y):
        """first_magnet eq:dipole-plane, summed."""
        return sum(-mu / math.hypot(x - mx, y - my) ** 3 for mx, my, mu in self.magnets)

    def efield(self, x, y):
        ex = ey = 0.0
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            ex += Q * dx / r3
            ey += Q * dy / r3
        return ex, ey

    def rhs(self, s):
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        vx, vy = px / (g * self.m), py / (g * self.m)
        ex, ey = self.efield(x, y)
        b = self.bz(x, y)
        return (vx, vy, self.q * (ex + vy * b), self.q * (ey - vx * b))

    def energy(self, s):
        x, y, px, py = s
        V = sum(self.q * Q / math.hypot(x - cx, y - cy) for cx, cy, Q in self.charges)
        return math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2) + V


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


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


def fly(st, h):
    """From the launch to the first event: a charge, a magnet, the world's edge, the detector
    (obstacles first; charges, then magnets, in the game's order); at the detector the cone
    decides between Arrived and Rejected (PHYSICS.md 6.1)."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    nq = len(st.charges)
    preds = [(f"Collided({i})", 0, lambda p, cx=cx, cy=cy: math.hypot(p[0] - cx, p[1] - cy) <= st.R)
             for i, (cx, cy, _) in enumerate(st.charges)]
    for j, (mx, my, _) in enumerate(st.magnets):
        preds.append((f"Collided({nq + j})", 0, lambda p, mx=mx, my=my: math.hypot(p[0] - mx, p[1] - my) <= st.Rm))
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])))
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
            ang = math.degrees(math.atan2(se[3], se[2]))
            face = ""
            if outcome == "Arrived":
                face = face_of(pe, st.dmin, st.dmax)
                if st.cone is not None and abs(ang - st.cone[0]) > st.cone[1]:
                    outcome = "Rejected"
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "face": face, "ang": ang,
                    "steps": n + 1, "wdev": wdev, "path": path}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "face": "",
                    "ang": 0.0, "steps": n, "wdev": wdev, "path": path}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1])) / 15


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    st0 = Setup(lvl, ref, 1)
    doc.check(len(lvl["shots"]) == 3 and len(st0.charges) == 2 and len(st0.magnets) == 1 and not lvl["coils"],
              "three ions, two placed charges and one placed magnet, no coil")
    doc.value("q", st0.q, 3, sci=True)
    doc.value("m", st0.m, 3)
    doc.value("c", st0.c, 3)
    doc.value("R", st0.R, 3)
    doc.value("Rm", st0.Rm, 3)
    doc.value("xZero", st0.x0[0], 3)
    doc.value("yZero", st0.x0[1], 3)
    doc.value("xWorld", st0.wmax[0], 3)
    doc.value("yWorld", st0.wmax[1], 3)
    doc.value("coneHalf", st0.cone[1], 3)
    doc.value("H", H, 4)
    reg = lvl["limits"]["region"]
    for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                  ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
        doc.value(k_, v, 3)
    for k, (cx, cy, Q) in enumerate(st0.charges):
        tag = ["One", "Two"][k]
        doc.value(f"Q{tag}", Q, 3, sci=True)
        doc.value(f"x{tag}", cx, 3)
        doc.value(f"y{tag}", cy, 3)
    (mx, my, mu), = st0.magnets
    doc.value("xMag", mx, 3)
    doc.value("yMag", my, 3)
    doc.value("muMag", mu, 3, sci=True)

    # The balance along the axis: E_y / B_z against the middle ion's speed (wien_filter eq:balance).
    rows = []
    for x in (6.0, 9.0, 11.0, 13.0, 15.0, 18.0, 22.0):
        y = st0.x0[1]
        b = st0.bz(x, y)
        ex, ey = st0.efield(x, y)
        rows.append([f"{x:g}", b, ey, ey / b])
    doc.table("axis", [r"$x$ (on $y = y_0$)", r"$B_z$", r"$E_y$", r"$E_y/B_z$"], rows,
              digits=[None, 6, 6, 4])

    # ---- The three ions ---------------------------------------------------------------------
    tags = ["Slow", "Mid", "Fast"]
    worst_err = 0.0
    for i, tag in enumerate(tags):
        st = Setup(lvl, ref, i)
        v0 = st.p0 / (st.gamma0 * st.m)
        doc.value(f"T{tag}", st.T0, 3)
        doc.value(f"v{tag}", v0, 6)
        doc.value(f"p{tag}", st.p0, 8)
        for k_, v in (("xDetLo", st.dmin[0]), ("yDetLo", st.dmin[1]), ("xDetHi", st.dmax[0]),
                      ("yDetHi", st.dmax[1])):
            doc.value(k_ + tag, v, 3)
        run, half = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", tag + ": the reference flight arrives")
        err = err_of(run, half)
        worst_err = max(worst_err, err)
        doc.value(f"t{tag}", run["t"], 10)
        doc.value(f"x{tag}", run["x"][0], 10)
        doc.value(f"y{tag}", run["x"][1], 10)
        doc.value(f"ang{tag}", run["ang"], 5)
        doc.text(f"face{tag}", run["face"])
        doc.value(f"wDev{tag}", run["wdev"], 2, sci=True)
        doc.data("path" + tag.lower(), run["path"])
        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"{['slow', 'middle', 'fast'][i]} ion: "
        doc.compare(f"Out{tag}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"Tc{tag}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"Xc{tag}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Yc{tag}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Pyc{tag}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        bare = Setup(lvl, [], i)
        rb = fly(bare, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{tag}", "empty board, " + title + "outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BareX{tag}", "empty board, " + title + r"$x_\mathrm{e}$", rb["x"][0], gb["end"]["x"][0], TOL)
        if i == 1:
            doc.value("bareXMid", rb["x"][0], 5)
            doc.value("bareYMid", rb["x"][1], 4)
            doc.text("bareOutMid", rb["outcome"])
    doc.value("errMax", worst_err, 2, sci=True)
    doc.check(worst_err < 1e-9, "the step-halving estimates are below 1e-9")
