"""Companion of docs/levels/mass_spectrometer/doc.tex (Level 37, Mass spectrometer).

Every number of the document, from the level file: the six ions' flights (three masses, two
launch angles) by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) in the
placed charges' Coulomb fields and the placed magnets' in-plane dipole fields
(first_magnet eq:dipole-plane), the gate passed in order before the detector (PHYSICS.md 6.2),
the entries by cubic Hermite interpolation and bisection, the errors by halving the step, and
each exit angle checked by the curvature relation of the magnetic part (calutron
eq:curvature does not hold here: the charges also act; the energy is the check). Compared with
the game: the six reference flights and the empty board's.
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
        self.gates = [([g["min"][0] / s, g["min"][1] / s], [g["max"][0] / s, g["max"][1] / s])
                      for g in lvl.get("gates") or []]
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
    passed = 0                      # the gates passed, in order
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, hs)
        f2 = st.rhs(s2)
        pts = [hermite(s, f, s2, f2, hs, u) for u in us]
        if passed < len(st.gates):
            lo_g, hi_g = st.gates[passed]
            if any(box_distance(p, lo_g, hi_g) <= 0 for p in pts):
                passed += 1
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
            if outcome == "Arrived" and passed < len(st.gates):
                outcome = f"SkippedGate({passed})"
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
    n = len(lvl["shots"])
    st0 = Setup(lvl, ref, 0)
    doc.check(n == 6 and len(st0.gates) == 1, "six ions and one gate")
    for k_, v in (("q", st0.q), ("c", st0.c), ("Tzero", st0.T0), ("R", st0.R), ("Rm", st0.Rm),
                  ("xZero", st0.x0[0]), ("yZero", st0.x0[1]), ("xWorld", st0.wmax[0]),
                  ("yWorld", st0.wmax[1]), ("H", H), ("tmax", st0.tmax)):
        doc.value(k_, v, 4, sci=(k_ == "q"))
    doc.value("angDeg", abs(math.degrees(math.atan2(st0.d[1], st0.d[0]))), 4)
    (g0, g1), = st0.gates
    for k_, v in (("gx0", g0[0]), ("gy0", g0[1]), ("gx1", g1[0]), ("gy1", g1[1])):
        doc.value(k_, v, 3)
    reg = lvl["limits"]["region"]
    for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                  ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
        doc.value(k_, v, 3)
    els = [["charge", f"({cx:g}, {cy:g})", Q] for cx, cy, Q in st0.charges] + \
          [["magnet", f"({mx:g}, {my:g})", mu] for mx, my, mu in st0.magnets]
    doc.table("elements", ["kind", "at", r"$Q$ or $\mu$"], els, align="llr", digits=[None, None, 3],
              sci=[False, False, True])
    rows = []
    worst = worst_w = 0.0
    names = "abcdef"
    for i in range(n):
        st = Setup(lvl, ref, i)
        # The finer run is reported: the step-halving estimate is its error (the coarser
        # one, 16 times larger, is only the estimate's partner).
        half, run = fly(st, H), fly(st, H / 2)
        doc.check(run["outcome"] == "Arrived", f"ion {i + 1}: the reference flight arrives")
        err = err_of(run, half)
        worst, worst_w = max(worst, err), max(worst_w, run["wdev"])
        ang = math.degrees(math.atan2(st.d[1], st.d[0]))
        rows.append([f"{i + 1}", f"{st.m:g}", f"{ang:+.0f}", run["x"][1], run["t"],
                     f"[{st.dmin[1]:g}, {st.dmax[1]:g}]"])
        doc.data("path" + names[i], run["path"])
        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        if i == 1:
            doc.game_value("gapTwo", gap, 1, sci=True)
        title = f"ion {i + 1}: "
        doc.compare(f"Out{names[i]}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{names[i]}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{names[i]}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{names[i]}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        bare = Setup(lvl, [], i)
        rb = fly(bare, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{names[i]}", f"empty board, ion {i + 1}: outcome", rb["outcome"], gb["outcome"])
        if i == 4:
            doc.text("bareOutFive", rb["outcome"])
    doc.table("results", ["ion", r"$m$", r"angle ($^\circ$)", r"$y_\mathrm{e}$", r"$t_\mathrm{e}$",
                          "detector's $y$"], rows, align="rrrrrl", digits=[None, None, None, 8, 8, None])
    doc.value("xDet", st0.dmin[0], 3)
    doc.value("errMax", worst, 2, sci=True)
    doc.value("wDevMax", worst_w, 2, sci=True)
    doc.check(worst < 1e-7, "the step-halving estimates are below 1e-7")
