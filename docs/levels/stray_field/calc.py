"""Companion of docs/levels/stray_field/doc.tex (Level 38, Stray field).

Every number of the document, from the level file. Newtonian mechanics (c = infinity). Under
the first disturbance (no outside field) the flight about the one placed charge is Level 1's
closed form in its Newtonian limit (first_bend eq:general with Gamma = 1): the entry by stepping
the angle and bisecting, the time by Simpson's rule on first_bend eq:dtdth. Under each
disturbance the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) in the
charge's field plus the uniform field, the entry by cubic Hermite interpolation and bisection,
the error against the step 2h, the energy (eq:energy-stray) the check. Compared with the game:
the reference under both disturbances, the empty board, and two other one-charge solutions
(placements.json).
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.004
N_SIMPSON = 4000


class Setup:
    def __init__(self, lvl, placement, dist):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        els = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in els]
        self.E = lvl["disturbances"][dist]["e"]
        self.p0 = math.sqrt(2 * self.m * self.T0)             # Newtonian: p = m v

    def rhs(self, s):
        x, y, px, py = s
        fx, fy = self.q * self.E[0], self.q * self.E[1]
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += self.q * Q * dx / r3
            fy += self.q * Q * dy / r3
        return (px / self.m, py / self.m, fx, fy)

    def energy(self, s):
        """eq:energy-stray: p^2/2m + V - q E.x."""
        x, y, px, py = s
        V = sum(self.q * Q / math.hypot(x - cx, y - cy) for cx, cy, Q in self.charges)
        return (px * px + py * py) / (2 * self.m) + V - self.q * (self.E[0] * x + self.E[1] * y)


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
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    preds = [(f"Collided({i})", 0, lambda p, cx=cx, cy=cy: math.hypot(p[0] - cx, p[1] - cy) <= st.R)
             for i, (cx, cy, _) in enumerate(st.charges)]
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
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else (
                face_of(pe, st.wmin, st.wmax) if outcome == "LeftBounds" else "")
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "face": face, "steps": n + 1,
                    "wdev": wdev, "path": path}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "face": "",
                    "steps": n, "wdev": wdev, "path": path}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1])) / 15


def latex_q(v):
    """A placed charge for the table's labels in LaTeX: Q = -10^{6}, Q = 2.5 x 10^{5}."""
    mant, exp = f"{v:.1e}".split("e")
    mant = mant.rstrip("0").rstrip(".")
    lead = {"1": "", "-1": "-"}.get(mant, mant + r" \times ")
    return f"Q = {lead}10^{{{int(exp)}}}"


class Kepler:
    """first_bend eq:general in Newton's limit (Gamma = 1, W/c^2 -> m), about the one charge."""

    def __init__(self, st):
        (cx, cy, Q), = st.charges
        self.cx, self.cy, self.m = cx, cy, st.m
        self.kappa = st.q * Q
        rx, ry = st.x0[0] - cx, st.x0[1] - cy
        px, py = st.p0 * st.d[0], st.p0 * st.d[1]
        self.r0 = math.hypot(rx, ry)
        self.phi0 = math.atan2(ry, rx)
        Lz = rx * py - ry * px
        self.sigma = 1.0 if Lz > 0 else -1.0
        self.L = abs(Lz)
        self.B = -self.kappa * st.m / self.L ** 2
        self.u0 = 1 / self.r0
        self.du0 = -((rx * px + ry * py) / self.r0) / self.L

    def u(self, th):
        return self.B + (self.u0 - self.B) * math.cos(th) + self.du0 * math.sin(th)

    def du(self, th):
        return -(self.u0 - self.B) * math.sin(th) + self.du0 * math.cos(th)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.cx + r * math.cos(phi), self.cy + r * math.sin(phi))

    def time(self, th, n):
        f = lambda a: self.m / (self.L * self.u(a) ** 2)          # first_bend eq:dtdth, c -> inf
        h = th / n
        s = f(0) + f(th) + 4 * sum(f((2 * k - 1) * h) for k in range(1, n // 2 + 1)) \
            + 2 * sum(f(2 * k * h) for k in range(1, n // 2))
        return s * h / 3

    def momentum(self, th):
        phi = self.phi0 + self.sigma * th
        pr, pphi = -self.L * self.du(th), self.sigma * self.L * self.u(th)
        return (pr * math.cos(phi) - pphi * math.sin(phi), pr * math.sin(phi) + pphi * math.cos(phi))


def kepler_entry(st, kp):
    th, step = 0.0, 1e-3
    while th < 10:
        nxt = th + step
        if box_distance(kp.pos(nxt), st.dmin, st.dmax) <= 0:
            lo, hi = th, nxt
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if box_distance(kp.pos(mid), st.dmin, st.dmax) <= 0:
                    hi = mid
                else:
                    lo = mid
            return hi
        th = nxt
    return None


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(lvl["physics"]["c"] is None and len(lvl["disturbances"]) == 2 and len(ref) == 1,
              "Newtonian, two disturbances, one placed charge")
    st0 = Setup(lvl, ref, 1)
    (cx, cy, Q), = st0.charges
    for k_, v in (("q", st0.q), ("m", st0.m), ("Tzero", st0.T0), ("R", st0.R), ("xZero", st0.x0[0]),
                  ("yZero", st0.x0[1]), ("xDetLo", st0.dmin[0]), ("yDetLo", st0.dmin[1]),
                  ("xDetHi", st0.dmax[0]), ("yDetHi", st0.dmax[1]), ("xWorld", st0.wmax[0]),
                  ("yWorld", st0.wmax[1]), ("Q", Q), ("xC", cx), ("yC", cy), ("H", H)):
        doc.value(k_, v, 4, sci=k_ in ("q", "Q"))
    doc.value("Ey", st0.E[1], 4)
    doc.value("pZero", st0.p0, 10)
    doc.value("vZero", st0.p0 / st0.m, 10)
    doc.value("accel", st0.q * st0.E[1] / st0.m, 4, sci=True)
    tcross = (st0.dmin[0] - st0.x0[0]) / (st0.p0 / st0.m)
    doc.value("tCross", tcross, 4)
    doc.value("dropEst", 0.5 * abs(st0.q * st0.E[1] / st0.m) * tcross ** 2, 3)
    doc.text("distOne", lvl["disturbances"][0]["name"])
    doc.text("distTwo", lvl["disturbances"][1]["name"])

    # The first disturbance in closed form (no outside field).
    kp = Kepler(Setup(lvl, ref, 0))
    doc.value("kappa", kp.kappa, 4)
    doc.value("L", kp.L, 10)
    doc.value("Bk", kp.B, 10)
    th = kepler_entry(Setup(lvl, ref, 0), kp)
    doc.check(th is not None, "the closed form reaches the detector")
    pe, tk = kp.pos(th), kp.time(th, N_SIMPSON)
    tk2 = kp.time(th, N_SIMPSON // 2)
    pk = kp.momentum(th)
    doc.value("thK", math.degrees(th), 6)
    doc.value("tK", tk, 12)
    doc.value("yK", pe[1], 12)
    doc.value("tKerr", abs(tk - tk2) / 15, 2, sci=True)

    tags = ["Off", "On"]
    largest = 0.0
    for i, tag in enumerate(tags):
        st = Setup(lvl, ref, i)
        run, coarse = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", tag + ": the reference flight arrives")
        err = err_of(run, coarse)
        doc.value(f"t{tag}", run["t"], 10)
        doc.value(f"x{tag}", run["x"][0], 10)
        doc.value(f"y{tag}", run["x"][1], 10)
        doc.text(f"face{tag}", run["face"])
        doc.value(f"err{tag}", err, 2, sci=True)
        doc.value(f"wDev{tag}", run["wdev"], 2, sci=True)
        doc.data("path" + tag.lower(), run["path"])
        if i == 0:
            dk = max(abs(run["t"] - tk), abs(run["x"][0] - pe[0]), abs(run["x"][1] - pe[1]),
                     abs(run["p"][0] - pk[0]), abs(run["p"][1] - pk[1]))
            doc.value("kDiff", dk, 2, sci=True)
            doc.check(dk < 1e-8, "the march agrees with the closed form")
        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"{lvl['disturbances'][i]['name']}: "
        doc.compare(f"Out{tag}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{tag}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{tag}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{tag}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Px{tag}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
        doc.compare(f"Py{tag}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        largest = max(largest, abs(run["t"] - g["end"]["t"]), abs(run["x"][1] - g["end"]["x"][1]),
                      abs(run["p"][0] - g["end"]["p"][0]), abs(run["p"][1] - g["end"]["p"][1]))
        rb = fly(Setup(lvl, [], i), H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{tag}", "empty board, " + title + "outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BareY{tag}", "empty board, " + title + r"$y_\mathrm{e}$", rb["x"][1], gb["end"]["x"][1], TOL)
        doc.value(f"bareY{tag}", rb["x"][1], 6)
    doc.game_value("largest", largest, 2, sci=True)

    # Two other one-charge solutions (placements.json's second and third): a negative charge
    # above the line pulls the ion up, a small positive one just below the launch pushes it.
    # Same march; tolerance as for the reference (ten times our step estimate, 1e-9 at least,
    # or the game's preview-verification gap).
    for k, alt in ((2, "Pull"), (3, "Push")):
        placement = doc.result(k)["placement"]
        (e,) = placement
        doc.value(f"x{alt}", e["node"][0], 4)
        doc.value(f"y{alt}", e["node"][1], 4)
        doc.value(f"Q{alt}", e["value"], 4, sci=True)
        for i, tag in enumerate(tags):
            st = Setup(lvl, placement, i)
            run, coarse = fly(st, H), fly(st, 2 * H)
            doc.check(run["outcome"] == "Arrived", f"{alt}, {tag}: the flight arrives")
            err = err_of(run, coarse)
            doc.value(f"y{alt}{tag}", run["x"][1], 6)
            g = doc.flight(k, flight=i)
            gp = doc.flight(k, flight=i, run="preview")
            gap = abs(gp["end"]["x"][1] - g["end"]["x"][1])
            tol = max(TOL, 10 * err, gap)
            title = f"${latex_q(e['value'])}$ at ({e['node'][0]}, {e['node'][1]}), {lvl['disturbances'][i]['name']}: "
            doc.compare(f"{alt}Out{tag}", title + "outcome", run["outcome"], g["outcome"])
            doc.compare(f"{alt}Y{tag}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
