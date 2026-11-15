"""Companion of docs/levels/crt_earth_field/doc.tex (Level 44, CRT in the Earth's field).

Every number of the document, from the level file. Relativistic mechanics (c finite) in the
charges' Coulomb fields and the uniform stray B_z: Foundations eq:eom marched by the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4), with Lienard's radiated energy as one
more state; the entry by cubic Hermite interpolation and bisection; the error against the step
2h; the energy gamma m c^2 + q phi conserved as the check (eq:energy-crt). The field's effect and
the correction: each beam's height on the line x = x_spot without the field, with it, and with
the added charge (marched), against the first-order sag of a circle (eq:sag-crt). Compared with
the game: the reference and the tube as built, both shots, in the field.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-6      # Lienard's energy: the game's three-point Gauss-Legendre per step, sixth order,
                    # ~1e-7 relative at its steps (as in rf_kick); ten times
H = 0.004


class Setup:
    def __init__(self, lvl, placement, shot_i, bz=None):
        shot = lvl["shots"][shot_i]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.c = lvl["physics"]["c"]
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
        if bz is None:
            (dz,) = lvl["disturbances"]
            assert dz["e"] == [0.0, 0.0] and not dz.get("waves")
            bz = dz["bz"]
        self.bz = bz
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)                 # Foundations eq:launch
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)

    def rhs(self, s):
        x, y, px, py = s[:4]
        m, c, q = self.m, self.c, self.q
        g = math.sqrt(1 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        ex = ey = 0.0
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            ex += Q * dx / r3
            ey += Q * dy / r3
        fx, fy = q * (ex + vy * self.bz), q * (ey - vx * self.bz)
        vF = vx * fx + vy * fy
        ax, ay = (fx - vx * vF / c ** 2) / (g * m), (fy - vy * vF / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        P = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        return (vx, vy, fx, fy, P)

    def energy(self, s):
        """eq:energy-crt: gamma m c^2 + sum q Q / r."""
        x, y, px, py = s[:4]
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.c) ** 2)
        return g * self.m * self.c ** 2 + sum(self.q * Q / math.hypot(x - cx, y - cy) for cx, cy, Q in self.charges)


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


def fly(st, h, x_stop=None):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0)
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    if x_stop is not None:
        preds = [("Crossed", 0, lambda p: p[0] >= x_stop),
                 ("LeftBounds", 1, lambda p: not (st.wmin[1] <= p[1] <= st.wmax[1]))]
    else:
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
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "path": path,
                    "rad": se[4], "wdev": max(wdev, abs(st.energy(se) - w0))}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "path": path,
                    "rad": s[4], "wdev": wdev}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    (dz,) = lvl["disturbances"]
    doc.check(lvl["physics"]["c"] is not None and len(lvl["shots"]) == 2 and len(ref) == 1
              and len(lvl["elements"]) == 1, "relativistic, two shots, one built-in and one placed charge")
    sts = [Setup(lvl, ref, i) for i in range(2)]
    st0 = sts[0]
    (cbx, cby, Qb), (crx, cry, Qr) = st0.charges
    for k_, v, dg in (("m", st0.m, 4), ("c", st0.c, 4), ("xZero", st0.x0[0], 4), ("yZero", st0.x0[1], 4),
                      ("xWorld", st0.wmax[0], 4), ("yWorld", st0.wmax[1], 4), ("R", st0.R, 4),
                      ("H", H, 4), ("B", st0.bz, 4), ("xB", cbx, 4), ("yB", cby, 4), ("xQ", crx, 4), ("yQ", cry, 4),
                      ("tMax", st0.tmax, 4)):
        doc.value(k_, v, dg)
    doc.value("q", st0.q, 4, sci=True)
    doc.value("QB", Qb, 4, sci=True)
    doc.value("QR", Qr, 4, sci=True)
    region = lvl["limits"]["region"]
    for k_, v in (("xRegLo", region["min"][0]), ("yRegLo", region["min"][1]),
                  ("xRegHi", region["max"][0]), ("yRegHi", region["max"][1])):
        doc.value(k_, v, 4)
    doc.text("distName", dz["name"])

    names = ["Slow", "Fast"]
    largest = 0.0
    lines = []
    for i, st in enumerate(sts):
        nm = names[i]
        for k_, v, dg in (("T", st.T0, 4), ("gamma", st.gamma0, 8), ("p", st.p0, 8),
                          ("v", st.p0 / (st.gamma0 * st.m), 8), ("xS", st.dmin[0], 4), ("ySLo", st.dmin[1], 4),
                          ("xSHi", st.dmax[0], 4), ("ySHi", st.dmax[1], 4)):
            doc.value(f"{k_}{nm}", v, dg)
        Rb = st.p0 / abs(st.q * st.bz)
        L = st.dmin[0] - st.x0[0]
        doc.value(f"RB{nm}", Rb, 5)
        doc.value(f"sag{nm}", L * L / (2 * Rb), 4)
        # The reference's flight, its error and checks, against the game.
        run, coarse = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", f"the {nm.lower()} beam arrives")
        err = err_of(run, coarse)
        doc.value(f"t{nm}", run["t"], 10)
        doc.value(f"x{nm}", run["x"][0], 10)
        doc.value(f"y{nm}", run["x"][1], 10)
        doc.value(f"px{nm}", run["p"][0], 8)
        doc.value(f"py{nm}", run["p"][1], 8)
        doc.text(f"face{nm}", run["face"])
        doc.value(f"err{nm}", err, 2, sci=True)
        doc.value(f"wDev{nm}", run["wdev"], 2, sci=True)
        doc.value(f"rad{nm}", run["rad"], 2, sci=True)
        doc.value(f"radFrac{nm}", run["rad"] / st.T0, 2, sci=True)
        doc.check(run["wdev"] < 1e-11 and run["rad"] < 1e-10 * st.T0, "energy conserved, radiation negligible")
        doc.data("path" + nm.lower(), run["path"])
        g, gp = doc.flight(0, flight=i), doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]), abs(gp["end"]["p"][0] - g["end"]["p"][0]),
                  abs(gp["end"]["p"][1] - g["end"]["p"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"$T = {st.T0:g}$: "
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
        # The tube as built, in the field (the game's second placement).
        sb = Setup(lvl, [], i)
        rb, cb = fly(sb, H), fly(sb, 2 * H)
        eb = err_of(rb, cb)
        gb, gbp = doc.flight(1, flight=i), doc.flight(1, flight=i, run="preview")
        gapb = max(abs(gbp["end"]["t"] - gb["end"]["t"]), abs(gbp["end"]["x"][0] - gb["end"]["x"][0]),
                   abs(gbp["end"]["x"][1] - gb["end"]["x"][1]))
        tolb = max(TOL, 10 * eb, gapb)
        title = f"as built, $T = {st.T0:g}$: "
        doc.compare(f"BuiltOut{nm}", title + "outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BuiltX{nm}", title + r"$x_\mathrm{e}$", rb["x"][0], gb["end"]["x"][0], tolb)
        doc.compare(f"BuiltY{nm}", title + r"$y_\mathrm{e}$", rb["x"][1], gb["end"]["x"][1], tolb)
        largest = max(largest, abs(rb["x"][0] - gb["end"]["x"][0]), abs(rb["x"][1] - gb["end"]["x"][1]))
        doc.text(f"builtOut{nm}", rb["outcome"])
        doc.text(f"builtFace{nm}", rb["face"])
        doc.value(f"builtX{nm}", rb["x"][0], 5)
        doc.value(f"builtY{nm}", rb["x"][1], 5)
        doc.data("built" + nm.lower(), rb["path"])
        # Heights on the spot's line x = x_S: as built without field (level 7), as built in the
        # field, and with the added charge in the field.
        xs = st.dmin[0]
        y_plain = fly(Setup(lvl, [], i, bz=0.0), H, x_stop=xs)
        y_field = fly(Setup(lvl, [], i), H, x_stop=xs)
        y_fixed = fly(Setup(lvl, ref, i), H, x_stop=xs)
        y_added = fly(Setup(lvl, ref, i, bz=0.0), H, x_stop=xs)
        ys = []
        for r in (y_plain, y_field, y_fixed, y_added):
            ys.append(r["x"][1] if r["outcome"] == "Crossed" else float("nan"))
        lines.append((st.T0, ys[0], ys[1], ys[1] - ys[0], ys[3] - ys[0], ys[2]))
        doc.value(f"shiftB{nm}", ys[1] - ys[0], 4)
        doc.value(f"shiftQ{nm}", ys[3] - ys[0], 4)
        doc.value(f"yPlain{nm}", ys[0], 6)
        doc.value(f"yFixed{nm}", ys[2], 6)
        doc.value(f"sum{nm}", ys[2] - ys[0], 4)
    doc.table("lines", [r"$T$", r"no field", r"field", r"field's shift", r"added charge's shift", r"both"],
              lines, digits=[3, 6, 6, 4, 4, 6])
    doc.game_value("largest", largest, 2, sci=True)
    pv = [st.p0 ** 2 / (st.gamma0 * st.m) for st in sts]
    doc.value("pRatio", sts[1].p0 / sts[0].p0, 4)
    doc.value("pvRatio", pv[1] / pv[0], 4)
    doc.value("TRatio", sts[1].T0 / sts[0].T0, 4)
