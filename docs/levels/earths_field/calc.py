"""Companion of docs/levels/earths_field/doc.tex (Level 43, Earth's field).

Every number of the document, from the level file. Newtonian mechanics (c = infinity). The empty
board in closed form: a straight line, or a circle of radius p/(|q| B) (eq:circle). With the
charges: Newton's equation with the Coulomb force and the uniform B_z (eq:eom-b), marched by the
classical Runge-Kutta method at a fixed step (Foundations eq:rk4), the entry by cubic Hermite
interpolation and bisection, the error against the step 2h, the energy (eq:energy-b) conserved
as the check. The stray field's effect: the arrival height on the line x = x_det as a function
of B_z (marched flights), its slope at B_z = 0 against the empty board's (eq:sag). Compared with
the game: the reference and the empty board under the three disturbances.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.004
H_CURVE = 0.01      # the curve's flights: step error ~1e-9, far below what the curve shows
N_CURVE = 41


class Setup:
    def __init__(self, lvl, placement, dist=None, bz=None):
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
        if dist is not None:
            dz = lvl["disturbances"][dist]
            assert dz["e"] == [0.0, 0.0] and not dz.get("waves")
            self.bz = dz["bz"]
        else:
            self.bz = bz
        self.p0 = math.sqrt(2 * self.m * self.T0)

    def rhs(self, s):
        x, y, px, py = s
        vx, vy = px / self.m, py / self.m
        fx, fy = self.q * vy * self.bz, -self.q * vx * self.bz
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += self.q * Q * dx / r3
            fy += self.q * Q * dy / r3
        return (vx, vy, fx, fy)

    def energy(self, s):
        """eq:energy-b: p^2/2m + sum q Q / r (B does no work)."""
        x, y, px, py = s
        return (px * px + py * py) / (2 * self.m) + sum(self.q * Q / math.hypot(x - cx, y - cy)
                                                        for cx, cy, Q in self.charges)


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
    """The march to the first event; with x_stop, to the crossing of the line x = x_stop."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    if x_stop is not None:
        preds = [("Crossed", 0, lambda p: p[0] >= x_stop)]
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
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else ""
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "face": face, "path": path,
                    "wdev": max(wdev, abs(st.energy(se) - w0))}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "face": "", "path": path,
                    "wdev": wdev}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def circle(st):
    """eq:circle: the empty board. B = 0: a straight line; otherwise a circle of radius
    R = p/(|q| |B|) tangent to the launch line, curving towards sign(-q B) in y. Returns the
    outcome, the end time, the end point and the momentum."""
    v0 = st.p0 / st.m
    if st.bz == 0:
        x_end = st.dmin[0]
        y_end = st.x0[1]
        if st.dmin[1] <= y_end <= st.dmax[1]:
            return "Arrived", (x_end - st.x0[0]) / v0, (x_end, y_end), (st.p0, 0.0)
        return "LeftBounds", (st.wmax[0] - st.x0[0]) / v0, (st.wmax[0], y_end), (st.p0, 0.0)
    R = st.p0 / abs(st.q * st.bz)
    sgn = 1.0 if -st.q * st.bz > 0 else -1.0
    y_of = lambda x: st.x0[1] + sgn * (R - math.sqrt(R * R - (x - st.x0[0]) ** 2))
    y_det = y_of(st.dmin[0])
    # the circle rises (or falls) monotonically over the board: it meets the detector's left
    # face if y_det is in its span, else it may cross a horizontal face between x_det and the
    # detector's far side; neither here (checked below), so it leaves the world's right edge
    doc_y = y_det
    assert not (st.dmin[1] <= doc_y <= st.dmax[1])
    x_end = st.wmax[0]
    theta = math.asin((x_end - st.x0[0]) / R)
    return ("LeftBounds", R * theta / v0, (x_end, y_of(x_end)),
            (st.p0 * math.cos(theta), sgn * st.p0 * math.sin(theta)))


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    dists = lvl["disturbances"]
    doc.check(lvl["physics"]["c"] is None and len(dists) == 3 and len(ref) == 2, "Newtonian, three fields, two charges")
    st0 = Setup(lvl, ref, 0)
    v0 = st0.p0 / st0.m
    for k_, v, dg in (("Tzero", st0.T0, 4), ("m", st0.m, 4), ("xZero", st0.x0[0], 4), ("yZero", st0.x0[1], 4),
                      ("xDetLo", st0.dmin[0], 4), ("yDetLo", st0.dmin[1], 4), ("xDetHi", st0.dmax[0], 4),
                      ("yDetHi", st0.dmax[1], 4), ("xWorld", st0.wmax[0], 4), ("yWorld", st0.wmax[1], 4),
                      ("R", st0.R, 4), ("H", H, 4), ("vZero", v0, 10)):
        doc.value(k_, v, dg)
    doc.value("q", st0.q, 4, sci=True)
    (c1x, c1y, Q1), (c2x, c2y, Q2) = st0.charges
    for k_, v in (("xOne", c1x), ("yOne", c1y), ("xTwo", c2x), ("yTwo", c2y)):
        doc.value(k_, v, 4)
    doc.value("QOne", Q1, 4, sci=True)
    doc.value("QTwo", Q2, 4, sci=True)
    doc.value("kOne", st0.q * Q1, 4)
    doc.value("kTwo", st0.q * Q2, 4)
    region = lvl["limits"]["region"]
    for k_, v in (("xRegLo", region["min"][0]), ("yRegLo", region["min"][1]),
                  ("xRegHi", region["max"][0]), ("yRegHi", region["max"][1])):
        doc.value(k_, v, 4)
    B = max(abs(d["bz"]) for d in dists)
    doc.check([d["bz"] for d in dists] == [0.0, B, -B], "no field, then +B (north), then -B (south)")
    doc.text("distNorth", dists[1]["name"])
    doc.text("distSouth", dists[2]["name"])
    doc.value("B", B, 4)
    Rl = st0.p0 / abs(st0.q * B)
    doc.value("Rlarmor", Rl, 6)
    L = st0.dmin[0] - st0.x0[0]
    doc.value("L", L, 4)
    doc.value("sag", Rl - math.sqrt(Rl * Rl - L * L), 4)
    doc.value("sagApprox", L * L / (2 * Rl), 4)
    doc.value("forceB", abs(st0.q) * v0 * B, 4, sci=True)
    slope_free = L * L * abs(st0.q) / (2 * st0.p0)
    doc.value("slopeFree", slope_free, 4, sci=True)
    doc.value("halfHeight", (st0.dmax[1] - st0.dmin[1]) / 2, 4)

    largest = 0.0
    tags = ["Off", "North", "South"]
    # The empty board in closed form.
    for i, dz in enumerate(dists):
        st = Setup(lvl, [], i)
        out, tc, xc, pc = circle(st)
        run = fly(st, H)
        doc.check(run["outcome"] == out and abs(run["x"][1] - xc[1]) < 1e-9 and abs(run["t"] - tc) < 1e-9,
                  f"the empty board's closed form, {dz['name']}")
        g = doc.flight(1, flight=i)
        title = f"empty board, {dz['name']}: "
        doc.compare(f"FreeOut{i}", title + "outcome", out, g["outcome"])
        doc.compare(f"FreeT{i}", title + r"$t_\mathrm{e}$", tc, g["end"]["t"], TOL)
        doc.compare(f"FreeY{i}", title + r"$y_\mathrm{e}$", xc[1], g["end"]["x"][1], TOL)
        doc.compare(f"FreePy{i}", title + r"$p_{y,\mathrm{e}}$", pc[1], g["end"]["p"][1], TOL)
        doc.value(f"freeY{tags[i]}", xc[1], 6)
        doc.value(f"freeT{tags[i]}", tc, 6)
        largest = max(largest, abs(tc - g["end"]["t"]), abs(xc[1] - g["end"]["x"][1]), abs(pc[1] - g["end"]["p"][1]))

    # The reference under the three fields.
    for i, dz in enumerate(dists):
        st = Setup(lvl, ref, i)
        run, coarse = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", f"the reference arrives, {dz['name']}")
        err = err_of(run, coarse)
        tag = tags[i]
        doc.value(f"t{tag}", run["t"], 10)
        doc.value(f"x{tag}", run["x"][0], 10)
        doc.value(f"y{tag}", run["x"][1], 10)
        doc.value(f"px{tag}", run["p"][0], 8)
        doc.value(f"py{tag}", run["p"][1], 8)
        doc.text(f"face{tag}", run["face"])
        doc.value(f"err{tag}", err, 2, sci=True)
        doc.value(f"wDev{tag}", run["wdev"], 2, sci=True)
        doc.check(run["wdev"] < 1e-11, "the energy is conserved")
        doc.data("path" + tag.lower(), run["path"])
        g, gp = doc.flight(0, flight=i), doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]), abs(gp["end"]["p"][0] - g["end"]["p"][0]),
                  abs(gp["end"]["p"][1] - g["end"]["p"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"{dz['name']}: "
        doc.compare(f"Out{i}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{i}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{i}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{i}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Px{i}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
        doc.compare(f"Py{i}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        largest = max(largest, abs(run["t"] - g["end"]["t"]), abs(run["x"][0] - g["end"]["x"][0]),
                      abs(run["x"][1] - g["end"]["x"][1]), abs(run["p"][0] - g["end"]["p"][0]),
                      abs(run["p"][1] - g["end"]["p"][1]))
    doc.game_value("largest", largest, 2, sci=True)

    # How the arrival height on x = x_det depends on B_z, with the reference's charges.
    curve = []
    for k in range(N_CURVE):
        b = -1.5 * B + 3.0 * B * k / (N_CURVE - 1)
        r = fly(Setup(lvl, ref, bz=b), H_CURVE, x_stop=st0.dmin[0])
        curve.append((b, r["x"][1] if r["outcome"] == "Crossed" else float("nan")))
    doc.data("curve", [(b, y) for b, y in curve if y == y])
    doc.data("curvefree", [(b, st0.x0[1] + (1 if -st0.q * b > 0 else -1) *
                            (st0.p0 / abs(st0.q * b) - math.sqrt((st0.p0 / abs(st0.q * b)) ** 2 - L * L))
                            if b != 0 else st0.x0[1]) for b, _ in curve])
    eps = B * 1e-3
    y_at = lambda b: fly(Setup(lvl, ref, bz=b), H, x_stop=st0.dmin[0])["x"][1]
    slope = (y_at(eps) - y_at(-eps)) / (2 * eps)
    doc.value("slopeRef", slope, 4, sci=True)
    doc.value("slopeRatio", slope_free / slope, 3)
    y0line = y_at(0.0)
    doc.value("yLineZero", y0line, 6)
    doc.value("yLineNorth", y_at(B), 6)
    doc.value("yLineSouth", y_at(-B), 6)
    doc.value("spreadRef", y_at(B) - y_at(-B), 4)
    doc.value("spreadFree", 2 * (Rl - math.sqrt(Rl * Rl - L * L)), 4)
    # the curve's maximum: golden-section search between the samples around the largest
    k_max = max(range(len(curve)), key=lambda k: curve[k][1] if curve[k][1] == curve[k][1] else -1e9)
    a_, b_ = curve[max(k_max - 1, 0)][0], curve[min(k_max + 1, len(curve) - 1)][0]
    gr = (math.sqrt(5) - 1) / 2
    yc = lambda b: fly(Setup(lvl, ref, bz=b), H_CURVE, x_stop=st0.dmin[0])["x"][1]
    c_, d_ = b_ - gr * (b_ - a_), a_ + gr * (b_ - a_)
    fc, fd = yc(c_), yc(d_)
    while b_ - a_ > 1.0:
        if fc > fd:
            b_, d_, fd = d_, c_, fc
            c_ = b_ - gr * (b_ - a_)
            fc = yc(c_)
        else:
            a_, c_, fc = c_, d_, fd
            d_ = a_ + gr * (b_ - a_)
            fd = yc(d_)
    b_top = 0.5 * (a_ + b_)
    doc.check(0 < b_top < B, "the maximum lies between no field and facing north")
    doc.value("foldB", b_top, 3)
    doc.value("foldY", yc(b_top), 4)
    # the closest approach to the upper charge, field off and facing north and south
    for tag, b in (("Off", 0.0), ("North", B), ("South", -B)):
        r = fly(Setup(lvl, ref, bz=b), H_CURVE, x_stop=st0.dmin[0])
        doc.value(f"near{tag}", min(math.hypot(x - c2x, y - c2y) for x, y in r["path"]), 3)
