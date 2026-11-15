"""Companion of docs/levels/mains_hum/doc.tex (Level 42, Mains hum).

Every number of the document, from the level file. Newtonian mechanics (c = infinity), where
each disturbance's plane wave is the uniform field E0 e cos(omega t + phi) (eq:hum). The empty
board in closed form (eq:free). With the charges: the exact split of the motion into the quiver
and the guiding motion (eq:split), the guiding motion's start per phase (eq:effective); the full
equations marched by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4), the
entry by cubic Hermite interpolation and bisection, the error against the step 2h. The
approximation that drops the quiver from the charges' force (eq:guiding): the static flights from
the effective starts, measured against the full march. The lens: the arrival height's
sensitivity to the start's angle and height (eq:lens), by central differences of static flights.
The phase sweep: the reference at 72 phases. Compared with the game: the reference and the
empty board under the four disturbances.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.004
H_SWEEP = 0.01      # the sweep's flights: step error ~1e-9, far below the margins it reports
N_SWEEP = 72


class Setup:
    def __init__(self, lvl, placement, dist, phase_deg=None):
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
        dz = lvl["disturbances"][dist]
        (w,) = dz["waves"]
        assert dz["e"] == [0.0, 0.0] and dz["bz"] == 0.0
        k = math.radians(w["direction_deg"])
        self.e_hat = (-math.sin(k), math.cos(k))                      # z x k-hat
        self.E0, self.omega = w["amplitude"], w["omega"]
        self.phi = math.radians(w["phase_deg"] if phase_deg is None else phase_deg)
        self.p0 = math.sqrt(2 * self.m * self.T0)
        self.static = False                                          # True: the charges alone

    def hum(self, t):
        if self.static:
            return 0.0, 0.0
        a = self.E0 * math.cos(self.omega * t + self.phi)
        return a * self.e_hat[0], a * self.e_hat[1]

    def rhs(self, s, t):
        x, y, px, py = s
        hx, hy = self.hum(t)
        fx, fy = self.q * hx, self.q * hy
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += self.q * Q * dx / r3
            fy += self.q * Q * dy / r3
        return (px / self.m, py / self.m, fx, fy)


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


def fly(st, h, start=None, x_stop=None):
    """The march from the launch (or `start` = (x, y, px, py)) to the first event; with x_stop,
    to the crossing of the line x = x_stop instead (the guiding flights)."""
    s = start or (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
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
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "face": face, "path": path}
        s, f, t, n = s2, f2, t + hs, n + 1
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "face": "", "path": path}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def tidy(x):
    """For the tables: the rounding of cos 90 deg and sin 180 deg (~1e-17) shown as 0."""
    return 0.0 if abs(x) < 1e-12 else x


def free_exit(st):
    """eq:free: the empty board, x = x0 + v0 t and y from the quiver and drift; the flight ends
    where the detector or the world's right edge is first reached (both found here by bisection
    on the closed form)."""
    v0 = st.p0 / st.m
    A = st.q * st.E0 / st.m
    w, phi = st.omega, st.phi
    ey = st.e_hat[1]
    y = lambda t: st.x0[1] + ey * (A / (w * w) * (math.cos(phi) - math.cos(w * t + phi)) - A / w * math.sin(phi) * t)
    t_det = (st.dmin[0] - st.x0[0]) / v0
    # the detector's left face, or its top or bottom face between x = 28 and 30, or the edge
    ts = [k * 1e-3 for k in range(int(((st.wmax[0] - st.x0[0]) / v0) * 1000) + 1)]
    inside = lambda t: box_distance((st.x0[0] + v0 * t, y(t)), st.dmin, st.dmax) <= 0
    for a_, b_ in zip(ts, ts[1:]):
        if b_ >= t_det - 1e-3 and inside(b_):
            lo, hi = a_, b_
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if inside(mid):
                    hi = mid
                else:
                    lo = mid
            return "Arrived", hi, st.x0[0] + v0 * hi, y(hi)
    t_edge = (st.wmax[0] - st.x0[0]) / v0
    return "LeftBounds", t_edge, st.wmax[0], y(t_edge)


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    dists = lvl["disturbances"]
    doc.check(lvl["physics"]["c"] is None and len(dists) == 4 and len(ref) == 2, "Newtonian, four phases, two charges")
    st0 = Setup(lvl, ref, 0)
    A = st0.q * st0.E0 / st0.m
    w = st0.omega
    v0 = st0.p0 / st0.m
    for k_, v, dg in (("q", st0.q, 4), ("m", st0.m, 4), ("Tzero", st0.T0, 4), ("xZero", st0.x0[0], 4),
                      ("yZero", st0.x0[1], 4), ("xDetLo", st0.dmin[0], 4), ("yDetLo", st0.dmin[1], 4),
                      ("xDetHi", st0.dmax[0], 4), ("yDetHi", st0.dmax[1], 4), ("xWorld", st0.wmax[0], 4),
                      ("yWorld", st0.wmax[1], 4), ("R", st0.R, 4), ("omega", w, 4), ("H", H, 4),
                      ("vZero", v0, 10)):
        doc.value(k_, v, dg)
    doc.value("Ezero", st0.E0, 4, sci=True)
    (c1x, c1y, Q1), (c2x, c2y, Q2) = st0.charges
    for k_, v in (("xOne", c1x), ("yOne", c1y), ("xTwo", c2x), ("yTwo", c2y)):
        doc.value(k_, v, 4)
    doc.value("QOne", Q1, 4, sci=True)
    doc.value("QTwo", Q2, 4, sci=True)
    region = lvl["limits"]["region"]
    doc.value("xRegLo", region["min"][0], 4)
    doc.value("yRegLo", region["min"][1], 4)
    doc.value("xRegHi", region["max"][0], 4)
    doc.value("yRegHi", region["max"][1], 4)
    doc.value("period", 2 * math.pi / w, 4)
    doc.value("accel", A, 4)
    doc.value("quiverV", A / w, 4)
    doc.value("quiverX", A / (w * w), 4)
    doc.value("tCross", (st0.dmin[0] - st0.x0[0]) / v0, 4)
    doc.value("driftAt", A / w * (st0.dmin[0] - st0.x0[0]) / v0, 4)
    doc.value("halfHeight", (st0.dmax[1] - st0.dmin[1]) / 2, 4)
    doc.value("periods", (st0.dmin[0] - st0.x0[0]) / v0 / (2 * math.pi / w), 3)

    # The empty board in closed form, against the march and the game.
    rows_free = []
    largest = 0.0
    for i, dz in enumerate(dists):
        st = Setup(lvl, [], i)
        out, tf, xf, yf = free_exit(st)
        run = fly(st, H)
        doc.check(run["outcome"] == out and abs(run["x"][1] - yf) < 1e-9 and abs(run["t"] - tf) < 1e-9,
                  f"the empty board's closed form, phase {i}")
        rows_free.append((dz["waves"][0]["phase_deg"], tidy(A / (w * w) * math.cos(st.phi)), tidy(-A / w * math.sin(st.phi)), yf, out))
        g = doc.flight(1, flight=i)
        title = f"empty board, {dz['name']}: "
        doc.compare(f"FreeOut{i}", title + "outcome", out, g["outcome"])
        doc.compare(f"FreeY{i}", title + r"$y_\mathrm{e}$", yf, g["end"]["x"][1], TOL)
        largest = max(largest, abs(yf - g["end"]["x"][1]))
    doc.table("free", [r"$\varphi$ (deg)", r"$\delta y_0$", r"$\delta v_y$", r"$y_\mathrm{e}$", "outcome"],
              rows_free, digits=[3, 6, 6, 8, 0])

    # The reference under the four phases: the full march, and the guiding approximation.
    rows = []
    guide_worst = 0.0
    tags = ["A", "B", "C", "D"]
    for i, dz in enumerate(dists):
        st = Setup(lvl, ref, i)
        run, coarse = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", f"the reference arrives under {dz['name']}")
        err = err_of(run, coarse)
        # eq:effective, eq:guiding: the static flight from the shifted start to x = x_e, plus the
        # quiver there
        dy0, dvy = A / (w * w) * math.cos(st.phi) * st.e_hat[1], -A / w * math.sin(st.phi) * st.e_hat[1]
        sg = Setup(lvl, ref, i)
        sg.static = True
        gs = fly(sg, H, start=(st.x0[0], st.x0[1] + dy0, st.p0 * st.d[0], st.p0 * st.d[1] + st.m * dvy),
                 x_stop=run["x"][0])
        y_guess = gs["x"][1] - A / (w * w) * math.cos(w * gs["t"] + st.phi) * st.e_hat[1]
        guide_worst = max(guide_worst, abs(y_guess - run["x"][1]))
        closest = min(math.hypot(x - c2x, y - c2y) for x, y in run["path"])
        rows.append((dz["waves"][0]["phase_deg"], tidy(dy0), tidy(dvy), run["t"], run["x"][1], y_guess, closest))
        doc.value(f"closest{tags[i]}", closest, 3)
        doc.value(f"guide{tags[i]}", y_guess, 4)
        doc.value(f"t{tags[i]}", run["t"], 10)
        doc.value(f"y{tags[i]}", run["x"][1], 10)
        doc.value(f"err{tags[i]}", err, 2, sci=True)
        doc.text(f"face{tags[i]}", run["face"])
        doc.data("path" + tags[i].lower(), run["path"])
        g, gp = doc.flight(0, flight=i), doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][1] - g["end"]["x"][1]),
                  abs(gp["end"]["p"][0] - g["end"]["p"][0]), abs(gp["end"]["p"][1] - g["end"]["p"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"{dz['name']}: "
        doc.compare(f"Out{i}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{i}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{i}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{i}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Px{i}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
        doc.compare(f"Py{i}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        largest = max(largest, abs(run["t"] - g["end"]["t"]), abs(run["x"][1] - g["end"]["x"][1]),
                      abs(run["p"][0] - g["end"]["p"][0]), abs(run["p"][1] - g["end"]["p"][1]))
    doc.table("ref", [r"$\varphi$ (deg)", r"$\delta y_0$", r"$\delta v_y$", r"$t_\mathrm{e}$", r"$y_\mathrm{e}$",
                      r"guiding $y_\mathrm{e}$", r"$d_2$"], rows, digits=[3, 4, 3, 9, 9, 5, 4])
    doc.value("guideWorst", guide_worst, 2)
    faces = {doc.values[f"face{t}"] for t in tags}
    doc.check(len(faces) == 1, "every phase enters through the same face")
    doc.text("faceAll", faces.pop())
    doc.game_value("largest", largest, 2, sci=True)
    margin = min(min(r[4] - st0.dmin[1], st0.dmax[1] - r[4]) for r in rows)
    doc.check(margin > 0, "every phase arrives inside the detector's span")
    doc.value("margin", margin, 2)

    # The lens: how the arrival height on the line x = x_det depends on the start's angle and height.
    sl = Setup(lvl, ref, 0)
    sl.static = True
    x_det = sl.dmin[0]

    def y_at(dy, dvy):
        r = fly(sl, H, start=(sl.x0[0], sl.x0[1] + dy, sl.p0, sl.m * dvy), x_stop=x_det)
        return r["x"][1]
    eps = 1e-3
    dydv = (y_at(0, eps) - y_at(0, -eps)) / (2 * eps)
    dydy = (y_at(eps, 0) - y_at(-eps, 0)) / (2 * eps)
    doc.value("yStatic", y_at(0, 0), 8)
    doc.value("dydv", dydv, 4)
    doc.value("dydy", dydy, 4)
    doc.value("freeDydv", x_det / v0, 4)
    doc.value("spreadV", abs(dydv) * A / w, 3)
    doc.value("spreadY", abs(dydy) * A / (w * w), 3)
    # the static arrival height on x = x_det against the start's transverse velocity
    curve = []
    for k in range(-24, 25):
        dv = k / 200
        r = fly(sl, H_SWEEP, start=(sl.x0[0], sl.x0[1], sl.p0, sl.m * dv), x_stop=x_det)
        if r["outcome"] == "Crossed":
            curve.append((dv, r["x"][1], min(math.hypot(x - c2x, y - c2y) for x, y in r["path"])))
    doc.data("curve", [(dv, y) for dv, y, _ in curve if y < 14])
    lowest = min(curve, key=lambda c: c[1])
    doc.value("foldV", lowest[0], 3)
    doc.value("foldY", lowest[1], 4)
    at = {round(dv, 3): (y, d) for dv, y, d in curve}
    doc.value("yAtMinus", at[-0.1][0], 4)
    doc.value("yAtPlus", at[0.1][0], 4)
    doc.value("closestMinus", at[-0.1][1], 3)
    doc.value("yAtSteep", at[-0.12][0], 4)
    doc.value("closestSteep", at[-0.12][1], 3)

    # Every phase: the reference at N_SWEEP phases.
    sweep = []
    for k in range(N_SWEEP):
        ph = 360.0 * k / N_SWEEP
        r = fly(Setup(lvl, ref, 0, phase_deg=ph), H_SWEEP)
        sweep.append((ph, r["x"][1] if r["outcome"] == "Arrived" else float("nan"), r["outcome"]))
    ok = [s for s in sweep if s[2] == "Arrived"]
    doc.value("nSweep", N_SWEEP, 4)
    doc.value("nSweepOk", len(ok), 4)
    doc.data("sweep", [(p, y) for p, y, o in sweep if o == "Arrived"])
    missed = [s for s in sweep if s[2] != "Arrived"]
    doc.text("missedPhases", ", ".join(f"{p:g}" for p, _, _ in missed) if missed else "none")
    doc.value("sweepMin", min(s[1] for s in ok), 4)
    doc.value("sweepMax", max(s[1] for s in ok), 4)
