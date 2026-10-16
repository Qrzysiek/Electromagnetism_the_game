"""Companion of docs/levels/reflectron/doc.tex (Level 10, Reflectron).

Two fixed centres, so no closed orbit (slingshot). Each ion is the classical Runge-Kutta
method at a fixed step (Foundations eq:rk4); the face crossing is the cubic Hermite
interpolant of that step. The truncation error is estimated by halving the step, and the
tolerance is set from that estimate before any comparison with the game.
"""
import math


# Reported step. Halving is the coarser run below: the error of the reported step is
# |y_h - y_{2h}|/15 (Foundations, the numerics section). The tolerance is ten times the
# largest such estimate and at least 1e-9, ten times the game's verified accuracy of
# about 1e-10. Computed in compute() before doc.compare, not from the game's numbers.
H = 0.005
# Simpson on the radiated power. 1e-4 still separates a wrong formula from the game's
# quadrature of the same power (about 3e-9 relative on an exact test). Set beforehand.
TOL_RAD = 1e-4


class Setup:
    """One shot and the charges of one placement, in the game's order."""

    def __init__(self, lvl, placement, shot_index):
        shot = lvl["shots"][shot_index]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.a = shot["particle"]["radius"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(d[0], d[1])
        self.d = [d[0] / n, d[1] / n]
        self.dir_len = n
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        self.kappas = tuple((self.q * Q, cx, cy) for cx, cy, Q in self.charges)
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)
        self.angle = math.degrees(math.atan2(self.d[1], self.d[0]))

    def V(self, x, y):
        total = 0.0
        for k, cx, cy in self.kappas:
            total += k / math.hypot(x - cx, y - cy)
        return total

    def rhs(self, state):
        x, y, px, py = state
        m = self.m
        c = self.c
        p2 = px * px + py * py
        g = math.sqrt(1.0 + p2 / (m * m * c * c))
        inv = 1.0 / (g * m)
        fx = fy = 0.0
        for k, cx, cy in self.kappas:
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            inv_r3 = k / (r2 * math.sqrt(r2))
            fx += dx * inv_r3
            fy += dy * inv_r3
        return (px * inv, py * inv, fx, fy)

    def energy(self, state):
        x, y, px, py = state
        rest = math.sqrt(self.m * self.m * self.c ** 4 + (px * px + py * py) * self.c * self.c)
        return rest + self.V(x, y)

    def kinetic(self, state):
        """(gamma - 1) m c^2 from the conserved energy and V at the position."""
        x, y, px, py = state
        return self.energy(state) - self.m * self.c * self.c - self.V(x, y)

    def launch_state(self):
        return (self.x0[0], self.x0[1], self.p0 * self.d[0], self.p0 * self.d[1])


def rk4_step(st, s, f, h):
    x, y, px, py = s
    k1x, k1y, k1px, k1py = f
    k2x, k2y, k2px, k2py = st.rhs((x + 0.5 * h * k1x, y + 0.5 * h * k1y,
                                   px + 0.5 * h * k1px, py + 0.5 * h * k1py))
    k3x, k3y, k3px, k3py = st.rhs((x + 0.5 * h * k2x, y + 0.5 * h * k2y,
                                   px + 0.5 * h * k2px, py + 0.5 * h * k2py))
    k4x, k4y, k4px, k4py = st.rhs((x + h * k3x, y + h * k3y, px + h * k3px, py + h * k3py))
    return (x + h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
            y + h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y),
            px + h / 6.0 * (k1px + 2.0 * k2px + 2.0 * k3px + k4px),
            py + h / 6.0 * (k1py + 2.0 * k2py + 2.0 * k3py + k4py))


def hermite(s, f, s2, f2, h, u):
    u2 = u * u
    u3 = u2 * u
    h00 = 2.0 * u3 - 3.0 * u2 + 1.0
    h10 = u3 - 2.0 * u2 + u
    h01 = -2.0 * u3 + 3.0 * u2
    h11 = u3 - u2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def _first(pred, pts, us):
    for u, p in zip(us, pts):
        if pred(p):
            return u
    return None


def _bisect(pred, s, f, s2, f2, h, lo, hi):
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pred(hermite(s, f, s2, f2, h, mid)):
            hi = mid
        else:
            lo = mid
    return hi


def assign_face(pe, lo, hi):
    x, y = pe
    segs = (
        (abs(x - lo[0]), "left", 0, lo[0], lo[1] <= y <= hi[1]),
        (abs(x - hi[0]), "right", 0, hi[0], lo[1] <= y <= hi[1]),
        (abs(y - lo[1]), "bottom", 1, lo[1], lo[0] <= x <= hi[0]),
        (abs(y - hi[1]), "top", 1, hi[1], lo[0] <= x <= hi[0]),
    )
    inside = [seg for seg in segs if seg[4]]
    res, name, axis, val, _ = min(inside or segs, key=lambda seg: seg[0])
    if res < 1e-6:
        out = [x, y]
        out[axis] = val
        return (out[0], out[1]), name, res
    return (x, y), name, res


def fly(st, h, path_every=0, table_every=0.0):
    """Integrate until the first event. Board order: collision, world, detector."""
    s = st.launch_state()
    f = st.rhs(s)
    t = 0.0
    n = 0
    w0 = st.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    rows = [(0.0, s[0], s[1], s[2], s[3], 0.0)]
    next_table = table_every
    us = [k / 64.0 for k in range(1, 65)]
    limit = st.R + st.a
    t_min = st.T0
    t_min_at = (0.0, s[0], s[1])
    # Closest approach to a centre, on the same nine samples as the kinetic dip.
    r_close = math.inf
    r_at = (0.0, s[0], s[1])
    r_idx = 0
    # Rightmost sample: the turnaround of a reflected ion, or the exit of one that is not.
    x_turn = s[0]
    turn_at = (0.0, s[0], s[1])
    while t <= st.tmax:
        s2 = rk4_step(st, s, f, h)
        f2 = st.rhs(s2)
        margin = st.c * h
        lo_x, hi_x = (s[0], s2[0]) if s[0] <= s2[0] else (s2[0], s[0])
        lo_y, hi_y = (s[1], s2[1]) if s[1] <= s2[1] else (s2[1], s[1])
        box_lo = (lo_x - margin, lo_y - margin)
        box_hi = (hi_x + margin, hi_y + margin)
        hits = []
        pts = None

        def samples():
            nonlocal pts
            if pts is None:
                pts = [hermite(s, f, s2, f2, h, u) for u in us]
            return pts

        for idx, (_, cx, cy) in enumerate(st.kappas):
            dx = 0.0 if box_lo[0] <= cx <= box_hi[0] else (box_lo[0] - cx if cx < box_lo[0] else cx - box_hi[0])
            dy = 0.0 if box_lo[1] <= cy <= box_hi[1] else (box_lo[1] - cy if cy < box_lo[1] else cy - box_hi[1])
            if dx * dx + dy * dy > limit * limit:
                continue
            u = _first(lambda p, cx=cx, cy=cy: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit,
                       samples(), us)
            if u is not None:
                hits.append((u, 0, f"Collided({idx})", idx))
        if box_lo[0] < st.wmin[0] or box_hi[0] > st.wmax[0] or box_lo[1] < st.wmin[1] or box_hi[1] > st.wmax[1]:
            u = _first(lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1]),
                       samples(), us)
            if u is not None:
                hits.append((u, 1, "LeftBounds", None))
        if not (box_hi[0] < st.dmin[0] or box_lo[0] > st.dmax[0] or box_hi[1] < st.dmin[1] or box_lo[1] > st.dmax[1]):
            u = _first(lambda p: st.dmin[0] <= p[0] <= st.dmax[0] and st.dmin[1] <= p[1] <= st.dmax[1],
                       samples(), us)
            if u is not None:
                hits.append((u, 2, "Arrived", None))
        # Minimum kinetic energy on the step, from the conserved W and V on the interpolant.
        # V is smooth; eight samples locate the dip well inside the step error.
        if st.kappas:
            for k in range(0, 9):
                p = hermite(s, f, s2, f2, h, k / 8.0)
                tk = (w0 - st.m * st.c * st.c) - st.V(p[0], p[1])
                if tk < t_min:
                    t_min = tk
                    t_min_at = (t + k / 8.0 * h, p[0], p[1])
                for idx, (_, cx, cy) in enumerate(st.kappas):
                    rr = math.hypot(p[0] - cx, p[1] - cy)
                    if rr < r_close:
                        r_close = rr
                        r_at = (t + k / 8.0 * h, p[0], p[1])
                        r_idx = idx
                if p[0] > x_turn:
                    x_turn = p[0]
                    turn_at = (t + k / 8.0 * h, p[0], p[1])
        if hits:
            u_hit, _rank, outcome, idx = min(hits)
            if outcome == "Arrived":
                pred = lambda p: st.dmin[0] <= p[0] <= st.dmax[0] and st.dmin[1] <= p[1] <= st.dmax[1]
                box = (st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                pred = lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])
                box = (st.wmin, st.wmax)
            else:
                cx, cy = st.kappas[idx][1], st.kappas[idx][2]
                pred = lambda p, cx=cx, cy=cy: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit
                box = None
            hi = _bisect(pred, s, f, s2, f2, h, u_hit - 1.0 / 64.0, u_hit)
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            face, residual = "", 0.0
            if box is not None:
                pe, face, residual = assign_face(pe, box[0], box[1])
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n + 1,
                    "face": face, "residual": residual, "u": hi, "t_before": t,
                    "before": (s[0], s[1]), "after": (s2[0], s2[1]),
                    "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                    "t_min": t_min, "t_min_at": t_min_at,
                    "r_close": r_close, "r_at": r_at, "r_idx": r_idx,
                    "turn_at": turn_at}
        s, f = s2, f2
        t += h
        n += 1
        drift = abs(st.energy(s) - w0)
        if drift > wdev:
            wdev = drift
        if table_every and t + 1e-12 >= next_table:
            rows.append((t, s[0], s[1], s[2], s[3], st.energy(s) - w0))
            next_table += table_every
        if path_every and n % path_every == 0:
            path.append((s[0], s[1]))
    return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": (s[2], s[3]), "steps": n,
            "face": "", "residual": 0.0, "u": 1.0, "t_before": t, "before": (s[0], s[1]),
            "after": (s[0], s[1]), "rows": rows, "path": path, "wdev": wdev, "w0": w0,
            "t_min": t_min, "t_min_at": t_min_at,
            "r_close": r_close, "r_at": r_at, "r_idx": r_idx,
            "turn_at": turn_at}


def radiated(st, h, t_end):
    m = st.m
    c = st.c
    q = st.q

    def power(state):
        x, y, px, py = state
        g = math.sqrt(1.0 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        _, _, fx, fy = st.rhs(state)
        vF = vx * fx + vy * fy
        ax = (fx - vx * vF / (c * c)) / (g * m)
        ay = (fy - vy * vF / (c * c)) / (g * m)
        vxa = vx * ay - vy * ax
        return (2.0 / 3.0) * q * q / (c ** 3) * g ** 6 * (ax * ax + ay * ay - vxa * vxa / (c * c))

    n = int(t_end / h)
    n -= n % 2
    state = st.launch_state()
    total = power(state)
    for i in range(1, n + 1):
        state = rk4_step(st, state, st.rhs(state), h)
        w = 1.0 if i == n else (4.0 if i % 2 else 2.0)
        total += w * power(state)
    simpson = h / 3.0 * total
    rest = t_end - n * h
    state2 = rk4_step(st, state, st.rhs(state), rest) if rest > 0.0 else state
    return simpson + rest * 0.5 * (power(state) + power(state2))


def step_gap(fine, coarse):
    return max(abs(fine["t"] - coarse["t"]),
               abs(fine["x"][0] - coarse["x"][0]),
               abs(fine["x"][1] - coarse["x"][1]),
               abs(fine["p"][0] - coarse["p"][0]),
               abs(fine["p"][1] - coarse["p"][1])) / 15.0


def pair(st, h, path_every=0, table_every=0.0):
    fine = fly(st, h, path_every=path_every, table_every=table_every)
    coarse = fly(st, 2.0 * h)
    return fine, coarse


def compute(doc):
    lvl = doc.level
    places = [lvl["reference_solution"],
              doc.result(1)["placement"],
              doc.result(2)["placement"]]
    names = ("Ref", "Lo", "Hi")
    # shot 0 is the slower ion, shot 1 the faster one. Both launch along +x.
    shot_tags = ("Slow", "Fast")
    runs = {}
    estimate = 0.0
    for pi, place in enumerate(places):
        for si in range(2):
            st = Setup(lvl, place, si)
            path_every = 8 if (pi == 0 or si == 0) else 0
            table_every = 5.0 if pi == 0 and si == 0 else 0.0
            fine, coarse = pair(st, H, path_every=path_every, table_every=table_every)
            doc.check(fine["outcome"] == coarse["outcome"],
                      f"{names[pi]} shot {si} has the same outcome at h and at 2h")
            doc.check(fine["residual"] < 1e-9, f"{names[pi]} shot {si} crosses a face cleanly")
            estimate = max(estimate, step_gap(fine, coarse))
            runs[(pi, si)] = (st, fine, coarse)
    tol = max(1e-9, 10.0 * estimate)

    slow = runs[(0, 0)][0]
    fast = runs[(0, 1)][0]
    doc.check(slow.x0 == fast.x0, "both ions leave the same point")
    doc.check(slow.dmin == fast.dmin and slow.dmax == fast.dmax,
              "both ions share the detector")
    doc.value("q", slow.q, 3, sci=True)
    doc.value("m", slow.m, 3)
    doc.value("c", slow.c, 3)
    doc.value("Tslow", slow.T0, 3)
    doc.value("Tfast", fast.T0, 3)
    doc.value("R", slow.R, 3)
    doc.value("aRad", slow.a, 3)
    doc.value("tmax", slow.tmax, 3)
    doc.value("xLaunch", slow.x0[0], 3)
    doc.value("yLaunch", slow.x0[1], 3)
    doc.value("gammaSlow", slow.gamma0, 8)
    doc.value("gammaFast", fast.gamma0, 8)
    doc.value("pSlow", slow.p0, 12)
    doc.value("pFast", fast.p0, 12)
    doc.value("vSlow", slow.p0 / (slow.gamma0 * slow.m), 8)
    doc.value("vFast", fast.p0 / (fast.gamma0 * fast.m), 8)
    doc.value("xWorld", slow.wmax[0], 3)
    doc.value("yWorld", slow.wmax[1], 3)
    doc.value("xDlo", slow.dmin[0], 3)
    doc.value("xDhi", slow.dmax[0], 3)
    doc.value("yDlo", slow.dmin[1], 3)
    doc.value("yDhi", slow.dmax[1], 3)
    doc.value("H", H, 4)
    doc.value("tolFlight", tol, 2, sci=True)
    # Largest of the six step-halving estimates, recorded before any comparison.
    doc.value("errMax", estimate, 2, sci=True)
    doc.value("maxCharges", lvl["limits"]["max_charges"])
    doc.value("nSample", 64)
    doc.value("nBisect", 80)
    region = lvl["limits"]["region"]
    doc.value("xRegLo", region["min"][0], 3)
    doc.value("xRegHi", region["max"][0], 3)
    doc.value("yRegLo", region["min"][1], 3)
    doc.value("yRegHi", region["max"][1], 3)
    (x_lo, y_lo, q_lo), (x_hi, y_hi, q_hi) = slow.charges
    doc.value("Qlo", q_lo, 3, sci=True)
    doc.value("Qhi", q_hi, 3, sci=True)
    doc.value("kLo", slow.q * q_lo, 4)
    doc.value("kHi", slow.q * q_hi, 4)
    doc.value("xLo", x_lo, 3)
    doc.value("yLo", y_lo, 3)
    doc.value("xHi", x_hi, 3)
    doc.value("yHi", y_hi, 3)
    doc.value("nCharges", len(slow.charges))

    V0 = slow.V(*slow.x0)
    doc.value("VZero", V0, 8)
    for tag, st in (("Slow", slow), ("Fast", fast)):
        W = st.gamma0 * st.m * st.c ** 2 + V0
        doc.value("W" + tag, W, 10)
        doc.value("Wbind" + tag, st.T0 + V0, 8)
    f0 = slow.rhs(slow.launch_state())
    doc.value("FxZero", f0[2], 6)
    doc.value("FyZero", f0[3], 6)

    st, fine, coarse = runs[(0, 0)]
    doc.value("steps", fine["steps"])
    doc.value("tE", fine["t"], 10)
    doc.value("xE", fine["x"][0], 10)
    doc.value("yE", fine["x"][1], 10)
    doc.value("pxE", fine["p"][0], 10)
    doc.value("pyE", fine["p"][1], 10)
    doc.text("faceE", fine["face"])
    doc.value("faceRes", fine["residual"], 2, sci=True)
    doc.value("uE", fine["u"], 6)
    doc.value("tBefore", fine["t_before"], 6)
    doc.value("xBefore", fine["before"][0], 6)
    doc.value("yBefore", fine["before"][1], 6)
    doc.value("xAfter", fine["after"][0], 6)
    doc.value("yAfter", fine["after"][1], 6)
    doc.value("wDev", fine["wdev"], 2, sci=True)
    doc.value("TminSlow", fine["t_min"], 6)
    doc.value("tTmin", fine["t_min_at"][0], 6)
    doc.value("xTmin", fine["t_min_at"][1], 6)
    doc.value("yTmin", fine["t_min_at"][2], 6)
    doc.value("tTurnSlow", fine["turn_at"][0], 6)
    doc.value("xTurnSlow", fine["turn_at"][1], 6)
    doc.value("yTurnSlow", fine["turn_at"][2], 6)
    ang_e = math.degrees(math.atan2(fine["p"][1], fine["p"][0]))
    doc.value("angE", ang_e, 6)
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$\Delta W$"],
              list(fine["rows"]),
              digits=[5, 5, 5, 5, 5, 2],
              sci=[False, False, False, False, False, True])
    fast_fine = runs[(0, 1)][1]
    doc.value("TminFast", fast_fine["t_min"], 6)
    doc.value("tTurnFast", fast_fine["turn_at"][0], 6)
    doc.value("xTurnFast", fast_fine["turn_at"][1], 6)
    doc.value("yTurnFast", fast_fine["turn_at"][2], 6)
    doc.data("pathslow", fine["path"])
    doc.data("pathfast", fast_fine["path"])
    doc.data("pathlo", runs[(1, 0)][1]["path"])
    doc.data("pathhi", runs[(2, 0)][1]["path"])

    def publish(pi, si):
        st, fine, coarse = runs[(pi, si)]
        prefix = names[pi] + shot_tags[si]
        doc.text(prefix + "Outcome", fine["outcome"])
        doc.text(prefix + "Face", fine["face"])
        doc.value(prefix + "T", fine["t"], 8)
        doc.value(prefix + "X", fine["x"][0], 8)
        doc.value(prefix + "Y", fine["x"][1], 8)
        doc.value(prefix + "Px", fine["p"][0], 8)
        doc.value(prefix + "Py", fine["p"][1], 8)
        doc.value(prefix + "Err", step_gap(fine, coarse), 2, sci=True)
        doc.value(prefix + "Ang", math.degrees(math.atan2(fine["p"][1], fine["p"][0])), 6)
        doc.value(prefix + "Tmin", fine["t_min"], 6)
        e1 = radiated(st, H / 2.0, fine["t"])
        e2 = radiated(st, H, fine["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
        doc.check(rel < TOL_RAD, prefix + " radiation quadrature is inside its tolerance")
        doc.value(prefix + "RadErr", rel, 2, sci=True)
        doc.value(prefix + "Rad", e1, 4, sci=True)
        return e1

    rads = {}
    for pi in range(3):
        for si in range(2):
            rads[(pi, si)] = publish(pi, si)

    for si, tag in enumerate(shot_tags):
        st, fine, _c = runs[(0, si)]
        W_i = st.gamma0 * st.m * st.c ** 2 + st.V(*st.x0)
        doc.value("Tarr" + tag, W_i - st.m * st.c ** 2 - st.V(*fine["x"]), 6)

    best_r = math.inf
    best_si = 0
    best_idx = 0
    for si in range(2):
        fine_i = runs[(0, si)][1]
        if fine_i["r_close"] < best_r:
            best_r = fine_i["r_close"]
            best_si = si
            best_idx = fine_i["r_idx"]
    doc.value("rClose", best_r, 6)
    doc.value("clearClose", best_r - slow.R, 6)
    doc.text("rCloseShot", ("slower", "faster")[best_si])
    doc.text("rCloseCharge", ("lower", "upper")[best_idx])
    doc.value("sampleArc", slow.c * H / 8.0, 3)
    doc.check(best_r > slow.R + slow.a, "both reference ions stay clear of both spheres")
    rad_frac_max = 0.0
    for pi in range(3):
        for si in range(2):
            st_i = runs[(pi, si)][0]
            rad_frac_max = max(rad_frac_max, rads[(pi, si)] / st_i.T0)
    doc.value("radFracMax", rad_frac_max, 2, sci=True)
    doc.value("arrowX", slow.x0[0] + 1.6 * slow.d[0], 4)
    doc.value("arrowY", slow.x0[1] + 1.6 * slow.d[1], 4)

    def put(pi, si):
        _st, fine, _c = runs[(pi, si)]
        g = doc.flight(pi, flight=si)
        ge = g["end"]
        tag = names[pi] + shot_tags[si]
        title = {"Ref": "ref.", "Lo": "lower", "Hi": "upper"}[names[pi]]
        title += " " + {"Slow": "slow", "Fast": "fast"}[shot_tags[si]]
        doc.compare(tag + "Outcome", title + " outcome", fine["outcome"], g["outcome"])
        doc.compare(tag + "T", title + " $t$", fine["t"], ge["t"], tol)
        doc.compare(tag + "X", title + " $x$", fine["x"][0], ge["x"][0], tol)
        doc.compare(tag + "Y", title + " $y$", fine["x"][1], ge["x"][1], tol)
        doc.compare(tag + "Px", title + " $p_x$", fine["p"][0], ge["p"][0], tol)
        doc.compare(tag + "Py", title + " $p_y$", fine["p"][1], ge["p"][1], tol)
        doc.compare(tag + "Rad", title + " $E$ (rel.)", rads[(pi, si)],
                    g["radiated_energy"], TOL_RAD, rel=True)

    for pi in range(3):
        for si in range(2):
            put(pi, si)
