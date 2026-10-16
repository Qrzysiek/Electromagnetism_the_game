"""Companion of docs/levels/soft_landing/doc.tex (Level 13, Soft landing).

Three fixed centres, so no closed orbit (slingshot). The flight is the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4); the face crossing is the cubic
Hermite interpolant of that step (Around the wall, sec:method). An entry counts only
when the kinetic energy (gamma - 1) m c^2 of the entry momentum lies in the file's
window, inclusive at both ends. That energy is T0 + V(x0) - V(xe) when the energy is
conserved (Einzel lens, eq:kinetic). The truncation error is estimated by halving the
step, and the tolerance is set from that estimate before any comparison with the game.
"""
import math


# Reported step. Halving is the coarser run below: the error of the reported step is
# |y_h - y_{2h}|/15 (Foundations, the numerics section). The tolerance is ten times the
# largest such estimate and at least 1e-9, ten times the game's verified accuracy of
# about 1e-10. Computed in compute() before doc.compare, not from the game's numbers.
H = 0.005
# Simpson on the radiated power. 1e-4 still separates a wrong formula from the game's
# quadrature of the same power. Set beforehand, and used for every flight.
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
        kinetic = (shot["detector"].get("acceptance") or {}).get("kinetic")
        self.t_lo = float(kinetic[0])
        self.t_hi = float(kinetic[1])

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

    def kinetic_p(self, px, py):
        """(gamma - 1) m c^2 from the momentum. This is the energy the window tests."""
        p2 = px * px + py * py
        g = math.sqrt(1.0 + p2 / (self.m * self.m * self.c * self.c))
        return (g - 1.0) * self.m * self.c * self.c

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


def hermite_scalar(s, f, s2, f2, h, u):
    u2 = u * u
    u3 = u2 * u
    h00 = 2.0 * u3 - 3.0 * u2 + 1.0
    h10 = u3 - 2.0 * u2 + u
    h01 = -2.0 * u3 + 3.0 * u2
    h11 = u3 - u2
    return h00 * s + h10 * h * f + h01 * s2 + h11 * h * f2


def hermite(s, f, s2, f2, h, u):
    return (hermite_scalar(s[0], f[0], s2[0], f2[0], h, u),
            hermite_scalar(s[1], f[1], s2[1], f2[1], h, u))


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
    r_close = math.inf
    r_at = (0.0, s[0], s[1])
    r_idx = 0
    px_min = s[2]
    px_min_t = 0.0
    y_dev = 0.0
    y_max = s[1]
    y_min = s[1]
    y0 = st.x0[1]
    crossed = {}
    for idx, (_, cx, cy) in enumerate(st.kappas):
        rr = math.hypot(s[0] - cx, s[1] - cy)
        if rr < r_close:
            r_close = rr
            r_at = (0.0, s[0], s[1])
            r_idx = idx

    def absorb(u, s0, f0, s1, f1, t0):
        nonlocal t_min, t_min_at, r_close, r_at, r_idx, px_min, px_min_t, y_dev, y_max, y_min
        p = hermite(s0, f0, s1, f1, h, u)
        px = hermite_scalar(s0[2], f0[2], s1[2], f1[2], h, u)
        tk = (w0 - st.m * st.c * st.c) - st.V(p[0], p[1])
        tu = t0 + u * h
        if tk < t_min:
            t_min = tk
            t_min_at = (tu, p[0], p[1])
        if px < px_min:
            px_min = px
            px_min_t = tu
        dy = abs(p[1] - y0)
        if dy > y_dev:
            y_dev = dy
        if p[1] > y_max:
            y_max = p[1]
        if p[1] < y_min:
            y_min = p[1]
        for j, (_, cx, cy) in enumerate(st.kappas):
            rr = math.hypot(p[0] - cx, p[1] - cy)
            if rr < r_close:
                r_close = rr
                r_at = (tu, p[0], p[1])
                r_idx = j

    def note_cross(s0, f0, s1, f1, t0, u_max):
        """First crossing of each vertical detector edge, on this step, up to u_max."""
        x_end = hermite(s0, f0, s1, f1, h, u_max)[0]
        for xt in (st.dmin[0], st.dmax[0]):
            if xt in crossed:
                continue
            gap0 = s0[0] - xt
            gap1 = x_end - xt
            if gap0 * gap1 > 0.0:
                continue
            lo, hi_u = 0.0, u_max
            if gap0 == 0.0:
                crossed[xt] = (t0, s0[1])
                continue
            for _ in range(60):
                mid = 0.5 * (lo + hi_u)
                if gap0 * (hermite(s0, f0, s1, f1, h, mid)[0] - xt) <= 0.0:
                    hi_u = mid
                else:
                    lo = mid
            p = hermite(s0, f0, s1, f1, h, hi_u)
            crossed[xt] = (t0 + hi_u * h, p[1])

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
            for k in range(0, 9):
                u = k / 8.0
                if u <= hi:
                    absorb(u, s, f, s2, f2, t)
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            face, residual = "", 0.0
            if box is not None:
                pe, face, residual = assign_face(pe, box[0], box[1])
            if se[2] < px_min:
                px_min = se[2]
                px_min_t = te
            dy = abs(pe[1] - y0)
            if dy > y_dev:
                y_dev = dy
            if pe[1] > y_max:
                y_max = pe[1]
            if pe[1] < y_min:
                y_min = pe[1]
            note_cross(s, f, s2, f2, t, hi)
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n + 1,
                    "face": face, "residual": residual, "u": hi, "t_before": t,
                    "before": (s[0], s[1]), "after": (s2[0], s2[1]),
                    "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                    "t_min": t_min, "t_min_at": t_min_at,
                    "r_close": r_close, "r_at": r_at, "r_idx": r_idx,
                    "px_min": px_min, "px_min_t": px_min_t, "y_dev": y_dev,
                    "y_max": y_max, "y_min": y_min, "crossed": crossed}
        for k in range(0, 9):
            absorb(k / 8.0, s, f, s2, f2, t)
        note_cross(s, f, s2, f2, t, 1.0)
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
            "px_min": px_min, "px_min_t": px_min_t, "y_dev": y_dev,
            "y_max": y_max, "y_min": y_min, "crossed": crossed}


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


def apply_window(st, result):
    """The window is closed at both ends: T == Tmin or T == Tmax still counts.
    A miss is Rejected and the reported state stays the entry. A flight that never
    enters the box is not tested. T is (gamma - 1) m c^2 of the entry momentum."""
    px, py = result["p"]
    result["Tmom"] = st.kinetic_p(px, py)
    result["Tarr"] = None
    result["margin"] = None
    if result["outcome"] != "Arrived":
        return result
    Tarr = result["Tmom"]
    result["Tarr"] = Tarr
    span = max(abs(st.t_hi), 1e-300)
    result["margin"] = min(Tarr - st.t_lo, st.t_hi - Tarr) / span
    if Tarr < st.t_lo or Tarr > st.t_hi:
        result["outcome"] = "Rejected"
    return result


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
    names = ("Ref", "Pair", "Half")
    charge_names = ("lower electrode", "upper electrode", "added charge")
    runs = {}
    estimate = 0.0
    for pi, place in enumerate(places):
        st = Setup(lvl, place, 0)
        table_every = 5.0 if pi == 0 else 0.0
        fine, coarse = pair(st, H, path_every=40, table_every=table_every)
        apply_window(st, fine)
        apply_window(st, coarse)
        doc.check(fine["outcome"] == coarse["outcome"],
                  f"{names[pi]} has the same outcome at h and at 2h")
        doc.check(fine["residual"] < 1e-9, f"{names[pi]} crosses a face cleanly")
        estimate = max(estimate, step_gap(fine, coarse))
        runs[pi] = (st, fine, coarse)
    # Ten times the step-halving estimate, and at least the exact-result floor.
    # The arrival kinetic energy uses this same tolerance: it is computed from p.
    tol = max(1e-9, 10.0 * estimate)

    ref = runs[0][0]
    doc.value("q", ref.q, 3, sci=True)
    doc.value("m", ref.m, 3)
    doc.value("c", ref.c, 3)
    doc.value("Tzero", ref.T0, 3)
    doc.value("R", ref.R, 3)
    doc.value("aRad", ref.a, 3)
    doc.value("tmax", ref.tmax, 3)
    doc.value("xLaunch", ref.x0[0], 3)
    doc.value("yLaunch", ref.x0[1], 3)
    doc.value("gammaZero", ref.gamma0, 8)
    doc.value("pZero", ref.p0, 12)
    doc.value("vZero", ref.p0 / (ref.gamma0 * ref.m), 8)
    doc.value("angLaunch", ref.angle, 3)
    doc.value("dirGap", abs(ref.dir_len - 1.0), 2, sci=True)
    doc.value("xWorld", ref.wmax[0], 3)
    doc.value("yWorld", ref.wmax[1], 3)
    doc.value("xDlo", ref.dmin[0], 3)
    doc.value("xDhi", ref.dmax[0], 3)
    doc.value("yDlo", ref.dmin[1], 3)
    doc.value("yDhi", ref.dmax[1], 3)
    doc.value("Tlo", ref.t_lo, 4)
    doc.value("Thi", ref.t_hi, 4)
    doc.value("riseMin", ref.T0 - ref.t_hi, 4)
    doc.value("riseMax", ref.T0 - ref.t_lo, 4)
    doc.value("Tquarter", ref.T0 / 4.0, 6)
    doc.value("quarterGap", ref.T0 / 4.0 - ref.t_hi, 4)
    doc.value("H", H, 4)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("errMax", estimate, 2, sci=True)
    doc.value("maxCharges", lvl["limits"]["max_charges"])
    doc.value("nSample", 64)
    doc.value("nBisect", 80)
    region = lvl["limits"]["region"]
    doc.value("xRegLo", region["min"][0], 3)
    doc.value("xRegHi", region["max"][0], 3)
    doc.value("yRegLo", region["min"][1], 3)
    doc.value("yRegHi", region["max"][1], 3)
    (x_lo, y_lo, q_lo), (x_hi, y_hi, q_hi) = ref.charges[0], ref.charges[1]
    doc.value("Qlo", q_lo, 3, sci=True)
    doc.value("Qhi", q_hi, 3, sci=True)
    doc.value("kLo", ref.q * q_lo, 4)
    doc.value("kHi", ref.q * q_hi, 4)
    doc.value("xLo", x_lo, 3)
    doc.value("yLo", y_lo, 3)
    doc.value("xHi", x_hi, 3)
    doc.value("yHi", y_hi, 3)
    doc.value("yMid", 0.5 * (y_lo + y_hi), 3)
    x_r, y_r, q_r = ref.charges[2]
    doc.value("Qref", q_r, 3, sci=True)
    doc.value("kRef", ref.q * q_r, 4)
    doc.value("xRef", x_r, 3)
    doc.value("yRef", y_r, 3)
    half = runs[2][0]
    x_h, y_h, q_h = half.charges[2]
    doc.value("Qhalf", q_h, 3, sci=True)
    doc.value("kHalf", half.q * q_h, 4)
    doc.value("xHalf", x_h, 3)
    doc.value("yHalf", y_h, 3)
    doc.value("nRef", len(ref.charges))
    doc.value("nPair", len(runs[1][0].charges))
    doc.value("nHalf", len(half.charges))
    doc.value("sampleArc", ref.c * H / 8.0, 4)
    doc.value("arrowX", ref.x0[0] + 1.6 * ref.d[0], 4)
    doc.value("arrowY", ref.x0[1] + 1.6 * ref.d[1], 4)

    rads = {}
    rad_errs = {}

    def publish(pi):
        st, fine, coarse = runs[pi]
        prefix = names[pi]
        doc.text(prefix + "Outcome", fine["outcome"])
        doc.text(prefix + "Face", fine["face"])
        doc.value(prefix + "Steps", fine["steps"])
        doc.value(prefix + "T", fine["t"], 8)
        doc.value(prefix + "X", fine["x"][0], 8)
        doc.value(prefix + "Y", fine["x"][1], 8)
        doc.value(prefix + "Px", fine["p"][0], 8)
        doc.value(prefix + "Py", fine["p"][1], 8)
        doc.value(prefix + "Ang", math.degrees(math.atan2(fine["p"][1], fine["p"][0])), 6)
        doc.value(prefix + "Err", step_gap(fine, coarse), 2, sci=True)
        doc.value(prefix + "Wdev", fine["wdev"], 2, sci=True)
        doc.value(prefix + "FaceRes", fine["residual"], 2, sci=True)
        doc.value(prefix + "U", fine["u"], 6)
        doc.value(prefix + "Tbefore", fine["t_before"], 6)
        doc.value(prefix + "Xbefore", fine["before"][0], 6)
        doc.value(prefix + "Ybefore", fine["before"][1], 6)
        doc.value(prefix + "Xafter", fine["after"][0], 6)
        doc.value(prefix + "Yafter", fine["after"][1], 6)
        doc.value(prefix + "Tmin", fine["t_min"], 6)
        doc.value(prefix + "TTmin", fine["t_min_at"][0], 6)
        doc.value(prefix + "XTmin", fine["t_min_at"][1], 6)
        doc.value(prefix + "YTmin", fine["t_min_at"][2], 6)
        doc.value(prefix + "PxMin", fine["px_min"], 8)
        doc.value(prefix + "TPxMin", fine["px_min_t"], 6)
        doc.value(prefix + "Ydev", fine["y_dev"], 6)
        doc.value(prefix + "Ymax", fine["y_max"], 6)
        doc.value(prefix + "Ymin", fine["y_min"], 6)
        for edge, xt in (("Dlo", st.dmin[0]), ("Dhi", st.dmax[0])):
            if xt in fine["crossed"]:
                tc, yc = fine["crossed"][xt]
                doc.value(prefix + "T" + edge, tc, 6)
                doc.value(prefix + "Y" + edge, yc, 6)
                doc.value(prefix + "Above" + edge, yc - st.dmax[1], 6)
        doc.value(prefix + "Rclose", fine["r_close"], 6)
        doc.value(prefix + "Clear", fine["r_close"] - st.R, 6)
        doc.text(prefix + "Rname", charge_names[fine["r_idx"]])
        V0 = st.V(*st.x0)
        Ve = st.V(*fine["x"])
        TarrV = st.T0 + V0 - Ve
        f0 = st.rhs(st.launch_state())
        doc.value(prefix + "Vzero", V0, 8)
        doc.value(prefix + "Vend", Ve, 8)
        doc.value(prefix + "Rise", Ve - V0, 8)
        doc.value(prefix + "W", st.gamma0 * st.m * st.c ** 2 + V0, 10)
        doc.value(prefix + "Wbind", st.T0 + V0, 8)
        doc.value(prefix + "Fx", f0[2], 6)
        doc.value(prefix + "Fy", f0[3], 6)
        doc.value(prefix + "Tmom", fine["Tmom"], 8)
        doc.value(prefix + "TarrV", TarrV, 8)
        doc.value(prefix + "Tgap", fine["Tmom"] - TarrV, 2, sci=True)
        if fine["Tarr"] is not None:
            doc.value(prefix + "Tarr", fine["Tarr"], 8)
            doc.value(prefix + "Margin", fine["margin"], 6)
            doc.value(prefix + "RoomHi", st.t_hi - fine["Tarr"], 6)
            doc.value(prefix + "RoomLo", fine["Tarr"] - st.t_lo, 6)
        e1 = radiated(st, H / 2.0, fine["t"])
        e2 = radiated(st, H, fine["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
        doc.check(rel < TOL_RAD, prefix + " radiation quadrature is inside its tolerance")
        rads[pi] = e1
        rad_errs[pi] = rel
        doc.value(prefix + "Rad", e1, 4, sci=True)
        doc.value(prefix + "RadErr", rel, 2, sci=True)
        doc.value(prefix + "RadFrac", e1 / st.T0, 2, sci=True)
        doc.data("path" + prefix.lower(), fine["path"])
        return fine

    for pi in range(3):
        publish(pi)

    fine = runs[0][1]
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$\Delta W$"],
              list(fine["rows"]),
              digits=[5, 5, 5, 5, 5, 2],
              sci=[False, False, False, False, False, True])
    rad_frac_max = max(rads[pi] / runs[pi][0].T0 for pi in range(3))
    doc.value("radFracMax", rad_frac_max, 2, sci=True)
    doc.check(min(fine["r_close"] for _st, fine, _c in runs.values()) > ref.R + ref.a,
              "every flight stays clear of the spheres")

    titles = {"Ref": "ref.", "Pair": "pair", "Half": "half"}

    def put(pi):
        st, fine, _c = runs[pi]
        g = doc.flight(pi)
        ge = g["end"]
        tag = names[pi]
        title = titles[tag]
        doc.compare(tag + "Outcome", title + " outcome", fine["outcome"], g["outcome"])
        doc.compare(tag + "T", title + " $t$", fine["t"], ge["t"], tol)
        doc.compare(tag + "X", title + " $x$", fine["x"][0], ge["x"][0], tol)
        doc.compare(tag + "Y", title + " $y$", fine["x"][1], ge["x"][1], tol)
        doc.compare(tag + "Px", title + " $p_x$", fine["p"][0], ge["p"][0], tol)
        doc.compare(tag + "Py", title + " $p_y$", fine["p"][1], ge["p"][1], tol)
        if fine["Tarr"] is not None:
            doc.compare(tag + "Tarr", title + " $T$", fine["Tarr"],
                        st.kinetic_p(ge["p"][0], ge["p"][1]), tol)
        doc.compare(tag + "Rad", title + " $E$ (rel.)", rads[pi],
                    g["radiated_energy"], TOL_RAD, rel=True)

    for pi in range(3):
        put(pi)
