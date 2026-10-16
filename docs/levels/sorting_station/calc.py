"""Companion of docs/levels/sorting_station/doc.tex (Level 14, Sorting station).

Six charges, so no closed orbit (slingshot). Each ray is the classical Runge-Kutta
method at a fixed step (Foundations eq:rk4); a face crossing is the cubic Hermite
interpolant of that step (Around the wall, sec:method). The flight must enter the
gate before the detector counts (Two stages, eq:gate). The slower rays bend more,
which is the 1/T kick of Thomson (eq:impulse); the paths themselves are integrated.
The truncation error is estimated by halving the step, and the tolerance is set
from that estimate before any comparison with the game.
"""
import math


# Reported step. Halving is the coarser run: the error of the reported step is
# |y_h - y_{2h}|/15 (Foundations, the numerics section). The tolerance is ten times
# the largest such estimate and at least 1e-9. Computed in compute() before
# doc.compare, not from the game's numbers.
H = 0.005
# Simpson on the radiated power. 1e-4 separates a wrong formula from the game's
# quadrature of the same power. Set beforehand.
TOL_RAD = 1e-4
# The angles named in the level text, in shot order. The aim used is the file's
# direction; the gap is reported.
NOMINAL_DEG = (-8.0, 0.0, 8.0, -8.0, 0.0, 8.0)
TAGS = ("SlowDown", "SlowAxis", "SlowUp", "FastDown", "FastAxis", "FastUp")


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
        g = lvl["gates"][0]
        self.gmin = [v / s for v in g["min"][:2]]
        self.gmax = [v / s for v in g["max"][:2]]
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
        _x, _y, px, py = state
        rest = math.sqrt(self.m * self.m * self.c ** 4 + (px * px + py * py) * self.c * self.c)
        return rest + self.V(state[0], state[1])

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


def _in_box(p, lo, hi):
    return lo[0] <= p[0] <= hi[0] and lo[1] <= p[1] <= hi[1]


def fly(st, h, path_every=0, table_every=0.0):
    """Integrate until the first ending event. A gate is entered and the flight goes on.
    Board order of endings: collision, world, detector. A detector entry with the gate
    still ahead is SkippedGate(0)."""
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
    r_close = math.inf
    r_idx = 0
    y_max = s[1]
    y_min = s[1]
    px_min = s[2]
    gate_passed = _in_box((s[0], s[1]), st.gmin, st.gmax)
    gate_info = None
    for idx, (_, cx, cy) in enumerate(st.kappas):
        rr = math.hypot(s[0] - cx, s[1] - cy)
        if rr < r_close:
            r_close = rr
            r_idx = idx

    def absorb(u, s0, f0, s1, f1):
        nonlocal r_close, r_idx, y_max, y_min, px_min
        p = hermite(s0, f0, s1, f1, h, u)
        px = hermite_scalar(s0[2], f0[2], s1[2], f1[2], h, u)
        if p[1] > y_max:
            y_max = p[1]
        if p[1] < y_min:
            y_min = p[1]
        if px < px_min:
            px_min = px
        for j, (_, cx, cy) in enumerate(st.kappas):
            rr = math.hypot(p[0] - cx, p[1] - cy)
            if rr < r_close:
                r_close = rr
                r_idx = j

    def mark_gate(s0, f0, s1, f1, t0, u_hit):
        nonlocal gate_passed, gate_info
        if gate_passed:
            return
        pts = [hermite(s0, f0, s1, f1, h, u) for u in us]
        u = _first(lambda p: _in_box(p, st.gmin, st.gmax), pts, us)
        if u is None or (u_hit is not None and u >= u_hit):
            return
        lo = max(0.0, u - 1.0 / 64.0)
        hi = _bisect(lambda p: _in_box(p, st.gmin, st.gmax), s0, f0, s1, f1, h, lo, u)
        if u_hit is not None and hi >= u_hit:
            return
        te = t0 + hi * h
        se = rk4_step(st, s0, f0, te - t0)
        pe, face, residual = assign_face(hermite(s0, f0, s1, f1, h, hi), st.gmin, st.gmax)
        gate_passed = True
        gate_info = {"t": te, "x": pe, "p": (se[2], se[3]), "face": face,
                     "residual": residual, "u": hi}

    def pack(outcome, te, pe, se, n_steps, face, residual, u, t_before, before, after):
        if outcome == "Arrived" and not gate_passed:
            outcome = "SkippedGate(0)"
        return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n_steps,
                "face": face, "residual": residual, "u": u, "t_before": t_before,
                "before": before, "after": after, "rows": rows, "path": path,
                "wdev": wdev, "w0": w0, "r_close": r_close, "r_idx": r_idx,
                "y_max": y_max, "y_min": y_min, "px_min": px_min,
                "gate_passed": gate_passed, "gate": gate_info}

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
        overlaps_world = (box_lo[0] < st.wmin[0] or box_hi[0] > st.wmax[0]
                          or box_lo[1] < st.wmin[1] or box_hi[1] > st.wmax[1])
        if overlaps_world:
            u = _first(lambda p: not _in_box(p, st.wmin, st.wmax), samples(), us)
            if u is not None:
                hits.append((u, 1, "LeftBounds", None))
        overlaps_det = not (box_hi[0] < st.dmin[0] or box_lo[0] > st.dmax[0]
                            or box_hi[1] < st.dmin[1] or box_lo[1] > st.dmax[1])
        if overlaps_det:
            u = _first(lambda p: _in_box(p, st.dmin, st.dmax), samples(), us)
            if u is not None:
                hits.append((u, 2, "Arrived", None))
        u_hit = min(hits)[0] if hits else None
        overlaps_gate = not (box_hi[0] < st.gmin[0] or box_lo[0] > st.gmax[0]
                             or box_hi[1] < st.gmin[1] or box_lo[1] > st.gmax[1])
        if overlaps_gate:
            mark_gate(s, f, s2, f2, t, u_hit)
        if hits:
            u_hit, _rank, outcome, idx = min(hits)
            if outcome == "Arrived":
                pred = lambda p: _in_box(p, st.dmin, st.dmax)
                box = (st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                pred = lambda p: not _in_box(p, st.wmin, st.wmax)
                box = (st.wmin, st.wmax)
            else:
                cx, cy = st.kappas[idx][1], st.kappas[idx][2]
                pred = lambda p, cx=cx, cy=cy: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit
                box = None
            hi = _bisect(pred, s, f, s2, f2, h, max(0.0, u_hit - 1.0 / 64.0), u_hit)
            for k in range(0, 9):
                if k / 8.0 <= hi:
                    absorb(k / 8.0, s, f, s2, f2)
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            face, residual = "", 0.0
            if box is not None:
                pe, face, residual = assign_face(pe, box[0], box[1])
            if pe[1] > y_max:
                y_max = pe[1]
            if pe[1] < y_min:
                y_min = pe[1]
            if se[2] < px_min:
                px_min = se[2]
            path.append(pe)
            return pack(outcome, te, pe, se, n + 1, face, residual, hi, t, (s[0], s[1]), (s2[0], s2[1]))
        for k in range(0, 9):
            absorb(k / 8.0, s, f, s2, f2)
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
    return pack("Timeout", t, (s[0], s[1]), s, n, "", 0.0, 1.0, t, (s[0], s[1]), (s[0], s[1]))


def radiated(st, h, t_end):
    m = st.m
    c = st.c
    q = st.q

    def power(state):
        _x, _y, px, py = state
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
    places = [lvl["reference_solution"], doc.result(1)["placement"]]
    names = ("Ref", "Lens")
    runs = {}
    estimate = 0.0
    for pi, place in enumerate(places):
        for si in range(6):
            st = Setup(lvl, place, si)
            table_every = 5.0 if pi == 0 and si == 1 else 0.0
            path_every = 40 if pi == 0 else 0
            fine, coarse = pair(st, H, path_every=path_every, table_every=table_every)
            doc.check(fine["outcome"] == coarse["outcome"],
                      f"{names[pi]} shot {si} has the same outcome at h and at 2h")
            doc.check(fine["residual"] < 1e-9, f"{names[pi]} shot {si} crosses a face cleanly")
            estimate = max(estimate, step_gap(fine, coarse))
            runs[(pi, si)] = (st, fine, coarse)
    tol = max(1e-9, 10.0 * estimate)

    ref0 = runs[(0, 1)][0]
    doc.value("q", ref0.q, 3, sci=True)
    doc.value("m", ref0.m, 3)
    doc.value("c", ref0.c, 3)
    doc.value("R", ref0.R, 3)
    doc.value("aRad", ref0.a, 3)
    doc.value("tmax", ref0.tmax, 3)
    doc.value("xLaunch", ref0.x0[0], 3)
    doc.value("yLaunch", ref0.x0[1], 3)
    doc.value("xWorld", ref0.wmax[0], 3)
    doc.value("yWorld", ref0.wmax[1], 3)
    doc.value("H", H, 4)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("errMax", estimate, 2, sci=True)
    doc.value("maxCharges", lvl["limits"]["max_charges"])
    doc.value("nSample", 64)
    doc.value("nBisect", 80)
    doc.value("nShots", 6)
    doc.value("nCharges", len(ref0.charges))
    doc.value("nLens", len(runs[(1, 0)][0].charges))
    region = lvl["limits"]["region"]
    doc.value("xRegLo", region["min"][0], 3)
    doc.value("xRegHi", region["max"][0], 3)
    doc.value("yRegLo", region["min"][1], 3)
    doc.value("yRegHi", region["max"][1], 3)
    doc.value("xGlo", ref0.gmin[0], 3)
    doc.value("xGhi", ref0.gmax[0], 3)
    doc.value("yGlo", ref0.gmin[1], 3)
    doc.value("yGhi", ref0.gmax[1], 3)
    slow = runs[(0, 0)][0]
    fast = runs[(0, 3)][0]
    doc.value("Tslow", slow.T0, 3)
    doc.value("Tfast", fast.T0, 3)
    doc.value("ratioT", fast.T0 / slow.T0, 4)
    doc.value("xSlowLo", slow.dmin[0], 3)
    doc.value("xSlowHi", slow.dmax[0], 3)
    doc.value("ySlowLo", slow.dmin[1], 3)
    doc.value("ySlowHi", slow.dmax[1], 3)
    doc.value("xFastLo", fast.dmin[0], 3)
    doc.value("xFastHi", fast.dmax[0], 3)
    doc.value("yFastLo", fast.dmin[1], 3)
    doc.value("yFastHi", fast.dmax[1], 3)
    doc.value("gammaSlow", slow.gamma0, 8)
    doc.value("pSlow", slow.p0, 12)
    doc.value("vSlow", slow.p0 / (slow.gamma0 * slow.m), 8)
    doc.value("gammaFast", fast.gamma0, 8)
    doc.value("pFast", fast.p0, 12)
    doc.value("vFast", fast.p0 / (fast.gamma0 * fast.m), 8)
    doc.value("sampleArc", ref0.c * H / 8.0, 4)
    for i, (cx, cy, Q) in enumerate(ref0.charges):
        doc.value("Q" + str(i + 1), Q, 3, sci=True)
        doc.value("k" + str(i + 1), ref0.q * Q, 4)
        doc.value("x" + str(i + 1), cx, 3)
        doc.value("y" + str(i + 1), cy, 3)
    doc.value("arrowX", ref0.x0[0] + 1.6 * ref0.d[0], 4)
    doc.value("arrowY", ref0.x0[1] + 1.6 * ref0.d[1], 4)

    same_slow = all(runs[(0, si)][0].dmin == slow.dmin and runs[(0, si)][0].dmax == slow.dmax
                    for si in range(3))
    same_fast = all(runs[(0, si)][0].dmin == fast.dmin and runs[(0, si)][0].dmax == fast.dmax
                    for si in range(3, 6))
    doc.check(same_slow, "the three slow shots share one detector")
    doc.check(same_fast, "the three fast shots share one detector")

    rads = {}

    def publish(pi, si):
        st, fine, coarse = runs[(pi, si)]
        prefix = names[pi] + TAGS[si]
        doc.text(prefix + "Outcome", fine["outcome"])
        doc.text(prefix + "Face", fine["face"])
        doc.value(prefix + "T", fine["t"], 8)
        doc.value(prefix + "X", fine["x"][0], 8)
        doc.value(prefix + "Y", fine["x"][1], 8)
        doc.value(prefix + "Px", fine["p"][0], 8)
        doc.value(prefix + "Py", fine["p"][1], 8)
        doc.value(prefix + "Ang", math.degrees(math.atan2(fine["p"][1], fine["p"][0])), 6)
        doc.value(prefix + "Err", step_gap(fine, coarse), 2, sci=True)
        doc.value(prefix + "Wdev", fine["wdev"], 2, sci=True)
        doc.value(prefix + "FaceRes", fine["residual"], 2, sci=True)
        doc.value(prefix + "Steps", fine["steps"])
        doc.value(prefix + "Ymax", fine["y_max"], 6)
        doc.value(prefix + "Ymin", fine["y_min"], 6)
        doc.value(prefix + "PxMin", fine["px_min"], 6)
        doc.value(prefix + "Rclose", fine["r_close"], 6)
        doc.value(prefix + "Clear", fine["r_close"] - st.R, 6)
        doc.text(prefix + "Passed", "yes" if fine["gate_passed"] else "no")
        if fine["gate"] is not None:
            g = fine["gate"]
            doc.text(prefix + "Gface", g["face"])
            doc.value(prefix + "Gt", g["t"], 6)
            doc.value(prefix + "Gx", g["x"][0], 6)
            doc.value(prefix + "Gy", g["x"][1], 6)
            doc.value(prefix + "Gpx", g["p"][0], 6)
            doc.value(prefix + "Gpy", g["p"][1], 6)
            doc.value(prefix + "Gang", math.degrees(math.atan2(g["p"][1], g["p"][0])), 6)
            doc.value(prefix + "Gres", g["residual"], 2, sci=True)
        if pi == 0:
            doc.value("ang" + TAGS[si], st.angle, 6)
            doc.value("angGap" + TAGS[si], st.angle - NOMINAL_DEG[si], 2, sci=True)
            doc.value("dirGap" + TAGS[si], abs(st.dir_len - 1.0), 2, sci=True)
            V0 = st.V(*st.x0)
            doc.value(prefix + "Vzero", V0, 6)
            doc.value(prefix + "Wbind", st.T0 + V0, 6)
        e1 = radiated(st, H / 2.0, fine["t"])
        e2 = radiated(st, H, fine["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
        doc.check(rel < TOL_RAD, prefix + " radiation quadrature is inside its tolerance")
        rads[(pi, si)] = e1
        doc.value(prefix + "Rad", e1, 4, sci=True)
        doc.value(prefix + "RadErr", rel, 2, sci=True)
        doc.value(prefix + "RadFrac", e1 / st.T0, 2, sci=True)
        if pi == 0:
            doc.data("path" + TAGS[si].lower(), fine["path"])

    for pi in range(2):
        for si in range(6):
            publish(pi, si)

    axis = runs[(0, 1)]
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$\Delta W$"],
              list(axis[1]["rows"]),
              digits=[5, 5, 5, 5, 5, 2],
              sci=[False, False, False, False, False, True])
    # The straight-line kick scales as 1/T. Compare the turn of the two axis rays
    # from the launch direction to the gate, as a check of that scaling, not as the flight.
    slow_turn = runs[(0, 1)][1]["gate"]
    fast_turn = runs[(0, 4)][1]["gate"]
    if slow_turn is not None and fast_turn is not None:
        d_slow = math.degrees(math.atan2(slow_turn["p"][1], slow_turn["p"][0])) - runs[(0, 1)][0].angle
        d_fast = math.degrees(math.atan2(fast_turn["p"][1], fast_turn["p"][0])) - runs[(0, 4)][0].angle
        doc.value("turnSlow", d_slow, 6)
        doc.value("turnFast", d_fast, 6)
        doc.value("turnRatio", d_slow / d_fast if d_fast else 0.0, 4)
        doc.value("kickSlow", d_slow * runs[(0, 1)][0].T0, 6)
        doc.value("kickFast", d_fast * runs[(0, 4)][0].T0, 6)
    clears = [runs[(pi, si)][1]["r_close"] for pi in range(2) for si in range(6)]
    doc.value("rCloseMin", min(clears), 6)
    doc.value("clearMin", min(clears) - ref0.R, 6)
    fracs = [rads[(pi, si)] / runs[(pi, si)][0].T0 for pi in range(2) for si in range(6)]
    doc.value("radFracMax", max(fracs), 2, sci=True)
    doc.check(min(clears) > ref0.R + ref0.a, "every flight stays clear of the spheres")
    doc.check(runs[(1, 0)][0].charges == ref0.charges[:2],
              "the lens placement is the first two reference charges")

    titles = ("slow $-$", "slow $0$", "slow $+$", "fast $-$", "fast $0$", "fast $+$")

    def put(pi, si, with_state, with_rad):
        st, fine, _c = runs[(pi, si)]
        g = doc.flight(pi, flight=si)
        ge = g["end"]
        tag = names[pi] + TAGS[si]
        title = ("L " if pi else "") + titles[si]
        doc.compare(tag + "Outcome", title + " outcome", fine["outcome"], g["outcome"])
        doc.compare(tag + "T", title + " $t$", fine["t"], ge["t"], tol)
        if with_state:
            doc.compare(tag + "X", title + " $x$", fine["x"][0], ge["x"][0], tol)
            doc.compare(tag + "Y", title + " $y$", fine["x"][1], ge["x"][1], tol)
        if with_rad:
            doc.compare(tag + "Rad", title + " $E$ (rel.)", rads[(pi, si)],
                        g["radiated_energy"], TOL_RAD, rel=True)

    for si in range(6):
        put(0, si, True, si in (1, 4))
    for si in range(6):
        put(1, si, False, False)
