"""Companion of docs/levels/around_the_wall/doc.tex (Level 8, Around the wall).

Eleven fixed centres, so there is no closed orbit (slingshot). The flight is the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4), the crossing of a face is the cubic
Hermite interpolant of that step (slingshot eq:hermite), and the detector's direction test is
the cone of Injection. The truncation error is estimated by halving the step, and the
tolerance is set from that estimate before any comparison with the game.
"""
import math


# The step that is reported. Halving it estimates the error as |y_h - y_{h/2}|/15
# (Foundations, the numerics section). The tolerance is ten times that estimate and at
# least 1e-9, which is ten times the game's verified accuracy of about 1e-10. It is
# computed in compute() from these flights, and it is not taken from the game's numbers.
H = 0.001
# Simpson's rule on the radiated power, against the game's quadrature of the same power
# (about 3e-9 relative on an exact test). 1e-4 still separates a wrong formula. Set
# before comparing, and not widened to fit a difference.
TOL_RAD = 1e-4


class Setup:
    """The level's numbers and the charges of one placement, wall first, then the player."""

    def __init__(self, lvl, placement):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.a = shot["particle"]["radius"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.sub = s
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(d[0], d[1])
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        acc = shot["detector"]["acceptance"]["direction"]
        self.axis_deg = acc[0]
        self.half_deg = acc[1]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        self.kappas = tuple((self.q * Q, cx, cy) for cx, cy, Q in self.charges)
        self.n_wall = sum(1 for e in lvl["elements"] if e["kind"] == "charge")
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)

    def V(self, x, y):
        """Potential energy q*phi (Foundations eq:coulomb), summed in file order."""
        total = 0.0
        for k, cx, cy in self.kappas:
            total += k / math.hypot(x - cx, y - cy)
        return total

    def rhs(self, state):
        """Right-hand side of the planar equation of motion, B = 0."""
        x, y, px, py = state
        m = self.m
        c = self.c
        p2 = px * px + py * py
        g = math.sqrt(1.0 + p2 / (m * m * c * c))
        inv = 1.0 / (g * m)
        fx = 0.0
        fy = 0.0
        for k, cx, cy in self.kappas:
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            inv_r3 = k / (r2 * math.sqrt(r2))
            fx += dx * inv_r3
            fy += dy * inv_r3
        return (px * inv, py * inv, fx, fy)

    def energy(self, state):
        """W = gamma m c^2 + V."""
        x, y, px, py = state
        rest = math.sqrt(self.m * self.m * self.c ** 4 + (px * px + py * py) * self.c * self.c)
        return rest + self.V(x, y)

    def launch_state(self):
        return (self.x0[0], self.x0[1], self.p0 * self.d[0], self.p0 * self.d[1])


def rk4_step(st, s, f, h):
    """One step of the classical Runge-Kutta method (Foundations eq:rk4)."""
    x, y, px, py = s
    k1x, k1y, k1px, k1py = f
    k2x, k2y, k2px, k2py = st.rhs((x + 0.5 * h * k1x, y + 0.5 * h * k1y,
                                   px + 0.5 * h * k1px, py + 0.5 * h * k1py))
    k3x, k3y, k3px, k3py = st.rhs((x + 0.5 * h * k2x, y + 0.5 * h * k2y,
                                   px + 0.5 * h * k2px, py + 0.5 * h * k2py))
    k4x, k4y, k4px, k4py = st.rhs((x + h * k3x, y + h * k3y,
                                   px + h * k3px, py + h * k3py))
    return (x + h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
            y + h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y),
            px + h / 6.0 * (k1px + 2.0 * k2px + 2.0 * k3px + k4px),
            py + h / 6.0 * (k1py + 2.0 * k2py + 2.0 * k3py + k4py))


def hermite(s, f, s2, f2, h, u):
    """Cubic Hermite interpolant of the position on a step, at fraction u."""
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
    """The face the free coordinate selects, with that coordinate set to the face value.

    The residual is the interpolant's distance to the face before the assignment. A residual
    below 1e-6 is the crossing of that face; anything larger is left as the interpolant gave it.
    """
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


def beta_of(p, axis_deg):
    """The cone angle (injection eq:cone), from the momentum. v is parallel to p."""
    a = math.radians(axis_deg)
    ex, ey = math.cos(a), math.sin(a)
    cross = abs(p[0] * ey - p[1] * ex)
    dot = p[0] * ex + p[1] * ey
    return math.atan2(cross, dot)


def fly(st, h, path_every=0, table_every=0.0):
    """Integrate from the launch until the first event.

    Events in one step are ordered as the board orders them: a collision, then the world's
    edge, then the detector. The position at the event is the Hermite interpolant; the
    momentum is one Runge-Kutta step of that length. A box event whose free coordinate lies
    in a face has that face's coordinate set to the face value.
    """
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
    # Closest approach to each charge, refined on the step that improves the endpoint minimum.
    closest = [math.inf] * len(st.kappas)
    close_xy = [(0.0, 0.0, 0.0)] * len(st.kappas)
    y_min = s[1]
    y_min_at = (0.0, s[0], s[1])
    while t <= st.tmax:
        s2 = rk4_step(st, s, f, h)
        f2 = st.rhs(s2)
        # The interpolant stays within c*h of the box of the two ends: the speed is less
        # than c, and the Hermite velocity terms add less than that. A boundary the enlarged
        # box misses cannot meet this step.
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
        leaves = (box_lo[0] < st.wmin[0] or box_hi[0] > st.wmax[0]
                  or box_lo[1] < st.wmin[1] or box_hi[1] > st.wmax[1])
        if leaves:
            u = _first(lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                      and st.wmin[1] <= p[1] <= st.wmax[1]),
                       samples(), us)
            if u is not None:
                hits.append((u, 1, "LeftBounds", None))
        enters = not (box_hi[0] < st.dmin[0] or box_lo[0] > st.dmax[0]
                      or box_hi[1] < st.dmin[1] or box_lo[1] > st.dmax[1])
        if enters:
            u = _first(lambda p: st.dmin[0] <= p[0] <= st.dmax[0] and st.dmin[1] <= p[1] <= st.dmax[1],
                       samples(), us)
            if u is not None:
                hits.append((u, 2, "Arrived", None))
        u_stop = min(hits)[0] if hits else 1.0
        # Distance to a charge falls by less than c*h on a step, plus the interpolant's bulge.
        span = 2.0 * margin
        active = [idx for idx, (_, cx, cy) in enumerate(st.kappas)
                  if closest[idx] == math.inf
                  or (s[0] - cx) ** 2 + (s[1] - cy) ** 2 <= (closest[idx] + span) ** 2]
        if box_lo[1] < y_min or active:
            for u, p in zip(us, samples()):
                if u > u_stop:
                    break
                tu = t + u * h
                if p[1] < y_min:
                    y_min = p[1]
                    y_min_at = (tu, p[0], p[1])
                for idx in active:
                    _, cx, cy = st.kappas[idx]
                    r = math.hypot(p[0] - cx, p[1] - cy)
                    if r < closest[idx]:
                        closest[idx] = r
                        close_xy[idx] = (tu, p[0], p[1])
        if hits:
            u_hit, rank, outcome, idx = min(hits)
            if outcome == "Arrived":
                pred = lambda p: st.dmin[0] <= p[0] <= st.dmax[0] and st.dmin[1] <= p[1] <= st.dmax[1]
                box = (st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                pred = lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                       and st.wmin[1] <= p[1] <= st.wmax[1])
                box = (st.wmin, st.wmax)
            else:
                cx, cy = st.kappas[idx][1], st.kappas[idx][2]
                pred = lambda p, cx=cx, cy=cy: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit
                box = None
            lo = u_hit - 1.0 / 64.0
            hi = _bisect(pred, s, f, s2, f2, h, lo, u_hit)
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            face, residual = "", 0.0
            if box is not None:
                pe, face, residual = assign_face(pe, box[0], box[1])
            path.append(pe)
            return {
                "outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n + 1,
                "face": face, "residual": residual, "u": hi, "t_before": t,
                "before": (s[0], s[1]), "after": (s2[0], s2[1]),
                "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                "closest": closest, "close_xy": close_xy, "y_min": y_min_at,
            }
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
    return {
        "outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": (s[2], s[3]), "steps": n,
        "face": "", "residual": 0.0, "u": 1.0, "t_before": t,
        "before": (s[0], s[1]), "after": (s[0], s[1]),
        "rows": rows, "path": path, "wdev": wdev, "w0": w0,
        "closest": closest, "close_xy": close_xy, "y_min": y_min_at,
    }


def radiated(st, h, t_end):
    """Radiated energy along the RK4 path (Foundations eq:larmor, first_bend eq:acceleration).

    Simpson's rule on an even number of steps, then the trapezoidal rule on the remainder.
    """
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
        return (2.0 / 3.0) * q * q / (c * c * c) * g ** 6 * (ax * ax + ay * ay - vxa * vxa / (c * c))

    n = int(t_end / h)
    n -= n % 2
    s = st.launch_state()
    acc = power(s)
    total = acc
    # Simpson weights 1, 4, 2, 4, ..., 2, 4, 1 on n+1 samples (n even).
    for i in range(1, n + 1):
        s = rk4_step(st, s, st.rhs(s), h)
        w = 4.0 if i % 2 else 2.0
        if i == n:
            w = 1.0
        total += w * power(s)
    simpson = h / 3.0 * total
    rest = t_end - n * h
    s2 = rk4_step(st, s, st.rhs(s), rest) if rest > 0.0 else s
    return simpson + rest * 0.5 * (power(s) + power(s2))


def finish(st, run):
    """Apply the cone. Entering outside it rejects the flight; the entry state stands."""
    if run["outcome"] == "Arrived":
        b = beta_of(run["p"], st.axis_deg)
        run["beta"] = b
        if b > math.radians(st.half_deg):
            run["outcome"] = "Rejected"
    else:
        run["beta"] = None
    return run


def step_gap(fine, coarse):
    """|y_h - y_{h/2}|/15 for the reported (finer) state. The face coordinate is assigned, so
    its gap is the free coordinate's and the momentum's, not a residual of the face value."""
    gaps = [abs(fine["t"] - coarse["t"]) / 15.0,
            abs(fine["x"][0] - coarse["x"][0]) / 15.0,
            abs(fine["x"][1] - coarse["x"][1]) / 15.0,
            abs(fine["p"][0] - coarse["p"][0]) / 15.0,
            abs(fine["p"][1] - coarse["p"][1]) / 15.0]
    if fine["beta"] is not None and coarse["beta"] is not None:
        gaps.append(abs(fine["beta"] - coarse["beta"]) / 15.0)
    return max(gaps)


def compute(doc):
    lvl = doc.level
    ref_place = lvl["reference_solution"]
    up_place = doc.result(1)["placement"]
    at_place = doc.result(2)["placement"]
    ref = Setup(lvl, ref_place)
    up = Setup(lvl, up_place)
    at = Setup(lvl, at_place)

    ref_f = finish(ref, fly(ref, H, path_every=40, table_every=5.0))
    ref_c = finish(ref, fly(ref, 2.0 * H))
    up_f = finish(up, fly(up, H, path_every=40))
    up_c = finish(up, fly(up, 2.0 * H))
    at_f = finish(at, fly(at, H, path_every=40))
    at_c = finish(at, fly(at, 2.0 * H))
    doc.check(ref_f["outcome"] == ref_c["outcome"] == "Arrived",
              "the reference arrives at h and at h/2, inside the cone")
    doc.check(up_f["outcome"] == up_c["outcome"] == "LeftBounds",
              "the upper charge alone leaves the world at both steps")
    doc.check(at_f["outcome"] == at_c["outcome"] == "LeftBounds",
              "the attracting charge alone leaves the world at both steps")
    doc.check(up_f["p"][0] < 0.0 and up_f["x"][0] < ref.charges[0][0],
              "the upper charge turns the particle back before the wall")
    doc.check(at_f["p"][0] > 0.0 and at_f["x"][0] < ref.charges[0][0],
              "the attractor alone reaches the floor before the wall")
    doc.check(ref_f["residual"] < 1e-9 and up_f["residual"] < 1e-9 and at_f["residual"] < 1e-9,
              "each event is a clean face crossing")

    # Ten times the largest step-halving estimate, and at least 1e-9. Before any comparison.
    estimate = max(step_gap(ref_f, ref_c), step_gap(up_f, up_c), step_gap(at_f, at_c))
    tol = max(1e-9, 10.0 * estimate)

    # ---- Symbols -----------------------------------------------------------------------------
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
    v0 = ref.p0 / (ref.gamma0 * ref.m)
    doc.value("vZero", v0, 8)
    doc.value("xWorld", ref.wmax[0], 3)
    doc.value("yWorld", ref.wmax[1], 3)
    doc.value("nWall", ref.n_wall)
    wall = ref.charges[:ref.n_wall]
    doc.value("Qwall", wall[0][2], 3, sci=True)
    doc.value("kWall", ref.q * wall[0][2], 4)
    doc.value("xWall", wall[0][0], 3)
    doc.value("yWallLo", wall[0][1], 3)
    doc.value("yWallHi", wall[-1][1], 3)
    doc.value("nCharges", len(ref.charges))
    # The reference's two charges are the last two, in file order: the upper one, then the attractor.
    (x_up, y_up, q_up), (x_at, y_at, q_at) = ref.charges[ref.n_wall:]
    doc.value("Qup", q_up, 3, sci=True)
    doc.value("kUp", ref.q * q_up, 4)
    doc.value("xUp", x_up, 3)
    doc.value("yUp", y_up, 3)
    doc.value("Qat", q_at, 3, sci=True)
    doc.value("kAt", ref.q * q_at, 4)
    doc.value("xAt", x_at, 3)
    doc.value("yAt", y_at, 3)
    doc.value("xDlo", ref.dmin[0], 3)
    doc.value("xDhi", ref.dmax[0], 3)
    doc.value("yDlo", ref.dmin[1], 3)
    doc.value("yDhi", ref.dmax[1], 3)
    doc.value("axisDeg", ref.axis_deg, 3)
    doc.value("halfDeg", ref.half_deg, 3)
    half_rad = math.radians(ref.half_deg)
    doc.value("halfRad", half_rad, 8)
    doc.value("H", H, 3)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("maxCharges", lvl["limits"]["max_charges"])
    doc.value("nSample", 64)
    doc.value("nBisect", 80)

    # ---- Launch, energy, the first push ------------------------------------------------------
    r_up = math.hypot(ref.x0[0] - x_up, ref.x0[1] - y_up)
    r_at = math.hypot(ref.x0[0] - x_at, ref.x0[1] - y_at)
    doc.value("rUp", r_up, 8)
    doc.value("rAt", r_at, 8)
    V0 = ref.V(*ref.x0)
    W = ref.gamma0 * ref.m * ref.c ** 2 + V0
    Wbind = ref.T0 + V0
    doc.value("VZero", V0, 8)
    doc.value("W", W, 10)
    doc.value("Wbind", Wbind, 8)
    doc.check(Wbind > 0.0, "W - mc^2 is positive, so energy does not confine the particle")
    f0 = ref.rhs(ref.launch_state())
    doc.value("FxZero", f0[2], 6)
    doc.value("FyZero", f0[3], 6)

    # ---- The reference flight ----------------------------------------------------------------
    doc.value("steps", ref_f["steps"])
    doc.value("tE", ref_f["t"], 12)
    doc.value("xE", ref_f["x"][0], 12)
    doc.value("yE", ref_f["x"][1], 12)
    doc.value("pxE", ref_f["p"][0], 12)
    doc.value("pyE", ref_f["p"][1], 12)
    doc.text("faceE", ref_f["face"])
    doc.value("faceRes", ref_f["residual"], 2, sci=True)
    doc.value("uE", ref_f["u"], 6)
    doc.value("tBefore", ref_f["t_before"], 8)
    doc.value("xBefore", ref_f["before"][0], 8)
    doc.value("yBefore", ref_f["before"][1], 8)
    doc.value("xAfter", ref_f["after"][0], 8)
    doc.value("yAfter", ref_f["after"][1], 8)
    beta = ref_f["beta"]
    doc.value("betaE", beta, 8)
    ang = math.degrees(beta)
    doc.value("angE", ang, 6)
    doc.value("marginDeg", ref.half_deg - ang, 6)
    doc.value("errRef", step_gap(ref_f, ref_c), 2, sci=True)
    doc.value("wDev", ref_f["wdev"], 2, sci=True)
    Ve = ref.V(*ref_f["x"])
    gamma_e = (W - Ve) / (ref.m * ref.c ** 2)
    p_energy = ref.m * ref.c * math.sqrt(gamma_e * gamma_e - 1.0)
    p_int = math.hypot(*ref_f["p"])
    doc.value("gammaE", gamma_e, 8)
    doc.value("pEnergy", p_energy, 10)
    doc.value("pDiff", abs(p_energy - p_int), 2, sci=True)
    # Closest approach: the attractor is the last charge; the wall is the first n_wall.
    i_at = len(ref.charges) - 1
    doc.value("rClose", ref_f["closest"][i_at], 6)
    doc.value("tClose", ref_f["close_xy"][i_at][0], 6)
    doc.value("xClose", ref_f["close_xy"][i_at][1], 6)
    doc.value("yClose", ref_f["close_xy"][i_at][2], 6)
    doc.value("clearAt", ref_f["closest"][i_at] - ref.R, 4)
    doc.value("rOverR", ref_f["closest"][i_at] / ref.R, 4)
    f_close = ref.rhs((ref_f["close_xy"][i_at][1], ref_f["close_xy"][i_at][2], 0.0, 0.0))
    doc.value("FxClose", f_close[2], 5)
    doc.value("FyClose", f_close[3], 5)
    doc.check(f_close[3] > 0.0, "at the closest approach the force points upward")
    doc.check(ref_f["closest"][i_at] > ref.R, "the reference does not hit the attractor")
    i_wall = min(range(ref.n_wall), key=lambda i: ref_f["closest"][i])
    doc.value("rWall", ref_f["closest"][i_wall], 6)
    doc.value("yWallHit", ref.charges[i_wall][1], 3)
    doc.value("clearWall", ref_f["closest"][i_wall] - ref.R, 4)
    doc.value("tDip", ref_f["y_min"][0], 6)
    doc.value("xDip", ref_f["y_min"][1], 6)
    doc.value("yDip", ref_f["y_min"][2], 6)
    # Where the rising path crosses the two vertical lines used in the account.
    def crossing_x(run_unused, target):
        prev = None
        for x, y in ref_f["path"]:
            if prev is not None and prev[0] < target <= x:
                span = x - prev[0]
                frac = 0.0 if span == 0.0 else (target - prev[0]) / span
                return target, prev[1] + frac * (y - prev[1])
            prev = (x, y)
        return None

    x15, y15 = crossing_x(None, ref.charges[0][0])
    x20, y20 = crossing_x(None, 20.0)
    doc.value("yPass", y15, 6)
    doc.value("yRise", y20, 6)
    doc.value("xRise", x20, 3)
    doc.check(y15 < wall[0][1] - ref.R, "the reference passes under the wall")

    def part_force(kappa, cx, cy):
        dx, dy = ref.x0[0] - cx, ref.x0[1] - cy
        r2 = dx * dx + dy * dy
        inv = kappa / (r2 * math.sqrt(r2))
        return dx * inv, dy * inv

    _, fy_up = part_force(ref.q * q_up, x_up, y_up)
    _, fy_at = part_force(ref.q * q_at, x_at, y_at)
    doc.value("FyUp", fy_up, 5)
    doc.value("FyAtLaunch", fy_at, 5)
    doc.check(fy_up < 0.0 and fy_at < 0.0, "both player charges push or pull downward at launch")

    def nearest(target):
        return min(ref_f["rows"], key=lambda r: abs(r[0] - target))

    for label, target in (("Five", 5.0), ("Fifteen", 15.0), ("TwentyFive", 25.0), ("Thirty", 30.0)):
        tt, xx, yy, _px, py, _dw = nearest(target)
        doc.value("t" + label, tt, 5)
        doc.value("x" + label, xx, 6)
        doc.value("y" + label, yy, 6)
        doc.value("py" + label, py, 6)
    doc.check(nearest(15.0)[4] < 0.0, "the particle is still descending at 15 time units")
    doc.check(nearest(25.0)[4] > 0.0, "the particle is climbing at 25 time units")

    e_fine = radiated(ref, H / 2.0, ref_f["t"])
    e_coarse = radiated(ref, H, ref_f["t"])
    rad_err = abs(e_fine - e_coarse) / abs(e_fine) / 15.0
    doc.value("Erad", e_fine, 4, sci=True)
    doc.value("radErr", rad_err, 2, sci=True)
    doc.value("radFrac", e_fine / ref.T0, 2, sci=True)
    doc.check(rad_err < TOL_RAD, "the radiation quadrature is inside its tolerance")

    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$\Delta W$"],
              [(t, x, y, px, py, dw) for t, x, y, px, py, dw in ref_f["rows"]],
              digits=[5, 5, 5, 5, 5, 2],
              sci=[False, False, False, False, False, True])
    doc.data("path", ref_f["path"])
    doc.data("wall", [(cx, cy) for cx, cy, _ in wall])

    # Coarser steps, for the convergence table only. Not used to set the tolerance.
    conv_rows = []
    for h in (0.008, 0.004, 0.002, H):
        run = finish(ref, fly(ref, h))
        ang_h = "" if run["beta"] is None else run["beta"] * 180.0 / math.pi
        conv_rows.append((h, run["outcome"], run["t"], run["x"][1], ang_h if ang_h != "" else "---"))
    doc.table("convergence",
              [r"$h$", "outcome", r"$t_\mathrm{e}$", r"$y_\mathrm{e}$", r"angle ($^\circ$)"],
              conv_rows, digits=[4, 5, 8, 8, 5])

    # ---- The two charges separately ----------------------------------------------------------
    def one(name, st, fine, coarse):
        doc.text(name + "Outcome", fine["outcome"])
        doc.text(name + "Face", fine["face"])
        doc.value(name + "T", fine["t"], 8)
        doc.value(name + "X", fine["x"][0], 8)
        doc.value(name + "Y", fine["x"][1], 8)
        doc.value(name + "Px", fine["p"][0], 8)
        doc.value(name + "Py", fine["p"][1], 8)
        doc.value(name + "Err", step_gap(fine, coarse), 2, sci=True)
        doc.value(name + "Res", fine["residual"], 2, sci=True)
        e1 = radiated(st, H / 2.0, fine["t"])
        e2 = radiated(st, H, fine["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 != 0.0 else 0.0
        doc.value(name + "Rad", e1, 4, sci=True)
        doc.value(name + "RadErr", rel, 2, sci=True)
        doc.check(rel < TOL_RAD, name + " radiation quadrature is inside its tolerance")
        doc.data("path" + name.lower(), fine["path"])
        return e1

    up_rad = one("Up", up, up_f, up_c)
    at_rad = one("At", at, at_f, at_c)

    # The cone drawn from the middle of the entrance, length 2.5, for the figure only.
    mouth_y = 0.5 * (ref.dmin[1] + ref.dmax[1])
    length = 2.5
    doc.value("coneOx", ref.dmin[0], 3)
    doc.value("coneOy", mouth_y, 3)
    doc.value("coneXp", ref.dmin[0] + length * math.cos(half_rad), 6)
    doc.value("coneYp", mouth_y + length * math.sin(half_rad), 6)
    doc.value("coneXm", ref.dmin[0] + length * math.cos(-half_rad), 6)
    doc.value("coneYm", mouth_y + length * math.sin(-half_rad), 6)
    doc.value("arrowX", ref.x0[0] + 1.6, 3)
    doc.value("arrowY", ref.x0[1], 3)

    # ---- Comparison, after the tolerance is fixed --------------------------------------------
    def put(i, tag, title, fine, rad):
        g = doc.flight(i)
        ge = g["end"]
        doc.compare(tag + "Outcome", title + " outcome", fine["outcome"], g["outcome"])
        doc.compare(tag + "T", title + " $t$", fine["t"], ge["t"], tol)
        doc.compare(tag + "X", title + " $x$", fine["x"][0], ge["x"][0], tol)
        doc.compare(tag + "Y", title + " $y$", fine["x"][1], ge["x"][1], tol)
        doc.compare(tag + "Px", title + " $p_x$", fine["p"][0], ge["p"][0], tol)
        doc.compare(tag + "Py", title + " $p_y$", fine["p"][1], ge["p"][1], tol)
        if fine["beta"] is not None:
            gb = beta_of(ge["p"], ref.axis_deg)
            doc.compare(tag + "Beta", title + r" $\beta$ (rad)", fine["beta"], gb, tol)
        doc.compare(tag + "Rad", title + " radiated energy (rel.)", rad, g["radiated_energy"],
                    TOL_RAD, rel=True)

    put(0, "Ref", "reference", ref_f, e_fine)
    put(1, "Up", "upper charge", up_f, up_rad)
    put(2, "At", "attractor", at_f, at_rad)
