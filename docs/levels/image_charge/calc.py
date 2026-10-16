"""Companion of docs/levels/image_charge/doc.tex (Level 17, Image charge).

One grounded sphere. The image is the one derived for the dome
(high_voltage_dome eq:kelvin). The potential is held at zero, so there is no
centre charge (high_voltage_dome eq:centre with V = 0), and the ion's own image
is that grounded image alone (high_voltage_dome eq:self). The flight is the
classical Runge-Kutta method at a fixed step (Foundations eq:rk4). A face or a
contact is the cubic Hermite interpolant of that step. The truncation error is
estimated by halving the step, and the tolerance is set from that estimate
before any comparison with the game.
"""
import math


# Reported step. Halving is the coarser run: the error of the reported step is
# |y_h - y_{2h}|/15 (Foundations, the numerics section). The tolerance is ten
# times the largest such estimate and at least 1e-9. Computed in compute()
# before doc.compare, not from the game's numbers.
H = 0.005
# Simpson on the radiated power. 1e-4 separates a wrong formula from the game's
# quadrature of the same power. Set beforehand.
TOL_RAD = 1e-4
# The obstacle is the sphere plus this contact shell (PHYSICS.md, conductors).
CONTACT = 0.02


class Setup:
    """One shot, one sphere, and the charges of one placement."""

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
        nrm = math.hypot(d[0], d[1])
        self.d = [d[0] / nrm, d[1] / nrm]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        cond = lvl["conductors"][0]
        bias = cond["bias"]
        self.cx = cond["center"][0] / s
        self.cy = cond["center"][1] / s
        self.radius = cond["radius"] / s
        self.kind = bias["kind"]
        self.bias_value = bias["value"] if "value" in bias else 0.0
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)
        self.fixed, self.alpha, self.induced_net = self._fixed_induced()
        self.hit_r = self.radius + CONTACT + self.a
        self.charge_r = self.R + self.a

    def kelvin(self, x, y, q):
        """Image of a charge q at (x, y) in the grounded sphere (eq:kelvin)."""
        dx = x - self.cx
        dy = y - self.cy
        r2 = dx * dx + dy * dy
        r = math.sqrt(r2)
        qi = -q * self.radius / r
        scale = self.radius * self.radius / r2
        return (self.cx + dx * scale, self.cy + dy * scale, qi)

    def _fixed_induced(self):
        images = []
        induced = 0.0
        for x, y, q in self.charges:
            ix, iy, qi = self.kelvin(x, y, q)
            images.append((ix, iy, qi))
            induced += qi
        # Capacitance of one sphere is its radius. Fixed potential: alpha = V.
        # Isolated net charge Q: alpha = (Q - induced) / a. Grounded: alpha = 0.
        if self.kind == "potential":
            alpha = self.bias_value
        elif self.kind == "charge":
            alpha = (self.bias_value - induced) / self.radius
        else:
            alpha = 0.0
        if alpha != 0.0:
            images.append((self.cx, self.cy, alpha * self.radius))
        return images, alpha, induced

    def self_charges(self, x, y):
        """Charges the particle induces (eq:self). A held potential does not add
        a centre charge that follows the particle; an isolated sphere does."""
        ix, iy, qi = self.kelvin(x, y, self.q)
        out = [(ix, iy, qi)]
        if self.kind == "charge":
            r = math.hypot(x - self.cx, y - self.cy)
            out.append((self.cx, self.cy, self.q * self.radius / r))
        return out

    def _sum(self, x, y, charges):
        phi = ex = ey = 0.0
        for cx, cy, q in charges:
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            r = math.sqrt(r2)
            phi += q / r
            inv_r3 = q / (r2 * r)
            ex += dx * inv_r3
            ey += dy * inv_r3
        return phi, ex, ey

    def self_force_mag(self, x, y):
        _phi, ex, ey = self._sum(x, y, self.self_charges(x, y))
        return abs(self.q) * math.hypot(ex, ey)

    def fields(self, x, y):
        phi_c, ex_c, ey_c = self._sum(x, y, self.charges)
        phi_i, ex_i, ey_i = self._sum(x, y, self.fixed)
        phi_s, ex_s, ey_s = self._sum(x, y, self.self_charges(x, y))
        return (phi_c + phi_i, phi_s, ex_c + ex_i + ex_s, ey_c + ey_i + ey_s)

    def surface_residual(self):
        worst = 0.0
        for k in range(16):
            ang = 2.0 * math.pi * k / 16.0
            x = self.cx + self.radius * math.cos(ang)
            y = self.cy + self.radius * math.sin(ang)
            phi, _, _ = self._sum(x, y, self.charges)
            phi_i, _, _ = self._sum(x, y, self.fixed)
            worst = max(worst, abs(phi + phi_i - self.alpha))
        return worst

    def rhs(self, state):
        x, y, px, py = state
        m = self.m
        c = self.c
        p2 = px * px + py * py
        g = math.sqrt(1.0 + p2 / (m * m * c * c))
        inv = 1.0 / (g * m)
        _phi, _phis, ex, ey = self.fields(x, y)
        return (px * inv, py * inv, self.q * ex, self.q * ey)

    def energy(self, state):
        x, y, px, py = state
        kinetic = (math.sqrt(1.0 + (px * px + py * py) / (self.m * self.m * self.c * self.c)) - 1.0)
        kinetic *= self.m * self.c * self.c
        phi, phis, _ex, _ey = self.fields(x, y)
        return kinetic + self.q * phi + 0.5 * self.q * phis

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
    """Integrate until the first ending event. Obstacles, then the world, then
    the detector (the engine's order)."""
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
    sphere_close = math.hypot(s[0] - st.cx, s[1] - st.cy)
    y_max = s[1]
    y_min = s[1]
    self_peak = st.self_force_mag(s[0], s[1])

    def note(x, y):
        nonlocal sphere_close, y_max, y_min, self_peak
        if y > y_max:
            y_max = y
        if y < y_min:
            y_min = y
        rr = math.hypot(x - st.cx, y - st.cy)
        if rr < sphere_close:
            sphere_close = rr
        ff = st.self_force_mag(x, y)
        if ff > self_peak:
            self_peak = ff

    def absorb(u, s0, f0, s1, f1):
        p = hermite(s0, f0, s1, f1, h, u)
        note(p[0], p[1])

    def pack(outcome, te, pe, se, n_steps, face, residual):
        return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n_steps,
                "face": face, "residual": residual, "rows": rows, "path": path,
                "wdev": wdev, "w0": w0, "sphere_close": sphere_close,
                "y_max": y_max, "y_min": y_min, "self_peak": self_peak}

    n_charges = len(st.charges)
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

        for idx, (cx, cy, _q) in enumerate(st.charges):
            dx = 0.0 if box_lo[0] <= cx <= box_hi[0] else (box_lo[0] - cx if cx < box_lo[0] else cx - box_hi[0])
            dy = 0.0 if box_lo[1] <= cy <= box_hi[1] else (box_lo[1] - cy if cy < box_lo[1] else cy - box_hi[1])
            if dx * dx + dy * dy > st.charge_r * st.charge_r:
                continue
            limit = st.charge_r
            u = _first(lambda p, cx=cx, cy=cy, limit=limit:
                       (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit, samples(), us)
            if u is not None:
                hits.append((u, 0, idx, "Collided(%d)" % idx))
        dx = 0.0 if box_lo[0] <= st.cx <= box_hi[0] else (box_lo[0] - st.cx if st.cx < box_lo[0] else st.cx - box_hi[0])
        dy = 0.0 if box_lo[1] <= st.cy <= box_hi[1] else (box_lo[1] - st.cy if st.cy < box_lo[1] else st.cy - box_hi[1])
        if dx * dx + dy * dy <= st.hit_r * st.hit_r:
            u = _first(lambda p: (p[0] - st.cx) ** 2 + (p[1] - st.cy) ** 2 <= st.hit_r * st.hit_r,
                       samples(), us)
            if u is not None:
                hits.append((u, 0, n_charges, "Collided(%d)" % n_charges))
        overlaps_world = (box_lo[0] < st.wmin[0] or box_hi[0] > st.wmax[0]
                          or box_lo[1] < st.wmin[1] or box_hi[1] > st.wmax[1])
        if overlaps_world:
            u = _first(lambda p: not _in_box(p, st.wmin, st.wmax), samples(), us)
            if u is not None:
                hits.append((u, 1, -1, "LeftBounds"))
        overlaps_det = not (box_hi[0] < st.dmin[0] or box_lo[0] > st.dmax[0]
                            or box_hi[1] < st.dmin[1] or box_lo[1] > st.dmax[1])
        if overlaps_det:
            u = _first(lambda p: _in_box(p, st.dmin, st.dmax), samples(), us)
            if u is not None:
                hits.append((u, 2, -1, "Arrived"))
        if hits:
            u_hit, _rank, idx, outcome = min(hits)
            if outcome == "Arrived":
                pred = lambda p: _in_box(p, st.dmin, st.dmax)
                box = (st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                pred = lambda p: not _in_box(p, st.wmin, st.wmax)
                box = (st.wmin, st.wmax)
            elif idx < n_charges:
                cx, cy, _q = st.charges[idx]
                limit = st.charge_r
                pred = lambda p, cx=cx, cy=cy, limit=limit: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit
                box = None
            else:
                pred = lambda p: (p[0] - st.cx) ** 2 + (p[1] - st.cy) ** 2 <= st.hit_r * st.hit_r
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
            else:
                if idx < n_charges:
                    cx, cy, _q = st.charges[idx]
                    residual = abs(math.hypot(pe[0] - cx, pe[1] - cy) - st.charge_r)
                    face = "charge"
                else:
                    residual = abs(math.hypot(pe[0] - st.cx, pe[1] - st.cy) - st.hit_r)
                    face = "sphere"
            if pe[1] > y_max:
                y_max = pe[1]
            if pe[1] < y_min:
                y_min = pe[1]
            path.append(pe)
            return pack(outcome, te, pe, se, n + 1, face, residual)
        for k in range(0, 9):
            absorb(k / 8.0, s, f, s2, f2)
        note(s2[0], s2[1])
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
    return pack("Timeout", t, (s[0], s[1]), s, n, "", 0.0)


def radiated(st, h, t_end):
    m = st.m
    c = st.c
    q = st.q

    def power(state):
        _x, _y, px, py = state
        g = math.sqrt(1.0 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        _vx, _vy, fx, fy = st.rhs(state)
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
    fine = fly(st, h, path_every, table_every)
    coarse = fly(st, 2.0 * h)
    return fine, coarse


def compute(doc):
    lvl = doc.level
    ref_place = lvl["reference_solution"]
    bare_place = doc.result(1)["placement"]
    ref = Setup(lvl, ref_place, 0)
    bare = Setup(lvl, bare_place, 0)
    runs = {
        0: (ref, *pair(ref, H, path_every=25, table_every=5.0)),
        1: (bare, *pair(bare, H, path_every=25)),
    }
    # One sphere: the image construction meets the boundary condition exactly.
    # The tolerance below is the flight's truncation, not a model residual.
    surf = ref.surface_residual()
    doc.check(ref.kind == "grounded" and ref.alpha == 0.0, "the sphere is not grounded")
    doc.check(len(ref.fixed) == len(ref.charges), "a grounded sphere has a centre charge")
    doc.check(surf < 1e-9, "the image potential on the sphere is not zero")
    doc.check(bare.surface_residual() < 1e-9, "the bare sphere is not at zero potential")
    doc.check(runs[0][1]["outcome"] == runs[0][2]["outcome"], "the reference outcome changes when the step is halved")
    doc.check(runs[1][1]["outcome"] == runs[1][2]["outcome"], "the bare outcome changes when the step is halved")
    err_max = max(step_gap(runs[i][1], runs[i][2]) for i in (0, 1))
    # Ten times the step-halving estimate, and at least the verified-flight floor.
    tol = max(1e-9, 10.0 * err_max)

    doc.value("H", H, 4)
    doc.value("contact", CONTACT, 2)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("errMax", err_max, 2, sci=True)
    doc.value("gamma0", ref.gamma0, 6)
    doc.value("p0", ref.p0, 8)
    doc.value("v0", ref.p0 / (ref.gamma0 * ref.m), 6)
    doc.value("vOverC", ref.p0 / (ref.gamma0 * ref.m * ref.c), 4)
    doc.value("sphereX", ref.cx, 4)
    doc.value("sphereY", ref.cy, 4)
    doc.value("sphereA", ref.radius, 4)
    doc.value("alpha", ref.alpha, 2)
    doc.value("inducedG", ref.induced_net, 6, sci=True)
    doc.value("qPart", ref.q, 3, sci=True)
    doc.value("mPart", ref.m, 2)
    doc.value("T0", ref.T0, 2)
    doc.value("cLight", ref.c, 2)
    doc.value("Rcharge", ref.R, 2)
    doc.value("xLaunch", ref.x0[0], 3)
    doc.value("yLaunch", ref.x0[1], 3)
    doc.value("xDetLo", ref.dmin[0], 3)
    doc.value("yDetLo", ref.dmin[1], 3)
    doc.value("xDetHi", ref.dmax[0], 3)
    doc.value("yDetHi", ref.dmax[1], 3)
    doc.value("xWorld", ref.wmax[0], 3)
    doc.value("yWorld", ref.wmax[1], 3)
    doc.value("sampleArc", ref.c * H / 8.0, 4)
    doc.value("arrowX", ref.x0[0] + 1.6 * ref.d[0], 4)
    doc.value("arrowY", ref.x0[1] + 1.6 * ref.d[1], 4)
    doc.value("surfRes", surf, 2, sci=True)
    doc.check(len(ref.charges) == 1, "the reference does not have one charge")
    cx, cy, q = ref.charges[0]
    ix, iy, qi = ref.kelvin(cx, cy, q)
    doc.value("xPlay", cx, 3)
    doc.value("yPlay", cy, 3)
    doc.value("QPlay", q, 3, sci=True)
    doc.value("imageX", ix, 6)
    doc.value("imageY", iy, 6)
    doc.value("imageQ", qi, 6, sci=True)
    doc.value("sourceR", math.hypot(cx - ref.cx, cy - ref.cy), 6)
    own = ref.self_charges(ref.x0[0], ref.x0[1])
    doc.check(len(own) == 1, "the grounded image has a centre charge")
    doc.value("selfX", own[0][0], 6)
    doc.value("selfY", own[0][1], 6)
    doc.value("selfQ", own[0][2], 4, sci=True)
    phi0, phis0, _ex, _ey = ref.fields(ref.x0[0], ref.x0[1])
    doc.value("phi0", phi0, 6, sci=True)
    doc.value("phiSelf0", phis0, 4, sci=True)
    doc.value("qPhi0", ref.q * phi0, 6)
    doc.value("halfSelf0", 0.5 * ref.q * phis0, 4, sci=True)

    def publish(i, prefix):
        st, fine, coarse = runs[i]
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
        doc.value(prefix + "Rclose", fine["sphere_close"], 6)
        doc.value(prefix + "Clear", fine["sphere_close"] - st.radius, 6)
        e1 = radiated(st, H / 2.0, fine["t"])
        e2 = radiated(st, H, fine["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
        doc.check(rel < TOL_RAD, prefix + " radiation quadrature is inside its tolerance")
        frac = e1 / st.T0
        doc.value(prefix + "Rad", e1, 4, sci=True)
        doc.value(prefix + "RadErr", rel, 2, sci=True)
        doc.value(prefix + "RadFrac", frac, 2, sci=True)
        doc.data("path" + prefix.lower(), fine["path"])
        return e1, frac

    rad_ref, frac_ref = publish(0, "Ref")
    _rad_bare, frac_bare = publish(1, "Bare")
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$\Delta W$"],
              list(runs[0][1]["rows"]),
              digits=[5, 5, 5, 5, 5, 2],
              sci=[False, False, False, False, False, True])
    doc.value("clearMin", min(runs[i][1]["sphere_close"] - runs[i][0].radius for i in (0, 1)), 6)
    doc.value("radFracMax", max(frac_ref, frac_bare), 2, sci=True)
    # Peak |q E_self| on the sampled flight, not the grounded image-only formula.
    # The samples are the step ends and nine Hermite stations on each step.
    self_force = max(runs[i][1]["self_peak"] for i in (0, 1))
    r_near = min(runs[i][1]["sphere_close"] for i in (0, 1))
    doc.value("rNear", r_near, 6)
    doc.value("selfForce", self_force, 2, sci=True)
    doc.value("selfBound", self_force / ref.T0, 2, sci=True)
    doc.check(runs[1][0].charges == [], "the bare placement adds a charge")
    doc.check(runs[0][1]["sphere_close"] > runs[0][0].hit_r, "the reference meets the sphere")

    def put(i, prefix, title, with_rad, ours_rad):
        _st, fine, _coarse = runs[i]
        g = doc.flight(i)
        ge = g["end"]
        doc.compare(prefix + "Outcome", title + "outcome", fine["outcome"], g["outcome"])
        doc.compare(prefix + "T", title + "$t$", fine["t"], ge["t"], tol)
        doc.compare(prefix + "X", title + "$x$", fine["x"][0], ge["x"][0], tol)
        doc.compare(prefix + "Y", title + "$y$", fine["x"][1], ge["x"][1], tol)
        if with_rad and g["radiated_energy"] != 0.0:
            doc.compare(prefix + "Rad", title + "$E$ (rel.)", ours_rad, g["radiated_energy"],
                        TOL_RAD, rel=True)

    put(0, "Ref", "", True, rad_ref)
    put(1, "Bare", "bare ", False, 0.0)
