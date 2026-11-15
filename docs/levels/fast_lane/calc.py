"""Companion of docs/levels/fast_lane/doc.tex (Level 27, Fast lane).

Every number of the document, from the level file, by the document's method: the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4) on the relativistic equations of
motion in the fixed charges' Coulomb field (eq:eom-plane), the entry into the detector by
cubic Hermite interpolation of the last step and bisection, the error estimated by halving
the step. Beside it, the same launch at the same speed in Newton's mechanics (eq:newton),
and the impulse estimate of a passing charge's deflection (eq:impulse), which show why the
fast particle is hard to bend. Compared with the game: the reference flight and the empty
board.
"""
import math


# The method's own error, estimated by doubling the step, is about 2e-11 at the step used
# (section 4.3); the game's verified flights are accurate to about 1e-10 (Foundations,
# section 6). Ten times the larger: 1e-9, in cells and time units. The reference flight
# grazes a charge (0.05 outside its radius), which magnifies every integration error: there
# the game's own estimate is larger, and a flight's tolerance is at least the game's gap
# between its preview and its verification (tolerances 1e-10 and 1e-12: the verified run's
# error is at most a tenth of that gap, so the gap is ten times it).
TOL = 1e-9
# The radiated energy: Simpson's rule on the RK4 path against the game's own quadrature of
# the same power (3e-9 relative on an exact test, PHYSICS.md 8); 1e-4 still tells a wrong
# formula apart.
TOL_RAD = 1e-4
H = 0.00125        # the step of the reference computation
H_TABLE = 0.01     # the step of the table of states (every 5 time units)


class Setup:
    """The level's numbers (Foundations 5) and the charges of a placement. With
    `newton=True` the momentum is m v (c taken to infinity in the mechanics only)."""

    def __init__(self, lvl, placement, newton=False):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.newton = newton
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.a = shot["particle"]["radius"]
        self.c = lvl["physics"]["c"]
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
        # The fixed charges of the level, then the placement's (the game's order).
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        self.kappas = [(self.q * Q, cx, cy) for cx, cy, Q in self.charges]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)                  # eq:launch (F)
        self.p_rel = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)
        self.v0 = self.p_rel / (self.gamma0 * self.m)
        # The Newtonian twin starts at the same speed v0 (its momentum m v0).
        self.p0 = self.m * self.v0 if newton else self.p_rel

    def V(self, x, y):
        """The potential energy q*phi (Foundations eq:coulomb, superposed)."""
        return sum(k / math.hypot(x - cx, y - cy) for k, cx, cy in self.kappas)

    def rhs(self, s):
        """The right-hand side of eq:eom-plane (or eq:newton)."""
        x, y, px, py = s
        g = 1.0 if self.newton else math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        fx = fy = 0.0
        for k, cx, cy in self.kappas:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += k * dx / r3
            fy += k * dy / r3
        return (px / (g * self.m), py / (g * self.m), fx, fy)

    def energy(self, s):
        """W = gamma m c^2 + V (eq:energy-plane), or p^2/2m + V for the twin."""
        x, y, px, py = s
        if self.newton:
            return (px * px + py * py) / (2 * self.m) + self.V(x, y)
        return math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2) + self.V(x, y)


def rk4_step(st, s, f, h):
    """One step of the classical Runge-Kutta method (Foundations eq:rk4)."""
    k1 = f
    k2 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)))
    k3 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)))
    k4 = st.rhs(tuple(a + h * b for a, b in zip(s, k3)))
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    """The cubic Hermite interpolant of the position on a step, at the fraction u."""
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def box_distance(p, lo, hi):
    """Signed distance to a box, measured per axis (negative inside): its sign is all that
    the bisection needs."""
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def outside_world(st, p):
    """Positive outside the world: the bisection's predicate for LeftBounds."""
    return max(st.wmin[0] - p[0], p[0] - st.wmax[0], st.wmin[1] - p[1], p[1] - st.wmax[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


def fly(st, h, record_every=None):
    """Integrates from the launch until the first event (Foundations 5): the detector, a
    charge, or the world's edge, located on the step's interpolant. Returns a dict with the
    outcome, the state at the event, the steps, and (optionally) a table of states."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    rows, path, wdev = [], [(s[0], s[1])], 0.0
    rmin, rarg = min((math.hypot(s[0] - cx, s[1] - cy), i) for i, (_, cx, cy) in enumerate(st.kappas))
    while True:
        s2 = rk4_step(st, s, f, h)
        f2 = st.rhs(s2)
        us = [k / 64 for k in range(1, 65)]
        pts = [hermite(s, f, s2, f2, h, u) for u in us]
        events = []
        for k, p in enumerate(pts):
            if box_distance(p, st.dmin, st.dmax) <= 0:
                events.append((us[k], 2, "Arrived", None))
                break
        for idx, (_, cx, cy) in enumerate(st.kappas):
            for k, p in enumerate(pts):
                if math.hypot(p[0] - cx, p[1] - cy) <= st.R + st.a:
                    events.append((us[k], 0, f"Collided({idx})", idx))
                    break
        for k, p in enumerate(pts):
            if outside_world(st, p) > 0:
                events.append((us[k], 1, "LeftBounds", None))
                break
        for p in pts:
            rmin, rarg = min((rmin, rarg), min((math.hypot(p[0] - cx, p[1] - cy), i)
                                               for i, (_, cx, cy) in enumerate(st.kappas)))
        if events:
            u_hit, _rank, outcome, idx = min(events)
            if outcome == "Arrived":
                pred = lambda p: box_distance(p, st.dmin, st.dmax) <= 0
            elif outcome == "LeftBounds":
                pred = lambda p: outside_world(st, p) > 0
            else:
                cx, cy = st.kappas[idx][1], st.kappas[idx][2]
                pred = lambda p: math.hypot(p[0] - cx, p[1] - cy) <= st.R + st.a
            lo, hi = u_hit - 1 / 64, u_hit
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if pred(hermite(s, f, s2, f2, h, mid)):
                    hi = mid
                else:
                    lo = mid
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            if outcome == "Arrived":
                face = face_of(pe, st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                face = face_of(pe, st.wmin, st.wmax)
            else:
                face = ""
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "steps": n + 1,
                    "face": face, "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                    "rmin": rmin, "rarg": rarg}
        s, f, t, n = s2, f2, t + h, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if record_every and n % record_every == 0:
            rows.append((t, s, st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if t > st.tmax:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "rows": rows,
                    "path": path, "steps": n, "wdev": wdev, "w0": w0, "face": "", "rmin": rmin,
                    "rarg": rarg}


def radiated(st, h, t_end):
    """The radiated energy along the RK4 path (Foundations eq:larmor) with the acceleration of
    Level 1 eq:acceleration, by Simpson's rule on the steps (Foundations eq:simpson)."""
    def power(s):
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (st.m ** 2 * st.c ** 2))
        v = (px / (g * st.m), py / (g * st.m))
        _, _, fx, fy = st.rhs(s)
        vF = v[0] * fx + v[1] * fy
        ax = (fx - v[0] * vF / st.c ** 2) / (g * st.m)
        ay = (fy - v[1] * vF / st.c ** 2) / (g * st.m)
        vxa = v[0] * ay - v[1] * ax
        return 2 / 3 * st.q ** 2 / st.c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / st.c ** 2)

    n = int(t_end / h)
    n -= n % 2
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    vals = [power(s)]
    for _ in range(n):
        s = rk4_step(st, s, st.rhs(s), h)
        vals.append(power(s))
    simpson = h / 3 * (vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2]))
    rest = t_end - n * h
    s2 = rk4_step(st, s, st.rhs(s), rest)
    return simpson + rest * 0.5 * (vals[-1] + power(s2))


def err_of(a, b):
    """The step-halving estimate of an end state (Foundations 6): |a - b| / 15."""
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1])) / 15


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    st = Setup(lvl, ref)
    doc.value("q", st.q, 3, sci=True)
    doc.value("m", st.m, 3)
    doc.value("c", st.c, 3)
    doc.value("Tzero", st.T0, 3)
    doc.value("R", st.R, 3)
    doc.value("xZero", st.x0[0], 3)
    doc.value("yZero", st.x0[1], 3)
    doc.value("dirX", st.d[0], 10)
    doc.value("dirY", st.d[1], 10)
    doc.value("launchAngle", math.degrees(math.atan2(st.d[1], st.d[0])), 6)
    doc.value("xDetLo", st.dmin[0], 3)
    doc.value("yDetLo", st.dmin[1], 3)
    doc.value("xDetHi", st.dmax[0], 3)
    doc.value("yDetHi", st.dmax[1], 3)
    doc.value("xWorld", st.wmax[0], 3)
    doc.value("yWorld", st.wmax[1], 3)
    names = ["One", "Two", "Three", "Four"]
    nfixed = len([e for e in lvl["elements"] if e["kind"] == "charge"])
    for k, (cx, cy, Q) in enumerate(st.charges):
        doc.value(f"Q{names[k]}", Q, 3, sci=True)
        doc.value(f"x{names[k]}", cx, 3)
        doc.value(f"y{names[k]}", cy, 3)
        doc.value(f"k{names[k]}", st.q * Q, 4)
    doc.check(nfixed == 3 and len(ref) == 1, "three fixed charges and one placed")

    # ---- The launch (Foundations eq:launch) -------------------------------------------------
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p_rel, 12)
    doc.value("vZero", st.v0, 10)
    doc.value("betaZero", st.v0 / st.c, 6)
    doc.value("pNewton", st.m * st.v0, 10)
    doc.value("TNewton", 0.5 * st.m * st.v0 ** 2, 10)
    V0 = st.V(*st.x0)
    W = st.gamma0 * st.m * st.c ** 2 + V0
    doc.value("VZero", V0, 10)
    doc.value("W", W, 12)

    # ---- The impulse estimate (eq:impulse): each fixed charge seen from the launch line -----
    # The launch line x0 + s d; its impact parameter with each charge, and the deflection
    # angle 2|kappa|/(b p v) relativistically and 2|kappa|/(b m v^2) for Newton's twin.
    rows = []
    for k, (cx, cy, Q) in enumerate(st.charges[:nfixed]):
        rx, ry = cx - st.x0[0], cy - st.x0[1]
        b = abs(rx * st.d[1] - ry * st.d[0])
        kap = st.q * Q
        th_rel = 2 * abs(kap) / (b * st.p_rel * st.v0)
        th_new = 2 * abs(kap) / (b * st.m * st.v0 ** 2)
        side = "left" if rx * st.d[1] - ry * st.d[0] < 0 else "right"
        sense = "towards" if kap < 0 else "away from"
        rows.append([f"{k + 1}", f"({cx:g}, {cy:g})", Q, b, math.degrees(th_rel), math.degrees(th_new),
                     f"{sense}, {side}"])
        if k == 0:
            doc.value("bOne", b, 4)
            doc.value("thRelOne", math.degrees(th_rel), 3)
            doc.value("thNewOne", math.degrees(th_new), 3)
    doc.table("impulse", ["", "charge at", r"$Q$", r"$b$", r"$\theta$ (deg)", r"$\theta_\mathrm{N}$ (deg)",
                          "bends it"],
              rows, align="llrrrrl", digits=[None, None, 3, 4, 3, 3, None],
              sci=[False, False, True, False, False, False, False])
    doc.value("gammaRatio", st.p_rel / (st.m * st.v0), 6)
    # eq:pv: the product p v that sets the deflection, against 2 T0 (Newton at the same T0).
    pv = st.p_rel * st.v0
    pv_formula = st.T0 * (st.T0 + 2 * st.m * st.c ** 2) / (st.T0 + st.m * st.c ** 2)
    doc.check(abs(pv - pv_formula) < 1e-12 * pv, "p v from the launch agrees with eq:pv")
    doc.value("pvZero", pv, 6)
    doc.value("twoT", 2 * st.T0, 3)
    doc.value("pvRatio", 2 * st.T0 / pv, 4)

    # ---- The reference flight by RK4 (eq:rk4) ------------------------------------------------
    run = fly(st, H)
    half = fly(st, 2 * H)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    doc.value("H", H, 4)
    doc.value("steps", run["steps"])
    te, xe, pe = run["t"], run["x"], run["p"]
    doc.value("tE", te, 12)
    doc.value("xE", xe[0], 12)
    doc.value("yE", xe[1], 12)
    doc.text("faceE", run["face"])
    doc.value("pxE", pe[0], 10)
    doc.value("pyE", pe[1], 10)
    doc.value("angE", math.degrees(math.atan2(pe[1], pe[0])), 5)
    err = err_of(run, half)
    doc.value("err", err, 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    doc.value("rMin", run["rmin"], 4)
    doc.value("rMinClear", run["rmin"] - st.R - st.a, 2)
    doc.value("rMinWhich", run["rarg"] + 1)
    # |p| at the entry from the energy alone (eq:energy-plane): an independent check.
    ge = (W - st.V(*xe)) / (st.m * st.c ** 2)
    pe_energy = st.m * st.c * math.sqrt(ge * ge - 1)
    doc.value("gammaE", ge, 10)
    doc.value("pEnergy", pe_energy, 10)
    doc.value("pDiff", abs(math.hypot(*pe) - pe_energy), 2, sci=True)
    doc.check(abs(math.hypot(*pe) - pe_energy) < 1e-9, "|p| at the entry agrees with the energy")
    doc.check(err < 1e-9, "the step-halving estimate is below 1e-9")

    tab = fly(st, H_TABLE, record_every=int(round(5 / H_TABLE)))
    doc.value("HTable", H_TABLE, 3)
    trows = [[f"{t:.0f}", s[0], s[1], s[2], s[3], dw] for t, s, dw in tab["rows"]]
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$W(t) - W$"], trows,
              align="rrrrrr", digits=[None, 8, 8, 8, 8, 2], sci=[False] * 5 + [True])

    # ---- The radiated energy -----------------------------------------------------------------
    E_rad = radiated(st, H, te)
    E_rad2 = radiated(st, H / 2, te)
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / st.T0, 3, sci=True)
    doc.value("EradErr", abs(E_rad - E_rad2) / abs(E_rad2) / 15, 2, sci=True)

    # ---- The empty board and the Newtonian twin (not compared with the game) ------------------
    bare = Setup(lvl, [])
    rb, rb2 = fly(bare, H), fly(bare, 2 * H)
    doc.text("bareOut", rb["outcome"])
    doc.text("bareFace", rb["face"])
    doc.value("bareT", rb["t"], 10)
    doc.value("bareX", rb["x"][0], 10)
    doc.value("bareY", rb["x"][1], 10)
    doc.value("bareErr", err_of(rb, rb2), 2, sci=True)
    twin = fly(Setup(lvl, ref, newton=True), H)
    doc.text("twinOut", twin["outcome"])
    doc.value("twinT", twin["t"], 6)
    doc.value("twinX", twin["x"][0], 6)
    doc.value("twinY", twin["x"][1], 6)
    doc.text("twinFace", twin["face"])
    doc.value("twinWdev", twin["wdev"], 2, sci=True)
    twin_bare = fly(Setup(lvl, [], newton=True), H)
    doc.text("twinBareOut", twin_bare["outcome"])
    doc.value("twinBareX", twin_bare["x"][0], 6)
    doc.value("twinBareY", twin_bare["x"][1], 6)
    doc.data("path", run["path"])
    doc.data("pathbare", rb["path"])
    doc.data("pathtwin", twin["path"])

    # ---- The comparison with the game ----------------------------------------------------------
    def game_gap(i):
        a, b = doc.flight(i, run="preview")["end"], doc.flight(i)["end"]
        return max(abs(a["t"] - b["t"]), *(abs(u - v) for u, v in zip(a["x"], b["x"])),
                   *(abs(u - v) for u, v in zip(a["p"], b["p"])))
    tol_ref = max(TOL, 10 * err, game_gap(0))
    tol_bare = max(TOL, 10 * err_of(rb, rb2), game_gap(1))
    doc.game_value("gapRef", game_gap(0), 2, sci=True)
    doc.game_value("gapBare", game_gap(1), 2, sci=True)
    g = doc.flight(0)
    doc.compare("Outcome", "outcome", run["outcome"], g["outcome"])
    doc.compare("T", r"time of arrival $t_\mathrm{e}$", te, g["end"]["t"], tol_ref)
    doc.compare("X", r"$x_\mathrm{e}$", xe[0], g["end"]["x"][0], tol_ref)
    doc.compare("Y", r"$y_\mathrm{e}$", xe[1], g["end"]["x"][1], tol_ref)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", pe[0], g["end"]["p"][0], tol_ref)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", pe[1], g["end"]["p"][1], tol_ref)
    doc.compare("Rad", "radiated energy (rel.)", E_rad, g["radiated_energy"], TOL_RAD, rel=True, digits=8)
    doc.game_value("gameSteps", g["steps"])
    gb = doc.flight(1)
    doc.compare("BareOut", "empty board: outcome", rb["outcome"], gb["outcome"])
    doc.compare("BareT", r"empty board: $t_\mathrm{e}$", rb["t"], gb["end"]["t"], tol_bare)
    doc.compare("BareX", r"empty board: $x_\mathrm{e}$", rb["x"][0], gb["end"]["x"][0], tol_bare)
    doc.compare("BareY", r"empty board: $y_\mathrm{e}$", rb["x"][1], gb["end"]["x"][1], tol_bare)
