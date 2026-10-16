"""Companion of docs/levels/slingshot/doc.tex (Level 2, Slingshot).

Every number of the document, from the level file, by the document's method: the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4) on the equations of motion
(eq:eom-two), the entry into the detector by cubic Hermite interpolation of the last step,
and the error estimated by halving the step. Compared with the game: the reference flight
and three more placements in detail.
"""
import math


# The method's own error, estimated by halving the step, is about 1e-11 at the steps used
# (section 5); the game's verified flights are accurate to about 1e-10 (Foundations,
# section 6). Ten times the larger: 1e-9, in cells and time units.
TOL = 1e-9
# The radiated energy: Simpson's rule on the RK4 path against the game's own quadrature of
# the same power (3e-9 relative on an exact test, PHYSICS.md 8); 1e-4 still tells a wrong
# formula apart.
TOL_RAD = 1e-4
H = 0.005          # the step of the reference computation
H_TABLE = 0.01     # the step of the table of states (every 5 time units)


class Setup:
    """The level's numbers (Foundations 5) and the charges of a placement."""

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
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def V(self, x, y):
        """The potential energy q*phi (Foundations eq:coulomb, superposed)."""
        return sum(k / math.hypot(x - cx, y - cy) for k, cx, cy in self.kappas)

    def rhs(self, s):
        """The right-hand side of eq:eom-two."""
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        fx = fy = 0.0
        for k, cx, cy in self.kappas:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += k * dx / r3
            fy += k * dy / r3
        return (px / (g * self.m), py / (g * self.m), fx, fy)

    def energy(self, s):
        """W = gamma m c^2 + V (eq:energy-two)."""
        x, y, px, py = s
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


def fly(st, h, record_every=None):
    """Integrates from the launch until the first event (Foundations 5). Returns a dict with
    the outcome, the state at the event, the steps, and (optionally) a table of states."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    rows, path, wdev = [], [(s[0], s[1])], 0.0
    while True:
        s2 = rk4_step(st, s, f, h)
        f2 = st.rhs(s2)
        # The events of this step, by the interpolant's samples: detector, charges, world.
        events = []
        us = [k / 64 for k in range(1, 65)]
        pts = [hermite(s, f, s2, f2, h, u) for u in us]
        for k, p in enumerate(pts):
            if box_distance(p, st.dmin, st.dmax) <= 0:
                events.append((us[k], "Arrived"))
                break
        for idx, (_, cx, cy) in enumerate(st.kappas):
            for k, p in enumerate(pts):
                if math.hypot(p[0] - cx, p[1] - cy) <= st.R + st.a:
                    events.append((us[k], f"Collided({idx})"))
                    break
        for k, p in enumerate(pts):
            if not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1]):
                events.append((us[k], "LeftBounds"))
                break
        if events:
            u_hit, outcome = min(events)
            if outcome == "Arrived":
                # Bisection on the interpolant: the first instant inside the box.
                lo, hi = u_hit - 1 / 64, u_hit
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    if box_distance(hermite(s, f, s2, f2, h, mid), st.dmin, st.dmax) <= 0:
                        hi = mid
                    else:
                        lo = mid
                te = t + hi * h
                # The state at the entry: one step of the method from t to te.
                se = rk4_step(st, s, f, te - t)
                pe = hermite(s, f, s2, f2, h, hi)
                face = ("left" if abs(pe[0] - st.dmin[0]) < 1e-9 else "right" if abs(pe[0] - st.dmax[0]) < 1e-9
                        else "bottom" if abs(pe[1] - st.dmin[1]) < 1e-9 else "top")
                path.append(pe)
                return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "steps": n + 1,
                        "face": face, "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                        "state_before": s, "t_before": t, "state_after": s2}
            return {"outcome": outcome, "t": t + u_hit * h, "rows": rows, "path": path,
                    "steps": n + 1, "wdev": wdev, "w0": w0}
        s, f, t, n = s2, f2, t + h, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if record_every and n % record_every == 0:
            rows.append((t, s, st.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if t > st.tmax:
            return {"outcome": "Timeout", "t": t, "rows": rows, "path": path, "steps": n,
                    "wdev": wdev, "w0": w0}


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
    # The rest of the flight, past the last even step, by the trapezoidal rule (tiny).
    rest = t_end - n * h
    s2 = rk4_step(st, s, st.rhs(s), rest)
    return simpson + rest * 0.5 * (vals[-1] + power(s2))


def zero_velocity_curve(st, w_bind, origin, n=360):
    """The curve V(x) = W - m c^2 (eq:allowed): along rays from `origin`, the outermost
    point where V rises to the bound energy, by bisection."""
    pts = []
    for k in range(n + 1):
        a = 2 * math.pi * k / n
        dx, dy = math.cos(a), math.sin(a)
        lo, hi = 0.5, 60.0
        if st.V(origin[0] + hi * dx, origin[1] + hi * dy) < w_bind:
            continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if st.V(origin[0] + mid * dx, origin[1] + mid * dy) <= w_bind:
                lo = mid
            else:
                hi = mid
        pts.append((origin[0] + lo * dx, origin[1] + lo * dy))
    return pts


def compute(doc):
    lvl = doc.level
    st = Setup(lvl, lvl["reference_solution"])
    (c1x, c1y, Q1), (c2x, c2y, Q2) = st.charges
    doc.value("q", st.q, 3, sci=True)
    doc.value("m", st.m, 3)
    doc.value("c", st.c, 3)
    doc.value("Tzero", st.T0, 3)
    doc.value("R", st.R, 3)
    doc.value("tmax", st.tmax, 4)
    doc.value("QOne", Q1, 3, sci=True)
    doc.value("QTwo", Q2, 3, sci=True)
    doc.value("kOne", st.q * Q1, 4)
    doc.value("kTwo", st.q * Q2, 4)
    doc.text("rOne", f"({c1x:g}, {c1y:g})")
    doc.text("rTwo", f"({c2x:g}, {c2y:g})")
    doc.value("dirX", st.d[0], 2, sci=True)

    # ---- The launch and the energy (eq:energy-two, eq:allowed) ---------------------------
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p0, 15)
    r10 = math.hypot(st.x0[0] - c1x, st.x0[1] - c1y)
    r20 = math.hypot(st.x0[0] - c2x, st.x0[1] - c2y)
    doc.value("rOneZero", r10, 10)
    doc.value("rTwoZero", r20, 10)
    V0 = st.V(*st.x0)
    doc.value("VZero", V0, 10)
    Wb = st.T0 + V0
    doc.value("Wbind", Wb, 10)
    doc.check(Wb < 0, "the particle is bound: W - mc^2 < 0")
    W = st.gamma0 * st.m * st.c ** 2 + V0
    doc.value("W", W, 15)

    # ---- The reference flight by RK4 (eq:rk4) ------------------------------------------------
    run = fly(st, H)
    half = fly(st, H / 2)
    coarse = fly(st, 2 * H)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    doc.value("H", H, 3)
    doc.value("steps", run["steps"])
    te, xe = run["t"], run["x"]
    doc.value("tE", te, 13)
    doc.value("xE", xe[0], 13)
    doc.value("yE", xe[1], 6)
    doc.text("faceE", run["face"])
    doc.value("errT", abs(half["t"] - run["t"]) / 15, 2, sci=True)
    doc.value("errX", abs(half["x"][0] - run["x"][0]) / 15, 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    pe = run["p"]
    doc.value("pxE", pe[0], 12)
    doc.value("pyE", pe[1], 12)
    # |p| at the entry from the energy alone (eq:energy-two): an independent check.
    ge = (W - st.V(*xe)) / (st.m * st.c ** 2)
    pe_energy = st.m * st.c * math.sqrt(ge * ge - 1)
    doc.value("gammaE", ge, 12)
    doc.value("pEnergy", pe_energy, 12)
    doc.value("pDiff", abs(math.hypot(*pe) - pe_energy), 2, sci=True)
    doc.check(abs(math.hypot(*pe) - pe_energy) < 1e-9, "|p| at the entry agrees with the energy")

    # The table of states every 5 time units (step H_TABLE).
    tab = fly(st, H_TABLE, record_every=int(round(5 / H_TABLE)))
    doc.value("HTable", H_TABLE, 3)
    rows = [[f"{t:.0f}", s[0], s[1], s[2], s[3], dw] for t, s, dw in tab["rows"]]
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$W(t) - W$"], rows,
              align="rrrrrr", digits=[None, 10, 10, 10, 10, 2], sci=[False] * 5 + [True])
    # The last step and the interpolation.
    sb, sa = run["state_before"], run["state_after"]
    doc.value("tBefore", run["t_before"], 6)
    doc.value("xBefore", sb[0], 10)
    doc.value("yBefore", sb[1], 10)
    doc.value("xAfter", sa[0], 10)
    doc.value("yAfter", sa[1], 10)
    doc.value("uE", (te - run["t_before"]) / H, 10)

    # The convergence with the step (eq:rk4's error, of order h^4).
    conv = []
    for h in (0.04, 0.02, 0.01, 0.005, 0.0025):
        r = fly(st, h)
        conv.append([f"{h:g}", r["steps"], r["t"], r["x"][0]])
    doc.table("convergence", [r"$h$", "steps", r"$t_\mathrm{e}$", r"$x_\mathrm{e}$"], conv,
              align="rrll", digits=[None, None, 13, 13])
    doc.value("coarseErrT", abs(coarse["t"] - run["t"]), 2, sci=True)
    rough = fly(st, 0.04)
    doc.value("roughSteps", rough["steps"])
    doc.value("roughErr", max(abs(rough["t"] - run["t"]), abs(rough["x"][0] - run["x"][0])), 1, sci=True)
    doc.value("tableSteps", fly(st, H_TABLE)["steps"])

    # The apex: where the particle comes closest to standing still.
    st_apex = min(tab["rows"], key=lambda r: math.hypot(r[1][2], r[1][3]))

    # ---- The radiated energy -------------------------------------------------------------------
    E_rad = radiated(st, H, te)
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / st.T0, 3, sci=True)

    # ---- The figure: the path and the zero-velocity curve (eq:allowed) ---------------------------
    doc.data("path", run["path"])
    doc.data("zvc", zero_velocity_curve(st, Wb, (6.0, 6.0)))

    # ---- The comparison with the game -----------------------------------------------------------
    g = doc.flight(0)
    doc.compare("Outcome", "outcome", "Arrived", g["outcome"])
    doc.compare("T", r"time of arrival $t_\mathrm{e}$", te, g["end"]["t"], TOL)
    doc.compare("X", r"$x_\mathrm{e}$", xe[0], g["end"]["x"][0], TOL)
    doc.compare("Y", r"$y_\mathrm{e}$", xe[1], g["end"]["x"][1], TOL)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", pe[0], g["end"]["p"][0], TOL)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", pe[1], g["end"]["p"][1], TOL)
    doc.compare("Rad", "radiated energy (rel.)", E_rad, g["radiated_energy"], TOL_RAD, rel=True, digits=8)
    doc.game_value("gameSteps", g["steps"])

    # ---- Three more placements in detail (placements.json) ----------------------------------------
    for k in range(1, doc.n_results()):
        placement = doc.result(k)["placement"]
        e = placement[0]
        tag = ["", "One", "Two", "Three", "Four"][k]
        sk = Setup(lvl, placement)
        rk, rk2 = fly(sk, H), fly(sk, H / 2)
        gk = doc.flight(k)
        doc.text(f"det{tag}Node", f"({e['node'][0]}, {e['node'][1]})")
        doc.value(f"det{tag}Q", e["value"], 3, sci=True)
        doc.compare(f"Det{tag}Out", rf"placement {k}: outcome", rk["outcome"], gk["outcome"])
        if rk["outcome"] != "Arrived":
            continue
        doc.text(f"det{tag}Face", rk["face"])
        doc.value(f"det{tag}T", rk["t"], 12)
        doc.value(f"det{tag}X", rk["x"][0], 12)
        doc.value(f"det{tag}Y", rk["x"][1], 12)
        doc.value(f"det{tag}Err", max(abs(rk2["t"] - rk["t"]), abs(rk2["x"][0] - rk["x"][0]),
                                      abs(rk2["x"][1] - rk["x"][1])) / 15, 2, sci=True)
        doc.compare(f"Det{tag}T", rf"placement {k}: $t_\mathrm{{e}}$", rk["t"], gk["end"]["t"], TOL)
        doc.compare(f"Det{tag}X", rf"placement {k}: $x_\mathrm{{e}}$", rk["x"][0], gk["end"]["x"][0], TOL)
        doc.compare(f"Det{tag}Y", rf"placement {k}: $y_\mathrm{{e}}$", rk["x"][1], gk["end"]["x"][1], TOL)

    # ---- The game's map of single-charge solutions (for orientation; the game's verdicts) ------------
    verdicts = doc.verdicts()
    fam = doc.family()
    names = {1e6: "plusone", 2e6: "plustwo", 4e6: "plusfour",
             -1e6: "minusone", -2e6: "minustwo", -4e6: "minusfour"}
    counts = {v: [] for v in names}
    for p, code in zip(fam, verdicts):
        if code == "A":
            e = p[0]
            counts[e["value"]].append((e["node"][0], e["node"][1]))
    for v, name in names.items():
        doc.game_value(f"n{name}", len(counts[v]))
        doc.data(f"map_{name}", counts[v] or [(-5, -5)])
    doc.game_value("nSolving", sum(len(v) for v in counts.values()))
    doc.game_value("nFamily", len(fam))
