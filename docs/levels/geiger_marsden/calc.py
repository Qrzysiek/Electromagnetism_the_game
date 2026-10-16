"""Companion of docs/levels/geiger_marsden/doc.tex (Level 3, Geiger–Marsden).

The reference flight has two centres, so it has no closed form (slingshot). It is integrated
by the classical Runge–Kutta method at a fixed step (Foundations eq:rk4) on eq:eom-two, and
the entry into the detector is located by cubic Hermite interpolation of the last step. The
error is estimated by halving the step, before any comparison with the game.

The level's Rutherford formula is derived in the document from the Newtonian limit of
Level 1's orbit and evaluated here two ways (eq:chi-tan and eq:chi-geom), as a check of that
algebra. It is not the number compared with the game: the flight is relativistic and has two
centres. A circle about the nucleus separates the encounter (eq:sphere) so the formula, the
one-centre orbit of Level 1, and the two-centre flight can be set side by side.
"""
import math


# Step-halving puts the truncation error near 1e-12 (section 5). The game's verified
# flights are accurate to about 1e-10 (Foundations, section 6). The tolerance is ten
# times the larger of those, and it is not taken from the difference with the game.
TOL = 1e-9
# Simpson's rule on the radiated power against the game's own quadrature of the same
# power (3e-9 relative on an exact test, PHYSICS.md 8). 1e-4 still separates a wrong
# formula. Set before comparing.
TOL_RAD = 1e-4
H = 0.01
# A circle about the nucleus on which the encounter is read off (eq:sphere). Several
# times the distance of closest approach, and not a number from the level file.
R_SPHERE = 8.0


class Setup:
    """The level's numbers and the charges of one placement, in the game's order."""

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
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        self.kappas = [(self.q * Q, cx, cy) for cx, cy, Q in self.charges]
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2.0 * self.m * self.T0)

    def V(self, x, y):
        """Potential energy q*phi, both charges (Foundations eq:coulomb)."""
        return sum(k / math.hypot(x - cx, y - cy) for k, cx, cy in self.kappas)

    def rhs(self, state):
        """Right-hand side of eq:eom-two."""
        x, y, px, py = state
        g = math.sqrt(1.0 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        fx = fy = 0.0
        for k, cx, cy in self.kappas:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += k * dx / r3
            fy += k * dy / r3
        return (px / (g * self.m), py / (g * self.m), fx, fy)

    def energy(self, state):
        """W = gamma m c^2 + V (eq:energy-two)."""
        x, y, px, py = state
        kinetic = math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2)
        return kinetic + self.V(x, y)

    def launch_state(self):
        return (self.x0[0], self.x0[1], self.p0 * self.d[0], self.p0 * self.d[1])


def rk4_step(st, s, f, h):
    """One step of the classical Runge–Kutta method (Foundations eq:rk4)."""
    k1 = f
    k2 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)))
    k3 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)))
    k4 = st.rhs(tuple(a + h * b for a, b in zip(s, k3)))
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
                 for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    """Cubic Hermite interpolant of the position on a step, at fraction u (eq:hermite)."""
    h00 = 2 * u ** 3 - 3 * u ** 2 + 1
    h10 = u ** 3 - 2 * u ** 2 + u
    h01 = -2 * u ** 3 + 3 * u ** 2
    h11 = u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def box_distance(p, lo, hi):
    """Signed distance to a box, per axis (negative inside)."""
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def _bisect(pred, s, f, s2, f2, h, lo, hi):
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pred(hermite(s, f, s2, f2, h, mid)):
            hi = mid
        else:
            lo = mid
    return hi


def _first_u(pred, pts, us):
    for u, p in zip(us, pts):
        if pred(p):
            return u
    return None


def fly(st, h, record_every=None, sphere=None):
    """Integrate from the launch until the first event (Foundations, section 5).

    Events are ranked as the board ranks them: a collision, then the world's edge, then
    the detector. `sphere`, if given, is (cx, cy, radius): the first inward and outward
    crossings of that circle are recorded with the same interpolant.
    """
    s = st.launch_state()
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    rows = []
    path = [(s[0], s[1])]
    wdev = 0.0
    closest = None
    sphere_in = sphere_out = None
    us = [k / 64 for k in range(1, 65)]
    while True:
        s2 = rk4_step(st, s, f, h)
        f2 = st.rhs(s2)
        pts = [hermite(s, f, s2, f2, h, u) for u in us]
        hits = []
        for idx, (_, cx, cy) in enumerate(st.kappas):
            limit = st.R + st.a
            u = _first_u(lambda p, cx=cx, cy=cy, limit=limit:
                         math.hypot(p[0] - cx, p[1] - cy) <= limit, pts, us)
            if u is not None:
                hits.append((u, 0, f"Collided({idx})"))
        u = _first_u(lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                    and st.wmin[1] <= p[1] <= st.wmax[1]), pts, us)
        if u is not None:
            hits.append((u, 1, "LeftBounds"))
        u = _first_u(lambda p: box_distance(p, st.dmin, st.dmax) <= 0, pts, us)
        if u is not None:
            hits.append((u, 2, "Arrived"))
        # Closest approach to the nucleus (the first charge), refined on this step.
        cx, cy = st.charges[0][0], st.charges[0][1]
        ru, rr, rp = _min_radius(s, f, s2, f2, h, cx, cy)
        if closest is None or rr < closest["r"]:
            te = t + ru * h
            se = s if ru == 0 else rk4_step(st, s, f, te - t)
            closest = {"t": te, "x": rp[0], "y": rp[1], "r": rr, "p": se[2:]}
        if sphere is not None and (sphere_in is None or sphere_out is None):
            scx, scy, srad = sphere
            def rad(p, scx=scx, scy=scy):
                return math.hypot(p[0] - scx, p[1] - scy)
            r_start = rad((s[0], s[1]))
            if sphere_in is None and r_start > srad:
                u_in = _first_u(lambda p: rad(p) <= srad, pts, us)
                if u_in is not None:
                    lo = u_in - 1 / 64
                    hi = _bisect(lambda p: rad(p) <= srad, s, f, s2, f2, h, lo, u_in)
                    sphere_in = _crossing(st, s, f, s2, f2, h, t, hi, scx, scy)
            if sphere_in is not None and sphere_out is None and r_start < srad:
                u_out = _first_u(lambda p: rad(p) >= srad, pts, us)
                if u_out is not None:
                    lo = u_out - 1 / 64
                    hi = _bisect(lambda p: rad(p) >= srad, s, f, s2, f2, h, lo, u_out)
                    sphere_out = _crossing(st, s, f, s2, f2, h, t, hi, scx, scy)
        if hits:
            u_hit, _, outcome = min(hits)
            if outcome == "Arrived":
                lo = u_hit - 1 / 64
                hi = _bisect(lambda p: box_distance(p, st.dmin, st.dmax) <= 0,
                             s, f, s2, f2, h, lo, u_hit)
            elif outcome == "LeftBounds":
                lo = u_hit - 1 / 64
                hi = _bisect(lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                            and st.wmin[1] <= p[1] <= st.wmax[1]),
                             s, f, s2, f2, h, lo, u_hit)
            else:
                hi = u_hit
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "steps": n + 1,
                    "rows": rows, "path": path, "wdev": wdev, "w0": w0,
                    "state_before": s, "t_before": t, "state_after": s2,
                    "closest": closest, "sphere_in": sphere_in, "sphere_out": sphere_out}
        s, f, t, n = s2, f2, t + h, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if record_every and n % record_every == 0:
            rows.append((t, s, st.energy(s) - w0))
        if n % 15 == 0:
            path.append((s[0], s[1]))
        if t > st.tmax:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:],
                    "rows": rows, "path": path, "steps": n, "wdev": wdev, "w0": w0,
                    "state_before": s, "t_before": t, "state_after": s,
                    "closest": closest, "sphere_in": sphere_in, "sphere_out": sphere_out}


def _min_radius(s, f, s2, f2, h, cx, cy):
    def at(u):
        p = hermite(s, f, s2, f2, h, u)
        return math.hypot(p[0] - cx, p[1] - cy), p
    best_u, best_r, _ = 0.0, at(0.0)[0], None
    for k in range(1, 65):
        u = k / 64
        r, _ = at(u)
        if r < best_r:
            best_u, best_r = u, r
    lo = max(0.0, best_u - 1 / 64)
    hi = min(1.0, best_u + 1 / 64)
    gr = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    for _ in range(40):
        u1 = b - gr * (b - a)
        u2 = a + gr * (b - a)
        if at(u1)[0] < at(u2)[0]:
            b = u2
        else:
            a = u1
    u = 0.5 * (a + b)
    r, p = at(u)
    return u, r, p


def _crossing(st, s, f, s2, f2, h, t, u, cx, cy):
    te = t + u * h
    se = rk4_step(st, s, f, te - t)
    pe = hermite(s, f, s2, f2, h, u)
    return {"t": te, "x": pe[0], "y": pe[1], "p": se[2:]}


def kinetic(st, p):
    """T = (gamma - 1) m c^2."""
    g = math.sqrt(1.0 + (p[0] ** 2 + p[1] ** 2) / (st.m ** 2 * st.c ** 2))
    return (g - 1.0) * st.m * st.c ** 2, g


def angle_deg(p):
    return math.degrees(math.atan2(p[1], p[0]))


def radiated(st, h, t_end):
    """Radiated energy along the RK4 path (Foundations eq:larmor, Level 1 eq:acceleration),
    by Simpson's rule on the steps (Foundations eq:simpson)."""
    def power(state):
        x, y, px, py = state
        g = math.sqrt(1.0 + (px * px + py * py) / (st.m ** 2 * st.c ** 2))
        v = (px / (g * st.m), py / (g * st.m))
        _, _, fx, fy = st.rhs(state)
        vF = v[0] * fx + v[1] * fy
        ax = (fx - v[0] * vF / st.c ** 2) / (g * st.m)
        ay = (fy - v[1] * vF / st.c ** 2) / (g * st.m)
        vxa = v[0] * ay - v[1] * ax
        return (2.0 / 3.0) * st.q ** 2 / st.c ** 3 * g ** 6 * (
            ax * ax + ay * ay - vxa * vxa / st.c ** 2)

    n = int(t_end / h)
    n -= n % 2
    state = st.launch_state()
    vals = [power(state)]
    for _ in range(n):
        state = rk4_step(st, state, st.rhs(state), h)
        vals.append(power(state))
    simpson = h / 3 * (vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2]))
    rest = t_end - n * h
    state2 = rk4_step(st, state, st.rhs(state), rest)
    return simpson + rest * 0.5 * (vals[-1] + power(state2))


def rutherford(kappa, T, b):
    """The Newtonian scattering angle from infinity, two ways (eq:chi-tan, eq:chi-geom).

    Both return chi in radians. `b` is positive.
    """
    alpha = kappa / (2.0 * T * b * b)
    A = math.sqrt(1.0 / (b * b) + alpha * alpha)
    theta_inf = math.acos(max(-1.0, min(1.0, alpha / A)))
    chi_geom = math.pi - 2.0 * theta_inf
    chi_tan = 2.0 * math.atan(kappa / (2.0 * T * b))
    return chi_geom, chi_tan, alpha, A, theta_inf


def one_centre(st, cross, kappa, cx, cy):
    """Level 1's orbit (first_bend eq:general) continued from a state on the circle.

    Returns the angle swept to the other crossing of the same radius, the change in the
    direction of p between the crossings, and the scattering angle from asymptote to
    asymptote of this same orbit.
    """
    rx, ry = cross["x"] - cx, cross["y"] - cy
    r0 = math.hypot(rx, ry)
    phi0 = math.atan2(ry, rx)
    px, py = cross["p"]
    Lz = rx * py - ry * px
    sigma = 1.0 if Lz >= 0.0 else -1.0
    L = abs(Lz)
    pr = (px * rx + py * ry) / r0
    u0 = 1.0 / r0
    up0 = -pr / L
    T, g = kinetic(st, cross["p"])
    W = g * st.m * st.c ** 2 + kappa / r0
    gamma2 = 1.0 - kappa ** 2 / (st.c ** 2 * L * L)
    gamma = math.sqrt(gamma2)
    B = -kappa * W / (st.c ** 2 * L * L * gamma2)

    def u(theta):
        return (B + (u0 - B) * math.cos(gamma * theta)
                + (up0 / gamma) * math.sin(gamma * theta))

    def uprime(theta):
        return (-(u0 - B) * gamma * math.sin(gamma * theta)
                + up0 * math.cos(gamma * theta))

    def root_down(th0, step, level):
        """First theta past th0, moving by `step`, where u falls through `level`.

        The bracket is then bisected. `step` may be negative: the incoming asymptote
        lies at negative theta.
        """
        th = th0
        up = u(th)
        guard = 0
        while abs(th) < 6.0 and guard < 20000:
            guard += 1
            th2 = th + step
            u2 = u(th2)
            if up > level >= u2:
                inside, outside = th, th2
                for _ in range(70):
                    mid = 0.5 * (inside + outside)
                    if u(mid) > level:
                        inside = mid
                    else:
                        outside = mid
                return 0.5 * (inside + outside)
            th, up = th2, u2
        return None

    # Out of the circle: u rises above u0, then falls back through it.
    th_exit = root_down(1e-3, 1e-3, u0)
    # The two asymptotes of the same orbit, where u falls through 0.
    th_out = root_down(0.0, 1e-3, 0.0)
    th_in = root_down(0.0, -1e-3, 0.0)

    def momentum(theta):
        pr_ = -L * uprime(theta)
        pp_ = sigma * L * u(theta)
        phi = phi0 + sigma * theta
        rhat = (math.cos(phi), math.sin(phi))
        phihat = (-math.sin(phi), math.cos(phi))
        return (pr_ * rhat[0] + pp_ * phihat[0], pr_ * rhat[1] + pp_ * phihat[1])

    p_exit = momentum(th_exit)
    chi = math.pi - (th_out - th_in)
    return {
        "theta_exit": th_exit,
        "ang_exit": angle_deg(p_exit),
        "turn": angle_deg(p_exit) - angle_deg(cross["p"]),
        "chi": chi,
        "gamma2": gamma2,
        "L": L,
        "W": W,
        "T": T,
        "th_in": th_in,
        "th_out": th_out,
    }


def circle(cx, cy, radius, n=180):
    return [(cx + radius * math.cos(2 * math.pi * k / n),
             cy + radius * math.sin(2 * math.pi * k / n)) for k in range(n + 1)]


def compute(doc):
    lvl = doc.level
    st = Setup(lvl, lvl["reference_solution"])
    (nx, ny, Qn), (sx, sy, Qs) = st.charges
    kN, kS = st.q * Qn, st.q * Qs
    doc.value("q", st.q, 3, sci=True)
    doc.value("m", st.m, 3)
    doc.value("c", st.c, 3)
    doc.value("Tzero", st.T0, 3)
    doc.value("R", st.R, 3)
    doc.value("tmax", st.tmax, 4)
    doc.value("Qn", Qn, 3, sci=True)
    doc.value("Qs", Qs, 3, sci=True)
    doc.value("kN", kN, 4)
    doc.value("kS", kS, 4)
    doc.text("rN", f"({nx:g}, {ny:g})")
    doc.text("rS", f"({sx:g}, {sy:g})")
    doc.value("dirX", st.d[0], 2)
    doc.value("gammaZero", st.gamma0, 12)
    doc.value("pZero", st.p0, 12)

    rN0 = math.hypot(st.x0[0] - nx, st.x0[1] - ny)
    rS0 = math.hypot(st.x0[0] - sx, st.x0[1] - sy)
    doc.value("rNZero", rN0, 10)
    doc.value("rSZero", rS0, 10)
    V0 = st.V(*st.x0)
    doc.value("VZero", V0, 10)
    doc.value("Wbind", st.T0 + V0, 10)
    W = st.gamma0 * st.m * st.c ** 2 + V0
    doc.value("W", W, 12)
    # A head-on Newtonian alpha of energy T0 stops where kappa/r = T0 (eq:headon).
    doc.value("headOn", kN / st.T0, 6)

    # ---- The reference flight (eq:rk4, the method of slingshot) ----------------------------
    sphere = (nx, ny, R_SPHERE)
    run = fly(st, H, record_every=int(round(2.0 / H)), sphere=sphere)
    half = fly(st, H / 2, sphere=sphere)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    doc.check(run["sphere_in"] is not None and run["sphere_out"] is not None,
              "the flight crosses the circle about the nucleus twice")
    doc.value("H", H, 3)
    doc.value("steps", run["steps"])
    doc.value("Rsphere", R_SPHERE, 3)
    te, xe, pe = run["t"], run["x"], run["p"]
    doc.value("tE", te, 12)
    doc.value("xE", xe[0], 12)
    doc.value("yE", xe[1], 12)
    doc.value("pxE", pe[0], 12)
    doc.value("pyE", pe[1], 12)
    doc.value("angE", angle_deg(pe), 6)
    TE, gE = kinetic(st, pe)
    doc.value("TE", TE, 8)
    err = max(abs(half["t"] - run["t"]),
              abs(half["x"][0] - run["x"][0]),
              abs(half["x"][1] - run["x"][1]),
              abs(half["p"][0] - run["p"][0]),
              abs(half["p"][1] - run["p"][1])) / 15
    doc.value("errStep", err, 2, sci=True)
    doc.check(err < 1e-10, "halving the step moves the entry by far less than 1e-9")
    doc.value("errT", abs(half["t"] - run["t"]) / 15, 2, sci=True)
    doc.value("errX", abs(half["x"][0] - run["x"][0]) / 15, 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    ge = (W - st.V(*xe)) / (st.m * st.c ** 2)
    p_energy = st.m * st.c * math.sqrt(ge * ge - 1)
    doc.value("gammaE", ge, 12)
    doc.value("pEnergy", p_energy, 12)
    doc.value("pDiff", abs(math.hypot(*pe) - p_energy), 2, sci=True)
    doc.check(abs(math.hypot(*pe) - p_energy) < 1e-8,
              "the momentum at the entry agrees with the energy")
    sb, sa = run["state_before"], run["state_after"]
    doc.value("tBefore", run["t_before"], 6)
    doc.value("xBefore", sb[0], 4)
    doc.value("yBefore", sb[1], 4)
    doc.value("pxBefore", sb[2], 4)
    doc.value("pyBefore", sb[3], 4)
    doc.value("xAfter", sa[0], 4)
    doc.value("yAfter", sa[1], 4)
    doc.value("uE", (te - run["t_before"]) / H, 8)

    rows = [[f"{t:.0f}", s[0], s[1], s[2], s[3], dw] for t, s, dw in run["rows"]]
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$W(t)-W$"],
              rows, align="rrrrrr", digits=[None, 6, 6, 6, 6, 2], sci=[False] * 5 + [True])

    conv = []
    for h in (0.04, 0.02, 0.01, 0.005):
        r = fly(st, h)
        conv.append([f"{h:g}", r["steps"], r["t"], r["x"][1]])
    doc.table("convergence", [r"$h$", "steps", r"$t_\mathrm{e}$", r"$y_\mathrm{e}$"],
              conv, align="rrll", digits=[None, None, 10, 10])
    rough = fly(st, 0.04)
    doc.value("roughSteps", rough["steps"])
    doc.value("roughErr", max(abs(rough["t"] - run["t"]), abs(rough["x"][1] - run["x"][1])),
              2, sci=True)

    # ---- Closest approach and the two forces there (eq:forces) -----------------------------
    cl = run["closest"]
    doc.value("rMin", cl["r"], 8)
    doc.value("tMin", cl["t"], 8)
    doc.value("xMin", cl["x"], 8)
    doc.value("yMin", cl["y"], 8)
    Tmin, _ = kinetic(st, cl["p"])
    doc.value("Tmin", Tmin, 8)
    fN = abs(kN) / cl["r"] ** 2
    rS = math.hypot(cl["x"] - sx, cl["y"] - sy)
    fS = abs(kS) / rS ** 2
    doc.value("rSteerMin", rS, 6)
    doc.value("forceRatio", fS / fN, 4, sci=True)
    # Energy in the nucleus's potential alone, W - mc^2 - kappa_S/r_S (eq:channel).
    e_launch = (st.T0 + V0) - kS / rS0
    e_min = Tmin + kN / cl["r"]
    doc.value("EnucLaunch", e_launch, 6)
    doc.value("EnucMin", e_min, 6)
    doc.value("rStop", kN / e_min, 6)
    doc.check(cl["r"] > st.R + st.a, "the alpha misses the nucleus")
    doc.check(fS / fN < 0.02, "at the turn the nucleus dominates the force")

    # ---- The circle: local Rutherford numbers, Level 1's orbit, the two-centre turn --------
    inn, out = run["sphere_in"], run["sphere_out"]
    doc.value("tIn", inn["t"], 8)
    doc.value("xIn", inn["x"], 8)
    doc.value("yIn", inn["y"], 8)
    Tin, _ = kinetic(st, inn["p"])
    doc.value("Tin", Tin, 8)
    doc.value("angIn", angle_deg(inn["p"]), 6)
    doc.value("pxIn", inn["p"][0], 6)
    doc.value("pyIn", inn["p"][1], 6)
    rx, ry = inn["x"] - nx, inn["y"] - ny
    Lz = rx * inn["p"][1] - ry * inn["p"][0]
    b_local = abs(Lz) / math.hypot(*inn["p"])
    doc.value("bIn", b_local, 8)
    doc.value("Lin", abs(Lz), 8)
    # The level's formula with the kinetic energy and the perpendicular distance read off
    # on the circle (eq:chi-tan). Both expressions of the derivation, as a check.
    chi_geom, chi_tan, alpha, A, theta_inf = rutherford(kN, Tin, b_local)
    doc.value("chiTan", math.degrees(chi_tan), 6)
    doc.value("chiGeom", math.degrees(chi_geom), 6)
    doc.value("chiAlg", abs(chi_tan - chi_geom), 2, sci=True)
    doc.check(abs(chi_tan - chi_geom) < 1e-12,
              "the two forms of the Rutherford angle agree")
    doc.value("alphaIn", alpha, 8)
    doc.value("AIn", A, 8)
    doc.value("thetaInf", math.degrees(theta_inf), 6)
    # How far the circle is from the formula's assumption of infinity.
    doc.value("potIn", kN / R_SPHERE, 6)
    doc.value("potFrac", (kN / R_SPHERE) / Tin, 4)

    oc = one_centre(st, inn, kN, nx, ny)
    doc.check(oc["theta_exit"] is not None and oc["th_in"] is not None
              and oc["th_out"] is not None, "the one-centre orbit reaches both asymptotes")
    doc.value("thetaExit", oc["theta_exit"], 8)
    doc.value("angOne", oc["ang_exit"], 6)
    doc.value("turnOne", oc["turn"], 6)
    doc.value("chiRel", math.degrees(oc["chi"]), 6)
    doc.value("gammaTwo", oc["gamma2"], 8)
    doc.value("eta", kN / (st.c * oc["L"]), 6)
    turn_two = angle_deg(out["p"]) - angle_deg(inn["p"])
    doc.value("tOut", out["t"], 8)
    doc.value("xOut", out["x"], 6)
    doc.value("yOut", out["y"], 6)
    doc.value("angOut", angle_deg(out["p"]), 6)
    doc.value("turnTwo", turn_two, 6)
    doc.value("turnGap", abs(turn_two - oc["turn"]), 6)
    # The same Newtonian formula with the energy and impact parameter the one-centre
    # problem actually carries to infinity (eq:chi-inf): T_inf = T + kappa/r, and
    # b_inf from the Newtonian angular momentum m r x v = L_rel/gamma.
    Tinf = Tin + kN / R_SPHERE
    g_in = kinetic(st, inn["p"])[1]
    L_newt = abs(Lz) / g_in
    p_inf = math.sqrt(2.0 * st.m * Tinf)
    b_inf = L_newt / p_inf
    _, chi_inf, _, _, _ = rutherford(kN, Tinf, b_inf)
    doc.value("Tinf", Tinf, 8)
    doc.value("bInf", b_inf, 8)
    doc.value("chiInf", math.degrees(chi_inf), 6)
    doc.value("chiRelGap", abs(math.degrees(oc["chi"]) - math.degrees(chi_inf)), 6)

    # ---- Radiation (as:radiation) ----------------------------------------------------------
    E_rad = radiated(st, H, te)
    E_rad_fine = radiated(st, H / 2, te)
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / st.T0, 3, sci=True)
    doc.value("radErr", abs(E_rad_fine - E_rad) / E_rad, 2, sci=True)
    doc.check(abs(E_rad_fine - E_rad) / E_rad < TOL_RAD / 10,
              "the radiation quadrature is finer than its comparison tolerance")

    doc.data("path", run["path"])
    doc.data("sphere", circle(nx, ny, R_SPHERE))

    # ---- Comparison with the game ----------------------------------------------------------
    g = doc.flight(0)
    doc.compare("Outcome", "outcome", "Arrived", g["outcome"])
    doc.compare("T", r"time of arrival $t_\mathrm{e}$", te, g["end"]["t"], TOL)
    doc.compare("X", r"$x_\mathrm{e}$", xe[0], g["end"]["x"][0], TOL)
    doc.compare("Y", r"$y_\mathrm{e}$", xe[1], g["end"]["x"][1], TOL)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", pe[0], g["end"]["p"][0], TOL)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", pe[1], g["end"]["p"][1], TOL)
    doc.compare("Rad", "radiated energy (rel.)", E_rad, g["radiated_energy"], TOL_RAD, rel=True)
    doc.game_value("gameSteps", g["steps"])

    # ---- The weaker steering charge (placements.json): a backscatter that misses -----------
    placement = doc.result(1)["placement"]
    e = placement[0]
    sk = Setup(lvl, placement)
    rk = fly(sk, H)
    half_k = fly(sk, H / 2)
    gk = doc.flight(1)
    doc.text("detOneNode", f"({e['node'][0]}, {e['node'][1]})")
    doc.value("detOneQ", e["value"], 3, sci=True)
    doc.value("detOneK", st.q * e["value"], 4)
    doc.compare("DetOneOut", "weaker charge: outcome", rk["outcome"], gk["outcome"])
    doc.value("weakT", rk["t"], 8)
    doc.value("weakX", rk["x"][0], 8)
    doc.value("weakY", rk["x"][1], 8)
    doc.value("weakAng", angle_deg(rk["p"]), 6)
    doc.value("weakRmin", rk["closest"]["r"], 6)
    err_w = max(abs(half_k["t"] - rk["t"]),
                abs(half_k["x"][0] - rk["x"][0]),
                abs(half_k["x"][1] - rk["x"][1])) / 15
    doc.value("weakErr", err_w, 2, sci=True)
    doc.check(rk["outcome"] == "LeftBounds", "the weaker charge leaves the world")
    if rk["outcome"] == "LeftBounds":
        doc.compare("WeakT", r"weaker charge: $t$", rk["t"], gk["end"]["t"], TOL)
        doc.compare("WeakX", r"weaker charge: $x$", rk["x"][0], gk["end"]["x"][0], TOL)
        doc.compare("WeakY", r"weaker charge: $y$", rk["x"][1], gk["end"]["x"][1], TOL)
    doc.data("pathweak", rk["path"])

    # ---- The nucleus alone: the formula's own problem, from this launch (eq:chi-tan) ----------
    bare = Setup(lvl, [])
    rb = fly(bare, H)
    half_b = fly(bare, H / 2)
    # Perpendicular distance from the nucleus to the launch velocity.
    rx0, ry0 = st.x0[0] - nx, st.x0[1] - ny
    p_launch = (st.p0 * st.d[0], st.p0 * st.d[1])
    b_launch = abs(rx0 * p_launch[1] - ry0 * p_launch[0]) / st.p0
    _, chi_naive, _, _, _ = rutherford(kN, st.T0, b_launch)
    oc_launch = one_centre(st, {"x": st.x0[0], "y": st.x0[1], "p": p_launch}, kN, nx, ny)
    T_loc, g_loc = kinetic(st, p_launch)
    Tinf_b = T_loc + kN / rN0
    b_inf_b = (b_launch * st.p0 / g_loc) / math.sqrt(2.0 * st.m * Tinf_b)
    _, chi_inf_b, _, _, _ = rutherford(kN, Tinf_b, b_inf_b)
    doc.value("bLaunch", b_launch, 6)
    doc.value("chiNaive", math.degrees(chi_naive), 6)
    doc.value("TbareInf", Tinf_b, 6)
    doc.value("bBareInf", b_inf_b, 6)
    doc.value("chiBareInf", math.degrees(chi_inf_b), 6)
    doc.value("chiBare", math.degrees(oc_launch["chi"]), 6)
    doc.value("chiBareGap", abs(math.degrees(oc_launch["chi"]) - math.degrees(chi_inf_b)), 4)
    doc.value("chiReadGap", abs(math.degrees(chi_naive) - math.degrees(chi_inf_b)), 4)
    doc.value("bareT", rb["t"], 6)
    doc.value("bareX", rb["x"][0], 6)
    doc.value("bareY", rb["x"][1], 6)
    doc.value("bareAng", angle_deg(rb["p"]), 6)
    doc.value("bareRmin", rb["closest"]["r"], 6)
    doc.text("bareOut", rb["outcome"])
    err_b = max(abs(half_b["t"] - rb["t"]), abs(half_b["x"][0] - rb["x"][0]),
                abs(half_b["x"][1] - rb["x"][1])) / 15
    doc.value("bareErr", err_b, 2, sci=True)
    doc.check(rb["outcome"] == "LeftBounds", "with only the nucleus the alpha leaves the world")
    doc.check(oc_launch["th_out"] is not None, "the unsteered orbit has an outgoing asymptote")
    gb = doc.flight(2)
    doc.compare("BareOut", "nucleus only: outcome", rb["outcome"], gb["outcome"])
    doc.compare("BareT", r"nucleus only: $t$", rb["t"], gb["end"]["t"], TOL)
    doc.compare("BareX", r"nucleus only: $x$", rb["x"][0], gb["end"]["x"][0], TOL)
    doc.compare("BareY", r"nucleus only: $y$", rb["x"][1], gb["end"]["x"][1], TOL)
    doc.data("pathbare", rb["path"])
