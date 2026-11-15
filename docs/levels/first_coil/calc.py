"""Companion of docs/levels/first_coil/doc.tex (Level 28, First coil).

Every number of the document, from the level file, by the document's method:
- the coil's field in its plane (eq:biot-savart) by the trapezoidal rule on the
  periodic Biot-Savart integral, checked against the closed form with complete elliptic
  integrals (eq:loop-field), whose K and E come from the arithmetic-geometric mean
  (eq:agm);
- the cyclotron radius and period at the launch (eq:larmor-radius, eq:cyclotron);
- the flight by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) on
  eq:eom-b, the entry into the detector by cubic Hermite interpolation and bisection, the
  error by halving the step.
Compared with the game: the reference flight and the empty board.
"""
import math


# The method's own error, by halving the step, is about 1e-12 at the step used; the game's
# verified flights are accurate to about 1e-10 (Foundations, section 6). Ten times the
# larger: 1e-9, in cells and time units. A flight's tolerance is at least the game's own gap
# between its preview and its verification (tolerances 1e-10 and 1e-12: the verified run's
# error is at most a tenth of that gap).
TOL = 1e-9
TOL_RAD = 1e-4
H = 0.005          # the step of the reference computation
H_TABLE = 0.01     # the table of states
N_TRAP = 4000      # nodes of the trapezoidal rule on the Biot-Savart integral


def agm_ke(m):
    """K(m) and E(m), m = k^2, by the arithmetic-geometric mean (eq:agm)."""
    a, b = 1.0, math.sqrt(1.0 - m)
    c2_sum = 0.5 * m          # the sum of 2^(n-1) c_n^2, from c_0^2 = m
    p = 1.0
    for _ in range(40):
        c = 0.5 * (a - b)
        a, b = 0.5 * (a + b), math.sqrt(a * b)
        p *= 2.0
        c2_sum += 0.5 * p * c * c
        if abs(c) < 1e-17:
            break
    K = math.pi / (2.0 * a)
    E = K * (1.0 - c2_sum)
    return K, E


class Setup:
    def __init__(self, lvl, placement):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.a_p = shot["particle"]["radius"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.wire = lvl["physics"]["wire_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        coil = lvl["coils"][0]
        self.cc = (coil["center"][0] / s, coil["center"][1] / s)
        self.ca = coil["radius"] / s
        self.kappa = coil["kappa"]
        elements = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in elements]
        self.kappas = [(self.q * Q, cx, cy) for cx, cy, Q in self.charges]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def bz(self, x, y):
        """eq:loop-field: the coil's B_z in its plane."""
        rho = math.hypot(x - self.cc[0], y - self.cc[1])
        a = self.ca
        alpha2 = (a - rho) ** 2
        beta = a + rho
        m = 4 * a * rho / (beta * beta)
        K, E = agm_ke(m)
        return 2 * self.kappa / (alpha2 * beta) * ((a * a - rho * rho) * E + alpha2 * K)

    def bz_biot_savart(self, rho):
        """eq:biot-savart by the trapezoidal rule (exponentially convergent: periodic)."""
        a = self.ca
        total = 0.0
        for j in range(N_TRAP):
            phi = 2 * math.pi * j / N_TRAP
            cphi = math.cos(phi)
            total += a * (a - rho * cphi) / (a * a + rho * rho - 2 * a * rho * cphi) ** 1.5
        return self.kappa * total * 2 * math.pi / N_TRAP

    def V(self, x, y):
        return sum(k / math.hypot(x - cx, y - cy) for k, cx, cy in self.kappas)

    def rhs(self, s):
        """eq:eom-b."""
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        vx, vy = px / (g * self.m), py / (g * self.m)
        fx = fy = 0.0
        for k, cx, cy in self.kappas:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += k * dx / r3
            fy += k * dy / r3
        b = self.bz(x, y)
        return (vx, vy, fx + self.q * vy * b, fy - self.q * vx * b)

    def energy(self, s):
        x, y, px, py = s
        return math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2) + self.V(x, y)


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


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


def fly(st, h, record_every=None):
    """From the launch to the first event: the detector, a charge's disc, the coil's wire
    (a distance from the wire's circle of at most its radius), the world's edge."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    w0 = st.energy(s)
    rows, path, wdev = [], [(s[0], s[1])], 0.0
    preds = []
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    for idx, (_, cx, cy) in enumerate(st.kappas):
        preds.append((f"Collided({idx})", 0,
                      lambda p, cx=cx, cy=cy: math.hypot(p[0] - cx, p[1] - cy) <= st.R + st.a_p))
    preds.append(("Collided(coil)", 0,
                  lambda p: abs(math.hypot(p[0] - st.cc[0], p[1] - st.cc[1]) - st.ca) <= st.wire + st.a_p))
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                                  and st.wmin[1] <= p[1] <= st.wmax[1])))
    us = [k / 64 for k in range(1, 65)]
    while True:
        # The last step ends exactly at t_max, where the game stops a flight (Timeout).
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
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else (
                face_of(pe, st.wmin, st.wmax) if outcome == "LeftBounds" else "")
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "steps": n + 1, "face": face,
                    "rows": rows, "path": path, "wdev": wdev, "w0": w0}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(st.energy(s) - w0))
        if record_every and n % record_every == 0:
            rows.append((t, s, st.energy(s) - w0))
        if n % 10 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            path.append((s[0], s[1]))
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "rows": rows,
                    "path": path, "steps": n, "wdev": wdev, "w0": w0, "face": ""}


def radiated(st, h, t_end):
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
    doc.value("wire", st.wire, 3)
    doc.value("xZero", st.x0[0], 3)
    doc.value("yZero", st.x0[1], 3)
    doc.value("dirX", st.d[0], 2, sci=True)
    doc.value("xDetLo", st.dmin[0], 3)
    doc.value("yDetLo", st.dmin[1], 3)
    doc.value("xDetHi", st.dmax[0], 3)
    doc.value("yDetHi", st.dmax[1], 3)
    doc.value("xWorld", st.wmax[0], 3)
    doc.value("yWorld", st.wmax[1], 3)
    doc.value("xCoil", st.cc[0], 3)
    doc.value("yCoil", st.cc[1], 3)
    doc.value("aCoil", st.ca, 3)
    doc.value("kappa", st.kappa, 3, sci=True)
    doc.check(len(ref) == 1 and len(lvl["coils"]) == 1, "one coil and one placed charge")
    (qx, qy, Q) = st.charges[0]
    doc.value("QOne", Q, 3, sci=True)
    doc.value("xOne", qx, 3)
    doc.value("yOne", qy, 3)
    doc.value("kOne", st.q * Q, 4)

    # ---- The coil's field (eq:loop-field against eq:biot-savart) -----------------------------
    b_centre = 2 * math.pi * st.kappa / st.ca
    doc.value("Bcentre", b_centre, 6)
    rows = []
    worst = 0.0
    for frac in (0.0, 0.25, 0.5, 0.737, 0.9, 0.97):
        rho = frac * st.ca
        b_closed = st.bz(st.cc[0] + rho, st.cc[1]) if rho > 0 else b_centre
        b_bs = st.bz_biot_savart(rho)
        worst = max(worst, abs(b_closed - b_bs) / abs(b_bs))
        rows.append([rho, b_bs, b_closed, b_bs / b_centre])
    doc.table("fieldtab", [r"$\rho$", r"$B_z$ by \eqref{eq:biot-savart}", r"$B_z$ by \eqref{eq:loop-field}",
                           r"$B_z/B_z(0)$"], rows, digits=[4, 10, 10, 4])
    doc.value("fieldAgree", worst, 2, sci=True)
    doc.check(worst < 1e-12, "the closed form agrees with the Biot-Savart quadrature")
    doc.value("nTrap", N_TRAP)

    # ---- The launch and the circle it would make (eq:larmor-radius, eq:cyclotron) -------------
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p0, 12)
    v0 = st.p0 / (st.gamma0 * st.m)
    doc.value("vZero", v0, 10)
    rho0 = math.hypot(st.x0[0] - st.cc[0], st.x0[1] - st.cc[1])
    b0 = st.bz(*st.x0)
    doc.value("rhoZero", rho0, 4)
    doc.value("Bzero", b0, 8)
    r_l = st.p0 / (abs(st.q) * b0)
    t_c = 2 * math.pi * st.gamma0 * st.m / (abs(st.q) * b0)
    doc.value("rL", r_l, 6)
    doc.value("Tc", t_c, 6)
    doc.value("omegaC", 2 * math.pi / t_c, 6)
    # The sense of rotation: the force q v x B at the launch (v along d, B along z).
    fx, fy = st.q * st.d[1] * b0, -st.q * st.d[0] * b0
    doc.text("sense", "clockwise" if st.q * b0 > 0 else "anticlockwise")
    doc.value("fxZero", fx, 3, sci=True)
    # The centre of the launch circle: from x0 a distance r_L along the force.
    nf = math.hypot(fx, fy)
    gx, gy = st.x0[0] + r_l * fx / nf, st.x0[1] + r_l * fy / nf
    doc.value("xGuide", gx, 5)
    doc.value("yGuide", gy, 5)

    # ---- The empty board ---------------------------------------------------------------------
    bare = Setup(lvl, [])
    rb, rb2 = fly(bare, H), fly(bare, H / 2)
    doc.text("bareOut", rb["outcome"])
    doc.value("bareT", rb["t"], 10)
    doc.value("bareX", rb["x"][0], 10)
    doc.value("bareY", rb["x"][1], 10)
    doc.value("bareErr", err_of(rb, rb2), 2, sci=True)
    doc.value("bareWdev", rb["wdev"], 2, sci=True)
    doc.value("bareTurns", rb["t"] / t_c, 4)
    doc.data("pathbare", rb["path"])

    # ---- The reference flight ----------------------------------------------------------------
    run, half = fly(st, H), fly(st, H / 2)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    err = err_of(run, half)
    doc.value("H", H, 4)
    doc.value("steps", run["steps"])
    te, xe, pe = run["t"], run["x"], run["p"]
    doc.value("tE", te, 12)
    doc.value("xE", xe[0], 12)
    doc.value("yE", xe[1], 12)
    doc.text("faceE", run["face"])
    doc.value("pxE", pe[0], 10)
    doc.value("pyE", pe[1], 10)
    doc.value("err", err, 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    W = run["w0"]
    doc.value("VZero", st.V(*st.x0), 10)
    doc.value("W", W, 12)
    ge = (W - st.V(*xe)) / (st.m * st.c ** 2)
    pe_energy = st.m * st.c * math.sqrt(ge * ge - 1)
    doc.value("pEnergy", pe_energy, 10)
    doc.value("pDiff", abs(math.hypot(*pe) - pe_energy), 2, sci=True)
    doc.check(abs(math.hypot(*pe) - pe_energy) < 1e-9, "|p| at the entry agrees with the energy")
    doc.check(err < 1e-10, "the step-halving estimate is below 1e-10")
    tab = fly(st, H_TABLE, record_every=int(round(5 / H_TABLE)))
    doc.value("HTable", H_TABLE, 3)
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$W(t) - W$"],
              [[f"{t:.0f}", s[0], s[1], s[2], s[3], dw] for t, s, dw in tab["rows"]],
              align="rrrrrr", digits=[None, 8, 8, 8, 8, 2], sci=[False] * 5 + [True])
    E_rad = radiated(st, H, te)
    E_rad2 = radiated(st, 2 * H, te)  # the coarser step: its error is 16 times ours
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / st.T0, 3, sci=True)
    doc.value("EradErr", abs(E_rad - E_rad2) / abs(E_rad2) / 15, 2, sci=True)
    doc.data("path", run["path"])

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
    gb = doc.flight(1)
    doc.compare("BareOut", "empty board: outcome", rb["outcome"], gb["outcome"])
    doc.compare("BareT", r"empty board: $t_\mathrm{e}$", rb["t"], gb["end"]["t"], tol_bare)
    doc.compare("BareX", r"empty board: $x_\mathrm{e}$", rb["x"][0], gb["end"]["x"][0], tol_bare)
    doc.compare("BareY", r"empty board: $y_\mathrm{e}$", rb["x"][1], gb["end"]["x"][1], tol_bare)
