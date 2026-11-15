"""Companion of docs/levels/dempster/doc.tex (Level 33, Dempster's mass spectrometer).

Every number of the document, from the level file: the six ions' flights (three masses, two
launch angles) by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) in the
coil's field (first_coil eq:loop-field, K and E by the arithmetic-geometric mean, eq:agm) and
the placed charge's, the entries by cubic Hermite interpolation and bisection, the errors by
halving the step; the uniform-field estimates (first_coil eq:larmor-radius) at each orbit's
top, and the 180-degree focusing estimate. Compared with the game: the six reference flights
and the empty board's outcomes.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
TOL_RAD = 1e-4
H = 0.005
H_TABLE = 0.01
N_TRAP = 4000


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
    def __init__(self, lvl, placement, i=0):
        shot = lvl["shots"][i]
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


def apex(st, run):
    """The flight's highest sampled point: there the momentum is horizontal, and the circle
    of eq:larmor-radius through it has its centre straight below or above it."""
    return max(run["path"], key=lambda p: p[1])


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    n_shots = len(lvl["shots"])
    doc.check(n_shots == 6 and len(ref) == 1 and len(lvl["coils"]) == 1,
              "six ions, one coil, one placed charge")
    st0 = Setup(lvl, ref, 0)
    (qx, qy, Q) = st0.charges[0]
    doc.value("q", st0.q, 3, sci=True)
    doc.value("c", st0.c, 3)
    doc.value("Tzero", st0.T0, 3)
    doc.value("R", st0.R, 3)
    doc.value("xZero", st0.x0[0], 3)
    doc.value("yZero", st0.x0[1], 3)
    doc.value("angDeg", abs(math.degrees(math.atan2(st0.d[0], st0.d[1]))), 4)
    doc.value("xCoil", st0.cc[0], 3)
    doc.value("yCoil", st0.cc[1], 3)
    doc.value("aCoil", st0.ca, 3)
    doc.value("kappa", st0.kappa, 3, sci=True)
    doc.value("Bcentre", 2 * math.pi * st0.kappa / st0.ca, 6)
    doc.value("QOne", Q, 3, sci=True)
    doc.value("xOne", qx, 3)
    doc.value("yOne", qy, 3)
    doc.value("xWorld", st0.wmax[0], 3)
    doc.value("yWorld", st0.wmax[1], 3)
    reg = lvl["limits"]["region"]
    for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                  ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
        doc.value(k_, v, 3)
    r_launch = math.hypot(st0.x0[0] - qx, st0.x0[1] - qy)
    doc.value("rLaunch", r_launch, 4)
    doc.value("VLaunch", st0.q * Q / r_launch, 6)
    doc.value("H", H, 4)

    rows, est_rows = [], []
    worst_err = worst_w = 0.0
    names = ["a", "b", "c", "d", "e", "f"]
    for i in range(n_shots):
        st = Setup(lvl, ref, i)
        run, half = fly(st, H), fly(st, H / 2)
        doc.check(run["outcome"] == "Arrived", f"shot {i + 1}: the reference flight arrives")
        err = err_of(run, half)
        worst_err = max(worst_err, err)
        worst_w = max(worst_w, run["wdev"])
        ang = math.degrees(math.atan2(st.d[0], st.d[1]))
        rows.append([f"{i + 1}", f"{st.m:g}", f"{ang:+.0f}", run["x"][0], run["t"],
                     f"[{st.dmin[0]:g}, {st.dmax[0]:g}]"])
        doc.data("path" + names[i], run["path"])
        # The uniform-field estimate at the orbit's top (eq:diameter).
        ax, ay = apex(st, run)
        W = run["w0"]
        g_top = (W - st.V(ax, ay)) / (st.m * st.c ** 2)
        p_top = st.m * st.c * math.sqrt(g_top * g_top - 1)
        b_top = st.bz(ax, ay)
        r_top = p_top / abs(st.q * b_top)
        if ang > 0:
            est_rows.append([f"{st.m:g}", p_top, b_top, r_top, 2 * r_top])
        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"ion {i + 1}: "
        doc.compare(f"Out{names[i]}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{names[i]}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{names[i]}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{names[i]}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        if i == 0:
            doc.value("topY", ay, 4)
            doc.value("TtopOne", (g_top - 1) * st.m * st.c ** 2, 5)
    doc.table("results", ["ion", r"$m$", r"angle ($^\circ$)", r"$x_\mathrm{e}$", r"$t_\mathrm{e}$",
                          "detector"], rows, align="rrrrrl", digits=[None, None, None, 8, 8, None])
    doc.table("estimates", [r"$m$", r"$p$ at the top", r"$B_z$ at the top",
                            r"$r = p/\abs{qB_z}$", r"$2r$"], est_rows, digits=[None, 5, 5, 4, 4])
    doc.value("errMax", worst_err, 2, sci=True)
    doc.value("wDevMax", worst_w, 2, sci=True)
    doc.check(worst_err < 1e-10, "the step-halving estimates are below 1e-10")
    # The 180-degree focusing estimate (eq:focus): in a uniform field the two angles land at
    # the same point, both 2 r cos(alpha) from the source; against the central ray, r alpha^2.
    alpha = math.atan2(abs(st0.d[0]), st0.d[1])
    doc.value("alphaRad", alpha, 6)
    doc.value("focusOne", est_rows[0][3] * alpha * alpha, 3)
    spreads = [abs(rows[2 * k][3] - rows[2 * k + 1][3]) for k in range(3)]
    doc.value("spreadOne", spreads[0], 3)
    doc.value("spreadTwo", spreads[1], 3)
    doc.value("spreadFour", spreads[2], 3)

    # The empty board (both outcomes as the game's, compared).
    for i in range(n_shots):
        st = Setup(lvl, [], i)
        rb = fly(st, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{names[i]}", f"empty board, ion {i + 1}: outcome", rb["outcome"], gb["outcome"])
        if i == 4:
            doc.text("bareOutFive", rb["outcome"])
