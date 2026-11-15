"""Companion of docs/levels/calutron/doc.tex (Level 35, Calutron).

Every number of the document, from the level file: the three isotopes' flights by the
classical Runge-Kutta method at a fixed step (Foundations eq:rk4) in the placed magnets'
in-plane dipole fields (first_magnet eq:dipole-plane), the entries by cubic Hermite
interpolation and bisection, the errors by doubling the step; and the rigidity estimate of each
isotope's bend (eq:kick). Compared with the game: the three reference flights and the empty
board's.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
H = 0.005
H_TABLE = 0.01


class Setup:
    def __init__(self, lvl, placement, i=0):
        shot = lvl["shots"][i]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.a_p = shot["particle"]["radius"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.Rm = lvl["physics"]["magnet_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        mags = [e for e in lvl["elements"] if e["kind"] == "magnet"] + \
               [e for e in placement if e["kind"] == "magnet"]
        self.magnets = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in mags]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def bz(self, x, y):
        """eq:dipole-plane: B_z = -mu / r^3 in the plane, summed over the magnets."""
        return sum(-mu / math.hypot(x - mx, y - my) ** 3 for mx, my, mu in self.magnets)

    def rhs(self, s):
        """eq:eom-dipole."""
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        vx, vy = px / (g * self.m), py / (g * self.m)
        b = self.bz(x, y)
        return (vx, vy, self.q * vy * b, -self.q * vx * b)

    def pabs2(self, s):
        return s[2] * s[2] + s[3] * s[3]

    def canonical(self, s):
        """eq:canonical about the first magnet: (r x p)_z + q mu / r."""
        mx, my, mu = self.magnets[0]
        rx, ry = s[0] - mx, s[1] - my
        return rx * s[3] - ry * s[2] + self.q * mu / math.hypot(rx, ry)


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
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = st.rhs(s)
    t, n = 0.0, 0
    p2_0 = st.pabs2(s)
    P0 = st.canonical(s) if st.magnets else 0.0
    rows, path, pdev, Pdev = [], [(s[0], s[1])], 0.0, 0.0
    preds = [("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0)]
    for idx, (mx, my, _) in enumerate(st.magnets):
        preds.append((f"Collided({idx})", 0,
                      lambda p, mx=mx, my=my: math.hypot(p[0] - mx, p[1] - my) <= st.Rm + st.a_p))
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                                  and st.wmin[1] <= p[1] <= st.wmax[1])))
    us = [k / 64 for k in range(1, 65)]
    while True:
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
                    "rows": rows, "path": path, "pdev": pdev, "Pdev": Pdev}
        s, f, t, n = s2, f2, t + hs, n + 1
        pdev = max(pdev, abs(math.sqrt(st.pabs2(s)) - math.sqrt(p2_0)))
        if st.magnets:
            Pdev = max(Pdev, abs(st.canonical(s) - P0))
        if record_every and n % record_every == 0:
            rows.append((t, s))
        if n % 10 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "rows": rows,
                    "path": path, "steps": n, "pdev": pdev, "Pdev": Pdev, "face": ""}


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




def bend_integral(st, t_end, h):
    """The integral of B_z along the path, |v| dt, by Simpson's rule on the march's steps up
    to the last even step and the trapezoidal rule on the rest (eq:curvature)."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    g = math.sqrt(1 + st.p0 ** 2 / (st.m * st.c) ** 2)
    v = st.p0 / (g * st.m)
    n = int(t_end / h)
    n -= n % 2
    vals = []
    for k in range(n + 1):
        vals.append(st.bz(s[0], s[1]) * v)
        if k < n:
            s = rk4_step(st, s, st.rhs(s), h)
    simpson = h / 3 * (vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2]))
    rest = t_end - n * h
    s2 = rk4_step(st, s, st.rhs(s), rest)
    return simpson + rest * 0.5 * (vals[-1] + st.bz(s2[0], s2[1]) * v)


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1])) / 15


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(len(lvl["shots"]) == 3 and len(ref) == 2 and all(e["kind"] == "magnet" for e in ref),
              "three isotopes and two placed magnets")
    st0 = Setup(lvl, ref, 0)
    doc.value("q", st0.q, 3, sci=True)
    doc.value("c", st0.c, 3)
    doc.value("Tzero", st0.T0, 3)
    doc.value("Rm", st0.Rm, 3)
    doc.value("xZero", st0.x0[0], 3)
    doc.value("yZero", st0.x0[1], 3)
    doc.value("xWorld", st0.wmax[0], 3)
    doc.value("yWorld", st0.wmax[1], 3)
    doc.value("H", H, 4)
    reg = lvl["limits"]["region"]
    for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                  ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
        doc.value(k_, v, 3)
    for k, (mx, my, mu) in enumerate(st0.magnets):
        tag = ["One", "Two"][k]
        doc.value(f"xMag{tag}", mx, 3)
        doc.value(f"yMag{tag}", my, 3)
        doc.value(f"mu{tag}", mu, 3, sci=True)
        doc.value(f"b{tag}", abs(my - st0.x0[1]), 3)

    # eq:curvature: the direction turns at the rate -q B_z / p per unit length, so the exit
    # angle is the launch angle minus (q/p) times the integral of B_z along the actual path.
    rows = []
    tags = ["One", "Two", "Four"]
    worst = 0.0
    worst_bend = 0.0
    for i, tag in enumerate(tags):
        st = Setup(lvl, ref, i)
        p = st.p0
        run, half = fly(st, H), fly(st, 2 * H)
        doc.check(run["outcome"] == "Arrived", tag + ": the reference flight arrives")
        err = err_of(run, half)
        worst = max(worst, err)
        ang = math.degrees(math.atan2(run["p"][1], run["p"][0]))
        I = bend_integral(st, run["t"], H)
        pred = math.degrees(math.atan2(st.d[1], st.d[0]) - st.q * I / p)
        worst_bend = max(worst_bend, abs(pred - ang))
        rows.append([f"{st.m:g}", p, I, pred, ang, run["x"][1],
                     f"[{st.dmin[1]:g}, {st.dmax[1]:g}]"])
        doc.value(f"m{tag}", st.m, 3)
        doc.value(f"p{tag}", p, 8)
        doc.value(f"bendInt{tag}", I, 6)
        doc.value(f"ang{tag}", ang, 5)
        doc.value(f"t{tag}", run["t"], 10)
        doc.value(f"x{tag}", run["x"][0], 10)
        doc.value(f"y{tag}", run["x"][1], 10)
        doc.text(f"face{tag}", run["face"])
        doc.value(f"pDev{tag}", run["pdev"], 2, sci=True)
        doc.data("path" + tag.lower(), run["path"])
        for k_, v in (("yDetLo", st.dmin[1]), ("yDetHi", st.dmax[1])):
            doc.value(k_ + tag, v, 3)
        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * err, gap)
        title = f"mass {st.m:g}: "
        doc.compare(f"Out{tag}", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(f"T{tag}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
        doc.compare(f"X{tag}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
        doc.compare(f"Y{tag}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
        doc.compare(f"Py{tag}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
        bare = Setup(lvl, [], i)
        rb = fly(bare, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{tag}", f"empty board, mass {st.m:g}: outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BareT{tag}", f"empty board, mass {st.m:g}: " + r"$t_\mathrm{e}$", rb["t"], gb["end"]["t"], TOL)
        doc.text(f"bareOut{tag}", rb["outcome"])
    doc.table("bends", [r"$m$", r"$p$", r"$\int B_z\,\dd s$", r"by \eqref{eq:curvature} (deg)",
                        r"exit angle (deg)", r"$y_\mathrm{e}$", "detector's $y$"], rows,
              align="rrrrrrl", digits=[None, 6, 6, 6, 6, 6, None])
    doc.value("errMax", worst, 2, sci=True)
    doc.check(worst < 1e-9, "the step-halving estimates are below 1e-9")
    doc.value("bendAgree", worst_bend, 2, sci=True)
    doc.check(worst_bend < 1e-6, "the exit angles agree with the curvature integral")
