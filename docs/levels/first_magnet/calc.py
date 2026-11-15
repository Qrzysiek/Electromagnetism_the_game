"""Companion of docs/levels/first_magnet/doc.tex (Level 29, First magnet).

Every number of the document, from the level file, by the document's two methods:
- the flight by the classical Runge-Kutta method at a fixed step (Foundations eq:rk4) on
  eq:eom-dipole, the entry into the detector by cubic Hermite interpolation and bisection,
  the error by halving the step;
- the same flight reduced to quadratures by its two conserved quantities, |p| (the magnetic
  force does no work) and the canonical angular momentum about the magnet (eq:canonical):
  the closest approach from a quadratic (eq:rmin), the angle and the time as integrals over
  r (eq:phi-of-r, eq:t-of-r) by mpmath's quadrature, the entry by bisection on r.
Compared with the game: the reference flight and the empty board.
"""
import math

import mpmath


TOL = 1e-9          # ten times the larger of our estimate (~1e-12) and the game's (1e-10)
TOL_RAD = 1e-4
H = 0.005
H_TABLE = 0.01


class Setup:
    def __init__(self, lvl, placement):
        shot = lvl["shots"][0]
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


def quadrature_flight(st, doc):
    """The reference flight from its conserved quantities (section 4.3), in mpmath."""
    mpmath.mp.dps = 30
    mx, my, mu = (mpmath.mpf(v) for v in st.magnets[0])
    q, m, c = mpmath.mpf(st.q), mpmath.mpf(st.m), mpmath.mpf(st.c)
    x0, y0 = mpmath.mpf(st.x0[0]), mpmath.mpf(st.x0[1])
    p = mpmath.sqrt(mpmath.mpf(st.T0) ** 2 / c ** 2 + 2 * m * mpmath.mpf(st.T0))
    px0, py0 = p * mpmath.mpf(st.d[0]), p * mpmath.mpf(st.d[1])
    gamma = mpmath.sqrt(1 + p ** 2 / (m * c) ** 2)
    rx0, ry0 = x0 - mx, y0 - my
    r0 = mpmath.hypot(rx0, ry0)
    phi0 = mpmath.atan2(ry0, rx0)
    qmu = q * mu
    P = rx0 * py0 - ry0 * px0 + qmu / r0              # eq:canonical
    L = lambda r: P - qmu / r                          # (r x p)_z on the orbit
    pr0 = (rx0 * px0 + ry0 * py0) / r0
    # eq:rmin: |L(r)| = p r, with L < 0 here: p r^2 + P r - q mu = 0 ... solved for the sign.
    sgn = -1 if L(r0) < 0 else 1
    # L(r) = sgn * p * r  <=>  p r^2 - sgn P r + sgn q mu = 0
    a_, b_, c_ = p, -sgn * P, sgn * qmu
    disc = b_ * b_ - 4 * a_ * c_
    roots = [(-b_ + mpmath.sqrt(disc)) / (2 * a_), (-b_ - mpmath.sqrt(disc)) / (2 * a_)]
    rmin = max(rt for rt in roots if rt < r0)
    pr2 = lambda r: p ** 2 - (L(r) / r) ** 2

    def integral(f, R):
        # r = rmin + u^2 removes the turning point's inverse square root.
        return mpmath.quad(lambda u: f(rmin + u * u) * 2 * u, [0, mpmath.sqrt(R - rmin)])

    dphi = lambda r: L(r) / (r * r * mpmath.sqrt(pr2(r)))   # eq:phi-of-r
    dt = lambda r: gamma * m / mpmath.sqrt(pr2(r))           # eq:t-of-r
    I0, J0 = integral(dphi, r0), integral(dt, r0)
    inward = pr0 < 0

    def point(r):
        """Position, time and momentum on the outgoing branch at the radius r."""
        phi = phi0 + (I0 if inward else -I0) + integral(dphi, r)
        t = (J0 if inward else -J0) + integral(dt, r)
        x, y = mx + r * mpmath.cos(phi), my + r * mpmath.sin(phi)
        prr, pphi = mpmath.sqrt(pr2(r)), L(r) / r
        pxv = prr * mpmath.cos(phi) - pphi * mpmath.sin(phi)
        pyv = prr * mpmath.sin(phi) + pphi * mpmath.cos(phi)
        return x, y, t, pxv, pyv

    lo_x, lo_y = mpmath.mpf(st.dmin[0]), mpmath.mpf(st.dmin[1])
    hi_x, hi_y = mpmath.mpf(st.dmax[0]), mpmath.mpf(st.dmax[1])

    def inside(r):
        x, y, *_ = point(r)
        return lo_x <= x <= hi_x and lo_y <= y <= hi_y

    # March out in r until the point is inside the detector, then bisect on r.
    r_a, r_b = rmin, None
    step = mpmath.mpf("0.25")
    r = rmin + step
    while r < 60:
        if inside(r):
            r_b = r
            break
        r_a, r = r, r + step
    doc.check(r_b is not None, "the quadrature flight reaches the detector")
    for _ in range(100):
        mid = (r_a + r_b) / 2
        if inside(mid):
            r_b = mid
        else:
            r_a = mid
    x, y, t, pxv, pyv = point(r_b)
    return {"P": P, "rmin": rmin, "r0": r0, "pr0": pr0, "x": (x, y), "t": t, "p": (pxv, pyv),
            "rE": r_b}


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    st = Setup(lvl, ref)
    doc.check(len(st.magnets) == 1 and not lvl["elements"], "one magnet, placed by the reference")
    mx, my, mu = st.magnets[0]
    for k, v in (("q", st.q), ("m", st.m), ("c", st.c), ("Tzero", st.T0), ("Rm", st.Rm)):
        doc.value(k, v, 3, sci=(k == "q"))
    doc.value("xZero", st.x0[0], 3)
    doc.value("yZero", st.x0[1], 3)
    doc.value("xDetLo", st.dmin[0], 3)
    doc.value("yDetLo", st.dmin[1], 3)
    doc.value("xDetHi", st.dmax[0], 3)
    doc.value("yDetHi", st.dmax[1], 3)
    doc.value("xWorld", st.wmax[0], 3)
    doc.value("yWorld", st.wmax[1], 3)
    doc.value("xMag", mx, 3)
    doc.value("yMag", my, 3)
    doc.value("muMag", mu, 3, sci=True)
    doc.value("qmu", st.q * mu, 4)
    reg = lvl["limits"]["region"]
    doc.value("regX0", reg["min"][0], 3)
    doc.value("regY0", reg["min"][1], 3)
    doc.value("regX1", reg["max"][0], 3)
    doc.value("regY1", reg["max"][1], 3)
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p0, 12)

    # ---- The field and the circle at the launch (eq:dipole-plane, eq:larmor-radius) --------
    r0 = math.hypot(st.x0[0] - mx, st.x0[1] - my)
    b0 = st.bz(*st.x0)
    doc.value("rZero", r0, 6)
    doc.value("Bzero", b0, 6)
    doc.value("rLzero", st.p0 / abs(st.q * b0), 4)

    # ---- The reference flight by RK4 ---------------------------------------------------------
    run, half = fly(st, H), fly(st, H / 2)
    doc.check(run["outcome"] == "Arrived", "the reference flight arrives")
    err = err_of(run, half)
    te, xe, pe = run["t"], run["x"], run["p"]
    doc.value("H", H, 4)
    doc.value("steps", run["steps"])
    doc.value("tE", te, 12)
    doc.value("xE", xe[0], 12)
    doc.value("yE", xe[1], 12)
    doc.text("faceE", run["face"])
    doc.value("pxE", pe[0], 10)
    doc.value("pyE", pe[1], 10)
    doc.value("angE", math.degrees(math.atan2(pe[1], pe[0])), 5)
    doc.value("err", err, 2, sci=True)
    doc.value("pDev", run["pdev"], 2, sci=True)
    doc.value("PDev", run["Pdev"], 2, sci=True)
    doc.check(err < 1e-10, "the step-halving estimate is below 1e-10")
    tab = fly(st, H_TABLE, record_every=int(round(5 / H_TABLE)))
    doc.value("HTable", H_TABLE, 3)
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$P(t) - P$"],
              [[f"{t:.0f}", s[0], s[1], s[2], s[3], st.canonical(s) - st.canonical(
                  (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1]))] for t, s in tab["rows"]],
              align="rrrrrr", digits=[None, 8, 8, 8, 8, 2], sci=[False] * 5 + [True])
    E_rad = radiated(st, H, te)
    E_rad2 = radiated(st, 2 * H, te)
    doc.value("Erad", E_rad, 6, sci=True)
    doc.value("EradFrac", E_rad / st.T0, 3, sci=True)
    doc.value("EradErr", abs(E_rad - E_rad2) / abs(E_rad2) / 15, 2, sci=True)
    doc.data("path", run["path"])

    # ---- The same flight by quadratures ------------------------------------------------------
    qf = quadrature_flight(st, doc)
    doc.value("Pcan", float(qf["P"]), 12)
    doc.value("Lzero", float(qf["P"] - st.q * mu / qf["r0"]), 12)
    doc.value("prZero", float(qf["pr0"]), 10)
    doc.value("rMin", float(qf["rmin"]), 12)
    doc.value("rEntry", float(qf["rE"]), 10)
    dq = max(abs(float(qf["t"]) - te), abs(float(qf["x"][0]) - xe[0]), abs(float(qf["x"][1]) - xe[1]),
             abs(float(qf["p"][0]) - pe[0]), abs(float(qf["p"][1]) - pe[1]))
    doc.value("tQ", float(qf["t"]), 12)
    doc.value("yQ", float(qf["x"][1]), 12)
    doc.value("xQ", float(qf["x"][0]), 12)
    doc.value("pxQ", float(qf["p"][0]), 10)
    doc.value("pyQ", float(qf["p"][1]), 10)
    doc.value("quadDiff", dq, 2, sci=True)
    doc.check(dq < 1e-9, "the quadratures agree with the march")
    # The closest approach in the march, against eq:rmin.
    rmin_march = min(math.hypot(x - mx, y - my) for x, y in run["path"])
    doc.value("rMinMarch", rmin_march, 6)

    # ---- The empty board ----------------------------------------------------------------------
    bare = Setup(lvl, [])
    rb, rb2 = fly(bare, H), fly(bare, H / 2)
    doc.text("bareOut", rb["outcome"])
    doc.text("bareFace", rb["face"])
    doc.value("bareT", rb["t"], 10)
    doc.value("bareX", rb["x"][0], 10)
    doc.value("bareY", rb["x"][1], 10)

    # ---- The comparison with the game ----------------------------------------------------------
    def game_gap(i):
        a, b = doc.flight(i, run="preview")["end"], doc.flight(i)["end"]
        return max(abs(a["t"] - b["t"]), *(abs(u - v) for u, v in zip(a["x"], b["x"])),
                   *(abs(u - v) for u, v in zip(a["p"], b["p"])))
    tol_ref = max(TOL, 10 * err, game_gap(0))
    tol_bare = max(TOL, 10 * err_of(rb, rb2), game_gap(1))
    doc.game_value("gapRef", game_gap(0), 2, sci=True)
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
