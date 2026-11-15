"""Companion of docs/levels/beta_spectrometer/doc.tex (Level 32, Beta-ray spectrometer).

Every number of the document, from the level file. Each electron's flight about the one
attracting charge is Level 1's closed form (first_bend eq:orbit, eq:general): the orbit
u(theta) = A cos(Gamma theta) + B from an apsis, the entry into the detector by stepping theta
and bisecting, the time by Simpson's rule on eq:dtdth (its error by halving the step), the
momentum by first_bend eq:momentum. Independently, the classical Runge-Kutta method at a
fixed step on the equations of motion (Foundations eq:rk4) checks every end state. Compared
with the game: both shots of the reference and of the empty board.
"""
import math


TOL = 1e-9          # ten times the larger of our estimates (~1e-11) and the game's (1e-10)
H = 0.002           # the march's step
N_SIMPSON = 4000    # Simpson's panels on the time integral (an even number)


class Shot:
    def __init__(self, lvl, placement, i):
        shot = lvl["shots"][i]
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
        els = [e for e in lvl["elements"] if e["kind"] == "charge"] + list(placement)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in els]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = math.sqrt(self.T0 ** 2 / self.c ** 2 + 2 * self.m * self.T0)

    def rhs(self, s):
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (self.m * self.m * self.c * self.c))
        fx = fy = 0.0
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            fx += self.q * Q * dx / r3
            fy += self.q * Q * dy / r3
        return (px / (g * self.m), py / (g * self.m), fx, fy)

    def energy(self, s):
        x, y, px, py = s
        V = sum(self.q * Q / math.hypot(x - cx, y - cy) for cx, cy, Q in self.charges)
        return math.sqrt(self.m ** 2 * self.c ** 4 + (px * px + py * py) * self.c ** 2) + V


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


# ---- The closed form (first_bend) -----------------------------------------------------------

class Orbit:
    """first_bend eq:general about one charge, with theta the angle swept since the launch."""

    def __init__(self, sh):
        (cx, cy, Q), = sh.charges
        self.cx, self.cy = cx, cy
        self.kappa = sh.q * Q
        rx, ry = sh.x0[0] - cx, sh.x0[1] - cy
        self.r0 = math.hypot(rx, ry)
        self.phi0 = math.atan2(ry, rx)
        px, py = sh.p0 * sh.d[0], sh.p0 * sh.d[1]
        Lz = rx * py - ry * px
        self.sigma = 1.0 if Lz > 0 else -1.0
        self.L = abs(Lz)
        self.m, self.c = sh.m, sh.c
        self.W = sh.gamma0 * sh.m * sh.c ** 2 + self.kappa / self.r0      # first_bend eq:energy
        self.G2 = 1 - self.kappa ** 2 / (sh.c ** 2 * self.L ** 2)          # eq:orbit
        # Gamma^2 < 0: first_bend's hyperbolic case, cos -> cosh (eq:spiral).
        self.hyp = self.G2 < 0
        self.G = math.sqrt(abs(self.G2))
        self.B = -self.kappa * self.W / (sh.c ** 2 * self.L ** 2 * self.G2)
        self.u0 = 1 / self.r0
        pr0 = (rx * px + ry * py) / self.r0
        self.du0 = -pr0 / self.L                                            # eq:radial
        self.A = self.u0 - self.B

    def u(self, th):
        g = self.G * th
        if self.hyp:
            return self.B + self.A * math.cosh(g) + self.du0 / self.G * math.sinh(g)
        return self.B + self.A * math.cos(g) + self.du0 / self.G * math.sin(g)

    def du(self, th):
        g = self.G * th
        if self.hyp:
            return self.A * self.G * math.sinh(g) + self.du0 * math.cosh(g)
        return -self.A * self.G * math.sin(g) + self.du0 * math.cos(g)

    def pos(self, th):
        phi = self.phi0 + self.sigma * th
        r = 1 / self.u(th)
        return (self.cx + r * math.cos(phi), self.cy + r * math.sin(phi))

    def dtdth(self, th):
        u = self.u(th)
        return (self.W - self.kappa * u) / (self.c ** 2 * self.L * u * u)   # first_bend eq:dtdth

    def time(self, th, n):
        h = th / n
        s = self.dtdth(0) + self.dtdth(th)
        s += 4 * sum(self.dtdth((2 * k - 1) * h) for k in range(1, n // 2 + 1))
        s += 2 * sum(self.dtdth(2 * k * h) for k in range(1, n // 2))
        return s * h / 3

    def momentum(self, th):
        phi = self.phi0 + self.sigma * th
        pr, pphi = -self.L * self.du(th), self.sigma * self.L * self.u(th)  # eq:momentum
        return (pr * math.cos(phi) - pphi * math.sin(phi), pr * math.sin(phi) + pphi * math.cos(phi))


def closed_form_entry(sh, orb):
    """The first theta at which the orbit is inside the detector, by steps of 1e-3 and
    bisection; the events before it (the charge, the world) are checked on the same steps."""
    th, step = 0.0, 1e-3
    while th < 50:
        nxt = th + step
        p = orb.pos(nxt)
        r = 1 / orb.u(nxt)
        if r <= sh.R + sh.a or not (sh.wmin[0] <= p[0] <= sh.wmax[0] and sh.wmin[1] <= p[1] <= sh.wmax[1]):
            return None
        if box_distance(p, sh.dmin, sh.dmax) <= 0:
            lo, hi = th, nxt
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if box_distance(orb.pos(mid), sh.dmin, sh.dmax) <= 0:
                    hi = mid
                else:
                    lo = mid
            return hi
        th = nxt
    return None


# ---- The march (independent check) ------------------------------------------------------------

def rk4_step(sh, s, f, h):
    k1 = f
    k2 = sh.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)))
    k3 = sh.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)))
    k4 = sh.rhs(tuple(a + h * b for a, b in zip(s, k3)))
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def fly(sh, h):
    s = (sh.x0[0], sh.x0[1], sh.p0 * sh.d[0], sh.p0 * sh.d[1])
    f = sh.rhs(s)
    t, n = 0.0, 0
    w0 = sh.energy(s)
    wdev = 0.0
    path = [(s[0], s[1])]
    preds = [(f"Collided({i})", 0, lambda p, cx=cx, cy=cy: math.hypot(p[0] - cx, p[1] - cy) <= sh.R + sh.a)
             for i, (cx, cy, _) in enumerate(sh.charges)]
    preds.append(("LeftBounds", 1, lambda p: not (sh.wmin[0] <= p[0] <= sh.wmax[0] and sh.wmin[1] <= p[1] <= sh.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, sh.dmin, sh.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
    while True:
        hs = min(h, sh.tmax - t)
        s2 = rk4_step(sh, s, f, hs)
        f2 = sh.rhs(s2)
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
            se = rk4_step(sh, s, f, te - t)
            pe = hermite(s, f, s2, f2, hs, hi)
            face = face_of(pe, sh.dmin, sh.dmax) if outcome == "Arrived" else (
                face_of(pe, sh.wmin, sh.wmax) if outcome == "LeftBounds" else "")
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:], "face": face, "steps": n + 1,
                    "wdev": wdev, "path": path}
        s, f, t, n = s2, f2, t + hs, n + 1
        wdev = max(wdev, abs(sh.energy(s) - w0))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "face": "",
                    "steps": n, "wdev": wdev, "path": path}


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(len(lvl["shots"]) == 2 and len(ref) == 1 and not lvl["elements"],
              "two electrons and the reference's one charge")
    tags = ["One", "Two"]
    for i, tag in enumerate(tags):
        sh = Shot(lvl, ref, i)
        orb = Orbit(sh)
        if i == 0:
            (cx, cy, Q), = sh.charges
            doc.value("q", sh.q, 3, sci=True)
            doc.value("m", sh.m, 3)
            doc.value("c", sh.c, 3)
            doc.value("R", sh.R, 3)
            doc.value("Q", Q, 3, sci=True)
            doc.value("xC", cx, 3)
            doc.value("yC", cy, 3)
            doc.value("kappa", orb.kappa, 6)
            doc.value("xZero", sh.x0[0], 3)
            doc.value("yZero", sh.x0[1], 3)
            doc.value("rZero", orb.r0, 6)
            doc.value("xWorld", sh.wmax[0], 3)
            doc.value("yWorld", sh.wmax[1], 3)
            doc.value("H", H, 4)
            doc.value("nSimpson", N_SIMPSON)
            reg = lvl["limits"]["region"]
            for k_, v in (("regX0", reg["min"][0]), ("regY0", reg["min"][1]),
                          ("regX1", reg["max"][0]), ("regY1", reg["max"][1])):
                doc.value(k_, v, 3)
        v0 = sh.p0 / (sh.gamma0 * sh.m)
        doc.value(f"T{tag}", sh.T0, 10)
        doc.value(f"gamma{tag}", sh.gamma0, 10)
        doc.value(f"p{tag}", sh.p0, 12)
        doc.value(f"beta{tag}", v0 / sh.c, 6)
        for k_, v in (("xDetLo", sh.dmin[0]), ("yDetLo", sh.dmin[1]), ("xDetHi", sh.dmax[0]),
                      ("yDetHi", sh.dmax[1])):
            doc.value(k_ + tag, v, 3)
        # The circular-orbit condition (eq:circle): |kappa| = gamma m v^2 r0 = p v r0.
        doc.value(f"pvr{tag}", sh.p0 * v0 * orb.r0, 10)
        doc.value(f"newton{tag}", 2 * sh.T0 * orb.r0, 6)
        doc.value(f"W{tag}", orb.W, 12)
        doc.value(f"L{tag}", orb.L, 12)
        doc.value(f"Gamma{tag}", orb.G, 10)
        doc.value(f"GammaSq{tag}", orb.G2, 10)
        doc.text(f"kind{tag}", "hyperbolic" if orb.hyp else "trigonometric")
        doc.value(f"B{tag}", orb.B, 12)
        doc.value(f"A{tag}", orb.A, 6, sci=abs(orb.A) < 1e-3)
        doc.value(f"du{tag}", orb.du0, 2, sci=True)
        th = closed_form_entry(sh, orb)
        doc.check(th is not None, tag + ": the closed form reaches the detector")
        pe = orb.pos(th)
        t1, t2 = orb.time(th, N_SIMPSON), orb.time(th, N_SIMPSON // 2)
        terr = abs(t1 - t2) / 15
        pm = orb.momentum(th)
        doc.value(f"thE{tag}", th, 10)
        doc.value(f"thEdeg{tag}", math.degrees(th), 6)
        doc.value(f"tE{tag}", t1, 12)
        doc.value(f"tErr{tag}", terr, 2, sci=True)
        doc.value(f"xE{tag}", pe[0], 12)
        doc.value(f"yE{tag}", pe[1], 12)
        doc.value(f"pxE{tag}", pm[0], 10)
        doc.value(f"pyE{tag}", pm[1], 10)
        face = face_of(pe, sh.dmin, sh.dmax)
        doc.text(f"face{tag}", face)
        # The march, an independent check.
        run, half = fly(sh, H), fly(sh, H / 2)
        doc.check(run["outcome"] == "Arrived", tag + ": the march arrives")
        merr = max(abs(run["t"] - half["t"]), abs(run["x"][0] - half["x"][0]), abs(run["x"][1] - half["x"][1])) / 15
        dmarch = max(abs(run["t"] - t1), abs(run["x"][0] - pe[0]), abs(run["x"][1] - pe[1]),
                     abs(run["p"][0] - pm[0]), abs(run["p"][1] - pm[1]))
        doc.value(f"mErr{tag}", merr, 2, sci=True)
        doc.value(f"march{tag}", dmarch, 2, sci=True)
        doc.value(f"wDev{tag}", run["wdev"], 2, sci=True)
        doc.check(dmarch < 1e-9, tag + ": the march agrees with the closed form")
        doc.data("path" + tag.lower(), run["path"])

        g = doc.flight(0, flight=i)
        gp = doc.flight(0, flight=i, run="preview")
        gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
                  abs(gp["end"]["x"][1] - g["end"]["x"][1]))
        tol = max(TOL, 10 * max(terr, merr), gap)
        doc.game_value(f"gap{tag}", gap, 2, sci=True)
        title = ("electron 1: " if i == 0 else "electron 2: ")
        doc.compare(f"Out{tag}", title + "outcome", "Arrived", g["outcome"])
        doc.compare(f"Tc{tag}", title + r"$t_\mathrm{e}$", t1, g["end"]["t"], tol)
        doc.compare(f"Xc{tag}", title + r"$x_\mathrm{e}$", pe[0], g["end"]["x"][0], tol)
        doc.compare(f"Yc{tag}", title + r"$y_\mathrm{e}$", pe[1], g["end"]["x"][1], tol)
        doc.compare(f"Pxc{tag}", title + r"$p_{x,\mathrm{e}}$", pm[0], g["end"]["p"][0], tol)
        doc.compare(f"Pyc{tag}", title + r"$p_{y,\mathrm{e}}$", pm[1], g["end"]["p"][1], tol)

        bare = Shot(lvl, [], i)
        rb = fly(bare, H)
        gb = doc.flight(1, flight=i)
        doc.compare(f"BareOut{tag}", "empty board, " + title + "outcome", rb["outcome"], gb["outcome"])
        doc.compare(f"BareT{tag}", "empty board, " + title + r"$t_\mathrm{e}$", rb["t"], gb["end"]["t"], TOL)
        if i == 0:
            doc.text("bareFace", rb["face"])
            doc.value("bareX", rb["x"][0], 4)
            doc.value("bareY", rb["x"][1], 4)
