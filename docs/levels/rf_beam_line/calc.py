"""Companion of docs/levels/rf_beam_line/doc.tex (Level 48, RF beam line).

Every number of the document, from the level file. The charges' Coulomb fields, the uniform
stray field of each disturbance, the antennas' exact retarded fields (rf_kick eq:fields) at the
level's frequency; each ion's flight from its launch time, the fields at t0 + t (rf_separator
eq:offset); the relativistic equations of motion (Foundations eq:eom) marched by the classical
Runge-Kutta method at a fixed step (Foundations eq:rk4) with Lienard's radiated energy as one more
state; the gate (two_stages eq:gate) and its direction cone (injection eq:cone) tested at every
entry; the entries by cubic Hermite interpolation and bisection; the error against the step 2h.
The empty board in closed form: a straight line, and the hyperbolic motion in the uniform field
(eq:hyperbolic). Half a period late is the antennas reversed (eq:late): checked by marching.
Compared with the game: the reference, the empty board and the reference's charges alone, under
both disturbances, every ion.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-6      # the radiated energy, as in rf_kick: the game's Gauss-Legendre per step, ~1e-7; ten times
H = 0.004
SHOTS = ("One", "Two")
DISTS = ("Off", "On")


class Antenna:
    """rf_kick eq:fields: the retarded fields of p(t) = p0 cos(omega t)."""

    def __init__(self, e, s, omega, c, sign=1.0):
        a = math.radians(e["angle_deg"])
        self.pos = (e["node"][0] / s, e["node"][1] / s, 0.0)
        self.amp = (sign * math.cos(a) * e["value"], sign * math.sin(a) * e["value"], 0.0)
        self.omega, self.c = omega, c

    def moment(self, t):
        co, si, w = math.cos(self.omega * t), math.sin(self.omega * t), self.omega
        return ([a * co for a in self.amp], [-w * a * si for a in self.amp],
                [-w * w * a * co for a in self.amp])

    def fields(self, x, t):
        d = [x[k] - self.pos[k] for k in range(3)]
        r = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        n = [v / r for v in d]
        c = self.c
        p, pd, pdd = self.moment(t - r / c)
        dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        cross = lambda a, b: (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        npd, npp, nq = dot(n, p), dot(n, pd), dot(n, pdd)
        e = [(3 * n[k] * npd - p[k]) / r ** 3 + (3 * n[k] * npp - pd[k]) / (c * r * r)
             + (n[k] * nq - pdd[k]) / (c * c * r) for k in range(3)]
        b1, b2 = cross(pd, n), cross(pdd, n)
        b = [(b1[k] / (r * r) + b2[k] / (c * r)) / (c * c) for k in range(3)]
        return e, b


class Setup:
    def __init__(self, lvl, placement, shot_i, dist, t0=None, antenna_sign=1.0):
        shot = lvl["shots"][shot_i]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.t0 = shot["launch"].get("time", 0.0) if t0 is None else t0
        self.c = ph["c"]
        self.omega = ph["rf_omega"]
        self.Ra = ph["antenna_radius"]
        self.Rq = ph["charge_radius"]
        self.tmax = ph["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        (g,) = lvl["gates"]
        self.gmin = [v / s for v in g["min"][:2]]
        self.gmax = [v / s for v in g["max"][:2]]
        axis, half = g["acceptance"]["direction"]
        self.gaxis = (math.cos(math.radians(axis)), math.sin(math.radians(axis)))
        self.ghalf = math.radians(half)
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        els = list(lvl["elements"]) + list(placement)
        assert all(e["kind"] in ("charge", "antenna") for e in els)
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in els if e["kind"] == "charge"]
        self.antennas = [Antenna(e, s, e.get("omega") or self.omega, self.c, antenna_sign)
                         for e in els if e["kind"] == "antenna"]
        self.E = lvl["disturbances"][dist]["e"]
        self.Bz = lvl["disturbances"][dist]["bz"]
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)
        self.v0 = self.p0 / (self.gamma0 * self.m)

    def field(self, x, y, t):
        """The stray field, the charges' Coulomb fields, the antennas at t0 + t (rf_separator
        eq:offset)."""
        ex, ey, bz = self.E[0], self.E[1], self.Bz
        for cx, cy, Q in self.charges:
            dx, dy = x - cx, y - cy
            r3 = (dx * dx + dy * dy) ** 1.5
            ex += Q * dx / r3
            ey += Q * dy / r3
        for a in self.antennas:
            e, b = a.fields((x, y, 0.0), self.t0 + t)
            ex, ey, bz = ex + e[0], ey + e[1], bz + b[2]
        return ex, ey, bz

    def rhs(self, s, t):
        x, y, px, py = s[:4]
        m, c, q = self.m, self.c, self.q
        g = math.sqrt(1 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        ex, ey, bz = self.field(x, y, t)
        fx, fy = q * (ex + vy * bz), q * (ey - vx * bz)
        vf = vx * fx + vy * fy
        ax, ay = (fx - vx * vf / c ** 2) / (g * m), (fy - vy * vf / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        power = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        return (vx, vy, fx, fy, power)


def box_distance(p, lo, hi):
    return max(lo[0] - p[0], p[0] - hi[0], lo[1] - p[1], p[1] - hi[1])


def face_of(pe, lo, hi):
    for name, val, k in (("left", lo[0], 0), ("right", hi[0], 0), ("bottom", lo[1], 1), ("top", hi[1], 1)):
        if abs(pe[k] - val) < 1e-9:
            return name
    return ""


def rk4_step(st, s, f, t, h):
    k1 = f
    k2 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k1)), t + h / 2)
    k3 = st.rhs(tuple(a + h / 2 * b for a, b in zip(s, k2)), t + h / 2)
    k4 = st.rhs(tuple(a + h * b for a, b in zip(s, k3)), t + h)
    return tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4))


def hermite(s, f, s2, f2, h, u):
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    return (h00 * s[0] + h10 * h * f[0] + h01 * s2[0] + h11 * h * f2[0],
            h00 * s[1] + h10 * h * f[1] + h01 * s2[1] + h11 * h * f2[1])


def locate(s, f, s2, f2, hs, pred, u_hit):
    """The first root of pred on the interpolant, between the sample before u_hit and u_hit."""
    lo, hi = u_hit - 1 / 64, u_hit
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pred(hermite(s, f, s2, f2, hs, mid)):
            hi = mid
        else:
            lo = mid
    return hi


def fly(st, h):
    """The march to the first event that ends the flight (Foundations sec:board). The gate
    (two_stages eq:gate) is tested at each entry against its cone (injection eq:cone) and does not
    end the flight; a detector entered before the gate is passed is a skipped gate."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    path = [(s[0], s[1])]
    gate = {"passed": False, "inside": False, "entries": []}
    in_gate = lambda p: box_distance(p, st.gmin, st.gmax) <= 0
    preds = [(f"Collided({i})", 0, lambda p, c=c: math.hypot(p[0] - c[0], p[1] - c[1]) <= st.Rq)
             for i, c in enumerate(st.charges)]
    preds += [(f"Collided({len(st.charges) + i})", 0,
               lambda p, a=a: math.hypot(p[0] - a.pos[0], p[1] - a.pos[1]) <= st.Ra)
              for i, a in enumerate(st.antennas)]
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 64 for k in range(1, 65)]
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, t, hs)
        f2 = st.rhs(s2, t + hs)
        pts = [hermite(s, f, s2, f2, hs, u) for u in us]
        if not gate["passed"]:
            if gate["inside"]:
                if any(not in_gate(p) for p in pts):
                    gate["inside"] = False      # left the gate: the next entry is tested again
            else:
                k_in = next((k for k, p in enumerate(pts) if in_gate(p)), None)
                if k_in is not None:
                    u = locate(s, f, s2, f2, hs, in_gate, us[k_in])
                    se = rk4_step(st, s, f, t, u * hs)
                    pe = hermite(s, f, s2, f2, hs, u)
                    ax, ay = st.gaxis
                    beta = math.atan2(abs(se[2] * ay - se[3] * ax), se[2] * ax + se[3] * ay)
                    theta = math.atan2(se[3] * ax - se[2] * ay, se[2] * ax + se[3] * ay)
                    ok = beta <= st.ghalf
                    gate["entries"].append({"t": t + u * hs, "x": pe, "beta": beta, "theta": theta, "ok": ok,
                                            "face": face_of(pe, st.gmin, st.gmax)})
                    gate["passed"], gate["inside"] = ok, not ok
        events = []
        for name, rank, pred in preds:
            for k, p in enumerate(pts):
                if pred(p):
                    events.append((us[k], rank, name, pred))
                    break
        if events:
            u_hit, _r, outcome, pred = min(events, key=lambda e: (e[0], e[1]))
            u = locate(s, f, s2, f2, hs, pred, u_hit)
            te = t + u * hs
            se = rk4_step(st, s, f, t, te - t)
            pe = hermite(s, f, s2, f2, hs, u)
            if outcome == "Arrived" and not gate["passed"]:
                outcome = "SkippedGate(0)"
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else ""
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "path": path,
                    "rad": se[4], "gate": gate}
        s, f, t, n = s2, f2, t + hs, n + 1
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "path": path,
                    "rad": s[4], "gate": gate}


def err_of(a, b):
    """Foundations sec:numerics: the h run's error from the 2h run, |y_2h - y_h| / 15."""
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def end_gap(a, b):
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1]))


def bare(st):
    """The empty board in closed form. Off: the straight line y = y0 to the right edge. On:
    eq:hyperbolic, p_x = p_i, p_y = F t, eps = sqrt(eps0^2 + c^2 F^2 t^2), y - y0 = (eps - eps0)/F,
    x - x0 = (p_i c/|F|) asinh(c |F| t/eps0); the first of the bottom edge and the right edge."""
    F = st.q * st.E[1]
    c, pi = st.c, st.p0
    if F == 0.0:
        return "LeftBounds", (st.wmax[0] - st.x0[0]) / st.v0, (st.wmax[0], st.x0[1]), (pi, 0.0)
    eps0 = st.gamma0 * st.m * c * c
    a = abs(F)
    y_edge = st.wmin[1] if F < 0 else st.wmax[1]
    eps_e = eps0 + F * (y_edge - st.x0[1])
    t_y = math.sqrt(eps_e ** 2 - eps0 ** 2) / (c * a)
    x_y = st.x0[0] + pi * c / a * math.asinh(c * a * t_y / eps0)
    if x_y <= st.wmax[0]:
        return "LeftBounds", t_y, (x_y, y_edge), (pi, F * t_y)
    t_x = eps0 / (c * a) * math.sinh(a * (st.wmax[0] - st.x0[0]) / (pi * c))
    eps = math.sqrt(eps0 ** 2 + (c * F * t_x) ** 2)
    return "LeftBounds", t_x, (st.wmax[0], st.x0[1] + (eps - eps0) / F), (pi, F * t_x)


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    lim = lvl["limits"]
    shots, dists = lvl["shots"], lvl["disturbances"]
    doc.check(lvl["physics"]["c"] is not None and len(shots) == 2 and len(dists) == 2
              and len(lvl["gates"]) == 1 and not lvl["elements"] and not lvl["physics"].get("radiation_reaction")
              and all(e["kind"] in ("charge", "antenna") and "omega" not in e for e in ref),
              "relativistic, two shots, two disturbances, one gate, charges and antennas at the level's frequency")
    doc.check(shots[0]["particle"] == shots[1]["particle"]
              and all(shots[0]["launch"][k] == shots[1]["launch"][k] for k in ("node", "direction", "kinetic_energy")),
              "the two shots differ only in launch time and detector")
    doc.check(all(d["e"][0] == 0.0 and d["bz"] == 0.0 for d in dists) and dists[0]["e"][1] == 0.0,
              "the disturbances: no field, and a uniform field along y")
    a = Setup(lvl, ref, 0, 0)
    b = Setup(lvl, ref, 1, 0)
    w = a.omega
    period = 2 * math.pi / w
    doc.check(abs(b.t0 - period / 2) < 1e-12, "the second ion leaves half a period after the first")
    charges = [e for e in ref if e["kind"] == "charge"]
    antennas = [e for e in ref if e["kind"] == "antenna"]
    doc.check(len(charges) == 2 and len(antennas) == 2, "the reference: two charges, two antennas")
    for k_, v, dg in (("Tzero", a.T0, 4), ("m", a.m, 4), ("c", a.c, 4), ("Ra", a.Ra, 4), ("Rq", a.Rq, 4),
                      ("xZero", a.x0[0], 4), ("yZero", a.x0[1], 4), ("xWorld", a.wmax[0], 4),
                      ("yWorld", a.wmax[1], 4), ("omega", w, 4), ("H", H, 4), ("tMax", a.tmax, 4),
                      ("dt", b.t0, 6), ("period", period, 6), ("vZero", a.v0, 10),
                      ("gammaZero", a.gamma0, 10), ("pZero", a.p0, 10),
                      ("xGLo", a.gmin[0], 4), ("yGLo", a.gmin[1], 4), ("xGHi", a.gmax[0], 4), ("yGHi", a.gmax[1], 4),
                      ("gAxis", lvl["gates"][0]["acceptance"]["direction"][0], 4),
                      ("gHalf", lvl["gates"][0]["acceptance"]["direction"][1], 4)):
        doc.value(k_, v, dg)
    doc.value("q", a.q, 4, sci=True)
    doc.value("Ey", dists[1]["e"][1], 4, sci=True)
    doc.text("distOff", dists[0]["name"])
    doc.text("distOn", dists[1]["name"])
    for st, nm in ((a, "One"), (b, "Two")):
        for k_, v in (("xLo", st.dmin[0]), ("yLo", st.dmin[1]), ("xHi", st.dmax[0]), ("yHi", st.dmax[1])):
            doc.value(f"{k_}{nm}", v, 4)
    for k_, v in (("xRegLo", lim["region"]["min"][0]), ("yRegLo", lim["region"]["min"][1]),
                  ("xRegHi", lim["region"]["max"][0]), ("yRegHi", lim["region"]["max"][1]),
                  ("maxQ", lim["max_charges"]), ("maxA", lim["max_antennas"])):
        doc.value(k_, v, 4)
    doc.text("magnitudes", ", ".join(rf"\num{{{x:g}}}" for x in lim["magnitudes"]))
    doc.text("amplitudes", ", ".join(rf"\num{{{x:g}}}" for x in lim["antenna_amplitudes"]))
    for i, e in enumerate(charges, 1):
        doc.value(f"xQ{i}", e["node"][0], 4)
        doc.value(f"yQ{i}", e["node"][1], 4)
        doc.value(f"Q{i}", e["value"], 4, sci=True)
    for i, e in enumerate(antennas, 1):
        doc.value(f"xA{i}", e["node"][0], 4)
        doc.value(f"yA{i}", e["node"][1], 4)
        doc.value(f"pA{i}", e["value"], 4, sci=True)
        doc.value(f"angA{i}", e["angle_deg"], 4)

    # The stray field's size (sec:stray): the transverse push and the drop it causes to the gate
    # and to the detectors' face along the undeflected line.
    on = Setup(lvl, [], 0, 1)
    F = on.q * on.E[1]
    t_gate = (on.gmin[0] - on.x0[0]) / on.v0
    t_det = (on.dmin[0] - on.x0[0]) / on.v0
    acc = F / (on.gamma0 * on.m)
    doc.value("force", F, 3, sci=True)
    doc.value("accel", acc, 3, sci=True)
    doc.value("tGate", t_gate, 4)
    doc.value("tDet", t_det, 4)
    doc.value("dropGate", abs(0.5 * acc * t_gate ** 2), 3)
    doc.value("dropDet", abs(0.5 * acc * t_det ** 2), 3)
    doc.value("angleGate", math.degrees(math.atan2(abs(F * t_gate), on.p0)), 3)

    # The four reference flights.
    largest = 0.0
    gates_rows = []
    flights = {}
    for si, sn in enumerate(SHOTS):
        for di, dn in enumerate(DISTS):
            st = Setup(lvl, ref, si, di)
            run, coarse = fly(st, H), fly(st, 2 * H)
            flights[(si, di)] = run
            nm = sn + dn
            doc.check(run["outcome"] == "Arrived", f"ion {si + 1}, {dists[di]['name']}: arrives")
            err = err_of(run, coarse)
            doc.value(f"t{nm}", run["t"], 10)
            doc.value(f"x{nm}", run["x"][0], 10)
            doc.value(f"y{nm}", run["x"][1], 10)
            doc.value(f"px{nm}", run["p"][0], 8)
            doc.value(f"py{nm}", run["p"][1], 8)
            doc.text(f"face{nm}", run["face"])
            doc.value(f"err{nm}", err, 2, sci=True)
            doc.value(f"rad{nm}", run["rad"], 2, sci=True)
            doc.data("path" + nm.lower(), run["path"])
            ents = run["gate"]["entries"]
            doc.check(run["gate"]["passed"] and len(ents) >= 1, f"ion {si + 1}, {dn}: passes the gate")
            ge = ents[-1]
            doc.check(len(ents) == 1, f"ion {si + 1}, {dn}: passes the gate at its first entry")
            doc.check(abs(abs(ge["theta"]) - ge["beta"]) < 1e-12, "the signed angle's size is the cone's deviation")
            gates_rows.append((rf"ion {si + 1}, {dists[di]['name']}", ge["t"], ge["x"][0], ge["x"][1], ge["face"],
                               math.degrees(ge["theta"]), math.degrees(st.ghalf - ge["beta"])))
            fi = 2 * si + di
            g, gp = doc.flight(0, flight=fi), doc.flight(0, flight=fi, run="preview")
            tol = max(TOL, 10 * err, end_gap(g["end"], gp["end"]))
            title = f"ion {si + 1}, {dists[di]['name']}: "
            doc.compare(f"Out{nm}", title + "outcome", run["outcome"], g["outcome"])
            doc.compare(f"T{nm}", title + r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
            doc.compare(f"X{nm}", title + r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
            doc.compare(f"Y{nm}", title + r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
            doc.compare(f"Px{nm}", title + r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
            doc.compare(f"Py{nm}", title + r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
            doc.compare(f"Rad{nm}", title + "radiated energy (relative)", run["rad"], g["radiated_energy"], RAD_REL, rel=True)
            largest = max(largest, end_gap(run, g["end"]))
    doc.table("gates", ["flight", r"$t_G$", r"$x_G$", r"$y_G$", "face", r"$\theta$ (deg)",
                        r"$\delta - \abs{\theta}$ (deg)"], gates_rows, digits=[4, 6, 6, 6, 4, 4, 3])
    # sec:stray: what the stray field changes at the gate, along the reference's paths.
    for si, sn in enumerate(SHOTS):
        off, on_ = gates_rows[2 * si], gates_rows[2 * si + 1]
        doc.value(f"dTheta{sn}", off[5] - on_[5], 4)
        doc.value(f"dY{sn}", off[3] - on_[3], 3)
        doc.value(f"thetaOff{sn}", off[5], 3)
        doc.value(f"thetaOn{sn}", on_[5], 3)
        doc.value(f"yGOff{sn}", off[3], 4)
        doc.value(f"yGOn{sn}", on_[3], 4)
    doc.value("coneWidth", 2 * math.degrees(a.ghalf), 3)
    doc.value("gateTop", min(a.gmax[1] - r[3] for r in gates_rows[0::2]), 2)
    doc.value("gateBottom", min(r[3] - a.gmin[1] for r in gates_rows[1::2]), 2)
    # The closest approach of ion 1, field on, to the negative charge (over the path's samples,
    # every 20 steps: about h v0 20 = 0.08 apart, so the minimum is a few 1e-3 high at most).
    qx, qy, _ = a.charges[1]
    doc.value("closeQ", min(math.hypot(p[0] - qx, p[1] - qy) for p in flights[(0, 1)]["path"]), 2)
    doc.value("gateHeight", a.gmax[1] - a.gmin[1], 3)
    doc.game_value("largest", largest, 2, sci=True)
    doc.value("minGateMargin", min(r[6] for r in gates_rows), 3)

    # eq:late: half a period late is the antennas reversed: ion 2 against itself launched at 0 with -p0.
    late = 0.0
    for di in range(2):
        rev = fly(Setup(lvl, ref, 1, di, t0=0.0, antenna_sign=-1.0), H)
        late = max(late, end_gap(rev, flights[(1, di)]))
    doc.check(late < 1e-9, "ion 2 flies as ion 1 would with the antennas reversed")
    doc.value("lateGap", late, 1, sci=True)

    # The empty board, in closed form; the reference's charges alone, marched (both ions alike).
    for di, dn in enumerate(DISTS):
        for si, sn in enumerate(SHOTS):
            st = Setup(lvl, [], si, di)
            out, tb, xb, pb = bare(st)
            gb = doc.flight(1, flight=2 * si + di)
            title = f"empty board, ion {si + 1}, {dists[di]['name']}: "
            doc.compare(f"BareOut{sn}{dn}", title + "outcome", out, gb["outcome"])
            doc.compare(f"BareT{sn}{dn}", title + r"$t_\mathrm{e}$", tb, gb["end"]["t"], TOL)
            if dn == "On":
                doc.compare(f"BareX{sn}{dn}", title + r"$x_\mathrm{e}$", xb[0], gb["end"]["x"][0], TOL)
                doc.compare(f"BarePy{sn}{dn}", title + r"$p_{y,\mathrm{e}}$", pb[1], gb["end"]["p"][1], TOL)
            if si == 0:
                doc.value(f"bareT{dn}", tb, 8)
                doc.value(f"bareX{dn}", xb[0], 6)
                doc.value(f"bareY{dn}", xb[1], 6)
        st = Setup(lvl, charges, 0, di)
        cr, cc = fly(st, H), fly(st, 2 * H)
        cerr = err_of(cr, cc)
        doc.text(f"chargesOut{dn}", {"LeftBounds": "leaves the world", "Arrived": "arrives"}.get(cr["outcome"], cr["outcome"]))
        doc.value(f"chargesY{dn}", cr["x"][1], 6)
        doc.text(f"chargesGate{dn}", "passes" if cr["gate"]["passed"] else "misses")
        for si, sn in enumerate(SHOTS):
            gc, gcp = doc.flight(2, flight=2 * si + di), doc.flight(2, flight=2 * si + di, run="preview")
            ctol = max(TOL, 10 * cerr, abs(gc["end"]["t"] - gcp["end"]["t"]), abs(gc["end"]["x"][1] - gcp["end"]["x"][1]))
            title = f"charges alone, ion {si + 1}, {dists[di]['name']}: "
            doc.compare(f"QOut{sn}{dn}", title + "outcome", cr["outcome"], gc["outcome"])
            doc.compare(f"QT{sn}{dn}", title + r"$t_\mathrm{e}$", cr["t"], gc["end"]["t"], ctol)
            doc.compare(f"QY{sn}{dn}", title + r"$y_\mathrm{e}$", cr["x"][1], gc["end"]["x"][1], ctol)
