"""Companion of docs/levels/synchrotron_light/doc.tex (Level 41, Synchrotron light).

Every number of the document, from the level file. The field is static and axially symmetric:
the coil's B_z (first_coil eq:loop-field) plus the dipole on the axis, -mu/rho^3; its gradient
analytically (eq:gradient, checked by central differences); the vector potential for the
canonical angular momentum (ramp_the_coil eq:vector-potential plus mu/rho^2). The motion: the
relativistic equation of motion with the Landau-Lifshitz force in this field (eq:ll-static),
marched by the classical Runge-Kutta method at a fixed step with three more states: the work of
the reaction (the radiated energy), Lienard's radiated energy (Foundations eq:larmor) and the
reaction's torque about the axis; the entry by cubic Hermite interpolation and bisection; the
error against the step 2h. Checks: kinetic energy plus the reaction's work constant
(eq:balance), the canonical angular momentum's change equal to the torque's integral
(eq:torque). Without radiation: the bound on how near the axis the ion can come (eq:bound), for
every allowed magnet. By hand: the adiabatic inspiral through concentric circular orbits
(eq:adiabatic), its time by Simpson's rule. Compared with the game: the reference, no magnet and
the five other magnets of placements.json.
"""
import math


TOL = 1e-9          # ten times the larger of our estimate and the game's (1e-10), or the game's own gap
RAD_REL = 1e-7      # Lienard's energy: the game's three-point Gauss-Legendre per step (steps ~1.6,
                    # the integrand varying at twice the orbit's frequency, ~0.3): (0.5)^6/2e6 ~ 1e-8; ten times
H = 0.02
N_SIMPSON = 2000


def agm_ke(m):
    """K(m), E(m) by the arithmetic-geometric mean (first_coil eq:agm)."""
    a, b = 1.0, math.sqrt(1.0 - m)
    csum, p = 0.5 * m, 1.0
    for _ in range(60):
        c = 0.5 * (a - b)
        a, b = 0.5 * (a + b), math.sqrt(a * b)
        p *= 2.0
        csum += 0.5 * p * c * c
        if abs(c) < 1e-17:
            break
    K = math.pi / (2.0 * a)
    return K, K * (1.0 - csum)


class Field:
    """B_z(rho), dB_z/drho and A_phi(rho) of the coil (kappa, radius a) and the axis dipole mu."""

    def __init__(self, kappa, a, mu):
        self.kappa, self.a, self.mu = kappa, a, mu

    def coil(self, rho):
        a, k = self.a, self.kappa
        al2, be = (a - rho) ** 2, a + rho
        m = 4 * a * rho / (be * be)
        K, E = agm_ke(m)
        N = (a * a - rho * rho) * E + al2 * K
        D = al2 * be
        b = 2 * k * N / D                                       # first_coil eq:loop-field
        # eq:gradient: dK/dm = (E - (1-m)K)/(2m(1-m)), dE/dm = (E - K)/(2m), dm/drho = 4a(a-rho)/beta^3
        dm = 4 * a * (a - rho) / be ** 3
        dK = (E - (1 - m) * K) / (2 * m * (1 - m)) * dm
        dE = (E - K) / (2 * m) * dm
        dN = -2 * rho * E + (a * a - rho * rho) * dE + 2 * (rho - a) * K + al2 * dK
        dD = 2 * (rho - a) * be + al2
        db = 2 * k * (dN * D - N * dD) / (D * D)
        h = (1 - m / 2) * K - E
        aphi = 2 * k * be * h / rho                              # ramp_the_coil eq:vector-potential
        return b, db, aphi

    def at(self, rho):
        b, db, aphi = self.coil(rho)
        mu = self.mu
        return b - mu / rho ** 3, db + 3 * mu / rho ** 4, aphi + mu / rho ** 2


class Setup:
    def __init__(self, lvl, placement, reaction=True):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        ph = lvl["physics"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.c = ph["c"]
        self.reaction = reaction and ph["radiation_reaction"]
        self.tmax = ph["t_max"]
        self.wire = ph["wire_radius"]
        self.Rm = ph["magnet_radius"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        n = math.hypot(*d)
        self.d = [d[0] / n, d[1] / n]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1, lvl["grid"]["ny"] / s + 1]
        (coil,) = lvl["coils"]
        self.cc = (coil["center"][0] / s, coil["center"][1] / s)
        mags = [e for e in list(lvl["elements"]) + list(placement) if e["kind"] == "magnet"]
        for e in mags:
            assert (e["node"][0] / s, e["node"][1] / s) == self.cc, "the magnet on the coil's axis"
        self.has_magnet = bool(mags)
        self.field = Field(coil["kappa"], coil["radius"] / s, sum(e["value"] for e in mags))
        self.gamma0 = 1 + self.T0 / (self.m * self.c ** 2)                # Foundations eq:launch
        self.p0 = self.m * self.c * math.sqrt(self.gamma0 ** 2 - 1)
        c, q, m = self.c, self.q, self.m
        self.k1 = 2 * q ** 3 / (3 * m * c ** 3)                            # eq:ll-static
        self.k2 = 2 * q ** 4 / (3 * m * m * c ** 3)

    def rhs(self, s, t):
        """State (x, y, px, py, W_LL, W_Lienard, torque integral, W_Lienard with the Lorentz
        force's acceleration alone, as the game reports it)."""
        x, y, px, py = s[:4]
        m, c, q = self.m, self.c, self.q
        dx, dy = x - self.cc[0], y - self.cc[1]
        rho = math.hypot(dx, dy)
        B, dB, _ = self.field.at(rho)
        g = math.sqrt(1 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        fLx, fLy = q * vy * B, -q * vx * B
        if self.reaction:
            vg = (vx * dx + vy * dy) / rho * dB                          # v . grad B_z
            fx = self.k1 * g * vg * vy - self.k2 * B * B * g * g * vx
            fy = -self.k1 * g * vg * vx - self.k2 * B * B * g * g * vy
        else:
            fx = fy = 0.0
        Fx, Fy = fLx + fx, fLy + fy
        vF = vx * Fx + vy * Fy
        ax, ay = (Fx - vx * vF / c ** 2) / (g * m), (Fy - vy * vF / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        P = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        vF = vx * fLx + vy * fLy
        ax, ay = (fLx - vx * vF / c ** 2) / (g * m), (fLy - vy * vF / c ** 2) / (g * m)
        vxa = vx * ay - vy * ax
        PL = 2 / 3 * q * q / c ** 3 * g ** 6 * (ax * ax + ay * ay - vxa * vxa / c ** 2)
        return (vx, vy, Fx, Fy, -(vx * fx + vy * fy), P, dx * fy - dy * fx, PL)

    def kinetic(self, px, py):
        p2 = px * px + py * py
        return p2 / (self.m * (math.sqrt(1 + p2 / (self.m * self.c) ** 2) + 1))

    def canonical(self, s):
        """eq:canonical: (r x p)_z + q rho A_phi about the axis."""
        dx, dy = s[0] - self.cc[0], s[1] - self.cc[1]
        rho = math.hypot(dx, dy)
        return dx * s[3] - dy * s[2] + self.q * rho * self.field.at(rho)[2]


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


def fly(st, h, record=0.0, checks=True):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1], 0.0, 0.0, 0.0, 0.0)
    f = st.rhs(s, 0.0)
    t, n = 0.0, 0
    j0 = st.canonical(s)
    wdev = jdev = 0.0
    rho_min = math.hypot(s[0] - st.cc[0], s[1] - st.cc[1])
    path = [(s[0], s[1])]
    rows = [(0.0,) + s[:4] + (st.T0,)] if record else []
    every = round(record / h) if record else 0
    cc, a = st.cc, st.field.a
    preds = []
    if st.has_magnet:
        preds.append(("Collided(magnet)", 0, lambda p: math.hypot(p[0] - cc[0], p[1] - cc[1]) <= st.Rm))
    preds += [("Collided(coil)", 0, lambda p: abs(math.hypot(p[0] - cc[0], p[1] - cc[1]) - a) <= st.wire),
              ("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0] and st.wmin[1] <= p[1] <= st.wmax[1])),
              ("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0)]
    us = [k / 64 for k in range(1, 65)]

    def check(s):
        nonlocal wdev, jdev
        if checks:
            wdev = max(wdev, abs(st.kinetic(s[2], s[3]) + s[4] - st.T0))
            jdev = max(jdev, abs(st.canonical(s) - j0 - s[6]))

    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, t, hs)
        f2 = st.rhs(s2, t + hs)
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
            se = rk4_step(st, s, f, t, te - t)
            pe = hermite(s, f, s2, f2, hs, hi)
            check(se)
            face = face_of(pe, st.dmin, st.dmax) if outcome == "Arrived" else ""
            path.append(pe)
            if record:
                rows.append((te, pe[0], pe[1], se[2], se[3], st.kinetic(se[2], se[3])))
            return {"outcome": outcome, "t": te, "x": pe, "p": se[2:4], "face": face, "steps": n + 1,
                    "wdev": wdev, "jdev": jdev, "j0": j0, "j": st.canonical(se), "wll": se[4], "wl": se[5],
                    "torque": se[6], "wlL": se[7], "path": path, "rows": rows, "rho_min": rho_min}
        s, f, t, n = s2, f2, t + hs, n + 1
        check(s)
        rho_min = min(rho_min, math.hypot(s[0] - st.cc[0], s[1] - st.cc[1]))
        if n % 25 == 0:
            path.append((s[0], s[1]))
        if record and n % every == 0:
            rows.append((t,) + s[:4] + (st.kinetic(s[2], s[3]),))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:4], "face": "", "steps": n,
                    "wdev": wdev, "jdev": jdev, "j0": j0, "j": st.canonical(s), "wll": s[4], "wl": s[5],
                    "torque": s[6], "wlL": s[7], "path": path, "rows": rows, "rho_min": rho_min}


def err_of(a, b):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / 15


def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rho_min_bound(st, j0):
    """eq:bound: without radiation |p| = p0 and J = j0, so a point at rho is reachable only if
    |j0 - q rho A_phi(rho)| <= rho p0. The smallest such rho, by stepping in from the launch
    radius and bisecting."""
    g = lambda r: abs(j0 - st.q * r * st.field.at(r)[2]) - r * st.p0
    r = math.hypot(st.x0[0] - st.cc[0], st.x0[1] - st.cc[1])
    step = 1e-3
    while r - step > 0.05 and g(r - step) <= 0:
        r -= step
    return bisect(g, r - step, r)


def adiabatic(st, r_from, r_to, n):
    """eq:adiabatic: concentric circular orbits p(r) = q r B(r); the energy E = sqrt(m^2c^4 + c^2p^2)
    falls at the rate P(r) = (2 q^4 B^2/3 m^4 c^3) p^2 (eq:power-circle), so the time from r_from to
    r_to is int dr (dE/dr)/P, by Simpson's rule with n intervals."""
    m, c, q = st.m, st.c, st.q

    def integrand(r):
        B, dB, _ = st.field.at(r)
        p = q * r * B
        E = math.sqrt(m * m * c ** 4 + c * c * p * p)
        dEdr = c * c * p * q * (B + r * dB) / E
        P = st.k2 * B * B * p * p / (m * m)
        return dEdr / P

    hh = (r_from - r_to) / n
    s = integrand(r_to) + integrand(r_from)
    s += 4 * sum(integrand(r_to + (2 * k - 1) * hh) for k in range(1, n // 2 + 1))
    s += 2 * sum(integrand(r_to + 2 * k * hh) for k in range(1, n // 2))
    return s * hh / 3


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    doc.check(lvl["physics"]["c"] is not None and lvl["physics"]["radiation_reaction"] and len(ref) == 1
              and ref[0]["kind"] == "magnet" and not lvl["elements"], "relativistic, with radiation reaction, one magnet")
    st = Setup(lvl, ref)
    fld = st.field
    for k_, v, dg in (("q", st.q, 4), ("m", st.m, 4), ("Tzero", st.T0, 4), ("c", st.c, 4), ("xZero", st.x0[0], 4),
                      ("yZero", st.x0[1], 4), ("xDetLo", st.dmin[0], 4), ("yDetLo", st.dmin[1], 4),
                      ("xDetHi", st.dmax[0], 4), ("yDetHi", st.dmax[1], 4), ("xC", st.cc[0], 4), ("yC", st.cc[1], 4),
                      ("a", fld.a, 4), ("kappa", fld.kappa, 4), ("mu", fld.mu, 4), ("tMax", st.tmax, 4),
                      ("wire", st.wire, 4), ("Rm", st.Rm, 4), ("H", H, 4)):
        doc.value(k_, v, dg)
    doc.value("gammaZero", st.gamma0, 10)
    doc.value("pZero", st.p0, 10)
    doc.value("vZero", st.p0 / (st.gamma0 * st.m), 10)
    tau0 = 2 * st.q ** 2 / (3 * st.m * st.c ** 3)
    doc.value("tauZero", tau0, 6)
    doc.value("kOne", st.k1, 6)
    r0 = math.hypot(st.x0[0] - st.cc[0], st.x0[1] - st.cc[1])
    doc.value("rZero", r0, 6)
    doc.value("rCorner", math.hypot(st.dmax[0] - st.cc[0], st.dmax[1] - st.cc[1]), 6)
    doc.value("rEdge", min(st.dmax[0] - st.cc[0], st.dmax[1] - st.cc[1]), 6)

    # The field: values and the gradient's check.
    worst = 0.0
    rows = []
    for rho in (2.0, 3.0, 4.0, 5.0, 6.0):
        B, dB, A = fld.at(rho)
        d = 1e-5
        num = (fld.at(rho + d)[0] - fld.at(rho - d)[0]) / (2 * d)
        worst = max(worst, abs(num - dB) / abs(dB))
        rows.append((rho, B, dB, A))
    doc.check(worst < 1e-8, "the gradient's closed form against central differences")
    doc.value("gradCheck", worst, 1, sci=True)
    doc.table("fieldtab", [r"$\rho$", r"$B_z$", r"$\dd B_z/\dd\rho$", r"$A_\phi$"], rows, digits=[3, 8, 8, 8])
    B6 = fld.at(r0)[0]
    doc.value("BZero", B6, 8)
    doc.value("pCircle", st.q * r0 * B6, 8)
    doc.value("rGuide", bisect(lambda r: st.q * r * fld.at(r)[0] - st.p0, 3.0, 8.0), 6)
    omega_c = st.q * B6 / (st.gamma0 * st.m)
    doc.value("omegaC", omega_c, 6)
    doc.value("period", 2 * math.pi / omega_c, 4)
    # the reaction against the Lorentz force: |f| / |q v x B| = (2 q^3 B gamma^2)/(3 m^2 c^3) at the launch
    doc.value("ratio", st.k2 * B6 * st.gamma0 ** 2 / st.q, 2)
    P0 = st.k2 * B6 ** 2 * st.p0 ** 2 / st.m ** 2
    doc.value("powerZero", P0, 6, sci=True)
    doc.value("lossTurn", P0 * 2 * math.pi / omega_c / st.T0 * 100, 3)

    # The reference flight.
    run, coarse = fly(st, H, record=100.0), fly(st, 2 * H, checks=False)
    doc.check(run["outcome"] == "Arrived", "the reference arrives")
    err = err_of(run, coarse)
    for k_, v, dg in (("tRef", run["t"], 10), ("xRef", run["x"][0], 10), ("yRef", run["x"][1], 10),
                      ("pxRef", run["p"][0], 10), ("pyRef", run["p"][1], 10), ("wLL", run["wll"], 8),
                      ("wLienard", run["wl"], 8), ("jZero", run["j0"], 8), ("jEnd", run["j"], 8),
                      ("torque", run["torque"], 8)):
        doc.value(k_, v, dg)
    doc.text("faceRef", run["face"])
    doc.value("errRef", err, 2, sci=True)
    doc.value("wDev", run["wdev"], 2, sci=True)
    doc.value("jDev", run["jdev"], 2, sci=True)
    doc.value("radFrac", run["wll"] / st.T0, 3)
    doc.value("tEnd", st.kinetic(*run["p"]), 6)
    doc.value("wLienardL", run["wlL"], 8)
    doc.value("lorentzOnly", abs(run["wlL"] - run["wl"]) / run["wl"], 2, sci=True)
    doc.value("schott", abs(run["wl"] - run["wll"]), 2, sci=True)
    doc.value("schottRel", abs(run["wl"] - run["wll"]) / run["wll"], 2, sci=True)
    doc.check(run["wdev"] < 1e-9 and run["jdev"] < 1e-9, "the energy and torque balances hold")
    doc.data("pathref", run["path"])
    doc.table("march", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$", r"$T$"],
              run["rows"],
              digits=[5, 7, 7, 7, 7, 7])

    # Without radiation: the bound for every allowed magnet and none.
    doc.value("rBoundRef", rho_min_bound(st, run["j0"]), 5)
    allowed = sorted(lvl["limits"]["magnet_strengths"])
    bounds = []
    for v in [0.0] + [s * v for v in allowed for s in (1, -1)]:
        sb = Setup(lvl, [] if v == 0 else [{"node": ref[0]["node"], "kind": "magnet", "value": v}])
        j0 = sb.canonical((sb.x0[0], sb.x0[1], sb.p0 * sb.d[0], sb.p0 * sb.d[1]))
        bounds.append((v, j0, rho_min_bound(sb, j0)))
    doc.check(all(b[2] > math.hypot(2, 2) for b in bounds), "without radiation no allowed magnet reaches the box")
    doc.table("bounds", [r"$\mu$", r"$J$", r"$\rho_\mathrm{min}$"], bounds, digits=[3, 7, 5])
    doc.value("rBoundWorst", min(b[2] for b in bounds), 5)
    # and the march without radiation agrees with the bound
    nr = fly(Setup(lvl, ref, reaction=False), 0.05, checks=False)
    doc.text("noRadOut", nr["outcome"])
    doc.value("noRadRhoMin", nr["rho_min"], 4)
    doc.check(nr["outcome"] == "Timeout" and nr["rho_min"] >= rho_min_bound(st, run["j0"]) - 1e-3,
              "the radiation-free march stays outside the bound")

    # By hand: the adiabatic inspiral.
    r_g = bisect(lambda r: st.q * r * fld.at(r)[0] - st.p0, 3.0, 8.0)
    for tag, r_end in (("Corner", math.hypot(2.0, 2.0)), ("Edge", 2.0)):
        ta = adiabatic(st, r_g, r_end, N_SIMPSON)
        tb = adiabatic(st, r_g, r_end, N_SIMPSON // 2)
        doc.check(abs(ta - tb) < 1e-6 * ta, "Simpson converged")
        doc.value(f"tAd{tag}", ta, 4)
    ta = adiabatic(st, r_g, math.hypot(2.0, 2.0), N_SIMPSON)
    doc.value("adRel", (run["t"] - ta) / run["t"] * 100, 2)

    # Compared with the game.
    g, gp = doc.flight(0), doc.flight(0, run="preview")
    gap = max(abs(gp["end"]["t"] - g["end"]["t"]), abs(gp["end"]["x"][0] - g["end"]["x"][0]),
              abs(gp["end"]["x"][1] - g["end"]["x"][1]), abs(gp["end"]["p"][0] - g["end"]["p"][0]),
              abs(gp["end"]["p"][1] - g["end"]["p"][1]))
    tol = max(TOL, 10 * err, gap)
    doc.compare("Out", "outcome", run["outcome"], g["outcome"])
    doc.compare("T", r"$t_\mathrm{e}$", run["t"], g["end"]["t"], tol)
    doc.compare("X", r"$x_\mathrm{e}$", run["x"][0], g["end"]["x"][0], tol)
    doc.compare("Y", r"$y_\mathrm{e}$", run["x"][1], g["end"]["x"][1], tol)
    doc.compare("Px", r"$p_{x,\mathrm{e}}$", run["p"][0], g["end"]["p"][0], tol)
    doc.compare("Py", r"$p_{y,\mathrm{e}}$", run["p"][1], g["end"]["p"][1], tol)
    doc.compare("Rad", "radiated energy (Liénard, Lorentz acceleration; relative)", run["wlL"],
                g["radiated_energy"], RAD_REL, rel=True)
    doc.game_value("tolRef", tol, 2, sci=True)
    doc.value("errRefTen", 10 * err, 2, sci=True)

    for k in range(1, doc.n_results()):
        placement = doc.result(k)["placement"]
        sk = Setup(lvl, placement)
        mu = placement[0]["value"] if placement else 0.0
        tag = "None" if not placement else ("M" if mu < 0 else "P") + str(abs(mu)).replace(".", "")
        rk, ck = fly(sk, H, checks=False), fly(sk, 2 * H, checks=False)
        ek = err_of(rk, ck)
        gk, gkp = doc.flight(k), doc.flight(k, run="preview")
        gapk = max(abs(gkp["end"]["t"] - gk["end"]["t"]), abs(gkp["end"]["x"][0] - gk["end"]["x"][0]),
                   abs(gkp["end"]["x"][1] - gk["end"]["x"][1]))
        tolk = max(TOL, 10 * ek, gapk)
        title = ("no magnet" if not placement else f"$\\mu = {mu:g}$") + ": "
        doc.compare(f"{tag}Out", title + "outcome", rk["outcome"], gk["outcome"])
        doc.compare(f"{tag}T", title + r"$t_\mathrm{e}$", rk["t"], gk["end"]["t"], tolk)
        doc.compare(f"{tag}X", title + r"$x_\mathrm{e}$", rk["x"][0], gk["end"]["x"][0], tolk)
        doc.compare(f"{tag}Y", title + r"$y_\mathrm{e}$", rk["x"][1], gk["end"]["x"][1], tolk)
        doc.value(f"t{tag}", rk["t"], 5)
        doc.value(f"rad{tag}", rk["wll"] / st.T0, 3)
        doc.value(f"rhoEnd{tag}", math.hypot(rk["x"][0] - st.cc[0], rk["x"][1] - st.cc[1]), 4)
        doc.text(f"out{tag}", rk["outcome"])
