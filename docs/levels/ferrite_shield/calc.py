"""Companion of docs/levels/ferrite_shield/doc.tex (Level 30, Ferrite shield).

One ferrite box (a magnetically soft insulator, mu = 1000) and magnets standing on the plane.
The magnetic problem is the dielectric one with H for E and mu for epsilon (the Glass screen
document): a bound magnetic surface charge sigma = (lambda/2 pi) H_n averaged over each
rectangle at four Gauss points, lambda = (mu - 1)/(mu + 1), here odd under the mirror z -> -z
(eq:odd), so no neutrality unknown. The rectangles' integral is the deflection plates'
(eq:rectangle). The flights by the classical Runge-Kutta method at a fixed step in B_z of the
magnets and the ferrite's charge. The tolerance: ten times the larger of the step-halving
estimate and the change when the rectangles' edge doubles, at least the game's own gap between
its preview and its verification, set before any comparison.
"""
import math

import numpy as np


H = 0.01                 # the reported step; the coarser flight uses twice it
PANEL_FINE = 0.25        # rectangles' edge, the reported mesh
PANEL_COARSE = 0.5       # its double, for the discretization estimate
TOL_RAD = 1e-4
CONTACT = 0.02           # the obstacle's contact shell (PHYSICS.md 2.6, 2.7)


def _edges(a, b, target):
    n = max(1, int(round((b - a) / target)))
    return np.linspace(a, b, n + 1)


def _pack(kind, u0, u1, v0, v1, c, electrode, outward):
    return {
        "kind": kind,
        "u0": np.asarray(u0, dtype=float),
        "u1": np.asarray(u1, dtype=float),
        "v0": np.asarray(v0, dtype=float),
        "v1": np.asarray(v1, dtype=float),
        "c": np.asarray(c, dtype=float),
        "electrode": np.asarray(electrode, dtype=int),
        "outward": np.asarray(outward, dtype=float),
    }


def _cat(parts):
    if not parts:
        return None
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "kind"}
    out["kind"] = parts[0]["kind"]
    return out


def build_panels(electrodes, target, subdivision):
    """Upper-half rectangles. kind 0, 1, 2: the normal is x, y, z.
    Tangent axes are (y, z), (x, z) and (x, y). Length is along x at angle 0."""
    groups = {0: [], 1: [], 2: []}
    boxes = []
    for e, el in enumerate(electrodes):
        s = subdivision
        cx = el["center"][0] / s
        cy = el["center"][1] / s
        ang = el["angle_deg"]
        along_x = int(round(ang / 90.0)) % 2 == 0
        L = el["length"]
        T = el["thickness"]
        hz = 0.5 * el["height"]
        hx, hy = (0.5 * L, 0.5 * T) if along_x else (0.5 * T, 0.5 * L)
        x0, x1 = cx - hx, cx + hx
        y0, y1 = cy - hy, cy + hy
        boxes.append((x0, x1, y0, y1))
        # Side faces. The upper half is z from 0 to hz; the mirror is added later.
        for sign, x in ((-1.0, x0), (1.0, x1)):
            ye, ze = _edges(y0, y1, target), _edges(0.0, hz, target)
            yy, zz = np.meshgrid(0.5 * (ye[:-1] + ye[1:]), 0.5 * (ze[:-1] + ze[1:]), indexing="ij")
            yb, zb = np.meshgrid(ye[:-1], ze[:-1], indexing="ij")
            yt, zt = np.meshgrid(ye[1:], ze[1:], indexing="ij")
            n = yy.size
            groups[0].append(_pack(0, yb.ravel(), yt.ravel(), zb.ravel(), zt.ravel(),
                                    np.full(n, x), np.full(n, e), np.full(n, sign)))
        for sign, y in ((-1.0, y0), (1.0, y1)):
            xe, ze = _edges(x0, x1, target), _edges(0.0, hz, target)
            xb, zb = np.meshgrid(xe[:-1], ze[:-1], indexing="ij")
            xt, zt = np.meshgrid(xe[1:], ze[1:], indexing="ij")
            n = xb.size
            groups[1].append(_pack(1, xb.ravel(), xt.ravel(), zb.ravel(), zt.ravel(),
                                    np.full(n, y), np.full(n, e), np.full(n, sign)))
        xe, ye = _edges(x0, x1, target), _edges(y0, y1, target)
        xb, yb = np.meshgrid(xe[:-1], ye[:-1], indexing="ij")
        xt, yt = np.meshgrid(xe[1:], ye[1:], indexing="ij")
        n = xb.size
        groups[2].append(_pack(2, xb.ravel(), xt.ravel(), yb.ravel(), yt.ravel(),
                                np.full(n, hz), np.full(n, e), np.full(n, 1.0)))
    packs = {}
    for k, v in groups.items():
        c = _cat(v)
        if c is not None:
            packs[k] = c
    return packs, boxes


def ferrite_bodies(lvl):
    """Level ferrite boxes, in the shape build_panels already meshes."""
    out = []
    for d in lvl.get("ferrites") or []:
        out.append({
            "center": d["center"],
            "length": float(d["length"]),
            "thickness": float(d["thickness"]),
            "height": float(d["height"]),
            "angle_deg": float(d.get("angle_deg", 0.0)),
            "permeability": float(d["permeability"]),
        })
    return out


def _accum(u, v, w):
    """Corner primitive P (eq:primitive). The potential uses |w|, so it is even;
    the normal field is odd. eq:rectangle is the corner sum built on this."""
    r = np.sqrt(u * u + v * v + w * w)
    ur = u + r
    vr = v + r
    lu = np.where(ur > 0.0, np.log(np.where(ur > 0.0, ur, 1.0)), 0.0)
    lv = np.where(vr > 0.0, np.log(np.where(vr > 0.0, vr, 1.0)), 0.0)
    ww = np.abs(w)
    denom = ww * r
    ang = np.where(
        denom > 0.0,
        np.arctan2(u * v, np.where(denom > 0.0, denom, 1.0)),
        np.where(u * v > 0.0, 0.5 * np.pi, np.where(u * v < 0.0, -0.5 * np.pi, 0.0)),
    )
    pot = u * lv + v * lu - ww * ang
    ez = np.where(w < 0.0, -ang, ang)
    return pot, -lv, -lu, ez


def _rect(x, y, z, x0, x1, y0, y1, z0):
    """Potential and field of unit density on rectangles in the plane z0."""
    w = z - z0
    pot = ex = ey = ez = 0.0
    for su, u in ((1.0, x - x0), (-1.0, x - x1)):
        for sv, v in ((1.0, y - y0), (-1.0, y - y1)):
            p, a, b, c = _accum(u, v, w)
            s = su * sv
            pot = pot + s * p
            ex = ex + s * a
            ey = ey + s * b
            ez = ez + s * c
    return pot, ex, ey, ez


def _panel_contrib(x, y, z, pack):
    """Per-panel potential and field, including the mirror. kind 0: normal x,
    tangents y, z. kind 1: normal y, tangents x, z. kind 2: normal z, tangents x, y.
    The mirror is the panel reflected through z = 0, carrying the opposite density:
    the bound magnetic charge is odd under the mirror (eq:odd)."""
    k = int(pack["kind"])
    u0, u1, v0, v1, c = pack["u0"], pack["u1"], pack["v0"], pack["v1"], pack["c"]
    if k == 0:
        pot, ey, ez, ex = _rect(y, z, x, u0, u1, v0, v1, c)
        pot2, ey2, ez2, ex2 = _rect(y, z, x, u0, u1, -v1, -v0, c)
    elif k == 1:
        pot, ex, ez, ey = _rect(x, z, y, u0, u1, v0, v1, c)
        pot2, ex2, ez2, ey2 = _rect(x, z, y, u0, u1, -v1, -v0, c)
    else:
        pot, ex, ey, ez = _rect(x, y, z, u0, u1, v0, v1, c)
        pot2, ex2, ey2, ez2 = _rect(x, y, z, u0, u1, v0, v1, -c)
    return pot - pot2, ex - ex2, ey - ey2, ez - ez2


def _panel_fields(x, y, z, pack, sigma):
    pot, ex, ey, ez = _panel_contrib(x, y, z, pack)
    return (pot * sigma).sum(), (ex * sigma).sum(), (ey * sigma).sum(), (ez * sigma).sum()


def centroids(pack):
    u = 0.5 * (pack["u0"] + pack["u1"])
    v = 0.5 * (pack["v0"] + pack["v1"])
    c = pack["c"]
    k = int(pack["kind"])
    if k == 0:
        return c, u, v
    if k == 1:
        return u, c, v
    return u, v, c




def dipole_h(mx, my, mu, x, y, z):
    """eq:dipole: the field of a moment mu z-hat at (mx, my, 0), at (x, y, z)."""
    rx, ry, rz = x - mx, y - my, z
    r2 = rx * rx + ry * ry + rz * rz
    r = math.sqrt(r2)
    r5 = r2 * r2 * r
    return mu * 3 * rz * rx / r5, mu * 3 * rz * ry / r5, mu * (3 * rz * rz - r2) / r5


def _rect_only(x, y, z, pack, j):
    """Rectangle j of a pack alone, without its mirror, at unit density."""
    k = int(pack["kind"])
    u0, u1, v0, v1, c = (pack[n][j:j + 1] for n in ("u0", "u1", "v0", "v1", "c"))
    if k == 0:
        pot, ey, ez, ex = _rect(y, z, x, u0, u1, v0, v1, c)
    elif k == 1:
        pot, ex, ez, ey = _rect(x, z, y, u0, u1, v0, v1, c)
    else:
        pot, ex, ey, ez = _rect(x, y, z, u0, u1, v0, v1, c)
    return float(ex[0]), float(ey[0]), float(ez[0])


class Shield:
    def __init__(self, lvl, placement, target, ferrite=True):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
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
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        mags = [e for e in lvl["elements"] if e["kind"] == "magnet"] + \
               [e for e in placement if e["kind"] == "magnet"]
        self.magnets = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in mags]
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)
        self.bodies = ferrite_bodies(lvl) if ferrite else []
        self.mu_r = self.bodies[0]["permeability"] if self.bodies else 1.0
        self.lam = (self.mu_r - 1.0) / (self.mu_r + 1.0)
        self.packs, self.boxes = build_panels(self.bodies, target, s) if ferrite else ({}, [])
        self.n_panels = sum(p["u0"].size for p in self.packs.values())
        self._solve()
        self._prepare_fast()

    def h_ext(self, x, y, z):
        hx = hy = hz = 0.0
        for mx, my, mu in self.magnets:
            a, b, c = dipole_h(mx, my, mu, x, y, z)
            hx, hy, hz = hx + a, hy + b, hz + c
        return hx, hy, hz

    def _solve(self):
        """eq:system: (I - lambda/2pi K) sigma = lambda/2pi <H_ext . n>, each row averaged at
        four Gauss points and projected on the receiving rectangle's outward normal."""
        n = self.n_panels
        self._parts = []
        self.sigma = np.zeros(n)
        self.resid = 0.0
        self._areas = np.zeros(n)
        if n == 0:
            return
        spans, owners, cur = [], [], 0
        for pack in self.packs.values():
            m = pack["u0"].size
            spans.append((pack, cur, cur + m))
            owners += [(pack, j) for j in range(m)]
            cur += m
        g = 1.0 / math.sqrt(3.0)
        K = np.zeros((n, n))
        rhs_v = np.zeros(n)
        for i, (own, j) in enumerate(owners):
            kind, sign = int(own["kind"]), float(own["outward"][j])
            u0, u1, v0, v1, c = (float(own[k][j]) for k in ("u0", "u1", "v0", "v1", "c"))
            uc, vc, hu, hv = 0.5 * (u0 + u1), 0.5 * (v0 + v1), 0.5 * (u1 - u0), 0.5 * (v1 - v0)
            row = np.zeros(n)
            for su in (-g, g):
                for sv in (-g, g):
                    u, v = uc + su * hu, vc + sv * hv
                    x, y, z = ((c, u, v) if kind == 0 else (u, c, v) if kind == 1 else (u, v, c))
                    hx, hy, hz = self.h_ext(x, y, z)
                    rhs_v[i] += sign * (hx if kind == 0 else hy if kind == 1 else hz) / 4.0
                    for pack, a, b in spans:
                        _p, px, py, pz = _panel_contrib(x, y, z, pack)
                        row[a:b] += sign * (px if kind == 0 else py if kind == 1 else pz) / 4.0
                    # The rectangle's own normal field averages to zero over its two sides:
                    # take its own part out (its mirror, carrying -sigma, stays).
                    ox, oy, oz = _rect_only(x, y, z, own, j)
                    row[i] -= sign * (ox if kind == 0 else oy if kind == 1 else oz) / 4.0
            K[i] = row
        scale = self.lam / (2.0 * math.pi)
        A = np.eye(n) - scale * K
        b = scale * rhs_v
        self.sigma = np.linalg.solve(A, b)
        self.resid = float(np.max(np.abs(A @ self.sigma - b)))
        self._parts = [self.sigma[a:b2] for _p, a, b2 in spans]
        self._areas = np.concatenate([(p["u1"] - p["u0"]) * (p["v1"] - p["v0"]) for p, _a, _b in spans])

    def _prepare_fast(self):
        """For the flights: B_z at z = 0 of every rectangle and its mirror, in one pass per
        orientation (the mirror's in-plane field cancels the rectangle's there)."""
        self._fast = []
        for pack, sig in zip(self.packs.values(), self._parts):
            k = int(pack["kind"])
            u0, u1, v0, v1, c = pack["u0"], pack["u1"], pack["v0"], pack["v1"], pack["c"]
            if k == 2:
                U0, U1, V0, V1 = (np.concatenate([a, a]) for a in (u0, u1, v0, v1))
                C = np.concatenate([c, -c])
            else:
                U0, U1 = np.concatenate([u0, u0]), np.concatenate([u1, u1])
                V0, V1 = np.concatenate([v0, -v1]), np.concatenate([v1, -v0])
                C = np.concatenate([c, c])
            S = np.concatenate([sig, -sig])
            self._fast.append((k, U0, U1, C, np.concatenate([V0, V1, V0, V1]),
                               np.concatenate([S, -S, -S, S])))

    def bz_ferrite(self, x, y):
        bz = 0.0
        for k, U0, U1, C, VV, SS in self._fast:
            if k == 0:
                U, V, W = y, 0.0, x
            elif k == 1:
                U, V, W = x, 0.0, y
            else:
                U, V, W = x, y, 0.0
            ua, ub = U - U0, U - U1
            w = W - C
            _pot, _eu, ev, ew = _accum(np.concatenate([ua, ua, ub, ub]), V - VV,
                                       np.concatenate([w, w, w, w]))
            bz += float(SS @ (ew if k == 2 else ev))
        return bz

    def bz_by_panels(self, x, y):
        bz = 0.0
        for pack, sig in zip(self.packs.values(), self._parts):
            _p, _x, _y, pz = _panel_contrib(x, y, 0.0, pack)
            bz += float((pz * sig).sum())
        return bz

    def bz(self, x, y):
        b = sum(-mu / math.hypot(x - mx, y - my) ** 3 for mx, my, mu in self.magnets)
        return b + (self.bz_ferrite(x, y) if self.n_panels else 0.0)


def rhs(st, s):
    x, y, px, py = s
    g = math.sqrt(1 + (px * px + py * py) / (st.m * st.m * st.c * st.c))
    vx, vy = px / (g * st.m), py / (g * st.m)
    b = st.bz(x, y)
    return (vx, vy, st.q * vy * b, -st.q * vx * b)


def rk4_step(st, s, f, h):
    k1 = f
    k2 = rhs(st, tuple(a + h / 2 * b for a, b in zip(s, k1)))
    k3 = rhs(st, tuple(a + h / 2 * b for a, b in zip(s, k2)))
    k4 = rhs(st, tuple(a + h * b for a, b in zip(s, k3)))
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
    """From the launch to the first event: a magnet's disc, the ferrite's box with its contact
    shell, the world's edge, the detector (obstacles first)."""
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = rhs(st, s)
    t, n = 0.0, 0
    p0 = math.hypot(s[2], s[3])
    rows, path, pdev = [], [(s[0], s[1])], 0.0
    dmin_box = float("inf")
    preds = []
    for idx, (mx, my, _) in enumerate(st.magnets):
        preds.append((f"Collided({idx})", 0, lambda p, mx=mx, my=my: math.hypot(p[0] - mx, p[1] - my) <= st.Rm))
    for idx, (x0, x1, y0, y1) in enumerate(st.boxes):
        hb = (x0 - CONTACT, x1 + CONTACT, y0 - CONTACT, y1 + CONTACT)
        preds.append((f"Collided({len(st.magnets) + idx})", 0,
                      lambda p, hb=hb: hb[0] <= p[0] <= hb[1] and hb[2] <= p[1] <= hb[3]))
    preds.append(("LeftBounds", 1, lambda p: not (st.wmin[0] <= p[0] <= st.wmax[0]
                                                  and st.wmin[1] <= p[1] <= st.wmax[1])))
    preds.append(("Arrived", 2, lambda p: box_distance(p, st.dmin, st.dmax) <= 0))
    us = [k / 32 for k in range(1, 33)]
    while True:
        hs = min(h, st.tmax - t)
        s2 = rk4_step(st, s, f, hs)
        f2 = rhs(st, s2)
        pts = [hermite(s, f, s2, f2, hs, u) for u in us]
        for (x0, x1, y0, y1) in st.boxes:
            for p in (pts[0], pts[15], pts[-1]):
                dx = max(x0 - p[0], 0.0, p[0] - x1)
                dy = max(y0 - p[1], 0.0, p[1] - y1)
                dmin_box = min(dmin_box, math.hypot(dx, dy))
        events = []
        for name, rank, pred in preds:
            for k, p in enumerate(pts):
                if pred(p):
                    events.append((us[k], rank, name, pred))
                    break
        if events:
            u_hit, _r, outcome, pred = min(events, key=lambda e: (e[0], e[1]))
            lo, hi = u_hit - 1 / 32, u_hit
            for _ in range(60):
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
                    "rows": rows, "path": path, "pdev": pdev, "dmin": dmin_box}
        s, f, t, n = s2, f2, t + hs, n + 1
        pdev = max(pdev, abs(math.hypot(s[2], s[3]) - p0))
        if record_every and n % record_every == 0:
            rows.append((t, s))
        if n % 20 == 0:
            path.append((s[0], s[1]))
        if hs < h:
            return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": s[2:], "rows": rows,
                    "path": path, "steps": n, "pdev": pdev, "face": "", "dmin": dmin_box}


def radiated(st, h, t_end):
    def power(s):
        x, y, px, py = s
        g = math.sqrt(1 + (px * px + py * py) / (st.m ** 2 * st.c ** 2))
        v = (px / (g * st.m), py / (g * st.m))
        _, _, fx, fy = rhs(st, s)
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
        s = rk4_step(st, s, rhs(st, s), h)
        vals.append(power(s))
    simpson = h / 3 * (vals[0] + vals[-1] + 4 * sum(vals[1:-1:2]) + 2 * sum(vals[2:-1:2]))
    rest = t_end - n * h
    s2 = rk4_step(st, s, rhs(st, s), rest)
    return simpson + rest * 0.5 * (vals[-1] + power(s2))


def gap_of(a, b, div=1.0):
    if a["outcome"] != b["outcome"]:
        return float("inf")
    return max(abs(a["t"] - b["t"]), abs(a["x"][0] - b["x"][0]), abs(a["x"][1] - b["x"][1]),
               abs(a["p"][0] - b["p"][0]), abs(a["p"][1] - b["p"][1])) / div


def compute(doc):
    lvl = doc.level
    ref = lvl["reference_solution"]
    empty = doc.result(1)["placement"]
    fine_ref = Shield(lvl, ref, PANEL_FINE)
    coarse_ref = Shield(lvl, ref, PANEL_COARSE)
    fine_bare = Shield(lvl, empty, PANEL_FINE)
    coarse_bare = Shield(lvl, empty, PANEL_COARSE)
    doc.check(len(fine_bare.bodies) == 1 and len(fine_bare.magnets) == 1, "one ferrite, one level magnet")
    doc.check(len(fine_ref.magnets) == 2, "the reference places one magnet")

    # The flights' field is the panel sum rearranged: the same numbers to rounding.
    worst = 0.0
    for st in (fine_ref, coarse_bare):
        for x, y in (tuple(st.x0), (st.magnets[0][0], st.x0[1]), (st.dmin[0], st.dmin[1])):
            fast, slow = st.bz_ferrite(x, y), st.bz_by_panels(x, y)
            worst = max(worst, abs(fast - slow) / max(1.0, abs(slow)))
    doc.check(worst < 1e-10, "the flights' field differs from the panel sum")
    # Odd mirror: the upper half's magnetic charge, and the whole body's (zero by symmetry).
    q_half = float((fine_bare.sigma * fine_bare._areas).sum())
    doc.value("qHalf", q_half, 4, sci=True)

    runs = {}
    for key, fine, coarse in (("Ref", fine_ref, coarse_ref), ("Bare", fine_bare, coarse_bare)):
        a, b, c_ = fly(fine, H), fly(fine, 2 * H), fly(coarse, H)
        runs[key] = (fine, a, b, c_)
        doc.check(a["outcome"] == b["outcome"] == c_["outcome"],
                  key + ": the outcome changes with the step or the mesh")
    err_step = max(gap_of(r[1], r[2], 15.0) for r in runs.values())
    err_mesh = max(gap_of(r[1], r[3]) for r in runs.values())

    def game_gap(i):
        a, b = doc.flight(i, run="preview")["end"], doc.flight(i)["end"]
        return max(abs(a["t"] - b["t"]), *(abs(u - v) for u, v in zip(a["x"], b["x"])),
                   *(abs(u - v) for u, v in zip(a["p"], b["p"])))
    # Ten times the larger estimate, at least the verified-flight floor and the game's own gap.
    tol = max(1e-9, 10.0 * max(err_step, err_mesh), game_gap(0), game_gap(1))
    doc.game_value("gapGame", max(game_gap(0), game_gap(1)), 2, sci=True)
    doc.value("errStep", err_step, 2, sci=True)
    doc.value("errMesh", err_mesh, 2, sci=True)
    doc.value("tolFlight", tol, 2, sci=True)

    st = fine_ref
    doc.value("q", st.q, 3, sci=True)
    doc.value("m", st.m, 3)
    doc.value("c", st.c, 3)
    doc.value("Tzero", st.T0, 3)
    doc.value("gammaZero", st.gamma0, 6)
    doc.value("pZero", st.p0, 10)
    doc.value("Rm", st.Rm, 3)
    doc.value("contact", CONTACT, 2)
    doc.value("xZero", st.x0[0], 3)
    doc.value("yZero", st.x0[1], 3)
    for k, v in (("xDetLo", st.dmin[0]), ("yDetLo", st.dmin[1]), ("xDetHi", st.dmax[0]),
                 ("yDetHi", st.dmax[1]), ("xWorld", st.wmax[0]), ("yWorld", st.wmax[1])):
        doc.value(k, v, 3)
    (m1x, m1y, mu1), (m2x, m2y, mu2) = st.magnets
    doc.value("xMagOne", m1x, 3)
    doc.value("yMagOne", m1y, 3)
    doc.value("muOne", mu1, 3, sci=True)
    doc.value("xMagTwo", m2x, 3)
    doc.value("yMagTwo", m2y, 3)
    doc.value("muTwo", mu2, 3, sci=True)
    body = st.bodies[0]
    x0b, x1b, y0b, y1b = st.boxes[0]
    doc.value("xF0", x0b, 4)
    doc.value("xF1", x1b, 4)
    doc.value("yF0", y0b, 4)
    doc.value("yF1", y1b, 4)
    doc.value("fL", body["length"], 3)
    doc.value("fT", body["thickness"], 3)
    doc.value("fH", body["height"], 3)
    doc.value("muR", st.mu_r, 4)
    doc.value("lam", st.lam, 6)
    doc.value("panelFine", PANEL_FINE, 2)
    doc.value("panelCoarse", PANEL_COARSE, 2)
    doc.value("nPanels", st.n_panels)
    doc.value("solveRes", st.resid, 2, sci=True)
    doc.value("H", H, 3)
    doc.value("stepReach", st.c * H, 3)

    # The field on the launch line and beyond: the level magnet with and without the ferrite.
    alone = Shield(lvl, empty, PANEL_FINE, ferrite=False)
    rows = []
    pts = [("launch line, under the magnet", m1x, st.x0[1]),
           ("2 lower", m1x, st.x0[1] - 2.0),
           ("4 lower", m1x, st.x0[1] - 4.0),
           ("under the ferrite's left end", x0b, st.x0[1]),
           ("3 beyond its right end", x1b + 3.0, st.x0[1])]
    for name, x, y in pts:
        b_with, b_without = fine_bare.bz(x, y), alone.bz(x, y)
        rows.append([name, f"({x:g}, {y:g})", b_without, b_with, b_with / b_without])
    doc.table("fieldtab", ["point", r"$(x, y)$", r"$B_z$, magnet alone", r"$B_z$, with the ferrite",
                           "ratio"], rows, digits=[None, None, 4, 4, 3])
    doc.value("shieldRatio", rows[0][4], 3)
    doc.value("bAlone", rows[0][2], 4)
    doc.value("bWith", rows[0][3], 4)
    doc.value("bInside", fine_bare.bz(m1x, 0.5 * (y0b + y1b)), 4)
    doc.value("bInsideAlone", alone.bz(m1x, 0.5 * (y0b + y1b)), 4)

    # The reference without the ferrite: not compared (the level always has it).
    nof = fly(Shield(lvl, ref, PANEL_FINE, ferrite=False), H)
    doc.text("noFerriteOut", nof["outcome"])
    doc.value("noFerriteX", nof["x"][0], 5)
    doc.value("noFerriteY", nof["x"][1], 5)
    doc.text("noFerriteFace", nof["face"])
    nof_bare = fly(alone, H)
    doc.value("bareAloneY", nof_bare["x"][1], 4)
    doc.value("bareAloneX", nof_bare["x"][0], 4)

    for key in ("Ref", "Bare"):
        s_, fine_run, half_run, coarse_run = runs[key]
        doc.text(key + "Outcome", fine_run["outcome"])
        doc.text(key + "Face", fine_run["face"])
        doc.value(key + "T", fine_run["t"], 8)
        doc.value(key + "X", fine_run["x"][0], 8)
        doc.value(key + "Y", fine_run["x"][1], 8)
        doc.value(key + "Px", fine_run["p"][0], 8)
        doc.value(key + "Py", fine_run["p"][1], 8)
        doc.value(key + "Ang", math.degrees(math.atan2(fine_run["p"][1], fine_run["p"][0])), 4)
        doc.value(key + "Steps", fine_run["steps"])
        doc.value(key + "Pdev", fine_run["pdev"], 2, sci=True)
        doc.value(key + "Dmin", fine_run["dmin"], 4)
        doc.data("path" + key.lower(), fine_run["path"])
    tab = fly(fine_ref, H, record_every=int(round(5 / H)))
    doc.table("states", [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$"],
              [[f"{t:.0f}", s[0], s[1], s[2], s[3]] for t, s in tab["rows"]],
              align="rrrrr", digits=[None, 7, 7, 7, 7])
    ref_run = runs["Ref"][1]
    e1 = radiated(fine_ref, H / 2, ref_run["t"])
    e2 = radiated(fine_ref, H, ref_run["t"])
    rel = abs(e1 - e2) / abs(e1) / 15.0
    doc.check(rel < TOL_RAD, "the radiation quadrature is inside its tolerance")
    doc.value("Erad", e1, 4, sci=True)
    doc.value("EradFrac", e1 / fine_ref.T0, 2, sci=True)

    for i, key, title in ((0, "Ref", ""), (1, "Bare", "empty board: ")):
        r = runs[key][1]
        g = doc.flight(i)
        doc.compare(key + "Out", title + "outcome", r["outcome"], g["outcome"])
        doc.compare(key + "Tc", title + r"$t_\mathrm{e}$", r["t"], g["end"]["t"], tol)
        doc.compare(key + "Xc", title + r"$x_\mathrm{e}$", r["x"][0], g["end"]["x"][0], tol)
        doc.compare(key + "Yc", title + r"$y_\mathrm{e}$", r["x"][1], g["end"]["x"][1], tol)
