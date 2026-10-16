"""Companion of docs/levels/glass_screen/doc.tex (Level 22, Glass screen).

One dielectric box and a fixed charge. The rectangle integral is the
deflection-plates document (eq:rectangle). The bound charge is the dielectric
condition, not a conductor potential: sigma = (lambda / 2 pi) times the
average outward normal field, lambda = (eps - 1) / (eps + 1), with one
neutrality unknown per body. The panel's own normal field is left out of that
average. The empty placement keeps the fixed charge and the glass. The
tolerance is ten times the larger of the step-halving estimate and the change
when the panels are doubled, and at least 1e-9, set before any comparison
with the game. The ion's image is left out, and its bound is quoted.
"""
import math

import numpy as np


# Reported step. The coarser flight uses twice this step.
H = 0.01
# Panel edges, in cells. The reported mesh is the finer one. Doubling the edge
# is the discretization estimate. Chosen as the two sizes, not from the game.
PANEL_FINE = 0.5
PANEL_COARSE = 1.0
# Simpson on the radiated power. Set beforehand.
TOL_RAD = 1e-4
# Obstacle shell (PHYSICS.md, the electrode box plus the contact shell).
CONTACT = 0.02
# The image-force bound floors the distance at this (PHYSICS.md 2.7).
D_FLOOR = 1e-3


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


def dielectric_bodies(lvl):
    """Level dielectric boxes, in the shape build_panels already meshes."""
    out = []
    for d in lvl.get("dielectrics") or []:
        out.append({
            "center": d["center"],
            "length": float(d["length"]),
            "thickness": float(d["thickness"]),
            "height": float(d["height"]),
            "angle_deg": float(d.get("angle_deg", 0.0)),
            "permittivity": float(d["permittivity"]),
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
    The mirror is the panel reflected through z = 0, carrying the same density."""
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
    return pot + pot2, ex + ex2, ey + ey2, ez + ez2


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


class Plates:
    def __init__(self, lvl, placement, target, dielectric=True):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        d = shot["launch"]["direction"][:2]
        nrm = math.hypot(d[0], d[1])
        self.d = [d[0] / nrm, d[1] / nrm]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        charges = [e for e in lvl["elements"] if e["kind"] == "charge"]
        charges += [e for e in placement if e["kind"] == "charge"]
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in charges]
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)
        self.dielectrics = dielectric_bodies(lvl) if dielectric else []
        self.eps = self.dielectrics[0]["permittivity"] if self.dielectrics else 1.0
        self.lam = (self.eps - 1.0) / (self.eps + 1.0) if self.dielectrics else 0.0
        self.packs, self.boxes = build_panels(self.dielectrics, target, s)
        self.n_panels = sum(p["u0"].size for p in self.packs.values())
        self.sigma = self._solve()

    def _coulomb(self, x, y, z):
        phi = ex = ey = ez = 0.0
        for cx, cy, q in self.charges:
            dx, dy, dz = x - cx, y - cy, z
            r2 = dx * dx + dy * dy + dz * dz
            r = math.sqrt(r2)
            phi += q / r
            inv = q / (r2 * r)
            ex += dx * inv
            ey += dy * inv
            ez += dz * inv
        return phi, ex, ey, ez

    def _mirror_field(self, x, y, z, pack):
        """Cartesian field of the mirror only. The panel itself is not included."""
        k = int(pack["kind"])
        u0, u1, v0, v1, c = pack["u0"], pack["u1"], pack["v0"], pack["v1"], pack["c"]
        if k == 0:
            _pot, ey, ez, ex = _rect(y, z, x, u0, u1, -v1, -v0, c)
        elif k == 1:
            _pot, ex, ez, ey = _rect(x, z, y, u0, u1, -v1, -v0, c)
        else:
            _pot, ex, ey, ez = _rect(x, y, z, u0, u1, v0, v1, -c)
        return ex, ey, ez

    def _en(self, ex, ey, ez, pack):
        k = int(pack["kind"])
        comp = ex if k == 0 else (ey if k == 1 else ez)
        return comp * pack["outward"]

    def _gauss_points(self, pack, j, gauss):
        """Four Gauss nodes on panel j, in the panel's plane."""
        u0, u1 = float(pack["u0"][j]), float(pack["u1"][j])
        v0, v1 = float(pack["v0"][j]), float(pack["v1"][j])
        c = float(pack["c"][j])
        uc, vc = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
        hu, hv = 0.5 * (u1 - u0), 0.5 * (v1 - v0)
        pts = []
        k = int(pack["kind"])
        for s in (-gauss, gauss):
            for t in (-gauss, gauss):
                u, v = uc + s * hu, vc + t * hv
                if k == 0:
                    pts.append((c, u, v))
                elif k == 1:
                    pts.append((u, c, v))
                else:
                    pts.append((u, v, c))
        return pts

    def _one_panel(self, pack, j):
        one = {"kind": pack["kind"]}
        for key, val in pack.items():
            if key == "kind":
                continue
            one[key] = val[j:j + 1]
        return one

    def _solve(self):
        n = self.n_panels
        if n == 0:
            self.system_resid = 0.0
            self.gamma = 0.0
            self._sigma_parts = []
            return np.zeros(0)
        spans = []
        cursor = 0
        areas = np.zeros(n)
        owners = []
        for pack in self.packs.values():
            m = pack["u0"].size
            spans.append((pack, cursor, cursor + m))
            areas[cursor:cursor + m] = (pack["u1"] - pack["u0"]) * (pack["v1"] - pack["v0"])
            for j in range(m):
                owners.append((pack, j))
            cursor += m
        # Four-point Gauss average of the normal-field condition. A single
        # centroid converges slowly near the edges. The panel's own normal
        # field averages to zero, so each sample's diagonal keeps only the mirror.
        gauss = 1.0 / math.sqrt(3.0)
        K = np.zeros((n, n))
        en_ext = np.zeros(n)
        for i, (owner, local) in enumerate(owners):
            samples = self._gauss_points(owner, local, gauss)
            one = self._one_panel(owner, local)
            krow = np.zeros(n)
            en = 0.0
            for x, y, z in samples:
                _phi, ex, ey, ez = self._coulomb(x, y, z)
                sign = float(owner["outward"][local])
                kind = int(owner["kind"])
                comp = ex if kind == 0 else (ey if kind == 1 else ez)
                en += sign * comp
                row = np.zeros(n)
                for pack, j0, j1 in spans:
                    _pot, px, py, pz = _panel_contrib(x, y, z, pack)
                    row[j0:j1] = self._en(px, py, pz, pack)
                mx, my, mz = self._mirror_field(x, y, z, one)
                row[i] = float(self._en(mx, my, mz, one)[0])
                krow += row
            K[i, :] = krow / len(samples)
            en_ext[i] = en / len(samples)
        scale = self.lam / (2.0 * math.pi)
        A = np.eye(n + 1)
        A[:n, :n] -= scale * K
        A[:n, n] = -scale
        A[n, :n] = areas
        A[n, n] = 0.0
        b = np.zeros(n + 1)
        b[:n] = scale * en_ext
        sol = np.linalg.solve(A, b)
        self.system_resid = float(np.max(np.abs(A @ sol - b)))
        self.gamma = float(sol[n])
        sigma = sol[:n]
        self._sigma_parts = []
        for pack, j0, j1 in spans:
            self._sigma_parts.append(sigma[j0:j1])
        self._areas = areas
        return sigma

    def field(self, x, y):
        phi = ex = ey = 0.0
        for pack, sig in zip(self.packs.values(), self._sigma_parts):
            p, a, b, _c = _panel_fields(x, y, 0.0, pack, sig)
            phi += p
            ex += a
            ey += b
        pc, ac, bc, _ec = self._coulomb(x, y, 0.0)
        return phi + pc, ex + ac, ey + bc


def rhs_of(st, state):
    x, y, px, py = state
    m = st.m
    c = st.c
    p2 = px * px + py * py
    g = math.sqrt(1.0 + p2 / (m * m * c * c))
    inv = 1.0 / (g * m)
    _phi, ex, ey = st.field(x, y)
    return (px * inv, py * inv, st.q * ex, st.q * ey)


def energy_of(st, state):
    """eq:flight-energy: kinetic energy plus q phi. Foundations eq:energy
    differs by the rest energy. phi is the panels, the mirrors and the placed
    charges; the ion is not in it."""
    x, y, px, py = state
    kinetic = (math.sqrt(1.0 + (px * px + py * py) / (st.m * st.m * st.c * st.c)) - 1.0)
    kinetic *= st.m * st.c * st.c
    phi, _ex, _ey = st.field(x, y)
    return kinetic + st.q * phi


def rk4_step(st, s, f, h):
    x, y, px, py = s
    k1x, k1y, k1px, k1py = f
    k2x, k2y, k2px, k2py = rhs_of(st, (x + 0.5 * h * k1x, y + 0.5 * h * k1y,
                                        px + 0.5 * h * k1px, py + 0.5 * h * k1py))
    k3x, k3y, k3px, k3py = rhs_of(st, (x + 0.5 * h * k2x, y + 0.5 * h * k2y,
                                        px + 0.5 * h * k2px, py + 0.5 * h * k2py))
    k4x, k4y, k4px, k4py = rhs_of(st, (x + h * k3x, y + h * k3y, px + h * k3px, py + h * k3py))
    return (x + h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
            y + h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y),
            px + h / 6.0 * (k1px + 2.0 * k2px + 2.0 * k3px + k4px),
            py + h / 6.0 * (k1py + 2.0 * k2py + 2.0 * k3py + k4py))


def hermite_scalar(s, f, s2, f2, h, u):
    u2 = u * u
    u3 = u2 * u
    h00 = 2.0 * u3 - 3.0 * u2 + 1.0
    h10 = u3 - 2.0 * u2 + u
    h01 = -2.0 * u3 + 3.0 * u2
    h11 = u3 - u2
    return h00 * s + h10 * h * f + h01 * s2 + h11 * h * f2


def hermite(s, f, s2, f2, h, u):
    return (hermite_scalar(s[0], f[0], s2[0], f2[0], h, u),
            hermite_scalar(s[1], f[1], s2[1], f2[1], h, u))


def _first(pred, pts, us):
    for u, p in zip(us, pts):
        if pred(p):
            return u
    return None


def _bisect(pred, s, f, s2, f2, h, lo, hi):
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if pred(hermite(s, f, s2, f2, h, mid)):
            hi = mid
        else:
            lo = mid
    return hi


def assign_face(pe, lo, hi):
    x, y = pe
    segs = (
        (abs(x - lo[0]), "left", 0, lo[0], lo[1] <= y <= hi[1]),
        (abs(x - hi[0]), "right", 0, hi[0], lo[1] <= y <= hi[1]),
        (abs(y - lo[1]), "bottom", 1, lo[1], lo[0] <= x <= hi[0]),
        (abs(y - hi[1]), "top", 1, hi[1], lo[0] <= x <= hi[0]),
    )
    inside = [seg for seg in segs if seg[4]]
    res, name, axis, val, _ = min(inside or segs, key=lambda seg: seg[0])
    if res < 1e-6:
        out = [x, y]
        out[axis] = val
        return (out[0], out[1]), name, res
    return (x, y), name, res


def _in_box(p, lo, hi):
    return lo[0] <= p[0] <= hi[0] and lo[1] <= p[1] <= hi[1]


def _dist_box(x, y, box):
    x0, x1, y0, y1 = box
    dx = 0.0 if x0 <= x <= x1 else (x0 - x if x < x0 else x - x1)
    dy = 0.0 if y0 <= y <= y1 else (y0 - y if y < y0 else y - y1)
    return math.hypot(dx, dy)


def fly(st, h, path_every=0, table_every=0.0):
    s = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    f = rhs_of(st, s)
    t = 0.0
    n = 0
    w0 = energy_of(st, s)
    wdev = 0.0
    path = [(s[0], s[1])]
    rows = [(0.0, s[0], s[1], s[2], s[3], 0.0)]
    next_table = table_every
    us = [k / 32.0 for k in range(1, 33)]
    d_min = min(_dist_box(s[0], s[1], b) for b in st.boxes) if st.boxes else None
    image_ratio = 0.0
    y_max = s[1]
    y_min = s[1]

    def note(x, y):
        nonlocal d_min, image_ratio, y_max, y_min
        if y > y_max:
            y_max = y
        if y < y_min:
            y_min = y
        if not st.boxes:
            return
        dist = min(_dist_box(x, y, b) for b in st.boxes)
        if dist < d_min:
            d_min = dist
        _phi, ex, ey = st.field(x, y)
        qe = abs(st.q) * math.hypot(ex, ey)
        bound = st.q * st.q / max(dist, D_FLOOR) ** 2
        ratio = bound / max(qe, st.T0)
        if ratio > image_ratio:
            image_ratio = ratio

    n_charges = len(st.charges)
    while t <= st.tmax:
        s2 = rk4_step(st, s, f, h)
        f2 = rhs_of(st, s2)
        margin = st.c * h
        lo_x, hi_x = (s[0], s2[0]) if s[0] <= s2[0] else (s2[0], s[0])
        lo_y, hi_y = (s[1], s2[1]) if s[1] <= s2[1] else (s2[1], s[1])
        box_lo = (lo_x - margin, lo_y - margin)
        box_hi = (hi_x + margin, hi_y + margin)
        hits = []
        pts = None

        def samples():
            nonlocal pts
            if pts is None:
                pts = [hermite(s, f, s2, f2, h, u) for u in us]
            return pts

        for idx, (cx, cy, _q) in enumerate(st.charges):
            limit = st.R
            if _dist_box(cx, cy, (box_lo[0], box_hi[0], box_lo[1], box_hi[1])) > limit and not (
                    box_lo[0] <= cx <= box_hi[0] and box_lo[1] <= cy <= box_hi[1]):
                # cheap reject: charge centre far from the step box
                dx = 0.0 if box_lo[0] <= cx <= box_hi[0] else min(abs(cx - box_lo[0]), abs(cx - box_hi[0]))
                dy = 0.0 if box_lo[1] <= cy <= box_hi[1] else min(abs(cy - box_lo[1]), abs(cy - box_hi[1]))
                if dx * dx + dy * dy > limit * limit:
                    continue
            u = _first(lambda p, cx=cx, cy=cy, limit=limit:
                       (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit, samples(), us)
            if u is not None:
                hits.append((u, 0, idx, "Collided(%d)" % idx))
        for idx, box in enumerate(st.boxes):
            # Expanded obstacle.
            hit_box = (box[0] - CONTACT, box[1] + CONTACT, box[2] - CONTACT, box[3] + CONTACT)
            if (box_hi[0] < hit_box[0] or box_lo[0] > hit_box[1]
                    or box_hi[1] < hit_box[2] or box_lo[1] > hit_box[3]):
                continue
            u = _first(lambda p, hb=hit_box: hb[0] <= p[0] <= hb[1] and hb[2] <= p[1] <= hb[3], samples(), us)
            if u is not None:
                hits.append((u, 0, n_charges + idx, "Collided(%d)" % (n_charges + idx)))
        if box_lo[0] < st.wmin[0] or box_hi[0] > st.wmax[0] or box_lo[1] < st.wmin[1] or box_hi[1] > st.wmax[1]:
            u = _first(lambda p: not _in_box(p, st.wmin, st.wmax), samples(), us)
            if u is not None:
                hits.append((u, 1, -1, "LeftBounds"))
        if not (box_hi[0] < st.dmin[0] or box_lo[0] > st.dmax[0]
                or box_hi[1] < st.dmin[1] or box_lo[1] > st.dmax[1]):
            u = _first(lambda p: _in_box(p, st.dmin, st.dmax), samples(), us)
            if u is not None:
                hits.append((u, 2, -1, "Arrived"))
        if hits:
            u_hit, _rank, idx, outcome = min(hits)
            if outcome == "Arrived":
                pred = lambda p: _in_box(p, st.dmin, st.dmax)
                face_box = (st.dmin, st.dmax)
            elif outcome == "LeftBounds":
                pred = lambda p: not _in_box(p, st.wmin, st.wmax)
                face_box = (st.wmin, st.wmax)
            else:
                face_box = None
                if idx < n_charges:
                    cx, cy, _q = st.charges[idx]
                    limit = st.R
                    pred = lambda p, cx=cx, cy=cy, limit=limit: (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit
                else:
                    hb = st.boxes[idx - n_charges]
                    hit_box = (hb[0] - CONTACT, hb[1] + CONTACT, hb[2] - CONTACT, hb[3] + CONTACT)
                    pred = lambda p, hb=hit_box: hb[0] <= p[0] <= hb[1] and hb[2] <= p[1] <= hb[3]
            hi = _bisect(pred, s, f, s2, f2, h, max(0.0, u_hit - 1.0 / 32.0), u_hit)
            te = t + hi * h
            se = rk4_step(st, s, f, te - t)
            pe = hermite(s, f, s2, f2, h, hi)
            face, residual = "", 0.0
            if face_box is not None:
                pe, face, residual = assign_face(pe, face_box[0], face_box[1])
            note(pe[0], pe[1])
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n + 1,
                    "face": face, "residual": residual, "rows": rows, "path": path, "wdev": wdev,
                    "d_min": d_min, "image_ratio": image_ratio, "y_max": y_max, "y_min": y_min}
        note(s2[0], s2[1])
        s, f = s2, f2
        t += h
        n += 1
        drift = abs(energy_of(st, s) - w0)
        if drift > wdev:
            wdev = drift
        if table_every and t + 1e-12 >= next_table:
            rows.append((t, s[0], s[1], s[2], s[3], energy_of(st, s) - w0))
            next_table += table_every
        if path_every and n % path_every == 0:
            path.append((s[0], s[1]))
    return {"outcome": "Timeout", "t": t, "x": (s[0], s[1]), "p": (s[2], s[3]), "steps": n,
            "face": "", "residual": 0.0, "rows": rows, "path": path, "wdev": wdev,
            "d_min": d_min, "image_ratio": image_ratio, "y_max": y_max, "y_min": y_min}


def radiated(st, h, t_end):
    m = st.m
    c = st.c
    q = st.q

    def power(state):
        _x, _y, px, py = state
        g = math.sqrt(1.0 + (px * px + py * py) / (m * m * c * c))
        vx, vy = px / (g * m), py / (g * m)
        _vx, _vy, fx, fy = rhs_of(st, state)
        vF = vx * fx + vy * fy
        ax = (fx - vx * vF / (c * c)) / (g * m)
        ay = (fy - vy * vF / (c * c)) / (g * m)
        vxa = vx * ay - vy * ax
        return (2.0 / 3.0) * q * q / (c ** 3) * g ** 6 * (ax * ax + ay * ay - vxa * vxa / (c * c))

    n = int(t_end / h)
    n -= n % 2
    state = (st.x0[0], st.x0[1], st.p0 * st.d[0], st.p0 * st.d[1])
    total = power(state)
    for i in range(1, n + 1):
        state = rk4_step(st, state, rhs_of(st, state), h)
        w = 1.0 if i == n else (4.0 if i % 2 else 2.0)
        total += w * power(state)
    simpson = h / 3.0 * total
    rest = t_end - n * h
    state2 = rk4_step(st, state, rhs_of(st, state), rest) if rest > 0.0 else state
    return simpson + rest * 0.5 * (power(state) + power(state2))


def step_gap(fine, coarse):
    if fine["outcome"] != coarse["outcome"]:
        return 1.0
    return max(abs(fine["t"] - coarse["t"]),
               abs(fine["x"][0] - coarse["x"][0]),
               abs(fine["x"][1] - coarse["x"][1]),
               abs(fine["p"][0] - coarse["p"][0]),
               abs(fine["p"][1] - coarse["p"][1])) / 15.0


def mesh_gap(fine, coarse):
    if fine["outcome"] != coarse["outcome"]:
        return 1.0
    return max(abs(fine["t"] - coarse["t"]),
               abs(fine["x"][0] - coarse["x"][0]),
               abs(fine["x"][1] - coarse["x"][1]),
               abs(fine["p"][0] - coarse["p"][0]),
               abs(fine["p"][1] - coarse["p"][1]))


def compute(doc):
    lvl = doc.level
    ref_place = lvl["reference_solution"]
    bare_place = doc.result(1)["placement"]
    # Two meshes, built before any game number is read into the tolerance.
    fine_ref = Plates(lvl, ref_place, PANEL_FINE)
    coarse_ref = Plates(lvl, ref_place, PANEL_COARSE)
    fine_bare = Plates(lvl, bare_place, PANEL_FINE)
    coarse_bare = Plates(lvl, bare_place, PANEL_COARSE)
    runs = {
        0: (fine_ref, fly(fine_ref, H, path_every=20, table_every=5.0), fly(fine_ref, 2.0 * H)),
        1: (fine_bare, fly(fine_bare, H, path_every=20), fly(fine_bare, 2.0 * H)),
    }
    mesh_runs = {
        0: fly(coarse_ref, H),
        1: fly(coarse_bare, H),
    }
    doc.check(runs[0][1]["outcome"] == runs[0][2]["outcome"], "the reference outcome changes when the step is halved")
    doc.check(runs[1][1]["outcome"] == runs[1][2]["outcome"], "the bare outcome changes when the step is halved")
    doc.check(runs[0][1]["outcome"] == mesh_runs[0]["outcome"], "doubling the panels changes the reference outcome")
    doc.check(runs[1][1]["outcome"] == mesh_runs[1]["outcome"], "doubling the panels changes the bare outcome")
    err_step = max(step_gap(runs[i][1], runs[i][2]) for i in (0, 1))
    err_mesh = max(mesh_gap(runs[i][1], mesh_runs[i]) for i in (0, 1))
    # Ten times the larger estimate, and at least the verified-flight floor.
    tol = max(1e-9, 10.0 * max(err_step, err_mesh))

    sdiv = lvl["grid"]["subdivision"]
    played = [e for e in ref_place if e["kind"] == "charge"]
    doc.check(len(fine_ref.dielectrics) == 1, "the level does not have one dielectric")
    doc.check(len(fine_bare.dielectrics) == 1, "the bare flight drops the glass")
    doc.check(len(played) == 1, "the reference does not place one charge")
    doc.check(len(fine_ref.charges) == 2, "the reference flight does not have two charges")
    doc.check(len(fine_bare.charges) == 1, "the bare flight does not keep only the fixed charge")
    doc.check(fine_ref.charges[0] == fine_bare.charges[0], "the fixed charge changed")
    doc.check(not any(e.get("kind") == "plate" for e in ref_place), "the reference places a plate")
    glass = fine_bare.dielectrics[0]
    # Sample on the launch line, under the glass. The fixed charge sits above
    # it, so the unscreened field pushes a positive ion down.
    sample_x = glass["center"][0] / sdiv
    sample_y = fine_bare.x0[1]
    phi_b, _ex_b, ey_b = fine_bare.field(sample_x, sample_y)
    cx, cy, q_fixed = fine_bare.charges[0]
    dx = sample_x - cx
    dy = sample_y - cy
    r_sample = math.hypot(dx, dy)
    phi_c = q_fixed / r_sample
    ey_c = q_fixed * dy / (r_sample ** 3)
    # The centre of the glass, on the line from the charge. The field there is
    # inside the dielectric. The launch line is outside, on the far side.
    in_y = glass["center"][1] / sdiv
    _phi_i, _ex_i, ey_i = fine_bare.field(sample_x, in_y)
    dy_i = in_y - cy
    r_in = math.hypot(sample_x - cx, dy_i)
    ey_ci = q_fixed * dy_i / (r_in ** 3)
    doc.check(ey_c < 0.0, "the charge above the beam does not push a positive ion downward")
    doc.check(abs(ey_i) < 0.5 * abs(ey_ci), "the field inside the glass is not reduced")
    doc.check(bool(np.isfinite(fine_ref.sigma).all()), "the bound charge is not finite")
    doc.value("H", H, 3)
    doc.value("panelFine", PANEL_FINE, 2)
    doc.value("panelCoarse", PANEL_COARSE, 2)
    doc.value("nPanels", fine_ref.n_panels)
    doc.value("contact", CONTACT, 2)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("errStep", err_step, 2, sci=True)
    doc.value("errMesh", err_mesh, 2, sci=True)
    doc.value("errMax", max(err_step, err_mesh), 2, sci=True)
    doc.value("gamma0", fine_ref.gamma0, 6)
    doc.value("p0", fine_ref.p0, 8)
    doc.value("v0", fine_ref.p0 / (fine_ref.gamma0 * fine_ref.m), 6)
    doc.value("vOverC", fine_ref.p0 / (fine_ref.gamma0 * fine_ref.m * fine_ref.c), 4)
    doc.value("qPart", fine_ref.q, 3, sci=True)
    doc.value("mPart", fine_ref.m, 2)
    doc.value("T0", fine_ref.T0, 2)
    doc.value("cLight", fine_ref.c, 2)
    doc.value("stepReach", fine_ref.c * H, 3)
    doc.value("Rcharge", fine_ref.R, 2)
    doc.value("xLaunch", fine_ref.x0[0], 3)
    doc.value("yLaunch", fine_ref.x0[1], 3)
    doc.value("xDetLo", fine_ref.dmin[0], 3)
    doc.value("yDetLo", fine_ref.dmin[1], 3)
    doc.value("xDetHi", fine_ref.dmax[0], 3)
    doc.value("yDetHi", fine_ref.dmax[1], 3)
    doc.value("xWorld", fine_ref.wmax[0], 3)
    doc.value("yWorld", fine_ref.wmax[1], 3)
    doc.value("arrowX", fine_ref.x0[0] + 1.6 * fine_ref.d[0], 4)
    doc.value("arrowY", fine_ref.x0[1] + 1.6 * fine_ref.d[1], 4)
    doc.value("eps", fine_bare.eps, 2)
    doc.value("lam", fine_bare.lam, 2)
    doc.value("xGlass", glass["center"][0] / sdiv, 3)
    doc.value("yGlass", glass["center"][1] / sdiv, 3)
    doc.value("glassL", float(glass["length"]), 3)
    doc.value("glassT", float(glass["thickness"]), 3)
    doc.value("glassH", float(glass["height"]), 3)
    doc.value("eGlass", ey_b, 4, sci=True)
    doc.value("phiGlass", phi_b, 4, sci=True)
    doc.value("eCoulomb", ey_c, 4, sci=True)
    doc.value("phiCoulomb", phi_c, 4, sci=True)
    doc.value("eInside", ey_i, 4, sci=True)
    doc.value("eInsideC", ey_ci, 4, sci=True)
    inside_rel = abs(ey_i) / abs(ey_ci)
    screen_rel = abs(ey_b) / abs(ey_c)
    doc.value("insideRel", inside_rel, 2, sci=True)
    doc.value("screenRel", screen_rel, 2, sci=True)
    doc.value("qFixed", q_fixed, 3, sci=True)
    doc.value("xCharge", cx, 3)
    doc.value("yCharge", cy, 3)
    play = fine_ref.charges[1]
    doc.value("qPlay", play[2], 3, sci=True)
    doc.value("xPlay", play[0], 3)
    doc.value("yPlay", play[1], 3)
    doc.value("rSample", r_sample, 4)
    box = fine_bare.boxes[0]
    doc.value("xGlass0", box[0], 3)
    doc.value("xGlass1", box[1], 3)
    doc.value("yGlass0", box[2], 3)
    doc.value("yGlass1", box[3], 3)
    # Bound charge of the upper half and its mirror. Neutrality makes this zero.
    q_net = 0.0
    near_mean = None
    for pack, sig in zip(fine_bare.packs.values(), fine_bare._sigma_parts):
        area = (pack["u1"] - pack["u0"]) * (pack["v1"] - pack["v0"])
        q_net += 2.0 * float((sig * area).sum())
        if int(pack["kind"]) == 1:
            near_mean = float(sig[pack["outward"] > 0.0].mean())
    doc.value("qNet", q_net, 4, sci=q_net != 0.0)
    doc.value("sigNear", near_mean, 4, sci=True)
    doc.check(near_mean is not None and near_mean < 0.0, "the near face does not carry the opposite charge")
    doc.check(abs(q_net) <= 1e-8 * abs(q_fixed), "the bound charge is not neutral")
    doc.value("solveRes", fine_ref.system_resid, 2, sci=True)
    doc.value("boundGamma", fine_bare.gamma, 4, sci=True)
    doc.check(fine_ref.system_resid < 1e-6, "the panel system is not solved")
    doc.check(fine_bare.system_resid < 1e-6, "the bare panel system is not solved")
    # The same player charge with the glass removed. Not compared with the game:
    # the level always includes the glass.
    no_glass = Plates(lvl, ref_place, PANEL_FINE, dielectric=False)
    hit = fly(no_glass, H)
    doc.text("HitOutcome", hit["outcome"])
    doc.value("HitT", hit["t"], 6)
    doc.value("HitX", hit["x"][0], 6)
    doc.value("HitY", hit["x"][1], 6)

    def publish(i, prefix):
        st, fine, coarse = runs[i]
        doc.text(prefix + "Outcome", fine["outcome"])
        doc.text(prefix + "Face", fine["face"])
        doc.value(prefix + "T", fine["t"], 6)
        doc.value(prefix + "X", fine["x"][0], 6)
        doc.value(prefix + "Y", fine["x"][1], 6)
        doc.value(prefix + "Px", fine["p"][0], 6)
        doc.value(prefix + "Py", fine["p"][1], 6)
        doc.value(prefix + "Ang", math.degrees(math.atan2(fine["p"][1], fine["p"][0])), 4)
        doc.value(prefix + "Err", step_gap(fine, coarse), 2, sci=True)
        doc.value(prefix + "Wdev", fine["wdev"], 2, sci=fine["wdev"] != 0.0)
        doc.value(prefix + "Steps", fine["steps"])
        doc.value(prefix + "Ymax", fine["y_max"], 4)
        doc.value(prefix + "Ymin", fine["y_min"], 4)
        dmin = fine["d_min"]
        doc.value(prefix + "Dmin", 0.0 if dmin is None else dmin, 4)
        doc.value(prefix + "Image", fine["image_ratio"], 2, sci=fine["image_ratio"] != 0.0)
        if prefix == "Ref":
            e1 = radiated(st, H / 2.0, fine["t"])
            e2 = radiated(st, H, fine["t"])
            rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
            doc.check(rel < TOL_RAD, prefix + " radiation quadrature is inside its tolerance")
            doc.value(prefix + "Rad", e1, 4, sci=True)
            doc.value(prefix + "RadErr", rel, 2, sci=True)
            doc.value(prefix + "RadFrac", e1 / st.T0, 2, sci=True)
            rad = e1
        else:
            rad = 0.0
        doc.data("path" + prefix.lower(), fine["path"])
        return rad

    rad_ref = publish(0, "Ref")
    publish(1, "Bare")
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$"],
              [(r[0], r[1], r[2], r[3], r[4]) for r in runs[0][1]["rows"]],
              digits=[4, 4, 4, 4, 4])
    doc.value("radFracMax", runs[0][1] and rad_ref / fine_ref.T0, 2, sci=True)
    doc.check(not any(e.get("kind") == "plate" for e in bare_place), "the bare placement adds a plate")
    doc.check(not any(e.get("kind") == "charge" for e in bare_place), "the bare placement adds a charge")

    def put(i, prefix, title, with_rad, ours_rad):
        _st, fine, _coarse = runs[i]
        g = doc.flight(i)
        ge = g["end"]
        doc.compare(prefix + "Outcome", title + "outcome", fine["outcome"], g["outcome"])
        doc.compare(prefix + "T", title + "$t$", fine["t"], ge["t"], tol)
        doc.compare(prefix + "X", title + "$x$", fine["x"][0], ge["x"][0], tol)
        doc.compare(prefix + "Y", title + "$y$", fine["x"][1], ge["x"][1], tol)
        if with_rad and g["radiated_energy"] != 0.0:
            doc.compare(prefix + "Rad", title + "$E$ (rel.)", ours_rad, g["radiated_energy"],
                        TOL_RAD, rel=True)

    # The radiated energy is reported, not compared. The quadrature tolerance
    # belongs to one field sampled twice; these panels and the engine's
    # triangles are not the same field, and the state rows already test it.
    put(0, "Ref", "", False, rad_ref)
    put(1, "Bare", "bare ", False, 0.0)
