"""Companion of docs/levels/microscope_column/doc.tex (Level 26, Microscope column).

Seven electrodes and the rectangle collocation of the deflection plates
(eq:rectangle, eq:collocation, eq:flight-energy). The condenser is three
apertures; the reference holds the middle pair at one potential. A deflector
supply, a fixed charge and one placed charge come after the gate. The gate
test is Two stages' rule, and the detector's direction test is the injection
cone. This file does not copy the engine.

The reported step is 0.05 rather than the deflection plates' 0.01 so that
three rays, on both placements and both meshes, finish under the calc time
limit. The mesh edges stay 0.5 and 1.0. The tolerance is ten times the larger
of the measured step-halving estimate and the measured panel-size estimate,
and at least 1e-9. It is fixed before any comparison with the game and is
not widened to match the game. The ion's image in the electrodes is left out,
as the game leaves it out, and its bound is quoted.
"""
import math

import numpy as np


# Reported step. The coarser flight uses twice this step.
H = 0.05
# Panel edges, in cells. The reported mesh is the finer one. Doubling the edge
# is the discretization estimate. The same two sizes as the deflection plates.
PANEL_FINE = 0.5
PANEL_COARSE = 1.0
# Simpson on the radiated power. Set beforehand, as on the deflection plates.
TOL_RAD = 1e-4
# Obstacle shell (PHYSICS.md, the electrode box plus the contact shell).
CONTACT = 0.02
# The image-force bound floors the distance at this (PHYSICS.md 2.7).
D_FLOOR = 1e-3


def _edges(a, b, target):
    n = max(1, int(round((b - a) / target)))
    return np.linspace(a, b, n + 1)


def _pack(kind, u0, u1, v0, v1, c, electrode):
    return {
        "kind": kind,
        "u0": np.asarray(u0, dtype=float),
        "u1": np.asarray(u1, dtype=float),
        "v0": np.asarray(v0, dtype=float),
        "v1": np.asarray(v1, dtype=float),
        "c": np.asarray(c, dtype=float),
        "electrode": np.asarray(electrode, dtype=int),
    }


def _cat(parts):
    if not parts:
        return None
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "kind"}
    out["kind"] = parts[0]["kind"]
    return out


def build_panels(electrodes, target, subdivision):
    """Upper-half rectangles. kind 0, 1, 2: the normal is x, y, z.
    Tangent axes are (y, z), (x, z) and (x, y). Length is along x at angle 0
    and along y when angle/90 is odd."""
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
        for _sign, x in ((-1.0, x0), (1.0, x1)):
            ye, ze = _edges(y0, y1, target), _edges(0.0, hz, target)
            yb, zb = np.meshgrid(ye[:-1], ze[:-1], indexing="ij")
            yt, zt = np.meshgrid(ye[1:], ze[1:], indexing="ij")
            n = yb.size
            groups[0].append(_pack(0, yb.ravel(), yt.ravel(), zb.ravel(), zt.ravel(),
                                    np.full(n, x), np.full(n, e)))
        for _sign, y in ((-1.0, y0), (1.0, y1)):
            xe, ze = _edges(x0, x1, target), _edges(0.0, hz, target)
            xb, zb = np.meshgrid(xe[:-1], ze[:-1], indexing="ij")
            xt, zt = np.meshgrid(xe[1:], ze[1:], indexing="ij")
            n = xb.size
            groups[1].append(_pack(1, xb.ravel(), xt.ravel(), zb.ravel(), zt.ravel(),
                                    np.full(n, y), np.full(n, e)))
        xe, ye = _edges(x0, x1, target), _edges(y0, y1, target)
        xb, yb = np.meshgrid(xe[:-1], ye[:-1], indexing="ij")
        xt, yt = np.meshgrid(xe[1:], ye[1:], indexing="ij")
        n = xb.size
        groups[2].append(_pack(2, xb.ravel(), xt.ravel(), yb.ravel(), yt.ravel(),
                                np.full(n, hz), np.full(n, e)))
    packs = {k: _cat(v) for k, v in groups.items()}
    return packs, boxes


def _accum(u, v, w):
    """Corner primitive P (the deflection plates' eq:primitive)."""
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
    """Per-panel potential and field, including the mirror."""
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
    def __init__(self, lvl, placement, target):
        shot = lvl["shots"][0]
        s = lvl["grid"]["subdivision"]
        self.lvl = lvl
        self.q = shot["particle"]["charge"]
        self.m = shot["particle"]["mass"]
        self.c = lvl["physics"]["c"]
        self.T0 = shot["launch"]["kinetic_energy"]
        self.R = lvl["physics"]["charge_radius"]
        self.tmax = lvl["physics"]["t_max"]
        self.x0 = [v / s for v in shot["launch"]["node"][:2]]
        self.d = [1.0, 0.0]
        self.dmin = [v / s for v in shot["detector"]["min"][:2]]
        self.dmax = [v / s for v in shot["detector"]["max"][:2]]
        self.wmin = [-1.0, -1.0]
        self.wmax = [lvl["grid"]["nx"] / s + 1.0, lvl["grid"]["ny"] / s + 1.0]
        direction = shot["detector"].get("acceptance", {}).get("direction")
        self.cone = direction is not None
        self.axis_deg = float(direction[0]) if self.cone else 0.0
        self.half_deg = float(direction[1]) if self.cone else 0.0
        gates = lvl.get("gates") or []
        self.gate = len(gates) == 1
        if self.gate:
            g = gates[0]
            self.gmin = [g["min"][0] / s, g["min"][1] / s]
            self.gmax = [g["max"][0] / s, g["max"][1] / s]
        else:
            self.gmin = None
            self.gmax = None
        charges = [e for e in lvl["elements"] if e["kind"] == "charge"]
        charges += [e for e in placement if e["kind"] == "charge"]
        self.charges = [(e["node"][0] / s, e["node"][1] / s, e["value"]) for e in charges]
        self.gamma0 = 1.0 + self.T0 / (self.m * self.c * self.c)
        self.p0 = math.sqrt(self.T0 * self.T0 / (self.c * self.c) + 2.0 * self.m * self.T0)
        self.packs, self.boxes = build_panels(lvl["electrodes"], target, s)
        volts = []
        for e in lvl["electrodes"]:
            kind = e["bias"]["kind"]
            volts.append(float(e["bias"]["value"]) if kind == "potential" else 0.0)
        for el in placement:
            if el["kind"] != "supply":
                continue
            hit = None
            for i, e in enumerate(lvl["electrodes"]):
                if e["center"][0] == el["node"][0] and e["center"][1] == el["node"][1]:
                    hit = i
            if hit is None:
                raise RuntimeError("supply is not on an electrode")
            volts[hit] = float(el["value"])
        self.voltages = np.array(volts, dtype=float)
        self.n_panels = sum(p["u0"].size for p in self.packs.values())
        self.sigma = self._solve()
        self.vacuum = len(self.charges) == 0 and float(np.max(np.abs(self.sigma))) == 0.0

    def set_direction(self, index):
        d = self.lvl["shots"][index]["launch"]["direction"][:2]
        nrm = math.hypot(d[0], d[1])
        self.d = [d[0] / nrm, d[1] / nrm]

    def _coulomb(self, x, y, z):
        phi = ex = ey = 0.0
        for cx, cy, q in self.charges:
            dx, dy, dz = x - cx, y - cy, z
            r2 = dx * dx + dy * dy + dz * dz
            r = math.sqrt(r2)
            phi += q / r
            inv = q / (r2 * r)
            ex += dx * inv
            ey += dy * inv
        return phi, ex, ey

    def _solve(self):
        n = self.n_panels
        G = np.zeros((n, n))
        b = np.zeros(n)
        points = []
        owners = []
        for pack in self.packs.values():
            xs, ys, zs = centroids(pack)
            for i in range(xs.size):
                points.append((float(xs[i]), float(ys[i]), float(zs[i])))
                owners.append(int(pack["electrode"][i]))
        spans = []
        cursor = 0
        for pack in self.packs.values():
            m = pack["u0"].size
            spans.append((pack, cursor, cursor + m))
            cursor += m
        for i, (x, y, z) in enumerate(points):
            phi_c, _, _ = self._coulomb(x, y, z)
            b[i] = self.voltages[owners[i]] - phi_c
            for pack, j0, j1 in spans:
                pot, _, _, _ = _panel_contrib(x, y, z, pack)
                G[i, j0:j1] = pot
        sigma = np.linalg.solve(G, b)
        self.system_resid = float(np.max(np.abs(G @ sigma - b)))
        self._sigma_parts = []
        for pack, j0, j1 in spans:
            self._sigma_parts.append(sigma[j0:j1])
        return sigma

    def field(self, x, y):
        if self.vacuum:
            return 0.0, 0.0, 0.0
        phi = ex = ey = 0.0
        for pack, sig in zip(self.packs.values(), self._sigma_parts):
            p, a, b, _c = _panel_fields(x, y, 0.0, pack, sig)
            phi += p
            ex += a
            ey += b
        pc, ac, bc = self._coulomb(x, y, 0.0)
        return phi + pc, ex + ac, ey + bc

    def electrode_charge(self):
        """Charge on each electrode, upper half plus the mirror."""
        totals = np.zeros(self.voltages.size)
        for pack, sig in zip(self.packs.values(), self._sigma_parts):
            area = (pack["u1"] - pack["u0"]) * (pack["v1"] - pack["v0"])
            np.add.at(totals, pack["electrode"], sig * area)
        return 2.0 * totals

    def surface_sample(self):
        """Potential on the long faces, in the plane z = 0, off the centroids.

        Aperture plates are long in y, so the large faces are x = x0 and
        x = x1. The deflector is long in x, so its large faces are the
        y faces. The sample stations are interior fractions of the long edge.
        """
        worst = 0.0
        for (x0, x1, y0, y1), v in zip(self.boxes, self.voltages):
            if (y1 - y0) > (x1 - x0):
                span = y1 - y0
                stations = [(x, y0 + frac * span) for frac in (0.17, 0.41, 0.73) for x in (x0, x1)]
            else:
                span = x1 - x0
                stations = [(x0 + frac * span, y) for frac in (0.17, 0.41, 0.73) for y in (y0, y1)]
            for x, y in stations:
                phi, _, _ = self.field(x, y)
                worst = max(worst, abs(phi - float(v)))
        return worst


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
    """The deflection plates' flight energy: kinetic energy plus q phi."""
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
    d_min = min(_dist_box(s[0], s[1], b) for b in st.boxes)
    image_ratio = 0.0
    y_max = s[1]
    y_min = s[1]
    gate_passed = bool(st.gate) and _in_box((s[0], s[1]), st.gmin, st.gmax)
    gate_info = None

    def mark_gate(s0, f0, s1, f1, t0, u_gate, u_hit):
        nonlocal gate_passed, gate_info
        if not st.gate or gate_passed or u_gate is None:
            return
        if u_hit is not None and u_gate >= u_hit:
            return
        lo = max(0.0, u_gate - 1.0 / 32.0)
        hi_g = _bisect(lambda p: _in_box(p, st.gmin, st.gmax), s0, f0, s1, f1, h, lo, u_gate)
        if u_hit is not None and hi_g >= u_hit:
            return
        te_g = t0 + hi_g * h
        se_g = rk4_step(st, s0, f0, te_g - t0)
        pe_g, face_g, res_g = assign_face(hermite(s0, f0, s1, f1, h, hi_g), st.gmin, st.gmax)
        gate_passed = True
        gate_info = {"t": te_g, "x": pe_g, "p": (se_g[2], se_g[3]), "face": face_g, "residual": res_g}

    def note(x, y):
        nonlocal d_min, image_ratio, y_max, y_min
        if y > y_max:
            y_max = y
        if y < y_min:
            y_min = y
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
            dx = 0.0 if box_lo[0] <= cx <= box_hi[0] else min(abs(cx - box_lo[0]), abs(cx - box_hi[0]))
            dy = 0.0 if box_lo[1] <= cy <= box_hi[1] else min(abs(cy - box_lo[1]), abs(cy - box_hi[1]))
            if dx * dx + dy * dy > limit * limit:
                continue
            u = _first(lambda p, cx=cx, cy=cy, limit=limit:
                       (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit, samples(), us)
            if u is not None:
                hits.append((u, 0, idx, "Collided(%d)" % idx))
        for idx, box in enumerate(st.boxes):
            hit_box = (box[0] - CONTACT, box[1] + CONTACT, box[2] - CONTACT, box[3] + CONTACT)
            if (box_hi[0] < hit_box[0] or box_lo[0] > hit_box[1]
                    or box_hi[1] < hit_box[2] or box_lo[1] > hit_box[3]):
                continue
            u = _first(lambda p, hb=hit_box: hb[0] <= p[0] <= hb[1] and hb[2] <= p[1] <= hb[3],
                       samples(), us)
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
        u_gate = None
        if st.gate and not gate_passed:
            if not (box_hi[0] < st.gmin[0] or box_lo[0] > st.gmax[0]
                    or box_hi[1] < st.gmin[1] or box_lo[1] > st.gmax[1]):
                u_gate = _first(lambda p: _in_box(p, st.gmin, st.gmax), samples(), us)
        if hits:
            u_hit, _rank, idx, outcome = min(hits)
            mark_gate(s, f, s2, f2, t, u_gate, u_hit)
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
                    pred = lambda p, cx=cx, cy=cy, limit=limit: (
                        (p[0] - cx) ** 2 + (p[1] - cy) ** 2 <= limit * limit)
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
            beta = None
            if outcome == "Arrived" and st.gate and not gate_passed:
                outcome = "SkippedGate(0)"
            elif outcome == "Arrived" and st.cone:
                axis = math.radians(st.axis_deg)
                ex, ey = math.cos(axis), math.sin(axis)
                px, py = se[2], se[3]
                beta = math.atan2(abs(px * ey - py * ex), px * ex + py * ey)
                if beta > math.radians(st.half_deg):
                    outcome = "Rejected"
            note(pe[0], pe[1])
            path.append(pe)
            return {"outcome": outcome, "t": te, "x": pe, "p": (se[2], se[3]), "steps": n + 1,
                    "face": face, "residual": residual, "rows": rows, "path": path, "wdev": wdev,
                    "d_min": d_min, "image_ratio": image_ratio, "y_max": y_max, "y_min": y_min,
                    "beta": beta, "gate_passed": gate_passed, "gate": gate_info}
        mark_gate(s, f, s2, f2, t, u_gate, None)
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
            "d_min": d_min, "image_ratio": image_ratio, "y_max": y_max, "y_min": y_min,
            "beta": None, "gate_passed": gate_passed, "gate": gate_info}


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


def _bodies_match(lvl):
    shots = lvl["shots"]
    s0 = shots[0]
    for s in shots[1:]:
        if s["particle"] != s0["particle"]:
            return False
        if s["launch"]["node"] != s0["launch"]["node"]:
            return False
        if s["launch"]["kinetic_energy"] != s0["launch"]["kinetic_energy"]:
            return False
        if s["detector"] != s0["detector"]:
            return False
    return True


def _angle(shot):
    d = shot["launch"]["direction"]
    return math.degrees(math.atan2(d[1], d[0]))

def _coulomb_one(x, y, cx, cy, q):
    dx, dy = x - cx, y - cy
    r2 = dx * dx + dy * dy
    r = math.sqrt(r2)
    inv = q / (r2 * r)
    return q / r, dx * inv, dy * inv


def compute(doc):
    lvl = doc.level
    ref_place = lvl["reference_solution"]
    bare_place = doc.result(1)["placement"]
    built = [e for e in lvl["elements"] if e["kind"] == "charge"]
    steer = [e for e in ref_place if e["kind"] == "charge"]
    supplies = [e for e in ref_place if e["kind"] == "supply"]
    doc.check(len(lvl["shots"]) == 3, "the level does not launch three rays")
    doc.check(_bodies_match(lvl), "the three shots do not share the particle, launch, energy and detector")
    doc.check(len(lvl["electrodes"]) == 7, "the level does not have seven electrodes")
    doc.check(len(lvl.get("gates") or []) == 1, "the level does not list one gate")
    doc.check(len(built) == 1, "the level does not build in one charge")
    doc.check(len(steer) == 1, "the reference does not place one charge")
    doc.check(len(supplies) == 3, "the reference does not set three supplies")
    doc.check(not bare_place, "the empty placement is not empty")
    acc = lvl["shots"][0]["detector"]["acceptance"]["direction"]
    doc.check(float(acc[0]) == 0.0, "the cone axis is not +x")

    fine = Plates(lvl, ref_place, PANEL_FINE)
    coarse = Plates(lvl, ref_place, PANEL_COARSE)
    bare = Plates(lvl, bare_place, PANEL_FINE)
    bare_c = Plates(lvl, bare_place, PANEL_COARSE)
    doc.check(len(fine.charges) == 2, "the reference flight does not carry two charges")
    doc.check(len(bare.charges) == 1, "the empty flight does not keep the fixed charge")
    doc.check(fine.gate and fine.cone, "the flight has no gate or no cone")
    doc.check(fine.system_resid < 1e-6 and bare.system_resid < 1e-6, "a panel system is not solved")

    mids = []
    defs = []
    for el in supplies:
        hit = None
        for i, electrode in enumerate(lvl["electrodes"]):
            if electrode["center"][0] == el["node"][0] and electrode["center"][1] == el["node"][1]:
                hit = i
        doc.check(hit is not None, "a supply is not on an electrode")
        doc.check(float(fine.voltages[hit]) == float(el["value"]), "a supply did not set its electrode")
        if int(round(lvl["electrodes"][hit]["angle_deg"] / 90.0)) % 2 == 1:
            mids.append((hit, float(el["value"])))
        else:
            defs.append((hit, float(el["value"])))
    doc.check(len(mids) == 2 and len(defs) == 1, "the supplies are not two apertures and one deflector")
    doc.check(mids[0][1] == mids[1][1] and mids[0][1] > 0.0, "the middle aperture is not one positive potential")
    doc.check(defs[0][1] > 0.0, "the deflector supply is not positive")
    raised = {mids[0][0], mids[1][0], defs[0][0]}
    doc.check(all(float(fine.voltages[i]) == 0.0 for i in range(7) if i not in raised),
              "an outer electrode is not grounded")
    doc.check(float(np.max(np.abs(bare.voltages))) == 0.0, "the empty placement raises an electrode")

    ref = []
    ref_step = []
    ref_mesh = []
    for k in range(3):
        fine.set_direction(k)
        ref.append(fly(fine, H, path_every=10, table_every=(5.0 if k == 1 else 0.0)))
        ref_step.append(fly(fine, 2.0 * H))
        coarse.set_direction(k)
        ref_mesh.append(fly(coarse, H))
    off = []
    off_step = []
    off_mesh = []
    for k in range(3):
        bare.set_direction(k)
        off.append(fly(bare, H, path_every=10))
        off_step.append(fly(bare, 2.0 * H))
        bare_c.set_direction(k)
        off_mesh.append(fly(bare_c, H))

    names = ("lower", "axis", "upper")
    half = math.radians(fine.half_deg)

    def _ended(label, run):
        beta = run["beta"]
        ang = math.degrees(beta) if beta is not None else float("nan")
        return "%s ends %s %s at %.6g, %.6g, cone %.4g deg, gate %s" % (
            label, run["outcome"], run["face"], run["x"][0], run["x"][1], ang,
            "yes" if run["gate_passed"] else "no")

    for k, name in enumerate(names):
        doc.check(ref[k]["outcome"] == ref_step[k]["outcome"],
                  "the " + name + " reference outcome changes when the step is doubled")
        doc.check(ref[k]["outcome"] == ref_mesh[k]["outcome"],
                  "doubling the panels changes the " + name + " reference outcome")
        doc.check(off[k]["outcome"] == off_step[k]["outcome"],
                  "the " + name + " empty outcome changes when the step is doubled")
        doc.check(off[k]["outcome"] == off_mesh[k]["outcome"],
                  "doubling the panels changes the " + name + " empty outcome")
        doc.check(ref[k]["outcome"] == "Arrived" and ref[k]["face"] == "left", _ended(name + " reference", ref[k]))
        doc.check(ref[k]["gate_passed"] and ref[k]["gate"] is not None, _ended(name + " gate", ref[k]))
        doc.check(ref[k]["beta"] is not None and ref[k]["beta"] <= half, _ended(name + " cone", ref[k]))
    doc.check(off[0]["outcome"] == "LeftBounds" and off[0]["face"] == "right", _ended("lower empty", off[0]))
    doc.check(off[1]["outcome"] == "LeftBounds" and off[1]["face"] == "right", _ended("axis empty", off[1]))
    collided = "Collided(%d)" % (len(bare.charges) + defs[0][0])
    doc.check(off[2]["outcome"] == collided, _ended("upper empty", off[2]))

    err_step = max([step_gap(ref[k], ref_step[k]) for k in range(3)]
                   + [step_gap(off[k], off_step[k]) for k in range(3)])
    # Both placements carry electrode charge, so every ray enters the panel estimate.
    err_mesh = max([mesh_gap(ref[k], ref_mesh[k]) for k in range(3)]
                   + [mesh_gap(off[k], off_mesh[k]) for k in range(3)])
    # Ten times the larger estimate, and at least the verified-flight floor.
    # Fixed here, before any game number is compared.
    tol = max(1e-9, 10.0 * max(err_step, err_mesh))

    sdiv = lvl["grid"]["subdivision"]
    y_axis = fine.x0[1]
    mid_boxes = [fine.boxes[i] for i, _v in mids]
    upper = [b for b in mid_boxes if b[2] > y_axis]
    lower = [b for b in mid_boxes if b[3] < y_axis]
    doc.check(len(upper) == 1 and len(lower) == 1, "the middle electrodes are not one aperture")
    y_rim_hi = upper[0][2]
    y_rim_lo = lower[0][3]
    opening = y_rim_hi - y_rim_lo
    x_mid = 0.5 * (upper[0][0] + upper[0][1])
    y_above = 0.5 * (y_axis + y_rim_hi)
    y_below = 0.5 * (y_axis + y_rim_lo)
    phi_launch, _ex_l, ey_launch = fine.field(fine.x0[0], y_axis)
    phi_mid, _ex_m, ey_axis = fine.field(x_mid, y_axis)
    _phi_a, _ex_a, ey_above = fine.field(x_mid, y_above)
    _phi_b, _ex_b, ey_below = fine.field(x_mid, y_below)
    outer_x = []
    for i, voltage in enumerate(fine.voltages):
        if float(voltage) != 0.0:
            continue
        x = 0.5 * (fine.boxes[i][0] + fine.boxes[i][1])
        if all(abs(x - u) > 1e-9 for u in outer_x):
            outer_x.append(x)
    outer_x.sort()
    doc.check(len(outer_x) == 2, "the grounded electrodes are not two apertures")
    phi_in, _ex_i, _ey_i = fine.field(outer_x[0], y_axis)
    phi_out, _ex_o, _ey_o = fine.field(outer_x[1], y_axis)
    t_mid = fine.T0 - fine.q * (phi_mid - phi_launch)
    dbox = fine.boxes[defs[0][0]]
    x_def = 0.5 * (dbox[0] + dbox[1])
    y_under = dbox[2] - 0.5
    _phi_d, _ex_d, ey_def = fine.field(x_def, y_under)
    totals = fine.electrode_charge()
    q_mid = float(totals[mids[0][0]] + totals[mids[1][0]])
    q_def = float(totals[defs[0][0]])
    q_outer = float(sum(totals[i] for i in range(7) if i not in raised))
    bx = built[0]["node"][0] / sdiv
    by = built[0]["node"][1] / sdiv
    bq = float(built[0]["value"])
    sx = steer[0]["node"][0] / sdiv
    sy = steer[0]["node"][1] / sdiv
    sq = float(steer[0]["value"])
    phi_pump, _ex_p, ey_pump = _coulomb_one(bx, by + 1.0, bx, by, bq)
    _phi_tot, _ex_tot, ey_over = fine.field(bx, by + 1.0)
    surf = fine.surface_sample()
    doc.check(math.isfinite(surf), "the surface potential is not finite")
    doc.check(ey_above < 0.0 and ey_below > 0.0,
              "the condenser does not point toward the axis: %.6g %.6g" % (ey_above, ey_below))
    doc.check(t_mid > 0.0, "the on-axis ion is reflected in the condenser")
    doc.check(ey_def < 0.0, "under the deflector the field does not point down: %.6g" % ey_def)
    doc.check(q_mid > 0.0 and q_def > 0.0, "a raised electrode does not carry a positive charge")

    doc.value("H", H, 3)
    doc.value("panelFine", PANEL_FINE, 2)
    doc.value("panelCoarse", PANEL_COARSE, 2)
    doc.value("nPanels", fine.n_panels)
    doc.value("contact", CONTACT, 2)
    doc.value("tolFlight", tol, 2, sci=True)
    doc.value("errStep", err_step, 2, sci=True)
    doc.value("errMesh", err_mesh, 2, sci=True)
    doc.value("gamma0", fine.gamma0, 6)
    doc.value("p0", fine.p0, 8)
    doc.value("v0", fine.p0 / (fine.gamma0 * fine.m), 6)
    doc.value("vOverC", fine.p0 / (fine.gamma0 * fine.m * fine.c), 4)
    doc.value("qPart", fine.q, 3, sci=True)
    doc.value("mPart", fine.m, 2)
    doc.value("T0", fine.T0, 2)
    doc.value("cLight", fine.c, 2)
    doc.value("stepReach", fine.c * H, 3)
    doc.value("xLaunch", fine.x0[0], 3)
    doc.value("yLaunch", fine.x0[1], 3)
    doc.value("xDetLo", fine.dmin[0], 3)
    doc.value("yDetLo", fine.dmin[1], 3)
    doc.value("xDetHi", fine.dmax[0], 3)
    doc.value("yDetHi", fine.dmax[1], 3)
    doc.value("xWorld", fine.wmax[0], 3)
    doc.value("yWorld", fine.wmax[1], 3)
    doc.value("axisDeg", fine.axis_deg, 4, sci=False)
    doc.value("halfDeg", fine.half_deg, 4)
    doc.value("xGate0", fine.gmin[0], 3)
    doc.value("yGate0", fine.gmin[1], 3)
    doc.value("xGate1", fine.gmax[0], 3)
    doc.value("yGate1", fine.gmax[1], 3)
    doc.value("opening", opening, 3)
    doc.value("yRimHi", y_rim_hi, 3)
    doc.value("yRimLo", y_rim_lo, 3)
    doc.value("xMid", x_mid, 3)
    doc.value("xIn", outer_x[0], 3)
    doc.value("xOut", outer_x[1], 3)
    doc.value("yAbove", y_above, 3)
    doc.value("yBelow", y_below, 3)
    doc.value("vMid", mids[0][1], 3, sci=True)
    doc.value("vDef", defs[0][1], 3, sci=True)
    volts = [abs(float(v)) for v in lvl["limits"]["supply_voltages"]]
    doc.value("vHalf", max(volts), 3, sci=True)
    doc.value("xDef", x_def, 4)
    doc.value("yDef0", dbox[2], 4)
    doc.value("yDef1", dbox[3], 4)
    doc.value("xDef0", dbox[0], 4)
    doc.value("xDef1", dbox[1], 4)
    doc.value("yUnder", y_under, 4)
    doc.value("xPump", bx, 3)
    doc.value("yPump", by, 3)
    doc.value("qPump", bq, 3, sci=True)
    doc.value("xSteer", sx, 3)
    doc.value("ySteer", sy, 3)
    doc.value("qSteer", sq, 3, sci=True)
    doc.value("Rcharge", fine.R, 2)
    mags = [abs(float(v)) for v in lvl["limits"]["magnitudes"]]
    doc.value("qMag", max(mags), 3, sci=True)
    doc.value("phiLaunch", phi_launch, 4, sci=True)
    doc.value("phiMid", phi_mid, 4, sci=True)
    doc.value("phiIn", phi_in, 4, sci=True)
    doc.value("phiOut", phi_out, 4, sci=True)
    doc.value("eAbove", ey_above, 4, sci=True)
    doc.value("eBelow", ey_below, 4, sci=True)
    doc.value("eAxis", ey_axis, 4, sci=ey_axis != 0.0)
    doc.value("eLaunch", ey_launch, 4, sci=ey_launch != 0.0)
    doc.value("Tmid", t_mid, 4)
    doc.value("slowFrac", fine.q * (phi_mid - phi_launch) / fine.T0, 3)
    doc.value("eDef", ey_def, 4, sci=True)
    doc.value("eyPump", ey_pump, 4, sci=True)
    doc.value("eyOver", ey_over, 4, sci=True)
    doc.value("phiPump", phi_pump, 4, sci=True)
    doc.value("qMid", q_mid, 4, sci=True)
    doc.value("qDef", q_def, 4, sci=True)
    doc.value("qOuter", q_outer, 4, sci=True)
    doc.value("surfRes", surf, 3, sci=True)
    doc.value("solveRes", fine.system_resid, 2, sci=True)
    doc.value("nShots", 3)
    for k, tag in enumerate(("Low", "Axis", "High")):
        doc.value("a" + tag, _angle(lvl["shots"][k]), 6)
    for i, box in enumerate(fine.boxes):
        doc.value("xP%da" % i, box[0], 4)
        doc.value("xP%db" % i, box[1], 4)
        doc.value("yP%da" % i, box[2], 4)
        doc.value("yP%db" % i, box[3], 4)
    direction = lvl["shots"][1]["launch"]["direction"]
    dnorm = math.hypot(direction[0], direction[1])
    doc.value("arrowX", fine.x0[0] + 1.6 * direction[0] / dnorm, 4)
    doc.value("arrowY", fine.x0[1] + 1.6 * direction[1] / dnorm, 4)
    cone_o = (fine.dmin[0], 0.5 * (fine.dmin[1] + fine.dmax[1]))
    cone_reach = 2.2
    doc.value("coneOx", cone_o[0], 4)
    doc.value("coneOy", cone_o[1], 4)
    doc.value("coneXp", cone_o[0] + cone_reach * math.cos(half), 4)
    doc.value("coneYp", cone_o[1] + cone_reach * math.sin(half), 4)
    doc.value("coneXm", cone_o[0] + cone_reach * math.cos(-half), 4)
    doc.value("coneYm", cone_o[1] - cone_reach * math.sin(half), 4)

    tags = ("Low", "Axis", "High")
    rad_frac = []
    for k, tag in enumerate(tags):
        fine.set_direction(k)
        e1 = radiated(fine, H / 2.0, ref[k]["t"])
        e2 = radiated(fine, H, ref[k]["t"])
        rel = abs(e1 - e2) / abs(e1) / 15.0 if e1 else 0.0
        doc.check(rel < TOL_RAD, tag + " radiation quadrature is inside its tolerance")
        doc.value(tag + "Rad", e1, 4, sci=True)
        doc.value(tag + "RadErr", rel, 2, sci=True)
        doc.value(tag + "RadFrac", e1 / fine.T0, 2, sci=True)
        rad_frac.append(e1 / fine.T0)
    doc.value("radFracMax", max(rad_frac), 2, sci=True)

    def publish(prefix, run, step_run):
        doc.text(prefix + "Outcome", run["outcome"])
        doc.text(prefix + "Face", run["face"])
        doc.text(prefix + "Passed", "yes" if run["gate_passed"] else "no")
        doc.value(prefix + "T", run["t"], 6)
        doc.value(prefix + "X", run["x"][0], 6)
        doc.value(prefix + "Y", run["x"][1], 6)
        doc.value(prefix + "Px", run["p"][0], 6)
        doc.value(prefix + "Py", run["p"][1], 6)
        doc.value(prefix + "Ang", math.degrees(math.atan2(run["p"][1], run["p"][0])), 4)
        if run["beta"] is not None:
            doc.value(prefix + "Beta", math.degrees(run["beta"]), 4)
        doc.value(prefix + "Err", step_gap(run, step_run), 2, sci=True)
        doc.value(prefix + "Wdev", run["wdev"], 2, sci=run["wdev"] != 0.0)
        doc.value(prefix + "Steps", run["steps"])
        doc.value(prefix + "Ymax", run["y_max"], 4)
        doc.value(prefix + "Ymin", run["y_min"], 4)
        doc.value(prefix + "Dmin", run["d_min"], 4)
        doc.value(prefix + "Image", run["image_ratio"], 2, sci=True)
        if run["gate"] is not None:
            gate = run["gate"]
            doc.value(prefix + "GT", gate["t"], 4)
            doc.value(prefix + "GX", gate["x"][0], 4)
            doc.value(prefix + "GY", gate["x"][1], 4)
            doc.text(prefix + "GFace", gate["face"])

    for k, tag in enumerate(tags):
        publish(tag, ref[k], ref_step[k])
        publish("Off" + tag, off[k], off_step[k])
    doc.value("refImage", max(ref[k]["image_ratio"] for k in range(3)), 2, sci=True)
    doc.value("offImage", max(off[k]["image_ratio"] for k in range(3)), 2, sci=True)
    doc.value("refDmin", min(ref[k]["d_min"] for k in range(3)), 4)
    doc.data("pathlow", ref[0]["path"])
    doc.data("pathaxis", ref[1]["path"])
    doc.data("pathhigh", ref[2]["path"])
    doc.data("pathofflow", off[0]["path"])
    doc.data("pathoffaxis", off[1]["path"])
    doc.data("pathoffhigh", off[2]["path"])
    doc.table("states",
              [r"$t$", r"$x$", r"$y$", r"$p_x$", r"$p_y$"],
              [(r[0], r[1], r[2], r[3], r[4]) for r in ref[1]["rows"]],
              digits=[4, 4, 4, 4, 4])

    def put(prefix, title, run, placement, flight):
        g = doc.flight(placement, flight=flight)
        ge = g["end"]
        doc.compare(prefix + "Outcome", title + "outcome", run["outcome"], g["outcome"])
        doc.compare(prefix + "T", title + "$t$", run["t"], ge["t"], tol)
        doc.compare(prefix + "X", title + "$x$", run["x"][0], ge["x"][0], tol)
        doc.compare(prefix + "Y", title + "$y$", run["x"][1], ge["x"][1], tol)

    # State rows only. Radiation is reported above and is not a comparison.
    put("Low", "low ", ref[0], 0, 0)
    put("Axis", "axis ", ref[1], 0, 1)
    put("High", "up ", ref[2], 0, 2)
    put("OffLow", "off low ", off[0], 1, 0)
    put("OffAxis", "off axis ", off[1], 1, 1)
    put("OffHigh", "off up ", off[2], 1, 2)
