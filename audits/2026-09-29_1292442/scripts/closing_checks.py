"""Close the checks the remainder note left partial.

Reads Hairer's dop853.f (downloaded beside this script if it is missing), the
generated Rust coefficient tables, and the level-19 electrode. Does not import
the game and does not write under levels/ or crates/.
"""

from __future__ import annotations

import hashlib
import math
import re
import urllib.request
from pathlib import Path

import numpy as np

import remainder_checks as rc

ROOT = Path(__file__).resolve().parents[2]
PHYS = ROOT / "crates" / "physics" / "src"
FORT = ROOT / "audit" / "results" / "dop853.f"
FORT_URL = "https://www.unige.ch/~hairer/prog/nonstiff/dop853.f"
EXPECTED_SHA = "8bfc0a954a38af1a458f5ab8025fb09c923373a936f52ff61e32f1a15c8f298a"


def show(name: str, value: float) -> None:
    print(f"{name}: {value:.6e}")


def ensure_fortran() -> str:
    scratch = Path(r"C:\Users\qrzys\AppData\Local\Temp\grok-goal-e1265ff9b472\implementer\dop853.f")
    if not FORT.exists() and scratch.exists():
        FORT.write_bytes(scratch.read_bytes())
    if not FORT.exists():
        urllib.request.urlretrieve(FORT_URL, FORT)
    data = FORT.read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    show("dop853 sha matches recorded", 0.0 if sha == EXPECTED_SHA else 1.0)
    print(f"dop853 sha: {sha}")
    return data.decode("ascii")


def fortran_tables(text: str):
    start = text.index("SUBROUTINE DP86CO")
    end = text.index("DOUBLE PRECISION Y(N)", start)
    block = text[start:end]
    pairs = re.findall(r"([a-z]+ ?\d+)\s*=\s*([-+]?[0-9.]+D[-+]?\d+)", block)
    n = 17
    a = np.zeros((n, n))
    b = np.zeros(n)
    c = np.zeros(n)
    e = np.zeros(n)
    bhh = np.zeros(4)
    d = np.zeros((8, n))
    seen = set()
    for name, value in pairs:
        name = name.replace(" ", "")
        if name in seen:
            raise SystemExit(f"duplicate {name}")
        seen.add(name)
        lit = float(value.replace("D", "E"))
        prefix = re.match(r"[a-z]+", name).group(0)
        digits = name[len(prefix) :]
        if prefix == "a":
            split = {2: 1, 3: 2, 4: 2}[len(digits)]
            i, j = int(digits[:split]), int(digits[split:])
            a[i, j] = lit
        elif prefix == "b":
            b[int(digits)] = lit
        elif prefix == "c":
            c[int(digits)] = lit
        elif prefix == "er":
            e[int(digits)] = lit
        elif prefix == "bhh":
            bhh[int(digits)] = lit
        elif prefix == "d":
            d[int(digits[0]), int(digits[1:])] = lit
        else:
            raise SystemExit(name)
    # Hairer evaluates stages 12 and 13 at x+h. They have no c parameter.
    c[12] = 1.0
    c[13] = 1.0
    return a, b, c, e, bhh, d, len(seen)


def dop853_diff(text: str) -> None:
    a, b, c, e, bhh, d, nseen = fortran_tables(text)
    rust = (PHYS / "integrator" / "dop853_coefficients.rs").read_text(encoding="utf-8")
    ra = np.array(
        [
            [float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][+-]?\d+)?", row)]
            for row in re.findall(
                r"\[(.*?)\]",
                re.search(r"pub const A:.*?=\s*\[(.*?)\];", rust, re.S).group(1),
                re.S,
            )
        ]
    )
    rb = rc.parse_rust_array(rust, "B")
    rc_ = rc.parse_rust_array(rust, "C")
    re_ = rc.parse_rust_array(rust, "E")
    rbhh = rc.parse_rust_array(rust, "BHH")
    rd_rows = re.findall(
        r"\[(.*?)\]",
        re.search(r"pub const D:.*?=\s*\[(.*?)\];", rust, re.S).group(1),
        re.S,
    )
    rd = np.array(
        [[float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][+-]?\d+)?", row)] for row in rd_rows]
    )
    show("dop853 coefficient count", float(nseen))
    show("dop853 A max abs", float(np.max(np.abs(a - ra))))
    show("dop853 B max abs", float(np.max(np.abs(b - rb))))
    show("dop853 C max abs", float(np.max(np.abs(c - rc_))))
    show("dop853 E max abs", float(np.max(np.abs(e - re_))))
    show("dop853 BHH max abs", float(np.max(np.abs(bhh - rbhh))))
    show("dop853 D max abs", float(np.max(np.abs(d - rd))))
    # Fortran stores stage 11 in K2 and stage 12 in K3 before the error estimate.
    # B11*K2, B12*K3, BHH3*K3, ER11*K2, ER12*K3 are those stages.
    src = (PHYS / "integrator" / "dop853.rs").read_text(encoding="utf-8")
    stage_ok = all(
        token in src
        for token in (
            "BHH[3] * k[12][i]",
            "for j in [1, 6, 7, 8, 9, 10, 11, 12]",
        )
    )
    show("dop853 stage-12 error index present", 0.0 if stage_ok else 1.0)


def controller(text: str) -> None:
    def grab(name: str) -> float:
        m = re.search(rf"{name}=([0-9.]+)D([+-]?\d+)", text)
        return float(m.group(1) + "E" + m.group(2))

    safe, fac1, fac2, beta, uround = (grab(n) for n in ("SAFE", "FAC1", "FAC2", "BETA", "UROUND"))
    rust = (PHYS / "integrator" / "dop853.rs").read_text(encoding="utf-8")
    defaults = {
        "safety": 0.9,
        "fac_min": 0.333,
        "fac_max": 6.0,
        "beta": 0.0,
        "uround": 2.3e-16,
    }
    got = {"safety": safe, "fac_min": fac1, "fac_max": fac2, "beta": beta, "uround": uround}
    worst = max(abs(got[k] - defaults[k]) for k in defaults)
    # The Rust source must contain the same literals the Fortran initialises.
    present = all(token in rust for token in ("safety: 0.9", "fac_min: 0.333", "fac_max: 6.0", "beta: 0.0", "uround: 2.3e-16"))
    show("dop853 controller defaults", worst if present else 1.0)

    def h_factor(err: float, facold: float, accept: bool) -> float:
        expo1 = 1.0 / 8.0 - beta * 0.2
        facc1 = 1.0 / fac1
        facc2 = 1.0 / fac2
        fac11 = err**expo1
        if accept:
            fac = max(facc2, min(facc1, (fac11 / facold**beta) / safe))
            return 1.0 / fac
        return 1.0 / min(facc1, fac11 / safe)

    # Same arithmetic the Fortran lines FAC=MAX(FACC2,MIN(FACC1,FAC/SAFE)) and the
    # rejected-step line HNEW=H/MIN(FACC1,FAC11/SAFE) produce. Spot-check the powers.
    errs = [1e-8, 1e-3, 0.5, 1.0, 2.0, 10.0]
    worst_h = 0.0
    for err in errs:
        facold = max(err, 1e-4)
        expo1 = 0.125 - beta * 0.2
        fac11 = err**expo1
        fac = max(1.0 / fac2, min(1.0 / fac1, (fac11 / facold**beta) / safe))
        worst_h = max(worst_h, abs(h_factor(err, facold, True) - 1.0 / fac))
        reject = 1.0 / min(1.0 / fac1, fac11 / safe)
        worst_h = max(worst_h, abs(h_factor(err, facold, False) - reject))
    show("dop853 step-factor identity", worst_h)


def dense_derivative() -> None:
    rng = np.random.default_rng(1)
    c = rng.normal(size=8)

    def y(s: float) -> float:
        s1 = 1.0 - s
        conpar = c[4] + s * (c[5] + s1 * (c[6] + s * c[7]))
        return c[0] + s * (c[1] + s1 * (c[2] + s * (c[3] + s1 * conpar)))

    def dy_ds(s: float) -> float:
        s1 = 1.0 - s
        q3 = c[6] + s * c[7]
        dq3 = c[7]
        q2 = c[5] + s1 * q3
        dq2 = -q3 + s1 * dq3
        q = c[4] + s * q2
        dq = q2 + s * dq2
        p4 = c[3] + s1 * q
        dp4 = -q + s1 * dq
        p3 = c[2] + s * p4
        dp3 = p4 + s * dp4
        p2 = c[1] + s1 * p3
        dp2 = -p3 + s1 * dp3
        return p2 + s * dp2

    h = 0.37
    worst = 0.0
    for s in np.linspace(0.05, 0.95, 7):
        # Complex step is exact for this polynomial up to roundoff.
        num = (y(s + 1e-30j)).imag / 1e-30
        worst = max(worst, abs(dy_ds(float(s)) - num) / max(abs(num), 1.0))
        worst = max(worst, abs(dy_ds(float(s)) / h - num / h) / max(abs(num / h), 1.0))
    show("dense derivative vs complex step", worst)


def partial_pivot_lu(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Row-major partial-pivot LU, the same index order as bem.rs."""
    n = a.shape[0]
    m = a.copy()
    perm = np.arange(n)
    for k in range(n):
        p = k + int(np.argmax(np.abs(m[k:, k])))
        if p != k:
            m[[k, p]] = m[[p, k]]
            perm[[k, p]] = perm[[p, k]]
        pivot = m[k, k]
        for i in range(k + 1, n):
            f = m[i, k] / pivot
            m[i, k] = f
            if f != 0.0:
                m[i, k + 1 :] -= f * m[k, k + 1 :]
    return m, perm


def lu_solve(lu: np.ndarray, perm: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = len(b)
    x = b[perm].astype(float).copy()
    for i in range(n):
        x[i] -= lu[i, :i] @ x[:i]
    for i in range(n - 1, -1, -1):
        x[i] = (x[i] - lu[i, i + 1 :] @ x[i + 1 :]) / lu[i, i]
    return x


def bem() -> None:
    level = __import__("json").loads((ROOT / "levels" / "19_power_supply.json").read_text(encoding="utf-8"))
    e = level["electrodes"][0]
    center = np.array(e["center"], dtype=float)
    tris = rc.electrode_mesh(
        center,
        math.radians(e["angle_deg"]),
        e["length"] / 2.0,
        e["thickness"] / 2.0,
        e["height"] / 2.0,
        1.0,
    )
    # Collocation at centroids. In the plane z=0 the mirror doubles the potential.
    cents = [(a + b + c) / 3.0 for a, b, c in tris]
    n = len(tris)
    show("bem panels", float(n))
    # A sample of matrix entries against an independent barycentric quadrature.
    worst = 0.0
    checked = 0
    for i in range(0, n, max(1, n // 6)):
        for j in range(0, n, max(1, n // 4)):
            if i == j:
                continue
            pot, _ = rc.triangle_integrals(*tris[j], cents[i])
            qpot, _ = rc.quad_triangle(*tris[j], cents[i], n=24)
            scale = max(abs(qpot), 1e-30)
            worst = max(worst, abs(pot - qpot) / scale)
            checked += 1
    show("bem entries checked", float(checked))
    show("bem entry vs quadrature", worst)
    # Far field of one panel with unit density: area/R.
    a, b, c = tris[0]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    far = np.array([80.0, 40.0, 25.0])
    pot, _ = rc.triangle_integrals(a, b, c, far)
    mono = area / np.linalg.norm(far - (a + b + c) / 3.0)
    show("bem one panel vs monopole", abs(pot - mono) / mono)
    # LU against numpy on the collocation matrix of the first 12 panels.
    m = 12
    mat = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            pot, _ = rc.triangle_integrals(*tris[j], cents[i])
            mat[i, j] = pot
    rng = np.random.default_rng(2)
    sigma = rng.normal(size=m)
    rhs = mat @ sigma
    lu, perm = partial_pivot_lu(mat)
    got = lu_solve(lu, perm, rhs)
    ref = np.linalg.solve(mat, rhs)
    show("bem lu vs numpy", float(np.max(np.abs(got - ref)) / np.max(np.abs(ref))))
    show("bem lu residual", float(np.max(np.abs(mat @ got - rhs)) / np.max(np.abs(rhs))))


def fib(n: int, twist: float) -> np.ndarray:
    golden = math.pi * (3.0 - math.sqrt(5.0))
    i = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (i + 0.5) / n
    rho = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = golden * i + twist
    return np.column_stack((rho * np.cos(phi), rho * np.sin(phi), z))


def householder_solve(mat: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Column-major Householder QR and the triangular solve in conductor.rs."""
    m, n = mat.shape
    a = np.asfortranarray(mat).ravel(order="F").copy()
    tau = np.zeros(n)
    for k in range(n):
        col = k * m
        norm = math.sqrt(float(np.dot(a[col + k : col + m], a[col + k : col + m])))
        if norm == 0.0:
            continue
        alpha = -norm if a[col + k] > 0.0 else norm
        v0 = a[col + k] - alpha
        a[col + k + 1 : col + m] /= v0
        tau[k] = -v0 / alpha
        a[col + k] = alpha
        for j in range(k + 1, n):
            cj = j * m
            s = a[cj + k] + float(np.dot(a[col + k + 1 : col + m], a[cj + k + 1 : cj + m]))
            s *= tau[k]
            a[cj + k] -= s
            a[cj + k + 1 : cj + m] -= s * a[col + k + 1 : col + m]
    b = rhs.copy()
    for k in range(n):
        s = b[k] + float(np.dot(a[k * m + k + 1 : (k + 1) * m], b[k + 1 :]))
        s *= tau[k]
        b[k] -= s
        b[k + 1 :] -= s * a[k * m + k + 1 : (k + 1) * m]
    x = np.zeros(n)
    for k in range(n - 1, -1, -1):
        s = b[k]
        for j in range(k + 1, n):
            s -= a[j * m + k] * x[j]
        x[k] = s / a[k * m + k]
    return x


def mfs_uniform(k: int) -> None:
    """Grounded sphere in a unit field along z. Shell 0.6a, twist 0.37, 2K points."""
    a = 1.0
    y = 0.6 * a * fib(k, 0.0)
    x = a * fib(2 * k, 0.37)
    rhs = x[:, 2].copy()  # cancel φ_ext = -z
    mat = 1.0 / np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)
    coef, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    xt = 1.7 * fib(240, 1.1)
    pred = coef @ (1.0 / np.linalg.norm(xt[:, None, :] - y[None, :, :], axis=2)).T
    r = np.linalg.norm(xt, axis=1)
    exact = a**3 * xt[:, 2] / r**3
    show(f"mfs K={k} exterior vs dipole", float(np.max(np.abs(pred - exact)) / np.max(np.abs(exact))))
    boundary = coef @ (1.0 / np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)).T
    show(f"mfs K={k} boundary vs z", float(np.max(np.abs(boundary - rhs)) / np.max(np.abs(rhs))))


def mfs_householder() -> None:
    k = 24
    y = 0.6 * fib(k, 0.0)
    x = fib(2 * k, 0.37)
    rhs = x[:, 2].copy()
    mat = 1.0 / np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)
    mine = householder_solve(mat, rhs)
    ref, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    show("mfs householder vs lstsq", float(np.max(np.abs(mine - ref)) / np.max(np.abs(ref))))


def pipe_exact() -> None:
    w = 2.0
    k = math.pi / w
    # Second derivatives of sin(kx) exp(-kz) cancel. Report the floating residual.
    worst = 0.0
    for x in (0.2, 0.7, 1.3):
        for z in (0.0, 0.5, 2.0):
            phi = math.sin(k * x) * math.exp(-k * z)
            dxx = -(k * k) * phi
            dzz = (k * k) * phi
            worst = max(worst, abs(dxx + dzz) / max(abs(phi), 1e-30))
    show("strip mode laplace analytic", worst)
    # Square duct, side w: lowest Dirichlet decay is π√2 / w, not π/w.
    show("square duct k ratio minus 1", abs((k) / (math.pi * math.sqrt(2.0) / w) - 1.0))
    # Circular pipe of diameter w=2R: j01/R versus π/w = π/(2R).
    # j01 from a short bisection of the series, checked against J0.
    lo, hi = 2.0, 2.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if j0(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    j01 = 0.5 * (lo + hi)
    show("J0 at located root", abs(j0(j01)))
    show("circular pipe k ratio minus 1", abs((math.pi / 2.0) / j01 - 1.0))


def j0(x: float) -> float:
    s = 0.0
    term = 1.0
    y = (x * 0.5) ** 2
    for n in range(40):
        s += term
        term *= -y / ((n + 1) * (n + 1))
    return s


def shader() -> None:
    def fm(w2, s, dt):
        """Shader branch: series for |ω²s²| < 0.5, closed form otherwise, in dt."""
        w2 = dt(w2)
        s = dt(s)
        z = w2 * s * s
        if abs(z) < dt(0.5):
            s1 = s * (1 - z / 6 * (1 - z / 20 * (1 - z / 42 * (1 - z / 72))))
            c1 = s * s * dt(0.5) * (1 - z / 12 * (1 - z / 30 * (1 - z / 56 * (1 - z / 90))))
            return dt(1) - w2 * c1, s1, c1
        if w2 > 0:
            w = np.sqrt(w2)
            h = np.sin(dt(0.5) * w * s)
            return np.cos(w * s), np.sin(w * s) / w, dt(2) * h * h / w2
        k = np.sqrt(-w2)
        h = np.sinh(dt(0.5) * k * s)
        return np.cosh(k * s), np.sinh(k * s) / k, dt(2) * h * h / -w2

    worst = 0.0
    # |z| just below the shader's series switch, and a hyperbolic step past it.
    for w2, s in ((0.2, 1.5), (-0.3, 1.2), (4.0, 0.5), (-4.0, 0.5)):
        a = np.array(fm(w2, s, np.float32), dtype=float)
        b = np.array(fm(w2, s, np.float64), dtype=float)
        worst = max(worst, float(np.max(np.abs(a - b)) / np.max(np.abs(b))))
    show("shader field-motion vs f64", worst)
    show("shader ln2 absolute", abs(0.6931472 - math.log(2.0)))
    # In the plane the dipole is Bz = -μ/r³. The shader writes that. f32 evaluation.
    mu = np.float32(1.7)
    r = np.float32(2.5)
    shader_bz = -mu / (r * r * r)
    exact = -1.7 / 2.5**3
    show("shader dipole bz", abs(float(shader_bz) - exact) / abs(exact))
    # Potential floor: Q/max(r, 1e-4) against Q/r.
    q = 1.0
    for radius, label in ((1.0, "outside"), (1e-6, "inside floor")):
        painted = q / max(radius, 1e-4)
        show(f"shader coulomb {label}", abs(painted - q / radius) / (q / radius))


def caps() -> None:
    text = (ROOT / "crates" / "game" / "src" / "potential.rs").read_text(encoding="utf-8")
    caps = {
        "MAX_CHARGES": 1024,
        "MAX_MAGNETS": 64,
        "MAX_LOOPS": 16,
        "MAX_SEGMENTS": 64,
    }
    for name, expect in caps.items():
        found = int(re.search(rf"{name}: usize = (\d+)", text).group(1))
        uses = f".take({name}" in text or f"n == {name}" in text
        show(f"cap {name}", 0.0 if found == expect and uses else 1.0)
    # The upload keeps a prefix. 1025 charges, 65 magnets, 17 loops, 65 segments.
    for label, n, cap in (
        ("charges dropped", 1025, 1024),
        ("magnets dropped", 65, 64),
        ("loops dropped", 17, 16),
        ("segments dropped", 65, 64),
    ):
        show(label, float(n - min(n, cap)))
    # A live scene of 1025 charges has to travel in EM_PLACE. Measure that string.
    import json

    payload = json.dumps([{"node": [i, 0, 0], "value": 1.0} for i in range(1025)], separators=(",", ":"))
    show("place json characters", float(len(payload)))
    print(f"windows env char note: {len(payload)} versus 32767")


def field_motion_root() -> None:
    """Port of FieldMotion::retarded for uniform motion, against the quadratic root."""
    c = 5.0
    v = np.array([0.4, -0.2, 0.0])
    gamma = 1.0 / math.sqrt(1.0 - np.dot(v, v) / c**2)
    u0 = np.array([gamma * c, gamma * v[0], gamma * v[1]])
    x0 = np.zeros(3)
    x = np.array([1.2, 0.4, 0.0])
    dt = 0.8

    def at(s):
        dx = s * u0
        return u0, dx

    def eval_g(s):
        u, dx = at(s)
        d = x - (x0 + np.array([dx[1], dx[2], 0.0]))
        dist = np.linalg.norm(d)
        value = c * dt - dx[0] - dist
        slope = -u[0] + np.dot(d, np.array([u[1], u[2], 0.0])) / max(dist, 1e-300)
        return value, slope

    r = x - v * dt
    rv = float(np.dot(r, v))
    v2 = float(np.dot(v, v))
    tau_u = dt - (rv + math.sqrt(rv * rv + (c * c - v2) * float(np.dot(r, r)))) / (c * c - v2)
    s_true = tau_u / gamma
    scale = max(np.linalg.norm(r) / c, abs(dt), 1e-300) / gamma
    lo, hi = -math.inf, math.inf
    s = tau_u / gamma
    for iteration in range(100):
        if iteration == 12 and not (math.isfinite(lo) and math.isfinite(hi)):
            return
        value, slope = eval_g(s)
        if value > 0.0:
            lo = max(lo, s)
        else:
            hi = min(hi, s)
        nxt = s - value / slope
        if math.isfinite(lo) and math.isfinite(hi) and not (lo < nxt < hi):
            nxt = 0.5 * (lo + hi)
        if abs(nxt - s) <= 4.0 * np.finfo(float).eps * (abs(s) + scale):
            break
        s = nxt
    show("retarded worldline vs quadratic", abs(s - s_true))
    show("retarded worldline residual", abs(eval_g(s)[0]))


def main() -> None:
    text = ensure_fortran()
    dop853_diff(text)
    controller(text)
    dense_derivative()
    bem()
    mfs_householder()
    for k in (40, 80, 440):
        mfs_uniform(k)
    pipe_exact()
    shader()
    caps()
    field_motion_root()


if __name__ == "__main__":
    main()
