"""Recompute the audit's numeric checks from level files and from the formulas in the source.

Does not import the game. The cloud and stray-field sections read ``levels/*.json`` and apply
the same expressions as ``crates/physics/src/field.rs`` (cloud interior) and
``crates/game/src/potential.wgsl`` (point-charge potential). The wave section transcribes
``PlaneWave::fields`` in ``crates/physics/src/external.rs``. The elliptic section transcribes
``elliptic_ke`` and ``loop_bracket`` in ``crates/physics/src/magnetic.rs``. Later sections
transcribe the detector box (``Level::box_region`` plus ``Aabb::signed_distance``), ``β`` from
the kinetic energy, ``beam::tapered``, and the Faraday-cup force after absorption is booked.

Run from anywhere:

    py -3 audit/scripts/audit_checks.py

Writes ``audit/results/checks.txt``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEVELS = ROOT / "levels"
OUT = ROOT / "audit" / "results" / "checks.txt"


def phi_cloud(q: float, radius: float, dist: float) -> float:
    """Potential of a uniform sphere, k = 1. Inside: Q (3R^2 - d^2) / (2 R^3)."""
    if dist < radius:
        return q * (3.0 * radius * radius - dist * dist) / (2.0 * radius**3)
    return q / dist


def phi_shader(q: float, dist: float) -> float:
    """Potential used by potential.wgsl for every Coulomb source, including clouds."""
    return q / max(math.sqrt(dist * dist), 1e-4)


def elliptic_ke(m: float, m1: float) -> tuple[float, float]:
    a = 1.0
    b = math.sqrt(m1)
    c = math.sqrt(m)
    acc = 0.5 * c * c
    pow2 = 0.5
    for _ in range(40):
        if abs(c) <= math.ulp(1.0) * a:
            break
        an = 0.5 * (a + b)
        bn = math.sqrt(a * b)
        c = 0.5 * (a - b)
        a, b = an, bn
        pow2 *= 2.0
        acc += pow2 * c * c
    k = math.pi / (2.0 * a)
    return k, k * (1.0 - acc)


def loop_bracket(m: float, m1: float) -> float:
    if m < 0.05:
        c = 1.0
        e_prev = 1.0
        k_prev = 1.0
        total = 0.0
        mp = 1.0
        for n in range(1, 17):
            nf = float(n)
            c *= ((2.0 * nf - 1.0) / (2.0 * nf)) ** 2
            e_n = c / (1.0 - 2.0 * nf)
            k_n = c
            mp *= m
            if n >= 2:
                total += mp * (e_n - 0.5 * e_prev - k_n + k_prev)
            e_prev, k_prev = e_n, k_n
        return 0.5 * math.pi * total
    k, e = elliptic_ke(m, m1)
    return (1.0 - 0.5 * m) * e - m1 * k


def plane_wave_fields(amplitude: float, angle: float, omega: float, phase: float, c: float, x, t: float):
    """Transcription of PlaneWave::in_plane + fields. c = inf is represented by math.inf."""
    direction = (math.cos(angle), math.sin(angle), 0.0)
    polarization = (-math.sin(angle), math.cos(angle), 0.0)
    if math.isfinite(c):
        retard = (direction[0] * x[0] + direction[1] * x[1] + direction[2] * x[2]) / c
    else:
        retard = 0.0
    ph = omega * (t - retard) + phase
    scale = amplitude * math.cos(ph)
    e = tuple(p * scale for p in polarization)
    if math.isfinite(c):
        # k̂ × E / c
        b = (
            (direction[1] * e[2] - direction[2] * e[1]) / c,
            (direction[2] * e[0] - direction[0] * e[2]) / c,
            (direction[0] * e[1] - direction[1] * e[0]) / c,
        )
    else:
        b = (0.0, 0.0, 0.0)
    return e, b


def momentum(mass: float, kinetic: float, c: float) -> float:
    """|p| from kinetic energy. c = inf is the Newtonian p = sqrt(2 m T)."""
    if math.isfinite(c):
        return math.sqrt(kinetic * kinetic / (c * c) + 2.0 * mass * kinetic)
    return math.sqrt(2.0 * mass * kinetic)


def load_levels() -> list[tuple[str, dict]]:
    files = sorted(LEVELS.glob("[0-9]*.json"))
    out = []
    for path in files:
        out.append((path.name, json.loads(path.read_text(encoding="utf-8"))))
    return out


def section_clouds(levels, lines: list[str]) -> None:
    lines.append("== Charge clouds vs the potential-map formula ==")
    lines.append(
        "The trajectory uses the uniform-sphere potential. The potential map evaluates Q/r "
        "for every Coulomb source (potential.wgsl), and paints a solid disk of charge_radius "
        "at the centre."
    )
    any_cloud = False
    for name, level in levels:
        clouds = level.get("clouds") or []
        if not clouds:
            continue
        any_cloud = True
        subdiv = float(level["grid"]["subdivision"])
        charge_radius = float(level["physics"]["charge_radius"])
        for cloud in clouds:
            center = [c / subdiv for c in cloud["center"]]
            radius = float(cloud["radius"])
            q = float(cloud["charge"])
            lines.append(f"\n{name}: cloud Q={q:g} R={radius:g} at {center}, charge_radius={charge_radius:g}")
            for frac, label in [(0.0, "centre"), (0.25, "R/4"), (0.5, "R/2"), (0.75, "3R/4"), (1.0, "surface")]:
                dist = frac * radius
                true = phi_cloud(q, radius, dist)
                shader = phi_shader(q, dist if dist > 0.0 else 0.0)
                ratio = shader / true if true else float("nan")
                lines.append(
                    f"  {label:8} d={dist:6.3f}  phi_true={true:.8g}  phi_map={shader:.8g}  map/true={ratio:.4f}"
                )
            # Turning line of each shot launched inside this cloud, mixing the CPU energy
            # (true potential at the launch) with the map's 1/r potential. That is what
            # potential.rs + potential.wgsl do: the limit is sampled with FieldSolver::sample,
            # the colour with Q/r.
            for shot in level.get("shots", []):
                node = shot["launch"]["node"]
                launch = [node[0] / subdiv, node[1] / subdiv, node[2] / subdiv]
                dist = math.dist(launch, center)
                if dist >= radius:
                    continue
                particle_q = float(shot["particle"]["charge"])
                kinetic = float(shot["launch"]["kinetic_energy"])
                phi0 = phi_cloud(q, radius, dist)
                energy = kinetic + particle_q * phi0
                # Allowed where particle_q * phi <= energy. For particle_q < 0 that is phi >= energy/particle_q.
                if particle_q == 0.0:
                    continue
                phi_turn = energy / particle_q

                def radius_at(phi_of):
                    # Binary search the radius where phi_of(d) = phi_turn, inside the cloud.
                    lo, hi = 0.0, radius
                    for _ in range(60):
                        mid = 0.5 * (lo + hi)
                        value = phi_of(mid)
                        # Positive cloud: phi falls with radius. particle may be negative.
                        if (q > 0 and value > phi_turn) or (q < 0 and value < phi_turn):
                            lo = mid
                        else:
                            hi = mid
                    return 0.5 * (lo + hi)

                r_true = radius_at(lambda d: phi_cloud(q, radius, d))
                r_map = radius_at(lambda d: phi_shader(q, d))
                lines.append(
                    f"  shot launch d={dist:.4f} T={kinetic:g} q={particle_q:g} "
                    f"total energy={energy:.8g}"
                )
                lines.append(
                    f"  turning radius  true={r_true:.4f}  map={r_map:.4f}  "
                    f"(map paints a solid disk of radius {charge_radius:g} at the centre)"
                )
                if kinetic > 0.0:
                    # Colour at the launch point: u_a uses the CPU potential, the shader
                    # uses Q/r, so u = (q/T) (phi_map - phi_true) instead of 0.
                    phi_m = phi_shader(q, dist)
                    u = particle_q * (phi_m - phi0) / kinetic
                    lines.append(
                        f"  colour at launch (q/T)(phi_map-phi_true) = {u:.6g} T0 "
                        f"(0 if the map used the trajectory potential)"
                    )
                    edge = phi_shader(q, charge_radius)
                    edge_true = phi_cloud(q, radius, charge_radius)
                    u_edge = particle_q * (edge - edge_true) / kinetic
                    lines.append(
                        f"  at the painted disk edge d={charge_radius:g}: "
                        f"phi_map={edge:.8g} phi_true={edge_true:.8g} "
                        f"colour offset={u_edge:.6g} T0"
                    )
    if not any_cloud:
        lines.append("No charge clouds in levels/*.json.")


def section_stray(levels, lines: list[str]) -> None:
    lines.append("\n== Stray fields absent from the potential and magnetic maps ==")
    lines.append(
        "potential.wgsl sums fixed charges, dipoles, circular coils and polygon segments. "
        "Uniform stray E and B_z are stored on the flight (Level::scenario_with) and are "
        "sampled by the Total field view, and they are not uploaded by potential::params. "
        "The potential and magnetic maps are built from Level::display_scenario, which always "
        "passes disturbance index 0, and the map cache key does not include the selected disturbance."
    )
    for name, level in levels:
        disturbances = level.get("disturbances") or []
        interesting = []
        for i, d in enumerate(disturbances):
            e = d.get("e") or [0.0, 0.0]
            bz = float(d.get("bz") or 0.0)
            waves = d.get("waves") or []
            if e != [0.0, 0.0] or bz != 0.0 or waves:
                interesting.append((i, d.get("name", ""), e, bz, len(waves)))
        if not interesting:
            continue
        lines.append(f"\n{name}  ({len(disturbances)} realization(s); map uses index 0)")
        nx = int(level["grid"]["nx"])
        ny = int(level["grid"]["ny"])
        for i, title, e, bz, n_waves in interesting:
            lines.append(
                f"  [{i}] {title!r}  E=({e[0]:g}, {e[1]:g})  Bz={bz:g}  waves={n_waves}"
            )
            if e != [0.0, 0.0]:
                # Potential of the uniform field across the arena, -E·x.
                span = abs(e[0]) * nx + abs(e[1]) * ny
                lines.append(f"      |Δφ| of -E·x across the grid ≈ {span:g}")
            if bz != 0.0 and i != 0:
                lines.append(
                    "      magnetic map stays on realization 0 while this flight is selected"
                )
            if i == 0 and (e != [0.0, 0.0] or bz != 0.0):
                # display_scenario uses disturbance 0, so this flight's uniform field is
                # the one the maps would show if they uploaded it. They upload neither.
                subdiv = float(level["grid"]["subdivision"])
                c_raw = level.get("physics", {}).get("c")
                c = math.inf if c_raw is None else float(c_raw)
                for s_i, shot in enumerate(level.get("shots", [])):
                    particle = shot["particle"]
                    q_p = float(particle["charge"])
                    mass = float(particle["mass"])
                    kinetic = float(shot["launch"]["kinetic_energy"])
                    node = shot["launch"]["node"]
                    x = node[0] / subdiv
                    y = node[1] / subdiv
                    phi0 = -(float(e[0]) * x + float(e[1]) * y)
                    lines.append(
                        f"      shot {s_i + 1} at ({x:g},{y:g}) q={q_p:g} T={kinetic:g}"
                    )
                    if kinetic > 0.0 and q_p != 0.0 and e != [0.0, 0.0]:
                        u = -(q_p * phi0) / kinetic
                        energy = kinetic + q_p * phi0
                        g = -energy / kinetic
                        lines.append(
                            f"        empty potential map: flat colour u={u:.6g} T0, "
                            f"forbidden excess g={g:.6g} (true u at launch is 0)"
                        )
                        e_vec = math.hypot(float(e[0]), float(e[1]))
                        # q (phi - phi0) = T on the energy boundary. phi - phi0 = -E·Δx.
                        climb = kinetic / abs(q_p * e_vec) if e_vec else float("nan")
                        lines.append(
                            f"        true energy boundary is {climb:.6g} cells against E "
                            f"(|q E|={abs(q_p * e_vec):.6g} per cell)"
                        )
                    if bz != 0.0 and q_p != 0.0 and kinetic > 0.0:
                        p = momentum(mass, kinetic, c)
                        b_ref = p / (abs(q_p) * 5.0)
                        lines.append(
                            f"        uniform Bz/b_ref={bz / b_ref:.6g} "
                            f"(b_ref={b_ref:.6g} for a 5-cell gyroradius; the magnetic map adds 0)"
                        )


def section_inventory(levels, lines: list[str]) -> None:
    lines.append("\n== Moments, free particles, static waves, metal with radiation reaction ==")
    moment_shots = []
    free = []
    static_waves = []
    rr_metal = []
    for name, level in levels:
        if any(s.get("particle", {}).get("moment", 0.0) for s in level.get("shots", [])):
            moment_shots.append(name)
        if level.get("free_particles"):
            free.append(name)
        for d in level.get("disturbances") or []:
            for w in d.get("waves") or []:
                if float(w.get("omega") or 0.0) == 0.0:
                    static_waves.append(name)
        if level.get("physics", {}).get("radiation_reaction") and (
            level.get("conductors") or level.get("electrodes")
        ):
            rr_metal.append(name)
    lines.append(f"  shots with a magnetic moment: {moment_shots or 'none'}")
    lines.append(f"  levels with free_particles: {free or 'none'}")
    lines.append(f"  levels with an ω=0 wave: {static_waves or 'none'}")
    lines.append(f"  radiation reaction together with metal: {rr_metal or 'none'}")


def section_wave(lines: list[str]) -> None:
    lines.append("\n== Static plane-wave term, finite c ==")
    lines.append(
        "PlaneWave::fields sets B = k̂ × E / c whenever c is finite, including ω = 0. "
        "External::sample's comment says that term has B = 0. W4 checks energy, and a "
        "uniform B does no work, so that test still passes."
    )
    amplitude, angle, phase, c = 0.03, 0.4, 0.8, 4.0
    for omega in (0.0, 1.3):
        e, b = plane_wave_fields(amplitude, angle, omega, phase, c, (0.0, 0.0, 0.0), 0.0)
        bl = math.sqrt(sum(v * v for v in b))
        el = math.sqrt(sum(v * v for v in e))
        lines.append(
            f"  ω={omega:g} c={c:g}: |E|={el:.6g} |B|={bl:.6g} |B|/( |E|/c )={bl / (el / c):.6g}"
        )
    e_inf, b_inf = plane_wave_fields(amplitude, angle, 0.0, phase, math.inf, (0.0, 0.0, 0.0), 0.0)
    lines.append(
        f"  ω=0 c=∞: |E|={math.sqrt(sum(v*v for v in e_inf)):.6g} "
        f"|B|={math.sqrt(sum(v*v for v in b_inf)):.6g}"
    )


def section_static_dipole(lines: list[str]) -> None:
    lines.append("\n== Static electric dipole: field versus −∇(n·p / r²) ==")
    lines.append(
        "antenna.rs::fields returns the Jackson E below and phi = 0 for every ω, including 0. "
        "LevelField::is_static is true when every antenna has ω = 0, and the trajectory energy "
        "diagnostic then uses q φ."
    )
    p = (0.8, -0.3)
    x = (1.7, 0.4)
    h = 1e-6

    def phi(pt):
        r2 = pt[0] * pt[0] + pt[1] * pt[1]
        r = math.sqrt(r2)
        return (p[0] * pt[0] + p[1] * pt[1]) / (r * r2)

    def dipole_e(pt):
        r = math.hypot(pt[0], pt[1])
        n = (pt[0] / r, pt[1] / r)
        nd = n[0] * p[0] + n[1] * p[1]
        scale = 1.0 / (r * r * r)
        return ((3.0 * n[0] * nd - p[0]) * scale, (3.0 * n[1] * nd - p[1]) * scale)

    e = dipole_e(x)
    e_phi = (
        -(phi((x[0] + h, x[1])) - phi((x[0] - h, x[1]))) / (2.0 * h),
        -(phi((x[0], x[1] + h)) - phi((x[0], x[1] - h))) / (2.0 * h),
    )
    err = math.hypot(e[0] - e_phi[0], e[1] - e_phi[1])
    lines.append(
        f"  at {x}, p={p}: |E|={math.hypot(*e):.8g}  φ=n·p/r²={phi(x):.8g}  "
        f"|E − (−∇φ)|={err:.3e}  (the code stores φ = 0)"
    )


def node_position(grid: dict, node: list) -> list[float]:
    s = float(grid.get("subdivision") or 1)
    z = float(node[2]) if len(node) > 2 else 0.0
    return [float(node[0]) / s, float(node[1]) / s, z / s]


def detector_box(grid: dict, det: dict) -> tuple[list[float], list[float]]:
    """Level::box_region: node box, and in 2D a slab z ∈ [-1, 1]."""
    a = node_position(grid, det["min"])
    b = node_position(grid, det["max"])
    mn = [min(a[i], b[i]) for i in range(3)]
    mx = [max(a[i], b[i]) for i in range(3)]
    if int(grid.get("nz") or 0) == 0:
        mn[2] = -1.0
        mx[2] = 1.0
    return mn, mx


def aabb_signed_distance(mn: list[float], mx: list[float], x: list[float]) -> float:
    """Aabb::signed_distance. Positive outside, negative inside."""
    q = [abs(x[i] - 0.5 * (mn[i] + mx[i])) - 0.5 * (mx[i] - mn[i]) for i in range(3)]
    outside = math.sqrt(sum(max(v, 0.0) ** 2 for v in q))
    inside = min(max(q), 0.0)
    return outside + inside


def beta_of(mass: float, kinetic: float, c: float) -> float:
    """|v|/c from T = (γ − 1) m c². ε = T/(m c²), β = sqrt(ε(ε+2)) / (1+ε)."""
    eps = kinetic / (mass * c * c)
    return math.sqrt(eps * (eps + 2.0)) / (1.0 + eps)


def section_launch_detector(levels, lines: list[str]) -> None:
    lines.append("\n== Launch node against its detector (and gates) ==")
    lines.append(
        "Signed distance is Aabb::signed_distance of Level::box_region. "
        "A 2D detector is a slab z ∈ [−1, 1]. The runner treats g ≤ 0 at launch as "
        "an immediate event, before gates and acceptance."
    )
    inside = []
    near_spread = []
    closest = []
    gate_hits = []
    for name, level in levels:
        grid = level["grid"]
        gates = []
        for g in level.get("gates") or []:
            gates.append(detector_box(grid, g))
        for i, shot in enumerate(level.get("shots") or []):
            x = node_position(grid, shot["launch"]["node"])
            sd = aabb_signed_distance(*detector_box(grid, shot["detector"]), x)
            acc = shot["detector"].get("acceptance") or {}
            has_acc = bool(acc.get("direction") or acc.get("kinetic") or acc.get("radiation"))
            closest.append((sd, name, i, has_acc))
            if sd <= 0.0:
                inside.append(f"{name} shot {i} sd={sd:.6g} acceptance={has_acc}")
            beam = shot.get("beam") or {}
            if beam:
                width = float(beam.get("width") or 0.0)
                length = float(beam.get("length") or 0.0)
                dist = (beam.get("distribution") or "gaussian").lower()
                factor = 3.0 if dist != "uniform" else 1.0
                reach = factor * math.hypot(width, length)
                if sd < reach:
                    near_spread.append(
                        f"{name} shot {i} sd={sd:.6g} spread_reach={reach:.6g} ({dist})"
                    )
            for gi, (mn, mx) in enumerate(gates):
                gsd = aabb_signed_distance(mn, mx, x)
                if gsd <= 0.0:
                    gate_hits.append(f"{name} shot {i} gate {gi} sd={gsd:.6g}")
        for i, fp in enumerate(level.get("free_particles") or []):
            if not fp.get("detector"):
                continue
            x = node_position(grid, fp["node"])
            sd = aabb_signed_distance(*detector_box(grid, fp["detector"]), x)
            if sd <= 0.0:
                inside.append(f"{name} free_particle {i} sd={sd:.6g}")
    closest.sort()
    lines.append(f"  launch nodes with g ≤ 0 on their detector: {inside or 'none'}")
    lines.append(f"  launch nodes inside a gate: {gate_hits or 'none'}")
    lines.append(
        "  beam samples whose box clearance is inside the distribution's reach "
        f"(Gaussian 3σ of width and length): {near_spread or 'none'}"
    )
    lines.append("  closest launch nodes (signed distance, file, shot, has acceptance):")
    for sd, name, i, has_acc in closest[:8]:
        lines.append(f"    {sd:.6g}  {name} shot {i}  acceptance={has_acc}")


def section_betas(levels, lines: list[str]) -> None:
    lines.append("\n== Launch speed in units of c ==")
    lines.append(
        "β = sqrt(ε(ε+2))/(1+ε), ε = T/(m c²), from Kinematics. "
        "A beam's energy spread reaches ±3σ (Gaussian) or ± the spread (uniform) of T₀. "
        "tapered adds up to 0.1c when the continuation points against the velocity, "
        "so β > 0.9 can make the past or the future superluminal."
    )
    rows = []
    for name, level in levels:
        c = level.get("physics", {}).get("c")
        if c is None:
            continue
        c = float(c)
        for i, shot in enumerate(level.get("shots") or []):
            mass = float(shot["particle"]["mass"])
            kinetic = float(shot["launch"]["kinetic_energy"])
            b0 = beta_of(mass, kinetic, c)
            beam = shot.get("beam") or {}
            spread = float(beam.get("energy_spread") or 0.0)
            dist = (beam.get("distribution") or "gaussian").lower() if beam else ""
            factor = 3.0 if beam and dist != "uniform" else (1.0 if beam else 0.0)
            t_hi = kinetic * (1.0 + factor * spread)
            b_hi = beta_of(mass, t_hi, c) if t_hi > 0.0 else b0
            rows.append((b_hi, b0, name, i, c, kinetic, mass, spread, dist or "shot"))
    rows.sort(reverse=True)
    over = [r for r in rows if r[0] > 0.9]
    lines.append(f"  shots whose upper-edge β exceeds 0.9: {len(over)}")
    for b_hi, b0, name, i, c, kinetic, mass, spread, dist in over:
        lines.append(
            f"    β_edge={b_hi:.6f} β₀={b0:.6f}  {name} shot {i}  "
            f"c={c:g} T={kinetic:g} m={mass:g} spread={spread:g} {dist}"
        )
    lines.append("  highest upper-edge speeds:")
    for b_hi, b0, name, i, c, kinetic, mass, spread, dist in rows[:12]:
        lines.append(
            f"    β_edge={b_hi:.6f} β₀={b0:.6f}  {name} shot {i}  c={c:g} T={kinetic:g}"
        )


def _tapered(r, v, a, tau: float, c: float):
    """beam::tapered, in the plane. Returns (position, velocity, acceleration)."""
    alen = math.hypot(a[0], a[1])
    t_scale = 0.1 * c / max(alen, 1e-300)
    u1 = 0.7
    w = 1.0 - u1
    u = tau / t_scale
    if abs(u) <= u1:
        return (
            (r[0] + v[0] * tau + a[0] * 0.5 * tau * tau, r[1] + v[1] * tau + a[1] * 0.5 * tau * tau),
            (v[0] + a[0] * tau, v[1] + a[1] * tau),
            a,
        )
    s = (abs(u) - u1) / w
    ln_cosh = s + math.log1p(math.exp(-2.0 * s)) - math.log(2.0)
    th = math.tanh(s)
    sign = 1.0 if u > 0.0 else -1.0
    disp = t_scale * t_scale * (0.5 * u1 * u1 + u1 * w * s + w * w * ln_cosh)
    vel = t_scale * sign * (u1 + w * th)
    return (
        (r[0] + v[0] * tau + a[0] * disp, r[1] + v[1] * tau + a[1] * disp),
        (v[0] + a[0] * vel, v[1] + a[1] * vel),
        (a[0] * (1.0 - th * th), a[1] * (1.0 - th * th)),
    )


def _accel_of(a, tau: float, c: float):
    alen = math.hypot(a[0], a[1])
    t_scale = 0.1 * c / max(alen, 1e-300)
    u = tau / t_scale
    u1 = 0.7
    w = 0.3
    if abs(u) <= u1:
        return a
    s = (abs(u) - u1) / w
    sech2 = 1.0 - math.tanh(s) ** 2
    return (a[0] * sech2, a[1] * sech2)


def section_tapered(lines: list[str]) -> None:
    lines.append("\n== tapered continuation and the speed of light ==")
    lines.append(
        "Transcription of beam::tapered. T = 0.1 c/|a|, u = τ/T, u₁ = 0.7. "
        "As |u| → ∞ the velocity is v + sign(u) â · 0.1 c. "
        "Braking (a opposite v) makes the past faster; speeding up makes the future faster."
    )
    c = 1.0
    cases = [
        ("brake, |v|=0.95c, a=−c", (0.95, 0.0), (-1.0, 0.0), -1.0),
        ("brake, |v|=0.95c, a=−c, future", (0.95, 0.0), (-1.0, 0.0), 1.0),
        ("speed up, |v|=0.95c, a=+c, future", (0.95, 0.0), (1.0, 0.0), 1.0),
        ("across, |v|=0.995c, a=+y", (0.995, 0.0), (0.0, 1.0), -1.0),
    ]
    for label, v, a, direction in cases:
        # Far into the taper: |u| = 20.
        alen = math.hypot(*a)
        t_scale = 0.1 * c / alen
        tau = direction * 20.0 * t_scale
        _, vel, _ = _tapered((0.0, 0.0), v, a, tau, c)
        speed = math.hypot(*vel)
        # Edge of the constant-acceleration piece.
        tau_edge = direction * 0.7 * t_scale
        _, vel_e, _ = _tapered((0.0, 0.0), v, a, tau_edge, c)
        lines.append(
            f"  {label}: |v| at |u|=0.7 is {math.hypot(*vel_e):.6f} c, "
            f"at |u|=20 is {speed:.6f} c"
        )

    # Closed form against integrating dv/dτ = a sech², dx/dτ = v.
    v0 = (0.95, 0.0)
    a0 = (-1.0, 0.0)
    tau_end = -0.5
    n = 8000
    h = tau_end / n
    state = (0.0, 0.0, v0[0], v0[1])

    def deriv_taper(tt, st):
        ax, ay = _accel_of(a0, tt, c)
        return (st[2], st[3], ax, ay)

    for i in range(n):
        state = _rk4(state, h * i, h, deriv_taper)
    closed_r, closed_v, _ = _tapered((0.0, 0.0), v0, a0, tau_end, c)
    lines.append(
        f"  closed form vs RK4 of a sech² over τ={tau_end}: "
        f"|Δr|={math.hypot(closed_r[0] - state[0], closed_r[1] - state[1]):.3e} "
        f"|Δv|={math.hypot(closed_v[0] - state[2], closed_v[1] - state[3]):.3e}"
    )

    # Newton's method as in accelerated_fields, 50 steps, no residual test.
    def newton(x, dt, r, v, a, c):
        tau = dt - math.hypot(x[0] - r[0], x[1] - r[1]) / c
        max_beta = 0.0
        for _ in range(50):
            rp, vp = _tapered(r, v, a, tau, c)[:2]
            max_beta = max(max_beta, math.hypot(*vp) / c)
            dx, dy = x[0] - rp[0], x[1] - rp[1]
            dist = math.hypot(dx, dy)
            g = c * (dt - tau) - dist
            slope = -c + (dx * vp[0] + dy * vp[1]) / dist
            next_tau = tau - g / slope
            done = abs(next_tau - tau) <= 4.0 * 2.220446049250313e-16 * max(abs(dt) + dist / c, 1e-300)
            tau = next_tau
            if done:
                break
        rp, vp, ap = _tapered(r, v, a, tau, c)
        dist = math.hypot(x[0] - rp[0], x[1] - rp[1])
        g = c * (dt - tau) - dist
        return tau, g, math.hypot(*vp) / c, max_beta, ap

    lines.append("  accelerated_fields' Newton loop (50 steps, stops on a step-size test, not on g = 0):")
    observers = [
        ("beside a braking charge", (0.0, 0.4), 0.0, (0.0, 0.0), (0.95, 0.0), (-1.0, 0.0)),
        ("ahead of a braking charge", (0.2, 0.0), 0.0, (0.0, 0.0), (0.95, 0.0), (-1.0, 0.0)),
        ("behind a speeding charge", (-0.3, 0.2), 0.0, (0.0, 0.0), (0.95, 0.0), (1.0, 0.0)),
    ]
    for label, x, dt, r, v, a in observers:
        tau, g, beta, max_beta, _ = newton(x, dt, r, v, a, 1.0)
        lines.append(
            f"    {label}: τ={tau:.6g} g={g:.3e} |v|/c at τ={beta:.6f} "
            f"max |v|/c during the iteration={max_beta:.6f}"
        )


def _rk4(state, t: float, h: float, deriv):
    def add(a, b, s):
        return tuple(ai + s * bi for ai, bi in zip(a, b))

    k1 = deriv(t, state)
    k2 = deriv(t + 0.5 * h, add(state, k1, 0.5 * h))
    k3 = deriv(t + 0.5 * h, add(state, k2, 0.5 * h))
    k4 = deriv(t + h, add(state, k3, h))
    return tuple(si + (h / 6.0) * (a + 2.0 * b + 2.0 * c + d) for si, a, b, c, d in zip(state, k1, k2, k3, k4))


def _cup_drift() -> tuple[float, float, float]:
    """Work a fading charge does after its energy has been booked as absorbed.

    Survivor starts at rest at the origin. The absorbed charge is at (2, 0), flies at
    speed 1 along +x, and its Coulomb factor is e^{−(π/2) t} (mouth width 2, entry speed 1:
    k = π/w, rate = k v_n). c = ∞, both charges 1, survivor mass 1. No other field.
    The force is the one in beam.rs: d * (q q_f factor / r³), d = x − x_fade.
    """
    rate = math.pi / 2.0

    def deriv(tt, st):
        x1x, x1y = 2.0 + tt, 0.0
        dx, dy = st[0] - x1x, st[1] - x1y
        r2 = dx * dx + dy * dy
        scale = math.exp(-rate * tt) / (r2 * math.sqrt(r2))
        return (st[2], st[3], dx * scale, dy * scale)

    state = (0.0, 0.0, 0.0, 0.0)
    t_end = 8.0
    n = 40000
    h = t_end / n
    for i in range(n):
        state = _rk4(state, h * i, h, deriv)
    ke = 0.5 * (state[2] * state[2] + state[3] * state[3])
    u0 = 1.0 / 2.0
    return ke, u0, t_end


def section_cup_budget(levels, lines: list[str]) -> None:
    lines.append("\n== Faraday cup: who can show the energy panel the fade disagrees with ==")
    lines.append(
        "instant_drain defaults off, so a detector is Fate::Cup. "
        "c = ∞ and beam_interaction (or free particles, or max_free > 0) make the "
        "panel say the total is conserved. Listed levels have at least two charged "
        "particles in the file (a beam count, or more than one charged shot)."
    )
    shown = []
    for name, level in levels:
        physics = level.get("physics") or {}
        if physics.get("c") is not None:
            continue
        if physics.get("instant_drain"):
            continue
        interacts = bool(physics.get("beam_interaction")) or bool(level.get("free_particles")) or int(
            (level.get("limits") or {}).get("max_free") or 0
        ) > 0
        if not interacts:
            continue
        charged = 0
        for shot in level.get("shots") or []:
            q = float(shot["particle"]["charge"])
            n = int((shot.get("beam") or {}).get("count") or 1)
            if q != 0.0:
                charged += n
        if charged >= 2:
            shown.append(f"{name} ({charged} charged particles)")
    lines.append(f"  {shown or 'none'}")
    ke, u0, t = _cup_drift()
    lines.append(
        "  Illustration, not a shipped flight: after booking, the survivor feels "
        "q1 q2 factor(t) R̂/R² and nothing in the budget changes except its kinetic energy. "
        f"U at entry = {u0:.6g}. At t={t:g}, the survivor's kinetic energy is {ke:.6g} "
        f"({ke / u0:.4f} of that U). An instant drain leaves it at 0."
    )


def section_elliptic(lines: list[str]) -> None:
    lines.append("\n== Coil bracket g(m): series against the direct elliptic combination ==")
    lines.append("Transcribed from magnetic.rs. The m^0 and m^1 coefficients are omitted in the source because they cancel.")
    for m in (1e-8, 1e-4, 0.02, 0.049999):
        direct_k, direct_e = elliptic_ke(m, 1.0 - m)
        direct = (1.0 - 0.5 * m) * direct_e - (1.0 - m) * direct_k
        series = loop_bracket(m, 1.0 - m)
        lead = 0.5 * math.pi * (3.0 / 16.0) * m * m
        rel = (series - direct) / direct if direct else float("nan")
        lines.append(
            f"  m={m:.6g}  series={series:.8g}  direct={direct:.8g}  "
            f"(series-direct)/direct={rel:.3e}  series/lead={series / lead:.6f}"
        )


def main() -> None:
    levels = load_levels()
    lines = [
        f"Levels read: {len(levels)} from {LEVELS}",
        "",
    ]
    section_clouds(levels, lines)
    section_stray(levels, lines)
    section_inventory(levels, lines)
    section_wave(lines)
    section_static_dipole(lines)
    section_launch_detector(levels, lines)
    section_betas(levels, lines)
    section_tapered(lines)
    section_cup_budget(levels, lines)
    section_elliptic(lines)
    text = "\n".join(lines) + "\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding="utf-8")
    print(text)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
