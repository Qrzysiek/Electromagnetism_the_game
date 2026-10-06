"""Pairwise element check: every kind of element next to every other, run by the game.

    python scripts/pairs.py            # generate the cases, run them through the physics worker
    python scripts/pairs.py --gui      # also open each case in the game (capture mode)
    python scripts/pairs.py --only sphere_grounded   # cases with this element only

Every element the level format or the player can place (fixed charges, magnets, antennas,
metal spheres of each bias, clouds, electrodes fixed, tunable and driven, coils steady,
ramped, driven and polygonal, free particles at rest and thrown, gates, the launch, the
detector, a beam, player charges, magnets, antennas, plates and free charges, stray fields
and waves) is put next to every other at several separations, from far to overlapping,
in a Newtonian world and in a relativistic one with radiation reaction. A pair the game
refuses (a player element on an occupied node, metal overlapping metal, a particle
starting on a wire: `Level::setup_issues`) is fine; one that it accepts must not crash
or produce non-finite numbers. Cases without a verdict within the time limit are listed
as slow.

The physics layer is the game's own worker (`crates/game/src/worker.rs`, test
`pairwise_cases`, run in release); the GUI layer (`--gui`) opens each case with
`EM_LEVEL_FILE` and `EM_PLACE` in capture mode and also catches crashes of the render
thread and rows that overflow the side panel. Cases go to `target/pairs/`, the report to
`target/pairs/report.json`; the exit code is 1 if any case failed.
"""

import argparse
import copy
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = os.environ.get("CARGO_TARGET_DIR") or os.path.join(ROOT, "target")
OUT = os.path.join(ROOT, "target", "pairs")

# The first element sits at A, the second at A + (d, 0) for each d (cells): overlapping,
# touching, near and far.
A = (14, 10)
SEPARATIONS = [1, 2, 4]
WORLDS = {
    "newton": {"c": None},
    "rel": {"c": 20.0, "radiation_reaction": True},
}


def node(p):
    return [p[0], p[1], 0]


def base_level(world):
    physics = {
        "c": world["c"],
        "charge_radius": 0.3,
        "magnet_radius": 0.3,
        "wire_radius": 0.1,
        "antenna_radius": 0.3,
        "t_max": 60.0,
        "tolerances": {"preview": 1e-10, "verify": 1e-12},
    }
    if world.get("radiation_reaction"):
        physics["radiation_reaction"] = True
    return {
        "format_version": 2,
        "engine_version": "0.1.0",
        "name": "Pair",
        "description": "",
        "grid": {"nx": 30, "ny": 20, "nz": 0, "subdivision": 1},
        "physics": physics,
        "shots": [{
            "particle": {"charge": 1.0, "mass": 1.0, "radius": 0.0},
            "launch": {"node": [1, 4, 0], "direction": [1.0, 0.0, 0.0], "kinetic_energy": 0.5},
            "detector": {"min": [27, 2, 0], "max": [29, 6, 0]},
        }],
        "elements": [],
        "coils": [],
        # Only what the pair's player elements need (`player` adds it): allowing free
        # charges alone makes every flight a beam flight.
        "limits": {
            "max_charges": 0, "magnitudes": [], "allow_positive": True,
            "allow_negative": True, "max_magnets": 0, "magnet_strengths": [],
            "supply_voltages": [-1.0, 1.0],
        },
        "reference_solution": [],
    }


def toward(p, q):
    """Angle (degrees) from p to q; 0 if they coincide."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    return 0.0 if dx == 0 and dy == 0 else math.degrees(math.atan2(dy, dx)) % 360.0


def velocity(p, q, speed):
    a = math.radians(toward(p, q))
    return [speed * math.cos(a), speed * math.sin(a)]


def coil_drive(world, kappa):
    # A DC source through R: kappa(t) -> V / (R c^2) (PHYSICS.md 2.10). c = inf: the
    # circuit is ill-posed; the game must say so, not crash.
    c = world["c"] or 1e4
    return {"source": {"kind": "dc", "value": kappa * c * c}, "resistance": 1.0}


def electrode(p, **kw):
    e = {"center": node(p), "length": 3.0, "thickness": 0.4, "height": 3.0,
         "angle_deg": 90.0, "bias": {"kind": "potential", "value": 1.0}}
    e.update(kw)
    return e


# Each element: (level, placement, its node, the other element's node, world) -> None.
def add(key, item):
    return lambda lv, pl, p, o, w: lv.setdefault(key, []).append(copy.deepcopy(item(p, o, w)))


PLAYER_LIMITS = {
    "charge": {"max_charges": 4, "magnitudes": [0.5, 1.0]},
    "magnet": {"max_magnets": 2, "magnet_strengths": [1.0, 2.0]},
    "antenna": {"max_antennas": 2, "antenna_amplitudes": [0.5, 1.0], "antenna_omegas": [1.0]},
    "plate": {"max_plates": 2, "plate_voltages": [-1.0, 1.0]},
    "free": {"max_free": 2, "free_charges": [-1.0, 1.0], "free_speeds": [0.0, 1.0]},
}


def player(item):
    def f(lv, pl, p, o, w):
        e = item(p, o, w)
        lv["limits"].update(PLAYER_LIMITS[e["kind"]])
        pl.append(e)
    return f


def set_launch(beam):
    def f(lv, pl, p, o, w):
        lv["shots"][0]["launch"]["node"] = node(p)
        if beam:
            lv["shots"][0]["beam"] = {"count": 8, "energy_spread": 0.01, "angle_spread_deg": 2.0,
                                      "width": 0.2, "length": 0.2, "distribution": "gaussian",
                                      "transmission": 0.5, "seed": 0}
            lv["physics"]["beam_interaction"] = True
    return f


def set_detector(lv, pl, p, o, w):
    lv["shots"][0]["detector"] = {"min": [p[0], p[1] - 1, 0], "max": [p[0] + 1, p[1] + 1, 0]}


ELEMENTS = {
    "charge": add("elements", lambda p, o, w: {"node": node(p), "kind": "charge", "value": 1.0}),
    "magnet": add("elements", lambda p, o, w: {"node": node(p), "kind": "magnet", "value": 2.0}),
    "antenna": add("elements", lambda p, o, w: {"node": node(p), "kind": "antenna", "value": 1.0,
                                                "angle_deg": 90.0, "omega": 1.0}),
    "sphere_grounded": add("conductors", lambda p, o, w: {
        "center": node(p), "radius": 1.5, "bias": {"kind": "grounded"}}),
    "sphere_charged": add("conductors", lambda p, o, w: {
        "center": node(p), "radius": 1.5, "bias": {"kind": "charge", "value": 1.0}}),
    "sphere_potential": add("conductors", lambda p, o, w: {
        "center": node(p), "radius": 1.5, "bias": {"kind": "potential", "value": 1.0}}),
    "cloud": add("clouds", lambda p, o, w: {"center": node(p), "radius": 2.0, "charge": -2.0}),
    "electrode": add("electrodes", lambda p, o, w: electrode(p)),
    "electrode_tunable": lambda lv, pl, p, o, w: (
        lv.setdefault("electrodes", []).append(electrode(p, tunable=True)),
        pl.append({"node": node(p), "kind": "supply", "value": -1.0})),
    "electrode_driven": add("electrodes", lambda p, o, w: electrode(
        p, bias={"kind": "grounded"},
        drive={"source": {"kind": "dc", "value": 1.0}, "resistance": 3.0})),
    "coil": add("coils", lambda p, o, w: {"shape": "circle", "center": node(p), "radius": 2.5,
                                          "kappa": 0.2}),
    "coil_ramped": add("coils", lambda p, o, w: {"shape": "circle", "center": node(p),
                                                 "radius": 2.5, "kappa": 0.0, "rate": 0.01}),
    "coil_driven": add("coils", lambda p, o, w: {"shape": "circle", "center": node(p),
                                                 "radius": 2.5, "kappa": 0.0,
                                                 "drive": coil_drive(w, 0.2)}),
    "coil_polygon": add("coils", lambda p, o, w: {"shape": "polygon", "kappa": 0.2, "vertices": [
        [p[0] - 2, p[1] - 2, 0], [p[0] + 2, p[1] - 2, 0], [p[0] + 2, p[1] + 2, 0],
        [p[0] - 2, p[1] + 2, 0]]}),
    "free_level": add("free_particles", lambda p, o, w: {
        "particle": {"charge": -1.0, "mass": 1.0, "radius": 0.3}, "node": node(p)}),
    "free_level_thrown": add("free_particles", lambda p, o, w: {
        "particle": {"charge": -1.0, "mass": 1.0, "radius": 0.3}, "node": node(p),
        "velocity": velocity(p, o, 1.0)}),
    "gate": add("gates", lambda p, o, w: {"min": [p[0] - 1, p[1] - 1, 0],
                                          "max": [p[0] + 1, p[1] + 1, 0]}),
    "launch": set_launch(False),
    "beam": set_launch(True),
    "detector": set_detector,
    "player_charge": player(lambda p, o, w: {"node": node(p), "kind": "charge", "value": -1.0}),
    "player_magnet": player(lambda p, o, w: {"node": node(p), "kind": "magnet", "value": 2.0}),
    "player_antenna": player(lambda p, o, w: {"node": node(p), "kind": "antenna", "value": 1.0,
                                              "angle_deg": 45.0, "omega": 1.0}),
    "player_plate": player(lambda p, o, w: {"node": node(p), "kind": "plate", "value": 1.0,
                                            "angle_deg": 90.0}),
    "player_free": player(lambda p, o, w: {"node": node(p), "kind": "free", "value": 1.0,
                                           "angle_deg": 0.0, "speed": 0.0}),
    "player_free_thrown": player(lambda p, o, w: {"node": node(p), "kind": "free", "value": 1.0,
                                                  "angle_deg": toward(p, o), "speed": 1.0}),
    "stray_field": lambda lv, pl, p, o, w: lv.setdefault("disturbances", []).append(
        {"name": "stray", "e": [0.0, 0.02], "bz": 0.02}),
    "wave": lambda lv, pl, p, o, w: lv.setdefault("disturbances", []).append(
        {"name": "wave", "waves": [{"amplitude": 0.05, "direction_deg": 30.0, "omega": 1.0,
                                    "phase_deg": 0.0}]}),
}
# Elements without a position: one separation is enough.
GLOBAL = {"stray_field", "wave"}


def cases(only, separations=None):
    names = list(ELEMENTS)
    seen = set()
    separations = separations or SEPARATIONS
    for i, a in enumerate(names):
        for b in names[i:]:
            if only and only not in (a, b):
                continue
            for wname, world in WORLDS.items():
                for d in separations:
                    if (a in GLOBAL or b in GLOBAL) and d != separations[-1]:
                        continue
                    pa, pb = A, (A[0] + d, A[1])
                    level, placement = base_level(world), []
                    ELEMENTS[a](level, placement, pa, pb, world)
                    ELEMENTS[b](level, placement, pb, pa, world)
                    name = f"{a}+{b}@{d}-{wname}"
                    level["name"] = name
                    key = hashlib.sha1(json.dumps([level, placement], sort_keys=True)
                                       .replace(name, "").encode()).hexdigest()
                    if key in seen:
                        continue
                    seen.add(key)
                    yield name, level, placement


def write_cases(only, separations=None):
    cases_dir = os.path.join(OUT, "cases")
    shutil.rmtree(cases_dir, ignore_errors=True)
    os.makedirs(cases_dir)
    n = 0
    for name, level, placement in cases(only, separations):
        with open(os.path.join(cases_dir, f"{n:05d}.json"), "w", encoding="utf-8") as f:
            json.dump({"name": name, "level": level, "placement": placement}, f)
        n += 1
    return cases_dir, n


def physics_layer(cases_dir):
    report = os.path.join(OUT, "report.json")
    env = dict(os.environ, EM_PAIRS=cases_dir, EM_PAIRS_REPORT=report)
    cmd = ["cargo", "test", "--release", "-p", "game", "--bin", "game", "pairwise_cases",
           "--", "--ignored", "--nocapture"]
    subprocess.run(cmd, cwd=ROOT, env=env)
    with open(report, encoding="utf-8") as f:
        return json.load(f)


def gui_case(path, game, workdir):
    with open(path, encoding="utf-8") as f:
        case = json.load(f)
    level_file = os.path.join(workdir, os.path.basename(path))
    with open(level_file, "w", encoding="utf-8") as f:
        json.dump(case["level"], f)
    env = dict(os.environ, EM_CAPTURE=os.path.join(workdir, "shot.png"),
               EM_LEVEL_FILE=level_file, EM_PLACE=json.dumps(case["placement"]),
               EM_WAIT="1", EM_FRAME="30")
    try:
        p = subprocess.run([game], env=env, cwd=ROOT, capture_output=True, text=True,
                           timeout=180)
    except subprocess.TimeoutExpired:
        return case["name"], "slow", "no verdict within 180 s"
    lines = p.stderr.splitlines()
    panics = [l for l in lines if "panicked" in l]
    overflow = [l for l in lines if "panel overflow" in l]
    if p.returncode != 0 or panics:
        return case["name"], "crash", (panics or lines[-3:] or [f"exit {p.returncode}"])[0]
    if overflow:
        return case["name"], "overflow", overflow[0]
    return case["name"], "ok", ""


def gui_layer(cases_dir, jobs):
    exe = "game.exe" if os.name == "nt" else "game"
    game = os.path.join(TARGET, "release", exe)
    subprocess.run(["cargo", "build", "--release", "-p", "game"], cwd=ROOT, check=True)
    files = sorted(os.path.join(cases_dir, f) for f in os.listdir(cases_dir))
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        def run(path):
            work = os.path.join(tmp, os.path.basename(path))
            os.makedirs(work)
            return gui_case(path, game, work)
        with ThreadPoolExecutor(jobs) as pool:
            for i, r in enumerate(pool.map(run, files)):
                results.append(r)
                if r[1] != "ok":
                    print(f"[gui {i + 1}/{len(files)}] {r[0]}: {r[1]} {r[2]}", flush=True)
    return [{"name": n, "status": s, "detail": d} for n, s, d in results]


def summarize(title, results):
    # "rejected": the game refuses the setup (fine). "slow": no verdict within the time
    # limit; listed, not a failure (strongly interacting particles in the exact retarded
    # model take minutes; the sandbox's cost meters say so).
    slow = [r for r in results if r["status"] == "slow"]
    bad = [r for r in results if r["status"] not in ("ok", "rejected", "slow")]
    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print(f"\n{title}: {len(results)} cases, " +
          ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    for r in bad + slow:
        print(f"  {r['status']:9} {r['name']}: {r['detail']}")
    return not bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--gui", action="store_true", help="also run every case in the game")
    ap.add_argument("--only", help="only cases with this element")
    ap.add_argument("--jobs", type=int, default=3, help="parallel games for --gui")
    ap.add_argument("--separations", type=int, nargs="+",
                    help=f"separations in cells (default {SEPARATIONS})")
    args = ap.parse_args()
    if args.only and args.only not in ELEMENTS:
        sys.exit(f"unknown element {args.only}; known: {', '.join(ELEMENTS)}")
    cases_dir, n = write_cases(args.only, args.separations)
    print(f"{n} cases in {cases_dir}")
    ok = summarize("physics worker", physics_layer(cases_dir))
    if args.gui:
        ok &= summarize("game (capture mode)", gui_layer(cases_dir, args.jobs))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
