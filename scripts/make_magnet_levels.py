"""Multi-shot levels with magnetic fields, modelled on real instruments.

Writes levels/2x_*.json (format v2). Detector positions come from simulating a reference
setup (the coil field is not perfectly uniform), see the DETECTORS table; the reference
solutions are stored and verified by the level tests.

Units: k = 1, cell = 1. Magnetic source strengths in field units (PHYSICS.md 2.2): a coil
of radius a with kappa = mu0 I / 4 pi has B = 2 pi kappa / a at its centre. A particle of
momentum p and charge q circles with radius p / (q B).
"""

import json
import math
import sys

M = 1e6


def shot(q, m, node, angle_deg, ke, det):
    a = math.radians(angle_deg)
    return {
        "particle": {"charge": q, "mass": m, "radius": 0.0},
        "launch": {"node": [node[0], node[1], 0], "direction": [math.cos(a), math.sin(a), 0.0],
                   "kinetic_energy": ke},
        "detector": {"min": [det[0], det[1], 0], "max": [det[2], det[3], 0]},
    }


def level(name, desc, grid, shots, elements, coils, limits, reference, c=5.0, t_max=400.0):
    return {
        "format_version": 2,
        "engine_version": "0.1.0",
        "name": name,
        "description": desc,
        "grid": {"nx": grid[0], "ny": grid[1], "nz": 0, "subdivision": 1},
        "physics": {"c": c, "charge_radius": 0.3, "magnet_radius": 0.3, "wire_radius": 0.1,
                    "t_max": t_max, "tolerances": {"preview": 1e-10, "verify": 1e-12}},
        "shots": shots,
        "elements": elements,
        "coils": coils,
        "limits": limits,
        "reference_solution": reference,
    }


def charge(x, y, q):
    return {"node": [x, y, 0], "kind": "charge", "value": q}


# Detector boxes per level and shot group; "probe" means a wide strip used to find where
# the reference setup lands (python make_magnet_levels.py probe).
PROBE = "probe" in sys.argv[1:]


def dempster():
    src = (9, 5)
    strip = (11, 2, 31, 5)
    det_light = strip if PROBE else DETECTORS["dempster_light"]
    det_heavy = strip if PROBE else DETECTORS["dempster_heavy"]
    shots = [
        shot(1e-6, 1.0, src, 85.0, 0.05, det_light),
        shot(1e-6, 1.0, src, 95.0, 0.05, det_light),
        shot(1e-6, 4.0, src, 85.0, 0.05, det_heavy),
        shot(1e-6, 4.0, src, 95.0, 0.05, det_heavy),
    ]
    return level(
        "Dempster's mass spectrometer (1918)",
        "Ions from the source are accelerated by your electrode, then bent by the magnetic "
        "field of the large coil. In a uniform field the radius is r = √(2mT)/(qB), so the "
        "heavy ions (m = 4) land twice as far out as the light ones (m = 1). Ions leaving "
        "at slightly different angles re-converge after half a turn (180° focusing). One "
        "accelerating electrode must bring both masses, at both angles, to their detectors.",
        (32, 22),
        shots,
        [],
        [{"shape": "circle", "center": [16, 11, 0], "radius": 10.5, "kappa": 5e5}],
        {"max_charges": 1, "magnitudes": [0.5 * M, 0.75 * M, 1 * M, 1.5 * M, 2 * M, 3 * M, 4 * M],
         "allow_positive": True, "allow_negative": False,
         "region": {"min": [8, 1, 0], "max": [10, 3, 0]}},
        [charge(9, 3, 2 * M)],
    )


def wien():
    # Source and detectors lie inside the coil: in the plane a closed coil can only be
    # entered by crossing its wire.
    src = (2, 10)
    strip = (26, 6, 28, 14)
    dets = [strip] * 3 if PROBE else [DETECTORS["wien_slow"], DETECTORS["wien_mid"],
                                      DETECTORS["wien_fast"]]
    shots = [
        shot(1e-6, 1.0, src, 0.0, 0.3, dets[0]),
        shot(1e-6, 1.0, src, 0.0, 0.5, dets[1]),
        shot(1e-6, 1.0, src, 0.0, 1.0, dets[2]),
    ]
    return level(
        "Wien filter",
        "Crossed electric and magnetic fields pass exactly one speed straight through: the "
        "electric force qE balances the magnetic force qvB when v = E/B. Slower ions are "
        "pushed one way, faster ones the other. The coil provides B; place charges to "
        "supply E so that each of the three ions reaches its own detector.",
        (30, 20),
        shots,
        [],
        # Counter-clockwise current: B out of the plane inside the coil, so for q > 0
        # moving along +x the magnetic force q v × B points along −y; E must point along +y.
        [{"shape": "polygon", "vertices": [[1, 5, 0], [29, 5, 0], [29, 15, 0], [1, 15, 0]],
          "kappa": 2.5e4}],
        {"max_charges": 4, "magnitudes": [0.1 * M, 0.2 * M, 0.4 * M],
         "allow_positive": True, "allow_negative": True,
         "region": {"min": [8, 7, 0], "max": [22, 13, 0]}},
        [charge(11, 7, 0.2 * M), charge(19, 7, 0.2 * M),
         charge(11, 13, -0.2 * M), charge(19, 13, -0.2 * M)],
    )


def calutron():
    src = (3, 10)
    shots = [
        shot(1e-6, 1.0, src, 0.0, 0.5, DETECTORS["calutron_light"]),
        shot(1e-6, 2.0, src, 0.0, 0.5, DETECTORS["calutron_heavy"]),
    ]
    return level(
        "Calutron",
        "Lawrence's calutron separated uranium isotopes on an industrial scale by bending "
        "ions in a magnetic field: at equal energy the radius r = √(2mT)/(qB) grows with "
        "the mass. Two isotopes leave the source together. Place magnets (uniformly "
        "magnetized spheres standing out of the plane) so that each reaches its own "
        "detector.",
        (30, 20),
        shots,
        [],
        [],
        {"max_charges": 0, "magnitudes": [], "allow_positive": True, "allow_negative": True,
         "max_magnets": 2, "magnet_strengths": [1 * M, 2 * M, 4 * M],
         "region": {"min": [8, 3, 0], "max": [22, 17, 0]}},
        [],
    )


DETECTORS = {
    # Reference landing points (probe): light 13.37, 13.56; heavy 22.66, 22.95 (at y = 5).
    "dempster_light": (13, 3, 14, 5),
    "dempster_heavy": (22, 3, 23, 5),
    # Reference exit heights at x = 26 (probe): slow 11.79, design speed 9.38, fast 8.12.
    "wien_slow": (26, 11, 28, 13),
    "wien_mid": (26, 9, 28, 10),
    "wien_fast": (26, 8, 28, 9),
    "calutron_light": (27, 13, 30, 17),
    "calutron_heavy": (27, 3, 30, 7),
}

LEVELS = {"21_dempster": dempster, "22_wien_filter": wien, "23_calutron": calutron}

if __name__ == "__main__":
    for key, make in LEVELS.items():
        with open(f"levels/{key}.json", "w", newline="\n", encoding="utf-8") as f:
            f.write(json.dumps(make(), indent=2, ensure_ascii=False) + "\n")
    print("wrote", len(LEVELS), "levels", "(probe detectors)" if PROBE else "")
