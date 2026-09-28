"""All shipped levels, as a curriculum (single source of truth).

    python scripts/levels.py            # write levels/NN_*.json
    python scripts/levels.py --only 07  # rebuild one level

Order: by how hard the phenomenon is to understand and use. Every new element or
concept first gets an easy introduction level; within a chapter the possibility space
and the rarity of solutions grow on average (measured with `generator analyze`), and
so does the number of elements needed.

Each level is described in physical terms. Where a level has a designed reference setup
(`reference=`), detectors can be placed automatically around where the reference
flights actually land (`Auto` detectors): the script writes a probe version with wide
detector strips, runs `generator check` on it, and puts a detector of the requested size
around each landing point. Levels without a designed reference get one from the solver
(`generator solve --write`).

Units: k = 1, cell = 1. Particles are weakly charged (q ~ 1e-6) and electrodes strongly
charged (Q ~ 1e6) so that neglected radiation stays far below the numerical accuracy
(PHYSICS.md section 8). Magnetic strengths are in field units (PHYSICS.md 2.2): a coil of
radius a with kappa = mu0 I / 4 pi has B = 2 pi kappa / a at its centre, and a particle of
momentum p circles with radius p / (q B).
"""

import json
import math
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The generator binary; EM_GENERATOR overrides it (e.g. a build in another target dir).
GENERATOR = os.environ.get("EM_GENERATOR") or os.path.join(ROOT, "target", "release", "generator")
M = 1e6


def charge(x, y, q):
    return {"node": [x, y, 0], "kind": "charge", "value": q}


def magnet(x, y, mu):
    return {"node": [x, y, 0], "kind": "magnet", "value": mu}


def circle_coil(cx, cy, radius, kappa):
    return {"shape": "circle", "center": [cx, cy, 0], "radius": radius, "kappa": kappa}


def rect_coil(x0, y0, x1, y1, kappa):
    return {"shape": "polygon", "vertices": [[x0, y0, 0], [x1, y0, 0], [x1, y1, 0], [x0, y1, 0]],
            "kappa": kappa}


class Auto:
    """Detector found by simulating the reference setup.

    `strip` is the probe box (x0, y0, x1, y1); `axis` is the coordinate the particle
    crosses when entering it ('x' for a vertical strip, 'y' for a horizontal one). The
    final detector keeps the strip's extent along that axis and is `size` nodes wide
    across it, around the landing point.
    """

    def __init__(self, strip, axis, size=1):
        self.strip, self.axis, self.size = strip, axis, size


def box(b, direction=None, kinetic=None, radiation=None):
    """Detector box; optional acceptance: direction = (axis_deg, half_angle_deg),
    kinetic = (min, max), radiation = a radiation goal (`radiation_goal`)."""
    d = {"min": [b[0], b[1], 0], "max": [b[2], b[3], 0]}
    acc = {}
    if direction is not None:
        acc["direction"] = list(direction)
    if kinetic is not None:
        acc["kinetic"] = list(kinetic)
    if radiation is not None:
        acc["radiation"] = radiation
    if acc:
        d["acceptance"] = acc
    return d


def beam(count, transmission, energy=0.01, angle_deg=0.5, width=0.2, length=0.2,
         distribution="gaussian", seed=0):
    """A shot fired as a beam (level::beam): spreads are Gaussian sigmas (truncated at
    3 sigma) or uniform half-widths; energy relative to T0."""
    return {"count": count, "energy_spread": energy, "angle_spread_deg": angle_deg,
            "width": width, "length": length, "distribution": distribution,
            "transmission": transmission, "seed": seed}


def radiation_goal(axis_deg, half_deg, energy, band=None):
    """Radiation goal of a shot's detector (PHYSICS.md §3.4): the energy per steradian the
    flight radiates into the directions axis ± half (degrees), in all frequencies or in
    `band` = (omega_min, omega_max), must lie in `energy` = (min, max)."""
    g = {"direction": [axis_deg, half_deg], "energy": list(energy)}
    if band is not None:
        g["band"] = list(band)
    return g


def shot(q, m, node, angle_deg, ke, detector, time=0.0, moment=0.0, beam=None):
    """`moment`: magnetic moment along z (spin up > 0, spin down < 0). `beam`: fire as a
    beam (see `beam`)."""
    a = math.radians(angle_deg)
    launch = {"node": [node[0], node[1], 0], "direction": [math.cos(a), math.sin(a), 0.0],
              "kinetic_energy": ke}
    if time:
        launch["time"] = time
    particle = {"charge": q, "mass": m, "radius": 0.0}
    if moment:
        particle["moment"] = moment
    extra = {"beam": beam} if beam else {}
    return {
        **extra,
        "particle": particle,
        "launch": launch,
        "detector": detector,
    }


def metal(x, y, r, kind="grounded", value=None):
    """A metal sphere: kind "grounded", "charge" (net charge value) or "potential"."""
    bias = {"kind": kind}
    if value is not None:
        bias["value"] = value
    return {"center": [x, y, 0], "radius": r, "bias": bias}


def cloud(x, y, r, q):
    """A charge cloud: a uniformly charged sphere particles fly through (Thomson's atom)."""
    return {"center": [x, y, 0], "radius": r, "charge": q}


def free_particle(q, m, node, velocity=(0.0, 0.0), radius=0.0, detector=None, moment=0.0):
    """A dynamic particle placed by the level (it moves and interacts; with a detector it
    is a goal that must arrive)."""
    d = {"particle": {"charge": q, "mass": m, "radius": radius}, "node": [node[0], node[1], 0],
         "velocity": list(velocity)}
    if moment:
        d["particle"]["moment"] = moment
    if detector is not None:
        d["detector"] = detector
    return d


def free_charge(x, y, q, angle_deg, speed):
    """The player's free charge (for references)."""
    return {"node": [x, y, 0], "kind": "free", "value": q, "angle_deg": angle_deg,
            "speed": speed}


def plate(x, y, length, thickness=0.4, height=4.0, angle_deg=0.0, kind="grounded",
          value=None, tunable=False):
    """A box electrode (plate, slab or wall) standing on the plane; `tunable`: the player
    sets its potential with a power supply."""
    bias = {"kind": kind}
    if value is not None:
        bias["value"] = value
    e = {"center": [x, y, 0], "length": length, "thickness": thickness,
         "height": height, "angle_deg": angle_deg, "bias": bias}
    if tunable:
        e["tunable"] = True
    return e


def player_plate(x, y, v, angle_deg=0.0):
    """A plate placed by the player (limits.plate size) at potential v."""
    e = {"node": [x, y, 0], "kind": "plate", "value": v}
    if angle_deg:
        e["angle_deg"] = angle_deg
    return e


def supply(x, y, v):
    """The player's power supply of the tunable electrode centred at (x, y)."""
    return {"node": [x, y, 0], "kind": "supply", "value": v}


def antenna(x, y, p0, angle_deg=0.0):
    e = {"node": [x, y, 0], "kind": "antenna", "value": p0}
    if angle_deg:
        e["angle_deg"] = angle_deg
    return e


def stray(name, e=(0.0, 0.0), bz=0.0, waves=()):
    """One disturbance: uniform stray fields plus plane waves (PHYSICS.md 2.3)."""
    d = {"name": name, "e": list(e), "bz": bz}
    if waves:
        d["waves"] = list(waves)
    return d


def wave(amplitude, omega, phase_deg, direction_deg=0.0):
    """A plane wave travelling towards `direction_deg`, polarized along z x k. For c = inf
    it is a uniform field E0 cos(omega t + phase) along that polarization."""
    return {"amplitude": amplitude, "direction_deg": direction_deg, "omega": omega,
            "phase_deg": phase_deg}


def level(name, desc, grid=(30, 20), shots=(), elements=(), coils=(), max_charges=0,
          magnitudes=(), signs=(True, True), max_magnets=0, strengths=(), region=None,
          reference=None, c=5.0, t_max=400.0, disturbances=(), max_antennas=0,
          amplitudes=(), rf_omega=0.0, radiation_reaction=False, omegas=(), conductors=(),
          electrodes=(), max_plates=0, plate_voltages=(), plate_size=None, supplies=(),
          beam_interaction=False, gates=(), clouds=(), free_particles=(), max_free=0,
          free_charges=(), free_speeds=(), free_mass=1.0, free_radius=0.3):
    limits = {"max_charges": max_charges, "magnitudes": list(magnitudes),
              "allow_positive": signs[0], "allow_negative": signs[1],
              "max_magnets": max_magnets, "magnet_strengths": list(strengths)}
    if max_antennas:
        limits["max_antennas"] = max_antennas
        limits["antenna_amplitudes"] = list(amplitudes)
        if omegas:
            limits["antenna_omegas"] = list(omegas)
    if region:
        limits["region"] = box(region)
    if max_plates:
        limits["max_plates"] = max_plates
        limits["plate_voltages"] = list(plate_voltages)
        if plate_size:
            limits["plate"] = dict(zip(("length", "thickness", "height"), plate_size))
    if supplies:
        limits["supply_voltages"] = list(supplies)
    if max_free:
        limits["max_free"] = max_free
        limits["free_charges"] = list(free_charges)
        limits["free_speeds"] = list(free_speeds)
        if free_mass != 1.0:
            limits["free_mass"] = free_mass
        if free_radius != 0.3:
            limits["free_radius"] = free_radius
    return {
        "format_version": 2, "engine_version": "0.1.0", "name": name, "description": desc,
        "grid": {"nx": grid[0], "ny": grid[1], "nz": 0, "subdivision": 1},
        "physics": {"c": c, "charge_radius": 0.3, "magnet_radius": 0.3, "wire_radius": 0.1,
                    "antenna_radius": 0.3,
                    "t_max": t_max, "tolerances": {"preview": 1e-10, "verify": 1e-12},
                    **({"rf_omega": rf_omega} if rf_omega else {}),
                    **({"radiation_reaction": True} if radiation_reaction else {}),
                    **({"beam_interaction": True} if beam_interaction else {})},
        "shots": list(shots), "elements": list(elements), "coils": list(coils),
        "limits": limits, "reference_solution": list(reference or []),
        **({"disturbances": list(disturbances)} if disturbances else {}),
        **({"conductors": list(conductors)} if conductors else {}),
        **({"clouds": list(clouds)} if clouds else {}),
        **({"free_particles": list(free_particles)} if free_particles else {}),
        **({"electrodes": list(electrodes)} if electrodes else {}),
        **({"gates": list(gates)} if gates else {}),
    }


# =======================================================================================
# charges.

def first_bend():
    return level(
        "First bend",
        "Place a charge to bend the beam into the detector (one is enough). Like charges "
        "repel, unlike "
        "charges attract.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 16, 30, 20)))],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M])


def slingshot():
    return level(
        "Slingshot",
        "The detector is behind you. Attraction can turn the particle around; the fixed "
        "negative charge is a hint, but not enough on its own.",
        shots=[shot(1e-6, 1.0, (6, 3), 90.0, 0.5, box((0, 0, 3, 4)))],
        elements=[charge(15, 12, -2 * M)],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M])


def geiger_marsden():
    return level(
        "Geiger–Marsden (1909)",
        "Alpha particles fired at gold foil sometimes bounced straight back: the atom has a "
        "tiny, heavy, positive nucleus. Steer the alpha near the source so that it "
        "backscatters into the detector. The deflection follows $\\tan(\\theta/2) = qQ/(2T_0 b)$: "
        "the impact parameter b decides everything.",
        shots=[shot(2e-6, 4.0, (0, 11), 0.0, 4.0, box((0, 18, 1, 20)))],
        elements=[charge(18, 10, 8 * M)],
        max_charges=3,
        magnitudes=[m * M for m in (0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5)],
        region=(2, 2, 10, 18), c=20.0)


def around_the_wall():
    # The wall joined with injection (docs/CURRICULUM.md): one charge gets the beam
    # around the wall, a second one straightens it into the entrance. The search finds
    # no one-charge solution and many with two.
    return level(
        "Around the wall",
        "A wall of positive charges blocks the direct path, and the detector behind it is "
        "the entrance of the next stage: the beam must enter it moving along the axis "
        "(within ±8°). Bend the beam around the wall, then straighten it.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((28, 8, 30, 12), direction=(0.0, 8.0)))],
        elements=[charge(15, y, 0.5 * M) for y in range(6, 15)],
        max_charges=4, magnitudes=[0.5 * M, 1 * M, 2 * M, 4 * M])

def thomson_crt():
    return level(
        "Thomson's cathode-ray tube (1897)",
        "J. J. Thomson deflected cathode rays with charged plates and showed they are "
        "negative particles: electrons. Here an electron at 0.2 c must hit a one-cell spot "
        "on the screen. Two electron energies must hit their own one-cell spots: the "
        "deflection falls as 1/T, like in an oscilloscope. Build deflection plates from "
        "point charges.",
        shots=[shot(-1e-6, 1.0, (0, 10), 0.0, e, Auto((29, 0, 30, 20), "x")) for e in (0.4, 0.6)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 0.75, 1, 1.5)],
        region=(8, 3, 17, 17),
        reference=[charge(8, 4, -1 * M)])


# =======================================================================================
# several shots, one setup.

def twin_beams():
    return level(
        "Twin beams",
        "Two particles, two detectors, one setup: every shot must reach its own detector. "
        "Switch between the shots with [ and ], or show them all.",
        shots=[shot(1e-6, 1.0, (0, 13), 0.0, 0.5, box((27, 15, 30, 19))),
               shot(1e-6, 1.0, (0, 7), 0.0, 0.5, box((27, 1, 30, 5)))],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M])


def hemispherical_analyzer():
    return level(
        "Hemispherical analyzer (XPS)",
        "Photoelectron spectrometers sort particles by energy between two concentric "
        "hemispheres, where the field is exactly that of a point charge at the centre. "
        "Three energies leave the source; each must reach its own detector after bending "
        "through 180°. A circular orbit of radius R needs $|qQ| = p v R$; the others become "
        "ellipses.",
        shots=[shot(1e-6, 1.0, (4, 17), 0.0, e, Auto((1, 0, 3, 11), "x"))
               for e in (0.45, 0.5, 0.55)],
        max_charges=3, magnitudes=[m * M for m in (1, 2, 3, 4, 5, 6, 7, 8)],
        region=(1, 3, 12, 16),
        reference=[charge(4, 10, -7 * M)])


def reflectron():
    return level(
        "Reflectron",
        "Time-of-flight mass spectrometers turn ions around with an electrostatic mirror, "
        "the reflectron. Two ions with different energies must both come back to the "
        "detector beside the source.",
        shots=[shot(1e-6, 1.0, (0, 13), 0.0, e, box((0, 4, 2, 7))) for e in (0.4, 0.6)],
        max_charges=4, magnitudes=[1 * M, 2 * M, 4 * M], signs=(True, False),
        region=(18, 0, 29, 20))


def einzel_lens():
    return level(
        "Einzel lens",
        "An Einzel lens focuses charged-particle beams in electron microscopes and ion "
        "sources. Three rays leave the source at −12°, 0° and +12°; all three must meet in "
        "the small detector on the axis. Only charges of the particle's sign are "
        "available (a decelerating lens).",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, box((28, 9, 30, 11)))
               for a in (-12.0, 0.0, 12.0)],
        max_charges=4, magnitudes=[0.25 * M, 0.5 * M, 1 * M, 2 * M], signs=(True, False),
        region=(8, 1, 18, 19))


# =======================================================================================
# relativity.

def fast_lane():
    return level(
        "Fast lane",
        "A relativistic particle (v = 0.85 c) is hard to bend: its momentum grows as $\\gamma m v$. "
        "Charges that would easily steer a slow particle barely move this one.",
        # q = 1e-7 with 10x larger fixed charges: the same trajectories as q = 1e-6, but
        # 100x less (neglected) radiation, which scales as q²/(m c² r) at fixed qQ.
        shots=[shot(1e-7, 1.0, (0, 4), math.degrees(math.atan2(0.2, 1.0)), 2.0,
                    box((27, 14, 30, 18)))],
        elements=[charge(12, 8, -20 * M), charge(18, 12, 20 * M), charge(20, 5, -20 * M)],
        max_charges=3, magnitudes=[m * M for m in (20, 30, 40, 60, 80)], c=1.5)


def beta_spectrometer():
    return level(
        "Beta-ray spectrometer",
        "Beta decay electrons are fast. Two electrons, at v = 0.91 c (γ = 2.41) and "
        "v = 0.87 c (γ = 2.0), must each reach their own detector after bending around "
        "one charge. The circular-orbit condition is $|qQ| = \\gamma m v^2 R$, not the Newtonian "
        "$2T_0 R$. Trust relativity, not intuition. (The electron charge is small and the "
        "fixed charges large, so that radiation stays negligible.)",
        shots=[shot(-1e-7, 1.0, (4, 17), 0.0, e, Auto((1, 0, 3, 13), "x"))
               for e in (math.sqrt(2.0), 1.0)],
        max_charges=3, magnitudes=[m * M for m in (20, 40, 70, 100, 120, 140, 170, 200)],
        region=(1, 3, 12, 16), c=1.0,
        reference=[charge(4, 10, 140 * M)])


# =======================================================================================
# magnetic fields (coils placed by the level).

def first_coil():
    return level(
        "First coil",
        "The coil's magnetic field bends moving charges into circles: the force $q\\,\\mathbf{v}\\times\\mathbf{B}$ is "
        "perpendicular to the velocity, so it changes the direction but never the speed "
        "(watch the energy bars). Use a charge to steer the circling particle into the "
        "detector.",
        shots=[shot(1e-6, 1.0, (8, 10), 90.0, 0.5, box((14, 8, 17, 11)))],
        coils=[circle_coil(15, 10, 9.5, 1.9e5)],
        max_charges=3, magnitudes=[0.5 * M, 1 * M, 2 * M])


def dempster():
    shots = [shot(1e-6, m, (9, 5), a, 0.05, Auto((11, 2, 31, 5), "y"))
             for m in (1.0, 2.0, 4.0) for a in (85.0, 95.0)]
    return level(
        "Dempster's mass spectrometer (1918)",
        "Ions from the source are accelerated by your electrodes, then bent by the magnetic "
        "field of the large coil. In a uniform field $r = \\sqrt{2mT}/(qB)$, so masses 1, 2 and 4 "
        "land at different places; ions leaving at slightly different angles re-converge "
        "after half a turn (180° focusing). One setup must bring all six ions to the "
        "detectors of their masses.",
        grid=(32, 22), shots=shots, coils=[circle_coil(16, 11, 10.5, 5e5)],
        max_charges=3,
        magnitudes=[m * M for m in (0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4)],
        region=(2, 1, 15, 4),
        reference=[charge(9, 3, 2 * M)])


# The Wien filters' detectors: the selected speed (T0 = 0.5) must leave straight along
# the axis (within ±3°), as it does in crossed fields that balance; the slower and faster
# ions to either side. With a free exit one charge's non-uniform field could sort the
# three speeds (measured, `analyze --fewest`); a straight exit needs E that balances
# qvB all along the path. (At ±4° one placement of a single charge still passed.)
WIEN_DETECTORS = {0.3: box((26, 12, 28, 15)), 0.5: box((26, 9, 28, 11), direction=(0.0, 3.0)),
                  1.0: box((26, 5, 28, 8))}


def wien_filter():
    return level(
        "Wien filter",
        "Crossed electric and magnetic fields pass exactly one speed straight through: $qE$ "
        "balances $qvB$ when $v = E/B$. Slower ions are pushed one way, faster ones the other. "
        "The coil provides B; place charges to supply E so that the middle ion leaves "
        "straight along the axis (within ±3°) and the slower and faster ones reach their "
        "own detectors.",
        shots=[shot(1e-6, 1.0, (2, 10), 0.0, e, WIEN_DETECTORS[e]) for e in (0.3, 0.5, 1.0)],
        coils=[rect_coil(1, 5, 29, 15, 2.5e4)],
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.15, 0.2, 0.3, 0.4)],
        region=(6, 6, 24, 14))

def first_magnet():
    return level(
        "First magnet",
        "A magnet here is a uniformly magnetized sphere standing out of the plane (⊙) or "
        "into it (⊗). In the plane its field is perpendicular to the plane and falls as "
        "$1/r^3$, so it bends the particle sideways, strongly only nearby. Place a magnet to "
        "bend the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        max_magnets=3, strengths=[1 * M, 2 * M, 4 * M], region=(6, 2, 24, 18))


def calutron():
    return level(
        "Calutron",
        "Lawrence's calutron separated uranium isotopes by bending ions in a magnetic "
        "field: at equal energy the radius $r = \\sqrt{2mT}/(qB)$ grows with the mass. Three "
        "isotopes leave the source together; place magnets so that each reaches its own "
        "detector.",
        # Separated probe strips: in the same field the lightest isotope bends most.
        shots=[shot(1e-6, m, (3, 10), 0.0, 0.5, Auto((27, y0, 30, y1), "x", 3))
               for m, (y0, y1) in ((1.0, (15, 19)), (2.0, (10, 15)), (4.0, (4, 10)))],
        max_magnets=4, strengths=[m * M for m in (1, 2, 3, 4)], region=(6, 2, 24, 18),
        reference=None)


def build_wien_filter():
    return level(
        "Build a Wien filter",
        "No coil this time: build the whole velocity selector yourself from charges (for "
        "E) and magnets (for B). Again the middle ion must leave straight along the axis "
        "(within ±3°), the slower and faster ones in their own detectors.",
        shots=[shot(1e-6, 1.0, (2, 10), 0.0, e, WIEN_DETECTORS[e]) for e in (0.3, 0.5, 1.0)],
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.4)],
        max_magnets=2, strengths=[m * M for m in (1, 2, 3)], region=(6, 3, 24, 17))


# =======================================================================================
# Magnetic moments (PHYSICS.md 3.2). A neutral particle with a magnetic moment m along z
# (a spin state) feels no Lorentz force but the force m grad B_z. For a magnet mu the
# potential is U = m mu / r^3; flying past at distance b it is kicked by 4 m mu / (v b^3).
# With m = 1e-6 and mu = 1e6 a pass at 3 cells deflects by about 0.15 rad.

def stern_gerlach():
    return level(
        "Stern–Gerlach (1922)",
        "Stern and Gerlach sent silver atoms through a strongly non-uniform magnetic field. "
        "The atoms are neutral, but each carries a tiny magnetic moment, whose energy $-\\mathbf{m}\\cdot\\mathbf{B}$ "
        "changes where the field changes: the beam split in two, one part for each spin "
        "state. Here the two spin states fly as two shots. Place a magnet so that each "
        "reaches its own detector: one state is pushed towards stronger field, the other "
        "away (see the magnetic map).",
        shots=[shot(0.0, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2), moment=m)
               for m in (1e-6, -1e-6)],
        max_magnets=3, strengths=[m * M for m in (0.5, 1, 2)], region=(6, 3, 22, 17),
        reference=[magnet(14, 13, M)])


# =======================================================================================
# noise. Fields from outside the arena; one setup must work under each of them.

def stray_field():
    # A 0.5 T0 particle crosses the arena in t ~ 30; a stray E = 6.7e3 gives a = qE/m =
    # 6.7e-3 and pushes it ~3 cells off course.
    return level(
        "Stray field",
        "Somewhere in the building a high-voltage supply is switched on and off, and with "
        "it a weak uniform field across your beam line. The beam must reach the detector "
        "in both cases. Compare the two flights and aim between them.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 12, 30, 16)))],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2, 4)],
        disturbances=[stray("supply off"), stray("supply on", e=(0.0, -6.7e3))],
        c=None)


def mains_hum():
    # Uniform field E0 cos(wt + phi) along y (c = inf). A particle launched at phase phi
    # drifts with v_y = -(qE0 / m w) sin(phi): the hum acts like a random launch angle.
    omega, e0 = 0.6, 6e4
    return level(
        "Mains hum",
        "An alternating field from the mains wiring shakes the beam. Launched at different "
        "moments of the cycle, the particle drifts off at different angles, so no single "
        "aim works. A lens that images the source onto the detector does not care about "
        "the launch angle.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((28, 9, 30, 11)))],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(6, 1, 22, 19),
        disturbances=[stray(f"phase {p}°", waves=[wave(e0, omega, p, 0.0)])
                      for p in (0, 90, 180, 270)],
        c=None)


def earths_field():
    # Electron with p = 1: in B_z = 7.7e3 it circles with R = p / (qB) ~ 130 cells and
    # drifts ~3 cells sideways over the arena. The lab can be turned around, so the
    # stray field has either sign.
    b = 7.7e3
    return level(
        "Earth's field",
        "Cathode-ray tubes had to be adjusted to how they stood in the Earth's magnetic "
        "field. Your electron beam must hit the detector whichever way the apparatus faces: "
        "no stray field, or a weak field in or out of the plane. Magnetic bending grows "
        "with the time spent in the field.",
        shots=[shot(-1e-6, 1.0, (0, 10), 0.0, 0.5, box((28, 9, 30, 11)))],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(4, 2, 24, 18),
        disturbances=[stray("no field"), stray("facing north", bz=b),
                      stray("facing south", bz=-b)],
        c=None)


# =======================================================================================
# radio frequency. Antennas: oscillating dipoles with their exact retarded
# fields (PHYSICS.md 2.4). c = 5; with omega = 0.6 the wavelength 2 pi c / omega = 52
# cells: the arena is in the near and induction zones, where the field is strongest.
# The RF period is 10.5; a particle (v = 1) needs ~30 to cross the arena.

RF = 0.6
RF_PERIOD = 2 * math.pi / RF


def rf_kick():
    return level(
        "RF kick",
        "An antenna is a dipole driven by a radio-frequency generator: its field reverses "
        "every half period. The kick it gives a passing particle depends on the moment it "
        "passes. Orient and place an antenna to steer the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        max_antennas=3, amplitudes=[m * M for m in (1, 2, 4, 8)], rf_omega=RF)


def rf_separator():
    return level(
        "RF separator",
        "Identical particles, one launched half an RF period after the other. No static "
        "field can tell them apart; an oscillating one can. RF separators at CERN sorted "
        "kaons from pions this way. Send each bunch to its own detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, d, 30, d + 3)), time=t)
               for t, d in ((0.0, 14), (RF_PERIOD / 2, 3))],
        max_antennas=2, amplitudes=[m * M for m in (1, 2, 4, 8)],
        max_charges=1, magnitudes=[m * M for m in (0.5, 1, 2)],
        region=(4, 3, 24, 17), rf_omega=RF)


def streak_camera():
    return level(
        "Streak camera",
        "A streak camera turns time into position: a sweeping field deflects what arrives "
        "early one way and what arrives late the other. Three bunches leave the source a "
        "third of an RF period apart; each must hit its own spot on the screen.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((28, d, 30, d + 2)), time=t)
               for t, d in ((0.0, 15), (RF_PERIOD / 3, 9), (2 * RF_PERIOD / 3, 3))],
        max_antennas=2, amplitudes=[m * M for m in (1, 2, 4, 8)],
        max_charges=2, magnitudes=[m * M for m in (0.5, 1, 2)],
        region=(4, 2, 24, 18), rf_omega=RF)


# =======================================================================================
# radiation. The particle's own radiation is part of the physics
# (Landau-Lifshitz radiation reaction, PHYSICS.md 3.1). A strongly charged particle
# (q = m = 1, c = 2: classical radius q^2/(m c^2) = 0.25 cells, tau0 = 2q^2/(3mc^3) =
# 0.083) in a coil field B ~ 0.1 radiates about 0.1 % of its energy per turn and
# spirals inwards over ~30 turns. |F_rad|/|F_Lorentz| ~ 1e-2, where the LL
# treatment is accurate.

T_SYNC = 600.0


def synchrotron_light():
    # T0 = 0.43 at c = 2 gives p = 0.95 = r q B(r) at r = 6: the orbit circles the coil
    # axis. The field is axisymmetric about the axis (coil + a magnet on it), so without
    # radiation the canonical angular momentum is conserved and the particle cannot
    # reach the axis; radiation damping (rate ~ B^2) shrinks the orbit onto it.
    return level(
        "Synchrotron light",
        "An accelerated charge radiates, and the particle loses the energy it radiates. "
        "In a magnetic field it circles, radiates on every turn and slowly spirals "
        "inwards: synchrotron radiation, the reason electron rings need constant "
        "re-acceleration. Without radiation this particle could never reach the "
        "centre. Put a magnet on the axis so that it radiates fast enough to get there "
        "in time: the radiated power grows with the square of the field. Switch the map "
        "to 'particle field' to watch the radiation leave.",
        shots=[shot(1.0, 1.0, (15, 4), 180.0, 0.43, box((13, 8, 17, 12)))],
        coils=[circle_coil(15, 10, 9.5, 0.162)],
        # Checked: with radiation reaction switched off no allowed magnet solves it
        # (a magnet of 8 would, by field geometry alone, so it is not offered).
        max_magnets=1, strengths=[0.5, 1.0, 2.0, 4.0],
        region=(15, 10, 15, 10), c=2.0, t_max=T_SYNC, radiation_reaction=True)


def tune_the_rf():
    # Bunches 4 time units apart. An antenna at omega sees them at a phase difference of
    # 4 omega: they are kicked oppositely when that is near pi (omega ~ 0.785) and alike
    # when it is near 2 pi (omega ~ 1.57).
    return level(
        "Tune the RF",
        "Now you choose the frequency. Two identical bunches leave the source 4 time "
        "units apart and must reach different detectors. An antenna kicks them "
        "differently only if they meet it at different phases of its oscillation: the "
        "phase difference is $\\omega\\,\\Delta t$. Pick the frequency, then the orientation and place.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, d, 30, d + 3)), time=t)
               for t, d in ((0.0, 14), (4.0, 3))],
        max_antennas=2, amplitudes=[m * M for m in (2, 4, 8)],
        omegas=[0.2, 0.4, 0.8, 1.2, 1.6],
        max_charges=1, magnitudes=[m * M for m in (0.5, 1, 2)],
        region=(4, 3, 24, 17), rf_omega=RF)


# =======================================================================================
# Metals (PHYSICS.md 2.6). A metal sphere is an equipotential: held at a potential V it
# acts from outside like a charge V a at its centre, plus the charges that nearby charges
# induce on it (their images). Values: V = 4e5 on a sphere of radius 2.5 is 1e6 of charge.

def high_voltage_dome():
    return level(
        "High-voltage dome",
        "A metal sphere held at high voltage, like the dome of a Van de Graaff generator. "
        "From outside it pushes like a big charge at its centre, but it is a conductor: a "
        "charge you put near it pulls the dome's charge towards itself. The dome lifts the "
        "beam; add a charge to bring it into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 16, 30, 19)))],
        conductors=[metal(15, 4, 2.5, "potential", 3e5)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)])


def image_charge():
    return level(
        "Image charge",
        "A large grounded sphere. Whatever charge you place near it, the metal answers with "
        "an opposite image charge inside: near the surface your charge is almost cancelled, "
        "far away the pair acts like a dipole. Steer the beam around the sphere.",
        shots=[shot(1e-6, 1.0, (0, 6), 0.0, 0.5, box((27, 13, 30, 17)))],
        conductors=[metal(14, 10, 3.5, "grounded")],
        max_charges=3, magnitudes=[m * M for m in (0.5, 1, 2, 4)])


def polarized_sphere():
    return level(
        "Polarised sphere",
        "An isolated, uncharged metal sphere next to a strong charge. It stays neutral, but "
        "its charges separate: the side facing the charge takes the opposite sign, the far "
        "side the same sign. The sphere becomes a dipole that moves with every charge you "
        "add. Bring the beam to the detector.",
        shots=[shot(1e-6, 1.0, (0, 4), 0.0, 0.5, box((27, 15, 30, 18)))],
        elements=[charge(8, 16, 4 * M)],
        conductors=[metal(16, 10, 2.5, "charge", 0.0)],
        max_charges=3, magnitudes=[m * M for m in (0.5, 1, 2, 4)], region=(4, 2, 24, 18))


# =======================================================================================
# Delivering beams: detectors that also require a direction or an energy on arrival, as
# the entrance of a next stage does (PHYSICS.md 6.1).

def injection():
    return level(
        "Injection",
        "The detector is the entrance of the next accelerator stage: the beam must not only "
        "hit it, it must enter moving along the axis (+x, within ±8°), or the next stage "
        "loses it. The beam already reaches the entrance, but at 20°: straighten it. The "
        "cone shows what is accepted.",
        shots=[shot(1e-6, 1.0, (0, 5), 20.0, 0.5, box((28, 12, 30, 17), direction=(0.0, 8.0)))],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)])


def soft_landing():
    # T0 = 0.5; the ions must arrive with 0.02 <= T <= 0.08: a potential hill of
    # 0.42-0.48 T0 at the detector, without reflecting them.
    return level(
        "Soft landing",
        "Ions for surface science must land gently, or they would smash the surface: they "
        "may arrive with at most a quarter of their launch energy. The two electrodes at "
        "the target already brake them; add charges so that they land slowly, on the "
        "target, without turning back.",
        # A retarding electrode pair at the target (as deceleration optics have); the
        # player sets the final braking and keeps the beam on the target.
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((29, 9, 30, 11), kinetic=(0.02, 0.12)))],
        elements=[charge(28, 6, 0.6 * M), charge(28, 14, 0.6 * M)],
        max_charges=3, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.5, 0.75, 1)],
        region=(12, 2, 26, 18))


def collimator():
    return level(
        "Collimator",
        "Three rays leave the source at −10°, 0° and +10°. All three must reach the "
        "detector travelling parallel to the axis (within ±3°): a parallel beam, as optics "
        "downstream expect. Only charges of the particle's sign are available.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, box((28, 6, 30, 14), direction=(0.0, 3.0)))
               for a in (-10.0, 0.0, 10.0)],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], signs=(True, False),
        region=(4, 1, 20, 19))


# =======================================================================================
# Electrodes (PHYSICS.md 2.7): real plates and apertures with their fringe fields,
# computed by the boundary element method. V is the electrode potential; for the
# particle q = 1e-6, a potential difference dV costs q dV = 1e-6 dV of kinetic energy
# (T0 = 0.5 for the usual beams).

def deflection_plates():
    return level(
        "Deflection plates",
        "Real deflection plates, as in a cathode-ray tube: two metal plates at +V and −V. "
        "Between them the field is nearly uniform, but it bulges out at the ends (the "
        "fringe field) and the plates' charge rearranges when you bring charges near. Add "
        "a charge to bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        electrodes=[plate(11, 13, 8, kind="potential", value=-3e4),
                    plate(11, 7, 8, kind="potential", value=3e4)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], region=(16, 1, 26, 19))


def power_supply():
    # The top plate's supply: -120k, -60k, 60k or 120k lands the beam about 4, 2, -2 or
    # -4 cells from the axis; off (grounded) it flies straight.
    return level(
        "Power supply",
        "Instead of placing charges, turn a knob: the top deflection plate is connected to "
        "a power supply (yellow frame). Click the plate, or choose in the panel, to set "
        "its potential. A positive plate pushes positive ions away, a negative one pulls "
        "them. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2))],
        electrodes=[plate(11, 13, 8, tunable=True), plate(11, 7, 8)],
        supplies=[-1.2e5, -6e4, 6e4, 1.2e5], reference=[supply(11, 13, -1.2e5)])


def tune_the_lens():
    # The middle aperture's halves have a supply each; only both at 300k focus the three
    # rays into the detector (100k and 200k focus too weakly, 400k too strongly).
    apertures = []
    for x, tunable in ((10, False), (14, True), (18, False)):
        apertures += [plate(x, 15, 6, height=3.0, angle_deg=90.0, tunable=tunable),
                      plate(x, 5, 6, height=3.0, angle_deg=90.0, tunable=tunable)]
    return level(
        "Tune the lens",
        "A real Einzel lens: three metal apertures, the outer ones grounded. The two halves "
        "of the middle aperture have a power supply each. At the right voltage the lens "
        "focuses the three diverging rays into the small detector; too little and they "
        "stay apart, too much and they cross before it.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, Auto((28, 0, 30, 20), "x", 1))
               for a in (-6.0, 0.0, 6.0)],
        electrodes=apertures, supplies=[1e5, 2e5, 3e5, 4e5],
        reference=[supply(14, 15, 3e5), supply(14, 5, 3e5)])


def build_a_deflector():
    return level(
        "Build a deflector",
        "Now you build the electrodes. Place a metal plate at a potential of your choice "
        "near the beam (R turns it, Q/E changes its potential). A plate at a positive "
        "potential pushes positive ions away. Plates keep a cell away from other metal "
        "and from the other elements.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 1))],
        max_plates=3, plate_voltages=[-1e5, -5e4, 5e4, 1e5], region=(5, 3, 22, 17),
        reference=[player_plate(10, 12, 1e5)])


def shielding():
    return level(
        "Shielding",
        "A strong charge next to the beam line throws the beam out of the arena. A "
        "grounded metal plate between them screens its field: the charge it induces on the "
        "plate cancels much of the field on the far side. Place a grounded plate.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2))],
        elements=[charge(15, 14, 2 * M)],
        max_plates=3, plate_voltages=[0.0], region=(5, 3, 25, 17),
        reference=[player_plate(15, 12, 0.0)])


def real_einzel_lens():
    # Three apertures (plate pairs with a gap of 4 cells around y = 10) at x = 10, 14,
    # 18; the outer ones grounded, the middle one at V. For positive ions a positive V
    # decelerates and focuses.
    apertures = []
    for x, kind, value in ((10, "grounded", None), (14, "potential", 2.5e5),
                           (18, "grounded", None)):
        apertures += [plate(x, 15, 6, height=3.0, angle_deg=90.0, kind=kind, value=value),
                      plate(x, 5, 6, height=3.0, angle_deg=90.0, kind=kind, value=value)]
    return level(
        "Real Einzel lens",
        "The Einzel lens of arc 1, as it is really built: three metal apertures, the "
        "middle one at high voltage. Its focus is not quite on the detector, and the rays "
        "leave at ±6° now. Add charges to bring all three rays into the small "
        "detector.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, box((28, 9, 30, 11))) for a in (-6.0, 0.0, 6.0)],
        electrodes=apertures,
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.5)], region=(20, 2, 27, 18))


# =======================================================================================
# Multi-stage instruments (PHYSICS.md 6.2): gates are stages every flight must pass, in
# order, before its detector counts, e.g. first prepare the beam, then do the experiment.

def two_stages():
    return level(
        "Two stages",
        "Real instruments work in stages, and a particle only counts if it went through "
        "every stage in order. Here the particle must first pass the gate (the violet "
        "dashed box) and only then enter the detector. Reaching the detector without "
        "passing the gate does not count.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 20)))],
        gates=[box((10, 11, 13, 14))],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M])


# =======================================================================================
# Beams (PHYSICS.md 3.3): many particles at once, with spreads in position, direction and
# energy; the goal is a verified transmission. Harder, less forgiving versions of earlier
# levels. Interaction between the particles is exact for c = inf, so beam levels are
# Newtonian. To make the particles' repulsion (space charge) matter, their charge is
# scaled up by K, with mass and energy scaled alike: each particle's path in the fixed
# fields is unchanged (same q/m and T/q), while the repulsion grows with K.

K = 5e7


def space_charge():
    return level(
        "Space charge",
        "The first bend again, but now a whole beam of 16 particles, with a small spread in "
        "position, direction and energy. The particles repel each other, so the beam spreads "
        "out on the way (space charge). Bend it into the detector: at least 90 % of it must "
        "arrive, verified.",
        c=None, t_max=80.0, beam_interaction=True,
        shots=[shot(1e-6 * K, K, (0, 10), 0.0, 0.5 * K, box((27, 15, 30, 19)),
                    beam=beam(16, 0.9))],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M],
        reference=[charge(0, 0, 4 * M)])


def stern_gerlach_beam():
    return level(
        "Stern–Gerlach beam",
        "The experiment as it was done: a beam of both spin states at once, 12 atoms each, "
        "with a spread in speed and direction. Each spin state must reach its own detector, "
        "at least 90 % of it. The state pulled towards the magnet passes closer to it and "
        "fans out more. (Neutral atoms hardly interact; in this Newtonian level, not at "
        "all.)",
        c=None, t_max=80.0,
        shots=[shot(0.0, 1.0, (0, 10), 0.0, 0.5, box((27, 7, 30, 9)), moment=1e-6,
                    beam=beam(12, 0.9, energy=0.02)),
               shot(0.0, 1.0, (0, 10), 0.0, 0.5, box((27, 10, 30, 15)), moment=-1e-6,
                    beam=beam(12, 0.9, energy=0.02))],
        max_magnets=3, strengths=[m * M for m in (0.5, 1, 2)], region=(6, 3, 22, 17),
        reference=[magnet(14, 13, M)])


def collimated_beam():
    return level(
        "Collimated beam",
        "The collimator again, for a real beam: 16 particles leave the source with a spread "
        "of directions (σ = 4°) and repel each other on the way. At least 90 % must reach "
        "the detector travelling parallel to the axis (within ±3°). Only charges of the "
        "particle's sign are available.",
        c=None, t_max=80.0, beam_interaction=True,
        shots=[shot(1e-6 * K, K, (0, 10), 0.0, 0.5 * K, box((28, 6, 30, 14), direction=(0.0, 3.0)),
                    beam=beam(16, 0.9, angle_deg=4.0))],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], signs=(True, False),
        region=(4, 1, 20, 19), reference=None)


def velocity_selector():
    # The Wien filter of arc 3 for a beam: three speeds (T0 = 0.3, 0.5 and 1.0), 8 ions
    # each with a 2 % energy spread, repelling each other. The middle speed must leave
    # straight (within ±5°, wider than the single ions' ±3°: the beam spreads). With two
    # speeds and free exits one charge sufficed (`analyze --fewest`); now the search finds
    # no one-charge solution, and two-charge ones in 2 of 32 runs (the reference). The
    # Wien filter's own two charges fail here: space charge blows the slow beams apart.
    det = dict(WIEN_DETECTORS)
    det[0.5] = box((26, 9, 28, 11), direction=(0.0, 5.0))
    return level(
        "Velocity selector",
        "The Wien filter for a real beam: three groups of ions, slow, middle and fast, "
        "mixed and repelling each other. The coil provides B; place charges for E so that "
        "the middle group leaves straight along the axis (within ±5°) and the slower and "
        "faster groups reach their own detectors: at least 80 % of each.",
        c=None, t_max=80.0, beam_interaction=True,
        shots=[shot(1e-6 * K, K, (2, 10), 0.0, e * K, det[e], beam=beam(8, 0.8, energy=0.02))
               for e in (0.3, 0.5, 1.0)],
        coils=[rect_coil(1, 5, 29, 15, 2.5e4)],
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.15, 0.2, 0.3, 0.4)],
        region=(6, 6, 24, 14),
        reference=[charge(10, 6, 0.4 * M), charge(15, 6, 0.4 * M)])

def beam_preparation():
    # The collimated beam of the previous level, now as the first stage of an
    # experiment: collimated at the gate, then steered into the target.
    return level(
        "Beam preparation",
        "Experiments need a prepared beam. Stage 1: collimate it, so that at least 90 % of "
        "the particles cross the gate travelling parallel to the axis (within ±4°). Stage "
        "2: steer the prepared beam down into the target. A particle counts only if it "
        "went through both stages in order. Only charges of the particle's sign are "
        "available.",
        c=None, t_max=80.0, beam_interaction=True,
        shots=[shot(1e-6 * K, K, (0, 10), 0.0, 0.5 * K, box((27, 1, 30, 6)),
                    beam=beam(16, 0.9, angle_deg=4.0))],
        gates=[box((16, 6, 18, 15), direction=(0.0, 4.0))],
        max_charges=5, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], signs=(True, False),
        region=(3, 1, 26, 19), reference=None)


def relativistic_beam():
    # v = 0.8c (c = 5, gamma = 5/3, T0 = (gamma - 1) m c^2). Moving side by side, the
    # particles attract magnetically: the net repulsion is 1/gamma^2 = 0.36 of the Coulomb
    # force (PHYSICS.md 3.3, tests B6/B11). Quasi-static interaction (exact in the
    # velocities; checked against the retarded fields in the level tests), with radiation
    # reaction: the particles radiate about 1e-8 of their energy, too much to neglect.
    return level(
        "Relativistic beam",
        "A beam at 0.8 c. Its particles repel each other electrically, but charges moving "
        "side by side also attract magnetically: at this speed the net repulsion is only "
        "36 % ($1/\\gamma^2$) of the Coulomb force, which is why fast beams hold together. They are "
        "also harder to bend (momentum $\\gamma m v$). Bring at least 90 % of the beam into the "
        "detector.",
        c=5.0, t_max=20.0, beam_interaction=True, radiation_reaction=True,
        shots=[shot(0.06, 1.0, (0, 10), 0.0, (5 / 3 - 1) * 25.0, box((27, 14, 30, 19)),
                    beam=beam(16, 0.9))],
        max_charges=3, magnitudes=[200.0, 400.0, 800.0, 1600.0], reference=None)


# =======================================================================================
# Real instruments (SPEC 3, harder variants): an earlier level's idealised solution is
# built in (its reference setup becomes fixed elements), and a real effect the ideal
# design ignored breaks it. The player adds a few elements to make it work again. The
# build checks that the idealised design alone now fails (`MUST_FAIL_ALONE`).

def idealised(slug):
    """The shipped level `slug` with its reference setup built in as fixed elements."""
    i = [k for k, _ in LEVELS].index(slug) + 1
    path = os.path.join(ROOT, "levels", f"{i:02d}_{slug}.json")
    lvl = json.load(open(path, encoding="utf-8"))
    lvl["elements"] = lvl["elements"] + lvl["reference_solution"]
    lvl["reference_solution"] = []
    return lvl


def realistic(slug, name, desc, change, **player):
    """A realistic iteration of level `slug`: `change(lvl)` adds the real effect;
    `player` are the player's elements (arguments of `level`)."""
    lvl = idealised(slug)
    lvl["name"] = name
    lvl["description"] = desc
    change(lvl)
    lvl["limits"] = level("", "", **player)["limits"]
    return lvl


def as_beams(lvl, **spec):
    """Every shot fired as a beam (`beam` arguments)."""
    for s in lvl["shots"]:
        s["beam"] = beam(**spec)


def chromatic_aberration():
    # The Einzel lens of level 8, whose three rays now each come with an 8 % energy
    # spread: a lens focuses faster particles further away.
    def change(lvl):
        as_beams(lvl, count=8, transmission=0.9, energy=0.08, angle_deg=0.0, width=0.05,
                 length=0.05)
        lvl["physics"]["t_max"] = 80.0
    return realistic(
        "einzel_lens", "Chromatic aberration",
        "Your Einzel lens from arc 1 (Charges) is built in. Real ion sources are not "
        "monochromatic: every ray now carries an 8 % spread in energy, and a lens focuses "
        "faster ions further away (chromatic aberration). At least 90 % of each ray must "
        "still reach the small detector. Add a few charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.1, 0.25, 0.5, 1)],
        region=(4, 1, 26, 19))


def real_analyzer():
    # The hemispherical analyzer of level 9 with a source that has an angular spread.
    def change(lvl):
        as_beams(lvl, count=8, transmission=0.9, energy=0.0, angle_deg=4.0, width=0.05,
                 length=0.05)
        lvl["physics"]["t_max"] = 80.0
    return realistic(
        "hemispherical_analyzer", "Real analyser",
        "Your electron-energy analyser from arc 1 (Charges) is built in. A real source emits "
        "into a cone: the rays now leave within a few degrees of the axis (σ = 4°), and "
        "each energy must still land in its own detector, 90 % of it. Half a turn in a "
        "central field focuses directions only approximately. Add a few charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(1, 3, 12, 16))


def crt_earth_field():
    # The CRT of level 3 installed in the lab, facing north, in the Earth's magnetic field
    # (as in level 33). A tube is adjusted where it stands: one orientation. (Working in
    # every orientation, as in level 33, with the one-cell spots of level 3: at most one
    # search in 32 solved it, even at 40 % of the field.)
    b = 7.7e3

    def change(lvl):
        lvl["disturbances"] = [stray("facing north", bz=b)]
    return realistic(
        "thomson_crt", "CRT in the Earth's field",
        "Your cathode-ray tube from arc 1 (Charges) is built in and installed in the lab, facing "
        "north. In the Earth's magnetic field the electrons drift sideways on their way, "
        "and both beams miss their spots. Adjust the tube where it stands: add a few "
        "charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(4, 2, 26, 18))


def calutron_space_charge():
    # The calutron of level 28 at production currents: the isotope beams repel each
    # other (scaled charges, as in chapter 14; c = 5 with radiation reaction, since
    # the scaled charges would radiate more than may be neglected).
    # A quarter of the usual scaling and 5 ions per isotope: enough space charge to
    # spoil the ideal design, still correctable, and affordable with radiation reaction.
    def change(lvl):
        k = K / 4
        for s in lvl["shots"]:
            s["particle"]["charge"] *= k
            s["particle"]["mass"] *= k
            s["launch"]["kinetic_energy"] *= k
        as_beams(lvl, count=5, transmission=0.8, energy=0.01, angle_deg=0.5, width=0.2,
                 length=0.2)
        lvl["physics"]["beam_interaction"] = True
        lvl["physics"]["radiation_reaction"] = True
        lvl["physics"]["t_max"] = 80.0
    return realistic(
        "calutron", "Calutron at full current",
        "Your calutron from arc 3 (Relativity and magnetism) is built in. The wartime calutrons ran intense "
        "beams, and the ions repel each other: the isotope beams spread out and their "
        "spots grow into each other. At least 80 % of each isotope must still reach its "
        "own collector. Add a few charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(4, 2, 26, 18))


def beam_pipe():
    # Level 10's injection, but the beam runs in a grounded metal pipe: the wall next to
    # the deflecting charge screens it.
    def change(lvl):
        lvl["electrodes"] = [plate(15, 2, 26, height=4.0)]
    return realistic(
        "injection", "Beam pipe",
        "Your injection setup from arc 1 (Charges) is built in, but beams run in grounded metal "
        "pipes, and the pipe wall lies between the beam and your steering charge: the "
        "charge it induces on the wall cancels much of the field. Straighten the beam into "
        "the next stage again. Add a few charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(4, 4, 26, 19))


def soft_landing_current():
    # The soft landing of level 12 with an intense beam: braked to low energy, the ions
    # crowd and repel each other (space charge grows as the beam slows).
    # 8 ions at half the usual charge scaling (12 at full scaling: no search solved it).
    def change(lvl):
        k = K / 2
        for s in lvl["shots"]:
            s["particle"]["charge"] *= k
            s["particle"]["mass"] *= k
            s["launch"]["kinetic_energy"] *= k
            acc = s["detector"].get("acceptance", {})
            if "kinetic" in acc:
                acc["kinetic"] = [e * k for e in acc["kinetic"]]
        as_beams(lvl, count=8, transmission=0.9, energy=0.01, angle_deg=0.5, width=0.2,
                 length=0.2)
        lvl["physics"]["beam_interaction"] = True
        lvl["physics"]["radiation_reaction"] = True
        lvl["physics"]["t_max"] = 80.0
    return realistic(
        "soft_landing", "Soft landing, full current",
        "Your soft-landing optics from arc 1 (Charges) are built in. With a real ion current the "
        "braked ions crowd together and repel each other: the slower they get, the "
        "stronger the space charge, and the beam blows up just before the target. At "
        "least 90 % must still land gently on the target. Add a few charges.",
        change, max_charges=3, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.5)],
        region=(12, 2, 26, 18))


# =======================================================================================
# Finales (docs/CURRICULUM.md): multi-stage levels that need several modules. Their
# references are composed stage by stage, as a player would build them: each stage is
# solved as a sub-level whose detector is the next gate, with the earlier stages fixed.

def solve_stage(lvl, fixed, detectors, gates, prepare=None, **player):
    """Solve one stage of `lvl`: the shots aim at `detectors` (one per shot), after
    passing `gates`, with `fixed` elements built in and the player's elements `player`.
    `prepare(sub)` may adjust the sub-level further (e.g. fix the earlier stages' power
    supplies). Returns the solver's elements."""
    sub = json.loads(json.dumps(lvl))
    sub["elements"] = sub["elements"] + fixed
    for shot_, d in zip(sub["shots"], detectors):
        shot_["detector"] = d
    sub["gates"] = gates
    sub["limits"] = level("", "", **player)["limits"]
    sub["reference_solution"] = []
    if prepare:
        prepare(sub)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump(sub, f)
        path = f.name
    try:
        run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "600")
        found = json.load(open(path, encoding="utf-8"))["reference_solution"]
    finally:
        os.unlink(path)
    if not found:
        sys.exit(f"{lvl['name']}: no solution for the stage with detectors {detectors[:1]}")
    print(f"{lvl['name']}: stage solved with {len(found)} elements", flush=True)
    return found


def sorting_station():
    # Arc 1 finale. Stage 1 (lens): the six rays (two energies x three angles) pass a
    # two-cell gate. Stage 2 (deflector and lens): each energy is refocused onto its own
    # two-cell spot. (A version asking for parallel arrival, within 6-15 degrees, was not
    # solvable by the staged search.)
    energies, angles = (0.4, 0.6), (-8.0, 0.0, 8.0)
    gate = box((18, 9, 20, 11))
    slow = box((37, 13, 40, 15))
    fast = box((37, 5, 40, 7))
    shots = [shot(1e-6, 1.0, (0, 10), a, e, slow if e < 0.5 else fast)
             for e in energies for a in angles]
    mags = [m * M for m in (0.25, 0.5, 1, 2, 4)]
    lvl = level(
        "Sorting station",
        "Finale of the first arc. Six rays leave the source: two energies, each at −8°, 0° "
        "and +8°. First bring all six through the gate (a lens). Behind it they spread "
        "out again: sort them by energy and focus each energy onto its own small spot. "
        "Lens, deflector, lens: everything this arc taught, in one setup.",
        grid=(40, 20), shots=shots, gates=[gate], t_max=200.0,
        max_charges=8, magnitudes=mags, region=(2, 1, 37, 19))
    stage1 = solve_stage(lvl, [], [gate] * len(shots), [],
                         max_charges=3, magnitudes=mags, region=(2, 1, 17, 19))
    stage2 = solve_stage(lvl, stage1, [s_["detector"] for s_ in shots], [gate],
                         max_charges=4, magnitudes=mags, region=(21, 1, 36, 19))
    lvl["reference_solution"] = stage1 + stage2
    return lvl


def microscope_column():
    # Arc 2 finale. Stage 1 (condenser): a three-aperture Einzel lens, its middle
    # aperture's halves on power supplies, images the source through a small crossover
    # gate. Stage 2 (deflector and projector): a pair of deflection plates, the top one
    # on a supply, pushes the re-diverging beam down towards the spot beside the ion
    # pump, whose high voltage (a fixed charge) repels it; charges bring it there moving
    # along the axis (within ±20°). Small, low electrodes: a column of full-size ones
    # (8 plates, 6 x 3 cells) took 8.5 s to set up, over the sandbox budget; these take
    # 0.75 s. The search solves stage 2 with 2 elements in 4 of 48 runs.
    apertures = []
    for x, tunable in ((6, False), (10, True), (14, False)):
        apertures += [plate(x, 14, 4, height=2.0, angle_deg=90.0, tunable=tunable),
                      plate(x, 6, 4, height=2.0, angle_deg=90.0, tunable=tunable)]
    deflector = [plate(27, 13, 4, height=2.0, tunable=True)]
    gate = box((19, 9, 21, 11))
    spot = box((38, 3, 40, 5), direction=(0.0, 20.0))
    shots = [shot(1e-6, 1.0, (0, 10), a, 0.5, spot) for a in (-6.0, 0.0, 6.0)]
    volts = [-2e5, -1e5, 1e5, 2e5, 3e5, 4e5, 5e5]
    mags = [m * M for m in (0.1, 0.2, 0.3, 0.5, 1)]
    lvl = level(
        "Microscope column",
        "Finale of the second arc: the column of an electron microscope, in stages. The "
        "condenser lens (three apertures, the middle one on two power supplies) must bring "
        "all three rays through the small crossover gate. Behind it, the deflection plate "
        "(on its own supply) must push the beam down to the small spot beside the ion "
        "pump, whose high voltage repels it, and your charges must bring it in moving "
        "along the axis (within ±20°).",
        grid=(40, 20), shots=shots, gates=[gate], t_max=200.0,
        electrodes=apertures + deflector, elements=[charge(34, 1, 1 * M)],
        supplies=volts, max_charges=5, magnitudes=mags, region=(31, 1, 37, 19))

    def condenser_only(sub):
        # The deflector stays grounded while the condenser is tuned.
        for e in sub["electrodes"][6:]:
            e.pop("tunable", None)

    stage1 = solve_stage(lvl, [], [gate] * len(shots), [], prepare=condenser_only,
                         supplies=volts)

    def condenser_fixed(sub):
        # The condenser at the voltages of stage 1: fixed potentials.
        for s_ in stage1:
            for e in sub["electrodes"]:
                if e.get("tunable") and e["center"][:2] == s_["node"][:2]:
                    e.pop("tunable")
                    e["bias"] = {"kind": "potential", "value": s_["value"]}

    stage2 = solve_stage(lvl, [], [spot] * len(shots), [gate], prepare=condenser_fixed,
                         supplies=volts, max_charges=3, magnitudes=mags,
                         region=(31, 1, 37, 19))
    lvl["reference_solution"] = stage1 + stage2
    return lvl


def mass_spectrometer():
    # Arc 3 finale. Three masses (1, 2, 4) at the same energy, each at -5 and +5 degrees.
    # Stage 1 (lens of charges): electrostatic optics depend only on T/q, not on the
    # mass, so one lens focuses every mass through the gate. Stage 2 (magnetic sector of
    # magnets): the momentum sqrt(2mT) differs, so magnets separate the masses onto their
    # own collectors, found in probe strips and then shrunk to the landing points.
    masses, angles = (1.0, 2.0, 4.0), (-5.0, 5.0)
    gate = box((12, 9, 14, 11))
    strips = {1.0: box((37, 14, 40, 20)), 2.0: box((37, 7, 40, 13)), 4.0: box((37, 0, 40, 6))}
    shots = [shot(1e-6, m, (0, 10), a, 0.5, Auto((37, 0, 40, 20), "x", 2))
             for m in masses for a in angles]
    mags = [m * M for m in (0.25, 0.5, 1, 2)]
    strengths = [m * M for m in (1, 2, 3, 4)]
    lvl = level(
        "Mass spectrometer from parts",
        "Finale of the third arc. Ions of three masses (1, 2 and 4) leave the source with "
        "the same energy, each at −5° and +5°. First focus them all through the gate: an "
        "electrostatic lens bends every mass alike at the same energy, since only T/q "
        "matters. Then separate them: in a magnetic field the momentum $\\sqrt{2mT}$ decides, so "
        "each mass lands on its own collector.",
        grid=(40, 20), shots=shots, gates=[gate], t_max=300.0,
        max_charges=5, magnitudes=mags, max_magnets=4, strengths=strengths,
        region=(2, 1, 36, 19))
    probe = json.loads(json.dumps(lvl, default=lambda d: box(d.strip)))
    stage1 = solve_stage(probe, [], [gate] * len(shots), [],
                         max_charges=3, magnitudes=mags, region=(2, 1, 11, 19))
    stage2 = solve_stage(probe, stage1, [strips[m] for m in masses for _ in angles], [gate],
                         max_magnets=3, strengths=strengths, max_charges=2, magnitudes=mags,
                         region=(15, 1, 36, 19))
    lvl["reference_solution"] = stage1 + stage2
    return lvl


def rf_beam_line():
    # Arc 4 finale. Two identical bunches half an RF period apart, while a supply in the
    # building switches a weak stray field on and off. Stage 1 (steering that works in
    # both cases, charges): both bunches pass the gate, with and without the stray field.
    # Stage 2 (RF separator, antennas): each bunch reaches its own detector.
    # The gate also asks for the direction (the injection module): with a free gate one
    # charge sufficed for stage 1, and with ±8° and half this stray field too.
    gate = box((16, 9, 18, 11), direction=(0.0, 6.0))
    up, down = box((37, 14, 40, 16)), box((37, 4, 40, 6))
    shots = [shot(1e-6, 1.0, (0, 10), 0.0, 0.5, d, time=t)
             for t, d in ((0.0, up), (RF_PERIOD / 2, down))]
    mags = [m * M for m in (0.25, 0.5, 1, 2)]
    amps = [m * M for m in (1, 2, 4, 8)]
    lvl = level(
        "RF beam line",
        "Finale of the fourth arc. Two identical bunches leave the source half an RF period "
        "apart, while a supply in the building switches a weak stray field on and off. "
        "First bring both bunches through the gate, moving along the axis (within ±6°), "
        "whether the stray field is on or off. "
        "Then send each bunch to its own detector: only a field that changes in time can "
        "tell them apart.",
        grid=(40, 20), shots=shots, gates=[gate], t_max=200.0, rf_omega=RF,
        disturbances=[stray("supply off"), stray("supply on", e=(0.0, -1.3e4))],
        max_charges=5, magnitudes=mags, max_antennas=3, amplitudes=amps,
        region=(3, 1, 36, 19))
    stage1 = solve_stage(lvl, [], [gate] * len(shots), [],
                         max_charges=3, magnitudes=mags, region=(3, 1, 15, 19))
    stage2 = solve_stage(lvl, stage1, [up, down], [gate],
                         max_antennas=2, amplitudes=amps, max_charges=1, magnitudes=mags,
                         region=(19, 1, 36, 19))
    lvl["reference_solution"] = stage1 + stage2
    return lvl


def isotope_separator():
    # Arc 5 finale. A beam of two isotopes (masses 1 and 2, same energy), 8 ions each, with
    # a spread of directions, repelling each other (scaled charges, as in the beam arc).
    # Stage 1 (collimation, charges of the ions' sign): 80 % of each isotope crosses the
    # gate moving along the axis (within ±5°). Stage 2 (magnetic sector, magnets): each
    # isotope reaches its own collector.
    gate = box((14, 7, 16, 13), direction=(0.0, 5.0))
    # Both collectors on the side the field bends to: a magnetic field bends both
    # isotopes the same way, the lighter (less momentum) more. (Collectors on opposite
    # sides of the axis had no solution.)
    light, heavy = box((37, 13, 40, 19)), box((37, 5, 40, 11))
    shots = [shot(1e-6 * K, m * K, (0, 10), 0.0, 0.5 * K, d,
                  beam=beam(8, 0.8, energy=0.01, angle_deg=3.0))
             for m, d in ((1.0, light), (2.0, heavy))]
    mags = [m * M for m in (0.25, 0.5, 1, 2)]
    strengths = [m * M for m in (1, 2, 3, 4)]
    lvl = level(
        "Isotope separator",
        "Finale of the last arc. A beam of two isotopes, 8 ions each, leaves the source "
        "with a spread of directions, and the ions repel each other. First collimate it: "
        "at least 80 % of each isotope must cross the gate moving along the axis (within "
        "±5°). Then separate the isotopes with magnets, each onto its own collector. Only "
        "charges of the ions' sign are available.",
        grid=(40, 20), shots=shots, gates=[gate], t_max=120.0, c=None,
        beam_interaction=True, max_charges=5, magnitudes=mags, signs=(True, False),
        max_magnets=4, strengths=strengths, region=(2, 1, 36, 19))
    stage1 = solve_stage(lvl, [], [gate] * len(shots), [], max_charges=3,
                         magnitudes=mags, signs=(True, False), region=(2, 1, 13, 19))
    stage2 = solve_stage(lvl, stage1, [light, heavy], [gate], max_magnets=3,
                         strengths=strengths, max_charges=2, magnitudes=mags,
                         signs=(True, False), region=(17, 1, 36, 19))
    lvl["reference_solution"] = stage1 + stage2
    return lvl


# =======================================================================================
# Jackson series (docs/JACKSON.md): effects from J. D. Jackson, Classical Electrodynamics
# (3rd ed., 1999), each level citing its section or problem. Arc: charges in fields
# (Ch. 12). Uniform fields across the arena are given as a stray field (a single
# "disturbance"), which is exactly uniform and static, as the book assumes.

def jackson_exb_drift():
    # E = 1e5, B = 3.3e5: drift speed E/B = 0.30 cells per time unit; a particle starting
    # at rest cycloids with an amplitude 2mE/(qB^2) = 1.8 cells. Both signs drift the same
    # way; the guiding centres follow the equipotentials, so a charge steers both alike.
    return level(
        "Jackson §12.3: E×B drift",
        "Jackson §12.3 (motion in combined, uniform, static electric and magnetic fields). "
        "In crossed fields a particle starting at rest does not follow E: it rolls along a "
        "cycloid and drifts with the velocity $\\mathbf{E}\\times\\mathbf{B}/B^2$, whatever its charge and mass. A "
        "positive and a negative ion start here; both drift to the right, rolling in "
        "opposite senses. Their guiding centres follow the equipotentials (see the "
        "potential map): place a charge to bend the drift of both into the detector.",
        shots=[shot(q, 1.0, (2, 10), 0.0, 1e-4, box((27, 14, 30, 18))) for q in (1e-6, -1e-6)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], c=None, t_max=300.0,
        disturbances=[stray("crossed fields", e=(0.0, 1e5), bz=3.3e5)])


def jackson_van_allen():
    # A dipole "Earth" (magnet 100M at the centre). At r = 6, B = 4.6e5 and a particle
    # with v = 1 gyrates with a = 2.2 cells (a/R ~ 0.4, not << 1): its guiding centre
    # drifts around the Earth in about 95 time units, protons one way and electrons the
    # other (Pr. 12.9b: dphi/dt = -(3/2)(a/R)^2 omega_B for a << R, ~70 here).
    return level(
        "Jackson Pr. 12.9: Van Allen equator",
        "Jackson Problem 12.9 (the Van Allen belts). In the equatorial plane of a dipole "
        "Earth the field grows towards the Earth, and a gyrating particle drifts around it: "
        "protons one way, electrons the other (the gradient drift of §12.4). A proton and "
        "an electron start at the same place; the proton must reach the northern satellite, "
        "the electron the southern one. Shift their drift shells with charges or magnets.",
        shots=[shot(1e-6, 1.0, (15, 4), 0.0, 0.5, box((14, 17, 16, 19))),
               shot(-1e-6, 1.0, (15, 4), 180.0, 0.5, box((14, 0, 16, 2)))],
        elements=[magnet(15, 10, 100 * M)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        max_magnets=3, strengths=[m * M for m in (2, 4, 8, 16)], c=None, t_max=300.0)


def jackson_gradient_drift():
    # A uniform B (3.3e5) plus a row of fixed magnets along the bottom edge (3M each):
    # |B| grows towards the bottom and the gyrating ions drift along the edge, the
    # positive one left, the negative one right, leaving the arena at y ~ 3 after ~100
    # time units. The detectors are higher: no one-element solution (search), 27 of 32
    # runs with two.
    return level(
        "Jackson §12.4: gradient drift",
        "Jackson §12.4 (particle drifts in nonuniform, static magnetic fields). A row of "
        "magnets makes the field stronger towards the bottom. A gyrating ion's orbit is "
        "tighter where the field is stronger, so it drifts sideways, along the lines of "
        "equal |B|: the positive ion to the left, the negative one to the right. Bring each "
        "into its detector. A charge adds an E×B drift that is the same for both signs; a "
        "magnet changes the gradient, which moves them in opposite directions.",
        shots=[shot(1e-6, 1.0, (15, 8), 90.0, 0.5, box((0, 14, 2, 16))),
               shot(-1e-6, 1.0, (15, 8), 90.0, 0.5, box((28, 14, 30, 16)))],
        elements=[magnet(x, 0, 3 * M) for x in range(0, 31, 2)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        max_magnets=3, strengths=[m * M for m in (1, 2, 3, 4)], region=(1, 3, 29, 19),
        c=None, t_max=300.0, disturbances=[stray("uniform B", bz=3.3e5)])


def jackson_runaway():
    # c = 2, B0 = 1e5 (c B0 = 2e5) and E = 3e5 > c B0: there is no frame in which the
    # electric field vanishes, and the particle runs away along E (out of the top of the
    # arena within ~10 time units). Near a magnet |B| exceeds E/c, and there the particle
    # drifts (E x B, to the right). The search finds no one-magnet solution and two-magnet
    # ones in 14 of 32 runs. The particle's charge is 1e-8 and every field 100 times the
    # usual (the trajectories depend only on qE and qB): with q = 1e-6 it radiated 1.0e-8
    # of its (tiny, starting from rest) launch energy, above the 1e-10 that may be
    # neglected; radiation scales as q^2.
    return level(
        "Jackson Pr. 12.5: E×B runaway",
        "Jackson Problem 12.5 and §12.3. Crossed fields make a particle drift only if "
        "$|E| < c|B|$: then a frame moving with $\\mathbf{E}\\times\\mathbf{B}/B^2$ sees no electric field. Here the "
        "electric field is stronger than c times the magnetic one (c = 2 cells per time "
        "unit), no such frame exists, and the particle runs away along E. Near a magnet "
        "the field is strong enough: place magnets so that the particle drifts into the "
        "detector instead.",
        shots=[shot(1e-8, 1.0, (4, 10), 0.0, 1e-4, box((27, 8, 30, 12)))],
        max_magnets=6, strengths=[m * 100 * M for m in (1, 2, 4, 8)], c=2.0, t_max=300.0,
        region=(2, 2, 28, 18), disturbances=[stray("crossed fields", e=(0.0, 3e7), bz=1e7)])


def throw_a_charge():
    # Dynamic particles, introduction: the player's free charge (like sign, 0.5 or 1)
    # pushes the particle (charge 1, Newtonian) off course with its Coulomb field as it
    # flies past or towards it, and recoils itself. The search finds 507 one-charge
    # launches. (A first version with charges of 1e-6 showed only the contact collision:
    # the owner's review.)
    return level(
        "Throw a charge",
        "A new element: a free charge. It is not fixed: pull the handle behind it back like "
        "a slingshot to give it a velocity, and it flies off with the particle. Its field "
        "pushes the particle, and the particle's field pushes it back: momentum passes "
        "between them without any contact. Push the particle into the detector.",
        shots=[shot(1.0, 1.0, (0, 5), 0.0, 0.5, box((14, 17, 17, 19)))], c=None, t_max=100.0,
        max_free=3, free_charges=[0.5, 1.0], free_speeds=[0.0, 0.25, 0.5, 1.0])

def jackson_recoil():
    # Jackson §13.1: Coulomb scattering with a target free to recoil. Equal masses, target
    # at rest, elastic: momentum and energy conservation alone make the two leave at
    # right angles (for any force between them). The player throws the projectile (a free
    # charge, charge 1, mass 1); the target (the shot, at rest) must reach its detector
    # with T >= 0.1 (a distant push gives at most ~0.05). Nothing else is placed.
    return level(
        "Jackson §13.1: recoil at right angles",
        "Jackson §13.1 (energy transfer in a Coulomb collision), with a target that is free "
        "to recoil. Throw a particle of the same mass at the one at rest: whatever the "
        "force between them, conservation of momentum and energy sends them off at right "
        "angles to each other. Choose the aim so that the target recoils into its "
        "detector, fast enough to count; watch the thrown one leave at 90°.",
        grid=(40, 20), c=None, t_max=150.0,
        shots=[shot(1.0, 1.0, (28, 9), 90.0, 1e-9, box((30, 0, 36, 2), kinetic=(0.1, 10.0)))],
        max_free=3, free_charges=[1.0], free_speeds=[0.5, 1.0, 1.5], region=(1, 1, 10, 19))

def jackson_faraday():
    # Jackson §5.15: a coil whose current rises linearly (kappa = 1e4 t) induces
    # E = -dA/dt, which drives a charge at rest around the coil's axis; as B grows the
    # orbit also contracts (the flux through it is an adiabatic invariant, §12.5). The
    # detector asks for T >= 0.3, which only the induced field supplies: placed charges
    # are at most 0.2M and at least 4 cells from the launch point, so three of them give
    # at most 3 x 1e-6 x 2e5 / 4 = 0.15. The search finds 5 one-charge solutions.
    coil = circle_coil(15, 10, 8.0, 0.0)
    coil["rate"] = 1e4
    return level(
        "Jackson §5.15: Faraday's law",
        "Jackson §5.15: a changing magnetic flux induces an electric field that circulates "
        "around it, $\\oint \\mathbf{E}\\cdot d\\mathbf{l} = -d\\Phi/dt$. The current in this coil rises steadily; no charge is "
        "anywhere near the particle at rest, yet the induced field drives it around the "
        "axis, faster and faster, while the growing field tightens its orbit. Bring it into "
        "the detector with at least 0.3 of kinetic energy: only the induction can give it "
        "that much.",
        shots=[shot(1e-6, 1.0, (15, 6), 0.0, 1e-6, box((18, 11, 22, 15), kinetic=(0.3, 10.0)))],
        coils=[coil], c=None, t_max=200.0,
        max_charges=3, magnitudes=[m * M for m in (0.05, 0.1, 0.2)], region=(8, 10, 22, 17))


def jackson_knock_on():
    # Jackson Pr. 13.1: the player throws a heavy ion (a free charge: charge 1, mass 40)
    # past a light particle at rest (the shot: charge 1, mass 1), which takes the energy
    # T(b) = T_max / (1 + (b/b_min)^2) and recoils at (pi - theta)/2 from the ion's
    # direction (theta: its scattering angle in the ion's frame): forward for close
    # passes, sideways for distant ones. Nothing else is placed, so the only interaction
    # is the ion's field (earlier versions steered the ion with charges, then magnets,
    # which also acted on the target: the owner's reviews). The detector asks for
    # T >= 0.08: a distant push gives at most qQ/r ~ 0.03, a kick up to 2 m v^2 = 0.7.
    return level(
        "Jackson Pr. 13.1: knock-on",
        "Jackson Problem 13.1 (energy transfer in a Coulomb collision): throw a heavy ion "
        "(pull its slingshot handle back) past the light particle at rest. The ion's field "
        "kicks it, the harder the closer the ion passes, $T(b) = T_\\text{max}/\\left(1 + (b/b_\\text{min})^2\\right)$: "
        "a close pass throws it forward, a distant one sideways. Kick the light particle "
        "into its detector, fast enough to count (a gentle push from afar is not a kick).",
        grid=(40, 20), c=None, t_max=120.0,
        shots=[shot(1.0, 1.0, (30, 12), 90.0, 1e-9, box((32, 17, 37, 20), kinetic=(0.08, 10.0)))],
        max_free=3, free_charges=[1.0], free_speeds=[0.3, 0.45, 0.6], free_mass=40.0,
        region=(1, 1, 10, 19))

def jackson_stormer():
    # Jackson §12.1 (canonical momentum): in the equatorial plane of a dipole Earth
    # (magnet 100M) the field is axially symmetric, so L = r p_phi + q mu / r is
    # conserved (A_phi = mu / r^2). Launched straight at the Earth from r = 15 with p = 1,
    # L = q mu / 15 and the particle cannot come closer than r = 7.2 (Stormer's forbidden
    # region). Placed charges break the symmetry. No one-charge solution; two-charge ones
    # in 7 of 24 search runs.
    return level(
        "Jackson §12.1: Störmer's forbidden region",
        "Jackson §12.1 (canonical momentum). The Earth's dipole field is symmetric about "
        "its axis, so the particle's canonical angular momentum $r p_\\varphi + q r A_\\varphi$ cannot "
        "change: aimed straight at the Earth, it is turned away more than seven cells "
        "out, and no aim of the launch helps (Störmer's forbidden region, which keeps "
        "slow cosmic rays away from the equator). Your charges break the symmetry: bring "
        "the particle down to the Earth.",
        c=None, t_max=300.0,
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((13, 8, 17, 12)))],
        elements=[magnet(15, 10, 100 * M)],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        reference=[charge(1, 11, 0.5 * M), charge(11, 10, 2 * M)])


def jackson_magnetosphere():
    # Arc finale. A uniform "solar wind" (E = 1e5, B = 3.3e5) carries a proton and an
    # electron from rest on the left towards a dipole Earth (magnet 30M at (30, 10));
    # near the Earth the gradient drift splits them, the proton around the north, the
    # electron around the south. Stages 1 and 2 (E x B steering, charges): both pass the
    # upper gate, then the lower one. Stage 3 (drift shells, charges and magnets): each
    # reaches its own satellite. (With one gate it took 3 elements in all.)
    gate = box((12, 16, 14, 18))
    low = box((23, 5, 25, 7))
    # Satellites behind the Earth, one cell tall. Beside the Earth (two cells tall, with
    # the gate at y 14-16) the search needed only 3 elements in all; with the satellites
    # swapped (each particle three quarters around) it found none with up to 3 charges
    # and 4 magnets in stage 2: the solar wind sweeps the particles off their shells.
    north, south = box((35, 13, 37, 14)), box((35, 6, 37, 7))
    shots = [shot(1e-6, 1.0, (2, 10), 0.0, 1e-4, north), shot(-1e-6, 1.0, (2, 10), 0.0, 1e-4, south)]
    mags = [m * M for m in (0.25, 0.5, 1, 2)]
    strengths = [m * M for m in (2, 4, 8, 16)]
    lvl = level(
        "Jackson Ch. 12: magnetosphere",
        "Finale of the Jackson arc on charges in fields (§12.3, §12.4, Pr. 12.9). The solar "
        "wind's crossed fields carry a proton and an electron from rest towards the Earth, "
        "both drifting the same way. First lift both through the upper gate, then bring "
        "them down through the lower one: steer their E×B drift along the "
        "equipotentials. Near the Earth the gradient drift splits them, "
        "protons one way and electrons the other: bring the proton around to the small "
        "northern satellite behind the Earth, and the electron to the southern one.",
        grid=(40, 20), shots=shots, gates=[gate, low], t_max=400.0, c=None,
        elements=[magnet(30, 10, 30 * M)],
        max_charges=6, magnitudes=mags, max_magnets=5, strengths=strengths,
        region=(3, 1, 37, 19), disturbances=[stray("solar wind", e=(0.0, 1e5), bz=3.3e5)])
    stage1 = solve_stage(lvl, [], [gate] * 2, [], max_charges=3, magnitudes=mags,
                         region=(3, 1, 11, 19))
    stage2 = solve_stage(lvl, stage1, [low] * 2, [gate], max_charges=3, magnitudes=mags,
                         region=(15, 1, 22, 19))
    stage3 = solve_stage(lvl, stage1 + stage2, [north, south], [gate, low], max_charges=3,
                         magnitudes=mags, max_magnets=4, strengths=strengths,
                         region=(24, 1, 37, 19))
    lvl["reference_solution"] = stage1 + stage2 + stage3
    return lvl


# Jackson arc on conductors (Ch. 2-3). The golden-ratio levels use a particle with charge
# 1 (Newtonian, c = inf, so nothing radiates): the image force on it is then comparable
# to the sphere's own field, as in Problem 2.4. Test K5 checks the engine against the
# problem's answer.

def jackson_own_image():
    # A particle with charge 1 passes a grounded sphere (R = 3): its image -q R/d at
    # R^2/d pulls it in with the force q^2 R d / (d^2 - R^2)^2 (test K1).
    return level(
        "Jackson §2.2: its own image",
        "Jackson §2.2 (a point charge near a grounded conducting sphere). The sphere "
        "carries no charge of its own, yet it attracts every charge that passes: the "
        "charge induces an opposite image charge inside it, $-qR/d$ at the distance $R^2/d$ "
        "from the centre, and the image pulls. The pull grows steeply near the surface. "
        "Bring the particle into the detector.",
        shots=[shot(1.0, 1.0, (0, 6), 0.0, 0.3, box((27, 12, 30, 16)))],
        conductors=[metal(15, 11, 3, "grounded")],
        max_charges=3, magnitudes=[0.25, 0.5, 1.0, 2.0], c=None, t_max=300.0)


def jackson_two_spheres():
    return level(
        "Jackson Pr. 2.6: two spheres",
        "Jackson Problem 2.6: two conducting spheres, one charged, one neutral. Each sphere's "
        "charge induces an image in the other, which induces an image in the first, and so "
        "on: the neutral sphere becomes a dipole facing the charged one, and both respond "
        "to your charges and to the particle itself. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        conductors=[metal(11, 5, 2.5, "charge", 3e6), metal(19, 13, 2.5, "charge", 0.0)],
        max_charges=3, magnitudes=[m * M for m in (0.5, 1, 2, 4)], region=(2, 1, 26, 19))


def jackson_slot():
    # A grounded wall across the arena with a slot on the axis; a strong charge behind
    # the wall. Jackson §3.13 treats a circular hole in a plane; in the slice the hole is
    # a slot (the wall's plates are 4 cells high). The wall screens the charge except
    # through the slot, where its field leaks out and falls off quickly.
    return level(
        "Jackson §3.13: field through a hole",
        "Jackson §3.13 (a conducting plane with a hole): a grounded wall screens what is "
        "behind it, except through the opening, where the field leaks through and falls "
        "off within a few widths of the hole (in this slice the hole is a slot). The strong "
        "charge behind the wall is screened; the beam passes the slot. Your charges on this "
        "side are screened from the far side too. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 2, 30, 6)))],
        electrodes=[plate(15, 4, 8, angle_deg=90.0), plate(15, 16, 8, angle_deg=90.0)],
        elements=[charge(19, 13, 2 * M)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], region=(2, 1, 13, 19))


def jackson_golden_ratio():
    # Sphere R = 3 with Q = 1, particle q = 1: the force vanishes at d = 1.618 R, the
    # golden ratio (Pr. 2.4a). Launched at T0 = 0.3 straight at the sphere, the particle
    # passes the repulsive maximum and is pulled onto the sphere. The search finds no
    # one-charge solution and two-charge ones in 28 of 32 runs.
    return level(
        "Jackson Pr. 2.4: golden-ratio capture",
        "Jackson Problem 2.4: a charge near an isolated conducting sphere carrying the same "
        "charge is repelled far away, but attracted close to it: its image in the sphere "
        "wins inside 1.618 radii from the centre, the golden ratio (for equal charges). "
        "This particle has enough energy to get that close, and the sphere would catch "
        "it. Bring it around the sphere into the detector behind.",
        shots=[shot(1.0, 1.0, (0, 10), 0.0, 0.3, box((27, 8, 30, 12)))],
        conductors=[metal(15, 10, 3, "charge", 1.0)],
        max_charges=4, magnitudes=[0.25, 0.5, 1.0, 2.0], c=None, t_max=300.0)


def jackson_sphere_slalom():
    # Arc finale: three isolated spheres with the particle's charge (Q = q = 1, R = 2.5),
    # alternately below and above the axis; the particle must weave between them without
    # being captured (Pr. 2.4), through a gate after each of the first two spheres, then
    # into the detector behind the third. (With two spheres and one gate the stages took
    # 3 elements in all.)
    gates = [box((15, 9, 17, 11)), box((24, 9, 26, 11))]
    target = box((37, 12, 40, 16))
    shots = [shot(1.0, 1.0, (0, 4), 0.0, 0.3, target)]
    mags = [0.25, 0.5, 1.0, 2.0]
    lvl = level(
        "Jackson Ch. 2: sphere slalom",
        "Finale of the Jackson arc on conductors (Problems 2.4 and 2.6). Three isolated "
        "spheres carry the particle's own charge: far away they repel it, close by their "
        "images pull it in (inside 1.618 radii for equal charges), and each polarizes the "
        "others. Weave between them: through the two gates, then around the last sphere "
        "into the detector.",
        grid=(40, 20), shots=shots, gates=gates, t_max=400.0, c=None,
        conductors=[metal(10, 5, 2.5, "charge", 1.0), metal(21, 15, 2.5, "charge", 1.0),
                    metal(31, 5, 2.5, "charge", 1.0)],
        max_charges=8, magnitudes=mags, region=(2, 1, 37, 19))
    stage1 = solve_stage(lvl, [], [gates[0]], [], max_charges=3, magnitudes=mags,
                         region=(2, 1, 14, 19))
    stage2 = solve_stage(lvl, stage1, [gates[1]], gates[:1], max_charges=3, magnitudes=mags,
                         region=(18, 1, 23, 19))
    stage3 = solve_stage(lvl, stage1 + stage2, [target], gates, max_charges=3, magnitudes=mags,
                         region=(27, 1, 37, 19))
    lvl["reference_solution"] = stage1 + stage2 + stage3
    return lvl

# Jackson arc on radiation damping (Ch. 16). A classical atom: an electron (charge -1,
# mass 1) around a fixed nucleus (+1), c = 2, radiation reaction (Landau-Lifshitz):
# tau = 2q^2/(3mc^3) = 1/12. A circular orbit shrinks as r^3 = r0^3 - 6 tau t (Pr. 16.2;
# test R4). The detector around the nucleus accepts kinetic energy up to 0.26: an electron
# that has not radiated keeps its energy and arrives faster (at r <= 2.83 at least 0.27
# for the orbits here), so every solution needs radiation. Its 2 x 2 box ends the flight
# before |F_RR|/|F_L| exceeds the Landau-Lifshitz range.

CLASSICAL_ATOM_DETECTOR = box((13, 8, 17, 12), kinetic=(0.0, 0.26))


def atom_orbit(r, node, angle_deg, energy=None):
    """An electron launched at `node` (distance r from the nucleus at (15, 10)) with the
    circular speed for r, or the speed for the orbital `energy`."""
    v2 = 1.0 / r if energy is None else 2.0 * (energy + 1.0 / r)
    return shot(-1.0, 1.0, node, angle_deg, 0.5 * v2, CLASSICAL_ATOM_DETECTOR)


def jackson_classical_atom():
    # From r0 = 6 the circular orbit needs ~400 time units to shrink to the detector; the
    # level allows 250. A small charge far away makes the orbit eccentric (its nearly
    # uniform field changes the eccentricity secularly), and an eccentric orbit radiates
    # most at the perihelion (Pr. 16.3): the electron arrives at t = 125. Without
    # radiation the same setup is rejected (too fast).
    return level(
        "Jackson Pr. 16.2: the classical atom",
        "Jackson Problem 16.2 (and §16.2): an accelerated charge radiates, so a classical "
        "electron circling a nucleus loses energy and spirals in, $r^3 = r_0^3 - 9Z(c\\tau)^2 c t$. "
        "This one would take too long. The detector around the nucleus only counts an "
        "electron that has lost energy by radiating (arriving slowly). Make it radiate "
        "faster: an orbit that dips close to the nucleus radiates most there.",
        shots=[atom_orbit(6.0, (15, 4), 0.0)],
        elements=[charge(15, 10, 1.0)], max_charges=3, magnitudes=[0.1, 0.25, 0.5, 1.0],
        c=2.0, t_max=250.0, radiation_reaction=True, region=(1, 1, 29, 19),
        reference=[charge(2, 1, 0.25)])


def jackson_circularization():
    # Two electrons of the same orbital energy (-1/12): circular at r = 6, and an ellipse
    # launched at its aphelion r = 9. The ellipse reaches the detector early but too fast
    # (it must radiate first); radiation also circularizes it (Pr. 16.3). The search
    # finds two one-charge solutions (the reference) and two-charge ones in 12 of 24 runs.
    return level(
        "Jackson Pr. 16.3: orbits circularize",
        "Jackson Problem 16.3: an electron on an elliptic orbit radiates most near the "
        "nucleus, and its orbit becomes rounder as it decays. Two electrons with the same "
        "energy, one circling, one on an ellipse: the ellipse swings close to the nucleus "
        "early, but too fast to count. Bring both into the detector slowed by their "
        "radiation, in time.",
        shots=[atom_orbit(6.0, (15, 4), 0.0),
               atom_orbit(9.0, (15, 19), 180.0, energy=-1.0 / 12.0)],
        elements=[charge(15, 10, 1.0)], max_charges=3, magnitudes=[0.1, 0.25, 0.5, 1.0],
        c=2.0, t_max=150.0, radiation_reaction=True, region=(1, 1, 29, 19),
        reference=[charge(21, 7, -0.25)])


def jackson_three_orbits():
    # Arc finale: three electrons, the two of the previous level and a wide circular orbit
    # (r = 8, ~1000 time units to collapse alone), all within 200. The search finds no
    # solution with fewer than three charges, and three-charge ones in 2 of 24 runs (the
    # reference).
    return level(
        "Jackson Ch. 16: three orbits",
        "Finale of the Jackson arc on radiation damping (Problems 16.2 and 16.3). Three "
        "electrons circle the nucleus: on a small circle, on an ellipse, and on a wide "
        "circle that alone would take five times too long to collapse. Shape all three "
        "orbits so that each radiates enough, soon enough, to arrive slowly in the "
        "detector.",
        shots=[atom_orbit(6.0, (15, 4), 0.0),
               atom_orbit(9.0, (15, 19), 180.0, energy=-1.0 / 12.0),
               atom_orbit(8.0, (7, 10), 270.0)],
        elements=[charge(15, 10, 1.0)], max_charges=6, magnitudes=[0.1, 0.25, 0.5, 1.0],
        c=2.0, t_max=200.0, radiation_reaction=True, region=(1, 1, 29, 19),
        reference=[charge(12, 6, 0.25), charge(21, 11, -1.0), charge(22, 12, 0.5)])


# Jackson arc on bound charges (§16.7-16.8, Pr. 13.2): an electron bound harmonically
# inside a charge cloud (Thomson's atom): omega_0^2 = |qQ|/(m R^3). Tests C1 (the cloud's
# field) and C2 (Pr. 16.1: the radiating oscillator decays at Gamma = omega_0^2 tau).

def cloud_omega(q, big_q, r, m=1.0):
    return math.sqrt(abs(q * big_q) / (m * r ** 3))


def jackson_bound_charge():
    # Electron (q = -1e-6) bound in a cloud (Q = 1e6, R = 4, omega_0 = 0.125), circling at
    # r = 2 (v = omega_0 r = 0.25): an atom's electron is not at rest at the centre (the
    # owner's review). The restoring field grows to Q/R^2 = 6.3e4 at the edge. Pull it out
    # with charges.
    return level(
        "Jackson §16.7: a bound charge",
        "Jackson §16.7 models an atom's electron as a charge bound by a spring. Here the "
        "spring is real electrostatics: inside a sphere of uniform positive charge (J. J. "
        "Thomson's atom) the field grows linearly from the centre, so the electron is "
        "pulled back in proportion to how far it strays: it circles the centre with the "
        "same period at any radius. Only inside: outside the sphere the pull is "
        "Coulomb's $1/r^2$, no spring at all, so an electron with more energy than the "
        "spring can hold leaves. Pull it out of the atom into the detector.",
        shots=[shot(-1e-6, 1.0, (10, 8), 0.0, 0.5 * 0.25 ** 2, box((24, 8, 28, 12)))],
        clouds=[cloud(10, 10, 4.0, 1e6)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], c=None, t_max=300.0,
        region=(2, 1, 28, 19))


def jackson_resonance():
    # Electron oscillating (amplitude 1) in a cloud with omega_0 = 0.125; antennas only on
    # the left, far enough for a nearly uniform drive: a linear driven oscillator, whose
    # amplitude grows only at resonance. With omega_0 removed from the frequency list no
    # search finds a solution (even with two antennas); with it, 25 one-antenna ones.
    w0 = cloud_omega(1e-6, 1e6, 4.0)
    return level(
        "Jackson §16.8: resonance",
        "Jackson §16.8 (scattering and absorption of radiation by an oscillator). The bound "
        "electron oscillates at its natural frequency ω₀ = 0.125. A weak field that "
        "oscillates at ω₀ pushes it in step every cycle and its swing grows steadily; at any "
        "other frequency the pushes cancel out. Tune an antenna (frequency, orientation, "
        "place) so that the electron swings out of the atom into the detector. "
        "The atom is Thomson's: a sphere of uniform positive charge. Inside it the field "
        "grows in proportion to the distance from the centre, so an electron inside feels "
        "the pull of a spring ($\\omega_0^2 = |qQ|/(mR^3)$), the effective description of Jackson's "
        "bound charge. Outside the sphere the pull is Coulomb's $1/r^2$, no spring at all: an "
        "electron with more energy than the spring can hold leaves.",
        shots=[shot(-1e-6, 1.0, (15, 10), 90.0, 0.5 * w0 ** 2, box((20, 8, 22, 12)))],
        clouds=[cloud(15, 10, 4.0, 1e6)], rf_omega=0.3, t_max=400.0,
        max_antennas=2, amplitudes=[m * M for m in (0.1, 0.2, 0.5, 1)],
        omegas=[round(w0 / 2, 4), round(w0, 4), round(2 * w0, 4), 0.3],
        region=(2, 1, 7, 19))


def jackson_bound_knock():
    # Pr. 13.2: a passing charge transfers energy to a bound one. Unit charges, c = inf:
    # the electron (the shot: q = -1, m = 1) circles at r = 2 inside a cloud (Q = 1, R = 4,
    # omega_0 = 0.125, v = omega_0 r = 0.25); the player throws a heavy negative ion (a
    # free charge: q = -1, m = 40), whose field kicks the electron out. Nothing else is
    # placed (earlier versions steered the ion with elements that also acted on the
    # electron: the owner's reviews).
    w0 = cloud_omega(1.0, 1.0, 4.0)
    return level(
        "Jackson Pr. 13.2: a kick for a bound charge",
        "Jackson Problem 13.2: a charged particle flying past an atom gives its bound "
        "electron a kick. Throw a heavy negative ion (pull its slingshot handle back) past "
        "the atom, whose electron circles inside it. A quick pass is a sharp kick; a slow "
        "one lets the electron follow and hand the energy back. Kick the electron out of "
        "the atom into the detector above. "
        "The atom is Thomson's: a sphere of uniform positive charge. Inside it the field "
        "grows in proportion to the distance from the centre, so an electron inside feels "
        "the pull of a spring ($\\omega_0^2 = |qQ|/(mR^3)$), the effective description of Jackson's "
        "bound charge. Outside the sphere the pull is Coulomb's $1/r^2$, no spring at all: an "
        "electron with more energy than the spring can hold leaves.",
        grid=(40, 20), c=None, t_max=150.0,
        shots=[shot(-1.0, 1.0, (30, 8), 0.0, 0.5 * (w0 * 2.0) ** 2, box((27, 17, 33, 20)))],
        clouds=[cloud(30, 10, 4.0, 1.0)],
        max_free=3, free_charges=[-1.0], free_speeds=[0.3, 0.45, 0.6], free_mass=40.0,
        region=(1, 1, 10, 19))

def jackson_spectroscopy():
    # Arc finale: two atoms of different sizes (R = 4 and 3: omega_0 = 0.125 and 0.192),
    # an electron in each. Each must be driven out at its own resonance, without driving
    # the other out too early or the wrong way.
    w1, w2 = cloud_omega(1e-6, 1e6, 4.0), cloud_omega(1e-6, 1e6, 3.0)
    return level(
        "Jackson Ch. 16: spectroscopy",
        "Finale of the Jackson arc on bound charges (§16.7-16.8). Two atoms of different "
        "sizes, so two different natural frequencies: the big one rings at 0.125, the small "
        "one at 0.192. Each electron must leave its atom into its own detector: drive each "
        "at its own resonance, and keep each drive from upsetting the other atom. "
        "Each atom is Thomson's: a sphere of uniform positive charge. Inside it the field "
        "grows in proportion to the distance from the centre, so an electron inside feels "
        "the pull of a spring ($\\omega_0^2 = |qQ|/(mR^3)$), the effective description of Jackson's "
        "bound charge. Outside the sphere the pull is Coulomb's $1/r^2$, no spring at all: an "
        "electron with more energy than the spring can hold leaves.",
        grid=(40, 20), t_max=500.0, rf_omega=0.3,
        shots=[shot(-1e-6, 1.0, (13, 10), 90.0, 0.5 * w1 ** 2, box((18, 13, 21, 16))),
               shot(-1e-6, 1.0, (27, 10), 90.0, 0.5 * w2 ** 2, box((19, 4, 22, 7)))],
        clouds=[cloud(13, 10, 4.0, 1e6), cloud(27, 10, 3.0, 1e6)],
        max_antennas=4, amplitudes=[m * M for m in (0.1, 0.2, 0.5, 1)],
        omegas=[round(w1, 4), round(w2, 4), 0.25, 0.3],
        max_charges=2, magnitudes=[m * M for m in (0.25, 0.5, 1)],
        region=(1, 1, 39, 19))


def check_indirect(lvl):
    """Levels whose goal particles must be moved by a shot (the free particles with a
    detector). Two conditions, with the shots taken out of play (launched inside their
    own detectors, so they end at once and act on nothing):

    1. The solver finds no placement within the limits that brings every goal home (a
       heuristic search, but over every placement).
    2. No real influence: the strongest allowed element of each kind and sign at the
       placement node closest to each goal shifts it by less than 0.1 cells, against the
       same flight without the element, over the time the shot needs to reach it (for a
       goal that moves on its own, e.g. an orbiting electron, the comparison removes its
       own motion).

    Returns a list of problems."""
    goals = [f for f in lvl.get("free_particles", []) if f.get("detector")]
    if not goals:
        return []
    probe = json.loads(json.dumps(lvl))
    for sh in probe["shots"]:
        n = sh["launch"]["node"]
        sh["detector"] = {"min": [n[0], n[1], 0], "max": [n[0] + 1, n[1] + 1, 0]}
    probe["reference_solution"] = []

    def run(args, level):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                         encoding="utf-8") as f:
            json.dump(level, f)
            path = f.name
        try:
            return run_generator(*args[:1], path, *args[1:])
        finally:
            os.unlink(path)

    problems = []
    text = run(["solve", "--restarts", "16", "--iterations", "300"], probe)
    found = re.findall(r"verified (\d+)-charge solutions[^:]*: (\d+)", text)
    if not found:
        sys.exit("check_indirect: unexpected solver output")
    problems += [f"{int(n)} element(s): {k} solutions without the shots"
                 for n, k in found if int(k) > 0]

    L = lvl["limits"]
    region = L.get("region")
    lo = region["min"][:2] if region else [0, 0]
    hi = region["max"][:2] if region else [lvl["grid"]["nx"], lvl["grid"]["ny"]]
    kinds = []
    if L.get("max_charges"):
        kinds += [("charge", q) for q in (max(L["magnitudes"]), -max(L["magnitudes"]))]
    if L.get("max_magnets"):
        kinds += [("magnet", m) for m in (max(L["magnet_strengths"]),
                                          -max(L["magnet_strengths"]))]
    sh = lvl["shots"][0]
    v0 = math.sqrt(2 * sh["launch"]["kinetic_energy"] / sh["particle"]["mass"])
    x0 = sh["launch"]["node"]
    probe["physics"]["t_max"] = max(
        math.hypot(g["node"][0] - x0[0], g["node"][1] - x0[1]) for g in goals) / v0

    def ends(placement):
        probe["reference_solution"] = placement
        text = run(["check"], probe)
        out = re.findall(r"shot goal-\d+ .*ends at \(([-\d.]+), ([-\d.]+)\)", text)
        if len(out) != len(goals):
            sys.exit("the generator does not report goal particles: rebuild it")
        return [(float(a), float(b)) for a, b in out]

    # After the kick: the magnets turn a moving goal (charge q, mass m) at the rate
    # |q| B / m (non-relativistically, whatever its speed); over the flight time T it
    # turns by at most (|q|/m) B_max T, which over the distance L to its detector moves
    # its arrival by at most that times L. B_max: the strongest allowed magnets, all at
    # the region's nearest nodes, anywhere in the box spanned by the goal's start and its
    # detector (the in-plane dipole field is |mu|/r^3).
    if L.get("max_magnets"):
        mu = max(L["magnet_strengths"]) * L["max_magnets"]
        t_flight = lvl["physics"]["t_max"]
        for k, g in enumerate(goals):
            d = g["detector"]
            xs = [g["node"][0], d["min"][0], d["max"][0]]
            ys = [g["node"][1], d["min"][1], d["max"][1]]
            box_lo, box_hi = (min(xs), min(ys)), (max(xs), max(ys))
            # Nearest distance between the placement region and that box.
            dx = max(0, box_lo[0] - hi[0], lo[0] - box_hi[0])
            dy = max(0, box_lo[1] - hi[1], lo[1] - box_hi[1])
            r = max(math.hypot(dx, dy), 1e-9)
            b_max = mu / r ** 3
            part = g["particle"]
            reach = math.hypot(box_hi[0] - box_lo[0], box_hi[1] - box_lo[1])
            shift = abs(part["charge"]) / part["mass"] * b_max * t_flight * reach
            if shift > 0.3:
                problems.append(f"magnets could turn goal {k + 1} by up to {shift:.2f} cells "
                                f"after the kick (nearest distance {r:.1f})")

    base = ends([])
    for g in goals:
        gx, gy = g["node"][:2]
        node = [min(max(gx, lo[0]), hi[0]), min(max(gy, lo[1]), hi[1]), 0]
        for kind, value in kinds:
            for k, ((ax, ay), (bx, by)) in enumerate(
                    zip(base, ends([{"node": node, "kind": kind, "value": value}]))):
                shift = math.hypot(ax - bx, ay - by)
                if shift > 0.1:
                    problems.append(f"{kind} {value:g} at {node[:2]} shifts goal {k + 1} by "
                                    f"{shift:.2f} cells before the shot arrives")
    return problems


def needs_elements(lvl, fewer):
    """Check that the level is not solved with only `fewer` charges anywhere (the
    solver with that cap must fail): a finale needs its modules. A heuristic check: a
    failed search is not a proof."""
    probe = json.loads(json.dumps(lvl))
    probe["limits"]["max_charges"] = fewer
    probe["reference_solution"] = []
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump(probe, f)
        path = f.name
    try:
        run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "600")
        found = json.load(open(path, encoding="utf-8"))["reference_solution"]
    finally:
        os.unlink(path)
    return not found


# Finales and the fewest charges that must not suffice.
FINALES = {"sorting_station": 4}


# Realistic iterations: their built-in idealised design must fail on its own.
MUST_FAIL_ALONE = {"chromatic_aberration", "real_analyzer", "crt_earth_field",
                   "calutron_space_charge", "beam_pipe", "soft_landing_current"}

# --- Jackson Ch. 14: radiation goals (PHYSICS.md §3.4) -------------------------------
# A world at c = 2 with a charge of 1/40 at gamma = 3 (T0 = 8): it radiates noticeably (a
# few per cent of its energy in a tight bend), so radiation reaction is included
# (Landau-Lifshitz). With a charge of 0.1 the tight bends radiated half the energy and
# max |F_RR| / max |F_L| reached 0.19: the charge was lowered (the ratio scales as q^2 at
# fixed q x magnet strength, the trajectories do not change) to keep it near 0.012. A
# radiation goal is a receiver far away, covering an arc of directions in the plane; it
# measures the energy per steradian of the whole flight (in all frequencies or in a band).
RAD_C = 2.0
RAD_T0 = (3.0 - 1.0) * RAD_C * RAD_C
RAD_Q = 1.0 / 40.0
RAD_MAGNETS = [40.0, 80.0, 160.0, 320.0]


def jackson_beaming():
    # Jackson §14.3: a relativistic charge radiates into a cone of ~1/gamma (19 degrees)
    # around its velocity. Straight on (no magnet) it radiates nothing; a bend that sweeps
    # its velocity through 30 degrees lights the receiver, one that bends it down barely
    # does (from 0 degrees, emission at 30 degrees is ~1e4 weaker). 106 single-magnet
    # placements pass 6e-4 per steradian. The reference is a moderate bend (a magnet of
    # 160 at (22, 4): 4.3e-3 per steradian) that does not pass the next level (only
    # 1.4e-4 of it in the band 80-160), so the two show the difference.
    return level(
        "Jackson §14.3: forward beaming",
        "Jackson §14.3: a fast charge radiates when it is accelerated, and almost all of it "
        "goes forward, into a narrow cone around its velocity (half-angle about $1/\\gamma$; here "
        "γ = 3, so about 19°). A receiver far away, at 30° (the band outside the arena), "
        "must collect at least 0.0006 per steradian from the flight. Flying straight the "
        "particle does not radiate at all: bend it so that it heads for the receiver "
        "while it is being bent. It may end anywhere on the right edge. (c = 2 cells per "
        "time unit; radiation reaction is included.)",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 5), 0.0, RAD_T0,
                    box((38, 0, 40, 20), radiation=radiation_goal(30.0, 5.0, (6e-4, 1e3))))],
        reference=[magnet(22, 4, 160.0)],
        max_magnets=3, strengths=RAD_MAGNETS, region=(5, 2, 34, 18))


def jackson_critical_frequency():
    # Jackson §14.6: a bend of radius R flashes the receiver for a time ~R/(gamma^3 c), so
    # its spectrum reaches up to the critical frequency ~ (3/2) gamma^3 c / R (81/R here).
    # The same arena and receiver as forward beaming, but only the band 80-160 counts:
    # of forward beaming's 106 single-magnet solutions, 26 deliver the 3e-3 per steradian
    # asked here (tight bends: e.g. a magnet of 320 at (20, 7), 9.8e-3). Forward beaming's
    # reference (160 at (22, 4)) puts 4.3e-3 into the receiver, but only 1.4e-4 in the
    # band. (A first version mirrored the arena and asked 6e-5 in the band: most of the
    # previous level's bends passed it, so it taught nothing new: the owner's review.)
    return level(
        "Jackson §14.6: the critical frequency",
        "Jackson §14.6: the same receiver as before, but now it only counts angular "
        "frequencies from 80 to 160 (the shaded band of the spectrum in the panel). A "
        "bend's radiation reaches the receiver as a short flash, the shorter the tighter "
        "the bend, and a short flash contains high frequencies: up to about $\\tfrac{3}{2}\\gamma^3 c/R$ "
        "for a bend of radius R. A gentle bend that lit the receiver before may put almost "
        "nothing into this band: bend the particle harder, while it heads for the "
        "receiver, to deliver at least 0.003 per steradian in the band. Watch the spectrum "
        "move up as the bend gets tighter.",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 5), 0.0, RAD_T0,
                    box((38, 0, 40, 20), radiation=radiation_goal(30.0, 5.0, (3e-3, 1e3),
                                                                  band=(80.0, 160.0))))],
        reference=[magnet(20, 7, 320.0)],
        max_magnets=3, strengths=RAD_MAGNETS, region=(5, 2, 34, 18))


def jackson_quiet_turn():
    # Jackson §14.2-14.3: turning from 0 to 90 degrees sweeps the velocity through 45
    # degrees, and the receiver there gets energy per steradian ~ a^2 x (time in the
    # cone) ~ 1/R: gentle turns are quiet. The detector wants the particle heading up
    # (90 +- 30 degrees) and the receiver at 45 degrees at most 6e-3 per steradian. A
    # second way: turn the other way round (clockwise, through -90 and 180 degrees), and
    # the velocity never points at the receiver (the search's one-magnet solutions loop
    # round a magnet: 2.7e-4 per steradian).
    return level(
        "Jackson §14.2: a quiet turn",
        "Jackson §14.2–14.3: bring the particle to the detector at the top, heading up "
        "(within 30°), without lighting the receiver at 45°: it may collect at most 0.006 "
        "per steradian. Any turn from right to up sweeps the velocity through 45°, and "
        "the radiation follows the velocity. The energy that reaches the receiver falls "
        "with the radius of the turn (the power grows with the acceleration squared, the "
        "time in the cone only with the radius). Or is there a way to turn without ever "
        "heading for the receiver?",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 3), 0.0, RAD_T0,
                    box((24, 17, 40, 20), direction=(90.0, 30.0),
                        radiation=radiation_goal(45.0, 5.0, (0.0, 6e-3))))],
        max_magnets=5, strengths=RAD_MAGNETS, region=(5, 2, 36, 16))


def jackson_thomson():
    # Jackson §14.8 with the Doppler shifts of §11.3: a plane wave (omega = 1, travelling
    # towards -x, amplitude 24: a0 = qE0/(m omega c) = 0.3, nearly linear Thomson
    # scattering) shakes the particle, which re-radiates. Heading at phi against the wave,
    # the particle meets it at omega (1 + beta cos phi); a receiver at theta sees
    # omega_s = omega (1 + beta cos phi) / (1 - beta cos(theta - phi)): 32 for phi = theta
    # = 30 degrees, 34 = ~4 gamma^2 head on. The band 27-33 at 30 degrees: the particle must
    # head for the receiver while in the wave. Measured on the reference (80 at (6, 4)):
    # the line peaks at 26-32 (4.5e-3 per steradian in 29-32), with the second harmonic
    # (nonlinear Thomson, a0^2 ~ 0.1) at 50-70. 22 single-magnet solutions reach 3e-3; with
    # the wave off none do (the bends' own flash: only gentle magnets, 40 and 80, are
    # offered; a first try with 8 as amplitude and 3e-4 let 17 bends pass without it).
    return level(
        "Jackson §14.8: Thomson scattering",
        "Jackson §14.8: a light wave (coming from the right, ω = 1) shakes a charge, and the "
        "shaken charge radiates: it scatters the light. A charge flying into the wave meets "
        "its crests faster, and moving towards the receiver it squeezes what it sends "
        r"there: heading at $\varphi$ against the wave, it is seen at an angle $\theta$ to "
        r"scatter $\omega_s = \omega\,\frac{1 + \beta\cos\varphi}{1 - \beta\cos(\theta - \varphi)}$, "
        r"up to $\omega_s \approx 4\gamma^2\omega$ head on: how laser light is turned into "
        "X-rays (inverse Compton sources). The receiver at 30° counts only 27–33: steer the "
        "particle so that the light it scatters towards the receiver lands in that band, at "
        "least 0.003 per steradian.",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 5), 0.0, RAD_T0,
                    box((38, 0, 40, 20), radiation=radiation_goal(30.0, 5.0, (3e-3, 1e3),
                                                                  band=(27.0, 33.0))))],
        reference=[magnet(6, 4, 80.0)],
        disturbances=[stray("light wave", waves=[wave(24.0, 1.0, 0.0, 180.0)])],
        max_magnets=3, strengths=[40.0, 80.0], region=(5, 2, 34, 18))


def jackson_undulator():
    # Jackson §14.7: magnets of alternating sign every s cells wiggle the particle with
    # period lambda_u = 2s; on the axis the wiggles add up coherently at
    # omega_1 = 2 gamma^2 omega_u / (1 + K^2/2), omega_u = 2 pi v / lambda_u: ~107/s here
    # (K << 1). The band [34, 38] wants s = 3 (a line at 36; s = 2 puts it at ~50, s = 4
    # at ~27). Ten magnets of 40 at spacing 3, two cells from the path (half strength at
    # the ends, so that the particle leaves straight) deliver 2.0e-4 per steradian into
    # the band; the same row at spacing 2 or 4, 1e-5 and 3e-6; the line grows as the
    # number of periods squared. The particle must arrive on the axis heading straight
    # (0 +- 3 degrees), as an undulator must not steer the beam: without that, a pair of
    # magnets that steered it close past a magnet near the end made a hard flash of 4.7e-4
    # (the best of 300 random pairs); with it, the best random pair or triple delivers
    # 4e-5 and 6.7e-5. Minimum 1e-4. The search also finds three-magnet solutions (e.g.
    # -40 at (5, 8), 20 at (18, 8), -40 at (32, 7)): kicks far apart whose flashes
    # interfere, fringes ~2 pi / (delay) ~ 8 apart in omega, one in the band. That is an
    # undulator's principle with few periods; the random solve rate is 2e-3.
    return level(
        "Jackson §14.7: undulator",
        "Jackson §14.7: a row of magnets of alternating sign wiggles a fast particle, and "
        "seen from straight ahead the wiggles' radiation adds up at one frequency, "
        "$\\omega_1 = 2\\gamma^2\\omega_u/(1 + K^2/2)$, where $\\omega_u = 2\\pi v/\\lambda_u$ is the wiggle frequency and K the "
        "wiggle strength: the Doppler effect squeezes a slow wiggle into a fast wave. "
        "Build an undulator that delivers at least 0.0001 per steradian to the receiver "
        "straight ahead, in the band 34–38, and lets the particle go on straight along the "
        "axis into the detector. The spacing sets the frequency, the number of periods the "
        "strength of the line. (Half-strength magnets at the ends keep the particle on "
        "course.)",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 10), 0.0, RAD_T0,
                    box((38, 9, 40, 11), direction=(0.0, 3.0),
                        radiation=radiation_goal(0.0, 1.0, (1e-4, 1e3), band=(34.0, 38.0))))],
        reference=[magnet(6 + 3 * k, 8, 40.0 * (-1) ** k * (0.5 if k in (0, 9) else 1.0))
                   for k in range(10)],
        max_magnets=12, strengths=[10.0, 20.0, 40.0], region=(4, 6, 36, 8))



# Arcs: (name, [(tier, [(slug, function)])]). The level files are numbered in this order;
# `levels/curriculum.json` tells the game each level's arc and tier.
ARCS = [
    ("Charges: steering and optics", [
        ("Introduction", [
            ("first_bend", first_bend),
            ("slingshot", slingshot),
            ("geiger_marsden", geiger_marsden),
            ("twin_beams", twin_beams),
            ("two_stages", two_stages),
            ("injection", injection),
        ]),
        ("Intermediate", [
            ("thomson_crt", thomson_crt),
            ("around_the_wall", around_the_wall),
            ("einzel_lens", einzel_lens),
            ("reflectron", reflectron),
            ("collimator", collimator),
            ("hemispherical_analyzer", hemispherical_analyzer),
            ("soft_landing", soft_landing),
        ]),
        ("Master", [
            ("sorting_station", sorting_station),
        ]),
    ]),
    ("Metal and electrodes", [
        ("Introduction", [
            ("high_voltage_dome", high_voltage_dome),
            ("polarized_sphere", polarized_sphere),
            ("image_charge", image_charge),
            ("deflection_plates", deflection_plates),
            ("power_supply", power_supply),
            ("build_a_deflector", build_a_deflector),
            ("shielding", shielding),
        ]),
        ("Intermediate", [
            ("tune_the_lens", tune_the_lens),
            ("real_einzel_lens", real_einzel_lens),
            ("beam_pipe", beam_pipe),
        ]),
        ("Master", [
            ("microscope_column", microscope_column),
        ]),
    ]),
    ("Relativity and magnetism", [
        ("Introduction", [
            ("fast_lane", fast_lane),
            ("first_coil", first_coil),
            ("first_magnet", first_magnet),
            ("stern_gerlach", stern_gerlach),
        ]),
        ("Intermediate", [
            ("beta_spectrometer", beta_spectrometer),
            ("dempster", dempster),
            ("wien_filter", wien_filter),
            ("calutron", calutron),
            ("build_wien_filter", build_wien_filter),
        ]),
        ("Master", [
            ("mass_spectrometer", mass_spectrometer),
        ]),
    ]),
    ("Time: noise, radio frequency and radiation", [
        ("Introduction", [
            ("stray_field", stray_field),
            ("rf_kick", rf_kick),
            ("synchrotron_light", synchrotron_light),
        ]),
        ("Intermediate", [
            ("mains_hum", mains_hum),
            ("earths_field", earths_field),
            ("crt_earth_field", crt_earth_field),
            ("rf_separator", rf_separator),
            ("tune_the_rf", tune_the_rf),
            ("streak_camera", streak_camera),
        ]),
        ("Master", [
            ("rf_beam_line", rf_beam_line),
        ]),
    ]),
    ("Beams", [
        ("Introduction", [
            ("space_charge", space_charge),
            ("stern_gerlach_beam", stern_gerlach_beam),
            ("relativistic_beam", relativistic_beam),
        ]),
        ("Intermediate", [
            ("collimated_beam", collimated_beam),
            ("velocity_selector", velocity_selector),
            ("beam_preparation", beam_preparation),
            ("chromatic_aberration", chromatic_aberration),
            ("real_analyzer", real_analyzer),
            ("calutron_space_charge", calutron_space_charge),
            ("soft_landing_current", soft_landing_current),
        ]),
        ("Master", [
            ("isotope_separator", isotope_separator),
        ]),
    ]),
    ("Jackson: charges in fields and collisions", [
        ("Introduction", [
            ("jackson_exb_drift", jackson_exb_drift),
            ("jackson_van_allen", jackson_van_allen),
            ("throw_a_charge", throw_a_charge),
            ("jackson_faraday", jackson_faraday),
        ]),
        ("Intermediate", [
            ("jackson_knock_on", jackson_knock_on),
            ("jackson_gradient_drift", jackson_gradient_drift),
            ("jackson_runaway", jackson_runaway),
            ("jackson_recoil", jackson_recoil),
            ("jackson_stormer", jackson_stormer),
        ]),
        ("Master", [
            ("jackson_magnetosphere", jackson_magnetosphere),
        ]),
    ]),
    ("Jackson: conductors", [
        ("Introduction", [
            ("jackson_own_image", jackson_own_image),
            ("jackson_two_spheres", jackson_two_spheres),
        ]),
        ("Intermediate", [
            ("jackson_slot", jackson_slot),
            ("jackson_golden_ratio", jackson_golden_ratio),
        ]),
        ("Master", [
            ("jackson_sphere_slalom", jackson_sphere_slalom),
        ]),
    ]),
    ("Jackson: radiation damping", [
        ("Introduction", [
            ("jackson_classical_atom", jackson_classical_atom),
        ]),
        ("Intermediate", [
            ("jackson_circularization", jackson_circularization),
        ]),
        ("Master", [
            ("jackson_three_orbits", jackson_three_orbits),
        ]),
    ]),
    ("Jackson: bound charges", [
        ("Introduction", [
            ("jackson_bound_charge", jackson_bound_charge),
            ("jackson_resonance", jackson_resonance),
        ]),
        ("Intermediate", [
            ("jackson_bound_knock", jackson_bound_knock),
        ]),
        ("Master", [
            ("jackson_spectroscopy", jackson_spectroscopy),
        ]),
    ]),
    ("Jackson: radiation", [
        ("Introduction", [
            ("jackson_beaming", jackson_beaming),
            ("jackson_critical_frequency", jackson_critical_frequency),
        ]),
        ("Intermediate", [
            ("jackson_quiet_turn", jackson_quiet_turn),
            ("jackson_thomson", jackson_thomson),
        ]),
        ("Master", [
            ("jackson_undulator", jackson_undulator),
        ]),
    ]),
]

LEVELS = [lv for _, tiers in ARCS for _, lvls in tiers for lv in lvls]


def write_curriculum():
    """`levels/curriculum.json`: every level's arc and tier, by slug, in order."""
    arcs = [{"name": name, "tiers": [{"name": t, "levels": [slug for slug, _ in lvls]}
                                     for t, lvls in tiers]}
            for name, tiers in ARCS]
    path = os.path.join(ROOT, "levels", "curriculum.json")
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write(json.dumps({"arcs": arcs}, indent=2, ensure_ascii=False) + "\n")



# =======================================================================================
# Machinery.

def run_generator(*args):
    out = subprocess.run([GENERATOR, *args], capture_output=True, text=True, cwd=ROOT,
                         encoding="utf-8")
    if out.returncode != 0:
        sys.exit(f"generator {' '.join(args)} failed:\n{out.stdout}\n{out.stderr}")
    return out.stdout


def landing_points(lvl):
    """Where the reference flights end, from `generator check` on a probe level."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump(lvl, f)
        path = f.name
    try:
        text = run_generator("check", path)
    finally:
        os.unlink(path)
    # One list of (outcome, x, y) per shot: one entry per disturbance.
    ends = {}
    for line in text.splitlines():
        m = re.search(r"shot (\d+)(?: disturbance \d+)?: outcome (\w+).*"
                      r"ends at \(([-\d.]+), ([-\d.]+)\)", line)
        if m:
            ends.setdefault(int(m.group(1)), []).append(
                (m.group(2), float(m.group(3)), float(m.group(4))))
    return [ends[k] for k in sorted(ends)]


def resolve_auto(lvl):
    autos = [s["detector"] for s in lvl["shots"]]
    if not any(isinstance(d, Auto) for d in autos):
        return lvl
    probe = json.loads(json.dumps(lvl, default=lambda d: box(d.strip)))
    if not probe["reference_solution"]:
        # No designed reference: let the solver find a setup that brings every shot into
        # its probe strip, and use it as the reference.
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                         encoding="utf-8") as f:
            json.dump(probe, f)
            path = f.name
        try:
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "600")
            probe = json.load(open(path, encoding="utf-8"))
        finally:
            os.unlink(path)
        if not probe["reference_solution"]:
            sys.exit(f"{lvl['name']}: no reference setup found for the probe strips")
        lvl["reference_solution"] = probe["reference_solution"]
    ends = landing_points(probe)
    placed = []
    for s, d, flights in zip(lvl["shots"], autos, ends):
        if not isinstance(d, Auto):
            continue
        for outcome, _, _ in flights:
            if outcome != "Arrived":
                sys.exit(f"{lvl['name']}: a reference flight misses its probe strip ({outcome})")
        x0, y0, x1, y1 = d.strip
        acrosses = [y if d.axis == "x" else x for _, x, y in flights]
        a_lo, a_hi = min(acrosses), max(acrosses)
        size = d.size
        if len(acrosses) == 1:
            across = a_lo
            lo = math.floor(across - d.size / 2 + 0.5)
            # Keep the landing point at least 0.2 cells inside the detector.
            while across - lo < 0.2:
                lo -= 1
            while lo + d.size - across < 0.2:
                lo += 1
        else:
            # Under disturbances the landing points spread: the detector covers all of
            # them with 0.2 cells to spare, growing beyond `size` if it must.
            while True:
                lo = math.floor(a_lo - 0.2)
                if lo + size - a_hi >= 0.2:
                    break
                size += 1
            # Centre the spread in the detector.
            lo += max(0, math.floor(((lo + size - a_hi) - (a_lo - lo)) / 2))
        # Stay inside the arena.
        n = lvl["grid"]["ny" if d.axis == "x" else "nx"] * lvl["grid"]["subdivision"]
        lo = max(min(lo, n - size), 0)
        across = (a_lo, a_hi)
        # Shots of the same species (particle and energy, differing only in direction)
        # share one detector; different species get detectors of their own.
        screen = (x0, x1) if d.axis == "x" else (y0, y1)
        species = (d.axis, screen, json.dumps(s["particle"]), s["launch"]["kinetic_energy"])
        placed.append([species, [s], across[0], across[1], lo, lo + size])
    groups = {}
    for p in placed:
        g = groups.setdefault(p[0], p)
        if g is not p:
            g[1] += p[1]
            g[2], g[3] = min(g[2], p[2]), max(g[3], p[3])
            g[4], g[5] = min(g[4], p[4]), max(g[5], p[5])
    # Split overlapping detectors on the same screen at the node halfway between the
    # landing points.
    placed = sorted(groups.values(), key=lambda p: (p[0][:2], p[2]))
    for a, b in zip(placed, placed[1:]):
        if a[0][:2] == b[0][:2] and a[5] > b[4]:
            m = round((a[3] + b[2]) / 2)
            if m - a[3] < 0.2 or b[2] - m < 0.2 or m - a[4] < 1 or b[5] - m < 1:
                sys.exit(f"{lvl['name']}: two landing points are too close to separate")
            a[5], b[4] = m, m
    for (axis, (c0, c1), _, _), shots, _, _, lo, hi in placed:
        for s in shots:
            s["detector"] = box((c0, lo, c1, hi) if axis == "x" else (lo, c0, hi, c1))
    return lvl


def previous_build(slug):
    """The level file of `slug` from the previous build (whatever its number), or None."""
    for name in os.listdir(os.path.join(ROOT, "levels")):
        if re.fullmatch(rf"\d+_{re.escape(slug)}\.json", name):
            return os.path.join(ROOT, "levels", name)
    return None


def same_physics(a, b):
    """Whether two levels differ at most in their name, description, limits and
    reference: then a reference of one is a solution of the other (if within limits)."""
    automatic = any(isinstance(sh["detector"], Auto) for sh in a["shots"] + b["shots"])

    def strip(l):
        l = {k: v for k, v in l.items()
             if k not in ("name", "description", "limits", "reference_solution")}
        if automatic:
            # Automatic detectors are placed from the reference: compare without them.
            l["shots"] = [{k: v for k, v in sh.items() if k != "detector"} for sh in l["shots"]]
        return json.dumps(l, sort_keys=True)
    return strip(a) == strip(b)


def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    # Renumber: move each slug's previous file to its place in the curriculum first.
    wanted = {slug: f"{i:02d}_{slug}.json" for i, (slug, _) in enumerate(LEVELS, start=1)}
    for slug, name in wanted.items():
        old = previous_build(slug)
        new = os.path.join(ROOT, "levels", name)
        if old and os.path.normcase(old) != os.path.normcase(new):
            os.replace(old, new)
    for name in os.listdir(os.path.join(ROOT, "levels")):
        m = re.fullmatch(r"\d+_(\w+)\.json", name)
        if m and m.group(1) not in wanted:
            print(f"removing {name} (no longer in the curriculum)")
            os.remove(os.path.join(ROOT, "levels", name))
    write_curriculum()
    for i, (slug, make) in enumerate(LEVELS, start=1):
        key = f"{i:02d}_{slug}"
        if only and not (key.startswith(only) or slug in only.split(",")):
            continue
        path = os.path.join(ROOT, "levels", f"{key}.json")
        previous = json.load(open(path, encoding="utf-8")) if os.path.exists(path) else None
        kept = False
        lvl = make()
        if (previous and previous["reference_solution"] and not lvl["reference_solution"]
                and same_physics(lvl, previous)):
            # Only the text or the limits changed: keep the verified reference (checked
            # again below) instead of searching anew.
            lvl["reference_solution"] = previous["reference_solution"]
            print(f"{key}: reference kept from the previous build")
            kept = True
        lvl = resolve_auto(lvl)
        if slug in MUST_FAIL_ALONE:
            alone = dict(lvl, reference_solution=[])
            with open(path, "w", newline="\n", encoding="utf-8") as f:
                f.write(json.dumps(alone, indent=2, ensure_ascii=False) + "\n")
            report = run_generator("check", path)
            lines = [l for l in report.splitlines() if "shot " in l]
            works = all("Arrived" in l and "Verified" in l for l in lines)
            print(f"{key}: idealised design alone: {'WORKS (not a real iteration)' if works else 'fails, as it should'}")
            for l in lines:
                print("   ", l.strip())
            if works:
                sys.exit(f"{key}: the idealised design still works on its own")
        with open(path, "w", newline="\n", encoding="utf-8") as f:
            f.write(json.dumps(lvl, indent=2, ensure_ascii=False) + "\n")
        indirect = check_indirect(lvl)
        if indirect:
            sys.exit(f"{key}: goal particles are reached directly: " + "; ".join(indirect))
        if slug in FINALES:
            fewer = FINALES[slug]
            ok_fewer = needs_elements(lvl, fewer)
            print(f"{key}: not solved with {fewer} charges: {ok_fewer}", flush=True)
            if not ok_fewer:
                sys.exit(f"{key}: solvable with only {fewer} charges")
        if not lvl["reference_solution"]:
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "400")
        report = run_generator("check", path)
        if "placement Ok" not in report and kept:
            # A kept reference outside the new limits (e.g. a smaller region): solve anew.
            print(f"{key}: kept reference no longer allowed; solving anew", flush=True)
            lvl["reference_solution"] = []
            with open(path, "w", newline="\n", encoding="utf-8") as f:
                f.write(json.dumps(lvl, indent=2, ensure_ascii=False) + "\n")
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "400")
            report = run_generator("check", path)
        ok = "placement Ok" in report and all(
            "Arrived" in l and "Verified" in l for l in report.splitlines() if "shot " in l)
        print(f"{key}: reference {'verified' if ok else 'NOT VERIFIED'}", flush=True)
        if not ok:
            print(report)

    render_math()


def render_math():
    """Renders the LaTeX formulas ($...$) of every level's description (shipped and
    custom) to SVG with MathJax (scripts/math/render.mjs; `npm install` in scripts/math
    once), into levels/math/ with a manifest the game reads (crates/game/src/math.rs)."""
    formulas = []
    for folder in (os.path.join(ROOT, "levels"), os.path.join(ROOT, "levels", "custom")):
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            if not name.endswith(".json") or name in ("curriculum.json", "golden_hashes.json"):
                continue
            text = json.load(open(os.path.join(folder, name), encoding="utf-8")).get(
                "description", "")
            if text.count("$") % 2:
                sys.exit(f"{name}: unmatched $ in the description")
            formulas += re.findall(r"\$([^$]+)\$", text)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                     encoding="utf-8") as f:
        json.dump(sorted(set(formulas)), f)
        path = f.name
    try:
        out = subprocess.run(["node", os.path.join(ROOT, "scripts", "math", "render.mjs"), path,
                              os.path.join(ROOT, "levels", "math")],
                             capture_output=True, text=True, encoding="utf-8")
    finally:
        os.unlink(path)
    if out.returncode != 0:
        sys.exit("rendering the formulas failed: " + out.stdout + out.stderr)
    print(out.stdout.strip())
    # Formulas are not broken across lines (MathJax 3's SVG output has no line breaking);
    # the game shrinks one wider than the panel to fit. Warn well before that: 30 ex is
    # about 60 % of the side panel (340 points, 1 ex = 0.52 x the 12.5-point body font).
    manifest = json.load(open(os.path.join(ROOT, "levels", "math", "manifest.json"),
                              encoding="utf-8"))
    for tex, m in sorted(manifest.items()):
        if m["width_ex"] > MAX_FORMULA_EX:
            print(f"warning: formula {m['width_ex']:.1f} ex wide (over {MAX_FORMULA_EX}): "
                  f"${tex}$ -- split it into several $...$ at its relation signs, so that it "
                  "can wrap between them")


MAX_FORMULA_EX = 30.0


if __name__ == "__main__":
    main()
