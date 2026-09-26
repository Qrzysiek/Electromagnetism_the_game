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
GENERATOR = os.path.join(ROOT, "target", "release", "generator")
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


def box(b):
    return {"min": [b[0], b[1], 0], "max": [b[2], b[3], 0]}


def shot(q, m, node, angle_deg, ke, detector, time=0.0):
    a = math.radians(angle_deg)
    launch = {"node": [node[0], node[1], 0], "direction": [math.cos(a), math.sin(a), 0.0],
              "kinetic_energy": ke}
    if time:
        launch["time"] = time
    return {
        "particle": {"charge": q, "mass": m, "radius": 0.0},
        "launch": launch,
        "detector": detector,
    }


def metal(x, y, r, kind="grounded", value=None):
    """A metal sphere: kind "grounded", "charge" (net charge value) or "potential"."""
    bias = {"kind": kind}
    if value is not None:
        bias["value"] = value
    return {"center": [x, y, 0], "radius": r, "bias": bias}


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
          amplitudes=(), rf_omega=0.0, radiation_reaction=False, omegas=(), conductors=()):
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
    return {
        "format_version": 2, "engine_version": "0.1.0", "name": name, "description": desc,
        "grid": {"nx": grid[0], "ny": grid[1], "nz": 0, "subdivision": 1},
        "physics": {"c": c, "charge_radius": 0.3, "magnet_radius": 0.3, "wire_radius": 0.1,
                    "antenna_radius": 0.3,
                    "t_max": t_max, "tolerances": {"preview": 1e-10, "verify": 1e-12},
                    **({"rf_omega": rf_omega} if rf_omega else {}),
                    **({"radiation_reaction": True} if radiation_reaction else {})},
        "shots": list(shots), "elements": list(elements), "coils": list(coils),
        "limits": limits, "reference_solution": list(reference or []),
        **({"disturbances": list(disturbances)} if disturbances else {}),
        **({"conductors": list(conductors)} if conductors else {}),
    }


# =======================================================================================
# Chapter 1: charges.

def first_bend():
    return level(
        "First bend",
        "Place one charge to bend the beam into the detector. Like charges repel, unlike "
        "charges attract.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 16, 30, 20)))],
        max_charges=1, magnitudes=[1 * M, 2 * M, 4 * M])


def slingshot():
    return level(
        "Slingshot",
        "The detector is behind you. Attraction can turn the particle around; the fixed "
        "negative charge is a hint, but not enough on its own.",
        shots=[shot(1e-6, 1.0, (6, 3), 90.0, 0.5, box((0, 0, 3, 4)))],
        elements=[charge(15, 12, -2 * M)],
        max_charges=1, magnitudes=[1 * M, 2 * M, 4 * M])


def geiger_marsden():
    return level(
        "Geiger–Marsden (1909)",
        "Alpha particles fired at gold foil sometimes bounced straight back: the atom has a "
        "tiny, heavy, positive nucleus. Steer the alpha near the source so that it "
        "backscatters into the detector. The deflection follows tan(θ/2) = qQ / (2 T₀ b): "
        "the impact parameter b decides everything.",
        shots=[shot(2e-6, 4.0, (0, 11), 0.0, 4.0, box((0, 18, 1, 20)))],
        elements=[charge(18, 10, 8 * M)],
        max_charges=1,
        magnitudes=[m * M for m in (0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5)],
        region=(2, 2, 10, 18), c=20.0)


def the_wall():
    return level(
        "The wall",
        "A wall of positive charges blocks the direct path. Several charges together can "
        "guide the particle around it.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 8, 30, 12)))],
        elements=[charge(15, y, 0.5 * M) for y in range(7, 14)],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M])


def thomson_crt():
    return level(
        "Thomson's cathode-ray tube (1897)",
        "J. J. Thomson deflected cathode rays with charged plates and showed they are "
        "negative particles: electrons. Here an electron at 0.2 c must hit a one-cell spot "
        "on the screen. Two electron energies must hit their own one-cell spots: the "
        "deflection falls as 1/T, like in an oscilloscope. Build deflection plates from "
        "point charges.",
        shots=[shot(-1e-6, 1.0, (0, 10), 0.0, e, Auto((29, 0, 30, 20), "x")) for e in (0.4, 0.6)],
        max_charges=2, magnitudes=[m * M for m in (0.25, 0.5, 0.75, 1, 1.5)],
        region=(8, 3, 17, 17),
        reference=[charge(8, 4, -1 * M)])


# =======================================================================================
# Chapter 2: several shots, one setup.

def twin_beams():
    return level(
        "Twin beams",
        "Two particles, two detectors, one setup: every shot must reach its own detector. "
        "Switch between the shots with [ and ], or show them all.",
        shots=[shot(1e-6, 1.0, (0, 13), 0.0, 0.5, box((27, 15, 30, 19))),
               shot(1e-6, 1.0, (0, 7), 0.0, 0.5, box((27, 1, 30, 5)))],
        max_charges=2, magnitudes=[1 * M, 2 * M, 4 * M])


def hemispherical_analyzer():
    return level(
        "Hemispherical analyzer (XPS)",
        "Photoelectron spectrometers sort particles by energy between two concentric "
        "hemispheres, where the field is exactly that of a point charge at the centre. "
        "Three energies leave the source; each must reach its own detector after bending "
        "through 180°. A circular orbit of radius R needs |qQ| = p·v·R; the others become "
        "ellipses.",
        shots=[shot(1e-6, 1.0, (4, 17), 0.0, e, Auto((1, 0, 3, 11), "x"))
               for e in (0.45, 0.5, 0.55)],
        max_charges=2, magnitudes=[m * M for m in (1, 2, 3, 4, 5, 6, 7, 8)],
        region=(1, 3, 12, 16),
        reference=[charge(4, 10, -7 * M)])


def reflectron():
    return level(
        "Reflectron",
        "Time-of-flight mass spectrometers turn ions around with an electrostatic mirror, "
        "the reflectron. Two ions with different energies must both come back to the "
        "detector beside the source.",
        shots=[shot(1e-6, 1.0, (0, 13), 0.0, e, box((0, 4, 2, 7))) for e in (0.4, 0.6)],
        max_charges=3, magnitudes=[1 * M, 2 * M, 4 * M], signs=(True, False),
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
# Chapter 3: relativity.

def fast_lane():
    return level(
        "Fast lane",
        "A relativistic particle (v = 0.85 c) is hard to bend: its momentum grows as γmv. "
        "Charges that would easily steer a slow particle barely move this one.",
        # q = 1e-7 with 10x larger fixed charges: the same trajectories as q = 1e-6, but
        # 100x less (neglected) radiation, which scales as q²/(m c² r) at fixed qQ.
        shots=[shot(1e-7, 1.0, (0, 4), math.degrees(math.atan2(0.2, 1.0)), 2.0,
                    box((27, 14, 30, 18)))],
        elements=[charge(12, 8, -20 * M), charge(18, 12, 20 * M), charge(20, 5, -20 * M)],
        max_charges=2, magnitudes=[m * M for m in (20, 30, 40, 60, 80)], c=1.5)


def beta_spectrometer():
    return level(
        "Beta-ray spectrometer",
        "Beta decay electrons are fast. Two electrons, at v = 0.91 c (γ = 2.41) and "
        "v = 0.87 c (γ = 2.0), must each reach their own detector after bending around "
        "one charge. The circular-orbit condition is |qQ| = γ m v² R, not the Newtonian "
        "2 T₀ R. Trust relativity, not intuition. (The electron charge is small and the "
        "fixed charges large, so that radiation stays negligible.)",
        shots=[shot(-1e-7, 1.0, (4, 17), 0.0, e, Auto((1, 0, 3, 13), "x"))
               for e in (math.sqrt(2.0), 1.0)],
        max_charges=2, magnitudes=[m * M for m in (20, 40, 70, 100, 120, 140, 170, 200)],
        region=(1, 3, 12, 16), c=1.0,
        reference=[charge(4, 10, 140 * M)])


# =======================================================================================
# Chapter 4: magnetic fields (coils placed by the level).

def first_coil():
    return level(
        "First coil",
        "The coil's magnetic field bends moving charges into circles: the force q v × B is "
        "perpendicular to the velocity, so it changes the direction but never the speed "
        "(watch the energy bars). Use one charge to steer the circling particle into the "
        "detector.",
        shots=[shot(1e-6, 1.0, (8, 10), 90.0, 0.5, box((14, 8, 17, 11)))],
        coils=[circle_coil(15, 10, 9.5, 1.9e5)],
        max_charges=1, magnitudes=[0.5 * M, 1 * M, 2 * M])


def dempster():
    shots = [shot(1e-6, m, (9, 5), a, 0.05, Auto((11, 2, 31, 5), "y"))
             for m in (1.0, 2.0, 4.0) for a in (85.0, 95.0)]
    return level(
        "Dempster's mass spectrometer (1918)",
        "Ions from the source are accelerated by your electrodes, then bent by the magnetic "
        "field of the large coil. In a uniform field r = √(2mT)/(qB), so masses 1, 2 and 4 "
        "land at different places; ions leaving at slightly different angles re-converge "
        "after half a turn (180° focusing). One setup must bring all six ions to the "
        "detectors of their masses.",
        grid=(32, 22), shots=shots, coils=[circle_coil(16, 11, 10.5, 5e5)],
        max_charges=2,
        magnitudes=[m * M for m in (0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4)],
        region=(2, 1, 15, 4),
        reference=[charge(9, 3, 2 * M)])


def wien_filter():
    return level(
        "Wien filter",
        "Crossed electric and magnetic fields pass exactly one speed straight through: qE "
        "balances qvB when v = E/B. Slower ions are pushed one way, faster ones the other. "
        "The coil provides B; place charges to supply E so that each of the three ions "
        "reaches its own detector.",
        shots=[shot(1e-6, 1.0, (2, 10), 0.0, e, Auto((26, 6, 28, 14), "x"))
               for e in (0.3, 0.5, 1.0)],
        coils=[rect_coil(1, 5, 29, 15, 2.5e4)],
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.15, 0.2, 0.3, 0.4)],
        region=(6, 6, 24, 14),
        reference=[charge(11, 7, 0.2 * M), charge(19, 7, 0.2 * M),
                   charge(11, 13, -0.2 * M), charge(19, 13, -0.2 * M)])


# =======================================================================================
# Chapter 5: your own magnets.

def first_magnet():
    return level(
        "First magnet",
        "A magnet here is a uniformly magnetized sphere standing out of the plane (⊙) or "
        "into it (⊗). In the plane its field is perpendicular to the plane and falls as "
        "1/r³, so it bends the particle sideways, strongly only nearby. Place one magnet to "
        "bend the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        max_magnets=1, strengths=[1 * M, 2 * M, 4 * M], region=(6, 2, 24, 18))


def calutron():
    return level(
        "Calutron",
        "Lawrence's calutron separated uranium isotopes by bending ions in a magnetic "
        "field: at equal energy the radius r = √(2mT)/(qB) grows with the mass. Three "
        "isotopes leave the source together; place magnets so that each reaches its own "
        "detector.",
        # Separated probe strips: in the same field the lightest isotope bends most.
        shots=[shot(1e-6, m, (3, 10), 0.0, 0.5, Auto((27, y0, 30, y1), "x", 3))
               for m, (y0, y1) in ((1.0, (15, 19)), (2.0, (10, 15)), (4.0, (4, 10)))],
        max_magnets=2, strengths=[m * M for m in (1, 2, 3, 4)], region=(6, 2, 24, 18),
        reference=None)


def build_wien_filter():
    return level(
        "Build a Wien filter",
        "No coil this time: build the whole velocity selector yourself from charges (for "
        "E) and magnets (for B). Three ions with different energies must each reach their "
        "own detector.",
        shots=[shot(1e-6, 1.0, (2, 10), 0.0, e, Auto((26, 2, 28, 18), "x"))
               for e in (0.3, 0.5, 1.0)],
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.4)],
        max_magnets=2, strengths=[m * M for m in (1, 2, 3)], region=(6, 3, 24, 17),
        reference=[magnet(15, 6, -2 * M), magnet(15, 14, -2 * M),
                   charge(11, 7, 0.2 * M), charge(19, 7, 0.2 * M),
                   charge(11, 13, -0.2 * M), charge(19, 13, -0.2 * M)])


# =======================================================================================
# Chapter: noise. Fields from outside the arena; one setup must work under each of them.

def stray_field():
    # A 0.5 T0 particle crosses the arena in t ~ 30; a stray E = 6.7e3 gives a = qE/m =
    # 6.7e-3 and pushes it ~3 cells off course.
    return level(
        "Stray field",
        "Somewhere in the building a high-voltage supply is switched on and off, and with "
        "it a weak uniform field across your beam line. The beam must reach the detector "
        "in both cases. Compare the two flights and aim between them.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 12, 30, 16)))],
        max_charges=1, magnitudes=[m * M for m in (0.25, 0.5, 1, 2, 4)],
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
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)],
        region=(4, 2, 24, 18),
        disturbances=[stray("no field"), stray("facing north", bz=b),
                      stray("facing south", bz=-b)],
        c=None)


# =======================================================================================
# Chapter 7: radio frequency. Antennas: oscillating dipoles with their exact retarded
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
        "passes. Orient and place one antenna to steer the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        max_antennas=1, amplitudes=[m * M for m in (1, 2, 4, 8)], rf_omega=RF)


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
# Chapter 8: radiation. The particle's own radiation is part of the physics
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
        "phase difference is ω Δt. Pick the frequency, then the orientation and place.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, d, 30, d + 3)), time=t)
               for t, d in ((0.0, 14), (4.0, 3))],
        max_antennas=1, amplitudes=[m * M for m in (2, 4, 8)],
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
        "beam; add one charge to bring it into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 16, 30, 19)))],
        conductors=[metal(15, 4, 2.5, "potential", 3e5)],
        max_charges=1, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)])


def image_charge():
    return level(
        "Image charge",
        "A large grounded sphere. Whatever charge you place near it, the metal answers with "
        "an opposite image charge inside: near the surface your charge is almost cancelled, "
        "far away the pair acts like a dipole. Steer the beam around the sphere.",
        shots=[shot(1e-6, 1.0, (0, 6), 0.0, 0.5, box((27, 13, 30, 17)))],
        conductors=[metal(14, 10, 3.5, "grounded")],
        max_charges=2, magnitudes=[m * M for m in (0.5, 1, 2, 4)])


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
        max_charges=2, magnitudes=[m * M for m in (0.5, 1, 2, 4)], region=(4, 2, 24, 18))


LEVELS = [
    # Chapter 1: charges (intro, then rising difficulty).
    ("first_bend", first_bend),
    ("geiger_marsden", geiger_marsden),
    ("thomson_crt", thomson_crt),
    ("slingshot", slingshot),
    ("the_wall", the_wall),
    # Chapter 2: several shots, one setup.
    ("twin_beams", twin_beams),
    ("reflectron", reflectron),
    ("einzel_lens", einzel_lens),
    ("hemispherical_analyzer", hemispherical_analyzer),
    # Chapter 3: metals (induced charge).
    ("high_voltage_dome", high_voltage_dome),
    ("polarized_sphere", polarized_sphere),
    ("image_charge", image_charge),
    # Chapter 4: relativity.
    ("fast_lane", fast_lane),
    ("beta_spectrometer", beta_spectrometer),
    # Chapter 5: magnetic fields (coils placed by the level).
    ("first_coil", first_coil),
    ("dempster", dempster),
    ("wien_filter", wien_filter),
    # Chapter 6: your own magnets.
    ("first_magnet", first_magnet),
    ("calutron", calutron),
    ("build_wien_filter", build_wien_filter),
    # Chapter 7: noise (outside fields; one setup for every disturbance).
    ("stray_field", stray_field),
    ("mains_hum", mains_hum),
    ("earths_field", earths_field),
    # Chapter 8: radio frequency.
    ("rf_kick", rf_kick),
    ("rf_separator", rf_separator),
    ("tune_the_rf", tune_the_rf),
    ("streak_camera", streak_camera),
    # Chapter 9: radiation.
    ("synchrotron_light", synchrotron_light),
]


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


def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    for i, (slug, make) in enumerate(LEVELS, start=1):
        key = f"{i:02d}_{slug}"
        if only and not key.startswith(only):
            continue
        lvl = resolve_auto(make())
        path = os.path.join(ROOT, "levels", f"{key}.json")
        with open(path, "w", newline="\n", encoding="utf-8") as f:
            f.write(json.dumps(lvl, indent=2, ensure_ascii=False) + "\n")
        if not lvl["reference_solution"]:
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "400")
        report = run_generator("check", path)
        ok = all("Arrived" in l and "Verified" in l for l in report.splitlines() if "shot " in l)
        print(f"{key}: reference {'verified' if ok else 'NOT VERIFIED'}", flush=True)
        if not ok:
            print(report)


if __name__ == "__main__":
    main()
