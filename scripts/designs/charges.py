"""Arcs: charges, several shots with one setup, relativity."""

import math
from designs.common import Auto, M, box, charge, level, shot

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
