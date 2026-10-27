"""Arcs: magnetic fields, magnetic moments, noise (fields from outside), radio frequency, radiation."""

import math
from designs.common import (
    Auto,
    M,
    box,
    charge,
    circle_coil,
    level,
    magnet,
    rect_coil,
    shot,
    stray,
    wave,
)

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
