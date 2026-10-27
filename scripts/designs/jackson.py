"""The Jackson series (docs/JACKSON.md): effects from J. D. Jackson, Classical Electrodynamics."""

import json
import math
import os
import re
import sys
import tempfile
from designs.common import (
    Auto,
    M,
    box,
    charge,
    circle_coil,
    cloud,
    free_charge,
    level,
    magnet,
    metal_sphere,
    plate,
    radiation_goal,
    run_generator,
    shot,
    slider,
    stray,
    supply,
    wave,
)
from designs.instruments import solve_stage

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


def ramp_the_coil(rate=-5.3e3, rates=(-1e4, -5e3, 5e3, 1e4), vmax=2.1e4, detector=None):
    # Coil supplies' introduction (the player sets a coil's ramp rate), as "Power supply"
    # introduced the electrodes' supplies: a coil (radius 6) whose current rises from 0
    # at the rate the player sets, and an ion crossing it from the left (T0 = 0.5, v = 1:
    # inside the coil from t = 9 to 21). It meets the field the ramp has built by then
    # (and the induced field -rate A), so the rate sets how far it is bent; the detector
    # is placed around the reference's landing, and `generator window` keeps the slider's
    # guessable points out of the solving interval. The supply sits on the coil's
    # terminal (19, 5), its wire's node due east of the centre (`Coil::supply_node`).
    coil = circle_coil(15, 5, 4.0, 0.0)
    coil["tunable"] = True
    detector = detector or Auto((28, 0, 30, 20), "x", 2)
    return level(
        "Ramp the coil",
        "The coil now has a power supply too (the yellow ring: click the coil, then set "
        "the slider). It does not set a fixed current: it sets how fast the current "
        "rises from zero. The ion crosses the coil a while after the start and meets "
        "the field the current has built by then: the faster the ramp, the more it is "
        "bent. (A rising field also induces an electric field around the axis, "
        "Faraday's law, which pushes the ion along a little.) Find the ramp that brings "
        "the ion into the detector.",
        shots=[shot(1e-6, 1.0, (0, 12), 0.0, 0.5, detector)],
        coils=[coil], c=None, t_max=100.0,
        coil_rates=slider(vmax, *rates, rate),
        reference=[supply(19, 5, rate)])



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



# Jackson arc on electrostatics (Ch. 2-4). The golden-ratio levels use a particle with charge
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
        conductors=[metal_sphere(15, 11, 3, "grounded")],
        max_charges=3, magnitudes=[0.25, 0.5, 1.0, 2.0], c=None, t_max=300.0)



def jackson_two_spheres():
    return level(
        "Jackson Pr. 2.6: two spheres",
        "Jackson Problem 2.6: two conducting spheres, one charged, one neutral. Each sphere's "
        "charge induces an image in the other, which induces an image in the first, and so "
        "on: the neutral sphere becomes a dipole facing the charged one, and both respond "
        "to your charges and to the particle itself. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        conductors=[metal_sphere(11, 5, 2.5, "charge", 3e6), metal_sphere(19, 13, 2.5, "charge", 0.0)],
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
        conductors=[metal_sphere(15, 10, 3, "charge", 1.0)],
        max_charges=4, magnitudes=[0.25, 0.5, 1.0, 2.0], c=None, t_max=300.0)



def jackson_multipoles():
    # Jackson §4.1 (the multipole expansion). Seen from a beam passing at distance b, a
    # group of charges kicks it, to first order, by (2q/(m v^2)) sum_i Q_i/(b - y_i) (y_i
    # across the beam): ~Q/b for a net charge, ~p/b^2 for a dipole, ~Qd^2/b^3 for a
    # quadrupole. Two far beams (A along x, 11-16 cells above the region; C along y, 15-21
    # cells to its right) must arrive within 0.05 degrees of straight; the near beam B, one
    # cell below the region, must be bent down about 10 degrees into its detector. A net
    # charge or a dipole strong enough to bend B turns A by 0.2 degrees or more; a compact
    # quadrupole (+q, -2q, +q, 1 cell apart) by 0.01-0.03. C catches what A cannot see (a
    # structure along x). The region is centred on both far flights, where the part of a
    # charge's field that a flight sees is stationary: a first version with the region near
    # the far beams' launch points let the search hide a net charge of -0.25 M by putting
    # the positive charges further along the flight.
    return level(
        "Jackson §4.1: multipoles",
        "Jackson §4.1 (the multipole expansion): far from a group of charges, only its "
        "lowest nonvanishing moment counts. A beam passing at a distance b is turned by a "
        "net charge Q by an angle $\\propto Q/b$, by a dipole p (charges ±q a distance d "
        "apart, p = qd) only $\\propto p/b^2$, by a quadrupole (such as +q, −2q, +q in a "
        "row) only $\\propto qd^2/b^3$. Bend the lower beam down into its detector, but "
        "leave the two far beams untouched: they must arrive within 0.05° of straight.",
        grid=(40, 30),
        shots=[shot(1e-6, 1.0, (0, 28), 0.0, 0.5, box((38, 26, 40, 30), direction=(0.0, 0.05))),
               shot(1e-6, 1.0, (37, 0), 90.0, 0.5, box((36, 29, 38, 30), direction=(90.0, 0.05))),
               shot(1e-6, 1.0, (0, 11), 0.0, 0.5, Auto((38, 0, 40, 30), "x", 2))],
        reference=[charge(19, 12, 0.25 * M), charge(19, 13, -0.5 * M),
                   charge(19, 14, 0.25 * M)],
        max_charges=5, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], region=(16, 12, 22, 17))



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
        conductors=[metal_sphere(10, 5, 2.5, "charge", 1.0), metal_sphere(21, 15, 2.5, "charge", 1.0),
                    metal_sphere(31, 5, 2.5, "charge", 1.0)],
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


# Levels with insulating bodies (dielectrics, ferrites): their reference must fail
# without them (they matter).
NEEDS_BODY = {"glass_screen", "ferrite_shield"}


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



def bound_scattered_power(q, m, w0sq, w, c, e0, axis_deg, half_deg, pol_deg):
    """Mean power per steradian a charge bound with omega_0^2 = w0sq scatters from a plane
    wave (omega w, amplitude e0, polarization at pol_deg), averaged over the arc axis ±
    half: the Landau-Lifshitz steady state (test S10), chi = (q/m)(1 - i w tau)/(w0^2 -
    w^2 - i w Gamma), Gamma = w0^2 tau, radiating (q^2/8 pi c^3) w^4 |chi|^2 e0^2 sin^2 of
    the angle from the polarization."""
    tau = 2.0 * q * q / (3.0 * m * c ** 3)
    num = abs(complex(1.0, -w * tau)) ** 2
    den = abs(complex(w0sq - w * w, -w * w0sq * tau)) ** 2
    chi2 = (q / m) ** 2 * num / den
    a, h = math.radians(axis_deg - pol_deg), math.radians(half_deg)
    # mean of sin^2(theta) over a +- h
    mean = 0.5 - math.cos(2 * a) * (math.sin(2 * h) / (4 * h) if h > 0 else 0.5)
    return q * q / (8 * math.pi * c ** 3) * w ** 4 * chi2 * e0 * e0 * mean



def jackson_cross_section():
    # Jackson §16.8, the cross section of a bound charge, with the steady receiver. A
    # Thomson atom (Q = 16, R = 4: omega_0 = 1/2) and its electron (q = -1, m = 1), c = 2
    # (tau = 1/12, Gamma = omega_0^2 tau = 0.021: a transient decays as e^{-Gamma t/2}),
    # lit from the left by a plane wave (omega = 0.35, polarized along y, amplitude 0.01:
    # the electron swings ~0.1). The receiver behind the light (180 +- 10 degrees: the
    # scattered light only, the wave goes the other way; in the plane the dipole radiates
    # the most along x) counts the band 0.32-0.38 over t = 400-1100 (the start's transient,
    # at the atom's own frequency, outside the band and decayed). The electron starts with a
    # small swing along x, which sends nothing back along the axis. Bare atom: 0.92 sigma_T
    # (Rayleigh's side). The goal, 3-6 times that, wants the spring along y softened: a
    # pair of positive charges at +-d along y gives omega_y^2 = omega_0^2 - 4Q/d^3 (no shift
    # of the electron): (2, 5), (4, 6), (8, 8), (16, 10), and (16, 7) past the resonance;
    # single charges also pull the electron aside.
    q, m, big_q, r, c = -1.0, 1.0, 16.0, 4.0, 2.0
    w0 = cloud_omega(q, big_q, r)
    w, e0 = 0.35, 0.01
    bare = bound_scattered_power(q, m, w0 * w0, w, c, e0, 180.0, 10.0, 90.0)
    return level(
        "Jackson §16.8: why the sky is blue",
        "Jackson §16.8: light scattered by a bound electron. The light (ω = 0.35, from the "
        "left, its field along y) shakes the atom's electron, and the shaken electron sends "
        "light back to the receiver. How much depends on how near the light's frequency is "
        "to the atom's own: "
        r"$\sigma = \sigma_T\,\frac{\omega^4}{(\omega_0^2 - \omega^2)^2 + (\omega\Gamma)^2}$. "
        "Far below the resonance it grows as $\\omega^4$ (Rayleigh: blue light, nearer the "
        "ultraviolet resonances of air, scatters more than red, so the sky is blue); at the "
        "resonance it is enormous; far above it the electron scatters like a free one "
        "(Thomson). This atom rings at ω₀ = 0.5, too stiff for this light: soften its spring "
        "along the light's field until the receiver counts 3 to 6 times what the bare atom "
        "sends. A charge pulls the swinging electron the harder the nearer it comes: on the "
        "line of the swing it weakens the spring, across it it stiffens it. Mind its steady "
        "pull too. The receiver averages over a window, after the start has died away; the "
        "electron must stay in its atom.",
        grid=(40, 20), c=c, t_max=1200.0, radiation_reaction=True,
        shots=[shot(q, m, (26, 10), 0.0, 0.5 * w0 ** 2 * 0.2 ** 2,
                    box((0, 8, 1, 12), radiation=radiation_goal(
                        180.0, 10.0, (3.0 * bare, 6.0 * bare), band=(0.32, 0.38),
                        window=(400.0, 1100.0))))],
        clouds=[cloud(26, 10, r, big_q)],
        disturbances=[stray("light", waves=[wave(e0, w, 0.0, 0.0)])],
        reference=[charge(26, 2, 8.0), charge(26, 18, 8.0)],
        max_charges=4, magnitudes=[1.0, 2.0, 4.0, 8.0, 16.0], region=(16, 0, 36, 20))



def jackson_braking():
    # Jackson §15.2 (radiation in a collision; the sudden stop): the target stops the
    # particle at once (the goal's abrupt_stop), which radiates the flat spectrum
    # (q^2/4 pi^2 c) beta^2 sin^2 theta / (1 - beta cos theta)^2 around its final velocity:
    # nothing straight ahead, the most at cos theta = beta (19.5 degrees here, gamma^2
    # times the broadside value). The receiver is straight ahead (0 +- 5 degrees), the band
    # 100-200 (where the gentle bends offered, 40 and 80, put almost nothing: their
    # critical frequency is ~81/R): at 6.3e-3 per steradian for the best heading, the
    # minimum 5e-3 wants the particle to strike the target ~12-32 degrees off the axis.
    # Flying straight in delivers nothing; 80 single-magnet solutions; without the stop's
    # radiation (the control) no one- or two-magnet placement reaches even 3e-3. The
    # reference (80 at (28, 9)) delivers 6.06e-3, near the analytic best.
    return level(
        "Jackson §15.2: braking radiation",
        "Jackson §15.2: a charge that is stopped abruptly radiates a flash of every "
        "frequency (bremsstrahlung, braking radiation). This target stops the particle "
        "at once. The flash goes mostly forward, but not straight ahead: seen at an angle "
        r"$\theta$ to the particle's last velocity it is "
        r"$\propto \sin^2\theta/(1 - \beta\cos\theta)^2$, zero along the velocity and "
        r"largest at $\cos\theta = \beta$ (here about 20°). The receiver straight ahead "
        "counts the band 100–200 and wants at least 0.005 per steradian: fired straight "
        "into the target the particle lights it not at all.",
        grid=(40, 20), c=RAD_C, t_max=100.0, radiation_reaction=True,
        shots=[shot(RAD_Q, 1.0, (2, 10), 0.0, RAD_T0,
                    box((38, 4, 40, 16),
                        radiation=dict(radiation_goal(0.0, 5.0, (5e-3, 1e3), band=(100.0, 200.0)),
                                       abrupt_stop=True)))],
        reference=[magnet(28, 9, 80.0)],
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




# --- Jackson Pr. 14.23: coherent radiation of several charges --------------------------
# Crossed fields make a charge circle (radius 2, speed 1.6 = 0.4 c in the frame drifting
# with E x B / B^2 = 0.2 along x; omega_0 = qB / (gamma m) = 0.8) while its circle drifts
# into the detector (about 11 turns). The receiver at 90 degrees sees every particle's
# field: a radiation goal in a level with free charges measures them all, their far
# fields added (PHYSICS.md §3.4). Measured (the charge alone: 4.84e-5 per sr in the band
# 0.6-1.0, 8.94e-5 in 0.6-2.6): a second charge in step one cell below: 3.96 times; an
# antiphase partner on the far side at the exact (relativistic) mirror speed 1.738: 0.0035
# times (at the Galilean 1.8: 0.107; at 1.7: 0.045); four charges evenly spaced: 0.0036
# times in 0.6-2.6 (only 4 omega_0 = 3.2 survives, at 15.1 times), 0.083 with Galilean
# velocities, while the pair gives 1.24 there (its 2 omega_0 line, 3.85 times).
RING_C = 4.0

RING_VG, RING_VD, RING_RHO = 1.6, 0.2, 2.0



def ring_gamma(v):
    return 1.0 / math.sqrt(1.0 - (v / RING_C) ** 2)



RING_B = ring_gamma(RING_VG) * RING_VG / (RAD_Q * RING_RHO)

RING_E = RING_VD * RING_B

RING_T0 = (ring_gamma(RING_VG - RING_VD) - 1.0) * RING_C * RING_C

# The shot's speed in the drifting frame (it starts at the bottom of its circle, moving
# along -x in the lab at 1.4).
RING_U = (RING_VG - RING_VD + RING_VD) / (1.0 + (RING_VG - RING_VD) * RING_VD / RING_C ** 2)



def ring_partner(x, y, ux, uy):
    """A free charge at (x, y) circling with velocity (ux, uy) in the drifting frame: the
    drift added relativistically."""
    den = 1.0 + ux * RING_VD / RING_C ** 2
    vx = (ux + RING_VD) / den
    vy = uy / (ring_gamma(RING_VD) * den)
    return free_charge(x, y, RAD_Q, math.degrees(math.atan2(vy, vx)), math.hypot(vx, vy))



def ring_level(name, desc, band, energy, reference, max_free):
    return level(
        name, desc,
        grid=(30, 20), c=RING_C, t_max=200.0, radiation_reaction=True, beam_interaction=True,
        shots=[shot(RAD_Q, 1.0, (6, 8), 180.0, RING_T0,
                    box((25, 5, 27, 15),
                        radiation=radiation_goal(90.0, 10.0, energy, band=band)))],
        reference=reference,
        disturbances=[stray("crossed fields", e=(0.0, RING_E), bz=RING_B)],
        max_free=max_free, free_charges=[RAD_Q], free_speeds=[1.0, 1.4, 1.8, 2.2],
        region=(2, 2, 28, 18))



def jackson_in_step():
    return ring_level(
        "Jackson Pr. 14.23: in step",
        "Jackson Pr. 14.23: in crossed fields a charge circles (its cyclotron motion, "
        "$\\omega_0 = 0.8$) while its circle drifts to the right into the detector (the E×B "
        "drift). Circling, it radiates at $\\omega_0$. The receiver far away at 90° counts "
        "the band 0.6–1.0, and it sees every particle: their fields add. Two charges "
        "circling in step send twice the field, so four times the energy. Add a charge "
        "that circles in step with this one, so that the receiver collects at least "
        "0.00014 per steradian, about three times what the charge sends alone. (c = 4; "
        "radiation reaction is included.)",
        (0.6, 1.0), (1.4e-4, 1e3),
        [free_charge(6, 7, RAD_Q, 180.0, RING_VG - RING_VD)], 3)



def jackson_quiet_ring():
    return ring_level(
        "Jackson Pr. 14.23: a quiet ring",
        "Jackson Pr. 14.23: the same drifting circle, but now the receiver must stay "
        "quiet: at most 0.000007 per steradian in the band 0.6–1.0, about a seventh of "
        "what the charge sends alone. Two charges on opposite sides of one circle, "
        "circling together, form a dipole that never turns: their radiation at "
        "$\\omega_0$ cancels. Put a second charge on the far side of the circle, circling "
        "about the same drifting centre: in the lab it moves at the drift plus its "
        "circling velocity.",
        (0.6, 1.0), (0.0, 7e-6),
        [ring_partner(6, 12, RING_U, 0.0)], 3)



def jackson_ring_of_four():
    return ring_level(
        "Jackson Pr. 14.23: a ring of four",
        "Jackson Pr. 14.23: N equal charges spaced evenly around one circle radiate only "
        "at the multiples of $N\\omega_0$ (N² times as strongly as one charge there), and a "
        "ring of many radiates almost nothing: a steady current does not radiate. The "
        "receiver now counts the band 0.6–2.6, which holds the first three harmonics "
        "($\\omega_0$, $2\\omega_0$, $3\\omega_0$); keep it below 0.000013 per steradian. An "
        "opposite pair cancels $\\omega_0$, but its fields add at $2\\omega_0$. Space "
        "charges around the circle so that only a harmonic outside the band survives.",
        (0.6, 2.6), (0.0, 1.3e-5),
        [ring_partner(8, 10, 0.0, -RING_U), ring_partner(6, 12, RING_U, 0.0),
         ring_partner(4, 10, 0.0, RING_U)], 5)
