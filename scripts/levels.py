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
    across it, around the landing point. Shots of one species share a detector unless
    `own` (e.g. the same ion launched at different times, sorted apart).
    """

    def __init__(self, strip, axis, size=1, own=False):
        self.strip, self.axis, self.size, self.own = strip, axis, size, own


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


def radiation_goal(axis_deg, half_deg, energy, band=None, window=None):
    """Radiation goal of a shot's detector (PHYSICS.md §3.4): the energy per steradian the
    flight radiates into the directions axis ± half (degrees), in all frequencies or in
    `band` = (omega_min, omega_max), must lie in `energy` = (min, max). With `window` =
    (t1, t2) a steady receiver: the mean power per steradian received over that window
    of its time (then `energy` bounds the power), and the particle need not arrive."""
    g = {"direction": [axis_deg, half_deg], "energy": list(energy)}
    if band is not None:
        g["band"] = list(band)
    if window is not None:
        g["window"] = list(window)
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
          value=None, tunable=False, drive=None, shape=None):
    """A box electrode (plate, slab or wall) standing on the plane; `tunable`: the player
    sets its potential with a power supply (a driven plate's source level); `drive`: a
    circuit (see `drive`); `shape` (tube levels): "rod" (a circle of diameter `length`) or
    "ring" (an annulus of mid-diameter `length`, wall `thickness`)."""
    bias = {"kind": kind}
    if value is not None:
        bias["value"] = value
    e = {"center": [x, y, 0], "length": length, "thickness": thickness,
         "height": height, "angle_deg": angle_deg, "bias": bias}
    if tunable:
        e["tunable"] = True
    if drive is not None:
        e["drive"] = drive
    if shape is not None:
        e["shape"] = shape
    return e


def dc(v=0.0):
    """A constant source (its level the player's power supply sets on a tunable plate)."""
    return {"kind": "dc", "value": v}


def drive(source, r, l=0.0, c=0.0, switch=None):
    """A circuit driving a plate or a circular coil (PHYSICS.md 2.10): `source` through
    the resistance r, with a series inductance l (plates) or capacitance c (coils);
    `switch` = (closed at first, time it toggles) makes r a switch (plates without l)."""
    d = {"source": source, "resistance": r}
    if l:
        d["inductance"] = l
    if c:
        d["capacitance"] = c
    if switch is not None:
        d["switch"] = {"closed": True, "at": switch[1]} if switch[0] else {"at": switch[1]}
    return d


def slider(vmax, *candidates):
    """The potentials of a power supply or a player plate: a slider over [-vmax, vmax].
    The owner (2026-10-06): voltages are tuned with sliders whose range is symmetric about
    0 and whose ends do not give the answer away (no answer near 0, the ends or the
    middle of either half: `guessable`, also checked by the level tests). The
    candidates are the values the solver and the difficulty analysis try."""
    return sorted({-vmax, vmax, *candidates})


def guessable(v, vmax):
    """Whether a potential is near a point of its slider that a player would try first:
    0, either end, or the middle of either half (within 5 % of the range's half)."""
    return any(abs(abs(v) - p) < 0.05 * vmax for p in (0.0, 0.5 * vmax, vmax))


def slider_problems(lvl):
    """The slider rules (see `slider`) on a built level: each list symmetric about 0, and
    no reference value at a guessable point of its own slider (a plate's potentials, or
    what its power supply operates: an electrode's potentials, a coil's ramp rates)."""
    out = []
    limits = lvl["limits"]
    lists = {"supply voltages": limits.get("supply_voltages", []),
             "coil rates": limits.get("coil_rates", []),
             "plate voltages": limits.get("plate_voltages", []),
             "solenoid fields": limits.get("solenoid_fields", [])}
    for what, vals in lists.items():
        if vals and min(vals) != -max(vals):
            out.append(f"{what} slider [{min(vals)}, {max(vals)}] not symmetric about 0")
    # A tunable coil's supply sits on its terminal (`Coil::supply_node`); an undriven
    # coil's is a ramp rate, a driven coil's its source's voltage.
    coil_centres = [[c["center"][0] + math.floor(c["radius"] + 0.5), c["center"][1], c["center"][2]]
                    for c in lvl.get("coils", [])
                    if c.get("shape") == "circle" and c.get("tunable") and not c.get("drive")]
    for e in lvl["reference_solution"]:
        if e["kind"] == "plate":
            vals = lists["plate voltages"]
        elif e["kind"] == "supply" and e["node"] == lvl.get("tube", {}).get("solenoid"):
            vals = lists["solenoid fields"]
        elif e["kind"] == "supply":
            vals = lists["coil rates" if e["node"] in coil_centres else "supply voltages"]
        else:
            continue
        if vals and max(vals) > 0 and guessable(e["value"], max(vals)):
            out.append(f"{e['kind']} at {e['node'][:2]}: {e['value']:g} is guessable on "
                       f"the slider +-{max(vals):g}")
    return out


def heat_goal(metal, lo, hi, start, end):
    """A tube level's goal on metal `metal`'s mean heat (the field's work on its carriers per
    time) over t = start..end, between lo and hi."""
    return {"metal": metal, "min": lo, "max": hi, "start": start, "end": end}


def cut(x0, y0, x1, y1):
    """The player's insulating cut from the node (x0, y0) to (x1, y1) (as drawn in the game:
    `Level::cut_between`)."""
    dx, dy = x1 - x0, y1 - y0
    e = {"node": [x0, y0, 0], "kind": "cut", "value": math.hypot(dx, dy)}
    angle = math.degrees(math.atan2(dy, dx))
    if angle:
        e["angle_deg"] = angle
    return e


def current_goal(electrode, lo, hi, start, end, carriers=None):
    """A tube level's current goal: the current into `electrode`, averaged over
    t = start..end, between lo and hi; `carriers`: None for the cathode's, "ions" for the
    ion source's (a bipolar level)."""
    g = {"electrode": electrode, "min": lo, "max": hi, "start": start, "end": end}
    if carriers is not None:
        g["carriers"] = carriers
    return g


def tube(cathode, electrode, lo, hi, start, end, emit_toward_deg=None, q_over_m=-1.0,
         saturation=None, projectiles=(), b_z=0.0, solenoid=None, circuits=(),
         potentials=(), goals=(), ions=None, step=None, particle_rate=None, particles=None,
         emitter_mesh=None, output=None, metals=(), solenoid_drive=None, heat_goals=()):
    """A tube level (z-invariant, docs/TUBES.md): electrode `cathode` emits electrons
    (from its faces toward `emit_toward_deg` only, if given; at most `saturation` per unit
    length of its surface, if given: a temperature-limited cathode); the goal is the
    current into `electrode`, averaged over t = start..end, between lo and hi; with
    `projectiles` (see `projectile`) also each one with a detector arriving there. `b_z`:
    a uniform field along z; `solenoid` = (x, y): the terminal of a solenoid the player
    powers (a supply there sets B along z, over `level(solenoid_fields=...)`). `goals`:
    further current goals (`current_goal`); `ions` = (electrode, q/m[, coated face's
    direction]): an ion source (the bipolar diode). Tube levels are Newtonian (pass
    c=None to `level`). `solenoid_drive`: the solenoid's field over time (a source, see
    `dc`, `sine`; its level set by the supply on the terminal), which induces an electric
    field circulating about the board's centre. `heat_goals`: goals on a metal's mean heat
    (`heat_goal`)."""
    t = {"cathode": cathode, "charge_per_mass": q_over_m}
    if heat_goals:
        t["heat_goals"] = list(heat_goals)
    if solenoid_drive is not None:
        t["solenoid_drive"] = solenoid_drive
    all_goals = ([current_goal(electrode, lo, hi, start, end)] if electrode is not None
                 else []) + list(goals)
    if all_goals:
        t["goals"] = all_goals
    if step is not None:
        t["step"] = step
    if particle_rate is not None:
        t["particle_rate"] = particle_rate
    if particles is not None:
        t["particles"] = particles
    if emitter_mesh is not None:
        t["emitter_mesh"] = emitter_mesh
    if output is not None:
        t["output"] = {"circuit": output[0], "c": output[1]}
    if ions is not None:
        t["ions"] = {"electrode": ions[0], "charge_per_mass": ions[1]}
        if len(ions) > 2:
            t["ions"]["emit_toward_deg"] = ions[2]
    if circuits:
        t["circuits"] = list(circuits)
    if potentials:
        t["potential_goals"] = list(potentials)
    if emit_toward_deg is not None:
        t["emit_toward_deg"] = emit_toward_deg
    if saturation is not None:
        t["saturation"] = saturation
    if projectiles:
        t["projectiles"] = list(projectiles)
    if b_z:
        t["b_z"] = b_z
    if solenoid is not None:
        t["solenoid"] = [solenoid[0], solenoid[1], 0]
    if metals:
        t["metals"] = list(metals)
    return t


def metal(lo, hi, contacts, density, thermal_speed, mean_free_path, core_radius,
          carriers_per_area, seed=0, q_over_m=-1.0):
    """A metal in a tube level (PHYSICS.md 2.11, the Drude metal): the rectangle from the
    grid corner `lo` to `hi`, its `contacts` (electrodes: boxes along its sides), its
    carriers' charge density (magnitude), q/m and thermal speed, its mean free path and
    its ion cores' radius (their number per area follows: a/(2l) at most 0.005), its
    macroparticles per area in the preview."""
    return {"corners": [[lo[0], lo[1], 0], [hi[0], hi[1], 0]], "contacts": list(contacts),
            "density": density, "charge_per_mass": q_over_m, "thermal_speed": thermal_speed,
            "mean_free_path": mean_free_path, "core_radius": core_radius,
            "carriers_per_area": carriers_per_area, "seed": seed}


def contact(x, y, length, vertical=True, tunable=False):
    """A metal's contact: a box 2 thick along one of its sides (its centre on the grid, its
    face on the metal's edge)."""
    return plate(x, y, length, thickness=2.0, angle_deg=90.0 if vertical else 0.0,
                 tunable=tunable)


# The metals' common material (docs/DRUDE.md, PHYSICS.md 2.11): carriers of charge density
# 2, q/m = -1, thermal speed 2 (kT/q = 4, small against the supplies' tens: a contact's own
# voltage stays a few per cent), mean free path 1 (a quarter of the strips' width), cores
# of radius 0.01 (a/2l = 0.005), 50 macroparticles per area; the step resolves the plasma
# oscillation (omega_p dt = 0.18).
METAL = dict(density=2.0, thermal_speed=2.0, mean_free_path=1.0, core_radius=0.01,
             carriers_per_area=50.0)
METAL_STEP = 0.035


def tube_circuit(electrode, source, resistance=None, capacitance=0.0, supply_offset=False,
                 linked=(), inductance=0.0, return_to=None):
    """An electrode's circuit in a tube level (per unit length): `source` (see `dc`,
    `sine`) through `resistance` (0: an ideal source; None: floating) and `inductance` (a
    tank's coil) in series, and `capacitance` to ground; the branch ends on ground, or on
    electrode `return_to`; a tunable electrode's supply sets the source's level, or with
    `supply_offset` a sine's offset; `linked`: further electrodes wired to it (on an ideal
    source a grid's wires; on a tank or a resistance one node with it: a magnetron's
    strapped vanes)."""
    c = {"electrode": electrode, "source": source}
    if inductance:
        c["inductance"] = inductance
    if return_to is not None:
        c["return_to"] = return_to
    if linked:
        c["linked"] = list(linked)
    if resistance is not None:
        c["resistance"] = resistance
    if capacitance:
        c["capacitance"] = capacitance
    if supply_offset:
        c["supply_offset"] = True
    return c


def sine(amplitude, omega, offset=0.0, phase=0.0):
    """A sine source: offset + amplitude sin(omega t + phase)."""
    s = {"kind": "sine", "amplitude": amplitude, "omega": omega}
    if offset:
        s["offset"] = offset
    if phase:
        s["phase"] = phase
    return s


def potential_goal(electrode, measure, lo, hi, start, end, omega=None, reference=None):
    """A tube level's goal on an electrode's potential over t = start..end: its "mean", its
    "swing" (half the peak-to-peak) or its oscillation's "amplitude" at `omega` (a drift
    removed) between lo and hi; against electrode `reference`'s potential if given (a
    voltmeter across the two)."""
    g = {"electrode": electrode, "measure": measure, "min": lo, "max": hi,
         "start": start, "end": end}
    if omega is not None:
        g["omega"] = omega
    if reference is not None:
        g["reference"] = reference
    return g


def projectile(x, y, direction_deg, speed, charge, mass, launch=0.0, detector=None):
    """A line charge fired through a tube level at time `launch`; it must arrive in the
    box `detector` = ((x0, y0), (x1, y1)) if given."""
    p = {"node": [x, y, 0], "direction_deg": direction_deg, "speed": speed,
         "charge": charge, "mass": mass}
    if launch:
        p["launch"] = launch
    if detector is not None:
        (x0, y0), (x1, y1) = detector
        p["detector"] = [[x0, y0, 0], [x1, y1, 0]]
    return p


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
          free_charges=(), free_speeds=(), free_mass=1.0, free_radius=0.3, dielectrics=(),
          coil_rates=(), ferrites=(), tube=None, solenoid_fields=(), max_cuts=0):
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
    if coil_rates:
        limits["coil_rates"] = list(coil_rates)
    if solenoid_fields:
        limits["solenoid_fields"] = list(solenoid_fields)
    if max_free:
        limits["max_free"] = max_free
    if max_cuts:
        limits["max_cuts"] = max_cuts
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
        **({"dielectrics": list(dielectrics)} if dielectrics else {}),
        **({"ferrites": list(ferrites)} if ferrites else {}),
        **({"gates": list(gates)} if gates else {}),
        **({"tube": tube} if tube else {}),
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
        "a power supply (yellow frame). Click the plate to switch it on, and set its "
        "potential with the slider in the panel. A positive plate pushes positive ions away, a negative one pulls "
        "them. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2))],
        electrodes=[plate(11, 13, 8, tunable=True), plate(11, 7, 8)],
        supplies=slider(3.4e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, -1.2e5)])


# Circuits (PHYSICS.md 2.9-2.10): a source through a resistor (with an inductor or a
# switch) drives a plate; the plates of "Power supply" (the top one's capacitance against
# the grounded lower plate is C = 2.63, from the boundary elements). The plate levels are
# Newtonian (c = infinity), where the circuits' quasi-static model is exact; the ions
# (T0 = 0.5, m = 1) fly at 1 cell per time unit and cross the plates (x = 7 to 15) at t
# = 7 to 15 after their launch.

def charging_a_plate():
    # Through R = 3 the top plate charges with tau = RC = 7.9: as the ions cross it holds
    # 59 to 85 % of its supply's voltage. -120k lands the beam at y = 12.49, -60k at 11.22
    # (in "Power supply", charged at once, -60k gave 12). The supply solves from -140k
    # to -85k (generator window); the slider's ends (+-310k) and their middles lie
    # outside.
    return level(
        "Charging a plate",
        "The top deflection plate is connected to its power supply through a resistor R. "
        "A plate is a capacitor: its charge has to flow in through R, so its potential "
        "rises as $V(t) = V_s(1 - e^{-t/RC})$ towards the supply's $V_s$, with the time "
        "constant $\\tau = RC$, here about 8 (C comes from the plates' shape). The ions "
        "cross the plates while it is still charging: the Circuits panel shows V(t), with "
        "the launch marked. Set the supply so that the beam lands in the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 12, 30, 13)))],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 3.0)), plate(11, 7, 8)],
        supplies=slider(3.1e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, -1.2e5)],
        c=None)


def chopper():
    # The top plate carries no charge until a switch connects it to its supply through
    # R = 0.4 at t = 30 (tau = 1.05): the ion launched at t = 0 crosses the plates at
    # t = 7 to 15 and flies straight on; the one launched at t = 40 meets a fully
    # charged plate: 120k lands it at y = 5.54, 60k at 8.23 (none: 10). The supply
    # solves from 90k to 165k; the slider's ends (+-360k) and their middles lie outside.
    return level(
        "Chopper",
        "A switch connects the top plate to its power supply at t = 30, through a small "
        "resistor: the plate charges within about one time unit, so the switch turns the "
        "deflector on. Two ions are launched, at t = 0 and at t = 40. Set the supply so "
        "that each reaches its own detector: the first flies past the plates before the "
        "switch closes, the second after.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 9, 30, 11))),
               shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 5, 30, 7)), time=40.0)],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 0.4, switch=(False, 30.0))),
                    plate(11, 7, 8)],
        supplies=slider(3.6e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, 1.2e5)],
        c=None)


def ringing_plate():
    # With L = 4.66 in series (and R = 0.2) the top plate (C = 2.63) rings: omega =
    # 1/sqrt(LC) = 0.286, the first maximum at t = pi/omega = 11, overshooting the
    # supply's voltage by 79 %. The ions cross the plates at t = 7 to 15, around it:
    # 120k lands the beam at y = 3.98 (6 cells down; "Power supply"'s static plate gave
    # 4), 60k at 7.60. The supply solves from 107.5k to 170k (generator window); the
    # slider's ends (+-205k) and their middles (+-102.5k) lie outside.
    return level(
        "Ringing plate",
        "Now an inductor L sits between the supply and the top plate. With the plate's "
        "capacitance it makes an LC circuit: the charge rushes in, overshoots and swings "
        "around the supply's voltage at $\\omega = 1/\\sqrt{LC}$, slowly damped by the "
        "small resistance. The first swing peaks just as the ions pass, at almost twice "
        "the supply's voltage (see the Circuits panel). Set the supply so that the beam "
        "lands in the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 3, 30, 5)))],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 0.2, l=4.66)),
                    plate(11, 7, 8)],
        supplies=slider(2.05e5, -1.4e5, -7e4, 7e4, 1.4e5), reference=[supply(11, 13, 1.4e5)],
        c=None)


def induction_kick():
    # Jackson's Faraday level (5.15) with the coil on a circuit: its current rises
    # through R from 0, kappa(t) = kappa_inf (1 - e^{-t/tau}), tau = L/R = 20, with L =
    # N/c^2 its self-inductance (N = 448.5 for radius 8 and wire 0.1). A driven coil
    # needs a finite c: c = 1e4, where the circuit's quasi-static parameter (light time
    # across the arena / tau) is 1.5e-4, its corrections of order its square (the
    # first-order retardation of a closed current vanishes), and the ion (v ~ 1)
    # Newtonian to 1e-8. The
    # canonical angular momentum holds p_phi = -q kappa(t) A_phi(rho) for the ion
    # starting at rest: the kick follows the flux, not its rate. kappa_inf = 1e6 (V/R =
    # kappa_inf c^2 = 1e14, R = L/tau = 2.24e-7, V = 2.24e7): at 4 cells from the axis,
    # where A_phi ~ 1.6 per unit kappa, a kick p ~ 1.6. As in Faraday's level the
    # detector asks for T >= 0.3, which only the induced field supplies: three charges
    # of at most 0.2M at least 4 cells away give at most 0.15.
    tau, n_coil, c = 20.0, 448.5, 1e4
    r = n_coil / c ** 2 / tau
    coil = circle_coil(15, 10, 8.0, 0.0)
    coil["drive"] = drive(dc(1e6 * c ** 2 * r), r)
    return level(
        "Induction kick",
        "The coil's current now comes from a circuit: switched on at t = 0 through a "
        "resistor, it rises as $I(t) = (V/R)(1 - e^{-t/\\tau})$ with $\\tau = L/R$, L the "
        "coil's own inductance. While it rises, the changing flux induces an electric "
        "field around the axis that drives the particle at rest; once the current "
        "settles, the induction stops. Faraday's law in integral form: the kick depends "
        "only on how much the flux changed, not on how fast. Bring the particle into the "
        "detector with at least 0.3 of kinetic energy.",
        shots=[shot(1e-6, 1.0, (15, 6), 0.0, 1e-6, box((18, 11, 22, 15), kinetic=(0.3, 10.0)))],
        coils=[coil], c=c, t_max=200.0,
        max_charges=3, magnitudes=[m * M for m in (0.05, 0.1, 0.2)], region=(8, 10, 22, 17))


def pulse(high, rise, delay=0.0, width=1e4, fall=0.0, low=0.0):
    """A trapezoid source (PHYSICS.md 2.9): `low` until `delay`, up to `high` over `rise`,
    held for `width`, down over `fall` (its top the player's power supply sets)."""
    return {"kind": "pulse", "low": low, "high": high, "delay": delay, "rise": rise,
            "width": width, "fall": fall}


def betatron(core=1.34e8, guide=2.3e3, kinetic=(8.6, 11.3), core_slider=3.1e8,
             guide_slider=8.3e3, rise=300.0, detector=None):
    # Wideroe's betatron: an ion at rest at R = 8 between a core coil (radius 4) and a
    # guide coil (radius 12), concentric. Its canonical angular momentum is 0, so
    # p = q A_phi(R, t); a circle at R needs p = q R B_z(R, t). With both coils' strengths
    # rising in proportion that holds at every t iff A_phi(R) = R B_z(R) for the mixture:
    # the core must rise 2.057 times as fast as the guide (exact loop integrals; field
    # index n = -3.8 < 1, radially stable). The guide ramps at the player's rate g; the
    # core follows a trapezoid source through a fast circuit (L/R = 2, c = 1e4): its
    # current rises in step with the source for `rise`, then holds. The rise needs
    # kappa_top / rise = 2.057 g, i.e. a source top V = 2.057 g rise N / tau (N = 189.5,
    # its self-flux). When the core stops, the flux through the orbit stops growing while
    # the guide's field rises on: the orbit shrinks onto the detector between the coils,
    # with the energy the betatron gave it. A wrong ratio moves the orbit during the rise:
    # inward onto the detector too slow, or outward into the guide. The core's supply is
    # its source's top voltage (`supply_voltages`), the guide's its ramp rate. Measured
    # (generator trace): the orbit holds at r = 8.00 +- 0.04 through the rise while T
    # climbs to 4.6, then shrinks onto the detector (r = 6 to 7) at t = 567, T ~ 10.
    # Solving intervals (generator window): the core 130M to 142.5M (the energy is its
    # flux's), the guide 1.48k to 3.39k (and islands up to 3.59k: a wrong ratio settles
    # the orbit at another radius, which the extraction still reaches); the sliders
    # (+-310M, +-8.3k) keep 0, their ends and middles outside.
    c = 1e4
    n_core = 189.5
    tau = 2.0
    r = n_core / c ** 2 / tau
    core_coil = circle_coil(15, 15, 4.0, 0.0)
    core_coil["drive"] = drive(pulse(0.0, rise), r)
    core_coil["tunable"] = True
    guide_coil = circle_coil(15, 15, 12.0, 0.0)
    guide_coil["tunable"] = True
    detector = detector or box((8, 14, 9, 16), kinetic=kinetic)
    return level(
        "Betatron",
        "Wideroe's betatron accelerates without any electrode: a changing flux induces an "
        "electric field around it, and a magnetic field holds the particle on its circle. "
        "The particle starts at rest between two coils. The inner one (the core, on a "
        "circuit whose source rises and then holds) supplies most of the flux, the outer "
        "one (the guide, ramped) the field at the orbit. The orbit keeps its radius only if "
        "the flux inside it grows twice as fast as the guide field times the orbit's area, "
        "$\\Phi(R) = 2\\pi R^2 B(R)$: then the particle gains energy turn after turn on "
        "the same circle. When the core's source stops rising, the flux stops growing but "
        "the guide's field rises on, and the orbit shrinks onto the detector. Set the two "
        "supplies so that the particle arrives there with its energy in the window.",
        grid=(30, 30), c=c, t_max=900.0,
        shots=[shot(1e-6, 1.0, (23, 15), 90.0, 1e-9, detector)],
        coils=[core_coil, guide_coil],
        supplies=slider(core_slider, core), coil_rates=slider(guide_slider, guide),
        reference=[supply(19, 15, core), supply(27, 15, guide)])


def pulse_sorter(refs=(1.46e5, 9.6e4, -3.36e5), detectors=None):
    # Arc 5 finale: three ions, launched at t = 0, 30 and 60, pass three deflection plate
    # pairs in a row (x = 5-9, 16-20, 27-31; gaps of 8 cells), each top plate on its own
    # circuit and its own supply slider. A charges through R = 20 (tau ~ 25 with these
    # small plates' C ~ 1.25), B is switched on at t = 40, C rings (L = 130: period ~ 80).
    # Measured deflections per volt (generator check, linear at 30 kV), each ion relative
    # to the largest: A 0.28, 0.82, 1; B 0, 0.96, 1; C 1, 0.89, 0.13. Each ion weighs the
    # three circuits differently, so the three landing points fix the three supplies and
    # no slider is redundant; the player must read the circuits' plots. The reference
    # (+146k, +96k, -336k: 0.32, 0.21 and 0.73 of the slider's +-460k, away from its
    # guessable points) lands the ions at y = 12.3, 9.3 and 6.5, top to bottom in launch
    # order, every flight clear of the plates; the detectors are placed around them.
    # Small, low plates (4 x 2 cells, as in the microscope column): six full-size ones
    # took 5.5 s to set up, over the sandbox budget.
    def pair(x, d):
        return [plate(x, 14, 4, height=2.0, tunable=True, drive=d), plate(x, 6, 4, height=2.0)]
    electrodes = (pair(7, drive(dc(), 20.0)) + pair(18, drive(dc(), 0.4, switch=(False, 40.0)))
                  + pair(29, drive(dc(), 0.5, l=130.0)))
    detectors = detectors or [Auto((38, 0, 40, 20), "x", 1, own=True) for _ in range(3)]
    shots = [shot(1e-6, 1.0, (0, 10), 0.0, 0.5, d, time=t)
             for t, d in zip((0.0, 30.0, 60.0), detectors)]
    return level(
        "Pulse sorter",
        "Finale of the circuits: three ions, launched at t = 0, 30 and 60, fly through "
        "three pairs of deflection plates, and each top plate hangs on its own circuit "
        "with its own power supply. The first charges through a resistor, the second is "
        "switched on at t = 40, the third rings with an inductor. Each ion meets the "
        "three plates at a different stage of their circuits (the Circuits panel shows "
        "each plate's V(t), with the launches marked). Set the three supplies so that "
        "every ion lands in its own detector.",
        grid=(40, 20), shots=shots, electrodes=electrodes, c=None, t_max=160.0,
        supplies=slider(4.6e5, -3e5, -2e5, -1e5, 1e5, 2e5, 3e5),
        reference=[supply(x, 14, v) for x, v in zip((7, 18, 29), refs) if v])


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
        electrodes=apertures, supplies=slider(4.9e5, 1e5, 2e5, 3e5, 4e5),
        reference=[supply(14, 15, 3e5), supply(14, 5, 3e5)])


def build_a_deflector():
    return level(
        "Build a deflector",
        "Now you build the electrodes. Place a metal plate at a potential of your choice "
        "near the beam (R turns it, Q/E changes its potential). A plate at a positive "
        "potential pushes positive ions away. Plates keep a cell away from other metal "
        "and from the other elements.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 1))],
        max_plates=3, plate_voltages=slider(2.3e5, -1e5, -7.5e4, -5e4, 5e4, 7.5e4, 1e5),
        region=(5, 3, 22, 17), reference=[player_plate(10, 12, 7.5e4)])


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


def glass(x, y, length, thickness, eps, height=4.0, angle_deg=0.0):
    """A dielectric box of relative permittivity eps (PHYSICS.md 2.7)."""
    return {"center": [x, y, 0], "length": length, "thickness": thickness, "height": height,
            "angle_deg": angle_deg, "permittivity": eps}


def ferrite(x, y, length, thickness, mu, height=4.0, angle_deg=0.0):
    """A ferrite box (high-mu insulator) of relative permeability mu (PHYSICS.md 2.7)."""
    return {"center": [x, y, 0], "length": length, "thickness": thickness, "height": height,
            "angle_deg": angle_deg, "permeability": mu}


def ferrite_shield():
    # Ferrites' introduction: a strong magnet (2M at (15, 15)) beside the beam line,
    # behind a ferrite plate (mu = 1000, 1.5 cells thick, 2 high, from x = 10 to 20).
    # Measured (generator check): without the plate the magnet bends the beam 2.27 cells
    # up at the right edge; with it 0.98 (a taller plate shields more: 0.50 at height 4;
    # mu = 5 still 1.50). The detector (one cell) is placed around the reference's landing
    # with the plate; without it the beam lands ~1.3 cells higher, which the build checks
    # (NEEDS_BODY: the reference fails without the ferrite).
    return level(
        "Ferrite shield",
        "A ferrite is a magnetically soft insulator: it carries no current, but its "
        "magnetic domains line up with a field and draw the field lines into it, so "
        "behind a ferrite plate a magnet's field is much weaker (a magnetic shield; "
        "electric fields pass, as through any insulator). The strong magnet behind the "
        "plate still nudges the beam a little. Place magnets to bring it into the "
        "detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((28, 0, 30, 20), "x", 1))],
        elements=[magnet(15, 15, 2 * M)],
        ferrites=[ferrite(15, 12, 10, 1.5, 1000.0, height=2.0)],
        max_magnets=3, strengths=[m * M for m in (0.5, 1, 2)], region=(4, 1, 26, 9),
        reference=[magnet(17, 7, -1 * M)])


def glass_screen():
    # Shielding's charge (2M at (15, 14)) behind a glass slab (eps = 4, 1.5 cells thick,
    # 2 high, from x = 10 to 20) instead of a metal plate. The glass's bound charge
    # weakens the field beyond it but does not cancel it, as metal does: the beam is still
    # thrown off, less, and the player's charges must make up the rest. The build checks
    # that the reference fails without the glass (NEEDS_BODY): the glass matters
    # (measured: the solver's one charge, -1M at (22, 2), brings the beam in with it; the
    # beam crashes into that charge without it). Small and low (setup 0.8 s; a 12 x 2 x 4
    # slab took 3.7 s).
    return level(
        "Glass screen",
        "Shielding with glass instead of metal. A dielectric is an insulator: its charges "
        "cannot flow, but its molecules polarize, and the bound charge on its surfaces "
        "weakens the field that crosses it, here by roughly the permittivity ε = 4, never "
        "to zero as metal does. The strong charge behind the glass still pushes the beam "
        "off; place charges to bring it into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 8, 30, 12)))],
        elements=[charge(15, 14, 2 * M)],
        dielectrics=[glass(15, 12, 10, 1.5, 4.0, height=2.0)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1)], region=(4, 1, 26, 9))


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
    volts = slider(7.1e5, -2e5, -1e5, 1e5, 2e5, 3e5, 4e5, 5e5)
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
            ("glass_screen", glass_screen),
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
            ("ferrite_shield", ferrite_shield),
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
            ("ramp_the_coil", ramp_the_coil),
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
    ("Circuits", [
        ("Introduction", [
            ("charging_a_plate", charging_a_plate),
            ("chopper", chopper),
            ("ringing_plate", ringing_plate),
        ]),
        ("Intermediate", [
            ("induction_kick", induction_kick),
            ("betatron", betatron),
        ]),
        ("Master", [
            ("pulse_sorter", pulse_sorter),
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
    ("Jackson: electrostatics", [
        ("Introduction", [
            ("jackson_own_image", jackson_own_image),
            ("jackson_two_spheres", jackson_two_spheres),
        ]),
        ("Intermediate", [
            ("jackson_slot", jackson_slot),
            ("jackson_golden_ratio", jackson_golden_ratio),
            ("jackson_multipoles", jackson_multipoles),
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
            ("jackson_in_step", jackson_in_step),
        ]),
        ("Intermediate", [
            ("jackson_quiet_turn", jackson_quiet_turn),
            ("jackson_thomson", jackson_thomson),
            ("jackson_cross_section", jackson_cross_section),
            ("jackson_braking", jackson_braking),
            ("jackson_quiet_ring", jackson_quiet_ring),
            ("jackson_ring_of_four", jackson_ring_of_four),
        ]),
        ("Master", [
            ("jackson_undulator", jackson_undulator),
        ]),
    ]),
]

# Inside the components: vacuum tubes (docs/TUBES.md, PHYSICS.md 2.11). z-invariant:
# the electrodes are long prisms, the electrons long lines of charge, the currents per
# unit length; Newtonian (tube voltages give v ~ 0.03 c in SI; here c is infinite).

def vacuum_diode():
    # The planar diode of crates/level/tests/tubes.rs: plates 10 long, 1 thick, faces 7
    # apart, inside the tube's grounded envelope (PHYSICS.md 2.11). Preview currents
    # (`generator tube-scan`): 1.00 at 22, 1.15 at 24, 1.29 at 26, 1.45 at 28; the
    # verification's are ~7 % lower (I ~ V^1.5) with an error of ~10 % (the envelope made
    # the convergence first order, PHYSICS.md 2.11). The goal 0.9 to 1.5 is verified for
    # about 23.6 to 28 (1.0 to 1.4 left only 25.3 to 26.9); the reference 26: 1.203 +-
    # 0.124.
    return level(
        "Vacuum diode",
        "Inside a vacuum tube. The left plate is a hot cathode: its right face emits "
        "electrons (blue dots), as many as its field pulls off. The right plate, the "
        "anode, is on your power supply. Set its voltage so that a current between 0.9 "
        "and 1.5 reaches the anode (averaged over $t = 10$ to $20$). The electrons in "
        "flight repel the ones behind them: this space charge chokes the current, which "
        "grows as $I \\propto V^{3/2}$ (Child and Langmuir), faster than in proportion "
        "to the voltage. The potential map shows the space charge as a flat region in "
        "front of the cathode, where the field is pulled down to zero.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(8, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(70.0, 13.0, 26.0, 52.0),
        tube=tube(0, 1, 0.9, 1.5, 10.0, 20.0, emit_toward_deg=0.0),
        reference=[supply(16, 8, 26.0)])


def sideways():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11, V8): a diode (30 across a gap of
    # 9: Hull's planar cutoff sqrt(2V)/d = 0.86) in the player's solenoid field, and a
    # collector at 10 above the top end of the gap; inside the grounded envelope. Preview
    # collector current: 0.47 without field (the anode 0.96), nothing with B negative (the
    # drift goes down); rising as the electrons curve upwards, 0.83 at B = 0.8, 1.28 at
    # 0.95, 1.34 at 1.1, 1.27 at 1.2, 1.15 at 1.3, while the anode's falls to 0 from 1.1
    # (the cutoff; below it, 0.13 at 0.95, comes in at the gap's ends, where the
    # translation symmetry behind Hull's theorem is broken); above ~1.35 the electrons
    # curl back onto the cathode, which takes them (0.87 at 1.4, 0.50 at 1.5). The goal
    # >= 1.1 asks for B in ~0.95-1.3: it solved on 0.925-1.225, 3.8 % of the solenoid's +-4
    # (the showcase rule's review, 2026-10-20). At least 0.9 (still nearly twice the 0.47
    # without field) asks for ~0.82-1.38, between Hull's cutoff and the electrons' curling
    # back onto the cathode.
    return level(
        "Sideways",
        "Around this tube is a solenoid: switch on its supply (the yellow coil at the top "
        "right) to fill the tube with a magnetic field along its length, out of the board "
        "when positive. Electrons crossing it curve, and above a critical field (Hull's "
        "cutoff) they cannot reach the anode at all: they drift along the cathode instead, "
        "sideways, perpendicular to both the electric and the magnetic field. Send the "
        "electrons into the collector above the gap: at least 0.9 must land there "
        "(averaged over $t = 20$ to $40$).",
        grid=(24, 22), c=None, t_max=60.0,
        electrodes=[plate(6, 11, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 11, 10, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(9, 19, 6, thickness=1.0, kind="potential", value=10.0)],
        solenoid_fields=slider(4.0, -1.1, 0.6, 1.1),
        tube=tube(0, 2, 0.9, 10.0, 20.0, 40.0, emit_toward_deg=0.0, solenoid=(23, 21)),
        reference=[supply(23, 21, 1.1)])


def bipolar_diode():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11): Vacuum diode's tube with the
    # anode's face emitting ions under the space-charge limit. The ions' q/m = 1/16 (16
    # times heavier per charge than the electrons, not 1836 or more: in a static field the
    # mass changes neither the paths nor the electrons' gain, only the time scale and the
    # ion current, sqrt(m_e/m_i) of the electrons' in a planar diode; heavier ions would
    # take longer to settle and longer runs). Preview over t = 20-30 (both currents settled:
    # at q/m = 1/25 the electrons' was 2.540 there and 2.553 over 30-50): electrons into the
    # anode 0.613 at 10 V, 1.723 at 20, 2.558 at 26, 3.496 at 32 (V^1.5; Vacuum diode's
    # 1.29 at 26 without ions: the ions' charge doubles it, 1.86 in Langmuir's infinite
    # planar diode), ions into the cathode 0.148, 0.410, 0.619, 0.830 (0.24 of the
    # electrons'; the planar law: 0.25). The goals, electrons >= 2.0 and ions <= 0.75, ask
    # for about 22 to 29.5 V in the preview; without the ions the electron goal would need
    # 41 V, far beyond the wear limit.
    return level(
        "Bipolar diode",
        "The same diode, but now the anode's face is coated with a salt that, hot, emits "
        "positive ions (orange), as many as its field pulls off, just as the cathode emits "
        "electrons. The ions fly the other way, and their positive charge cancels part of "
        "the electrons' cloud: more electrons get through than in the vacuum diode at the "
        "same voltage (1.86 times as many in an infinitely wide diode: Langmuir, 1929). But "
        "the ions strike the cathode and wear it away. Set the anode's voltage so that at "
        "least 2.0 of electrons reach the anode while at most 0.75 of ions strike the "
        "cathode (both averaged over $t = 20$ to $30$). Both currents grow as "
        "$I \\propto V^{3/2}$.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(8, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(40.0, 13.0, 26.0, 33.0),
        tube=tube(0, 1, 2.0, 10.0, 20.0, 30.0, emit_toward_deg=0.0,
                  ions=(1, 0.0625, 180.0),
                  goals=[current_goal(0, 0.0, 0.75, 20.0, 30.0, carriers="ions")]),
        reference=[supply(16, 8, 26.0)])


def through_the_cloud():
    # Designed with `generator tube-scan` and verified with `check` (PHYSICS.md 2.11),
    # inside the tube's grounded envelope. The owner (2026-10-12) found the flow at 30 V
    # too weak for the ion to cross a cloud: the anode is now at 60 (a current of 2.40
    # without plates, 2.5 times the former, its space charge twice), and the ion (q/m =
    # 0.02, launched at t = 10 when the flow is established) at 1.8: slower, at 60-70 V it
    # ends on the cathode whatever single plate helps. Without a plate it leaves the top at
    # x = 4.8, left of the detector (8 to 14). Of 90 single plates, 24 bring it in:
    # positive ones below the gap, left of its path, kick it rightwards (the original
    # idea) but draw electrons (2.7-3.7 at +30 to +75), and one above the cathode's top
    # end (6, 18) at +65/+75 pushes it right as it leaves the gap (2.85). Verified with
    # the anode at 60: (9, 2) vertical at +30, current 2.898 +- 0.036, the ion arrives at
    # t = 21.23 and 21.20 (the reference); (8, 2) horizontal at +45, 2.956 +- 0.018, t =
    # 21.33 and 21.32; (6, 18) vertical at +65, 2.830 +- 0.073, t = 23.80 and 23.85. Too
    # strong a plate below draws too much current (3.2 and more at +55 to +75).
    # The showcase rule's review (2026-10-20): with the reference plate the anode's current
    # is 1.75, 2.26, 2.58, 2.90, 3.23, 3.59, 4.14, 4.88, 5.56 at 40, 50, 55, 60, 65, 70, 80,
    # 90, 100 V, the ion leaving the top at x = 9.66, 8.42, 7.78, 7.15, 6.53, 5.91, 4.59, 3.41
    # (the detector 8 to 14) and at 100 turned back onto the cathode: the supply's +-100
    # shows it all. With the current goal 2.5-3.1 the supply solved on 54-61 (3.5 %); at
    # 2.3-3.3 the ion bounds it above (~62) and the current below (~51).
    # Measured: the supply 51-61 (5.0 %), the plate 30-49, next to +-100's half's middle (50):
    # the plates span +-120 (its middle 60; the interval 7.9 %).
    return level(
        "Through the cloud",
        "An ion (a long line of positive charge, the large dot) is fired up through the gap "
        "of a working diode, after its electrons have started to flow. Keep the diode "
        "running, with a current between 2.3 and 3.3 into the anode (averaged over "
        "$t = 15$ to $30$), and bring the ion into the detector at the top. The anode's "
        "field pushes positive charges toward the cathode, and the electrons' cloud pulls "
        "on the ion as it passes: at the voltage the current needs, the ion leaves the "
        "tube far left of the detector. Your two plates can push it back on course, but "
        "they act on the electrons as well, and the ion's own charge disturbs the flow it "
        "crosses.",
        grid=(24, 20), c=None, t_max=40.0,
        electrodes=[plate(6, 10, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 10, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(100.0, 40.0, 60.0, 80.0),
        max_plates=2, plate_voltages=slider(120.0, -40.0, 20.0, 30.0, 45.0, 65.0),
        plate_size=(3.0, 0.5, 4.0),
        tube=tube(0, 1, 2.3, 3.3, 15.0, 30.0, emit_toward_deg=0.0,
                  projectiles=[projectile(11, 0, 90.0, 1.8, 0.2, 10.0, launch=10.0,
                                          detector=((8, 18), (14, 20)))]),
        reference=[supply(16, 10, 60.0), player_plate(9, 2, 30.0, angle_deg=90.0)])


def electron_gun():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11). The
    # cathode is temperature limited (0.2 per unit length of its 8-long face: at most 1.6),
    # so pulling harder gains nothing. Inside the tube's grounded envelope the drift
    # region's potential peaks on the axis, which focuses the beam a little: with the
    # former 0.1 it reached the collector unaided (0.113 +- 0.042 for the goal 0.08), so
    # the cathode was made hotter, its beam's space charge spreading it harder. Without
    # plates now 0.155 +- 0.026 (the goal 0.24: missed). Two horizontal plates beside the
    # beam after the slit focus it (preview 0.25-0.31 at +2 to +8, x = 17 to 23).
    # Verified: (21, 13/7) at +6 0.294 +- 0.008 (the reference); (17, 13/7) at +6 0.288
    # +- 0.023; (19, 13/7) at +4 0.294 +- 0.019. At 0 the lens is too strong (0.192 +-
    # 0.025) and at +12 too weak (0.161 +- 0.012): both missed.
    # The plates' range (the showcase rule's review, 2026-10-20): the collector's current
    # (preview) with both plates at V: 0.180, 0.253, 0.281, 0.293, 0.273, 0.196, 0.155 at
    # V = 0, 2, 4, 6, 8, 10, 14; nothing at -10 to -40 (the beam turned back), 0.125, 0.101,
    # 0.098 at +20, +30, +40 (the plates, above the anode's 30, collect it). One plate alone
    # steers the beam off: each slider solves on ~3-7.5 with the other at 6, 2.8 % of +-80.
    # Every behaviour lies within +-40 (beyond it only more of the same), so the plates span
    # +-40 (that interval 5.6 %).
    return level(
        "Electron gun",
        "An electron gun: the hot cathode (left) gives at most 0.2 per unit length of its "
        "face (its temperature limits it: a stronger field pulls no more), the anode's slit "
        "lets a thin beam through, and the small collector on the right should catch it. "
        "But the beam's electrons repel one another: it spreads and misses. Place two plates "
        "beside the beam to focus it, so that at least 0.24 reaches the collector (averaged "
        "over $t = 20$ to $40$). A plate at a lower potential than the beam's surroundings "
        "pushes the electrons back toward the axis; too low, and it throws them back.",
        grid=(32, 20), c=None, t_max=60.0,
        electrodes=[plate(4, 10, 8, thickness=1.0, angle_deg=90.0),
                    plate(12, 14, 6, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(12, 6, 6, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(28, 10, 2, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0)],
        max_plates=2, plate_voltages=slider(40.0, -10.0, 3.0, 6.0, 12.0, 20.0),
        plate_size=(3.0, 0.5, 4.0),
        tube=tube(0, 3, 0.24, 5.0, 20.0, 40.0, emit_toward_deg=0.0, saturation=0.2),
        reference=[player_plate(21, 13, 6.0), player_plate(21, 7, 6.0)])


def ion_diode():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11): an inverted coaxial diode, the
    # cathode the ring (mid-diameter 10, wall 1: inner radius b = 4.5), the anode the rod
    # (radius a = 1) emitting ions (q/m = 1/16, as in Bipolar diode) on every face, in the
    # player's solenoid field. Hull's cutoff from r0 to r1 for a charge released at rest:
    # B = sqrt(2 m V/q) 2 r1/|r0^2 - r1^2|: for the electrons (b to a) 0.80 at 30 V, for the
    # ions (a to b) 4 * 2b/2a = 18 times that: the ions cross, the electrons are insulated
    # in a sheath along the ring, whose charge draws more ions (the magnetically insulated
    # ion diode). The sheath keeps spreading towards the rod (the impedance collapse of real
    # ion diodes), so the goal's window is a fixed time after switch-on. A first design (with
    # the game's usual 150 macroparticles per emitting face) was shelved: the sheath's ion
    # current depended mostly on the number of macroparticles (x1 to x8: 8.21, 7.42, 6.90,
    # 6.66 at 28 V, B = 1.8) and the verdicts' bars were 15-25 %. With 8 times as many (the
    # tree code makes them affordable: PHYSICS.md 2.11) it converges in their number (x16:
    # 8.30 against 8.19 at 30 V, B = 1.7), and the rest of the preview's error was the
    # emitters' release layer, which the sheath magnifies: emitters meshed twice as finely
    # bring the preview within 7.3 % of the verification at B = 1.7 (4.8 % at 2.4): bars of
    # ~10 %, a preview of ~3.4 s and a verification of ~13 s. (V, B) grid in the preview
    # (ions onto the ring / electrons leaking into the rod, t = 20-30):
    #   V \ B  1.5         1.7        1.9        2.2        2.6        3.0
    #   24     5.55/0      4.82/0     4.39/0     3.97/0     3.64/0     3.46/0
    #   27     7.09/0      6.12/0     5.52/0     4.92/0     4.47/0     4.21/0
    #   30     8.87/0.006  7.78/0     6.80/0     6.01/0     5.40/0     5.05/0
    #   33    11.84/0.094  9.88/0.002 8.31/0     7.20/0     6.42/0     5.95/0
    #   36    14.81/0.171 11.32/0    10.06/0     8.55/0     7.53/0     6.95/0
    # The boost is largest just above the cutoff (1.76 times the strong field's at 30 V).
    # The showcase rule's review (2026-10-20): with ions 5.5 to 8 the supply solved on 27-32
    # (3.1 % of +-80; the +-80 shows the leak above ~50). Widened to 5.0-8.5 (not 9: the ion
    # current climbs as B falls towards the electrons' cutoff, 8.87 at B = 1.5, and the
    # field's interval would then take in 1.5, its half's middle). Measured: the supply
    # 26-33.5 (4.7 %), the field 1.55-3.0, up to +-3's end (5.05 at B = 3): the field spans
    # +-8 (its half's middle 4; the ions bent harder towards their own cutoff, ~14).
    # The goals (ions 5.5 to 8, leak <= 0.05) leave a diagonal strip, B roughly in proportion
    # to sqrt(V): at the reference (30, 1.9) V ~ 27-32 and B ~ 1.66-2.55; full voltage and
    # field (~8.1) overshoot it, the sliders' middles fall short. The leak the description
    # promises is within both sliders (the owner, 2026-10-13: whatever a level says can
    # happen must be possible): at 30 V the electrons pour into the rod with the field off
    # (17.4 at B = 0, 0.76 at 1.2, 0.006 at 1.5), and at B = 1.9 above ~50 V (0.066 at 50,
    # 0.25 at 60, 3.6 at 80): hence the supply's +-80. A grounded envelope (a thin ring
    # around the cathode) closes the tube: without it the open 2D boundary left a field
    # outside the cathode, which drew electrons off its outer face into an arc drifting
    # round it (E x B) off the board.
    return level(
        "Ion diode",
        "An ion diode turned inside out: the ring is the cathode, emitting electrons from "
        "its inner face, and the rod in its middle, on your supply, emits ions. Without a "
        "magnetic field the electrons rush to the rod and waste the supply's current. Your "
        "solenoid (the yellow coil, top right) curls them: above Hull's critical field they "
        "cannot reach the rod, and stay in a sheath along the ring. The heavy ions still "
        "cross, and the sheath's negative charge pulls more of them across. Deliver an ion "
        "current between 5 and 8.5 to the ring while at most 0.05 of electrons leak to the "
        "rod (averaged over $t = 20$ to $30$). The critical field grows with the voltage, "
        "as $B \\propto \\sqrt{V}$: too weak a field for your voltage lets the electrons "
        "leak; the sheath helps the ions most just above it. The thin outer ring is the "
        "tube's grounded envelope.",
        grid=(24, 16), c=None, t_max=100.0,
        electrodes=[plate(12, 8, 10, thickness=1.0, shape="ring"),
                    plate(12, 8, 2, thickness=2.0, shape="rod", tunable=True),
                    plate(12, 8, 14.5, thickness=0.5, shape="ring")],
        supplies=slider(80.0, 27.0, 30.0, 33.0),
        solenoid_fields=slider(8.0, 1.7, 1.9, 2.2),
        tube=tube(0, 1, 0.0, 0.05, 20.0, 30.0, solenoid=(23, 15), step=0.05,
                  ions=(1, 0.0625), particles=8.0, emitter_mesh=2.0,
                  goals=[current_goal(0, 5.0, 8.5, 20.0, 30.0, carriers="ions")]),
        reference=[supply(12, 8, 30.0), supply(23, 15, 1.9)])


def rectifier():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11): a
    # cathode between two anodes 6 away, each on an AC supply (period 20) whose setting is
    # its amplitude, the right one's half a period behind; the cathode is the output,
    # through R = 30 and C = 1 to ground, and charges positive as it loses electrons.
    # Output over t = 60-100 (preview): both anodes at 20/30/40: mean 6.0/9.9/13.9, swing
    # 0.34/0.55/0.77 (full wave); one at 20/30/40: 3.7/6.2/9.0, swing 1.7/2.6/3.5 (half
    # wave); both, but of opposite signs (in phase): 4.5/7.4/10.5, swing 3.3/5.0/6.7. The
    # goal (mean 8.5 to 11.5, ripple below 1.2, over t = 40-80, after the output has
    # settled) needs both anodes in opposite phase: each supply alone, the other at the
    # reference, solves on 24 to 36 (preview; `generator window`), the verification's
    # error (+-0.7 on the mean) narrowing it (at 30/30: 9.26 +- 0.66, at the edge). The
    # reference 32 sits in the middle.
    omega = 2.0 * math.pi / 20.0
    return level(
        "Rectifier",
        "A rectifier tube: a hot cathode in the middle, an anode on each side. Each anode is "
        "on its own AC power supply, a sine of period 20 whose amplitude you set; the right "
        "one's runs half a period behind the left one's. Electrons flow from the cathode to "
        "an anode only while that anode is positive: the tube conducts one way. The cathode "
        "is connected to ground through a load resistor (30) and a capacitor (1): losing "
        "electrons, it charges positive, and its potential is the rectifier's output, DC "
        "from AC. Make a smooth supply: the cathode's mean potential between 8.5 and 11.5, "
        "and its ripple (half its peak-to-peak) below 1.2, over $t = 40$ to $80$.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(12, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(5, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True),
                    plate(19, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(80.0, 16.0, 32.0, 48.0),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0,
                  circuits=[tube_circuit(1, sine(30.0, omega), resistance=0.0),
                            tube_circuit(2, sine(30.0, omega, phase=math.pi), resistance=0.0),
                            tube_circuit(0, dc(0.0), resistance=30.0, capacitance=1.0)],
                  potentials=[potential_goal(0, "mean", 8.5, 11.5, 40.0, 80.0),
                              potential_goal(0, "swing", 0.0, 1.2, 40.0, 80.0)]),
        reference=[supply(5, 8, 32.0), supply(19, 8, 32.0)])


def triode():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11). A
    # cathode 6 long, five grid wires (0.6 x 0.5, 2 apart) 1.25 before its face on one
    # sine source (amplitude 1, period 40) whose offset, the bias, the player sets, and an
    # anode 6 long through R = 60 to B+ = 100. The anode's own capacitance (it faces the
    # grounded envelope and the grid) shunts the load: g_m/(omega C) capped the gain at ~1.3
    # with a 14-long anode at period 20, so the anode is short and the signal slow. Output
    # amplitude (verification, over one period t = 40-80, the particles softened at an
    # eighth of a segment): 1.96 at bias -1.2, 2.37 at -1, 2.76 at -0.8, 3.12 at -0.6,
    # 3.49 at -0.4, 3.79 at -0.2, 4.09 at 0, the preview ~7 % above (bars +-0.3 to
    # +-0.36): the goal 2.2 to 3.8 is verified for a bias of about -0.9 to -0.4 (at 0 and
    # -0.2 the gain is too high, at -1 too low); cut off below about -2 (only the grid's
    # capacitive feedthrough, 0.32, is left).
    omega = 2.0 * math.pi / 40.0
    wires = [plate(7, y, 0.6, thickness=0.5, angle_deg=90.0, tunable=(y == 8))
             for y in (8, 4, 6, 10, 12)]
    return level(
        "Triode",
        "A triode: between the hot cathode (left) and the anode (right) a grid of five thin "
        "wires, all on one supply: a small signal, a sine of amplitude 1 and period 40, "
        "around a steady potential, the bias, which you set. Close to the cathode, where the "
        "electrons are slow, the grid controls the current far more strongly than the "
        "anode does. The anode takes its current through a resistor (60) from a fixed 100, "
        "so as the current swings, the anode's potential swings too, and much more than "
        "the grid's: an amplifier. Bias the grid so that the anode's potential oscillates "
        "with an amplitude between 2.2 and 3.8 (over $t = 40$ to $80$). Too negative a "
        "bias cuts the current off.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(5, 8, 6, thickness=1.0, angle_deg=90.0),
                    plate(12, 8, 6, thickness=1.0, angle_deg=90.0)] + wires,
        supplies=slider(2.6, -1.4, -0.7, 0.7),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, emit_toward_deg=0.0,
                  circuits=[tube_circuit(2, sine(1.0, omega), resistance=0.0,
                                         supply_offset=True, linked=(3, 4, 5, 6)),
                            tube_circuit(1, dc(100.0), resistance=60.0)],
                  potentials=[potential_goal(1, "amplitude", 2.2, 3.8, 40.0, 80.0,
                                             omega=omega)]),
        reference=[supply(7, 8, -0.7)])


def klystron():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11). A weak beam (a temperature-
    # limited cathode, 0.02 per unit length of its 4-long face, 400 macroparticles per unit
    # time: 1/400 of a period apart) from a Wehnelt gun (the cathode between two plates at
    # its potential, the accelerating gap into the grounded drift tube's mouth a
    # converging lens: most of the beam reaches the collector) flies past a buncher
    # (plates 14-20 on a sine of amplitude V1, period 1: its entrance and exit gaps
    # modulate the electrons' speeds) and a catcher 21 further (one plate, a tank: L = 0.051
    # in series with R = 0.04 to ground, resonant at omega = 2 pi with the plate's own
    # capacitance, 0.497 from the boundary elements; Q ~ 8, rung up in 2.5). The player
    # sets the cathode's potential -V0 (the beam's energy) and V1. The beam's space charge
    # is negligible: the bunching is ballistic. Output amplitude against the drive at
    # V0 = 50 (preview): 0.073 at V1 = 8, 0.111 at 14, 0.121 at 20, 0.109 at 26
    # (overbunching); with the buncher at a fixed offset the J1 shape was measured to 2 %
    # (J1(X)/J1max with X = 0.092 V1 fitted at the peak: 0.306, 0.711, 0.956 at V1 = 4,
    # 10, 24 for the measured 0.300, 0.694, 0.959). The output's phase moves with the
    # drift angle omega L/v0, 7.5 degrees per volt near 50 (theta0/(2 V0) = 7.6: phase
    # pushing): 122.7 at V0 = 40, 169.6 at 45, -151.4 at 50, -114.4 at 55, -82.7 at 60
    # (V1 = 20); the amplitude hardly depends on V0. At V0 = 54: the amplitude 0.072 at
    # V1 = 8, 0.088 at 10, 0.095 at 11, 0.102 at 12, 0.111 at 14 (phase -125); past the
    # optimum 0.103 at 28, 0.096 at 30, 0.089 at 32, 0.079 at 34 (phase -113 to -103).
    # Refine 1 to 2: phases within 2.6 degrees, amplitudes 5 %. A flat optimum and those
    # bars leave no narrow window at saturation, so the goal is a klystron backed off below
    # it (as one feeding an accelerator is, so that its output follows its drive):
    # amplitude 0.085-0.105 and phase -135 to -105, which ask for V0 ~ 52-56 and V1 ~
    # 9.7-12.5, or overbunched ~27.5-33; the next island of phase lies beyond the slider.
    omega = 2.0 * math.pi
    P = plate
    return level(
        "Klystron",
        "A klystron amplifier. The gun on the left shoots a weak electron beam into a "
        "grounded tube; you set the cathode's voltage, and with it the beam's speed. On its "
        "way the beam passes the buncher, plates on a small sine (period 1) whose amplitude "
        "you set: crossing its gaps speeds some electrons up and slows others down. Further "
        "on the fast ones catch up with the slow ones and the beam gathers into bunches, one "
        "per period, which induce a current in the catcher plate, a tuned circuit ringing at "
        "the same frequency. Run it as a linear amplifier, below saturation: the catcher "
        "must oscillate with an amplitude between 0.085 and 0.105, and with a phase, "
        "relative to the buncher's signal, between $-135°$ and $-105°$ (over $t = 10$ to "
        "$20$). A stronger signal bunches the beam tighter, up to a point: too strong a one "
        "overbunches it, the fast electrons overtaking the slow ones before the catcher. "
        "The phase depends on how long the bunches take to arrive.",
        grid=(48, 20), c=None, t_max=40.0,
        electrodes=[
            P(3, 10, 4, thickness=1.0, angle_deg=90.0, tunable=True),
            P(3, 13, 4, thickness=1.0), P(3, 7, 4, thickness=1.0),
            P(10, 13, 6, thickness=1.0), P(10, 7, 6, thickness=1.0),
            P(17, 13, 6, thickness=1.0, tunable=True), P(17, 7, 6, thickness=1.0),
            P(28, 13, 14, thickness=1.0), P(28, 7, 14, thickness=1.0),
            P(38, 13, 4, thickness=1.0),
            P(38, 7, 4, thickness=1.0),
            P(43, 13, 4, thickness=1.0), P(43, 7, 4, thickness=1.0),
            P(47, 10, 4, thickness=1.0, angle_deg=90.0)],
        supplies=slider(80.0, -30.0, -54.0, -60.0, 11.0, 30.0),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, emit_toward_deg=0.0, saturation=0.02,
                  particle_rate=400.0,
                  circuits=[tube_circuit(0, dc(-50.0), resistance=0.0, linked=(1, 2)),
                            tube_circuit(5, sine(11.0, omega), resistance=0.0, linked=(6,)),
                            tube_circuit(10, dc(0.0), resistance=0.04, inductance=0.051)],
                  potentials=[potential_goal(10, "amplitude", 0.085, 0.105, 10.0, 20.0,
                                             omega=omega),
                              potential_goal(10, "phase", -135.0, -105.0, 10.0, 20.0,
                                             omega=omega)]),
        reference=[supply(3, 10, -54.0), supply(17, 13, 11.0)])



def magnetron():
    # Designed with `generator tube-scan --trace` (PHYSICS.md 2.11; docs/TUBES.md). A strapped
    # eight-vane magnetron in the z-invariant model: the cathode a rod (radius a = 2) on the
    # player's supply, eight vanes (tips at b = 4.5; centres on grid nodes, so the diagonal
    # ones are shorter) and the anode block, a ring (inner radius 8.5). The vanes at 0, 90,
    # 180 and 270 degrees and the block are grounded; the diagonal ones are strapped onto one
    # node and a tank to ground (L = 0.036 in series with R = 0.002, and 1.34 to ground: with
    # the vanes' own 1.337 a cold resonance at 3.22, Q ~ 55; the added capacitance halves the
    # electrons' pull on the frequency, which without it hopped between branches near 3.0 and
    # 3.6). Its oscillation is the pi mode: four crests around the vanes, turning at omega/4.
    # The oscillation starts from the switch-on transient (the cathode's step rings the tank
    # to ~4) and the particles' noise, so its start-up time depends on the resolution (t ~ 20
    # to 80); its saturated swing does not (45 V, B = 6: preview 14.6, verification 15.6, x4
    # particles 16.0): the goal's window is late, after every resolution has started. An
    # injected priming signal did not make the start reliable (below threshold the hub
    # absorbs it) and is not used. Swing over t = 90-120 (preview; frequency from the
    # up-crossings 3.03-3.11 where it oscillates):
    #   V \ B   5.5    6      6.5
    #   37.5    0.57   0.55   0.32
    #   40      8.50   0.68   0.36
    #   42.5    3.13   9.03   0.41
    #   45      0.80  14.86   1.01
    #   47.5    5.16* 14.55  15.99
    #   50      4.20* 3.47   21.51
    #   52.5    1.53   1.61  22.84
    #   55      0.82   1.50   5.47   (* the other branch, near 3.6)
    # The oscillating band follows the Buneman-Hartree condition, the electrons' drift in
    # step with the wave, V = Omega B (b'^2 - a^2)/2 - Omega^2 b'^2/2 with Omega = omega/4 =
    # 0.77 and an effective anode radius b' ~ 5 (between the vanes' tips and their roots):
    # 37, 41 and 45 at B = 5.5, 6 and 6.5, the band's lower edges. Far above it the electrons
    # outrun the wave and the oscillation stops (B = 6: none at 50, below Hull's cutoff 58.7);
    # below Hull's critical field they cross to the vanes (at 45 V: 34 into the vanes at
    # B = 4.4; the cutoff is ~4.5, Hull's with b' = 5). Just above it a second solution:
    # at 45 V, B = 4.58-4.67 the tank swings 13.6 (omega 3.29) while 8.4 flows into the
    # vanes, seven times the main band's current (B = 4.8-5.25: the other branch, 5.8-7.8,
    # short of the goal). `generator window` at the reference: the supply 43.3-48.0 (2.9 %
    # of the range), the solenoid 5.83-6.25 and that pocket (and their mirrors: the field's
    # sign only turns the hub the other way); every middle verified.
    # Its output (the owner, 2026-10-13: the radiation on the map, from all the charges): the
    # tube moved to the middle of a 96 x 96 board (about a wavelength, 2 pi c/omega ~ 82 at
    # c = 40, across) and an antenna, a rod of diameter 2 twenty cells from the axis, is
    # strapped to the tank's vanes; the field solver (level::tube_em) shows its waves. Its
    # radiation resistance, measured with `generator tube-field`, is the tank's load: at R =
    # 0.002 the board radiated 0.91 of the load's R<I^2>, at R = 0.0018 the two agree to ~3 %
    # (test magnetron_antenna_radiates_its_load). Redesigned windows: the supply 42.0-47.3,
    # the solenoid 5.83-6.33 (the near-cutoff pocket shrank to a point at 4.58).
    return level(
        "Magnetron",
        "A magnetron, the heart of a microwave oven and of wartime radar. Electrons leave "
        "the cathode in the middle (your supply sets its voltage below the vanes) and, in "
        "your solenoid's field (the yellow coil, top right), curl around it instead of "
        "crossing to the vanes. The vanes are strapped in two sets: one grounded with the "
        "block, the other through a tuned circuit ringing near $\\omega = 3$. Its "
        "oscillation is a wave of four crests running around the vanes at $\\omega/4$; "
        "electrons drifting in step with it give it energy and gather into four spokes "
        "that turn with it. Make the tuned circuit oscillate with a swing of at least 10 "
        "(half its peak-to-peak, over $t = 90$ to $120$). The electrons drift at $E/B$: "
        "to keep pace with the wave the voltage must rise with the field. Too weak a field "
        "lets them straight onto the vanes; too fast a drift outruns the wave. The rod "
        "outside the block is the antenna, strapped to the tuned set of vanes: it radiates "
        "the microwaves: the magnetic or energy flow map shows them spreading over the "
        "board once the oscillation has grown out of the noise, after about $t = 50$ (the "
        "mouse wheel zooms in on the spokes, a drag pans).",
        grid=(96, 96), c=None, t_max=150.0,
        electrodes=[plate(48, 48, 4, thickness=4.0, shape="rod", tunable=True),
                    plate(54, 48, 3.0, thickness=0.6, angle_deg=0.0),
                    plate(52, 52, 2.313708, thickness=0.6, angle_deg=45.0),
                    plate(48, 54, 3.0, thickness=0.6, angle_deg=90.0),
                    plate(44, 52, 2.313708, thickness=0.6, angle_deg=135.0),
                    plate(42, 48, 3.0, thickness=0.6, angle_deg=180.0),
                    plate(44, 44, 2.313708, thickness=0.6, angle_deg=225.0),
                    plate(48, 42, 3.0, thickness=0.6, angle_deg=270.0),
                    plate(52, 44, 2.313708, thickness=0.6, angle_deg=315.0),
                    plate(48, 48, 18, thickness=1.0, shape="ring"),
                    plate(68, 48, 2, thickness=2.0, shape="rod")],
        supplies=slider(80.0, -38.0, -45.0, -52.0),
        solenoid_fields=slider(10.0, 5.5, 6.0, 6.5),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, solenoid=(95, 95), step=0.05,
                  circuits=[tube_circuit(1, dc(0.0), resistance=0.0, linked=(3, 5, 7, 9)),
                            tube_circuit(2, dc(0.0), resistance=0.0018, inductance=0.036,
                                         capacitance=1.34, linked=(4, 6, 8, 10))],
                  potentials=[potential_goal(2, "swing", 10.0, 1000.0, 90.0, 120.0)],
                  output=(1, 40.0)),
        reference=[supply(48, 48, -45.0), supply(95, 95, 6.0)])

def ohms_law():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11, the Drude metal): a strip
    # 16 x 4 between two contacts. Preview currents into the right contact over t = 5-35
    # (the contacts emitting as thermal walls): 0.91 at 8, 1.72 at 16, 2.43 at 24, 3.17 at
    # 32, 3.84 at 40, -2.61 at -24: in proportion to within a slow fall of I/V (0.1135 to
    # 0.096: the carriers warm in the field); the bulk theory gives 0.1175 V. The goal
    # 2.6-3.2: the supply solves on 26-32 (`generator window`).
    return level(
        "Ohm's law",
        "Inside a metal. The rectangle is a conductor seen end-on: its electrons (blue; a "
        "few traced with their paths) fly at their thermal speed and bounce off the "
        "metal's ions (the dots), a zigzag. Your supply, on the contact at the right end, "
        "adds a slow drift along the field to that random motion: a current, in "
        "proportion to the voltage, $I = V/R$ (Ohm's law). Make a current between 2.6 and "
        "3.2 flow into the right contact (averaged over $t = 10$ to $40$). At every "
        "collision the electrons hand the ions the energy the field gave them: the Joule "
        "heat, shown below. A classical model of a metal (its real electrons are a "
        "quantum gas): see the physics notes.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(3, 8, 4), contact(21, 8, 4, tunable=True)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  goals=[current_goal(1, 2.6, 3.2, 10.0, 40.0, carriers="net")],
                  metals=[metal((4, 6), (20, 10), [0, 1], seed=1, **METAL)]),
        reference=[supply(21, 8, 30.0)])


def long_and_thin():
    # Designed with `generator tube-scan`: two strips of the same metal in series, the upper
    # one 16 x 4 (L/W = 4), the lower one 8 x 8 (L/W = 1). Supply 1 on the upper strip's
    # left end, supply 2 on the junction (the upper strip's right end, strapped to the lower
    # strip's left end), the lower strip's right end grounded. Each strip alone (the
    # contacts emitting as thermal walls): the upper 2.63 at 26.5, the lower 2.32 at 5.5:
    # R = 10.1 and 2.37, a ratio of 4.3 (L/W gives 4; the thin strip's stronger field warms
    # its carriers). The goals 2.0-2.8 each: at the reference (32, 5.5) supply 1 solves on
    # 28-38, the junction on about 4.7-6.6 (the textbook ratio, V2 = V1/5 = 6.4, passes).
    return level(
        "Long and thin",
        "Two strips of the same metal: the upper one long and thin, the lower one short and "
        "wide. A strip's resistance grows with its length and falls with its width, "
        "$R = L/(\\sigma W)$ (per unit length along it): the upper one has about four times "
        "the lower one's. They are in series: your first supply holds the upper strip's "
        "left end, your second the junction (the upper strip's right end, wired to the lower "
        "strip's left end), and the lower strip's right end is grounded. Make a current "
        "between 2.0 and 2.8 flow through each strip (into its left end, averaged over "
        "$t = 10$ to $40$).",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(17, 5, 8), contact(3, 12, 4, tunable=True),
                    contact(21, 12, 4, tunable=True), contact(7, 5, 8)],
        supplies=slider(40.0, 5.5, 15.0, 32.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  circuits=[tube_circuit(2, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 2.0, 2.8, 10.0, 40.0, carriers="net"),
                         current_goal(3, 2.0, 2.8, 10.0, 40.0, carriers="net")],
                  metals=[metal((4, 10), (20, 14), [1, 2], seed=4, **METAL),
                          metal((8, 1), (16, 9), [3, 0], seed=5, **METAL)]),
        reference=[supply(3, 12, 32.0), supply(21, 12, 5.5)])


def potentiometer():
    # Designed with `generator tube-scan`: a strip 18 x 4, its left end grounded, a side
    # contact on its top a third of the way along (x 8-10). Preview with the right end at
    # 30: the side contact's current (positive in) -0.41 at 8, -0.16 at 9, 0.17 at 10,
    # 0.45 at 11, 0.65 at 12: none at about 9.5 (a uniform field gives 10; the ohmic
    # contacts' layers take a little). The goals: 2.4-3.0 along the strip (the right end
    # at about 27-33) and -0.35..0.35 through the side contact (about 8.3-10.6 at 30;
    # +-0.25 left the verification uncertain: the side contact's current carries the
    # strip's thermal noise, an error bar of ~0.18 with 50 macroparticles per area over
    # t = 10-40; now 100 over 10-80). With the contacts emitting as a thermal wall (they
    # had been heat sinks, PHYSICS.md 2.11) the side contact's current at 30 is -0.22 at 8,
    # -0.02 at 9, 0.31 at 10: none at about 9.05 (V1/3 = 10: the contacts take ~1 V).
    # The side contact's current verified with a bar of ~0.3 (+-0.35 left it uncertain):
    # the window +-0.5, about 7.5-10.6 at 30. The showcase rule's review (2026-10-20): it solved
    # on 6.75-10.25, 4.4 % of +-40: the window +-0.7 (the side current changes by ~0.27 per
    # volt, still a null to find).
    return level(
        "Potentiometer",
        "A longer strip of metal with a third contact on its side, a third of the way "
        "along, on your second supply. Inside a conductor carrying a steady current the "
        "field is uniform: the potential falls evenly from one end to the other (the "
        "contacts take a little). Make a "
        "current between 2.4 and 3.0 flow along the strip into its right end, and set the "
        "side contact so that no current flows through it (between $-0.7$ and $0.7$, "
        "averaged over $t = 10$ to $80$): it must sit at the potential the strip has "
        "there. A potentiometer measures a voltage this way, drawing no current from it.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(2, 8, 4), contact(22, 8, 4, tunable=True),
                    contact(9, 11, 2, vertical=False, tunable=True)],
        supplies=slider(40.0, 10.0, 15.0, 30.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  goals=[current_goal(1, 2.4, 3.0, 10.0, 80.0, carriers="net"),
                         current_goal(2, -0.7, 0.7, 10.0, 80.0, carriers="net")],
                  metals=[metal((3, 6), (21, 10), [0, 1, 2], seed=2,
                                **dict(METAL, carriers_per_area=100.0))]),
        reference=[supply(22, 8, 30.0), supply(9, 11, 9.0)])


def hall_probe():
    # Designed with `generator tube-scan` (100 macroparticles per area: the floating
    # contacts' voltage carries the resistor's thermal noise, ~0.1 over the window): a
    # strip 16 x 4, the right end on the supply, two floating contacts on its sides at
    # mid-length (C = 1 to ground). At 30 (I ~3.1) the top's potential against the
    # bottom's: -0.58 at B = 0.4, -1.15 at 0.7, -1.44 at 1.0, -2.1 at 1.2, -2.33 at 1.6,
    # +0.75 at -0.7 (electrons pushed to the top; about -0.15 at B = 0). Below the weak-field
    # R_H I B = (4/pi) I B/rho as omega_c tau grows (0.4 B: the Hall factor falls towards 1)
    # and the contacts short the strip's ends. The goals: 2.7-3.5 into the right end (about
    # 26-34) and the top 1.0-1.9 below the bottom (B about 0.6-1.15 at 30). With the
    # contacts emitting as a thermal wall: -0.90 at B = 0.9, -1.44 at 1.2, -1.68 at 1.5,
    # +0.22 at 0. `generator window` at the reference (30, 1.3): the supply solves on
    # 26-34, the solenoid on 1.0-1.4 (its slider +-4, the middles at +-2), over t = 10-80
    # (the floating contacts' noise).
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at 30 (target/scratch/b_range_scan.py; preview): current 3.14, 3.12, 3.06,
    # 3.00, 2.91, 2.73, 2.55, 2.21, 1.96 and the top against the bottom -0.45, -1.05, -2.43,
    # -3.68, -4.67, -6.61, -8.12, -10.57, -12.52 at B = 0.5, 1, 2, 3, 4, 6, 8, 12, 16 (the
    # Hall voltage ever more sublinear). The goal moved to -6.0 to -3.5: B about 2.9-5.5.
    hall = dict(METAL, carriers_per_area=100.0)
    return level(
        "Hall probe",
        "A strip of metal carrying a current, in your solenoid's field (B along the "
        "strip's length, out of the board when positive; the yellow coil, top right). The "
        "field pushes the moving electrons sideways (the Lorentz force) until the charge "
        "they pile on one side holds the rest back: a voltage across the strip, the Hall "
        "voltage, $V_H \\approx I B/(n q)$ (per unit length). The two contacts on the "
        "strip's sides read it, floating: no current flows through them. Try the field's "
        "whole range: weak, the electrons' paths barely bend; strong, they curl into small "
        "circles between collisions. Make a current between 2.7 and 3.5 flow into the right "
        "end and the top contact sit between 3.5 and 6 below the bottom one (averaged over "
        "$t = 10$ to $80$). Which way must the field point? Electrons are negative.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(3, 8, 4), contact(21, 8, 4, tunable=True),
                    contact(12, 11, 2, vertical=False), contact(12, 5, 2, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -4.2, 3.0, 4.2, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 15),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0)],
                  goals=[current_goal(1, 2.7, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", -6.0, -3.5, 10.0, 80.0,
                                             reference=3)],
                  metals=[metal((4, 6), (20, 10), [0, 1, 2, 3], seed=3, **hall)]),
        reference=[supply(21, 8, 30.0), supply(23, 15, 4.2)])


def positive_carriers():
    # Designed with `generator tube-scan`: two strips 16 x 4 in the player's solenoid, the
    # upper of the levels' metal (electrons), the lower the same but its carriers positive
    # (q/m = +1), each with its left end grounded, its right end on a supply and floating
    # side contacts at mid-length (C = 1). The Hall voltage's sign follows the carriers' for
    # a given current: with the carriers moving the same way the two strips' Hall voltages
    # agree in sign, so the tops both sit above the bottoms only with opposite supplies.
    # Preview at (30, B = -1.3, -15) over t = 10-80: upper current 3.10, top above bottom
    # 1.52; lower current 1.66, top above bottom 1.21 (both 0.3 lower at B = 0 for the upper
    # strip, 0 for the lower; the lower's potential shifts the upper's reading by ~0.2). The
    # floating contacts' noise over the window: 0.11 (upper) and 0.15 (lower) between seeds
    # and resolutions (50 macroparticles per area: two strips at 100 would pass the
    # verification's carrier limit). At B = -1.3 the lower Hall voltage at -12 (0.65-0.84)
    # did not verify against its noise; at B = -1.8: upper 2.23 +- 0.10 at 30, lower 0.89
    # at -10 and 1.23 at -14 (two seeds, both resolutions). `generator window`: the upper
    # supply solves on 28-34, the lower on -16 to -8, the solenoid on -2.5 to -1.5: its
    # slider is +-7 (at +-5 the half's middle, -2.5, solved).
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at (30, B, -12) (preview): upper current 3.11, 3.06, 3.06, 2.98, 2.85, 2.68,
    # 2.50, 2.22, 1.95; upper Hall 0.40, 1.12, 2.31, 3.57, 4.69, 6.12, 7.63, 10.47, 11.99;
    # lower 0.57, 0.73, 1.14, 1.47, 2.05, 2.80, 3.15, 4.44, 4.92 at B = -0.5, -1, -2, -3,
    # -4, -6, -8, -12, -16. The goals moved: current 2.6-3.5, upper 3.4-6.0, lower 1.3-2.6:
    # B about -5.6 to -2.9.
    material = dict(METAL, carriers_per_area=50.0)
    return level(
        "Positive carriers",
        "Two strips of different metals in your solenoid's field, each with two floating "
        "contacts on its sides that read its Hall voltage. The upper one conducts by "
        "electrons; the lower one's carriers are positive (orange). Some real metals, such "
        "as zinc and beryllium, behave so: band theory explains it, and in this classical "
        "model their carriers simply carry a positive charge. The field reaches strengths "
        "that curl the carriers into small circles. Make a current between 2.6 and 3.5 flow "
        "into the upper strip's right end, and make both strips' top contacts sit above "
        "their bottom ones: the upper by 3.4 to 6, the lower by 1.3 to 2.6 (averaged over "
        "$t = 10$ to $80$). Which way must each current flow?",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[contact(3, 14, 4), contact(21, 14, 4, tunable=True),
                    contact(12, 17, 2, vertical=False), contact(12, 11, 2, vertical=False),
                    contact(3, 5, 4), contact(21, 5, 4, tunable=True),
                    contact(12, 8, 2, vertical=False), contact(12, 2, 2, vertical=False)],
        supplies=slider(40.0, -12.0, 15.0, 30.0),
        solenoid_fields=slider(16.0, -4.0, 3.0, 4.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(6, dc(0.0), capacitance=1.0),
                            tube_circuit(7, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), resistance=0.0)],
                  goals=[current_goal(1, 2.6, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", 3.4, 6.0, 10.0, 80.0, reference=3),
                              potential_goal(6, "mean", 1.3, 2.6, 10.0, 80.0, reference=7)],
                  metals=[metal((4, 12), (20, 16), [0, 1, 2, 3], seed=6, **material),
                          metal((4, 3), (20, 7), [4, 5, 6, 7], seed=7, q_over_m=1.0,
                                **material)]),
        reference=[supply(21, 14, 30.0), supply(23, 19, -4.0), supply(21, 5, -12.0)])


# Why Hall probes are semiconductors (the owner, 2026-10-20): two strips of the same size and
# the same conductance in one field, the upper the levels' metal (few carriers that fly far: it
# resembles n-type germanium at 1e16 per cm^3 in these fields, game scales.rs), the lower four
# times as many carriers colliding four times as often (density 8, mean free path 0.25: sigma =
# rho tau the same; n-type silicon at 1e17 per cm^3, within x1.21). At the same current the
# lower's carriers drift four times slower, the magnetic push on each is four times weaker, and
# so is its Hall voltage (V_H = IB/(nq)); at the same field they turn a quarter as far between
# collisions (omega_c tau). Copper is the trend a million times further. The step follows the
# dense strip's plasma oscillation (omega_p = 10.0: 0.019, omega_p dt = 0.19).
DENSE = dict(density=8.0, thermal_speed=2.0, mean_free_path=0.25, core_radius=0.002,
             carriers_per_area=50.0)
DENSE_STEP = 0.019


def few_or_many():
    # Designed with `generator tube-scan` (target/scratch/few_scratch.py): at (30, B, 30)
    # (preview, t = 10-80) over B = 0, 1, 2, 3, 4.2, 6, 9, 12, 16: upper current 3.18, 3.16,
    # 3.14, 3.03, 2.93, 2.77, 2.53, 2.29, 2.04; lower 3.25, 3.40, 3.33, 3.29, 3.25, 3.10, 3.13,
    # 2.94, 2.68; upper Hall (top less bottom) -0.19, -1.58, -3.13, -3.91, -5.22, -6.73, -9.27,
    # -10.7, -13.3; lower +0.76, +0.02, -0.81, -0.85, -1.45, -1.80, -2.85, -4.00, -4.97 (its
    # floating contacts scatter ~0.4: its macroparticles are four times heavier). Per unit IB
    # the Hall voltages are 0.42 and 0.11: four times apart, as the densities. Verified at
    # (30, 3.8, 30): currents 3.00 +- 0.05 and 3.32 +- 0.16, Hall -4.69 +- 0.08 and -1.37 +-
    # 0.19. The Hall readings' error bars grow with B (upper 0.08, 0.51, 0.74 and lower 0.29,
    # 0.36, 0.70 at B = 2.5, 5, 5.5): with goals -6 to -3.5 and -2 to 0 the solenoid solved on
    # 3.0-4.5 only (4.7 %), so -7 to -3 and -2.5 to 0. `generator window`: the solenoid
    # solves on 2.0-6.0 (12.5 % of +-16; 1.5 and 6.5-8.5 fail), the upper supply on 28-35
    # (8.8 % of +-40), the lower on 25-31 (7.5 %; both measured with the narrower Hall goals,
    # their edges set by the current goals). The solenoid's +-16 spans straight paths to
    # circles (omega_c tau 6.4 in the upper strip, 1.6 in the lower).
    material = dict(METAL, carriers_per_area=50.0)
    return level(
        "Few carriers or many",
        "Two strips of the same size and the same resistance in your solenoid's field, each "
        "with two floating contacts on its sides that read its Hall voltage. The upper one's "
        "carriers are few and fly far between collisions (like germanium's); the lower one has "
        "four times as many, colliding four times as often (like silicon's). Make a current "
        "between 2.7 and 3.5 flow into each strip's right end, and the upper strip's top "
        "contact sit between 3 and 7 below its bottom one, while the lower strip's sits at "
        "most 2.5 below (averaged over $t = 10$ to $80$). Carrying the same current, the many "
        "carriers drift slowly, the field pushes each one gently, and the Hall voltage is "
        "small: real Hall probes are semiconductors with few, fast carriers (InSb, GaAs), "
        "while copper's carriers are about a million times as many, and its Hall voltage a "
        "million times smaller: microvolts.",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[contact(3, 14, 4), contact(21, 14, 4, tunable=True),
                    contact(12, 17, 2, vertical=False), contact(12, 11, 2, vertical=False),
                    contact(3, 5, 4), contact(21, 5, 4, tunable=True),
                    contact(12, 8, 2, vertical=False), contact(12, 2, 2, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -3.8, 3.0, 3.8, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=DENSE_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(6, dc(0.0), capacitance=1.0),
                            tube_circuit(7, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), resistance=0.0)],
                  goals=[current_goal(1, 2.7, 3.5, 10.0, 80.0, carriers="net"),
                         current_goal(5, 2.7, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", -7.0, -3.0, 10.0, 80.0, reference=3),
                              potential_goal(6, "mean", -2.5, 0.0, 10.0, 80.0, reference=7)],
                  metals=[metal((4, 12), (20, 16), [0, 1, 2, 3], seed=6, **material),
                          metal((4, 3), (20, 7), [4, 5, 6, 7], seed=7, **DENSE)]),
        reference=[supply(21, 14, 30.0), supply(23, 19, 3.8), supply(21, 5, 30.0)])


def magnetoresistance():
    # Designed with `generator tube-scan` (target/scratch/design_mr_intro.py): a strip 16 x 4
    # (current along it) and one 4 x 10 (short along the current, wide across it), both on
    # one supply, in the player's solenoid. In the long strip the Hall field holds the
    # carriers back from the sides and they flow straight; the short wide strip's contacts
    # short that field: the carriers' paths bend and its resistance grows with |B|
    # (geometric magnetoresistance; with L/W = 0.4 short of the Corbino limit 1 + (0.4 B)^2).
    # Preview over t = 10-80, the currents into the two right contacts: at 7 V the long
    # strip 0.73, 0.69, 0.68, 0.63, 0.63 at B = 3..7 (its own, intrinsic magnetoresistance:
    # tau ~ 1/v), the short one 5.43, 4.62, 4.08, 3.64, 3.04 (about 7.2 at B = 0); at 6 V
    # and 8 V in proportion. A run takes ~16 s, at B = 6 as at 0. `generator window`: B solves on
    # 4.5-5.5 and -5.5 to -4.5 (its mirror's middle, -5, not verified: near the goals' edges);
    # with the long strip's goal 0.62-0.76 the supply solved on 7.0-7.5 only (1.2 % of +-20):
    # widened to 0.58-0.78, the slider +-10.
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at 7 V (preview): long strip 0.795, 0.772, 0.738, 0.683, 0.679, 0.626, 0.573,
    # 0.547, 0.522, 0.452; short 6.96, 6.09, 5.41, 4.63, 4.06, 3.61, 2.68, 2.09, 1.56, 1.18 at
    # B = 1, 2, 3, 4, 5, 6, 8, 10, 13, 16. At B = 5 both currents go as the supply (long
    # 0.096 V, short 0.58 V), so its window is the overlap of the two goals' ratios: with
    # 0.58-0.78 and 3.2-4.6 it solved on 6.5-7.5 only (5 %; 6.0 gave 0.550, 8.0 gave 4.68).
    # Widened to 0.56-0.80 and 3.2-4.8: the supply about 5.8-8.2, B about 3.8-6.6 (the half's
    # middle, 8, clear).
    return level(
        "Magnetoresistance",
        "Two strips of the same metal on one supply, in your solenoid's field: a long one and "
        "a short, wide one. In the long strip the charge the field piles on its sides (the "
        "Hall voltage) pushes the carriers back, and they flow straight on. The short strip's "
        "contacts short that charge out: its carriers' paths bend and lengthen, and it "
        "resists more, whichever way the field points; at full strength its carriers circle "
        "on the spot between collisions. Make between 0.56 and 0.8 flow into the long "
        "strip's right end and between 3.2 and 4.8 into the short strip's (averaged over "
        "$t = 10$ to $80$).",
        grid=(24, 20), c=None, t_max=200.0,
        # 0, 1: the long strip's contacts; 2, 3: the short wide strip's.
        electrodes=[contact(3, 15, 4), contact(21, 15, 4, tunable=True),
                    contact(9, 6, 10), contact(15, 6, 10)],
        supplies=slider(10.0, 7.0, 4.0),
        solenoid_fields=slider(16.0, -5.0, 3.0, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(1, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 0.56, 0.80, 10.0, 80.0, carriers="net"),
                         current_goal(3, 3.2, 4.8, 10.0, 80.0, carriers="net")],
                  metals=[metal((4, 13), (20, 17), [0, 1], seed=11, **METAL),
                          metal((10, 1), (14, 11), [2, 3], seed=12, **METAL)]),
        reference=[supply(21, 15, 7.0), supply(23, 19, 5.0)])


def strength_and_sign():
    # The metals' master (task 6): a magnetic-field sensor of two parts in series on one
    # supply, joined at a floating junction. The short wide strip of Magnetoresistance (4 x
    # 10, electrons) resists more as |B| grows, whichever way B points; the Hall bar 16 x 4
    # carries positive carriers (Positive carriers' lower strip), its Hall voltage odd in B.
    # The current pins the supply, the junction's potential (the wide strip's share of it)
    # |B|, and the Hall voltage's sign B's sign.
    # Designed with `generator tube-scan` (target/scratch/level_scratch.py): at 30 over B =
    # -16, -8, -5, -3, 0, 3, 5, 8, 12, 16 (preview, t = 10-80) the current 1.20, 2.02, 2.48,
    # 2.72, 2.98, 2.75, 2.43, 2.05, 1.58, 1.16, the junction 7.44, 5.10, 4.16, 3.37, 2.86, 3.34,
    # 4.06, 5.00, 5.94, 6.48 (both even in B: the two strips' magnetoresistance) and the Hall
    # voltage -10.9, -7.25, -5.43, -3.69, -0.11, 3.67, 5.36, 6.97, 9.34, 10.9 (odd: the top
    # above the bottom for B out of the board, the current then flowing leftwards through
    # positive carriers). Then mirrored left to right (the wide strip under the Hall bar's
    # right end, the wire between them short and clear of the board's legend): the current
    # flows rightwards, and the top sits above the bottom for B into the board; at B = 5,
    # before the mirror, 2.40, 3.98, 5.15. Verified at (30, -5): 2.47 +- 0.06, 4.065 +-
    # 0.005, 5.43 +- 0.19. `generator window`: the supply solves on 27-32 (6.2 % of +-40; 26
    # and 33 fail on the junction), the solenoid on -6.0 to -3.5 (7.8 % of +-16; -6.5 to -9
    # and -3 to -2 fail). The junction over the current is the wide strip's resistance, whatever the
    # supply: 0.96, 1.22, 1.67, 2.44, 3.75, 5.59 at |B| = 0, 3, 5, 8, 12, 16. Goals: the
    # current 2.1-2.9, the junction 3.6-4.5 (their ratio 1.24-2.14: |B| about 3-7, clear of
    # the half's middle, 8), the Hall voltage 3-8 (the sign).
    return level(
        "Strength and sign",
        "A magnetic field sensor in two parts, in series on your supply. The current enters "
        "the Hall bar's left end (its carriers are positive), leaves its right end into the "
        "junction, and crosses a short, wide strip of the electron metal to ground; that "
        "strip resists more as the field grows, whichever way the field points. Make a "
        "current between 2.1 and 2.9 flow through the wide strip, the junction sit between "
        "3.6 and 4.5, and the Hall bar's top contact sit 3 to 8 above its bottom one "
        "(averaged over $t = 10$ to $80$). The wide strip reads the field's strength but not "
        "its sign; the Hall bar reads its sign.",
        grid=(24, 20), c=None, t_max=200.0,
        # 0: the wide strip's bottom contact (grounded); 1: the Hall bar's left end (the
        # supply); 2, 3: its side contacts; 4, 5: the junction (the Hall bar's right end
        # wired to the wide strip's top contact, straight below it).
        electrodes=[contact(17, 2, 10, vertical=False), contact(3, 15, 4, tunable=True),
                    contact(12, 18, 2, vertical=False), contact(12, 12, 2, vertical=False),
                    contact(21, 15, 4), contact(17, 8, 10, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -5.0, 3.0, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), capacitance=1.0, linked=[5])],
                  goals=[current_goal(5, 2.1, 2.9, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(4, "mean", 3.6, 4.5, 10.0, 80.0),
                              potential_goal(2, "mean", 3.0, 8.0, 10.0, 80.0, reference=3)],
                  metals=[metal((4, 13), (20, 17), [1, 2, 3, 4], seed=13, q_over_m=1.0,
                                **METAL),
                          metal((12, 3), (22, 7), [5, 0], seed=14, **METAL)]),
        reference=[supply(3, 15, 30.0), supply(23, 19, -5.0)])


# Inside the components: transistors (docs/DRUDE.md, the owner approved the plan
# 2026-10-20). The field effect: a gate beside a thin channel of few carriers pushes them out
# or draws more in, and the channel's current follows. A metal cannot be gated (it screens
# within an atom); the model's lightly doped channel can (its Debye length 0.56, a quarter of
# its width). The channel's material: density 1, mean free path 0.5 (n-type silicon at 1e16
# per cm^3 by its own ratios, within x1.1: game scales.rs, `resemblance`).
CHANNEL = dict(density=1.0, thermal_speed=2.0, mean_free_path=0.5, core_radius=0.005,
               carriers_per_area=400.0)
CHANNEL_STEP = 0.04


def gate():
    # Designed with `generator tube-scan` (target/scratch/gate_design*.py): a channel 8 x 2
    # between a grounded source and a drain held at 10, gates 1.5 above and below it on the
    # player's supply. The channel charges through its own resistance (R C ~ L^2: a channel 16
    # long had not settled by t = 30; 8 long settles within ~10). Convergence: with the gates
    # 0.5 away and 100 macroparticles per area the preview's current was 0.342 at gate 0
    # against the verification's 0.220 (and 0.243 at refine 3); the macroparticles were the
    # cause (twice as many: 0.262; four times: 0.266; the verification with twice as many:
    # 0.256; the softening and the contacts' segments hardly mattered), and the gates moved
    # 1.5 away (three preview segments). Depletion needs more: the current at -10, -15, -25,
    # -35 was 0.176, 0.163, 0.150, 0.159 at 200 per area (a plateau, an artefact), 0.182, 0.136,
    # 0.071, 0.050 at 400, 0.181, 0.139, 0.069, 0.027 at 800, 0.171, 0.134, 0.065, 0.023 at the
    # verification with 800: below threshold the current crosses the gate's barrier
    # thermally, and coarse macroparticles' grainy field rocks the barrier (at -46: 0.078,
    # 0.052, 0.007 at 200, 400, 800). So 400 per area: converged within 10 % to about -25;
    # pinch-off beyond -40 is shown as it is, its current not converged at the preview. At
    # 200 per area over the slider: 0.957, 0.598, 0.262 at +100, 40, 0, the channel's carriers
    # 13575, 6917, 2969 (doubled at 400: at +60 the verification's ~36000 pass the particle
    # cap and pause emission, where nothing solves). Goal 0.08-0.16: the gate about -13 to
    # -24.
    return level(
        "Gate",
        "The principle of the transistor. A thin strip of lightly doped silicon carries a "
        "current from the drain, held at 10, to the grounded source. Two gate electrodes beside "
        "it, on your supply, never touch it: a negative gate pushes the strip's carriers out "
        "and the current fades, until none flows; a positive one draws more in from the "
        "contacts. Set the gate so that between 0.08 and 0.16 flows into the drain, a third to "
        "two thirds of what flows with the gate at 0 (averaged over $t = 20$ to $60$). A metal screens any field within an atom, so no gate reaches "
        "inside copper; billions of these gates, each over a thin oxide on a silicon channel, "
        "switch in every chip.",
        grid=(20, 16), c=None, t_max=100.0,
        # 0: the source (grounded); 1: the drain (held at 10); 2, 3: the gates (your supply).
        electrodes=[contact(5, 8, 2), contact(15, 8, 2),
                    plate(10, 11, 6, thickness=1.0, tunable=True),
                    plate(10, 5, 6, thickness=1.0)],
        supplies=slider(60.0, -18.0, 0.0, 20.0),
        tube=tube(0, None, 0, 0, 0, 0, step=CHANNEL_STEP,
                  circuits=[tube_circuit(1, dc(10.0), resistance=0.0),
                            tube_circuit(2, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 0.08, 0.16, 20.0, 60.0, carriers="net")],
                  metals=[metal((6, 7), (14, 9), [0, 1], seed=3, **CHANNEL)]),
        reference=[supply(10, 11, -18.0)])


# The eddy-current levels' drive: the solenoid's field swings as B0 sin(wt), period ~50: long
# against the metal's settling (a few tenths) and its slow polarizations (a few units), and
# slow enough that a cut's two faces, a Debye length apart, carry little current across as
# a capacitor (at period 25 a core cut along its length took 61 % more heat than Ohm's law
# gives: the cut's 1/(wC) was comparable to a strip's resistance; at period 50, 14 %).
EDDY_OMEGA = 0.125


def cold_electrode(x, y):
    """A small floating electrode (zero charge) that emits nothing (saturation 0): the cathode
    a tube level with only isolated metals needs. Floating, it cannot draw charge from this
    open world's "infinity", which made a metal beside a grounded electrode run away
    (docs/OPEN_PROBLEMS.md, OP-1)."""
    return plate(x, y, 0.5, thickness=0.5, kind="charge", value=0.0)


def eddy_currents():
    # Designed with `generator tube-scan` (target/scratch/eddy_window.py): an isolated plate
    # 12 x 12 on the solenoid's axis, its field B0 sin(0.125 t). Local Ohm's law (PHYSICS.md
    # 2.11, D8): heat sigma <Bdot^2> int u dA = 0.47 x 728.75 x (0.125 B0)^2 / 2 = 2.68 B0^2;
    # the carriers' linear response holds up to Bdot ~ 0.1 (qEl/kT = 0.1 at the edges'
    # middles): B0 up to ~0.8. Measured (preview / verification): 0.744 / 0.753 at B0 = 0.55,
    # 1.153 / 1.181 at 0.68, 1.504 / 1.595 at 0.8 (5-7 % below local Ohm's law: the walls,
    # PHYSICS.md 2.11): about 2.55 B0^2. With the goal 1.0-1.4 the slider solved on 0.66-0.72
    # only (`generator window`, every 0.0225; the verdict's error straddles the goal's edges
    # beyond): widened to 0.9-1.5. Then the swing's full range (the owner, 2026-10-20): +-3;
    # heat over B0^2 2.31, 2.51, 2.47, 2.48, 2.39 at B0 = 0.5, 1, 1.5, 2, 3 (previews; the
    # carriers barely warm), so the goal 1.2-4.0 solves on B0 ~ 0.69-1.26.
    return level(
        "Eddy currents",
        "A metal plate alone in your solenoid, whose field now swings back and forth, "
        "$B = B_0 \\sin \\omega t$. A changing magnetic field induces an electric field "
        "that circulates about the solenoid's axis (Faraday's law), and in the metal it drives "
        "the carriers round in loops: eddy currents, which heat it, the faster the field "
        "changes the more (as $\\dot B^2$). Set the swing $B_0$ so that the plate takes "
        "between 1.2 and 4 of heat per unit time (averaged over $t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        solenoid_fields=slider(3.0, 0.95),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, solenoid=(23, 19),
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 1.2, 4.0, 25.0, 125.0)],
                  metals=[metal((6, 4), (18, 16), [], seed=21, **METAL)]),
        reference=[supply(23, 19, 0.95)])


# Eddy currents in Ohm's regime (the owner, 2026-10-20: Eddy currents' plate reads as a doped
# semiconductor, omega_c tau 0.38 at its answer; make one like copper's). The same plate with
# carriers that collide far more often, in a field that swings faster, so that every ratio
# deciding its eddy currents is small, as in copper: omega tau = 0.049 (the swing slow against
# the collisions), omega_c tau <= 0.18 over the whole slider (0.052 at the answer; copper's
# 0.0044 per tesla), omega/(4 pi sigma) = 0.10 (the edge charge keeps up with the swing), the
# plate 40 mean free paths across (the other's 12). The metal's own ratios, though, are a lightly
# doped semiconductor's: omega_p tau 0.70 and a screening length of 1.14 mean free paths match
# n-type silicon at 1e16 per cm^3 within x1.08 (game scales.rs, `resemblance`), where copper's are
# 400 and 0.0014 (real metals' omega_p tau 50-500: the step must follow the plasma oscillation,
# so the model affords about 1). Named for silicon (the owner's choice). The mean free path 0.3
# (tau = 0.12), cores of radius 0.0025 (a/2l = 0.004, the model's limit 0.005; 96 000 of them),
# density 2.7 (sigma = 0.32), the step 0.033 (omega_p dt = 0.19).
SILICON = dict(density=2.7, thermal_speed=2.0, mean_free_path=0.3, core_radius=0.0025,
               carriers_per_area=50.0)
SILICON_STEP = 0.033
# 13 whole cycles of the heat (it goes as Bdot^2, period pi/omega) in the goal's window.
SILICON_OMEGA = 0.13 * math.pi


def eddy_silicon():
    # Designed with `generator tube-scan` (target/scratch/copper_scratch.py): heat (preview)
    # 0.114, 0.471, 1.02, 1.11, 1.19, 1.67, 2.06, 4.42, 7.35, 16.4 at B0 = 0.1, 0.2, 0.3, 0.32,
    # 0.35, 0.38, 0.43, 0.6, 0.8, 1.2: 11.5 B0^2, local Ohm's law's 11.1 (sigma <Bdot^2> int u
    # dA; 40 mean free paths across, its walls matter less than the other plate's, 5-7 %
    # below), the preview scattering ~10 % at small swings. With the goal 1.2-4.0 B0 solved on
    # 0.40-0.55 (5.0 % of +-1.5): from 1.0. Verified 2.067 +- 0.012 at 0.43 (the other plate's
    # +- 0.28: the shorter collision time averages the thermal noise out). The swing's +-1.5
    # spans heat 0 to ~26 with omega_c tau <= 0.18.
    return level(
        "Eddy currents in silicon",
        "The last plate's carriers flew far between collisions, so that the field turned their "
        "paths by a good part of a radian on the way ($\\omega_c \\tau \\approx 0.4$; the Scales "
        "section). This plate is lightly doped silicon's: its carriers collide far more often, "
        "in a field that swings faster. The field barely turns them between collisions "
        "($\\omega_c \\tau \\approx 0.05$), so the eddy currents follow Ohm's law alone, as they "
        "do in copper at any laboratory field, and their heat grows as $B_0^2$. Set the swing "
        "$B_0$ so that the plate takes between 1 and 4 of heat per unit time (averaged over "
        "$t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        solenoid_fields=slider(1.5, 0.43),
        tube=tube(0, None, 0, 0, 0, 0, step=SILICON_STEP, saturation=0.0, solenoid=(23, 19),
                  solenoid_drive=sine(1.0, SILICON_OMEGA),
                  heat_goals=[heat_goal(0, 1.0, 4.0, 25.0, 125.0)],
                  metals=[metal((6, 4), (18, 16), [], seed=21, **SILICON)]),
        reference=[supply(23, 19, 0.43)])


def laminations():
    # Designed with `generator tube-scan` (target/scratch/lam20_*.py): a transformer's core
    # 20 x 12 on the solenoid's axis (30 carriers per area: 14 400 in the verification), the
    # field fixed (B0 = 0.68 at omega 0.125), three cuts of any length and angle (the owner:
    # a bigger plate, cuts not only across the whole metal). Verification (local Ohm's law):
    # uncut 2.833 (3.060); three full cuts along, evenly 0.516 (0.277), at y = 7, 11, 14
    # 0.597, at 6, 10, 14 0.751; two along and one across 0.407; two along 0.755 (0.475);
    # one along and one across 0.791 (0.765); one along 1.260 (0.991); corner to corner 1.435;
    # one across 1.633 (1.692); two half-length cuts along 1.745, two staggered 16 long 1.737.
    # The long cuts' faces carry current across as capacitors (PHYSICS.md 2.11): along them
    # the gain falls short of Ohm's law. The goal, at most 0.65, takes all three cuts across
    # the whole core, along it, about evenly (the best preview of three along read 0.609).
    return level(
        "Laminations",
        "A transformer's core in a field that swings back and forth: the eddy currents waste "
        "its energy as heat. They run in loops across the whole core; cut it, and each loop "
        "must close within a narrower piece, so the currents are weaker and the heat falls, "
        "the narrower the pieces, the more. Drag across the core to cut it, any length and "
        "angle; you have three cuts. Make the core take at most 0.65 of heat per unit time "
        "(averaged over $t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        max_cuts=3,
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, b_z=0.68,
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 0.0, 0.65, 25.0, 125.0)],
                  metals=[metal((2, 4), (22, 16), [], seed=122,
                                **dict(METAL, carriers_per_area=30.0))]),
        reference=[cut(2, 8, 22, 8), cut(2, 12, 22, 12), cut(12, 4, 12, 16)])


ARCS.append(("Inside the components: vacuum tubes", [
    ("Introduction", [
        ("vacuum_diode", vacuum_diode),
        ("sideways", sideways),
        ("bipolar_diode", bipolar_diode),
    ]),
    ("Intermediate", [
        ("through_the_cloud", through_the_cloud),
        ("electron_gun", electron_gun),
        ("ion_diode", ion_diode),
        ("rectifier", rectifier),
        ("triode", triode),
        ("klystron", klystron),
    ]),
    ("Master", [
        ("magnetron", magnetron),
    ]),
]))

ARCS.append(("Inside the components: metals", [
    ("Introduction", [
        ("ohms_law", ohms_law),
        ("long_and_thin", long_and_thin),
        ("potentiometer", potentiometer),
    ]),
    ("Intermediate", [
        ("hall_probe", hall_probe),
        ("positive_carriers", positive_carriers),
        ("few_or_many", few_or_many),
        ("magnetoresistance", magnetoresistance),
        ("eddy_currents", eddy_currents),
        ("eddy_silicon", eddy_silicon),
        ("laminations", laminations),
    ]),
    ("Master", [
        ("strength_and_sign", strength_and_sign),
    ]),
]))

ARCS.append(("Inside the components: transistors", [
    ("Introduction", [
        ("gate", gate),
    ]),
]))

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
        species = (d.axis, screen, json.dumps(s["particle"]), s["launch"]["kinetic_energy"],
                   id(d) if d.own else None)
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
    for (axis, (c0, c1), _, _, _), shots, _, _, lo, hi in placed:
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
        if slug in NEEDS_BODY:
            built = json.load(open(path, encoding="utf-8"))
            vacuum = dict(built, dielectrics=[], dielectric_spheres=[], ferrites=[])
            with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                             encoding="utf-8") as f:
                json.dump(vacuum, f)
                probe_path = f.name
            try:
                alone = run_generator("check", probe_path)
            finally:
                os.unlink(probe_path)
            works = all("Arrived" in l and "Verified" in l
                        for l in alone.splitlines() if "shot " in l)
            print(f"{key}: reference without the dielectrics and ferrites: "
                  f"{'WORKS (the body does not matter)' if works else 'fails, as it should'}")
            if works:
                sys.exit(f"{key}: the reference works without its dielectrics and ferrites")
        if not ok:
            print(report)
        for problem in slider_problems(json.load(open(path, encoding="utf-8"))):
            sys.exit(f"{key}: {problem}")

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
