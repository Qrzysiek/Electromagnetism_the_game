"""Shared by every level design: units and constants, the particles, the board's elements, automatic detectors (`Auto`), helpers, and running the generator."""

import math
import os
import subprocess
import sys

# The repository (this file is scripts/designs/common.py).
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# `--text`: rewrite the built levels' texts only, running nothing (`text_only`).
TEXT_ONLY = "--text" in sys.argv

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



def metal_sphere(x, y, r, kind="grounded", value=None):
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
        limits["free_charges"] = list(free_charges)
        limits["free_speeds"] = list(free_speeds)
        if free_mass != 1.0:
            limits["free_mass"] = free_mass
        if free_radius != 0.3:
            limits["free_radius"] = free_radius
    if max_cuts:
        limits["max_cuts"] = max_cuts
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
# Machinery.

def run_generator(*args):
    out = subprocess.run([GENERATOR, *args], capture_output=True, text=True, cwd=ROOT,
                         encoding="utf-8")
    if out.returncode != 0:
        sys.exit(f"generator {' '.join(args)} failed:\n{out.stdout}\n{out.stderr}")
    return out.stdout
