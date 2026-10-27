"""Arcs: multi-stage instruments, beams, real instruments (idealised designs that must fail), finales."""

import json
import os
import sys
import tempfile
from designs.common import (
    Auto,
    M,
    ROOT,
    TEXT_ONLY,
    beam,
    box,
    charge,
    level,
    magnet,
    plate,
    rect_coil,
    run_generator,
    shot,
    slider,
    stray,
)
from designs.fields import RF, RF_PERIOD, WIEN_DETECTORS

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
    from designs.curriculum import LEVELS  # at run time: the curriculum imports this module

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
    supplies). Returns the solver's elements (none with `--text`: the built reference is
    kept)."""
    if TEXT_ONLY:
        return []
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
