# Curriculum: modules, practice, finales

The owner's review (2026-09-28): most levels are solved with one or two elements, which gets
boring outside the introductions. Hard levels should need several steps: the player builds
mental models of what geometric arrangements do, recognises that something similar worked
before, and modifies or combines earlier solutions. Inspiration: *Turing Complete*, where
each set of elements is introduced in small steps and ends in a grand finale that combines
everything learned; then a new element starts small again and interacts with the old ones.

The existing level ideas stay; they are regrouped, improved or joined.

## Principles

1. **Modules.** Each arc teaches named building blocks, one per short introduction level: the
   deflector, the attractor turn, the lens (a symmetric pair), the mirror, the energy
   analyser, the collimator; later the plate deflector, the aperture lens, the shield, the
   magnetic sector, the velocity filter, the RF kicker. The level text names the module, so
   that the player builds a vocabulary.
2. **Practice.** After the introductions, levels need two or three modules or a module with
   a twist (several shots, an acceptance condition, a real effect).
3. **Combination.** Later levels of an arc combine its modules with those of earlier arcs.
4. **Finale.** Each arc ends in a multi-stage level (gates) that needs several modules at
   once: 5–10 elements, each stage a sub-problem the player recognises.
5. **Realistic iterations** (the former chapter 15) start from an idealised design and add a
   real effect; they sit late in the arc of the instrument they improve.
6. **Measured.** `generator analyze` reports the fewest elements that solve a level (a
   sweep with the element counts capped); a finale must not be solvable with fewer than its
   stated minimum. Staged levels are checked stage by stage: each stage must be solvable
   given the previous ones.
7. **References by composition.** The solver cannot find 8-element designs from scratch;
   finale references are built stage by stage in `scripts/levels.py` (each stage solved as
   a sub-level whose detector is the next gate, with the earlier stages fixed), as a
   player would build them.
8. Later: par scores (fewest elements, least total charge, largest margins) and
   blueprints (SPEC Stage 8).

## Arcs (existing levels regrouped; new finales)

| Arc | Introductions (modules) | Practice and combination | Finale (new) |
|---|---|---|---|
| 1 Charges: steering and optics | First bend (deflector), Slingshot (attractor turn), Geiger–Marsden (scattering), Twin beams (several shots), Einzel lens (lens), Reflectron (mirror), Hemispherical analyser (analyser), Collimator (collimator), Two stages (gates) | The wall, Thomson CRT, Injection, Soft landing | **Sorting station**: two energies × three angles; focus through a gate, then sort the energies |
| 2 Metal and electrodes | High-voltage dome, Polarised sphere, Image charge, Deflection plates, Power supply, Build a deflector | Shielding, Tune the lens, Real Einzel lens, Beam pipe | **Microscope column**: condenser lens, deflector and shield from plates, in stages |
| 3 Relativity and magnetism | Fast lane, First coil, First magnet, Stern–Gerlach | Beta spectrometer, Dempster, Wien filter, Calutron, Build a Wien filter, CRT in the Earth's field | **Mass spectrometer from parts**: velocity filter, magnetic sector, collectors |
| 4 Time: noise and RF | Stray field, RF kick | Mains hum, Earth's field, RF separator, Tune the RF, Streak camera, Synchrotron light | **RF beam line**: kick, separate and streak bunches under a stray field |
| 5 Beams | Space charge, Stern–Gerlach beam | Collimated beam, Velocity selector, Relativistic beam, Beam preparation, Chromatic aberration, Real analyser, Calutron at full current, Soft landing at full current | **Isotope separator**: an interacting two-isotope beam through collimation, velocity selection and a sector |

The arc 1 finale is built first as the prototype of the method: level 52 (see LEVELS.md).
Its reference needs 2 charges in stage 1 and 4 in stage 2; the build checks that no search
finds a solution with 4 charges or fewer. The other arcs follow after the owner's review of
the prototype.

**Colours show destinations.** Shots aimed at the same detector share its colour and differ
only in lightness (`draw::shot_color`), so the player sees which ray must go where. Before,
each shot had its own colour and a shared detector was drawn in all of them at once (the
owner's review of the prototype).

**Measured starting point** (`analyze --fewest`, LEVELS.md): before the restructure, 40 of
51 levels were solved with a single element, 10 with two, and one (44) with three. The
Wien filter (26) and the velocity selector (43), meant as crossed-field levels, are solved
with one charge: they are the first to tighten.
