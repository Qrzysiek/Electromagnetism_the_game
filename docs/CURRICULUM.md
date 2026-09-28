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
8. **Margin.** A level's limits are never the fewest elements that solve it: they allow at
   least two more (and a choice of magnitudes), so that players can find less-than-optimal
   solutions, or just play around with the system (the owner, 2026-09-28).
9. Later: par scores (fewest elements, least total charge, largest margins) and
   blueprints (SPEC Stage 8).

## Arcs as built

The Jackson series (docs/JACKSON.md) continues the curriculum with arcs of the same
shape, each citing the sections and problems of the book it follows. Every arc has three tiers (`ARCS` in `scripts/levels.py`, written to
`levels/curriculum.json`; the game groups its level list by arc and shows each level's arc
and tier):

- **Introduction:** one new element or condition per level, solvable with one element (the
  limits still allow more).
- **Intermediate:** the arc's modules interact, or an earlier design meets a real effect:
  usually two or more elements.
- **Master:** a multi-stage finale where all bets are off, but every stage is a module the
  player has used before.

| Arc | Introduction | Intermediate | Master |
|---|---|---|---|
| 1 Charges: steering and optics | First bend, Slingshot, Geiger–Marsden, Twin beams, Two stages, Injection | Thomson's CRT, Around the wall (the wall joined with injection), Einzel lens, Reflectron, Collimator, Hemispherical analyser, Soft landing | **Sorting station**: lens through a gate, then sort two energies onto their own spots |
| 2 Metal and electrodes | High-voltage dome, Polarised sphere, Image charge, Deflection plates, Power supply, Build a deflector, Shielding | Tune the lens, Real Einzel lens, Beam pipe | **Microscope column**: tune the condenser through a crossover gate, then deflector plate and projector charges onto a spot beside an ion pump, arriving along the axis |
| 3 Relativity and magnetism | Fast lane, First coil, First magnet, Stern–Gerlach | Beta spectrometer, Dempster, Wien filter, Calutron, Build a Wien filter | **Mass spectrometer from parts**: one electrostatic lens for three masses (only T/q matters), then a magnetic sector sorts them by momentum |
| 4 Time: noise, RF and radiation | Stray field, RF kick, Synchrotron light | Mains hum, Earth's field, CRT in the Earth's field, RF separator, Tune the RF, Streak camera | **RF beam line**: steer two bunches through a gate with a stray field on and off, then separate them with RF |
| 5 Beams | Space charge, Stern–Gerlach beam, Relativistic beam | Collimated beam, Velocity selector, Beam preparation, Chromatic aberration, Real analyser, Calutron at full current, Soft landing at full current | **Isotope separator**: collimate an interacting two-isotope beam through a gate, then separate the isotopes with magnets |
| 6 Jackson: charges in fields and collisions (Ch. 12–13) | §12.3 E×B drift, Pr. 12.9 Van Allen equator, Collision course (the free charge), Pr. 13.1 knock-on | §12.4 gradient drift, Pr. 12.5 E×B runaway, §12.1 Störmer's forbidden region | **Magnetosphere**: steer the solar wind's E×B drift through two gates, then the gradient drift splits proton and electron around a dipole Earth |
| 7 Jackson: conductors (Ch. 2–3) | §2.2 its own image, Pr. 2.6 two spheres | §3.13 field through a hole, Pr. 2.4 golden-ratio capture | **Sphere slalom**: weave between three spheres carrying the particle's charge |
| 8 Jackson: radiation damping (Ch. 16) | Pr. 16.2 the classical atom | Pr. 16.3 orbits circularize | **Three orbits**: three radiating electrons, each shaped to collapse in time |
| 9 Jackson: bound charges (§16.7–16.8, Pr. 13.2) | §16.7 a bound charge, §16.8 resonance | Pr. 13.2 a kick for a bound charge | **Spectroscopy**: two atoms, two natural frequencies, each driven at its own |

Changes to existing levels: The wall and Injection were joined into Around the wall (the
beam must enter the detector along the axis after going around the wall: two charges).
The Wien filters now require the selected speed to leave straight along the axis (±4°):
before, one charge's non-uniform field could sort the three speeds; now the search needs
a crossed-field pair. Every level's limits leave a margin of at least two elements above
its fewest solution (rule 8). The real-instrument levels moved into the arcs of the
effects they add (CRT in the Earth's field to arc 4, Beam pipe to arc 2, the beam ones to
arc 5), after the levels they build on.

The arc 1 finale is built first as the prototype of the method (see LEVELS.md).
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
