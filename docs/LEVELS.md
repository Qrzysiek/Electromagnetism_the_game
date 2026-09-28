# Levels: curriculum and difficulty

`scripts/levels.py` is the single source of truth for the shipped levels: their order,
designs, reference solutions and automatic detector placement. Rebuild with
`python scripts/levels.py` (or `--only NN` for one level). Every level must have a
verified reference solution for every shot and negligible radiation
(`crates/level/tests/levels.rs`).

## Curriculum rules

- The levels form five arcs (docs/CURRICULUM.md; `ARCS` in `scripts/levels.py`, written
  to `levels/curriculum.json`). Each arc has introductions (one new element or condition
  each), intermediate levels (the arc's modules interact, or an earlier design meets a
  real effect) and a master level (a multi-stage finale built from the arc's modules).
- Every new element or concept is introduced by an easy level first.
- Limits leave a margin of at least two elements above a level's fewest solution, so that
  players can find less-than-optimal solutions or just play around.
- Notes below name levels rather than numbers: the numbers change when levels are added.
- A good puzzle has a large configuration space, a small solution set, and a "distance to
  solution" that changes piecewise continuously. A player can then learn from attempts,
  but guessing or brute force is slow.

| # | Level | Arc, tier | Idea | Player elements | Reference |
|---|---|---|---|---|---|
| 1 | First bend | 1 Introduction | one charge bends a beam | ≤ 3 charges | 1 |
| 2 | Slingshot | 1 Introduction | bending around a level charge | ≤ 3 charges | 1 |
| 3 | Geiger–Marsden | 1 Introduction | Coulomb scattering, sign and strength | ≤ 3 charges | 1 |
| 4 | Twin beams | 1 Introduction | several shots, one setup | ≤ 3 charges | 1 |
| 5 | Two stages | 1 Introduction | a gate to pass before the detector counts | ≤ 3 charges | 1 |
| 6 | Injection | 1 Introduction | the detector also requires a direction (±8°) | ≤ 3 charges | 1 |
| 7 | Thomson's cathode-ray tube | 1 Intermediate | two energies onto one-cell spots; deflection ∝ 1/T | ≤ 3 charges | 1 |
| 8 | Around the wall | 1 Intermediate | the wall joined with injection: around the obstacle, then enter along the axis (±8°) | ≤ 4 charges | 2 |
| 9 | Einzel lens | 1 Intermediate | focusing of an angular spread | ≤ 4 charges | 2 |
| 10 | Reflectron | 1 Intermediate | reflection; energy-dependent turning point | ≤ 4 charges | 2 |
| 11 | Collimator | 1 Intermediate | three rays must arrive parallel (±3°) | ≤ 4 charges | 2 |
| 12 | Hemispherical analyzer (XPS) | 1 Intermediate | energy dispersion on a circular orbit | ≤ 3 charges | 1 |
| 13 | Soft landing | 1 Intermediate | arrival energy window: brake without reflecting | ≤ 3 charges | 1 |
| 14 | Sorting station | 1 Master | two energies × three angles: a lens through a gate, then sort the energies onto their own spots | ≤ 8 charges | 6 |
| 15 | High-voltage dome | 2 Introduction | a sphere held at a potential acts like a charge at its centre, but it is a conductor | ≤ 3 charges | 1 |
| 16 | Polarised sphere | 2 Introduction | an isolated neutral sphere becomes a dipole near a charge | ≤ 3 charges | 1 |
| 17 | Image charge | 2 Introduction | a grounded sphere answers every charge with an opposite image | ≤ 3 charges | 1 |
| 18 | Deflection plates | 2 Introduction | real plates with fringe fields (BEM) | ≤ 3 charges | 1 |
| 19 | Power supply | 2 Introduction | power supplies: set an electrode's potential instead of placing charges | 1 supply | 1 |
| 20 | Build a deflector | 2 Introduction | placing your own plates (electrodes) | ≤ 3 plates | 1 |
| 21 | Shielding | 2 Introduction | a grounded plate screens a charge's field | ≤ 3 plates | 1 |
| 22 | Tune the lens | 2 Intermediate | focus a real Einzel lens with its voltage | 2 supplies | 2 |
| 23 | Real Einzel lens | 2 Intermediate | three apertures at high voltage focus three rays | ≤ 4 charges | 2 |
| 24 | Beam pipe | 2 Intermediate | the injection (6) built in; a grounded pipe wall between the beam and the steering charge screens it | ≤ 3 charges | 1 |
| 25 | Microscope column | 2 Master | tune the condenser through a crossover gate; deflector plate and charges onto a spot beside an ion pump, along the axis (±20°) | ≤ 5 charges, 3 supplies | 4 |
| 26 | Fast lane | 3 Introduction | γ changes the bending | ≤ 3 charges | 1 |
| 27 | First coil | 3 Introduction | magnetic fields from level coils | ≤ 3 charges | 1 |
| 28 | First magnet | 3 Introduction | placing your own magnets | ≤ 3 magnets | 1 |
| 29 | Stern–Gerlach | 3 Introduction | neutral atoms with spin up or down, force m grad B_z | ≤ 3 magnets | 1 |
| 30 | Beta-ray spectrometer | 3 Intermediate | relativistic circular orbits: \|qQ\| = γmv²R, two electron energies | ≤ 3 charges | 1 |
| 31 | Dempster's mass spectrometer | 3 Intermediate | 180° focusing and mass separation | ≤ 3 charges | 1 |
| 32 | Wien filter | 3 Intermediate | crossed E (your charges) and B (a coil); the selected speed leaves straight (±3°) | ≤ 4 charges | 2 |
| 33 | Calutron | 3 Intermediate | isotope separation with magnets only | ≤ 4 magnets | 2 |
| 34 | Build a Wien filter | 3 Intermediate | the whole velocity selector from charges and magnets; straight exit (±3°) | ≤ 4 charges, ≤ 2 magnets | 3 |
| 35 | Mass spectrometer from parts | 3 Master | one electrostatic lens for three masses (only T/q matters), then a magnetic sector sorts them by momentum | ≤ 5 charges, ≤ 4 magnets | 5 |
| 36 | Stray field | 4 Introduction | a stray field switched on and off; aim between | ≤ 3 charges | 1 |
| 37 | RF kick | 4 Introduction | an antenna's kick depends on the passing time | ≤ 3 antennas | 1 |
| 38 | Synchrotron light | 4 Introduction | radiative damping spirals the particle to the axis; only possible with radiation (canonical angular momentum) | ≤ 1 magnet | 1 |
| 39 | Mains hum | 4 Intermediate | uniform AC field at 4 phases acts like a random launch angle; imaging | ≤ 4 charges | 2 |
| 40 | Earth's field | 4 Intermediate | stray B_z of either sign on an electron beam | ≤ 4 charges | 2 |
| 41 | CRT in the Earth's field | 4 Intermediate | the CRT (7) built in, installed facing north in the Earth's field (stray B_z); adjusted where it stands | ≤ 3 charges | 1 |
| 42 | RF separator | 4 Intermediate | identical bunches half a period apart to different detectors | ≤ 1 charge, ≤ 2 antennas | 1 |
| 43 | Tune the RF | 4 Intermediate | choose the antenna frequency: bunches 4 apart meet it at phase difference ωΔt | ≤ 1 charge, ≤ 2 antennas | 1 |
| 44 | Streak camera | 4 Intermediate | three bunches a third of a period apart to three spots | ≤ 2 charges, ≤ 2 antennas | 1 |
| 45 | RF beam line | 4 Master | steer two bunches through a gate (±6°) with a stray field on and off, then separate them with RF | ≤ 5 charges, ≤ 3 antennas | 4 |
| 46 | Space charge | 5 Introduction | 16 particles repel each other; ≥ 90 % must arrive, verified | ≤ 3 charges | 1 |
| 47 | Stern–Gerlach beam | 5 Introduction | both spin states as spread beams, each to its own detector | ≤ 3 magnets | 1 |
| 48 | Relativistic beam | 5 Introduction | a beam at 0.8c: magnetic attraction weakens space charge to 1/γ²; quasi-static interaction with radiation reaction | ≤ 3 charges | 1 |
| 49 | Collimated beam | 5 Intermediate | the collimator for a spread, interacting beam (±3°) | ≤ 4 charges | 2 |
| 50 | Velocity selector | 5 Intermediate | the Wien filter for a beam: three speeds, 8 ions each; the middle one leaves straight (±5°) | ≤ 4 charges | 2 |
| 51 | Beam preparation | 5 Intermediate | two stages for a beam: collimate it through a gate (±4°), then steer it into the target | ≤ 5 charges | 3 |
| 52 | Chromatic aberration | 5 Intermediate | the Einzel lens (9) built in; each ray becomes a beam with an 8 % energy spread | ≤ 3 charges | 1 |
| 53 | Real analyser | 5 Intermediate | the hemispherical analyser (12) built in; the source emits into a cone (σ = 4°) | ≤ 3 charges | 1 |
| 54 | Calutron at full current | 5 Intermediate | the calutron (33) built in; the isotope beams repel each other (5 ions each, radiation reaction) | ≤ 3 charges | 1 |
| 55 | Soft landing, full current | 5 Intermediate | the soft landing (13) built in; space charge grows as the ions are braked | ≤ 3 charges | 1 |
| 56 | Isotope separator | 5 Master | collimate an interacting two-isotope beam through a gate (±5°), then separate the isotopes with magnets | ≤ 5 charges, ≤ 4 magnets | 5 |

## Automatic detectors

A shot's detector may be given as `Auto(strip, axis, size)`. The script then runs the
reference solution (or, if the level has none, first asks the solver for a setup that
brings every shot into its probe strip). It places a `size`-cell detector around the
landing point, clamped to the arena. Shots of the same species (same particle and energy,
differing only in direction, as in Dempster's 180° focusing) share one detector.
Overlapping detectors of different species on the same screen are split at the node
halfway between their landing points. If the points are too close to split, the build
fails.

## Difficulty measurements

`cargo run --release -p generator -- analyze [--fewest] levels/[0-9]*.json` (by default 2000 random samples,
32 search runs with a budget of 400 evaluations each; the table below used 1000 and 16).
Columns:

- **log10 configs:** number of distinct player placements within the limits.
- **random solve rate:** fraction of uniformly random placements that solve every shot,
  verified. When none were found, the upper estimate 3/n (rule of three) is shown.
- **search:** a local search that only sees the distance objective, like a player watching
  the preview (greedy moves with occasional worse moves accepted). Shows its success rate
  and the mean evaluations of successful runs.
- **expected effort:** expected number of tried placements for the search (restarting after
  each failed run) versus for random guessing.
- **smoothness:** Spearman rank correlation of the distance objective between a placement
  and a one-move neighbour. High means that attempts carry information.

Measured 2026-09-28 after the curriculum restructure, with 1000 random samples and 16
search runs of 400 evaluations per level:

| level | log10 configs | random solve rate | search success | mean evals | expected effort: search / guessing | smoothness |
|---|---|---|---|---|---|---|
| 01_first_bend | 10.0 | 7.1e-2 | 16/16 | 10 | 10 / 14 | 0.77 |
| 02_slingshot | 10.0 | 2.1e-2 | 16/16 | 57 | 57 / 48 | 0.71 |
| 03_geiger_marsden | 9.5 | 1.7e-2 | 15/16 | 42 | 69 / 59 | 0.41 |
| 04_twin_beams | 10.0 | 6.0e-3 | 16/16 | 56 | 56 / 167 | 0.87 |
| 05_two_stages | 10.0 | 3.5e-2 | 16/16 | 58 | 58 / 29 | 0.79 |
| 06_injection | 10.4 | 3.2e-2 | 16/16 | 41 | 41 / 31 | 0.80 |
| 07_thomson_crt | 8.8 | 1.3e-2 | 15/16 | 82 | 109 / 77 | 0.75 |
| 08_around_the_wall | 13.5 | < 3.0e-3 | 3/16 | 61 | 1794 / 333 | 0.78 |
| 09_einzel_lens | 10.3 | < 3.0e-3 | 15/16 | 125 | 151 / 333 | 0.84 |
| 10_reflectron | 10.1 | 2.0e-3 | 16/16 | 108 | 108 / 500 | 0.84 |
| 11_collimator | 11.1 | < 3.0e-3 | 16/16 | 84 | 84 / 333 | 0.89 |
| 12_hemispherical_analyzer | 9.5 | 2.0e-3 | 9/16 | 97 | 408 / 500 | 0.72 |
| 13_soft_landing | 9.7 | 1.0e-3 | 2/16 | 174 | 2974 / 1000 | 0.81 |
| 14_sorting_station | 26.1 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.93 |
| 15_high_voltage_dome | 10.3 | 1.2e-1 | 16/16 | 12 | 12 / 8 | 0.76 |
| 16_polarized_sphere | 9.5 | 1.4e-2 | 16/16 | 30 | 30 / 71 | 0.65 |
| 17_image_charge | 10.3 | 2.1e-2 | 15/16 | 58 | 85 / 48 | 0.83 |
| 18_deflection_plates | 8.9 | 1.6e-1 | 16/16 | 5 | 5 / 6 | 0.75 |
| 19_power_supply | 0.7 | 2.7e-1 | 16/16 | 3 | 3 / 4 | -0.43 |
| 20_build_a_deflector | 8.7 | 1.2e-1 | 16/16 | 12 | 12 / 8 | 0.65 |
| 21_shielding | 7.1 | 3.7e-2 | 16/16 | 16 | 16 / 27 | 0.61 |
| 22_tune_the_lens | 1.4 | 4.2e-2 | 11/16 | 23 | 205 / 24 | 0.30 |
| 23_real_einzel_lens | 10.8 | 1.6e-2 | 16/16 | 42 | 42 / 62 | 0.84 |
| 24_beam_pipe | 9.6 | 3.2e-2 | 16/16 | 20 | 20 / 31 | 0.81 |
| 25_microscope_column | 16.2 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.80 |
| 26_fast_lane | 10.7 | 1.1e-2 | 16/16 | 79 | 79 / 91 | 0.67 |
| 27_first_coil | 10.0 | 4.3e-1 | 16/16 | 4 | 4 / 2 | 0.44 |
| 28_first_magnet | 9.1 | 7.7e-2 | 16/16 | 18 | 18 / 13 | 0.68 |
| 29_stern_gerlach | 8.8 | 6.8e-2 | 16/16 | 15 | 15 / 15 | 0.80 |
| 30_beta_spectrometer | 9.5 | 2.0e-3 | 8/16 | 117 | 517 / 500 | 0.62 |
| 31_dempster | 8.4 | 8.0e-2 | 16/16 | 12 | 12 / 12 | 0.56 |
| 32_wien_filter | 11.6 | < 3.0e-3 | 5/16 | 188 | 1068 / 333 | 0.66 |
| 33_calutron | 12.3 | < 3.0e-3 | 5/16 | 101 | 981 / 333 | 0.82 |
| 34_build_wien_filter | 18.2 | < 3.0e-3 | 2/16 | 332 | 3132 / 333 | 0.82 |
| 35_mass_spectrometer | 30.1 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.89 |
| 36_stray_field | 10.7 | 5.1e-2 | 16/16 | 32 | 32 / 20 | 0.78 |
| 37_rf_kick | 12.2 | 6.1e-2 | 16/16 | 19 | 19 / 16 | 0.65 |
| 38_synchrotron_light | 1.0 | 5.3e-1 | 16/16 | 2 | 2 / 2 | -0.08 |
| 39_mains_hum | 12.3 | < 3.0e-3 | 9/16 | 131 | 442 / 333 | 0.88 |
| 40_earths_field | 12.4 | 1.0e-3 | 5/16 | 218 | 1098 / 1000 | 0.81 |
| 41_crt_earth_field | 9.7 | 5.0e-3 | 15/16 | 68 | 95 / 200 | 0.80 |
| 42_rf_separator | 11.0 | 4.0e-3 | 15/16 | 88 | 115 / 250 | 0.76 |
| 43_tune_the_rf | 12.1 | 2.0e-3 | 15/16 | 63 | 89 / 500 | 0.72 |
| 44_streak_camera | 14.2 | < 3.0e-3 | 7/16 | 194 | 708 / 333 | 0.85 |
| 45_rf_beam_line | 28.7 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.90 |
| 46_space_charge | 10.0 | 2.1e-2 | 15/16 | 59 | 86 / 48 | 0.80 |
| 47_stern_gerlach_beam | 8.8 | 3.3e-2 | 14/16 | 17 | 74 / 30 | 0.82 |
| 48_relativistic_beam | 10.4 | 2.8e-2 | 15/16 | 23 | 50 / 36 | 0.75 |
| 49_collimated_beam | 11.1 | 1.0e-3 | 16/16 | 96 | 96 / 1000 | 0.88 |
| 50_velocity_selector | 11.6 | < 3.0e-3 | 1/16 | 306 | 6306 / 333 | 0.68 |
| 51_beam_preparation | 14.2 | < 3.0e-3 | 7/16 | 274 | 789 / 333 | 0.91 |
| 52_chromatic_aberration | 9.8 | 6.0e-3 | 14/16 | 43 | 101 / 167 | 0.87 |
| 53_real_analyzer | 8.6 | 1.3e-2 | 12/16 | 17 | 151 / 77 | 0.74 |
| 54_calutron_space_charge | 9.7 | 1.0e-3 | 9/16 | 99 | 410 / 1000 | 0.81 |
| 55_soft_landing_current | 9.1 | < 3.0e-3 | 12/16 | 200 | 334 / 333 | 0.74 |
| 56_isotope_separator | 28.6 | < 3.0e-3 | 1/16 | 187 | 6187 / 333 | 0.91 |

Reading the table by tier: the introductions are solved by the search in 15–16 of 16 runs
within 3–80 evaluations. The intermediate levels need more (Around the wall 3/16, Soft
landing 2/16, Build a Wien filter 2/16, Velocity selector 1/16). The master levels are not
found by a search over the whole level (0/16, the Isotope separator 1/16), since they
are meant to be built stage by stage, as their references were; random guessing of the
whole level is hopeless (log10 configs 16–30). Their stages are each within reach of the
modules the arc taught (see the finale notes below).

Notes:

- The introductions are meant to be easy. Synchrotron light is a concept level (8
  choices).
- Under disturbances a placement counts as a solution only if every flight arrives,
  verified. The distance objective sums over all flights.
- Dempster has in effect one physical parameter (the accelerating voltage, which sets all
  radii). It is kept as the concept level for 180° focusing.
- The Calutron is among the hardest intermediate levels: no random placement out of 2000
  solved it, and the search succeeded in 4 of 32 runs. Its difficulty comes from the
  magnets' strongly non-uniform dipole fields.
- Power supply and Tune the lens are knob levels with tiny configuration spaces (5 and 25
  settings). The analysis caches the result of every distinct placement, so such levels
  cost only their distinct placements; before the cache (and before it knew about
  supplies), analysing them did not terminate.
- Real instruments (Beam pipe, CRT in the Earth's field, Chromatic aberration, Real
  analyser, Calutron at full current, Soft landing at full current) build an earlier
  level's reference setup in as fixed elements and add one real effect the ideal design
  ignored. They sit in the arc of the effect they add, after the level they build on.
  `scripts/levels.py` checks at build time that the built-in design alone fails, and the
  test `realistic_levels_break_the_idealised_design` keeps checking it. Measured with the
  built-in design alone: Chromatic aberration one ray 7/8 arriving; Real analyser the
  highest energy 7/8; CRT in the Earth's field both beams lost; Beam pipe rejected (wrong
  direction); Calutron at full current the light isotope 3/5; Soft landing at full current
  0/8. Tuning for playability (no search solved the first versions): the CRT first had to
  work facing north, south and shielded with one-cell spots, and at most 1 in 32 searches
  solved it even at 40 % of the field; a tube is adjusted where it stands, so the level now
  has one orientation. Soft landing at full current had 12 ions at the full charge scaling
  (0/32); now 8 at half.
- Beam levels are harder versions of single-particle levels: a placement solves them only
  if the required share of the beam arrives, verified. Their objective for the search sums
  the detector distances (plus gate shortfalls) of the closest missing particles.
  Relativistic beam includes radiation reaction because its particles would otherwise
  neglect 5.3e-6 of their energy (above the 1e-10 allowed); the test
  `radiation_levels_need_radiation` accepts that reason for beams. Beam preparation
  combines the collimated beam with a second stage; the search found its reference in 44 s
  and solves it in 9 of 32 runs.
- The analysis logs its progress to stderr (`[analyze]` lines: phases with timing, distinct
  placements flown, the slowest evaluation, and every evaluation slower than 1 s with its
  placement). The log found why the velocity selector's analysis ran for hours (placements
  trapping the ions for 1e6 steps, and a beam bug); searches now stay within the step
  budget of the cost meters.
- Finales (master levels) have references built stage by stage (`solve_stage` in
  `scripts/levels.py`): each stage is solved as a sub-level whose detector is the next
  gate, with the earlier stages fixed, as a player would build it.
  - Sorting station: the lens (stage 1, at most 3 charges left of the gate) needs 2
    charges; sorting and refocusing (stage 2, at most 4 charges) needs 4. The build checks
    that the search finds no solution with 4 charges anywhere (`FINALES`), and the
    whole-level search with up to 8 charges finds none either (0/4 runs). Versions that
    asked for parallel arrival (within ±6°, ±10° or ±15°) had no staged solution;
    five-cell spots were solved with 3 charges in stage 2 (too easy).
  - Microscope column: the condenser (two supplies, both at 400 kV) and then the deflector
    supply with one charge: 4 elements. The spot beside the ion pump asks for arrival
    along the axis (±20°): stage 2 is solved in 4 of 48 search runs; at ±10° in none of
    16. The first version (8 full-size electrodes) took 8.5 s to set up, over the sandbox
    budget, and was solved with 3 elements; smaller, lower electrodes take 0.75 s.
  - Mass spectrometer from parts: one charge makes the lens for all three masses (only T/q
    matters for electrostatic optics); the magnetic sector takes 4 elements. The
    collectors are placed automatically at the reference's landing points.
  - RF beam line: stage 1 (through a two-cell gate within ±6°, with the stray field on and
    off) takes one charge; the RF separation onto two-cell spots takes 3 elements. Looser
    versions (a four-cell gate without direction, ±8°, half the stray field, four-cell
    spots) were solved with 3 elements in all.
  - Isotope separator: collimation (stage 1) takes one charge; the magnetic separation 4
    elements. A first version with the collectors on opposite sides of the axis had no
    solution: a magnetic field bends both isotopes the same way, the lighter more.
- `analyze --fewest` adds the fewest elements that solve a level: the reference solution's
  count, unless the search finds a verified solution with fewer (exhaustive for one
  element, annealing with 64 restarts of 600 iterations for 2 up to the reference count
  minus one). A failed search is not a proof, so the number is an upper bound.
  Measured 2026-09-28, before the restructure (old numbering): 2 elements for Reflectron,
  Einzel lens, Collimator, Tune the lens, Real Einzel lens, Calutron, Build a Wien filter,
  Mains hum, Earth's field and Collimated beam; 3 for Beam preparation; 6 for Sorting
  station; 1 for every other level. Three references used more elements than needed: the
  Wien filter (reference 4) and the velocity selector (reference 4) were solved with one
  charge, and Build a Wien filter (reference 6) with two elements. The Wien filters now
  require a straight exit (±3°; at ±4° one placement of a single charge still passed) and
  need 2 and 3 elements. The velocity selector became the Wien filter for a beam (three
  speeds, straight exit within ±5°): no one-charge solution, two-charge ones in 2 of 32
  search runs; the Wien filter's own design fails there (space charge). The Reference column of the table above gives the current
  counts.
- The search is a rough stand-in for a player. A player who understands the physics (for
  example the Wien condition v = E/B) should do much better than the search on the later
  levels, and a player who doesn't should do much worse.
