# Levels: curriculum and difficulty

`scripts/levels.py` is the single source of truth for the shipped levels: their order,
designs, reference solutions and automatic detector placement. Rebuild with
`python scripts/levels.py` (or `--only NN` for one level). Every level must have a
verified reference solution for every shot and negligible radiation
(`crates/level/tests/levels.rs`).

## Curriculum rules

- Levels are ordered mostly by how hard the phenomenon is to understand and use.
- Every new element or concept is introduced by an easy level first.
- Within a chapter, difficulty and the number of elements needed rise on average. This
  is a trend, not a strict rule.
- A good puzzle has a large configuration space, a small solution set, and a "distance to
  solution" that changes piecewise continuously. A player can then learn from attempts,
  but guessing or brute force is slow.

| # | Level | Chapter / new idea | Player elements |
|---|---|---|---|
| 1 | First bend | 1 Charges: one charge bends a beam | ≤ 1 charge |
| 2 | Geiger–Marsden | Coulomb scattering, sign and strength | ≤ 1 charge |
| 3 | Thomson's CRT | two energies; deflection ∝ 1/T | ≤ 2 charges |
| 4 | Slingshot | bending around a level charge | ≤ 1 charge |
| 5 | The wall | steering around an obstacle | ≤ 3 charges |
| 6 | Twin beams | 2 Several shots, one setup | ≤ 2 charges |
| 7 | Reflectron | reflection; energy-dependent turning point | ≤ 3 charges |
| 8 | Einzel lens | focusing of an angular spread | ≤ 4 charges |
| 9 | Hemispherical analyzer | energy dispersion on a circular orbit | ≤ 2 charges |
| 10 | Injection | 3 Delivering beams: the detector also requires a direction (±8°) | ≤ 1 charge |
| 11 | Collimator | three rays must arrive parallel (±3°) | ≤ 4 charges |
| 12 | Soft landing | arrival energy window: brake without reflecting | ≤ 2 charges |
| 13 | High-voltage dome | 4 Metals: a sphere held at a potential acts like a charge at its centre, but it is a conductor | ≤ 1 charge |
| 14 | Polarised sphere | an isolated neutral sphere becomes a dipole near a charge | ≤ 2 charges |
| 15 | Image charge | a grounded sphere answers every charge with an opposite image | ≤ 2 charges |
| 16 | Deflection plates | 5 Electrodes: real plates with fringe fields (BEM) | ≤ 1 charge |
| 17 | Power supply | power supplies: set an electrode's potential instead of placing charges | 1 supply (4 potentials) |
| 18 | Tune the lens | focus a real Einzel lens with its voltage | 2 supplies (4 potentials) |
| 19 | Real Einzel lens | three apertures at high voltage focus three rays | ≤ 2 charges |
| 20 | Build a deflector | placing your own plates (electrodes) | ≤ 1 plate (4 potentials) |
| 21 | Shielding | a grounded plate screens a charge's field | ≤ 1 grounded plate |
| 22 | Fast lane | 6 Relativity: γ changes the bending | ≤ 2 charges |
| 23 | Beta spectrometer | relativistic circular orbits of two energies | ≤ 2 charges |
| 24 | First coil | 7 Magnetic fields from level coils | ≤ 1 charge |
| 25 | Dempster | 180° focusing and mass separation | ≤ 2 charges |
| 26 | Wien filter | crossed E and B select one speed | ≤ 4 charges |
| 27 | First magnet | 8 Placing your own magnets | ≤ 1 magnet |
| 28 | Calutron | isotope separation with magnets only | ≤ 2 magnets |
| 29 | Build a Wien filter | crossed fields from charges and magnets | ≤ 4 charges, ≤ 2 magnets |
| 30 | Stern–Gerlach | 9 Magnetic moments: neutral atoms with spin up or down, force m grad B_z | ≤ 1 magnet |
| 31 | Stray field | 10 Noise: a stray field switched on and off; aim between | ≤ 1 charge |
| 32 | Mains hum | uniform AC field at 4 phases acts like a random launch angle; imaging | ≤ 4 charges |
| 33 | Earth's field | stray B_z of either sign on an electron beam | ≤ 3 charges |
| 34 | RF kick | 11 Radio frequency: an antenna's kick depends on the passing time | ≤ 1 antenna |
| 35 | RF separator | identical bunches half a period apart to different detectors | ≤ 2 antennas, ≤ 1 charge |
| 36 | Tune the RF | choose the antenna frequency: the phase difference ω Δt decides | ≤ 1 antenna (5 frequencies), ≤ 1 charge |
| 37 | Streak camera | three bunches a third of a period apart to three spots | ≤ 2 antennas, ≤ 2 charges |
| 38 | Synchrotron light | 12 Radiation: radiative damping spirals the particle to the axis; only possible with radiation (canonical angular momentum) | 1 magnet on the axis |
| 39 | Two stages | 13 Multi-stage instruments: a gate to pass before the detector counts | ≤ 1 charge |
| 40 | Space charge | 14 Beams: 16 particles repel each other; ≥ 90 % must arrive, verified | ≤ 1 charge |
| 41 | Stern–Gerlach beam | both spin states as spread beams, each to its own detector | ≤ 1 magnet |
| 42 | Collimated beam | the collimator for a spread, interacting beam (±3°) | ≤ 4 charges |
| 43 | Velocity selector | two speeds in one interacting beam sorted by crossed fields | ≤ 4 charges |
| 44 | Beam preparation | two stages for a beam: collimate it through a gate (±4°), then steer it into the target | ≤ 4 charges |
| 45 | Relativistic beam | a beam at 0.8c: magnetic attraction weakens space charge to 1/γ²; quasi-static interaction with radiation reaction | ≤ 2 charges |
| 46 | Chromatic aberration | 15 Real instruments: level 8's Einzel lens built in; each ray becomes a beam with an 8 % energy spread | ≤ 2 charges |
| 47 | Real analyser | level 9's hemispherical analyser built in; the source emits into a cone (σ = 4°) | ≤ 2 charges |
| 48 | CRT in the Earth's field | level 3's CRT built in, installed facing north in the Earth's field (stray B_z); adjusted where it stands | ≤ 2 charges |
| 49 | Beam pipe | level 10's injection built in; a grounded pipe wall between the beam and the steering charge screens it | ≤ 2 charges |
| 50 | Calutron at full current | level 28's calutron built in; the isotope beams repel each other (5 ions each, a quarter of the usual charge scaling, radiation reaction) | ≤ 2 charges |
| 51 | Soft landing, full current | level 12's soft landing built in; an interacting beam: space charge grows as the ions are braked | ≤ 2 charges |
| 52 | Sorting station | Finale prototype (docs/CURRICULUM.md): six rays, two energies × three angles; a lens brings them through a two-cell gate, then a deflector and a lens refocus each energy onto its own two-cell spot | ≤ 8 charges (reference 6: 2 + 4) |

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

`cargo run --release -p generator -- analyze [--fewest] levels/[0-9]*.json` (2000 random samples,
32 search runs with a budget of 400 evaluations each; the electrode levels, whose flights
cost more, with 600 samples and 16 runs). Columns:

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

| level | log10 configs | random solve rate | search success | mean evals | expected effort: search / guessing | smoothness |
|---|---|---|---|---|---|---|
| 01_first_bend | 3.6 | 1.1e-1 | 32/32 | 8 | 8 / 9 | 0.61 |
| 02_geiger_marsden | 3.4 | 2.8e-2 | 32/32 | 33 | 33 / 36 | 0.32 |
| 03_thomson_crt | 6.1 | 1.9e-2 | 32/32 | 35 | 35 / 53 | 0.73 |
| 04_slingshot | 3.6 | 9.5e-3 | 27/32 | 155 | 229 / 105 | 0.69 |
| 05_the_wall | 10.0 | 3.5e-3 | 26/32 | 101 | 194 / 286 | 0.66 |
| 06_twin_beams | 6.9 | 7.5e-3 | 29/32 | 87 | 128 / 133 | 0.86 |
| 07_reflectron | 7.9 | 4.5e-3 | 31/32 | 88 | 101 / 222 | 0.83 |
| 08_einzel_lens | 10.3 | 1.0e-3 | 31/32 | 98 | 111 / 1000 | 0.90 |
| 09_hemispherical_analyzer | 6.6 | 1.5e-3 | 19/32 | 120 | 394 / 667 | 0.75 |
| 10_injection | 3.7 | 4.8e-2 | 31/32 | 85 | 98 / 21 | 0.80 |
| 11_collimator | 11.1 | 5.0e-4 | 32/32 | 100 | 100 / 2000 | 0.87 |
| 12_soft_landing | 6.7 | 5.0e-4 | 9/32 | 171 | 1194 / 2000 | 0.74 |
| 13_high_voltage_dome | 3.7 | 1.0e-1 | 32/32 | 11 | 11 / 10 | 0.65 |
| 14_polarized_sphere | 6.6 | 1.2e-2 | 31/32 | 28 | 41 / 83 | 0.61 |
| 15_image_charge | 7.1 | 1.8e-2 | 25/32 | 68 | 180 / 56 | 0.80 |
| 16_deflection_plates | 3.2 | 1.7e-1 | 16/16 | 12 | 12 / 6 | 0.44 |
| 17_power_supply | 0.7 | 2.8e-1 | 16/16 | 3 | 3 / 4 | -0.43 |
| 18_tune_the_lens | 1.4 | 3.7e-2 | 11/16 | 23 | 205 / 27 | 0.30 |
| 19_real_einzel_lens | 5.8 | 1.0e-2 | 16/16 | 78 | 78 / 100 | 0.80 |
| 20_build_a_deflector | 3.2 | 1.7e-1 | 16/16 | 8 | 8 / 6 | 0.47 |
| 21_shielding | 2.6 | 3.2e-2 | 16/16 | 36 | 36 / 32 | 0.68 |
| 22_fast_lane | 7.3 | 1.3e-2 | 32/32 | 65 | 65 / 80 | 0.60 |
| 23_beta_spectrometer | 6.6 | 1.5e-3 | 20/32 | 131 | 371 / 667 | 0.77 |
| 24_first_coil | 3.6 | 3.3e-1 | 32/32 | 17 | 17 / 3 | 0.44 |
| 25_dempster | 5.8 | 7.7e-2 | 32/32 | 21 | 21 / 13 | 0.56 |
| 26_wien_filter | 11.6 | 1.5e-3 | 28/32 | 77 | 134 / 667 | 0.63 |
| 27_first_magnet | 3.3 | 8.1e-2 | 32/32 | 14 | 14 / 12 | 0.40 |
| 28_calutron | 6.5 | < 1.5e-3 | 4/32 | 124 | 2924 / 667 | 0.70 |
| 29_build_wien_filter | 18.2 | 2.5e-3 | 14/32 | 130 | 644 / 400 | 0.85 |
| 30_stern_gerlach | 3.2 | 8.2e-2 | 32/32 | 10 | 10 / 12 | 0.66 |
| 31_stray_field | 3.8 | 5.0e-2 | 32/32 | 19 | 19 / 20 | 0.62 |
| 32_mains_hum | 12.3 | < 1.5e-3 | 18/32 | 96 | 407 / 667 | 0.87 |
| 33_earths_field | 9.6 | 5.0e-4 | 17/32 | 184 | 537 / 2000 | 0.82 |
| 34_rf_kick | 4.3 | 4.3e-2 | 32/32 | 30 | 30 / 24 | 0.27 |
| 35_rf_separator | 11.0 | 6.0e-3 | 29/32 | 66 | 107 / 167 | 0.76 |
| 36_tune_the_rf | 7.9 | 1.5e-3 | 26/32 | 130 | 223 / 667 | 0.76 |
| 37_streak_camera | 14.2 | < 1.5e-3 | 16/32 | 206 | 606 / 667 | 0.84 |
| 38_synchrotron_light | 1.0 | 5.0e-1 | 32/32 | 2 | 2 / 2 | -0.05 |
| 39_two_stages | 3.6 | 6.0e-2 | 32/32 | 74 | 74 / 17 | 0.70 |
| 40_space_charge | 3.6 | 5.2e-2 | 32/32 | 40 | 40 / 19 | 0.71 |
| 41_stern_gerlach_beam | 3.2 | 3.3e-2 | 27/32 | 38 | 113 / 30 | 0.67 |
| 42_collimated_beam | 11.1 | < 1.5e-3 | 31/32 | 95 | 107 / 667 | 0.88 |
| 43_velocity_selector | 11.6 | < 1.5e-3 | 17/32 | 143 | 496 / 667 | 0.68 |
| 44_beam_preparation | 11.7 | < 1.5e-3 | 9/32 | 215 | 1237 / 667 | 0.86 |
| 45_relativistic_beam | 7.1 | 5.1e-2 | 31/32 | 28 | 41 / 19 | 0.69 |
| 46_chromatic_aberration | 6.8 | 9.5e-3 | 30/32 | 76 | 103 / 105 | 0.88 |
| 47_real_analyzer | 6.0 | 1.9e-2 | 32/32 | 50 | 50 / 53 | 0.73 |
| 48_crt_earth_field | 6.7 | 7.0e-3 | 26/32 | 72 | 164 / 143 | 0.74 |
| 49_beam_pipe | 6.6 | 4.3e-2 | 32/32 | 44 | 44 / 23 | 0.80 |
| 50_calutron_space_charge | 6.7 | 2.5e-3 | 15/32 | 123 | 576 / 400 | 0.82 |
| 51_soft_landing_current | 6.3 | 1.0e-3 | 24/32 | 195 | 328 / 1000 | 0.76 |

Notes:

- The introduction levels (1, 10, 13, 16, 17, 20, 24, 27, 30, 31, 34, 38, 39, 40, 45) are meant to be easy. Level 38 is a
  concept level (8 choices); harder radiation levels are to come.
- Under disturbances (chapter 10) a placement counts as a solution only if every flight
  arrives, verified. The distance objective sums over all flights.
- Dempster (25) has in effect one physical parameter (the accelerating voltage, which sets
  all radii). It is kept as the concept level for 180° focusing.
- The Calutron (28) is the hardest: no random placement out of 2000 solved it, and the
  search succeeded in 4 of 32 runs. Its difficulty comes from the magnets' strongly
  non-uniform dipole fields.
- Power supply (17) and Tune the lens (18) are knob levels with tiny configuration spaces
  (5 and 25 settings). The analysis caches the result of every distinct placement, so
  such levels cost only their distinct placements; before the cache (and before it knew
  about supplies), analysing them did not terminate.
- Real instruments (46–51) build an earlier level's reference setup in as fixed elements and
  add one real effect the ideal design ignored. `scripts/levels.py` checks at build time that
  the built-in design alone fails, and the test `realistic_levels_break_the_idealised_design`
  keeps checking it. Measured with the built-in design alone: level 46 one ray 7/8 arriving;
  47 the highest energy 7/8; 48 both beams lost; 49 rejected
  (wrong direction); 50 the light isotope 3/5; 51 0/8. Tuning for playability (no search
  solved the first versions): 48 first had to work facing north, south and shielded with
  one-cell spots, and at most 1 in 32 searches solved it even at 40 % of the field; a tube
  is adjusted where it stands, so the level now has one orientation. 51 had 12 ions at the
  full charge scaling (0/32); now 8 at half.
- Beam levels (40–44) are harder versions of single-particle levels (1, 30, 11, 26, 39): a
  placement solves them only if the required share of the beam arrives, verified. Their
  objective for the search sums the detector distances (plus gate shortfalls) of the
  closest missing particles. Relativistic beam (45) includes radiation reaction because
  its particles would otherwise neglect 5.3e-6 of their energy (above the 1e-10 allowed);
  the test `radiation_levels_need_radiation` accepts that reason for beams. Beam preparation (44) combines the collimated beam (42) with
  a second stage; the search found its reference in 44 s and solves it in 9 of 32 runs.
- The analysis logs its progress to stderr (`[analyze]` lines: phases with timing, distinct
  placements flown, the slowest evaluation, and every evaluation slower than 1 s with its
  placement). The log found why the velocity selector's analysis ran for hours (placements
  trapping the ions for 1e6 steps, and a beam bug); searches now stay within the step
  budget of the cost meters.
- Sorting station (52) is the first finale (docs/CURRICULUM.md). Its reference is built stage
  by stage (`solve_stage` in `scripts/levels.py`): stage 1 (the lens, at most 3 charges left
  of the gate) needs 2 charges; stage 2 (at most 4 charges right of the gate, with stage 1
  fixed) needs 4. The build checks that the search finds no solution with 4 charges
  anywhere (`FINALES`), and the whole-level search with up to 8 charges finds none either
  (0/4 runs): the level is meant to be solved in stages, as it was built. Versions that asked
  for parallel arrival (within ±6°, ±10° or ±15°) had no staged solution; five-cell spots
  were solved with 3 charges in stage 2 (too easy).
- `analyze --fewest` adds the fewest elements that solve a level: the reference solution's
  count, unless the search finds a verified solution with fewer (exhaustive for one
  element, annealing with 64 restarts of 600 iterations for 2 up to the reference count
  minus one). A failed search is not a proof, so the number is an upper bound.
  Measured 2026-09-28: 2 elements for 7, 8, 11, 18, 19, 28, 29, 32, 33 and 42; 3 for 44; 6
  for 52; 1 for every other level. Three references use more elements than needed: the
  Wien filter (26, reference 4) and the velocity selector (43, reference 4) are solved
  with one charge, and Build a Wien filter (29, reference 6) with two elements. These
  numbers are the starting point of the curriculum restructure (docs/CURRICULUM.md).
- The search is a rough stand-in for a player. A player who understands the physics (for
  example the Wien condition v = E/B) should do much better than the search on the later
  levels, and a player who doesn't should do much worse.
