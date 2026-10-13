# Inside the components, step 1: vacuum tubes — plan of the stage

`SPEC.md` ("Inside the components"): thermionic emission from a cathode, space-charge-limited
current validated against the Child–Langmuir law, the vacuum diode (rectification), the
triode (a grid controls the current: amplification), later magnetron and klystron. This
note is the plan, written before building anything (as `docs/CIRCUITS.md` was). Questions
for the owner are marked **(owner)**; nothing is built until they are settled.

## What the engine has and lacks

- **Has:** beams of interacting particles (`beam.rs`: pairwise, exact Coulomb for c = ∞,
  quasi-static or retarded for finite c; collisions; per-particle events); metal by
  boundary elements (`bem.rs`: box electrodes, the Maxwell capacitance matrix, factorized
  once per geometry); spheres with exact Kelvin images, including a particle's own image
  (`conductor.rs`, `self_field`); circuits driving electrodes, one way (§2.10).
- **Lacks:**
  1. *Particles' induced charges on box electrodes.* Today they are neglected with a
     bound (`IMAGE_FORCE_LIMIT`, 1e-10, §2.7). In a tube they are the physics: the
     cathode's induced charge is what limits the current.
  2. *Continuous emission:* particles appear over time at a cathode, at a rate the space
     charge limits. Shots launch at fixed times only.
  3. *Many particles:* a tube's steady state holds hundreds to thousands of electrons in
     flight. Pairwise forces cost O(N²) per evaluation.
  4. *The particles' back-action on circuits:* the anode current through a load (Shockley–
     Ramo), which the triode's gain needs (the circuits are one-way now).
  5. *Goals on currents and collected charge* rather than single particles in detectors.

## Model

- **Induced charges (two-way BEM).** At each evaluation: the particles' potential at the
  panels (N·P), the surface charge that cancels it from the existing factorization (a
  back-substitution, P²), its field at the particles (P·N). Exact for the BEM's
  discretization, conservative (energy with `½ q φ_induced` conserved), and the same
  machinery gives the Ramo currents (the rate of change of each electrode's induced
  charge). Spheres keep their exact images.
- **Emission, space-charge-limited first.** An emitter: a face of an electrode (or a
  sphere's cap) split into patches. Each step, each patch emits the charge that brings the
  normal field at its surface to zero (Gauss's law: the standard space-charge-limited
  emission of PIC codes), as macroparticles of weight w: a macroparticle carries charge
  `w q` and mass `w m`, so it follows an electron's trajectory exactly; only the
  graininess of the space charge depends on w (a stated approximation, measured by
  halving w). Thermionic emission (Richardson–Dushman, Maxwellian launch velocities, the
  virtual cathode) later, as a second emitter kind.
- **Many particles.** First pairwise (exact, the existing beam engine, c = ∞ for the
  tubes: electrons at a few hundred volts move at ~0.03 c, so retardation and magnetic
  forces are ~1e-3 corrections, measured and stated). Then, for interactivity, the mean
  field of SPEC's plan (charge deposited on a grid; the 3D Green's function convolved by
  FFT over sources confined to the plane), validated against pairwise on small N.
- **Two-way circuits.** The Ramo current into each electrode's node of the circuit;
  coupled by splitting (the circuit advanced over each particle step with the step's
  current, RADAU5 as now), valid while the circuit's time constants are long against the
  step (checked per level). Validated by the energy balance: the sources' work = Joule
  heat + the particles' kinetic energy gained + the change of field energy.

## The slice **(owner)**

The game is a slice z = 0 of a 3D world: electrodes extend in z, particles fly in the
plane. Electrons confined to the plane form a *sheet* between the electrodes, not the
volume-filling flow of the textbook diode, so the textbook prefactors do not apply to the
game's tubes. Proposal:

- validate the engine in 3D (it is 3D already: particles may leave the plane) against the
  textbook laws, below;
- build the levels in the slice, with the sheet tube stated as the model ("Physics model
  and its limits"), and its own current law measured and checked: the exponent 3/2 is
  exact in any geometry (the steady space-charge equations are invariant under
  `V → λV, J → λ^{3/2} J`), the prefactor against an independent steady-state solution
  of the sheet (an integral equation, by Wolfram).

The alternative, volume-filling electrons in the slice levels, needs 3D views; it would
wait for the 3D stage.

## Validation (planned tests, references fixed before measuring)

- **V1 induced charges:** a charge before a large grounded plate: the image force
  `q²/(2d)²` as the plate grows (the BEM's measured discretization error, E1–E2); Green
  reciprocity as Z1; energy conservation of a charge oscillating between two electrodes.
- **V2 Ramo:** a charge crossing a parallel-plate gap at speed v induces `I = q v/d`
  (infinite plates; finite plates to their measured edge correction).
- **V3 Child–Langmuir (3D, planar):** wide plates, a circular emitting patch of radius R
  at gap D: the 1D law `J = (√2/9π) √(q/m) V^{3/2}/D²` (in k = 1 units) at the centre
  as R/D grows, and the total current against Lau's 2D result `I/I_CL = 1 + D/(4R)`
  (Lau 2001; the strip emitter's coefficient, about 0.31, to check against the paper).
- **V4 Langmuir–Blodgett (3D, cylindrical):** a cathode wire inside a cylinder: the
  current per length with β² from its series `β = μ − 2μ²/5 + 11μ³/120 − 47μ⁴/3300 + …`,
  μ = ln(r/r_c) (Langmuir & Blodgett 1923).
- **V5 scaling:** the steady current ∝ V^{3/2} in the slice (any geometry) to the
  statistical error of the macroparticles; and the sheet's prefactor against the
  Wolfram solution.
- **V6 circuit energy balance** with a load resistor, as above.

## Levels (sketch)

- *Introduction:* space charge chokes the current (a diode, the anode's slider: the
  current grows as V^{3/2}, not V); the diode rectifies (AC on the anode: current flows
  in one half-cycle only).
- *Intermediate:* the triode: a grid's potential (a slider, symmetric range) controls the
  anode current; amplification: a small AC on the grid swings the anode voltage across a
  load by more.
- *Master:* later (magnetron or klystron), with the circuits' and beams' pieces.
- Goals needed: the current into an electrode (averaged over a window, as the radiation
  receiver), the charge collected, or an amplitude at a node of the circuit **(owner)**.

## Questions for the owner

1. **(owner)** The slice: sheet tubes in the levels, validated in 3D (proposed), or wait
   for 3D views?
  - depends if that changes physics. At this stage I demand that there should be only 3d physics, and only physics that have symmetry that allows to be 
2. **(owner)** Pairwise first (exact, slower, a few hundred electrons interactive?) and
   the mean field after, or the mean field from the start?
  - Whatever makes it interactive. If exact is good enough to make it interactive then okay, if not start with mean field as preview.
3. **(owner)** Space-charge-limited emission first, thermionic later?
  - sure
4. **(owner)** Macroparticles (weight w) as a stated approximation?
  - sure
5. **(owner)** Two-way circuit coupling by splitting (checked per level), or the fully
   implicit system (particles and circuit in one RADAU5 system, several times slower)?
  - again, if it is fast enough for real time, make it more accurate, if not then make it approximate as preview, maybe automatic check after creating in sandbox, so that those levels can be also done easily by hand.
6. **(owner)** The new goal kinds (current, charge, circuit amplitude).
  - sure, whatever is fitting the framework or just need relatively small extension.

## Decisions (owner, 2026-10-07)

1. **Geometry: z-invariant tubes.** Only real 3D physics, justified by a symmetry. The
   tubes are translation-invariant along z (the owner: every thin strip of the cathode
   behaves the same, and the strips glued together form the full 3D cathode), the
   symmetry the textbook diode theory itself assumes. Consequences:
   - a separate geometry mode for such levels: a particle is a line of charge along z
     (charge per unit length), fields of line charges (potential `−2λ ln r`, force
     `∝ 1/r`), electrodes as 2D cross-sections solved by 2D boundary elements;
   - everything in such a level must be z-invariant (electrodes, wires along z, uniform
     or solenoid B along z); point charges, magnets and in-plane coils break the
     symmetry and are refused there;
   - the references then hold exactly: Child–Langmuir (planar), Langmuir–Blodgett
     (coaxial), Lau's 2D law for a strip cathode of finite width;
   - retardation and magnetic forces between electrons (~v²/c², 1e-3 at tube energies)
     neglected, measured and stated, as for beams.
   Finite 3D tubes wait for the 3D stage.
2. **Speed decides the method:** exact pairwise interaction where it is interactive,
   otherwise a mean field as the preview (with the exact model in the verification).
3. **Emission:** space-charge-limited first, thermionic later.
4. **Macroparticles** (weight w: charge and mass per unit length) as a stated
   approximation.
5. **Circuits:** the accurate (fully coupled) model if it runs in real time; otherwise an
   approximate preview with the accurate model in the verification (an automatic check
   after a sandbox edit), so that such levels can be built by hand too.
6. **Goals:** current, collected charge and circuit amplitude, as far as they fit the
   framework or need only a small extension.

Sources: Lau, "Simple theory for the two-dimensional Child–Langmuir law", PRL 87,
278301 (2001), and Lau et al., "On the Child–Langmuir law in one, two, and three
dimensions" (2023, arXiv:2307.14552); Langmuir & Blodgett, Phys. Rev. 22, 347 (1923),
the series as quoted by Greenwood et al. (2016) and Zhang et al., "100 years of the
physics of diodes" (2017).

## Progress

- **Step 1 done:** z-invariant electrostatics (`physics::zinv`; tests V1–V3: segment
  integrals, coaxial capacitance, image force, second order).
- **Step 2 done (physics):** space-charge-limited emission with line-charge macroparticles
  (`physics::tube`; test V4: the coaxial diode against Langmuir–Blodgett, first-order
  convergence, Richardson-extrapolated 0.5 %). PHYSICS.md §2.11.
- **Step 3 done:** Ramo currents (test V5, second order).
- **Step 4, first part, done:** the level mode (`level::tube`: format, setup checks,
  preview and verified runs, verdict, cost; game: worker, panel, particles, potential map
  with space charge; sandbox: the Tube section; generator: check, solve, window) and the
  first level, Vacuum diode (V^{3/2}). PHYSICS.md §2.11.
- **The owner's review of Vacuum diode (2026-10-07):** too short to watch, no speed
  control, "more like a demo than a level". Done: runs streamed to t_max (200 for the
  diode) with a speed slider; live potential and electric maps; ion emitters (the bipolar
  diode, V6), projectiles that interact with the electrons and land on the electrodes
  (V7), the player's free charges thrown into tubes, a uniform B along z (crossed fields,
  V8), temperature-limited cathodes (V9); and two puzzles where placing plates matters:
  Through the cloud (an ion across a working diode, two goals that pull against each
  other) and Electron gun (focusing a space-charge-spreading beam onto a small
  collector). Magnetic and Poynting maps wait for the circuit coupling (B_z depends on
  the in-plane return path of the tube's current).
- **The player's solenoid** (a supply sets B along z) and Sideways (level 100): E × B
  drift above Hull's cutoff into a collector beyond the gap's end.
- **Circuits in tubes** (PHYSICS.md §2.11, V10–V11): electrodes as circuit nodes (a
  DC, sine or pulse source through R, C to ground; ideal and floating), solved implicitly
  each step with the particles' induced charge, so the tube and its circuit act on each
  other; one source may drive several electrodes (a grid's wires); goals on electrodes'
  potentials (mean, ripple, oscillation amplitude), plotted in the panel. Level:
  Rectifier (103), a full-wave rectifier whose two AC anodes the player drives.
- **Physics found on the way** (PHYSICS.md §2.11): carriers coming back to their emitter
  are absorbed on reaching their release layer (V12; they used to bounce there for
  ever); every tube level is enclosed by a grounded envelope at the board's edge (V13;
  open boards carried a net charge that depended on the logarithm's reference, and with
  it the far field: a uniform shift of the potentials changed a triode's losses 2.7-fold);
  the push's far field from a multipole series (V14) and the induced charge solved once
  per step, which pay for the envelope. The existing tube levels were re-tuned.
- **The triode** (level 104): a grid of five wires on one biased sine source before the
  cathode, the anode through a load; the player biases the grid for a gain between 2.2 and
  3.8. Two obstacles on the way, both now handled (PHYSICS.md §2.11): the anode's own
  capacitance caps the gain at `g_m/(ωC)` (a short anode, a slow signal), and the preview
  erred by +29 % at the narrow cathode–grid gap (the release distance's first-order
  error: the emitters are now meshed at a quarter, V15). On the way the macroparticles'
  weight was found to depend on the slider's range through its floor (a grid-screened
  cathode got 27 macroparticles instead of 152): the weight now comes from the emitting
  face's charge over the sources' period, with a far lower floor.
- **The owner's review (2026-10-12):** Through the cloud's flow at 30 V was too thin for
  the ion to cross a real cloud (now 60 V, the ion at 1.8, a plate below the gap: the
  original idea works again); the potential map showed blue and red with a sharp wall
  (its scale came from the potentials at t = 0, all 0 with AC anodes: fixed, and the
  tube maps have gain and range controls now). The triode's remaining preview error was
  the particles' softening (half a segment: the narrow gap's space charge smoothed): now
  an eighth (PHYSICS.md §2.11).
- **State (2026-10-12):** the stage's plan is built through the triode: levels 99–104
  (diode, crossed fields, ion through a diode, gun, rectifier, triode), tests V1–V16,
  all merged to main.

## Next, and questions for the owner

1. **Magnetic and Poynting maps** (requested together with the electric map, which is
   done). The tube's currents flow in the plane, so its magnetic field is B_z, a stream
   function of the total current (conduction plus displacement). It therefore depends on
   the path by which the current returns through the external circuit, and a z-invariant
   model must draw that path in the plane. **(owner)** Choose:
   (a) straight leads from each electrode to the envelope;
   (b) leads the level designer draws;
   (c) only the particles' own Biot–Savart field (first order in v/c, without the
   electrodes' surface currents), labelled as such.
2. **A bipolar level.** Ions emitted from the anode neutralize the electrons' space charge
   and raise their current (Langmuir: ×1.86 in a planar diode; V6 validates the coaxial
   bipolar flow). A puzzle needs a trade-off, e.g. more electron current while the ion
   current onto the cathode (its wear) stays limited. That probably needs a current goal
   per species, an extension of the format. **(owner)** Is that goal kind wanted?
3. **The stage's master level** (the plan above: a magnetron or klystron). Both fit the
   z-invariant model: a magnetron (a cathode, an anode of vanes, B along z) and a
   sheet-beam klystron (a beam crossing gaps). A magnetron needs resonators: in the
   quasi-static model, LC tanks between neighbouring vanes, driven by the electrons' Ramo
   currents. A klystron's output gap needs one to be efficient: without it, the gap's
   capacitance shunts the signal, as the anode's did in the triode. The tube circuits have
   no inductors yet. Validation: an LC node against its exact solution, the energy
   balance (V16), and for a magnetron the oscillation's onset against the Buneman–Hartree
   condition. **(owner)** Approve before building: which device, and inductors as a new
   circuit element.
4. **Tube previews take 0.7–2.3 s.** Five of the six levels exceed the preview budget
   (`PREVIEW_LIMIT`, 1 s, in `crates/level/src/cost.rs`); the display streams, so only
   the preview's numbers wait. Not profiled since the envelope and the finer softening.
   **(owner)** Make them faster, or give tubes a budget of their own?

## Decisions (owner, 2026-10-12)

1. **Maps:** straight leads from each electrode to the nearest point of the envelope,
   drawn on the map, the convention stated in the model notes. B_z from the total current
   (conduction plus displacement), zero outside the envelope; S = E × B_z ẑ/μ₀. Not (c):
   the particles' own Biot–Savart field is not divergence-free beside electrodes (their
   surface currents are of the beam current's size), so it is wrong there at first order.
   Leads drawn by the designer (b) can come later as an optional override.
2. **Per-carrier current goals:** yes (`carriers`: all, electrons or ions), and a bipolar
   level where the player places ion-emitting strips on the anode: an electron current
   above the unipolar Child–Langmuir value with the ion current onto a protected cathode
   segment limited.
3. **Inductors and both devices.** Inductors as a circuit element, with an
   energy-conserving scheme for the L and C states (backward Euler would damp an LC tank
   by ~ω²dt/2 per unit time, a numerical loss that would set the resonator's Q); tests: an
   LC node against its exact solution over many periods, V16 with ½LI². Then a klystron
   (an intermediate level; references: ballistic bunching, the optimum X = 1.841, and the
   reduced plasma frequency's debunching) and a magnetron as the stage's finale (crossed
   fields, space charge, LC tanks, the solenoid; Buneman–Hartree). The owner saw no point
   in choosing one. The smooth-anode magnetron's instability is seeded by the
   macroparticles' noise (V8): if the oscillation's onset does not converge between
   preview and verification, the tanks are primed (a small injected signal, as real
   magnetrons are).
4. **Preview cost:** profile and speed up first; decide on a tube budget with the
   measured numbers (if one is needed: time to the first streamed frame within the
   current limit, and a stated separate limit for the preview's numbers).

Order: 4, 2, inductors and the klystron, the magnetron, 1.

- **4 done (2026-10-12):** profiled per phase, the particles' fields (50 %) and the
  induced density (45 %) took the step's time. Both now run in parallel (each sum whole by
  one thread, collected in order: bit-identical), and the induced density uses the
  collocation matrix's inverse (a row-parallel product) instead of the serial triangular
  solves. Previews 0.24–0.63 s (were 0.7–2.3 s), verifications 0.9–2.6 s (were
  4.4–13.9 s): all within the existing budget, so tubes need no budget of their own. The
  inverse changed the rounding, and with it the tube levels' results by up to 0.74 % (the
  macroparticles' noise, inside every verdict's bar; PHYSICS.md §2.11, Speed). Further
  speed, if the magnetron needs it: a 2D fast multipole (treecode) for the particle sums.
- **2, the goal kind, done (2026-10-13):** a tube level has a list of current goals
  (`goals`; older files' single `goal` still loads), each counting the cathode's carriers
  or the ions (`carriers`), with the verdict, the panel (one plot per goal), the sandbox
  editor (add or remove goals, a carriers choice), `generator check`/`tube-scan` and
  `levels.py` (`current_goal`, `tube(..., goals=, ions=)`). The golden hashes did not
  change. First numbers for the bipolar levels (the game-scale diode with its anode
  emitting ions): at q/m = 0.01 the electron current rose ×2.2 (finite plates, not
  Langmuir's infinite 1.86) but the run must reach t ≈ 80 for the ions to settle, a
  preview of 2.6 s; with q/m = 0.04 (the mass enters only the timing and the ions' space
  charge, not the paths or the factor 1.86) the ions cross five times faster than the
  electrons, and the ion current into the cathode came to 0.21 of the electrons' (the
  planar law: √(m_e/m_i) = 0.2).
- **Inductors done (2026-10-13):** a tube circuit's branch may have a coil in series with
  its resistance and may return to another electrode instead of ground (`inductance`,
  `return_to`; editor, setup checks, pairs.py `tube_tank`). Coil branches by the
  trapezoidal rule (with the nodes' charge balance the implicit midpoint rule: an LC tank
  keeps its energy and is only detuned by (ω dt)²/12), resistive ones by backward Euler as
  before (no level changed). V17 (an LC tank to rounding), V18 (a tank between two
  electrodes), V16 with a coil (PHYSICS.md §2.11).
- **The bipolar levels (owner, 2026-10-13):** the approved ion strips with a protected
  cathode segment had a flaw (the ions that neutralize the cloud by the cathode land on
  it beneath the cloud: no real decision), so the owner chose an introduction plus a
  magnetically insulated ion diode. Done: rods and rings in tube levels (PHYSICS.md
  §2.11, test `rods_and_rings`) and **Bipolar diode** (level 101: the anode's face emits
  ions, q/m = 1/16; electrons ≥ 2.0 with ion wear ≤ 0.75 solve for 22.3–29.7 V; without
  the ions the electron goal would need 41 V).
- **The Ion diode, not shipped (2026-10-13; a question for the owner):** an inverted
  coaxial diode (the ring the cathode, the rod the anode emitting ions) in the player's
  solenoid. Its physics comes out as theory says (preview, 26 V, ring 10, rod 2, t =
  20–30): the electrons leak 14 into the rod below Hull's cutoff (0.75), the leak stops at
  ~2.3 times it, and the electron sheath raises the ion current from 3.8 (strong field) to
  11 (just above the cutoff); the sheath keeps spreading towards the rod (still rising at
  t = 80–90, the impedance collapse of real ion diodes). But the boost is dominated by the
  macroparticles' noise, far above a real tube's shot noise: four times as many changed
  the ion current at 28 V, B = 1.8 from 8.04 to 6.93 (the verification's 6.22), and the
  preview runs 7–17 % above the verification where the sheath acts (4 % in the strong-field
  limit; first-order convergence, refine 1, 2, 4: 4.94, 4.36, 4.08). With verdict bars of
  15–25 % no placement of the goals verified cleanly at its reference while keeping the
  slider rule. The two Hull cutoffs are robust, but with a realistic mass ratio their
  window (√(m_i/m_e) b/a apart) is far wider than the rule allows (the solving interval
  must stay between 0.55 and 0.95 of the range). Options: (a) reduce the noise first
  (quieter emission, more macroparticles with a faster field sum: a treecode), which the
  magnetron needs anyway, then return; (b) an artificial mass ratio of ~4 so that the two
  cutoffs bound the field; (c) relax the slider rule for this level's field. The level's
  design is kept in `scripts/levels.py` (`ion_diode`, not in the curriculum).
- **Klystron done (2026-10-13; level 106):** a weak beam (ballistic bunching) from a
  Wehnelt gun, a buncher on a sine and a tuned catcher (a tank) 21 further; the player sets
  the cathode's −V₀ (the beam's speed) and the drive V₁. Against Webster's small-signal
  theory from the level's own fields: the optimum drive within 13.7 % of the prediction,
  the rising side's shape to 0.006; past the optimum the two-gap buncher and the finite
  catcher depart from the thin-gap J₁ (0.822 against 0.726 at 1.5 times the optimum's
  drive), measured and stated. The amplitude alone hardly depends on V₀ (the buncher's
  coupling grows with the speed, offsetting the drift angle), so the level asks for the
  output's phase too (a new potential goal, `phase`): 7.5° per volt, phase pushing. A flat
  optimum and the verdicts' bars (amplitudes 5 %) leave no narrow window at saturation, so
  the goal is a klystron backed off below it, with an overbunched island as well. On the
  way: a per-level macroparticle rate for temperature-limited cathodes (`particle_rate`),
  and the amplitude measure corrected (fitted together with the trend; the old fit took
  part of the sine with the line).
- **The owner's review of the Klystron (2026-10-13):** "noise" in the picture, and no
  indication of which electrodes a slider drives. The maps themselves were clean on this
  machine (checked at gain 0 and 2 and in the electric map); the noise was the beam: every
  segment of the cathode started its emission account at zero and, emitting at the same
  rate, all fired together in volleys, which drew regular arcs (and made a coherent
  fluctuation of the current: V9's quantization error was −0.9 %/+0.3 %, now 1e-4). Fixed
  by a quiet start (PHYSICS.md §2.11), and many macroparticles are drawn fainter. The
  supplies are numbered: each slider's row starts with a badge whose number is drawn on
  the map beside every electrode it drives (linked ones too) or at a coil's or solenoid's
  terminal; hovering a slider (or moving it), or an electrode on the map, outlines its
  electrodes in cyan (`EM_FOCUS=<n>` holds one in a capture).
- **Noise, stage 1 (owner's choice (a), 2026-10-13).** Measured on the ion diode's
  scratch level (28 V over t = 20–30, ions onto the ring at B = 1.8 / the strong-field
  B = 3): the macroparticles' number dominates: ×1, ×2, ×4, ×8 gave 8.21/5.05, 7.42/4.91,
  6.90/4.75, 6.66/4.66 at the preview's mesh (converging toward ~6.4); the softening alone
  ×2, ×4 only 7.78/5.00, 7.68/4.88; the mesh and step another first-order ~6 % (refine 2:
  6.88/4.68; with ×4 particles 6.27/4.48). The quiet start alone did not change it (8.21
  against 8.04). So the sheath needs ~8 times as many macroparticles, ~30 times the cost
  with direct sums: a tree code. Its first step, done: the macroparticles are uniformly
  charged rods (radius an eighth of a segment), exactly Coulomb outside, so that far
  groups can be summed as multipoles exactly (PHYSICS.md §2.11; the calibration diode now
  converges at first order, ratio 0.51). `generator tube-scan --particles F --softening F`
  vary the two alone (`TubeResolution`).
- **The owner's review, continued (2026-10-13):** the "light blue noise over the blue
  potential" on the right of the Klystron was its beam's macroparticles, spread to the
  whole aperture (physically: the gun's lens and the buncher's fringe fields; a slow beam
  spreads most), drawn as dots. The particles are now drawn by default as a glow of their
  charge's density (a View choice: density, dots, off). And a general rule from the owner:
  at the lowest playback speed the animation must stay smooth, not a slideshow. Flights and
  beams already interpolated along their paths; tube levels jumped from step to step (at
  speed 0.05, 1.5 pictures per second): their particles now carry identities and the
  display interpolates between steps (particles, electrode charges, potentials, the
  projectiles along their tracks).
- **Noise, stage 2 (2026-10-13):** the tree code (PHYSICS.md §2.11) makes 8 times as many
  macroparticles affordable at about 3 s per preview (6.9 s with direct sums): the ion
  diode's sheath then differs by 8.3 % between preview and verification instead of 16 %,
  a verdict bar of ~12 % instead of ~24 %. The remaining difference is the mesh's and the
  step's (first order); halving the preview's mesh would cost ~4 times more. So magnetized
  sheaths stay the noisiest, costliest physics of the stage: levels built on them need
  goals that tolerate ~12 % and previews of ~3 s, or observables less sensitive than the
  sheath's current (a magnetron's frequency, set by its tanks; its onset above the
  Buneman–Hartree condition, with priming).
- **Ion diode done (2026-10-13; level 104, the owner's B before A):** shipped with 8 times
  the macroparticles and emitters meshed twice as finely (two new per-level fields,
  `particles` and `emitter_mesh`): the sheath converged in the particles' number, and the
  emitters' release layer was the rest of the preview's error. Bars ~9 % at the reference,
  a preview of 3.4 s (the stage's costliest level), a verification of 13.6 s.
- **Magnetron done (2026-10-13; level 108, the stage's Master tier):** straps (one node
  over several electrodes, test V19) make the strapped vanes; the goal is the tank's
  saturated swing (at least 10 over t = 90-120), robust where the start-up time is not
  (PHYSICS.md §2.11, the magnetron). Frequency goals were not needed: the strapped π mode
  rings at 3.03-3.11 wherever it oscillates. Left for later: a frequency measure (pushing,
  a tuning plunger), the output's power as a goal, and the spokes' picture in the owner's
  review.
- **Owner review 2026-10-13 (Ion diode glitchy and slow, its leak unreachable; the
  magnetron's emission must be seen):** done. Ion diode: a grounded envelope removes the
  arc of electrons drawn off the cathode's outer face by the open 2D boundary's field; the
  leak is reachable with either slider (field off: 17.4 at 30 V; at B = 1.9 above ~50 V,
  the supply now ±80) and visible (each species' glow scaled to its own density, the
  goal labelled "the electron current"); the map merges particles above 1500 (~116 frames
  per second instead of ~73 at ~6000 particles); the tree code's bit-identical speed-ups
  and leaves of 24 bring a preview from 3.4 s to 2.9 s (a fast multipole method was built,
  validated and found slower at these sizes). Magnetron: `tube.output` (the tank's load
  radiates through an antenna on the axis), drawn as the exact 2D retarded wave of the load
  current beside the tube (PHYSICS.md §2.11). New owner rule: whatever a level says can
  happen must be possible within its sliders.
