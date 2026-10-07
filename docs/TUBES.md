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
- **Next:** the rectifier (AC anode: a driven electrode in a tube), the triode (a grid of
  small electrodes), then the circuit coupling (the anode current through a load).
