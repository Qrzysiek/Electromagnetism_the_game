# Specification: Electromagnetism – the game

Status: **draft v2** (revised from `SPEC_original_pl.md`). Before writing code, Claude Code presents a Stage 1 plan (repository layout, modules, order of work, tests) and waits for approval.

Companion document: `PHYSICS.md` describes the physics that is implemented and how. It must be kept up to date with the code.

## 1. Goal

A puzzle game built on accurate classical electrodynamics. The player places fixed charges on a grid so that test particles fly from start point A into detector region B, avoiding obstacles (charges placed by the level).

Priorities, in order:
1. **Physical accuracy.** No simplifications of the equations except explicitly documented and controlled approximations. Numerical error must be small enough that it never decides victory or defeat (see §2.3).
2. **Playability.** Puzzles are hard but solvable by a human. A solution tolerates small imprecision, and the live preview gives visible feedback as the player gets closer.
3. **Performance.** An interactive live preview, and a generator fast enough to search many candidate levels.

Audience and commercial success are not goals. Classical electrodynamics is the backbone. Quantum effects are out of scope for the foreseeable future.

## 2. Physics

### 2.1 Model
- The engine is **three-dimensional from the start**. 2D mode uses the same engine (§4).
- Field sources: charges with a finite radius. The field outside the sphere is exact Coulomb. Touching the sphere means the particle is lost. There is no artificial softening.
- Test particles: defined per level by charge `q`, mass `m` and radius. They don't have to be elementary particles. "Scaling up the particles" (larger `q`, e.g. charged microspheres) is how the level controls how strongly particles affect each other (§2.4).
- Equation of motion, relativistic: `dp/dt = q(E + v×B)`, `p = γmv`. The state is `(x, p)`, so `|v| < c` holds by construction.
- Static field: exact Coulomb superposition in f64 with a fixed summation order.
- Units: dimensionless internally, with an explicit, documented conversion to SI.

### 2.2 Numerical integration
- A high-order adaptive integrator with error control and dense output (baseline: Dormand–Prince 8(5,3), "DOP853", Hairer et al.). Boris is **not** used while the fields are purely electric, because an adaptive step destroys its structure-preserving properties. It may be reconsidered with magnetic fields.
- Collisions with charge spheres, detector entry and leaving the map are found as **events** on the dense-output interpolant: root finding within the step, so a step can never jump over a sphere. The step size is additionally limited by the distance to the nearest surface.
- Conserved quantities (energy `γmc² + qφ` in static fields, and angular momentum where applicable) are monitored as diagnostics. They are **never** used to correct the solution.

### 2.3 Outcome verification ("no numerical luck")
- Every outcome (hit / miss / lost) is computed at a working tolerance and re-checked at a tolerance at least 100× tighter. If the outcome differs, or the trajectory passes a boundary (sphere, detector edge) closer than the estimated error, the result is **marginal**.
- The live preview uses a fast tolerance and is refined in the background. The UI shows when a result is verified.
- The generator rejects levels whose reference solutions are marginal.
- Sensitivity that is physical (chaotic, near-grazing trajectories) is not numerical error. It is handled by the robustness requirements for levels (§7).

### 2.4 Interaction model ladder
Each level declares its model. The model does not switch automatically mid-simulation. The active model is shown to the player.

A single particle in the static field of fixed charges is treated **exactly**, including relativity, at any speed. Approximations enter only through the interaction between moving particles. Each model has a validity range in `v/c`. A run (level + launch speed) outside the declared model's range is rejected, never silently simulated. Example: for co-moving particles, magnetic attraction cancels part of the Coulomb repulsion, and the net force scales as `1/γ²`. Pure Coulomb waves at relativistic speed would therefore overstate space charge. Relativistic waves require at least the Darwin model, and highly relativistic ones require Liénard–Wiechert.
1. Coulomb interaction between test particles (Stage 2): exact for `c = ∞`.
2. Darwin approximation: magnetic interaction to order `v²/c²` (Stage 6). Beams do not use it: they use 1 for `c = ∞` and 3 otherwise (§3, Beams).
3. Liénard–Wiechert retarded fields of point charges, with Landau–Lifshitz radiation reaction. Exact classical electrodynamics in vacuum (later).
4. FDTD / PIC on a grid, only where materials, cavities or waveguides need it (optional, far future).

### 2.5 Other
- B is a region (detector). An electrostatic field has no stable equilibrium (Earnshaw's theorem).
- Each level defines the particle species and the energy and direction at launch from A.
- Game rules (not physics): leaving the map bounds, or exceeding a maximum flight time, counts as a loss.

## 3. Gameplay

- Charges of the player and the level are fixed and sit on grid nodes. Particle trajectories are continuous.
- Each level has a recommended grid resolution (the one the generator used). The player may refine it by an integer factor (2×, 3×…), which keeps all old nodes. Arbitrary resolution is available in custom mode.
- Each level limits the player's charges: count, sign, allowed magnitudes, and optionally a **placement region** (the electrode region of a real instrument).
- Besides generated levels, a set of levels reproduces real experiments and instruments that point charges can model: Geiger–Marsden scattering, Thomson's cathode-ray tube, an Einzel lens, a hemispherical analyzer (exact: a point charge's field is exactly the field between concentric spheres), a reflectron, and a relativistic beta-ray spectrometer, and the magnetic instruments below. `scripts/levels.py` defines all shipped levels.
- **Curriculum** (`docs/LEVELS.md`): levels are grouped in chapters by how hard the phenomenon is to understand. Every new element gets an easy introduction level first, and within a chapter difficulty and the number of elements needed rise on average.
- **The level is solved** when wave 1 (a single particle) reaches B.
- **Waves** are the endurance and high-score layer: 1 particle, pause, 2, pause, 4, pause, … Later waves feel the particles' interaction more strongly (and, in the future, the response of materials). Pauses let materials relax.
- **High score (open, candidates below; may be combined):**
  - *Wave score:* the highest wave number that arrives completely, with no losses (collision with a charge or another particle, leaving the map, timeout). Tie-breaker: the fraction of the first failed wave that arrived.
  - *Speed score:* the highest launch speed (or energy) for which the configuration still delivers the particle(s) to B. Launch direction, A and B stay fixed. Note: in non-relativistic electrostatics, multiplying the launch energy by `s` is equivalent to dividing all charges by `s`. The challenge therefore comes only from the limits on charge magnitude and, at high speed, from relativity, which breaks this scale invariance. Relativity then becomes part of the gameplay.
  - *Mixed:* e.g. waves at a chosen launch speed, or a two-dimensional score (wave number × speed).
  - To be decided after playtesting.
- Wave particles start with a defined, small spread in position, direction and energy, as a real beam does. Identical initial conditions would be singular.
- Fast feedback: the wave-1 trajectory recomputes live on every change. There is no start button.
- **Multi-shot levels:** a level may contain several shots (particle species, launch and detector each). One setup must deliver every shot to its own detector. The UI switches between shots or shows all at once. This is the natural format for instruments that sort particles by mass, energy or angle.
- **Disturbances** (PHYSICS.md §2.3): fields from outside the arena, i.e. uniform stray E and B and plane waves travelling in the plane (for `c = ∞`, uniform fields oscillating in time). A level lists a finite set of realizations, and every shot must arrive under each of them. The player learns robust design: compromise aiming, imaging (arrival independent of the launch angle), fast transit. Shielding needs conductors and comes with them (§8).
- **Radiation** (PHYSICS.md §2.5, §3.1, §10): levels may include radiation reaction (Landau–Lifshitz); the particle's own field and the antennas' waves can be shown as animated maps.
- **Antennas** (PHYSICS.md §2.4): oscillating electric dipoles in the plane with exact retarded fields, at the level's RF frequency or at a frequency the player chooses from the level's list. Shots may have a launch time, so identical particles launched at different moments can be sorted (RF separator, streak camera).
- **Magnetic elements** (SPEC §9, PHYSICS.md §2.2): magnets (uniformly magnetized spheres standing out of the plane, exact dipole field) as level or player elements, and coils lying in the plane (circle or polygon, exact Biot–Savart) as level elements. Magnet levels: Dempster's mass spectrometer (180° focusing and mass separation), a Wien filter (velocity selector), and the calutron (magnetic isotope separation with player-placed magnets).
- **Harder variants of instruments (started: the real instruments, placed in the arcs of the effects they add).** A level shows the idealised solution of an earlier level, then adds a real effect that breaks it: fringe fields, finite electrode size, stray fields, radiation losses, space charge. The player adds elements to make the whole setup work again. Where patching cannot work, the level instead gives a realistic starting point or asks for a design from scratch.
- **Beams (planned, with Stage 2).** Many particles that interact through their fields. The goal is to control the whole beam, not just one particle. Design (agreed with the owner, details open to playtesting):
  - *Role:* beam levels are harder, less forgiving evolutions of single-particle levels. Examples: a Stern–Gerlach intro with two shots (one per spin state), then a mixed beam of both states that also interact; the Collimator, then the whole beam concentrated; the Wien filter, then a beam with an energy spread from which one speed band is selected.
  - *Species:* charge, mass, magnetic moment and fraction of the beam. A magnetic moment is fixed along ±z (the two spin states measured along z; an unpolarized beam is a 50/50 mixture). In the plane every level source gives B along z, so μ ∥ B: no torque, and the states stay fixed, while the force is ∇(μ·B). Neutral atoms with a moment give Stern–Gerlach. A moving moment in an electric field needs the motional (v×E/c²) and hidden-momentum terms; until those are modelled and tested, moments combined with E fields are rejected by `model_issues`.
  - *Beam parameters:* particle count, energy spread, angular spread, position spread (together the emittance), and bunch length (spread of launch times). Distributions: truncated Gaussian (±3σ) or uniform. Samples are deterministic (seeded quasi-random), so every level has one fixed, verifiable beam. Levels fix the parameters, and the sandbox edits them.
  - *Interaction:* all particles form one ODE system. For `c = ∞` (Newtonian levels) the instantaneous Coulomb interaction, plus the static dipole–dipole interaction of moments, is exact within the model. For finite `c` the exact classical choice is retarded Liénard–Wiechert fields between the particles, with a delay: each needs the others' past. Keeping every step shorter than the light travel time between the closest pair (`h < r_min / c`) makes every retarded time fall in history that is already computed, so the system stays explicit (now: delays shorter than a step are taken from the last step's extrapolated polynomial, PHYSICS.md §3.3). Before launch, particles are taken to keep their launch acceleration (tapered so that their velocity changes by at most 0.1c). Approximations such as Coulomb-only or Darwin at finite `c` would be off by O(v²/c²) (4 % at v = 0.2c) and are not used. The cost grows as N² and is shown by the resource meters, which limit N.
  - *Status:* implemented (PHYSICS.md §3.3): exact Coulomb interaction for `c = ∞`, retarded Liénard–Wiechert interaction for finite `c`, per-particle events, gates and verification, deterministic truncated-Gaussian or uniform sampling; the beam levels. Radiation reaction for beam particles (Landau–Lifshitz in the total field). For finite c the default interaction is quasi-static, an approximation marked in the game: Liénard–Wiechert fields of each source's past continued along the motion it would have in the fields it feels now, taken as uniform (exact and relativistic, in closed form; exact in the velocities and for sources turning in uniform fields, so the beam's radiation fields act in the dynamics; how the fields change along the sources' paths is left out, estimated per particle). The exact retarded interaction gives the verdict (the preview stays quasi-static), and every view switches to the verdict's flight when it arrives; it is also a sandbox option for the preview. It costs about 15× the preview (Relativistic beam: 0.9 s against 61 ms). The particle-field view shows the beam's full retarded field, and optionally the part the quasi-static model leaves out, from the recorded world lines (visual, CPU, about 20 ms per frame for 16 particles). A relativistic space-charge level (Relativistic beam) needs radiation reaction: its particles radiate 1e-9 to 1e-7 of their energy when space charge visibly matters (the fraction grows as (δ/L)³ (r/L)² (v/c)⁵ γ⁷/(γ−1) for a deflection δ over a length L, beam size r), above the 1e-10 allowed for neglected radiation.
  - *Goal:* a verified transmission fraction (for example ≥ 95 % of the beam in the detector); only particles whose arrival is verified count. A beam that becomes chaotic makes the preview and the verification disagree, so laminar designs are enforced by the existing verification.
  - *Collisions:* particles with a radius collide as rigid spheres (elastic, along the line of centres; PHYSICS.md §3.3, tests B18–B20). Opposite point charges (radius 0) would fall into each other, where the point-charge theory breaks down; the level validation rejects them.
  - *Multi-stage levels (beam preparation for a further experiment):* for example, first focus the beam, then separate it. The mechanism is ordered **gates**: pass-through detectors, each with optional direction and energy acceptance, that every flight must cross in order before the final detector. They reuse the detector-acceptance boundaries and their verified margins (PHYSICS.md §6.1). A campaign variant freezes the player's solution of one stage as level elements and adds the next stage.
- **Dynamic particles (implemented).** Particles that move and interact with the shots, with the intent of relative motion and collisions (the owner, 2026-09-28): the level's *free particles* (targets, partners: any species, a launch velocity, an optional detector: a goal that must arrive, verified) and the player's *free charges* (a new element: a node, a charge, a speed and a direction from the level's lists; hardcore: continuous). A free charge shows its velocity as an arrow whose tip is a handle: drag it to set the direction and the speed. The velocity is continuous in every mode (any direction, any speed in the range of the level's speeds; the owner: snapping to the lists was too clunky); the solver still searches the listed speeds and 15° directions. In relativistic levels the arrow's length measures rapidity `c·artanh(v/c)` (additive under boosts along a line; the speed itself for slow particles, unbounded as v approaches c), an option in the panel, which also shows v/c and φ. A level with dynamic particles flies as one interacting system (as beams do).
  - *Later:* annihilation and pair creation as an explicit, marked model option (a rule applied at contact; energy including rest energy leaves as a radiation pulse). Not classical electrodynamics; Bopp's work (as suggested by the owner) to be researched before designing it.
  - *Later (F:polar, docs/JACKSON.md):* neutral particles with higher moments. Induced electric dipoles (polarizability α: force α∇(E²)/2, energy −αE²/2; simple, conservative) and permanent electric dipoles (polar molecules: orientation and moment of inertia, torque p×E; the Stark deflector, Jackson Pr. 4.5). Quadrupoles are physical (CO₂, N₂) but their forces, from field second derivatives, are tiny in practice: listed, not planned.
- **Circuits (planned stage, after beams and gates).** Lumped electronic components coupled to the field world, SPICE-like (modified nodal analysis with an implicit, stiff integrator coupled to the particle dynamics):
  - Electrodes are circuit nodes with the capacitance matrix of their geometry (already computed by the BEM, PHYSICS.md §2.7). Time-varying potentials act through the cached unit solutions: E(x, t) = Σ Vᵢ(t) Eᵢ(x).
  - Coils are inductors with self and mutual inductances from their geometry (Neumann's formula); a changing current adds the induced field −∂A/∂t.
  - Particles are currents: absorbed charge flows into the electrode's node; passing charges induce currents by the Shockley–Ramo theorem, from the same unit solutions.
  - Components: R, C, L and sources exact; diodes (Shockley equation) and transistors (Ebers–Moll or square-law) as stated semiclassical device models with fitted parameters.
  - Validity: quasi-static (circuit much smaller than the wavelength), measured. Tests: RC/RLC transients, diode I–V, energy balance of circuit + field + particles, Ramo (induced current integrates to the delivered charge).
  - *Status (2026-10-06):* the stiff integrator is in place and validated: a port of Hairer & Wanner's RADAU5 (PHYSICS.md §5.4; the same decisions as the Fortran code on three classical stiff problems, accuracy against high-precision references, a circuit as an index-1 DAE). So is the circuit engine on it (PHYSICS.md §2.9: modified nodal analysis of R, L, C, sources and switches; tests Z4–Z9 against closed forms), and its one-way coupling to plates and circular coils (§2.10, tests Z10–Z13), in the level format (`drive` on electrodes and circular coils, with editor controls and plots in the panel) and in the curriculum's arc 5 "Circuits" (Charging a plate, Chopper, Ringing plate, Induction kick). The geometry's capacitances and inductances are validated too (tests Z1–Z3). Agreed design (owner, 2026-10-05): RADAU5 only in circuit levels; linear components first (DC, AC and pulse sources, R, L, C, timed switches); template circuits attached to plates and coils rather than a schematic editor. First stage planned with one-way coupling (the circuit drives the plates' potentials and the coils' currents; the particles' back-action, absorbed charge and induced currents, left out as a marked approximation with its bound), full coupling with the vacuum tubes.
  - First-principles devices (drift–diffusion semiconductors, resistive bodies with J = σE) stay far-future (§10, Stage 7).
- **Order of the next work** (owner, 2026-09-27): realistic iterations of the instrument levels, then the Jackson series, then circuits, then the inside of the components.
- **Inside the components: electrons (planned, in this order).**
  1. *Vacuum tubes* (exact classical electrodynamics): thermionic emission from a cathode, space-charge-limited current (validated against the Child–Langmuir law), the vacuum diode (rectification), the triode (a grid controls the current: amplification), later magnetron and klystron. Needs beams interacting near metal (the charges each electron induces on the electrodes) and continuous emission (particles launched over time).
  2. *Metals and resistors*, the Drude picture: classical electrons scattering off the lattice with a mean free time; Ohm's law and Joule heating emerge. Stated as a model: real conduction electrons are a quantum Fermi gas. Random collisions (fixed seed) make outcomes statistical.
  - *Mean field for many particles (planned, with the vacuum tubes).* Pairwise interaction (PHYSICS.md §3.3) is cheaper and more accurate up to about 10³ particles; beyond that a particle-in-cell mean field. Close-collision physics is then lost, which is acceptable only when each particle feels many others (a large plasma parameter), and must be stated. Steps:
     1. One electrostatic class (slow electrons in a tube). The particles' charge is deposited on a grid. Their field is the 3D Green's function convolved over sources confined to the plane (FFT), not a 2D Poisson solution: the world is a slice of 3D.
     2. Classes by velocity for relativistic clouds with a spread of speeds. "Relative speed non-relativistic" is not transitive, so the classes are fixed bins in rapidity (4-velocity), narrow enough (Δη ≪ 1) for slow motion within a bin. Each bin's field is solved in its rest frame and boosted (energy-binned quasi-static space charge, as in IMPACT-T, OPAL and ASTRA).
     3. Retardation between classes only where accelerations make it matter: a reduced radiation model (e.g. coherent synchrotron radiation along a bend). A full retarded mean field would need each class's charge history on the grid, integrated over the past light cone, which is costly.
     - *Validation:* against the exact pairwise retarded model on small N; the difference must shrink as the grid is refined and the bins narrowed.
     - *Worth it:* for the vacuum tubes (step 1, thousands of electrons) and for dense relativistic beams with an energy spread (bunch compression). Not for the game's current beams (≤ 32 particles).
  3. *Semiconductor devices*, ensemble Monte Carlo: electrons and holes as classical particles with effective masses, scattering at quantum-derived rates, in the field of a Poisson solution for the doped crystal (p–n junction, depletion layer, transistor action). Semiclassical: the band structure and scattering come from quantum theory, stated in the game. Far future.
- **Jackson series (planned).** Levels showcasing every effect and instrument of J. D. Jackson, *Classical Electrodynamics* (3rd ed.) that the simulation and hardware can compute. Each is named distinctively (e.g. "Jackson §12.3: E×B drift", "Jackson Pr. 2.4: The golden ratio") and cites its section or problem in the description (planned level field `reference`). The catalogue, with section and problem numbers checked against the book, the features each item needs, cross-topic highlights from the problems, and the items that wait for 3D, is `docs/JACKSON.md`.
- **Stated limits:** the game lists, per level, what is exact and what is approximated or left out ("Physics model and its limits", `Level::model_notes`), and shows the measured size of each approximation in the flight details (neglected radiation, electrode image force, radiation-reaction ratio).
- **Detector conditions:** a detector can also require the arrival direction (cone) and kinetic energy (window). They are verified with margins like every boundary (PHYSICS.md §6.1).
- **Radiation goals (implemented):** a detector can require the flight's radiation into a far receiver: an arc of in-plane directions, all frequencies or a frequency band, an energy per steradian in a window (a minimum, or a maximum for "quiet" goals). The game draws the receiver outside the arena's edge and shows the measured value and the spectrum over the receiver's directions with the band shaded (PHYSICS.md §3.4; the owner's choice, 2026-09-28: direction and frequency band). With other particles flying (free particles, free charges) the receiver sees them all and their far fields add, coherently (Pr. 14.23, a ring of charges); a beam shot has no radiation goal (its particles arrive one by one).
- **Hardcore mode:** continuous values (sliders, linear or log) within the ranges of the level's lists instead of the discrete values, and any antenna orientation. A level can also be hardcore by design (`limits.continuous`).
- **Moving elements:** mouse drag, or keyboard grab (G) and arrows.
- **Sandbox mode** (level editor): place level charges freely; set the launch point, direction and energy, the particle, c, the detector box, the grid and the player limits. Check solvability with the solver, store a reference solution, and save levels as JSON (`levels/custom/`) alongside generated ones. Resource meters show the measured cost of the setup (preview and verification time, steps, metal setup time, matrix memory, metal model error) against budgets, so levels stay interactive and inside the tested range of the solvers; the level tests enforce the machine-independent budgets (steps, memory) on shipped levels.

## 4. 2D and 3D modes

- One physics core, two presentation and control layers. 2D mode is its own experience: its own look, levels, generator settings and simpler controls.
- **Decision:** physics is always truly 3D, including in 2D mode. 2D mode is a cross-section of the 3D world. All charges lie in a plane of symmetry, and particles start in it with velocity parallel to it, so they never leave it. The field is the real 3D field (`1/r²`). "Flat" physics with a `1/r` field is excluded.
- In 2D, field lines show direction only. They are drawn evenly spaced (Jobard–Lefer style), start and end only on charges or the arena edge, and carry arrowheads along **E**. Their density in the plane does not represent field strength for a 3D field, and the UI says so. Strength is shown by the potential map and the arrows.
- 3D mode: editing by layers (the current grid slice is active, the others are semi-transparent).

## 5. Controls

- Mouse: place, remove and edit charges in 2D. In 3D, a layer mode with the scroll wheel changing layers.
- Keyboard: the cursor jumps between grid nodes. Arrow keys move it within a layer, PageUp/PageDown change layer (3D only). Keys to place, remove, flip sign and change magnitude.
- Later: group selection, copy, paste at the cursor, mirroring a group about a chosen plane, repeat counts (e.g. "paste 5 times every 2 cells"), configurable key bindings.

## 6. Visualisation

- Field lines start on a small sphere around each charge, with the number of lines proportional to the charge (3D). Traced on the CPU in f64 (cheap and exact enough).
- A vector field (arrows) on a chosen slice. A potential map in 2D (a fragment shader evaluating the potential per pixel, antialiased contours, forbidden region shaded: points no particle shown can reach, `U > E` for every shot, where a beam counts with its highest total energy `E = T + U(x₀)` over its particles; with one shot selected ("show all" off) only that shot's particles) and equipotential surfaces in 3D, computed on the GPU in f32 (visual only, never used for gameplay).
- Learning aids: slow motion, the force vector on the particle, a kinetic/potential energy bar during flight (relativistic kinetic energy `(γ−1)mc²`), the active physics model indicator, and the verified/marginal status of the result.

## 7. Level generator

1. Randomises the level's charges, A, B and launch parameters.
2. Searches the grid for a placement of k player charges for which the particle reaches B (simulated annealing or beam search). The objective is shaped to be continuous even when a trajectory is lost: distance to B, plus a penalty for how early the particle was lost.
3. Rejects trivial levels: without player charges the particle misses, and solutions with k−1 charges are not found in N independent searches.
4. Requires several distinct solutions (repeated searches from different starting points), so the level can also be solved in ways the generator didn't anticipate.
5. Robustness: moving any single charge by one cell must not always destroy the solution. The fraction of one-cell perturbations that survive is measured.
6. Computes a difficulty measure from points 3–5. Numeric thresholds are part of the acceptance criteria. Implemented as `generator analyze`: configuration-space size, random solve rate, the success and effort of a local search that sees only the distance objective, and smoothness of the objective across one-move neighbours (`docs/LEVELS.md`).
7. All reference solutions must be verified, not marginal (§2.3).
8. From Stage 2: verifies the level by simulating the full wave sequence (expensive, so done precisely only for the best candidates).
9. A native command-line tool, parallel over candidates. Saves levels as JSON: level charges, A, B, particle species, launch parameters, recommended grid, limits, reference solution, metrics, integrator tolerances and **physics engine version**.

## 8. Materials (future stages, accounted for in the architecture now)

- Metals: equipotential. Relaxation takes about `ε₀/σ ≈ 1e-19 s`, which is instantaneous here. They carry induced charge. Their electrostatic (instantaneous) response is a stated approximation: waves are not reflected or scattered by metal, and there are no waveguide or cavity modes. A full-wave solver (e.g. FDTD on a grid, possibly on the GPU) would lift it and unlock Jackson Ch. 8 (waveguides, resonant cavities) and scattering (Ch. 10): a possible future stage. In vacuum the fields are already exact (retarded Green's-function solutions of Maxwell's equations: Liénard–Wiechert, antennas, plane waves).
- Dielectrics: `∇·(ε∇φ) = −ρ`. Charged by stopped particles, discharged during pauses through leakage.
- **Optics: wave optics (planned, owner 2026-10-05; after dielectrics).** Mirrors, prisms, lenses, apertures and gratings, shown as waves: classical electrodynamics at the scale of the wavelength, neither geometric nor quantum, so that what geometric optics drops is visible and playable: diffraction-limited foci, interference, evanescent waves at total internal reflection (frustrated TIR, the Goos–Hänchen shift), the Gouy phase, the Arago spot, Talbot self-imaging, Brewster's angle, phase and group velocity in dispersive media, chromatic aberration, waveguide modes.
  - *Elements.* Mirrors (metal: a perfect conductor, or Drude) and dielectric bodies (lenses, prisms, slabs) with a dispersive `ε(ω)` (Sellmeier or Lorentz, so that a prism separates colours), apertures and gratings. Shapes are free curves (splines), customisable by the player: a lens can be morphed continuously into a prism, and a mirror can take any profile (flat, spherical, parabolic, freeform).
  - *Goals.* Light of a band into a detector region with a required intensity or spot size; designing a lens or mirror (the goal is the focus); separating colours; coupling into a waveguide.
  - *Model.* In the 2D slice optical elements are taken as uniform along z (cylindrical lenses, troughs), where the 2D Maxwell equations for each polarisation (`E_z`, or `H_z`) are exact; a stated difference from the game's other bodies, which are 3D.
  - *Method* (chosen by accuracy, then cost): a frequency-domain boundary integral solver for piecewise-homogeneous bodies (Helmholtz kernel, Hankel functions; Müller or PMCHWT formulation), exact up to the boundary discretisation and an extension of the electrodes' BEM; validated against Mie series for cylinders, the Fresnel coefficients and Sommerfeld's half-plane diffraction. Pulses and animation: Fourier synthesis over a band, or FDTD (general, with its numerical dispersion measured).
  - *Coupling with particles.* A charge crossing a dielectric radiates (Cherenkov and transition radiation, Jackson §13.4, §13.7; Smith–Purcell radiation from gratings), which also opens those catalogue items (`docs/JACKSON.md`) and many of Brau's.
- Semiconductors: drift-diffusion (Poisson with continuity equations for electrons and holes, mobility, diffusion, recombination). The most expensive stage, possibly never.
- Superconductors: London equations and the Meissner effect. Only meaningful once magnetic fields exist.
- Test particles fly in vacuum or channels. Materials respond as the environment.
- Methods: boundary element method for metals and dielectrics, and a 3D multigrid Poisson solver for semiconductors. Established native libraries (e.g. PETSc, hypre) are allowed, since the game is native.

## 9. Architecture and technology

- **Native desktop application.** Windows and Linux required, macOS desirable. No web version.
- **Rust** throughout.
- `physics` crate: pure library with no graphics dependency. It contains the field solvers, integrator, events and diagnostics. The game and the generator use the same code, so a generator solution always works in the game. Deterministic.
- A `FieldSolver` interface separates field sources from particle motion. Implementations: Coulomb; test-only analytic fields (e.g. uniform E); later Darwin, Liénard–Wiechert, BEM, grid Poisson, FDTD.
- Particles are stored as arrays (structure of arrays), not objects, so waves and streams scale.
- Libraries preferred over custom code wherever a mature, tested one exists:
  - Game and rendering: **Bevy** (pinned to one version for the whole project; built on wgpu, so Vulkan, DX12 or Metal) with **bevy_egui** for editor panels.
  - Integrator (decided in M1, see PHYSICS.md §5.1): our own step-at-a-time port of Hairer's DOP853. The coefficients are generated from `dop853.f`, and the result is cross-checked against `ode_solvers`. No existing crate exposed the per-step interpolant without heavy dependencies.
  - N-body (Stage 2): IAS15 (Rein & Spiegel 2015, REBOUND) and individual/block time steps are evaluated against the global-step DOP853 approach.
  - Parallelism: `rayon`. Serialisation: `serde`/`serde_json`. Vector math: `glam` (`DVec3`) or `nalgebra`. Float comparison in tests: `approx`.
- The live preview runs off the main thread. Rendering never waits for physics.

## 10. Stages and acceptance criteria

1. **Stage 1: vacuum, single particle.** 3D core (Coulomb, DOP853, events, outcome verification), 2D mode as a cross-section, editor (mouse and basic keyboard), live trajectory, field lines and potential map, generator points 1–7, saving and loading levels, `PHYSICS.md`. Criteria: all physics tests pass; the trajectory recomputes within one frame (16 ms) with 50 charges at preview tolerance; the generator produces levels meeting points 3–5 and 7.
2. **Stage 2: waves.** Coulomb interaction between particles, particle–particle collisions, the wave sequence with pauses, score, and a generator that verifies waves.
3. **Stage 3: full 3D mode.** Layer editing, 3D visualisation, group editing tools.
4. **Stage 4: static metals and dielectrics.** Started: metal spheres (grounded, fixed potential, floating) with image forces, PHYSICS.md §2.6, the metal levels of arc 2. Box electrodes (plates, slabs, walls) by BEM with measured error (§2.7), the electrode levels of arc 2, including player-placed plates and power supplies for tunable electrodes (Power supply, Tune the lens, Build a deflector, Shielding). Next: electrodes and spheres together, dielectrics.
5. **Stage 5: materials charged by particles, relaxation during pauses.**
6. **Stage 6: magnetic fields, Darwin approximation, superconductors.**
7. **Stage 7: Liénard–Wiechert with radiation reaction.** Optionally FDTD/PIC, semiconductors. (Radiation reaction, antennas and plane waves already exist for single particles; see PHYSICS.md §2.3–2.5, §3.1.)
8. **Stage 8: part design (§13).** Part editor and library, multipole abstraction with certified error, nesting, scaling, then responsive parts (T-matrix) once materials exist.

## 11. Physics tests (mandatory from Stage 1)

Each test compares against an analytic result, with thresholds stated in `PHYSICS.md`.
- Energy conservation in a static field: `γmc² + qφ` relative drift below a threshold.
- Rutherford scattering: deflection angle versus the analytic formula (non-relativistic limit), and versus the exact relativistic Coulomb scattering angle.
- Kepler orbit around an opposite charge: period, closure of the orbit, angular momentum conservation (non-relativistic limit).
- Relativistic Coulomb orbit: perihelion precession versus the analytic (Sommerfeld) result.
- Hyperbolic motion in a uniform field (test field solver): analytic `x(t)`, `v(t)`.
- Speed limit: `|v| < c` for extreme energies and fields.
- Symmetry: a particle starting in a plane of symmetry never leaves it.
- Event location: a particle aimed to graze a sphere at a known distance is classified correctly, just inside and just outside, down to the error tolerance.
- Convergence: the error scales with tolerance as expected for the method's order.
- Determinism: the same level gives bit-identical results across runs and across Windows and Linux builds.

## 12. Open decisions

- The physical scale of the reference world (cell size, reference particle species, typical energies). Non-relativistic electrostatics is scale-invariant, so this matters from Stage 2 (particle charge scale) and for materials.
- The final scoring formula (wave-based, speed-based, or a mix; see §3), after playtesting. The architecture must support both: launch speed is a per-run parameter, not baked into the level.
- Part design (§13): whether particles may pass through parts (opaque, transparent, or per-part apertures); how drive parameters and scaling appear in the UI; how the part library interacts with level progression (unlocks, as in Turing Complete).

## 13. Part design (future stage, after materials)

Inspired by Turing Complete's "build a component, then use it as a component". Here the component is a physical device, and what abstracts it is its field.

### 13.1 What a part is

- **Definition.** A part is a small design made in a part editor (a sandbox at its own scale). It holds elements, coils, antennas, materials (metals, dielectrics, magnetic materials, superconductors) and other parts. It is saved in a part library.
- **Parameters.** A part declares drive parameters, such as an electrode voltage, a coil current or an RF amplitude and phase, and possibly geometric ones.
- **Linear parts.** For parts made of linear materials, the field is linear in the drives: `F = Σ_k d_k F_k`. So a part is stored as a few basis fields, one per drive, and an instance is just a set of drive values.
- **Instances.** An instance has a position, an in-plane rotation, a mirror flag, a size scale `s` and its drive values.
- **Scaling.** Scaling is exact for linear static parts:
  - fixed charges give `E ∝ Q/s²`;
  - electrodes at fixed voltage give `E ∝ V/s`;
  - coils at fixed current give `B ∝ I/s`.

  So "the same device, twice as big" needs no re-solve.
- **The 2D slice.** It must stay exact. A part used in 2D must have the plane as its symmetry plane: sources symmetric about z = 0, magnetic moments and currents arranged as in PHYSICS.md §2.2. Rotations about z and mirrors in the plane preserve this, and the editor checks it.
- **Nesting.** Parts can contain parts (stacks of technology: electrode → lens → column → instrument). A level can require building a part, and later levels can offer it from the player's library, as in Turing Complete's progression.

### 13.2 What a part does to the world (minimal and extended)

1. **Minimum: the external field.** Outside the part's bounding sphere, its field is represented by a multipole expansion: solid harmonics for electrostatic and magnetostatic parts, vector spherical harmonics for time-harmonic, radiating parts. The expansion converges outside the circumscribing sphere, and its truncation error at distance `r` falls as `(R/r)^(p+1)` with a computable bound. The order `p` is chosen per evaluation to meet the tolerance. So the abstraction's error is **certified**, like everything else (§2.3), and verification (a tighter tolerance) automatically uses a higher order.
2. **Particles through a part: open question, with a natural default.** A part is a sub-scene, so inside its bounding sphere its exact field can always be evaluated from its children. Options:
   - **Opaque part:** its body is an obstacle, and only the external field matters. This is cheapest and needs only (1).
   - **Transparent part:** particles fly through. Near and inside the part the field is evaluated from the children (recursively), and far away from the expansion. This is exactly the near/far split of the fast multipole method, so a design's cost stays close to its number of visible parts, not the size of its contents.

   A per-part choice (opaque, or with defined apertures) is likely best. It is decided when the first real parts exist.
3. **Optional, later: transfer maps.** For parts built as beam optics (lenses, deflectors, analyzers), a design aid can show the part's transfer map from an entry plane to an exit plane: first-order matrices, or higher-order maps as in differential-algebra beam codes such as COSY INFINITY. This is for design hints only. Flights are always tracked exactly through the field.

### 13.3 Materials inside parts, and parts responding to each other

- **Linear response.** A part with conductors or dielectrics responds to external fields: induced charge, polarization. Around the part this response is linear. It is described by a response operator that maps the multipoles of the incoming field to the multipoles of the induced field (the T-matrix of multiple-scattering theory, Waterman). This operator is precomputed once per part definition, e.g. with the BEM solver of §8.
- **Scene solve.** A scene of several responsive parts then solves a small linear system for their mutual induction: parts × multipole coefficients. It is iterated to a tolerance instead of re-solving the whole geometry.
- **Nonlinear materials.** Ferromagnets with saturation and some superconductor regimes break superposition. Parts containing them fall back to a full solve of the enclosing scene, or are restricted to their linear range. This is stated per material.
- **Time-dependent parts.** Antennas and RF structures use frequency-domain expansions (exact outside the bounding sphere, including the radiation zone). Quasi-static parts use static expansions with time-dependent drives.

### 13.4 Performance and accuracy plan

- **Precompute per part definition:** its basis fields, i.e. multipole moments to a high order per drive, plus the response operator if the part is responsive. These are cached in the library with a hash of the definition.
- **Composition is cheap and exact.** Expansions translate and combine with the fast multipole method's translation operators (multipole-to-multipole, multipole-to-local). These are exact at the truncation order. Rotations about z and mirrors act exactly on the coefficients. So parts of parts cost only their own coefficients.
- **Evaluation** walks the part tree with an opening criterion, near → children, far → expansion. The order and the criterion come from the error bound, so every field value carries a certified bound that feeds the verification margin (§2.3). Evaluation order is fixed, so results stay deterministic.
- **Libraries.** Our source counts are small (tens of parts, a few levels of nesting), so an own implementation of solid harmonics and translations is simple and testable. Mature FMM libraries (e.g. exafmm, PVFMM) are references and fallbacks if scenes grow.
- **Tests,** mandatory when the stage starts:
  - The expansion converges to the direct sum at the predicted rate `(R/r)^(p+1)`.
  - Translation and rotation operators match a direct re-expansion.
  - A nested part matches the flat scene.
  - Scaling laws hold.
  - The response operator matches a full BEM solve of two nearby parts.
  - A particle through a transparent part matches the flat scene.
