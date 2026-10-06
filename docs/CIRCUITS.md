# Circuits: plan of the stage

The design agreed with the owner on 2026-10-05 (`SPEC.md`, "Circuits"): lumped
components coupled to the field world; the stiff integrator RADAU5 only in circuit
levels; linear components first (DC, AC and pulse sources, R, L, C, timed switches);
*template circuits attached to plates and coils* rather than a schematic editor. This
note is the plan of the first stage, written before building it. Questions the owner may
want to settle are marked **(owner)**.

## Done

- **The stiff integrator.** A port of Hairer & Wanner's RADAU5 (PHYSICS.md §5.4, tests
  I1–I5): the same decisions as the Fortran code on three classical stiff problems,
  accuracy against high-precision references, and a series RLC circuit written as an
  index-1 DAE against its closed form.
- **The geometry's circuit parameters** (PHYSICS.md §2.8, tests Z1–Z3): the electrodes'
  capacitance matrix by the boundary element method, mutual inductances by Neumann's
  formula, a ring's self-inductance.
- **The circuit engine** (`circuit.rs`, PHYSICS.md §2.9, tests Z4–Z9): modified nodal
  analysis of R, L (with mutual inductances), C, capacitance blocks, DC, AC and pulse
  sources and timed switches, integrated by RADAU5 between breakpoints with consistent
  restarts; the inductors' current rates from their law, for the induced fields. Against
  closed forms (RC, LC, charge sharing through a switch with its heat, a pulse, coupled
  inductors, a driven tank and its energy balance): within 1e-9.
- **Time-dependent coils.** Coils already take a linear ramp of their strength with the
  exact induced field `−κ̇ A_unit` (circular and polygonal coils, test M10).

## Units

The game's units are SI in form with `k = 1/(4πε₀) = 1` and `μ₀/4π = 1/c²` (PHYSICS.md
§1, §2.2). Hence, in circuits:

- a capacitor `Q = C V`, with C in lengths (a sphere of radius R has C = R);
- a resistor `V = R I`; an inductor `V = L dI/dt`; RC and L/R are times;
- a coil is specified in field units, `κ = μ₀ I/4π = I/c²`; its self- and mutual
  inductances are `N/c²` with N the geometric Neumann integral (a length). Coils in
  circuits therefore need a finite c (for c = ∞ a finite field needs an infinite
  current); plates work at any c.

## Model of the first stage

- **Components:** R, C, L; voltage sources DC, AC `V₀ sin(ωt + φ)` and pulses
  (trapezoids: rise, flat top, fall); switches that close or open at set times.
- **Elements as circuit parts.** A plate (electrode) is a node: the Maxwell capacitance
  matrix of all electrodes (from the BEM) is a block of capacitors between them and to
  ground (the electrodes not in a circuit keep their bias: fixed potential, grounded, or
  floating with a fixed charge, which then follows the others through the capacitance
  matrix). A coil is an inductor with its self-inductance and the mutual inductances to
  the other driven coils.
- **Equations.** Modified nodal analysis: unknowns the node potentials, the inductor
  currents and the voltage sources' currents; `M y' = A y + b(t)` with M from the
  capacitances and inductances, an index-1 DAE (checked: no loop of capacitors and
  voltage sources, no cut set of inductors and current sources). Integrated by RADAU5
  with the exact (constant) Jacobian, piecewise between breakpoints (switching times,
  pulse corners), restarting at each.
- **Coupling, one way.** The circuit is solved first, independently of the particles; its
  dense output gives each driven plate's potential `V(t)` and each coil's current. The
  particles then fly (DOP853, as now) in the fields of these sources: a plate's surface
  charge `σ(t) = σ_static + Σⱼ Vⱼ(t) wⱼ` (the unit responses of the BEM, including the
  floating electrodes' response), a coil's field `κ(t) B_unit` and induced field
  `−κ̇(t) A_unit` (generalizing the linear ramp). The particle integrator stops at the
  circuit's breakpoints, where the fields' time derivatives jump.
- **Left out, marked and bounded (first stage):** the particles' back-action on the
  circuit, the charge a particle deposits on a plate and the currents it induces
  (Shockley–Ramo); bound `q/(C V)` per particle, checked small in the levels, which use
  test charges. The circuit's own radiation and retardation: quasi-static, valid while
  the setup is small against `c/ω` of the circuit's frequencies (the parameter shown in
  the model notes). Full coupling (one stiff system of particles and circuit) belongs to
  the vacuum tubes (`SPEC.md`, inside the components), where the space charge sets the
  current.
- **(owner)** Is the one-way coupling acceptable for the first circuit levels? Its error
  is the back-action above; the alternative, full coupling from the start, puts every
  flight of a circuit level into RADAU5 (an implicit solve with the field Jacobian per
  step), several times slower.

## Validation (planned tests)

- RC charging of a plate pair through a resistor (C from the BEM): `V(t) = V₀(1 − e^{−t/RC})`.
- A coil's current through R from a DC source (L from the geometry): `I(t) = (V₀/R)(1 − e^{−Rt/L})`, and its induced field against the ramp's.
- An LC circuit made of a plate pair and an inductor: the frequency `1/√(LC)` and the
  energy exchange; with R, the damping.
- Energy balance: the sources' work = Joule heat + `½ C V²` + `½ L I²`, to the
  integration's accuracy.
- A switch: the charge redistribution when it closes (two capacitors in parallel: the
  energy lost is the classic `½ C V²/2` into the resistor, independent of R).
- Particles: a charge between large plates whose voltage rises as RC: its deflection
  against the closed form; a charge at rest in a coil ramped by an RL circuit: the
  induced field `−κ̇ A` (Faraday's law around the coil) and the betatron condition.

## Game

- Level format: an optional `drive` on plates and coils, a template with parameters, e.g.
  `{"rc": {"source": {"dc": 1e5}, "r": 2.0}}`, `{"lc": {"l": 3.0}}`,
  `{"switch": {"source": ..., "r": ..., "at": 12.0}}`. The level editor gets the controls
  (the exhaustive destructuring of `level_editor.rs` enforces it); the player's plates and
  coils can take the templates the level offers, with parameters from lists.
- The panel shows each driven element's `V(t)` or `I(t)` over the flight, with the
  particle's times marked.
- Levels: a timed deflector (an RC delay sets when a plate is charged), an LC kicker
  (the particle must meet the right phase), a chopper (a switch passes one bunch), a
  betatron driven by an RL circuit (the induced field accelerates, the field bends; the
  betatron condition), each with an introduction first (`docs/CURRICULUM.md`).

## Later stages

Nonlinear devices (diodes: Shockley; transistors: Ebers–Moll or square-law), stated as
semiclassical device models with fitted parameters; current sources; the full coupling
with particles (absorbed charge into the plate's node, Ramo currents), which the vacuum
tubes need.
