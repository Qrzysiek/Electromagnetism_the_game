# The carriers' own magnetic field (Phase 3.1): design

The owner, 2026-10-26: "I would like to see how the whole circuit works at once and what
field it emits externally", and, asked how to bring circuits onto the board, chose this
stage first ("Full Phase 3.1 first"): every current's magnetic field acting on everything,
so that inductors, transformers and the skin effect come out of geometry, before the circuit
levels are rebuilt with real parts (docs/DRUDE.md, Phase 3; docs/DECISIONS.md, DEC-13).

This document is the plan; PHYSICS.md gets the model, method and measured tests as they
land.

## 1. What changes, physically

The tube world is z-invariant: charges are lines along z moving in the plane, and today
they feel the instantaneous electrostatic field (electrodes by BEM, particles directly) and,
optionally, a uniform B_z from the player's solenoid with its induced field. Currents in the
plane make no field of their own.

In the Darwin (magnetoquasistatic) limit of Maxwell's equations, valid while the system is
small against a wavelength (L ≪ c/ω) and the speeds small against c, a current in the plane
makes a magnetic field along z, and its change an induced (transverse) electric field:

- Ampère: ∇ × B = μ₀ J_T, so for currents in the plane B = B_z ẑ with B_z = μ₀ ψ, ψ the
  stream function of the current's divergence-free part (J_T = ∇ × ψẑ).
- Faraday: E_T = −∂A_T/∂t, ∇ · A_T = 0, B_z = ∂_x A_y − ∂_y A_x.
- Force: F = q (E_L + E_T + v × B), B the self-field plus the player's solenoid.

A metal ring is then a one-turn solenoid along z: its current makes a uniform field inside
and none outside, its inductance per unit length is μ₀ times its hole's area. Nested rings
share flux (a transformer), a slab in a swinging field carries its currents within a skin
depth δ = √(2/(μ₀σω)) of its faces, a ring keeps a swinging field out of its hole
(shielding), and a strip's current loop through its leads has the inductance of the loop.

Units: the game's Coulomb constant is 1 (ε₀ = 1/4π), so μ₀ = 4π/c². The metals'
collisionless skin depth is c/ω_p (ω_p = √(4πnq²/m)); for the levels' metal (ω_p ≈ 5) a
board-sized effect needs c of order 10–50, a "slow light" world like the 3D levels' (the
`SLOW_LIGHT` remark). The speed of light is a level setting (`TubeSpec::speed_of_light`);
without it (the default) nothing changes: every existing level and golden hash stays.

## 2. Method

### 2.1 The canonical momentum (magnetoinductive) formulation

Explicit Darwin (compute E_T from the last steps' ∂J/∂t) is unstable once a metal is larger
than its skin depth c/ω_p: a current's change induces a field that changes the current more
than it started with (loop gain (L/λ)², λ = c/ω_p), the regime where the effects are
visible. The standard cure solves for the transverse field implicitly with the medium's
inertial response (Nielson and Lewis 1976; Busnardo-Neto et al. 1977; Hewett and Nielson
1978; Gibbons and Hewett 1995). In the canonical-momentum form each particle carries
P = m v + q A_T(x); from the Darwin Lagrangian

- dP/dt = q (E_L + ∇(v · A_T)) + (the solenoid's force, as now),
- v = (P − q A_T(x)) / m,

and A_T is fixed at each instant by the particles' canonical momenta: with β(x) = Σ q²/m
per area (the carriers' density times q²/m: ω_p²/4π inside a metal, 0 outside),

  (−∇² + μ₀ β) A_T = μ₀ [J_P]_T,   J_P = Σ q P/m (per area),

projected on divergence-free fields. No time derivative of J appears: the induced field
comes out of P's evolution and A's dependence on it, implicitly and stably. In the plane,
A_T = ∇ × χẑ (B_z = −∇²χ) turns it into one scalar unknown, χ.

The operator is symmetric positive definite (it is the second variation of the magnetic
energy plus the carriers' kinetic energy), so a conjugate-gradient solve on a grid
converges, warm-started from the last step's χ.

Collisions (lattice.rs) redraw the kinetic velocity v; P follows as m v + q A_T(x). A
contact gives off carriers with their thermal v and P = m v + q A_T. Beams and the 3D
levels are untouched.

### 2.2 The grid

A uniform grid over the board (as `tube_em`'s, 4 per cell: ~100 × 80 nodes), charges'
P deposited by area weighting (CIC), χ solved there, A_T and B_z interpolated back at the
particles with the same weights (momentum-conserving pairing). The electrostatic field stays
as it is (BEM and the particles' direct field).

### 2.3 Boundaries, electrodes, circuits

- The board's frame: the world outside the board is empty, so the field must fall off as
  in free space. First version: the grid padded by a margin with χ = 0 on its edge
  (a distant perfectly conducting box: flux it holds is conserved, zero); the margin's
  effect on the tests is measured and bounded. Later, if needed, an open boundary
  (a boundary-integral closure or a FFT free-space Green's function on a doubled grid).
- Electrodes: perfect conductors, today equipotentials. With an induced field their
  surfaces' tangential E must vanish too, so their potential varies along them by the
  induced EMF. First version: electrodes are transparent to E_T (they keep the
  electrostatic condition alone); the error is the EMF along an electrode, ~ E_T × its
  size, measured in the tests and stated in the model notes.
- Circuits: a circuit's current runs along its lead (`tube_em`'s convention: from the
  enclosing ring or the board's edge), so the leads become sources of the field like any
  current, and the induced field's EMF along the lead and through the circuit enters the
  circuit's branch equation (the loop's self-inductance from geometry). The leads' paths
  move from `tube_em` into the level (one function both use), drawn in every view.

### 2.4 Cost

One CG solve per step on ~8000 nodes (warm-started: a few iterations), deposits and
interpolations linear in the particles: well below the particles' direct field. Measured
on the levels as they come.

## 3. Validation (written before the measurements, thresholds from the target accuracy)

Field solver alone (prescribed currents, no particles):

- V-D1 a thin ring's field: a current I per length around a square loop makes
  B_z = μ₀ I inside, 0 outside (1 % at the grid's resolution, away from the wire).
- V-D2 inductance: the loop's magnetic energy ½ L I² with L = μ₀ × area (per length).
- V-D3 mutual inductance of nested loops: M = μ₀ × the inner area.
- V-D4 the screened equation: a disc of carriers (uniform β) in a round loop's field
  screens it as London's B₀ I₀(r/λ)/I₀(a/λ), λ = 1/√(μ₀β) (the collisionless skin depth).
  (A slab that ends was tried first and failed: its screening current had nowhere to close,
  nothing in the field solver alone conserving charge; a round disc closes it by symmetry.)

With carriers:

- V-D5 a metal ring's current after a step of EMF rises as I = (V/R)(1 − e^{−tR/L}), L the
  geometric inductance, R the ring's resistance from Drude (2 % on L/R). (Measured so first;
  since 2026-10-26 as the impedance R + iωL by lock-in on a thin ring, the rise's integral
  being dominated by the ring's thermal current: PHYSICS.md §2.12.)
- V-D6 the skin effect: a disc in the swinging solenoid, its field's profile against a
  cylinder's I₀(kr)/I₀(kR), k² = iωμ₀σ(ω) with the model's AC conductivity (0.05 of the
  drive). (A slab was planned and tried first: at ω = 1 this metal is too poor a conductor
  for the slab formula, the charges closing its eddy loops making fields of their own;
  PHYSICS.md §2.12, DEC-13.)
- V-D7 nested rings as a transformer: the secondary's EMF against M dI₁/dt.
- Conservation: the electrostatic, kinetic and magnetic energy (½∫B²/μ₀) plus the
  lattices' heat and the circuits' work, constant to the method's accuracy; the total
  canonical momentum's balance.
- Symmetries (CLAUDE.md): mirror (B_z pseudoscalar: odd under a mirror in the plane),
  time reversal (B odd, the run reversed with v → −v and B → −B returning), rotation and
  translation of a setup, charge conjugation (all charges and B flipped), and the
  scaling c → ∞ recovering the electrostatic model (the self-field's effects ∝ 1/c²).

## 4. Order of work

1. `physics::darwin2d`: the grid, deposit, CG solve, interpolation; V-D1–V-D4.
2. The tube engine: canonical momentum for particles when `speed_of_light` is set, the
   self-field's force, collisions and emission; the level field and its editor control;
   V-D5–V-D7, energy, symmetries; PHYSICS.md §2.12.
3. Ring-shaped metals (Phase 2.3, needed for rings and transformers) unless square frames
   built from a rectangle with a hole are enough for the tests.
4. Leads as level geometry, drawn in every view; circuits' EMF from the induced field.
5. The display: the physics' own B_z on the magnetic map (today's map is `tube_em`'s FDTD
   field of the run, one-way).
6. Levels: Inductor, Transformer, Skin effect, Shielding, Induction heating; then the
   circuit levels rebuilt with real parts, the strip's properties as the player's controls
   (the owner: "consider also using properties of the strip as a parameter").

## 5. Progress

- 2026-10-26, step 1 done: `physics::darwin2d` (grid, sine-transform preconditioner,
  CG); V-D1–V-D3 exact, V-D4 (a disc's London screening) 0.61 %.
- 2026-10-26, step 2 (engine) in: `Tube::darwin`, canonical velocities in `TubeState::u`,
  the solve after the push (`Tube::darwin_step`), collisions and walls changing `u` by
  their Δv; `Lattice::with_hole` (rings). V-D5 (a ring's current rising at L/R) passes over
  three seeds (2.095 against 2.125); `crates/physics/tests/darwin.rs`, `--ignored`.
- 2026-10-26, V-D5 redone as the thin ring's impedance by lock-in: L/R from the phases within
  0.3 % of the field's energy over the drive's work (the rise time's integral could not decide
  between the kick's orders over nine seeds).
- 2026-10-26, V-D6 passes on a disc (round metals, `Lattice::with_disc`): the field's
  profile within 0.041 of a cylinder's skin effect (required 0.05; 0.027 since the kick follows
  the solve). The slab tried first
  failed, most likely for its charges rather than the coupling (PHYSICS.md §2.12).
- 2026-10-26, symmetries and energy Y21–Y25 pass (PHYSICS.md § Symmetries). They found the
  backgrounds counted twice below 2000 particles (a metals bug, fixed) and the kick taken from
  the last step's field (a self-drag, fixed: the kick now follows the solve). Known limits:
  first order in time (levels with the self-field: dt ≲ 0.015), and impulses (collisions,
  walls) lose the carrier's electromagnetic mass, so heat is not measurable under the
  self-field yet.
- 2026-10-27, impulses carry the carrier's own inductance (DEC-24): `Δu = Δv + q(q/m) g Δv`,
  `g` the solver's self-response at the carrier from two solves per metal (`SelfResponse`).
  Y25's walled disc 14.6 % → 4.6 % of U per time unit; the closed disc 3.5 % → 2.1 %, the rest
  fitting `α dt + β/N` (the step's first order now dominating at dt 0.02). V-D5 and V-D6 pass
  as before (L/R +3.9 %, R ratio 1.056; V-D6 0.025).
- 2026-10-27, V-D7 tried twice, not passed yet. (1) A slotted secondary round a closed
  primary, the slot's voltage against the primary's flux: the slot charges through the ring's
  resistance with `RC` ≈ 24, longer than the run (its voltage reached 5 of an EMF of ~30).
  (2) Two closed thin rings by lock-in, `G − G′ = Φ₁G′/(DS₂)` from the secondary's circuit law,
  `S₂` measured as `F′/G′` (357–363, the uniform current's 361): 15 % apart over two seeds on
  the engine before the impulses' correction, 25.5 % after (required 20 %): the seeds scatter
  by ±25 % (the runs with the primary noisier) and the measured change lags the prediction by
  6–10° in all four, as a ring that is not one lumped impedance would: at ω = 0.8 this metal is
  a poor conductor (`4πσ/ω` ≈ 3, its surface charges lagging the drive, as in V-D6's slab), and
  the primary's EMF and the solenoid's are spread differently along it. Next: ω ≈ 0.3
  (`4πσ/ω` ≈ 9) and more seeds (`target/scratch/vd7_rings_wip.rs`).
- 2026-10-27, a level's speed of light (`TubeSpec::speed_of_light`), its grid
  (`self_field_grid`), its issues (step ≤ 0.015, speed 0.5–10⁴, grid ≤ 511 a side) and the
  sandbox's control.
- 2026-10-27, step 5: the field drawn. Each frame samples the self-field's `B_z` on the board
  (0.5 apart); the magnetic map shows it for a level with a speed of light (`maps.rs`), its
  caption saying so. A strip with its supply off, the self-field on, booked −11.4 of Joule
  heat per time unit in the panel: the heat under the self-field is the next thing to measure.
- Still to do in step 2: V-D7 (above); the lattice's heat under the self-field against the Joule
  heat; perhaps an implicit-midpoint step; then the field drawn on the board (step 5), leads
  (task 43), electrodes' induced EMF.
