# Design decisions

Every design choice made in building the game: the model, algorithm or design taken, its
purpose, what else was considered and why not, and where it lives. When a choice is changed,
its entry is not rewritten: the change is appended to its history with the reasons, so that
the log shows how the game got where it is (the owner, 2026-10-26: "make a note of every design
choice and every time you change it ... append this note with the same info and reasons for
change").

Scope: features, models and algorithms (physics, numerics, display, tools), and the main
choices behind each level's design (the full measurements stay in the designs' comments,
`scripts/designs/`, and `docs/LEVELS.md`). The physics' models and their validation are in `PHYSICS.md`,
whose sections also record their own earlier versions; plans in `docs/TUBES.md` and
`docs/DRUDE.md`. Choices made before this log began (2026-10-26) are recorded there and are
added here when they are next touched.

Entry format:

```markdown
## DEC-<n>: <short title>

- **Status:** in use | replaced by DEC-<m> | deferred | dropped (since <date>)
- **Area:** physics | numerics | display | interface | tools | level design | process
- **Choice:** <the model, algorithm or design, precisely>
- **Purpose:** <what it is for>
- **Alternatives:** <what else was considered, and why not>
- **Where:** <files, functions, tests, docs>
- **History:**
  - <date>: chosen. <the reasons, the measurements that decided it>
  - <date>: changed from <a> to <b>. <why: what was found inappropriate, the evidence>
```

Entries are numbered in order and never renumbered.

## DEC-1: Carriers taken and given off by a metal's contacts between display frames

- **Status:** in use (since 2026-10-26)
- **Area:** display
- **Choice:** between two stored frames (`TubeRun::frame_between`), a metal's carrier that a
  contact takes before the next frame goes into it, and one a contact gives off after the
  frame comes out of it. A carrier in the per-step sample (`TubeFrame::paths`, every traced
  one) follows its recorded steps up to its last step before being taken, or from its first
  step after being given off. Any other moves in a straight line to (or from) the nearest
  point of its metal's contact faces (`MetalView::contacts`, `Level::contact_faces`) at the
  carriers' mean speed `√(π/2)·σ` (a 2D Maxwellian's), arriving by the next frame at the
  latest. Particles outside metals keep the old rule (their frames are steps).
- **Purpose:** a smooth picture by the contacts. Before, a taken carrier stood still at its
  last frame's position until the next frame and then vanished, and the given-off ones
  appeared all at once with the next frame; where carriers come and go constantly (a
  contact's neutral layer) the dots looked "locked" and moved "in bursts" (the owner,
  2026-10-26, in the Inverter and the Amplifier, whose frames are 0.4 and 0.74 apart within
  `FRAME_BUDGET`; Gate's 0.16 hid it).
- **Alternatives:** (1) record every absorption and emission (time, place) in the frames:
  exact, but ~2000 carriers a time unit cross the Amplifier's contacts, 15–30 MB over a run;
  (2) store frames more often: the budget ×5 (~340 MB) for Gate's spacing; (3) fade the
  carriers in and out where they are: they would still stand still. The sampled carriers'
  paths already hold the exact answer for the traced ones (the ones with trails, which the
  eye follows); the others' straight-line guess is the same kind of approximation as their
  straight-line motion between frames.
- **Where:** `crates/level/src/tube/run.rs` (`frame_between`, `taken`, `given_off`,
  `nearest_contact`), `crates/level/src/tube/geometry.rs` (`contact_faces`); test `contacts_take_and_give_off_carriers_between_frames`.
- **History:**
  - 2026-10-26: chosen, replacing "stays where it was / appears with the next frame".

## DEC-2: The field driving a metal's carriers in the Scales (qEℓ/kT)

- **Status:** in use (since 2026-10-26)
- **Area:** interface
- **Choice:** a metal's driving field is the largest voltage between two of its contacts
  over the distance between their faces' middles, from the run's mean electrode potentials
  over the later half of its frames so far; for each material the largest over its metals.
  Before a run, the player's largest supply over the metal's longer side; a swinging
  solenoid's induced field at the edge as before.
- **Purpose:** the qEℓ/kT ratio (whether the carriers warm, Ohm's law holds) must use the
  voltage actually across the metal.
- **Alternatives:** the circuits' source levels (static, an upper bound: 120 for the
  Amplifier's drain, which sits near 45 behind its load); the latest frame's potentials
  (they flicker with the noise and the signal); the mean over all frames (biased by a
  supply's start-up, the Amplifier's drain settling from 120).
- **Where:** `crates/game/src/scales.rs` (`drive_field`), `crates/game/src/ui/results.rs` (the mean
  potentials).
- **History:**
  - before 2026-10-26: the player's largest supply over the metal's longer side. Wrong where
    a supply is fixed in a circuit (the Wheatstone bridge read 0), across a wide strip (its
    10 instead of the 4 between its contacts), and on the transistors, where the player's
    supply is the gate: Gate read 0.28 from its gate's 18 instead of 0.156 from its drain's
    10.
  - 2026-10-26: changed to the contacts' potentials from the run, for those reasons.

## DEC-3: The tube runs' particle caps

- **Status:** in use, unchanged (`physics::tube::MAX_PARTICLES` 20 000, the metals' share
  `MAX_METAL_CARRIERS` 15 000 in the verification)
- **Area:** numerics
- **Choice:** keep the caps. A run stops emitting at the cap (its verdict uncertain, never a
  crash); a level whose metals hold more in the verification is refused.
- **Purpose:** a guard against a degenerate macroparticle weight (one once asked for ~1e300
  particles) and a bound on memory and run time.
- **Alternatives:** doubling both (tried for the Flip-flop, DEC-4: it ran, but the cell
  latched too slowly to make a level, so there was no reason to keep it).
- **Where:** `crates/physics/src/tube/`, `crates/level/src/tube/`; PHYSICS.md §2.9, §2.11.
- **History:**
  - 2026-10-26: doubling considered and reverted (DEC-4).

## DEC-4: Flip-flop (the transistors' master)

- **Status:** deferred (since 2026-10-26)
- **Area:** level design
- **Choice:** two of Gate's transistors, each pair of gates wired to the other's drain
  through a level-shifting source −S (a branch ending on the other drain: `return_to`; a
  depletion gate needs ~−70 to switch off while a drain sits at 20–80), each drain fed from
  120 through 200; reset and set by dropping a drain's supply to 0 for ~40.
- **Purpose:** a bit of memory: a bistable cell holding its state after a pulse.
- **Alternatives:** 4-long channels (fit the particle cap but do not switch: at V_DD 60,
  0.73 off against 1.32 on); fewer carriers per area (the depletion is then an artefact, at
  200 per area Gate's current plateaued from −10); gates shorter than the channel by the
  source (less gate capacitance, a drift region shielding them) — not yet tried.
- **Where:** `target/scratch/ff_proto.py`; docs/DRUDE.md (Phase 1, item 3).
- **History:**
  - 2026-10-26: deferred. Two 8-long channels need 25 600 carriers in the verification
    (DEC-3); with the caps doubled the cell latched too slowly: each drain carries the other
    transistor's gates (~0.6, more through the Miller effect) through its load, a time
    constant of ~200–400, so after a reset at S = 70 the drains crept back toward the
    symmetric point (A 50, B 54) and a set separated them by only ~7 by t = 300. A level
    would need runs of ~600–1000 (previews of minutes).

## DEC-5: Amplifier (level 123): channel, signal and supply

- **Status:** in use (since 2026-10-26)
- **Area:** level design
- **Choice:** Gate's 8-long channel; a sine of amplitude 3 and period 300 on both gates
  around the player's bias; the drain fed from 120 through 200; the goal the drain's
  oscillation amplitude 4.5–7 (a gain of at least 1.5) over one period from t = 150.
- **Purpose:** voltage gain, and the bias point's role (pinched: no gain; strongly on: the
  channel a plain resistor).
- **Alternatives:** shorter channels (faster but less gain: L 4 gained 0.7); a period of 100
  or 200 (the drain's capacitance through the load, g_m/C ~ 0.06, left gains of 1.1–1.4);
  400 (gain ~2, runs ~50 % longer); a load of 300 (no more gain at period 300: the
  gain-bandwidth limit); V_DD 100 (the window −20 to −13, 6 %, against −20 to −5, 12.5 %, at
  120: a higher drain stays saturated to higher biases).
- **Where:** `scripts/designs/components.py` (`amplifier`, its comment holds the scans), docs/LEVELS.md.
- **History:**
  - 2026-10-26: chosen. The window's lower edge is set by the preview's unconverged
    pinch-off (the verdict's error is 1.5 × the runs' difference: at −22 the verification
    gains 1.61, the preview 0.92).

## DEC-6: Wheatstone bridge (level 116): a half bridge of two strips and two resistors

- **Status:** in use (since 2026-10-26)
- **Area:** level design
- **Choice:** Magnetoresistance's two strips on one fixed supply of 12, each far end to
  ground through an ideal resistor (5 and 1.5), a voltmeter between the far ends; the
  player's solenoid; the goal the voltmeter within ±0.3 over t = 10–150.
- **Purpose:** a magnetoresistor read by a bridge: a balance on a ratio, a field sensor.
- **Alternatives:** four strips (over the particle cap, DEC-3; two arms ideal resistors is
  also how real sensor half bridges are built); resistors 10 and 2.8 with a supply of 10 and
  t = 10–80 (the output moved ~0.17 per unit of B against error bars up to 0.19: noisy; the
  noise of a strip's end grows with the resistor behind it, so the resistors were halved,
  the supply raised and the window lengthened).
- **Where:** `scripts/designs/components.py` (`wheatstone_bridge`), docs/LEVELS.md.
- **History:**
  - 2026-10-26: chosen.

## DEC-7: Ringing tank (level 124, the AC and induction arc's first)

- **Status:** in use (since 2026-10-26)
- **Area:** level design
- **Choice:** a parallel tank: a strip's left contact fed by the player's step at t = 5
  through a coil L = 7.5 with no series resistance, a capacitor C = 1.5 from it to ground,
  the strip's right contact grounded; the goal the ringing's amplitude 4.5–9 at its observed
  frequency (2π/23) over three whole periods (t = 5–74); the step's slider ±50.
- **Purpose:** the strip as a circuit's resistor: its collisions set the tank's Q. An
  introduction (the slider scales the ringing in proportion); the arc's next level gives the
  player the Q itself (the solenoid on a magnetoresistive strip).
- **Alternatives:** a series resistance 0.5 in the coil's branch (it shared the damping: Q
  2.4 instead of ~6; removed so that the strip alone damps); the goal at the ideal 1/√(LC)
  (the strip and its contact add ~0.3 to C: the tank rings at period ~23, not 21); a
  non-whole number of periods (leaks the step into the Fourier amplitude); a ±40 slider (the
  window would have reached the guessable 18–22).
- **Where:** `scripts/designs/components.py` (`ringing_tank`), `target/scratch/rlc_proto.py`.
- **History:**
  - 2026-10-26: designed. The reference 14 verified 6.57 ± 0.29; the step solves on about
    10.5–18.5 (9.5–10 too, 10.25 uncertain at the goal's edge; an isolated 19.5) and −18 to
    −14 at least (−10 not), 8–9 % of ±50.

## DEC-8: The design-decision log itself

- **Status:** in use (since 2026-10-26)
- **Area:** process
- **Choice:** this file, append-only per entry; CLAUDE.md requires an entry with every
  choice and a history line with every change.
- **Purpose:** the owner wants to see each model and algorithm chosen, why, and why it was
  changed when it was.
- **Where:** `docs/DECISIONS.md`, CLAUDE.md (Workflow).
- **History:**
  - 2026-10-26: started, with this session's choices; older ones are added when touched.

## DEC-9: Field-tuned damping (level 125, the AC arc's first intermediate)

- **Status:** in use (since 2026-10-26)
- **Area:** level design
- **Choice:** Ringing tank's parallel tank (L 6, C 1.5, no series resistance) with
  Magnetoresistance's short wide strip as its resistor, in the player's solenoid; a fixed
  step of 10 at t = 5; the goal the ringing's amplitude 2.15–3 at its observed frequency
  (2π/20) over three periods (t = 5–65); the strip at 100 carriers per area (twice the
  levels' metal).
- **Purpose:** the player sets a tank's Q with a field: from shorted out (Q ~0.5) to ringing
  (~3), a magnetoresistor as a contactless variable resistor.
- **Alternatives:** the step 6 or 8 at 50 per area (error bars up to 0.48 against a slope of
  ~0.17 per unit of B: the strip's thermal noise rings the tank, the more the higher its Q;
  more carriers and a larger step halved it); a goal from 2.0 (its window would have reached
  the slider's guessable 7.2–8.8).
- **Where:** `scripts/designs/components.py` (`field_tuned_damping`), `target/scratch/ftd_proto.py`.
- **History:**
  - 2026-10-26: designed. The reference 11 verified 2.42 ± 0.04; B solves on 9.5–13
    (10.9 % of ±16) and −13 to −11 (−9 not).

## DEC-10: Lock-in Hall (level 126, the AC arc's second intermediate)

- **Status:** in use (since 2026-10-26)
- **Area:** level design
- **Choice:** Hall probe's strip (100 carriers per area, floating side contacts) with a
  steady current from the player's supply, in the player's solenoid swinging as
  B₀ sin(0.125 t); two goals on the top contact against the bottom one at the field's
  frequency over two whole periods from t = 10: its amplitude 4.9–7 and its phase −60° to
  30° (in step with the field). The reference (supply −30, B₀ 5).
- **Purpose:** lock-in detection: a signal's amplitude and phase at one frequency, out of
  the noise and the strip's other responses (the magnetoresistance at 2ω); the phase gives
  the current's (and the carriers') sign, so the supply's sign matters.
- **Alternatives:** an amplitude goal alone (any sign of the current would solve: the phase
  is the lock-in's other half and makes the sign a choice); goals of 3.5–6 or 3.2–6.3 (the
  voltage goes as I·B₀, so one ratio sets both sliders' windows: the supply's would have
  reached its end, 38–40, or its middle, 18–22, and the field's 7.2–8.8, where the Hall
  voltage turns sublinear); the supply in its lower half (a smaller signal against the
  same noise); a slower swing (twice the run time).
- **Where:** `scripts/designs/components.py` (`lock_in_hall`), `target/scratch/lockin_proto.py`.
- **History:**
  - 2026-10-26: designed. The phase lags the ideal by ~17° for either sign of the current:
    the contacts' capacitance charging through the strip (an inductive pickup would shift
    the two signs opposite ways); kept and explained, as in a real lock-in. The reference
    verified 5.81 ± 0.24 at −18.0° ± 2.5°.
  - 2026-10-26: changed from two periods and the goal 4.9–7 to three periods and 4.7–7: B₀
    had solved on 4.5–6.0 only (4.7 %, under the 5 % the windows need; the supply −36 to
    −26). Three periods shrink the error bars (~0.8 times) for runs 1.4 times longer; the
    lower bound gives the field's window its lower part while the supply's stays above the
    guessable 22. Now B₀ 4.0–6.5 (7.8 %), the supply −36 to −24 (15 %).

## DEC-11: Resistive anode (Phase 1.6, a combination of the tubes and the metals)

- **Status:** in design (prototype, 2026-10-26)
- **Area:** level design
- **Choice:** a temperature-limited cathode's beam lands on a metal strip whose two ends are
  grounded contacts; the beam's current divides between the ends as the strip's resistances
  either side of the landing point (charge division); the player steers the beam (a pair of
  deflection plates, one on the supply) to make the current divide in a given ratio.
- **Purpose:** the first level joining the tubes' beams and the metals' carriers, as a real
  device: the resistive-anode position-sensitive detector (imaging detectors read where a
  particle lands from how its charge divides).
- **Alternatives:** the plan's first idea, the strip's heat and IR drop where a beam lands
  (the heat goal measures the field's work on the carriers inside the metal, not the energy
  a beam brings in: it read −3.6, the wrong measure); steering with the player's solenoid
  (it also acts on the strip: Hall and magnetoresistance in the anode).
- **Where:** `target/scratch/anode_proto.py`.
- **History:**
  - 2026-10-26: prototype. The engine joins them as it is: a beam electron (q/m −1, the
    metal's carriers' too) inside a metal at a step's start scatters off its cores and walls
    like the metal's own carriers and leaves by a contact. A beam of 0.4 at the middle of a
    16 × 4 strip divided 0.217 / 0.200 (52 %). Steering: plates 4 apart took part of the
    beam (the strip collected ~0.2 of 0.4); 8 apart around a 40 V beam the left end's share
    went only 0.51 → 0.63 over the supply 0 → 40 (a shift of ~2 of the strip's 16): a slower
    beam or another geometry needed before a level.

## DEC-12: A switched-off supply drives nothing

- **Status:** in use (since 2026-10-26)
- **Area:** level design, tests
- **Choice:** a player-tunable electrode's circuit (or driven plate's chain) carries, in
  the level, a source at rest: a DC source at 0, a sine with no swing (no offset where the
  supply sets the offset), a pulse whose top is its rest; the player's supply sets its
  level. Switched on, a supply starts at 0, which the slider rule keeps out of every
  solving window. Test `switched_off_supplies_drive_nothing`.
- **Purpose:** a level must not solve, or half solve, with nothing done (the owner,
  2026-10-26: "the winning conditions are trivial as it solves with the initial setting").
- **Alternatives:** a source's level anywhere outside the solving window (an arbitrary
  number the "off" hover would still show, and a test could not tell from a giveaway).
- **Where:** `scripts/designs/components.py` (rectifier, klystron, inverter, ringing_tank),
  `crates/level/tests/levels.rs`.
- **History:**
  - 2026-10-26: chosen. Ringing tank's step (14), the Inverter's pulse top (−10) and the
    Rectifier's anodes' swing (30) were their references, the Klystron's buncher swing
    (11) inside its window: they solved, or half solved, with the supplies off. Now each
    fails with everything off (`target/scratch/check_off.py`).

## DEC-13: The carriers' own magnetic field: the Darwin limit, canonical momenta

- **Status:** in design (2026-10-26; docs/DARWIN.md)
- **Area:** physics, numerics
- **Choice:** the magnetoquasistatic (Darwin) limit in the z-invariant tube world: currents
  in the plane make B_z, its change the transverse E_T. Particles carry canonical momenta
  P = mv + qA_T; A_T (as a stream function χ, B_z = −∇²χ) solves
  (−∇² + μ₀β)A_T = μ₀[J_P]_T on a board grid each step by conjugate gradients, β the
  carriers' q²/m per area. Collisions redraw v (P follows). A finite speed of light per
  level, none by default (existing levels unchanged). Electrodes transparent to E_T at
  first; circuits' leads as current sources and their EMF in the branch equations.
- **Purpose:** the owner chose this stage first (2026-10-26) to see whole circuits and the
  fields they make: inductance from geometry, transformers, skin effect, shielding.
- **Alternatives:** explicit Darwin (E_T from the last steps' dJ/dt: unstable once a metal
  is larger than its skin depth c/ω_p, the regime that shows the effects); full
  electrodynamics (FDTD acting on the particles, Phase 3.3: a time step bound by c/h, far
  costlier, and radiation not needed at these sizes); lumped coils driving local solenoids
  (the owner preferred the full field); drawing the circuits as symbols only (no fields).
- **Where:** docs/DARWIN.md (plan), `physics::darwin2d` (to come).
- **History:**
  - 2026-10-26: designed (docs/DARWIN.md), tests written first (V-D1–V-D7, energy,
    symmetries).
  - 2026-10-26: the field solver (`physics::darwin2d`): V-D1–V-D3 exact to the solver's
    tolerance (a loop is exactly a curl on the staggered grid), V-D4 0.61 % in 9 CG
    iterations. V-D4 changed from a slab to a disc: in a slab that ends the screening
    current −βA had nowhere to close (the field at its faces came out 0.4 B₀), as nothing
    in the field solver alone conserves charge; the electrostatic field will close it in
    the full model.
  - 2026-10-26: in the tube engine (canonical velocities `u = v + (q/m)A`, the solve after
    the push, collisions changing `u` by their `Δv`), metals with holes (rings). V-D5, a
    ring's current rising at L/R: the measured rise (less the carriers' own inertia) 2.095
    against L/R 2.125 over three seeds, the final current unchanged. The measurement
    changed twice: from one arm's carriers' velocities (noisy, and `dt/τ` below the drift the
    field sees, the splitting's first-order lag) to the field in the hole, `B/μ₀`; and R from
    the external field's work instead of the Joule heat (10 % scatter against 1 %).
  - 2026-10-26: V-D6 (the skin effect) changed from a slab to a disc (round metals, DEC-14).
    The slab 40 × 6 at ω = 1 did not follow the slab formula even at weak coupling (its
    middle's field 1.03 of the drive's at c = 10 against 0.96): the formula assumes a good
    conductor (`4πσ/ω` ≫ 1; here 2.7) without charges, and in a long slab the charges that
    close the eddy loops, and the Coulomb displacement current the Darwin limit keeps, make
    fields of their own (PHYSICS.md §2.12). A disc needs no charge (its eddy current runs
    along the induced field) and has the exact cylinder solution: worst 0.041 against 0.05.
  - 2026-10-26: the canonical kick `(q/m) dt ∇(v·A)` moved after the solve, from the field
    solved at the same positions: taken from the last step's field (solved a move before),
    each particle's own vector potential dragged it (the energy test Y25: a closed
    collisionless disc lost 0.35 per time unit at dt 0.02, 0.77 on a grid twice as fine;
    0.18 and 0.21 after, the rest the walls' impulses below and the macroparticles' own
    field; V-D6's worst 0.041 became 0.027). Two limits found and kept for now: the scheme is first order (time reversal 6.1 %
    of the self-field's effect at dt 0.02, 3.1 % at 0.01: levels with the self-field step at
    most ~0.015 until an implicit-midpoint scheme); collisions and walls lose the carrier's
    electromagnetic mass (`½q²|Δv|²Γ` per impulse; the lattice's heat is not the Joule heat
    under the self-field until impulses carry it: `Δu = Δv/(1 − q²Γ/m)`, to come). The same
    test found the backgrounds counted twice below 2000 particles (PHYSICS.md §2.11), a bug
    of the metals from before this stage, fixed. It had shaped one shipped level: Inverter's
    contacts empty its metals below 2000 particles while its gate is off (457 in the
    preview); its golden hash changed, its reference still solves (the drain off 24.97
    against 24.47 before, on 17.45 unchanged).
  - 2026-10-26: V-D5 changed from a ramp's rise time to the impedance at one frequency
    (lock-in), on a thin ring (wall 1, a lumped inductor). The rise-time integral carried the
    ring's thermal current long after the rise: ±1 a seed on τ ≈ 2, so nine seeds under either
    kick order could not decide (1.18 against 2.00 with the new, a three-seed 2.89 against 2.35
    with the old). The impedance gives L/R from the phases within 0.3 % of L from the field's
    energy over R from the drive's work. Lesson kept in the test: a mean square of a measured
    current holds its thermal part (`kT/L`), so powers and energies divide by the driven
    amplitude's `|I_ω|²/2`.
  - 2026-10-27: a tube level turns the self-field on by its speed of light
    (`TubeSpec::speed_of_light`; absent, the electrostatic model and every level file
    unchanged). Its grid: 0.25 apart as validated (V-D4–V-D6), the board and 4 round it, sides
    rounded up to powers of two, centred (`self_field_grid`; a 40 × 20 board 255 × 127 nodes);
    the verification's refinement keeps it. The level's issues hold its step to 0.015 (first
    order: Y24's 5 % from there), its speed to 0.5–10⁴ (the editor's range) and its grid to 511
    nodes a side. The sandbox edits it (the tube's section).
  - 2026-10-27: drawn on the magnetic map from the engine itself: each frame samples `B_z`
    on the board 0.5 apart (`TubeFrame::self_b`, 13 KB a frame on a 40 × 20 board), which the
    map interpolates and averages as the frames' other fields, converted to the maps' Gaussian
    units (`c` times the engine's). The alternative, `tube_em`'s one-way FDTD field, would show
    a second model of the same field; the energy flow, which needs E with B, stays `tube_em`'s.
  - 2026-10-27: the solver's preconditioner faster, bit for bit the same: its factor per sine
    mode computed once with the grid, and the rows' and the columns' transforms in parallel
    (rayon; each row the same arithmetic on any thread). Measured: the transforms were 78 % of
    a self-field run; their time halved (Ohm's law's strip, 6 time units: 42 s → 25 s on a
    loaded machine). The operator's stencils in parallel were slower (rows too short); two
    rows packed in one complex transform (another 2×, the rounding changed) left for later.
  - 2026-10-27: the step limit 0.015 → 0.012 (`SELF_FIELD_MAX_STEP`). With the impulses'
    currents fixed (DEC-25) the start takes the drawn velocities as the carriers' currents, and
    Y24's return error became 4.0 % at dt 0.01 and 7.8 % at 0.02 (3.1 % and 6.1 % before): the
    5 % target holds up to dt ≈ 0.0125.

## DEC-14: Round metals: a polygon of walls with the disc's area, the disc's exact background

- **Status:** in use (2026-10-26; tests only so far: V-D6)
- **Area:** physics
- **Choice:** `Lattice::with_disc(center, radius)` cuts a rectangular lattice (its bounding
  square) to a disc: the outline a regular polygon of 96 walls (`DISC_SIDES`) whose area
  is the disc's (vertices at `R √(θ/sin θ)`, θ = 2π/96: 0.036 % beyond the radius), the
  cores centred outside it removed, no carriers put outside it, and the background the
  uniform disc's in closed form (`disc_integrals`), not the polygon's.
- **Purpose:** a conductor with no surface charge under a field about its centre: the skin
  effect against an exact solution (V-D6); later round metals in levels (Corbino disc,
  rotating fields; task 34).
- **Alternatives:** a true circular wall (a new kind of wall in the engine's reflection loop,
  which handles straight segments only); the polygon's own background (the integrals of a
  uniformly charged polygon, a sum over its sides: exact for the walls, but the disc's is
  simpler and the two differ only in slivers 0.002 thick at R = 6, alternating in sign);
  a staircase of rectangles (holes) cut from the square (rough edges at any cost).
- **Where:** `crates/physics/src/lattice.rs` (`with_disc`, `outside`, `disc_integrals`);
  tests `disc_background_is_the_disc`, `disc_lattice_stays_inside`, V-D6 in
  `crates/physics/tests/darwin.rs`.
- **History:**
  - 2026-10-26: introduced for V-D6. With it, a carrier drawn into a core that straddles a
    hole's edge or the outline is put on the core's inner side (it had been pushed outward,
    into the hole or out of the metal); levels have neither, so none changes.

## DEC-15: Build speed: the default linker and no compiler cache (both measured)

- **Status:** decided (2026-10-26)
- **Area:** tooling
- **Choice:** keep the MSVC linker and build without a compiler cache.
- **Purpose:** the owner asked for faster builds (2026-10-26); two established tools were
  measured first.
- **Alternatives (measured on this machine):**
  - `rust-lld` as the linker (`target.x86_64-pc-windows-msvc.linker`): the generator's
    release rebuild after a touched file 12.7 s with either (thin LTO and one codegen unit:
    the time is code generation, not linking); the game's development rebuild 4.8 s with
    MSVC's linker, 4.5 s with `rust-lld`. Changing the linker also recompiles everything once
    (4 min 20 s). Not worth it.
  - `sccache` as `rustc-wrapper`, to share compiled dependencies between worktrees and build
    directories: two release builds of the generator in fresh directories, cold then warm,
    40 s and 37 s with 0 hits out of 92 compiles (the arguments it hashes carry each build
    directory's absolute paths, so another directory never hits without path remapping), and
    the generator's dependencies build in under a minute anyway. Not adopted; uninstalled.
- **Where:** nothing configured (no `.cargo/config.toml`).
- **History:**
  - 2026-10-26: measured and decided. What costs time is the heavy per-level tests, which
    test selection (DEC-16, to come) addresses.

## DEC-16: Tests in two tiers; the heavy ones one per check and level, selected by what they executed

- **Status:** in use (2026-10-27)
- **Area:** tooling, tests
- **Choice:** two tiers. The fast tier (every unit and integration test that runs in seconds
  to minutes) runs whole before every commit and in CI (`cargo test --workspace`).
  The heavy tier is one test per check and level in `crates/level/tests/per_level.rs`
  (`golden::`, `repeated::`, `reference::`, `radiation::`, `needs_radiation::<level>`: 521
  tests), data-driven through `libtest-mimic` (`harness = false`) and ignored by default, so
  that any of them runs by name (`-- --ignored <part of a name>`, or nextest's filters), each
  reported as it finishes, four at a time (they share the process's rayon threads;
  `.config/nextest.toml` reserves four threads each). The pass records
  (`level::run_cache::TestPasses`) stay, so an interrupted run resumes; the golden check
  compares a recorded hash with the golden file all the same. `UPDATE_GOLDEN=1` is a mode of
  the same program, writing the file once. Selection: `scripts/affected.py` keeps, for each
  heavy test, the source files its last run executed (LLVM source-based coverage of an
  instrumented build, one process per test) with their hashes, and the inputs that are not
  code (the level's physics, its golden entry, `Cargo.lock`, the toolchain); it runs every
  test whose recorded inputs differ from the working tree's (instrumented, refreshing its
  record) and writes the commit's change note (files changed, what they affect, tests run,
  results). The runs are deterministic, so a test whose inputs are bit for bit unchanged can
  only give its recorded result.
- **Purpose:** the owner (2026-10-26): the per-level checks had grown to ~1.5 h and reran in
  full after any change to the engine's or the levels' code (the pass records are keyed on the
  test program as a whole); test only what a change affects, note per change what it
  affects, and run everything when a bug escapes the selection. Their progress was counted
  per test binary, one of which held an hour of levels ("hangs on one item for several
  hours"): one test per check and level makes it linear.
- **Alternatives:** static selection from the module graph (STARTS-like: cheaper, but trait
  impls, macros and generics make Rust's graph unsafe to under-approximate); a hand-written
  map from level kinds to modules (goes stale silently); selection per crate (cargo's
  granularity: any change to `physics` reruns everything); predictive (learned) selection
  (needs a history this project does not have); Bazel (a migration out of proportion);
  per-level loops inside one test each (the status quo before: correct but coarse and opaque);
  nextest alone without `libtest-mimic` (it runs tests, it cannot make one per level).
- **Full runs:** when a global input changes (the toolchain, `Cargo.lock`, the build profile:
  every record then differs), on request, and when a bug escapes the selected tests (logged
  with the dependency the records missed). No fixed cadence (the owner, 2026-10-26: "every
  few levels" may be superseded by a better alternative). CI checks the golden hashes of one
  level per feature on Linux (T11 across platforms).
- **Where:** `crates/level/tests/per_level.rs`, `.config/nextest.toml`,
  `.github/workflows/ci.yml`, `scripts/affected.py`, `docs/TESTING.md`; CLAUDE.md § Workflow.
- **History:**
  - 2026-10-26: chosen with the owner (Ekstazi, Gligoric et al., ISSTA 2015; Microsoft's Test
    Impact Analysis).
  - 2026-10-27: the per-level tests (from `levels.rs`'s and `determinism.rs`'s loops, bodies
    unchanged); found while testing them that skipping a golden check on its pass record would
    no longer compare the hash with the golden file (a planted wrong entry passed): the
    record now stands in for the computation only.
  - 2026-10-27: the selection (`scripts/affected.py`). Its first rule (a test depends on the
    files it executed, plus those without any code region) missed constants: `reference::`
    reads `level::cost::STEPS_LIMIT` without executing a function of `cost.rs`, as constants
    are folded in at compile time. Replaced by declarations against bodies: a change outside
    the functions' bodies (and in `const fn` bodies) selects every test, a body only the
    tests that executed its file; the unit tests' files are left out (not in the program).
    Measured on two levels (`docs/TESTING.md`): each kind of change selected what it should.
    Coverage costs: with every crate instrumented, 6.5 times the plain build's time on one
    thread (a tube level 38 s against 5.8 s); with the workspace's crates alone
    (`RUSTC_WORKSPACE_WRAPPER`: the dependencies' inner loops, glam's and libm's, without
    counters) 2.5 times (14 s); more on several threads, the counters being shared memory
    (23 s at 4), so each instrumented test gets one thread and the machine's threads run as
    many tests. Two
    platform details: the merge-pool profile pattern `%m` writes nothing on Windows (`%p`
    used), and `libtest-mimic`'s `exit()` skips the exit handlers that write the profile (the
    program returns its exit code instead).
  - 2026-10-27: verdicts and coverage taken apart. The first recording, in chunks by check,
    would have taken about ten hours instead of the one or two estimated: each chunk ended
    with a few metal tube levels running on one instrumented thread each (tens of minutes)
    while the other cores waited, and every selected run would have cost 2.5 times the plain
    one. Now `run` gives the verdicts with the plain program at full speed, then takes the
    coverage, the longest tests first over every core (`--no-coverage` defers it to
    `refresh`); a test whose coverage is older than its verdict depends on every engine file
    until refreshed, which keeps the selection sound.
  - 2026-10-27: `repeated::` takes `golden::`'s coverage. The first finished test of the
    recording, `repeated::109_ohms_law`, took 34 minutes on one instrumented thread (57 times
    its plain run at full speed), putting the whole recording at 3.5 to 4 hours. It computes
    `golden::`'s hash twice on the same level, the runs deterministic, so it executes the
    same files: measured identical on 39 and 99, whose two checks had been recorded
    independently. It now needs only its plain verdict, and its coverage follows its
    source's when both are of the same snapshot (else it runs instrumented).
  - 2026-10-27: the records shared by every checkout, in the main one's
    `target/affected/records.json` (`EM_AFFECTED_RECORDS` to name another): each worktree
    had its own, so a recording made in one (10 h for all the heavy tests' coverage) was
    unknown to the next, which selected every test again. Builds stay per checkout (each
    compiles its own sources). A checkout's own records are read once where the shared file
    is missing. One heavy run at a time (two would lose one's writes). `refresh_coverage.py`
    (and `.cmd`) runs `affected.py run` unattended, keeping Windows awake.

## DEC-17: A pinned toolchain and CI's hygiene checks

- **Status:** in use (2026-10-26)
- **Area:** tooling
- **Choice:** `rust-toolchain.toml` pins Rust 1.98.1 (with rustfmt, clippy and the LLVM tools
  for coverage), and CI installs that toolchain instead of the moving `stable`. A CI job
  `hygiene` (Linux) runs cargo-deny (`deny.toml`: permissive licences only, the security
  advisories, the sources), cargo-machete (unused dependencies) and the documentation with
  warnings as errors.
- **Purpose:** the owner (2026-10-26): reproducible builds and a CI that checks what matters.
  The golden hashes, the run cache and test selection (DEC-16) assume the compiler does not
  change under them; an upgrade is a deliberate bump followed by every test.
- **Alternatives:** `stable` in CI (results could change with a release, unnoticed until a
  golden hash moved); a nightly pin (needed by nothing here).
- **Where:** `rust-toolchain.toml`, `deny.toml`, `.github/workflows/ci.yml`.
- **History:**
  - 2026-10-26: introduced. cargo-deny found, all through Bevy and egui, a decompression bomb
    in `ruzstd` 0.8 (RUSTSEC-2026-0344; the fix, 0.9.1, needs a Bevy upgrade; the game decodes
    only its own bundled assets) and three unmaintained crates (`paste`, `rustybuzz`,
    `ttf-parser`), ignored with those reasons, and two yanked ones, updated (`yoke-derive`
    0.8.4, `clipboard-win` 5.4.3). cargo-machete found no unused dependency; the documentation
    builds without warnings. CI's test job still runs every per-level check (hours on its
    runners); the per-level harness (DEC-16) makes them selectable and its sampled golden
    hashes on Linux the cross-platform check.

## DEC-18: Engine structure for safety (one force pipeline, a Tube builder, a named Material)

- **Status:** in use (2026-10-26)
- **Area:** code structure (physics, tests)
- **Choice:** (1) the tube step's fields come from one place: `charges_fields` (the trees or
  the direct sums, neither with the backgrounds) and `with_external` (the backgrounds and a
  changing field's induced part) for the push and the display alike, with test
  `trees_and_direct_sums_agree`; (2) `Tube::new(electrodes, softening, dt)` with every
  optional part empty, literals naming only what they set (`..Tube::new(...)`); (3)
  `lattice::Material` names the four numbers of a conductor's make-up in `Lattice::new`;
  (4) the metal tests share `tests/common` (their metal, the small floating electrode, the
  Boltzmann mobility, libm complex helpers on num-complex).
- **Purpose:** the owner (2026-10-26) approved the refactors after this session's bugs:
  the doubled background came from the fields being assembled in two places; the Darwin
  field needed 18 literals edited; nine positional arguments of four alike floats invited
  swaps; the mobility had been copied between test files.
- **Alternatives:** leaving them (each bug class stays possible); a full builder type with a
  method per field for `Tube` (more code than the struct update syntax for the same safety).
- **Where:** `crates/physics/src/tube/` (`fields.rs`: the pipeline; `mod.rs`: `Tube::new`),
  `crates/physics/src/lattice.rs`, `crates/physics/tests/common/mod.rs`.
- **History:**
  - 2026-10-26: done in four commits, each bit for bit (golden hashes of representative
    levels, the shared helpers against the old ones bit for bit). The consistency test's
    threshold came from `tree2d`'s documented truncation (5e-7 a node), after a comment in the
    step claiming ~1e-11 proved wrong (measured 1.2e-8 of the largest field).

## DEC-19: The largest engine files split into modules by concern

- **Status:** in use (2026-10-27)
- **Area:** code structure (physics, level)
- **Choice:** the three largest engine files became directories of modules: the types and
  the main loop stay in the parent (`mod.rs`, `lib.rs`), each concern's methods and helpers
  move to a child module. `physics/src/tube/`: `circuits`, `fields`, `emission`, `metals`,
  `darwin` (`mod.rs` keeps the types, `Tube::new`, `state` and `step`). `level/src/tube/`:
  `frames`, `run`, `sim`, `verdict`, `geometry`, `build`, `verify`, `issues`, `tests`.
  `level/src/`: `geometry`, `placement`, `issues`, `drives`, `fields`, `scenarios`, and the
  unit tests in `tests`, `flight_tests`, `model_tests`. Pure moves: a private item takes
  `pub(super)` (exactly its old visibility: the old module and its descendants), the parent
  re-exports what was public, so every path stays.
- **Purpose:** test selection by coverage (DEC-16) works per file: a test depends on the
  files it executes, so a 4000-line file made every test depend on any change in it. The
  files had also become hard to find one's way in. Approved by the owner (2026-10-26).
- **Alternatives:** selection per function (the coverage maps report line regions, which
  shift with every edit; per-file hashes are robust); new crates (separate compilation would
  help incremental builds, but the level crate's types and methods are interlocked, and a
  crate boundary changes privacy, which a pure move must not).
- **Where:** the directories above. `physics/src/beam.rs` (2850 lines) was left whole: half
  of it is one function, which only a change of code, not a move, could split, and only beam
  levels execute it.
- **History:**
  - 2026-10-27: moved bit for bit (golden hashes of representative levels, the tests); the
    level crate's split made by a mover that lexes Rust (strings, raw strings, chars,
    comments) and checks that the files after hold the lines before.
  - 2026-10-27: the game's three largest files too, by the same mover: `ui.rs` (3986 lines)
    into `ui/` (`format`, `tooltips`, `labels`, `results`, `controls`, `views`, `plots`,
    `tests`; the 1300-line `contents` stays whole in `mod.rs`: splitting a function is a
    change of code, not a move), `level_editor.rs` (3500) into `level_editor/` (`inputs`,
    `circuits`, `shots`, `bodies`, `tube`, `limits`, `check`, `tests`), `main.rs` (2739,
    now 983) with `capture`, `keys`, `maps`, `physics_sync`, `tube_views`, `camera`,
    `metal_maps` beside it.

## DEC-20: The level designs as a package by arc

- **Status:** in use (2026-10-27)
- **Area:** code structure (level designs)
- **Choice:** `scripts/levels.py` (4858 lines) became the package `scripts/designs/`:
  `common` (units, particles, elements, `Auto` detectors, helpers, running the generator),
  one module per group of arcs (`charges`, `fields`, `matter`, `instruments`, `jackson`,
  `components`), `curriculum` (`ARCS`, `LEVELS`) and `build` (detector placement, solving
  and checking references, writing levels and their texts, the formulas). `scripts/levels.py`
  stays the entry point, with its commands unchanged. Statements moved whole with their
  comments; each module's imports were computed from its global names (Python's `symtable`,
  exact about scopes). Three deliberate edits: `ROOT` one directory up, `run_generator` into
  `common` (two arcs call it), and the real instruments' one use of `LEVELS` imported in its
  function (the curriculum imports that module).
- **Purpose:** approved by the owner (2026-10-26) with the other refactors: one arc's design
  is found and edited in its own file, and a change to one arc reads as such in a diff.
- **Alternatives:** star imports between the modules (shorter, but hides where a helper comes
  from); one module per level (126 files for designs that share helpers by arc).
- **Where:** `scripts/designs/`, `scripts/levels.py`.
- **History:**
  - 2026-10-27: every level built by the old module and by the package (text-only, nothing
    run) and compared in full: 126 alike, the 15 arcs alike; `python scripts/levels.py --text`
    rewrote no level file and rendered the 95 formulas identically; every module's global
    names resolve.

## DEC-21: The scripts type-checked and tested

- **Status:** in use (2026-10-27)
- **Area:** tooling (Python scripts)
- **Choice:** `python -m mypy` (`mypy.ini`: mypy's defaults, which check the bodies of
  annotated functions; imports resolved from `scripts/`) and `python -m pytest` (`pytest.ini`,
  tests in `scripts/tests/`) on every script, locally and in CI's `python` job. The tests cover
  what could go wrong silently: `affected.py`'s view of Rust (a body against declarations,
  `const fn`, comments and layout, braces in strings and chars, lifetimes, nested comments,
  every engine file parsing) and the designs (each builds without running anything, the
  built files follow the curriculum's order, `curriculum.json` is current).
- **Purpose:** approved by the owner (2026-10-26) with the other improvements: the scripts
  had grown to 9300 lines (the levels' single source of truth among them) with no checks, and
  the test selection now depends on one of them.
- **Alternatives:** strict mypy (thousands of annotations to add to designs that are data);
  pyright (equivalent here, one more tool); ruff for style (not asked; the code is consistent).
- **Where:** `mypy.ini`, `pytest.ini`, `scripts/tests/`, `.github/workflows/ci.yml`.
- **History:**
  - 2026-10-27: mypy found five things, none a bug: three missing annotations in
    `gen_dop853_coefficients.py` (its output could not be regenerated to compare: the source's
    server refused the download; the edits are annotations and an assert on a match every
    coefficient name makes), a reused exception variable's name in `level_docs.py`, and
    `ARCS`'s type, now written out.

## DEC-22: Measured values written into the documents by the tests that measure them

- **Status:** in use (2026-10-27), for V-D5 and V-D6 first
- **Area:** documentation, tests
- **Choice:** a test records each number a document quotes (`measured(id, value, digits)`,
  `crates/physics/tests/common`: the value to its significant digits in `target/measured/<id>`);
  the document marks the quote, `<!--m:<id>-->0.027<!--/m-->` (invisible on GitHub and in the
  PDF); `python scripts/measured.py update` writes the recorded values in, `check` lists the
  quotes that differ. The prose around the numbers stays written by hand.
- **Purpose:** approved by the owner (2026-10-26) with the other improvements: the numbers in
  PHYSICS.md's tables were copied from test output by hand, and drifted (this session's V-D5
  and V-D6 rows were rewritten twice as the method changed).
- **Alternatives:** generating whole tables (the rows' setups and references are prose that
  no test knows); the tests asserting the documented numbers (a test would then fail on a
  documentation lag, not on physics); splitting PHYSICS.md into chapters (proposed with this;
  not done: 325 references cite its sections by number from the code and documents, and a
  generated PHYSICS.md would stop being the source the rules name).
- **Where:** `scripts/measured.py`, `crates/physics/tests/common/mod.rs` (`measured`,
  `significant`), the markers in PHYSICS.md (§2.12: V-D5, V-D6), `scripts/tests/test_measured.py`.
- **History:**
  - 2026-10-27: introduced on V-D5 and V-D6; the other tables move over as their tests are
    next run.

## DEC-23: Long runs of independent items under one controllable runner

- **Status:** in use (2026-10-27)
- **Area:** tooling
- **Choice:** `scripts/runner.py`, a library and a command line for any run of independent
  items (tests, scans, design runs): each item its own process, N at a time; a ledger per job
  (`target/progress/<job>.ledger.jsonl`) taking each outcome as it comes, keyed by the item's
  fingerprint (its command and the contents of its declared inputs); controls while it runs
  (`runner.py ctl <job> pause|resume|stop|skip ITEM|jobs N`); an item past four times its
  expected time is hung: by default the job stops there (the next items may hang too), the
  incident log names the item, its time against its limit and its output's end; run again
  after the fix, the stopping item runs first, finished items whose fingerprint still matches
  are not rerun, and those the fix changed are; `state --save/--load` keeps and restores a
  ledger. `affected.py` runs its tests through it (jobs `affected-verdicts`,
  `affected-coverage`; its coverage records, finer than file hashes, stay its validity);
  `progress.py` shows its jobs as its own.
- **Purpose:** the owner (2026-10-27): a run where one item hangs or fails had to be killed
  whole and redone (this session's recordings twice); skip, pause, resume, save and load at run
  time, nothing finished lost, a hang stopping the run with a record of where, and the reruns
  after a fix redoing only what the fix affects; reusable for anything new of this kind.
- **Alternatives:** nextest alone (per-test processes and timeouts, but no ledger across runs
  nor validity by inputs, and tests only); GNU parallel and job queues (not on this Windows
  machine, no fingerprints); a skip that goes on past a hang by default (the owner: the next
  items may hang too, wasted time).
- **Where:** `scripts/runner.py`, `scripts/tests/test_runner.py`, `scripts/affected.py`;
  docs/PROGRESS.md; CLAUDE.md (long computations).
- **History:**
  - 2026-10-27: introduced, a hang stopping the job by default (the owner's revision of a first
    design that skipped it and went on).

## DEC-24: A carrier's impulses carry its own inductance on the grid

- **Status:** superseded by DEC-25 (2026-10-27)
- **Area:** physics, numerics (the self-field, DEC-13)
- **Choice:** when a collision or a wall changes a carrier's velocity by Δv, its canonical
  velocity takes `Δu = Δv + q(q/m) g Δv`, `g` the 2×2 self-response of the field solver at the
  carrier (its own vector potential per unit of its canonical current, through the deposit's
  and the interpolation's stencil, screened by the other carriers): so its next velocity,
  `u − (q/m)A`, has the jump whole. `g` comes from two solves per metal at its carriers'
  response per area (a unit current on one x-edge and one y-edge at the grid's centre), read
  at any carrier by translation (`darwin2d::SelfResponse`, held in `TubeState::self_response`).
- **Purpose:** with `Δu = Δv` the solve realised only `(I + q(q/m)g)⁻¹ Δv`, and each jump lost
  `½q²Δv·gΔv` to the macroparticle's own field: a self-energy of the macroparticle (its charge
  that of many electrons, which jump independently) that the electrons do not have. The
  walled disc of Y25 lost 14.6 % of the field's energy per time unit, the carriers ~2 % of
  their kinetic energy, and the lattice's heat was not the Joule heat. In a one-carrier model of
  the solver's functional the energy is `(m + q²g)v²/2`, which the corrected reversal keeps.
- **Alternatives:** removing each carrier's own field from its velocity (`v = u − (q/m)(A −
  A_self)`: breaks the scheme's variational structure, whose energy the tests check); a solve
  per carrier for an exact `g` (thousands a step); the bulk `β` including the carrier's own
  (taken: its share of the screening over a skin depth is ~0.5 %).
- **Where:** `crates/physics/src/darwin2d.rs` (`SelfResponse`, `Darwin2d::self_response`),
  `crates/physics/src/tube/metals.rs` (`scatter`), `crates/physics/src/tube/mod.rs` (`state`);
  tests `self_response_matches_a_direct_solve`, Y25; PHYSICS.md §2.12.
- **History:**
  - 2026-10-27: Y25's walled disc 14.6 % → 4.6 % of U per time unit, the closed disc 3.5 % →
    2.1 % (its core collisions too). What remains fits `α dt + β/N` (α ≈ 82 % of U per time per
    unit dt, β ≈ 720 % from ~2800): the step's first-order error now dominates at dt 0.02, so
    Y25's criterion of 1.6× with four times the carriers became 1.6× with the step halved as
    well (measured 2.41×).
  - 2026-10-27: a particle given off (a contact's carrier, an emitter's electron) carries its
    own inductance from the start too, `u = v + q(q/m) g v + (q/m)A` (its metal's g, else the
    vacuum's, a third solve with β = 0; the kernels made before the metals are filled): its
    current appears at once as an impulse's. Found by the heat a strip with its supply off
    booked under the self-field (−6 to −12 per time unit against −3 without), which it
    reduced only partly: an excess of ~3–5, `∝ 1/c²`, not following the step, remains
    unexplained (PHYSICS.md §2.12). Y25 2.3 % and 0.95 % refined (2.7×), walled 4.5 %.
  - 2026-10-27: superseded by DEC-25. The kernel's uniform medium screens a carrier's field
    everywhere, but a metal ends at its edges: against direct solves with the carriers' own
    response, the true `g` is 1.13 times the kernel's on average in Y25's disc at c = 3 (1.25
    at its rim, 1.05 at its centre; 1.07 at c = 10). Each collision lost the difference: the
    disc with its cores lost 5.1e-3 of its kinetic energy per time unit, a loss linear in a
    factor on `g` and zero at 1.13, and Ohm's law's strip at rest booked −8.4 of heat per time
    unit at c = 3 against −3.3 without the self-field, the excess the last entry had left
    unexplained. A kernel for the actual medium would depend on the position (hundreds of
    solves per metal) and still miss the carriers' own fluctuations.

## DEC-25: A velocity set at once enters the next field solve as a fixed current

- **Status:** in use (2026-10-27)
- **Area:** physics, numerics (the self-field, DEC-13)
- **Choice:** each particle records what set its velocity since the last solve
  (`TubeState::set`, `VelocitySet`): nothing, given off (a contact's carrier, an emitter's
  electron, a metal's carriers at the start), or a jump Δv (a collision, a wall), its `u` kept
  from before it. The next solve takes such a particle's current fixed, without its response
  `q²/m`, the other carriers responding to it, and sets its `u` from the solution,
  `u = v + (q/m)A(x)`: a given-off one at the velocity it was given; a jumped one at the
  velocity it would have had without the jump, `u − (q/m)A⁽⁰⁾`, plus the jump, `A⁽⁰⁾` from a
  first solve with every carrier responding (`darwin_step`: two solves in a step with any
  jump).
- **Purpose:** a velocity jump realised exactly, the induction kept. On the grid a carrier's
  own field acts as an inductance: a change of `u` moves `v = u − (q/m)A` by only
  `(I + q(q/m)g)⁻¹` of it, and the part not realised is energy lost to the macroparticle's own
  field, a self-energy the electrons it stands for do not have. Fixing the current realises the
  jump to the solve's tolerance, with the real medium's screening (its edges, the carriers'
  own fluctuations); the first solve gives the jumped particle the induction of the other
  carriers' changing current (an inductor's back-EMF) as every carrier feels it. The disc of
  Y26 with its cores keeps its energy and booked heat within 7.4e-4 of its kinetic energy per
  time unit at dt 0.02 (5.1e-3 with DEC-24's kernel), Y25's drift no longer depends on the
  carriers' count (the walled disc 0.19 %, 4.5 % before), V-D5's L/R by the phases is +4.6 %
  of the one by the energy, its R the same with and without the self-field within 0.5 %
  (5.6 % before), and Ohm's law's strip at rest books the same heat with the self-field as
  without, within the runs' spread.
- **Alternatives:** DEC-24's uniform medium's kernel (13 % short at a metal's edges); a kernel
  for each metal's actual medium (position-dependent: hundreds of solves per metal, and still
  blind to the fluctuations); a first-order (Born) correction of the uniform kernel for the
  medium's edges (approximate where the correction is large); the fixed current at the drawn
  velocity alone, one solve (measured: the jumped particles miss a step of induction, a
  fraction `ν dt` at every step, and V-D5's inductance by the phases came out 13.1 % low);
  that with the induction of the particle's last step added, the change of A along its path
  (free, but the change holds its own field's, which varies as it crosses the cells: Y26
  failed, 1.5e-3). The cost: a second solve as long as the first (13 and 14 iterations to
  1e-10), the self-field's time doubled (Ohm's law's strip at c = 3: 92 s for 20 time units
  against 50 s with one solve, 7 s without the self-field); to be won back in the solver
  (its tolerance, its preconditioner).
- **Where:** `crates/physics/src/tube/darwin.rs` (`darwin_step`, `deposit`),
  `crates/physics/src/tube/metals.rs` (`scatter`), `crates/physics/src/tube/emission.rs`
  (`add`), `crates/physics/src/tube/mod.rs` (`VelocitySet`, `TubeState::set`); tests Y25, Y26,
  V-D5; PHYSICS.md §2.12.
- **History:**
  - 2026-10-27: chosen, replacing DEC-24; `darwin2d::SelfResponse` and its test removed. The
    one-solve form came first and V-D5 found its missing induction; the two-solve form
    replaced it before the commit.
  - 2026-10-27: the solve's tolerance kept at 1e-10. At 1e-6 the two solves take 8–9
    iterations against 13–14, and Y21–Y26 do not change (nor at 1e-4), but V-D5's driven ring
    fails: its R with the self-field 1.195 times the one without, its magnitudes 1.137 times
    its phases', its heat/work 1.095 times (bounds 15 %, 5 %, 5 %; 1.005, 0.959, 0.984 at
    1e-10). The relative residual weighs the ring's large-scale current little against its
    error. The second solve's time is to be won back by a preconditioner with the carriers'
    screening (#55).
  - 2026-10-27: nor a preconditioner with the carriers' screening: a uniform `β̄`
    (`μ₀/(λ(λ + μ₀β̄))`, from half to eight times the grid's mean β) takes 14.5 to 21.4
    iterations against the vacuum's 13.5 (Y26's disc). The time is in each iteration's own
    cost (#55).

## DEC-26: Leads and straps routed as circuit wiring, the level's geometry

- **Status:** in use (2026-10-27)
- **Area:** display, the tube field (`tube_em`), level geometry
- **Choice:** every lead (the path of an electrode's circuit current, from its enclosing
  grounded ring or the board's border) and every strap (between a floating circuit's
  electrodes) is routed once per placement by `level::leads` (`Level::tube_leads`, carried by
  `TubeRun::leads`): on a grid half a cell apart over the board, a shortest-path search in
  whole steps where a bend costs 4, a stub meeting its surface on a slant 2, each step off a
  box's face's middle 1 and a node another wire took 60 (a box met on its face, never at a
  corner); a node within 0.3 of an electrode or a metal is closed to every wire but those
  ending on it; straps first, then the leads, the shortest first; integer costs and ties by
  index (deterministic). `tube_em` drives the leads' currents along these polylines (each end
  snapped to a node of its metal); the game draws leads and straps in every view, beads on the
  leads in the magnetic and energy-flow maps.
- **Purpose:** the owner (2026-10-27): the cables were drawn only in the magnetic and energy
  views; they should be drawn throughout, not overlap, and run along x and y with right-angle
  bends unless there is a reason otherwise, as in circuit diagrams. One routing for the field
  and the drawing keeps the wire drawn the one whose current the field shows.
- **Alternatives:** the straight shortest segment (before: drawn only where the field ran,
  lines at any angle, crossing one another); a channel or visibility-graph router (finer on
  crowded boards, much more code; the grid search wires every shipped level without a shared
  node in under 10 ms); routing in the game alone (the field would carry the currents along
  other paths than those drawn).
- **Where:** `crates/level/src/leads.rs` (tests in `leads/tests.rs`:
  `every_shipped_tube_level_is_wired`, `a_free_electrode_gets_a_straight_lead`,
  `leads_go_round_what_is_in_the_way`), `crates/level/src/tube_em.rs`,
  `crates/game/src/draw.rs` (`draw_leads`, `draw_straps`), `crates/game/src/legend.rs`.
- **History:**
  - 2026-10-27: chosen. The first tie-break (by index) put Magnetron's antenna feed half a
    node off its axis, slanted onto the ring and the rod (`magnetron_feed_field_is_smooth`
    found the field's jump off the line it reads across): the slant's cost now keeps a wire
    on a round electrode's axis where it can. The first captures showed leads meeting boxes
    at their corners and running on along their edges (Positive carriers' end contacts'
    leads continued the strips' bottom edges across the board, one line with them): a wire
    now meets a box's face clear of its corners (by the clearance, or half a face shorter
    than twice that) and pays a step for every grid step off the face's middle.
  - 2026-10-28: the beads in every view (the owner: the leads were shown in every view but
    without the charge moving through them, as the magnetic and energy maps show it). The
    leads' charge accounting moved out of `TubeEm` into `LeadAccounts`, shared by the field
    (bit for bit as before) and `LeadFlow`, which reckons each lead's charge frame by frame
    without the field, as far as the playback has reached; the field maps keep the field's
    own currents, in step with the field they show.

## DEC-27: The leads' circuits under the self-field, closed along the envelope

- **Status:** in use (2026-10-27)
- **Area:** physics (the tube engine's circuits and the Darwin field), level geometry, display
- **Choice:** with the carriers' own field, every electrode with a lead from the border and a
  supply to ground is a circuit node whose branch runs through the supply (its circuit, or
  its fixed potential as an ideal source) and along the lead; the lead's current is a source
  of the field, and the field's EMF along the lead, `−dΛ/dt`, `Λ = ∫ A·dl`, is in the branch's
  trapezoidal equation: the leads' own inductances `Mᵢⱼ` (one solve per lead at the start,
  the metals' carriers screening as in the steps' solves) implicit with the circuits, the
  rest of the flux along each lead (the carriers' currents') from the last solve, a step
  behind (`physics::tube` `leads.rs`). Every circuit is closed: each lead returns along the
  grounded envelope's inner face, the shorter way round, to a ground terminal below the
  middle of the board's lower edge (`LeadPlan::ground`, `returns`, `circuit_paths`), drawn
  with a ground symbol and in the legend. A lead from a ring, a circuit between two
  electrodes and a strap are refused with the self-field until modelled.
- **Purpose:** DARWIN.md §2.3 and step 4 (the owner's "how the whole circuit works at once and
  what field it emits", 2026-10-26): a circuit's inductance from its geometry, its leads'
  field acting on the carriers and theirs on it. The closure is what makes that inductance
  physical: an open lead's current comes from nowhere, its field falls as 1/r and its
  energy grows with the log of the grid's frame's distance (16.5 in a frame 6 away, 148 at
  50, for the same lead), so a circuit closed only through infinity (Ohm's law's leads from
  opposite sides) would have had a numerical parameter for its inductance. The envelope is
  the electrostatics' ground already (every tube level's, `ENVELOPE_MARGIN`), so the return
  runs where the ground conductor is.
- **Alternatives:** open leads (the display's convention, `tube_em`; box-dependent, above);
  the EMF explicit from the last steps (unstable: a loop's L/R of ~100 against a step of
  0.01); the carriers' flux implicit too (a joint solve of the circuits and the field,
  coupling the dense node system to the conjugate gradients; Y27 shows the lag's cost
  unmeasurable, 0.004 % of the field's energy per time unit); the leads to one common terminal
  routed by the router (tidier loops, but the routing of every level changes; the returns
  keep DEC-26's leads as they are); the return along the board's edge rather than the
  envelope (half a cell inside the conductor that carries it); a perfectly conducting
  boundary for the field at the board's edge (the solver's sine transforms need the frame's
  `χ = 0, B_z = 0`).
- **Where:** `crates/physics/src/tube/leads.rs` (tests there: `two_leads_make_the_loop_inductance`
  against `scripts/wolfram/d8_two_leads.wls`, `leads_inductance_is_the_field_energy`,
  `a_lead_makes_an_rlc_circuit`, `a_lead_lc_keeps_its_energy_and_runs_back`,
  `a_lead_field_is_odd_under_a_mirror`), `tube/circuits.rs`, `tube/darwin.rs`;
  `crates/physics/tests/darwin.rs` (Y27); `crates/level/src/leads.rs` (`border_path`,
  `returns`, `circuit_paths`; tests `returns_go_round_the_border_the_short_way`,
  `self_field_circuits_return_to_the_ground_terminal`), `tube/build.rs`, `tube/issues.rs`
  (test `self_field_tube_takes_the_paths`); `crates/game/src/draw.rs` (`draw_returns`),
  `legend.rs`.
- **History:**
  - 2026-10-27: chosen. Y27's first setup had no envelope and ran away (both electrodes'
    potentials falling as their currents rose: the open world's negative capacitance for
    a common charge, OP-1); with the envelope every level has, it held.
