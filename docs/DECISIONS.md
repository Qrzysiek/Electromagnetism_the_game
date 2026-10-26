# Design decisions

Every design choice made in building the game: the model, algorithm or design taken, its
purpose, what else was considered and why not, and where it lives. When a choice is changed,
its entry is not rewritten: the change is appended to its history with the reasons, so that
the log shows how the game got where it is (the owner, 2026-10-26: "make a note of every design
choice and every time you change it ... append this note with the same info and reasons for
change").

Scope: features, models and algorithms (physics, numerics, display, tools), and the main
choices behind each level's design (the full measurements stay in `scripts/levels.py`'s
comments and `docs/LEVELS.md`). The physics' models and their validation are in `PHYSICS.md`,
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
- **Where:** `crates/level/src/tube.rs` (`frame_between`, `taken`, `given_off`,
  `nearest_contact`, `contact_faces`); test `contacts_take_and_give_off_carriers_between_frames`.
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
- **Where:** `crates/game/src/scales.rs` (`drive_field`), `crates/game/src/ui.rs` (the mean
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
- **Where:** `crates/physics/src/tube.rs`, `crates/level/src/tube.rs`; PHYSICS.md §2.9, §2.11.
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
- **Where:** `scripts/levels.py` (`amplifier`, its comment holds the scans), docs/LEVELS.md.
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
- **Where:** `scripts/levels.py` (`wheatstone_bridge`), docs/LEVELS.md.
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
- **Where:** `scripts/levels.py` (`ringing_tank`), `target/scratch/rlc_proto.py`.
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
- **Where:** `scripts/levels.py` (`field_tuned_damping`), `target/scratch/ftd_proto.py`.
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
- **Where:** `scripts/levels.py` (`lock_in_hall`), `target/scratch/lockin_proto.py`.
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
- **Where:** `scripts/levels.py` (rectifier, klystron, inverter, ringing_tank),
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
