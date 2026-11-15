# Handoff (2026-11-14): where the work stands and what comes next

Read this instead of the previous sessions' transcripts. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; the metals' plan: `docs/DRUDE.md`
"Implementation plan"; decisions: `docs/DECISIONS.md` (DEC-43 the latest).

Levels are named here, not numbered: the numbers move when a level is added (the Corbino
disc moved every later level up by one: Self-inductance from 125 to 126). The task numbers
quoted (#34, #61) are a past session's list, kept as labels only.

## Waiting on the owner (the session of 2026-11-15, 04:20–12:00)

The owner's instructions for that session: make progress on this plan until noon, run no
new test suites (the next night runs them), take the recommended choice where confident and
list the others here.

- **Tests owed by the night run.** The commits after `1d483cc` were made without the fast
  tier, pytest or the affected heavy tests (each commit message says what it skipped). Before
  pushing: `cargo test --workspace`, `py -3.14 -m pytest`, `python scripts/affected.py run`
  on the last commit.
- **Nothing pushed.** `1d483cc` (the thermal-wall contacts, its 529 heavy verdicts passed)
  and the later commits are on `worktree-tubes-next` only; pushing, fast-forwarding main
  and rebuilding the game there wait for the owner (as the session of 2026-11-14 left it).
- **The thermal-noise level's capacitance control** (Next, item 1;
  `prototypes/noise/LEVEL_CHECKS.md`): the plan's control, a grounded plate beside the
  player's floating contact, moves its capacitance by at most 15 % (the ratio by at most 7 %,
  less than its statistical error: the verdict's error bar is ±19 % on the ratio, ±14 % on the
  decay), so the two goals could not pin `C`; a control must move `C` by a factor of two or
  more. Options: (1) a capacitor as the player's control, a lumped capacitance on the contact
  chosen from a list or a log slider over a decade (the engine has `circuits[].capacitance`;
  a new control kind: recommended); (2) a geometry where a plate dominates (a small contact,
  a long plate, plates allowed nearer: unlikely to reach a factor of two); (3) longer windows
  (halving the bars costs four times the preview, past 30 s); (4) the pace alone, an
  introduction only. The rest is measured and fine (RC about 2, preview 24.4 s, verification
  112.5 s). Not chosen here: a new control kind is a design choice.
- **Glass screen's message changed** (Next, item 2): the level was the dielectrics'
  introduction as "a glass screen"; it now says the glass polarizes and does not screen,
  which is what it does. Chosen without asking, as the only text the physics allows. The
  alternative, a level where a dielectric visibly weakens a field, needs another design (a
  slab tall enough to approach the infinite one is too slow to solve, and even that leaves
  0.73 of the field beyond it); say if you want one.

## Next

Reviewed 2026-11-14 against the tree, the plans and the tools' reports. In order:

1. **The thermal-noise level's measures are written and not committed** (DRUDE.md Phase 2
   item 1). The worktree has held them since 2026-11-05, in nine files under `crates/`:
   `Measure::SpreadRatio` (a contact's noise's spread over a reference contact's),
   `Measure::Decay` (a kicked contact's time to 1/e), `TubeRun::noise_pace` (the noise's
   correlation time, shown beside the goal), the setup issue for a spread without a
   reference, the editor's two controls, the panel's text, `tube-scan`'s `pace`. Their unit
   tests pass (`decay_and_spread_measures`, `correlation_time_of_an_ar1_noise`; run
   2026-11-14).
   - Owed before a commit (CLAUDE.md): DEC-42 in `docs/DECISIONS.md` (the code cites it);
     PHYSICS.md (the measures, the tests' results); the plan's physics tests (equipartition
     on a floating contact against D7's kT; the correlation time against a kicked contact's
     relaxation); `scripts/pairs.py`; the pre-commit checks and the affected heavy tests.
   - The prototypes: `prototypes/noise/` (`proto2.py`: two strips grounded at one end, a
     floating contact on each, a kicker over the player's; `fold.py` and `acf.py`: a
     trace's decay, spread and pace; `runs.out`: each run's `tube-scan` result line). The
     traces themselves stayed in the worktree's scratch, which `cleanup.py run` prunes from
     2026-11-19; they are reproducible from the scripts.
   - Read from `runs.out`, not rerun: the last run with both goals (`n2`) gave a spread
     ratio of 0.622 in the preview and 0.580 at the verification's resolution, a decay of
     6.15 and 6.16, a pace of 8.99 and 9.03. The ratio holds between the resolutions to
     7 %. The pace is 1.46 times the decay, where the level's idea is that they are one
     time (fluctuation–dissipation); the pace's own error over that window is about
     `√(2τ/T)` ≈ 20 %, so it is not settled either way: settle it before the level's text
     says so. The last runs kicked the contact with a square wave (`sq1`, `sq2`).
   - Still to measure, as the plan asks: the capacitances' range, and how many `RC` the
     30 s preview holds.
   - **Found 2026-11-15, before writing the equipartition test: the model's floating contacts
     are not in equilibrium** (`prototypes/noise/EQUIPARTITION.md`, its table and runs). On
     a strip fed by D7's contacts, `⟨V²⟩ C_node/kT` (exact in the model: 1) is 1.18 ± 0.03
     at length 4 and 1.43 ± 0.03 at length 2, the same at half the step, four times the
     macroparticles and a lumped capacitance five times the contact's: resolution-free, and
     the size depends on the feeding resistance, against the level's idea. The same
     contacts pass Nyquist between fixed potentials on the short strip (1.01 ± 0.05). Two
     thermal-reservoir contacts tried (a scratch patch, reverted) move it either way (a
     regular train 0.78–0.92, Poisson 1.15–1.21) and neither gives 1.
     **The cause, found the same day (the owner: research it first):** two departures from
     detailed balance. The contacts give carriers off δ (half a segment, about a Debye
     length) out from the face but take them in at the face, so a contact's electrochemical
     potential follows its face's field (the excess falls linearly with δ); and the refill
     clamp is no detailed-balance reservoir (~15 % at small δ). **A fix, measured:** Poisson
     emission of the half-space's thermal flux, each release weighted by the Boltzmann factor
     of its potential step, to first order `exp(βq 4πσδ)` from the segment's surface density
     (no cost): 0.987 ± 0.018 and 0.990 ± 0.021 at lengths 2 and 4, Nyquist 1.06 ± 0.06.
     Its price: Sharvin's contact resistance (`kT/(q j₀ w)`, ~15 % on a Hall-probe strip),
     which the clamp was designed to avoid. **The owner (2026-11-15): a contact kind for the
     noise level.** Done: `MetalSpec::thermal_contacts` (`OhmicContact::thermal`,
     `Tube::emit_thermal`, an editor checkbox; DEC-43), every other level bit for bit;
     tests D13 (equipartition, two lengths and a lumped capacitance), D14 (a kick's decay
     against the noise's pace), Y31 (charge conjugation, bit for bit), D10 (c) (the energy
     account); PHYSICS.md §2.11 *Thermal noise*; DEC-42 written for the measures; `pairs.py`
     `tube_noise`. **Done 2026-11-15:** the prototype rerun with `thermal_contacts`
     (`proto2.py --thermal`) and the two checks (`prototypes/noise/LEVEL_CHECKS.md`): `RC`
     is about 2 time units, 200 in the window (fine); the grounded plate beside the
     player's contact changes its capacitance by at most 15 %, less than the measures'
     statistical error: **the plan's capacitance control fails**. Waiting on the owner (the
     options in LEVEL_CHECKS.md; recommended: a capacitor as the player's control). The
     preview to t = 420 takes 24.4 s and the verification 112.5 s (fine); the verdict's error
     bars (the ratio ±19 %, the decay ±14 %) ask for a control that moves `C` by a factor of
     two or more. A grand-canonical layer
     (equilibrium without the contact resistance) stays untried: the way to bring every
     contact into equilibrium, should the owner want it.
2. **Glass screen: resolved 2026-11-15** (`docs/levels/glass_screen/NOTES.md`). The
   disagreement was the document draft's bug (each panel's field projected on the source
   panel's normal, not the receiving one's); fixed, the draft's field and flights match the
   game's (the field to ~0.3 %, the reference's end 0.056 apart, the empty board's 0.002).
   It showed that **the level's text was wrong**: the glass does not screen the beam's line,
   it strengthens the push there (+8 % under the charge, +43 % under its ends; inside the
   glass the field falls to 0.39). Text, reality entry, LEVELS.md and CURRICULUM.md
   corrected; physics and golden hash unchanged. The level's document is written
   (`docs/levels/glass_screen/`: eight comparisons, the largest difference 7e-4 against a
   tolerance of 0.55; the infinite slab as the analytic limit).
3. **The two open numerical questions that block levels** (*State* and *The Transformer*
   below): the AC mobility's phase through a level (Carriers' inertia and every phase goal
   at ωτ ≈ 1 wait on it), and the dense coupled rings' dynamics (V-D10: an accurate
   Transformer in Self-inductance's metal waits on it).
4. **Levels the engine can build as it is.** With the frequency slider (DRUDE.md Phase 2
   item 2): Tuned tank, Plasma resonance. With the rings (item 3): the Generator, to check
   first against the cut-ring finding below (at the levels' frequencies the slit's own
   capacitance may carry the current past a load across it; not measured). Phase 1's
   remainder: Relaxation, Calibrated field, Resistive anode (DEC-11, in design: the steering
   too weak so far).
   Two arcs have no master level yet (docs/CURRICULUM.md): the transistors (Flip-flop
   deferred, DEC-4) and AC and induction (Crystal radio, after Phase 2 items 2 and 4).
5. **Phase 2's two extensions not started**, neither waiting on the above: a uniform AC
   field (item 4: Polarization, Faraday cage, Receiving antenna) and temperatures (item 5,
   opening with Hot carriers).
6. **Level documents:** the first 32 levels in curriculum order have one (2026-11-15 added
   Glass screen, Fast lane, First coil, First magnet, Ferrite shield, Stern–Gerlach and the
   beta-ray spectrometer: the first coil carries the magnetic arc's derivations, first
   magnet and Stern–Gerlach check their marches against quadratures of their conserved
   quantities); no tube, metal or AC level does (`eddy2d.py` would be Self-inductance's
   `calc.py`). Next in order: Dempster, Wien filter, Calutron. Seen while building them:
   `deflection_plates`' `calc.py` exceeds the 120 s limit when every core is busy (135–138 s
   during the coverage pass), and every later document that cites it rebuilds it.
   `scripts/level_docs.py` needs `py -3.14` here: the documents' `calc.py` import numpy and
   mpmath, which this machine's `python` (3.15) lacks, and under it `status` shows 14
   documents failing on the import alone. Under 3.14 all 25 pass their checks (2026-11-14;
   a full `status` runs every `calc.py`, some ten minutes).
7. **Later stages, as planned:** interludes (SPEC.md §14; the Transformer is the first
   candidate and waits on their format), optics (SPEC.md §8; what is in place and missing:
   below), DRUDE.md Phase 3 items 2–4.

Housekeeping, none urgent:

- The main checkout has two untracked files: a session's transcript at its root
  (`2026-09-29-…-in-electromagnetism_the_game.txt`, 6 MB) and `docs/CODE_TOUR.md` (a reading
  guide of 2026-10-09, before the tubes and the metals; no tracked file cites it): to commit,
  ignore or delete.
- Left from October: the branch `level-docs/glass-screen-draft` (item 2's draft) and a stash
  on `level-docs/metal-15-26`.
- Fixed 2026-11-15: the side panel overflowed on the metal levels (Hall probe by 9.6 px).
  The cause was the *Scales* header, whose title egui lays out on one line ("ω_c τ = … ·
  qEℓ/kT = … · like germanium"); it is now a galley wrapped to the panel. A capture sweep of
  levels 100–128 shows no overflow. A collapsing header's title never wraps: a long one needs
  the same treatment (`ui/mod.rs`, the Scales header).
- SPEC.md's status lines now cover the tube world (2026-11-15: §3 *Inside the
  components*, §10 stages 4 and 6).

## Waiting on the owner

- Amplifier's preview, 33 s for a 30 s rule (10 % over): whether the allowance for a much
  better level applies (docs/LEVELS.md, *Timing audit*).
- Three solving intervals under the 5 % rule (docs/LEVELS.md, *Showcase review*, with each
  one's options): Klystron (2.0 % and 1.5 %) and Magnetron (3.3 % and 2.5 %), resonances
  narrow by their physics; Long and thin's junction (1.6 %).
- Self-inductance's metal is about 20 times denser than the other levels' (density 41
  against 2; 30 times in its first design), where ten times was accepted (2026-11-02).
- The interludes' format (SPEC.md §14, *Open*).
- Whether electrodes may touch when both are grounded (it would close the Corbino disc's
  outer contact).
- DRUDE.md's list *For the owner, when each step starts*: which of its questions are still
  open now that the arcs and the self-field exist (the gate's settling, the later stages'
  order).

## State

`main` = `worktree-tubes-next` = `a6564e1` (worktree `.claude/worktrees/tubes-next`, which
also holds item 1's uncommitted changes).

- Heavy tests: 529, every one recorded with its coverage at `a6564e1` (`python
  scripts/affected.py status` on the clean main checkout, 2026-11-14: none selected). Item
  1's changes touch the tube engine's level code and will select many.
- Disk (`python scripts/cleanup.py`, 2026-11-14): 13 GB over the two checkouts, nothing over
  its cap.
- Ring-shaped metals in the level format (task #34, DEC-28): `MetalSpec::hole`, one
  rectangular hole with walls at least a cell wide; the builder (`Lattice::with_hole`,
  carriers from the net area), setup issues (the thinnest wall for the mean-free-path rule, a
  metal allowed wholly in a ring's hole), cuts (the metal whose conductor holds the middle;
  past a hole's edge), the editor (a hole checkbox, two node pickers; `check_editable`), the
  display (`MetalView::{hole, holds, pieces, cell_area, background_at}`: the hole's outline,
  the backgrounds wall by wall, map cells counting their conductor only, no drift chevrons on
  a ring), the solver's cuts across a ring's walls, `pairs.py` `tube_ring` (and its cut
  element now cuts a ring's wall). Tests: `rings_circulate` (lattice; `L/P = −2⟨Ḃ⟩/⟨Ḃ²⟩` to
  −0.75 %; a cut ring; the quarter turn, PHYSICS.md Y29), `ring_setup_rules`,
  `ring_holes_are_editable`. Levels without a hole are bit for bit unchanged (every number
  keeps its original arithmetic; the 521 verdicts on the commit say so).
- Contacts on a ring's inner side, DRUDE.md Phase 2 item 3's remainder (an electrode in the
  hole was refused: its lead would cross the ring; in the electrostatic model a lead is only
  drawn, so it is allowed there with the lead drawn across).
  **Done 2026-11-04 (the owner: "whatever makes the most sense"; DEC-40):** contacts in a
  ring's hole (a box flush with one wall, or filling the hole and emitting from every face),
  the lead drawn across the ring (refused with the self-field), ohmic layers shared where
  contacts meet at a corner; D12 (the Corbino disc against the Boltzmann equation's σ_xx,
  four seeds) and Y30 (its mirror and quarter turn); the level **Corbino disc** (after
  Magnetoresistance; the twelve levels after it renumbered, their golden keys renamed). Its
  outer contact stops a cell short of each corner (electrodes may not touch), so its current
  falls somewhat more slowly than σ_xx (LEVELS.md; D12 is the exact test). An idea left: let
  electrodes touch when both are grounded (the infinite field between touching electrodes
  needs different potentials), which would close the outer contact.
  **Found with it, fixed (DEC-41):** the panel's Joule heat (the field's work, DEC-33) fell
  short of `I V` by about a constant per contact (the Corbino disc's panel −14 per time at
  its reference against `I V` 15, its heat map negative along every contact): a contact's
  own loop through its depleted layer (given off at δ, taken back at 1.5 δ). A metal with
  contacts now shows the macroscopic `J̄·Ē` (PHYSICS.md §2.11, *The contacts' layers*):
  Ohm's law +1.85 %, the Corbino disc −2.0 %, at rest 0.36 % of `I V`.
- The thermal-noise level's design (Phase 2 item 1), decided by the owner 2026-11-05 in place
  of the RMS goal (DRUDE.md): two resolution-free goals, two floating contacts' noises as a
  ratio and the noise's pace (its correlation time, about `RC`). Where it stands: *Next*,
  item 1.
- README screenshots (DEC-29): `python scripts/readme_shots.py` retakes them.
- **A supply's frequency on a slider** (DRUDE.md Phase 2 item 2, DEC-30): `limits.supply_omegas`,
  a logarithmic slider on every supply whose target a sine drives, `Element::omega`;
  `generator window --omega`; the slider rule (no reference within 5 % of the log length from
  an end or the middle); tests `crates/level/tests/frequency.rs` (bit for bit the sine's own
  frequency), `slider_ranges_hide_the_answer`; `pairs.py` `tube_tuned`.
- **Open: the AC mobility through a level** misses its 5 % (5.2 % at ω/κs = 0.4, 7.6 % at
  1.2): the drift lags the field 3–4° less than the Boltzmann equation, 2.5–3° more than
  D2's box shows. Data, what was ruled out and what to try: DRUDE.md Progress, 2026-11-02;
  the test's draft: `prototypes/ac_mobility_test.rs`. Phase goals at
  ωτ ≈ 1 (Carriers' inertia) wait for it.

## The finding that changes the Transformer's design

**A cut does not open a ring at a level's frequencies** (PHYSICS.md §2.11 *Rings*, test
`rings_circulate` (c)). A cut ring is a capacitor in series with its resistance: the slit's
faces, screening lengths apart (~0.2 per unit length), in parallel with its two arms facing
each other across the hole and round the outside (a few tenths). At ω = 0.125 a cut ring 12 ×
12 (walls 2) carried 103 % of the closed ring's current though no carrier crossed the cut. A
slit opens a ring only where `ωP/(4πσd) ≪ 1`: the levels' metals have `4πσ ≈ 3`, conducting
only a few times faster than the drives swing (a real metal's `4πσ` is ~10¹⁸ s⁻¹), so a cut,
or a gap a unit or two wide, opens a turn only at ω of a few hundredths, where the self-field's
effects (`ωL/R`, ~0.45 ω for V-D7's secondary at c = 3: L ≈ 100, R ≈ 220 by its two-coil
fit) are gone. A more conductive metal hardly helps: a cut ring is a series RLC, open only
below its resonance `1/√(LC)` (~0.15–0.2 at c = 3) and inductive only above `R/L`; both at
once need `√(L/C) ≫ R` (√(L/C) ≈ 14 here), which the levels' metals miss by an order of
magnitude, and the densest carriers the self-field's step allows (`ω_p dt ≤ 0.2` at dt 0.012:
ρ ≲ 20) only bring R to about √(L/C). (That resonance is itself a level: a split-ring
resonator, the metamaterials' element.)

## Done: Self-inductance (DEC-32), the self-field's heat measure fixed

- The owner's answers (2026-11-02): find the residual's cause first; make the runs faster;
  a metal ten times denser is acceptable (the level's is 30 times: more metal-like, ω_p τ 5.5;
  say so if that is not what the owner meant).
- The residual's cause (DEC-13's history, PHYSICS.md §2.12, test V-D8): the heat measure
  booked the drive's work at the push's velocity, before the self-field's part of the step
  took back its impulse: ~`dt/2τ` of the unscreened heat too much, whatever the density or c.
  Fixed (`Tube::book_work`); V-D8 (a disc's screened heat against the cylinder's) +0.9 %, and
  +10.9 % on the old booking (it fails). Left: a discretisation residual ~0.6 % of the
  unscreened heat (dt 0.012, cells 0.25), felt where the screening leaves under ~10 %.
- With the self-field a heat goal averages over whole half-swings of the drive (DEC-31): the
  ring's stored energy swings by up to `√(1 + (ωL/R)²)` times the heat per period.
- The level (redesigned 2026-11-02 for the owner's 30 s preview, DEC-32's history): a ring
  20 × 20 (walls 2), density 41, c = 27.4, dt 0.0088 (`c/ω_p` 1.21, `ω_p dt` 0.2), B₀ 0.45
  fixed, ω on a log slider 0.05–5, heat 1000–1500 over t = 2.5–20 at the bend (answer ω 1.3;
  verified 1183 ± 84, the preview 1127; ω solves on 1.2–1.6 verified, 6.2 %). Preview 17 s,
  verification 54 s (docs/LEVELS.md holds the level's numbers; to t = 10 and before the
  solver's speedup of 2026-11-04 they were 1195 ± 136, 1.15–1.55, 24.4 s and 62 s). (The
  first design: density 60, c 24, the bend at 0.7, ~75 s.)
  Scales:
  `ωL/R`, the skin depth, the board against c/ω (`self_field_lines`). The exact curve:
  `prototypes/inductance/eddy2d.py`; the scans: `NOTES.md` there.
- Checked: the panel's "Joule heat" (then the lattice's take) read ~750 per time above the
  drive's work on this level: the macroparticles' stirring of each other (`∝ n²/N`; with
  the drive at rest 748 without the self-field, 742 with it, 400 at twice the carriers),
  not the self-field. The panel now shows the field's work, what the goal measures (DEC-33;
  PHYSICS.md §2.12 *The heat at rest*; the runs: `prototypes/joule/`).
- Not written: the level's document (`scripts/level_docs.py new self_inductance`; none of the
  metal or AC levels has one yet); `eddy2d.py` would be its `calc.py`.
- **The timing rules (the owner, 2026-11-02, clarified 2026-11-03; CLAUDE.md):** the
  preliminary answer (the preview to the goal window's end) within 30 s; the animation
  streamed without stutter at speed 1 (a tube level's run at least one time unit a second
  of computing), however long; the verification within a few minutes; each bound exceeded
  by at most 100 % where that makes the level a lot better. Audited since (below; at their
  designs LEVELS.md recorded Eddy currents in silicon at 38 s and Lock-in Hall at 30–70 s:
  the audit measured 18 s and 26 s). Self-inductance streamed at 0.41 units a
  second. The owner (2026-11-04): fewer carriers are fine if the physics barely suffers.
  Measured first: they are not the cost (at a tenth of them the preview took 20 s against
  27); 70 % of a step was the self-field's two solves, 16–17 conjugate-gradient iterations
  each, whatever the carriers. Not changed: the carriers. Then, the owner's choice of more
  solver work: the iteration 0.47 → 0.075 ms (DEC-13's history, PHYSICS.md §2.12), the
  scattering and the frames' middle field in parallel: it streams 1.2 units a second
  (preview 16.8 s to t = 20, the goal's window lengthened for the verdict's noise, DEC-32;
  the commit's 13.8 s was not reproduced).
  The owner (2026-11-04): the animation stood still while the verification ran. Fixed
  (DEC-39): the verification now runs beside the rest of the display run's stream, and a
  run whole at the goal's end is marked done first, so the playback loops; `EM_TIMING=1`
  logs the frame rate and the time shown. The audit of the other tube levels (`generator
  check` with a fresh run cache and `tube-field`, `prototypes/timing/`) is done: the table
  in docs/LEVELS.md, *Timing audit*. Every tube level keeps the rules; the closest are
  Amplifier's preview (33 s, its window one period of a 300-unit signal: within the 100 %
  allowance, to confirm with the owner) and Self-inductance's stream (1.2 units a second);
  the slowest field maps Klystron's (2.8) and Wheatstone bridge's (3.5). Not part of the
  rules but seen: Potentiometer draws ~47 frames a second with nothing computing.
- **Optics (asked 2026-11-02):** SPEC.md §8 plans wave optics. In place: `physics::fdtd2d`
  (Yee, the TE polarisation: E in the plane, B_z; perfect conductors, the CPML border,
  charges as sources), the Poynting tools; mirrors, apertures and slits as conductors work
  now, and a 2D run of ~10⁹ cell updates takes seconds (well inside 30 s). Missing: a
  permittivity ε(x) and its dispersion ε(ω) in the wave solver (lenses, prisms), the other
  polarisation (E_z), a light source element (plane wave or beam of a band), detectors and
  goals on intensity and spot size, free-curve shapes with their editor tools, and the
  validation (Fresnel and Snell, total internal reflection and the evanescent wave, a
  dielectric cylinder's exact scattering, the grid's numerical dispersion, the symmetry and
  energy tests).

## The Transformer (task #61): an interlude, not a level (the owner, 2026-11-02)

Designed as a level with B₀ and ω on sliders, then measured: no accurate version fits a
level's 30 s preview (the bullet *simulation material* below), so it waits on the interludes'
format (SPEC.md §14). The design and what was learned on the way:

The suggested controls (the drive's amplitude and a cut in the secondary) cannot work, and a
level with only the solenoid's amplitude is one slider scaling everything (the owner's rule:
an introduction at most). The owner approved (2026-11-02): the Transformer with two sliders
on the solenoid's supply, its amplitude B₀ and its frequency ω (now in place, DEC-30);
the frequency slider logarithmic; the run time decided after measuring the designed level:

- Physics shown: the secondary screens the primary as `ωL₂/R₂` grows (the field inside it
  `B_ext R₂/(R₂ + iωL₂)`): over `ωL₂/R₂` = 0.3–3 the primary's heat against the secondary's
  falls about tenfold. Goals on both rings' heats (or on the gap field's amplitude and phase):
  the ratio pins ω, the size then pins B₀. Interacting controls, both rings visible on the
  self-field map.
- Parameters to reach `ωL/R` ~ 1 at playable frequencies (ω ≲ 1, well below ω_p and with
  ωτ < 1): c = 2 (L ×2.25 against c = 3), walls 2 with ℓ = 0.5 (R halved against V-D7's walls
  1; ℓ ≤ a quarter of the wall, `metal_issues`): `L₂/R₂` ≈ 2 instead of 0.45.
- The preview's time, to measure and show the owner: V-D7 runs ~4 min per 9 periods at ω =
  0.8 with ~7000 carriers; at ω ≈ 0.5 a few periods is 30–40 time units, ~2 min a preview
  and ~8 a verification unless the carriers per area drop (more noise).
- Scales (`scales.rs`): `self_field_lines` shows the first ring's `ωL/R`; extend it to both
  rings (and `ωτ`).
- **The heat measure for two rings: done (DEC-35).** An isolated metal's heat goal measures
  the drive's work plus the power the other metals' macroscopic currents pass it
  (`TubeState::mutual_work`: window-averaged currents on the self-field's grid, each one's
  vacuum field, `−∫ J̄_k·∂Ā_j/∂t`); one metal: as before, bit for bit. V-D9 (light rings)
  within 5 %. The research (the owner's choice, 2026-11-02) and the measures tried and
  dropped: PHYSICS.md §2.12 *Two metals*, DEC-34/35.
- **The dense rings' dynamics is off (V-D10, recorded, not asserted):** in Self-inductance's
  metal two nested rings come out +20 % (inner) and +0.3 % against the exact heats at the
  levels' step, +9 % and +6 % at half of it, and the drive's work on the inner ring, no heat
  measure in it, is +12 % (+17.5 % at half the step): the outer ring screens too little.
  Tested (2026-11-02, scratch, one seed each): each ring alone within a few per cent, so
  not the σ; coupled, the gap's field 6.7 % high (the outer ring screens too little); a gap
  of four cells no better; the self-field's cell 0.125 (`c/ω_p` 7 cells) brings the inner
  ring to +0.3 % (outer +6.3 %, ~2 seeds' spread). Open: the inner ring's drive's work
  stays +21 % there (the split drive/exchange, not the heat).
- **The Transformer is simulation material, not a level (the owner, 2026-11-02: a level's
  preview at most 30 s, accurate, simpler if need be; otherwise an interlude, SPEC.md §14).**
  Measured: the steps a preview needs are set by the Darwin limit's quality at the outer
  ring's bend (c/ω_b against the rings' size) and the validated accuracy (ω_p dt ≈ 0.04,
  c/ω_p ≥ 3.4 cells), not by the density; the cost per step is the self-field's two solves
  (~7 ms on the 128-node grid every board up to 24 rounds up to; 64 nodes only for a board
  of 8 or less) plus ~1.5 µs per carrier. Exact heats (`prototypes/inductance/
  transformer2d.py`, nested square rings): thin walls (1) screen too little (the bend at
  ω ≈ 2–3, c/ω ≈ 2 against rings 10–12 wide); thick walls (3 in a ring 24 wide, density 1,
  thermal speed 1, c 3.2) work (the inner/outer heat ratio 0.19 → 0.02 over ω 0.03–1.2) but
  preview in ~3–9 min. A small dense prototype (rings 10 and 6 wide, walls 1, density 44,
  c 21, dt 0.0017; `prototypes/proto_transformer.json`): preview 73 s for 12 time units,
  the outer ring's heat ±0.4 % but the inner one's 44 ± 30: at two periods (~26 s) a goal on
  it is noise. As an interlude it can run the thick-walled rings or the dense metal on a
  finer cell (V-D11).
- **Earlier notes for the Transformer:** coupled dense rings need
  the self-field's cell at about `c/ω_p`/7: `SELF_FIELD_CELL` per level (from its density
  and c) or a finer constant, ~4 times the preview (Self-inductance's ~75 s would be ~5 min);
  or a metal whose `c/ω_p` spans more cells at 0.25 (a lower density, with c lowered so
  that c/ω stays about the board). Scratch tests: `prototypes/scratch_cell.rs`
  (copy into `crates/physics/tests/` to run), outputs `prototypes/scratch_*.out`.
- With the self-field the residual above (~0.6 % of the unscreened heat) and DEC-31's
  half-swings apply; Self-inductance's metal (density 60, c 24, dt 0.0072) is a starting
  point: two rings nested need the secondary's `L₂/R₂` near the primary's frequency.

Should a version fit a level after all, a level's checklist stands (CLAUDE.md):
`scripts/designs/` and `curriculum.py`; `scripts/reality.py` (headline: a transformer's mutual inductance;
simplified: z-invariant one-turn coils, slow light); scales; `python scripts/levels.py --only
<slug>`; golden hash with `GOLDEN_ONLY=<slug>`; reference and radiation tests; the level
document (`scripts/level_docs.py new`); `docs/LEVELS.md`; a capture.

## Working notes

- Finished work: `git fetch -q origin`, `git merge-base --is-ancestor origin/main HEAD`, then
  `git push -q origin worktree-tubes-next:main` and `git push -q origin worktree-tubes-next`.
  The main checkout is fast-forwarded (`git merge --ff-only worktree-tubes-next` there) and
  rebuilt (`cargo build --release -p game -p generator` there: the owner plays from it), then
  `python scripts/cleanup.py run` (CLAUDE.md, *Disk*; DEC-37).
- Before committing: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`, the affected heavy tests (`affected.py run`, coverage on: a
  verdict-only run leaves the records' coverage stale and every later change runs all 529);
  `python -m mypy` / `pytest` for scripts; `python scripts/measured.py update`,
  `python scripts/physics_pdf.py` after PHYSICS.md edits.
- Edit scripts written through a Bash heredoc lose doubled backslashes (`\\` became `\`, so a
  replacement matching Rust's string continuations silently missed): use the Edit tool for
  such text, or write the script with the Write tool.
- The Bash guard refuses compound commands with computed values: write small Python scripts in
  `target/scratch` and run them. Windows file locks (Errno 22): retry.
- Never send clicks or keys to the desktop; captures through the game's own `EM_CAPTURE`.
