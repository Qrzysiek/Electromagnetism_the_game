# Handoff (2026-11-02): where the work stands and what comes next

Read this instead of the previous sessions' transcripts. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; the metals' plan: `docs/DRUDE.md`
"Implementation plan"; decisions: `docs/DECISIONS.md` (DEC-28 to DEC-35 the latest).

## State

`main` = `worktree-tubes-next` (worktree `.claude/worktrees/tubes-next`).

- The heavy-test verification of `4597a78` (the self-field circuits, nested metals): all 521
  verdicts passed, golden hashes included (`scripts/refresh_coverage.py`, 02:50–07:24, 4.6 h
  against the 11.4 h estimated; the coverage of 395 tests recorded with it).
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
- Not done from DRUDE.md Phase 2 item 3: contacts on a ring's inner side (an electrode in the
  hole is refused: its lead would cross the ring; the Corbino disc needs one). In the
  electrostatic model a lead is only drawn, so it could be allowed there with the lead drawn
  across (refused with the self-field, as a lead from a ring is, DEC-27).
- Coverage of the heavy tests not yet recorded since the rings: `python scripts/affected.py
  refresh` (hours; or `scripts\refresh_coverage.cmd`) before the next engine change, else that
  change selects all 521 again.
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

## Done: Self-inductance (level 125, DEC-32), the self-field's heat measure fixed

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
  fixed, ω on a log slider 0.05–5, heat 1000–1500 over t = 2.5–10 (to 20 since 2026-11-04) at the bend (answer ω 1.3;
  verified 1195 ± 136, the preview 1105; ω solves on 1.15–1.55 verified, 6.5 %). Preview
  24.4 s, verification 62 s. (The first design: density 60, c 24, the bend at 0.7, ~75 s.)
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
  by at most 100 % where that makes the level a lot better. To audit: every tube level's
  preview time and streaming rate; LEVELS.md already records Eddy currents in silicon (118)
  at 38 s and Lock-in Hall (127) at 30–70 s. Self-inductance (125) streamed at 0.41 units a
  second. The owner (2026-11-04): fewer carriers are fine if the physics barely suffers.
  Measured first: they are not the cost (at a tenth of them the preview took 20 s against
  27); 70 % of a step was the self-field's two solves, 16–17 conjugate-gradient iterations
  each, whatever the carriers. Not changed: the carriers. Then, the owner's choice of more
  solver work: the iteration 0.47 → 0.075 ms (DEC-13's history, PHYSICS.md §2.12), the
  scattering and the frames' middle field in parallel: it streams 1.45 units a second
  (preview 13.8 s to t = 20, the goal's window lengthened for the verdict's noise, DEC-32).
  The owner (2026-11-04): the animation stood still while the verification ran. Fixed
  (DEC-39): the verification now runs beside the rest of the display run's stream, and a
  run whole at the goal's end is marked done first, so the playback loops; `EM_TIMING=1`
  logs the frame rate and the time shown. The audit of the other tube levels (`generator
  check` with a fresh run cache, `target/scratch/timing_audit.py`; the field maps' rate:
  `field_audit.py`, not run yet) has 99–112 so far: previews 0.25–17.8 s, all streaming
  far above a unit a second; 113–127 still to measure. Simplify what is over, or move it to
  the interludes.
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

## Then the Transformer (task #61), with B₀ and ω on sliders (the owner, 2026-11-02)

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

The rest of the Transformer's checklist stands (CLAUDE.md): `scripts/designs/` and
`curriculum.py`; `scripts/reality.py` (headline: a transformer's mutual inductance;
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
  verdict-only run leaves the records' coverage stale and every later change runs all 525);
  `python -m mypy` / `pytest` for scripts; `python scripts/measured.py update`,
  `python scripts/physics_pdf.py` after PHYSICS.md edits.
- Edit scripts written through a Bash heredoc lose doubled backslashes (`\\` became `\`, so a
  replacement matching Rust's string continuations silently missed): use the Edit tool for
  such text, or write the script with the Write tool.
- The Bash guard refuses compound commands with computed values: write small Python scripts in
  `target/scratch` and run them. Windows file locks (Errno 22): retry.
- Never send clicks or keys to the desktop; captures through the game's own `EM_CAPTURE`.
