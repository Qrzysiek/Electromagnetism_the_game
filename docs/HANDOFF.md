# Handoff (2026-11-02): where the work stands and what comes next

Read this instead of the previous sessions' transcripts. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; the metals' plan: `docs/DRUDE.md`
"Implementation plan"; decisions: `docs/DECISIONS.md` (DEC-28 to DEC-33 the latest).

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
  the test's draft: `target/scratch/ac_mobility_test.rs` (this worktree). Phase goals at
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
- The level: a ring 20 × 20 (walls 2), density 60, c = 24, dt 0.0072, B₀ 0.45 fixed, ω on a
  log slider 0.03–3, heat 360–650 over t = 4–24 at the bend (answer ω 0.67, ωL/R 0.93;
  verified 490 ± 46, the preview 521; ω solves on 0.58–0.8 verified, 7 %). Preview ~75 s,
  verification ~185 s. Scales:
  `ωL/R`, the skin depth, the board against c/ω (`self_field_lines`). The exact curve:
  `target/scratch/inductance/eddy2d.py` (this worktree); the scans: `NOTES.md` there.
- Checked: the panel's "Joule heat" (then the lattice's take) read ~750 per time above the
  drive's work on this level: the macroparticles' stirring of each other (`∝ n²/N`; with
  the drive at rest 748 without the self-field, 742 with it, 400 at twice the carriers),
  not the self-field. The panel now shows the field's work, what the goal measures (DEC-33;
  PHYSICS.md §2.12 *The heat at rest*; the runs: `target/scratch/joule/`).
- Not written: the level's document (`scripts/level_docs.py new self_inductance`; none of the
  metal or AC levels has one yet); `eddy2d.py` would be its `calc.py`.

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
- **The heat measure for two rings**: an isolated metal's heat is the drive's induced field's
  work (`induced_work`), which for one ring is its Joule heat on average (its own field's
  work averages out). A secondary also takes power from the primary through the self-field
  (mutual induction), which `induced_work` misses. Book the self-field's part of the step's
  work per metal too (`½m(v_end² − v_push²)` in `darwin_step`) and measure the total, or the
  lattice's heat (noisy).
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
  The owner pulls and builds `cargo build --release -p game -p generator`.
- Before committing: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`, the affected heavy tests (`affected.py run --no-coverage`, then
  `refresh`); `python -m mypy` / `pytest` for scripts; `python scripts/measured.py update`,
  `python scripts/physics_pdf.py` after PHYSICS.md edits.
- Edit scripts written through a Bash heredoc lose doubled backslashes (`\\` became `\`, so a
  replacement matching Rust's string continuations silently missed): use the Edit tool for
  such text, or write the script with the Write tool.
- The Bash guard refuses compound commands with computed values: write small Python scripts in
  `target/scratch` and run them. Windows file locks (Errno 22): retry.
- Never send clicks or keys to the desktop; captures through the game's own `EM_CAPTURE`.
