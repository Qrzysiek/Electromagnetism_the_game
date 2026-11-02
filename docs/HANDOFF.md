# Handoff (2026-11-02): where the work stands and what comes next

Read this instead of the previous sessions' transcripts. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; the metals' plan: `docs/DRUDE.md`
"Implementation plan"; decisions: `docs/DECISIONS.md` (DEC-28 to DEC-30 the latest).

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

## Next: first a self-inductance introduction (the owner, 2026-11-02)

The curriculum rule wants an introduction before the Transformer (the first level with the
carriers' own field and with a frequency slider). The owner chose: one ring alone in the
swinging solenoid, B₀ fixed, ω on the slider; its heat grows as ω² until `ωL/R ≈ 1`, where
its own field screens the drive and the heat levels off; a heat goal pins ω.

- Needs a supply that sets only a frequency (its target's level list empty: switched off the
  solenoid's field stays `b_z` = 0, DEC-12; on, its sine's own amplitude at the player's ω),
  a small extension of DEC-30: placement, the solver's candidates, `analyze`'s steps, the
  panel's row (no level slider), the solenoid's wave in `tube/build.rs`.
- First scans (ring 20 × 20, walls 2, ℓ 0.5, c = 2, B₀ 0.5, heat over t = 10–40; the
  worktree's `target/scratch/inductance/`): heat 0.11, 0.64, 1.50, 2.42, 6.5, 22.5, 73 at
  ω = 0.05–3.2 (doubling). Without the self-field (c = 100) 0.92, 13.6, 219 at 0.1, 0.4, 1.6
  (as ω²). The ratio 0.70 at ω = 0.1 is the loop's `1/(1 + (ωL/R)²)` with `L/R` = 6.55
  (estimated `4πA/c² · σw/P` = 6.6); at high ω ~10 % of the unscreened heat remains, the
  wall's own eddy currents (skin depth 1.3 at ω = 1.6, under the wall's 2). Next: walls 1
  with ℓ 0.25 (that share ~4× smaller), c ≈ 1.2 to keep the bend near ω ≈ 0.2. About 4 min
  per run of 40 time units alone.
- How real: at c = 2 the Darwin limit's own condition (the board small against c/ω) holds
  only below ω ≈ 0.08; the whole bend lies beyond it (slow light, as V-D7).

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
- Scales (`scales.rs`): with a frequency slider, show `ωτ` and the rings' `ωL/R` for the
  current setting.

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
