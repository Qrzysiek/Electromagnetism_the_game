# Handoff (2026-10-28, afternoon): where the work stands and what comes next

Read this instead of the previous sessions' transcripts. Project rules: `CLAUDE.md`; the
self-field plan and its log: `docs/DARWIN.md` §4–5; the metals' plan: `docs/DRUDE.md`
"Implementation plan"; decisions: `docs/DECISIONS.md` (DEC-27, DEC-28 the latest).

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
- Coverage of the commit's heavy tests not yet recorded: `python scripts/affected.py refresh`
  (hours; or `scripts\refresh_coverage.cmd`) before the next engine change, else that change
  selects all 521 again.

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

## Next: the Transformer (task #61) — a decision for the owner

The handoff's suggested controls (the drive's amplitude and a cut in the secondary) cannot
work, and a level with only the solenoid's amplitude is one slider scaling everything (the
owner's rule: an introduction at most). The other player controls of a swinging-field level
cannot act on the induction: plates and supplies are refused beside an isolated metal (every
electrode must float, OP-1), free charges do not couple.

Recommended: **first the drive's frequency on a slider** (DRUDE.md Phase 2 item 2, already
approved; 1–2 days: the format, the editor and its ranges, `generator window` over it, the
AC conductivity test through the level path against D2), then the Transformer with two
sliders, the solenoid's amplitude B₀ and its frequency ω:

- Physics shown: the secondary screens the primary as `ωL₂/R₂` grows (the field inside it
  `B_ext R₂/(R₂ + iωL₂)`): over `ωL₂/R₂` = 0.3–3 the primary's heat against the secondary's
  falls about tenfold. Goals on both rings' heats (or on the gap field's amplitude and phase):
  the ratio pins ω, the size then pins B₀. Interacting controls, both rings visible on the
  self-field map.
- Parameters to reach `ωL/R` ~ 1 at playable frequencies (ω ≲ 1, well below ω_p and with
  ωτ < 1): c = 2 (L ×2.25 against c = 3), walls 2 with ℓ = 0.5 (R halved against V-D7's walls
  1; ℓ ≤ a quarter of the wall, `metal_issues`): `L₂/R₂` ≈ 2 instead of 0.45.
- To ask the owner with it: the slider rule for a frequency (the potentials' "symmetric about
  0" does not apply; a range [ω_min, ω_max] with the answer clear of its ends and middle?),
  and the preview's time: V-D7 runs ~4 min per 9 periods at ω = 0.8 with ~7000 carriers; at
  ω ≈ 0.5 a few periods is 30–40 time units, ~2 min a preview and ~8 a verification unless the
  carriers per area drop (more noise).

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
