# Audit of Electromagnetism, the game

Read-only review of the Rust workspace, done on 2026-09-29. Nothing under `crates/`, `levels/`, `docs/`, `scripts/`, or `PHYSICS.md` was edited. The notes here are evidence and transcriptions, not a patch.

The ranked results are in [FINDINGS.md](FINDINGS.md). What each one does on a running level is in [IMPACT.md](IMPACT.md), with screenshots in `screenshots/`. The per-level arithmetic is in [results/checks.txt](results/checks.txt). A pattern count (`unwrap`, `expect`, `unsafe`, `mul_add`, `f32`, `todo!`) is in [results/static_scan.txt](results/static_scan.txt). Shorter write-ups:

- [FORMULAS.md](FORMULAS.md) — each engine formula against an independent reference, and whether the docs and the legends state its limit
- [REMAINDER.md](REMAINDER.md) — the areas that pass left open: integrator, panels, beams, geometry, tests, and the release-game launch
- [notes/display_maps.md](notes/display_maps.md) — why the potential and magnetic maps disagree with the flight
- [notes/what_held.md](notes/what_held.md) — the earlier short list of formulas and checks that matched

## How to reproduce the numbers

From the repository root, in PowerShell:

```powershell
py -3 audit/scripts/audit_checks.py
py -3 audit/scripts/static_scan.py
wolframscript -file audit/scripts/wolfram/checks.wls
wolframscript -file audit/scripts/wolfram/formulas.wls
py -3 audit/scripts/remainder_checks.py
py -3 audit/scripts/closing_checks.py
py -3 audit/scripts/ll_drain_check.py
py -3 audit/scripts/difficulty_diff.py
wolframscript -file audit/scripts/wolfram/closing.wls
```

The Wolfram scripts need Wolfram Engine 14.2 (`wolframscript`), one kernel at a time. `checks.wls` is kept in [results/wolfram.txt](results/wolfram.txt): the dipole identity, the uniform-sphere field, `β(ε)`, the taper integral, the coil brackets, and the cup illustration. `formulas.wls` is the formula audit; its transcript is [results/formulas_wolfram.txt](results/formulas_wolfram.txt). `formulas_rest.wls` is an earlier draft: its Landau–Lifshitz line printed 0.024 because it used one `dB` for both unit systems. That number is a script error. The corrected residual is 0.

The pictures in [IMPACT.md](IMPACT.md) come from the release game, run from the repository root. One example, the crossed-field potential map:

```powershell
$env:EM_CAPTURE = "audit\screenshots\57_potential.png"
$env:EM_LEVEL = "57"
$env:EM_MAP = "potential"
$env:EM_WAIT = "1"
& "target\release\game.exe"
```

The window opens off-screen and exits. `generator trace audit\scenarios\<name>.json` integrates the cases that are not shipped levels.

`audit_checks.py` reads `levels/[0-9]*.json` and applies the expressions from `crates/physics/src/field.rs` (uniform-sphere potential), `crates/game/src/potential.wgsl` (`Q / max(r, 1e-4)`), `crates/physics/src/external.rs` (`PlaneWave::fields`), `crates/physics/src/magnetic.rs` (`elliptic_ke`, `loop_bracket`), `Level::box_region` with `Aabb::signed_distance`, `β` from the kinetic energy, `beam::tapered`, and the Coulomb force of a fading cup charge. It does not link the crate. `static_scan.py` counts substrings in `crates/{physics,level,game,generator}/src`. Comments are included, so the counts are a scan, not a verdict.

## What was read

The force assembly (`dynamics.rs`, `beam.rs::force_in`), the cloud and external samples (`field.rs`, `external.rs`), the event search (`events.rs`, the launch test and the acceptance block in `trajectory.rs`, the same split in `beam.rs`), the Landau–Lifshitz expression against the formula in its own comment, the coil bracket against a direct elliptic evaluation, the panel jump test, the map upload (`potential.rs`, `potential.wgsl`, `main.rs`), the total-field view (`radiation.rs`), `Level::model_issues`, the level editor’s exhaustive destructure, `beam::tapered` against its integral, the cup booking against `budget` and the energy panel in `ui.rs`, and the antenna potential. `PHYSICS.md` §§2–4 and §6 were compared with that code. All 88 shipped level files were scanned by the script: launch nodes against detectors and gates, and launch `β`.

## What was not done

The follow-up is [REMAINDER.md](REMAINDER.md): the integrator against Hairer’s source, the panels, the sphere fit, the retarded-time search, the geometry, the workspace tests, one shipped flight, the seven fast launches, the charge-cap pictures, and the formula size against Ubuntu Light. Two debug beam tests panic; that log is part of the note. The cup number in `checks.txt` integrates the fade force on one invented survivor; it is not level 46. Where a statement depends on a test result printed in `PHYSICS.md` and not re-run in `REMAINDER.md`, this audit treats that table as the project’s own measurement.
