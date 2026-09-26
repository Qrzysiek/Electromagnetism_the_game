# Project rules

A puzzle game on accurate classical electrodynamics. `SPEC.md` describes the game, `PHYSICS.md` the implemented physics and its validation, `README.md` how to run it.

## Physics accuracy
- Any change to physics or numerics updates `PHYSICS.md` in the same commit: the model, the method, and the measured test results.
- Validation tests compare against analytic or independent results. Set their thresholds from the target accuracy before measuring. If a threshold must change, document why in `PHYSICS.md`; never loosen one just to make a test pass.
- The authoritative trajectory path is f64 on the CPU and deterministic: no fused multiply-add, and a fixed summation order. The GPU is for visuals only.
- Every shipped level must have a verified reference solution for every shot, and neglected radiation below 1e-10 of the launch energy (`crates/level/tests/levels.rs`).
- Golden hashes (`levels/golden_hashes.json`) change only when the physics or the levels change intentionally.

## Sandbox completeness
Every level that ships must be reproducible by hand in the sandbox.
- `crates/game/src/level_editor.rs` edits every field of the level format. Each level type is destructured exhaustively, with no `..`. A new field in the level format therefore fails to compile until it gets a control there, or is explicitly marked as not user-editable with a reason.
- The editor ranges are shared constants. `check_editable` validates each level against them, and the test `all_levels_are_reproducible_in_the_sandbox` runs it on every file in `levels/` and `levels/custom/`.
- A new element kind, coil shape or other level feature needs a sandbox tool or control in the same change.

## Workflow
- Before committing: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`.
- Level design scripts: `scripts/make_experiment_levels.py`, `scripts/make_magnet_levels.py`. Level tools: `cargo run --release -p generator -- {solve,check,normalize}`.
- Commit messages end with the co-author line used in the history.
