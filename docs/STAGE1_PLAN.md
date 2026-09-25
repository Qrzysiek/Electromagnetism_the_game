# Stage 1 plan: vacuum, single particle

Approved 2026-09-25. Review checkpoints with the project owner: **M1, M2, M4**.

## Crates

| Crate | Role |
|---|---|
| `physics` | Pure library. Units, `FieldSolver` (Coulomb; uniform field for tests), relativistic dynamics in `(x, p)`, DOP853 integrator wrapper, events, trajectory, verification, diagnostics. Validation tests T1–T11 in `tests/`. |
| `level` | Level format (serde JSON), integer grid indices, limits, engine version. |
| `generator` | CLI: random levels, simulated-annealing search, triviality / multiple solutions / robustness checks, verification, difficulty metric. |
| `game` | Bevy app (pinned version) + bevy_egui: 2D view, editor, preview worker thread, potential map shader, CPU field lines. |

## Design decisions
- Charge positions are integer grid indices. Refining by k multiplies the indices, so old nodes stay exact.
- `FieldSolver` returns `(E, φ)` (later B) at `(x, t)`. Obstacle geometry is a separate trait.
- Particles are stored as structure of arrays from day one.
- 2D mode: everything at z = 0, pz = 0. z stays exactly 0 (test T8 checks this bitwise).
- The preview worker drops stale requests and computes the preview first, then the verified run. Rendering never blocks.
- Visuals only: the potential map uses an f32 GPU shader. Field lines use f64 on the CPU, in parallel.

## Milestones
| # | Content | Checkpoint |
|---|---|---|
| M0 | Toolchain, workspace, git, fmt/clippy config, CI (Windows + Linux) | |
| M1 | Integrator selection: benchmark `russell_ode`, `ode_solvers` (fallback: own DOP853) on Hairer reference problems + T6. Record the decision in PHYSICS.md | **yes** |
| M2 | Physics core + tests T1–T10 with measured errors in PHYSICS.md; 50-charge preview under 16 ms | **yes** |
| M3 | Verification (verified/marginal), level format, save/load, T11 determinism with golden hashes across OSes | |
| M4 | Playable 2D game with hand-made levels | **yes** |
| M5 | Generator | |
| M6 | Stage 1 acceptance check (SPEC §10), PHYSICS.md statuses | |

## Testing
- Validation tests print their measured errors. The numbers are recorded in PHYSICS.md.
- `proptest` for the speed limit, symmetry and grazing events. `criterion` for the preview budget and generator throughput.
- Editor logic is plain Rust, separate from Bevy, and unit tested.
- CI compares trajectory hashes between Windows and Linux.
