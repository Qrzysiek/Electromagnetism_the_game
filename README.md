# Electromagnetism – the game

A puzzle game built on accurate classical electrodynamics. Place charges so that a test particle flies from the launch point into the detector.

- `SPEC.md`: what the game is and the development stages.
- `PHYSICS.md`: the implemented physics, numerics and measured validation results.

## Running

Requires Rust (stable) and, on Windows, the Visual Studio C++ build tools.

```
cargo run --release -p game          # the game (run from the repository root)
cargo test --workspace               # all tests, including physics validation
cargo run --release -p generator -- check levels/01_first_bend.json
cargo run --release -p generator -- solve levels/02_the_wall.json
```

## Controls

- Mouse: left click places a charge, right click removes it, the wheel changes the magnitude.
- Keyboard:
  - Arrows move the cursor (Shift: ×5); Space/Enter place; Del/X remove.
  - S flips the sign; Q/E change the magnitude; C clears; 1–4 set the grid refinement.
  - N/P switch levels; V toggles the potential map; F toggles field lines; A toggles the animation.
