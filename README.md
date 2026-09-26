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

## Sandbox (level editor)

Switch to **Sandbox** at the top of the panel.

- **Tools** (click on the map):
  - level charges (any value; right click removes one)
  - player charges (to test solutions against the limits)
  - launch point
  - aim (click where the launch direction should point)
  - detector (drag a box)
- **Panel:** name, grid, c, particle, launch energy and angle, charge radius, time limit and player limits.
- **Check solvability** runs the solver in the background. You can then store the solution it found, or your own charges, as the reference solution.
- **Save** writes to `levels/custom/`. The game loads those levels too. On saving, the panel warns if the reference solution is not verified, or if neglected radiation exceeds 10⁻¹⁰ of the launch energy.

## Controls

- Mouse: left click places a charge, right click removes it, the wheel changes the magnitude.
- Keyboard:
  - Arrows move the cursor (Shift: ×5); Space/Enter place; Del/X remove.
  - S flips the sign (charges) or orientation (magnets: ⊙ moment out of the plane, ⊗ into it).
  - Q/E change the magnitude; M switches between charge and magnet; C clears; 1–4 set the grid refinement.
  - Multi-shot levels: `[` / `]` switch the shot; H shows all shots.
  - N/P switch levels; V cycles the map (potential → magnetic B → off); F toggles field lines; A toggles the animation.

## Multi-shot levels

A level can have several shots: particles with their own species, launch and detector. **One setup must deliver every shot.** The panel lists the shots with their status (✔ arrived, ✗ lost, ⚠ marginal, … computing).
