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
cargo run --release -p generator -- solve levels/05_the_wall.json
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
  - S flips the sign (charges), the orientation (magnets: ⊙ moment out of the plane, ⊗ into it) or the phase (antennas).
  - R rotates an antenna by 45° (Shift+R: back); W changes its frequency where the level offers several.
  - Q/E change the magnitude; M cycles the element kind (charge, magnet, antenna) the level allows; C clears; 1–4 set the grid refinement.
  - Multi-shot levels: `[` / `]` switch the shot; H shows all shots.
  - N/P switch levels; V cycles the map (potential → magnetic B → off); F toggles field lines; A toggles the animation.

## Multi-shot levels

A level can have several shots: particles with their own species, launch and detector. **One setup must deliver every shot.** The panel lists the shots with their status (✔ arrived, ✖ lost, ⚠ marginal, … computing).

## Antennas

Antennas are small oscillating electric dipoles in the plane. They run at the level's RF generator frequency, or at a frequency you choose where the level offers a choice (W cycles it, Shift+W back; the sandbox can give any antenna its own ω). All start in phase at t = 0. Their fields are exact, including the induction and radiation terms (PHYSICS.md §2.4). The kick they give depends on when a particle passes, so shots launched at different times (shown in the shot's launch time) can be sorted: RF separators, streak cameras. The potential map and field lines show static sources only.

## Radiation

- **Maps:** V cycles potential → magnetic B → waves (antennas and plane waves at the animation time) → particle field (the particle's own Liénard–Wiechert field, optionally only its radiation part) → total (everything at once) → off.
- The time-dependent maps have optional E arrows, a colour choice (B_z or |E|) and an adjustable dynamic range. Raise the range to see weak contributions such as the particle's own field next to strong electrodes.
- **Radiation reaction:** levels can include it (chapter 8). The particle then loses the energy it radiates, and the energy bars show kinetic + potential + radiated = constant.

## Disturbances

Some levels have fields from outside the arena: a stray field that is sometimes switched on, mains hum at different phases, the Earth's field for different orientations. Every shot is flown under every listed disturbance, and the setup must work for all of them. The panel lists the disturbances with their status. Selecting one shows its flight in detail and draws it brightest; the other flights of the shot are drawn fainter, as a bundle. The potential map and the field lines show the level's own sources only. In the sandbox, disturbances are edited in the "Disturbances" section.
