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
- **Computational cost** (bars under the tools): how hard the setup is on the computer, measured while it is computed.
  - Meters: preview time (the delay after every change), verification time, steps of the longest flight, and with metal or electrodes the one-off setup time, matrix memory and the metal model's error.
  - Green is within budget, amber is high, red is over the limit (the marks on each bar). Hover a meter for what drives it.
  - Budgets and measurement: `crates/level/src/cost.rs`. `generator check` prints the same meters.
- **Check solvability** runs the solver in the background. You can then store the solution it found, or your own charges, as the reference solution.
- **Save** writes to `levels/custom/`. The game loads those levels too. On saving, the panel warns if the reference solution is not verified, if neglected radiation exceeds 10⁻¹⁰ of the launch energy, or if a cost meter is over its limit.

## Controls

- Mouse: left click places a charge, right click removes it, the wheel changes the magnitude.
- Keyboard:
  - Arrows move the cursor (Shift: ×5); Space/Enter place; Del/X remove.
  - S flips the sign (charges), the orientation (magnets: ⊙ moment out of the plane, ⊗ into it) or the phase (antennas).
  - R rotates an antenna by 45° (Shift+R: back); W changes its frequency where the level offers several.
  - Moving elements: drag them with the mouse, or G to grab / drop with the arrow keys in between, Esc to cancel. In the sandbox, level elements can be dragged too.
  - **Hardcore** (checkbox next to "Your elements"): sliders instead of the level's fixed values. Any magnitude or frequency between the smallest and largest listed value, and any antenna angle. Each slider has its own linear/log switch. With the cursor on one of your elements the sliders edit it live; Q/E and W then step by ×1.1.
  - Q/E change the magnitude; M cycles the element kind (charge, magnet, antenna) the level allows; C clears; 1–4 set the grid refinement.
  - Multi-shot levels: `[` / `]` switch the shot; H shows all shots.
  - N/P switch levels; V cycles the map (potential → magnetic B → off); F toggles field lines; A toggles the animation.

## Multi-shot levels

A level can have several shots: particles with their own species, launch and detector. **One setup must deliver every shot.** The panel lists the shots with their status (✔ arrived, ✖ lost, ⚠ marginal, … computing).

## Detectors with conditions

Some detectors also require how the particle arrives: moving in a direction (drawn as a cone at the detector) and/or with a kinetic energy in a window, as the entrance of a next stage does. A particle entering outside these conditions is "rejected". The panel shows the allowed values and how the particle arrived. In the sandbox, set them in the shot's "Accept direction" and "Accept energy" rows.

## Electrodes

Levels can contain real electrodes: metal plates, slabs and walls with finite size, grounded or held at a voltage. Their field, including the fringe field at their ends, is computed by the boundary element method (PHYSICS.md §2.7). Field lines end on them, and touching one loses the particle. In the sandbox they are edited in the "Electrodes" section.

Some levels let you build or tune electrodes yourself:
- **Plates:** select "plate" in the palette and click to place one. R turns it (along x or y), Q/E or the wheel changes its potential, and it can be dragged by its body. Plates keep 1 cell away from other metal and stay clear of elements, coils and detectors; an outline at the cursor shows whether a plate fits there.
- **Power supplies:** electrodes with a yellow frame are tunable. Click one to switch its supply on or to the next potential (S: the opposite one, right click: off), or choose in the "Power supplies" panel.
- The flight details show a bound on the electrodes' neglected image force, with a warning if a flight passes so close to metal that it is no longer negligible.
- In the sandbox: "Player plates", "Plate potentials", "Plate L×T×H" and "Supply voltages" in the limits, "tunable" on each electrode.

## Gates (multi-stage instruments)

Some levels have gates: violet dashed boxes, numbered by dots, that every flight must pass in order before its detector counts, sometimes with a required direction (a cone) or energy. They are the stages of an instrument, for example "first make the beam parallel, then bend it into the experiment". A flight that reaches the detector without passing a gate "entered the detector without passing gate N". In the sandbox, gates are edited in the "Gates" section (PHYSICS.md §6.2).

## Beams

Some levels fire a whole beam: many particles with a spread in position, direction and energy (Gaussian or uniform), always the same fixed sample so the result is verifiable. The goal is a share of the beam, for example 90 %, arriving verified. In Newtonian levels the particles can repel each other exactly (space charge), so the beam spreads out on its way. Every particle's path is drawn in its shot's colour; lost particles are fainter. The panel shows each beam shot's verified transmission against its requirement. In the sandbox a shot becomes a beam with "fire as a beam" (count, spreads, distribution, transmission, seed), and "Beam interaction" switches the repulsion on (PHYSICS.md §3.3).

## Magnetic moments (spin)

Particles can carry a magnetic moment perpendicular to the plane: spin up or spin down. Neutral atoms with a moment feel no Lorentz force but the force m ∇B_z, which pushes one spin state towards stronger field and the other away: the Stern–Gerlach experiment (PHYSICS.md §3.2). The potential map then shows the particle's magnetic energy −m B_z, and the magnetic map is scaled to it. In the sandbox the moment is set per shot ("Magnetic moment m_z").

## Metal spheres

Levels can contain metal spheres. They can be grounded, isolated with a net charge, or held at a potential (like a Van de Graaff dome).
- Their surface is an equipotential, and charges nearby induce opposite charges on them. So a charge placed next to a grounded sphere is partly cancelled by its "image".
- Field lines end on the metal at right angles, and the potential map includes the induced charges.
- Touching a sphere loses the particle.
- In the sandbox they are edited in the "Metal spheres" section.

## Antennas

Antennas are small oscillating electric dipoles in the plane. They run at the level's RF generator frequency, or at a frequency you choose where the level offers a choice (W cycles it, Shift+W back; the sandbox can give any antenna its own ω). All start in phase at t = 0. Their fields are exact, including the induction and radiation terms (PHYSICS.md §2.4). The kick they give depends on when a particle passes, so shots launched at different times (shown in the shot's launch time) can be sorted: RF separators, streak cameras. The potential map and field lines show static sources only.

## Radiation

- **Maps:** V cycles potential → magnetic B → waves (antennas and plane waves at the animation time) → particle field (the particle's own Liénard–Wiechert field, optionally only its radiation part) → total (everything at once) → off.
- The time-dependent maps have optional E arrows, a colour choice (B_z or |E|) and an adjustable dynamic range. Raise the range to see weak contributions such as the particle's own field next to strong electrodes.
- **Radiation reaction:** levels can include it (chapter 8). The particle then loses the energy it radiates, and the energy bars show kinetic + potential + radiated = constant.

## Disturbances

Some levels have fields from outside the arena: a stray field that is sometimes switched on, mains hum at different phases, the Earth's field for different orientations. Every shot is flown under every listed disturbance, and the setup must work for all of them. The panel lists the disturbances with their status. Selecting one shows its flight in detail and draws it brightest; the other flights of the shot are drawn fainter, as a bundle. The potential map and the field lines show the level's own sources only. In the sandbox, disturbances are edited in the "Disturbances" section.
