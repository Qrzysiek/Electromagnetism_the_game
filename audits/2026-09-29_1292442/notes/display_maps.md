# Display maps and the flight

The flight samples `LevelField` in f64. Two GPU paths draw a field, and they do not sample the same object.

## Potential and magnetic maps

`update_map` (`crates/game/src/main.rs`) builds `Level::display_scenario(active_shot, placement)`. That function always passes disturbance index 0. The redraw key does not contain `active_disturbance`.

`potential::params` then uploads:

- every Coulomb source as `(x, y, Q, charge_radius)`, clouds included, through `Coulomb::charges()`
- induced sphere charges, and each electrode panel as a point charge at its centroid and at its mirror image
- dipole `m_z / b_ref`, loop `κ / b_ref`, segment `κ / b_ref`

`b_ref` is `|p₀| / (|q| · 5)` for a charged particle, or `T₀ / |m|` for a neutral moment. Colour 1 on the magnetic map is that reference field.

`potential.wgsl` evaluates `Σ Q / max(√(d²+z²), 1e-4)`. For the first `params.solid` entries, `z` is forced to 0 and `w` is a solid radius. Clouds are in that prefix, so a cloud centre becomes a disk of `charge_radius` and the interior potential is `Q/r`. There is no slot for uniform `E` or uniform `B_z`.

The launch offset `u_a` and the dark-region limits *do* call `FieldSolver::sample`, which has the cloud interior and the stray fields. The shader then forms

```text
u = (q / T) φ_shader + (−m b_ref / T) B_shader − u_a
g = (q / T) φ_shader + (−m b_ref / T) B_shader − E_cpu / T
```

`u_a` and `E_cpu` use the trajectory’s potential. `φ_shader` does not. Where those potentials differ, the colour at the launch point is not zero and the turning line is not the energy boundary.

On level 57 with no placed charge, `φ_shader = 0` and `B_shader = 0`. The script’s result: colour `±10000 T` and a forbidden excess of about `±10000`, against a true energy boundary `0.001` cells from the launch line. `B_z / b_ref = 116.7`.

On level 76 the ion starts at `d = 2` inside a cloud of radius 4. The trajectory turns at `2.828`. The map turns at `3.200`. The colour at the launch point is `−5 T`.

## Field lines

`update_field_lines` integrates `sample(x, 0).e` on `display_scenario(0, …)`. Uniform `E` of disturbance 0 is included, because it is inside `sample`. Uniform `E` of any later disturbance is not. The lines are electric; they never draw `B`.

If the scenario has no Coulomb sources, no metal spheres, and no electrodes, the function returns an empty list before the grid sweep. Level 57 before the player places a charge is that case, so the crossed field draws no lines. After a charge is placed, the lines include disturbance 0’s `E` and the point-charge field of the electrodes if any, at display resolution.

## Total, waves, and particle field

`radiation.rs` builds `scenarios_at(…, Display)` and keeps the active flight, so the selected disturbance is the one on screen.

The total view rasters `static_part(field).sample` onto a grid and uploads f32. `static_part` drops oscillating antennas and waves with `ω ≠ 0`. It keeps uniform stray fields, `ω = 0` waves, clouds, magnets, coils, and metal. Cloud interiors and `−E·x` are therefore in this view.

Oscillating waves are uploaded separately and evaluated in `radiation.wgsl` as `E = E₀ ê cos(phase)`, `B_z = (k̂ × E)_z / c` when the uploaded `c` is positive. `c = ∞` is uploaded as 0, and both the antenna and the wave branches then skip retardation and `B`. That matches `PlaneWave::fields` and `antenna.rs` at infinite `c`.

Waves with `ω = 0` are not in that oscillating upload. They remain in `static_part`, so the total view still samples them, including `B = k̂ × E / c` at finite `c`.

The particle-field view draws the Liénard–Wiechert field of the stored world line. Acceleration on a single-particle preview is rebuilt from `PathPoint.force`, which omits the image force and the radiation reaction. See FINDINGS §7.
