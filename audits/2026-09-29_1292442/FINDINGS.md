# Findings

The flights are integrated on the CPU in f64 from `LevelField::sample`. That sample uses the uniform-sphere potential inside a charge cloud and adds uniform stray fields. The reference solutions ride on that path. The defects below are places where a picture, a gate, or a sentence does not.

What each one does to a running level is in [IMPACT.md](IMPACT.md): screenshots from the game's own capture, and `generator trace` on the cases no shipped file hits. The formula-by-formula comparison with independent references, for the CPU solvers and the GPU shaders, and whether the docs and the legends state the limits, is in [FORMULAS.md](FORMULAS.md). The areas left open by that comparison are checked in [REMAINDER.md](REMAINDER.md).

Shipped levels are affected by the two display defects and by the beam energy total on the five Newtonian interacting beams (finding 12). The moment gate, the static antenna, and a launch already inside its detector are sandbox issues: no shipped file hits them. The tapered past can exceed `c`; every shipped beam launch is still slow enough that adding `0.1c` stays under `c` (`results/checks.txt`).

## 1. Uniform stray fields are missing from the potential and magnetic maps

**Shipped. The trajectory is unaffected. The pictures the levels point at are not.**

`potential.wgsl` sums `Q/r` and the `B_z` of dipoles, circular coils, and polygon segments (`potential_colour`, `magnetic_field`). `potential::params` fills those arrays from `Coulomb::charges()`, the dipoles, the loops, and the segments. It never reads `LevelField::external`.

Uniform `E` has potential `−E·x` and uniform `B_z` is stored on the flight (`external.rs`, `External::Uniform`; added in `LevelField::sample`). The colour’s zero and the energy limit *are* taken from that full sample (`potential.rs`, `u_a` and `limits`), then combined in the shader with a potential that lacks `−E·x`. On an empty arena the shader potential is 0, so the map is one flat colour, shifted by the launch value of `−q E·x / T`.

The maps are also stuck on disturbance 0. `Level::display_scenario` always calls `scenario_with(shot, 0, …)`, and the redraw key in `update_map` is `(sent_revision, active_shot, show_all_shots, map, has_flight)`. Selecting another disturbance changes the flight and leaves the map. Field lines use the same `display_scenario(0, …)` (`main.rs`). With no charges, no metal, and no electrodes, `field_lines_cancellable` returns immediately, so a pure uniform `E` draws no lines until something else is placed.

The legend calls the potential map exact and says the dark region is the energetically forbidden set (`ui.rs`). `PHYSICS.md` §2.1 and §3.2 describe that map as the particle’s `U = q(φ − φ_A) − m(B_z − B_{z,A})`. The magnetic-map hover is narrower (“magnets and coils”); the legend still defines colour 1 as the field in which the particle gyrates at 5 cells, which is the total `B_z`.

The total-field view is the one that matches the flight. It rasters `static_part` of the active flight’s `LevelField` with `FieldSolver::sample` (`radiation.rs`), active disturbance included, clouds included.

### Level 57, Jackson §12.3, E×B drift

The description tells the player that the guiding centres follow the equipotentials and to look at the potential map. The only sources are `E = (0, 10^5)` and `B_z = 3.3×10^5`. There are no elements until the player places a charge. `c = ∞`.

For either ion (`q = ±10^{-6}`, `T = 10^{-4}`, launch at `(2, 10)`):

| quantity | value |
|---|---|
| shader potential | 0 |
| colour at launch | `±10000` times `T` (the true relative potential is 0) |
| forbidden excess `g` | `+9999` for the positive ion, `−10001` for the negative ion |
| true distance to the energy boundary, against `E` | `0.001` cells |
| `B_z / b_ref` | `116.7` (the magnetic map adds 0) |

`b_ref = |p| / (|q| · 5)` is the field of a 5-cell gyroradius. A value of 117 means the orbit is about `5/117 ≈ 0.043` cells across. The magnetic map is the zero field.

With both shots shown, the dark region is the minimum of the two excesses, so it disappears. Each shot on its own is either entirely dark or entirely allowed. The true allowed set for one sign is a half-plane starting `0.001` cells off the launch line. Placing the reference charge adds that charge’s `Q/r` and still omits `−E·y`, so the equipotentials stay the circles of the charge.

The same empty-map arithmetic is in `results/checks.txt` for levels 63 (runaway electrons: colour `30000 T`, `B_z/b_ref = 35`), 66 (magnetosphere: same `10000 T` and `117` as level 57), and 82, 86, 87 (crossed fields of the rings: colour about `1.3 T`, energy boundary about `6.2` cells, `B_z/b_ref ≈ 2.9`).

### Other shipped disturbances

| level | what the selected flight has | what the map shows |
|---|---|---|
| 36 Stray field, 45 RF beam line | uniform `E` on realization 1 (`\|Δφ\|` of `−E·x` about `1.3×10^5` and `2.6×10^5`) | realization 0 |
| 40 Earth’s field | `B_z = ±7700` on “facing north / south” | realization 0 |
| 41 CRT, Earth’s field | `B_z = 7700` on the only realization | 0; normalized field is only `0.04`, so the missing colour is faint, and the sign of `B_z` is still absent |
| 62 Gradient drift | a row of magnets, plus `B_z = 3.3×10^5` (`B_z/b_ref = 1.65`) | the magnets only; the equal-`\|B\|` contours miss that additive field |

Oscillating waves (mains hum, the light wave of level 84) are a different view. They are drawn by the waves and total views, which keep `External::Wave` when `ω ≠ 0`.

## 2. Inside a charge cloud the potential map is a point charge

**Shipped, levels 76–79. The trajectory uses the uniform sphere.**

`Coulomb::sample` inside radius `R` is

```text
φ = Q (3 R² − d²) / (2 R³),    E = Q d / R³
```

which is what `PHYSICS.md` §2.1 states, including the sentence that the potential map needs nothing new. The map does not use this. `potential::params` uploads every entry of `Coulomb::charges()`, and that iterator is documented as the point-charge reduction. The shader then evaluates `Q / max(r, 1e-4)` and paints a disk of `physics.charge_radius` (0.3) at the centre. A cloud is penetrable. The disk is not a body.

The turning contour compares that shader potential with a launch energy computed from the CPU potential. For a positive cloud and a negative ion the map’s turning radius is larger than the orbit’s.

| level | launch distance | turning radius, trajectory | turning radius, map | colour error at the launch |
|---|---|---|---|---|
| 76 bound charge | 2 (mid-radius) | 2.828 | 3.200 | `−5 T` |
| 77 resonance | 0 (the centre) | 1.000 | 2.723 | launch sits under the false disk; at the disk edge, `−379 T` |
| 78 bound knock | 2 | 2.828 | 3.200 | `−5 T` |
| 79 spectroscopy, `R = 4` | 0 | 1.000 | 2.723 | disk edge `−379 T` |
| 79 spectroscopy, `R = 3` | 0 | 1.000 | 2.077 | disk edge `−153 T` |

Outside the sphere the two potentials agree (ratio 1 on the surface in `checks.txt`). At `d = R/2` the map is high by a factor `1.45`; at `d = R/4`, by `2.72`. The centre value in the shader is the clamp `Q / 10^{-4}`, which the disk covers.

Level 77’s ion is launched at the cloud centre and oscillates out to 1 cell. The map draws the turning line at 2.72 cells and a solid dot of radius 0.3 on the launch point.

## 3. A static antenna stores the dipole field and a zero potential

**Sandbox. No shipped antenna has `ω = 0` (the RF levels use 0.6, 0.4, 0.3, 0.125, 0.1925, or a wave at `ω = 1`). The blank sandbox starts `rf_omega` at 0, and an antenna with no frequency of its own uses that value.**

`OscillatingDipole::fields` returns Jackson’s `E` and `B` and sets `phi: 0.0` on both the `c = ∞` branch and the retarded branch (`antenna.rs`). At `ω = 0`, `ṗ = p̈ = 0`, so

```text
E = [3 n (n·p) − p] / r³,    B = 0
```

which is `−∇φ` for `φ = n·p / r²`. Wolfram reduces `−∇φ − E` to `{0, 0}` (`audit/scripts/wolfram/checks.wls`). The Python central difference at `(1.7, 0.4)` with `p = (0.8, −0.3)` gives `|E − (−∇φ)| = 1.9×10^{-11}`, while the code’s `φ` is 0. `PHYSICS.md` §2.4 calls the `c = ∞` case the electrostatic dipole field.

`LevelField::is_static` is true when every antenna has `ω = 0` and nothing else depends on time. The energy diagnostic then integrates `W = T + q φ − m B_z` (`trajectory.rs`) and records `|ΔW|`. With `φ` missing, `W` changes by the work `q v·E` and the diagnostic reports that a static field does not conserve energy. The force itself uses `E`, so the trajectory is the dipole orbit. Oscillating antennas are unaffected: `is_static` is false and the diagnostic is `NaN`, which is the right treatment of a non-conservative field. Test A2 compares `E` with a charge pair and does not fly a static antenna.

The potential map has no antenna term either, so the uphill colour and the dark region omit `φ = n·p / r²`. The total view samples `E` and shows the dipole.

## 4. Magnetic-moment bans look only at shots

**Sandbox. No shipped level has `free_particles`, and the only moments in `levels/` are the shots of 29 and 47.**

`PHYSICS.md` §3.2 says `Level::model_issues` rejects a magnetic moment together with electric sources at finite `c`, with radiation reaction, and with ramped coils. The force that is integrated is `q(E + v×B) + m ∇B_z` (`dynamics.rs::force`, `beam.rs::force_in`). The Aharonov–Casher term is absent, which is why the rejection exists.

The rejection runs only when some shot has `particle.moment ≠ 0` (`lib.rs` around the `model_issues` moment block, and again for ramped coils). A moment on `free_particles[].particle.moment` is edited in the sandbox (`level_editor.rs`, “moment m_z”), copied onto the beam particle (`beam.rs::extra_particles`), and given `m ∇B_z`. It does not open the rejection. Player-placed free charges are hard-wired to moment 0, so they are not this hole.

The electric-source test, when it does run, also skips charge clouds. It looks at elements that are not magnets, at `limits.max_charges`, antennas, plates, conductors, electrodes, and disturbances. Clouds live in `Level::clouds` and contribute `E` through `Coulomb::sample`. A shot that carries a moment, at finite `c`, in a level whose only electric source is a cloud, passes `model_issues` and is integrated without the hidden-momentum term.

Level 29 is finite `c` and magnets only, which is the case §3.2 calls exact. Level 47 is Newtonian. Neither file is this hole.

## 5. Electrode pictures are point charges, and the legend says the map is exact

**Shipped on the electrode levels (19, 22–25, and the other files that contain `"electrodes"`). The integrated field is the panel integral.**

`Electrodes::sample` uses the triangle integrals except at display resolution, where `point_evaluation` replaces each panel by `σ · area` at the centroid and at its mirror image (`bem.rs`). `potential::params` does the same upload itself, so the potential map never evaluates a panel. Field lines are built from `display_scenario`, which is display resolution, so they follow the point charges too. Preview and verification resolutions keep the integrals. The code comment marks the point charges as pictures.

Far from a panel the monopole is the right leading term. On the panel the true potential is finite and the normal field jumps; a centroid in the plane is a `1/r` spike clamped at `10^{-4}`. This audit did not measure that difference on a shipped electrode. The legend’s “exact” covers this map as well as findings 1 and 2.

## 6. Sentences that the implementation has left behind

These do not change a flight. They will mislead the next edit.

- `PHYSICS.md` §2.1: the cloud “potential map need nothing new.” Finding 2.
- `PHYSICS.md` §2.4 calls the infinite-`c` antenna the electrostatic dipole and does not say that `φ` is stored as 0. Finding 3.
- `PHYSICS.md` §3.2, the map formula and the claim that `model_issues` rejects every moment-plus-electric combination. Findings 1 and 4.
- `PHYSICS.md` §3.3, the paragraph that ends the beam section: `model_issues` “rejects … radiation reaction,” and “particle–particle contact is not an event.” `beam_issues` does not mention radiation reaction. Radiation reaction on a beam is implemented and described earlier in the same section. Opposite charges are rejected only when `radius == 0`; spheres with a radius take `elastic_impulse` and the step ends at contact. The collisions paragraph a few lines later already says so. The same section says the energy books “balance by construction” and that the tapered past “changes by at most 0.1 c” in order to stay below `c`. Findings 12 and 11.
- `PHYSICS.md` §4 defines `energy_max_abs_error` as the largest `|W(t) − W(0)|` with `W = (γ−1)mc² + qφ`. The runner stores `|W − W(0) − ∫ v·F_RR|`, and `W` also contains `½ q φ_self − m B_z` (`trajectory.rs`).
- `PHYSICS.md` §6 names `Trajectory::closest_sampled`. The trajectory struct has no such field. Closest approach is the margin from `minimize_on`.
- `beam.rs` crate header: the default finite-`c` interaction is “constant acceleration (`accelerated_fields`)”. `continuations` builds `FieldMotion` for a charged particle with no moment, and mixes in `accelerated_fields` only as that motion’s weight falls. The comment above the step cap still says every retarded step is “shorter than a third of the light travel time”; the formula two lines later is `h ≤ (E L + r_min/c) / (1 + w/c)`, which is what runs.
- `External::sample`, on an `ω = 0` wave: “with `B = 0`, since `k·x/c` enters only multiplied by `ω`.” `PlaneWave::fields` still sets `B = k̂ × E / c` whenever `c` is finite. The phase is independent of position at `ω = 0`; the ratio `B/E` is not. `PHYSICS.md` §2.3 writes `B = k̂ × E / c` for the wave and, separately, that `ω = 0` is a static uniform field with its potential. The potential is implemented (`φ = −E·x`). The magnetic piece matches the wave formula and not the comment. At `c = ∞`, `B` is 0 either way. The script’s transcription gives `|B| = |E|/c` at `ω = 0` and at `ω = 1.3` (`c = 4`). Test W4 checks the energy, and a uniform `B` does no work. No shipped wave has `ω = 0`. The UI labels that term “static field.”
- `panel.rs` module note: the normal field “has the jump `±2πσ`.” The test `jump_across_the_surface` requires a jump of `4π` in `∫∇(1/R)`, and `E = −σ ∫∇(1/R)`. With `k = 1` the discontinuity of `E` is `4πσ`. `2πσ` is the field on one side of an infinite sheet. The test and the integral agree with the unit system. The note mixes the one-sided field with the jump.
- `Level::verify_flights` and the worker header: preview and verification “differ only with metal spheres.” `has_metal` is spheres, electrodes, or player plates. For an interacting beam at finite `c` the verification flight is also the retarded interaction.
- `level/src/beam.rs`: mixing a beam shot with a single shot “is not supported, `model_issues`.” `model_issues` has no such check. A non-beam shot in a beam level flies as one particle, which is what §3.3 describes.

## 7. Latent

Nothing here is reachable from a shipped file. Each one is a real disagreement in the source.

**`FieldMotion::Retarded.t` is the emission time minus the observer time.** `retarded()` stores `t − t0 − dt` (`field_motion.rs`). `fields()` passes that value into `lienard::fields_from` as the retarded time. The Liénard–Wiechert expression uses position, velocity, and acceleration; it does not use the time value to form `E` or `B`. A caller that reads `Retarded.t` as an absolute emission time gets the difference instead. `LwField.retarded_time` keeps the same number.

**Landau–Lifshitz does not see the image field.** `ParticleOde::force` adds `self_field`. `radiation_reaction_force` differentiates `field.sample` only, and the beam copy does the same with the other particles’ retarded fields added (`beam.rs`). The image field is smooth at the particle and is part of the field in the Landau–Lifshitz formula. No shipped radiation-reaction level has conductors or electrodes (inventory in `checks.txt`). `model_issues` does not reject the combination.

**A NaN event value is dropped in release.** `events::search` returns `None` when an endpoint is NaN, after a `debug_assert`. The comment records a zero-length wire that once produced one. In release that interval is treated as clear.

**The particle-field view rebuilds acceleration from `PathPoint.force`.** That force is `(E + v×B) q + m ∇B_z` (`worker.rs::build_path`). The integrator’s force also has the image field and, when the level asks, the radiation-reaction force. Positions on the track come from the integrator, so the drawn path is the integrated one. The `β̇` used for the Liénard picture is not. Shipped radiation-reaction levels are required to keep `|F_RR|/|F_L| < 0.05`, which bounds the missing reaction piece. The image piece is present on metal levels.

**Antenna and wave phase is reduced in f64, then `ω r / c` is subtracted in f32** (`radiation.rs` upload, `radiation.wgsl` `antennas_at` / `waves_at`). A long light-travel phase loses low bits in f32. Shipped arenas are tens of cells.

**A repeated polygon vertex makes a ramped coil’s vector potential NaN.** “+ vertex” copies the last vertex (`level_editor.rs`). `segment_potential` divides by `|b − a|`. A steady current never asks for that potential (`induced_e` returns 0 when `rate == 0`), and the segment field off the point is 0 because `r_a × r_b = 0`. A nonzero ramp calls `unit_vector_potential` at every sample, so one repeated vertex makes `E_induced` NaN everywhere. Shipped coils have distinct vertices.

**Body-overlap checks add a hardcoded 0.3.** Electrode and sphere containment in `model_issues` uses `CONTACT_DISTANCE + 0.3` for elements and `CONTACT_DISTANCE` for launches. `0.3` is `default_radius()`. The editor can change `charge_radius`, and a particle can have its own radius. Shipped files use 0.3, so the shipped levels match the constant. A sandbox body with another radius can overlap metal while `model_issues` is empty. Contact during a flight uses the particle radius (`trajectory.rs`).

## 8. Already stated in PHYSICS.md

These are limits of the model, measured there, and not re-opened here.

- Event speed bound `v_max` is `1.25` times the fastest of five samples on the step, capped at `c` (`trajectory.rs`). A faster excursion inside the step can pass a boundary. `PHYSICS.md` §6 describes this residual.
- The quasi-static beam interaction, its indicator, the Faraday-cup fade as a model of the pipe, and the neglected electrode image force are approximations with the bounds in §2.7 and §3.3. The energy panel’s total is a separate defect (finding 12).
- Metal together with antennas, waves, or ramped coils is rejected. Time-dependent `B` is left out of `grad_bz` for the same reason.
- `c = ∞` is Newtonian: `inv_mc_sq = 0`, wave `B = 0`, no radiation reaction.

## 9. Small things that are easy to misread

`potential::params` silently keeps the first `MAX_CHARGES` (1024), `MAX_MAGNETS` (64), `MAX_LOOPS` (16), and `MAX_SEGMENTS` (64). A sandbox past those caps draws a prefix of the sources. Shipped levels are far under the caps.

`PlaneWave::vector_potential` divides by `ω`. It is used from the wave tests at `ω ≠ 0`. An `ω = 0` call diverges. The flight does not call it; the static term uses the scalar potential.

The static scan (`results/static_scan.txt`) found no `unsafe`, no `mul_add`, no `todo!`, and no `unimplemented!` under the four `src` trees. `f32` appears in the game crate (the GPU upload) and not in `physics`, `level`, or `generator`. The `unwrap` and `expect` hits in `physics` are lock poisoning, non-empty geometry, test asserts, and “bounds/detector present,” not a `Result` ignored on the force evaluation.

## 10. A launch already inside the detector is an arrival

**Sandbox. No shipped launch node has signed distance ≤ 0 to its detector, and none lies inside a gate. The closest shipped node is 2 cells outside (`58_jackson_van_allen.json`).**

`run_cancellable` builds the events in order: obstacles, then the arena edge, then the detector. If any event function is ≤ 0 at `x0`, it sets `outcome` from that event and returns. The acceptance block sits after the integration loop. It rewrites `Arrived` to `SkippedGate` when a gate is still missing, and to `Rejected` when the direction, energy, or radiation window has a negative margin. The early return never reaches it.

The detector event is the box’s signed distance, negative inside. A 2D detector is a slab `z ∈ [−1, 1]` (`Level::box_region`), so a launch node in the rectangle is inside the slab. `model_issues` rejects a launch inside a gate, inside a metal sphere, and inside an electrode. It does not reject a launch inside the detector. The editor edits the detector box and the launch node.

The beam runner does the same test and marks the particle `Done` with `Outcome::Arrived` before the loop. Gate and acceptance checks run only when an event fires during a step. A cup entered that way sets `t_off = 0` and does not install a `Fade`, so the charge is absent from the force at once.

A flight that enters the detector during the integration does go through those checks. The hole is the launch state.

## 11. The tapered past can go faster than light

**The bound in the comment is false. Every shipped beam launch is still under it: the fastest upper edge is `β = 0.805` on level 48, and `0.805c + 0.1c = 0.905c`. Seven single-particle shots launch faster than `0.9c`; their integrated state is a momentum, which stays below `c`.**

`tapered` continues a world line from `(r, v, a)` by `τ`. With `T = 0.1 c / |a|` and `u = τ/T`, the acceleration is constant for `|u| ≤ 0.7` and `a sech²` after that. The velocity change along `â` is at most `0.1c`. Far out,

```text
v_∞ = v + sign(u) · â · 0.1 c
```

A particle that is braking has `â` opposite `v`, so the past (`u < 0`) is faster: `|v_∞| = |v| + 0.1c`. Inside the constant piece the change is at most `0.07c`, so a brake from `|v| > 0.93c` is already superluminal before the `sech²` starts. Speeding up does the same thing to the future. The comment on `accelerated_fields` says the velocity “changes by at most 0.1c, so the past never becomes superluminal.” `PHYSICS.md` §3.3 says the same thing as the reason the old `v + aτ` was replaced. Wolfram integrates `a sech²` and gets this velocity and the code’s position with difference 0 (`checks.wls`). The past limit is `v₀ − a T`, and `T = 0.1 c/|a|`, so a pure brake reaches `c` at `β = 9/10` and, inside the constant piece (`0.07c`), at `β = 93/100`. A Python RK4 of the same curve agrees to `1.4×10^{-14}` in position over `τ = −0.5`.

The script’s transcription, `c = 1`, `|v| = 0.95c`, `a = −c`:

| piece | speed |
|---|---|
| `|u| = 0.7`, braking into the past | `1.02 c` |
| `|u| = 20`, braking into the past | `1.05 c` |
| same acceleration, into the future | `0.85 c` |
| speeding up, into the future | `1.05 c` |

`accelerated_fields` then solves `g(τ) = c(dt − τ) − |x − r(τ)|` with 50 Newton steps and calls `fields_from` on whatever `τ` it has. There is no test that `g` is zero or that `|v| < c`. On a braking charge at `0.95c`, an observer at `(0, 0.4)` stops at `τ = −0.373`, `g = −0.183`, `|v| = 1.05c`. An observer ahead of it stops at `g = −0.211` after the iteration has already visited `1.05c`. An observer behind a charge that is speeding up, whose past is slower, stops on a root (`g = 0`). `lienard::retarded_time` expands its bracket while `g ≤ 0` and has no superluminal guard either; once `|v| > c`, `g` need not become positive.

Callers: the exact beam interaction’s pre-launch world line (`View::before_launch`); the quasi-static interaction, both as `Continuation::Accelerated` and as the blend inside `Continuation::Fields` when the uniform-field motion’s weight drops; the field view (`radiation.rs`, `with_past` and the neglected-field comparison). `FieldMotion` is a different continuation. This audit did not re-derive it.

Shipped launch speeds (`β = √(ε(ε+2)) / (1+ε)`, `ε = T/(mc²)`; a Gaussian beam’s upper edge is `T(1 + 3 σ)`):

| file | `β` | what it is |
|---|---|---|
| 80, 81, 83, 84, 85, 88 | 0.9428 (`γ = 3`, `c = 2`, `T = 8`) | one particle |
| 30 Beta-ray spectrometer, shot 0 | 0.9102 | one particle; the description’s `0.91c` |
| 48 Relativistic beam | 0.800, edge 0.805 | the fastest shipped beam |

A transverse `0.1c` added to `β = 0.943` stays below `c`. The anti-parallel case crosses `c`. Launch accelerations of those seven shots were not evaluated. Their integrated trajectories stay below `c`, because the state is a momentum.

## 12. The cup fade does work the energy total never receives

**Shipped on the Newtonian interacting beams: 46 Space charge, 49 Collimated beam, 50 Velocity selector, 51 Beam preparation, 56 Isotope separator. The particles still fly under the fading force. The panel’s total does not stay constant, and for `c = ∞` it says the total is conserved.**

The levels’ detector fate is `Fate::Cup` unless `instant_drain` is set. None of these files set it. On entry the runner adds to `absorbed` the drop in kinetic plus potential plus pairwise Coulomb energy, then removes the particle from the charges present. A `Fade` keeps pushing the others: `force += d · (q q_f factor(t) / r³)` with `factor = exp(−k v_n (t − t_off))` (`beam.rs`). That force is absent from `budget`, which sums only the charges still present and the charges stopped on bodies. Nothing afterwards adds the fade’s work to `absorbed`.

After the booking, the displayed total

```text
kinetic + (potential − potential₀) + (interaction − interaction₀) + absorbed − kinetic₀
```

(`ui.rs`) therefore changes at the rate `v · F_fade`. For `c = ∞` the note under the bar says “The total is conserved.”

The diagnostic knows. While any fade has `factor > 10^{-15}` it resets the baseline, and the comment says the energy is going into the cup. `PHYSICS.md` §3.3 says the same, and also that the books balance by construction. The construction balances a drain. `Fates::default` uses `Fate::Drain`, which is the fate in the description of test B16. A cup is the fate the levels use.

The script integrates one survivor, charge 1, mass 1, starting at rest at the origin, while the absorbed charge starts at `(2, 0)`, flies at speed 1, and fades at rate `π/2` (mouth width 2). Mutual energy at entry is `0.5`. At `t = 8` the survivor’s kinetic energy is `5.204356351439677×10^{-3}` in Wolfram’s `NDSolve`, which matches the Python integrator and is the whole drift of the total. An instant drain leaves it at 0. This geometry is not a shipped flight; it is the force and the booking.

The same match arm adds the neglected impulse `|q| / (r c)` for both `Drain` and `Cup`, and only in the quasi-static interaction. The comment says the field drops at once and lingers for a light time. A cup does not drop the field at once. The term still inflates `neglected_retardation`. At finite `c` the panel does not claim conservation: the note already leaves magnetic and radiation field energy out of the mutual column, and `energy_max_rel_error` is `NaN` whenever the particles interact. The booking of the fade is the same.
