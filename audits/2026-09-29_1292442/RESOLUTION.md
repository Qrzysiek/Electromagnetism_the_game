# Resolution of the 2026-09-29 audit

- **Audited code:** commit `25f821a` (the working tree had no other changes; `git status` in
  REMAINDER.md).
- **Fixed in:** commit `1292442` (2026-10-05).
- **Next audit:** review `git diff 1292442..HEAD` (plus anything this table marks as open).

Every remark was checked against the code before it was fixed. Verdicts: **valid** (the
remark holds and was fixed), **valid, documented** (it holds, and the fix is a stated limit
rather than a code change), **no action** (the audit itself classed it as already stated).
Paths are relative to the repository root.

## FINDINGS.md

| # | Remark | Verdict | What was done |
|---|---|---|---|
| 1 | Potential and magnetic maps omit the uniform stray fields, stay on disturbance 0; field lines likewise, none for a pure uniform E | valid | The maps draw the static part of the active flight's field, including the disturbance's uniform `E` (potential `−E·(x − x_launch)` in the shader) and `B_z`; `Level::display_scenario(shot, disturbance, …)`; the map and field-line keys include the disturbance; field lines draw a uniform E and static antennas (`crates/game/src/potential.rs`, `potential.wgsl`, `main.rs`, `visuals.rs`). Found while fixing: a particle launched at rest gave `b_ref = 1e-300`, so `B_z/b_ref` was infinite in f32 and the map black (latent with magnets before); `b_ref` now uses the momentum of the energy unit. The dark region was also drawn for oscillating antennas and waves (non-conservative): now cleared for any time-dependent field. Screenshots of levels 57 and 36 confirm. |
| 2 | Inside a charge cloud the map is a point charge with a solid disk | valid | Clouds are their own kind of source in the shader, with the uniform sphere's potential inside and no disk (`Coulomb::sources`). Level 76/77 screenshots: turning lines at the flight's radii. |
| 3 | A static antenna stores φ = 0 | valid | `OscillatingDipole::fields` returns the Lorenz-gauge `φ = n·[p]/r² + n·[ṗ]/(cr)`, `vector_potential` the `A`; `LevelField` adds `φ` for ω = 0 antennas; the map draws them. Tests: unit `potentials_generate_the_fields`, A2 (potential), A4 (a static antenna's flight conserves `T + qφ`; it reported half the launch energy as error before). |
| 4 | Moment bans look only at shots; clouds not counted as electric | valid | `model_issues` counts free particles' moments and clouds (and free charges) (`crates/level/src/lib.rs`); test `moments_with_electric_sources_are_rejected`. |
| 5 | Electrode pictures are centroid point charges; legend says exact | valid | Measured 4.7e-2 of the plate potential next to a plate (level 19). Pictures now integrate panels within two sizes exactly and the rest by a three-point rule, on the CPU (`Electrodes::picture`, also the total view and field lines, which the audit did not list) and per pixel in the shader. Test E6: 4.3e-6 of the plate potential, 3.5e-5 of the field. The display mesh's own error (4.2e-3 against the verification mesh) is stated in the legend. |
| 6 | Sentences left behind | valid (each) | PHYSICS.md §2.1, §2.3, §2.4, §2.7, §3.1–§3.3, §4, §6, §10 rewritten; `beam.rs` header and step-cap comment; `External::sample` (code changed: an ω = 0 wave term is now the uniform electric field the format and UI describe, `B = 0`); `panel.rs` jump `4πσ`; `verify_flights` and the worker header; `level::beam` mixing comment. |
| 7a | `FieldMotion::Retarded.t` is relative | valid | Now the absolute emission time (`field_motion.rs`). |
| 7b | Landau–Lifshitz does not see the image field | valid | The reaction's fields include `self_field` (single flights and beams). Radiation reaction next to metal is now rejected by `model_issues` (the electrostatic metal reflects no radiation, an O(1) change of the reaction near it); test `radiation_reaction_near_metal_is_rejected`. |
| 7c | A NaN event value is clear in release | valid | It fails the flight (`Outcome::Failed(NonFiniteEvent)`), both runners; unit test `nan_is_reported_not_cleared`. |
| 7d | The particle-field view's acceleration omits image force and radiation reaction | valid | `build_path` takes the flight's whole force (`ParticleOde::force` + reaction); the energy panel's potential gains `½ q φ_self` (found while fixing). |
| 7e | f32 phase of antennas and waves | valid, documented | Below 1e-5 rad for `ω r/c ≤ 100` (arenas are tens of cells); PHYSICS.md §10. |
| 7f | A repeated polygon vertex makes a ramped coil's induced field NaN | valid | `segment_potential` skips a zero-length side; unit test `repeated_vertex_adds_nothing`. |
| 7g | Containment checks add a hardcoded 0.3 | valid | Each element kind's radius, and the particle's for launch points (`Level::bodies`); test `containment_uses_the_bodies_radii`. |
| 8 | Already stated in PHYSICS.md | no action | |
| 9 | Silent caps of the map's arrays; `PlaneWave::vector_potential` at ω = 0 | valid | The map's sources are a storage buffer (no caps); `vector_potential` returns `−E t` at ω = 0 (unit test extended). |
| 10 | A launch inside the detector is an arrival without gates or acceptance | valid | Judged like any entry in both runners (a cup also fades the charge from launch); `model_issues` lists it. Tests D4, `containment_uses_the_bodies_radii`. |
| 11 | The tapered past can exceed c | valid | `taper_reach`: `Δv = min(0.1c, 0.9(c − |v|))` (unchanged for every shipped beam, all below 0.89c); the retarded-time solve on the taper is bracketed; `lienard::retarded_time` gives up after 200 doublings; the WGSL copies follow. Test B24. |
| 12 | The cup fade does work the energy total never receives | valid | The work of absorbed charges still acting (cup fade; retarded drains) is integrated per step and taken out of `absorbed`; the neglected impulse `|q|/(Rc)` is added for drains only. Test B16 extended: 8.0e-11 of T through a cup fade whose work is 1.1e-2 of T. |

## REMAINDER.md

| Remark | Verdict | What was done |
|---|---|---|
| Hyperbolic `FieldMotion::retarded` returns a "root" from overflowing terms | valid | `None` where the light-cone function's rounding exceeds 1e-10 of the observer's distance (the beam path gave that point weight 0 already). |
| Debug builds of B6 and B11 panic at `debug_assert!(ga > 0)` | valid | Two particles entering their detectors together: the second one's event is now at the step's start. This was also the Windows CI failure. |
| Generator exits 101 on level 59 | stale binary | `target/release/generator.exe` predated the free-charge element; the rows 57–88 were measured with a fresh build and added to `docs/LEVELS.md`. |
| Difficulty table ends at level 56 | valid | Rows 57–88 measured (1000 samples, 16 × 400) and added. |
| Four display names omit the year | valid | `docs/LEVELS.md` matches the JSON names (0 mismatches). |
| Electrode image-force denominator and range in PHYSICS.md | valid | §2.7: `max(|qE|, F₀)`, range 1.8e-13 to 1.8e-11. |
| The ghost continuation's 10 000-step cap is undocumented | valid | §6 states it (conservative). |
| Landau–Lifshitz across an instant drain's light cone | valid | `model_issues` rejects radiation reaction with the instant drain; §3.1. |
| Everything marked "match" | no action | |

## FORMULAS.md

| Remark | Verdict | What was done |
|---|---|---|
| Poynting prefactor vs Jackson's Gaussian form | valid | PHYSICS.md §10 states the conversion (`B_G = cB`). |
| Legends: potential "Exact", particle field with moments, total-field grid, waves in f32 | valid | All four legends state their limits (`crates/game/src/ui.rs`). |

## Found while fixing, beyond the audit

- **CI never passed on Linux**: the golden hashes differed (level 8 first) because `sin`,
  `exp`, `pow`, `atan2`, … come from the platform's C library. The authoritative path now
  uses the `libm` crate; `clippy.toml` disallows the standard library's versions there; the
  golden hashes were regenerated once (all 88 changed in the last bits).
- The black map for particles launched at rest with a magnetic field (finding 1's row).
- The dark region shown in non-conservative fields (finding 1's row).

## Files not archived

`results/dop853.f` (E. Hairer's Fortran DOP853, used to cross-check the coefficients) has
no licence statement and is not committed; its SHA-256 is recorded in
`crates/physics/src/integrator/dop853_coefficients.rs`. `scripts/__pycache__` is a
build artefact.
