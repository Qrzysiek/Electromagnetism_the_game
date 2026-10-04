# Checks that matched

These were compared with the source or recomputed by `audit/scripts/audit_checks.py`. They are not a substitute for `cargo test`. The later pass, with Wolfram residuals for the CPU formulas and the shaders, is [FORMULAS.md](../FORMULAS.md).

## Antennas, except the static potential

`OscillatingDipole::fields` matches Jackson §9.2 as printed in `PHYSICS.md` §2.4: the near, induction, and radiation terms, `B = (ṗ × n / r² + p̈ × n / (c r)) / c²`, and `B = 0` at `c = ∞`. The Larmor power is `p₀² ω⁴ / (3 c³)`. At `ω = 0` that `E` is `−∇(n·p / r²)`; the stored potential is 0. See FINDINGS §3.

## Kinematics and the force

`Kinematics::kinetic_energy` is `p² / (m (γ + 1))`, finite at `c = ∞`. `Kinematics::new` stores `inv_mc_sq = 0` when `c` is infinite, so `γ = 1` and `v = p/m`.

`ParticleOde::force` is `q (E + E_image + v×B) + m ∇B_z`. `grad_bz` sums magnets and coils, scaling a ramped coil by `κ(t)`, and omits antennas and waves. The omission is the one `model_issues` is meant to enforce for shots.

`landau_lifshitz` matches the expression printed in `dynamics.rs`: the coefficients are `2 q³ / (3 m c³)` and `2 q⁴ / (3 m² c⁴)`, the first term is `γ` times the convective derivative of the fields, and the last term subtracts `v` times `γ²` times the square of the Lorentz force with the ` (E·v)²/c² ` piece removed. The derivative is a central difference over `10^{-5}` cells. This audit did not re-run tests R1–R6.

## Clouds, outside the map

`Coulomb::sample` branches on `d² < R²` and uses `φ = Q (3 R² − d²) / (2 R³)`, `E = Q d / R³`. Wolfram’s gradient of that potential is `Q r / R³` inside and `Q r / r³` outside, and the two potentials agree on `r = R` (`audit/results/wolfram.txt`). Metal is given the same clouds as point charges (`field_at`), and `model_issues` keeps each cloud clear of spheres and electrodes. That split is what §2.1 describes. The map is the exception, FINDINGS §2.

## Coil bracket

`loop_bracket` drops the `m⁰` and `m¹` terms. The script evaluates the series and the direct combination `(1 − m/2) E − (1 − m) K` from the same AGM.

| `m` | `(series − direct) / direct` | `series / (π m² / 32)` |
|---|---|---|
| `10^{-8}` | the direct value is rounding (`~10^{-16}`); `series / lead = 1` | 1.000000 |
| `10^{-4}` | `−4.5×10^{-8}` | 1.000025 |
| `0.02` | `−9.9×10^{-13}` | 1.005047 |
| `0.049999` | `−1.8×10^{-13}` | 1.012802 |

The series is the stable evaluation. The direct form only has digits once `m` is large enough that the cancellation is mild, and there the two agree. Wolfram’s series of `g(m) = (1 − m/2) E(m) − (1 − m) K(m)` starts at `(3π/32) m²`; the coefficients of `m⁰` and `m¹` are absent. The vector-potential bracket `h(m) = (1 − m/2) K(m) − E(m)` starts at `(π/32) m²`, which is the leading term in the comment on `potential_bracket`.

## Panel jump

`jump_across_the_surface` expects `∫∇(1/R)` to change by `4π` across the triangle. With `E = −σ` times that integral and `k = 1`, the jump in the normal field is `4πσ`. That is the Gaussian discontinuity for this unit system (`σ/ε₀ = 4π k σ`). The module comment’s `±2πσ` is the one-sided infinite-sheet field, not this jump.

## Events

`first_crossing` certifies an interval when `g(a) + g(b) − v_max (b − a) > 0`, and otherwise bisects, left half first. A NaN endpoint returns no crossing (release) or asserts (debug). The speed bound used by the trajectory is the documented `1.25` times five samples, capped at `c`.

## Arrival during a flight

Once the integrator is running, a detector entry is rewritten from `Arrived` when a gate is still missing or the acceptance margin is negative (`trajectory.rs` after the loop; the beam copy at the event). The launch state returns before that rewrite. See FINDINGS §10.

## The taper’s integral

`tapered`’s closed form matches an RK4 of `dv/dτ = a sech²((|u|−0.7)/0.3)` to `1.4×10^{-14}` in position over `τ = −0.5` (`results/checks.txt`). Wolfram’s integral of the same acceleration differs from that closed form by 0, and the past limit is `v₀ − a T` (`results/wolfram.txt`). The speed bound claimed for that curve is FINDINGS §11. A finite momentum still has `|v| < c`: `Kinematics::velocity` is `p / (γ m)`.

## Beams

`elastic_impulse` at `c = ∞` returns `2 μ (v₁ − v₂)·n` with `μ = m₁ m₂ / (m₁ + m₂)`, and only when the particles approach. The relativistic root is a bisection of `T(p₁ − J n) + T(p₂ + J n) = T_before`.

The cup’s charge factor `e^{−k v_n (t−t_off)}` with `k = π/w` is the pipe mode `PHYSICS.md` §3.3 writes down, and the force loop applies it. The energy total’s columns leave that work out. See FINDINGS §12.

Opposite point charges are rejected using shots, free particles, and the player’s free-charge list (`beam_issues`). That gate is wider than the moment gate.

Rayon’s `into_par_iter().map().collect()` on a slice keeps index order. The beam right-hand side and the spectrum’s parallel directions rely on that. The conductor `HashMap` in `merge` is a first-seen index; the summed vector stays in input order.

## Editor

`edit_level` and `check_editable` destructure `Level` by field name. A new level field fails to compile until it has a control or is explicitly set aside. The `..` tokens in that file are numeric ranges.

## Static scan

Under `crates/physics/src`, `crates/level/src`, and `crates/generator/src`: no `unsafe`, no `mul_add`, no `f32`, no `todo!`. The game crate’s `f32` count is the GPU path, which the project rules keep off the authoritative trajectory.
