# Formula audit

The CPU force and field formulas that were re-derived here match the independent references inside the limits the model states. The GPU shaders copy those formulas, or their exact reduction to the plane `z = 0`, with floors that sit inside a painted body or below `f32` precision. No new disagreement turned up in a solver. The formulas that disagree with the picture, or with a comment, are the ones already filed in [FINDINGS.md](FINDINGS.md). Three legends omit a limit that `PHYSICS.md` already states, and the Poynting citation is easy to misread across unit systems.

The reference is Wolfram Engine 14.2, plus the closed forms in Jackson and Landau & Lifshitz written for this force law: `F = q(E + v × B)`, `k = 1/(4πε₀) = 1`, so `μ₀/4π = 1/c²`. Reproduce with

```powershell
wolframscript -file audit/scripts/wolfram/formulas.wls
```

The transcript is [results/formulas_wolfram.txt](results/formulas_wolfram.txt). Every line printed a number. Residuals below are copied from that file. The shaders were read against the CPU expressions; they were not executed in this pass. `cargo test` was not run, and the tables in `PHYSICS.md` are the project’s own measurements.

An earlier draft, `formulas_rest.wls`, printed `0.024` for the Landau–Lifshitz comparison because it used one `dB` for both unit systems. With `B_G = c B` and `∂B_G = c ∂B` the residual is 0. That `0.024` is not a physics result.

## Residuals

| Check | Residual |
|---|---|
| CODATA 2018 `k` from `μ₀` against the constant in `units.rs` | `4.36e-12` relative |
| That constant against `10⁻⁷ c²` | `5.49e-10` relative |
| `(γ − 1)mc² = p² / (m(γ + 1))` | 0 |
| Uniform-`E` parallel and perpendicular integrals (T6) | 0 |
| Relativistic Coulomb orbit, differentiated, and the reduced form | 0 |
| Periapsis energy, `m=1`, `c=20`, `v=4`, `b=3`, `κ=2` | `0.e-89` (0 at 40 digits) |
| Deflection against Rutherford at `c = 10⁶` | `1.93e-14` relative |
| Coil `B_ρ` bracket, `A_φ` prefactor, on-axis field | 0 |
| AGM against `EllipticK`, over `m` from `10⁻⁸` to `1 − 10⁻⁸` | `4.44e-16` |
| Series of `g` and `h`, coefficient of `m²` | 0 |
| Segment: the two log forms; `∇ × A − B` at one point | 0 |
| Dipole `∇B_z` against `Grad` of `B_z` | 0 |
| Square-polygon `∇B_z` against a central difference, `h = 10⁻⁶` | `1.27e-9` |
| Loop `∂B_z/∂ρ` from `g(m)` against a central difference, `h = 10⁻⁶` | `2.12e-12` |
| Kelvin potential on the sphere; image energy with the factor `½` | 0 |
| Heaviside field against Liénard–Wiechert at zero acceleration | 0 |
| Landau–Lifshitz, Gaussian rewritten with `B_G = c B` | 0 |
| Uniform-`B` magnetic power against the Liénard power | 0 |
| `(uc)² − \|S\|²` factor; equality when `\|E\| = c\|B\|` and `E ⊥ B` | 0 |
| Plane wave: `∇·E`, `∇×E + ∂B/∂t`, `∇×B − ∂E/∂t / c²` | 0, identically |
| Frozen wave, `ω = 0`: `\|B\| − \|E\|/c` | 0 |
| Dipole `E`, `B` reconstructed from the Lorenz-gauge potentials, step `10⁻⁶` | `1.40e-13` |
| Static dipole: `−∇(n·p/r²) − E` | 0 |
| Field motion: `Λ³ + ω² Λ` and `dU/ds − Λ U` | 0 |
| Same series as the shader, at the shader’s switch `\|ω² s²\| = 1/2`, against `sin` | `8.49e-10` relative |
| Filon recurrence against `Integrate`; series at `D = 1/5` | 0 and `8.06e-17` |
| Sudden-stop integrand is a derivative (Jackson 14.65) | 0 |
| Dipole Larmor average | 0 |
| Quintic blend: value and first two derivatives at the two ends | 0 |
| Shader `ln 2` constant `0.6931472` | `1.94e-8` absolute |
| Elastic impulse: slope at `J = 0` is minus the closing speed | 0 |
| Equal masses, target at rest: the positive root is the incoming momentum | 0 |
| That root against `2μ(v₁ − v₂)` at `c = 10⁶` | 0 at the printed digits |

`N::meprec` fires on the segment curl, the Heaviside comparison, and the field-motion algebra. Those expressions cancel as exact radicals and then hit `$MaxExtraPrecision`. Each of those lines still printed 0.

## Units

`units.rs` stores `K_SI = 8.9875517923e9`, the CODATA 2018 value from `μ₀ = 1.25663706212e-6`. It is the printed-digit rounding of `μ₀ c² / 4π` (`4.36e-12` relative) and sits `5.49e-10` relative from `10⁻⁷ c²`. Inside a level, `k = 1` and `μ₀ ε₀ c² = 1` are imposed exactly, and `c = ∞` is exact Newtonian mechanics (`inv_mc_sq = 0`). The CODATA difference is a conversion-table note. It does not enter a trajectory.

`PHYSICS.md` §1 states the internal choice and names CODATA 2018. That is the right limit to state.

## Motion

The state is `(x, p)`, with `v = p / (γ m)` and `γ = sqrt(1 + |p|²/(m² c²))`. Any finite `p` has `|v| < c`. The kinetic identity that avoids cancellation at low speed is exact (residual 0). `PHYSICS.md` §3 states both, and states that for `γ ≳ 10⁷` the rounded quotient can equal `c`. The integrator’s scaled state and the DOP853 tableau were not re-derived.

**Uniform electric field.** The parallel and perpendicular integrals printed in `PHYSICS.md` for T6 reduce identically to the claimed antiderivatives (residual 0).

**Relativistic Coulomb orbit.** Differentiating the energy and angular-momentum constraint gives the orbit equation in §T3/T5, and that equation’s reduced form is an identity (both residuals 0). At one periapsis the energy residual is 0 at 40 digits. The deflection tends to the Rutherford angle as `c → ∞`, to `1.93e-14` relative at `c = 10⁶`. The text states the condition `Lc > |κ|` and the two regimes (scattering and bound precession).

**Landau–Lifshitz.** The expression in `dynamics.rs` matches Landau’s Gaussian formula once `B_G = c B` and the same factor is applied to the field derivative (residual 0 at one numeric point). In a uniform `B`, with `v` perpendicular to `B`, the two velocity terms together equal the Liénard power, because `1 + γ² β² = γ²` (relative residual 0). `PHYSICS.md` §3.1 states the validity condition (small against the Lorentz force, ratio reported, shipped levels below 0.05) and that the derivative is a central difference over `10⁻⁵` cells. The model note says the same, and badges it approximate.

The derivative is taken of `field.sample` only. The image field is in the Lorentz force and is left out of the reaction. No shipped radiation-reaction level has conductors or electrodes. The model note does not mention that omission. It is already filed with the force assembly.

**Elastic impulse.** `elastic_impulse` returns `2μ (v₁ − v₂)·n` when `c = ∞`, and otherwise the positive root of kinetic-energy balance. The slope of that balance at zero impulse is minus the closing speed (residual 0), so an approaching pair has a positive root. For equal masses and a target at rest the root is the incoming momentum, which is the one-dimensional exchange (residual 0). At `c = 10⁶` that root agrees with the Newtonian value at the printed digits. `PHYSICS.md` §3.3 states the rule and that a relativistic rigid contact is a convention. The bisection loop itself was not executed.

## Static sources

**Charge clouds.** Inside, `E = Q d / R³` and `φ = Q(3R² − d²)/(2R³)`; outside, a point charge. The two potentials meet on the sphere. That identity is in [results/wolfram.txt](results/wolfram.txt) from `checks.wls`. The model note says the cloud is fixed and does not respond, and badges the field exact. That matches the premise in `PHYSICS.md` §8: the sources are prescribed. The potential *map* does not use the interior formula. See below, and finding 2.

**Magnets.** Outside a uniformly magnetized sphere the field is the point dipole. `∇B_z` as written in `magnetic.rs` is `Grad` of that `B_z` (residual 0). In the plane, `m ∥ ẑ` so `m·n = 0` and `B_z = −μ/r³`. The magnetic-map shader evaluates that reduction, with `r` floored at `10⁻³` cells, inside the magnet. The model note badges the exterior dipole exact. `PHYSICS.md` §3.2 states the further limit: with `E = 0` the force `m ∇B_z` is exact relativistically, and with `E ≠ 0` at finite `c` the Aharonov–Casher term is rejected by `model_issues`.

**Circular coils.** The brackets match Smythe’s form (`B_ρ` identity 0, `A_φ` identity 0). On axis, `K = E = π/2` and `B = μ₀ I / (2a)` at `z = 0` (identity 0). The AGM agrees with `EllipticK` to `4.44e-16`. `g` and `h` open at `m²`; the leading coefficients are `3π/32` and `π/32` (both residuals 0), which is what `PHYSICS.md` prints for `h`. The series is used for `m < 0.05`, where the direct difference cancels. `∂B_z/∂ρ` taken from `g(m)` agrees with a central difference of the elliptic `B_z` to `2.12e-12`, the step error at `h = 10⁻⁶`.

The shader evaluates the direct `B_z` form at `z = 0` (there `B_ρ` is identically zero, so the cancelling bracket `g` is not needed), 12 AGM iterations, and `α²` floored at `10⁻⁸`. Quadratic convergence of the AGM is past `f32` epsilon well before 12 iterations, and the floor is inside the wire. The model note badges coils approximate and says the closed-form field of a filament is exact, with thickness only an obstacle. The badge and the sentence together are the right claim. Ramped coils are badged approximate, and the note prints the light-time ratio `(size/c)·|rate|/peak`. That is the validity condition `PHYSICS.md` states.

**Polygon coils.** The closed-polygon gradient, per segment `−F (a − b) × ẑ`, agrees with a central difference of `B_z` to `1.27e-9` at one interior point. That is the finite-difference error. `PHYSICS.md` §3.2 is right that a single segment is not curl-free and the closed sum is. A repeated vertex makes `segment_potential` divide by a zero length. The editor can create one. No shipped coil has one (finding 7).

**Straight segments.** The two logarithms for `A` differ by 0, because both products equal the squared perpendicular distance. `∇ × A − B` is 0 at the test point. The shader’s in-plane `B_z` is the `z` component of the same cross-product formula, with the denominator floored at `10⁻¹²`.

## Conductors and electrodes

A grounded sphere’s Kelvin image puts the potential at 0 on the sphere (residual 0). The field energy of a point charge uses `½ q φ_image`; that factor gives the image force (residual 0 at one exterior point). Without the `½` the force would be twice the image force. `budget` uses `½ q φ_self`. The model note states the measured boundary error, the two resolutions, the neglected electrode image force, and that the metal response is electrostatic. Those are the right limits. The MFS least-squares fit and the BEM matrix were not re-derived.

The triangle integrals are the Wilton 1984 / Graglia 1993 formulas, and the jump test requires a jump of `4π` in `∫ ∇(1/R)` for unit `σ`. In these units the discontinuity of normal `E` is `4πσ` (an infinite sheet, seen from one side, is `2πσ`). The code and the test use `4πσ`. The module comment says the jump is `±2πσ`. That sentence is finding 6. On the triangle the normal component is set to the average of the two sides, which is 0. The model note badges electrodes approximate and names the mesh and the missing image force.

At display resolution each panel is replaced by `σ · area` at its centroid, and the potential map is built from those point charges. The legend still says the map is exact (finding 5).

## Waves, antennas, and the particle’s own field

**Plane wave.** With `E₀ ⊥ k̂`, `|k̂| = 1` and `B = k̂ × E / c`, the three vacuum Maxwell equations hold identically, including at `ω = 0` (symbolic residual 0, and 0 at a numeric point). A frozen wave is a uniform `E` together with a uniform `B` of magnitude `|E|/c`. Both `(E, 0)` and `(E, k̂ × E/c)` are static vacuum solutions. The implemented one is the second. `c = ∞` drops `B`, which is the Newtonian limit the shader also takes (`scales.x == 0`).

`PHYSICS.md` §2.3 prints `B = k̂ × E / c`, and the header formula matches `PlaneWave::fields`. The prose then says `ω = 0` “gives a static uniform field `E₀ ê cos φ`, with its potential,” and does not mention the accompanying `B`. `External::sample`’s comment says `B = 0` for `ω = 0`, which is not what `fields` returns. No shipped wave has `ω = 0`. The in-game model note says disturbances are applied exactly, which matches the Maxwell solution the code evaluates. `vector_potential` divides by `ω` and is test-only. This is finding 6.

**Antennas.** `E` and `B` are Jackson §9.2 in these units, including the static near field `[3n(n·p) − p]/r³` with `B = 0`. Rebuilding `E` and `B` from the Lorenz-gauge potentials by a central difference of step `10⁻⁶` differs by `1.40e-13`, the difference error. The static field is identically `−∇(n·p/r²)` (residual 0). Both CPU branches store `φ = 0`. For `ω = 0` the energy diagnostic then misses `q φ`. For `ω ≠ 0` the potential is not the energy diagnostic’s potential, and the diagnostic is NaN. Shipped antennas have `ω ≠ 0`. The shader copies `E` and `B` (the identity `n × (n × p̈) = (n·p̈)n − p̈` is the one it uses) and has no scalar potential, because the field view does not draw `φ`.

The model note badges antennas exact. That is fair for `E` and `B`. It does not mention `φ = 0`. Finding 3 is the static sandbox case: on `audit/scenarios/static_antenna.json` the kinetic energy swings by `0.48 T₀` while `T + qφ` stays at `10⁻⁵ T₀`, the size of printing positions to four decimals.

The time-averaged power `⟨P⟩ = p₀² ω⁴ / (3 c³)` equals `(2/3)⟨p̈²⟩/c³` with `⟨cos²⟩ = 1/2` (residual 0). `PHYSICS.md` states that the generator supplies it.

**Liénard–Wiechert.** `fields_from` is Jackson §14.1 with `β = v/c` and `β̇ = a/c`:

```
E = q (n − β)(1 − β²) / (κ³ R²) + (q/c) n × ((n − β) × β̇) / (κ³ R),    B = n × E / c
```

At zero acceleration this is the Heaviside present-position field (residual 0 at `v=2`, `c=5`, `X=4`, `Y=3`, `q=3`). The shader’s two-dimensional reduction, `n × (w ẑ) = (w n_y, −w n_x)`, is the cross product for an in-plane `n` and a `z` component. The GPU floors distance at `10⁻⁶` and leaves a hole of radius `0.15` cells on the charge, so the singularity does not fill a pixel. `PHYSICS.md` §2.5 and §10 state that the view continues the world line uniformly outside the sampled interval, and that a real launch or stop would radiate. The single-particle legend says the same for the motion before launch.

**Moving moments in that view.** `PHYSICS.md` §10 says a moment is drawn as `B_z = −m/r³` from the retarded position, and that the moving dipole’s velocity and radiation terms, of order `v/c`, are left out of the picture. The shader does that (`radiation.wgsl`, the `moment` branch). The single-particle legend says “Liénard–Wiechert, exact” and has no moment case. Level 29 has `c = 5` and a moment, so the omitted terms are not negligible on that picture. The dynamics of that level are magnets only; the omission is the view. The legend for a beam does say that the launch acceleration is kept, which matches the dynamics.

**Field motion.** For uniform `E` in the plane and `B` along `z`,

```
U = (γc, γ v),   dU/ds = Λ U,   Λ³ = −ω² Λ,   ω² = (q/m)² (B² − E²/c²)
```

and `U` and `X` are the linear combinations of `U₀`, `ΛU₀` and `Λ²U₀` printed in `field_motion.rs` and in `PHYSICS.md` §3.3. Both matrix identities are 0. The CPU switches to the series at `|ω² s²| < 0.01`. The shader uses the same truncated series and switches at `|z| < 1/2`, where `f32` would lose the closed form’s small differences. At that switch the series differs from `sin` by `8.49e-10` relative, under `f32` epsilon. The hyperbolic branch was read, not given a separate residual. The Newton loop that solves the retarded time was not re-derived.

The quasi-static beam note states the assumption (fields uniform along the path over the light time), that the `1/γ²` magnetic attraction is included, and that the estimated error is shown. That matches the formula’s domain.

**Faraday cup.** `k = π/w` is the decay constant of the lowest Dirichlet mode of a pipe. The model note says the faster modes, the rim charge, and the currents to ground are left out, and badges it approximate. The CPU and the shader use the same exponential, with the positive light delay, and apply it only to emission after entry. The three-dimensional mode was not re-solved here. The fade law is the right textbook term. The energy panel’s total is a separate bug (finding 12): for `c = ∞` the caption says “The total is conserved,” while a fading charge still does work that the total never receives.

**Tapered past.** The quintic blend is `1 − smoothstep`, and its value and first two derivatives vanish at both ends (residual 0). The shader’s `ln cosh` uses `0.6931472` for `ln 2`, absolute error `1.94e-8`. Times `w² = 0.09` and the position scale, that is about `10⁻¹¹ c²/|a|` in the continued position, invisible in `f32`. The comment that “the past never becomes superluminal” is false for a pure brake that is already fast: the speed reaches `c` at `β = 0.9`, inside the constant piece (finding 11). Shipped beam edges stay under that. The integrated trajectory uses momentum and cannot exceed `c`. The guard’s formula is the one `PHYSICS.md` describes; the comment overclaims what it guarantees.

## Radiation integrals and the energy flow

**Spectrum.** `dW/dΩ` and `d²I/dωdΩ` are Jackson 14.38 and 14.65 with `β̇ = a/c`. The Filon moments match `Integrate` for `k = 0…4` at `D = 7/5` (residual 0), and the series at `D = 1/5` to `8e-17`. The sudden-stop identity, that the integrand of 14.65 is the derivative of `n × (n × β)/κ`, is 0. The quartic divided differences were not re-derived. The model notes state the uniform-before-launch and uniform-after-arrival approximation, and that an abrupt stop is exact only for frequencies far below the inverse of a real stopping time. `PHYSICS.md` §3.4 states the arc spacing, the frequency grid, and the in-plane receiver.

**Poynting.** In these units `u = (E² + c² B²)/8π` and `S = (c²/4π) E × B`. Then `(uc)² − |S|²` reduces to a square that vanishes when `|E| = c|B|` and `E ⊥ B` (both residuals 0), which is the statement that the energy velocity is at most `c` and equals `c` in a radiation field. This is the SI expression with `k = 1` and `μ₀ = 4π/c²`. Jackson’s Gaussian formula in §6.7 is `(c/4π) E × B`, smaller by a factor `c`.

`poynting.rs` says “Jackson §6.7, in the game’s units” and then writes the converted factor. `PHYSICS.md` §10 and the on-screen line write `S = (c²/4π) E × B` and cite Jackson §6.7 without the unit sentence. The factor on the screen is the one that belongs to this force law. A reader who copies the Gaussian prefactor out of that section would be off by `c`.

The rest of the on-screen text is the right set of limits: only the flux through a closed surface is unique, the map is a slice of a three-dimensional flow, charges beside magnets circulate field momentum, and near a particle a third of the work at low speed comes from the stored exchange energy. The tracer hover states the speed bound that the identity proves.

## What the legends say

The hover on an exact model note is “Exact within classical electrodynamics, up to the integration accuracy.” The hover on an approximate note is “An approximation or omission; where it can matter, its size is measured.” Those two sentences match how the notes are badged. The notes that state a limit well are listed above: radiation neglected or Landau–Lifshitz, electrodes, metal, ramped coils, the cup, quasi-static beams, radiation-goal end conditions, fixed clouds, and the moment’s neglected radiation.

Four wordings do not carry the limit the formula actually has.

| Where | What it says | The limit that belongs with it |
|---|---|---|
| Potential-map legend | “Exact” | Uniform stray fields are absent (finding 1). Inside a cloud, outside the solid disk of radius `charge_radius`, the map is `Q/r` (finding 2). Display electrodes are point charges at the centroids (finding 5). The shader also floors `r` at `10⁻⁴`. |
| Particle-field legend, one particle | “Liénard–Wiechert, exact” | True for a charge whose world line is the sampled flight, continued uniformly outside it. A moment is `B_z = −m/r³` from the retarded position, without the moving-dipole terms of order `v/c` (`PHYSICS.md` §10). Level 29 is that case. |
| Total-field legend | “The total field” of the level sources plus the particle | The static part is computed once, in f64, on 6 points per cell, and interpolated bilinearly (`PHYSICS.md` §10). Antennas, waves, and the particle are evaluated per pixel. The legend does not say which part is the grid. |
| Waves legend | “exact retarded fields” | The analytic expressions match the CPU fields. The picture is `f32`, with phases reduced in f64 before upload (`PHYSICS.md` §10). At `c = ∞` the shader drops to the instantaneous dipole and `B = 0`, which is the Newtonian limit. |

The particle-field legend does say, in a following line, that near metal the induced charges are not drawn while the force on a sphere is included and the electrode image force is bounded. The field-line line says the lines show direction only. Both of those limits are stated where the player sees them.

## Display substitutions already filed

These are formula disagreements between the authoritative sample and a picture or a comment. They are not new solver bugs.

- `potential.wgsl` uses `Q / max(r, 10⁻⁴)` for every charge, including a cloud. `Coulomb::sample` uses the uniform-sphere field inside. Finding 2.
- `potential::params` never uploads `External`. Uniform `E` has `φ = −E·x` on the flight and 0 on the map. Finding 1.
- Display electrodes are the centroid point charges. Preview and verify keep the triangle integrals. Finding 5.
- A static antenna stores `φ = 0` while `E = −∇φ` for `φ = n·p/r²`. Finding 3.
- `External::sample` comments that a frozen wave has `B = 0`. `fields` sets `B = k̂ × E / c`, and that pair satisfies Maxwell. Finding 6.
- `panel.rs` comments that the normal field jumps by `±2πσ`. The jump is `4πσ`. Finding 6.
- `beam.rs` comments that the tapered past never becomes superluminal. A brake already at `β ≥ 0.9` does. Finding 11.
- For `c = ∞` the energy caption says the total is conserved. A cup fade does work the total does not book. Finding 12.

## Left unevaluated

Those items are checked in [REMAINDER.md](REMAINDER.md). The coefficient comparison, the metal-sphere fit, the pipe mode, and the radiation integral are in [results/closing_checks.txt](results/closing_checks.txt) and [results/closing_wolfram.txt](results/closing_wolfram.txt). Where a measured error is quoted from `PHYSICS.md` and the workspace tests did not print that table again, this audit treats the table as the project’s own measurement.
