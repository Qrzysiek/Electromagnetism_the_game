# Physics: model, numerics, validation

This document is the authoritative description of the physics that is implemented and how. Each section has a status:
- **design:** planned, not yet in code
- **implemented:** in code, with the file and function named
- **validated:** implemented, with the tests that check it and their measured results

Any change to the physics or numerics must update this file in the same commit.

---

## 1. Units — *implemented* (`crates/physics/src/units.rs`)

Internal units are dimensionless. A world is defined by four SI reference values:

| Quantity | Symbol | Meaning |
|---|---|---|
| length | `L₀` | grid cell size |
| charge | `Q₀` | reference charge |
| mass | `M₀` | reference mass |
| Coulomb constant | `k = 1/(4πε₀)` | set to 1 internally |

Derived units:
- time `T₀ = sqrt(M₀ L₀³ / (k Q₀²))`
- energy `E₀ = k Q₀² / L₀`
- velocity `V₀ = L₀ / T₀`
- the speed of light becomes a dimensionless parameter `c = c_SI / V₀`

`WorldScale` converts between the two. CODATA 2018 constants are used for `k`, `e`, `m_e`, `m_p`; `c` is exact.

In the non-relativistic limit, trajectories depend only on the dimensionless ratio `k q Q / (L E_kin)`. The physical scale enters through `c` (relativity) and, in later stages, through material properties and particle–particle interaction.

## 2. Field sources — *validated* (`crates/physics/src/field.rs`: `Coulomb`, `UniformElectric`)

Fixed charges `Q_i` at positions `r_i`, each a rigid sphere of radius `R_i` with a spherically symmetric charge distribution. Outside the sphere (the only region particles can reach), the shell theorem gives exactly:

```
φ(r) = Σ_i k Q_i / |r − r_i|
E(r) = Σ_i k Q_i (r − r_i) / |r − r_i|³
```

- The sum is evaluated in f64 in a fixed order (the order in the level file). No softening.
- Vector arithmetic uses `glam::DVec3`. Its `dot`, `cross` and `length` are plain multiply/add plus IEEE `sqrt`, with no fused multiply-add (checked in the source; `mul_add` is never called). Results are therefore deterministic across platforms.
- `UniformElectric` (potential `−E·x`) exists for tests with analytic solutions (T6).
- Test particles are rigid, spherically symmetric and **non-polarizable** (a stated approximation, since a conducting microsphere would feel image forces). Under these assumptions, the force on a test particle of radius `a` equals `q E(center)` exactly, because the force between two non-overlapping spherically symmetric distributions equals the point-charge force.
- Contact condition: `|x − r_i| ≤ R_i + a`, which means the particle is lost.

## 2.2 Magnetic sources — *validated* (`crates/physics/src/magnetic.rs`, `field.rs::LevelField`)

**Units.** Source strengths are given directly in field units, i.e. already multiplied by `μ₀/4π`:
- a dipole by `μ = μ₀ m / 4π`,
- a coil by `κ = μ₀ I / 4π`.

With `k = 1/(4πε₀) = 1` internally, `μ₀/4π = 1/c²`. Specifying the sources in field units keeps magnets meaningful in Newtonian levels (`c = ∞`), while the SI relations hold exactly. The unit of B is `M₀ / (Q₀ T₀)` (from `F = q v × B`).

**Sources.** All three are exact closed forms:

| Source | Field | Obstacle |
|---|---|---|
| Uniformly magnetized sphere | Outside the sphere, exactly the point-dipole field `B = μ (3(m̂·r̂) r̂ − m̂) / r³` | sphere |
| Circular coil | Simpson et al. (NASA/TM-2001-209961). With `α² = (a−ρ)² + z²`, `β² = (a+ρ)² + z²`, `m = k² = 4aρ/β²`: `B_z = 2κ/(α²β) [(a² − ρ² − z²) E + α² K]`, `B_ρ = 2κ z β / (α² ρ) g(m)`, with `g(m) = (1 − m/2) E − (1 − m) K` | torus (wire) |
| Polygonal coil | Sum of exact segment fields `B = κ (r_a × r_b)(|r_a| + |r_b|) / (|r_a| |r_b| (|r_a| |r_b| + r_a·r_b))` | capsules (wires) |

**Circular-coil numerics.**
- `K` and `E` come from the arithmetic–geometric mean (DLMF 19.8), converging quadratically to full precision.
- `1 − m = α²/β²` is passed separately, which avoids cancellation near the wire.
- `g(m)` is `O(m²)`. Evaluated directly it cancels catastrophically near the axis: the first version gave an 85 % error at a point 10⁻¹⁶ off the axis. For `m < 0.05` it is therefore summed as a power series in `m` (DLMF 19.5.1–2), whose `m⁰` and `m¹` terms cancel analytically.

**The 2D slice.** Charges, dipoles perpendicular to the plane (`m ∥ ẑ`) and coils lying in the plane give, at `z = 0`:
- **E** exactly in the plane,
- **B** exactly along `ẑ`, with the in-plane components exactly 0, not just small.

So `q v × B` keeps the particle in the plane (test M5). A straight wire perpendicular to the plane would make B in-plane and push the particle out; such wires are therefore for 3D levels only.

**Energy.** B does no work, so the conserved energy `W` of §4 is unchanged. The forbidden-region map stays exact with magnets.

## 2.3 External fields: stray fields and plane waves — *validated* (`crates/physics/src/external.rs`)

Fields whose sources are outside the arena (disturbances, SPEC §3). `LevelField::external` adds them to the level's own sources. Both kinds are exact solutions of the vacuum Maxwell equations, so they add no approximation:

- **Uniform stray fields** `E` (in the plane) and `B_z`. They are static, and the electric part has potential `−E·x`.
- **Plane waves** travelling in the plane, linearly polarized in the plane:
  ```
  E = E₀ ê cos(ω (t − k̂·x / c) + φ),   B = k̂ × E / c,   ê = ẑ × k̂
  ```
  - `B` is along `ẑ` and `E` lies in the plane, so the 2D slice stays exact (test W3).
  - `k̂·x/c` is multiplied by `ω`, so for `c = ∞` the wave becomes a spatially uniform field oscillating in time with `B = 0`. This is the exact Newtonian limit, used for "mains hum" levels.
  - `ω = 0` gives a static uniform field `E₀ ê cos φ`, with its potential.
- A wave travelling out of the plane (k̂ = ±ẑ) with `E` in the plane would have `B` in the plane and push particles out of it. It is therefore not offered in 2D levels.

**Unit consistency.** `B = k̂ × E / c` follows from Faraday's law `∇×E = −∂B/∂t`. That law has this form in every unit system with force `q(E + v×B)`, which is the one used here (§2.2). A unit test checks all four vacuum Maxwell equations for the implemented wave by central differences: `∇·E = ∇·B = 0`, `∇×E = −∂B/∂t`, `∇×B = ∂E/∂t / c²`.

**Levels.** A level lists `disturbances`, each one realization: uniform `E`, `B_z`, and any number of waves. Every shot is flown under every disturbance, and a setup solves the level only if every flight arrives, verified.
- A "random" disturbance, such as the phase of mains hum, is therefore represented by a finite, fixed set of realizations. The game stays deterministic, and every flight is verified individually.
- Robustness between the listed realizations (for example at intermediate phases) is not guaranteed. It is a design parameter of the level (how many realizations it lists).

**Diagnostics.** In time-dependent fields the energy of §4 is not conserved (the wave does work). The energy diagnostic is then reported as NaN rather than a misleading number (`FieldSolver::is_static`). Validation instead uses the exact invariants of motion in a plane wave (test W2).

## 2.4 Antennas: oscillating electric dipoles — *validated* (`crates/physics/src/antenna.rs`)

A small antenna is modelled as a point electric dipole `p(t) = p₀ cos(ωt + φ)` inside a solid sphere (an obstacle). Its fields are the exact retarded fields (Jackson §9.2), at `t_r = t − r/c`, with `n` the unit vector from the dipole and units `k = 1`, `μ₀/4π = 1/c²`:

```
E = [3n(n·p) − p]/r³ + [3n(n·ṗ) − ṗ]/(c r²) + [n × (n × p̈)]/(c² r)
B = (1/c²) [ṗ × n / r² + p̈ × n / (c r)]
```

- This is an exact solution of the vacuum Maxwell equations outside the source point, covering the near (quasi-static), induction and radiation zones. Nothing is truncated.
- For `c = ∞` it becomes the electrostatic dipole field of the instantaneous `p(t)`, with `B = 0`.
- A dipole in the plane (`p ⊥ ẑ`) gives, at `z = 0`, `E` exactly in the plane and `B` exactly along `ẑ`, so the 2D slice stays exact (test A3).
- The radiated power is `⟨P⟩ = p₀² ω⁴ / (3c³)` (Larmor). It is supplied by the generator that drives the antenna, which is outside the model.

**Levels.** Antennas run at the level's RF generator frequency `physics.rf_omega`, or at their own `omega` if they have one. Levels may offer the player a list of frequencies (`limits.antenna_omegas`). All antennas start in phase at t = 0 (φ = 0). An antenna element has amplitude `value` (a negative value is the opposite phase) and orientation `angle_deg`. The player can choose orientations of 0°, 45°, 90° and 135°.
- A shot may have a launch time `t₀` (`launch.time`). The flight then sees the time-dependent sources at `t₀ + t`, implemented as `LevelField::time_offset`, which is exact.
- Identical particles launched at different times follow different trajectories only in time-dependent fields. This is the basis of the RF levels.

## 2.5 The particle's own field: Liénard–Wiechert — *validated* (`crates/physics/src/lienard.rs`)

The exact field of a point charge on an arbitrary world line (Jackson §14.1). With `R = x − r(t_r)`, `n = R/|R|`, `β = v/c`, `κ = 1 − n·β` at the retarded time `t_r`, where `|x − r(t_r)| = c (t − t_r)`:

```
E = q (n − β)(1 − β²) / (κ³ R²)  +  (q/c) n × ((n − β) × β̇) / (κ³ R),     B = n × E / c
```

- The first term is the velocity (generalized Coulomb) field; the second, proportional to the acceleration, is the radiation field.
- The retarded time is the root of `g(t_r) = c(t − t_r) − |x − r(t_r)|`, which is strictly decreasing for `|v| < c`. It is found by bracketing and bisection to full precision.
- For a computed flight the world line is interpolated between the preview's dense-output points: cubic Hermite for position, linear for velocity and acceleration (the acceleration comes from the Lorentz force). Before launch and after the end, the particle is continued with uniform motion. A real launch or stop would itself radiate; this is a stated simplification of the view.
- **Use.** Visual only: the "particle field" map (§10). The particle's own field does not act on the particle except through radiation reaction (§3.1), the consistent classical treatment of that self-interaction.

## 3. Equation of motion — *validated* (`crates/physics/src/dynamics.rs`)

State `y = (x, p)`, with:

```
γ     = sqrt(1 + |p|² / (m² c²))
dx/dt = v = p / (γ m)
dp/dt = q (E(x, t) + v × B(x, t))
```

`|v| = |p| c / sqrt(m²c² + |p|²) < c` holds for any finite `p`. The speed limit is structural, not enforced after the fact.

Implementation details:
- `c = ∞` is allowed and gives exact Newtonian mechanics (`1/(mc)² = 0`, so `γ = 1`). The non-relativistic tests use this rather than a large finite `c`.
- The integrator works with the scaled state `(x, p/p_ref)`, where `p_ref = |p₀|` (or `m` if the particle starts at rest). One tolerance is then meaningful for positions (grid units) and momenta alike.
- Kinetic energy is evaluated as `(γ−1)mc² = p²/(m(γ+1))`, which avoids cancellation at low speed.
- **Speed limit in floating point:** `|v| = |p|c/sqrt(m²c² + p²)` never exceeds `c`. For `γ ≳ 10⁷`, the rounded quotient can equal `c` exactly. The strict inequality `|v| < c` is guaranteed and tested for `γ < 10⁷`. Nothing uses `v` as state, so this has no effect on the trajectory.

## 3.1 Radiation reaction — *validated* (`dynamics.rs::radiation_reaction_force`)

An accelerated charge radiates, and the energy comes from the particle. Levels can include this (`physics.radiation_reaction`) through the Landau–Lifshitz force (Classical Theory of Fields §76), written for force `q(E + v×B)`, `k = 1`, `μ₀/4π = 1/c²`:

```
f = (2q³/3mc³) γ [DE/Dt + v × DB/Dt]
  + (2q⁴/3m²c⁴) [c E×B + c B×(B×v) + E (v·E)/c]
  − (2q⁴/3m²c⁵) γ² v [(E + v×B)² − (E·v)²/c²]
```

- **Why LL.** The Lorentz–Abraham–Dirac equation has runaway and pre-acceleration solutions. LL is its reduction of order, consistent to first order in `τ₀ = 2q²/(3mc³)`, with no runaways. It is valid while the reaction force is small against the Lorentz force; each trajectory reports `max |F_RR| / |F_L|`. Shipped levels require it below 0.05 (`level/tests/levels.rs`), and the game flags it.
- **The derivative term.** `D/Dt = ∂/∂t + v·∇` is a central difference of the field along the world line over 1e-5 cells. Its relative error (about 1e-10) is far below the O(τ₀) accuracy of LL itself. It is deterministic.
- **Energy bookkeeping.** The work of `f` is integrated as a 7th ODE component. So `kinetic + potential + radiated = const` holds to integration accuracy, and the energy diagnostic of §4 checks exactly this.
- **Consistency.** `c = ∞` gives no reaction. With the flag off, the dynamics are bit-identical to before (6-component state).
- **Levels.** A level with radiation reaction must need it: with the flag off, its reference solution must fail (test `radiation_levels_need_radiation`).

## 4. Conserved quantities (diagnostics only) — *validated* (`crates/physics/src/trajectory.rs`)

Static fields (no plane waves with ω ≠ 0, §2.3):
- Total energy `W = γ m c² + q φ(x)` is conserved. The code tracks `W − mc² = (γ−1)mc² + qφ`, which is finite for `c = ∞`. `Trajectory::energy_max_abs_error` is the largest `|W(t) − W(0)|` over the accepted steps and the final event state.
- A single central charge also conserves angular momentum `L = x × p`, which holds relativistically too.

They are computed along every trajectory and reported. They are **never** used to project or correct the state.

## 5. Integrator — *validated* (`crates/physics/src/integrator/dop853.rs`)

- Dormand–Prince 8(5,3) (DOP853, Hairer, Nørsett & Wanner, *Solving ODEs I*, §II.10), with a 7th-order dense-output interpolant.
- Mixed error control on each component: `err ≤ atol + rtol·|y|`.
- Two tolerance levels:
  - **preview** (fast, live editing)
  - **verify** (at least 100× tighter)
- **Chosen values (M2):** preview `rtol = atol = 1e-10`, verify `1e-12`.
  - Benchmark (`cargo bench -p physics --bench preview`): one trajectory through a 50-charge, 40×30-cell level takes 175 µs at tol 1e-8 (41 steps) and 356 µs at 1e-11 (90 steps). Both are far inside the 16 ms frame budget, so even the preview can use a tight tolerance.
  - Calibration (T10): after one Kepler orbit the global position error is about 30–80 × tol. At 1e-12 that is below 1e-10 grid cells, just above the rounding floor found in M1.

### 5.1 Why our own port (M1 decision)

| Candidate | Result |
|---|---|
| `ode_solvers` 0.6.2 (DOP853) | Faithful port, but it offers no per-step access to the interpolant. Its callback only sees step endpoints, so events cannot be located on the dense output. Kept as a **test-only cross-check**. |
| `russell_ode` 3.3.1 | Needs OpenBLAS through MSYS2 on Windows, a heavy system dependency for a 6-variable non-stiff problem. Rejected. |
| `diffsol` 0.16.2 | Good step/event API, but its only explicit method is 5th order (Tsit5). At the 1e-12 tolerances we need, that means several times more steps. Its root finding looks only at sign changes at step endpoints, which misses grazing events. Rejected. |
| **Own port of Hairer's `DP86CO`** | Step-at-a-time API exposing the dense output; a hard per-step cap `h_cap`. **Adopted.** |

How the port is made trustworthy:
- **Coefficients** are generated mechanically from Hairer's original `dop853.f` (`scripts/gen_dop853_coefficients.py`; the source's SHA-256 is recorded in the generated file). There is no hand transcription.
- **Coefficient table vs. order conditions** (unit tests, independent of how the table was made): row sums `Σ_j a_sj = c_s`; quadrature `Σ b_i c_i^(q−1) = 1/q` for q = 1…8; the 3rd-order estimator is exact to degree 2; the 5th-order estimator annihilates polynomials to degree 4.
- **Step controller, error norm, HINIT and dense output** follow the Fortran operation by operation, including summation order.
- **Documented deviations:** forward time only; a step reaching `t_limit` ends exactly on it but never exceeds `h_cap`; no stiffness detector.

### 5.2 Measured results (`cargo test -p physics --test integrator -- --nocapture`, `cargo run --release -p integrator_eval`)

**Cross-check against `ode_solvers` on the Arenstorf orbit** (Hairer's standard test, periodic with `T = 17.0652165601579625588917206249`). The accepted and rejected step counts are **identical at every tolerance**:

| tol | error after one period (ours) | error (ode_solvers) | \|ours − theirs\| | steps (acc+rej) | f-evals |
|---|---|---|---|---|---|
| 1e-6 | 6.91e-3 | 6.91e-3 | 6.5e-10 | 87 | 1035 |
| 1e-8 | 8.43e-5 | 8.43e-5 | 7.5e-11 | 147 | 1736 |
| 1e-10 | 8.55e-7 | 8.56e-7 | 6.4e-10 | 235 | 2785 |
| 1e-12 | 7.70e-10 | 7.91e-10 | 9.9e-11 | 357 | 4249 |
| 1e-14 | 4.05e-10 | 9.93e-11 | 5.0e-10 | 535 | 6436 |

The two codes round differently (`ode_solvers` forms `y + Σ(h·a)k`, Hairer and we form `y + h·Σ a k`). The orbit's close lunar passage amplifies this to about 1e-10. Below tol ≈ 1e-12 both hit the same rounding floor (about 1e-10 to 1e-9 for this problem), so tighter tolerances buy nothing there. Speed is comparable (about 250 µs for the 1e-12 run on both).

**Convergence order (fixed steps):**

| problem | errors | observed order |
|---|---|---|
| Kepler, e = 0.5, one period, N = 128 → 256 | 1.68e-10 → 7.57e-13 | 7.80 |
| harmonic oscillator, t = 10, N = 20 → 40 | 2.34e-9 → 8.66e-12 | 8.08 |
| dense output, local error, h = 0.8 → 0.4 → 0.2 | 2.10e-7 → 8.14e-10 → 3.17e-12 | 8.01, 8.00 |

The theory predicts 8 for the global error and 8 for the local error of the 7th-order interpolant.

**Relativistic hyperbolic motion** (prototype of T6; `m = c = qE = 1`, from rest; exact `x = sqrt(1+t²) − 1`), tol 1e-12:

| t_end | γ | relative error in x | steps |
|---|---|---|---|
| 1 | 1.4 | 1.6e-13 | 14 |
| 10 | 10 | 4.7e-15 | 33 |
| 1000 | 1000 | 1.9e-14 | 57 |

### 5.3 Why not Boris

- Why not Boris: with `B = 0`, Boris is leapfrog. Its good long-term energy behaviour relies on a fixed step, and an adaptive step (unavoidable with close encounters) breaks that. It may be revisited for magnetic fields.

## 6. Events — *validated* (`crates/physics/src/events.rs`, `trajectory.rs`)

Events are located on the dense-output polynomial of each accepted step:
- `g_i(t) = |x(t) − r_i| − (R_i + a)` crosses zero: collision with charge `i`.
- `x(t)` enters region B: success.
- `x(t)` leaves the world bounds: lost (a game rule).
- `t > t_max`: timeout (a game rule).

All event functions are signed distances (sphere, torus, capsule, box, world bounds), which are 1-Lipschitz in position. Along the trajectory, `|dg/dt| ≤ |v| ≤ v_max`. On an interval `[a, b]` with `g(a), g(b) > 0`:

```
min g ≥ (g(a) + g(b) − v_max (b − a)) / 2
```

- If this bound is positive, the interval is **certified crossing-free**.
- Otherwise it is bisected, left half first, so the **earliest** crossing is found even when `g` dips below zero and comes back between samples.
- The crossing time is resolved to adjacent floating-point numbers: the result is the first representable time at which the interpolated `g ≤ 0`.
- In a typical step all distant obstacles are certified from the two endpoint values alone, so the check costs almost nothing.
- The earliest event over all functions ends the trajectory. Simultaneous events are ranked obstacle, then bounds, then detector.
- If an interval can no longer be subdivided and is still not certified (a graze within `v_max · ulp(t)` of the surface), it is treated as a miss. Such a result is marginal by construction (§7).

**Speed bound:** `v_max` is the largest speed sampled on the interpolant at 5 points of the step, times 1.25, capped at `c`. This is the one non-rigorous element: it assumes the speed does not rise by more than 25 % between those samples within one error-controlled step. A rigorous alternative is `v_max = c`, which is looser and slower. It remains available if a counterexample is ever found.

The step size is **not** otherwise limited near obstacles; the certification replaces the displacement cap considered in the design. `Trajectory::closest_sampled` holds, per obstacle, the smallest surface distance among the points evaluated. That is an upper bound on the true closest approach; exact margins are computed in M3.

## 7. Outcome verification — *validated* (`crates/physics/src/verify.rs`, margins in `trajectory.rs`)

Every flight that matters is computed at two tolerances: preview (1e-10) and verify (1e-12).

**Margins.** For each event boundary (every obstacle, the world bounds, the detector), the runner computes the minimum of its event function over the flight:
- **Boundaries not crossed:** the closest approach (positive). Wherever the certified lower bound of §6 drops below `2 × MARGIN_SAFE` (0.5 cells), it is refined on the dense output with 32 samples plus golden-section search. Larger margins are only lower-bounded.
- **The boundary that ended the flight:** the penetration depth (negative) of the *continued* trajectory, integrated past the event as if the boundary were absent. Integration continues until the minimum is reached, or the depth exceeds `MARGIN_SAFE` (0.25 cells), or it breaks down at a point charge (depth −∞).
- A first version used only the step containing the event. That made the depth depend on where steps happened to end, and flagged clear hits as marginal. The continued-trajectory definition is a property of the trajectory alone.
- **Time:** for flights ending in an event, the remaining time `t_max − t_event`.
- **Stopping the continuation:** it goes on while the minimum of the latest step lies at the step's right end, judged from *where* the minimum lies. A first version compared the event function at the step end, evaluated once from the integrator state and once from the dense output. Those differ by rounding, so the continuation sometimes stopped after one step. Clear exits through the map edge were then flagged marginal (depths −0.57 vs −0.11). Regression test: `clear_exit_through_bounds_is_verified`.

**Classification.**
- **Verified:** both runs give the same outcome, and every margin `m` below `MARGIN_SAFE` satisfies `|m| > SAFETY · |m_preview − m_verify| + FLOOR`, with `SAFETY = 10` and `FLOOR = 1e-9` cells (`FLOOR · t_max` for time).
- The difference between the runs is a direct estimate of the preview run's error. Since error ∝ tol (T10), it overestimates the verify run's error about 100×. The criterion is therefore conservative by roughly three orders of magnitude.
- **Marginal** otherwise, with the reason: outcome mismatch, which boundary, its margin, and its error estimate.

**Measured** (`cargo test -p physics --test verification -- --nocapture`). Relativistic (`c = 3`) Coulomb flyby of a sphere of radius `r_min(1 + δ)`:

| δ | gap / depth | outcome | status |
|---|---|---|---|
| ±1e-3 | ∓2.27e-3 | hit / miss | verified |
| ±1e-6 | ∓2.27e-6 | hit / miss | verified |
| ±1e-11 | ∓(1.8 … 2.7)e-11 | hit / miss | marginal (error estimate 2.6e-12) |
| ±1e-13 | 4.2e-12, 4.7e-12 | miss, miss (the δ = +1e-13 hit is misclassified) | marginal |

The last row shows the mechanism working as intended. At δ = 1e-13 the integration error exceeds the physical gap, so the outcome is not trustworthy, and the result is flagged instead of reported.

Accuracy of the margins against the analytic gap `r_min δ`: relative error 2e-10 (δ = 1e-2), 2e-8 (1e-4), 2e-6 (1e-6). The absolute error stays about 5e-12 cells.

## 8. Neglected effects and their control — *implemented*

| Effect | Status in Stage 1 | Control |
|---|---|---|
| Radiation (Larmor/Liénard) | neglected, unless the level includes radiation reaction (§3.1) | The radiated energy `∫P dt`, with `P = (2/3) k q² γ⁶ (a² − (v × a)²/c²) / c³`, is computed as a diagnostic (`Trajectory::radiated_energy`, trapezoidal rule over accepted steps; exactly 0 for `c = ∞`). Flagging levels above a threshold is part of the generator (M5). |
| | | **Consistency requirement (found in M4):** radiated fraction per close pass ≈ `r_cl / r`, with `r_cl = k q²/(mc²)` the particle's classical radius. The first demo levels used `q = m = 1`, `c = 1.5…5`, giving `r_cl` = 0.04–0.44 cells. Their neglected radiation was 1e-4 to 2 × T₀, so they were physically inconsistent. Fix, as in real accelerators: weakly charged particle (`q = 10⁻⁶`), strongly charged electrodes (`Q ~ 10⁶`). Trajectories depend only on `qQ/m` and `c`, so the puzzles are unchanged. Neglected radiation of the reference flights is now 1e-16 … 1.7e-11 × T₀. Test `level/tests/levels.rs` requires `< 1e-10` for every shipped level; the game flags player setups above that. |
| Magnetic field of the moving particle acting on others | n/a (single particle) | Stage 6 (Darwin) |
| Polarization of test particles | neglected (non-polarizable by assumption) | documented assumption |
| Recoil of fixed charges | none (held fixed by definition) | game rule |
| Radiation reaction and scattering in waves (Thomson scattering) | neglected | The same radiated-energy diagnostic and the same `< 1e-10 T₀` requirement apply to flights under disturbances (all flights of all shipped levels are checked). |

## 9. Validation tests — *validated*

Run with `cargo test -p physics --test validation --test properties -- --nocapture --test-threads=1`. Unless stated otherwise, tolerance is 1e-12.

**How the thresholds were set.** Criteria were fixed before measuring, from the target accuracy. The exception is T10: its first version used a 0.8–1.2 slope window guessed before measuring. That version fitted the pre-asymptotic regime of an e = 0.5 orbit (11–43 steps per orbit) and measured 0.78. The theory predicts slope → 1 asymptotically: Hairer's error estimate `err5²/sqrt(err5² + 0.01·err3²)` scales as h⁹, like the local error of the 8th-order solution. The test now measures in the asymptotic regime (e = 0.2, tol 1e-9…1e-14). It additionally guards `error < 100·tol` for both orbits over tol 1e-6…1e-13.

| # | Test | Analytic reference | Pass criterion | Measured |
|---|---|---|---|---|
| T1 | Energy conservation, 20 random charges in 3D, 4 layouts, `c = ∞` and `c = 5` | `W` constant | `max|ΔW| / T₀ < 1e-10` | 1.1e-12 … 1.4e-11 |
| T2 | Rutherford, non-relativistic, 4 geometries (θ from 0.22 to 2.75 rad) | `tan(θ/2) = κ/(m v∞² b∞)` | relative error < 1e-8 | 8e-13 … 4.9e-12 |
| T3 | Relativistic Coulomb scattering, repulsive and attractive, `v/c` = 0.45 … 0.98 | `θ = |π − (2/Γ) arccos(−B₀/A)|` | relative error < 1e-8 | 9e-14 … 2.5e-12 |
| T4 | Kepler orbit, e = 0.5, 10 orbits | `T = 2π sqrt(m a³/|κ|)`, closure, `L` const. | period < 1e-9, closure < 1e-8, `ΔL/L` < 1e-10 | period 3.2e-11, closure 3.5e-9, `ΔL/L` 1.9e-11, `ΔE/E` 5.5e-11 |
| T5 | Relativistic precession, `c` = 10 and 3 (`v/c` up to 0.5), 20 revolutions | `Δφ = 2π(1/Γ − 1)`; orbit `r(φ)` | relative error < 1e-6; `r(φ)` < 1e-9 | 1.4e-10, 5.9e-11; `r(φ)` 3.2e-10, 8.3e-10 |
| T6 | Uniform `E`, `p⊥` = 0, 0.5, 3, up to γ ≈ 1000 | `x∥ = (sqrt(ε₀² + (qEct)²) − ε₀)/qE`, `x⊥ = (p⊥c/qE) asinh(qEct/ε₀)` | x, y < 1e-11; p < 1e-12 | x ≤ 4.2e-13, y ≤ 6.2e-13, p ≤ 6.7e-16 |
| T7 | Speed limit, 64 random fields, `|p|` from 1e-3 to 1e7 | `|v| ≤ c`; `< c` for γ < 1e7 | always | holds |
| T8 | Plane of symmetry, 64 random in-plane configurations | `z = p_z = 0` | exactly ±0 | holds |
| T9a | Grazing a sphere on a straight line, closest approach `d` | hit iff `R > d` | correct for `R = d(1 ± δ)` | correct down to δ = 1e-14 |
| T9b | Grazing in a Coulomb orbit (`c = ∞`, 3; repulsive, attractive) | `r_min = 1/(A + B₀)` | `r_min` < 1e-10; classification at δ = 1e-6, 1e-8 | `r_min` 3e-13 … 2e-12; correct |
| T10 | Convergence with tolerance | error ∝ tol | asymptotic slope in [0.9, 1.1]; error < 100 tol | slope 0.957; error/tol 28 … 47 (e = 0.2), ≤ 80 (e = 0.5) |
| T11 | Determinism: FNV-1a hash of every bit of the reference trajectories of all `levels/*.json` (both tolerances, samples, margins, step statistics) | – | identical across runs; equal to `levels/golden_hashes.json` (generated on Windows) on Windows and Linux CI | holds |

### Magnetism tests (`cargo test -p physics --test magnetism -- --nocapture --test-threads=1`)

| # | Test | Reference | Criterion | Measured |
|---|---|---|---|---|
| M1a | Circular coil (tilted, off-centre) on its axis | `B = 2πκ a² / (a² + z²)^(3/2)` | < 1e-14 | ≤ 7.9e-16 |
| M1b | Circular coil off the axis, including its plane, inside and outside | Biot–Savart, trapezoidal rule with 20000 nodes | < 1e-11 | 1.1e-14 |
| M1c | Coil far field | Dipole `μ = κπa²`, correction `O((a/r)²)` | < 3(a/r)² | 1.1e-4 at r = 200, 1.1e-6 at r = 2000 (scales as (a/r)²) |
| M1d | Rectangular coil | Segment-wise Gauss–Legendre quadrature | < 1e-12 | 2.7e-15 |
| M1e | Regular N-gon converges to the circle | order 2 | slope 1.9–2.1 | 2.00 |
| M1f | Dipole: ∇·B = 0, ∇×B = 0 | central differences | < 1e-7 of scale | ≤ 3.7e-8 (difference error) |
| M1g | In-plane sources at z = 0 | `B_x = B_y = 0` | exactly ±0 | holds (200 random points) |
| M2 | Relativistic cyclotron, γ = 1, 1.12, 3.16, 10 turns | `r = p/(|q|B)`, `T = 2πγm/(|q|B)` | r < 1e-10, closure < 1e-9, |p| < 1e-10 | r ≤ 3.1e-11, closure ≤ 9.7e-10, |p| ≤ 3.1e-11 |
| M3a | E×B from rest, Newtonian | cycloid `x = (E/Bω)(ωt − sin ωt)`, `y = (E/Bω)(1 − cos ωt)` | < 1e-11 | 5.5e-13 |
| M3b | E×B from rest, relativistic, E/cB = 0.3, 0.6, 0.9 | drift `v_d = E/B`; cusp period `2π m γ_d³ / (qB)` (Lorentz boost to the frame where E' = 0, B' = B/γ_d) | < 1e-9 | period ≤ 6.3e-13, drift ≤ 7.4e-16 |
| M4 | Energy with charges, dipoles and a coil; c = ∞ and 5 | `W` constant | < 1e-10 | ≤ 2.7e-11 |
| M5 | Plane of symmetry with dipoles and a coil (32 random configurations) | `z = p_z = 0` | exactly ±0 | holds |
| M6 | Straight flight into a ring wire and a straight wire | known hit time | < 1e-12 | ≈ 1e-15 |

Note on M2: the first threshold for |p| drift (1e-12) was stricter than the integration tolerance (also 1e-12). DOP853 does not preserve |p| exactly, so the criterion is the same 1e-10 as for energy in T1.

### External-field tests (`cargo test -p physics --test waves -- --nocapture --test-threads=1`)

| # | Test | Reference | Criterion | Measured |
|---|---|---|---|---|
| W0 | Plane wave (`c` = 3, oblique): Maxwell equations; vector potential | central differences: `∇·E = ∇·B = 0`, `∇×E = −∂B/∂t`, `∇×B = ∂E/∂t/c²`; `E = −∂A/∂t`, `B = ∇×A` | < 1e-7 of the field scale (difference error) | holds |
| W1 | Newtonian particle in a uniform field `E₀ cos(ωt + φ)`, 4 phases, oblique launch | `v = v₀ + (a/ω)[sin(ωt+φ) − sin φ] ê`, `x = x₀ + v₀t − (a/ω)[(cos(ωt+φ) − cos φ)/ω + t sin φ] ê`, `a = qE₀/m` | `|Δx|/(a/ω²)`, `|Δv|/(a/ω)` < 1e-9 | ≤ 4.5e-12, ≤ 1.7e-13 |
| W2 | Relativistic particle in a plane wave (`c` = 2), a₀ = qE₀/(mcω) = 0.1, 1, 2, 3; 40 periods; one case with p_z ≠ 0 | Exact invariants (Landau–Lifshitz §47–48): light-front momentum `γmc − p·k̂` and canonical transverse momentum `p_⊥ + qA`. With the mass shell they fix `p` as a function of the wave phase, so this checks the analytic solution for `p`. | drift < 1e-9 of `max(mc, |p|max)` | 3.3e-13 … 3.6e-12 (|p| up to 6 mc) |
| W3 | Plane of symmetry: charge, dipole, in-plane wave and uniform stray fields | `z = p_z = 0` | exactly ±0 | holds |
| W4 | Energy with uniform stray E and B_z and a static (ω = 0) wave term, `c = ∞` and 4 | `W` constant, with potential `−E·x` | `max|ΔW|/T₀` < 1e-10 | 2.6e-12 |

The thresholds are the ones of the corresponding T and M tests (1e-9 for trajectory comparisons at tolerance 1e-12, 1e-10 for conserved quantities), fixed before measuring.

### Antenna tests (`cargo test -p physics --test antenna -- --nocapture --test-threads=1`, and unit tests in `antenna.rs`)

| # | Test | Reference | Criterion | Measured |
|---|---|---|---|---|
| A0 | Maxwell equations at r = 1.4, 5 and 42 (λ = 6.3): near, induction and radiation zones | central differences, step `2.5e-5 · min(r, λ)` | residual < 1e-6 of the compared terms | ≤ 6.6e-8; residuals scale as h² (16× per 4× step), so they are difference error |
| A0' | `c = ∞` limit | instantaneous static dipole field; `c = 1e9` | < 1e-8 | holds |
| A1 | Time-averaged Poynting flux `(c²/4π) E×B` through spheres R = 0.2 λ, λ, 20 λ (general p₀ direction, off-origin dipole) | Larmor `p₀²ω⁴/(3c³)`, the same at every R (near-field terms are reactive) | < 1e-10 | 8e-16, 2e-15, 5e-15 |
| A2 | Static limit ω = 0 | charges ±Q at ±d/2, `Qd = |p₀|`; the difference is the octupole term `O((d/r)²)` | < (d/1.8)² | 9.7e-6 (d = 1e-2), 9.7e-8 (d = 1e-3): scales as d² |
| A3 | Plane of symmetry: two in-plane antennas, a charge, a launch-time offset | `z = p_z = 0` | exactly ±0 | holds |

A1 is the strongest check. Its quadratures (Gauss–Legendre in cos θ, uniform in φ and in time) are exact for the trigonometric polynomials involved, so any error in the field formulas, including the near and induction terms, would show up directly.

### Radiation tests (`cargo test -p physics --test radiation_reaction --test lienard -- --nocapture --test-threads=1`)

| # | Test | Reference | Criterion | Measured |
|---|---|---|---|---|
| R1 | Synchrotron damping in uniform B, γ₀ = 1.05, 2, 10, 600 time units (\|p\| falls 8–40×) | LL reduces exactly to `du/dt = −(κ/m) u sqrt(1 + u²)`, `κ = 2q⁴B²/(3m²c³)`, so `u(t) = 1/sinh(asinh(1/u₀) + κt/m)` | \|Δu\|/u₀ < 1e-9; energy balance < 1e-10 | ≤ 8.5e-13; balance ≤ 3.2e-12 |
| R2 | Coulomb scattering at c = 5, 10, 20 | LL work = −∫P_Liénard dt + ΔE_Schott, `E_S = τ₀ m γ⁴ v·a`, up to O(τ₀ω) | relative difference < τ₀ω | 5.0e-9, 6.5e-9, 7.1e-9 (τ₀ω = 7e-3 … 1e-4) |
| R3 | `c = ∞` | no reaction | bit-identical to the flag off | holds |
| L1 | Uniformly moving charge, v/c = 0.05 … 0.95 | Heaviside field from the present position | < 1e-12 | 2.6e-15 |
| L2 | Slowly oscillating charge (Aω/c = 5e-4, 5e-5), radiation zone | oscillating-dipole field `p₀ = qA` (independent code, §2.4) | < 20 Aω/c | 2.8e-5, 2.8e-6: linear in A as expected |
| L3 | Circular motion at v = 0.8 c, near and far points | vacuum Maxwell equations, central differences | < 1e-6 | ≤ 1.7e-7 (difference error) |

Note on R2: the first version compared LL work with `∫P dt` alone. It found a difference of 1e-4 at every c, falling as 1/distance². That is the Schott energy of the finite start and end points; with it included the agreement is 5–7e-9.

### Analytic reference for T3 and T5 (relativistic Coulomb problem)

Take potential energy `κ/r` with `κ = k q Q`, conserved energy `W = γmc² + κ/r`, angular momentum `L`, and `u = 1/r`. The orbit equation is:

```
u'' + Γ² u = −κ W / (L² c²),        Γ² = 1 − κ² / (L² c²)   (requires Lc > |κ|)
```

Its solution, with `φ = 0` at the periapsis:

```
u(φ) = A cos(Γφ) + B₀
B₀   = −κ W / (L² c² Γ²)
A    = sqrt( (W² − m²c⁴) / (L² c² Γ²) + B₀² )
```

- **Scattering (`W > mc²`):** the asymptotes satisfy `cos(Γφ∞) = −B₀/A`, so the deflection angle is `θ = |π − (2/Γ) arccos(−B₀/A)|`. As `c → ∞`, this reduces to the Rutherford formula of T2.
- **Bound orbit (`κ < 0`, `W < mc²`):** successive periapsides are separated by `2π/Γ`, so the precession per revolution is `Δφ = 2π(1/Γ − 1)`.

## 10. Visualisation of time-dependent fields (visual only) — `crates/game/src/radiation.rs`

- **Waves** map: the fields of antennas and plane waves at the animation's lab time.
- **Particle field** map: the Liénard–Wiechert field of the active flight at the animation's flight time. Optionally the radiation part only (the acceleration term, which falls as 1/R).
- **Total** map: everything at once at the animation time. That is the level's static sources, antennas, waves and disturbances (the full level field of the flight), plus the particle's own field: retarded for finite c, the instantaneous Coulomb field for c = ∞. Inside source bodies nothing is drawn. With beams (SPEC §3) the particle term becomes the sum over the beam.
- **Colour and range.** The colour is B_z (signed) or |E|. The dynamic range is adjustable up to 14 decades. In most levels the particle's field is more than 10 orders of magnitude weaker than the electrodes' (q ~ 1e-6 against Q ~ 1e6, chosen to keep radiation negligible, §8). So in the total view it only shows at a large range. At about 12 decades the map becomes nearly sign-only, and the particle's contribution is visible mainly near the zero lines of the other sources. That is the true proportion, not a display artefact. In radiation levels (q = 1) and, later, with beams, the particle's field is comparable to the rest and shows directly.
- Both are computed in f64 on the CPU with the tested physics code (§2.3–2.5) on a 5-pixels-per-cell texture, which is bilinearly filtered, plus optional E arrows. Measured cost: at 5 px/cell, median 2.1 ms per frame for the particle field, 0.3 ms for waves and 4.4–5.0 ms for the total field (measured in release builds).
- **Anti-aliasing.**
  - **CPU maps.** The field is evaluated at the sample points, and the compressed signed value is interpolated bilinearly onto a 4× finer texture before colouring. Where neighbouring samples change sign or differ by more than 0.3 of the colour scale (sign changes, steep edges at large dynamic range), the field is instead evaluated exactly at every fine texel (adaptive refinement). This adds about 0.9 ms of the 7.3 ms per frame in "Synchrotron light" at t = 150. Most of that cost is the retarded-time search over the long spiral history.
  - **GPU magnetic map.** Where the colour changes by more than 5 % between neighbouring pixels (sign changes next to wires and magnets), it is supersampled 4 × 4 per pixel.
  - **Contours.** On the GPU maps, contour lines fade out smoothly where they would crowd closer than a few pixels, instead of being cut off (which aliased into moiré).
- **Particle.** The particle-field and total views show one particle, the selected shot and disturbance, and only its flight is drawn.
- **Scales.** The colour is B_z, which is the whole of B in the plane, on an asinh (log-like) scale. It is saturated at the 99th percentile of |B_z| sampled over one RF period or over the flight.
- **Wavelength.** The wavelength of the particle's cyclotron radiation, 2πc/ω_c, is often much larger than the arena (≈ 120 cells in "Synchrotron light"). The arena then lies in the near and induction zones, and the map shows the rotating near field rather than detached spiral wave fronts. That is the physically correct picture.
