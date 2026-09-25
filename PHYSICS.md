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

## 3. Equation of motion — *validated* (`crates/physics/src/dynamics.rs`)

State `y = (x, p)`, with:

```
γ     = sqrt(1 + |p|² / (m² c²))
dx/dt = v = p / (γ m)
dp/dt = q (E(x) + v × B(x))
```

`|v| = |p| c / sqrt(m²c² + |p|²) < c` holds for any finite `p`. The speed limit is structural, not enforced after the fact.

Stage 1 has `B = 0`.

Implementation details:
- `c = ∞` is allowed and gives exact Newtonian mechanics (`1/(mc)² = 0`, so `γ = 1`). The non-relativistic tests use this rather than a large finite `c`.
- The integrator works with the scaled state `(x, p/p_ref)`, where `p_ref = |p₀|` (or `m` if the particle starts at rest). One tolerance is then meaningful for positions (grid units) and momenta alike.
- Kinetic energy is evaluated as `(γ−1)mc² = p²/(m(γ+1))`, which avoids cancellation at low speed.
- **Speed limit in floating point:** `|v| = |p|c/sqrt(m²c² + p²)` never exceeds `c`. For `γ ≳ 10⁷`, the rounded quotient can equal `c` exactly. The strict inequality `|v| < c` is guaranteed and tested for `γ < 10⁷`. Nothing uses `v` as state, so this has no effect on the trajectory.

## 4. Conserved quantities (diagnostics only) — *validated* (`crates/physics/src/trajectory.rs`)

Static electric field:
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

All event functions are signed distances (sphere, box, world bounds), which are 1-Lipschitz in position. Along the trajectory, `|dg/dt| ≤ |v| ≤ v_max`. On an interval `[a, b]` with `g(a), g(b) > 0`:

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
| Radiation (Larmor/Liénard) | neglected | The radiated energy `∫P dt`, with `P = (2/3) k q² γ⁶ (a² − (v × a)²/c²) / c³`, is computed as a diagnostic (`Trajectory::radiated_energy`, trapezoidal rule over accepted steps; exactly 0 for `c = ∞`). Flagging levels above a threshold is part of the generator (M5). |
| | | **Consistency requirement (found in M4):** radiated fraction per close pass ≈ `r_cl / r`, with `r_cl = k q²/(mc²)` the particle's classical radius. The first demo levels used `q = m = 1`, `c = 1.5…5`, giving `r_cl` = 0.04–0.44 cells. Their neglected radiation was 1e-4 to 2 × T₀, so they were physically inconsistent. Fix, as in real accelerators: weakly charged particle (`q = 10⁻⁶`), strongly charged electrodes (`Q ~ 10⁶`). Trajectories depend only on `qQ/m` and `c`, so the puzzles are unchanged. Neglected radiation of the reference flights is now 1e-16 … 1.7e-11 × T₀. Test `level/tests/levels.rs` requires `< 1e-10` for every shipped level; the game flags player setups above that. |
| Magnetic field of the moving particle acting on others | n/a (single particle) | Stage 6 (Darwin) |
| Polarization of test particles | neglected (non-polarizable by assumption) | documented assumption |
| Recoil of fixed charges | none (held fixed by definition) | game rule |

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
