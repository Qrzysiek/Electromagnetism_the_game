# Physics: model, numerics, validation

This document is the authoritative description of the physics that is implemented and how. Each section has a status:
- **design:** planned, not yet in code
- **implemented:** in code, with the file and function named
- **validated:** implemented, with the tests that check it and their measured results

Any change to the physics or numerics must update this file in the same commit.

---

## 1. Units — *design*

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

In the non-relativistic limit, trajectories depend only on the dimensionless ratio `k q Q / (L E_kin)`. The physical scale enters through `c` (relativity) and, in later stages, through material properties and particle–particle interaction.

## 2. Field sources — *design*

Fixed charges `Q_i` at positions `r_i`, each a rigid sphere of radius `R_i` with a spherically symmetric charge distribution. Outside the sphere (the only region particles can reach), the shell theorem gives exactly:

```
φ(r) = Σ_i k Q_i / |r − r_i|
E(r) = Σ_i k Q_i (r − r_i) / |r − r_i|³
```

- The sum is evaluated in f64 in a fixed order (the order in the level file). No softening.
- Test particles are rigid, spherically symmetric and **non-polarizable** (a stated approximation, since a conducting microsphere would feel image forces). Under these assumptions, the force on a test particle of radius `a` equals `q E(center)` exactly, because the force between two non-overlapping spherically symmetric distributions equals the point-charge force.
- Contact condition: `|x − r_i| ≤ R_i + a`, which means the particle is lost.

## 3. Equation of motion — *design*

State `y = (x, p)`, with:

```
γ     = sqrt(1 + |p|² / (m² c²))
dx/dt = v = p / (γ m)
dp/dt = q (E(x) + v × B(x))
```

`|v| = |p| c / sqrt(m²c² + |p|²) < c` holds for any finite `p`. The speed limit is structural, not enforced after the fact.

Stage 1 has `B = 0`.

## 4. Conserved quantities (diagnostics only) — *design*

Static electric field:
- Total energy `W = γ m c² + q φ(x)` is conserved.
- A single central charge also conserves angular momentum `L = x × p`, which holds relativistically too.

They are computed along every trajectory and reported. They are **never** used to project or correct the state.

## 5. Integrator — *validated* (`crates/physics/src/integrator/dop853.rs`)

- Dormand–Prince 8(5,3) (DOP853, Hairer, Nørsett & Wanner, *Solving ODEs I*, §II.10), with a 7th-order dense-output interpolant.
- Mixed error control on each component: `err ≤ atol + rtol·|y|`.
- Two tolerance levels:
  - **preview** (fast, live editing)
  - **verify** (at least 100× tighter)
- The exact values are fixed in M2 after benchmarking the physics problem, and recorded here.

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
- The step size is also capped so that the particle's displacement per step is a fraction of its distance to the nearest charge surface (a safeguard; event location is the primary mechanism).
- Why not Boris: with `B = 0`, Boris is leapfrog. Its good long-term energy behaviour relies on a fixed step, and an adaptive step (unavoidable with close encounters) breaks that. It may be revisited for magnetic fields.

## 6. Events — *design*

Events are located on the dense-output polynomial of each accepted step:
- `g_i(t) = |x(t) − r_i| − (R_i + a)` crosses zero: collision with charge `i`.
- `x(t)` enters region B: success.
- `x(t)` leaves the world bounds: lost (a game rule).
- `t > t_max`: timeout (a game rule).

Within a step, each `g` is sampled at several points of the interpolant to catch grazing approaches that dip below zero and come back. Roots are refined with Brent's method. The earliest event ends the trajectory.

## 7. Outcome verification — *design*

For every outcome that matters (reference solutions, final judgement, and the background refinement of the preview):
1. Integrate at the verify tolerance.
2. Record the minimum margin to every event boundary that was **not** triggered, and the margin at the boundary that was.
3. The result is **verified** if the outcome matches the preview result and every margin exceeds a safety factor times the error estimate (from comparing against a run at the next tighter tolerance). Otherwise it is **marginal**.

This separates numerical uncertainty (which must never decide the outcome) from physical sensitivity (chaos and grazing orbits, which are real and are handled by level-design robustness requirements).

## 8. Neglected effects and their control — *design*

| Effect | Status in Stage 1 | Control |
|---|---|---|
| Radiation (Larmor/Liénard) | neglected | The radiated energy `∫P dt`, with `P = (2/3) k q² γ⁶ (a² − (v × a)²/c²) / c³`, is computed as a diagnostic. If it exceeds a set fraction of the kinetic energy, the level is flagged. |
| Magnetic field of the moving particle acting on others | n/a (single particle) | Stage 6 (Darwin) |
| Polarization of test particles | neglected (non-polarizable by assumption) | documented assumption |
| Recoil of fixed charges | none (held fixed by definition) | game rule |

## 9. Validation tests — *design*

Thresholds are provisional until the integrator is benchmarked. Measured results are filled in when each test is implemented.

| # | Test | Analytic reference | Pass criterion | Result |
|---|---|---|---|---|
| T1 | Energy conservation, many-charge field | `W` constant | relative drift < 1e-10 at verify tolerance | – |
| T2 | Rutherford, non-relativistic | `tan(θ/2) = k q Q / (m v∞² b)` | relative error < 1e-8 | – |
| T3 | Relativistic Coulomb scattering | see below | relative error < 1e-8 | – |
| T4 | Kepler orbit (`c → ∞`) | `T = 2π sqrt(m a³ / |kqQ|)`, orbit closes, `L` constant | period error < 1e-9, closure < 1e-8 | – |
| T5 | Relativistic orbit precession | `Δφ = 2π(1/Γ − 1)` per revolution | relative error < 1e-6 | – |
| T6 | Hyperbolic motion in uniform `E` (test solver) | `x(t) = (mc²/qE)(sqrt(1 + (qEt/mc)²) − 1)` | relative error < 1e-10 | – |
| T7 | Speed limit | `|v| < c` | always, for extreme fields and energies | – |
| T8 | Plane of symmetry | out-of-plane coordinate stays 0 | exactly 0 (bitwise) | – |
| T9 | Grazing events | known closest approach `d` versus `R ± δ` | correct classification down to `δ` at the tolerance level | – |
| T10 | Convergence | error ∝ tolerance | slope consistent with the method | – |
| T11 | Determinism | – | bit-identical across runs and Windows/Linux builds | – |

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
