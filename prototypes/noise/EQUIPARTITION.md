# A floating contact's noise against equipartition (2026-11-15)

The thermal-noise level (docs/DRUDE.md Phase 2 item 1) rests on equipartition: a floating
contact's potential spreads as `⟨V²⟩ = kT/C` whatever feeds it, so the resistance sets the
noise's pace and not its size. The plan's first physics test checks it. Measured here, before
writing the test: **the model's contacts do not give it.**

## The setup and what is exact

`equipartition.rs` (copy into `crates/physics/tests/` to run; the `RESERVOIR` modes need
`reservoir_contacts.patch` applied to `crates/physics/src/tube/emission.rs`). A strip `L × 4`
of the tests' metal (D1's: background 1, `q/m = −1`, thermal speed 1, `ℓ = 0.5`), its left
contact grounded (an ideal source), its right contact a floating circuit node with a lumped
capacitance `C_L` to ground, a grounded frame 3 units round. Both contacts ohmic, as in D7.
Sampled every step over t = 10 to the run's end; 16 or 32 seeds.

The reference is exact in the model, not an approximation: the node's potential is
`V = (Q − ind)/C_node` (`circuits.rs`), with `C_node` the node's own coefficient of the
boundary elements' Maxwell matrix plus `C_L`, so its energy is `(Q − ind)²/2C_node`. In any
Gibbs state in which `Q` exchanges with a reservoir, `Q − ind` is Gaussian with variance
`kT C_node` at every configuration of the carriers: `⟨V²⟩ C_node/kT = 1`, with kT a
macroparticle's (`|q_c| σ_v²/|q/m|`, D7's), whatever the strip, its resistance or its own
capacitance. (Equivalently, by fluctuation–dissipation: the response to a charge put on the
node is `V(0⁺) = δQ/C_node`, and `V(∞) = 0` once the strip has drained it.)

## Measured: `⟨V²⟩ C_node/kT` (statistics: the seeds' standard error)

| Contacts | L | C_L (C_node) | carriers per area, dt | result |
|---|---|---|---|---|
| the game's (ohmic: the layer refilled to its density) | 4 | 0 (0.410) | 31, 0.02 | **1.178 ± 0.025** (16 × 500) |
| the game's | 4 | 0 | 31, 0.01 | 1.178 ± 0.045 |
| the game's | 4 | 0 | 125, 0.02 | 1.179 ± 0.04 |
| the game's | 2 | 0 (0.483) | 31, 0.02 | **1.427 ± 0.03** |
| the game's | 2 | 2 (2.483) | 31, 0.02 | **1.386 ± 0.034** (32 × 500) |
| regular reservoir (the half-Maxwellian flux of density 1, emitted on a schedule) | 2 | 0 | 31, 0.02 | 0.884 ± 0.022 |
| regular reservoir | 4 | 0 | 31, 0.02 | 0.916 ± 0.032 |
| regular reservoir | 2 | 0 | segments 0.125 | 0.921 ± 0.026 |
| regular reservoir | 2 | 0 | softening 0.025 | 0.860 ± 0.023 |
| regular reservoir | 2 | 0 | dt 0.01 / 0.04 | 0.882 ± 0.020 / 0.844 ± 0.019 |
| regular reservoir | 2 | 0 | 125, 0.02 | 0.857 ± 0.023 |
| regular reservoir | 2 | 2 | 31, 0.02 | 0.775 ± 0.029 |
| Poisson reservoir (the same flux at random times) | 2 | 0 | 31, 0.02 | 1.211 ± 0.033 |
| Poisson reservoir | 2 | 2 | 31, 0.02 | 1.147 ± 0.036 |

And Nyquist on the short strip between two fixed contacts (D7's measure, the charge carried's
spread slope over T = 1–4 against `2kT/R`; L = 2 at 0.2, 16 × 400): the game's contacts
1.011 ± 0.050 (R 2.23), Poisson reservoirs 1.059 ± 0.055 (R 3.64), regular ones 0.699 ± 0.037
(R 3.50).

## What this says

- With the game's contacts the floating contact's noise is **18 % (L = 4) to 43 % (L = 2)
  above kT/C**, unchanged by the step, the macroparticles' weight and the lumped capacitance:
  a property of the continuum model, resolution-free, and growing as the strip shortens, so
  **the size depends on the resistance feeding the contact**, the opposite of the level's
  idea. Its spread *ratio* between two contacts stays resolution-free (the excess cancels only
  where both strips are alike).
- The same contacts pass Nyquist between fixed potentials even on the short strip (as D7 at
  L = 16): the excess appears only when the contact floats. So D7 does not establish that the
  contacts are in equilibrium with a floating node (PHYSICS.md §2.11's *Thermal noise*
  paragraph says it does).
- The excess comes from the contacts' emission rule: replacing it changes the answer by tens
  of per cent either way. The regular train is sub-Poissonian (its Nyquist 0.70, its floating
  noise below kT/C, falling further as C grows); Poisson emission passes Nyquist and still
  floats 15–21 % high at the game's release distance.

## The cause (found 2026-11-15, the owner's choice "research the cause first")

Two departures from detailed balance, each measured alone (strip 2 × 4, `C_L` = 0, 16 × 250):

| Release distance `δ` (`EMISSION_OFFSET`, segments; segments 0.25) | the game's contacts | Poisson reservoir |
|---|---|---|
| 0.5 (the game's) | 1.427 ± 0.03 | 1.211 ± 0.033 |
| 0.25 | 1.142 ± 0.024 | 1.080 ± 0.025 |
| 0.125 | 1.152 ± 0.036 | 1.035 ± 0.032 |
| 0.5, without the return absorption (`RETURN_LEAVE` 1e9) | 1.484 ± 0.031 | |

1. **The release at δ.** A contact gives a carrier off δ out from its face (half a segment,
   0.12–0.2 here: about the Debye length) but takes carriers in at the face. The potential at
   δ differs from the electrode's by the face's field over δ, so the contact's electrochemical
   potential follows the field at its face, which a floating contact's own charge sets: a
   feedback that heats the contact's mode. The excess falls linearly with δ for the Poisson
   reservoir, which tends to 1. The fix: weight each release by the Boltzmann factor of its
   potential step, `exp(−βq(φ(δ) − V))`, which puts the reservoir's electrochemical potential
   back at the electrode's. Measured with the game's δ: 0.982 ± 0.026 (L = 2) and
   1.018 ± 0.035 (L = 4) with the potential evaluated at each emitting segment's release
   point; 0.987 ± 0.018 and 0.990 ± 0.021 with the step to first order,
   `φ(δ) − V ≈ −4πσ δ` (σ the segment's surface density, already computed each step: no
   cost). Nyquist between fixed contacts with it: 1.062 ± 0.061.
2. **The refill clamp.** Filling the layer to its density whenever it lacks is no
   detailed-balance reservoir: it keeps ~15 % at small δ, where the Poisson wall reaches 1.
   (At finer segments the clamp also stops emitting from segments whose layer holds less
   than one macroparticle, `density × length × depth < weight`: at segments 0.0625 the
   contact never emits and the floating noise collapses to 0.12. The game's levels are far
   from it.)

The return absorption (a contact's own carrier back within 1.5 δ) is not a cause.

**The price of a detailed-balance wall:** its contact resistance, Sharvin's
(`kT/(q j₀ w)`, `j₀ = ρσ_v/√(2π)` the charge flux a half-space carries: 0.63 per contact of
width 4 here; the short strip's R rises from 2.23 to 3.54). The clamp has none, which is what
the metal levels' ohmic contacts were designed for. With the Hall probe's metal (density 2,
thermal speed 2, ℓ = 1) a 16 × 4 strip's 8.5 would gain about 1.3 (15 %). A grand-canonical
layer (carriers inserted and removed in the layer with the Boltzmann factor of each move)
would be a detailed-balance contact without that resistance, at the cost of a potential per
trial; not tried.

- The prototype's pace 1.46 times its decay (docs/HANDOFF.md) was with the clamp.

## In the engine (DEC-43, the owner's choice: a contact kind)

`MetalSpec::thermal_contacts` (`Tube::emit_thermal`, the first-order factor). Tests in
`crates/physics/tests/drude.rs`: D13 `⟨V²⟩ C/kT` = 1.004 ± 0.014 (L = 2), 0.986 ± 0.024
(L = 4), 0.961 ± 0.021 (L = 2 with `C_L` = 2, 32 seeds over 1000; the variance about each
run's own mean reads ~1 % low); D14 a kick's decay 1.092 ± 0.013 against the noise's pace
1.104 ± 0.017 (ratio 0.990 ± 0.019; the curves within 0.016 over three time constants);
Y31 charge conjugation bit for bit; D10 (c) the energy account to 2e-15.
