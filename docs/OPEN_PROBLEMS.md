# Open problems

Questions met in this project that the literature does not answer (as far as we searched)
and that we could not derive in reasonable time with the tools at hand. Each entry records
the question, the evidence, what was searched and tried, and the workaround the project uses
meanwhile. Resolved entries stay, marked resolved. Maintained with the `open-problems` skill.

## OP-1: Eddy currents in an elongated strip of the scaled Drude metal

- **Status:** resolved 2026-10-19: not the eddy currents, and not a bug of the particles. A
  grounded electrode in an open z-invariant world (no enclosing conductor) is held at 0 in
  the free-space gauge `φ = −2λ ln r`, whose zero lies at the unit radius; it can draw net
  charge from that "infinity", and with net charge the potential coefficients `−2 ln d`
  (negative beyond d = 1) leave the electrostatic energy unbounded below: a large floating
  conductor next to it polarizes along the line to it while the electrode charges up, without
  end. The strip 8 × 32 with a small grounded electrode 1.5 below it ran away with no field
  at all (B = 0, no wave: the electrode's charge 0.7 → 228 over t = 5–60, doubling every ~5;
  the lower half's carriers −128 → −214), as in a steady B = 1 (0.6 → 169); with the
  electrode floating (`Bias::Charge(0)`) it stayed at 0 ± 0.007, and with a grounded box
  around everything (a reference at a finite distance, all coefficients `−2 ln(d/R)`
  positive) at 0 ± 0.04 (temporary tests in drude.rs, `static_field_explore`, removed). The
  eddy physics itself agrees with the continuum and the analytic law when the reference is
  sound (D8, the square and the strip 8 × 32 with a floating electrode). The shipped levels'
  metals all sit between contacts held at potentials (none floats); a level with an isolated
  metal needs its electrodes floating or an enclosing grounded conductor.
- **Arose in:** the eddy currents and laminations (task of 2026-10-19): `physics::tube::
  SolenoidWave` (a uniform B_z(t) and its induced `E = −½ Ḃ ẑ × (x − x₀)`), validation D8 in
  `crates/physics/tests/drude.rs`; the explorations that gave the numbers below were
  temporary tests in that file (`eddy_param_explore`, `eddy_profile_explore`,
  `strip_profile_explore`, `strip_long_explore`; their code is in the session's scratch
  scripts `target/scratch/add_*_explore.py`).
- **Question:** an isolated rectangle of the tests' metal (macroparticles of `q/m = −1` over
  a background 1, thermal speed 1, random hard-disc cores with `ℓ = 0.5`, specular walls,
  z-invariant electrostatics with the particles as rods of radius 0.05; σ = 0.235 from the
  Boltzmann equation, D1; dielectric relaxation `1/(4πσ)` = 0.34) in a field ramping in a
  triangle between `±b0`. Ohm's law predicts the steady eddy current `ψ = −σ Ḃ u` (u the
  rectangle's torsion function) within each ramp. What sets the time the strip takes to
  reach it, and what is the large, Ḃ-signed, B-dependent current along an elongated strip
  that replaces it? An answer: the strip's charge and current modes (relaxation times,
  frequencies) as functions of the aspect ratio, `ω_c τ` and `4πσ`, from a continuum model
  (Ohm's law with the Hall tensor, charge conservation, 2D Poisson) that reproduces the
  particles' numbers below, or a reason why the particles must differ from it.
- **Why it matters:** laminations are elongated strips. Until this is understood, a level
  whose heat depends on strips' eddy currents rests on numbers that are not Ohm's law's.
- **Known:**
  - A square 16 × 16 at `Ḃ = 0.05`, `b0 = 0.25`, period 20 (2 time units skipped after each
    corner): power −1.74 ± 0.67 %, current through the left half +1.16 ± 1.72 % of the
    torsion series (D8); it settles within half a time unit of each corner. At `b0 = 1`,
    period 80: about right throughout (power ≈ 1.3 against 1.35, current 2.0–3.0 against
    2.43). At `Ḃ = 0.1` it comes out 5 % low at any `b0` (0.2 or 0.5): the carriers warm
    (`qEℓ/kT ≈ 0.2` at its edges), a known nonlinearity.
  - A strip 8 × 32 centred on the axis, `Ḃ = 0.05`, `b0 = 0.25`, period 20: power +5.1 ± 1.5 %,
    current through the left half −37 ± 2.5 %; the current's profile across the middle
    third linear but shifted (a net current along the strip, reversing with Ḃ), and it
    changed when the one grounded electrode the electrostatics needs moved from the
    strip's corner to below its middle (then symmetric, but 40–60 % of the theory at the
    edges): not a steady state, since an outside electrostatic field cannot change a steady
    eddy pattern. With period 80 and `b0 = 1` (same Ḃ): the left half's current swings over
    each ramp from −3.6 to −9.5 and back to +3.7 (theory +2.52), following B, while the
    power is only ~0.25 (theory 0.68): mostly not an Ohmic current.
  - Four touching strips 4 × 16 (the square cut in four), `Ḃ = 0.1`: each outer strip
    carries, at each corner, the uniform current of the induced field's uniform part
    across it (~5.4, σE_y A with E_y ≈ 0.3), decaying with τ ≈ 6 (expected
    `1/(4πσN)` ≈ 1.7 with N ≈ 0.2 for a single strip): no carrier crosses between the
    strips (0 of 20 000 in 10 time units).
- **Narrowed (later on 2026-10-19), continuum model** (`target/scratch/eddy_continuum.py`,
  dense, h 0.5/0.25; `eddy_continuum_fft.py`, FFT, h 0.1 resolving the Debye layers):
  Ohm's law with the model's Hall tensor (drude.rs `mobility`), diffusion, charge
  conservation and the z-invariant kernel reproduces the square's series to 4 digits and
  takes the strip 8 × 32 to the stream function within ~3 time units, with or without a
  grounded electrode (centred below: no effect; at the corner: a slow rise, the strip's
  longitudinal polarization, τ of several units), at `b0` 0.25 or 1; its electrode's charge
  stays ≤ 0.04. The particles agree with it when the electrode is far (40 away: the strip's
  left-half current ≈ 2.5, power ≈ 0.67, also at `b0 = 1`). With the electrode 1.5 below the
  strip at `b0 = 1` they do not: the electrode's charge (the solver's, consistent with the
  strip's potential at it, −0.4 × φ) grows 0.18 → 0.93 by t = 10 and to 150–220 later; the
  lower half's carriers go from −128 to −167…−201 (a quarter of all the carriers moving
  toward the electrode); the anomalous current along y flows in both halves' interiors with
  the same sign (−12 each), not in the wall layers (≤ 1): a runaway longitudinal
  polarization that the continuum (inertia-free, local) does not have. Background and
  carriers' potentials have consistent signs (they cancel to −0.45 beside the strip at
  t = 0). Candidates: the carriers' inertia (a magnetoplasmon-like instability with the
  electrode's image charge as feedback), something kinetic near the walls, or a numerical
  effect of the particle push with a nearby electrode; not yet tested.
- **Literature checked:** web searches for edge magnetoplasmons in conducting strips with
  Drude dissipation, and for eddy-current transients in thin strips with the Hall effect.
  Closest: edge magnetoplasmons in 2D sheets (classical, Drude: e.g. arXiv:1905.05105 on
  gate screening and dissipation, arXiv:2112.00892, arXiv:1812.08798): modes along the
  edges of a 2D conductor in B whose velocity and damping depend on σ_xy, σ_xx and screening
  by nearby conductors (as the electrode's move did here), but for a 2D sheet, not a
  z-invariant slab, and without a drive by Ḃ; eddy-current transients in slabs and
  laminations (Knoepfel's book, via a CERN accelerator school lecture): the skin and step
  responses of a conductor with an instant dielectric relaxation, not a scaled one; a 1997
  phys-l discussion of the surface charge that steers induced currents into a steady loop,
  the Hall effect named but not worked out.
- **Attempts:** in-session (about two hours): the walls (specular: no size effect), leakage
  between strips (none), the field's rate (nonlinearity found and removed for the square),
  B's size (no effect on the square), the electrode's place (changes the strip), ramps four
  times longer (the strip departs further). No continuum calculation yet.
- **Workaround in use:** D8 validates the induction on the square only; no level uses eddy
  currents in elongated pieces. Levels with eddy currents wait for this entry (or use
  compact pieces, validated, and slow ramps).
- **Next steps:** solve the continuum model (`j = σ̂(B(t)) E`, `∂ρ/∂t = −∇·j`, E from the
  charges in the z-invariant geometry plus the induced field) on the strip by finite
  differences and compare its current with the particles'; if they agree, the effect is the
  scaled metal's (a slow, Hall-coupled surface-charge mode, an edge-magnetoplasmon analogue)
  and the levels need a denser metal or slower fields; if not, look for the particles'
  departure (the softening against the Debye length, 0.05 against 0.28; the time step).
