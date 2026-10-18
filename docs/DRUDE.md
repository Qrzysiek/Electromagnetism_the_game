# Inside the components, step 2: metals and resistors (Drude) — plan of the stage

`SPEC.md` ("Inside the components", item 2): classical electrons scattering off the lattice
with a mean free time; Ohm's law and Joule heating emerge; stated as a model (real
conduction electrons are a quantum Fermi gas); random collisions with a fixed seed. This
note is the plan, written before building anything (as `docs/TUBES.md` was). Questions for
the owner are marked **(owner)**; nothing is built until they are settled.

## What the engine has and lacks

- **Has (from the tubes, PHYSICS.md §2.11):** the z-invariant world (electrodes as long
  prisms, carriers as lines of charge, exact 2D electrostatics by boundary elements, the
  carriers' induced charges and Ramo currents), thousands of macroparticles per step
  (the tree code), circuits with sources, resistors, capacitors and coils, the uniform
  field B_z of a solenoid, current and potential goals, the potential, electric, magnetic
  and energy-flow maps (the last two from the run's Maxwell field), smooth playback.
- **Lacks:**
  1. *A conductor's interior:* a region the carriers live in (its walls reflect them),
     with a neutralizing positive background (the ions).
  2. *Scattering:* collisions with the lattice that randomize the carriers' motion and take
     the energy the field gives them (Joule heat), with a seeded random generator.
  3. *Ohmic contacts:* electrodes that take carriers in and give them off so that the
     region next to them stays neutral (the standard boundary condition of device
     Monte Carlo codes), at the lattice's temperature.
  4. *Thermal carriers:* a Maxwellian start and a temperature.

## Model

- **The world:** z-invariant again (the owner's rule: only 3D physics a symmetry
  justifies): a resistor long along z, the current flowing in the plane. Carriers are
  macroparticles (lines of many electrons) as in the tubes.
- **The lattice, two options (owner, question 2):**
  - *(a) Explicit ions:* fixed discs (the ion cores) on a lattice or scattered at random;
    an electron that hits one bounces off it and leaves with a speed drawn from the
    lattice's temperature (a thermalizing collision: the lattice takes the energy the
    field gave it, the Joule heat). The mean free time is not a parameter: it emerges
    from the geometry, `ℓ = 1/(2 a n_ion)` in 2D for discs of radius a at random, `τ = ℓ/v`.
    Everything is visible: the zigzag paths, the drift superposed on the thermal motion,
    the collisions where the heat goes. Without the thermalizing rule the model is the
    Lorentz gas, deterministic chaos with normal diffusion (Bunimovich–Sinai), but the
    field would heat the electrons without bound.
  - *(b) Relaxation:* no ions drawn; each carrier, at a rate 1/τ (a Poisson process),
    takes a new velocity from the Maxwellian (Drude–Lorentz). τ is a parameter; the
    textbook formulas hold exactly for this model; the mechanism is not seen.
- **The background:** the ions' charge as a uniform positive density over the conductor
  (jellium), in the field sums as fixed charges, so the conductor is neutral at rest;
  screening emerges (the Debye length `λ_D = √(kT/4π n q²)`).
- **Scales:** a real metal (n ~ 1e29 m⁻³, τ ~ 1e-14 s, a degenerate Fermi gas) cannot be
  simulated particle by particle; the game's metal is a classical plasma with the
  orderings that make it behave as one: a screening length and a mean free path well
  below the resistor's size (fields inside a conductor carrying current are uniform and
  its charges sit at its surfaces; transport is diffusive), the plasma period resolved.
  Stated in the model notes, with what is classical (the thermal speed, the electronic
  specific heat, the Lorenz number) and what carries over (σ = n q² τ/m, the Hall
  coefficient 1/(nq), the relaxation times).
- **Heat (owner, question 4):** a lattice at a fixed temperature (a bath) first; the Joule
  heat it takes is measured. A lattice that warms (its temperature rising where the heat
  goes, τ changing with it) later, if wanted.

## Validation (planned tests, references fixed before measuring)

- **D1 conductivity:** a uniform field in a periodic box: the drift `v_d = qEτ/m` and
  `σ = n q² τ/m` (exact for (b); for (a) with τ measured from the collision count, and
  against `ℓ = 1/(2 a n_ion)` at low density).
- **D2 AC conductivity:** `σ(ω) = σ₀/(1 − iωτ)` (exact for (b)).
- **D3 Joule heating:** the power the lattice takes equals `J·E` per area, and the sources'
  work equals it in a circuit (the energy balance, as V16).
- **D4 Hall effect:** a strip in B_z carrying a current: the Hall field `E_H = J B/(n q)`
  (`R_H = 1/(nq)`, independent of τ), the sign telling the carriers' sign (Jackson Pr.
  6.15).
- **D5 relaxation:** an excess charge inside the conductor: `ρ̈ + ρ̇/τ + ω_p² ρ = 0`, the
  plasma oscillation damped, or the decay with `τ_d = 1/(4πσ)` when `ω_p τ ≪ 1`.
- **D6 Einstein relation:** carriers injected at a point spread with `D = μ kT/q`
  (`μ = qτ/m`).
- **D7 Johnson–Nyquist:** the open terminals' voltage fluctuates with `⟨V²⟩ = kT/C`
  (equipartition; the classical noise `4kTR` per unit bandwidth).

## Levels (sketch)

- *Introduction:* Ohm's law (a resistor between two contacts, the supply's slider: the
  current grows in proportion; the drift is slow against the thermal motion); the
  geometry sets the resistance (long and thin against short and wide, `R = L/(σ w)` per
  unit length along z).
- *Intermediate:* the Hall probe (a strip with side contacts in the player's solenoid: the
  Hall voltage tells the field, and the carriers' sign); the relaxation of a charge put
  inside; the thermal noise of a resistor (statistical goals).
- *Master:* later, with the tubes' and circuits' pieces (e.g. a Hall magnetometer reading
  a field the player must find, or a resistor network inside the tube world).
- Goals: the tubes' current and potential goals cover most; the noise needs an RMS measure
  (a small extension).

## Questions for the owner

1. **(owner)** The z-invariant world again (a resistor long along z)?
2. **(owner)** The lattice: explicit ions the electrons bounce off, the Joule heat taken at
   each collision and τ emerging from the geometry (proposed: the mechanism is seen), or
   Drude's relaxation at a rate 1/τ (exact textbook formulas, nothing to see)?
3. **(owner)** A scaled classical metal (screening length and mean free path below the
   size, stated as a model) — acceptable?
4. **(owner)** A lattice at a fixed temperature first (the heat measured), a warming
   lattice later?
5. **(owner)** The first levels: Ohm's law and geometry (introduction), the Hall probe
   (intermediate)?

## Decisions (owner, 2026-10-13)

1. **Explicit ions** (option a): fixed ion cores the electrons bounce off; each collision
   hands the lattice the energy the field gave (the Joule heat); the mean free time
   emerges from the geometry. Details fixed with the model: the cores' positions jittered
   off a lattice (a regular lattice has open corridors, infinite free paths, no Ohm's law);
   the bounce is specular off the disc, the new speed drawn from the lattice's temperature
   with the collision-weighted law `∝ v² e^{−mv²/2kT}` (in 2D, the 3D Maxwell speed law),
   so that a Maxwellian gas stays Maxwellian at rest; the ions' charge a uniform positive
   background over the conductor (jellium, its field in closed form for rectangles).
2. **Accepted:** the z-invariant world; a scaled classical metal, stated as a model; the
   lattice at a fixed temperature first; the first levels Ohm's law and geometry, then the
   Hall probe.
3. **References for the tests** from the Boltzmann equation of this model (exact for dilute
   cores, the Boltzmann–Grad limit), derived before measuring (Wolfram where needed).

## Progress

- **2026-10-13, the physics (PHYSICS.md §2.11, "Conductors: the Drude metal"):** the
  conductor's interior (`physics::lattice`: the jellium background in closed form, the
  cores, the collision step `Lattice::advance`), integrated in `Tube` (`lattices`, the
  Joule heat and the collision count in `TubeState`), validated against the Boltzmann
  equation of the model: D1 conductivity −1.8 %, collision rate +0.3 %, D3 energy
  balance to 2e-12, D4 the conductivity tensor in B_z 1.2 % (Hall factor 4/π), D6
  Einstein −0.4 %, and the self-consistent resistor between two contacts +0.1 %.
- **Changed from decision 1:** the cores are placed at random, not jittered off a
  lattice. A lattice jittered by less than a cell (to keep the cores apart) leaves open
  corridors, and even one jittered over whole cells has free paths more regular than the
  Boltzmann equation's: σ 10–22 % off. Random cores (an amorphous metal) have exactly
  exponential free paths; at `a² n_ion ≤ 0.005` the medium's own correction is ≤ 3 %.
- **Levels' constraint:** `a² n_ion = a/(2ℓ) ≤ 0.005`: the cores are small against
  their spacing (`a/s ≤ 0.07`); they are drawn as dots, larger than their radius.
- **2026-10-14, the levels' machinery:** metals in tube levels (`MetalSpec`: a rectangle,
  its contacts, the carriers' density, q/m and thermal speed, the mean free path and core
  radius, macroparticles per area), setup checks (dense cores, the plasma period and the
  mean free time resolved, nothing inside a metal), the editor's controls, a contact's
  lead current (`Carriers::Net`), potential goals against another electrode (a Hall
  voltmeter), the display (the outline, the cores, the carriers faint with a traced
  sample showing the zigzag and the drift, the Joule heat), the maps' and the Maxwell
  field's backgrounds, frames kept sparsely for long runs, the inspection tooltip.
- **Contacts, changed three times:** the cathodes' space-charge limit made a contact a
  cold reservoir (the Joule heat negative), then, releasing thermal carriers, a layer
  with its own voltage in front of an injecting contact (a potentiometer read 8 for 10).
  Ohmic contacts (the device Monte Carlo codes' condition: the layer in front of a
  contact kept neutral) give the bulk field V/L to 0.2 %; first they placed Maxwellian
  carriers in the layer (kT each, where those reaching the face bring 1.5 kT: every
  contact a heat sink, found when a wide strip's lattice heat came out negative), now
  they emit from the face as a thermal wall: the lattice takes 90 % of IV.
- **The levels' material:** kT/q = 4 (thermal speed 2), small against the supplies' tens
  so that the contacts' and walls' own voltages stay a few per cent; ℓ = 1, a = 0.01.
- **Levels:** Ohm's law and Potentiometer (introduction), Hall probe (intermediate), in
  the arc "Inside the components: metals".
- **The walls' layer is physics:** a metal at a potential against its surroundings
  carries a surface charge, spread here over the Debye length (0.4): depleted at the
  walls when positive, accumulated when negative, uniform at 0 (PHYSICS.md §2.11).
- **2026-10-14, D2 and D5:** the AC conductivity against the Boltzmann equation at
  `ωτ ≈ 1` (2.3 %); screening within the Debye length (+1.0 %, once the carriers'
  softening is well below it: at the levels' preview resolution the layers are ~8 % too
  thick). D5 as planned (an excess charge relaxing) has no clean reference here: with
  `τ ∝ 1/v` the slow carriers' long free times give the decay no single rate, and a
  finite strip's depolarization field is not uniform.
- **Levels:** Long and thin (introduction, after Ohm's law): two strips in series, the
  junction on a second supply; equal currents need V₂/V₁ ≈ 1/5.3 (L/W alone gives 1/5;
  the thin strip's stronger field warms its carriers). The Potentiometer and the Hall
  probe retuned for the thermal-wall contacts (the null at 9.05 for 30; B on 1.0–1.4),
  their goals averaged over t = 10–80 against the floating contacts' noise.
- **Open:** D7 (Johnson noise); a warming lattice; a Hall level with both signs of
  carriers.
- **2026-10-18, D7:** fluctuation and dissipation as Nyquist's: a shorted resistor's
  charge carried spreads as `2kT/R` (−1.4 %, statistics 1.7 %, 16 seeds; PHYSICS.md
  §2.11), with kT a macroparticle's. **No noise level**, as the levels' sketch had
  planned: a macroparticle is a line of charge infinite along z moving as one, so the
  noise is the model's graininess (it halves when the macroparticles double); a real
  conductor's carriers fluctuate independently along z and its noise vanishes in a
  z-invariant world. A noise goal would be neither verifiable (preview and verification
  disagree by design) nor real physics.
