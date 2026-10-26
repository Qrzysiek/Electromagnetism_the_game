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

## What the metal can still show (survey, 2026-10-20)

The owner asked (2026-10-20) what more the metal can do: AC, capacitors, other components,
RLC, receiving antennas, the response to general fields, and what it combines with. The
model as it stands: carriers bouncing off ion cores (Joule heat), ohmic and floating
contacts, the tube world's circuits (DC, sine and pulse sources, R, C, a tank's L, floating
nodes and wires between electrodes), the solenoid's field (steady or swinging, with its
induced E), cathodes and beams in the same board. Its metals are doped semiconductors by
their own ratios (ω_p τ 0.7–2, a screening length 0.4–1.1 mean free paths: `scales.rs`,
`resemblance`), which is what makes several of the items below possible at all. It lacks:
the carriers' own magnetic field (no inductance from geometry), two carrier species in one
metal, a lattice whose temperature varies, external fields other than the solenoid's, a
drive's frequency as a control, metal shapes other than rectangles with cuts.

**Now: levels only (the engine as it is)**

- *The field effect, the transistor's principle.* A gate electrode beside a thin strip
  pushes its carriers out (depletion) or draws more in (accumulation): the strip's current
  follows the gate's potential. A metal cannot be gated (it screens within an atom); the
  model's low density can. Prototype (`target/scratch/fet_proto*.py`; a strip 16 × 2 of
  density 0.5 between a source and a drain at 5, gates 1 above and below on one supply):
  the strip held 3158 carriers at gate 0, 1945, 1006, 655, 494 at −20, −40, −60, −80 (84 %
  depleted), 7761 at +20; the drain current over t = 10–30 still followed the channel's
  charging through its own resistance (a level needs it settled first). Levels: a gate
  that pinches the current off; an amplifier (a sine on the gate, a larger swing at the
  drain through a load strip); an inverter (a gated strip and a load strip in series).
- *A Wheatstone bridge of strips,* joined by the wires between floating contacts (Strength
  and sign's junction): four strips, a supply across, a voltmeter between the two
  midpoints (a potential goal with a reference). One or two arms magnetoresistive in the
  solenoid: the bridge balances at one |B|, a real magnetoresistive sensor's layout.
- *A charge's relaxation and the plasma's ringing* (the plan's "relaxation of a charge put
  inside"): a pulse on a contact or a floating plate; with ω_p τ ≈ 2 the metal rings for a
  couple of periods before it settles, the silicon-like one (0.7) creeps without ringing.
  Real: the dielectric relaxation time ε/σ (copper 1e-19 s, the reason a metal is neutral
  inside).
- *RC and RLC with a microscopic resistor:* a strip charging a floating contact (C); a
  tank (the circuits' L and C) damped by a strip, its Q set by the strip's resistance, and
  with a magnetoresistive strip by the solenoid: the field tunes the ringing (an amplitude
  goal).
- *The Hall effect on AC:* the swinging solenoid on a strip carrying a steady current, the
  Hall voltage's amplitude at ω (the amplitude goal; the maps' lock-in): lock-in detection.
- *Beams and metals together:* a cathode's beam landing on a metal (its electrons become
  carriers: the same species scatter alike) is a resistive anode, its heat where the beam
  lands (the heat map, task 11); a Hall probe reading the solenoid that must turn a beam
  onto a collector (Sideways with a magnetometer), a master across the tubes and metals.

**Small extensions: a field or a control, with its tests (days)**

- *An RMS goal* (the plan's "statistical goals"): thermal noise, Johnson–Nyquist's
  4kTRΔf (test D7 has it) and kT/C, a floating contact's voltage noise that does not depend
  on the resistance feeding it (task 20).
- *A drive's frequency on a slider:* resonance puzzles (tune a tank; find a floating
  plate's plasma resonance); the Drude roll-off at ωτ = 1, where the current lags the
  field by 45°: the carriers' inertia, a kinetic inductance m/(nq²) per square (test D2
  has the AC conductivity).
- *Contact temperatures* (a hot and a cold contact): thermoelectricity. Seebeck's voltage
  between hot and cold ends, Peltier's heat carried by the current (one junction cools),
  a thermocouple of Few or many's two densities. For metals Drude's classical estimate is
  famously a hundred times too large; for the nondegenerate carriers of a doped
  semiconductor it is right, and those are what thermoelectric generators are made of. A
  lattice that warms (the owner's decision 4, "later") adds a resistance rising with
  temperature and thermal runaway.
- *Metal shapes beyond rectangles* (a ring, an L or U): the Corbino disc (radial current,
  the Hall current circling: the largest geometric magnetoresistance), and a pickup loop
  whose ends feed a load: Faraday's law driving an external circuit, a generator with the
  solenoid as its primary. (A straight strip closed by an abstract wire will not do: the
  EMF of a loop needs the loop's area.)
- *A uniform AC electric field* (a long wave, quasi-static) on the carriers and the
  electrodes: a metal's response to a general field. An isolated plate polarizes (its
  charge piles on its faces within 1/ω_p); a closed metal box shields its inside below
  ω_p and leaks above (a Faraday cage); a strip with contacts on a tuned tank is a
  receiving antenna, and with a vacuum diode and an RC after it a crystal radio (antenna,
  tuned circuit, rectifier, smoothing: AM demodulation), a master across the metals,
  circuits and tubes.

**New stages: new physics with its validation (weeks)**

- *The carriers' own magnetic field* (magnetoquasistatics, the Darwin limit; particle
  codes since Nielson and Lewis, 1976). In the z-invariant world a current in the plane
  makes a field along z (a metal ring is a solenoid: uniform field inside, none outside),
  and its change an induced E. It brings inductance from geometry, transformers (nested
  rings sharing flux), the skin effect (a slab's current crowding to its faces at
  δ ∝ 1/√(σω)), AC shielding by a ring, induction heating at surfaces, and the eddy
  currents' own field (Eddy currents is valid for thin plates only). It needs a finite c;
  the references are analytic (a ring's inductance, a slab's skin depth, nested rings'
  mutual inductance).
- *Semiconductor devices:* electrons and holes in one lattice (Positive carriers has
  holes), a doping profile (the background's sign changing across a junction), generation
  and recombination (np = n_i² at rest). The p–n junction (its built-in potential, its
  depletion layer, the diode's exponential law), the LED and the photodiode (recombination
  as light, generation by light), the bipolar transistor. Device simulation's ensemble
  Monte Carlo (Jacoboni and Lugli, 1989) is this model plus band structure.
- *Full electrodynamics for the carriers* (the tube world's Maxwell field, today computed
  after the run for the views, acting on them: a 2D electromagnetic particle-in-cell):
  antennas with retardation (the half-wave resonance), a metal that reflects below ω_p and
  turns transparent above it (why metals shine and alkali metals pass ultraviolet), lossy
  waveguide walls, surface plasmons.

**Combining with what exists.** The metals live in the tube world, with its circuits,
cathodes, beams, the solenoid and the Maxwell field; they do not reach the particle levels'
3D slice. The pieces that join them: the wires between floating contacts (series and
parallel networks of strips, bridges), the tanks and sources of the tube circuits (RC, RLC,
receivers), the cathodes (beams into metals; the vacuum diode as a radio's detector), the
swinging solenoid (AC Hall, generators, transformers once the currents make their own
field).

**Order:** see the implementation plan below.

## Implementation plan (the owner approved every idea, 2026-10-20)

**Order.** Levels the engine can do now come first: they need no new physics, and they
exercise the pieces the later steps reuse. Then the small extensions in the order of what
they unlock: the RMS goal stands alone; the frequency control serves the resonances and the
receiver; ring shapes serve the Corbino disc, the generator, the Faraday cage and the
magnetic stage; the uniform AC field needs both. Temperatures stand alone. Then the three
stages: the carriers' own magnetic field first (the AC world the owner's list began with:
inductors, transformers, the skin effect), the semiconductor devices second, full
electrodynamics last (the largest, and it overlaps the first). The owner (2026-10-20):
"first we need the separate elements, then we can think of joining them": the levels that
join worlds (metals with the tubes' beams and diodes) wait for Phase 4, toward the goal
that at each scale every element cooperates with every other (the phenomena and
materials, not necessarily today's models), showcasing everything that can be shown.

**Every step, done** (CLAUDE.md): new physics with its PHYSICS.md section, validation tests
against references fixed before measuring, and the symmetry and conservation tests; new
level fields with editor controls and `check_editable`; new marks with legend rows; new
element kinds through `scripts/pairs.py`. Every level: designed by scans, its windows at
least 5 % (clear of 0, the ends and the halves' middles), golden hash, `docs/LEVELS.md`,
its scales (`scales.rs`), a capture; the full suite; fast-forwarded to main.

**Arcs.** The metals arc keeps the sensors (the bridge, the Corbino disc). Three new arcs,
each with its introductions, intermediates and a finale (docs/CURRICULUM.md): *Inside the
components: transistors* (the field effect, later the junctions), *AC and induction* (the
plasma's ringing, RLC with strips, lock-in, resonance, the generator, the receiver, later
inductors and transformers), *Heat and noise* (thermal noise, thermoelectricity, a warming
lattice). Finales across arcs join them with the tubes.

**Step 0, tooling: the run cache (the owner: "a must"; a day).** A cache of the
generator's run results, hash-locked: each result is stored under the hash of everything
that could change it, so a stale one can never be served:
- the level, canonically serialized (every field, not the file's name or number);
- the placement, the run's kind and parameters (preview or verification, `--refine`,
  `--particles`, `--softening`, the end time);
- the engine: a fingerprint of the generator's own executable (any change of the physics,
  the level code, the dependencies or the compiler makes another binary, another key).

The runs are deterministic (the golden hashes), so a cached result is bit for bit the run's.
`--no-cache` bypasses it; `--verify-cache` reruns a sample and fails loudly on any
difference. Writes are atomic (a temporary file renamed: parallel `window` runs share the
directory). With it, JSON output from `window` and `check`. Each level below takes 5–15
design runs, repeated today across `tube-scan`, `check`, `window` and captures.

### Phase 1: levels with the engine as it is (about 10 levels, a week of sessions)

1. *Transistors, introduction:* **Gate.** One strip, a drain supply, a gate slider from
   accumulation (more carriers, more current) to pinch-off. Scales: the gate's potential
   in kT/q, the channel's width in Debye lengths. Risk: the channel charges through its own
   resistance (the prototype had not settled by t = 30); if long windows cost too much, an
   option to start the carriers in equilibrium with the electrodes (a run before t = 0), an
   engine change decided then.
2. *Transistors, intermediate:* **Amplifier** (a sine on the gate with its bias on a
   slider, a load strip: the drain swings more than the gate where the bias puts the
   steepest transconductance); **Inverter** (a pulse on the gate, a load strip: the output
   high, then low, two potential goals over two windows).
3. *Transistors, master:* **Flip-flop.** Two inverters cross-coupled by wires: a bistable
   cell, set by a pulse and holding its state, a bit of memory. *Deferred (2026-10-26):* a
   prototype (`target/scratch/ff_proto.py`: two of Gate's transistors, each pair of gates
   wired to the other's drain through a level-shifting source −S, the depletion gates
   needing −S ~ −70 to switch off; reset and set by dropping a drain's supply) met two
   obstacles. (1) Two 8-long channels hold 25 600 carriers in the verification, over the
   engine's 15 000 (`MAX_METAL_CARRIERS`, a share of `physics::tube::MAX_PARTICLES` 20 000);
   4-long channels fit but do not switch (at V_DD 60 0.73 off against 1.32 on: the short
   gates lose the channel). With the caps doubled the cell ran, but (2) it latches too
   slowly: each drain carries the other transistor's gates (~0.6, more through the Miller
   effect) through its load of 200, a time constant of ~200–400, so after a reset at
   S = 70 both drains crept back towards the symmetric point (A 50, B 54 at t = 100–150)
   and a set separated them by only ~7 by t = 300. A level would need runs of ~600–1000
   (previews of minutes). Options: gates shorter than the channel, placed by the source
   (less gate capacitance and a drift region shielding them from the drain), or a faster
   engine.
4. *Metals, intermediate:* **Wheatstone bridge.** Four strips joined by wires,
   magnetoresistive arms: the bridge balances at one |B|.
5. *AC and induction, introduction:* **Relaxation** (a pulse on a floating contact: the
   plasma rings at ω_p and decays; the amplitude goal at ω_p); **RLC with a strip** (a tank
   damped by a strip). *Intermediate:* **Field-tuned damping** (the strip
   magnetoresistive: the solenoid sets the tank's Q); **Lock-in Hall** (the swinging
   solenoid, a steady current: the Hall voltage's amplitude).
6. *Easy-win combinations* (the owner: an easy combination is done right away, not left for
   Phase 4): **Resistive anode** (a beam landing on a strip: its heat and IR drop) and
   **Calibrated field** (a Hall probe reads the solenoid that must turn a beam onto a
   collector), each after a prototype shows the engine joins them as it is.

### Phase 2: small extensions (about three weeks)

1. **An RMS goal** (`Measure::Rms`, 1–2 days; task 20). Test: equipartition, ⟨Q²⟩/2C =
   kT/2 on a floating contact, against D7's spectrum. *Heat and noise:* **Thermal noise**
   (a floating contact's kT/C, whatever the resistance feeding it); **Hot carriers** (a
   field that warms the carriers raises their noise).
2. **A drive's frequency on a slider** (1–2 days). Test: the AC conductivity through the
   level path against D2. *AC and induction:* **Carriers' inertia** (at ωτ = 1 the current
   lags 45°: a phase goal), **Tuned tank**, **Plasma resonance** (a floating plate driven
   at ω_p).
3. **Ring shapes** (a rectangle with a rectangular hole, contacts on its outer and inner
   sides; 2–3 days). Tests: the ring's background field (two rectangles' closed forms), the
   Corbino magnetoresistance against Drude's 1 + (ω_c τ)², its quarter-turn symmetry.
   **Corbino disc** (metals); **Generator** (a cut ring in the swinging solenoid, contacts
   either side of the cut on a load: Faraday's law driving a circuit; AC and induction).
4. **A uniform AC electric field** (`TubeSpec` gains E(t), the electrodes' response a
   further unit system, an AM waveform; 3–4 days). Tests: an isolated plate's polarization
   against the electrostatic solution; a ring's shielding against the Drude dielectric
   function ε(ω) = 1 − ω_p²/(ω(ω + i/τ)) (shielded below ω_p, passed above); the
   symmetries of linear response. **Polarization**, **Faraday cage** (a ring keeps the
   field out: floating contacts inside read nothing), **Receiving antenna** (a strip on a
   tuned tank).
5. **Temperatures** (a contact's temperature; then a lattice that warms cell by cell, its
   heat flowing to a bath; about a week). Tests: the model's Seebeck coefficient and
   thermal conductivity from its Boltzmann equation (derived first, Wolfram), the Kelvin
   relation Π = ST (Onsager's reciprocity, a symmetry test), the Lorenz number
   (Wiedemann–Franz). *Heat and noise:* **Seebeck**, **Peltier**, **Thermocouple** (Few or
   many's two densities), **Self-heating** (the resistance rising with temperature).

### Phase 3: new stages (about two to three months)

1. **The carriers' own magnetic field** (the Darwin limit in the z-invariant world, a
   finite c; about two weeks). The field along z from the current in the plane (the
   current density's stream function on a grid), the induced E (Darwin's transverse field,
   solved implicitly), the magnetic force. Tests first: a ring's inductance (a solenoid's,
   per length), nested rings' mutual inductance, a slab's skin depth, the magnetic energy
   ½LI², the symmetries (B odd under time reversal and the mirror). *AC and induction:*
   **Inductor** (a ring's current rising at L/R), **Transformer** (nested rings), **Skin
   effect**, **Shielding** (a ring keeps an AC field out), **Induction heating**; the
   eddy-current levels revisited for thick plates.
2. **Semiconductor devices** (two species in one lattice, a doping profile, generation and
   recombination; about three weeks). Tests: a junction's built-in potential (kT/q)
   ln(N_A N_D/n_i²) and depletion width, the Shockley diode law, detailed balance
   (np = n_i² at rest, a symmetry test), the minority carriers' diffusion length √(Dτ).
   *Transistors:* **Diode**, **LED** and **Photodiode** (recombination as light, light as
   generation), **Bipolar transistor**; a junction-transistor finale.
3. **Full electrodynamics for the carriers** (the tube world's FDTD field acting on them;
   three to four weeks). Tests: the plasma's dispersion ω² = ω_p² + c²k², a Drude slab's
   reflection (Fresnel with ε(ω)), a half-wave antenna's resonance. **Shiny metal**
   (reflects below ω_p), **Transparent metal** (passes above), **Half-wave antenna**;
   perhaps surface plasmons.

### Phase 4: joining (the backlog)

The owner (2026-10-20): more worlds will come as they are thought of, and joining them is
reserved for when there is nothing better to do; a combination that is an easy win is done
right away, as soon as its pieces exist (Phase 1's item 6 now; the Faraday cage with a beam
after Phase 2's items 3 and 4, the crystal radio after 2 and 4, the bolometer after 5).
Toward every element cooperating with
every other at each scale. First a cooperation
matrix: every element and material per world (the particles' slice, the tube world) against
every other, each pair marked works, refused (with the reason) or not yet, kept as a test
the way `scripts/pairs.py` keeps the slice's pairs, and the gaps it shows ranked by what
they would showcase. The levels already in view (besides Phase 1's item 6): **Faraday cage with a beam** (a beam inside a ring flies
straight while one outside wobbles), the AC arc's finale **Crystal radio** (antenna, tank,
vacuum diode, RC: the AM signal recovered), a **bolometer** (a beam's power read by the
heat it leaves on a strip). Further: the same materials at every scale where their
phenomena exist (a resistive conductor, a dielectric, a magnetic material in both worlds).

**For the owner, when each step starts** (none blocks Phase 1): the gate's settling
(long windows or a pre-settled start); the three new arcs or one larger metals arc; the c
of the metal levels for the magnetic stage; the order of the two later stages (proposed:
the magnetic field first).

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
- **Owner 2026-10-20: can the current's fluctuations from the carriers' thermal motion be
  measured in a resistor?** Yes, in a real one: Johnson–Nyquist noise, `⟨I²⟩ = 4kTΔf/R` (a
  shorted resistor) or `⟨V²⟩ = kT/C` across its own capacitance; the owner asked to note it
  as a level idea. In this model it stands as decided above: the noise is real for the
  model's line charges (D7 verifies `2kT/R` and `kT/C` with kT a macroparticle's) but not for
  a metal's carriers, and a goal on it would not survive the verification's finer
  macroparticles. What the game shows of it now: the metal maps' noise (each cell's current
  and power with its standard error, the carriers' thermal motion: PHYSICS.md §2.11), the
  floating contacts' wandering potentials, and the beads on the leads moving both ways. A
  level could still make it a goal if the verification kept the preview's macroparticles
  (refining the step only) and the goal were stated in the model's own kT (the noise
  against `√(kT/C)`, resolution-free): the owner to decide (task list).
- **Owner 2026-10-19: the Hall probe does not show that the field pushes the electrons
  aside.** Each strip in a field now carries the Hall effect's arrows: the carriers'
  drift, the magnetic push on it and the electric push back across the strip
  (`TubeRun::hall_forces`; test `hall_forces_balance`: opposite, the electric 1.022 of the
  magnetic once the charge has piled up, within the dielectric relaxation time ≈ 0.2).
  A force map (the force on the carriers, its electric and magnetic parts over the board)
  was proposed and left: the arrows answer the request.
- **Owner 2026-10-19: Ohm's law looked AC from its wiring; show the electrons' drift on
  DC, the traced ones constant in number.** The level is DC (supply 30, lead currents
  ±3.0 ± 0.3 throughout); with the supply off, as the level opens, the leads carry the
  thermal noise's current, both ways, and their beads wander. Chevrons across each metal
  now move with its carriers' mean displacement, the drift (0.38 per time unit in Ohm's
  law, `I/(ρW)` within 1.7 %: test `metal_drift_carries_the_current`), faint when it is
  noise; a traced carrier a contact takes is replaced by another drawn uniformly from the
  metal (40 throughout; replaced by the newest, just given off, they gathered at the
  contacts).
- **2026-10-19, Positive carriers (intermediate, after the Hall probe):** two strips in
  the player's solenoid, the upper of the levels' metal, the lower's carriers positive
  (`q/m = +1`, as zinc's or beryllium's behave), each with floating side contacts. Both
  tops above their bottoms need the currents opposite: the carriers then move the same
  way, and the field pushes electrons down and positive carriers up. Building it showed
  three faults of the display, fixed (PHYSICS.md §2.11):
  - the frames' field on the CPU (the Hall arrows, the field map's arrows) softened each
    carrier as `r⃗/(r² + ε²)`, never quite Coulomb's: a strip's thousands of carriers
    failed to cancel its background by `2λε²/r³` each, a field across the other strip
    (0.6 between the lower strip's contacts) that took its electric arrow to 1.7 times
    the magnetic. The carriers are now rods of radius 0.5, Coulomb's beyond, like the
    physics' thinner ones;
  - the electric push, taken on one frame along the centre line, scattered by more than
    its size on the weakly driven lower strip: it is now the field across the strip's
    middle, recorded by the run at every frame, averaged over the drift's time unit;
  - the inner contacts' leads ran through the other strip: leads now avoid metals (test
    `leads_cross_no_metal`).

  `hall_forces_balance` now covers all three strips: 1.049, 0.957, 1.288. The lower
  strip's ratio follows its weaker drive (driven at −30 like the upper: 1.043), not its
  carriers' sign.
