# Levels: curriculum and difficulty

`scripts/levels.py` is the single source of truth for the shipped levels: their order,
designs, reference solutions and automatic detector placement. Rebuild with
`python scripts/levels.py` (or `--only NN` for one level). Every level must have a
verified reference solution for every shot and negligible radiation
(`crates/level/tests/levels.rs`).

## Curriculum rules

- The levels form arcs (docs/CURRICULUM.md; `ARCS` in `scripts/levels.py`, written
  to `levels/curriculum.json`). Each arc has introductions (one new element or condition
  each), intermediate levels (the arc's modules interact, or an earlier design meets a
  real effect) and a master level (a multi-stage finale built from the arc's modules).
- Every new element or concept is introduced by an easy level first.
- Limits leave a margin of at least two elements above a level's fewest solution, so that
  players can find less-than-optimal solutions or just play around.
- Notes below name levels rather than numbers: the numbers change when levels are added.
- A good puzzle has a large configuration space, a small solution set, and a "distance to
  solution" that changes piecewise continuously. A player can then learn from attempts,
  but guessing or brute force is slow.

| # | Level | Arc, tier | Idea | Player elements | Reference |
|---|---|---|---|---|---|
| 1 | First bend | 1 Introduction | one charge bends a beam | ≤ 3 charges | 1 |
| 2 | Slingshot | 1 Introduction | bending around a level charge | ≤ 3 charges | 1 |
| 3 | Geiger–Marsden (1909) | 1 Introduction | Coulomb scattering, sign and strength | ≤ 3 charges | 1 |
| 4 | Twin beams | 1 Introduction | several shots, one setup | ≤ 3 charges | 1 |
| 5 | Two stages | 1 Introduction | a gate to pass before the detector counts | ≤ 3 charges | 1 |
| 6 | Injection | 1 Introduction | the detector also requires a direction (±8°) | ≤ 3 charges | 1 |
| 7 | Thomson's cathode-ray tube (1897) | 1 Intermediate | two energies onto one-cell spots; deflection ∝ 1/T | ≤ 3 charges | 1 |
| 8 | Around the wall | 1 Intermediate | the wall joined with injection: around the obstacle, then enter along the axis (±8°) | ≤ 4 charges | 2 |
| 9 | Einzel lens | 1 Intermediate | focusing of an angular spread | ≤ 4 charges | 2 |
| 10 | Reflectron | 1 Intermediate | reflection; energy-dependent turning point | ≤ 4 charges | 2 |
| 11 | Collimator | 1 Intermediate | three rays must arrive parallel (±3°) | ≤ 4 charges | 2 |
| 12 | Hemispherical analyzer (XPS) | 1 Intermediate | energy dispersion on a circular orbit | ≤ 3 charges | 1 |
| 13 | Soft landing | 1 Intermediate | arrival energy window: brake without reflecting | ≤ 3 charges | 1 |
| 14 | Sorting station | 1 Master | two energies × three angles: a lens through a gate, then sort the energies onto their own spots | ≤ 8 charges | 6 |
| 15 | High-voltage dome | 2 Introduction | a sphere held at a potential acts like a charge at its centre, but it is a conductor | ≤ 3 charges | 1 |
| 16 | Polarised sphere | 2 Introduction | an isolated neutral sphere becomes a dipole near a charge | ≤ 3 charges | 1 |
| 17 | Image charge | 2 Introduction | a grounded sphere answers every charge with an opposite image | ≤ 3 charges | 1 |
| 18 | Deflection plates | 2 Introduction | real plates with fringe fields (BEM) | ≤ 3 charges | 1 |
| 19 | Power supply | 2 Introduction | power supplies: set an electrode's potential instead of placing charges | 1 supply | 1 |
| 20 | Build a deflector | 2 Introduction | placing your own plates (electrodes) | ≤ 3 plates | 1 |
| 21 | Shielding | 2 Introduction | a grounded plate screens a charge's field | ≤ 3 plates | 1 |
| 22 | Glass screen | 2 Introduction | a dielectric (ε = 4) between a strong charge and the beam weakens its field but, unlike metal, does not cancel it; charges make up the rest (the reference fails without the glass) | ≤ 3 charges | 1 |
| 23 | Tune the lens | 2 Intermediate | focus a real Einzel lens with its voltage | 2 supplies | 2 |
| 24 | Real Einzel lens | 2 Intermediate | three apertures at high voltage focus three rays | ≤ 4 charges | 2 |
| 25 | Beam pipe | 2 Intermediate | the injection (6) built in; a grounded pipe wall between the beam and the steering charge screens it | ≤ 3 charges | 1 |
| 26 | Microscope column | 2 Master | tune the condenser through a crossover gate; deflector plate and charges onto a spot beside an ion pump, along the axis (±20°) | ≤ 5 charges, 3 supplies | 4 |
| 27 | Fast lane | 3 Introduction | γ changes the bending | ≤ 3 charges | 1 |
| 28 | First coil | 3 Introduction | magnetic fields from level coils | ≤ 3 charges | 1 |
| 29 | First magnet | 3 Introduction | placing your own magnets | ≤ 3 magnets | 1 |
| 30 | Ferrite shield | 3 Introduction | a ferrite plate (μ = 1000) between a strong magnet and the beam draws in its field (deflection 2.27 → 0.98 cells); magnets bring the beam in (the reference fails without the ferrite) | ≤ 3 magnets | 1 |
| 31 | Stern–Gerlach (1922) | 3 Introduction | neutral atoms with spin up or down, force m grad B_z | ≤ 3 magnets | 1 |
| 32 | Beta-ray spectrometer | 3 Intermediate | relativistic circular orbits: \|qQ\| = γmv²R, two electron energies | ≤ 3 charges | 1 |
| 33 | Dempster's mass spectrometer (1918) | 3 Intermediate | 180° focusing and mass separation | ≤ 3 charges | 1 |
| 34 | Wien filter | 3 Intermediate | crossed E (your charges) and B (a coil); the selected speed leaves straight (±3°) | ≤ 4 charges | 2 |
| 35 | Calutron | 3 Intermediate | isotope separation with magnets only | ≤ 4 magnets | 2 |
| 36 | Build a Wien filter | 3 Intermediate | the whole velocity selector from charges and magnets; straight exit (±3°) | ≤ 4 charges, ≤ 2 magnets | 3 |
| 37 | Mass spectrometer from parts | 3 Master | one electrostatic lens for three masses (only T/q matters), then a magnetic sector sorts them by momentum | ≤ 5 charges, ≤ 4 magnets | 5 |
| 38 | Stray field | 4 Introduction | a stray field switched on and off; aim between | ≤ 3 charges | 1 |
| 39 | RF kick | 4 Introduction | an antenna's kick depends on the passing time | ≤ 3 antennas | 1 |
| 40 | Ramp the coil | 4 Introduction | coil supplies: the player sets a coil's ramp rate; the ion passing beside it meets the field built by then (and the induced field), so the rate sets its deflection | 1 supply | 1 |
| 41 | Synchrotron light | 4 Introduction | radiative damping spirals the particle to the axis; only possible with radiation (canonical angular momentum) | ≤ 1 magnet | 1 |
| 42 | Mains hum | 4 Intermediate | uniform AC field at 4 phases acts like a random launch angle; imaging | ≤ 4 charges | 2 |
| 43 | Earth's field | 4 Intermediate | stray B_z of either sign on an electron beam | ≤ 4 charges | 2 |
| 44 | CRT in the Earth's field | 4 Intermediate | the CRT (7) built in, installed facing north in the Earth's field (stray B_z); adjusted where it stands | ≤ 3 charges | 1 |
| 45 | RF separator | 4 Intermediate | identical bunches half a period apart to different detectors | ≤ 1 charge, ≤ 2 antennas | 1 |
| 46 | Tune the RF | 4 Intermediate | choose the antenna frequency: bunches 4 apart meet it at phase difference ωΔt | ≤ 1 charge, ≤ 2 antennas | 1 |
| 47 | Streak camera | 4 Intermediate | three bunches a third of a period apart to three spots | ≤ 2 charges, ≤ 2 antennas | 1 |
| 48 | RF beam line | 4 Master | steer two bunches through a gate (±6°) with a stray field on and off, then separate them with RF | ≤ 5 charges, ≤ 3 antennas | 4 |
| 49 | Charging a plate | 5 Introduction | a circuit: the plate charges through a resistor, `V(t) = V_s(1 − e^{−t/RC})` (τ = 7.9), while the ions pass; −60k, enough when charged at once, now falls short | 1 supply | 1 |
| 50 | Chopper | 5 Introduction | a switch turns the deflector on between two launches: the first ion flies straight, the second must reach its own detector | 1 supply | 1 |
| 51 | Ringing plate | 5 Introduction | an inductor in series: the plate rings (LC), its first swing at 1.8 times the supply as the ions pass; only +120k lands 6 cells down | 1 supply | 1 |
| 52 | Induction kick | 5 Intermediate | the Faraday coil on a circuit: its current rises through R (τ = L/R = 20, L from its own flux); the induced field kicks a charge at rest, by the flux change, not its rate; T ≥ 0.3 needs the induction | ≤ 3 charges | 1 |
| 53 | Betatron | 5 Intermediate | Wideröe's 2:1 condition: a core coil on a circuit (its source rises, then holds) and a ramped guide coil; the orbit holds at R = 8 (±0.04) while T climbs to 4.6, then shrinks onto the detector (T ≈ 10, window 8.6–11.3) | 2 supplies (a source voltage, a ramp rate) | 2 |
| 54 | Pulse sorter | 5 Master | three ions launched at t = 0, 30, 60 through three plate pairs on an RC, a switched and an LC circuit: each ion meets the circuits at a different stage, so the three landings fix the three supply sliders (deflection fractions A 0.28/0.82/1, B 0/0.96/1, C 1/0.89/0.13) | 3 supplies | 3 |
| 55 | Space charge | 6 Introduction | 16 particles repel each other; ≥ 90 % must arrive, verified | ≤ 3 charges | 1 |
| 56 | Stern–Gerlach beam | 6 Introduction | both spin states as spread beams, each to its own detector | ≤ 3 magnets | 1 |
| 57 | Relativistic beam | 6 Introduction | a beam at 0.8c: magnetic attraction weakens space charge to 1/γ²; quasi-static interaction with radiation reaction | ≤ 3 charges | 1 |
| 58 | Collimated beam | 6 Intermediate | the collimator for a spread, interacting beam (±3°) | ≤ 4 charges | 2 |
| 59 | Velocity selector | 6 Intermediate | the Wien filter for a beam: three speeds, 8 ions each; the middle one leaves straight (±5°) | ≤ 4 charges | 2 |
| 60 | Beam preparation | 6 Intermediate | two stages for a beam: collimate it through a gate (±4°), then steer it into the target | ≤ 5 charges | 3 |
| 61 | Chromatic aberration | 6 Intermediate | the Einzel lens (9) built in; each ray becomes a beam with an 8 % energy spread | ≤ 3 charges | 1 |
| 62 | Real analyser | 6 Intermediate | the hemispherical analyser (12) built in; the source emits into a cone (σ = 4°) | ≤ 3 charges | 1 |
| 63 | Calutron at full current | 6 Intermediate | the calutron (33) built in; the isotope beams repel each other (5 ions each, radiation reaction) | ≤ 3 charges | 1 |
| 64 | Soft landing, full current | 6 Intermediate | the soft landing (13) built in; space charge grows as the ions are braked | ≤ 3 charges | 1 |
| 65 | Isotope separator | 6 Master | collimate an interacting two-isotope beam through a gate (±5°), then separate the isotopes with magnets | ≤ 5 charges, ≤ 4 magnets | 5 |
| 66 | Jackson §12.3: E×B drift | 7 Introduction | crossed uniform fields: both signs drift with E×B/B² along the equipotentials | ≤ 3 charges | 1 |
| 67 | Jackson Pr. 12.9: Van Allen equator | 7 Introduction | gradient drift around a dipole Earth, protons and electrons in opposite directions | ≤ 3 charges, ≤ 3 magnets | 1 |
| 68 | Throw a charge | 7 Introduction | a new element: the free charge (pull its slingshot handle back); its field pushes the particle and recoils, momentum passing without contact | ≤ 3 free charges | 1 |
| 69 | Jackson §5.15: Faraday's law | 7 Introduction | a ramped coil's induced field drives a charge at rest around the axis; the detector asks for energy only induction supplies | ≤ 3 charges | 1 |
| 70 | Jackson Pr. 13.1: knock-on | 7 Intermediate | the player throws a heavy ion (a free charge) past a light particle at rest, which must reach its detector with T ≥ 0.08: forward for close passes, sideways for distant ones | ≤ 3 free charges (mass 40) | 1 |
| 71 | Jackson §12.4: gradient drift | 7 Intermediate | drift along lines of equal \|B\|, opposite for the two signs; charges move both alike, magnets oppositely | ≤ 3 charges, ≤ 3 magnets | 2 |
| 72 | Jackson Pr. 12.5: E×B runaway | 7 Intermediate | \|E\| > c\|B\|: no drift frame, the particle runs away; magnets make B strong enough to drift | ≤ 6 magnets | 2 |
| 73 | Jackson §13.1: recoil at right angles | 7 Intermediate | the player throws a particle of equal mass at one at rest, which must reach its detector with T ≥ 0.1; the two leave at right angles | ≤ 3 free charges | 1 |
| 74 | Jackson §12.1: Störmer's forbidden region | 7 Intermediate | canonical angular momentum in a dipole's equatorial plane keeps the particle 7 cells out; charges break the symmetry | ≤ 4 charges | 2 |
| 75 | Jackson Ch. 12: magnetosphere | 7 Master | the solar wind's E×B drift steered up and down through two gates, then the gradient drift splits proton and electron around a dipole Earth | ≤ 6 charges, ≤ 5 magnets | 5 |
| 76 | Jackson §2.2: its own image | 8 Introduction | a grounded sphere attracts every passing charge through its image −qR/d | ≤ 3 charges | 1 |
| 77 | Jackson Pr. 2.6: two spheres | 8 Introduction | a charged and a neutral sphere image each other; the neutral one becomes a dipole | ≤ 3 charges | 1 |
| 78 | Jackson §3.13: field through a hole | 8 Intermediate | a grounded wall screens a charge except through the slot, where its field leaks out | ≤ 3 charges | 2 |
| 79 | Jackson Pr. 2.4: golden-ratio capture | 8 Intermediate | like charges attract inside 1.618 radii of an equally charged isolated sphere; go around it | ≤ 4 charges | 2 |
| 80 | Jackson §4.1: multipoles | 8 Intermediate | bend a near beam into its detector while two far beams (11–21 cells away, along x and y) arrive within 0.05° of straight: a net charge or a dipole that bends the near beam turns them; a compact quadrupole (+q, −2q, +q) does not | ≤ 5 charges | 3 |
| 81 | Jackson Ch. 2: sphere slalom | 8 Master | weave between three spheres carrying the particle's charge, through two gates, without being captured | ≤ 8 charges | 4 |
| 82 | Jackson Pr. 16.2: the classical atom | 9 Introduction | a radiating electron spirals into the nucleus; the detector counts only an electron slowed by its radiation | ≤ 3 charges | 1 |
| 83 | Jackson Pr. 16.3: orbits circularize | 9 Intermediate | a circular and an elliptic orbit of the same energy; the ellipse radiates most near the nucleus | ≤ 3 charges | 1 |
| 84 | Jackson Ch. 16: three orbits | 9 Master | three electrons on different orbits, all slowed by radiation into the detector in time | ≤ 6 charges | 3 |
| 85 | Jackson §16.7: a bound charge | 10 Introduction | an electron bound harmonically inside a charge cloud (Thomson's atom); pull it out | ≤ 3 charges | 1 |
| 86 | Jackson §16.8: resonance | 10 Introduction | a weak drive grows the bound electron's swing only at ω₀; tune an antenna | ≤ 2 antennas | 1 |
| 87 | Jackson Pr. 13.2: a kick for a bound charge | 10 Intermediate | the player throws a heavy negative ion past an atom; its field kicks the orbiting electron out into the detector | ≤ 3 free charges (mass 40) | 1 |
| 88 | Jackson Ch. 16: spectroscopy | 10 Master | two atoms with different natural frequencies; drive each at its own resonance | ≤ 4 antennas, ≤ 2 charges | 2 |
| 89 | Jackson §14.3: forward beaming | 11 Introduction | a radiation goal: a far receiver at 30° must collect ≥ 6e-4 per steradian; the radiation is beamed along the velocity (1/γ = 19°), so bend the particle while it heads for the receiver | ≤ 3 magnets | 1 |
| 90 | Jackson §14.6: the critical frequency | 11 Introduction | the same arena and receiver as forward beaming, but only ω 80–160 counts (≥ 3e-3 per steradian): a tight bend's short flash reaches ~(3/2)γ³c/R; the previous level's moderate bend puts only 1.4e-4 there | ≤ 3 magnets | 1 |
| 91 | Jackson Pr. 14.23: in step | 11 Introduction | a charge circles in crossed fields (ω₀ = 0.8) while its circle drifts into the detector; the receiver at 90° sees every particle: add a charge circling in step, whose field adds (four times the energy; ≥ 1.4e-4 per steradian in 0.6–1.0, the charge alone sends 4.84e-5) | ≤ 3 free charges | 1 |
| 92 | Jackson §14.2: a quiet turn | 11 Intermediate | turn 90° into the top detector with at most 6e-3 per steradian at 45°: gentle turns are quiet (∝ 1/R), or loop the other way round so the velocity never points at the receiver | ≤ 5 magnets | 1 (loop) |
| 93 | Jackson §14.8: Thomson scattering | 11 Intermediate | a light wave from the right (ω = 1) shakes the particle, which scatters it; the receiver at 30° counts 27–33, the Doppler-shifted line ω(1 + β cos φ)/(1 − β cos(θ − φ)) of a particle heading for it (≈ 4γ²ω head on) | ≤ 3 magnets (40, 80) | 1 |
| 94 | Jackson §16.8: why the sky is blue | 11 Intermediate | a steady receiver: a plane wave (ω = 0.35) lights a Thomson atom (ω₀ = 0.5, c = 2: radiation damping Γ = ω₀²τ = 0.021); the receiver behind the light averages the backscattered power in 0.32–0.38 over t = 400–1100, and wants 3–6 times the bare atom's (σ = 0.92 σ_T on Rayleigh's side): soften the spring along the light's field (a charge on the line of the swing weakens it; ω_y² = ω₀² − 4Q/d³ for a pair) | ≤ 4 charges | 1 |
| 95 | Jackson §15.2: braking radiation | 11 Intermediate | the target stops the particle abruptly; its flash ∝ sin²θ/(1 − β cos θ)² is zero straight ahead and largest at cos θ = β (≈ 20°): strike the target ~12–32° off the receiver's axis (band 100–200, ≥ 5e-3 per steradian) | ≤ 3 magnets (40, 80) | 1 |
| 96 | Jackson Pr. 14.23: a quiet ring | 11 Intermediate | the same drifting circle, and the receiver must stay below 7e-6 per steradian in 0.6–1.0: a partner on the far side of the circle, circling in antiphase about the same drifting centre, cancels the fundamental | ≤ 3 free charges | 1 |
| 97 | Jackson Pr. 14.23: a ring of four | 11 Intermediate | the band 0.6–2.6 holds ω₀, 2ω₀, 3ω₀ (below 1.3e-5 per steradian): a pair leaves 2ω₀, three at 120° leave 3ω₀; four evenly spaced radiate only at 4ω₀, outside the band | ≤ 5 free charges | 3 |
| 98 | Jackson §14.7: undulator | 11 Master | build an alternating row whose on-axis line 2γ²ω_u lands in the band 34–38 (spacing 3), strong enough (the line grows as N²), without steering the beam | ≤ 12 magnets | 10 |
| 99 | Vacuum diode | 12 Introduction | a z-invariant tube (PHYSICS.md §2.11) in its grounded envelope: the cathode emits under the space-charge limit; set the anode's supply so that 0.9–1.5 reaches it (verified 1.10 ± 0.04 at 26): the current grows as V^{3/2} | 1 supply (anode) | 1 |
| 100 | Sideways | 12 Introduction | the player's solenoid: above Hull's cutoff the electrons cannot reach the anode and drift along the cathode (E × B) into the collector beyond the gap's end (≥ 1.1: B in 0.9–1.3 of the slider ±4 in the preview, positive: out of the board); the wrong sign sends them down, too strong a field curls them back onto the cathode, which takes them (verified 1.18 ± 0.03 at 1.1) | 1 supply (the solenoid) | 1 |
| 101 | Bipolar diode | 12 Introduction | ions (q/m = 1/16) from the anode's face: their space charge cancels part of the electrons' and the electron current doubles (Vacuum diode's tube: 2.56 at 26 V against 1.29 in the preview; Langmuir's infinite planar diode: 1.86 times), but the ions strike the cathode: electrons ≥ 2.0 with ions ≤ 0.75 (both ∝ V^{3/2}, over t = 20–30) solve for 22.3–29.7 V of ±40 (preview; the middle verified); without the ions the electron goal would need 41 V; verified at 26: electrons 2.467 ± 0.131, ions 0.592 ± 0.003 | 1 supply (the anode, ±40) | 1 |
| 102 | Through the cloud | 12 Intermediate | an ion (q/m = 0.02, launched at 1.8) fired up through a working diode (60 V: a dense flow, the owner's request) must reach the top detector (x 8–14) while the anode current stays in 2.5–3.1: the field and the electron cloud push it out of the tube left of the detector; a positive plate below the gap, left of its path, kicks it back on course but draws electrons (verified at 60 V with a vertical plate at (9, 2) at +30: current 2.908 ± 0.015, the ion arrives at t = 21.23 and 21.20); too strong a plate draws too much current | ≤ 2 plates (±100; ±80 until its plate's window, 29.3–41.3, was found to hold the half's middle), the anode's supply (±100) | 2 |
| 103 | Electron gun | 12 Intermediate | a temperature-limited cathode (at most 1.6) and a slit: the beam through it spreads under its own charge and misses the small collector (0.155 ± 0.026 without plates); two plates beside the beam after the slit focus ≥ 0.24 onto it (+6 at (21, 13/7): 0.282 ± 0.015; +2 to +8 at x = 17–23 in the preview); at 0 the lens is too strong, at +12 too weak | ≤ 2 plates (±80) | 2 |
| 104 | Ion diode | 12 Intermediate | an inverted coaxial diode (ring cathode, rod anode emitting ions, a grounded envelope) in the player's solenoid: above Hull's cutoff the electrons stay in a sheath along the ring and draw more ions across, most just above it; ions onto the ring 5.5–8 with at most 0.05 of electrons leaking to the rod (t = 20–30) leave a diagonal strip (B ∝ √V): at (30 V, B = 1.9) V ~27–32 and B ~1.66–2.55 (preview); full voltage and field overshoot it. The leak is within both sliders (field off, or above ~50 V at B = 1.9). 8 times the macroparticles and emitters meshed twice as finely (PHYSICS.md §2.11); verified 6.43 ± 0.55; preview 2.9 s | 2 supplies (anode ±80, solenoid ±3) | 2 |
| 105 | Rectifier | 12 Intermediate | circuits in the tube (PHYSICS.md §2.11): a cathode between two anodes on AC supplies (period 20) half a period apart; the cathode is the output, through a load (R = 30, C = 1) to ground: its mean in 8.5–11.5 and its ripple below 1.2 (over t = 40–80) need both anodes driven alike, full wave, at about 24–36 each (preview); one anode (half wave) leaves a ripple of 2.6, opposite signs put the anodes in phase (ripple 5); verified at 32/32: mean 9.47 ± 0.25, ripple 0.764 ± 0.020 | 2 supplies (the anodes, ±80) | 2 |
| 106 | Triode | 12 Intermediate | the grid (PHYSICS.md §2.11): five wires 1.25 before the cathode on one sine (amplitude 1, period 40) around the bias the player sets; the anode through R = 60 to a fixed 100 swings far more than the grid: bias it so the anode's oscillation has an amplitude of 2.2–3.8 (verified at −0.7: 2.91 ± 0.26; about −0.9 to −0.4 verifies: at 0 the gain is too high, at −1 too low; cut off below about −2) | 1 supply (the grid's bias, ±2.6; ±2 until its preview window, −1.15 to −0.35, was found to hold the half's middle) | 1 |
| 107 | Klystron | 12 Intermediate | a weak beam (no space charge: ballistic bunching) from a Wehnelt gun passes a buncher (V₁ sin 2πt) and a tuned catcher 21 further; the player sets the cathode's −V₀ and V₁. The output amplitude follows Webster's 2I₀J₁(X) (the optimum drive within 14 % of the prediction from the buncher's own field, the rising side to 0.006; test `klystron_against_bunching_theory`), its phase moves 7.5° per volt (phase pushing). Backed off below saturation: amplitude 0.085–0.105 and phase −135° to −105° solve for −V₀ in −56.7 to −52.7 and V₁ in 10–12, or overbunched 27.3–32 (of ±80; the middles verified on the first); verified at (−54, 11): 0.0967 ± 0.0022, −122.6° ± 3.5° | 2 supplies (cathode, buncher; ±80) | 2 |
| 108 | Magnetron | 12 Master | a strapped eight-vane magnetron (rod cathode on the player's supply, the diagonal vanes strapped onto a tank near ω = 3 and onto an antenna outside the block, the player's solenoid) on a 96 × 96 board: the tank must oscillate with a swing of at least 10 over t = 90–120. The π mode starts where the electrons' drift keeps step with its wave (Buneman–Hartree, V rising with B): at (45 V, B = 6) the supply solves on 42.0–47.3 and the solenoid on 5.83–6.33. The antenna radiates what the tank's load absorbs (its radiation resistance is the load: the field solver's outflow matches R⟨I²⟩ to ~3 %); the magnetic and energy-flow maps show the waves leaving. Verified 16.80 ± 0.43; preview 2.1 s, verification 8.1 s | 2 supplies (cathode ±80, solenoid ±10) | 2 |
| 109 | Ohm's law | 13 Introduction | a metal (PHYSICS.md §2.11, the Drude model): a strip 16 × 4 between two contacts, the right one on the player's supply. The electrons zigzag between the ion cores at the thermal speed; the supply adds a drift, a current in proportion to the voltage (I/V from 0.1135 at 8 to 0.096 at 40: the carriers warm). Goal 2.6–3.2 into the right contact over t = 10–40: the supply solves on 26–32. Verified 3.03 ± 0.04; the Joule heat shown (the lattice takes 81 of I V = 90; the carriers bring the rest into the contact) | 1 supply (±40) | 1 |
| 110 | Long and thin | 13 Introduction | the geometry of a resistor: two strips of the same metal in series, the upper 16 × 4, the lower 8 × 8 (R ∝ L/W: four times), the upper strip's left end on supply 1, the junction (the upper strip's right end strapped to the lower strip's left end) on supply 2, the lower strip's right end grounded: a current 2.0–2.8 through each. Measured R = 10.1 and 2.37 (a ratio of 4.3: the thin strip's stronger field warms its carriers); the textbook ratio (V₂ = V₁/5) passes too. At the reference (32, 5.5) supply 1 solves on 28–38, the junction on about 4.7–6.6 (on the 2 V grid only 6). Verified 2.41 ± 0.03 and 2.43 ± 0.08 | 2 supplies (±40) | 2 |
| 111 | Potentiometer | 13 Introduction | a strip 18 × 4, its left end grounded, the right end and a side contact a third of the way along on the player's supplies: a current 2.4–3.0 along the strip and none (±0.5) through the side contact over t = 10–80, which must sit at the strip's potential there (V₁/3 in a uniform field; the null measured at 9.05 for 30: the contacts take about a volt). At the reference (30, 9) the end solves on 28–32 and the side contact on 8–10. Verified 2.73 ± 0.06 and −0.09 ± 0.31 (the floating side's thermal noise) | 2 supplies (±40) | 2 |
| 112 | Hall probe | 13 Intermediate | a strip 16 × 4 carrying a current in the player's solenoid, floating side contacts at mid-length (C = 1): the Lorentz force piles electrons on one side until the Hall field holds them (V_H ≈ I B/nq, below it as ω_c τ = 0.4 B grows). The solenoid spans its full range (±16: from straight paths to circles 0.12 across, ω_c τ = 6.4; the owner's showcase rule, 2026-10-20); over it at 30 the current falls 3.14 → 1.96 and the Hall voltage grows ever more slowly (−1.05 at B = 1, −4.67 at 4, −8.12 at 8, −12.52 at 16). Goals: 2.7–3.5 into the right end, the top contact 3.5–6.0 below the bottom one over t = 10–80 (B out of the board). At the reference (30, B = 4.2) the solenoid solves on 3.2–5.2 (6.3 % of ±16) and the supply on 29–36 (8.8 % of ±40). With the slider at ±4 the goal had been 1.0–1.9 below (B 1.0–1.4) | 2 supplies (strip ±40, solenoid ±16) | 2 |
| 113 | Positive carriers | 13 Intermediate | two strips 16 × 4 in the player's solenoid, each with floating side contacts at mid-length (C = 1): the upper of the levels' metal (electrons), the lower's carriers positive (q/m = +1, as zinc's and beryllium's behave; both strips at 50 macroparticles per area). For a given current the Hall voltage's sign follows the carriers': both tops above their bottoms need the two currents opposite. The solenoid spans ±16 (the owner's showcase rule); over it at (30, B, −12) the upper Hall voltage 1.12, 4.69, 7.63, 11.99 and the lower 0.73, 2.05, 3.15, 4.92 at B = −1, −4, −8, −16. Goals: 2.6–3.5 into the upper strip's right end, its top contact 3.4–6.0 and the lower strip's 1.3–2.6 above the bottom ones over t = 10–80: B into the board, the lower supply negative. At the reference (30, B = −4, −12) the solenoid solves on −5.6 to −3.2 (7.5 % of ±16), the upper supply on 27–37 (12.5 % of ±40), the lower on −15 to −8 (8.8 %). (At ±7 the goals had been 2.7–3.5, 1.6–2.8, 0.6–1.6, B −2.5 to −1.5) | 3 supplies (strips ±40, solenoid ±16) | 3 |
| 114 | Magnetoresistance | 13 Intermediate | a strip 16 × 4 and one 4 × 10 (short along the current, wide across it) on one supply in the player's solenoid: in the long strip the Hall field holds the carriers back and they flow straight; the short strip's contacts short it, the carriers' paths bend and its resistance grows with |B|, either sign (geometric magnetoresistance). The solenoid spans ±16 (the owner's showcase rule): at 7 V the long strip carries 0.795, 0.683, 0.573, 0.452 and the short 6.96, 4.63, 2.68, 1.18 at B = 1, 4, 8, 16 (at full strength the carriers circle on the spot between collisions). Goals over t = 10–80: 0.56–0.80 into the long strip's right end, 3.2–4.8 into the short one's (3.7–4.5 at ±7). At B = 5 both currents go as the supply (long 0.096 V, short 0.58 V), so the supply's window is the overlap of the two goals' ratios: with 0.58–0.78 and 3.2–4.6 it solved on 6.5–7.5 only (5 %). At the reference (7, B = 5) the solenoid solves on 4.0–6.4 (7.5 % of ±16; its mirror too: the effect is even in B) and the supply on 6.5–8.0 (7.5 % of ±10) | 2 supplies (strips ±10, solenoid ±16) | 2 |
| 115 | Eddy currents | 13 Intermediate | an isolated plate 12 × 12 on the solenoid's axis, its field swinging as B₀ sin(0.125 t) (`TubeSpec::solenoid_drive`): the induced field drives the carriers round in loops, eddy currents, whose heat goes as `Ḃ²` (the heat goal measures the induced field's work, PHYSICS.md §2.11). The board's only electrode is a small floating cold one (it emits nothing; isolated metals need every electrode floating, OP-1). Local Ohm's law gives `σ ⟨Ḃ²⟩ ∫u dA` = 2.68 B₀²; measured heat over B₀² 2.31, 2.51, 2.47, 2.48, 2.39 at B₀ = 0.5, 1, 1.5, 2, 3 (5–7 % below Ohm's law: the walls; the carriers barely warm even at 3). The swing spans ±3 (the owner's showcase rule: heat over a factor of 36). Goal 1.2–4.0 over t = 25–125 (two periods): EDDY_WINDOW. Verified 2.30 ± 0.28 at the reference B₀ = 0.95 | 1 supply (solenoid swing ±3) | 1 |
| 116 | Laminations | 13 Intermediate | a transformer's core 20 × 12 (30 macroparticles per area) on the solenoid's axis in the level's own swinging field (B₀ = 0.68, ω = 0.125), three cuts of any length and angle, drawn by dragging (`ElementKind::Cut`: walls inside the metal). Each eddy loop must go round the cuts: measured (local Ohm's law) uncut 2.833 (3.060); three full cuts along, evenly 0.516 (0.277), at y = 7, 11, 14 0.597, at 6, 10, 14 0.751; two along and one across 0.407; two along 0.755 (0.475); one along and one across 0.791 (0.765); one along 1.260 (0.991); corner to corner 1.435; one across 1.633 (1.692); two half-length cuts along 1.745, two staggered 16 long 1.737. The long cuts' faces carry current across as capacitors, so the gain along them falls short of Ohm's law. Goal at most 0.65 over t = 25–125: all three cuts across the whole core, about evenly, mostly along it. Verified at the reference (two along, one across) | 3 cuts | 3 |

## Automatic detectors

A shot's detector may be given as `Auto(strip, axis, size)`. The script then runs the
reference solution (or, if the level has none, first asks the solver for a setup that
brings every shot into its probe strip). It places a `size`-cell detector around the
landing point, clamped to the arena. Shots of the same species (same particle and energy,
differing only in direction, as in Dempster's 180° focusing) share one detector.
Overlapping detectors of different species on the same screen are split at the node
halfway between their landing points. If the points are too close to split, the build
fails.

## Difficulty measurements

`cargo run --release -p generator -- analyze [--fewest] levels/[0-9]*.json` (by default 2000 random samples,
32 search runs with a budget of 400 evaluations each; the table below used 1000 and 16).
Columns:

- **log10 configs:** number of distinct player placements within the limits.
- **random solve rate:** fraction of uniformly random placements that solve every shot,
  verified. When none were found, the upper estimate 3/n (rule of three) is shown.
- **search:** a local search that only sees the distance objective, like a player watching
  the preview (greedy moves with occasional worse moves accepted). Shows its success rate
  and the mean evaluations of successful runs.
- **expected effort:** expected number of tried placements for the search (restarting after
  each failed run) versus for random guessing.
- **smoothness:** Spearman rank correlation of the distance objective between a placement
  and a one-move neighbour. High means that attempts carry information.

Measured 2026-09-28 after the curriculum restructure, with 1000 random samples and 16
search runs of 400 evaluations per level; levels 61–92 (the Jackson series, added later)
on 2026-10-05 and levels 46–49 (the circuits) on 2026-10-06 with the same settings:

| level | log10 configs | random solve rate | search success | mean evals | expected effort: search / guessing | smoothness |
|---|---|---|---|---|---|---|
| 01_first_bend | 10.0 | 7.1e-2 | 16/16 | 10 | 10 / 14 | 0.77 |
| 02_slingshot | 10.0 | 2.1e-2 | 16/16 | 57 | 57 / 48 | 0.71 |
| 03_geiger_marsden | 9.5 | 1.7e-2 | 15/16 | 42 | 69 / 59 | 0.41 |
| 04_twin_beams | 10.0 | 6.0e-3 | 16/16 | 56 | 56 / 167 | 0.87 |
| 05_two_stages | 10.0 | 3.5e-2 | 16/16 | 58 | 58 / 29 | 0.79 |
| 06_injection | 10.4 | 3.2e-2 | 16/16 | 41 | 41 / 31 | 0.80 |
| 07_thomson_crt | 8.8 | 1.3e-2 | 15/16 | 82 | 109 / 77 | 0.75 |
| 08_around_the_wall | 13.5 | < 3.0e-3 | 3/16 | 61 | 1794 / 333 | 0.78 |
| 09_einzel_lens | 10.3 | < 3.0e-3 | 15/16 | 125 | 151 / 333 | 0.84 |
| 10_reflectron | 10.1 | 2.0e-3 | 16/16 | 108 | 108 / 500 | 0.84 |
| 11_collimator | 11.1 | < 3.0e-3 | 16/16 | 84 | 84 / 333 | 0.89 |
| 12_hemispherical_analyzer | 9.5 | 2.0e-3 | 9/16 | 97 | 408 / 500 | 0.72 |
| 13_soft_landing | 9.7 | 1.0e-3 | 2/16 | 174 | 2974 / 1000 | 0.81 |
| 14_sorting_station | 26.1 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.93 |
| 15_high_voltage_dome | 10.3 | 1.2e-1 | 16/16 | 12 | 12 / 8 | 0.76 |
| 16_polarized_sphere | 9.5 | 1.4e-2 | 16/16 | 30 | 30 / 71 | 0.65 |
| 17_image_charge | 10.3 | 2.1e-2 | 15/16 | 58 | 85 / 48 | 0.83 |
| 18_deflection_plates | 8.9 | 1.6e-1 | 16/16 | 5 | 5 / 6 | 0.75 |
| 19_power_supply | 0.7 | 2.7e-1 | 16/16 | 3 | 3 / 4 | -0.43 |
| 20_build_a_deflector | 8.7 | 1.2e-1 | 16/16 | 12 | 12 / 8 | 0.65 |
| 21_shielding | 7.1 | 3.7e-2 | 16/16 | 16 | 16 / 27 | 0.61 |
| 22_tune_the_lens | 1.4 | 4.2e-2 | 11/16 | 23 | 205 / 24 | 0.30 |
| 23_real_einzel_lens | 10.8 | 1.6e-2 | 16/16 | 42 | 42 / 62 | 0.84 |
| 24_beam_pipe | 9.6 | 3.2e-2 | 16/16 | 20 | 20 / 31 | 0.81 |
| 25_microscope_column | 16.2 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.80 |
| 26_fast_lane | 10.7 | 1.1e-2 | 16/16 | 79 | 79 / 91 | 0.67 |
| 27_first_coil | 10.0 | 4.3e-1 | 16/16 | 4 | 4 / 2 | 0.44 |
| 28_first_magnet | 9.1 | 7.7e-2 | 16/16 | 18 | 18 / 13 | 0.68 |
| 29_stern_gerlach | 8.8 | 6.8e-2 | 16/16 | 15 | 15 / 15 | 0.80 |
| 30_beta_spectrometer | 9.5 | 2.0e-3 | 8/16 | 117 | 517 / 500 | 0.62 |
| 31_dempster | 8.4 | 8.0e-2 | 16/16 | 12 | 12 / 12 | 0.56 |
| 32_wien_filter | 11.6 | < 3.0e-3 | 5/16 | 188 | 1068 / 333 | 0.66 |
| 33_calutron | 12.3 | < 3.0e-3 | 5/16 | 101 | 981 / 333 | 0.82 |
| 34_build_wien_filter | 18.2 | < 3.0e-3 | 2/16 | 332 | 3132 / 333 | 0.82 |
| 35_mass_spectrometer | 30.1 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.89 |
| 36_stray_field | 10.7 | 5.1e-2 | 16/16 | 32 | 32 / 20 | 0.78 |
| 37_rf_kick | 12.2 | 6.1e-2 | 16/16 | 19 | 19 / 16 | 0.65 |
| 38_synchrotron_light | 1.0 | 5.3e-1 | 16/16 | 2 | 2 / 2 | -0.08 |
| 39_mains_hum | 12.3 | < 3.0e-3 | 9/16 | 131 | 442 / 333 | 0.88 |
| 40_earths_field | 12.4 | 1.0e-3 | 5/16 | 218 | 1098 / 1000 | 0.81 |
| 41_crt_earth_field | 9.7 | 5.0e-3 | 15/16 | 68 | 95 / 200 | 0.80 |
| 42_rf_separator | 11.0 | 4.0e-3 | 15/16 | 88 | 115 / 250 | 0.76 |
| 43_tune_the_rf | 12.1 | 2.0e-3 | 15/16 | 63 | 89 / 500 | 0.72 |
| 44_streak_camera | 14.2 | < 3.0e-3 | 7/16 | 194 | 708 / 333 | 0.85 |
| 45_rf_beam_line | 28.7 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.90 |
| 46_charging_a_plate | 0.7 | 2.4e-1 | 16/16 | 4 | 4 / 4 | -0.39 |
| 47_chopper | 0.7 | 2.7e-1 | 16/16 | 4 | 4 / 4 | -0.42 |
| 48_ringing_plate | 0.7 | 2.7e-1 | 16/16 | 4 | 4 / 4 | -0.42 |
| 49_induction_kick | 7.8 | 9.7e-2 | 16/16 | 32 | 32 / 10 | 0.62 |
| 50_space_charge | 10.0 | 2.1e-2 | 15/16 | 59 | 86 / 48 | 0.80 |
| 51_stern_gerlach_beam | 8.8 | 3.3e-2 | 14/16 | 17 | 74 / 30 | 0.82 |
| 52_relativistic_beam | 10.4 | 2.8e-2 | 15/16 | 23 | 50 / 36 | 0.75 |
| 53_collimated_beam | 11.1 | 1.0e-3 | 16/16 | 96 | 96 / 1000 | 0.88 |
| 54_velocity_selector | 11.6 | < 3.0e-3 | 1/16 | 306 | 6306 / 333 | 0.68 |
| 55_beam_preparation | 14.2 | < 3.0e-3 | 7/16 | 274 | 789 / 333 | 0.91 |
| 56_chromatic_aberration | 9.8 | 6.0e-3 | 14/16 | 43 | 101 / 167 | 0.87 |
| 57_real_analyzer | 8.6 | 1.3e-2 | 12/16 | 17 | 151 / 77 | 0.74 |
| 58_calutron_space_charge | 9.7 | 1.0e-3 | 9/16 | 99 | 410 / 1000 | 0.81 |
| 59_soft_landing_current | 9.1 | < 3.0e-3 | 12/16 | 200 | 334 / 333 | 0.74 |
| 60_isotope_separator | 28.6 | < 3.0e-3 | 1/16 | 187 | 6187 / 333 | 0.91 |
| 61_jackson_exb_drift | 10.4 | 3.0e-3 | 16/16 | 90 | 90 / 333 | 0.65 |
| 62_jackson_van_allen | 20.7 | 3.0e-3 | 16/16 | 70 | 70 / 333 | 0.62 |
| 63_throw_a_charge | 14.2 | 7.0e-3 | 16/16 | 92 | 92 / 143 | 0.84 |
| 64_jackson_faraday | 7.8 | 8.0e-3 | 16/16 | 60 | 60 / 125 | 0.62 |
| 65_jackson_knock_on | 11.6 | 2.0e-3 | 5/16 | 236 | 1116 / 500 | 0.90 |
| 66_jackson_gradient_drift | 20.0 | 1.0e-3 | 12/16 | 147 | 281 / 1000 | 0.65 |
| 67_jackson_runaway | 18.5 | < 3.0e-3 | 13/16 | 110 | 202 / 333 | 0.81 |
| 68_jackson_recoil | 11.6 | 1.4e-2 | 15/16 | 131 | 158 / 71 | 0.84 |
| 69_jackson_stormer | 13.5 | < 3.0e-3 | 9/16 | 229 | 540 / 333 | 0.71 |
| 70_jackson_magnetosphere | 36.0 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.80 |
| 71_jackson_own_image | 10.3 | 4.3e-2 | 16/16 | 18 | 18 / 23 | 0.84 |
| 72_jackson_two_spheres | 9.8 | 8.0e-3 | 15/16 | 68 | 94 / 125 | 0.77 |
| 73_jackson_slot | 9.0 | 2.0e-3 | 15/16 | 239 | 266 / 500 | 0.71 |
| 74_jackson_golden_ratio | 13.4 | 1.0e-3 | 10/16 | 87 | 327 / 1000 | 0.66 |
| 75_jackson_multipoles | 10.5 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.69 |
| 76_jackson_sphere_slalom | 25.0 | < 3.0e-3 | 2/16 | 157 | 2957 / 333 | 0.77 |
| 77_jackson_classical_atom | 10.1 | 5.1e-2 | 16/16 | 17 | 17 / 20 | 0.62 |
| 78_jackson_circularization | 10.1 | < 3.0e-3 | 7/16 | 210 | 724 / 333 | 0.67 |
| 79_jackson_three_orbits | 19.0 | < 3.0e-3 | 2/16 | 219 | 3019 / 333 | 0.71 |
| 80_jackson_bound_charge | 10.1 | 8.0e-3 | 16/16 | 101 | 101 / 125 | 0.60 |
| 81_jackson_resonance | 8.0 | 2.0e-3 | 16/16 | 97 | 97 / 500 | 0.86 |
| 82_jackson_bound_knock | 11.6 | 7.0e-3 | 15/16 | 93 | 120 / 143 | 0.76 |
| 83_jackson_spectroscopy | 25.5 | < 3.0e-3 | 14/16 | 194 | 251 / 333 | 0.75 |
| 84_jackson_beaming | 10.1 | 5.1e-2 | 16/16 | 81 | 81 / 20 | 0.59 |
| 85_jackson_critical_frequency | 10.1 | 1.4e-2 | 12/16 | 107 | 240 / 71 | 0.62 |
| 86_jackson_in_step | 13.2 | 4.1e-1 | 16/16 | 9 | 9 / 2 | 0.78 |
| 87_jackson_quiet_turn | 15.8 | 3.0e-3 | 12/16 | 165 | 298 / 333 | 0.79 |
| 88_jackson_thomson | 9.2 | 2.2e-2 | 15/16 | 94 | 120 / 45 | 0.60 |
| 89_jackson_cross_section | 13.2 | 2.8e-2 | 31/32 | 60 | 73 / 36 | 0.73 |
| 90_jackson_braking | 9.2 | 6.4e-2 | 16/16 | 30 | 30 / 16 | 0.63 |
| 91_jackson_quiet_ring | 13.2 | 3.0e-3 | 3/16 | 107 | 1840 / 333 | 0.82 |
| 92_jackson_ring_of_four | 21.1 | < 3.0e-3 | 0/16 | NaN | inf / 333 | 0.87 |
| 93_jackson_undulator | 24.6 | < 3.0e-3 | 5/16 | 123 | 1003 / 333 | 0.52 |
| 99_vacuum_diode | 0.8 | 1.9e-1 | 8/8 | 6 | 6 / 5 | 0.19 |
| 100_sideways | 0.8 | 2.1e-1 | 8/8 | 6 | 6 / 5 | 0.29 |
| 101_bipolar_diode | 0.8 | 1.9e-1 | 8/8 | 6 | 6 / 5 | 0.24 |
| 102_through_the_cloud | 7.7 | < 1.0e-2 | 5/8 | 145 | 385 / 100 | 0.66 |
| 103_electron_gun | 7.2 | 1.0e-2 | 8/8 | 133 | 133 / 100 | 0.78 |
| 104_ion_diode | 1.6 | 7.3e-2 | 8/8 | 23 | 23 / 14 | 0.43 |
| 105_rectifier | 1.6 | 8.0e-2 | 8/8 | 18 | 18 / 12 | 0.48 |
| 106_triode | 0.8 | 1.9e-1 | 8/8 | 3 | 3 / 5 | 0.33 |
| 107_klystron | 1.8 | 3.0e-2 | 8/8 | 36 | 36 / 33 | 0.21 |
| 108_magnetron | 1.6 | 1.1e-1 | 8/8 | 26 | 26 / 9 | 0.34 |
| 109_ohms_law | 0.8 | 1.9e-1 | 8/8 | 4 | 4 / 5 | -0.05 |
| 110_long_and_thin | 1.6 | 3.0e-2 | 8/8 | 34 | 34 / 33 | 0.58 |
| 111_potentiometer | 1.6 | 3.0e-2 | 6/8 | 66 | 200 / 33 | 0.30 |
| 112_hall_probe | 1.6 | 6.3e-2 | 8/8 | 19 | 19 / 16 | 0.28 |
| 113_positive_carriers | 2.3 | 1.0e-2 | 8/8 | 128 | 128 / 100 | 0.61 |
| 114_magnetoresistance | 1.5 | 6.3e-2 | 6/8 | 36 | 169 / 16 | 0.28 |

Reading the table by tier (levels 1–56): the introductions are solved by the search in
15–16 of 16 runs within 3–80 evaluations. The intermediate levels need more (Around the wall 3/16, Soft
landing 2/16, Build a Wien filter 2/16, Velocity selector 1/16). The master levels are not
found by a search over the whole level (0/16, the Isotope separator 1/16), since they
are meant to be built stage by stage, as their references were; random guessing of the
whole level is hopeless (log10 configs 16–30). Their stages are each within reach of the
modules the arc taught (see the finale notes below).

The Jackson series (57–89): the search solves most introductions in 15–16 of 16 runs
(Faraday's law, the classical atom and the bound charge among them); none of its 16 runs
solved Magnetosphere (66, the arc's master), Multipoles (71) or A ring of four (88), and
random guessing found none of them either. The search sees only the distance objective
(plus the deficits of the conditions), so these numbers also measure how much a goal's
conditions (Multipoles' 0.05° direction windows, a radiation band, gates) hide from it,
not only how hard the physics is.

Notes:

- The introductions are meant to be easy. Synchrotron light is a concept level (8
  choices).
- Under disturbances a placement counts as a solution only if every flight arrives,
  verified. The distance objective sums over all flights.
- Dempster has in effect one physical parameter (the accelerating voltage, which sets all
  radii). It is kept as the concept level for 180° focusing.
- The Calutron is among the hardest intermediate levels: no random placement out of 2000
  solved it, and the search succeeded in 4 of 32 runs. Its difficulty comes from the
  magnets' strongly non-uniform dipole fields.
- Power supply and Tune the lens are knob levels. Since 2026-10-06 every power supply and
  player plate is a slider over a range symmetric about 0 (the owner); the analysis still
  samples the listed candidates (7 and 49 settings). The analysis caches the result of every distinct placement, so such levels
  cost only their distinct placements; before the cache (and before it knew about
  supplies), analysing them did not terminate.
- The circuit levels (46–48) are knob levels on the plates of Power supply: Charging a
  plate punishes the static answer (−60k, which reaches 2 cells when charged at once, falls
  short at 1.2 while charging), Ringing plate rewards the overshoot, and Chopper sends the
  first ion straight in any setting, before the switch closes.
- Sliders (2026-10-06): `generator window <level>` scans each supply or plate of the
  reference over its range and prints the solving intervals. Each level's range is chosen
  so that no interval contains 0, either end or the middle of either half (a player's
  first tries); `levels.py` and the test `slider_ranges_hide_the_answer` check the
  references. Solving intervals: Power supply −155k…−85k on ±340k, Build a deflector
  +45k…+100k on ±230k (reference +75k), Tune the lens +280k…+350k per half on ±490k,
  Microscope column condenser +385k…+435k and deflector +415k…+620k on ±710k, Charging
  a plate −140k…−85k on ±310k, Chopper +90k…+165k on ±360k, Ringing plate +107.5k…+170k
  on ±205k (reference +140k), Pulse sorter +113k…+173k, +47k…+130k and −407k…−287k on
  ±460k. Ramp the coil (a coil supply: ramp rates) −6.33k…−4.17k on ±21k. Betatron: core source 130M…142.5M on ±310M (the energy is its flux's), guide rate 1.48k…3.39k (islands to 3.59k) on ±8.3k.
  Induction kick (49) is Jackson's Faraday level with the coil on a circuit (c = 1e4,
  where the circuit's quasi-static parameter is 1.5e-4): the reference's one charge
  brings the ion in with T = 0.37 against the 0.3 required.
  Pulse sorter (50) is the arc's master: three supplies, three ions at different launch
  times; designed by measuring the deflection matrix and solving for landings 12.5, 9.5
  and 6.5 (plate gaps of 8 cells keep every flight clear of the plates).
- Real instruments (Beam pipe, CRT in the Earth's field, Chromatic aberration, Real
  analyser, Calutron at full current, Soft landing at full current) build an earlier
  level's reference setup in as fixed elements and add one real effect the ideal design
  ignored. They sit in the arc of the effect they add, after the level they build on.
  `scripts/levels.py` checks at build time that the built-in design alone fails, and the
  test `realistic_levels_break_the_idealised_design` keeps checking it. Measured with the
  built-in design alone: Chromatic aberration one ray 7/8 arriving; Real analyser the
  highest energy 7/8; CRT in the Earth's field both beams lost; Beam pipe rejected (wrong
  direction); Calutron at full current the light isotope 3/5; Soft landing at full current
  0/8. Tuning for playability (no search solved the first versions): the CRT first had to
  work facing north, south and shielded with one-cell spots, and at most 1 in 32 searches
  solved it even at 40 % of the field; a tube is adjusted where it stands, so the level now
  has one orientation. Soft landing at full current had 12 ions at the full charge scaling
  (0/32); now 8 at half.
- Beam levels are harder versions of single-particle levels: a placement solves them only
  if the required share of the beam arrives, verified. Their objective for the search sums
  the detector distances (plus gate shortfalls) of the closest missing particles.
  Relativistic beam includes radiation reaction because its particles would otherwise
  neglect 5.3e-6 of their energy (above the 1e-10 allowed); the test
  `radiation_levels_need_radiation` accepts that reason for beams. Beam preparation
  combines the collimated beam with a second stage; the search found its reference in 44 s
  and solves it in 9 of 32 runs.
- The analysis logs its progress to stderr (`[analyze]` lines: phases with timing, distinct
  placements flown, the slowest evaluation, and every evaluation slower than 1 s with its
  placement). The log found why the velocity selector's analysis ran for hours (placements
  trapping the ions for 1e6 steps, and a beam bug); searches now stay within the step
  budget of the cost meters.
- Finales (master levels) have references built stage by stage (`solve_stage` in
  `scripts/levels.py`): each stage is solved as a sub-level whose detector is the next
  gate, with the earlier stages fixed, as a player would build it.
  - Sorting station: the lens (stage 1, at most 3 charges left of the gate) needs 2
    charges; sorting and refocusing (stage 2, at most 4 charges) needs 4. The build checks
    that the search finds no solution with 4 charges anywhere (`FINALES`), and the
    whole-level search with up to 8 charges finds none either (0/4 runs). Versions that
    asked for parallel arrival (within ±6°, ±10° or ±15°) had no staged solution;
    five-cell spots were solved with 3 charges in stage 2 (too easy).
  - Microscope column: the condenser (two supplies, both at 400 kV) and then the deflector
    supply with one charge: 4 elements. The spot beside the ion pump asks for arrival
    along the axis (±20°): stage 2 is solved in 4 of 48 search runs; at ±10° in none of
    16. The first version (8 full-size electrodes) took 8.5 s to set up, over the sandbox
    budget, and was solved with 3 elements; smaller, lower electrodes take 0.75 s.
  - Mass spectrometer from parts: one charge makes the lens for all three masses (only T/q
    matters for electrostatic optics); the magnetic sector takes 4 elements. The
    collectors are placed automatically at the reference's landing points.
  - RF beam line: stage 1 (through a two-cell gate within ±6°, with the stray field on and
    off) takes one charge; the RF separation onto two-cell spots takes 3 elements. Looser
    versions (a four-cell gate without direction, ±8°, half the stray field, four-cell
    spots) were solved with 3 elements in all.
  - Jackson Ch. 12: magnetosphere: lifting both particles through the upper gate takes 2
    charges, bringing them down through the lower one 2 more, and the drift shells to the
    satellites behind the Earth one element. With one gate the level took 3 elements; with
    the satellites swapped (each particle three quarters around the Earth) stage 2 had no
    solution with up to 3 charges and 4 magnets, since the solar wind sweeps the
    particles off their drift shells.
  - Jackson Ch. 2: sphere slalom: 2 charges to the first gate, one each for the second
    gate and the detector. With two spheres and one gate the level took 3 elements.
  - Jackson §4.1: multipoles: every one of the 3444 dipoles in the region (±q, q from
    0.25 to 2 M) fails the far beams' 0.05° (340 of them bring the near beam home); of 1644
    compact quadrupoles (triples +q, −2q, +q and rectangles), 1181 keep the far beams
    straight and 71 also bring the near beam home. The region is centred on both far
    flights: a first version with it near their launch points let the search hide a net
    charge of −0.25 M (a charge's kick counts with the part of its field that the flight
    sees, so a negative charge near the start was balanced by positive ones further along;
    in mid-flight that part is stationary). Measured (analyze --fewest, 32 runs): random
    solve rate < 1.5e-3, search 1/32, smoothness 0.69, fewest 3; the description names the
    shape of the answer.
  - Jackson Ch. 16: three orbits: not a staged level (the three electrons radiate at
    the same time), but no search finds a solution with fewer than 3 charges, and 3-charge
    ones in 2 of 24 runs. Every radiation-damping level needs radiation by construction:
    its detector accepts kinetic energy up to 0.26, and an electron that has not radiated
    arrives with at least 0.27 (test `radiation_levels_need_radiation` checks it).
  - Jackson Ch. 16: spectroscopy: the reference drives each atom with its own antenna at
    its own natural frequency (0.125 and 0.1925). In Jackson §16.8: resonance, with ω₀
    removed from the frequency list, no search finds a solution even with two antennas:
    the antennas stand far enough away for a nearly uniform drive, so the electron is a
    linear oscillator. (Close to the atom, a strong antenna at 2ω₀ also worked, by
    parametric resonance in its steep near field.)
  - Indirect goals (knock-on, recoil at right angles, a kick for a bound charge): the
    player throws the projectile as a free charge (slingshot) and places nothing else, so
    the projectile's field is the only thing that acts on the target. Earlier versions
    steered a level-launched projectile with charges, then with magnets far from the
    target; both still acted on it directly (charges shifted it by 14 cells; magnets
    turned a kicked target by several cells and precessed the orbiting electron: the
    owner's reviews). Knock-on and recoil detectors take only T ≥ 0.08 and ≥ 0.1: a
    distant slow push gives at most ~0.03–0.05, so a gentle nudge that drifts in does not
    count (it cut the knock-on solutions from 1149 to 40). `check_indirect` still guards
    levels whose goals are level free particles with elements placed by the player.
  - Bound electrons start on an orbit (radius 2, v = ω₀r), as an atom's electron does, not
    at rest at the atom's centre (the owner's review: at rest, the approaching ion set
    it moving before the kick, and the solutions used that).
  - Radiation goals (arc 11, PHYSICS.md §3.4): a world at c = 2 with a charge of 1/40 at
    γ = 3, radiation reaction included. With a charge of 0.1 the tight bends radiated half
    the particle's energy and max |F_RR| / max |F_L| reached 0.19; at 1/40 (magnets ×4, the
    same trajectories) it is 0.018–0.022 and they radiate 6–8 %. Measured (analyze
    --fewest, 16 runs): forward beaming random solve rate 4.4e-2, fewest 1; critical
    frequency (first version) 2.5e-2, fewest 1; quiet turn 3.0e-3, fewest 2; undulator 2.0e-3, fewest 3.
  - Forward beaming: 106 single-magnet placements pass the goal; flying straight
    radiates nothing. Its reference is a moderate bend (160 at (22, 4)) that fails the
    critical frequency, so the pair shows the difference.
  - The critical frequency: the first version mirrored forward beaming and asked 6e-5 per
    steradian in the band, which most of forward beaming's bends passed (71 of 106): it
    taught nothing new (the owner's review). Now it has forward beaming's arena and
    receiver and asks 3e-3 in the band: 26 single-magnet solutions, the tight bends.
  - Quiet turn: besides gentle turns, one magnet can loop the particle clockwise
    (0° → −90° → 180° → 90°), so its velocity never points at the receiver at 45°: the
    reference does that (6.1e-4 per steradian). Radiation goes where the particle heads
    while it is accelerated.
  - Thomson scattering: the wave's amplitude is 24 (a₀ = qE₀/(mωc) = 0.3, nearly linear;
    the second harmonic shows at 50–70), the magnets only 40 and 80, the minimum 3e-3 per
    steradian: 22 single-magnet solutions, none with the wave off (the bends' own flash).
    A first try (amplitude 8, minimum 3e-4) let 17 bends pass without the wave. Measured
    on the reference: the line peaks at 26–32, as the Doppler formula gives for headings
    near 30°.
  - Braking radiation: flying straight into the target delivers nothing (rejected); 80
    single-magnet solutions; without the stop's radiation (the control) no one- or
    two-magnet placement reaches even 3e-3 (the bends offered are gentle, their critical
    frequency far below the band). The reference delivers 6.06e-3 per steradian against
    an analytic best of 6.3e-3 for the ideal heading.
  - Undulator: the first version (band 31–41, any exit) was solved by two magnets and by
    6 % of random placements: a pair that steered the particle close past a magnet near
    the end made a hard, broadband flash (4.7e-4 per steradian in the band against 2e-4
    for the whole row). Now the particle must arrive on the axis heading straight
    (0° ± 3°), as an undulator must not steer the beam, and the band is 34–38: the best
    of 300 random pairs or triples then delivers 4e-5 and 6.7e-5 (minimum 1e-4). The
    search still finds three-magnet solutions: kicks far apart whose flashes interfere
    (fringes 2π/delay ≈ 8 apart in ω), an undulator of few periods.
  - Rings of charges (Pr. 14.23; PHYSICS.md §3.4, the system measure): c = 4, a charge
    of 1/40 circling at 0.4c (radius 2, ω₀ = 0.8) in the frame drifting with E×B at 0.2,
    about 11 turns before it drifts into its detector; the player's free charges circle
    with it. Measured (the charge alone: 4.84e-5 per steradian in 0.6–1.0, 8.94e-5 in
    0.6–2.6): a second charge in step one cell below, 3.96 times; an antiphase partner on
    the far side at the exact (relativistic) mirror speed 1.738, 0.0035 times, at the
    Galilean 1.8, 0.107 (the drift frame's speeds then differ: its circle turns 0.8 %
    slower and drifts out of antiphase over the flight), at 1.7, 0.045; four evenly
    spaced, 0.0036 times in 0.6–2.6 (Galilean velocities: 0.083), only 4ω₀ left at 15.1
    times, while the pair leaves 1.24 times there (its 2ω₀ line, 3.85 times). The goals
    (in step ≥ 1.4e-4, quiet ring ≤ 7e-6, ring of four ≤ 1.3e-5) pass the Galilean
    placements; the charge alone fails each.
  - Isotope separator: collimation (stage 1) takes one charge; the magnetic separation 4
    elements. A first version with the collectors on opposite sides of the axis had no
    solution: a magnetic field bends both isotopes the same way, the lighter more.
- The tube levels 99–104 (2026-10-12, after the tube's envelope, the emitters' finer mesh,
  the weights from the emitting face and the particles softened at an eighth of a
  segment, PHYSICS.md §2.11) with 300 samples and 8 search runs (a tube evaluation takes
  0.7–3.5 s). Until then `analyze` and
  `generator window` counted every tube placement as verified (the shots' verdicts, of
  which a tube level has none): their search figures and "middle verified" were
  optimistic, and the search ignored the potential goals; both are fixed and the rows
  re-measured. Through the cloud (at 60 V since the owner's review): the search succeeds
  in 5 of 8 runs, random placements rarely (below 1 in 100). Of 90 single plates with the
  anode at 60, 24 bring the ion in; the plates below the gap that do also draw electrons,
  and only the milder ones keep the current in range: the plate trades against the
  anode's voltage. Electron gun: the plates'
  window is +2 to +8 together, 4–6 with one of them held at the reference (the lens is
  symmetric). Ion diode (2026-10-13; after the owner's review a preview takes 2.9 s, a
  verification 12.4 s): the supply solves on 27.3–32.0 at B = 1.9 (2.9 % of ±80, its middle
  verified; before the envelope and the wider slider 27.0–32.3 of ±40), the solenoid
  on 1.65–2.50 at 30 V and on the mirror interval −2.50 to −1.68 (the insulation does not
  depend on the field's sign); the solenoid interval's middle, 2.08, is uncertain, not
  verified: the ion current flattens with B (6.0 at 2.2, 5.05 at 3), so the middle lies
  within the verdict's bar of the goal's lower bound. Magnetron (2026-10-13; a preview
  takes 2.1 s, a verification 8.1 s, after its redesign around the antenna on a larger
  board): the table counts the sliders' candidate values (the band's three points on
  each), so its random rate flatters it; on the continuous sliders the supply solves on
  42.0–47.3 at B = 6 (3.3 % of ±80) and the solenoid on 5.83–6.33 at 45 V (2.5 % of ±10),
  plus the mirrors and a near-cutoff point at 4.58, the main windows' middles verified. Rectifier: two supplies, the full wave on each 24–36 with the other at 32
  (preview), random choices of the five candidate values solve 8 %. Sideways, Vacuum
  diode and Triode are one-slider introductions (the search's five candidate values;
  Sideways solves on 0.9–1.3 of ±4 in the preview, Vacuum diode verifies on 23.5–30.5 of
  ±70, Triode on about −0.9 to −0.4 of ±2).
- Vacuum diode (2026-10-07; re-measured 2026-10-12 in the envelope, its goal widened to
  0.9–1.5 for the verdict's first-order error): the search tries the slider's candidate
  voltages only (±70, 13, 26, 52), so the table counts five configurations; on the
  continuous slider the verified interval is 23.5–30.5, 5.0 % of the range (`generator
  check` at its edges, uncertain at 23 and 31; with the emitters' finer mesh, whose
  smaller error bars widened it from 23.5–28). An introduction: the player learns the
  V^{3/2} law by trying.
- The metals' levels 109–113 (2026-10-14 to 19) with 300 samples and 8 search runs: an evaluation takes ~30 s (a strip's thousands of carriers), minutes far from a solution (Positive carriers at |B| ≥ 3: 4–8 min; a supply reversed or at its slider's end alike). Positive carriers (2026-10-19): the table counts the sliders' candidate values; `generator window` over the half of each slider that can solve (the other reverses a current or a Hall voltage; every 2 V and 0.25, landing on the slider rule's points): the upper supply solves on 28–34 of ±40, the lower on −16 to −8, the solenoid on −2.5 to −1.5 of ±7 (at ±5 the half's middle −2.5 solved). The search needs 128 evaluations against the Hall probe's 19: the two Hall voltages each depend on two sliders, and the lower one's sign on its supply's.
- `analyze --fewest` adds the fewest elements that solve a level: the reference solution's
  count, unless the search finds a verified solution with fewer (exhaustive for one
  element, annealing with 64 restarts of 600 iterations for 2 up to the reference count
  minus one). A failed search is not a proof, so the number is an upper bound.
  Measured 2026-09-28, before the restructure (old numbering): 2 elements for Reflectron,
  Einzel lens, Collimator, Tune the lens, Real Einzel lens, Calutron, Build a Wien filter,
  Mains hum, Earth's field and Collimated beam; 3 for Beam preparation; 6 for Sorting
  station; 1 for every other level. Three references used more elements than needed: the
  Wien filter (reference 4) and the velocity selector (reference 4) were solved with one
  charge, and Build a Wien filter (reference 6) with two elements. The Wien filters now
  require a straight exit (±3°; at ±4° one placement of a single charge still passed) and
  need 2 and 3 elements. The velocity selector became the Wien filter for a beam (three
  speeds, straight exit within ±5°): no one-charge solution, two-charge ones in 2 of 32
  search runs; the Wien filter's own design fails there (space charge). The Reference column of the table above gives the current
  counts.
- The search is a rough stand-in for a player. A player who understands the physics (for
  example the Wien condition v = E/B) should do much better than the search on the later
  levels, and a player who doesn't should do much worse.
