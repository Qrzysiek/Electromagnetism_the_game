# Jackson series: catalogue of implementable effects

Source: J. D. Jackson, *Classical Electrodynamics*, 3rd ed. (Wiley, 1999). Every section and
problem number below was checked against the book's table of contents and problem sets (the
book is kept locally in `references/`, which is gitignored and never committed). The solutions
manual will be consulted for reference outcomes when a level or test is built.

Level names follow the pattern "Jackson §12.3: E×B drift" or "Jackson Pr. 2.4: The golden
ratio", and the description cites the section or problem (planned level field `reference`).

## How to read this list

Each item is tagged with what it needs:

- **NOW**: the current engine can do it (point charges, conducting spheres, electrode plates
  and supplies, dipole magnets and coils with moments along z, antennas, plane waves, moments,
  beams, gates, Liénard–Wiechert fields, Landau–Lifshitz radiation reaction).
- **F:name**: needs a new feature from the table in the next section.
- **3D**: needs the full 3D simulation (listed separately at the end, for when we get there).
- **OUT**: out of reach for a particle simulation (media, wave optics, quantum).

Also:

- **Level** means a puzzle.
- **Test** means a validation test against the book's analytic result (PHYSICS.md). Many items are both.

The 2D game is a slice z = 0 of a 3D world. The slice is exact when the sources are mirror
symmetric under z → −z: E then lies in the plane and B is along z. So:

- coaxial coil pairs placed symmetrically about the plane are consistent with the slice;
- a coil whose axis lies in the plane is not;
- neither is any force with an out-of-plane component.

## New features, ranked by what they unlock

| Feature | What it is | Unlocks |
|---|---|---|
| **F:bound** (bound oscillator; *done* as the charge cloud) | A charge on an isotropic spring `−k x` (mass m, charge e), optionally with damping. Driven by fields and by LL radiation reaction. | §16.7 line breadth; §16.8 scattering by an oscillator; Pr. 13.2/13.3 energy transfer to a bound charge; Pr. 16.1 exponential decay; Pr. 16.12 collision broadening; §4.6 polarizability model; §7.5 dispersion (Lorentz model); Thomson → Rayleigh crossover (Pr. 16.13 sum rule) |
| **F:prescribed** (prescribed-motion sources) | Charges moving on given paths (circle, SHM, figure-8), as radiation sources. | Pr. 9.2 rotating charges (quadrupole at 2ω); Pr. 14.4/14.12/14.14 harmonics of SHM and circular motion; Pr. 14.23/14.24 N charges on a ring (steady current doesn't radiate); Pr. 14.19/14.20 rotating/flipping moment |
| **F:coils+** (more coil kinds) | Helmholtz pairs and solenoids with axis along z, symmetric about the plane (exact in the slice); cylindrical permanent magnets ⊥ plane. | Pr. 5.3/5.5/5.7 (solenoid, Helmholtz field uniformity), Pr. 5.19 (cylindrical magnet), uniform-B regions for §12.2–12.4 without an infinite field |
| **F:ramp** (time-dependent coils, induction) | Coil currents ramped in time, with the induced `E = −∂A/∂t` (quasi-static, §5.15). | §5.15 Faraday's law; betatron (§5.15 + §12.5); Pr. 6.24 (E outside a changing solenoid); §12.5 adiabatic invariance of flux (slowly rising B) |
| **F:cavity** (analytic cavity modes) | Given RF mode fields in a region: pillbox TM₀₁₀ with axis in the plane. | §8.7 resonant cavities: an RF accelerating gap; bunching (klystron idea); linac stages with gates |
| **F:spin** (BMT spin) | Spin vector as an internal variable, precessing by Thomas–BMT (§11.11). The Stern–Gerlach force is neglected, with a bound. A detector condition on spin. | §11.8 Thomas precession; §11.11 BMT; Pr. 12.11–12.12 muon g−2 (spin tune vs orbit) |
| **F:decay** (unstable particles) | A particle with a proper lifetime decays into given products (fixed seed); time dilation from proper time. | §11.2 (muon lifetime experiments); Pr. 11.19–11.21 (two-body decays; Pr. 11.20 Λ → pπ); Pr. 15.13/15.16 (radiation at the sudden creation of a charge) |
| **F:polar** (polarizable / polar neutral particles) | Induced dipole α: force `α∇(E²)/2`; or a permanent electric dipole. | §4.1–4.2 energy of a dipole in a field; Pr. 4.5 force on a dipole; Stark deflector / electrostatic trap for neutral molecules (the electric analogue of Stern–Gerlach); induced dipoles are simple (force α∇(E²)/2); permanent ones need orientation dynamics (torque p×E, moment of inertia) |
| **F:dielectric** (dielectric spheres) | Linear dielectric spheres by the method of fundamental solutions, like the conductors. | §4.4 dielectric sphere in a field; Pr. 4.9 point charge near a dielectric sphere; §4.7 energy (force on a dielectric) |
| **F:sphere+electrode** (spheres and electrodes together) | Conducting spheres in the field of plates and supplies. | §2.5 sphere in a uniform field (between two plates); Pr. 2.10 hemispherical boss on a plane; §2.7 hemispheres at different potentials |
| **F:screened** (world options) | Debye-screened Coulomb `e^{−r/λ}/r`; Proca (massive photon) Yukawa static fields. | Pr. 13.5 screened Coulomb scattering; §12.8 photon mass effects (orbits precess; a "what if" level); §12.9 London screening of B (superconductor as a field option) |
| **F:finite-c beams** (retarded interaction) | Darwin Lagrangian (§12.6), later full retardation. | §12.6 / Pr. 12.13 two-body Darwin dynamics; space charge at relativistic speed (the magnetic attraction partly cancels the electric repulsion, factor 1/γ²) |
| **Radiation goals** (*done*, PHYSICS.md §3.4) | A far receiver over an arc of in-plane directions: energy per steradian of the flight, in all frequencies (Liénard) or in a band (Jackson's radiation integral 14.65), in a window [min, max]. With other particles flying (free particles, free charges) it sees them all: their far fields add (*done*). | §14.3 beaming, §14.6 critical frequency, §14.7 undulators, a quiet turn, §15.2 braking radiation (a target's sudden stop, Pr. 14.10), Pr. 14.23 ring of charges; next: Pr. 14.6 collision radiation with both partners moving (two identical charges: no dipole radiation), Pr. 9.2 rotating charges |
| **Dynamic particles** (*done*) | Level free particles (targets, with optional detectors) and player free charges with a launch velocity (drag the arrow; rapidity for relativistic levels); rigid-sphere collisions. | Pr. 13.1 knock-on and Pr. 13.2 as proper two-body levels; two-body problems with reduced mass (§12.6 / Pr. 12.13 Darwin); Rutherford scattering with a recoiling target; Pr. 14.6 collision radiation with both partners moving |
| **F:plasma** (neutralizing background) | A beam or cloud with a uniform fixed neutralizing background. | Pr. 7.12 plasma oscillations; §7.5 plasma frequency; ion channels |
| **F:antennas+** (more antenna types) | Magnetic dipole (AC loop in plane, moment along z) and centre-fed linear antenna in the plane. | §9.3 magnetic dipole fields; Pr. 9.14; §9.4 centre-fed antenna, Pr. 9.16; §9.12 |
| **F:image-force** (electrode image force) | The force of a particle's own induced charges on electrodes (today a bound, `IMAGE_FORCE_LIMIT`). | Pr. 1.13, Pr. 3.19–3.20 (induced charge between grounded plates, Green reciprocity); vacuum tubes (SPEC "inside components") |

## Built levels

Arc 6, "Jackson: charges in fields and collisions" (Ch. 12–13), `scripts/levels.py`:

- Jackson §12.3: E×B drift (introduction). Uniform crossed fields as a stray field; both
  signs drift with E×B/B² and their guiding centres follow the equipotentials.
- Jackson Pr. 12.9: Van Allen equator (introduction). A dipole Earth; protons and electrons
  drift around it in opposite directions. The drift period measured in the level (about
  95 time units) is longer than the problem's small-gyroradius formula gives (about 70):
  here a/R ≈ 0.4.
- Jackson §12.4: gradient drift (intermediate). A row of magnets makes a gradient; the two
  signs drift in opposite directions.
- Jackson Pr. 12.5: E×B runaway (intermediate). |E| > c|B|: no drift frame; magnets make the
  field strong enough locally. The particle's charge is 1e-8 with fields 100 times the usual
  (the trajectory depends only on qE and qB): with q = 1e-6 it radiated 1.0e-8 of its
  launch energy (it starts from rest, so that energy is tiny), above the 1e-10 that may be
  neglected.
- Jackson Pr. 13.1: knock-on. The player throws a heavy ion (a free charge) past a light
  particle at rest, which must reach its detector with enough energy to count as a kick;
  without the interaction the reference fails.
- Jackson §12.1: Störmer's forbidden region (intermediate). The conserved canonical
  angular momentum keeps a particle aimed at a dipole Earth beyond r = 7.2 cells; placed
  charges break the axial symmetry.
- Jackson Ch. 12: magnetosphere (master). E×B steering through two gates, then gradient
  drift around a dipole Earth.

Arc 7, "Jackson: conductors" (Ch. 2–3):

- Jackson §2.2: its own image (introduction). A particle with charge 1 (Newtonian) passes a
  grounded sphere and is pulled by its image.
- Jackson Pr. 2.6: two spheres (introduction). A charged and a neutral floating sphere.
- Jackson §3.13: field through a hole (intermediate). In the slice the hole is a slot in a
  grounded wall of two plates.
- Jackson Pr. 2.4: golden-ratio capture (intermediate). Test K5 checks the engine against
  the problem's answers.
- Jackson Ch. 2: sphere slalom (master). Three spheres with the particle's charge, two
  gates.
- Not built: §2.5 (sphere in a uniform field) needs metal spheres in a stray field, which
  the validator rejects today (the sphere solver does not include external fields);
  Pr. 2.20 (quadrupole lens) focuses a beam moving across the plane, a 3D effect.

Arc 8, "Jackson: radiation damping" (Ch. 16):

- Jackson Pr. 16.2: the classical atom (introduction). Test R4 checks the r³ law (0.64 %
  at v/c = 0.056). The level makes radiation necessary with an energy window at the
  detector: without radiation the electron arrives too fast.
- Jackson Pr. 16.3: orbits circularize (intermediate). A circular and an elliptic orbit
  of the same energy.
- Jackson Ch. 16: three orbits (master).

Arc 9, "Jackson: bound charges" (§16.7–16.8, Pr. 13.2), built on the new charge cloud
(F:bound, realised as Thomson's atom: an electron inside a sphere of uniform charge is
bound harmonically, PHYSICS.md §2.1):

- Jackson §16.7: a bound charge (introduction).
- Jackson §16.8: resonance (introduction). Without the resonant frequency there is no
  solution.
- Jackson Pr. 13.2: a kick for a bound charge (intermediate). Unit charges, Newtonian;
  without the interaction the reference fails.
- Jackson Ch. 16: spectroscopy (master). Two atoms, two resonances.

Arc 10, "Jackson: radiation" (Ch. 14), on the new radiation goals (PHYSICS.md §3.4):

- Jackson §14.3: forward beaming (introduction). The radiation follows the velocity
  within 1/γ: bend the particle while it heads for the receiver.
- Jackson §14.6: the critical frequency (introduction). A band goal at high frequency: a
  tight bend's short flash.
- Jackson §14.2: a quiet turn (intermediate). A maximum: gentle turns radiate less into
  a direction the velocity sweeps; or turn the other way round.
- Jackson §14.8: Thomson scattering (intermediate). A plane wave shakes the particle; the
  scattered line is Doppler-shifted by the particle's heading (≈ 4γ²ω head on: inverse
  Compton); the band picks the heading.
- Jackson §15.2: braking radiation (intermediate). The goal's target stops the particle
  abruptly (the new `abrupt_stop` option, test S4); the flash is zero straight ahead, so
  the particle must strike the target off the receiver's axis.
- Jackson Pr. 14.23: in step (introduction). A charge circles in crossed fields (ω₀ = 0.8)
  while its circle drifts into the detector; a second charge circling in step doubles the
  field at the receiver: four times the energy (the first level where the receiver sees
  several particles, their far fields added).
- Jackson Pr. 14.23: a quiet ring (intermediate). A partner on the far side of the circle,
  in antiphase, silences the fundamental (0.0035 of one charge at the exact relativistic
  mirror speed, 0.107 at the Galilean one; the goal allows 0.14).
- Jackson Pr. 14.23: a ring of four (intermediate). The band holds ω₀, 2ω₀ and 3ω₀: an
  opposite pair leaves 2ω₀ (3.85 times one charge), three at 120° leave 3ω₀ (9 times);
  four evenly spaced leave only 4ω₀, outside the band (15.1 times one charge there).
- Jackson §14.7: undulator (master). The on-axis line 2γ²ω_u/(1 + K²/2) in a band.
- Tried and set aside first: braking radiation by stopping a γ = 3 particle with charges
  (it needs a large potential hill and the particle slips round the charges). Built
  instead with a target that stops it abruptly.

Validation tests from the book (PHYSICS.md): K5 (Pr. 2.4, golden ratio), M7 (Pr. 12.9b,
equatorial drift), M11 (§12.5, adiabatic invariance), S1 (Pr. 14.15, harmonics of circular
motion), S2 (Parseval, 14.60–14.65), S5 (Pr. 14.23, the form factor of charges on a ring), B17 (Pr. 13.1, knock-on energy transfer), R4 (Pr. 16.2, classical atom
collapse), R5 (Pr. 14.5b, radiation in a head-on collision), C2 (Pr. 16.1, the
radiating oscillator decays at Γ = ω₀²τ).

## Catalogue by chapter

### Ch. 1: Introduction to electrostatics

- **§1.1–1.5 Coulomb, Gauss, potential.** NOW. Already the basis of the early chapters.
- **§1.11 energy and capacitance.**
  - Test: capacitance of a cube, Pr. 1.20b. The book's value is 0.655(4) (variational); the modern value is C = 0.66068 × 4πε₀a. NOW (BEM).
  - Test: the square plate C = 0.3667874 (see Pr. 3.21).
- **§1.13 relaxation method.** Background for our BEM. Test only (compare with Pr. 1.22–1.24 examples).
- **Pr. 1.1b shielding.** No field inside a closed conductor. Level: a Faraday cage of plates around a detector region. NOW (plates), or F:sphere+electrode.
- **Pr. 1.13, induced charge between plates (Green reciprocity).** A charge between grounded plates. F:image-force for a level; a test is NOW with BEM.

### Ch. 2: Boundary-value problems I

- **§2.1–2.4 images.** Charge near a grounded sphere, an insulated sphere, and a sphere at fixed potential. NOW (images; test exists).
  - Levels: "Jackson §2.2: the image charge" (the particle is attracted to a neutral grounded sphere).
- **Pr. 2.4, the golden ratio.** A charge q near an insulated sphere carrying the *same* charge is attracted when it comes close enough: d/R − 1 < 0.618…
  - Level: a like-signed particle is captured by the sphere only inside the golden-ratio zone.
  - Test: the force changes sign at exactly this distance.
  - NOW (floating spheres).
- **§2.5 sphere in a uniform field.** A dipole image; the field enhancement is 3E₀ at the poles. F:sphere+electrode (uniform field from a plate pair); a test with an external uniform field is NOW.
- **Pr. 2.6 two spheres.** Image series and capacitance coefficients. Test exists (K2b); a level "Jackson Pr. 2.6: two spheres" is NOW.
- **§2.7 hemispheres at different potentials.** The split is ⊥ to the plane, so it is consistent with the slice. F:sphere+electrode (split spheres).
- **Pr. 2.10, hemispherical boss on a plane.** Field enhancement 3× at the top. F:sphere+electrode.
- **§2.10–2.11 2D problems, corners.** The field near a wedge/corner scales as ρ^{π/β−1}. Test with BEM plates forming a wedge; a level about field concentration at a sharp edge (§3.4 too). NOW (plates).
- **Pr. 2.20 electric quadrupole lens.** Four line charges give a quadrupole focusing field. As a slice: four point charges or four plates. The level "Jackson Pr. 2.20: quadrupole" focuses a beam in the plane. NOW; alternating-gradient focusing in both planes is 3D.

### Ch. 3: Boundary-value problems II

- **§3.4 conical hole or sharp point.** Field enhancement at a tip (lightning rod). The test is a sharp plate edge; the level is an electron emitted from a tip. NOW / F:image-force.
- **Pr. 3.3, Pr. 3.21: charged disc capacitance.**
  - Disc: C = 8ε₀a.
  - Square plate: C = 0.3667874 × 4πε₀ × side.
  - BEM tests. NOW.
- **§3.13 / Pr. 3.12: field through a circular hole in a conducting plane.** A slit in a plate lets the field leak through. Level: steer a particle through a field leaking through a gap in a grounded wall. NOW (plates).
- **Pr. 3.19–3.20 induced charge (reciprocity).** F:image-force.

### Ch. 4: Multipoles, dielectrics

- **§4.1–4.2 multipoles.**
  - Levels that build a dipole and a quadrupole from charges and fly through their fields. NOW.
  - The energy of a dipole in a field needs polar particles: F:polar.
- **Pr. 4.5 force and torque on a dipole.** Deflection of polar/polarizable neutral molecules (Stark deflector). F:polar.
- **§4.4 dielectric sphere; Pr. 4.9 charge near a dielectric sphere.** F:dielectric.
- **§4.5–4.6 molecular polarizability; the oscillator model.** F:bound (charge on a spring polarizes in a field).
- **§4.7 energy in dielectrics** (a dielectric is pulled into a capacitor). F:dielectric; the materials stage.

### Ch. 5: Magnetostatics, Faraday's law

- **§5.2 Biot–Savart; §5.5 circular loop.** NOW (loops, polygons).
- **§5.6 magnetic moment; §5.7 force on a moment.** NOW (Stern–Gerlach and Stern–Gerlach beam).
- **Pr. 5.3/5.5/5.7 solenoid and Helmholtz pair field uniformity.** A level asks for a uniform-field region from a coil pair. F:coils+.
- **Pr. 5.19 cylindrical magnet.** F:coils+.
- **§5.10–5.12 magnetized sphere, permanent magnets, magnetic shielding.** The materials stage (F:dielectric analogue for μ).
- **§5.15 Faraday's law.** Betatron (induction accelerator): the "2:1 rule", with the orbit flux twice the flux at the orbit field. F:ramp.
- **§5.17 inductances.** The circuits stage (SPEC).
- **§5.18 eddy currents.** OUT (conductors as media).

### Ch. 6: Maxwell equations, conservation laws

- **§6.4–6.5 retarded solutions (Jefimenko).** NOW (Liénard–Wiechert fields). Test: Pr. 6.2 Heaviside–Feynman field formula against LW.
- **§6.7 Poynting theorem.** A view: energy flow around a charge in a field. NOW (visual).
- **Pr. 6.15 Hall effect.** The Drude stage (SPEC "inside components").
- **Pr. 6.24 changing solenoid.** E outside, where B = 0 (induction without B). F:ramp.
- **§6.11 monopoles.** A charge–monopole system has out-of-plane forces. 3D.

### Ch. 7–8: Waves, waveguides, cavities

- **§7.5 plasma frequency; Pr. 7.12 plasma oscillations.** F:plasma.
- **§7.5 Lorentz dispersion model.** F:bound (driven oscillators by a plane wave).
- **§7.7 MHD waves.** OUT.
- **§8.7 resonant cavities (TM₀₁₀ pillbox with axis in the plane).** An RF accelerator gap: timing a particle to the RF phase. F:cavity.
- **§8.1–8.6 waveguides, §8.10–8.11 fibres.** OUT (they need FDTD or mode solvers, not particles).

### Ch. 9: Radiating systems

- **§9.2 electric dipole.** NOW (antennas).
- **§9.3 magnetic dipole and quadrupole; Pr. 9.14.** F:antennas+.
- **§9.4 centre-fed linear antenna; Pr. 9.16; §9.12.** F:antennas+.
- **Pr. 9.2 rotating charges:**
  - a rotating dipole radiates at ω;
  - a symmetric pair radiates at 2ω (quadrupole).
  - F:prescribed.
- **Pr. 9.9 + §16.2 the classical atom collapses.**
  - Level: an electron on a Coulomb orbit spirals in (LL radiation reaction, exaggerated constants in a "what if" level).
  - Test: fall time t = r₀³/(4 r_e² c) (Pr. 16.2a).
  - NOW.
- **Pr. 9.8 radiated angular momentum.** A diagnostic for circular motion. NOW (test).

### Ch. 10: Scattering and diffraction

- **OUT**: Rayleigh scattering, Mie scattering, diffraction and Babinet need wave solvers.
- Thomson scattering by free charges (§14.8) is covered in Ch. 14.

### Ch. 11: Special relativity

- **§11.2 experiments; Pr. 11.19–11.21: decays in flight.** Lifetime dilation; Λ → pπ kinematics (Pr. 11.20) (reconstruct the parent from the two tracks in B). F:decay.
- **§11.3–11.5 kinematics.** NOW (every relativistic level).
- **Pr. 11.13 moving line charge.** E and B transform. With LW fields: NOW.
- **Pr. 11.17–11.18 fields of a fast charge.** The flattened "pancake" field (γ). A view, NOW (particle field map).
- **§11.8 Thomas precession; §11.11 BMT.** F:spin.
- **§11.10 field transformations.** A level: in the frame moving with v = E×B/B² the electric field vanishes (the §12.3 drift explained). NOW (a text plus a drift level).

### Ch. 12: Dynamics of relativistic particles

- **§12.2 uniform B.** Gyration, and the gyroradius from momentum. NOW (the cyclotron and mass spectrometer levels).
- **§12.3 E×B drift; Pr. 12.4–12.5.**
  - Crossed fields with |E| < c|B| give a drift.
  - With |E| > c|B| the particle runs away (hyperbolic).
  - A level with both regimes; a test against the exact solution.
  - NOW (with plate pairs in coils, or F:coils+).
- **§12.4 gradient drift in the plane.** Drift along the contours of |B|. Level plus test. NOW. Curvature drift needs motion along B: 3D.
- **§12.5 adiabatic invariance of flux.** Slowly increasing B shrinks the orbit, keeping p⊥²/B. Test M11 (*done*): μ conserved to 1e-3 while B rises 4×, and the guiding centre drifts inward keeping ρ²B. The Faraday level shows the orbit contracting. A dedicated level waits for coil supplies (the player sets the ramp and compresses the orbit into a detector): with a fixed ramp and placed charges steering by E×B drift, the drift is weak against the energy the charges can give directly, and the guiding centre's inward drift into weaker field (in a loop's plane) gains only ~1.9× in energy for 4× in current.
- **Pr. 12.9 Van Allen belts, equatorial part.** A dipole "Earth" (moment along z): the gradient drift of trapped particles around it in the equatorial plane.
  - Level: guide a particle once round the Earth.
  - NOW. The bounce between mirror points is 3D.
- **Pr. 12.11–12.12 muon g−2.** Spin precession relative to momentum ∝ a = (g−2)/2 in a storage ring. F:spin.
- **§12.6 Darwin Lagrangian; Pr. 12.13.** F:finite-c beams.
- **§12.8 Proca; §12.9 London.** F:screened (world options).
- **§12.1 canonical momentum.** A level: the conserved canonical momentum in a symmetric field predicts where the particle turns (a Störmer-type forbidden region). NOW, plus a text.

### Ch. 13: Collisions and energy loss

- **§13.1 / Pr. 13.1 energy transfer in a Coulomb collision.** ΔE(b) for a heavy particle passing a free light charge.
  - A beam two-body level: knock the light charge into a detector.
  - A test against the impulse approximation.
  - NOW (beams with two species).
- **Pr. 13.2/13.3 energy transfer to a harmonically bound charge.** Adiabatic cut-off for b > v/ω₀: slow collisions transfer nothing. F:bound.
- **§13.2–13.3 energy loss in matter, density effect.** OUT (a statistical medium).
- **§13.4 Cherenkov; §13.7 transition radiation.** OUT (media).
- **§13.5 / Pr. 13.5 screened Coulomb scattering.** F:screened.
- **§13.6 multiple scattering.** A statistical beam through a field of random fixed charges (fixed seed): measure the rms angle against the Gaussian law. NOW (beams + many charges; cost limits).
- **Rutherford scattering (§13.5 context).** NOW (Geiger–Marsden level exists).

### Ch. 14: Radiation by moving charges

- **§14.1 LW fields; Pr. 14.2 near fields of an accelerated charge.** NOW.
- **§14.2 Larmor and Liénard power; Pr. 14.11 invariant form.** NOW (radiation diagnostics, test).
- **§14.3 angular distribution.** *Level built (forward beaming).*
  - Forward beaming (1/γ cone) for acceleration along and across v.
  - Level: place a detector for radiation (the waves view).
  - NOW (in-plane cut of the pattern).
- **Pr. 14.5–14.8 radiation in Coulomb collisions.** For example, the head-on nonrelativistic collision radiates ΔW = 8 z m v₀⁵ / (45 Z c³) (Pr. 14.5b). Test with LL/Larmor integration. NOW.
- **Pr. 14.9 synchrotron energy decay.** NOW (exact test exists).
- **Pr. 14.10 sudden stop.** Bremsstrahlung of a charge stopped by a hard wall (obstacle). NOW (a view: the radiation shell).
- **§14.6 synchrotron spectrum; Pr. 14.15, 14.17–14.18 circular and helical motion.** *Test S1 (Pr. 14.15); level built (the critical frequency).* The in-plane spectrum at harmonics of ω₀. NOW (spectrum diagnostic to add). Helical motion is 3D.
- **Pr. 14.12/14.14/14.22 harmonics of SHM, circular and elliptic orbits.** F:prescribed, or NOW for a Coulomb orbit (the ellipse's harmonics).
- **Pr. 14.21 correspondence principle.** Radiation from a classical hydrogen orbit versus Bohr transition rates. A test with the §16.2 decay. NOW.
- **Pr. 14.23 N charges on a ring.** A beam of N equally spaced charges on a circle: the radiation falls exponentially with N ("a steady current does not radiate"). *Levels built (in step, a quiet ring, a ring of four); test S5 (the form factor).* The player's free charges circle in crossed fields with the shot, whose circle drifts into its detector; the receiver sees them all.
- **Pr. 14.25–14.26 synchrotron polarization, Crab nebula.** Polarization is 3D; energetics as a text.
- **§14.7 undulators and wigglers; Pr. 14.27 second harmonic.** *Level built (undulator).*
  - The trajectory and its in-plane radiation. NOW (alternating magnets).
  - K parameter: undulator for K < 1, wiggler for K ≫ 1; the level asks for a given K.
- **§14.8 Thomson scattering.** *Level built (Doppler-shifted Thomson scattering).* A free charge driven by a plane wave re-radiates. NOW (plane waves + LW). Radiation pressure pushes the charge forward: the LL force. A test against σ_T.

### Ch. 15: Bremsstrahlung, virtual quanta

- **§15.1–15.2 radiation in collisions; Pr. 15.2–15.3 hard-sphere collision.**
  - Particle bouncing off a hard sphere (an obstacle).
  - The radiated spectrum is flat up to ω ~ 1/(collision time).
  - NOW (with a spectrum diagnostic).
- **§15.1 / Pr. 15.8 collisions of like particles: no dipole radiation.** Two charges interacting only with each other have `d̈ = F (q₁/m₁ − q₂/m₂)`: with equal charge-to-mass ratios the dipole moment stays with the centre of mass and only quadrupole radiation is left (Pr. 15.8 with a hard sphere; here the Coulomb collision). NOW, with the system measure (radiation goals see every particle). *Prototype (2026-09-29, not shipped):* c = 4, the shot q = m = 1 at β = 0.3 deflected ~40° into a detector by the player's free charge thrown at it: with an identical partner the pattern is inversion-symmetric, as a quadrupole's (5.2e-5 per steradian at 0° and at 180°, 2.3e-5 at 45° and 225°); with a partner of charge 2 (mass 1) dipole and quadrupole interfere (8e-6 at 0°, 2.1e-4 at 180°). Level idea tried and dropped (2026-09-29): free charges of mass 2 and charge 1 or 2 thrown at the shot, a quiet receiver. The partner of the shot's charge-to-mass ratio (charge 2) has no dipole radiation, but its collision is twice as strong, and its quadrupole radiation (3.5e-5 to 7e-5 per steradian behind the shot, 2.2e-5 to 6e-5 elsewhere) exceeds the unequal-ratio partner's dipole-plus-quadrupole radiation (charge 1: 8e-7 to 4.4e-6 behind, where the two interfere destructively; 1.7e-5 to 1e-4 ahead): at β = 0.3 "no dipole" does not mean "quieter", and a quiet goal would reward the wrong partner. The effect needs v ≪ c and equally strong collisions, which a player's free choice of partner does not give.
- **Pr. 15.9 screened; Pr. 15.10 hyperbolic Coulomb orbit.** The dipole-approximation spectrum with the Hankel-function formula. A test against the numerical Fourier transform of the acceleration. NOW (Coulomb); F:screened.
- **§15.4–15.5 virtual quanta (Weizsäcker–Williams).** The field of a fast charge as a pulse of radiation. A view: the particle field map of a γ ≫ 1 charge. NOW.
- **§15.6–15.7, Pr. 15.4, 15.13–15.16 radiation at sudden creation/disappearance of charges** (beta decay, π → μ, K decays). F:decay (instant creation gives the radiation shell).
- **Pr. 15.7 fission inner bremsstrahlung, Pr. 15.12 energy loss.** OUT (nuclear / statistical); a text only.

### Ch. 16: Radiation damping

- **§16.2 radiation reaction from energy conservation.** NOW (LL).
- **Pr. 16.2 circular orbit collapse; Pr. 16.3 elliptic orbits circularize.** Eccentricity ∝ (L/L₀)^{3/2}. Tests plus a "classical atom" level. NOW.
- **§16.3–16.6 Abraham–Lorentz, runaways, preacceleration, classical electron models.**
  - Text only (LL avoids them).
  - Pr. 16.10–16.11 (gap acceleration with damping) as a test of LL against the integro-differential solution: NOW.
- **Pr. 16.1, §16.7 radiating oscillator.** Level breadth: energy decays as e^{−Γt}, and the line width is Γ. F:bound.
- **§16.8 scattering and absorption by an oscillator.** The resonance cross section, Thomson at high ω, Rayleigh (ω⁴) at low ω. F:bound + plane waves.
- **Pr. 16.12 collision broadening.** F:bound + random interruptions.
- **Pr. 16.13 dipole sum rule.** A test for F:bound.

## Cross-topic highlights (most interesting levels)

These combine several chapters; most come from the problems.

1. **Golden-ratio capture** (Pr. 2.4): like charges attract near an insulated sphere. NOW.
2. **Classical atom** (Pr. 9.9, 14.21, 16.2–16.3): the orbit decays and circularizes; compare with Bohr rates. NOW.
3. **Van Allen equator** (Pr. 12.9, §12.4, §5.6): gradient drift around a dipole Earth. NOW.
4. **E×B runaway** (§12.3, §11.10, Pr. 12.4–12.5): drift below c|B|, escape above. NOW.
5. **Ring of charges** (Pr. 14.23): radiation vanishing with N. NOW.
6. **Undulator** (§14.7, Pr. 14.27): tune K for the radiation. NOW.
7. **Knock-on** (Pr. 13.1): a heavy projectile kicks a light charge into a detector. NOW.
8. **Quadrupole lens** (Pr. 2.20) + collimated beam (Collimated beam, Velocity selector). NOW.
9. **Betatron** (§5.15, §12.5): induction acceleration with the 2:1 condition. F:ramp.
10. **g−2 ring** (Pr. 12.11–12.12, §11.11): spin tune. F:spin.
11. **Λ decay reconstruction** (Pr. 11.20): decay in flight in a B field. F:decay.
12. **Stark deflector** (Pr. 4.5): the electric Stern–Gerlach for polar molecules. F:polar.
13. **RF cavity linac** (§8.7): phase stability with gates. F:cavity.
14. **Radiating oscillator** (§16.7–16.8): line width and resonance scattering of a plane wave. F:bound.

## For the 3D simulation (future reference)

These need motion or fields out of the plane, so they wait for the 3D stage.

- **Magnetic mirror and bottle** (§12.4–12.5): the adiabatic invariant p⊥²/B reflects particles; loss cone.
- **Van Allen belts, full** (Pr. 12.9): bounce between mirror points, the three periodic motions.
- **Curvature drift** (§12.4): needs motion along curved B lines.
- **Helical motion along B** (§12.2): pitch angle; helical undulators and radiation from helical orbits (Pr. 14.17–14.18).
- **Magnetic lenses** (the Busch solenoid lens, §5 fields + §12) and alternating-gradient quadrupole focusing in both transverse planes (Pr. 2.20 in 3D).
- **Coils with axes in the plane; arbitrary magnet orientations** (§5.5–5.7): B in the plane gives out-of-plane forces.
- **Monopoles** (§6.11–6.12): charge–monopole angular momentum, cone of motion.
- **Stern–Gerlach with a realistic magnet** (inhomogeneous field in 3D, arbitrary spin orientation; §5.7, §11.11).
- **Synchrotron radiation out of the orbital plane** (§14.6): vertical opening angle 1/γ, polarization (Pr. 14.25–14.26).
- **Radiation patterns in 3D** (§9.2–9.4, §14.3): full lobes, polarization, dipole antennas oriented freely.
- **Thomson scattering angular pattern and polarization** (§14.8).
- **Penning-type traps** (electric quadrupole + axial B): the axial motion.
- **Dielectric and permeable spheres in arbitrary fields; magnetic shielding by a spherical shell** (§4.4, §5.12). These are feasible in 2D slices only for symmetric setups.
- **Multiple scattering in a volume** (§13.6), energy loss geometry.
- **Rutherford differential cross section in 3D** (the 2D slice is exact per orbit, but the solid-angle statistics need 3D sampling).

## Out of reach (all dimensions)

- Media and wave optics: reflection and refraction (§7.3–7.4), dispersion in media as a continuum (§7.5–7.11), waveguides and fibres (Ch. 8 apart from analytic cavity modes), scattering and diffraction (Ch. 10), Cherenkov and transition radiation (§13.4, §13.7), energy loss in matter (§13.2–13.3), eddy currents (§5.18).
- Quantum topics: the Dirac quantization condition (§6.12), nuclear and atomic applications beyond classical estimates (§9.11, Ch. 15 decays except as classical sudden-creation radiation).
