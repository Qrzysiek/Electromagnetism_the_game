# Brau series: catalogue of implementable problems

Source: C. A. Brau, *Modern Problems in Classical Electrodynamics* (Oxford University Press,
2004). The book is kept locally in `references/` (gitignored, never committed). It is a
scan; its OCR text (Tesseract 5.5, local, also gitignored) was used to list the sections
and exercises, and the page images for the equations. Section and exercise numbers below
were checked against the pages (book page = PDF page − 16).

Brau covers much of Jackson's ground in a relativistic, field-theoretic order, adding
modern optics (Chapters 6–9) and accelerator radiation (Chapter 10). Most of what the
engine can do from it is therefore already built from Jackson (`docs/JACKSON.md`); this
list records the overlap and what Brau adds.

## Status (2026-10-05)

- **Built from Brau:** test C5 (§5.2.3, excitation by a fast charged particle, the same
  problem as Jackson Pr. 13.3); test C6 (§4.3.2 and §10.5, a cluster's collective
  radiation damping: it showed that the beam runner's quasi-static interaction loses the
  collective damping of bound charges, PHYSICS.md §3.3); test P5 (§11.1.2, the 4/3
  problem); test P6 (Ex. 10.15, the electron's shadow: the optical theorem through the
  exchange flux, and the axis dark above the resonance, bright below it); test W5
  (Ex. 4.6, Lawson–Woodward, and radiation reaction's violation of it); test S9 (§10.3.2,
  nonlinear Thomson harmonics); test R8 (Ex. 10.14, radiation pressure on a moving
  charge); test C7 (§7.3.1, Bohr's classical energy loss); test C8 (Ex. 3.6, the
  polarization force of an atom); test P7 (Ex. 2.24, magnets through the stress tensor).
- **Proposed (NOW, for the owner to choose):** the ten items in "Proposals" below, all
  validation tests and one level.
- **Needs a feature:** F:emission (space-charge flow, Ex. 3.1; planned with vacuum tubes),
  F:screened (Proca, Ex. 2.25), F:dielectric and F:magnetic (Chapter 6), F:optics
  (Chapters 7–9, transition and Cherenkov radiation), F:spin (§11.3, Ex. 11.9).
- **3D or OUT:** magnetic monopoles (§11.2: their forces on charges moving in the plane
  leave the plane), solenoids and spheroids (Ex. 3.12–3.13, 3.21–3.22), quantum and
  statistical topics (photons, blackbody radiation, thermal plasmas, Raman and nonlinear
  media).

## How to read this list

Tags as in `docs/JACKSON.md`: **NOW** (the engine can do it), **F:name** (a feature in
JACKSON.md's table, or named here), **3D**, **OUT**. *Covered* gives the test or level
that already does it. "Ex." is an exercise, "§" a section.

## Chapters

### Prologue, Ch. 1 Relativistic kinematics

- §0, §1.1–1.4: Coulomb to Maxwell, Lorentz transformations, 4-vectors, the field tensor.
  Covered by the engine's foundations (T-tests: relativistic dynamics; L-tests:
  Liénard–Wiechert fields).
- Ex. 1.5, 1.7 (the twin and the constantly accelerated rocket): hyperbolic motion is a
  charge in a uniform E (covered, T-tests); proper time is not shown by the game.
- Ex. 1.11 (Doppler shift): covered by level 84 (Thomson scattering with Doppler shifts).
- Ex. 1.18 (E ∥ B in some frame), Ex. 1.21 (a moving capacitor), Ex. 1.22 (a moving
  current): field transformations; moving conductors are OUT (electrodes are fixed).

### Ch. 2 Relativistic mechanics and field theory

- §2.1–2.2 (Lagrangian mechanics, canonical momentum): covered (M10: canonical angular
  momentum conserved in a ramped field).
- Ex. 2.8 (relativistic corrections to a near-circular Coulomb orbit; precession):
  covered (T3, T5 against the exact relativistic Coulomb problem).
- **Ex. 2.10 (the aperture lens of an electron gun):** an anode with a hole between
  fields of different strength focuses the beam. In the slice the hole is a slit between
  two plates. NOW (a level, and a test of the focal length; the slit's analogue of
  Davisson and Calbick's f = 4V/ΔE for a round hole).
- Ex. 2.13, 2.14, 2.21 (uniform E, uniform B): covered.
- Ex. 2.20 (gravity as a scalar field: the perihelion advance is −1/6 of general
  relativity's): F:scalar-field, a "what if" world option; low priority.
- **Ex. 2.24 (the attraction of two magnets from the stress tensor):** NOW, a test (the
  magnetostatic counterpart of P4).
- Ex. 2.25 (Proca mass term): F:screened. Ex. 2.17, 2.26 (nonlinear electrodynamics):
  OUT.
- §2.4.2 (the symmetric stress tensor): covered (P4).

### Ch. 3 Time-independent fields

- §3.1–3.2 (electrostatics, conductors, capacitance, images, numerical methods): covered
  (E1–E6 boundary elements, K1–K5 metal spheres, Z1 capacitance and reciprocation).
- Ex. 3.1 (space-charge-limited flow, Child–Langmuir): F:emission, planned with the vacuum
  tubes (SPEC "inside the components").
- Ex. 3.2, 3.16 (Cavendish's photon-mass experiment): F:screened and 3D.
- **Ex. 3.6 and §3.1.3 (multipoles; the induced dipole of an atom):** a charge z passing a
  charge-cloud atom feels its polarization, the potential `−α z²/(2r⁴)`, with the bound
  charge's polarizability `α = q²/(mω₀²)`. NOW, a test (slow passes; the adiabatic limit
  of C5).
- Ex. 3.10 (an isolated sphere and a point charge): covered (K-tests).
- Ex. 3.12–3.13 (spheroids), Ex. 3.21–3.22 (solenoids): 3D shapes.
- §3.3 (Biot–Savart, forces, multipoles): covered (M-tests, Z2–Z3 inductances).
- Ex. 3.28 (a superconducting sphere expels B): F:magnetic (a μ = 0 sphere, with the
  dielectric bodies, #30).

### Ch. 4 Electromagnetic waves

- §4.1.2 (a charge in a plane wave): covered (W2 invariants, R7 radiation pressure, level
  84).
- **Ex. 4.6 (a charge in a plane-wave pulse):** during the pulse the mean energy rises
  with the local intensity; after it the charge is at rest again, displaced (the
  Lawson–Woodward theorem: a plane wave in vacuum cannot accelerate a charge net). NOW, a
  test with a pulse envelope (a test-local field, as in C5), relativistic intensities
  included.
- Ex. 4.8 (the steady figure-of-eight with cusps): covered (W2).
- §4.2.2 (spontaneous emission by a harmonic oscillator): covered (C2).
- §4.2.4 (blackbody radiation, Einstein coefficients): OUT (statistical, quantum).
- **§4.3.2 (plasma oscillations), with §10.5 (coherent emission):** N electrons in a
  charge cloud. Their centre of mass oscillates at exactly the single-charge ω₀ whatever
  their mutual repulsion (Kohn's theorem, exact for c = ∞), so the mode is the cluster's
  dipole plasmon. At finite c each electron feels the others' radiation fields, and the
  mode should decay at N times the single-charge rate (classical superradiance). NOW, a
  test of the beam runner's retarded fields at order 1/c³, which no test reaches yet, and
  a measurement of how the default quasi-static interaction treats it.
- §4.3.1, Ex. 4.18–4.19 (waves in magnetized and thermal plasmas): F:plasma or OUT
  (continuum and statistical).

### Ch. 5 Fourier techniques and virtual quanta

- §5.1, Ex. 5.1–5.9 (Fourier theorems, pulse compression, chirped pulses): mathematics and
  optics; nothing for the engine beyond what the spectrum code already does (S-tests).
- §5.2.1, Ex. 5.13 (virtual quanta, the equivalent radiation field): the spectrum of a
  passing charge's field; the content of C5.
- §5.2.2 (bremsstrahlung): covered (S7, Coulomb bremsstrahlung).
- **§5.2.3 (excitation by a fast charged particle):** built, test C5.
- §5.2.4 (transition radiation): F:optics (the radiation of the charge and its image at a
  conducting surface).

### Ch. 6–9 Materials, dispersive media, nonlinear optics, diffraction

- Ch. 6 (polarization, magnetization, Ex. 6.5 and 6.10 iron-core magnets): F:dielectric,
  F:magnetic (#30).
- Ch. 7 (Kramers–Kronig, Lorentz–Drude, reflection, surface plasmons, Ex. 7.1–7.18):
  F:optics; the single Lorentz oscillator is built (C3, C4).
- **§7.3.1 (energy loss of a fast particle, Bohr's classical formula):** the energy
  given to bound charges summed over impact parameters, `2π ∫ ΔE(b) b db`, has the closed
  form `ξ K₀ K₁ − (β²/2) ξ² (K₁² − K₀²)` per unit length (Jackson §13.2's distant
  collisions). NOW, a test extending C5 (Jackson tags energy loss in matter OUT as a
  statistical medium, but the per-atom integral is deterministic).
- Ch. 8 (nonlinear optics): OUT (nonlinear media).
- Ch. 9 (geometrical optics, Gaussian beams, resonators, diffraction): F:optics. Ex. 9.1
  (the ray as a particle: the refractive index plays the momentum) is the bridge between
  the game's electron optics and the planned light optics.

### Ch. 10 Radiation by relativistic particles

- §10.1 (Liénard–Wiechert fields, multipoles, spectra; Ex. 10.1–10.8, 10.16): covered
  (L-tests, A-tests antennas, S1 circular harmonics, S2–S3).
- **Ex. 10.9 (the radiation of a head-on collision with any central potential):** covered
  for Coulomb (R5); NOW as a variant through a charge cloud (harmonic inside, Coulomb
  outside), a small test.
- Ex. 10.10 (Coulomb scattering), Ex. 10.12 (a hard-sphere collision): covered (S7, S4).
- Ex. 10.11 (the momentum radiated in a circular orbit): covered (R1).
- **Ex. 10.14 (radiation pressure on a charge crossing a wave at relativistic speed):**
  NOW, a test extending R7 (from rest) to relativistic crossings (the Doppler factors of
  inverse Compton recoil).
- **Ex. 10.15 (the electron's shadow):** behind a charge driven by a plane wave the
  scattered wave interferes destructively with the incident one; its flux deficit is the
  extinction (the optical theorem). NOW, a test joining C3 (the bound charge's steady
  state) and the energy-flow tests (P1–P3: the exchange flux through a large sphere).
- **§10.3.2 (nonlinear Thomson scattering):** the harmonics of a charge in a strong wave
  (a₀ ~ 1). Level 84 shows the second harmonic; NOW, a test of the harmonic spectrum
  against the exact Bessel-function result.
- §10.4 (synchrotron and undulator radiation, Ex. 10.6, 10.17): covered (levels 38, 81,
  88; S1, S8).
- §10.5 (coherence, form factors, Ex. 10.18–10.19): covered (S5, level 87); Ex. 10.20
  (coherent transition radiation): F:optics.
- §10.6 (Cherenkov radiation, Ex. 10.22): F:optics. Ex. 10.21 (tachyons): OUT.

### Ch. 11 Fundamental particles in classical electrodynamics

- **§11.1.2, Ex. 11.1 (electromagnetic mass and the 4/3 problem):** the field energy and
  momentum of a moving charged sphere, `U = U₀ γ (1 + β²/3)` and `p = (4/3)(U₀/c²) γ v`,
  do not form a 4-vector without Poincaré's stresses. NOW, a test with P4's machinery
  (the engine's Liénard–Wiechert field integrated outside the contracted sphere).
- §11.1.3 (point particles and radiation reaction): covered (R1–R7, Landau–Lifshitz).
- §11.1.4, Ex. 11.4–11.5 (extended particles, the Sommerfeld–Page equation): F:extended
  (a delay-differential self-force), OUT for now.
- §11.2, Ex. 11.6–11.7 (magnetic monopoles): 3D.
- §11.3, Ex. 11.9 (spin, Thomas precession, the muon's magic γ): F:spin.

## Proposals

Ranked by what they would add (all tests unless marked):

1. **Collective radiation damping of a cluster** (§4.3.2, §10.5): built, test C6. The
   exact retarded interaction gives it within 0.75 %; the quasi-static one loses it.
2. **The electron's shadow** (Ex. 10.15): built, test P6 (within 3.4e-12; the shadow
   on the axis is dark above the resonance and bright below it).
3. **The 4/3 problem** (§11.1.2, Ex. 11.1): built, test P5 (within 1.3e-14).
4. **Lawson–Woodward** (Ex. 4.6): built, test W5 (at rest after the pulse to 1.8e-13;
   with radiation reaction the exact push, 1.7e-13).
5. **Nonlinear Thomson harmonics** (§10.3.2): built, test S9 (≤ 3.8e-9 against the
   closed-form orbit; along the wave only the fundamental radiates).
6. **Relativistic radiation pressure** (Ex. 10.14): built, test R8 (≤ 4.9e-5 at
   β = 0.5, 0.9 and four angles).
7. **Bohr's classical energy loss** (§7.3.1): built, test C7 (7.5e-6).
8. **Polarization force of an atom** (Ex. 3.6): built, test C8 (2.8e-5 against the
   exact linear response; the adiabatic `−α z²/(2r⁴)` impulse approached as ~1.9/ξ²).
9. **Magnets attracting through the stress tensor** (Ex. 2.24): built, test P7 (2.8e-15).
10. **The aperture lens** (Ex. 2.10): a level and a test (electron optics).
