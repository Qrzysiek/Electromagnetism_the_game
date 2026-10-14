//! The physical model of a level in words: what is exact and what is approximated or
//! left out, with where the game measures the size of each approximation. Shown in the
//! game, so that the model's limits are stated where the physics is played (details in
//! PHYSICS.md).

use crate::{Element, ElementKind, Level};

/// One statement about the model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModelNote {
    /// Exact within classical electrodynamics (up to integration accuracy), or an
    /// approximation / omission.
    pub exact: bool,
    pub text: String,
}

fn note(exact: bool, text: &str) -> ModelNote {
    ModelNote {
        exact,
        text: text.into(),
    }
}

impl Level {
    /// The model of this level with a placement, exact parts first.
    pub fn model_notes(&self, player: &[Element]) -> Vec<ModelNote> {
        if let Some(spec) = &self.tube {
            return tube_notes(spec);
        }
        let all = || self.elements.iter().chain(player);
        let has = |k: ElementKind| all().any(|e| e.kind == k);
        let antennas = has(ElementKind::Antenna) || self.limits.max_antennas > 0;
        let magnets = has(ElementKind::Magnet) || self.limits.max_magnets > 0;
        let electrodes = !self.electrodes.is_empty() || self.limits.max_plates > 0;
        let mut out = vec![note(
            true,
            "Classical electrodynamics in the plane z = 0 of a 3D world: charges, magnets \
             and metal are 3D objects, and their full 3D fields are used in the plane.",
        )];
        out.push(if self.physics.c.is_some() {
            note(
                true,
                "Relativistic motion of the particle in the given fields (exact at any speed).",
            )
        } else {
            note(
                true,
                "Newtonian mechanics (c = ∞): the particle's motion is exact in this limit; \
                 there is no radiation and no retardation.",
            )
        });
        if magnets {
            out.push(note(
                true,
                "Magnets: uniformly magnetized spheres, exactly a dipole field outside them.",
            ));
        }
        if !self.coils.is_empty() {
            out.push(note(
                false,
                "Coils: the current flows in a thin filament (exact closed-form fields); the \
                 wire's thickness only makes it an obstacle.",
            ));
        }
        if antennas {
            out.push(note(
                true,
                "Antennas: point dipoles with their exact retarded fields (near, induction and \
                 radiation terms).",
            ));
        }
        if self.shots.iter().any(|s| s.particle.moment != 0.0)
            || self.free_particles.iter().any(|f| f.particle.moment != 0.0)
        {
            out.push(note(
                true,
                "Magnetic moments: fixed perpendicular to the plane (a spin state, up or \
                 down); B is perpendicular to the plane too, so there is no torque and the \
                 force is m grad B_z, exactly (also relativistically, without electric fields).",
            ));
            out.push(note(
                false,
                "Spin is classical here: each particle is in one of the two states measured \
                 along z (an unpolarised beam is a 50/50 mixture), not in a superposition.",
            ));
            out.push(note(
                false,
                "The radiation of the accelerated magnetic moment is neglected: its estimate \
                 m² ∫|da/dt|² dt / c⁷ (a moving moment carries an electric dipole v×m/c²) is \
                 checked to be below 1e-10 of the launch energy in the shipped levels.",
            ));
        }
        if self.coils.iter().any(crate::Coil::is_ramped) {
            // The coils' own retardation: the light time across the arena against the time
            // over which the current changes appreciably, κ / (dκ/dt).
            let size = f64::from(self.grid.nx.max(self.grid.ny));
            let ratio = self
                .coils
                .iter()
                .map(|c| match c {
                    crate::Coil::Circle { kappa, rate, .. }
                    | crate::Coil::Polygon { kappa, rate, .. } => {
                        let peak = kappa.abs().max((kappa + rate * self.physics.t_max).abs());
                        if peak > 0.0 { rate.abs() / peak } else { 0.0 }
                    }
                })
                .fold(0.0, f64::max)
                * size
                / self.c();
            out.push(note(
                false,
                &format!(
                    "Ramped coils: the current rises linearly; its induced field −∂A/∂t is \
                     quasi-static (the coil's retardation and radiation are left out; light \
                     time across the arena × relative ramp rate = {ratio:.1e})."
                ),
            ));
        }
        if self.has_drives() {
            out.push(note(
                true,
                "Circuits: ideal lumped components (Kirchhoff's laws), sources through \
                 resistors, inductors and capacitors into the plates and coils, integrated \
                 to about 1e-9 of their scale; a plate's capacitances come from the boundary \
                 elements, a coil's inductance from the flux of its own field.",
            ));
            let back = self.drive_back_action_bound(player);
            out.push(note(
                false,
                &format!(
                    "Circuits drive one way: the particles' back-action on them (the charge a \
                     particle induces on a plate, the voltage it induces in a coil) is left \
                     out, at most {back:.1e} of the drive here."
                ),
            ));
            let ratio = self.drive_quasi_static_ratio(player);
            out.push(note(
                false,
                &format!(
                    "Circuits are quasi-static: the plates' and coils' own retardation and \
                     radiation, and the magnetic field of the charging currents, are left out \
                     (light time across the arena × the drives' fastest relative rate = \
                     {ratio:.1e})."
                ),
            ));
        }
        if !self.clouds.is_empty() {
            out.push(note(
                true,
                "Charge clouds: spheres of uniform charge that particles fly through, with their exact field inside (linear: a harmonic restoring force) and outside (a point charge). They are fixed: the cloud itself does not respond.",
            ));
        }
        if !self.disturbances.is_empty() {
            out.push(note(
                true,
                "Disturbances: given uniform fields and plane waves from outside the arena, \
                 applied exactly.",
            ));
        }
        if self.physics.c.is_some() {
            out.push(if self.physics.radiation_reaction {
                note(
                    false,
                    "Radiation reaction by the Landau–Lifshitz force: valid while it is small \
                     compared with the Lorentz force (ratio in the flight details).",
                )
            } else {
                note(
                    false,
                    "Radiation is neglected: the energy the particle would radiate is shown in \
                     the flight details and must stay below 1e-10 of its launch energy.",
                )
            });
        }
        if self
            .shots
            .iter()
            .any(|s| s.detector.acceptance.is_some_and(|a| a.radiation.is_some()))
        {
            out.push(note(
                true,
                "Radiation goal: the energy per steradian of the particle's far-zone \
                 radiation (Liénard; the spectrum from Jackson's radiation integral), exact for \
                 the computed motion, averaged over the receiver's directions in the plane.",
            ));
            out.push(note(
                false,
                "Radiation goal: the particle is taken to fly uniformly before its launch and \
                 after it enters the detector, so no radiation is counted from its launch or \
                 from the detector stopping it (unless the goal says the target stops it).",
            ));
            if self.has_beams() {
                out.push(note(
                    true,
                    "Radiation goal with other particles flying: the receiver sees them all, so \
                     it measures every charge's radiation until the particle arrives, their \
                     far fields added (they interfere: in phase they reinforce, in antiphase \
                     they cancel). Each flies uniformly before its launch and after its end.",
                ));
            }
            if self.shots.iter().any(|s| {
                s.detector
                    .acceptance
                    .and_then(|a| a.radiation)
                    .is_some_and(|r| r.window.is_some())
            }) {
                out.push(note(
                    false,
                    "Steady receiver: the mean power per steradian over its window, the \
                     radiation weighted by a smooth window (sin⁴) so that with a band a line \
                     outside it leaks in by about 1e-6; what is left of the start's transient \
                     (damped by the radiation reaction, at the atom's own frequency) counts \
                     unless the band leaves it out.",
                ));
            }
            if self.shots.iter().any(|s| {
                s.detector
                    .acceptance
                    .and_then(|a| a.radiation)
                    .is_some_and(|r| r.abrupt_stop)
            }) {
                out.push(note(
                    false,
                    "Radiation goal: the target stops the particle instantly; its flash is \
                     exact for frequencies far below the inverse of a real stopping time.",
                ));
            }
        }
        if !self.conductors.is_empty() {
            out.push(note(
                false,
                "Metal spheres: image charges plus a fitted surface charge; the boundary error \
                 is measured (below 1e-10 in the shipped levels), and the preview and the \
                 verification use different resolutions.",
            ));
        }
        if electrodes || player.iter().any(|e| e.kind == ElementKind::Plate) {
            out.push(note(
                false,
                "Electrodes: surface charge on flat triangles (boundary element method); the \
                 discretization error enters the verification through two mesh resolutions.",
            ));
            out.push(note(
                false,
                "Electrodes: the charge the particle itself induces on them (its image force) \
                 is neglected; its bound is shown in the flight details.",
            ));
        }
        if !self.ferrites.is_empty() {
            out.push(note(
                false,
                "Ferrites: bound magnetic charge on flat triangles from the continuity of the \
                 normal B, as for the dielectrics (the same accuracy); magnetized by the static \
                 magnets and coils only (PHYSICS.md §2.7).",
            ));
        }
        if !self.dielectrics.is_empty() || !self.dielectric_spheres.is_empty() {
            out.push(note(
                false,
                "Dielectrics: bound charge on flat triangles from the continuity of the normal \
                 displacement, averaged over each triangle; about 1 % off at the verification's \
                 resolution for glass (ε = 4; a sphere 1.4 %), up to 5 % as ε grows large \
                 (PHYSICS.md §2.7). \
                 Static response (no dispersion), no charge the particle induces in them.",
            ));
        }
        if !self.conductors.is_empty() || electrodes {
            out.push(note(
                false,
                "Metal responds electrostatically (instantly): valid for particles slow \
                 compared with light crossing the setup.",
            ));
        }
        if self.has_beams() {
            let charged = self.shots.iter().any(|s| s.particle.charge != 0.0);
            out.push(
                if self.interacts() && self.physics.c.is_none() {
                    note(
                        true,
                        "Beam: all particles fly together and repel each other with the exact \
                     Coulomb force (for c = ∞ that is the whole interaction).",
                    )
                } else if self.interacts() && !self.physics.beam_retarded {
                    note(
                        false,
                        if self.physics.radiation_reaction {
                            "Beam: all particles fly together. Each feels the others' fields, \
                     radiation included, computed from where they are now as if they had \
                     always moved in the fields they feel now, those fields changing along \
                     their paths as they do now (exact in the velocities, so the magnetic \
                     attraction that weakens the repulsion of a fast beam by 1/γ² is \
                     included; exact to the rate of change of the fields, so the others' \
                     share of each one's radiation damping is included too; approximate \
                     where the fields change unevenly over the light time between them: \
                     the estimated error is shown per particle, and checked against the \
                     exact retarded fields in the level tests). Each particle \
                     feels its own radiation reaction (Landau–Lifshitz)."
                        } else {
                            "Beam: all particles fly together. Each feels the others' fields, \
                     radiation included, computed from where they are now as if they had \
                     always moved in the fields they feel now, those fields changing along \
                     their paths as they do now (exact in the velocities, so the magnetic \
                     attraction that weakens the repulsion of a fast beam by 1/γ² is \
                     included; exact to the rate of change of the fields, so the others' \
                     share of each one's radiation damping is included too; approximate \
                     where the fields change unevenly over the light time between them: \
                     the estimated error is shown per particle, and checked against the \
                     exact retarded fields in the level tests)."
                        },
                    )
                } else if self.interacts() && self.physics.radiation_reaction {
                    note(
                        true,
                        "Beam: all particles fly together and act on each other with their exact \
                     retarded (Liénard–Wiechert) fields, magnetic attraction and radiation \
                     included, and each feels its own radiation reaction (Landau–Lifshitz). \
                     Before launch they are taken to keep their launch acceleration.",
                    )
                } else if self.interacts() {
                    note(
                        true,
                        "Beam: all particles fly together and act on each other with their exact \
                     retarded (Liénard–Wiechert) fields, magnetic attraction and radiation \
                     included. Before launch they are taken to keep their launch acceleration. Each \
                     particle's own radiation reaction is left out (its radiated energy is \
                     checked to be negligible).",
                    )
                } else if self.physics.c.is_none() && !charged {
                    note(
                        true,
                        "Beam of neutral particles: for c = ∞ their moments do not interact (the \
                     interaction scales as 1/c²).",
                    )
                } else {
                    note(
                        false,
                        "Beam: the particles' interaction with each other is left out.",
                    )
                },
            );
            if self.interacts() {
                out.push(note(
                    false,
                    "Absorbed particles: one that hits a body stops there and its charge stays, at rest, still acting on the others; one that leaves the arena flies on and keeps acting (it only counts as lost).",
                ));
                out.push(note(
                    false,
                    if self.physics.instant_drain {
                        "Detectors (instant drain): a particle entering its detector is carried away at once; its field vanishes where the news of its absorption, travelling at c, has arrived (at once for c = ∞)."
                    } else {
                        "Detectors: each is the mouth of a deep grounded Faraday cup. A particle entering it flies on into the cup while the cup screens its charge, seen from outside, as e^(−k·depth) (k = π/width of the mouth: the slowest-decaying mode of a grounded pipe), so its field fades over the entry time instead of vanishing at once. Approximations: the faster-decaying modes and the charge induced on the rim are left out (they matter near the mouth only), and so are the fields of the currents carrying the charge away."
                    },
                ));
            }
            out.push(note(
                false,
                "The beam is a fixed sample (quasi-random) of its distribution: the same \
                 particles every time, so the result is verifiable.",
            ));
        } else if self.shots.len() > 1 {
            out.push(note(
                false,
                "Shots fly one at a time: particles do not interact with each other.",
            ));
        }
        out.sort_by_key(|n| !n.exact);
        out
    }
}

/// A tube level's model (PHYSICS.md §2.11): its own world, z-invariant, with its own
/// approximations.
fn tube_notes(spec: &crate::tube::TubeSpec) -> Vec<ModelNote> {
    let mut out = vec![
        note(
            true,
            "A z-invariant tube: every electrode, carrier and field is the same along z \
             (long prisms and lines of charge), the symmetry the textbook diode laws assume. \
             Its fields are exact 2D electrostatics: a line charge's potential is −2λ ln r.",
        ),
        note(
            false,
            "Newtonian motion (carriers below 0.1 c; at tube energies about 0.03 c): the \
             magnetic forces between carriers and retardation, of order v²/c², are left out.",
        ),
        note(
            false,
            "Electrodes: surface charge on straight segments (boundary elements), the \
             emitters meshed four times finer. The verification halves every scale \
             (segments, time step, macroparticle weight, softening); its error bar is 1.5 \
             times its difference to the preview.",
        ),
        note(
            false,
            "Carriers are macroparticles: each a line of many carriers following one \
             carrier's path. Only the graininess of the space charge depends on their weight \
             (1/150 of the charge the emitting face holds, in the preview); their mutual \
             force is softened within an eighth of a segment.",
        ),
        note(
            false,
            "Space-charge-limited emission: each step a segment emits the charge its surface \
             holds, released at rest half a segment out. Skipping that first layer makes the \
             current high by an error of first order in the resolution (a few per cent at \
             the verification's), which the error bar covers. A carrier coming back is \
             absorbed by its emitter on reaching the layer it was released from.",
        ),
        note(
            false,
            "A grounded envelope half a cell outside the board encloses the tube: what \
             leaves the board ends on it (counted as lost).",
        ),
        note(
            false,
            "Many interacting carriers make a run chaotic: a change at the level of rounding \
             moved verified results by up to 0.7 %, inside the error bars (the \
             macroparticles' noise).",
        ),
        note(
            false,
            "Metal responds electrostatically (instantly): valid for carriers slow compared \
             with light crossing the tube.",
        ),
    ];
    if spec.b_z != 0.0 || spec.solenoid.is_some() {
        out.push(note(
            true,
            "A uniform magnetic field along z (a long solenoid around the tube): the \
             carriers' gyration is exact (Boris rotation).",
        ));
    }
    if spec.saturation.is_some() {
        out.push(note(
            false,
            "Temperature-limited cathode: it emits at most its saturation current (the \
             Richardson current of its temperature), whatever the field.",
        ));
    }
    if spec.ions.is_some() {
        out.push(note(
            false,
            "Ions: lighter per charge than real ones (real ions are at least 1836 times \
             heavier than electrons per charge). In a static field the mass changes neither \
             the paths nor the electrons' gain from the ions' space charge, only the time \
             scale and the ion current (√(m_e/m_i) of the electrons' in a planar diode).",
        ));
    }
    if !spec.circuits.is_empty() {
        let coils = spec.circuits.iter().any(|c| c.inductance > 0.0);
        out.push(note(
            false,
            if coils {
                "Circuits: lumped, per unit length along z, acting on the carriers and acted \
                 on through the charge they induce (Ramo). Resistive branches are integrated \
                 by backward Euler (first order in the step); coils by the trapezoidal rule, \
                 so an LC tank keeps its energy and is only detuned, by (ω Δt)²/12 of its \
                 frequency."
            } else {
                "Circuits: lumped, per unit length along z, acting on the carriers and acted \
                 on through the charge they induce (Ramo); integrated by backward Euler \
                 (first order in the step)."
            },
        ));
    }
    if !spec.metals.is_empty() {
        out.push(note(
            false,
            "Metals (the Drude model): classical carriers over the ions' uniform charge, \
             bouncing off fixed ion cores placed at random (an amorphous metal). A real \
             metal's conduction electrons are a degenerate quantum gas, far denser and with \
             far shorter free times: this is a scaled classical plasma, whose screening \
             length and mean free path are well below the conductor's size, so that it \
             conducts as a metal does. What carries over: Ohm's law, σ = n q² τ/m, the Joule \
             heat, the Hall coefficient (1/(nq) times this model's Hall factor, 4/π). What \
             does not: the thermal speed, the carriers' heat capacity.",
        ));
        out.push(note(
            false,
            "A collision is a mirror bounce off a core with a new speed drawn from the \
             lattice's temperature: the lattice takes the energy the field gave (a bath at \
             a fixed temperature: it does not warm). The mean free path ℓ = 1/(2 a n) of \
             cores of radius a, n per area; the conductivity (3/8) √(π/2) ρ|q/m| ℓ/σ_v \
             of the Boltzmann equation holds within 3 % for the dilute cores allowed \
             (PHYSICS.md §2.11; the carriers warm a little in strong fields, lowering it).",
        ));
        out.push(note(
            false,
            "Contacts: a contact takes in the carriers that reach it and gives off, at \
             rest, as many as keep the metal next to it neutral (space-charge limited); \
             they gain the lattice's temperature at their first collision. A current goal \
             on a contact counts what it takes in less what it gives off (its lead's \
             current).",
        ));
        out.push(note(
            false,
            "The carriers' current fluctuates (the thermal, Johnson noise of a resistor): \
             averaged over the goal's window it is uncertain by a few per cent, which the \
             error bar covers.",
        ));
    }
    if let Some(o) = spec.output {
        out.push(note(
            false,
            &format!(
                "The output: circuit {}'s load is what the antenna strapped to it radiates. \
                 The field maps show the run's electromagnetic field at the speed of light \
                 c = {} (a display parameter, one-way coupled: the carriers' motion stays \
                 Newtonian), the board's edge absorbing the waves.",
                o.circuit + 1,
                o.c
            ),
        ));
    }
    out.sort_by_key(|n| !n.exact);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every level states its model, and the notes follow what the level contains.
    #[test]
    fn notes_follow_the_level() {
        let has =
            |l: &Level, p: &[Element], s: &str| l.model_notes(p).iter().any(|n| n.text.contains(s));
        let bend = crate::shipped("first_bend");
        assert!(!has(&bend, &[], "Electrodes"));
        // A tube is its own world: z-invariant, not the plane of a 3D one.
        let diode = crate::shipped("vacuum_diode");
        assert!(has(&diode, &[], "z-invariant") && !has(&diode, &[], "3D world"));
        assert!(!has(&diode, &[], "Circuits"));
        let rectifier = crate::shipped("rectifier");
        assert!(has(&rectifier, &[], "backward Euler") && !has(&rectifier, &[], "coils by"));
        let plates = crate::shipped("build_a_deflector");
        assert!(has(&plates, &[], "image force"));
        assert!(has(&plates, &[], "Radiation is neglected"));
        let sync = crate::shipped("synchrotron_light");
        assert!(has(&sync, &[], "Landau–Lifshitz"));
        assert!(has(&sync, &[], "Coils"));
        // Exact statements come first.
        let n = sync.model_notes(&[]);
        assert!(n.windows(2).all(|w| w[0].exact || !w[1].exact));
    }
}
