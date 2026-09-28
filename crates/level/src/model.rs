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
        if self.shots.iter().any(|s| s.particle.moment != 0.0) {
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
                if self.physics.beam_interaction && self.physics.c.is_none() {
                    note(
                        true,
                        "Beam: all particles fly together and repel each other with the exact \
                     Coulomb force (for c = ∞ that is the whole interaction).",
                    )
                } else if self.physics.beam_interaction && !self.physics.beam_retarded {
                    note(
                        false,
                        if self.physics.radiation_reaction {
                            "Beam: all particles fly together. Each feels the others' fields, \
                     radiation included, computed from where they are now as if they had \
                     kept their present acceleration (exact in the velocities, so the \
                     magnetic attraction that weakens the repulsion of a fast beam by 1/γ² \
                     is included; approximate in how the acceleration changes: the \
                     estimated error is shown per particle, and checked against the exact \
                     retarded fields in the level tests). Each particle \
                     feels its own radiation reaction (Landau–Lifshitz)."
                        } else {
                            "Beam: all particles fly together. Each feels the others' fields, \
                     radiation included, computed from where they are now as if they had \
                     kept their present acceleration (exact in the velocities, so the \
                     magnetic attraction that weakens the repulsion of a fast beam by 1/γ² \
                     is included; approximate in how the acceleration changes: the \
                     estimated error is shown per particle, and checked against the exact \
                     retarded fields in the level tests)."
                        },
                    )
                } else if self.physics.beam_interaction && self.physics.radiation_reaction {
                    note(
                        true,
                        "Beam: all particles fly together and act on each other with their exact \
                     retarded (Liénard–Wiechert) fields, magnetic attraction and radiation \
                     included, and each feels its own radiation reaction (Landau–Lifshitz). \
                     Before launch they are taken to move uniformly; a particle that is \
                     absorbed stops acting once the news has travelled at c.",
                    )
                } else if self.physics.beam_interaction {
                    note(
                        true,
                        "Beam: all particles fly together and act on each other with their exact \
                     retarded (Liénard–Wiechert) fields, magnetic attraction and radiation \
                     included. Before launch they are taken to move uniformly; a particle \
                     that is absorbed stops acting once the news has travelled at c. Each \
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
            if self.physics.beam_interaction {
                out.push(note(
                    false,
                    "Absorbed particles: one that hits a body stops there and its charge stays, at rest, still acting on the others; one that enters the detector (a grounded Faraday cup) is carried away; one that leaves the arena flies on and keeps acting (it only counts as lost).",
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
