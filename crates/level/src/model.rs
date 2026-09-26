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
        if self.shots.len() > 1 {
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

    fn shipped(file: &str) -> Level {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../levels")
            .join(file);
        Level::from_json(&std::fs::read_to_string(path).unwrap()).unwrap()
    }

    /// Every level states its model, and the notes follow what the level contains.
    #[test]
    fn notes_follow_the_level() {
        let has =
            |l: &Level, p: &[Element], s: &str| l.model_notes(p).iter().any(|n| n.text.contains(s));
        let bend = shipped("01_first_bend.json");
        assert!(!has(&bend, &[], "Electrodes"));
        let plates = shipped("20_build_a_deflector.json");
        assert!(has(&plates, &[], "image force"));
        assert!(has(&plates, &[], "Radiation is neglected"));
        let sync = shipped("37_synchrotron_light.json");
        assert!(has(&sync, &[], "Landau–Lifshitz"));
        assert!(has(&sync, &[], "Coils"));
        // Exact statements come first.
        let n = sync.model_notes(&[]);
        assert!(n.windows(2).all(|w| w[0].exact || !w[1].exact));
    }
}
