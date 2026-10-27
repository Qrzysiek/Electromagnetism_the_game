//! Property editors for every part of a level (sandbox).
//!
//! **Completeness invariant:** everything a level file can contain must be editable here
//! (or by a sandbox tool), so any shipped level can be rebuilt by hand. Every level type
//! is destructured *exhaustively* (no `..`): adding a field to the level format fails to
//! compile until it gets a control here, or is explicitly marked as not user-editable
//! with a reason.

use bevy_egui::egui;
use level::beam::{BeamSpec, Distribution};
use level::tube::{
    Carriers, CurrentGoal, DEFAULT_STEP, HeatGoal, IonSource, Measure, MetalSpec, PotentialGoal,
    ProjectileSpec, TL_RATE, TubeCircuit, TubeOutput, TubeSpec,
};
use level::{
    Cloud, Coil, Conductor, ConductorBias, Detector, DetectorAcceptance, Dielectric,
    DielectricSphere, Disturbance, Drive, Electrode, ElectrodeShape, Element, ElementKind,
    FerritePlate, FreeParticle, Grid, Launch, Level, Limits, Node, ParticleSpec, RadiationGoal,
    Reality, Region2, Remark, RemarkKind, Shot, Source, Switch, TolerancesSpec, Wave, WorldPhysics,
};

use crate::ui::{fmt_si, parse_si};

mod bodies;
mod check;
mod circuits;
mod inputs;
mod limits;
mod shots;
mod tube;
use bodies::{
    edit_clouds, edit_coils, edit_conductors, edit_dielectric_spheres, edit_dielectrics,
    edit_electrodes, edit_elements, edit_ferrites,
};
pub use check::check_editable;
use circuits::{drive_editable, drive_editor, edit_source, source_editable};
pub use inputs::list_to_text;
use inputs::{
    bias_editor, kind_combo, label_editor, node, non_negative, parse_list, parse_signed_list,
    positive, row, sci, si,
};
use limits::{edit_disturbances, edit_limits};
use shots::{
    acceptance_in_range, edit_free_particles, edit_gates, edit_grid, edit_physics, edit_shot,
};
use tube::edit_tube;

// Ranges the editors allow; `check_editable` uses the same constants, so every shipped
// level is checked to lie within reach of the editors.
const GRID_CELLS: std::ops::RangeInclusive<u32> = 4..=200;
const SUBDIVISION: std::ops::RangeInclusive<u32> = 1..=8;
const MAX_POSITIVE: f64 = 1e12;
const MIN_POSITIVE: f64 = 1e-12;
/// Largest relative permittivity of a dielectric (water is 80; beyond 1e4 a dielectric is
/// metal for every purpose here), and permeability of a ferrite (soft ferrites reach
/// thousands).
const MAX_PERMITTIVITY: f64 = 1e4;
const MAX_RADIUS: f64 = 5.0;
/// Largest beam (particles; the interaction costs grow as N²).
const MAX_BEAM: u32 = 200;
const MAX_COUNT: u32 = 20;
const TOLERANCE: std::ops::RangeInclusive<f64> = 1e-18..=1.0;
/// Width of list text fields, so that grid rows fit the panel.
const LIST_WIDTH: f32 = 100.0;

/// Text buffers for list fields (edited as comma-separated text).
#[derive(Default)]
pub struct EditTexts {
    pub magnitudes: String,
    pub magnet_strengths: String,
    pub antenna_amplitudes: String,
    pub antenna_omegas: String,
    pub plate_voltages: String,
    pub supply_voltages: String,
    pub coil_rates: String,
    pub solenoid_fields: String,
    pub free_charges: String,
    pub free_speeds: String,
}

/// Result of editing: whether a text field has focus, and a requested grid refinement.
pub struct EditResult {
    pub focus: bool,
    pub refine_by: Option<u32>,
}

/// Editors for every field of `level`; the reference solution is edited in the
/// sandbox's Solution section, placements by the map tools.
pub fn edit_level(
    ui: &mut egui::Ui,
    level: &mut Level,
    texts: &mut EditTexts,
    active_shot: usize,
) -> EditResult {
    let Level {
        // Written on save: format and engine version of the file.
        format_version: _,
        engine_version: _,
        name,
        description,
        reality,
        grid,
        physics,
        shots,
        elements,
        coils,
        limits,
        // Edited in the Solution section (store own elements or a solver result).
        reference_solution: _,
        disturbances,
        conductors,
        clouds,
        free_particles,
        electrodes,
        dielectrics,
        dielectric_spheres,
        ferrites,
        gates,
        tube,
    } = level;
    let mut focus = false;
    let mut refine_by = None;
    let g = *grid;
    egui::CollapsingHeader::new("Level and world")
        .default_open(true)
        .show(ui, |ui| {
            // Name and description take the full width (a multi-line field inside the
            // grid would overlap the rows below it).
            ui.label("Name");
            focus |= ui
                .add(egui::TextEdit::singleline(name).desired_width(f32::INFINITY))
                .has_focus();
            ui.label("Description");
            focus |= ui
                .add(
                    egui::TextEdit::multiline(description)
                        .desired_rows(4)
                        .desired_width(f32::INFINITY),
                )
                .has_focus();
            focus |= edit_reality(ui, reality);
            egui::Grid::new("level_props")
                .num_columns(2)
                .show(ui, |ui| {
                    focus |= edit_grid(ui, grid, &mut refine_by);
                    focus |= edit_physics(ui, physics);
                });
        });
    if let Some(s) = shots.get_mut(active_shot) {
        egui::CollapsingHeader::new(format!("Shot {}", active_shot + 1))
            .default_open(true)
            .show(ui, |ui| {
                egui::Grid::new("shot_props").num_columns(2).show(ui, |ui| {
                    focus |= edit_shot(ui, s, &g);
                });
            });
    }
    egui::CollapsingHeader::new(format!("Level elements ({})", elements.len()))
        .show(ui, |ui| focus |= edit_elements(ui, elements, &g));
    egui::CollapsingHeader::new(format!("Coils ({})", coils.len()))
        .show(ui, |ui| focus |= edit_coils(ui, coils, &g));
    egui::CollapsingHeader::new(format!("Electrodes ({})", electrodes.len()))
        .show(ui, |ui| focus |= edit_electrodes(ui, electrodes, &g));
    egui::CollapsingHeader::new(format!("Dielectrics ({})", dielectrics.len()))
        .show(ui, |ui| {
            focus |= edit_dielectrics(ui, dielectrics, &g);
            focus |= edit_dielectric_spheres(ui, dielectric_spheres, &g);
        })
        .header_response
        .on_hover_text(
            "Insulating boxes of permittivity ε: their bound charge screens a field partly \
             (metal screens it fully)",
        );
    egui::CollapsingHeader::new(format!("Ferrites ({})", ferrites.len()))
        .show(ui, |ui| focus |= edit_ferrites(ui, ferrites, &g))
        .header_response
        .on_hover_text(
            "Magnetically soft insulating boxes of permeability μ: they draw in and bend \
             magnetic fields (from magnets and steady coils), not electric ones",
        );
    egui::CollapsingHeader::new(format!("Metal spheres ({})", conductors.len()))
        .show(ui, |ui| focus |= edit_conductors(ui, conductors, &g));
    egui::CollapsingHeader::new(format!("Charge clouds ({})", clouds.len()))
        .show(ui, |ui| focus |= edit_clouds(ui, clouds, &g));
    egui::CollapsingHeader::new(format!("Free particles ({})", free_particles.len()))
        .show(ui, |ui| {
            focus |= edit_free_particles(ui, free_particles, &g)
        })
        .header_response
        .on_hover_text(
            "Particles that move and interact with the shots (targets, partners); with a \
             detector they must arrive too",
        );
    egui::CollapsingHeader::new(format!("Disturbances ({})", disturbances.len()))
        .show(ui, |ui| focus |= edit_disturbances(ui, disturbances));
    egui::CollapsingHeader::new(format!("Gates ({})", gates.len()))
        .show(ui, |ui| focus |= edit_gates(ui, gates, &g))
        .header_response
        .on_hover_text("Stages every flight must pass, in order, before its detector counts");
    let n_electrodes = electrodes.len();
    egui::CollapsingHeader::new(if tube.is_some() { "Tube (on)" } else { "Tube" })
        .show(ui, |ui| focus |= edit_tube(ui, tube, n_electrodes, &g))
        .header_response
        .on_hover_text(
            "A z-invariant vacuum tube: the electrodes are prisms along z, the cathode emits \
             electrons under the space-charge limit, the goal is an electrode's current",
        );
    egui::CollapsingHeader::new("Player limits")
        .default_open(true)
        .show(ui, |ui| {
            egui::Grid::new("limit_props")
                .num_columns(2)
                .show(ui, |ui| {
                    focus |= edit_limits(ui, limits, texts, &g);
                });
        });
    EditResult { focus, refine_by }
}

/// What the level represents and how it differs from the real thing (`level::Reality`):
/// a headline and remarks, each a kind, a short line and its explanation.
fn edit_reality(ui: &mut egui::Ui, reality: &mut Reality) -> bool {
    let mut focus = false;
    ui.label("Represents").on_hover_text(
        "The real thing the level shows, in a few words (the header of its Represents section)",
    );
    focus |= ui
        .add(egui::TextEdit::singleline(&mut reality.represents).desired_width(f32::INFINITY))
        .has_focus();
    let mut remove = None;
    for (i, r) in reality.remarks.iter_mut().enumerate() {
        let Remark { kind, line, detail } = r;
        ui.push_id(("remark", i), |ui| {
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("kind")
                    .selected_text(kind.label())
                    .show_ui(ui, |ui| {
                        for k in RemarkKind::ALL {
                            ui.selectable_value(kind, k, k.label());
                        }
                    });
                focus |= ui
                    .add(egui::TextEdit::singleline(line).desired_width(200.0))
                    .has_focus();
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            focus |= ui
                .add(
                    egui::TextEdit::multiline(detail)
                        .desired_rows(2)
                        .desired_width(f32::INFINITY),
                )
                .has_focus();
        });
    }
    if let Some(i) = remove {
        reality.remarks.remove(i);
    }
    if ui
        .small_button("+ remark")
        .on_hover_text(
            "What carries over from the real thing, what is simplified, what is left out",
        )
        .clicked()
    {
        reality.remarks.push(Remark {
            kind: RemarkKind::Simplified,
            line: String::new(),
            detail: String::new(),
        });
    }
    focus
}

#[cfg(test)]
mod tests;
