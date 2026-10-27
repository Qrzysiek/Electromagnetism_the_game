//! Labels on the board and in the panel: the supplies' badges, the Hall arrows'
//! labels, the level's reality remarks, the level documents' buttons.

use super::*;

/// Each supply's number on the map, beside every electrode it drives (its tunable electrode
/// and the electrodes linked to it), or at a coil's or solenoid's terminal: the same number
/// as its slider's row in the panel.
pub(super) fn supply_badges(
    ctx: &egui::Context,
    game: &Game,
    (cam, cam_transform): (&Camera, &GlobalTransform),
    panel_left: f32,
) {
    let level = &game.editor.level;
    let targets = level.supply_targets();
    if targets.is_empty() {
        return;
    }
    let focus = game.focused_supply();
    let painter = ctx.layer_painter(egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new("supply badges"),
    ));
    let outlines = level.electrode_outlines(&[]);
    for (k, (centre, target)) in targets.iter().enumerate() {
        // Where its badges go: the upper right of each electrode it drives (world).
        let anchors: Vec<physics::DVec3> = match target {
            level::SupplyTarget::Electrode(_) => level
                .supply_electrodes(*target)
                .iter()
                .filter_map(|&i| outlines.get(i))
                .map(|o| match *o {
                    level::Outline::Rect(c) => {
                        let x = c.iter().map(|p| p.x).fold(f64::MIN, f64::max);
                        let y = c.iter().map(|p| p.y).fold(f64::MIN, f64::max);
                        physics::DVec3::new(x + 0.3, y + 0.3, 0.0)
                    }
                    level::Outline::Circle { center, radius: r }
                    | level::Outline::Ring {
                        center, outer: r, ..
                    } => center + physics::DVec3::new(0.75 * r + 0.3, 0.75 * r + 0.3, 0.0),
                })
                .collect(),
            _ => {
                let p = level.grid.position(*centre);
                vec![physics::DVec3::new(p.x + 0.45, p.y + 0.45, 0.0)]
            }
        };
        let ring = if focus == Some(*centre) {
            SUPPLY_FOCUS
        } else {
            SUPPLY_FRAME
        };
        for a in anchors {
            #[allow(clippy::cast_possible_truncation)]
            let world = Vec3::new(a.x as f32, a.y as f32, 0.0);
            let Ok(v) = cam.world_to_viewport(cam_transform, world) else {
                continue;
            };
            if v.x > panel_left - 10.0 {
                continue;
            }
            let at = egui::pos2(v.x, v.y);
            painter.circle(
                at,
                8.5,
                egui::Color32::from_rgb(22, 24, 30),
                egui::Stroke::new(1.6, ring),
            );
            painter.text(
                at,
                egui::Align2::CENTER_CENTER,
                format!("{}", k + 1),
                egui::FontId::proportional(11.0),
                egui::Color32::WHITE,
            );
        }
    }
}

/// Each Hall diagram's labels (`draw::hall_diagram`), just past its arrows' tips: what each
/// arrow is, on the board (the legend says more).
pub(super) fn hall_labels(
    ctx: &egui::Context,
    game: &Game,
    (cam, cam_transform): (&Camera, &GlobalTransform),
    panel_left: f32,
) {
    if game.hall.is_empty() {
        return;
    }
    let painter = ctx.layer_painter(egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new("hall labels"),
    ));
    let to_screen = |p: Vec2| {
        cam.world_to_viewport(cam_transform, p.extend(0.0))
            .ok()
            .map(|v| egui::pos2(v.x, v.y))
    };
    for h in &game.hall {
        let (centre, arrows) = crate::draw::hall_diagram(h);
        for a in &arrows {
            let d = a.tip - centre;
            if d.length() < 0.2 {
                continue;
            }
            let Some(at) = to_screen(a.tip + d.normalize() * 0.3) else {
                continue;
            };
            if at.x > panel_left - 30.0 {
                continue;
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let [r, g, b] = a.rgb.map(|c| (c * 255.0).round() as u8);
            let galley = painter.layout_no_wrap(
                a.label.to_owned(),
                egui::FontId::proportional(12.0),
                egui::Color32::from_rgb(r, g, b),
            );
            let rect = egui::Align2::CENTER_CENTER.anchor_size(at, galley.size());
            painter.rect_filled(
                rect.expand(2.0),
                3.0,
                egui::Color32::from_rgba_unmultiplied(12, 13, 17, 200),
            );
            painter.galley(rect.min, galley, egui::Color32::WHITE);
        }
    }
}

/// What the level represents and how far it is the real thing (`level::Reality`): the
/// headline always in view, the remarks folded under it (the owner, 2026-10-21: leaving them
/// out gives a wrong impression; folded, they do not clutter the panel).
pub(super) fn reality(ui: &mut egui::Ui, reality: &level::Reality) {
    use level::RemarkKind;
    if reality.is_empty() {
        return;
    }
    // The headline wraps on its own line (a long collapsing header would widen the panel).
    ui.horizontal_wrapped(|ui| {
        ui.label(egui::RichText::new("Represents:").strong())
            .on_hover_text("The real thing this level shows; how far it is that thing, below");
        ui.label(&reality.represents);
    });
    let simplified = reality
        .remarks
        .iter()
        .filter(|r| r.kind != RemarkKind::Kept)
        .count();
    egui::CollapsingHeader::new(format!(
        "How real: {} kept, {simplified} simplified or left out",
        reality.remarks.len() - simplified
    ))
    .id_salt("reality")
    // Opened for developer captures with EM_REALITY=1.
    .default_open(std::env::var("EM_REALITY").is_ok_and(|v| v == "1"))
    .show(ui, |ui| {
        for r in &reality.remarks {
            ui.horizontal_wrapped(|ui| {
                let (colour, meaning) = match r.kind {
                    RemarkKind::Kept => (
                        egui::Color32::from_rgb(110, 210, 120),
                        "As in the real thing: this carries over",
                    ),
                    RemarkKind::Simplified => (
                        egui::Color32::from_rgb(230, 180, 70),
                        "Idealized, scaled or approximated against the real thing",
                    ),
                    RemarkKind::Omitted => (
                        egui::Color32::from_rgb(225, 120, 95),
                        "Absent from the model: the real thing has it",
                    ),
                };
                ui.label(egui::RichText::new(r.kind.label()).small().color(colour))
                    .on_hover_text(meaning);
                ui.label(egui::RichText::new(&r.line).strong());
            });
            // Full size: folded away, it takes no room until opened (no small print).
            if !r.detail.is_empty() {
                ui.label(&r.detail);
            }
            ui.add_space(4.0);
        }
    })
    .header_response
    .on_hover_text(
        "The real thing this level shows: what carries over from it, what the level \
         simplifies, and what it leaves out. The numerical model's own approximations are \
         under Physics model and its limits.",
    );
}

/// The level's document and the Foundations (docs/levels/README.md), opened in the system's PDF
/// viewer; nothing for a level without a document yet.
pub(super) fn level_doc_buttons(ui: &mut egui::Ui, game: &Game) {
    use crate::level_docs::{Doc, FOUNDATIONS};
    let Some(path) = game.level_paths.get(game.level_index) else {
        return;
    };
    let doc = crate::level_docs::for_level(path);
    if matches!(doc, Doc::Missing) {
        return;
    }
    ui.horizontal(|ui| {
        match &doc {
            Doc::Ready(pdf) => {
                if ui
                    .button("The physics of this level")
                    .on_hover_text(
                        "This level derived from first principles, down to the numbers the game \
                         computes (a PDF: docs/levels/)",
                    )
                    .clicked()
                {
                    crate::level_docs::open(pdf);
                }
            }
            Doc::NotBuilt(command) => {
                ui.add_enabled(false, egui::Button::new("The physics of this level"))
                    .on_disabled_hover_text(format!("Not typeset yet: {command}"));
            }
            Doc::Missing => {}
        }
        if let Doc::Ready(pdf) = crate::level_docs::find(FOUNDATIONS)
            && ui
                .button("Foundations")
                .on_hover_text(
                    "The units, principles and conventions every level's document builds on",
                )
                .clicked()
        {
            crate::level_docs::open(&pdf);
        }
    });
}
