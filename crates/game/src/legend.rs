//! The board's legend (the owner, 2026-10-19: every element of the interface is
//! self-explanatory or explained where the eye is): for what the board draws now, a small
//! picture of each mark in its own colours, a line saying what it is, and on hover the
//! longer explanation. Tube levels first: their maps' arrows, the metals' marks, the
//! Hall effect's arrows, the leads, the electrodes' glows.

use bevy_egui::egui;

use crate::Game;
use crate::potential::MapMode;

/// A mark's picture in the legend (its colours as `draw.rs` draws it).
#[derive(Clone, Copy)]
enum Swatch {
    /// An arrow pointing right.
    Arrow(egui::Color32),
    /// A filled dot.
    Dot(egui::Color32),
    /// Two dots side by side (the carriers' two signs).
    TwoDots(egui::Color32, egui::Color32),
    /// A zigzag path ending in a dot (a traced carrier).
    Traced(egui::Color32),
    /// Two chevrons pointing right.
    Chevrons(egui::Color32),
    /// A wire with beads on it.
    Beads(egui::Color32),
    /// Rims of an electrode in these colours, side by side.
    Rims([egui::Color32; 3]),
    /// A rectangle's outline drawn three times, growing (a glow or a focus).
    Glow(egui::Color32, egui::Color32),
    /// A circle with a dot (out of the board) and one with a cross (into it).
    OutIn(egui::Color32, egui::Color32),
    /// A track ending in a large dot, beside a box (a projectile and its detector).
    Track(egui::Color32),
    /// Short paths in these colours, the last ending in a cross (computing, not verified,
    /// lost).
    Paths([egui::Color32; 3]),
    /// A thin line (a field line).
    Line(egui::Color32),
}

fn rgb(r: f32, g: f32, b: f32) -> egui::Color32 {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let c = |v: f32| (v * 255.0).round() as u8;
    egui::Color32::from_rgb(c(r), c(g), c(b))
}

/// The colours of `draw.rs` (kept beside each mark there).
const WIDTH: f32 = 34.0;
const HEIGHT: f32 = 14.0;

fn paint(ui: &mut egui::Ui, swatch: Swatch) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(WIDTH, HEIGHT), egui::Sense::hover());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 2.0, egui::Color32::from_gray(22));
    let (l, r, t, b) = (
        rect.left() + 4.0,
        rect.right() - 4.0,
        rect.top() + 3.0,
        rect.bottom() - 3.0,
    );
    let cy = rect.center().y;
    let stroke = |c: egui::Color32| egui::Stroke::new(1.5, c);
    match swatch {
        Swatch::Arrow(c) => {
            p.arrow(egui::pos2(l, cy), egui::vec2(r - l, 0.0), stroke(c));
        }
        Swatch::Dot(c) => {
            p.circle_filled(rect.center(), 3.0, c);
        }
        Swatch::TwoDots(a, b) => {
            p.circle_filled(egui::pos2(l + 8.0, cy), 3.0, a);
            p.circle_filled(egui::pos2(r - 8.0, cy), 3.0, b);
        }
        Swatch::Traced(c) => {
            let pts = vec![
                egui::pos2(l, b),
                egui::pos2(l + 7.0, t),
                egui::pos2(l + 12.0, b - 2.0),
                egui::pos2(l + 19.0, t + 2.0),
                egui::pos2(r - 3.0, cy),
            ];
            p.add(egui::Shape::line(pts, egui::Stroke::new(1.0, c)));
            p.circle_filled(egui::pos2(r - 3.0, cy), 2.5, c);
        }
        Swatch::Chevrons(c) => {
            for x0 in [l + 4.0, l + 15.0] {
                p.add(egui::Shape::line(
                    vec![
                        egui::pos2(x0, t),
                        egui::pos2(x0 + 7.0, cy),
                        egui::pos2(x0, b),
                    ],
                    stroke(c),
                ));
            }
        }
        Swatch::Beads(c) => {
            p.line_segment(
                [egui::pos2(l, cy), egui::pos2(r, cy)],
                egui::Stroke::new(1.0, c.gamma_multiply(0.5)),
            );
            for x in [l + 3.0, l + 13.0, l + 23.0] {
                p.circle_filled(egui::pos2(x, cy), 1.8, c);
            }
        }
        Swatch::Rims(cs) => {
            let w = (r - l) / 3.0;
            for (k, c) in cs.into_iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let x0 = l + w * k as f32;
                let box_ =
                    egui::Rect::from_min_max(egui::pos2(x0 + 1.0, t), egui::pos2(x0 + w - 1.0, b));
                p.rect_filled(box_, 0.0, egui::Color32::from_rgb(107, 115, 128));
                p.add(egui::Shape::closed_line(
                    vec![
                        box_.left_top(),
                        box_.right_top(),
                        box_.right_bottom(),
                        box_.left_bottom(),
                    ],
                    egui::Stroke::new(1.5, c),
                ));
            }
        }
        Swatch::Glow(warm, cool) => {
            for (k, c) in [warm, cool].into_iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let x0 = l + 1.0 + 14.0 * k as f32;
                let core = egui::Rect::from_min_max(
                    egui::pos2(x0 + 3.0, t + 3.0),
                    egui::pos2(x0 + 10.0, b - 3.0),
                );
                p.rect_filled(core, 0.0, egui::Color32::from_rgb(107, 115, 128));
                for (g, a) in [(1.0, 0.9), (2.0, 0.6), (3.0, 0.35)] {
                    let rr = core.expand(g);
                    p.add(egui::Shape::closed_line(
                        vec![
                            rr.left_top(),
                            rr.right_top(),
                            rr.right_bottom(),
                            rr.left_bottom(),
                        ],
                        egui::Stroke::new(1.0, c.gamma_multiply(a)),
                    ));
                }
            }
        }
        Swatch::OutIn(out, into) => {
            let (a, bb) = (egui::pos2(l + 6.0, cy), egui::pos2(r - 6.0, cy));
            p.circle_stroke(a, 5.0, egui::Stroke::new(1.2, out));
            p.circle_filled(a, 1.3, out);
            p.circle_stroke(bb, 5.0, egui::Stroke::new(1.2, into));
            let d = 3.4;
            p.line_segment(
                [bb + egui::vec2(-d, -d), bb + egui::vec2(d, d)],
                egui::Stroke::new(1.2, into),
            );
            p.line_segment(
                [bb + egui::vec2(-d, d), bb + egui::vec2(d, -d)],
                egui::Stroke::new(1.2, into),
            );
        }
        Swatch::Paths(cs) => {
            let w = (r - l) / 3.0;
            for (k, c) in cs.into_iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let x0 = l + w * k as f32;
                p.add(egui::Shape::line(
                    vec![
                        egui::pos2(x0, b),
                        egui::pos2(x0 + w * 0.5, t + 2.0),
                        egui::pos2(x0 + w - 2.0, t),
                    ],
                    stroke(c),
                ));
                if k == 2 {
                    let e = egui::pos2(x0 + w - 2.0, t + 1.0);
                    p.line_segment(
                        [e + egui::vec2(-2.5, -2.5), e + egui::vec2(2.5, 2.5)],
                        stroke(c),
                    );
                    p.line_segment(
                        [e + egui::vec2(-2.5, 2.5), e + egui::vec2(2.5, -2.5)],
                        stroke(c),
                    );
                }
            }
        }
        Swatch::Line(c) => {
            p.add(egui::Shape::line(
                vec![egui::pos2(l, b), egui::pos2(l + 10.0, cy), egui::pos2(r, t)],
                egui::Stroke::new(1.0, c),
            ));
        }
        Swatch::Track(c) => {
            let det = egui::Rect::from_min_max(egui::pos2(r - 7.0, t), egui::pos2(r, b));
            p.add(egui::Shape::closed_line(
                vec![
                    det.left_top(),
                    det.right_top(),
                    det.right_bottom(),
                    det.left_bottom(),
                ],
                egui::Stroke::new(1.0, c.gamma_multiply(0.7)),
            ));
            p.add(egui::Shape::line(
                vec![
                    egui::pos2(l, b),
                    egui::pos2(l + 8.0, cy),
                    egui::pos2(l + 16.0, t + 3.0),
                ],
                egui::Stroke::new(1.0, c.gamma_multiply(0.7)),
            ));
            p.circle_filled(egui::pos2(l + 16.0, t + 3.0), 3.0, c);
        }
    }
}

/// One row: the mark, what it is, and on hover what it means.
fn row(ui: &mut egui::Ui, swatch: Swatch, label: &str, detail: &str) {
    ui.horizontal(|ui| {
        paint(ui, swatch);
        ui.label(label);
    })
    .response
    .on_hover_text(detail);
}

/// The legend on the board, at its lower left (where the eye is, beside the marks; in the
/// panel it sat below the results, out of sight): a frame that folds to its title. It
/// names the marks that cannot be hovered (arrows, chevrons, glows, beads, colours); the
/// objects themselves (electrodes, metals, the solenoid) explain themselves on hover.
pub fn board(ctx: &egui::Context, game: &Game) {
    let tube = game.editor.level.is_tube();
    if !tube && game.flights.is_empty() && game.beams.is_empty() && !game.show_field_lines {
        return;
    }
    egui::Area::new(egui::Id::new("board_legend"))
        // Starts at the board's lower left; drag it aside (or fold it) where it covers
        // something.
        .pivot(egui::Align2::LEFT_BOTTOM)
        .default_pos(egui::pos2(10.0, ctx.viewport_rect().bottom() - 10.0))
        .movable(true)
        .show(ctx, |ui| {
            egui::Frame::popup(ui.style())
                .fill(egui::Color32::from_rgba_unmultiplied(18, 18, 22, 215))
                .show(ui, |ui| {
                    let shown = egui::CollapsingHeader::new("Legend")
                        .id_salt("board_legend_open")
                        .default_open(true)
                        .show(ui, |ui| {
                            if tube {
                                tube_rows(ui, game);
                            } else {
                                particle_rows(ui, game);
                            }
                        });
                    shown.header_response.on_hover_text(
                        "The marks drawn on the board: hover a line for what it means, or an \
                         object on the board for what it is.",
                    );
                });
        });
}

/// A tube level's marks that are not objects: what its board draws now
/// (`draw::draw_tube`), each with a few words and a hover text.
#[allow(clippy::too_many_lines)]
fn tube_rows(ui: &mut egui::Ui, game: &Game) {
    let level = &game.editor.level;
    let Some(spec) = level.tube.as_ref() else {
        return;
    };
    let view = game.tube_preview.as_ref();
    let frame = game.tube_display_frame();
    let map_on = matches!(
        game.map,
        Some(
            MapMode::Potential
                | MapMode::Electric
                | MapMode::Magnetic
                | MapMode::EnergyFlow
                | MapMode::Current
                | MapMode::Heat
        )
    );
    let metals = !spec.metals.is_empty();
    let (mut negative, mut positive) = (false, false);
    if let Some(f) = &frame {
        for k in &f.kind {
            negative |= *k == level::tube::KIND_NEGATIVE;
            positive |= *k == level::tube::KIND_POSITIVE;
        }
    }
    let (red, blue, grey) = (
        rgb(1.0, 0.45, 0.35),
        rgb(0.4, 0.65, 1.0),
        rgb(0.85, 0.85, 0.85),
    );
    row(
        ui,
        Swatch::Rims([red, blue, grey]),
        "electrode: +, −, 0",
        "An electrode's rim shows its potential: red above 0, blue below, grey at 0 \
         (grounded). Hover an electrode for its potential and circuit.",
    );
    if spec
        .circuits
        .iter()
        .any(|c| c.resistance.is_none() && !c.linked.is_empty())
    {
        row(
            ui,
            Swatch::Line(from_bevy(crate::draw::STRAP)),
            "wire: one floating node",
            "Electrodes joined by a wire are one conductor: they share their potential, and \
             what one takes in the other can give off. No supply holds it: its potential \
             settles where the currents in and out balance.",
        );
    }
    if game.swing.as_ref().is_some_and(|s| s.largest > 0.0) {
        row(
            ui,
            Swatch::Glow(red, blue),
            "potential swinging",
            "An electrode whose potential swings over the run glows with its departure from \
             its mean: warm above it, cool below, the brighter the further. Electrodes that \
             hardly swing do not glow.",
        );
    }
    if game.tube_field_now().is_some() || spec.solenoid_drive.is_some() {
        row(
            ui,
            Swatch::OutIn(rgb(0.95, 0.6, 0.15), rgb(0.2, 0.75, 0.7)),
            "B out of / into board",
            "The solenoid's field, uniform over the board: dots along the edges when it \
             points out of the board, crosses when into it (the magnetic map's colours). A \
             field that swings flips them, and they fade as it passes through zero.",
        );
    }
    let non_metal_dots = match game.tube_particles {
        crate::TubeParticles::Dots => true,
        crate::TubeParticles::Density => !map_on,
        crate::TubeParticles::Hidden => false,
    };
    let metal_dots = !matches!(game.tube_particles, crate::TubeParticles::Hidden);
    let (neg, pos) = (rgb(0.6, 0.85, 1.0), rgb(1.0, 0.7, 0.4));
    let carriers = "Charges in flight, each dot a macroparticle: a long line of many charges \
                    along z.";
    if (metals && metal_dots) || non_metal_dots {
        match (negative, positive) {
            (true, true) => row(ui, Swatch::TwoDots(neg, pos), "charge: −, +", carriers),
            (true, false) => row(
                ui,
                Swatch::Dot(neg),
                if metals {
                    "carriers (−)"
                } else {
                    "electrons"
                },
                carriers,
            ),
            (false, true) => row(
                ui,
                Swatch::Dot(pos),
                if metals { "carriers (+)" } else { "ions" },
                carriers,
            ),
            (false, false) => {}
        }
    } else if matches!(game.tube_particles, crate::TubeParticles::Density) && map_on {
        row(
            ui,
            Swatch::TwoDots(neg, pos),
            "charge density: glow",
            "The particles show as a glow on the map, their charge's density: electrons \
             light blue, ions orange.",
        );
    }
    if metals {
        if metal_dots {
            row(
                ui,
                Swatch::Traced(neg),
                "traced carriers",
                "A sample of the metal's carriers with their paths: straight flights at the \
                 thermal speed, a turn at each collision with an ion core (the faint dots), a \
                 slow drift along the field over many of them. One a contact takes is \
                 replaced by another.",
            );
        }
        row(
            ui,
            Swatch::Chevrons(rgb(0.95, 0.95, 0.9)),
            "drift (the current)",
            "The chevrons move with the metal's carriers as a whole (their mean \
             displacement) and point the way they went over the last time unit: the drift, \
             slow against each carrier's thermal speed. They fade when it is only the thermal \
             motion's noise.",
        );
        if view.is_some_and(|v| v.run.has_field()) && !game.hall.is_empty() {
            let c = |k: [f32; 3]| rgb(k[0], k[1], k[2]);
            row(
                ui,
                Swatch::Arrow(c(crate::draw::HALL_DRIFT)),
                "v: drift of the average carrier",
                "At each metal's middle, its average carrier over the last time unit (the white \
                 dot) with the forces on it, all drawn from it. The thin arrow is the way the \
                 carriers drift there (a velocity, not a force).",
            );
            row(
                ui,
                Swatch::Arrow(c(crate::draw::HALL_MAGNETIC)),
                "qv×B: magnetic push",
                "The magnetic field's push on the drifting carriers, at right angles to their \
                 drift; its way depends on the carriers' sign as well as the current's. It \
                 piles carriers on one side: charge on the metal's sides.",
            );
            row(
                ui,
                Swatch::Arrow(c(crate::draw::HALL_ELECTRIC)),
                "qE: electric push",
                "The electric field's push: along the metal from the supply, and across it from \
                 the charge the magnetic push piled on its sides (the Hall field; side contacts \
                 read it as the Hall voltage). In a long strip the piling goes on until the \
                 push across cancels the magnetic one, and the carriers drift straight along; \
                 in a short wide one the contacts carry that charge off, and the drift turns.",
            );
            row(
                ui,
                Swatch::Arrow(c(crate::draw::HALL_COLLISIONS)),
                "collisions",
                "The collisions with the ion cores take the carriers' momentum away: their mean \
                 push holds the drift steady, so it balances the other two (minus their sum: \
                 the three arrows close). The forces are drawn to the largest one's scale.",
            );
        }
    }
    if !spec.projectiles.is_empty() {
        row(
            ui,
            Swatch::Track(rgb(1.0, 0.95, 0.6)),
            "projectile, detector",
            "A projectile (large dot) with its track so far, and the box it must reach \
             (bright once it has arrived).",
        );
    }
    match game.map {
        Some(MapMode::Electric) => row(
            ui,
            Swatch::Arrow(rgb(1.0, 0.9, 0.6)),
            "E's direction",
            "The electric field's direction at the time shown (the force on a positive \
             charge; electrons are pushed the other way), faded where it is weak. Its \
             strength is the map's colour.",
        ),
        Some(MapMode::EnergyFlow) => row(
            ui,
            Swatch::Arrow(rgb(1.0, 0.9, 0.6)),
            "energy flow's way",
            "The direction of the Poynting vector S = (c/4π) E × B, where the field's energy \
             goes; its size is the map's colour.",
        ),
        Some(MapMode::Current) => row(
            ui,
            Swatch::Arrow(from_bevy(crate::draw::CURRENT_ARROW)),
            "current's way",
            "The direction of the current density j in each cell of the metal: the way \
             positive charge goes (electrons drift against it), faded where it is weak; its \
             size is the map's colour. Averaged over the last few time units and the cells \
             about each one that no cut separates from it (the carriers' thermal noise).",
        ),
        _ => {}
    }
    if game
        .tube_preview
        .as_ref()
        .is_some_and(|v| v.run.leads.leads.iter().any(Option::is_some))
    {
        row(
            ui,
            Swatch::Line(rgb(0.95, 0.95, 0.85)),
            "lead",
            "Each electrode's lead: the wire its circuit current comes along, from the board's \
             edge or the grounded ring enclosing it, drawn as a circuit diagram draws wires: \
             along the axes with right-angle bends, round the electrodes and the metals. Its \
             path is a convention (in the model's world, uniform along z, a real lead would \
             leave along z); the field maps take its current along it.",
        );
    }
    if game
        .tube_preview
        .as_ref()
        .is_some_and(|v| v.run.leads.ground.is_some())
    {
        row(
            ui,
            Swatch::Line(rgb(0.95, 0.95, 0.85)),
            "ground and returns",
            "With the carriers' own field every circuit is a closed loop: each lead's current \
             returns along the board's edge to the ground terminal (the symbol below the lower \
             edge, where the supplies stand). The loop's area is what its current's own \
             magnetic field threads, so its inductance is the loop drawn; the field the loop's \
             changing current induces along its wires acts back on the supply's circuit.",
        );
    }
    if matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow))
        && game.em_snapshot.is_some()
    {
        row(
            ui,
            Swatch::Beads(rgb(0.95, 0.95, 0.85)),
            "lead current",
            "The beads on each lead move with the charge it has carried, steadily on a steady \
             current, to and fro on an oscillating one. The current makes the field around it \
             on the map.",
        );
    }
}

/// The colour `draw.rs` gives a mark, for its swatch.
fn from_bevy(c: bevy::color::Color) -> egui::Color32 {
    let c = c.to_srgba();
    rgb(c.red, c.green, c.blue)
}

/// A particle or beam level's marks that are not objects (charges, magnets, coils, plates
/// and detectors explain themselves on hover): the paths and their colours, the particle
/// and the force on it while animating, the field lines.
fn particle_rows(ui: &mut egui::Ui, game: &Game) {
    let level = &game.editor.level;
    if !game.flights.is_empty() || !game.beams.is_empty() {
        let shot = game.active_shot.min(level.shots.len().saturating_sub(1));
        row(
            ui,
            Swatch::Line(from_bevy(crate::draw::shot_color(level, shot))),
            "path: arrives (verified)",
            "A flight's path in its shot's colour: verified, it arrives in its detector (the \
             detector's colour). With several flights per shot (disturbances) the selected \
             one is brightest.",
        );
        row(
            ui,
            Swatch::Paths([
                rgb(0.95, 0.95, 0.95),
                rgb(1.0, 0.9, 0.2),
                rgb(1.0, 0.45, 0.2),
            ]),
            "computing · unchecked · lost",
            "White: the path is still being computed. Yellow: computed but its verdict not \
             verified (the finer run disagrees or is still running). Orange, ending in a \
             cross: verified, and the particle does not arrive (it hit something, left the \
             arena or ran out of time).",
        );
        if game.animate {
            row(
                ui,
                Swatch::Dot(rgb(1.0, 1.0, 0.6)),
                "the particle now",
                "Where each particle is at the animation time; the one whose energy the panel \
                 shows is ringed in its shot's colour.",
            );
            if game.flights.iter().any(|f| f.preview.is_some()) {
                row(
                    ui,
                    Swatch::Arrow(rgb(1.0, 0.8, 0.2)),
                    "force on it",
                    "The force on the selected particle now: its direction exact, its length \
                     growing as the logarithm of its size (1 + |F|/F₀, F₀ the force that \
                     changes the kinetic energy by the energy unit over a cell).",
                );
            }
        }
    }
    if game.show_field_lines {
        row(
            ui,
            Swatch::Line(rgb(0.95, 0.95, 0.75)),
            "field line (direction)",
            "A line along the static electric field (arrowheads: its way). In this 2D slice of \
             a 3D field the lines show direction only, not strength.",
        );
    }
}
