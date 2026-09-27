//! Side panel: level selection, charge palette, result and verification status, energy
//! bars and controls.

use bevy::prelude::*;
use bevy_egui::{EguiContexts, egui};
use level::ElementKind;
use physics::trajectory::Outcome;
use physics::verify::{Boundary, Status};

use crate::potential::MapMode;
use crate::worker::{PathPoint, Preview};
use crate::{Game, PANEL_WIDTH};

/// Linearly interpolated path point at time `t` (for animation only).
pub fn point_at(p: &Preview, t: f64) -> Option<PathPoint> {
    let path = &p.path;
    let last = *path.last()?;
    if t >= last.t {
        return Some(last);
    }
    let i = path.partition_point(|q| q.t <= t).max(1);
    let (a, b) = (path[i - 1], path[i]);
    let w = if b.t > a.t {
        (t - a.t) / (b.t - a.t)
    } else {
        0.0
    };
    let lerp = |x: f64, y: f64| x + (y - x) * w;
    Some(PathPoint {
        t,
        x: a.x.lerp(b.x, w),
        p: a.p.lerp(b.p, w),
        kinetic: lerp(a.kinetic, b.kinetic),
        potential: lerp(a.potential, b.potential),
        radiated: lerp(a.radiated, b.radiated),
        force: a.force.lerp(b.force, w),
        speed_over_c: lerp(a.speed_over_c, b.speed_over_c),
    })
}

/// Compact number with an SI-style suffix (2e6 → "2M", 1e-6 → "1µ").
pub fn fmt_si(v: f64) -> String {
    let a = v.abs();
    let (scale, suffix) = [
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "k"),
        (1.0, ""),
        (1e-3, "m"),
        (1e-6, "µ"),
        (1e-9, "n"),
    ]
    .into_iter()
    .find(|(s, _)| a >= *s)
    .unwrap_or((1.0, ""));
    let m = v / scale;
    let text = format!("{m:.3}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    format!("{text}{suffix}")
}

/// Parses a number with an optional SI-style suffix ("2M", "1.5k", "3µ" or "3u", "1e6").
pub fn parse_si(text: &str) -> Option<f64> {
    let t = text.trim();
    let (num, scale) = match t.chars().last()? {
        'G' => (&t[..t.len() - 1], 1e9),
        'M' => (&t[..t.len() - 1], 1e6),
        'k' => (&t[..t.len() - 1], 1e3),
        'm' => (&t[..t.len() - 1], 1e-3),
        'u' => (&t[..t.len() - 1], 1e-6),
        'µ' => (&t[..t.len() - 'µ'.len_utf8()], 1e-6),
        'n' => (&t[..t.len() - 1], 1e-9),
        _ => (t, 1.0),
    };
    num.trim().parse::<f64>().ok().map(|v| v * scale)
}

/// What obstacle `i` of a level's scenario is (see `Level::field_at` for the order).
fn obstacle_name(game: &Game, i: usize) -> &'static str {
    let all = game
        .editor
        .level
        .elements
        .iter()
        .chain(&game.editor.placement);
    let charges = all
        .clone()
        .filter(|e| e.kind == ElementKind::Charge)
        .count();
    let magnets = all
        .clone()
        .filter(|e| e.kind == ElementKind::Magnet)
        .count();
    let antennas = all.filter(|e| e.kind == ElementKind::Antenna).count();
    let wires: usize = game
        .editor
        .level
        .coils
        .iter()
        .map(|c| match c {
            level::Coil::Circle { .. } => 1,
            level::Coil::Polygon { vertices, .. } => vertices.len(),
        })
        .sum();
    let electrodes = game.editor.level.electrodes.len();
    let plates = game
        .editor
        .placement
        .iter()
        .filter(|e| e.kind == ElementKind::Plate)
        .count();
    // Obstacle order of `Level::field`: charges, magnets, antennas, coil wires, level
    // electrodes, player plates, metal spheres.
    if i < charges {
        "a charge"
    } else if i < charges + magnets {
        "a magnet"
    } else if i < charges + magnets + antennas {
        "an antenna"
    } else if i < charges + magnets + antennas + wires {
        "a coil wire"
    } else if i < charges + magnets + antennas + wires + electrodes {
        "an electrode"
    } else if i < charges + magnets + antennas + wires + electrodes + plates {
        "one of your plates"
    } else {
        "a metal sphere"
    }
}

/// Shortens `s` to at most `n` characters, with an ellipsis.
fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        s.chars().take(n - 1).collect::<String>() + "…"
    }
}

fn outcome_text(game: &Game, o: Outcome) -> String {
    match o {
        Outcome::Arrived => "reached the detector".into(),
        Outcome::Rejected => {
            "entered the detector outside its acceptance (direction or energy)".into()
        }
        Outcome::Collided(i) => format!("hit {}", obstacle_name(game, i)),
        Outcome::LeftBounds => "left the map".into(),
        Outcome::Timeout => "ran out of time".into(),
        Outcome::Failed(e) => format!("integration failed ({e})"),
        Outcome::SkippedGate(i) => format!("entered the detector without passing gate {}", i + 1),
    }
}

fn boundary_text(b: Boundary) -> &'static str {
    match b {
        Boundary::Obstacle(_) => "a charge",
        Boundary::Bounds => "the map edge",
        Boundary::Detector => "the detector edge",
        Boundary::Acceptance => "the detector's direction/energy window",
        Boundary::TimeLimit => "the time limit",
        Boundary::Gate(_) => "a gate's edge",
        Boundary::GateAcceptance(_) => "a gate's direction/energy window",
    }
}

/// A horizontal bar for an energy relative to T₀, on a scale of −2…2 T₀.
fn energy_bar(ui: &mut egui::Ui, label: &str, value: f64, color: egui::Color32) {
    ui.horizontal(|ui| {
        ui.add_sized([70.0, 16.0], egui::Label::new(label));
        let (rect, _) = ui.allocate_exact_size(egui::vec2(170.0, 14.0), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 2.0, egui::Color32::from_gray(35));
        let zero = rect.center().x;
        #[allow(clippy::cast_possible_truncation)]
        let x = zero + (value.clamp(-2.0, 2.0) as f32) * rect.width() / 4.0;
        let bar = egui::Rect::from_x_y_ranges(zero.min(x)..=zero.max(x), rect.y_range());
        painter.rect_filled(bar, 2.0, color);
        painter.line_segment(
            [
                egui::pos2(zero, rect.top()),
                egui::pos2(zero, rect.bottom()),
            ],
            egui::Stroke::new(1.0, egui::Color32::GRAY),
        );
        ui.label(format!("{value:+.3}"));
    });
}

pub fn panel(mut contexts: EguiContexts, mut game: ResMut<Game>) -> Result {
    let ctx = contexts.ctx_mut()?;
    let mut root = egui::Ui::new(
        ctx.clone(),
        "root".into(),
        egui::UiBuilder::new()
            .layer_id(egui::LayerId::background())
            .max_rect(ctx.viewport_rect()),
    );
    let game = &mut *game;
    let response = egui::Panel::right("side_panel")
        .exact_size(PANEL_WIDTH)
        .resizable(false)
        .show(&mut root, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| contents(ui, game));
        });
    game.panel_left_px = Some(response.response.rect.left() * ctx.pixels_per_point());
    Ok(())
}

fn contents(ui: &mut egui::Ui, game: &mut Game) {
    game.text_focus = false;
    ui.horizontal(|ui| {
        ui.heading("Electromagnetism");
        ui.add_space(20.0);
        let mut sandbox = game.sandbox.active;
        ui.selectable_value(&mut sandbox, false, "Play");
        ui.selectable_value(&mut sandbox, true, "Sandbox");
        if sandbox != game.sandbox.active {
            if sandbox {
                crate::sandbox::enter(game);
            } else {
                game.sandbox.active = false;
            }
        }
    });

    // Level selection.
    let names: Vec<String> = game.levels.iter().map(|l| l.name.clone()).collect();
    let mut selected = game.level_index;
    ui.horizontal(|ui| {
        if ui.button("◀").clicked() {
            game.next_level(-1);
            selected = game.level_index;
        }
        egui::ComboBox::from_id_salt("level")
            .selected_text(truncate(
                &format!("{}. {}", selected + 1, names[selected]),
                30,
            ))
            .width(230.0)
            .show_ui(ui, |ui| {
                for (i, n) in names.iter().enumerate() {
                    ui.selectable_value(&mut selected, i, format!("{}. {n}", i + 1));
                }
            });
        if ui.button("▶").clicked() {
            game.next_level(1);
            selected = game.level_index;
        }
    });
    game.select_level(selected);
    if game.sandbox.active {
        ui.separator();
        crate::sandbox::panel(ui, game);
        ui.separator();
    }
    let level = game.editor.level.clone();
    if !level.description.is_empty() {
        ui.label(egui::RichText::new(&level.description).italics());
    }

    // Shots.
    let n_shots = level.shots.len();
    if n_shots > 1 {
        ui.horizontal_wrapped(|ui| {
            ui.label("Shot:");
            for i in 0..n_shots {
                let mark = if level.has_beams() {
                    // Beams: the shot's verified transmission in every flight.
                    let need = level.shots[i].beam.map_or(1.0, |b| b.transmission);
                    let all: Option<Vec<bool>> = (0..game.beams.len().max(1))
                        .map(|d| {
                            #[allow(clippy::cast_precision_loss)]
                            game.beam_transmission(d, i)
                                .map(|(ok, n)| ok as f64 >= need * n as f64 - 1e-9)
                        })
                        .collect();
                    match all {
                        None => "…",
                        Some(v) if v.iter().all(|&x| x) => "✔",
                        Some(_) => "✖",
                    }
                } else {
                    match game.verdict(i) {
                        Some((Status::Verified, Outcome::Arrived)) => "✔",
                        Some((Status::Verified, _)) => "✖",
                        Some(_) => "⚠",
                        None => "…",
                    }
                };
                let c = crate::draw::shot_color(i).to_srgba();
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let color = egui::Color32::from_rgb(
                    (c.red * 255.0) as u8,
                    (c.green * 255.0) as u8,
                    (c.blue * 255.0) as u8,
                );
                let text = egui::RichText::new(format!("{} {mark}", i + 1)).color(color);
                if ui
                    .selectable_label(game.active_shot == i, text)
                    .on_hover_text("Switch shot ([ and ])")
                    .clicked()
                {
                    game.select_shot(i);
                }
            }
            ui.checkbox(&mut game.show_all_shots, "show all (H)");
        });
        ui.label(
            egui::RichText::new("One setup must deliver every shot to its own detector.").small(),
        );
    }
    // Disturbances.
    if !level.disturbances.is_empty() {
        ui.horizontal_wrapped(|ui| {
            ui.label("Disturbance:");
            for d in 0..level.disturbances.len() {
                let mark = match game.disturbance_verdict(d) {
                    Some((Status::Verified, Outcome::Arrived)) => "✔",
                    Some((Status::Verified, _)) => "✖",
                    Some(_) => "⚠",
                    None => "…",
                };
                let name = &level.disturbances[d].name;
                let label = if name.is_empty() {
                    format!("{} {mark}", d + 1)
                } else {
                    format!("{} {name} {mark}", d + 1)
                };
                if ui
                    .selectable_label(game.active_disturbance == d, label)
                    .on_hover_text(disturbance_text(&level.disturbances[d]))
                    .clicked()
                {
                    game.active_disturbance = d;
                }
            }
        });
        ui.label(
            egui::RichText::new(
                "Outside fields: the setup must work under every disturbance. All flights \
                 of a shot are drawn; the selected one is brightest.",
            )
            .small(),
        );
    }
    let shot_index = game.active_shot.min(n_shots.saturating_sub(1));
    let shot = level.shots[shot_index];

    // Launch of the active shot.
    let kin = physics::dynamics::Kinematics::new(shot.particle.mass, level.c());
    let p0 = level.launch_momentum(shot_index);
    let gamma0 = kin.gamma(p0);
    let v0 = kin.velocity(p0).length();
    let t0 = shot.launch.kinetic_energy;
    ui.label(match level.physics.c {
        Some(c) => format!(
            "Launch: T₀ = {}, v₀ = {:.3} c, γ₀ = {gamma0:.4}   (c = {c})",
            energy_text(t0),
            v0 / c
        ),
        None => format!("Launch: T₀ = {}, v₀ = {v0:.3} (Newtonian)", energy_text(t0)),
    });
    ui.label(if shot.particle.moment == 0.0 {
        format!(
            "Particle: q = {}, m = {}",
            fmt_si(shot.particle.charge),
            fmt_si(shot.particle.mass)
        )
    } else {
        format!(
            "Particle: q = {}, m = {}, magnetic moment along {} of {} (spin {})",
            fmt_si(shot.particle.charge),
            fmt_si(shot.particle.mass),
            if shot.particle.moment > 0.0 {
                "+z"
            } else {
                "−z"
            },
            fmt_si(shot.particle.moment.abs()),
            if shot.particle.moment > 0.0 {
                "up"
            } else {
                "down"
            }
        )
    });
    egui::CollapsingHeader::new("Physics model and its limits")
        .id_salt("model_notes")
        // Opened for developer captures with EM_MODEL=1 (see `dev_capture`).
        .default_open(std::env::var("EM_MODEL").is_ok_and(|v| v == "1"))
        .show(ui, |ui| {
            for n in level.model_notes(&game.editor.placement) {
                ui.horizontal_wrapped(|ui| {
                    if n.exact {
                        ui.colored_label(egui::Color32::from_rgb(110, 210, 120), "exact")
                    } else {
                        ui.colored_label(egui::Color32::from_rgb(230, 180, 70), "approx.")
                    }
                    .on_hover_text(if n.exact {
                        "Exact within classical electrodynamics, up to the integration accuracy"
                    } else {
                        "An approximation or omission; where it can matter, its size is measured"
                    });
                    ui.label(egui::RichText::new(n.text).small());
                });
            }
            ui.label(
                egui::RichText::new("Details and tests: PHYSICS.md.")
                    .small()
                    .weak(),
            );
        });
    ui.separator();

    // Palette.
    let limits = &level.limits;
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Your elements").strong());
        let mut hard = game.editor.continuous();
        if ui
            .checkbox(&mut hard, "hardcore")
            .on_hover_text(
                "Continuous values: any magnitude or frequency between the level's smallest \
                 and largest listed value, any antenna orientation, set with sliders",
            )
            .changed()
        {
            game.editor.set_continuous(hard);
        }
    });
    let kind_name = |k: ElementKind| match k {
        ElementKind::Charge => "charges",
        ElementKind::Magnet => "magnets",
        ElementKind::Antenna => "antennas",
        ElementKind::Plate => "plates",
        ElementKind::Supply => "power supplies",
    };
    let allowed: Vec<ElementKind> = crate::editor::KINDS
        .into_iter()
        .filter(|&k| crate::editor::max_of(&level, k) > 0)
        .collect();
    ui.horizontal_wrapped(|ui| {
        for &k in &allowed {
            ui.label(format!(
                "{} {} of {} left",
                kind_name(k),
                game.editor.left(k),
                crate::editor::max_of(&level, k)
            ));
        }
    });
    if allowed.len() > 1 {
        ui.horizontal(|ui| {
            ui.label("Place:");
            let mut kind = game.editor.kind;
            for &k in &allowed {
                ui.selectable_value(&mut kind, k, kind_name(k).trim_end_matches('s'))
                    .on_hover_text("Cycle with M");
            }
            game.editor.set_kind(kind);
        });
    }
    if allowed.is_empty() {
        // Nothing to place (e.g. only power supplies to set).
    } else if crate::editor::is_signed(game.editor.kind) {
        plate_palette(ui, game, &level);
    } else {
        ui.horizontal_wrapped(|ui| {
            let kind = game.editor.kind;
            let both_signs =
                kind != ElementKind::Charge || (limits.allow_positive && limits.allow_negative);
            let sign = match (kind, game.editor.positive) {
                (ElementKind::Charge, true) => "+",
                (ElementKind::Charge, false) => "−",
                (ElementKind::Magnet, true) => "out (+z)",
                (ElementKind::Magnet, false) => "in (−z)",
                (ElementKind::Antenna, true) => "phase 0°",
                (ElementKind::Antenna, false) => "phase 180°",
                (ElementKind::Plate | ElementKind::Supply, _) => "",
            };
            let hover = match kind {
                ElementKind::Charge => "Flip sign (S)",
                ElementKind::Magnet => {
                    "Flip orientation (S): moment out of the plane (+z) or into it (−z)"
                }
                ElementKind::Antenna => "Flip phase (S): opposite phase of the RF generator",
                ElementKind::Plate | ElementKind::Supply => "",
            };
            ui.label(match kind {
                ElementKind::Charge => "New charge:",
                ElementKind::Magnet => "New magnet μ:",
                ElementKind::Antenna => "New antenna p₀:",
                ElementKind::Plate | ElementKind::Supply => "",
            });
            if both_signs {
                if ui.button(sign).on_hover_text(hover).clicked() {
                    game.editor.positive = !game.editor.positive;
                }
            } else {
                ui.label(sign);
            }
            if !game.editor.continuous() {
                for (i, m) in crate::editor::magnitudes(&level, kind).iter().enumerate() {
                    ui.selectable_value(&mut game.editor.magnitude_index, i, fmt_si(*m));
                }
            }
        });
    }
    supplies(ui, game, &level);
    if game.editor.continuous() {
        hardcore_controls(ui, game, &level);
    } else if game.editor.kind == ElementKind::Antenna {
        ui.horizontal(|ui| {
            ui.label("Orientation:");
            for a in level::ANTENNA_ANGLES {
                ui.selectable_value(&mut game.editor.angle_deg, a, format!("{a:.0}°"))
                    .on_hover_text("Rotate with R");
            }
        });
        let omegas = &level.limits.antenna_omegas;
        if omegas.is_empty() {
            ui.label(
                egui::RichText::new(format!(
                    "Antennas oscillate at ω = {} (the level's RF generator).",
                    fmt_si(level.physics.rf_omega)
                ))
                .small(),
            );
        } else {
            ui.horizontal_wrapped(|ui| {
                ui.label("Frequency ω:");
                for (i, w) in omegas.iter().enumerate() {
                    ui.selectable_value(&mut game.editor.omega_index, i, fmt_si(*w))
                        .on_hover_text("Cycle with W (Shift+W: down)");
                }
            });
            ui.label(egui::RichText::new("All antennas start in phase at t = 0.").small());
        }
    }
    ui.horizontal(|ui| {
        ui.label("Grid:");
        for f in 1..=4u32 {
            let current = game.editor.refinement() == f;
            if ui.selectable_label(current, format!("{f}×")).clicked() {
                game.editor.set_refinement(f);
            }
        }
        if ui.button("Clear").clicked() {
            game.editor.clear();
        }
    });
    if let Some(msg) = &game.editor.message {
        ui.colored_label(egui::Color32::from_rgb(255, 170, 80), msg);
    }
    ui.separator();

    if level.has_beams() {
        beam_result(ui, game, &level);
    } else {
        // Result.
        ui.label(egui::RichText::new("Result").strong());
        if game.solved() {
            ui.label(
                egui::RichText::new(if n_shots > 1 {
                    "✔ SOLVED: every shot arrives (verified)"
                } else {
                    "✔ SOLVED (verified)"
                })
                .color(egui::Color32::from_rgb(90, 240, 110))
                .size(18.0),
            );
        }
        let view = game
            .flights
            .get(game.active_flight())
            .cloned()
            .unwrap_or_default();
        match (&view.preview, view.verdict) {
            (None, _) => {
                ui.label("Computing…");
            }
            (Some(p), verdict) => {
                let who = if n_shots > 1 {
                    format!("Shot {}", shot_index + 1)
                } else {
                    "The particle".into()
                };
                ui.label(format!(
                    "{who} {} after t = {:.3}.",
                    outcome_text(game, p.outcome),
                    p.flight_time
                ));
                // Detector acceptance: what is allowed, and how the particle arrived.
                if let Some(acc) = level.shots[shot_index].detector.acceptance
                    && let Some(last) = p.path.last()
                {
                    let v = last.p;
                    let dir = v.y.atan2(v.x).to_degrees();
                    let mut parts = Vec::new();
                    if let Some([axis, half]) = acc.direction {
                        parts.push(format!(
                            "direction {axis:.0}° ± {half:.1}° (arrives at {dir:.1}°)"
                        ));
                    }
                    if let Some([lo, hi]) = acc.kinetic {
                        parts.push(format!(
                            "energy {lo:.3}–{hi:.3} (arrives with {:.3})",
                            last.kinetic
                        ));
                    }
                    ui.label(
                        egui::RichText::new(format!("Detector accepts: {}", parts.join("; ")))
                            .small(),
                    );
                }
                match verdict {
                    None => {
                        ui.label("Verifying at 100× tighter tolerance…");
                    }
                    Some((Status::Verified, _)) => {
                        ui.colored_label(egui::Color32::LIGHT_GRAY, "Verified.");
                    }
                    Some((
                        Status::SmallMargin {
                            boundary,
                            margin,
                            error_estimate,
                        },
                        _,
                    )) => {
                        ui.colored_label(
                        egui::Color32::YELLOW,
                        format!(
                            "⚠ Marginal: passes {} by {margin:.2e} cells, numerical error ≈ {error_estimate:.1e}. Too close to call.",
                            boundary_text(boundary)
                        ),
                    );
                    }
                    Some((Status::OutcomeMismatch { .. }, _)) => {
                        ui.colored_label(
                            egui::Color32::YELLOW,
                            "⚠ Marginal: the outcome depends on numerical precision.",
                        );
                    }
                    Some((Status::Failed, _)) => {
                        ui.colored_label(egui::Color32::RED, "Integration failed.");
                    }
                }
                egui::Grid::new("diag").num_columns(2).show(ui, |ui| {
                if level.physics.c.is_some() {
                    ui.label("max v/c");
                    ui.label(format!("{:.4}", p.max_speed_over_c));
                    ui.end_row();
                }
                ui.label("energy error |ΔW|/T₀");
                if p.energy_rel_error.is_nan() {
                    // Time-dependent fields do work on the particle.
                    ui.label("n/a (fields vary in time)");
                } else {
                    ui.label(format!("{:.1e}", p.energy_rel_error));
                }
                ui.end_row();
                if p.radiation_reaction {
                    ui.label("radiated / T₀ (included)");
                    ui.label(format!("{:.3e}", p.radiation_loss_fraction));
                    ui.end_row();
                    ui.label("max |F_rad| / |F_Lorentz|");
                    let text = format!("{:.1e}", p.reaction_ratio_max);
                    if p.reaction_ratio_max > 0.05 {
                        ui.colored_label(
                            egui::Color32::YELLOW,
                            text + "  LL approximation strained",
                        );
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                } else if level.physics.c.is_some() {
                    ui.label("radiated / T₀ (neglected)");
                    let text = format!("{:.1e}", p.radiated_fraction);
                    if p.radiated_fraction > 1e-10 {
                        ui.colored_label(egui::Color32::YELLOW, text + "  not negligible!");
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                }
                if p.image_force_bound > 0.0 {
                    ui.label("electrode image force (neglected)")
                        .on_hover_text(
                            "Bound on the force of the charge the particle induces on the                              electrodes, relative to the force that matters (PHYSICS.md §2.7).                              It grows as the particle passes closer to metal.",
                        );
                    let text = format!("≤ {:.1e}", p.image_force_bound);
                    if p.image_force_bound > level::IMAGE_FORCE_LIMIT {
                        ui.colored_label(
                            egui::Color32::YELLOW,
                            text + "  not negligible: keep away from the metal",
                        );
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                }
            });
            }
        }
        ui.separator();

        // Energy bars at the animated point of the active shot.
        ui.label(egui::RichText::new("Energy along the flight (units of T₀)").strong());
        ui.horizontal(|ui| {
            ui.checkbox(&mut game.animate, "Animate (A)");
            ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"));
        });
        if let Some(p) = &view.preview
            && let Some(pt) = point_at(
                p,
                if game.animate {
                    game.anim_time
                } else {
                    p.flight_time
                },
            )
        {
            energy_bar(
                ui,
                "kinetic",
                pt.kinetic / t0,
                egui::Color32::from_rgb(240, 200, 60),
            );
            energy_bar(
                ui,
                "potential",
                pt.potential / t0,
                egui::Color32::from_rgb(200, 90, 230),
            );
            if p.radiation_reaction {
                energy_bar(
                    ui,
                    "radiated",
                    pt.radiated / t0,
                    egui::Color32::from_rgb(90, 200, 255),
                );
                energy_bar(
                    ui,
                    "total − T₀",
                    (pt.kinetic + pt.potential + pt.radiated - t0) / t0,
                    egui::Color32::from_rgb(120, 220, 120),
                );
                ui.label(
                    egui::RichText::new(
                        "total = kinetic + potential + radiated: the energy the particle radiates \
                     is carried away by its field.",
                    )
                    .small(),
                );
            } else {
                energy_bar(
                    ui,
                    "total − T₀",
                    (pt.kinetic + pt.potential - t0) / t0,
                    egui::Color32::from_rgb(120, 220, 120),
                );
            }
            if level.physics.c.is_some() {
                ui.label(format!("t = {:.2}   v/c = {:.4}", pt.t, pt.speed_over_c));
            } else {
                ui.label(format!("t = {:.2}", pt.t));
            }
        }
        ui.separator();
    }
    ui.label(egui::RichText::new("View").strong());
    ui.horizontal_wrapped(|ui| {
        ui.label("Map (V):");
        ui.selectable_value(&mut game.map, Some(MapMode::Potential), "potential");
        ui.selectable_value(&mut game.map, Some(MapMode::Magnetic), "magnetic B");
        if crate::radiation::waves_available(&level) {
            ui.selectable_value(&mut game.map, Some(MapMode::Waves), "waves");
        }
        if crate::radiation::particle_field_available(&level) {
            ui.selectable_value(
                &mut game.map,
                Some(MapMode::ParticleField),
                "particle field",
            );
        }
        ui.selectable_value(&mut game.map, Some(MapMode::Total), "total")
            .on_hover_text("Everything at once: level sources, antennas, waves and the particle");
        ui.selectable_value(&mut game.map, None, "off");
    });
    if matches!(game.map, Some(MapMode::ParticleField | MapMode::Total)) {
        let (shot, d) = level.flight_of(game.active_flight());
        let which = if level.disturbances.is_empty() {
            format!("shot {}", shot + 1)
        } else {
            format!("shot {}, disturbance {}", shot + 1, d + 1)
        };
        let multi = level.flight_count() > 1;
        ui.label(
            egui::RichText::new(format!(
                "Field of the particle of {which}{}; only its flight is drawn.",
                if multi {
                    " (choose with the shot/disturbance tabs or [ ])"
                } else {
                    ""
                }
            ))
            .small(),
        );
    }
    if matches!(
        game.map,
        Some(MapMode::Waves | MapMode::ParticleField | MapMode::Total)
    ) {
        ui.horizontal(|ui| {
            use crate::radiation::FieldQuantity;
            ui.label("Colour:");
            ui.selectable_value(&mut game.field_quantity, FieldQuantity::Bz, "B_z");
            ui.selectable_value(&mut game.field_quantity, FieldQuantity::E, "|E|");
        });
        ui.add(
            egui::Slider::new(&mut game.field_range_decades, 1.0..=14.0).text("range (decades)"),
        )
        .on_hover_text("How many decades below the strongest field are still visible");
        ui.horizontal(|ui| {
            ui.checkbox(&mut game.show_field_arrows, "E arrows");
            if game.map == Some(MapMode::ParticleField) {
                ui.checkbox(&mut game.radiation_only, "radiation part only")
                    .on_hover_text("Only the acceleration term of the field (falls as 1/R)");
            }
        });
    }
    let legend = match game.map {
        Some(MapMode::Potential) => {
            "Red: uphill for this shot's particle, blue: downhill; contours every T₀/4. \
             Dark: forbidden by energy conservation for every particle shown (all shots \
             with \"show all\"; a beam counts with its most energetic particle). Exact, \
             also with magnets; with interacting beam particles only a guide, since they \
             exchange energy."
        }
        Some(MapMode::Magnetic) if shot.particle.charge == 0.0 && shot.particle.moment != 0.0 => {
            "B perpendicular to the plane. Orange: out of the plane, teal: into it. \
             Value 1 = field in which this particle's magnetic energy |m B| equals T₀; \
             contours every 0.25."
        }
        Some(MapMode::Magnetic) => {
            "B perpendicular to the plane. Orange: out of the plane, teal: into it. \
             Value 1 = field in which this particle circles with a 5-cell radius; \
             contours every 0.25."
        }
        Some(MapMode::Waves) => {
            "Fields of the antennas and waves at the animation time (exact retarded fields). \
             Colour: B perpendicular to the plane (orange out, blue in), on a logarithmic \
             scale; arrows: E. Near an antenna the field is quasi-static; further out the \
             radiation travels outwards at c."
        }
        Some(MapMode::Total) => {
            "The total field at the animation time: level charges, magnets and coils, \
             antennas, waves and disturbances, and the particle's own field (retarded, \
             Liénard–Wiechert). The particle's field is usually far weaker than the \
             electrodes'; raise the range to see it. Colour: B_z or |E| on a logarithmic \
             scale; arrows: E."
        }
        Some(MapMode::ParticleField) => {
            "The field of the particle itself (Liénard–Wiechert, exact) at the animation \
             time. Every change of velocity sends out a radiation pulse at c; the energy it \
             carries is what the particle loses (radiation reaction, when included). \
             Colour: B perpendicular to the plane, logarithmic; arrows: E. Before launch the \
             particle is taken to move uniformly."
        }
        None => "",
    };
    if !legend.is_empty() {
        ui.label(egui::RichText::new(legend).small());
    }
    ui.checkbox(&mut game.show_field_lines, "Electric field lines (F)");
    ui.add(egui::Slider::new(&mut game.field_line_spacing, 0.5..=4.0).text("spacing (cells)"));
    ui.add(egui::Slider::new(&mut game.field_line_opacity, 0.05..=1.0).text("opacity"));
    ui.label(
        egui::RichText::new(
            "In this 2D slice of a 3D field, lines show direction only, not strength.",
        )
        .small(),
    );
    ui.separator();

    ui.collapsing("Controls", |ui| {
        ui.label("Mouse: left click place, right click remove, wheel changes magnitude.");
        ui.label("Arrows move the cursor (Shift: ×5). Space/Enter place, Del/X remove.");
        ui.label("S flip sign/phase, Q/E change magnitude, C clear.");
        ui.label("1–4 grid refinement, [ ] switch shot, H show all shots.");
        ui.label("N/P next/previous level, V map, F field lines, A animation.");
        ui.label("R rotate antenna (Shift+R back), W antenna frequency, M cycles kinds.");
        ui.label("Drag your elements with the mouse, or G to grab / drop and Esc to cancel.");
        ui.label("Hardcore (checkbox): sliders instead of fixed values; Q/E and W step ×1.1.");
    });
}

/// Beam levels: verified transmission of every beam shot in every flight against its
/// requirement, and the diagnostics of the selected flight.
fn beam_result(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    ui.label(egui::RichText::new("Result").strong());
    if game.solved() {
        ui.label(
            egui::RichText::new("✔ SOLVED: every beam delivers its share (verified)")
                .color(egui::Color32::from_rgb(90, 240, 110))
                .size(18.0),
        );
    }
    let flights = game.beams.len().max(1);
    for d in 0..flights {
        if flights > 1 {
            ui.label(egui::RichText::new(format!("Disturbance {}", d + 1)).small());
        }
        for (s, shot) in level.shots.iter().enumerate() {
            let need = shot.beam.map_or(1.0, |b| b.transmission);
            let text = match game.beam_transmission(d, s) {
                None => format!("Shot {}: computing…", s + 1),
                Some((ok, n)) => {
                    #[allow(clippy::cast_precision_loss)]
                    let enough = ok as f64 >= need * n as f64 - 1e-9;
                    format!(
                        "Shot {}: {ok} of {n} arrive, verified (need {:.0} %) {}",
                        s + 1,
                        need * 100.0,
                        if enough { "✔" } else { "✘" }
                    )
                }
            };
            ui.colored_label(shot_color32(s), text);
        }
    }
    let d = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    if let Some(p) = game.beams.get(d).and_then(|b| b.preview.as_ref()) {
        let lost = p
            .outcomes
            .iter()
            .filter(|o| **o != Outcome::Arrived)
            .count();
        let energy = if p.energy_rel_error.is_nan() {
            // Retarded interaction: the particles exchange energy with the field.
            format!(
                "largest radiated share of a particle {:.1e} (neglected)",
                p.radiated_max
            )
        } else {
            format!("energy drift of the whole beam {:.1e}", p.energy_rel_error)
        };
        ui.label(
            egui::RichText::new(format!(
                "{} particles{}; {lost} lost in the preview; {energy}.",
                p.paths.len(),
                match (level.physics.beam_interaction, level.physics.c) {
                    (false, _) => ", not interacting",
                    (true, None) => ", interacting",
                    (true, Some(_)) => ", interacting (retarded fields)",
                },
            ))
            .small(),
        );
    }
    ui.separator();
}

/// The drawing colour of a shot, for egui.
fn shot_color32(shot: usize) -> egui::Color32 {
    let c = crate::draw::shot_color(shot).to_srgba();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let b = |v: f32| (v.clamp(0.0, 1.0) * 255.0) as u8;
    egui::Color32::from_rgb(b(c.red), b(c.green), b(c.blue))
}

/// A launch energy: plain for ordinary values, with SI prefixes for the large ones of
/// scaled beams.
fn energy_text(t: f64) -> String {
    if t.abs() >= 1e4 {
        fmt_si(t)
    } else {
        format!("{t:.3}")
    }
}

/// A potential for display: "grounded" or its value.
pub fn fmt_potential(v: f64) -> String {
    if v == 0.0 {
        "grounded".into()
    } else {
        format!("V = {}", fmt_si(v))
    }
}

/// Palette row for new plates: potential and orientation.
fn plate_palette(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    ui.horizontal_wrapped(|ui| {
        ui.label("New plate:");
        if !game.editor.continuous() {
            for (i, v) in level.limits.plate_voltages.iter().enumerate() {
                ui.selectable_value(&mut game.editor.magnitude_index, i, fmt_potential(*v))
                    .on_hover_text("Q/E or the wheel: next potential; S: the opposite one");
            }
        }
    });
    if !game.editor.continuous() {
        ui.horizontal(|ui| {
            ui.label("Orientation:");
            for (a, name) in level::PLATE_ANGLES.into_iter().zip(["along x", "along y"]) {
                ui.selectable_value(&mut game.editor.plate_angle_deg, a, name)
                    .on_hover_text("Rotate with R");
            }
        });
    }
    let s = level.limits.plate;
    ui.label(
        egui::RichText::new(format!(
            "Plates are {} × {} cells in the plane, {} high, and keep {} cell from other \
             electrodes.",
            s.length,
            s.thickness,
            s.height,
            level::PLATE_CLEARANCE
        ))
        .small(),
    );
}

/// The power supplies of the level's tunable electrodes: off (the electrode keeps its
/// own bias) or one of the listed potentials. Clicking an electrode on the map steps it.
fn supplies(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    let tunable: Vec<(usize, &level::Electrode)> = level
        .electrodes
        .iter()
        .enumerate()
        .filter(|(_, e)| e.tunable)
        .collect();
    if tunable.is_empty() {
        return;
    }
    ui.label(egui::RichText::new("Power supplies").strong())
        .on_hover_text(
            "Click a tunable electrode (yellow frame) to switch its supply on or to the next \
             potential; S flips it, right click switches it off.",
        );
    let list = level.limits.supply_voltages.clone();
    for (i, e) in tunable {
        let own = match e.bias {
            level::ConductorBias::Grounded => "grounded".to_string(),
            level::ConductorBias::Potential(v) => fmt_potential(v),
            level::ConductorBias::Charge(q) => format!("charge {}", fmt_si(q)),
        };
        let current = game.editor.supply(e.center);
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Electrode {}:", i + 1));
            let mut choice = current;
            if game.editor.continuous() {
                let mut on = current.is_some();
                ui.checkbox(&mut on, "on");
                if let Some((lo, hi)) = level::value_range(&list) {
                    let mut v = current.unwrap_or(lo.max(0.0_f64.min(hi)));
                    let r = ui.add_enabled(
                        on,
                        egui::Slider::new(&mut v, lo..=hi)
                            .custom_formatter(|v, _| fmt_si(v))
                            .custom_parser(parse_si),
                    );
                    game.text_focus |= r.has_focus();
                    choice = on.then_some(v);
                }
            } else {
                ui.selectable_value(&mut choice, None, format!("off ({own})"));
                for v in &list {
                    ui.selectable_value(&mut choice, Some(*v), fmt_potential(*v));
                }
            }
            game.editor.set_supply(e.center, choice);
        });
    }
}

/// Hardcore sliders: the element under the cursor if there is one (edited live), else
/// the values of new elements.
fn hardcore_controls(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    use level::{Element, value_range};
    let target = game.editor.element_at_cursor();
    let mut e = match target {
        Some(i) => game.editor.placement[i],
        None => Element {
            node: game.editor.cursor,
            kind: game.editor.kind,
            value: game.editor.selected_value(),
            angle_deg: game.editor.angle_deg,
            omega: game.editor.selected_omega(),
        },
    };
    ui.label(
        egui::RichText::new(match target {
            Some(_) => "Sliders edit the element under the cursor.",
            None => "Sliders set new elements (move the cursor onto one of yours to edit it).",
        })
        .small(),
    );
    let mut focus = false;
    let name = match e.kind {
        ElementKind::Charge => "|Q|",
        ElementKind::Magnet => "|μ|",
        ElementKind::Antenna => "|p₀|",
        ElementKind::Plate | ElementKind::Supply => "V",
    };
    let signed = crate::editor::is_signed(e.kind);
    if signed {
        // Signed potentials on a linear scale (the range may include 0).
        if let Some((lo, hi)) = value_range(crate::editor::magnitudes(level, e.kind)) {
            ui.horizontal(|ui| {
                ui.label(name);
                if lo < hi {
                    let r = ui.add(
                        egui::Slider::new(&mut e.value, lo..=hi)
                            .custom_formatter(|v, _| fmt_si(v))
                            .custom_parser(parse_si),
                    );
                    focus |= r.has_focus();
                } else {
                    ui.label(fmt_si(lo));
                }
            });
        }
        if e.kind == ElementKind::Plate {
            ui.horizontal(|ui| {
                ui.label("angle");
                let r = ui.add(
                    egui::Slider::new(&mut e.angle_deg, 0.0..=180.0)
                        .suffix("°")
                        .step_by(0.5),
                );
                focus |= r.has_focus();
            });
        }
    } else if let Some((lo, hi)) = value_range(crate::editor::magnitudes(level, e.kind)) {
        let mut m = e.value.abs();
        ui.horizontal(|ui| {
            ui.label(name);
            if lo < hi {
                let r = ui.add(
                    egui::Slider::new(&mut m, lo..=hi)
                        .logarithmic(game.log_magnitude)
                        .custom_formatter(|v, _| fmt_si(v))
                        .custom_parser(parse_si),
                );
                focus |= r.has_focus();
                ui.checkbox(&mut game.log_magnitude, "log");
            } else {
                ui.label(fmt_si(lo));
            }
        });
        e.value = e.value.signum() * m;
    }
    if e.kind == ElementKind::Antenna {
        ui.horizontal(|ui| {
            ui.label("angle");
            let r = ui.add(
                egui::Slider::new(&mut e.angle_deg, 0.0..=360.0)
                    .suffix("°")
                    .step_by(0.5),
            );
            focus |= r.has_focus();
        });
        let omegas = &level.limits.antenna_omegas;
        match (value_range(omegas), e.omega.as_mut()) {
            (Some((lo, hi)), Some(w)) if lo < hi => {
                ui.horizontal(|ui| {
                    ui.label("ω");
                    let r = ui.add(
                        egui::Slider::new(w, lo..=hi)
                            .logarithmic(game.log_omega)
                            .custom_formatter(|v, _| fmt_si(v))
                            .custom_parser(parse_si),
                    );
                    focus |= r.has_focus();
                    ui.checkbox(&mut game.log_omega, "log");
                });
            }
            _ => {
                ui.label(
                    egui::RichText::new(format!(
                        "ω = {}",
                        fmt_si(e.omega.unwrap_or(level.physics.rf_omega))
                    ))
                    .small(),
                );
            }
        }
    }
    game.text_focus |= focus;
    match target {
        Some(i) => game.editor.set_element(i, e),
        None if signed => {
            let k = crate::editor::kind_index(e.kind);
            game.editor.continuous_magnitude[k] = e.value;
            game.editor.plate_angle_deg = e.angle_deg;
        }
        None => {
            let k = crate::editor::kind_index(e.kind);
            game.editor.continuous_magnitude[k] = e.value.abs();
            game.editor.angle_deg = e.angle_deg;
            if let Some(w) = e.omega {
                game.editor.continuous_omega = w;
            }
        }
    }
}

/// One-line description of a disturbance.
pub fn disturbance_text(d: &level::Disturbance) -> String {
    let mut parts = Vec::new();
    if d.e != [0.0, 0.0] {
        parts.push(format!(
            "stray E = ({}, {})",
            fmt_si(d.e[0]),
            fmt_si(d.e[1])
        ));
    }
    if d.bz != 0.0 {
        parts.push(format!("stray B_z = {}", fmt_si(d.bz)));
    }
    for w in &d.waves {
        parts.push(if w.omega == 0.0 {
            format!(
                "static field {} at {:.0}°",
                fmt_si(w.amplitude),
                w.direction_deg + 90.0
            )
        } else {
            format!(
                "wave E₀ = {}, ω = {}, towards {:.0}°, phase {:.0}°",
                fmt_si(w.amplitude),
                fmt_si(w.omega),
                w.direction_deg,
                w.phase_deg
            )
        });
    }
    if parts.is_empty() {
        "undisturbed".into()
    } else {
        parts.join("; ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn si_round_trip() {
        for v in [2e6, 1.5e3, 1e-6, -4e6, 0.5, 12.0] {
            let back = parse_si(&fmt_si(v)).unwrap();
            assert!(
                (back - v).abs() <= 1e-12 * v.abs(),
                "{v} -> {} -> {back}",
                fmt_si(v)
            );
        }
        assert_eq!(parse_si("3u").map(f64::to_bits), Some(3e-6f64.to_bits()));
        assert_eq!(parse_si("1e6").map(f64::to_bits), Some(1e6f64.to_bits()));
    }
}
