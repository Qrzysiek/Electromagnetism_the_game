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

/// Decimals that resolve a direction tolerance of `half` degrees, for it and the arrival
/// angle: one from 1° up, one more per decade below (0.05° shows as 0.050°).
fn angle_decimals(half: f64) -> usize {
    // At most 6: finer tolerances are not met by any flight anyway.
    let (mut digits, mut step) = (1, 1.0);
    while half > 0.0 && half < step && digits < 6 {
        step /= 10.0;
        digits += 1;
    }
    digits
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
        // Bounded width: a flight that falls onto a charge reaches 1e5 T₀, and a long
        // number widened the panel (and with it moved the board under the cursor).
        if value != 0.0 && !(1e-3..1e3).contains(&value.abs()) {
            ui.label(format!("{value:+.2e}"));
        } else {
            ui.label(format!("{value:+.3}"));
        }
    });
}

pub fn panel(
    mut contexts: EguiContexts,
    mut game: ResMut<Game>,
    radiation: Res<crate::radiation::RadiationView>,
) -> Result {
    let ctx = contexts.ctx_mut()?;
    let mut root = egui::Ui::new(
        ctx.clone(),
        "root".into(),
        egui::UiBuilder::new()
            .layer_id(egui::LayerId::background())
            .max_rect(ctx.viewport_rect()),
    );
    let game = &mut *game;
    egui::Panel::right("side_panel")
        .exact_size(PANEL_WIDTH)
        .resizable(false)
        .show(&mut root, |ui| {
            let mut area = egui::ScrollArea::vertical();
            // Screenshots (EM_CAPTURE) can scroll the panel: EM_SCROLL=<points>.
            if let Some(y) = std::env::var("EM_SCROLL")
                .ok()
                .and_then(|v| v.parse::<f32>().ok())
            {
                area = area.vertical_scroll_offset(y);
            }
            let available = ui.available_width();
            let out = area.show(ui, |ui| contents(ui, game, &radiation));
            // A row wider than the panel makes egui draw a stray full-height line: report
            // it in capture runs (the audit looks for it).
            let used = out.content_size.x;
            if std::env::var("EM_CAPTURE").is_ok() && used > available + 0.5 {
                eprintln!("panel overflow: content {used:.1} wider than {available:.1}");
            }
        });
    // The panel's nominal edge, not its measured rect: an overflowing row must not move
    // the board (and with it the node under the cursor).
    let left = ctx.viewport_rect().right() - PANEL_WIDTH;
    game.panel_left_px = Some(left * ctx.pixels_per_point());
    Ok(())
}

fn contents(ui: &mut egui::Ui, game: &mut Game, radiation: &crate::radiation::RadiationView) {
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
                // Grouped by the curriculum: a heading for every arc, the tier before
                // its first level.
                let mut last: Option<&crate::curriculum::Place> = None;
                for (i, n) in names.iter().enumerate() {
                    let place = game.places.get(i).and_then(Option::as_ref);
                    if let Some(p) = place {
                        if last.is_none_or(|l| l.arc != p.arc) {
                            ui.label(egui::RichText::new(p.arc_title()).strong());
                        }
                        if last.is_none_or(|l| l.arc != p.arc || l.tier != p.tier) {
                            ui.label(egui::RichText::new(format!("  {}", p.tier)).weak().small());
                        }
                    } else if last.is_some() {
                        ui.label(egui::RichText::new("Custom levels").strong());
                    }
                    last = place;
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
    if let Some(Some(p)) = game.places.get(game.level_index) {
        ui.label(
            egui::RichText::new(format!("{} — {}", p.arc_title(), p.tier))
                .small()
                .weak(),
        )
        .on_hover_text(
            "Each arc introduces its elements one at a time, then combines them in intermediate levels, and ends in a master level that needs everything the arc taught.",
        );
    }
    let level = game.editor.level.clone();
    if !level.description.is_empty() {
        crate::math::paragraph(ui, &level.description, egui::RichText::italics);
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
                let c = crate::draw::shot_color(&level, i).to_srgba();
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
            ui.checkbox(&mut game.show_all_shots, "show all (H)")
                .on_hover_text(
                    "Draw every shot's flight, not only the selected one ([ and ] select)",
                );
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
    // A level without shots (a tube level) has no launch; a neutral placeholder keeps
    // the sections below (which show nothing for it) simple.
    let shot = level.shots.get(shot_index).copied().unwrap_or(level::Shot {
        particle: level::ParticleSpec {
            charge: 0.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        launch: level::Launch {
            node: [0, 0, 0],
            direction: [1.0, 0.0, 0.0],
            kinetic_energy: 1.0,
            time: 0.0,
        },
        detector: level::Detector {
            min: [0, 0, 0],
            max: [0, 0, 0],
            acceptance: None,
        },
        beam: None,
    });
    let t0 = shot.launch.kinetic_energy;
    if !level.shots.is_empty() {
        // Launch of the active shot.
        let kin = physics::dynamics::Kinematics::new(shot.particle.mass, level.c());
        let p0 = level.launch_momentum(shot_index);
        let gamma0 = kin.gamma(p0);
        let v0 = kin.velocity(p0).length();
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
    } else if let Some(t) = level.tube {
        ui.label(format!(
            "Carriers: q/m = {} (electrons; Newtonian), emitted by electrode {}",
            fmt_si(t.charge_per_mass),
            t.cathode + 1
        ));
    }
    egui::CollapsingHeader::new("Physics model and its limits")
        .id_salt("model_notes")
        // Opened for developer captures with EM_MODEL=1 (see `dev_capture`).
        .default_open(std::env::var("EM_MODEL").is_ok_and(|v| v == "1"))
        .show(ui, |ui| {
            // With circuits the notes need the solved circuit: the worker makes them.
            let notes = if level.has_drives() {
                match &game.circuit {
                    Some((r, view)) if *r == game.sent_revision => view.notes.clone(),
                    _ => {
                        ui.label(egui::RichText::new("Solving the circuit…").small().weak());
                        Vec::new()
                    }
                }
            } else {
                level.model_notes(&game.editor.placement)
            };
            for n in notes {
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
    if level.has_drives() {
        egui::CollapsingHeader::new("Circuits")
            .id_salt("circuits")
            .default_open(true)
            .show(ui, |ui| {
                let launches: Vec<f64> = level.shots.iter().map(|s| s.launch.time).collect();
                // Up to the end of the last flight (with a margin), once the previews are
                // in: the flights are often a few % of the circuit's time.
                let ends: Vec<f64> = game
                    .flights
                    .iter()
                    .enumerate()
                    .filter_map(|(i, f)| {
                        let p = f.preview.as_ref()?;
                        Some(level.shots[level.flight_of(i).0].launch.time + p.flight_time)
                    })
                    .collect();
                let window = (ends.len() == game.flights.len() && !ends.is_empty())
                    .then(|| 1.15 * ends.iter().fold(0.0_f64, |m, t| m.max(*t)));
                match &game.circuit {
                    Some((r, view)) if *r == game.sent_revision => {
                        for (label, samples) in &view.plots {
                            ui.label(egui::RichText::new(label).small());
                            let shown: Vec<(f64, f64)> = match window {
                                Some(w) => samples.iter().copied().filter(|s| s.0 <= w).collect(),
                                None => samples.clone(),
                            };
                            circuit_plot(ui, &shown, &launches);
                        }
                        ui.label(
                            egui::RichText::new(
                                "Over lab time, up to the end of the last flight; the yellow \
                                 lines are the launches.",
                            )
                            .small()
                            .weak(),
                        );
                    }
                    _ => {
                        ui.label(egui::RichText::new("Solving the circuit…").small().weak());
                    }
                }
            });
    }
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
        ElementKind::Free => "free charges",
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
    } else if game.editor.kind == ElementKind::Free {
        free_palette(ui, game, &level);
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
                (ElementKind::Plate | ElementKind::Supply | ElementKind::Free, _) => "",
            };
            let hover = match kind {
                ElementKind::Charge => "Flip sign (S)",
                ElementKind::Magnet => {
                    "Flip orientation (S): moment out of the plane (+z) or into it (−z)"
                }
                ElementKind::Antenna => "Flip phase (S): opposite phase of the RF generator",
                ElementKind::Plate | ElementKind::Supply | ElementKind::Free => "",
            };
            ui.label(match kind {
                ElementKind::Charge => "New charge:",
                ElementKind::Magnet => "New magnet μ:",
                ElementKind::Antenna => "New antenna p₀:",
                ElementKind::Plate | ElementKind::Supply | ElementKind::Free => "",
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

    // What the physics thread is doing for this setup.
    match game.progress() {
        crate::Progress::Flying => {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(
                    egui::RichText::new("Computing the new flights…")
                        .color(egui::Color32::from_rgb(255, 210, 120)),
                );
            });
            ui.label(egui::RichText::new("Faded paths: before your change.").small());
        }
        crate::Progress::Verifying => {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(if exact_verdict(&level) {
                    "Verifying with the exact retarded fields…"
                } else {
                    "Verifying…"
                });
            });
        }
        crate::Progress::Done => {}
    }
    if let Some(issues) = game.invalid_setup().map(<[String]>::to_vec) {
        ui.label(egui::RichText::new("Result").strong());
        ui.colored_label(
            egui::Color32::from_rgb(255, 170, 80),
            "This setup cannot be computed:",
        );
        for issue in issues {
            ui.colored_label(egui::Color32::from_rgb(255, 170, 80), issue);
        }
        ui.label(
            egui::RichText::new(
                "Metal bodies must stay apart, and particles must not start in metal or on \
                 a wire: the fields are singular there.",
            )
            .small(),
        );
    } else if level.is_tube() {
        tube_result(ui, game, &level);
    } else if level.has_beams() {
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
                // A steady receiver measures a particle that need not arrive.
                let steady = level.shots[shot_index]
                    .detector
                    .acceptance
                    .and_then(|a| a.radiation)
                    .is_some_and(|r| r.window.is_some());
                match (steady, p.outcome) {
                    (true, Outcome::Arrived) => ui.label(format!(
                        "{who} kept scattering until t = {:.3}; the receiver's mean power is in \
                         its window.",
                        p.flight_time
                    )),
                    (true, Outcome::Rejected) => ui.label(format!(
                        "{who} scattered until t = {:.3}; the receiver's mean power is outside \
                         its window (or the flight ended before the window did).",
                        p.flight_time
                    )),
                    _ => ui.label(format!(
                        "{who} {} after t = {:.3}.",
                        outcome_text(game, p.outcome),
                        p.flight_time
                    )),
                };
                // Detector acceptance: what is allowed, and how the particle arrived.
                if let Some(acc) = level.shots[shot_index].detector.acceptance
                    && let Some(last) = p.path.last()
                {
                    let v = last.p;
                    let dir = v.y.atan2(v.x).to_degrees();
                    let mut parts = Vec::new();
                    if let Some([axis, half]) = acc.direction {
                        let digits = angle_decimals(half);
                        parts.push(format!(
                            "direction {axis:.0}° ± {half:.digits$}° (arrives at {dir:.digits$}°)"
                        ));
                    }
                    if let Some([lo, hi]) = acc.kinetic {
                        parts.push(format!(
                            "energy {lo:.3}–{hi:.3} (arrives with {:.3})",
                            last.kinetic
                        ));
                    }
                    if !parts.is_empty() {
                        ui.label(
                            egui::RichText::new(format!("Detector accepts: {}", parts.join("; ")))
                                .small(),
                        );
                    }
                    if let Some(r) = acc.radiation {
                        radiation_goal(ui, &r, p.radiation, &p.spectrum, p.spectrum_top, false);
                    }
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
                    // Time-dependent fields do work on the particle. Short, with the reason
                    // on hover: with the electrodes' long first column the row widened the
                    // panel.
                    ui.label("n/a").on_hover_text(
                        "The fields vary in time and do work on the particle: its energy is \
                         not conserved",
                    );
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
                        // Short, with the explanation on hover: a long row widens the panel.
                        ui.colored_label(egui::Color32::YELLOW, text + " (!)")
                            .on_hover_text("Above 0.05: the Landau–Lifshitz approximation is strained");
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                } else if level.physics.c.is_some() {
                    ui.label("radiated / T₀ (neglected)");
                    let text = format!("{:.1e}", p.radiated_fraction);
                    if p.radiated_fraction > 1e-10 {
                        ui.colored_label(egui::Color32::YELLOW, text + " (!)")
                            .on_hover_text("Above 1e-10: not negligible");
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                }
                if p.image_force_bound > 0.0 {
                    ui.label("electrode image force (neglected)")
                        .on_hover_text(
                            "Bound on the force of the charge the particle induces on the electrodes, relative to the force that matters (PHYSICS.md §2.7). It grows as the particle passes closer to metal.",
                        );
                    let text = format!("≤ {:.1e}", p.image_force_bound);
                    if p.image_force_bound > level::IMAGE_FORCE_LIMIT {
                        ui.colored_label(egui::Color32::YELLOW, text + " (!)")
                            .on_hover_text("Not negligible: keep away from the metal");
                    } else {
                        ui.label(text);
                    }
                    ui.end_row();
                }
            });
            }
        }
        ui.separator();

        // Energy bars at the animated point of the active shot: one particle, named.
        let (unit, unit_is_t0) = game.energy_unit(level.flight_of(game.active_flight()).0);
        // Wrapped: the long unit note (a particle that gains far more than T₀) overflowed the
        // panel in level 90 with a charge pulling the electron out of its atom.
        ui.horizontal_wrapped(|ui| {
            ui.label(egui::RichText::new("Energy of").strong());
            let (shot_i, d) = level.flight_of(game.active_flight());
            let who = if level.disturbances.is_empty() {
                format!("shot {}", shot_i + 1)
            } else {
                format!("shot {}, disturbance {}", shot_i + 1, d + 1)
            };
            ui.label(
                egui::RichText::new(who)
                    .strong()
                    .color(shot_color32(&level, shot_i)),
            )
            .on_hover_text(
                "The energy of this shot's particle along its flight, in units of its \
                     launch energy T₀; choose the shot with the tabs above or [ ]",
            );
            ui.label(
                egui::RichText::new(if unit_is_t0 {
                    "(units of its T₀)"
                } else {
                    "(units of its largest kinetic energy: it starts at rest)"
                })
                .small(),
            );
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut game.animate, "Animate (A)").on_hover_text(
                "Play the flights in time: the particles move along their paths and the \
             time-dependent maps (waves, particle fields) follow",
            );
            ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"))
                .on_hover_text("Playback speed: 1 plays 4 time units per second");
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
                pt.kinetic / unit,
                egui::Color32::from_rgb(240, 200, 60),
            );
            energy_bar(
                ui,
                "potential",
                pt.potential / unit,
                egui::Color32::from_rgb(200, 90, 230),
            );
            if p.radiation_reaction {
                energy_bar(
                    ui,
                    "radiated",
                    pt.radiated / unit,
                    egui::Color32::from_rgb(90, 200, 255),
                );
                energy_bar(
                    ui,
                    "total − T₀",
                    (pt.kinetic + pt.potential + pt.radiated - t0) / unit,
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
                    (pt.kinetic + pt.potential - t0) / unit,
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
    if level.is_tube() {
        // The z-invariant field: only its potential map (the field lines and the other
        // maps are the 3D slice's).
        ui.horizontal(|ui| {
            ui.label("Map:");
            ui.selectable_value(&mut game.map, Some(MapMode::Potential), "potential");
            ui.selectable_value(&mut game.map, Some(MapMode::Electric), "electric field");
            ui.selectable_value(&mut game.map, None, "off");
        });
        if game.map == Some(MapMode::Electric) {
            ui.label(
                egui::RichText::new(
                    "|E| at the time shown (live while animating), with the electrons' space charge: bright \
                     where strong, on a log scale over 1.5 decades below the field at the 99th \
                     percentile of the board; contours every quarter decade; arrows show its \
                     direction (a force on electrons the opposite way). Inside metal it is 0. \
                     Magnetic and Poynting maps need the circuit the tube's current returns \
                     through, which comes with the circuit coupling.",
                )
                .small(),
            );
            ui.separator();
            controls(ui);
            return;
        }
        ui.label(
            egui::RichText::new(
                "The electrons' potential energy at the time shown (live), relative to the \
                 cathode, in units of the largest voltage: blue downhill (towards the \
                 anode), red uphill; contours every quarter. It includes the electrons' \
                 own space charge and the charge it induces on the electrodes: the flat \
                 region in front of the cathode, where the field is pulled down to zero, \
                 is the space charge that limits the current.",
            )
            .small(),
        );
        ui.separator();
        controls(ui);
        return;
    }
    ui.horizontal_wrapped(|ui| {
        ui.label("Map (V):")
            .on_hover_text("What the background shows; V cycles through the maps");
        ui.selectable_value(&mut game.map, Some(MapMode::Potential), "potential")
            .on_hover_text(
                "The particle's potential energy (static fields): where it is pushed uphill \
                 or downhill, and where energy conservation forbids it to go",
            );
        ui.selectable_value(&mut game.map, Some(MapMode::Magnetic), "magnetic B")
            .on_hover_text(
                "The static magnetic field of magnets, coils and the stray field of the \
                 disturbance shown (perpendicular to the plane)",
            );
        if crate::radiation::waves_available(&level) {
            ui.selectable_value(&mut game.map, Some(MapMode::Waves), "waves")
                .on_hover_text(
                    "The oscillating fields of the antennas and plane waves, moving at c",
                );
        }
        if crate::radiation::particle_field_available(&level) {
            ui.selectable_value(
                &mut game.map,
                Some(MapMode::ParticleField),
                if level.has_beams() {
                    "beam field"
                } else {
                    "particle field"
                },
            )
            .on_hover_text(
                "The field the moving particles themselves carry and radiate (retarded, \
                 Liénard–Wiechert), from their computed flights",
            );
        }
        ui.selectable_value(&mut game.map, Some(MapMode::Total), "total")
            .on_hover_text("Everything at once: level sources, antennas, waves and the particles");
        ui.selectable_value(&mut game.map, None, "off")
            .on_hover_text("No background map");
    });
    if matches!(game.map, Some(MapMode::ParticleField | MapMode::Total)) {
        let (shot, d) = level.flight_of(game.active_flight());
        let which = if level.disturbances.is_empty() {
            format!("shot {}", shot + 1)
        } else {
            format!("shot {}, disturbance {}", shot + 1, d + 1)
        };
        let multi = level.flight_count() > 1;
        let text = if level.has_beams() {
            if level.disturbances.is_empty() {
                "Field of every particle of the beam.".to_string()
            } else {
                format!(
                    "Field of every particle of the beam, disturbance {}.",
                    d + 1
                )
            }
        } else {
            format!(
                "Field of the particle of {which}{}; only its flight is drawn.",
                if multi {
                    " (choose with the shot/disturbance tabs or [ ])"
                } else {
                    ""
                }
            )
        };
        ui.label(egui::RichText::new(text).small());
    }
    if matches!(
        game.map,
        Some(MapMode::Waves | MapMode::ParticleField | MapMode::Total)
    ) {
        field_view_controls(ui, game, &level, radiation);
    }
    // Fields that change in time do work on the particle: no energy limit then.
    let (_, disturbance) = level.flight_of(game.active_flight());
    let time_dependent = level.coils.iter().any(level::Coil::is_time_dependent)
        || level.has_drives()
        || level
            .elements
            .iter()
            .chain(&game.editor.placement)
            .any(|e| {
                e.kind == level::ElementKind::Antenna
                    && e.omega.unwrap_or(level.physics.rf_omega) != 0.0
            })
        || level
            .disturbances
            .get(disturbance)
            .is_some_and(|d| d.waves.iter().any(|w| w.omega != 0.0));
    let legend = match game.map {
        Some(MapMode::Potential) if time_dependent => {
            "Red: uphill for this shot's particle, blue: downhill (the static part of the \
             field; plates and coils driven by circuits as they are at t = 0). No dark \
             region: the time-dependent fields (a ramped or driven coil's induced field, \
             driven plates, antennas, waves) do work on the particle, so energy \
             conservation forbids nothing here."
        }
        Some(MapMode::Potential) => {
            "Red: uphill for this shot's particle, blue: downhill; contours every T₀/4. \
             Dark: forbidden by energy conservation for every particle shown (all shots \
             with \"show all\"; a beam counts with its most energetic particle). The \
             flight's own static field, with the disturbance shown and magnets: exact for \
             charges, clouds, coils, antennas and stray fields; metal from its picture \
             model (spheres to about 1e-4, electrodes on a coarser mesh, within about 0.5 % \
             of their potential). With interacting beam particles only a guide, since \
             they exchange energy."
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
            "Fields of the antennas and waves at the animation time (their exact retarded \
             fields, drawn in single precision). Colour: B perpendicular to the plane \
             (orange out, blue in), on the chosen scale (bar above); arrows: E. Near an \
             antenna the field is quasi-static; further out the radiation travels outwards \
             at c."
        }
        Some(MapMode::Total) => {
            "The total field at the animation time: level charges, magnets and coils, \
             antennas, waves and disturbances, and the particle's own field (retarded, \
             Liénard–Wiechert). The static part is computed in double precision on 6 \
             points per cell and interpolated between them (metal from its picture model); \
             antennas, waves and the particle are evaluated per pixel. The particle's field \
             is usually far weaker than the electrodes'; raise the range to see it. Colour: \
             B_z or |E| on the chosen scale (bar above); arrows: E."
        }
        Some(MapMode::ParticleField) if level.has_beams() && game.neglected_only => {
            "What the quasi-static beam interaction leaves out of the dynamics: the full \
             retarded field of all particles minus the fields it uses (each particle's past \
             continued along the motion it would have in the fields it feels now). Mostly \
             how those fields change along its path during the light travel time, and the \
             delay with which the news of an absorption spreads. Shown on the full field's \
             colour scale (linear by default), so its true size is seen; its size is stated \
             above."
        }
        Some(MapMode::ParticleField) if level.has_beams() => {
            "The retarded (Liénard–Wiechert, exact) field of every particle of the beam at \
             the animation time, from the flight shown (the exact one once the verdict is \
             in). Every change of velocity sends out radiation at c. Before launch each \
             particle is taken to move with its launch acceleration. One entering its \
             detector flies on into the screening cup, its charge fading as seen from \
             outside (with the instant drain its field disappears as the news of its \
             absorption spreads at c); one absorbed by a body stays there at rest. Colour: \
             B perpendicular to the plane, on the chosen scale; arrows: E."
        }
        Some(MapMode::ParticleField) if shot.particle.moment != 0.0 => {
            "The field of the particle itself at the animation time. Its charge's field is \
             Liénard–Wiechert (exact); its magnetic moment is drawn as the static dipole \
             field B_z = −m/r³ from its retarded position, without the terms of a moving \
             dipole (of order v/c). Colour: B perpendicular to the plane, on the chosen \
             scale; arrows: E. Before launch the particle is taken to move uniformly."
        }
        Some(MapMode::ParticleField) => {
            "The field of the particle itself (Liénard–Wiechert, exact) at the animation \
             time. Every change of velocity sends out a radiation pulse at c; the energy it \
             carries is what the particle loses (radiation reaction, when included). \
             Colour: B perpendicular to the plane, on the chosen scale; arrows: E. Before \
             launch the particle is taken to move uniformly."
        }
        // Tube levels only (their own View section).
        Some(MapMode::Electric) | None => "",
    };
    if !legend.is_empty() {
        ui.label(egui::RichText::new(legend).small());
    }
    if game.map == Some(MapMode::ParticleField)
        && (!level.conductors.is_empty() || !level.electrodes.is_empty())
    {
        ui.label(
            egui::RichText::new(
                "Near metal: the charges the particle induces on it reshape its field there \
                 (and screen it inside); that part is not drawn. (In the dynamics their \
                 force is included for spheres and bounded for electrodes, PHYSICS.md \
                 §2.6–2.7.)",
            )
            .small()
            .color(egui::Color32::from_rgb(255, 210, 120)),
        );
    }
    ui.checkbox(&mut game.show_field_lines, "Electric field lines (F)")
        .on_hover_text(
            "Lines along the static electric field (with the stray field of the disturbance \
             shown)",
        );
    ui.add(egui::Slider::new(&mut game.field_line_spacing, 0.5..=4.0).text("spacing (cells)"))
        .on_hover_text("Smallest distance between neighbouring field lines");
    ui.add(egui::Slider::new(&mut game.field_line_opacity, 0.05..=1.0).text("opacity"))
        .on_hover_text("How strongly the field lines are drawn");
    ui.label(
        egui::RichText::new(
            "In this 2D slice of a 3D field, lines show direction only, not strength. Screened regions (behind grounded metal) get few lines.",
        )
        .small(),
    );
    ui.separator();

    controls(ui);
}

fn controls(ui: &mut egui::Ui) {
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

/// Tube levels: the goal electrode's current, preview and verified, against the goal,
/// and its collected charge over time (the window's ends marked).
fn tube_result(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    use level::tube::TubeStatus;
    let Some(spec) = level.tube else {
        return;
    };
    let g = spec.goal;
    ui.horizontal(|ui| {
        ui.checkbox(&mut game.animate, "Animate (A)")
            .on_hover_text("Play the preview: the electrons leave the cathode and fly");
        ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"))
            .on_hover_text("Playback speed: 1 plays 1 time unit per second");
    });
    if let Some(v) = &game.tube_preview
        && let Some(last) = v.run.frames.last()
    {
        let end = level.tube_display_end();
        ui.label(
            egui::RichText::new(if v.done {
                format!("t = {:.1} of {end:.0}", game.anim_time.min(last.t))
            } else {
                format!(
                    "t = {:.1} of {end:.0} (computed to {:.0}; playback waits for it)",
                    game.anim_time.min(last.t),
                    last.t
                )
            })
            .small(),
        );
    }
    ui.label(egui::RichText::new("Result").strong());
    if game.solved() {
        ui.label(
            egui::RichText::new("✔ SOLVED: the current is in range (verified)")
                .color(egui::Color32::from_rgb(90, 240, 110))
                .size(18.0),
        );
    }
    ui.label(format!(
        "Goal: the current into electrode {} between {} and {}, averaged over t = {} to {}",
        g.electrode + 1,
        fmt_si(g.min),
        fmt_si(g.max),
        fmt_si(g.start),
        fmt_si(g.end)
    ));
    let current = game.sent_revision;
    match &game.tube_preview {
        Some(v) if v.revision == current => {
            let run = &v.run;
            let shown = game
                .tube_frame()
                .and_then(|k| run.frames.get(k))
                .map_or(0, |f| f.x.len());
            ui.label(match run.current {
                Some(c) => format!("Preview: {}   ({shown} particles in flight now)", fmt_si(c)),
                None => format!("Preview: computing…   ({shown} particles in flight now)"),
            });
            ui.label(egui::RichText::new("Charge collected by the goal electrode:").small());
            circuit_plot(ui, &run.collected, &[g.start, g.end]);
            let quiet = run.current.is_some()
                && run
                    .frames
                    .iter()
                    .take_while(|f| f.t <= g.end)
                    .all(|f| f.x.is_empty());
            if quiet {
                ui.label(
                    egui::RichText::new(
                        "No electrons leave the cathode: the field at its face pushes them \
                         back. A diode conducts one way only.",
                    )
                    .color(egui::Color32::from_rgb(255, 210, 120)),
                );
            }
        }
        _ => {
            ui.label("Preview: computing…");
        }
    }
    match &game.tube_verdict {
        Some((r, v)) if *r == current => {
            let (text, color) = match v.status {
                TubeStatus::Met => ("in range ✔", egui::Color32::from_rgb(90, 240, 110)),
                TubeStatus::Missed => ("out of range ✖", egui::Color32::from_rgb(255, 120, 90)),
                TubeStatus::Uncertain => (
                    "too close to a bound to tell ✖",
                    egui::Color32::from_rgb(255, 210, 90),
                ),
            };
            ui.colored_label(
                color,
                format!("Verified: {} ± {}: {text}", fmt_si(v.fine), fmt_si(v.error)),
            );
        }
        _ => {
            ui.label("Verified: computing…");
        }
    }
    ui.label(
        egui::RichText::new(
            "A z-invariant tube: every electrode is a long prism, the electrons are long \
             lines of charge, and the current is per unit length. The cathode emits as \
             much as its field allows: the electrons' own charge limits the current.",
        )
        .small(),
    );
}

/// Beam levels: verified transmission of every beam shot in every flight against its
/// requirement, and the diagnostics of the selected flight.
fn beam_result(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    ui.horizontal(|ui| {
        ui.checkbox(&mut game.animate, "Animate (A)").on_hover_text(
            "Play the flights in time: the particles move along their paths and the \
             time-dependent maps (waves, particle fields) follow",
        );
        ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"))
            .on_hover_text("Playback speed: 1 plays 4 time units per second");
    });
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
                        if enough { "✔" } else { "✖" }
                    )
                }
            };
            ui.colored_label(shot_color32(level, s), text);
        }
        // Goal particles: the level's free particles with a detector.
        let n = level.shots.len();
        for (k, f) in level.free_particles.iter().enumerate() {
            if f.detector.is_none() {
                continue;
            }
            let text = match game.beam_transmission(d, n + k) {
                None => format!("Particle {}: computing…", k + 1),
                Some((1, _)) => format!("Particle {}: arrives, verified ✔", k + 1),
                Some(_) => format!("Particle {}: does not arrive (verified) ✖", k + 1),
            };
            let c = crate::draw::free_goal_color(k).to_srgba();
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let color = egui::Color32::from_rgb(
                (c.red * 255.0) as u8,
                (c.green * 255.0) as u8,
                (c.blue * 255.0) as u8,
            );
            ui.colored_label(color, text);
        }
    }
    let d = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    let view = game.beams.get(d);
    if let Some(p) = view.and_then(crate::BeamView::shown) {
        // Lost: particles with a goal (the shots' and the level's free particles with a
        // detector) that do not arrive; the player's free charges have none.
        let n_shots = level.shots.len();
        let has_goal = |s: usize| {
            s < n_shots
                || (s != usize::MAX
                    && level
                        .free_particles
                        .get(s - n_shots)
                        .is_some_and(|f| f.detector.is_some()))
        };
        let lost = p
            .outcomes
            .iter()
            .zip(&p.shots)
            .filter(|(o, s)| has_goal(**s) && **o != Outcome::Arrived)
            .count();
        // Which flight the views show: the quick preview, or the verdict's own flight
        // once it has arrived (exact at finite c).
        let verdict_shown = view.is_some_and(|v| v.verdict_flight.is_some());
        let energy = if p.energy_rel_error.is_nan() {
            // Finite c: the particles exchange energy with the field.
            let mut s = format!(
                "largest radiated share of a particle {:.1e} ({})",
                p.radiated_max,
                if level.physics.radiation_reaction {
                    "included"
                } else {
                    "neglected"
                }
            );
            if p.retardation_max > 0.0 {
                s += &format!(
                    "; approximation of the fields (quasi-static), estimated error {:.1e} of the interaction",
                    p.retardation_max
                );
            }
            s
        } else {
            format!("energy drift of the whole beam {:.1e}", p.energy_rel_error)
        };
        // Radiation goals: every charge's radiation together.
        for r in &p.radiation {
            if let Some(goal) = level
                .shots
                .get(r.shot)
                .and_then(|s| s.detector.acceptance)
                .and_then(|a| a.radiation)
            {
                radiation_goal(ui, &goal, Some(r.value), &r.spectrum, r.spectrum_top, true);
            }
        }
        ui.label(
            egui::RichText::new(format!(
                "{} particles{}; {lost} lost{}; {energy}.",
                p.paths.len(),
                match (level.interacts(), level.physics.c) {
                    (false, _) => ", not interacting",
                    (true, None) => ", interacting",
                    (true, Some(_)) if p.retarded => ", interacting (exact retarded fields)",
                    (true, Some(_)) => ", interacting (quasi-static fields)",
                },
                if verdict_shown { "" } else { " in the preview" },
            ))
            .small(),
        );
        let shown_text = match view {
            Some(v) if verdict_shown => v.preview.as_ref().map(|quick| {
                let (worst, differ) = quick_vs_verdict(quick, p);
                let mut s = String::from(if p.retarded {
                    "Shown in every view: the verdict's own flight (exact retarded fields). "
                } else {
                    "Shown in every view: the verdict's own flight (metal at its verification \
                     resolution). "
                });
                // The sandbox's exact option flies the preview retarded too.
                s += if p.retarded && !quick.retarded {
                    "The quick preview (quasi-static)"
                } else {
                    "The preview (looser tolerance)"
                };
                s += &format!(" ended within {worst:.1e} cells of it");
                if quick.retardation_max > 0.0 {
                    s += &format!(
                        " (its estimated error: {:.1e} of the interaction)",
                        quick.retardation_max
                    );
                }
                if differ > 0 {
                    s += &format!("; {differ} particles ended differently");
                }
                s + "."
            }),
            Some(v) if v.verified.is_none() && exact_verdict(level) => Some(
                "Shown: the quick preview (quasi-static fields). The exact flight, the \
                 verdict's, replaces it in every view when it is ready."
                    .to_string(),
            ),
            _ => None,
        };
        if let Some(text) = shown_text {
            ui.label(egui::RichText::new(text).small());
        }
    }
    beam_energy(ui, game, level);
    ui.separator();
}

/// Whether the verdict of a beam level flies the exact retarded interaction while the
/// preview is quasi-static (`Level::verification_beam_scenarios`).
fn exact_verdict(level: &level::Level) -> bool {
    level.has_beams()
        && level.interacts()
        && level.physics.c.is_some()
        && !level.physics.beam_retarded
}

/// How far the quick preview's end points are from the verdict flight's (particles with
/// the same outcome), and how many particles ended differently.
fn quick_vs_verdict(
    quick: &crate::worker::BeamPreview,
    verdict: &crate::worker::BeamPreview,
) -> (f64, usize) {
    let mut worst: f64 = 0.0;
    let mut differ = 0;
    for (i, (a, b)) in quick.paths.iter().zip(&verdict.paths).enumerate() {
        if quick.outcomes.get(i) != verdict.outcomes.get(i) {
            differ += 1;
        } else if let (Some(x), Some(y)) = (a.last(), b.last()) {
            worst = worst.max((x.1 - y.1).length());
        }
    }
    (worst, differ)
}

/// The drawing colour of a shot, for egui.
fn shot_color32(level: &level::Level, shot: usize) -> egui::Color32 {
    let c = crate::draw::shot_color(level, shot).to_srgba();
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

/// A slider for a potential over `[lo, hi]` (power supplies, plates); the value can also
/// be typed (SI prefixes).
pub fn potential_slider(v: &mut f64, lo: f64, hi: f64) -> egui::Slider<'_> {
    egui::Slider::new(v, lo..=hi)
        .custom_formatter(|v, _| fmt_si(v))
        .custom_parser(parse_si)
}

/// Palette row for new plates: potential and orientation.
/// New free charges: their charge, launch speed and direction. A placed one's velocity
/// is set by dragging the handle at the tip of its arrow.
fn free_palette(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    if !game.editor.continuous() {
        ui.horizontal_wrapped(|ui| {
            ui.label("New free charge q:");
            for (i, v) in level.limits.free_charges.iter().enumerate() {
                ui.selectable_value(&mut game.editor.magnitude_index, i, fmt_si(*v))
                    .on_hover_text("Q/E or the wheel: next charge");
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("speed:");
            let mut k = game.editor.speed_index;
            for (i, v) in level.limits.free_speeds.iter().enumerate() {
                ui.selectable_value(&mut k, i, fmt_si(*v));
            }
            game.editor.set_speed_index(k);
            ui.label(format!("direction {:.0}°", game.editor.free_angle_deg))
                .on_hover_text("R / Shift+R: turn by 15°");
        });
    }
    // The placed ones: velocity as v/c and rapidity.
    let c = level.c();
    for (k, e) in game
        .editor
        .placement
        .iter()
        .filter(|e| e.kind == ElementKind::Free)
        .enumerate()
    {
        let v = e.speed.unwrap_or(0.0);
        let text = if c.is_finite() {
            format!(
                "Free charge {}: q = {}, v = {:.3} c, rapidity φ = {:.3}, towards {:.0}°",
                k + 1,
                fmt_si(e.value),
                v / c,
                (v / c).atanh(),
                e.angle_deg
            )
        } else {
            format!(
                "Free charge {}: q = {}, v = {}, towards {:.0}°",
                k + 1,
                fmt_si(e.value),
                fmt_si(v),
                e.angle_deg
            )
        };
        ui.label(egui::RichText::new(text).small());
    }
    if c.is_finite() {
        let mut rapidity = game.arrow_measure == crate::editor::ArrowMeasure::Rapidity;
        if ui
            .checkbox(&mut rapidity, "arrows show rapidity")
            .on_hover_text(
                "Arrow length ∝ c·artanh(v/c) instead of v: rapidities add under boosts along \
                 a line; for slow particles it is the speed itself, and it grows without \
                 bound as v approaches c.",
            )
            .changed()
        {
            game.arrow_measure = if rapidity {
                crate::editor::ArrowMeasure::Rapidity
            } else {
                crate::editor::ArrowMeasure::Speed
            };
        }
    }
    ui.label(
        egui::RichText::new(format!(
            "Free charges move and push and pull every other particle with their fields \
             (mass {}, radius {}; if two ever touch, they bounce as rigid spheres). To launch \
             one, pull the handle behind it back, like a slingshot: the arrow shows where \
             it goes and how fast (c = {}).",
            fmt_si(level.limits.free_mass),
            level.limits.free_radius,
            level.physics.c.map_or_else(|| "∞".to_string(), fmt_si)
        ))
        .small(),
    );
}

fn plate_palette(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    // The potential: a slider in every mode (hardcore shows its own), for the plate under
    // the cursor if there is one, else for new plates.
    if !game.editor.continuous()
        && let Some((lo, hi)) = level::value_range(&level.limits.plate_voltages)
    {
        let target = game
            .editor
            .element_at_cursor()
            .filter(|&i| game.editor.placement[i].kind == ElementKind::Plate);
        let k = crate::editor::kind_index(ElementKind::Plate);
        let mut v = target.map_or(game.editor.continuous_magnitude[k], |i| {
            game.editor.placement[i].value
        });
        ui.horizontal(|ui| {
            ui.label(if target.is_some() {
                "This plate:"
            } else {
                "New plate:"
            });
            let r = ui
                .add(potential_slider(&mut v, lo, hi))
                .on_hover_text("Q/E or the wheel: a step; S: the opposite potential");
            game.text_focus |= r.has_focus();
        });
        match target {
            Some(i) => {
                let mut e = game.editor.placement[i];
                if e.value.to_bits() != v.to_bits() {
                    e.value = v;
                    game.editor.set_element(i, e);
                }
            }
            None => game.editor.continuous_magnitude[k] = v,
        }
    }
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

/// The power supplies of the level's tunable electrodes and coils: off (the electrode
/// keeps its own bias, the coil its own ramp) or a value on its slider. Clicking an
/// electrode, or a coil's ring, on the map switches it on and steps it.
fn supplies(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    let targets = level.supply_targets();
    if targets.is_empty() {
        return;
    }
    ui.label(egui::RichText::new("Power supplies").strong())
        .on_hover_text(
            "Click a tunable electrode (yellow frame) or a tunable coil's ring (yellow inner \
             ring) to switch its supply on or to step it; S flips it, right click switches it \
             off.",
        );
    for (centre, target) in targets {
        let (name, own) = match target {
            level::SupplyTarget::Electrode(i) => {
                let own = match level.electrodes[i].bias {
                    level::ConductorBias::Grounded => "grounded".to_string(),
                    level::ConductorBias::Potential(v) => fmt_potential(v),
                    level::ConductorBias::Charge(q) => format!("charge {}", fmt_si(q)),
                };
                (
                    format!("Electrode {}:", i + 1),
                    format!("the electrode is {own}"),
                )
            }
            level::SupplyTarget::Coil(i) => {
                let own = match &level.coils[i] {
                    level::Coil::Circle { drive: Some(_), .. } => {
                        "its circuit keeps its own source".to_string()
                    }
                    level::Coil::Circle { rate, .. } => {
                        format!("the coil ramps at dκ/dt = {}", fmt_si(*rate))
                    }
                    level::Coil::Polygon { .. } => String::new(),
                };
                // A driven coil's supply is its circuit's source voltage, an undriven
                // one's its ramp rate.
                let what = if level.coils[i].is_driven() {
                    "source V"
                } else {
                    "dκ/dt"
                };
                (format!("Coil {} {what}:", i + 1), own)
            }
        };
        let list = level.supply_list(target).to_vec();
        let current = game.editor.supply(centre);
        ui.horizontal_wrapped(|ui| {
            ui.label(name);
            // A slider over the level's range, in every mode (the owner: values are tuned,
            // not picked from a list whose ends give the answer away).
            let mut choice = current;
            let mut on = current.is_some();
            ui.checkbox(&mut on, "on")
                .on_hover_text(format!("Off, {own}"));
            if let Some((lo, hi)) = level::value_range(&list) {
                let mut v = current.unwrap_or(lo.max(0.0_f64.min(hi)));
                let r = ui.add_enabled(on, potential_slider(&mut v, lo, hi));
                game.text_focus |= r.has_focus();
                choice = on.then_some(v);
            }
            game.editor.set_supply(centre, choice);
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
            angle_deg: if game.editor.kind == ElementKind::Free {
                game.editor.free_angle_deg
            } else {
                game.editor.angle_deg
            },
            omega: game.editor.selected_omega(),
            speed: (game.editor.kind == ElementKind::Free).then(|| game.editor.selected_speed()),
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
        ElementKind::Free => "q",
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
        if e.kind == ElementKind::Free {
            ui.horizontal(|ui| {
                ui.label("direction");
                let r = ui.add(
                    egui::Slider::new(&mut e.angle_deg, 0.0..=360.0)
                        .suffix("°")
                        .step_by(0.5),
                );
                focus |= r.has_focus();
            });
            if let (Some((lo, hi)), Some(v)) =
                (value_range(&level.limits.free_speeds), e.speed.as_mut())
            {
                ui.horizontal(|ui| {
                    ui.label("speed");
                    if lo < hi {
                        let r = ui.add(egui::Slider::new(v, lo..=hi));
                        focus |= r.has_focus();
                    } else {
                        ui.label(fmt_si(lo));
                    }
                });
            }
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
        None if e.kind == ElementKind::Free => {
            let k = crate::editor::kind_index(e.kind);
            game.editor.continuous_magnitude[k] = e.value;
            game.editor.free_angle_deg = e.angle_deg;
            if let Some(v) = e.speed {
                game.editor.continuous_speed = v;
            }
        }
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

/// Controls of the time-dependent field views: what to show (the full field, or a part of
/// it measured against the full field), the colour quantity and scale, a colour bar with
/// ticks, and the size of the part shown.
fn field_view_controls(
    ui: &mut egui::Ui,
    game: &mut Game,
    level: &level::Level,
    radiation: &crate::radiation::RadiationView,
) {
    use crate::radiation::FieldQuantity;
    // What to show: the full field (baseline) or a part of it.
    let quasi_static = level.has_beams()
        && level.interacts()
        && !level.physics.beam_retarded
        && level.physics.c.is_some();
    if game.map == Some(MapMode::ParticleField) {
        let before = (game.radiation_only, game.neglected_only);
        ui.horizontal_wrapped(|ui| {
            ui.label("Show:")
                .on_hover_text("The full field, or a part of it on the full field's colour scale");
            if ui
                .selectable_label(!game.radiation_only && !game.neglected_only, "full field")
                .on_hover_text("The whole field of the moving charges: the baseline")
                .clicked()
            {
                (game.radiation_only, game.neglected_only) = (false, false);
            }
            // At c = ∞ nothing radiates.
            if level.physics.c.is_some()
                && ui
                    .selectable_label(game.radiation_only, "radiation part")
                    .on_hover_text(
                        "Only the acceleration term of the field: what travels away at c and \
                     falls off as 1/R (the rest is the velocity field that moves with the \
                     charge)",
                    )
                    .clicked()
            {
                (game.radiation_only, game.neglected_only) = (true, false);
            }
            if quasi_static
                && ui
                    .selectable_label(game.neglected_only, "left out by the model")
                    .on_hover_text(
                        "The full retarded field minus the fields the quasi-static beam \
                         interaction uses (present states continued with constant \
                         acceleration): what the dynamics leaves out (PHYSICS.md §3.3)",
                    )
                    .clicked()
            {
                (game.radiation_only, game.neglected_only) = (false, true);
            }
        });
        // A part is best seen on a linear scale: switch when the choice changes.
        let now = (game.radiation_only, game.neglected_only);
        if now != before {
            game.field_linear = now != (false, false);
        }
    }
    // The energy flow needs a finite speed of light (the magnetic energy of a finite B is
    // infinite at c = ∞ in these units).
    if level.physics.c.is_none() && game.field_quantity == FieldQuantity::S {
        game.field_quantity = FieldQuantity::E;
    }
    let flow = game.field_quantity == FieldQuantity::S;
    ui.horizontal_wrapped(|ui| {
        ui.label("Colour:")
            .on_hover_text("Which quantity colours the map");
        ui.selectable_value(&mut game.field_quantity, FieldQuantity::Bz, "B_z")
            .on_hover_text("B, perpendicular to the plane (all of B in the plane); signed");
        ui.selectable_value(&mut game.field_quantity, FieldQuantity::E, "|E|")
            .on_hover_text("The magnitude of E");
        if level.physics.c.is_some() {
            ui.selectable_value(&mut game.field_quantity, FieldQuantity::S, "S")
                .on_hover_text(
                    "The flow of field energy, the Poynting vector S = (c²/4π) E × B \
                     (Jackson §6.7): its magnitude colours the map, arrows and tracers show \
                     where it goes",
                );
        }
        ui.checkbox(
            &mut game.show_field_arrows,
            if flow { "S arrows" } else { "E arrows" },
        )
        .on_hover_text(if flow {
            "Arrows of the energy flow S: direction exact, length on the colour scale"
        } else {
            "Arrows of E: direction exact, length on the colour scale"
        });
    });
    if flow {
        flow_controls(ui, game, radiation);
    }
    ui.horizontal(|ui| {
        ui.label("Scale:").on_hover_text("How field values map to colour");
        ui.selectable_value(&mut game.field_linear, false, "log")
            .on_hover_text("Logarithmic: shows fields over many decades at once, but makes small ones look big");
        ui.selectable_value(&mut game.field_linear, true, "linear")
            .on_hover_text("Linear: brightness proportional to the field; honest sizes, small parts faint");
    });
    if game.field_linear {
        ui.add(egui::Slider::new(&mut game.field_gain_decades, 0.0..=8.0).text("gain (decades)"))
            .on_hover_text(
                "Amplifies the colours: full colour at 10^−gain of the full field's scale. \
                 The bar below states the values",
            );
    } else {
        ui.add(
            egui::Slider::new(&mut game.field_range_decades, 1.0..=14.0).text("range (decades)"),
        )
        .on_hover_text("How many decades below the full field's scale are still visible");
    }
    colour_bar(ui, game, radiation);
    // How large the part shown is (for the energy flow, each part has its own scale).
    if game.map == Some(MapMode::ParticleField)
        && (game.radiation_only || game.neglected_only)
        && !flow
    {
        let part = match game.field_quantity {
            FieldQuantity::Bz => radiation.part_b,
            FieldQuantity::E | FieldQuantity::S => radiation.part_e,
        };
        ui.label(
            egui::RichText::new(format!(
                "Largest values of this part: {} of the full field's scale.",
                percent(part)
            ))
            .color(egui::Color32::from_rgb(255, 210, 120)),
        );
    }
}

/// Controls of the energy-flow view: the part shown, tracers, the average over a period,
/// and what the view means (with its caveats).
fn flow_controls(ui: &mut egui::Ui, game: &mut Game, radiation: &crate::radiation::RadiationView) {
    use crate::radiation::FlowPart;
    ui.horizontal_wrapped(|ui| {
        ui.label("Flow:")
            .on_hover_text("Which part of the energy flow: the view's field is the moving charges' field plus the rest");
        for (part, name, hint) in [
            (
                FlowPart::Total,
                "total",
                "All of it: the flow of the whole field shown",
            ),
            (
                FlowPart::Own,
                "own",
                "The moving charges' field alone: the energy moving with them, and their \
                 radiation leaving at c",
            ),
            (
                FlowPart::Exchange,
                "exchange",
                "The cross terms of the charges' field with the rest, (c²/4π)(E₁ × B₂ + \
                 E₂ × B₁): they carry the work the rest does on the charges, so they show \
                 where a particle's energy comes from",
            ),
            (
                FlowPart::External,
                "external",
                "The rest alone: the static sources (in the total view), antennas and waves",
            ),
        ] {
            ui.selectable_value(&mut game.flow_part, part, name)
                .on_hover_text(hint);
        }
    });
    ui.horizontal_wrapped(|ui| {
        ui.checkbox(&mut game.flow_tracers, "tracers")
            .on_hover_text(
                "Dots drifting with the energy's velocity S/u (u the energy density): never \
             faster than c, exactly c in a radiation field. The exchange terms have no \
             velocity of their own: there the dots move along S, at c at full colour",
            );
        if radiation.oscillates() {
            ui.checkbox(&mut game.flow_average, "average over a period")
                .on_hover_text(
                    "The antennas' and waves' flow averaged over their period (for several \
                     frequencies, each exactly), the charges held still: the energy that \
                     sloshes back and forth near an antenna drops out and the net outflow \
                     remains. Off: the flow at the moment shown",
                );
        }
    });
    ui.label(
        egui::RichText::new(
            "S = (c²/4π) E × B, the flow of field energy (Jackson §6.7). Only its flux \
             through closed surfaces is unique: where the energy goes in between is \
             Poynting's convention. The map is a slice of a 3D flow; energy also leaves \
             the plane. Charges beside magnets keep a steady circulation (field momentum) \
             that carries energy nowhere. Near a particle part of the work done on it comes \
             from the exchange energy stored around it (a third at low speed), the rest \
             flows in.",
        )
        .small(),
    );
}

/// A fraction as a percentage with sensible digits.
fn percent(f: f64) -> String {
    let p = f * 100.0;
    if p >= 10.0 {
        format!("{p:.0} %")
    } else if p >= 0.1 {
        format!("{p:.1} %")
    } else {
        format!("{p:.1e} %")
    }
}

/// The beam's energy budget at the animation time (the active disturbance's flight), in
/// units of the beam's launch kinetic energy: where the energy went.
fn beam_energy(ui: &mut egui::Ui, game: &Game, level: &level::Level) {
    let d = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    let Some(p) = game.beams.get(d).and_then(crate::BeamView::shown) else {
        return;
    };
    let (Some(first), Some(last)) = (p.energy.first(), p.energy.last()) else {
        return;
    };
    let t = if game.animate { game.anim_time } else { last.t };
    // Linear interpolation between the recorded samples.
    let k = p
        .energy
        .partition_point(|e| e.t <= t)
        .clamp(1, p.energy.len().max(2) - 1);
    let e = if p.energy.len() < 2 || t >= last.t {
        *last
    } else {
        let (a, b) = (p.energy[k - 1], p.energy[k]);
        let f = if b.t > a.t {
            ((t - a.t) / (b.t - a.t)).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let mix = |x: f64, y: f64| x + (y - x) * f;
        physics::beam::EnergySample {
            t,
            kinetic: mix(a.kinetic, b.kinetic),
            potential: mix(a.potential, b.potential),
            interaction: mix(a.interaction, b.interaction),
            absorbed: mix(a.absorbed, b.absorbed),
            radiated: mix(a.radiated, b.radiated),
        }
    };
    let t0 = first.kinetic.max(1e-300);
    ui.label(egui::RichText::new("Energy of the whole beam").strong())
        .on_hover_text(
            "Where the beam's energy went, in units of its kinetic energy at launch T₀ \
         (all particles together)",
        );
    energy_bar(
        ui,
        "kinetic",
        e.kinetic / t0,
        egui::Color32::from_rgb(240, 200, 60),
    );
    energy_bar(
        ui,
        "potential",
        (e.potential - first.potential) / t0,
        egui::Color32::from_rgb(200, 90, 230),
    );
    if level.interacts() {
        energy_bar(
            ui,
            "mutual",
            (e.interaction - first.interaction) / t0,
            egui::Color32::from_rgb(230, 130, 90),
        );
    }
    energy_bar(
        ui,
        "absorbed",
        e.absorbed / t0,
        egui::Color32::from_rgb(160, 160, 160),
    );
    let rr = level.physics.radiation_reaction;
    if level.physics.c.is_some() {
        energy_bar(
            ui,
            "radiated",
            e.radiated / t0,
            egui::Color32::from_rgb(90, 200, 255),
        );
    }
    let mut total = e.kinetic
        + (e.potential - first.potential)
        + (e.interaction - first.interaction)
        + e.absorbed
        - first.kinetic;
    if rr {
        total += e.radiated;
    }
    energy_bar(
        ui,
        "total − T₀",
        total / t0,
        egui::Color32::from_rgb(120, 220, 120),
    );
    let note = match (level.physics.c.is_some(), rr, level.interacts()) {
        (false, _, _) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector. The total is conserved."
        }
        (true, true, _) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: carried away by the field (radiation reaction included, so \
             it is part of the total). Mutual: the Coulomb part of the particles' \
             interaction only (their magnetic and radiation field energy is not counted)."
        }
        (true, false, true) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: an estimate of what the particles would radiate, neglected \
             in the dynamics (not in the total). Mutual: the Coulomb part of the \
             interaction only."
        }
        (true, false, false) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: an estimate, neglected in the dynamics (not in the total)."
        }
    };
    ui.label(egui::RichText::new(note).small());
    ui.label(format!("t = {:.2}", e.t));
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]
/// The colour bar of the field views, with tick marks: the value of full colour, and
/// where tenths (linear) or decades (logarithmic) of it fall.
fn colour_bar(ui: &mut egui::Ui, game: &Game, radiation: &crate::radiation::RadiationView) {
    use crate::radiation::{FieldQuantity, colour_value};
    let signed = game.field_quantity == FieldQuantity::Bz;
    let sat = match game.field_quantity {
        FieldQuantity::Bz => radiation.b_sat(),
        FieldQuantity::E => radiation.e_sat(),
        FieldQuantity::S => radiation.s_sat(game.flow_part),
    };
    let name = match game.field_quantity {
        FieldQuantity::Bz => "±B_z",
        FieldQuantity::E => "|E|",
        FieldQuantity::S => "|S|",
    };
    let linear = game.field_linear;
    let scale = if linear {
        10f64.powf(game.field_gain_decades)
    } else {
        10f64.powf(game.field_range_decades)
    };
    if sat <= 0.0 {
        ui.label(
            egui::RichText::new(format!(
                "{} = 0 everywhere: nothing to colour.",
                name.trim_start_matches('±')
            ))
            .small(),
        );
        return;
    }
    let width = ui.available_width().min(260.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 26.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let bar = egui::Rect::from_min_size(rect.min, egui::vec2(width, 10.0));
    let colour = |s: f64| {
        let a = (s.abs().powf(0.8) * 0.9) as f32;
        let base = if game.field_quantity == FieldQuantity::S {
            egui::Color32::from_rgb(140, 255, 158)
        } else if !signed {
            egui::Color32::from_rgb(255, 235, 140)
        } else if s >= 0.0 {
            egui::Color32::from_rgb(255, 150, 40)
        } else {
            egui::Color32::from_rgb(40, 190, 255)
        };
        egui::Color32::from_rgb(20, 22, 28).lerp_to_gamma(base, a)
    };
    let n = 64;
    for k in 0..n {
        let f = (f64::from(k) + 0.5) / f64::from(n);
        let s = if signed { 2.0 * f - 1.0 } else { f };
        let x0 = bar.left() + bar.width() * k as f32 / n as f32;
        let x1 = bar.left() + bar.width() * (k + 1) as f32 / n as f32;
        painter.rect_filled(
            egui::Rect::from_min_max(egui::pos2(x0, bar.top()), egui::pos2(x1, bar.bottom())),
            0.0,
            colour(s),
        );
    }
    // Ticks at values: position where the colour value of that field falls.
    let at = |v: f64| {
        let s = colour_value(v, sat, linear, scale);
        let f = if signed { 0.5 * (s + 1.0) } else { s };
        bar.left() + bar.width() * f as f32
    };
    let full = if linear { sat / scale } else { sat };
    let mut ticks: Vec<(f64, bool)> = Vec::new();
    if linear {
        for k in 1..=4 {
            ticks.push((full * f64::from(k) / 4.0, k == 4));
        }
    } else {
        let decades = game.field_range_decades.floor() as i32;
        for k in 0..=decades {
            ticks.push((sat * 10f64.powi(-k), k == 0));
        }
    }
    let text = ui.visuals().text_color();
    for (v, label) in ticks {
        for sign in if signed { vec![1.0, -1.0] } else { vec![1.0] } {
            let x = at(sign * v);
            painter.line_segment(
                [egui::pos2(x, bar.top()), egui::pos2(x, bar.bottom() + 3.0)],
                egui::Stroke::new(1.0, text),
            );
            if label && sign > 0.0 {
                painter.text(
                    egui::pos2(x.min(bar.right() - 30.0), bar.bottom() + 3.0),
                    egui::Align2::LEFT_TOP,
                    format!("{v:.1e}"),
                    egui::FontId::proportional(10.0),
                    text,
                );
            }
        }
    }
    let x0 = at(0.0);
    painter.line_segment(
        [
            egui::pos2(x0, bar.top()),
            egui::pos2(x0, bar.bottom() + 3.0),
        ],
        egui::Stroke::new(1.0, text),
    );
    ui.label(
        egui::RichText::new(if linear {
            format!(
                "Full colour at {} = {full:.2e} (the full field's scale ×10^−{:.1}); ticks every quarter.",
                name,
                game.field_gain_decades
            )
        } else {
            format!(
                "Full colour at {} = {sat:.2e}; ticks every decade below it.",
                name
            )
        })
        .small(),
    );
}

/// A radiation goal: what is required, what the flight radiates into it, and (with a
/// band) the spectrum over the goal's directions with the band marked. `system`: other
/// particles fly too, and the receiver sees all of them.
fn radiation_goal(
    ui: &mut egui::Ui,
    goal: &level::RadiationGoal,
    radiation: Option<f64>,
    spectrum: &[(f64, f64)],
    spectrum_top: f64,
    system: bool,
) {
    let [axis, half] = goal.direction;
    let [lo, hi] = goal.energy;
    let band = goal.band.map_or("all frequencies".to_string(), |[a, b]| {
        format!("ω {a:.3}–{b:.3}")
    });
    // A steady receiver measures a mean power (per unit time), the others an energy.
    let unit = if goal.window.is_some() {
        "per sr per unit time"
    } else {
        "per sr"
    };
    let need = if lo > 0.0 {
        format!("{} – {} {unit}", fmt_si(lo), fmt_si(hi))
    } else {
        format!("below {} {unit}", fmt_si(hi))
    };
    if let Some([t1, t2]) = goal.window {
        ui.label(
            egui::RichText::new(format!(
                "Steady receiver: into {axis:.0}° ± {half:.0}°, {band}, from t = {t1:.0} to {t2:.0}: {need}"
            ))
            .small(),
        )
        .on_hover_text(
            "The receiver (the band outside the arena) stands far away and collects the light the particle scatters while a wave drives it: the mean power it receives over its window, after the start's transient has died away (its time is the arrival time less the light time from the centre of the arena). The particle need not arrive anywhere; it must keep scattering, bound, until the window has passed.",
        );
        ui.label(
            egui::RichText::new(if system {
                "The receiver is far away and averages the scattered light over its window; the particle need not arrive. It sees every particle: their fields add."
            } else {
                "The receiver is far away and averages the scattered light over its window; the particle need not arrive."
            })
            .small()
            .italics(),
        );
    } else {
        ui.label(
            egui::RichText::new(format!(
                "Radiation goal: into {axis:.0}° ± {half:.0}°, {band}: {need}"
            ))
            .small(),
        )
        .on_hover_text(
            "The receiver (the band outside the arena) stands far away: it collects the particle's radiation, not the particle. The particle itself must still end in its detector. The radiation leaves at the speed of light whenever the particle is accelerated, mostly in the direction it is heading, and reaches the receiver after the flight: switch the map to 'particle field' to watch it go.",
        );
        ui.label(
            egui::RichText::new(if system {
                "The receiver is far away and collects radiation, not the particle; the particle must still end in its detector. It sees every particle: their fields add, in phase they reinforce, in antiphase they cancel."
            } else {
                "The receiver is far away and collects radiation, not the particle; the particle must still end in its detector."
            })
            .small()
            .italics(),
        );
    }
    if let Some(e) = radiation {
        let ok = e >= lo && e <= hi;
        let color = if ok {
            egui::Color32::from_rgb(90, 240, 110)
        } else {
            egui::Color32::from_rgb(255, 170, 80)
        };
        let text = if goal.window.is_some() {
            format!("Received: {} {unit} (mean power)", fmt_si(e))
        } else {
            format!("Radiated into it: {} {unit}", fmt_si(e))
        };
        ui.colored_label(color, egui::RichText::new(text).small());
    }
    if spectrum.is_empty() {
        return;
    }
    // The spectrum d²I/dωdΩ from 0 to twice the band's top (the band shaded), or without a
    // band to where the flight's spectrum has fallen off.
    let (w, h) = (ui.available_width().min(300.0), 70.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, h), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_gray(25));
    let omega_max = spectrum_top.max(1e-300);
    #[allow(clippy::cast_possible_truncation)]
    let x_of = |omega: f64| rect.left() + (omega / omega_max) as f32 * rect.width();
    if let Some([a, b]) = goal.band {
        painter.rect_filled(
            egui::Rect::from_x_y_ranges(x_of(a)..=x_of(b), rect.y_range()),
            0.0,
            egui::Color32::from_rgba_unmultiplied(90, 240, 110, 40),
        );
    }
    let peak = spectrum.iter().map(|s| s.1).fold(0.0, f64::max);
    if peak > 0.0 {
        #[allow(clippy::cast_possible_truncation)]
        let points: Vec<egui::Pos2> = spectrum
            .iter()
            .map(|&(omega, v)| {
                egui::pos2(
                    x_of(omega),
                    rect.bottom() - 2.0 - (v / peak) as f32 * (h - 6.0),
                )
            })
            .collect();
        painter.add(egui::Shape::line(
            points,
            egui::Stroke::new(1.5, egui::Color32::from_rgb(240, 220, 120)),
        ));
    }
    ui.label(
        egui::RichText::new(format!(
            "Spectrum into these directions, ω from 0 to {omega_max:.3}{}; peak {} {unit} per unit ω{}",
            if goal.band.is_some() { " (band shaded)" } else { " (all of it counts)" },
            fmt_si(peak),
            if goal.window.is_some() { ", over the window" } else { "" }
        ))
        .small(),
    );
}

/// A small plot of a circuit quantity over lab time (PHYSICS.md §2.10), its range
/// including 0, with lines at the times `marks` (the launches).
fn circuit_plot(ui: &mut egui::Ui, samples: &[(f64, f64)], marks: &[f64]) {
    let (Some(&(t0, _)), Some(&(t1, _))) = (samples.first(), samples.last()) else {
        return;
    };
    let width = ui.available_width().min(280.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 54.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let (lo, hi) = samples
        .iter()
        .fold((0.0_f64, 0.0_f64), |(a, b), s| (a.min(s.1), b.max(s.1)));
    let span = if hi > lo { hi - lo } else { 1.0 };
    let duration = if t1 > t0 { t1 - t0 } else { 1.0 };
    #[allow(clippy::cast_possible_truncation)]
    let at = |t: f64, v: f64| {
        egui::pos2(
            rect.left() + ((t - t0) / duration) as f32 * rect.width(),
            rect.bottom() - 2.0 - ((v - lo) / span) as f32 * (rect.height() - 14.0),
        )
    };
    let text = ui.visuals().text_color();
    painter.line_segment(
        [at(t0, 0.0), at(t1, 0.0)],
        egui::Stroke::new(1.0, egui::Color32::from_gray(90)),
    );
    for &m in marks {
        if (t0..=t1).contains(&m) {
            painter.line_segment(
                [at(m, lo), at(m, hi)],
                egui::Stroke::new(1.0, egui::Color32::from_rgb(230, 200, 80)),
            );
        }
    }
    let points: Vec<egui::Pos2> = samples.iter().map(|&(t, v)| at(t, v)).collect();
    painter.add(egui::Shape::line(
        points,
        egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255)),
    ));
    painter.text(
        rect.left_top(),
        egui::Align2::LEFT_TOP,
        format!(
            "max {}   min {}   t ≤ {}",
            fmt_si(hi),
            fmt_si(lo),
            fmt_si(t1)
        ),
        egui::FontId::proportional(10.0),
        text,
    );
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
