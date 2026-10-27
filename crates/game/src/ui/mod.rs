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

mod controls;
mod format;
mod labels;
mod plots;
mod results;
mod tooltips;
mod views;
use controls::{cut_palette, free_palette, hardcore_controls, plate_palette, supplies};
pub(crate) use format::fmt_num;
use format::{angle_decimals, boundary_text, energy_text, outcome_text, percent, truncate};
pub use format::{disturbance_text, fmt_potential, fmt_si, parse_si};
use labels::{hall_labels, level_doc_buttons, reality, supply_badges};
pub use plots::PlotLayout;
use plots::{circuit_plot, potential_plots};
use results::{beam_result, energy_bar, exact_verdict, shot_color32, tube_result};
use tooltips::{inspect_tooltip, obstacle_name};
use views::{beam_energy, field_view_controls, radiation_goal};

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

pub fn panel(
    mut contexts: EguiContexts,
    mut game: ResMut<Game>,
    radiation: Res<crate::radiation::RadiationView>,
    camera: Single<(&Camera, &GlobalTransform), With<Camera2d>>,
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
    supply_badges(ctx, game, *camera, left);
    hall_labels(ctx, game, *camera, left);
    crate::legend::board(ctx, game);
    inspect_tooltip(ctx, game, *camera);
    Ok(())
}

/// The colour of a supply's badge and of its electrodes' highlight while it is in focus
/// (none of the status colours: yellow unverified, orange lost, white computing, violet gate).
pub const SUPPLY_FOCUS: egui::Color32 = egui::Color32::from_rgb(90, 240, 255);
/// The tunable electrodes' frame (and an idle badge's ring).
const SUPPLY_FRAME: egui::Color32 = egui::Color32::from_rgb(255, 217, 77);

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
    reality(ui, &level.reality);

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
        explain(
            ui,
            "It must work under every disturbance",
            "Outside fields: the setup must work under every disturbance. All flights of a \
             shot are drawn; the selected one is brightest.",
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
    } else if let Some(t) = &level.tube {
        if t.saturation == Some(0.0) {
            ui.label(format!(
                "{}: floating and cold (it emits nothing)",
                level.electrode_name(t.cathode)
            ));
        } else if t.metal_of(t.cathode).is_none() {
            ui.label(format!(
                "Carriers: q/m = {} (electrons; Newtonian), emitted by {}",
                fmt_si(t.charge_per_mass),
                level.electrode_ref(t.cathode)
            ));
        }
        for (k, m) in t.metals.iter().enumerate() {
            let contacts: Vec<String> = m.contacts.iter().map(|e| (e + 1).to_string()).collect();
            let joined = match contacts.len() {
                0 => "isolated (no contacts)".to_owned(),
                1 => format!("contact: electrode {}", contacts[0]),
                _ => format!("contacts: electrodes {}", contacts.join(", ")),
            };
            // Game units, plainly (an SI prefix would read 0.25 as 250m).
            let plain = |v: f64| {
                let text = format!("{v:.3}");
                text.trim_end_matches('0').trim_end_matches('.').to_owned()
            };
            ui.label(format!(
                "{}: carriers q/m = {}, charge density {}, thermal speed {}, mean free \
                 path {}; {joined}",
                level.metal_name(k),
                plain(m.charge_per_mass),
                plain(m.density),
                plain(m.thermal_speed),
                plain(m.mean_free_path),
            ));
        }
    }
    // What the numbers mean (`scales.rs`): the decisive ratios in the header, each quantity's
    // real-world scale inside.
    // The run's mean electrode potentials over the later half of its frames so far (the
    // fields across the metals; past a supply's start-up, the Amplifier's drain settling
    // from 120 to ~45).
    let potentials = game
        .tube_preview
        .as_ref()
        .filter(|v| v.revision == game.sent_revision && !v.run.frames.is_empty())
        .map(|v| {
            let later = &v.run.frames[v.run.frames.len() / 2..];
            #[allow(clippy::cast_precision_loss)]
            let n = later.len() as f64;
            let mut sum = vec![0.0; later[0].potentials.len()];
            for f in later {
                for (s, p) in sum.iter_mut().zip(&f.potentials) {
                    *s += f64::from(*p);
                }
            }
            sum.into_iter().map(|s| s / n).collect::<Vec<f64>>()
        });
    if let Some(sc) = crate::scales::of(&level, &game.editor.placement, potentials.as_deref()) {
        egui::CollapsingHeader::new(format!("Scales: {}", sc.headline))
            .id_salt("scales")
            // Opened for developer captures with EM_SCALES=1.
            .default_open(std::env::var("EM_SCALES").is_ok_and(|v| v == "1"))
            .show(ui, |ui| {
                for l in &sc.lines {
                    ui.label(egui::RichText::new(l.quantity).strong().small())
                        .on_hover_text(l.detail);
                    ui.label(egui::RichText::new(format!("here: {}", l.here)).small());
                    ui.label(
                        egui::RichText::new(format!("real world: {}", l.real))
                            .small()
                            .weak(),
                    );
                }
            })
            .header_response
            .on_hover_text(
                "What the level's numbers mean: the ratios that decide what you see (the same \
                 at any scale), and where such things happen in the real world.",
            );
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
    if !game.sandbox.active {
        level_doc_buttons(ui, game);
    }
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
                            circuit_plot(ui, &shown, &launches, "");
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
        ElementKind::Cut => "cuts",
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
    } else if game.editor.kind == ElementKind::Cut {
        cut_palette(ui, game);
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
                (
                    ElementKind::Plate | ElementKind::Supply | ElementKind::Free | ElementKind::Cut,
                    _,
                ) => "",
            };
            let hover = match kind {
                ElementKind::Charge => "Flip sign (S)",
                ElementKind::Magnet => {
                    "Flip orientation (S): moment out of the plane (+z) or into it (−z)"
                }
                ElementKind::Antenna => "Flip phase (S): opposite phase of the RF generator",
                ElementKind::Plate | ElementKind::Supply | ElementKind::Free | ElementKind::Cut => {
                    ""
                }
            };
            ui.label(match kind {
                ElementKind::Charge => "New charge:",
                ElementKind::Magnet => "New magnet μ:",
                ElementKind::Antenna => "New antenna p₀:",
                ElementKind::Plate | ElementKind::Supply | ElementKind::Free | ElementKind::Cut => {
                    ""
                }
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
        explain(
            ui,
            "Why: the fields are singular there",
            "Metal bodies must stay apart, and particles must not start in metal or on a \
             wire: the fields are singular there.",
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
        // The z-invariant world's maps, in one row as every level's (V cycles them in this
        // order).
        ui.horizontal_wrapped(|ui| {
            ui.label("Map (V):")
                .on_hover_text("What the background shows; V cycles through the maps");
            ui.selectable_value(&mut game.map, Some(MapMode::Potential), "potential")
                .on_hover_text("The carriers' potential energy at the time shown");
            let self_field = level
                .tube
                .as_ref()
                .is_some_and(|t| t.speed_of_light.is_some());
            ui.selectable_value(&mut game.map, Some(MapMode::Magnetic), "magnetic B")
                .on_hover_text(if self_field {
                    "B_z of the carriers' own currents, the field that acts back on them in \
                     this level (the Darwin limit: their inductance, the skin effect; PHYSICS.md \
                     §2.12), as the run computed it. Your solenoid's field is not in it: uniform \
                     and much stronger, it would swamp the map; the ⊙/⊗ marks round the board \
                     show its way."
                } else {
                    "B_z that the run's own charges and currents make (the current through \
                     the metal and along its leads, the beams), with retardation, computed as \
                     the animation plays. Your solenoid's field is not in it: uniform and \
                     hundreds of times stronger, it would swamp the map; the ⊙/⊗ marks round \
                     the board show its way. So the map flips with the current, not with the \
                     solenoid."
                });
            ui.selectable_value(&mut game.map, Some(MapMode::Electric), "electric field")
                .on_hover_text("|E| and its direction at the time shown");
            ui.selectable_value(&mut game.map, Some(MapMode::EnergyFlow), "energy flow")
                .on_hover_text("The Poynting vector of that field: where its energy goes");
            if level.tube.as_ref().is_some_and(|t| !t.metals.is_empty()) {
                ui.selectable_value(&mut game.map, Some(MapMode::Current), "current")
                    .on_hover_text(
                        "The current density j inside the metals: where the charge flows and \
                         how strongly (colour |j|, arrows its way)",
                    );
                ui.selectable_value(&mut game.map, Some(MapMode::Heat), "heat")
                    .on_hover_text(
                        "The power per area the field gives the metals' carriers, E·j: where \
                         the heat is made",
                    );
            }
            ui.selectable_value(&mut game.map, None, "off")
                .on_hover_text("No background map");
        });
        ui.horizontal(|ui| {
            use crate::TubeParticles;
            ui.label("Particles:");
            ui.selectable_value(&mut game.tube_particles, TubeParticles::Density, "density")
                .on_hover_text(
                    "The density of their charge, a glow on the map (electrons light blue, \
                     ions orange): a beam of many macroparticles reads as one stream",
                );
            ui.selectable_value(&mut game.tube_particles, TubeParticles::Dots, "dots")
                .on_hover_text("Each macroparticle (a line of many carriers) as a dot");
            ui.selectable_value(&mut game.tube_particles, TubeParticles::Hidden, "off")
                .on_hover_text("The map alone");
        });
        ui.add(
            egui::Slider::new(&mut game.map_tau, 0.0..=20.0)
                .logarithmic(true)
                .smallest_positive(0.1)
                .text("average τ"),
        )
        .on_hover_text(
            "A meter's view: the maps averaged over time with this time constant (a \
             first-order low-pass, as an RC meter or a power meter integrates). 0: the \
             instantaneous field, every fluctuation of the real physics.",
        );
        if matches!(game.map, Some(MapMode::Current | MapMode::Heat)) {
            ui.add(
                egui::Slider::new(&mut game.tube_range_decades, 0.5..=4.0).text("range (decades)"),
            )
            .on_hover_text("How many decades below the run's strongest so far are still coloured");
            let span = crate::metal_map_span(game.map_tau);
            let (headline, what) = if game.map == Some(MapMode::Current) {
                (
                    "|j| in the metals; arrows: its way",
                    "the current density |j| (arrows: the way positive charge goes; electrons \
                     drift against them)",
                )
            } else {
                (
                    "E·j: orange heats the carriers, teal takes back",
                    "the power per area the field gives the carriers, E·j, orange (teal where \
                     it takes energy back): in a steady state the heat made there; in an \
                     isolated metal the induced field's part, which the heat goals measure",
                )
            };
            explain(
                ui,
                headline,
                &format!(
                    "Inside the metals: {what}, on a log scale over {:.1} decades below the \
                     strongest so far on this level (the 95th percentile of the metals' cells: \
                     one setup compares with another, a cut core shows dimmer). Each \
                     cell's mean over the last {} time units (twice the average τ when one is \
                     set) and over the 3 × 3 cells about it that no cut separates from it. \
                     Under a swinging field, over its last whole half periods (up to four): \
                     the heat as a mean over whole cycles of it, the current by lock-in (its \
                     part in step with the field's change, which a plain mean would cancel: \
                     the rms, its way while the field rises). The carriers' thermal motion \
                     makes every cell's value noisy (Johnson–Nyquist noise): the span is cut \
                     into parts whose spread gives each value's standard error, and a cell \
                     fades out as its value sinks into its noise (gone within one standard \
                     error, full beyond three). Hover a metal for the values.",
                    game.tube_range_decades,
                    fmt_num(span),
                ),
            );
            ui.separator();
            controls(ui);
            return;
        }
        if matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow)) {
            ui.add(
                egui::Slider::new(&mut game.tube_range_decades, 0.5..=4.0).text("range (decades)"),
            )
            .on_hover_text("How many decades below the run's strongest field are still coloured");
            let (c, t, delay) = game
                .em_snapshot
                .as_ref()
                .map_or((0.0, 0.0, 0.0), |s| (s.c, s.t, s.delay));
            let (headline, what) = if game.map == Some(MapMode::Magnetic) {
                (
                    "B_z: orange out of the plane, teal into it",
                    "B_z, out of the plane (orange) or into it (teal)",
                )
            } else {
                (
                    "|S|, the field's energy flow; arrows: its way",
                    "the energy flow |S| = (c/4π)|E × B| (arrows: its direction)",
                )
            };
            let open = level.tube.as_ref().is_some_and(|t| t.output.is_some());
            let smoothed = if delay <= 0.0 {
                String::new()
            } else if open {
                format!(
                    " The field is shown smoothed over about {}: the run's currents are known \
                     at its frames only and step from one to the next, and the steps would \
                     radiate waves at the frames' rate.",
                    fmt_num(delay)
                )
            } else {
                format!(
                    " The tube's envelope is a closed, lossless box: the field is shown \
                     smoothed over about {} (a few of the box's light-crossing times), the \
                     field of the run's currents; unsmoothed, the waves its macroparticles \
                     radiate (far stronger than a real tube's electrons') would stay in the \
                     box and cover the map.",
                    fmt_num(delay)
                )
            };
            explain(
                ui,
                headline,
                &format!(
                    "The electromagnetic field of the run's charges and currents, with \
                     retardation at the speed of light c = {} (a display parameter: the \
                     carriers' own motion stays Newtonian): {what}, on a log scale over {:.1} \
                     decades below the run's strongest so far (the board's 99th \
                     percentile); field at t = {}. Every \
                     electrode is a perfect conductor; each one's circuit current comes \
                     along a straight lead from the enclosing grounded ring or the board's \
                     edge (a convention: a z-invariant tube's circuit runs in the \
                     plane).{smoothed}",
                    fmt_num(c),
                    game.tube_range_decades,
                    fmt_num(t),
                ),
            );
            ui.separator();
            controls(ui);
            return;
        }
        if game.map == Some(MapMode::Electric) {
            ui.add(
                egui::Slider::new(&mut game.tube_range_decades, 0.5..=4.0).text("range (decades)"),
            )
            .on_hover_text("How many decades below the strong field are still coloured");
            explain(
                ui,
                "|E|, log scale; arrows: its direction",
                &format!(
                    "|E| at the time shown (live while animating), with the electrons' space \
                     charge: bright where strong, on a log scale over {:.1} decades below the \
                     field at the 99th percentile of the board; contours every quarter \
                     decade; arrows show its direction (a force on electrons the opposite \
                     way). Inside metal it is 0.{}",
                    game.tube_range_decades,
                    if level.tube.as_ref().is_some_and(|t| t.output.is_some()) {
                        " This is the run's own, instantaneous field: the waves the tube \
                         radiates, which take time to travel, show on the magnetic and energy \
                         flow maps."
                    } else {
                        ""
                    }
                ),
            );
            ui.separator();
            controls(ui);
            return;
        }
        if game.map == Some(MapMode::Potential) {
            ui.add(
                egui::Slider::new(&mut game.tube_gain_decades, 0.0..=2.0).text("gain (decades)"),
            )
            .on_hover_text(
                "Amplifies the colours: full colour at 10^−gain of the largest potential \
                 difference; small differences (the space charge's dip) become visible",
            );
        }
        let metals = level.tube.as_ref().is_some_and(|t| !t.metals.is_empty());
        if game.map == Some(MapMode::Potential) {
            explain(
                ui,
                if metals {
                    "Carriers' potential energy: blue downhill, red uphill"
                } else {
                    "Electrons' potential energy: blue downhill, red uphill"
                },
                if metals {
                    "The carriers' potential energy at the time shown (live), relative to the \
                 first contact then, in units of the largest potential difference the \
                 electrodes reach divided by the gain: blue downhill, red uphill; contours \
                 every quarter. It includes the carriers' charge and the ions' background: \
                 inside a metal carrying a current they cancel but for a thin layer at its \
                 surfaces, and the potential falls evenly from contact to contact (a uniform \
                 field: Ohm's law)."
                } else {
                    "The electrons' potential energy at the time shown (live), relative to the \
                 cathode then, in units of the largest potential difference to it the \
                 electrodes reach (their supplies over a period) divided by the gain: blue \
                 downhill (towards the anode), red uphill, smoothly graded and saturating \
                 beyond about twice the unit; contours every quarter. It includes the \
                 electrons' own space charge and the charge it induces on the electrodes: \
                 the flat region in front of the cathode, where the field is pulled down to \
                 zero, is the space charge that limits the current."
                },
            );
        }
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
    let legend: (&str, &str) = match game.map {
        Some(MapMode::Potential) if time_dependent => (
            "Red uphill, blue downhill (the static part)",
            "Red: uphill for this shot's particle, blue: downhill (the static part of the \
             field; plates and coils driven by circuits as they are at t = 0). No dark \
             region: the time-dependent fields (a ramped or driven coil's induced field, \
             driven plates, antennas, waves) do work on the particle, so energy \
             conservation forbids nothing here.",
        ),
        Some(MapMode::Potential) => (
            "Red uphill, blue downhill; dark: out of reach",
            "Red: uphill for this shot's particle, blue: downhill; contours every T₀/4. \
             Dark: forbidden by energy conservation for every particle shown (all shots \
             with \"show all\"; a beam counts with its most energetic particle). The \
             flight's own static field, with the disturbance shown and magnets: exact for \
             charges, clouds, coils, antennas and stray fields; metal from its picture \
             model (spheres to about 1e-4, electrodes on a coarser mesh, within about 0.5 % \
             of their potential). With interacting beam particles only a guide, since \
             they exchange energy.",
        ),
        Some(MapMode::Magnetic) if shot.particle.charge == 0.0 && shot.particle.moment != 0.0 => (
            "B_z: orange out of the plane, teal into it",
            "B perpendicular to the plane. Orange: out of the plane, teal: into it. \
             Value 1 = field in which this particle's magnetic energy |m B| equals T₀; \
             contours every 0.25.",
        ),
        Some(MapMode::Magnetic) => (
            "B_z: orange out of the plane, teal into it",
            "B perpendicular to the plane. Orange: out of the plane, teal: into it. \
             Value 1 = field in which this particle circles with a 5-cell radius; \
             contours every 0.25.",
        ),
        Some(MapMode::Waves) => (
            "Antennas' and waves' fields: colour B_z, arrows E",
            "Fields of the antennas and waves at the animation time (their exact retarded \
             fields, drawn in single precision). Colour: B perpendicular to the plane \
             (orange out, blue in), on the chosen scale (bar above); arrows: E. Near an \
             antenna the field is quasi-static; further out the radiation travels outwards \
             at c.",
        ),
        Some(MapMode::Total) => (
            "Everything's field: colour B_z or |E|, arrows E",
            "The total field at the animation time: level charges, magnets and coils, \
             antennas, waves and disturbances, and the particle's own field (retarded, \
             Liénard–Wiechert). The static part is computed in double precision on 6 \
             points per cell and interpolated between them (metal from its picture model); \
             antennas, waves and the particle are evaluated per pixel. The particle's field \
             is usually far weaker than the electrodes'; raise the range to see it. Colour: \
             B_z or |E| on the chosen scale (bar above); arrows: E.",
        ),
        Some(MapMode::ParticleField) if level.has_beams() && game.neglected_only => (
            "What the beam's interaction leaves out",
            "What the quasi-static beam interaction leaves out of the dynamics: the full \
             retarded field of all particles minus the fields it uses (each particle's past \
             continued along the motion it would have in the fields it feels now). Mostly \
             how those fields change along its path during the light travel time, and the \
             delay with which the news of an absorption spreads. Shown on the full field's \
             colour scale (linear by default), so its true size is seen; its size is stated \
             above.",
        ),
        Some(MapMode::ParticleField) if level.has_beams() => (
            "The beam's own field: colour B_z, arrows E",
            "The retarded (Liénard–Wiechert, exact) field of every particle of the beam at \
             the animation time, from the flight shown (the exact one once the verdict is \
             in). Every change of velocity sends out radiation at c. Before launch each \
             particle is taken to move with its launch acceleration. One entering its \
             detector flies on into the screening cup, its charge fading as seen from \
             outside (with the instant drain its field disappears as the news of its \
             absorption spreads at c); one absorbed by a body stays there at rest. Colour: \
             B perpendicular to the plane, on the chosen scale; arrows: E.",
        ),
        Some(MapMode::ParticleField) if shot.particle.moment != 0.0 => (
            "The particle's own field: colour B_z, arrows E",
            "The field of the particle itself at the animation time. Its charge's field is \
             Liénard–Wiechert (exact); its magnetic moment is drawn as the static dipole \
             field B_z = −m/r³ from its retarded position, without the terms of a moving \
             dipole (of order v/c). Colour: B perpendicular to the plane, on the chosen \
             scale; arrows: E. Before launch the particle is taken to move uniformly.",
        ),
        Some(MapMode::ParticleField) => (
            "The particle's own field: colour B_z, arrows E",
            "The field of the particle itself (Liénard–Wiechert, exact) at the animation \
             time. Every change of velocity sends out a radiation pulse at c; the energy it \
             carries is what the particle loses (radiation reaction, when included). \
             Colour: B perpendicular to the plane, on the chosen scale; arrows: E. Before \
             launch the particle is taken to move uniformly.",
        ),
        // Tube levels only (their own View section).
        Some(MapMode::Electric | MapMode::EnergyFlow | MapMode::Current | MapMode::Heat) | None => {
            ("", "")
        }
    };
    if !legend.0.is_empty() {
        explain(ui, legend.0, legend.1);
    }
    if game.map == Some(MapMode::ParticleField)
        && (!level.conductors.is_empty() || !level.electrodes.is_empty())
    {
        explain_coloured(
            ui,
            "Near metal: its induced charges not drawn",
            "Near metal: the charges the particle induces on it reshape its field there \
             (and screen it inside); that part is not drawn. (In the dynamics their force is \
             included for spheres and bounded for electrodes, PHYSICS.md §2.6–2.7.)",
            Some(egui::Color32::from_rgb(255, 210, 120)),
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
    explain(
        ui,
        "Lines show direction only",
        "In this 2D slice of a 3D field, lines show direction only, not strength. Screened \
         regions (behind grounded metal) get few lines.",
    );
    ui.separator();

    controls(ui);
}

/// An explanation that does not crowd the panel (the owner, 2026-10-19: no walls of small
/// print; every element explained where it is): a short line, marked ℹ, and the whole
/// explanation on hover.
pub(crate) fn explain(ui: &mut egui::Ui, headline: &str, detail: &str) {
    explain_coloured(ui, headline, detail, None);
}

/// `explain` in a colour (a warning).
fn explain_coloured(
    ui: &mut egui::Ui,
    headline: &str,
    detail: &str,
    colour: Option<egui::Color32>,
) {
    let mut text = egui::RichText::new(format!("{headline}  ℹ")).small();
    if let Some(c) = colour {
        text = text.color(c);
    }
    ui.label(text).on_hover_text(detail);
}

fn controls(ui: &mut egui::Ui) {
    ui.collapsing("Controls", |ui| {
        ui.label("Mouse: left click place, right click remove, wheel changes magnitude.");
        ui.label("Arrows move the cursor (Shift: ×5). Space/Enter place, Del/X remove.");
        ui.label("S flip sign/phase, Q/E change magnitude, C clear.");
        ui.label("1–4 grid refinement, [ ] switch shot, H show all shots.");
        ui.label("N/P next/previous level, V map, F field lines, A animation.");
        ui.label("Hold I: the probe, everything known at the pointer (fields, particles, …).");
        ui.label("R rotate antenna (Shift+R back), W antenna frequency, M cycles kinds.");
        ui.label("Drag your elements with the mouse, or G to grab / drop and Esc to cancel.");
        ui.label("Hardcore (checkbox): sliders instead of fixed values; Q/E and W step ×1.1.");
        ui.label(
            "Tube levels: the wheel zooms, Z or Home shows the whole board; drag with the \
             middle button, or the left one on empty board, or Ctrl + arrows to pan.",
        );
    });
}

/// A plot's line: its colour and its samples `(t, v)`.
type Line<'a> = ([u8; 3], &'a [(f64, f64)]);

/// The plot lines' colours, in the order of the electrodes plotted (the first the single
/// plots' blue).
const PLOT_COLOURS: [[u8; 3]; 7] = [
    [120, 200, 255],
    [255, 165, 80],
    [130, 225, 140],
    [235, 130, 220],
    [240, 220, 110],
    [110, 230, 210],
    [255, 120, 110],
];

#[cfg(test)]
mod tests;
