//! Sandbox mode: a level editor. Everything about the level can be changed; levels are
//! saved to `levels/custom/` (loaded by the game, but not part of the shipped-level
//! tests).

use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::mpsc::{Receiver, channel};

use bevy::prelude::Vec2;
use bevy_egui::egui;
use level::{
    Coil, Detector, ENGINE_VERSION, Element, ElementKind, FORMAT_VERSION, Grid, Launch, Level,
    Limits, Node, ParticleSpec, Region2, Shot, TolerancesSpec, WorldPhysics, solve,
};
use physics::trajectory::Outcome;

use crate::Game;
use crate::draw::to_vec2;
use crate::ui::{fmt_si, parse_si};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tool {
    LevelElement,
    PlayerElement,
    Launch,
    Aim,
    Detector,
    Region,
    CoilCircle,
    CoilRect,
    CoilPolygon,
}

pub struct SolverJob {
    rx: Mutex<Receiver<String>>,
    solution_rx: Mutex<Receiver<Vec<Element>>>,
}

pub struct Sandbox {
    pub active: bool,
    pub tool: Tool,
    /// Kind of level elements placed with the LevelElement tool.
    pub element_kind: ElementKind,
    /// Value of level charges placed with the LevelElement tool.
    pub charge_value: f64,
    /// Value μ of level magnets placed with the LevelElement tool.
    pub magnet_value: f64,
    /// Amplitude p₀ and orientation (degrees) of level antennas placed with the
    /// LevelElement tool.
    pub antenna_value: f64,
    pub antenna_angle: f64,
    /// Own frequency of level antennas placed with the tool (`None`: RF generator).
    pub antenna_omega: Option<f64>,
    /// Strength κ of coils drawn with the coil tools.
    pub coil_kappa: f64,
    pub file_name: String,
    pub overwrite: bool,
    /// First corner (or centre) of a box or coil being dragged.
    pub detector_drag: Option<Node>,
    /// Level element being dragged (its node when picked up).
    pub element_drag: Option<Node>,
    /// Text buffers of list fields.
    pub texts: crate::level_editor::EditTexts,
    /// Vertices of a polygon coil being drawn.
    pub pending_polygon: Vec<Node>,
    pub status: Vec<String>,
    pub solver: Option<SolverJob>,
    pub solver_report: Vec<String>,
    pub solver_solution: Option<Vec<Element>>,
}

impl Default for Sandbox {
    fn default() -> Self {
        Self {
            active: false,
            tool: Tool::LevelElement,
            element_kind: ElementKind::Charge,
            charge_value: 1e6,
            magnet_value: 10.0,
            antenna_value: 1e6,
            antenna_angle: 90.0,
            antenna_omega: None,
            coil_kappa: 1.0,
            file_name: "my_level".into(),
            overwrite: false,
            detector_drag: None,
            element_drag: None,
            texts: crate::level_editor::EditTexts::default(),
            pending_polygon: Vec::new(),
            status: Vec::new(),
            solver: None,
            solver_report: Vec::new(),
            solver_solution: None,
        }
    }
}

fn default_shot() -> Shot {
    Shot {
        particle: ParticleSpec {
            charge: 1e-6,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        launch: Launch {
            node: [0, 10, 0],
            direction: [1.0, 0.0, 0.0],
            kinetic_energy: 0.5,
            time: 0.0,
        },
        detector: Detector {
            min: [27, 8, 0],
            max: [30, 12, 0],
            acceptance: None,
        },
        beam: None,
    }
}

/// A blank level to start from.
pub fn empty_level() -> Level {
    Level {
        format_version: FORMAT_VERSION,
        engine_version: ENGINE_VERSION.to_string(),
        name: "New level".into(),
        description: String::new(),
        grid: Grid {
            nx: 30,
            ny: 20,
            nz: 0,
            subdivision: 1,
        },
        physics: WorldPhysics {
            c: Some(5.0),
            charge_radius: 0.3,
            magnet_radius: 0.3,
            wire_radius: 0.1,
            antenna_radius: 0.3,
            rf_omega: 0.0,
            radiation_reaction: false,
            beam_interaction: false,
            beam_retarded: false,
            t_max: 200.0,
            tolerances: TolerancesSpec {
                preview: 1e-10,
                verify: 1e-12,
            },
        },
        shots: vec![default_shot()],
        elements: Vec::new(),
        coils: Vec::new(),
        limits: Limits {
            max_charges: 2,
            magnitudes: vec![1e6, 2e6, 4e6],
            allow_positive: true,
            allow_negative: true,
            max_magnets: 0,
            magnet_strengths: vec![],
            max_antennas: 0,
            antenna_amplitudes: vec![],
            antenna_omegas: vec![],
            continuous: false,
            region: None,
            max_plates: 0,
            plate_voltages: vec![],
            plate: level::PlateSize::default(),
            supply_voltages: vec![],
        },
        reference_solution: Vec::new(),
        disturbances: Vec::new(),
        conductors: Vec::new(),
        electrodes: Vec::new(),
        gates: Vec::new(),
    }
}

pub fn custom_dir() -> PathBuf {
    crate::levels_dir().join("custom")
}

fn sync_texts(game: &mut Game) {
    use crate::level_editor::list_to_text;
    let limits = &game.editor.base().limits;
    game.sandbox.texts.magnitudes = list_to_text(&limits.magnitudes);
    game.sandbox.texts.magnet_strengths = list_to_text(&limits.magnet_strengths);
    game.sandbox.texts.antenna_amplitudes = list_to_text(&limits.antenna_amplitudes);
    game.sandbox.texts.antenna_omegas = list_to_text(&limits.antenna_omegas);
    game.sandbox.texts.plate_voltages = list_to_text(&limits.plate_voltages);
    game.sandbox.texts.supply_voltages = list_to_text(&limits.supply_voltages);
}

/// Enters sandbox mode, editing the level currently loaded.
pub fn enter(game: &mut Game) {
    game.sandbox.active = true;
    sync_texts(game);
    game.sandbox.file_name = slug(&game.editor.base().name);
    game.editor.edit_level(|_| {});
}

fn slug(name: &str) -> String {
    let s: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect();
    let s = s.trim_matches('_').to_string();
    if s.is_empty() { "my_level".into() } else { s }
}

/// Index of the coil whose wire is nearest to `p` (within 0.6 cells).
fn coil_near(level: &Level, p: Vec2) -> Option<usize> {
    let grid = level.grid;
    level
        .coils
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let d = match c {
                Coil::Circle { center, radius, .. } => {
                    #[allow(clippy::cast_possible_truncation)]
                    let r = *radius as f32;
                    (p.distance(to_vec2(grid.position(*center))) - r).abs()
                }
                Coil::Polygon { vertices, .. } => {
                    let v: Vec<Vec2> = vertices
                        .iter()
                        .map(|n| to_vec2(grid.position(*n)))
                        .collect();
                    (0..v.len())
                        .map(|k| {
                            let (a, b) = (v[k], v[(k + 1) % v.len()]);
                            let t = ((p - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                            p.distance(a + (b - a) * t)
                        })
                        .fold(f32::INFINITY, f32::min)
                }
            };
            (i, d)
        })
        .filter(|&(_, d)| d < 0.6)
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(i, _)| i)
}

/// Handles a click (or drag end) on grid node `node` (world point `world`) with the
/// current tool. Returns true if the click was consumed (player-element clicks are left to
/// the normal editor).
/// A sandbox click at world position `world`. Tools that edit the level work on its own
/// grid: a refined play grid is left first, so that the clicked node is in the level's
/// units (`Editor::edit_level` would reset it anyway, and a node of the finer grid would
/// then name a different place).
pub fn pointer_at(
    game: &mut Game,
    world: Vec2,
    pressed: bool,
    released: bool,
    right: bool,
) -> bool {
    if game.sandbox.tool != Tool::PlayerElement && game.editor.refinement() != 1 {
        game.editor.edit_level(|_| {});
    }
    let node = crate::nearest_node(game, world);
    if !game.editor.level.grid.contains(node) {
        return false;
    }
    game.editor.set_cursor(node);
    pointer(game, node, world, pressed, released, right)
}

fn pointer(
    game: &mut Game,
    node: Node,
    world: Vec2,
    pressed: bool,
    released: bool,
    right: bool,
) -> bool {
    let tool = game.sandbox.tool;
    let shot = game.active_shot;
    match tool {
        Tool::PlayerElement => false,
        Tool::LevelElement => {
            // Pressing on an existing level element picks it up; releasing it on another
            // free node moves it there (releasing where it was sets its value, as a click).
            let at = |n: Node| game.editor.base().elements.iter().position(|c| c.node == n);
            if pressed && at(node).is_some() {
                game.sandbox.element_drag = Some(node);
                return true;
            }
            if released && let Some(from) = game.sandbox.element_drag.take() {
                if from != node {
                    let blocked = at(node).is_some()
                        || game
                            .editor
                            .base()
                            .shots
                            .iter()
                            .any(|s| s.launch.node == node);
                    if !blocked {
                        game.editor.edit_level(|l| {
                            if let Some(e) = l.elements.iter_mut().find(|c| c.node == from) {
                                e.node = node;
                            }
                        });
                        game.editor.placement.retain(|c| c.node != node);
                    }
                    return true;
                }
            } else if !(pressed || right) {
                return true;
            }
            let kind = game.sandbox.element_kind;
            let (value, angle_deg, omega) = match kind {
                ElementKind::Charge => (game.sandbox.charge_value, 0.0, None),
                ElementKind::Magnet => (game.sandbox.magnet_value, 0.0, None),
                ElementKind::Antenna => (
                    game.sandbox.antenna_value,
                    game.sandbox.antenna_angle,
                    game.sandbox.antenna_omega,
                ),
                // Level electrodes are edited in the "Electrodes" section.
                ElementKind::Plate | ElementKind::Supply => return true,
            };
            game.editor.edit_level(|l| {
                let launch = l.shots.iter().any(|s| s.launch.node == node);
                let existing = l.elements.iter().position(|c| c.node == node);
                if right {
                    if let Some(i) = existing {
                        l.elements.remove(i);
                    }
                } else if !launch && value != 0.0 {
                    let e = Element {
                        node,
                        kind,
                        value,
                        angle_deg,
                        omega,
                    };
                    match existing {
                        Some(i) => l.elements[i] = e,
                        None => l.elements.push(e),
                    }
                }
            });
            // A player element on the same node would overlap.
            game.editor.placement.retain(|c| c.node != node);
            true
        }
        Tool::Launch => {
            if pressed {
                let occupied = game.editor.base().elements.iter().any(|c| c.node == node);
                if !occupied {
                    game.editor.edit_level(|l| l.shots[shot].launch.node = node);
                }
            }
            true
        }
        Tool::Aim => {
            if pressed {
                game.editor.edit_level(|l| {
                    let a = l.grid.position(l.shots[shot].launch.node);
                    let b = l.grid.position(node);
                    let d = b - a;
                    if d.length() > 0.0 {
                        let d = d.normalize();
                        l.shots[shot].launch.direction = [d.x, d.y, 0.0];
                    }
                });
            }
            true
        }
        Tool::CoilPolygon => {
            // A vertex clicked twice (or closing on the first) would be a side of zero
            // length: skipped.
            if pressed && game.sandbox.pending_polygon.last() != Some(&node) {
                game.sandbox.pending_polygon.push(node);
            }
            if right {
                let mut vertices = std::mem::take(&mut game.sandbox.pending_polygon);
                vertices.dedup();
                if vertices.len() > 1 && vertices.first() == vertices.last() {
                    vertices.pop();
                }
                if vertices.len() >= 3 {
                    let kappa = game.sandbox.coil_kappa;
                    game.editor
                        .edit_level(|l| l.coils.push(Coil::Polygon { vertices, kappa }));
                } else if let Some(i) = coil_near(game.editor.base(), world) {
                    game.editor.edit_level(|l| {
                        l.coils.remove(i);
                    });
                }
            }
            true
        }
        Tool::Detector | Tool::Region | Tool::CoilCircle | Tool::CoilRect => {
            if right && matches!(tool, Tool::CoilCircle | Tool::CoilRect) {
                if let Some(i) = coil_near(game.editor.base(), world) {
                    game.editor.edit_level(|l| {
                        l.coils.remove(i);
                    });
                }
                return true;
            }
            if pressed {
                game.sandbox.detector_drag = Some(node);
            }
            if released && let Some(a) = game.sandbox.detector_drag.take() {
                let min = [a[0].min(node[0]), a[1].min(node[1]), 0];
                let max = [a[0].max(node[0]), a[1].max(node[1]), 0];
                let non_empty = min[0] < max[0] && min[1] < max[1];
                let kappa = game.sandbox.coil_kappa;
                match tool {
                    Tool::Region => {
                        let region = Region2 { min, max };
                        game.editor.edit_level(|l| l.limits.region = Some(region));
                        game.editor.placement.retain(|c| region.contains(c.node));
                    }
                    Tool::Detector if non_empty => {
                        game.editor.edit_level(|l| {
                            let acceptance = l.shots[shot].detector.acceptance;
                            l.shots[shot].detector = Detector {
                                min,
                                max,
                                acceptance,
                            };
                        });
                    }
                    Tool::CoilCircle => {
                        let grid = game.editor.base().grid;
                        let radius = grid.position(a).distance(grid.position(node));
                        if radius > 0.5 {
                            game.editor.edit_level(|l| {
                                l.coils.push(Coil::Circle {
                                    center: a,
                                    radius,
                                    kappa,
                                });
                            });
                        }
                    }
                    Tool::CoilRect if non_empty => {
                        // Counter-clockwise, so positive κ means counter-clockwise current.
                        let vertices = vec![
                            [min[0], min[1], 0],
                            [max[0], min[1], 0],
                            [max[0], max[1], 0],
                            [min[0], max[1], 0],
                        ];
                        game.editor
                            .edit_level(|l| l.coils.push(Coil::Polygon { vertices, kappa }));
                    }
                    _ => {}
                }
            }
            true
        }
    }
}

/// Polls the background solver.
pub fn poll(game: &mut Game) {
    let Some(job) = &mut game.sandbox.solver else {
        return;
    };
    if let Ok(sol) = job.solution_rx.lock().unwrap().try_recv() {
        game.sandbox.solver_solution = Some(sol);
    }
    let mut done = false;
    let lines: Vec<String> = job.rx.lock().unwrap().try_iter().collect();
    for line in lines {
        if line == "__done__" {
            done = true;
            game.sandbox.solver_report.retain(|l| l != "Searching…");
        } else {
            game.sandbox.solver_report.push(line);
        }
    }
    if done
        && let Some(job) = game.sandbox.solver.take()
        && let Ok(sol) = job.solution_rx.lock().unwrap().try_recv()
    {
        game.sandbox.solver_solution = Some(sol);
    }
}

fn start_solver(game: &mut Game) {
    let level = game.editor.base().clone();
    let (tx, rx) = channel::<String>();
    let (sol_tx, sol_rx) = channel::<Vec<Element>>();
    game.sandbox.solver_report = vec!["Searching…".into()];
    game.sandbox.solver_solution = None;
    std::thread::spawn(move || {
        let (score, outcome) = solve::objective(&level, &[]);
        let _ = tx.send(format!(
            "Without player elements: {outcome:?} (summed distance to the detectors {score:.2})"
        ));
        let mut best: Option<Vec<Element>> = None;
        let singles = solve::single_element_solutions(&level);
        let _ = tx.send(format!("Verified 1-element solutions: {}", singles.len()));
        if let Some(c) = singles.first() {
            best = Some(vec![*c]);
        }
        let max = (level.limits.max_charges + level.limits.max_magnets).min(4) as usize;
        for k in 2..=max {
            let found = solve::anneal(&level, k, 32, 300, 0x5EED + k as u64);
            let _ = tx.send(format!(
                "Verified {k}-element solutions found: {}",
                found.len()
            ));
            if best.is_none() {
                best = found.into_iter().next();
            }
        }
        if let Some(b) = best {
            let _ = sol_tx.send(b);
        }
        let _ = tx.send("__done__".into());
    });
    game.sandbox.solver = Some(SolverJob {
        rx: Mutex::new(rx),
        solution_rx: Mutex::new(sol_rx),
    });
}

/// Checks and writes the level to `levels/custom/<file_name>.json`.
fn save(game: &mut Game) {
    let mut level = game.editor.base().clone();
    level.engine_version = ENGINE_VERSION.to_string();
    let mut status = Vec::new();
    if level.reference_solution.is_empty() {
        status.push("No reference solution stored (use a solver result or your elements).".into());
    } else {
        let mut all_ok = true;
        for (i, v) in level
            .verify_flights(&level.reference_solution)
            .iter()
            .enumerate()
        {
            if !(v.outcome() == Outcome::Arrived && v.status.is_verified()) {
                all_ok = false;
                status.push(format!(
                    "Warning: shot {}: reference solution not verified ({:?}, {:?}).",
                    i + 1,
                    v.outcome(),
                    v.status
                ));
            }
            let tr = &v.verified;
            let rad = tr.radiated_energy / tr.kinetic_initial;
            if rad > 1e-10 && !level.physics.radiation_reaction {
                status.push(format!(
                    "Warning: shot {}: neglected radiation is {rad:.1e} × T₀ (> 1e-10): \
                     physically inconsistent. Use a smaller particle charge with larger \
                     fixed charges.",
                    i + 1
                ));
            }
        }
        if all_ok {
            status.push("Reference solution: verified for every shot.".into());
        }
    }
    for issue in level.model_issues() {
        status.push(format!("Warning: {issue}."));
    }
    if let Some((_, cost)) = &game.cost {
        for m in cost.over_budget() {
            status.push(format!(
                "Warning: {} {} is over its limit ({}): see Computational cost.",
                m.label.to_lowercase(),
                m.text(),
                level::cost::format_value(m.limit, m.unit)
            ));
        }
    }
    let dir = custom_dir();
    let path = dir.join(format!("{}.json", slug(&game.sandbox.file_name)));
    if path.exists() && !game.sandbox.overwrite {
        status.push(format!(
            "{} exists; tick 'overwrite' to replace it.",
            path.display()
        ));
        game.sandbox.status = status;
        return;
    }
    let result =
        std::fs::create_dir_all(&dir).and_then(|()| std::fs::write(&path, level.to_json() + "\n"));
    match result {
        Ok(()) => {
            status.push(format!("Saved {}", path.display()));
            game.reload_levels(Some(&path));
        }
        Err(e) => status.push(format!("Could not save: {e}")),
    }
    game.sandbox.status = status;
}

fn tool_help(ui: &mut egui::Ui, game: &mut Game) {
    let small = |ui: &mut egui::Ui, t: &str| {
        ui.label(egui::RichText::new(t).small());
    };
    match game.sandbox.tool {
        Tool::LevelElement => {
            ui.horizontal(|ui| {
                ui.selectable_value(
                    &mut game.sandbox.element_kind,
                    ElementKind::Charge,
                    "charge",
                );
                ui.selectable_value(
                    &mut game.sandbox.element_kind,
                    ElementKind::Magnet,
                    "magnet",
                );
                ui.selectable_value(
                    &mut game.sandbox.element_kind,
                    ElementKind::Antenna,
                    "antenna",
                );
                let value = match game.sandbox.element_kind {
                    ElementKind::Charge => &mut game.sandbox.charge_value,
                    ElementKind::Magnet => &mut game.sandbox.magnet_value,
                    ElementKind::Antenna | ElementKind::Plate | ElementKind::Supply => {
                        &mut game.sandbox.antenna_value
                    }
                };
                let r = ui.add(
                    egui::DragValue::new(value)
                        .speed(1e4)
                        .custom_formatter(|v, _| fmt_si(v))
                        .custom_parser(parse_si),
                );
                game.text_focus |= r.has_focus();
                if ui.button("±").clicked() {
                    *value = -*value;
                }
            });
            small(
                ui,
                "Left click: place or set value; drag an element to move it. Right click: \
                 remove. Magnet value μ = μ₀m/4π \
                 (positive: moment out of the plane; in the plane B_z = −μ/r³). Antenna: \
                 dipole amplitude p₀ along the orientation, oscillating at the level's RF ω.",
            );
            if game.sandbox.element_kind == ElementKind::Antenna {
                ui.horizontal(|ui| {
                    ui.label("orientation");
                    let r = ui.add(
                        egui::DragValue::new(&mut game.sandbox.antenna_angle)
                            .speed(1.0)
                            .suffix("°"),
                    );
                    game.text_focus |= r.has_focus();
                    let mut own = game.sandbox.antenna_omega.is_some();
                    ui.checkbox(&mut own, "own ω")
                        .on_hover_text("Otherwise the antenna runs at the level's RF generator");
                    match (own, game.sandbox.antenna_omega) {
                        (true, None) => game.sandbox.antenna_omega = Some(1.0),
                        (false, Some(_)) => game.sandbox.antenna_omega = None,
                        _ => {}
                    }
                    if let Some(w) = &mut game.sandbox.antenna_omega {
                        let r = ui.add(egui::DragValue::new(w).speed(0.01).range(0.0..=1e12));
                        game.text_focus |= r.has_focus();
                    }
                });
            }
        }
        Tool::PlayerElement => {
            small(
                ui,
                "Places elements as a player would (with the limits below).",
            );
        }
        Tool::Launch => small(
            ui,
            "Click a node to move the launch point of the active shot.",
        ),
        Tool::Aim => small(
            ui,
            "Click a point: the active shot's launch direction points at it.",
        ),
        Tool::Detector => small(
            ui,
            "Active shot's detector: press on one corner, release on the opposite one.",
        ),
        Tool::Region => {
            small(ui, "Drag the box where players may place elements.");
            if ui.button("Remove region (place anywhere)").clicked() {
                game.editor.edit_level(|l| l.limits.region = None);
            }
        }
        Tool::CoilCircle | Tool::CoilRect | Tool::CoilPolygon => {
            ui.horizontal(|ui| {
                ui.label("κ = μ₀I/4π:");
                let r = ui.add(
                    egui::DragValue::new(&mut game.sandbox.coil_kappa)
                        .speed(0.01)
                        .custom_formatter(|v, _| fmt_si(v))
                        .custom_parser(parse_si),
                );
                game.text_focus |= r.has_focus();
                if ui.button("±").clicked() {
                    game.sandbox.coil_kappa = -game.sandbox.coil_kappa;
                }
            });
            small(
                ui,
                if game.sandbox.tool == Tool::CoilCircle {
                    "Press at the centre, release at the radius (exact radius: Coils list below). Right click on a wire: remove."
                } else if game.sandbox.tool == Tool::CoilPolygon {
                    "Click the vertices in order; right click closes the polygon (at least 3 vertices). Right click with no pending vertices removes a coil."
                } else {
                    "Drag the rectangle. Right click on a wire: remove. Positive κ: \
                     counter-clockwise current."
                },
            );
        }
    }
}

/// Sandbox section of the side panel.
#[allow(clippy::too_many_lines)]
/// Resource meters of the current setup (`level::cost`): measured times, steps, memory
/// and metal accuracy against their budgets, on logarithmic bars.
fn cost_meters(ui: &mut egui::Ui, game: &Game) {
    let Some((revision, cost)) = &game.cost else {
        ui.label("Measuring…");
        return;
    };
    if *revision != game.sent_revision {
        ui.label(
            egui::RichText::new("Updating… (values of the previous setup)")
                .small()
                .color(egui::Color32::YELLOW),
        );
    }
    for m in cost.meters() {
        ui.horizontal(|ui| {
            ui.add_sized([80.0, 16.0], egui::Label::new(m.label))
                .on_hover_text(m.hint);
            let color = match m.load() {
                level::cost::Load::Fine => egui::Color32::from_rgb(80, 170, 90),
                level::cost::Load::High => egui::Color32::from_rgb(220, 170, 50),
                level::cost::Load::Over => egui::Color32::from_rgb(220, 70, 60),
            };
            let (rect, response) =
                ui.allocate_exact_size(egui::vec2(150.0, 14.0), egui::Sense::hover());
            response.on_hover_text(m.hint);
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(35));
            #[allow(clippy::cast_possible_truncation)]
            let at = |f: f64| rect.left() + f as f32 * rect.width();
            if !m.value.is_nan() {
                let bar =
                    egui::Rect::from_x_y_ranges(rect.left()..=at(m.fraction()), rect.y_range());
                painter.rect_filled(bar, 2.0, color);
            }
            let (good, limit) = m.marks();
            for (f, gray) in [(good, 150), (limit, 230)] {
                painter.line_segment(
                    [
                        egui::pos2(at(f), rect.top()),
                        egui::pos2(at(f), rect.bottom()),
                    ],
                    egui::Stroke::new(1.0, egui::Color32::from_gray(gray)),
                );
            }
            ui.label(m.text());
        });
    }
    let per_eval =
        level::cost::format_value(cost.seconds_per_evaluation(), level::cost::Unit::Seconds);
    let mut info = format!(
        "{} flight(s), {per_eval} per force evaluation",
        cost.flights
    );
    if cost.setup.unknowns > 0 {
        info += &format!(", {} metal unknowns", cost.setup.unknowns);
    }
    ui.label(egui::RichText::new(info).small());
    ui.label(
        egui::RichText::new(
            "Bars are logarithmic; the marks are the budget (grey) and the limit (white). \
             Times are measured on this computer and extrapolated to all flights while \
             they are computed.",
        )
        .small()
        .weak(),
    );
}

/// Applies an edited copy `l` of the level (sandbox panel). A refinement of the level's
/// grid by `refine_by` keeps every position: the level (edited copy), and the player's
/// elements, whose nodes are scaled too.
fn commit_level_edit(game: &mut Game, mut l: Level, before: &Level, refine_by: Option<u32>) {
    if let Some(f) = refine_by {
        l.refine(f, &mut []);
    }
    if l == *before {
        return;
    }
    // Keep everything inside the (possibly resized) grid.
    let m = l.grid.max_node();
    let clamp = |n: Node| [n[0].clamp(0, m[0]), n[1].clamp(0, m[1]), 0];
    let grid = l.grid;
    l.elements.retain(|c| grid.contains(c.node));
    l.reference_solution.retain(|c| grid.contains(c.node));
    for s in &mut l.shots {
        s.launch.node = clamp(s.launch.node);
        s.detector.min = clamp(s.detector.min);
        s.detector.max = clamp(s.detector.max);
    }
    // `edit_level` first returns the player's elements to the level's own grid (as it
    // was before the edit); a refinement then scales them with the level.
    game.editor.edit_level(|target| *target = l);
    if let Some(f) = refine_by {
        let f = i64::from(f);
        for e in &mut game.editor.placement {
            e.node = e.node.map(|v| v * f);
        }
        game.editor.cursor = game.editor.cursor.map(|v| v * f);
    }
    game.editor.placement.retain(|c| grid.contains(c.node));
}

pub fn panel(ui: &mut egui::Ui, game: &mut Game) {
    ui.label(
        egui::RichText::new("Sandbox: level editor")
            .strong()
            .size(16.0),
    );
    ui.label(format!("Editing: {}", game.editor.base().name));
    if let Err(problem) = crate::level_editor::check_editable(game.editor.base()) {
        ui.colored_label(
            egui::Color32::YELLOW,
            format!("⚠ Not fully reproducible with the editor: {problem}"),
        );
    }
    ui.horizontal(|ui| {
        if ui.button("New empty level").clicked() {
            let l = empty_level();
            game.sandbox.file_name = slug(&l.name);
            game.sandbox.solver_report.clear();
            game.sandbox.solver_solution = None;
            game.sandbox.status.clear();
            game.editor = crate::editor::Editor::new(l);
            game.active_shot = 0;
            sync_texts(game);
        }
        if ui.button("Play-test").clicked() {
            game.sandbox.active = false;
        }
    });

    // Shots.
    let n = game.editor.base().shots.len();
    ui.horizontal(|ui| {
        ui.label(format!("Shot {} of {n}", game.active_shot + 1));
        if ui.button("◀").clicked() {
            game.cycle_shot(-1);
        }
        if ui.button("▶").clicked() {
            game.cycle_shot(1);
        }
        if ui
            .button("+ copy")
            .on_hover_text("Add a shot (copy of the active one)")
            .clicked()
        {
            let s = game.editor.base().shots[game.active_shot];
            game.editor.edit_level(|l| l.shots.push(s));
            game.active_shot = n;
        }
        if n > 1 && ui.button("− remove").clicked() {
            let i = game.active_shot;
            game.editor.edit_level(|l| {
                l.shots.remove(i);
            });
            game.active_shot = i.saturating_sub(1);
        }
    });

    ui.label("Tool (click on the map):");
    ui.horizontal_wrapped(|ui| {
        let t = &mut game.sandbox.tool;
        ui.selectable_value(t, Tool::LevelElement, "Level element");
        ui.selectable_value(t, Tool::PlayerElement, "Player element");
        ui.selectable_value(t, Tool::Launch, "Launch point");
        ui.selectable_value(t, Tool::Aim, "Aim");
        ui.selectable_value(t, Tool::Detector, "Detector (drag)");
        ui.selectable_value(t, Tool::Region, "Player region (drag)");
        ui.selectable_value(t, Tool::CoilCircle, "Coil ○ (drag)");
        ui.selectable_value(t, Tool::CoilRect, "Coil □ (drag)");
        ui.selectable_value(t, Tool::CoilPolygon, "Coil polygon (click)");
    });
    tool_help(ui, game);
    ui.separator();
    egui::CollapsingHeader::new("Computational cost")
        .default_open(true)
        .show(ui, |ui| cost_meters(ui, game));
    ui.separator();

    // Every field of the level (see level_editor: completeness invariant). Edit a copy,
    // then apply if anything changed.
    let mut l = game.editor.base().clone();
    let before = l.clone();
    let result =
        crate::level_editor::edit_level(ui, &mut l, &mut game.sandbox.texts, game.active_shot);
    game.text_focus |= result.focus;
    if !l.limits.allow_positive && !l.limits.allow_negative {
        l.limits.allow_positive = true;
    }
    commit_level_edit(game, l, &before, result.refine_by);
    ui.separator();

    // Solutions.
    ui.label(egui::RichText::new("Solution").strong());
    let reference = game.editor.base().reference_solution.len();
    ui.label(format!(
        "Stored reference solution: {reference} element(s)."
    ));
    ui.horizontal_wrapped(|ui| {
        if ui.button("Use my elements as reference").clicked() {
            let p = game.editor.placement.clone();
            game.editor.edit_level(|l| l.reference_solution = p);
        }
        if ui.button("Show reference").clicked() {
            game.editor.placement = game.editor.base().reference_solution.clone();
            game.editor.edit_level(|_| {});
        }
        let busy = game.sandbox.solver.is_some();
        if ui
            .add_enabled(
                !busy,
                egui::Button::new(if busy {
                    "Solving…"
                } else {
                    "Check solvability"
                }),
            )
            .clicked()
        {
            start_solver(game);
        }
    });
    for line in &game.sandbox.solver_report {
        ui.label(egui::RichText::new(line).small());
    }
    if let Some(sol) = game.sandbox.solver_solution.clone() {
        ui.horizontal(|ui| {
            ui.label(format!("Found: {} element(s)", sol.len()));
            if ui.button("Show").clicked() {
                game.editor.placement = sol.clone();
                game.editor.edit_level(|_| {});
            }
            if ui.button("Store as reference").clicked() {
                game.editor.edit_level(|l| l.reference_solution = sol);
            }
        });
    }
    ui.separator();

    // Saving.
    ui.label(egui::RichText::new("Save").strong());
    ui.horizontal(|ui| {
        ui.label("levels/custom/");
        let r =
            ui.add(egui::TextEdit::singleline(&mut game.sandbox.file_name).desired_width(120.0));
        game.text_focus |= r.has_focus();
        ui.label(".json");
    });
    ui.horizontal(|ui| {
        ui.checkbox(&mut game.sandbox.overwrite, "overwrite");
        if ui.button("Save level").clicked() {
            save(game);
        }
    });
    for s in &game.sandbox.status {
        ui.label(egui::RichText::new(s).small());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn game() -> Game {
        let level = Level::from_json(include_str!("../../../levels/01_first_bend.json")).unwrap();
        let mut g = Game::new(vec![level], vec![PathBuf::from("01_first_bend.json")]);
        enter(&mut g);
        g.sandbox.tool = Tool::LevelElement;
        g
    }

    fn click(g: &mut Game, x: f32, y: f32, right: bool) {
        pointer_at(g, Vec2::new(x, y), !right, false, right);
        pointer_at(g, Vec2::new(x, y), false, !right, false);
    }

    /// Regression: after refining the play grid, level elements could not be removed
    /// (the clicked node was in the finer grid's units, the level in its own), and new
    /// ones landed at twice their position.
    #[test]
    fn level_edits_after_a_grid_change_hit_the_clicked_place() {
        let mut g = game();
        let n0 = g.editor.base().elements.len();
        click(&mut g, 5.0, 5.0, false);
        assert_eq!(g.editor.base().elements.len(), n0 + 1);
        g.editor.set_refinement(2);
        click(&mut g, 5.0, 5.0, true);
        assert_eq!(
            g.editor.base().elements.len(),
            n0,
            "removed on the finer grid"
        );
        g.editor.set_refinement(3);
        click(&mut g, 7.0, 3.0, false);
        let e = g.editor.base().elements.last().unwrap();
        let pos = g.editor.base().grid.position(e.node);
        assert_eq!((pos.x, pos.y), (7.0, 3.0));
    }

    /// Changing "Nodes per cell" refines the level; the player's elements keep their
    /// positions too (they used to keep their node numbers and jump).
    #[test]
    fn refining_the_level_keeps_the_player_elements_in_place() {
        let mut g = game();
        g.editor.set_cursor([5, 5, 0]);
        g.editor.place().unwrap();
        let before_pos = g.editor.level.grid.position(g.editor.placement[0].node);
        let l = g.editor.base().clone();
        commit_level_edit(&mut g, l.clone(), &l, Some(2));
        assert_eq!(g.editor.base().grid.subdivision, 2);
        assert_eq!(
            g.editor.level.grid.position(g.editor.placement[0].node),
            before_pos
        );
        // And it can still be removed where it is.
        g.editor.set_cursor(g.editor.placement[0].node);
        g.editor.remove();
        assert!(g.editor.placement.is_empty());
    }
}
