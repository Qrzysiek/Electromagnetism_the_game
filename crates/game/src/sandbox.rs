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
use physics::trajectory::{Outcome, RunSettings, run};
use physics::verify::verify;

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
    /// Strength κ of coils drawn with the coil tools.
    pub coil_kappa: f64,
    pub file_name: String,
    pub overwrite: bool,
    /// First corner (or centre) of a box or coil being dragged.
    pub detector_drag: Option<Node>,
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
            coil_kappa: 1.0,
            file_name: "my_level".into(),
            overwrite: false,
            detector_drag: None,
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
        },
        launch: Launch {
            node: [0, 10, 0],
            direction: [1.0, 0.0, 0.0],
            kinetic_energy: 0.5,
        },
        detector: Detector {
            min: [27, 8, 0],
            max: [30, 12, 0],
        },
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
            region: None,
        },
        reference_solution: Vec::new(),
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
pub fn pointer(
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
            if !(pressed || right) {
                return true;
            }
            let kind = game.sandbox.element_kind;
            let value = match kind {
                ElementKind::Charge => game.sandbox.charge_value,
                ElementKind::Magnet => game.sandbox.magnet_value,
            };
            game.editor.edit_level(|l| {
                let launch = l.shots.iter().any(|s| s.launch.node == node);
                let existing = l.elements.iter().position(|c| c.node == node);
                if right {
                    if let Some(i) = existing {
                        l.elements.remove(i);
                    }
                } else if !launch && value != 0.0 {
                    let e = Element { node, kind, value };
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
            if pressed {
                game.sandbox.pending_polygon.push(node);
            }
            if right {
                let vertices = std::mem::take(&mut game.sandbox.pending_polygon);
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
                        game.editor
                            .edit_level(|l| l.shots[shot].detector = Detector { min, max });
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
        for (i, scn) in level
            .scenarios(&level.reference_solution)
            .iter()
            .enumerate()
        {
            let v = verify(scn, level.tolerances());
            if !(v.outcome() == Outcome::Arrived && v.status.is_verified()) {
                all_ok = false;
                status.push(format!(
                    "Warning: shot {}: reference solution not verified ({:?}, {:?}).",
                    i + 1,
                    v.outcome(),
                    v.status
                ));
            }
            let tr = run(
                scn,
                &RunSettings::with_tolerance(level.physics.tolerances.verify),
            );
            let rad = tr.radiated_energy / tr.kinetic_initial;
            if rad > 1e-10 {
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
                let value = match game.sandbox.element_kind {
                    ElementKind::Charge => &mut game.sandbox.charge_value,
                    ElementKind::Magnet => &mut game.sandbox.magnet_value,
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
                "Left click: place or set value. Right click: remove. Magnet value μ = μ₀m/4π \
                 (positive: moment out of the plane; in the plane B_z = −μ/r³).",
            );
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
                    "Press at the centre, release at the radius (exact radius: Coils list                      below). Right click on a wire: remove."
                } else if game.sandbox.tool == Tool::CoilPolygon {
                    "Click the vertices in order; right click closes the polygon (at least 3                      vertices). Right click with no pending vertices removes a coil."
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

    // Every field of the level (see level_editor: completeness invariant). Edit a copy,
    // then apply if anything changed.
    let mut l = game.editor.base().clone();
    let before = l.clone();
    let result =
        crate::level_editor::edit_level(ui, &mut l, &mut game.sandbox.texts, game.active_shot);
    game.text_focus |= result.focus;
    if let Some(f) = result.refine_by {
        l.refine(f, &mut []);
    }
    if !l.limits.allow_positive && !l.limits.allow_negative {
        l.limits.allow_positive = true;
    }
    if l != before {
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
        game.editor.edit_level(|target| *target = l);
        game.editor.placement.retain(|c| grid.contains(c.node));
    }
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
