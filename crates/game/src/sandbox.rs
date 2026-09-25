//! Sandbox mode: a level editor. Everything about the level can be changed; levels are
//! saved to `levels/custom/` (loaded by the game, but not part of the shipped-level
//! tests).

use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::mpsc::{Receiver, channel};

use bevy_egui::egui;
use level::{
    Charge, Detector, ENGINE_VERSION, FORMAT_VERSION, Grid, Launch, Level, Limits, Node,
    ParticleSpec, TolerancesSpec, WorldPhysics, solve,
};
use physics::trajectory::{Outcome, RunSettings, run};
use physics::verify::verify;

use crate::Game;
use crate::ui::{fmt_si, parse_si};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tool {
    LevelCharge,
    PlayerCharge,
    Launch,
    Aim,
    Detector,
}

pub struct SolverJob {
    rx: Mutex<Receiver<String>>,
    solution_rx: Mutex<Receiver<Vec<Charge>>>,
}

pub struct Sandbox {
    pub active: bool,
    pub tool: Tool,
    /// Value of level charges placed with the LevelCharge tool.
    pub charge_value: f64,
    pub file_name: String,
    pub overwrite: bool,
    /// First corner of a detector being dragged.
    pub detector_drag: Option<Node>,
    /// Allowed magnitudes as edited text.
    pub magnitudes_text: String,
    pub status: Vec<String>,
    pub solver: Option<SolverJob>,
    pub solver_report: Vec<String>,
    pub solver_solution: Option<Vec<Charge>>,
}

impl Default for Sandbox {
    fn default() -> Self {
        Self {
            active: false,
            tool: Tool::LevelCharge,
            charge_value: 1e6,
            file_name: "my_level".into(),
            overwrite: false,
            detector_drag: None,
            magnitudes_text: String::new(),
            status: Vec::new(),
            solver: None,
            solver_report: Vec::new(),
            solver_solution: None,
        }
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
            t_max: 200.0,
            tolerances: TolerancesSpec {
                preview: 1e-10,
                verify: 1e-12,
            },
        },
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
        level_charges: Vec::new(),
        limits: Limits {
            max_charges: 2,
            magnitudes: vec![1e6, 2e6, 4e6],
            allow_positive: true,
            allow_negative: true,
        },
        reference_solution: Vec::new(),
    }
}

fn magnitudes_to_text(m: &[f64]) -> String {
    m.iter()
        .map(|v| format!("{v:e}"))
        .collect::<Vec<_>>()
        .join(", ")
}

pub fn custom_dir() -> PathBuf {
    crate::levels_dir().join("custom")
}

/// Enters sandbox mode, editing the level currently loaded.
pub fn enter(game: &mut Game) {
    game.sandbox.active = true;
    game.sandbox.magnitudes_text = magnitudes_to_text(&game.editor.base().limits.magnitudes);
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

/// Handles a click (or drag end) on grid node `node` with the current tool. Returns true
/// if the click was consumed (player-charge clicks are left to the normal editor).
pub fn pointer(game: &mut Game, node: Node, pressed: bool, released: bool, right: bool) -> bool {
    let tool = game.sandbox.tool;
    match tool {
        Tool::PlayerCharge => false,
        Tool::LevelCharge => {
            if !(pressed || right) {
                return true;
            }
            let value = game.sandbox.charge_value;
            let launch = game.editor.base().launch.node;
            game.editor.edit_level(|l| {
                let existing = l.level_charges.iter().position(|c| c.node == node);
                if right {
                    if let Some(i) = existing {
                        l.level_charges.remove(i);
                    }
                } else if node != launch && value != 0.0 {
                    match existing {
                        Some(i) => l.level_charges[i].charge = value,
                        None => l.level_charges.push(Charge {
                            node,
                            charge: value,
                        }),
                    }
                }
            });
            // A player charge on the same node would overlap.
            game.editor.placement.retain(|c| c.node != node);
            true
        }
        Tool::Launch => {
            if pressed {
                let occupied = game
                    .editor
                    .base()
                    .level_charges
                    .iter()
                    .any(|c| c.node == node);
                if !occupied {
                    game.editor.edit_level(|l| l.launch.node = node);
                }
            }
            true
        }
        Tool::Aim => {
            if pressed {
                game.editor.edit_level(|l| {
                    let a = l.grid.position(l.launch.node);
                    let b = l.grid.position(node);
                    let d = b - a;
                    if d.length() > 0.0 {
                        let d = d.normalize();
                        l.launch.direction = [d.x, d.y, 0.0];
                    }
                });
            }
            true
        }
        Tool::Detector => {
            if pressed {
                game.sandbox.detector_drag = Some(node);
            }
            if released && let Some(a) = game.sandbox.detector_drag.take() {
                let min = [a[0].min(node[0]), a[1].min(node[1]), 0];
                let max = [a[0].max(node[0]), a[1].max(node[1]), 0];
                if min[0] < max[0] && min[1] < max[1] {
                    game.editor
                        .edit_level(|l| l.detector = Detector { min, max });
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
    let (sol_tx, sol_rx) = channel::<Vec<Charge>>();
    game.sandbox.solver_report = vec!["Searching…".into()];
    game.sandbox.solver_solution = None;
    std::thread::spawn(move || {
        let (score, outcome) = solve::objective(&level, &[]);
        let _ = tx.send(format!(
            "Without player charges: {outcome:?} (closest approach to the detector {score:.2})"
        ));
        let mut best: Option<Vec<Charge>> = None;
        let singles = solve::single_charge_solutions(&level);
        let _ = tx.send(format!("Verified 1-charge solutions: {}", singles.len()));
        if let Some(c) = singles.first() {
            best = Some(vec![*c]);
        }
        for k in 2..=level.limits.max_charges.min(4) as usize {
            let found = solve::anneal(&level, k, 32, 300, 0x5EED + k as u64);
            let _ = tx.send(format!(
                "Verified {k}-charge solutions found: {}",
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
        status.push("No reference solution stored (use a solver result or your charges).".into());
    } else {
        let v = verify(
            &level.scenario(&level.reference_solution),
            level.tolerances(),
        );
        if v.outcome() == Outcome::Arrived && v.status.is_verified() {
            status.push("Reference solution: verified.".into());
        } else {
            status.push(format!(
                "Warning: reference solution not verified ({:?}, {:?}).",
                v.outcome(),
                v.status
            ));
        }
        let tr = run(
            &level.scenario(&level.reference_solution),
            &RunSettings::with_tolerance(level.physics.tolerances.verify),
        );
        let rad = tr.radiated_energy / tr.kinetic_initial;
        if rad > 1e-10 {
            status.push(format!(
                "Warning: neglected radiation is {rad:.1e} × T₀ (> 1e-10): physically inconsistent. \
                 Use a smaller particle charge with larger fixed charges."
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

/// Sandbox section of the side panel.
#[allow(clippy::too_many_lines)]
pub fn panel(ui: &mut egui::Ui, game: &mut Game) {
    ui.label(
        egui::RichText::new("Sandbox: level editor")
            .strong()
            .size(16.0),
    );
    ui.label(format!("Editing: {}", game.editor.base().name));
    ui.horizontal(|ui| {
        if ui.button("New empty level").clicked() {
            let l = empty_level();
            game.sandbox.magnitudes_text = magnitudes_to_text(&l.limits.magnitudes);
            game.sandbox.file_name = slug(&l.name);
            game.sandbox.solver_report.clear();
            game.sandbox.solver_solution = None;
            game.sandbox.status.clear();
            game.editor = crate::editor::Editor::new(l);
        }
        if ui.button("Play-test").clicked() {
            game.sandbox.active = false;
        }
    });

    ui.label("Tool (click on the map):");
    ui.horizontal_wrapped(|ui| {
        let t = &mut game.sandbox.tool;
        ui.selectable_value(t, Tool::LevelCharge, "Level charge");
        ui.selectable_value(t, Tool::PlayerCharge, "Player charge");
        ui.selectable_value(t, Tool::Launch, "Launch point");
        ui.selectable_value(t, Tool::Aim, "Aim");
        ui.selectable_value(t, Tool::Detector, "Detector (drag)");
    });
    match game.sandbox.tool {
        Tool::LevelCharge => {
            ui.horizontal(|ui| {
                ui.label("Charge value:");
                let r = ui.add(
                    egui::DragValue::new(&mut game.sandbox.charge_value)
                        .speed(1e4)
                        .custom_formatter(|v, _| fmt_si(v))
                        .custom_parser(parse_si),
                );
                game.text_focus |= r.has_focus();
                if ui.button("±").clicked() {
                    game.sandbox.charge_value = -game.sandbox.charge_value;
                }
            });
            ui.label(
                egui::RichText::new("Left click: place or set value. Right click: remove.").small(),
            );
        }
        Tool::PlayerCharge => {
            ui.label(
                egui::RichText::new("Places charges as a player would (with the limits below).")
                    .small(),
            );
        }
        Tool::Launch => {
            ui.label(egui::RichText::new("Click a node to move the launch point.").small());
        }
        Tool::Aim => {
            ui.label(
                egui::RichText::new("Click a point: the launch direction points at it.").small(),
            );
        }
        Tool::Detector => {
            ui.label(
                egui::RichText::new("Press on one corner, release on the opposite corner.").small(),
            );
        }
    }
    ui.separator();

    // Numeric properties: edit a copy, then apply if anything changed.
    let mut l = game.editor.base().clone();
    let before = l.clone();
    let mut focus = false;
    egui::Grid::new("sandbox_props")
        .num_columns(2)
        .show(ui, |ui| {
            ui.label("Name");
            focus |= ui.text_edit_singleline(&mut l.name).has_focus();
            ui.end_row();
            ui.label("Description");
            focus |= ui.text_edit_singleline(&mut l.description).has_focus();
            ui.end_row();
            ui.label("Grid (cells)");
            ui.horizontal(|ui| {
                focus |= ui
                    .add(egui::DragValue::new(&mut l.grid.nx).range(4..=200))
                    .has_focus();
                ui.label("×");
                focus |= ui
                    .add(egui::DragValue::new(&mut l.grid.ny).range(4..=200))
                    .has_focus();
            });
            ui.end_row();
            ui.label("Nodes per cell");
            let mut sub = l.grid.subdivision;
            focus |= ui
                .add(egui::DragValue::new(&mut sub).range(1..=8))
                .has_focus();
            if sub != l.grid.subdivision && sub.is_multiple_of(l.grid.subdivision) {
                let f = sub / l.grid.subdivision;
                l.refine(f, &mut []);
            }
            ui.end_row();
            ui.label("Speed of light c");
            ui.horizontal(|ui| {
                let mut newtonian = l.physics.c.is_none();
                ui.checkbox(&mut newtonian, "Newtonian");
                if newtonian {
                    l.physics.c = None;
                } else {
                    let mut c = l.physics.c.unwrap_or(5.0);
                    focus |= ui
                        .add(egui::DragValue::new(&mut c).speed(0.05).range(0.01..=1e9))
                        .has_focus();
                    l.physics.c = Some(c);
                }
            });
            ui.end_row();
            ui.label("Particle q, m");
            ui.horizontal(|ui| {
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut l.particle.charge)
                            .speed(1e-8)
                            .custom_formatter(|v, _| fmt_si(v))
                            .custom_parser(parse_si),
                    )
                    .has_focus();
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut l.particle.mass)
                            .speed(0.01)
                            .range(1e-12..=1e12),
                    )
                    .has_focus();
            });
            ui.end_row();
            ui.label("Launch energy T₀");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut l.launch.kinetic_energy)
                        .speed(0.01)
                        .range(1e-9..=1e9),
                )
                .has_focus();
            ui.end_row();
            ui.label("Launch angle (°)");
            let d = l.launch.direction;
            let mut angle = d[1].atan2(d[0]).to_degrees();
            if ui
                .add(
                    egui::DragValue::new(&mut angle)
                        .speed(0.5)
                        .range(-180.0..=180.0),
                )
                .changed()
            {
                let a = angle.to_radians();
                l.launch.direction = [a.cos(), a.sin(), 0.0];
            }
            ui.end_row();
            ui.label("Charge radius");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut l.physics.charge_radius)
                        .speed(0.01)
                        .range(0.01..=2.0),
                )
                .has_focus();
            ui.end_row();
            ui.label("Time limit");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut l.physics.t_max)
                        .speed(1.0)
                        .range(1.0..=1e7),
                )
                .has_focus();
            ui.end_row();
            ui.label("Player charges (max)");
            focus |= ui
                .add(egui::DragValue::new(&mut l.limits.max_charges).range(0..=20))
                .has_focus();
            ui.end_row();
            ui.label("Allowed magnitudes");
            let r = ui.text_edit_singleline(&mut game.sandbox.magnitudes_text);
            focus |= r.has_focus();
            if r.lost_focus() {
                let parsed: Vec<f64> = game
                    .sandbox
                    .magnitudes_text
                    .split(',')
                    .filter_map(|t| t.trim().parse::<f64>().ok())
                    .filter(|v| *v > 0.0 && v.is_finite())
                    .collect();
                if !parsed.is_empty() {
                    l.limits.magnitudes = parsed;
                }
                game.sandbox.magnitudes_text = magnitudes_to_text(&l.limits.magnitudes);
            }
            ui.end_row();
            ui.label("Signs allowed");
            ui.horizontal(|ui| {
                ui.checkbox(&mut l.limits.allow_positive, "+");
                ui.checkbox(&mut l.limits.allow_negative, "−");
            });
            ui.end_row();
        });
    game.text_focus |= focus;
    if !l.limits.allow_positive && !l.limits.allow_negative {
        l.limits.allow_positive = true;
    }
    if l != before {
        // Keep everything inside the (possibly resized) grid.
        let m = l.grid.max_node();
        let clamp = |n: Node| [n[0].clamp(0, m[0]), n[1].clamp(0, m[1]), 0];
        l.level_charges.retain(|c| l.grid.contains(c.node));
        l.reference_solution.retain(|c| l.grid.contains(c.node));
        l.launch.node = clamp(l.launch.node);
        l.detector.min = clamp(l.detector.min);
        l.detector.max = clamp(l.detector.max);
        let grid = l.grid;
        game.editor.edit_level(|target| *target = l);
        game.editor.placement.retain(|c| grid.contains(c.node));
    }
    ui.separator();

    // Solutions.
    ui.label(egui::RichText::new("Solution").strong());
    let reference = game.editor.base().reference_solution.len();
    ui.label(format!("Stored reference solution: {reference} charge(s)."));
    ui.horizontal_wrapped(|ui| {
        if ui.button("Use my charges as reference").clicked() {
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
            ui.label(format!("Found: {} charge(s)", sol.len()));
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
