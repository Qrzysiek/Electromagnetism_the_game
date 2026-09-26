//! Property editors for every part of a level (sandbox).
//!
//! **Completeness invariant:** everything a level file can contain must be editable here
//! (or by a sandbox tool), so any shipped level can be rebuilt by hand. Every level type
//! is destructured *exhaustively* (no `..`): adding a field to the level format fails to
//! compile until it gets a control here, or is explicitly marked as not user-editable
//! with a reason.

use bevy_egui::egui;
use level::{
    Coil, Detector, Disturbance, Element, ElementKind, Grid, Launch, Level, Limits, Node,
    ParticleSpec, Region2, Shot, TolerancesSpec, Wave, WorldPhysics,
};

use crate::ui::{fmt_si, parse_si};

// Ranges the editors allow; `check_editable` uses the same constants, so every shipped
// level is checked to lie within reach of the editors.
const GRID_CELLS: std::ops::RangeInclusive<u32> = 4..=200;
const SUBDIVISION: std::ops::RangeInclusive<u32> = 1..=8;
const MAX_POSITIVE: f64 = 1e12;
const MIN_POSITIVE: f64 = 1e-12;
const MAX_RADIUS: f64 = 5.0;
const MAX_COUNT: u32 = 20;
const TOLERANCE: std::ops::RangeInclusive<f64> = 1e-18..=1.0;

/// Text buffers for list fields (edited as comma-separated text).
#[derive(Default)]
pub struct EditTexts {
    pub magnitudes: String,
    pub magnet_strengths: String,
}

pub fn list_to_text(m: &[f64]) -> String {
    m.iter()
        .map(|v| format!("{v:e}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn parse_list(text: &str) -> Vec<f64> {
    text.split(',')
        .filter_map(|t| parse_si(t.trim()))
        .filter(|v| *v > 0.0 && v.is_finite())
        .collect()
}

/// Any finite value, with SI suffixes accepted when typed.
fn si(ui: &mut egui::Ui, v: &mut f64, speed: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .custom_formatter(|v, _| fmt_si(v))
            .custom_parser(parse_si),
    )
    .has_focus()
}

/// A small positive number shown in scientific notation (tolerances).
fn sci(ui: &mut egui::Ui, v: &mut f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(0.0)
            .range(TOLERANCE)
            .custom_formatter(|v, _| format!("{v:.0e}"))
            .custom_parser(|t| t.trim().parse::<f64>().ok()),
    )
    .has_focus()
}

fn positive(ui: &mut egui::Ui, v: &mut f64, speed: f64, hi: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .range(MIN_POSITIVE..=hi),
    )
    .has_focus()
}

/// A grid node (z shown only for 3D grids).
fn node(ui: &mut egui::Ui, n: &mut Node, grid: &Grid) -> bool {
    let m = grid.max_node();
    let mut focus = false;
    ui.horizontal(|ui| {
        focus |= ui
            .add(egui::DragValue::new(&mut n[0]).range(0..=m[0]).prefix("x "))
            .has_focus();
        focus |= ui
            .add(egui::DragValue::new(&mut n[1]).range(0..=m[1]).prefix("y "))
            .has_focus();
        if !grid.is_2d() {
            focus |= ui
                .add(egui::DragValue::new(&mut n[2]).range(0..=m[2]).prefix("z "))
                .has_focus();
        }
    });
    focus
}

fn row(ui: &mut egui::Ui, label: &str, add: impl FnOnce(&mut egui::Ui) -> bool) -> bool {
    ui.label(label);
    let f = ui.horizontal(add).inner;
    ui.end_row();
    f
}

fn edit_grid(ui: &mut egui::Ui, g: &mut Grid, refine_by: &mut Option<u32>) -> bool {
    let Grid {
        nx,
        ny,
        // 3D levels are not playable yet (SPEC stage 3); 2D levels have nz = 0.
        nz: _,
        subdivision,
    } = g;
    let mut focus = false;
    focus |= row(ui, "Grid (cells)", |ui| {
        let a = ui
            .add(egui::DragValue::new(nx).range(GRID_CELLS))
            .has_focus();
        ui.label("×");
        a | ui
            .add(egui::DragValue::new(ny).range(GRID_CELLS))
            .has_focus()
    });
    let current = *subdivision;
    let mut sub = current;
    focus |= row(ui, "Nodes per cell", |ui| {
        ui.add(egui::DragValue::new(&mut sub).range(SUBDIVISION))
            .has_focus()
    });
    // Refinement keeps every node in place, so only integer multiples are possible.
    if sub != current && sub.is_multiple_of(current) {
        *refine_by = Some(sub / current);
    }
    focus
}

fn edit_physics(ui: &mut egui::Ui, p: &mut WorldPhysics) -> bool {
    let WorldPhysics {
        c,
        charge_radius,
        magnet_radius,
        wire_radius,
        t_max,
        tolerances,
    } = p;
    let TolerancesSpec { preview, verify } = tolerances;
    let mut focus = false;
    focus |= row(ui, "Speed of light c", |ui| {
        let mut newtonian = c.is_none();
        ui.checkbox(&mut newtonian, "Newtonian");
        if newtonian {
            *c = None;
            false
        } else {
            let mut v = c.unwrap_or(5.0);
            let f = positive(ui, &mut v, 0.05, MAX_POSITIVE);
            *c = Some(v);
            f
        }
    });
    focus |= row(ui, "Radii: charge, magnet, wire", |ui| {
        positive(ui, charge_radius, 0.01, MAX_RADIUS)
            | positive(ui, magnet_radius, 0.01, MAX_RADIUS)
            | positive(ui, wire_radius, 0.01, MAX_RADIUS)
    });
    focus |= row(ui, "Time limit", |ui| {
        positive(ui, t_max, 1.0, MAX_POSITIVE)
    });
    focus |= row(ui, "Tolerance: preview, verify", |ui| {
        sci(ui, preview) | sci(ui, verify)
    });
    focus
}

fn edit_shot(ui: &mut egui::Ui, s: &mut Shot, grid: &Grid) -> bool {
    let Shot {
        particle,
        launch,
        detector,
    } = s;
    let ParticleSpec {
        charge,
        mass,
        radius,
    } = particle;
    let Launch {
        node: launch_node,
        direction,
        kinetic_energy,
    } = launch;
    let Detector { min, max } = detector;
    let mut focus = false;
    focus |= row(ui, "Particle q, m, radius", |ui| {
        si(ui, charge, 1e-8) | positive(ui, mass, 0.01, MAX_POSITIVE) | {
            ui.add(
                egui::DragValue::new(radius)
                    .speed(0.01)
                    .range(0.0..=MAX_RADIUS),
            )
            .has_focus()
        }
    });
    focus |= row(ui, "Launch node", |ui| node(ui, launch_node, grid));
    focus |= row(ui, "Launch energy T₀", |ui| {
        positive(ui, kinetic_energy, 0.01, MAX_POSITIVE)
    });
    focus |= row(ui, "Launch angle (°)", |ui| {
        let mut angle = direction[1].atan2(direction[0]).to_degrees();
        let r = ui.add(
            egui::DragValue::new(&mut angle)
                .speed(0.5)
                .range(-180.0..=180.0),
        );
        if r.changed() {
            let a = angle.to_radians();
            *direction = [a.cos(), a.sin(), 0.0];
        }
        r.has_focus()
    });
    focus |= row(ui, "Launch direction (x, y)", |ui| {
        // Exact components, e.g. to reproduce a file; normalized by the physics.
        ui.add(egui::DragValue::new(&mut direction[0]).speed(0.01))
            .has_focus()
            | ui.add(egui::DragValue::new(&mut direction[1]).speed(0.01))
                .has_focus()
    });
    focus |= row(ui, "Detector corner 1", |ui| node(ui, min, grid));
    focus |= row(ui, "Detector corner 2", |ui| node(ui, max, grid));
    focus
}

fn kind_combo(ui: &mut egui::Ui, id: usize, kind: &mut ElementKind) {
    let text = |k: ElementKind| match k {
        ElementKind::Charge => "charge",
        ElementKind::Magnet => "magnet μ",
    };
    egui::ComboBox::from_id_salt(("element_kind", id))
        .selected_text(text(*kind))
        .width(80.0)
        .show_ui(ui, |ui| {
            for k in [ElementKind::Charge, ElementKind::Magnet] {
                ui.selectable_value(kind, k, text(k));
            }
        });
}

fn edit_elements(ui: &mut egui::Ui, elements: &mut Vec<Element>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, e) in elements.iter_mut().enumerate() {
        let Element {
            node: n,
            kind,
            value,
        } = e;
        ui.horizontal(|ui| {
            kind_combo(ui, i, kind);
            focus |= si(ui, value, 1e4);
            focus |= node(ui, n, grid);
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        elements.remove(i);
    }
    focus
}

fn edit_coils(ui: &mut egui::Ui, coils: &mut Vec<Coil>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, c) in coils.iter_mut().enumerate() {
        ui.push_id(("coil", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Coil {}", i + 1));
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            match c {
                Coil::Circle {
                    center,
                    radius,
                    kappa,
                } => {
                    ui.horizontal(|ui| {
                        ui.label("circle, centre");
                        focus |= node(ui, center, grid);
                    });
                    ui.horizontal(|ui| {
                        ui.label("radius");
                        focus |= positive(ui, radius, 0.05, MAX_POSITIVE);
                        ui.label("κ");
                        focus |= si(ui, kappa, 100.0);
                    });
                }
                Coil::Polygon { vertices, kappa } => {
                    ui.horizontal(|ui| {
                        ui.label("polygon, κ");
                        focus |= si(ui, kappa, 100.0);
                    });
                    let mut drop = None;
                    for (k, v) in vertices.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("  vertex {}", k + 1));
                            focus |= node(ui, v, grid);
                            if ui.small_button("×").clicked() {
                                drop = Some(k);
                            }
                        });
                    }
                    if let Some(k) = drop
                        && vertices.len() > 3
                    {
                        vertices.remove(k);
                    }
                    if ui.small_button("+ vertex").clicked() {
                        let last = *vertices.last().expect("polygon has vertices");
                        vertices.push(last);
                    }
                }
            }
        });
    }
    if let Some(i) = remove {
        coils.remove(i);
    }
    focus
}

fn edit_disturbances(ui: &mut egui::Ui, list: &mut Vec<Disturbance>) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, d) in list.iter_mut().enumerate() {
        let Disturbance { name, e, bz, waves } = d;
        ui.push_id(("disturbance", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("{}", i + 1));
                focus |= ui
                    .add(egui::TextEdit::singleline(name).desired_width(120.0))
                    .has_focus();
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("stray E");
                focus |= si(ui, &mut e[0], 1e-3);
                focus |= si(ui, &mut e[1], 1e-3);
                ui.label("B_z");
                focus |= si(ui, bz, 1e-3);
            });
            let mut drop = None;
            for (k, w) in waves.iter_mut().enumerate() {
                let Wave {
                    amplitude,
                    direction_deg,
                    omega,
                    phase_deg,
                } = w;
                ui.horizontal(|ui| {
                    ui.label("  wave E₀");
                    focus |= si(ui, amplitude, 1e-3);
                    ui.label("ω");
                    focus |= ui
                        .add(
                            egui::DragValue::new(omega)
                                .speed(0.01)
                                .range(0.0..=MAX_POSITIVE),
                        )
                        .has_focus();
                    if ui.small_button("×").clicked() {
                        drop = Some(k);
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("    towards");
                    focus |= ui
                        .add(egui::DragValue::new(direction_deg).speed(1.0).suffix("°"))
                        .has_focus();
                    ui.label("phase");
                    focus |= ui
                        .add(egui::DragValue::new(phase_deg).speed(1.0).suffix("°"))
                        .has_focus();
                });
            }
            if let Some(k) = drop {
                waves.remove(k);
            }
            if ui.small_button("+ wave").clicked() {
                waves.push(Wave {
                    amplitude: 0.01,
                    direction_deg: 0.0,
                    omega: 1.0,
                    phase_deg: 0.0,
                });
            }
        });
        ui.separator();
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ disturbance").clicked() {
        list.push(Disturbance::default());
    }
    focus
}

fn edit_limits(ui: &mut egui::Ui, l: &mut Limits, texts: &mut EditTexts, grid: &Grid) -> bool {
    let Limits {
        max_charges,
        magnitudes,
        allow_positive,
        allow_negative,
        max_magnets,
        magnet_strengths,
        region,
    } = l;
    let mut focus = false;
    focus |= row(ui, "Player charges (max)", |ui| {
        ui.add(egui::DragValue::new(max_charges).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Allowed charge magnitudes", |ui| {
        let r = ui.text_edit_singleline(&mut texts.magnitudes);
        if r.lost_focus() {
            let parsed = parse_list(&texts.magnitudes);
            if !parsed.is_empty() {
                *magnitudes = parsed;
            }
            texts.magnitudes = list_to_text(magnitudes);
        }
        r.has_focus()
    });
    row(ui, "Charge signs allowed", |ui| {
        ui.checkbox(allow_positive, "+");
        ui.checkbox(allow_negative, "−");
        false
    });
    focus |= row(ui, "Player magnets (max)", |ui| {
        ui.add(egui::DragValue::new(max_magnets).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Allowed magnet strengths μ", |ui| {
        let r = ui.text_edit_singleline(&mut texts.magnet_strengths);
        if r.lost_focus() {
            *magnet_strengths = parse_list(&texts.magnet_strengths);
            texts.magnet_strengths = list_to_text(magnet_strengths);
        }
        r.has_focus()
    });
    let mut restricted = region.is_some();
    row(ui, "Player region", |ui| {
        ui.checkbox(&mut restricted, "restrict placement");
        false
    });
    match (restricted, region.is_some()) {
        (true, false) => {
            *region = Some(Region2 {
                min: [0, 0, 0],
                max: grid.max_node(),
            });
        }
        (false, true) => *region = None,
        _ => {}
    }
    if let Some(Region2 { min, max }) = region {
        focus |= row(ui, "  region corner 1", |ui| node(ui, min, grid));
        focus |= row(ui, "  region corner 2", |ui| node(ui, max, grid));
    }
    focus
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
        grid,
        physics,
        shots,
        elements,
        coils,
        limits,
        // Edited in the Solution section (store own elements or a solver result).
        reference_solution: _,
        disturbances,
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
    egui::CollapsingHeader::new(format!("Disturbances ({})", disturbances.len()))
        .show(ui, |ui| focus |= edit_disturbances(ui, disturbances));
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

/// Checks that every value of `level` can be produced with the sandbox editors (the same
/// ranges and constraints as the controls above, and the same exhaustive destructuring).
/// Returns the first problem found.
pub fn check_editable(level: &Level) -> Result<(), String> {
    let Level {
        format_version: _,
        engine_version: _,
        name: _,
        description: _,
        grid,
        physics,
        shots,
        elements,
        coils,
        limits,
        reference_solution,
        disturbances,
    } = level;
    let Grid {
        nx,
        ny,
        nz,
        subdivision,
    } = grid;
    let fail = |what: &str| Err(what.to_string());
    if !GRID_CELLS.contains(nx) || !GRID_CELLS.contains(ny) {
        return fail("grid size outside the editor's range");
    }
    if *nz != 0 {
        return fail("3D grids cannot be edited yet");
    }
    if !SUBDIVISION.contains(subdivision) {
        return fail("subdivision outside the editor's range");
    }
    let positive = |v: f64| (MIN_POSITIVE..=MAX_POSITIVE).contains(&v);
    let radius = |v: f64| (MIN_POSITIVE..=MAX_RADIUS).contains(&v);
    let WorldPhysics {
        c,
        charge_radius,
        magnet_radius,
        wire_radius,
        t_max,
        tolerances: TolerancesSpec { preview, verify },
    } = physics;
    if c.is_some_and(|c| !positive(c))
        || !radius(*charge_radius)
        || !radius(*magnet_radius)
        || !radius(*wire_radius)
        || !positive(*t_max)
        || !TOLERANCE.contains(preview)
        || !TOLERANCE.contains(verify)
    {
        return fail("physics value outside the editor's range");
    }
    let on_grid = |n: &Node| grid.contains(*n);
    for (i, s) in shots.iter().enumerate() {
        let Shot {
            particle:
                ParticleSpec {
                    charge,
                    mass,
                    radius: r,
                },
            launch:
                Launch {
                    node,
                    direction,
                    kinetic_energy,
                },
            detector: Detector { min, max },
        } = s;
        if !charge.is_finite()
            || !positive(*mass)
            || !(0.0..=MAX_RADIUS).contains(r)
            || !on_grid(node)
            || !direction.iter().all(|d| d.is_finite())
            || direction[2] != 0.0
            || !positive(*kinetic_energy)
            || !on_grid(min)
            || !on_grid(max)
        {
            return Err(format!(
                "shot {} has a value outside the editor's range",
                i + 1
            ));
        }
    }
    for e in elements.iter().chain(reference_solution) {
        let Element { node, kind, value } = e;
        match kind {
            ElementKind::Charge | ElementKind::Magnet => {}
        }
        if !on_grid(node) || !value.is_finite() {
            return fail("element outside the grid or with an invalid value");
        }
    }
    for c in coils {
        match c {
            Coil::Circle {
                center,
                radius: r,
                kappa,
            } => {
                if !on_grid(center) || !positive(*r) || !kappa.is_finite() {
                    return fail("circular coil outside the editor's range");
                }
            }
            Coil::Polygon { vertices, kappa } => {
                if vertices.len() < 3 || !vertices.iter().all(on_grid) || !kappa.is_finite() {
                    return fail(
                        "polygon coil not reproducible (needs 3 or more vertices on the grid)",
                    );
                }
            }
        }
    }
    if disturbances.len() > MAX_COUNT as usize {
        return fail("too many disturbances");
    }
    for d in disturbances {
        let Disturbance {
            name: _,
            e,
            bz,
            waves,
        } = d;
        let waves_ok = waves.iter().all(|w| {
            let Wave {
                amplitude,
                direction_deg,
                omega,
                phase_deg,
            } = w;
            amplitude.is_finite()
                && direction_deg.is_finite()
                && (0.0..=MAX_POSITIVE).contains(omega)
                && phase_deg.is_finite()
        });
        if !e.iter().all(|v| v.is_finite()) || !bz.is_finite() || !waves_ok {
            return fail("disturbance outside the editor's range");
        }
    }
    let Limits {
        max_charges,
        magnitudes,
        allow_positive: _,
        allow_negative: _,
        max_magnets,
        magnet_strengths,
        region,
    } = limits;
    if *max_charges > MAX_COUNT
        || *max_magnets > MAX_COUNT
        || !magnitudes
            .iter()
            .chain(magnet_strengths)
            .all(|v| *v > 0.0 && v.is_finite())
    {
        return fail("limits outside the editor's range");
    }
    if let Some(Region2 { min, max }) = region
        && !(on_grid(min) && on_grid(max))
    {
        return fail("player region outside the grid");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every shipped (and custom) level can be rebuilt by hand in the sandbox: all its
    /// values are within reach of the editors. New fields are covered at compile time by
    /// the exhaustive destructuring above.
    #[test]
    fn all_levels_are_reproducible_in_the_sandbox() {
        let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
        let mut checked = 0;
        for d in [dir.clone(), dir.join("custom")] {
            let Ok(entries) = std::fs::read_dir(&d) else {
                continue;
            };
            for e in entries {
                let p = e.unwrap().path();
                if p.extension().is_none_or(|x| x != "json") || p.ends_with("golden_hashes.json") {
                    continue;
                }
                let level = Level::from_json(&std::fs::read_to_string(&p).unwrap()).unwrap();
                assert_eq!(check_editable(&level), Ok(()), "{}", p.display());
                checked += 1;
            }
        }
        assert!(checked > 0);
    }
}
