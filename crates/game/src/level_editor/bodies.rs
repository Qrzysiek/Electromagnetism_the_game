//! The level's bodies: elements, coils, metal spheres, clouds, electrodes,
//! dielectrics and ferrites.

use super::*;

pub(super) fn edit_elements(ui: &mut egui::Ui, elements: &mut Vec<Element>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, e) in elements.iter_mut().enumerate() {
        let Element {
            node: n,
            kind,
            value,
            angle_deg,
            omega,
            // Free charges are player elements (the level's are `free_particles`).
            speed: _,
        } = e;
        ui.horizontal(|ui| {
            kind_combo(ui, i, kind);
            focus |= si(ui, value, 1e4);
            if *kind == ElementKind::Antenna {
                focus |= ui
                    .add(egui::DragValue::new(angle_deg).speed(1.0).suffix("°"))
                    .has_focus();
                let mut own = omega.is_some();
                ui.checkbox(&mut own, "ω")
                    .on_hover_text("Own frequency (otherwise the level's RF generator)");
                match (own, *omega) {
                    (true, None) => *omega = Some(1.0),
                    (false, Some(_)) => *omega = None,
                    _ => {}
                }
                if let Some(w) = omega {
                    focus |= ui
                        .add(
                            egui::DragValue::new(w)
                                .speed(0.01)
                                .range(0.0..=MAX_POSITIVE),
                        )
                        .has_focus();
                }
            }
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

pub(super) fn edit_coils(ui: &mut egui::Ui, coils: &mut Vec<Coil>, grid: &Grid) -> bool {
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
                    rate,
                    drive,
                    tunable,
                } => {
                    ui.horizontal(|ui| {
                        ui.label("circle, centre");
                        focus |= node(ui, center, grid);
                        ui.checkbox(tunable, "tunable").on_hover_text(
                            "The player sets its ramp rate (a driven coil: its source's level) \
                             with a power supply (Limits: coil ramp rates)",
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("radius");
                        focus |= positive(ui, radius, 0.05, MAX_POSITIVE);
                        ui.label("κ");
                        focus |= si(ui, kappa, 100.0);
                    });
                    ui.horizontal(|ui| {
                        ui.label("ramp dκ/dt").on_hover_text(
                            "The current changes linearly in time, κ + rate·t, inducing an \
                             electric field −∂A/∂t (Faraday's law)",
                        );
                        focus |= si(ui, rate, 1.0);
                    });
                    focus |= drive_editor(ui, ("coil drive", i), drive, false);
                }
                Coil::Polygon {
                    vertices,
                    kappa,
                    rate,
                    drive,
                } => {
                    ui.horizontal(|ui| {
                        ui.label("polygon, κ");
                        focus |= si(ui, kappa, 100.0);
                        ui.label("ramp dκ/dt");
                        focus |= si(ui, rate, 1.0);
                    });
                    focus |= drive_editor(ui, ("coil drive", i), drive, false);
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

pub(super) fn edit_conductors(ui: &mut egui::Ui, list: &mut Vec<Conductor>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, c) in list.iter_mut().enumerate() {
        let Conductor {
            center,
            radius,
            bias,
        } = c;
        ui.push_id(("conductor", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Sphere {}", i + 1));
                focus |= node(ui, center, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("radius");
                focus |= positive(ui, radius, 0.05, MAX_POSITIVE);
                let kind = match bias {
                    ConductorBias::Grounded => 0,
                    ConductorBias::Charge(_) => 1,
                    ConductorBias::Potential(_) => 2,
                };
                let mut k = kind;
                egui::ComboBox::from_id_salt(("bias", i))
                    .selected_text(["grounded", "charge Q", "potential V"][k])
                    .width(90.0)
                    .show_ui(ui, |ui| {
                        for (j, t) in ["grounded", "charge Q", "potential V"].iter().enumerate() {
                            ui.selectable_value(&mut k, j, *t);
                        }
                    });
                if k != kind {
                    *bias = match k {
                        1 => ConductorBias::Charge(0.0),
                        2 => ConductorBias::Potential(0.0),
                        _ => ConductorBias::Grounded,
                    };
                }
                match bias {
                    ConductorBias::Grounded => {}
                    ConductorBias::Charge(v) | ConductorBias::Potential(v) => {
                        focus |= si(ui, v, 1e4);
                    }
                }
            });
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ metal sphere").clicked() {
        let m = grid.max_node();
        list.push(Conductor {
            center: [m[0] / 2, m[1] / 2, 0],
            radius: 2.0,
            bias: ConductorBias::Grounded,
        });
    }
    focus
}

pub(super) fn edit_clouds(ui: &mut egui::Ui, list: &mut Vec<Cloud>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, c) in list.iter_mut().enumerate() {
        let Cloud {
            center,
            radius,
            charge,
        } = c;
        ui.push_id(("cloud", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Cloud {}", i + 1));
                focus |= node(ui, center, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("radius");
                focus |= positive(ui, radius, 0.05, MAX_POSITIVE);
                ui.label("charge Q");
                focus |= si(ui, charge, 1e4);
            });
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize
        && ui
            .small_button("+ charge cloud")
            .on_hover_text(
                "A sphere of uniform charge that particles fly through (Thomson's atom): \
                 inside, an opposite charge is bound harmonically.",
            )
            .clicked()
    {
        let m = grid.max_node();
        list.push(Cloud {
            center: [m[0] / 2, m[1] / 2, 0],
            radius: 3.0,
            charge: 1.0,
        });
    }
    focus
}

pub(super) fn edit_electrodes(ui: &mut egui::Ui, list: &mut Vec<Electrode>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, e) in list.iter_mut().enumerate() {
        let Electrode {
            center,
            length,
            thickness,
            height,
            angle_deg,
            shape,
            bias,
            tunable,
            drive,
            label,
        } = e;
        ui.push_id(("electrode", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Electrode {}", i + 1));
                focus |= node(ui, center, grid);
                ui.checkbox(tunable, "tunable").on_hover_text(
                    "The player sets its potential with a power supply (Limits: supply \
                     voltages); otherwise it keeps the bias below",
                );
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("shape")
                    .selected_text(match shape {
                        ElectrodeShape::Box => "box",
                        ElectrodeShape::Rod => "rod",
                        ElectrodeShape::Ring => "ring",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(shape, ElectrodeShape::Box, "box");
                        ui.selectable_value(shape, ElectrodeShape::Rod, "rod")
                            .on_hover_text("Tube levels: a circle, diameter L");
                        ui.selectable_value(shape, ElectrodeShape::Ring, "ring")
                            .on_hover_text(
                                "Tube levels: an annulus, mid-diameter L, wall T, open inside",
                            );
                    });
                ui.label("L×T×H");
                focus |= positive(ui, length, 0.05, MAX_POSITIVE);
                focus |= positive(ui, thickness, 0.02, MAX_POSITIVE);
                focus |= positive(ui, height, 0.05, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                focus |= ui
                    .add(egui::DragValue::new(angle_deg).speed(1.0).suffix("°"))
                    .has_focus();
                focus |= bias_editor(ui, ("ebias", i), bias);
            });
            focus |= label_editor(ui, label, "cathode, grid, source, gate, deflection plate…");
            focus |= drive_editor(ui, ("electrode drive", i), drive, true);
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ electrode").clicked() {
        let m = grid.max_node();
        list.push(Electrode {
            center: [m[0] / 2, m[1] / 2, 0],
            length: 6.0,
            thickness: 0.4,
            height: 4.0,
            angle_deg: 0.0,
            shape: ElectrodeShape::Box,
            bias: ConductorBias::Grounded,
            tunable: false,
            drive: None,
            label: String::new(),
        });
    }
    focus
}

pub(super) fn edit_dielectrics(ui: &mut egui::Ui, list: &mut Vec<Dielectric>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, d) in list.iter_mut().enumerate() {
        let Dielectric {
            center,
            length,
            thickness,
            height,
            angle_deg,
            permittivity,
        } = d;
        ui.push_id(("dielectric", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Dielectric {}", i + 1));
                focus |= node(ui, center, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("L×T×H");
                focus |= positive(ui, length, 0.05, MAX_POSITIVE);
                focus |= positive(ui, thickness, 0.02, MAX_POSITIVE);
                focus |= positive(ui, height, 0.05, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                focus |= ui
                    .add(egui::DragValue::new(angle_deg).speed(1.0).suffix("°"))
                    .has_focus();
                ui.label("ε");
                focus |= ui
                    .add(
                        egui::DragValue::new(permittivity)
                            .speed(0.1)
                            .range(1.0..=MAX_PERMITTIVITY),
                    )
                    .has_focus();
            });
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ dielectric").clicked() {
        let m = grid.max_node();
        list.push(Dielectric {
            center: [m[0] / 2, m[1] / 2, 0],
            length: 4.0,
            thickness: 2.0,
            height: 4.0,
            angle_deg: 0.0,
            permittivity: 4.0,
        });
    }
    focus
}

pub(super) fn edit_dielectric_spheres(
    ui: &mut egui::Ui,
    list: &mut Vec<DielectricSphere>,
    grid: &Grid,
) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, d) in list.iter_mut().enumerate() {
        let DielectricSphere {
            center,
            radius,
            permittivity,
        } = d;
        ui.push_id(("dielectric sphere", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Sphere {}", i + 1));
                focus |= node(ui, center, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("radius");
                focus |= positive(ui, radius, 0.05, MAX_POSITIVE);
                ui.label("ε");
                focus |= ui
                    .add(
                        egui::DragValue::new(permittivity)
                            .speed(0.1)
                            .range(1.0..=MAX_PERMITTIVITY),
                    )
                    .has_focus();
            });
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ dielectric sphere").clicked() {
        let m = grid.max_node();
        list.push(DielectricSphere {
            center: [m[0] / 2, m[1] / 2, 0],
            radius: 2.0,
            permittivity: 4.0,
        });
    }
    focus
}

pub(super) fn edit_ferrites(ui: &mut egui::Ui, list: &mut Vec<FerritePlate>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, f) in list.iter_mut().enumerate() {
        let FerritePlate {
            center,
            length,
            thickness,
            height,
            angle_deg,
            permeability,
        } = f;
        ui.push_id(("ferrite", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Ferrite {}", i + 1));
                focus |= node(ui, center, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            ui.horizontal(|ui| {
                ui.label("L×T×H");
                focus |= positive(ui, length, 0.05, MAX_POSITIVE);
                focus |= positive(ui, thickness, 0.02, MAX_POSITIVE);
                focus |= positive(ui, height, 0.05, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                focus |= ui
                    .add(egui::DragValue::new(angle_deg).speed(1.0).suffix("°"))
                    .has_focus();
                ui.label("μ");
                focus |= ui
                    .add(
                        egui::DragValue::new(permeability)
                            .speed(1.0)
                            .range(1.0..=MAX_PERMITTIVITY),
                    )
                    .has_focus();
            });
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ ferrite").clicked() {
        let m = grid.max_node();
        list.push(FerritePlate {
            center: [m[0] / 2, m[1] / 2, 0],
            length: 6.0,
            thickness: 1.0,
            height: 4.0,
            angle_deg: 0.0,
            permeability: 1000.0,
        });
    }
    focus
}
