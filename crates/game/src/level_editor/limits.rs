//! The disturbances and the player's limits.

use super::*;

pub(super) fn edit_disturbances(ui: &mut egui::Ui, list: &mut Vec<Disturbance>) -> bool {
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

pub(super) fn edit_limits(
    ui: &mut egui::Ui,
    l: &mut Limits,
    texts: &mut EditTexts,
    grid: &Grid,
) -> bool {
    let Limits {
        max_charges,
        magnitudes,
        allow_positive,
        allow_negative,
        max_magnets,
        magnet_strengths,
        max_antennas,
        antenna_amplitudes,
        antenna_omegas,
        continuous,
        region,
        max_plates,
        plate_voltages,
        plate,
        supply_voltages,
        coil_rates,
        solenoid_fields,
        max_free,
        max_cuts,
        free_charges,
        free_speeds,
        free_mass,
        free_radius,
    } = l;
    let mut focus = false;
    focus |= row(ui, "Player charges (max)", |ui| {
        ui.add(egui::DragValue::new(max_charges).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Charge values |Q|", |ui| {
        let r = ui.add(egui::TextEdit::singleline(&mut texts.magnitudes).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            let parsed = parse_list(&texts.magnitudes);
            if !parsed.is_empty() {
                *magnitudes = parsed;
            }
            texts.magnitudes = list_to_text(magnitudes);
        }
        r.has_focus()
    });
    row(ui, "Charge signs", |ui| {
        ui.checkbox(allow_positive, "+");
        ui.checkbox(allow_negative, "−");
        false
    });
    focus |= row(ui, "Player magnets (max)", |ui| {
        ui.add(egui::DragValue::new(max_magnets).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Magnet values |μ|", |ui| {
        let r = ui
            .add(egui::TextEdit::singleline(&mut texts.magnet_strengths).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *magnet_strengths = parse_list(&texts.magnet_strengths);
            texts.magnet_strengths = list_to_text(magnet_strengths);
        }
        r.has_focus()
    });
    focus |= row(ui, "Player antennas (max)", |ui| {
        ui.add(egui::DragValue::new(max_antennas).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Antenna values |p₀|", |ui| {
        let r = ui.add(
            egui::TextEdit::singleline(&mut texts.antenna_amplitudes).desired_width(LIST_WIDTH),
        );
        if r.lost_focus() {
            *antenna_amplitudes = parse_list(&texts.antenna_amplitudes);
            texts.antenna_amplitudes = list_to_text(antenna_amplitudes);
        }
        r.has_focus()
    });
    focus |= row(ui, "Antenna ω values", |ui| {
        let r =
            ui.add(egui::TextEdit::singleline(&mut texts.antenna_omegas).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *antenna_omegas = parse_list(&texts.antenna_omegas);
            texts.antenna_omegas = list_to_text(antenna_omegas);
        }
        r.on_hover_text("Empty: player antennas use the level's RF generator")
            .has_focus()
    });
    focus |= row(ui, "Player plates (max)", |ui| {
        ui.add(egui::DragValue::new(max_plates).range(0..=MAX_COUNT))
            .has_focus()
    });
    focus |= row(ui, "Plate potentials V", |ui| {
        let r =
            ui.add(egui::TextEdit::singleline(&mut texts.plate_voltages).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *plate_voltages = parse_signed_list(&texts.plate_voltages);
            texts.plate_voltages = list_to_text(plate_voltages);
        }
        r.on_hover_text("Signed; 0 is a grounded plate").has_focus()
    });
    focus |= row(ui, "Plate L×T×H", |ui| {
        let mut f = positive(ui, &mut plate.length, 0.05, MAX_POSITIVE);
        f |= positive(ui, &mut plate.thickness, 0.02, MAX_POSITIVE);
        f |= positive(ui, &mut plate.height, 0.05, MAX_POSITIVE);
        f
    });
    focus |= row(ui, "Supply voltages V", |ui| {
        let r = ui
            .add(egui::TextEdit::singleline(&mut texts.supply_voltages).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *supply_voltages = parse_signed_list(&texts.supply_voltages);
            texts.supply_voltages = list_to_text(supply_voltages);
        }
        r.on_hover_text("Potentials of the power supplies of tunable electrodes (signed)")
            .has_focus()
    });
    focus |= row(ui, "Solenoid fields B", |ui| {
        let r = ui
            .add(egui::TextEdit::singleline(&mut texts.solenoid_fields).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *solenoid_fields = parse_signed_list(&texts.solenoid_fields);
            texts.solenoid_fields = list_to_text(solenoid_fields);
        }
        r.on_hover_text("Fields B along z of a tube level's solenoid supply (signed)")
            .has_focus()
    });
    focus |= row(ui, "Coil ramp rates dκ/dt", |ui| {
        let r = ui.add(egui::TextEdit::singleline(&mut texts.coil_rates).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *coil_rates = parse_signed_list(&texts.coil_rates);
            texts.coil_rates = list_to_text(coil_rates);
        }
        r.on_hover_text(
            "Ramp rates of the power supplies of tunable coils (signed; a driven coil's \
             supply sets its source's level from this list)",
        )
        .has_focus()
    });
    focus |= row(ui, "Player cuts (max)", |ui| {
        ui.add(egui::DragValue::new(max_cuts).range(0..=MAX_COUNT))
            .on_hover_text(
                "Insulating cuts through a tube level's metals (each a slit across a metal, \
                 along x or y)",
            )
            .has_focus()
    });
    focus |= row(ui, "Player free charges (max)", |ui| {
        ui.add(egui::DragValue::new(max_free).range(0..=MAX_COUNT))
            .on_hover_text(
                "Charges that move and interact (dynamic particles), launched with a velocity",
            )
            .has_focus()
    });
    focus |= row(ui, "Free charge values q", |ui| {
        let r =
            ui.add(egui::TextEdit::singleline(&mut texts.free_charges).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *free_charges = parse_signed_list(&texts.free_charges);
            texts.free_charges = list_to_text(free_charges);
        }
        r.on_hover_text("Signed charges a free charge may carry")
            .has_focus()
    });
    focus |= row(ui, "Free charge speeds", |ui| {
        let r =
            ui.add(egui::TextEdit::singleline(&mut texts.free_speeds).desired_width(LIST_WIDTH));
        if r.lost_focus() {
            *free_speeds = parse_signed_list(&texts.free_speeds)
                .into_iter()
                .filter(|v| *v >= 0.0)
                .collect();
            texts.free_speeds = list_to_text(free_speeds);
        }
        r.on_hover_text("Launch speeds, cells per time unit (0: released at rest); below c")
            .has_focus()
    });
    focus |= row(ui, "Free charge m, radius", |ui| {
        positive(ui, free_mass, 0.01, MAX_POSITIVE)
            | ui.add(
                egui::DragValue::new(free_radius)
                    .speed(0.01)
                    .range(0.0..=MAX_RADIUS),
            )
            .has_focus()
    });
    row(ui, "Continuous values", |ui| {
        ui.checkbox(continuous, "hardcore")
            .on_hover_text("Any value within the ranges of the lists above, any orientation");
        false
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
