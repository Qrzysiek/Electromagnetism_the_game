//! The player's controls: the potential sliders, the palettes of free charges,
//! cuts and plates, the power supplies, the hardcore mode's.

use super::*;

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
pub(super) fn free_palette(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
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
    explain(
        ui,
        "Free charges move and interact; pull the handle back to launch",
        &format!(
            "Free charges move and push and pull every other particle with their fields \
             (mass {}, radius {}; if two ever touch, they bounce as rigid spheres). To launch \
             one, pull the handle behind it back, like a slingshot: the arrow shows where \
             it goes and how fast (c = {}).",
            fmt_si(level.limits.free_mass),
            level.limits.free_radius,
            level.physics.c.map_or_else(|| "∞".to_string(), fmt_si)
        ),
    );
}

/// New cuts: how to draw one, and what a cut does.
pub(super) fn cut_palette(ui: &mut egui::Ui, _game: &mut Game) {
    explain(
        ui,
        "A cut: drag across a metal to draw an insulating slit",
        "Press on a point of a metal and release on another: a straight slit between them, any \
         length and angle within the metal, which its carriers cannot cross. A current loop \
         must go round it, so eddy currents weaken (laminations); a slit from edge to edge \
         splits the metal in two. Drag a cut by its start to move it, R turns it by 15°, the \
         right button removes it.",
    );
}

pub(super) fn plate_palette(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
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
pub(super) fn supplies(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
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
    let focus = game.focused_supply();
    for (k, (centre, target)) in targets.into_iter().enumerate() {
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
            level::SupplyTarget::Solenoid => {
                let own = level.tube.as_ref().map_or(0.0, |t| t.b_z);
                (
                    "Solenoid B:".to_string(),
                    format!("off, B = {} (out of the board positive)", fmt_si(own)),
                )
            }
        };
        let list = level.supply_list(target).to_vec();
        let current = game.editor.supply(centre);
        let row = ui.horizontal_wrapped(|ui| {
            // Its number, as on the map beside the electrodes it drives.
            let badge = if focus == Some(centre) {
                SUPPLY_FOCUS
            } else {
                SUPPLY_FRAME
            };
            ui.label(
                egui::RichText::new(format!(" {} ", k + 1))
                    .strong()
                    .color(egui::Color32::BLACK)
                    .background_color(badge),
            )
            .on_hover_text("Its number on the map, beside every electrode it drives");
            // The electrode's function on hover: in the row it would overflow the panel.
            let named = ui.label(name);
            if let level::SupplyTarget::Electrode(i) = target
                && !level.electrodes[i].label.is_empty()
            {
                named.on_hover_text(level.electrode_name(i));
            }
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
            if choice != current {
                game.focus_supply(centre, 1.5);
            }
            game.editor.set_supply(centre, choice);
        });
        if ui.rect_contains_pointer(row.response.rect) {
            game.focus_supply(centre, 0.15);
        }
    }
}

/// Hardcore sliders: the element under the cursor if there is one (edited live), else
/// the values of new elements.
pub(super) fn hardcore_controls(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
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
        // A cut has no value (its orientation is the palette's).
        ElementKind::Cut => return,
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
