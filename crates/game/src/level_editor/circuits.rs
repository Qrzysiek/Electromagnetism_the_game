//! The circuits' sources and drives, and whether the player may tune them.

use super::*;

/// A source's kind (DC, sine, pulse) and its parameters.
pub(super) fn edit_source(ui: &mut egui::Ui, id: (&str, usize), source: &mut Source) -> bool {
    let mut focus = false;
    let kind = match source {
        Source::Dc { .. } => 0,
        Source::Sine { .. } => 1,
        Source::Pulse { .. } => 2,
    };
    let mut k = kind;
    ui.horizontal(|ui| {
        ui.label("source");
        egui::ComboBox::from_id_salt(id)
            .selected_text(["DC", "sine", "pulse"][k])
            .width(60.0)
            .show_ui(ui, |ui| {
                for (j, t) in ["DC", "sine", "pulse"].iter().enumerate() {
                    ui.selectable_value(&mut k, j, *t);
                }
            });
    });
    if k != kind {
        *source = match k {
            1 => Source::Sine {
                offset: 0.0,
                amplitude: 1.0,
                omega: 1.0,
                phase: 0.0,
            },
            2 => Source::Pulse {
                low: 0.0,
                high: 1.0,
                delay: 1.0,
                rise: 0.5,
                width: 5.0,
                fall: 0.5,
                period: 0.0,
            },
            _ => Source::Dc { value: 1.0 },
        };
    }
    match source {
        Source::Dc { value } => {
            ui.horizontal(|ui| {
                ui.label("V");
                focus |= si(ui, value, 1e3);
            });
        }
        Source::Sine {
            offset,
            amplitude,
            omega,
            phase,
        } => {
            ui.horizontal(|ui| {
                ui.label("V₀ + V sin(ωt + φ): V₀");
                focus |= si(ui, offset, 1e3);
                ui.label("V");
                focus |= si(ui, amplitude, 1e3);
            });
            ui.horizontal(|ui| {
                ui.label("ω");
                focus |= non_negative(ui, omega, 0.01);
                ui.label("φ");
                focus |= si(ui, phase, 0.01);
            });
        }
        Source::Pulse {
            low,
            high,
            delay,
            rise,
            width,
            fall,
            period,
        } => {
            ui.horizontal(|ui| {
                ui.label("low");
                focus |= si(ui, low, 1e3);
                ui.label("high");
                focus |= si(ui, high, 1e3);
            });
            ui.horizontal(|ui| {
                ui.label("delay");
                focus |= non_negative(ui, delay, 0.1);
                ui.label("rise");
                focus |= non_negative(ui, rise, 0.1);
                ui.label("top");
                focus |= non_negative(ui, width, 0.1);
            });
            ui.horizontal(|ui| {
                ui.label("fall");
                focus |= non_negative(ui, fall, 0.1);
                ui.label("period").on_hover_text("0: a single pulse");
                focus |= non_negative(ui, period, 0.1);
            });
        }
    }
    focus
}

/// The circuit driving an electrode or a circular coil (PHYSICS.md §2.10): a source
/// through a resistance, with a series inductance or a switch (electrodes) or a series
/// capacitance (coils).
pub(super) fn drive_editor(
    ui: &mut egui::Ui,
    id: (&str, usize),
    drive: &mut Option<Drive>,
    electrode: bool,
) -> bool {
    let mut focus = false;
    let mut on = drive.is_some();
    if ui
        .checkbox(&mut on, "driven by a circuit")
        .on_hover_text(
            "A source through a resistance drives it (one way: the particles do not act \
             back); it starts from its static state",
        )
        .changed()
    {
        *drive = on.then_some(Drive {
            source: Source::Dc { value: 1.0 },
            resistance: 1.0,
            inductance: 0.0,
            capacitance: 0.0,
            switch: None,
        });
    }
    let Some(Drive {
        source,
        resistance,
        inductance,
        capacitance,
        switch,
    }) = drive
    else {
        return focus;
    };
    focus |= edit_source(ui, id, source);
    ui.horizontal(|ui| {
        ui.label("R");
        focus |= positive(ui, resistance, 0.01, MAX_POSITIVE);
        if electrode {
            ui.label("series L")
                .on_hover_text("0: none; with the plate's capacitance an LC");
            focus |= non_negative(ui, inductance, 0.01);
        } else {
            ui.label("series C")
                .on_hover_text("0: none; with the coil's inductance an LC");
            focus |= non_negative(ui, capacitance, 0.01);
        }
    });
    if electrode {
        let mut has = switch.is_some();
        ui.horizontal(|ui| {
            if ui
                .checkbox(&mut has, "switch")
                .on_hover_text(
                    "The resistance is a switch toggling once; open, the plate floats with \
                     its charge (not with a series L)",
                )
                .changed()
            {
                *switch = has.then_some(Switch {
                    closed: false,
                    at: 1.0,
                });
            }
            if let Some(Switch { closed, at }) = switch {
                ui.checkbox(closed, "closed first");
                ui.label("toggles at");
                focus |= non_negative(ui, at, 0.1);
            }
        });
    }
    focus
}

/// A source within the editors' ranges.
pub(super) fn source_editable(source: &Source) -> bool {
    let non_negative = |v: &f64| (0.0..=MAX_POSITIVE).contains(v);
    match source {
        Source::Dc { value } => value.is_finite(),
        Source::Sine {
            offset,
            amplitude,
            omega,
            phase,
        } => {
            offset.is_finite() && amplitude.is_finite() && non_negative(omega) && phase.is_finite()
        }
        Source::Pulse {
            low,
            high,
            delay,
            rise,
            width,
            fall,
            period,
        } => {
            low.is_finite()
                && high.is_finite()
                && [delay, rise, width, fall, period]
                    .iter()
                    .all(|v| non_negative(v))
        }
    }
}

/// A drive within the editors' ranges.
pub(super) fn drive_editable(d: &Drive) -> bool {
    let Drive {
        source,
        resistance,
        inductance,
        capacitance,
        switch,
    } = d;
    let non_negative = |v: &f64| (0.0..=MAX_POSITIVE).contains(v);
    let source_ok = source_editable(source);
    source_ok
        && (MIN_POSITIVE..=MAX_POSITIVE).contains(resistance)
        && non_negative(inductance)
        && non_negative(capacitance)
        && switch.is_none_or(|Switch { closed: _, at }| non_negative(&at))
}
