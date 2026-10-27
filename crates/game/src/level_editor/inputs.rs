//! The editor's input widgets: lists, numbers with SI prefixes, ranges, nodes,
//! rows, kinds, biases and labels.

use super::*;

pub fn list_to_text(m: &[f64]) -> String {
    m.iter()
        .map(|v| format!("{v:e}"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// A list of signed values (potentials; 0 allowed), without duplicates.
pub(super) fn parse_signed_list(text: &str) -> Vec<f64> {
    let mut out: Vec<f64> = text
        .split(',')
        .filter_map(|t| parse_si(t.trim()))
        .filter(|v| v.is_finite())
        .collect();
    out.sort_by(f64::total_cmp);
    out.dedup();
    out
}

pub(super) fn parse_list(text: &str) -> Vec<f64> {
    text.split(',')
        .filter_map(|t| parse_si(t.trim()))
        .filter(|v| *v > 0.0 && v.is_finite())
        .collect()
}

/// Any finite value, with SI suffixes accepted when typed.
pub(super) fn si(ui: &mut egui::Ui, v: &mut f64, speed: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .custom_formatter(|v, _| fmt_si(v))
            .custom_parser(parse_si),
    )
    .has_focus()
}

/// A small positive number shown in scientific notation (tolerances).
pub(super) fn sci(ui: &mut egui::Ui, v: &mut f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(0.0)
            .range(TOLERANCE)
            .custom_formatter(|v, _| format!("{v:.0e}"))
            .custom_parser(|t| t.trim().parse::<f64>().ok()),
    )
    .has_focus()
}

pub(super) fn positive(ui: &mut egui::Ui, v: &mut f64, speed: f64, hi: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .range(MIN_POSITIVE..=hi)
            // Large and tiny values in scientific notation: "50000000.00" made a row wider
            // than the panel (scaled beam particles).
            .custom_formatter(|v, range| {
                if v.abs() >= 1e4 || (v != 0.0 && v.abs() < 1e-3) {
                    format!("{v:.2e}")
                } else {
                    egui::emath::format_with_decimals_in_range(v, range)
                }
            })
            .custom_parser(|t| t.trim().parse::<f64>().ok()),
    )
    .has_focus()
}

/// A value not below 0 (times, frequencies, series inductance and capacitance).
pub(super) fn non_negative(ui: &mut egui::Ui, v: &mut f64, speed: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .range(0.0..=MAX_POSITIVE)
            .custom_formatter(|v, _| fmt_si(v))
            .custom_parser(parse_si),
    )
    .has_focus()
}

/// A grid node (z shown only for 3D grids).
pub(super) fn node(ui: &mut egui::Ui, n: &mut Node, grid: &Grid) -> bool {
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

pub(super) fn row(ui: &mut egui::Ui, label: &str, add: impl FnOnce(&mut egui::Ui) -> bool) -> bool {
    ui.label(label);
    let f = ui.horizontal(add).inner;
    ui.end_row();
    f
}

pub(super) fn kind_combo(ui: &mut egui::Ui, id: usize, kind: &mut ElementKind) {
    let text = |k: ElementKind| match k {
        ElementKind::Charge => "charge",
        ElementKind::Magnet => "magnet μ",
        ElementKind::Antenna => "antenna p₀",
        ElementKind::Plate => "plate V",
        ElementKind::Supply => "supply V",
        ElementKind::Free => "free q",
        ElementKind::Cut => "cut",
    };
    egui::ComboBox::from_id_salt(("element_kind", id))
        .selected_text(text(*kind))
        .width(80.0)
        .show_ui(ui, |ui| {
            for k in [
                ElementKind::Charge,
                ElementKind::Magnet,
                ElementKind::Antenna,
            ] {
                ui.selectable_value(kind, k, text(k));
            }
        });
}

/// Bias selector shared by metal spheres and electrodes.
pub(super) fn bias_editor(ui: &mut egui::Ui, id: (&str, usize), bias: &mut ConductorBias) -> bool {
    let kind = match bias {
        ConductorBias::Grounded => 0,
        ConductorBias::Charge(_) => 1,
        ConductorBias::Potential(_) => 2,
    };
    let mut k = kind;
    egui::ComboBox::from_id_salt(id)
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
        ConductorBias::Grounded => false,
        ConductorBias::Charge(v) | ConductorBias::Potential(v) => si(ui, v, 1e4),
    }
}

/// An object's function in its device (`Electrode::label`, `MetalSpec::label`): text only.
pub(super) fn label_editor(ui: &mut egui::Ui, label: &mut String, examples: &str) -> bool {
    ui.horizontal(|ui| {
        ui.label("function");
        ui.add(
            egui::TextEdit::singleline(label)
                .hint_text(examples)
                .desired_width(180.0),
        )
        .on_hover_text("Its function in the device, shown with its number wherever it is named")
        .has_focus()
    })
    .inner
}
