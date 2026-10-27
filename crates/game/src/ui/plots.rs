//! The panel's plots: the electrodes' potentials, the circuits, their ticks and
//! the overlaid quantities.

use super::*;

/// How the panel groups the electrodes' potential plots (`potential_plots`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PlotLayout {
    /// Each metal's contacts on one plot, a potential goal's electrode with its reference.
    #[default]
    Grouped,
    /// Every electrode on one plot.
    Together,
    /// Each electrode on its own plot (its own scale).
    Separate,
}

/// The potentials of the electrodes in circuits (`(electrode, samples)`), grouped on plots
/// (`Game::plot_layout`): the lines of a plot share its scale, so that they compare at a
/// glance. Under each plot a chip per electrode in its line's colour: a click shows or hides
/// the line (the scale follows the lines shown), hovering outlines the electrode on the
/// board in that colour.
pub(super) fn potential_plots(
    ui: &mut egui::Ui,
    game: &mut Game,
    spec: &level::tube::TubeSpec,
    potentials: &[(usize, Vec<(f64, f64)>)],
) {
    let marks: Vec<f64> = spec
        .potential_goals
        .iter()
        .flat_map(|g| [g.start, g.end])
        .collect();
    let colour = |k: usize| PLOT_COLOURS[k % PLOT_COLOURS.len()];
    if potentials.len() > 1 {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Potentials:").small());
            for (layout, name, hint) in [
                (
                    PlotLayout::Grouped,
                    "by metal",
                    "Each metal's contacts on one plot (and a voltage goal's two electrodes), \
                     on one scale",
                ),
                (
                    PlotLayout::Together,
                    "together",
                    "Every electrode on one plot, one scale",
                ),
                (
                    PlotLayout::Separate,
                    "separate",
                    "Each electrode on its own plot and scale",
                ),
            ] {
                ui.selectable_value(&mut game.plot_layout, layout, name)
                    .on_hover_text(hint);
            }
        });
    }
    // The groups: indices into `potentials`.
    let n = potentials.len();
    let mut group: Vec<usize> = (0..n).collect();
    let find = |g: &[usize], mut i: usize| {
        while g[i] != i {
            i = g[i];
        }
        i
    };
    let together = |a: usize, b: usize, g: &mut Vec<usize>| {
        let (ra, rb) = (find(g, a), find(g, b));
        g[ra.max(rb)] = ra.min(rb);
    };
    let at = |e: usize| potentials.iter().position(|(j, _)| *j == e);
    match game.plot_layout {
        PlotLayout::Separate => {}
        PlotLayout::Together => {
            for i in 1..n {
                together(0, i, &mut group);
            }
        }
        PlotLayout::Grouped => {
            for i in 0..n {
                for j in 0..i {
                    let (a, b) = (potentials[i].0, potentials[j].0);
                    if spec.metal_of(a).is_some() && spec.metal_of(a) == spec.metal_of(b) {
                        together(i, j, &mut group);
                    }
                }
            }
            for g in &spec.potential_goals {
                if let (Some(i), Some(j)) = (at(g.electrode), g.reference.and_then(at)) {
                    together(i, j, &mut group);
                }
            }
        }
    }
    for root in 0..n {
        let members: Vec<usize> = (0..n).filter(|&i| find(&group, i) == root).collect();
        if members.is_empty() {
            continue;
        }
        let electrodes: Vec<String> = members
            .iter()
            .map(|&i| format!("{}", potentials[i].0 + 1))
            .collect();
        let metal = members
            .iter()
            .map(|&i| spec.metal_of(potentials[i].0))
            .reduce(|a, b| if a == b { a } else { None })
            .flatten();
        let title = match (members.len(), metal) {
            (1, _) => format!(
                "{}: potential",
                game.editor.level.electrode_name(potentials[members[0]].0)
            ),
            (_, Some(m)) => format!(
                "{}: potentials of electrodes {}",
                game.editor.level.metal_name(m),
                electrodes.join(", ")
            ),
            _ => format!("Potentials of electrodes {}", electrodes.join(", ")),
        };
        ui.label(egui::RichText::new(title).small());
        let lines: Vec<Line> = members
            .iter()
            .filter(|&&i| !game.plot_hidden.contains(&potentials[i].0))
            .map(|&i| (colour(i), potentials[i].1.as_slice()))
            .collect();
        multi_plot(ui, &lines, &marks, "V");
        if members.len() > 1 {
            ui.horizontal_wrapped(|ui| {
                for &i in &members {
                    let e = potentials[i].0;
                    let shown = !game.plot_hidden.contains(&e);
                    let [r, g, b] = colour(i);
                    let text = egui::RichText::new(format!("— electrode {}", e + 1))
                        .small()
                        .color(if shown {
                            egui::Color32::from_rgb(r, g, b)
                        } else {
                            egui::Color32::from_gray(110)
                        });
                    let response = ui
                        .add(egui::Label::new(text).sense(egui::Sense::click()))
                        .on_hover_text(
                            "Click to show or hide its line (the scale follows the lines \
                             shown); it is outlined on the board in this colour.",
                        );
                    if response.clicked() {
                        if shown {
                            game.plot_hidden.push(e);
                        } else {
                            game.plot_hidden.retain(|&j| j != e);
                        }
                    }
                    if response.hovered() {
                        game.focus_electrode(e, [r, g, b], 0.15);
                    }
                }
            });
        }
    }
}

/// A small plot of a circuit quantity over lab time (PHYSICS.md §2.10), its range
/// including 0, with lines at the times `marks` (the launches), `y` the quantity's symbol.
pub(super) fn circuit_plot(ui: &mut egui::Ui, samples: &[(f64, f64)], marks: &[f64], y: &str) {
    multi_plot(ui, &[(PLOT_COLOURS[0], samples)], marks, y);
}

/// Round values over `lo..=hi` for an axis: steps of 1, 2 or 5 times a power of ten, at most
/// about `n` of them.
pub(super) fn ticks(lo: f64, hi: f64, n: u32) -> Vec<f64> {
    if !(lo.is_finite() && hi.is_finite()) || hi <= lo {
        return vec![lo];
    }
    let raw = (hi - lo) / f64::from(n.max(1));
    // The power of ten at or below `raw` (without logarithms: the platform's differ).
    let mut mag = 1.0;
    while mag * 10.0 <= raw {
        mag *= 10.0;
    }
    while mag > raw {
        mag /= 10.0;
    }
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|f| f * mag)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * mag);
    let mut out = Vec::new();
    let mut v = (lo / step).ceil() * step;
    while v <= hi + 1e-9 * step {
        out.push(if v.abs() < 1e-9 * step { 0.0 } else { v });
        v += step;
    }
    out
}

/// A small plot of several series over lab time on one scale (its range including 0),
/// each in its colour, with lines at the times `marks`: axes with round ticks and their
/// numbers, t along the bottom, `y` (the quantity's symbol, if any) at the top of the left
/// axis.
#[allow(clippy::cast_possible_truncation)]
pub(super) fn multi_plot(ui: &mut egui::Ui, lines: &[Line], marks: &[f64], y: &str) {
    let t0 = lines
        .iter()
        .filter_map(|(_, s)| s.first().map(|p| p.0))
        .fold(f64::INFINITY, f64::min);
    let t1 = lines
        .iter()
        .filter_map(|(_, s)| s.last().map(|p| p.0))
        .fold(f64::NEG_INFINITY, f64::max);
    let width = ui.available_width().min(280.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 84.0), egui::Sense::hover());
    if !(t0.is_finite() && t1.is_finite()) {
        return;
    }
    let painter = ui.painter_at(rect);
    let (lo, hi) = lines
        .iter()
        .flat_map(|(_, s)| s.iter())
        .fold((0.0_f64, 0.0_f64), |(a, b), s| (a.min(s.1), b.max(s.1)));
    let (lo, hi) = if hi > lo {
        (lo, hi)
    } else {
        (lo - 1.0, hi + 1.0)
    };
    let (t0, t1) = if t1 > t0 { (t0, t1) } else { (t0, t0 + 1.0) };
    // The plotting area: room for the numbers on the left and below, the symbol above.
    let area = egui::Rect::from_min_max(
        egui::pos2(rect.left() + 34.0, rect.top() + 9.0),
        egui::pos2(rect.right() - 4.0, rect.bottom() - 13.0),
    );
    let at = |t: f64, v: f64| {
        egui::pos2(
            area.left() + ((t - t0) / (t1 - t0)) as f32 * area.width(),
            area.bottom() - ((v - lo) / (hi - lo)) as f32 * area.height(),
        )
    };
    let weak = ui.visuals().weak_text_color();
    let font = egui::FontId::proportional(9.5);
    let grid = egui::Stroke::new(1.0, egui::Color32::from_gray(48));
    for v in ticks(lo, hi, 4) {
        let p = at(t0, v);
        painter.line_segment([p, egui::pos2(area.right(), p.y)], grid);
        painter.text(
            egui::pos2(area.left() - 3.0, p.y),
            egui::Align2::RIGHT_CENTER,
            fmt_si(v),
            font.clone(),
            weak,
        );
    }
    for t in ticks(t0, t1, 5) {
        let p = at(t, lo);
        painter.line_segment([egui::pos2(p.x, area.top()), p], grid);
        painter.text(
            egui::pos2(p.x, area.bottom() + 1.0),
            egui::Align2::CENTER_TOP,
            fmt_si(t),
            font.clone(),
            weak,
        );
    }
    let axis = egui::Stroke::new(1.0, egui::Color32::from_gray(120));
    painter.line_segment([area.left_bottom(), area.left_top()], axis);
    painter.line_segment([area.left_bottom(), area.right_bottom()], axis);
    if lo < 0.0 && hi > 0.0 {
        painter.line_segment(
            [at(t0, 0.0), at(t1, 0.0)],
            egui::Stroke::new(1.0, egui::Color32::from_gray(95)),
        );
    }
    painter.text(
        egui::pos2(area.left() - 3.0, rect.top()),
        egui::Align2::RIGHT_TOP,
        y,
        font.clone(),
        weak,
    );
    painter.text(
        egui::pos2(area.right(), area.bottom() + 1.0),
        egui::Align2::RIGHT_TOP,
        "t",
        font,
        weak,
    );
    for &m in marks {
        if (t0..=t1).contains(&m) {
            painter.line_segment(
                [at(m, lo), at(m, hi)],
                egui::Stroke::new(1.0, egui::Color32::from_rgb(230, 200, 80)),
            );
        }
    }
    for ([r, g, b], samples) in lines {
        let points: Vec<egui::Pos2> = samples.iter().map(|&(t, v)| at(t, v)).collect();
        painter.add(egui::Shape::line(
            points,
            egui::Stroke::new(1.5, egui::Color32::from_rgb(*r, *g, *b)),
        ));
    }
}
