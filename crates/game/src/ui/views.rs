//! The views' controls and keys: the field maps, the flow, the beams' energies,
//! the colour bars, the radiation goals.

use super::*;

/// Controls of the time-dependent field views: what to show (the full field, or a part of
/// it measured against the full field), the colour quantity and scale, a colour bar with
/// ticks, and the size of the part shown.
pub(super) fn field_view_controls(
    ui: &mut egui::Ui,
    game: &mut Game,
    level: &level::Level,
    radiation: &crate::radiation::RadiationView,
) {
    use crate::radiation::FieldQuantity;
    // What to show: the full field (baseline) or a part of it.
    let quasi_static = level.has_beams()
        && level.interacts()
        && !level.physics.beam_retarded
        && level.physics.c.is_some();
    if game.map == Some(MapMode::ParticleField) {
        let before = (game.radiation_only, game.neglected_only);
        ui.horizontal_wrapped(|ui| {
            ui.label("Show:")
                .on_hover_text("The full field, or a part of it on the full field's colour scale");
            if ui
                .selectable_label(!game.radiation_only && !game.neglected_only, "full field")
                .on_hover_text("The whole field of the moving charges: the baseline")
                .clicked()
            {
                (game.radiation_only, game.neglected_only) = (false, false);
            }
            // At c = ∞ nothing radiates.
            if level.physics.c.is_some()
                && ui
                    .selectable_label(game.radiation_only, "radiation part")
                    .on_hover_text(
                        "Only the acceleration term of the field: what travels away at c and \
                     falls off as 1/R (the rest is the velocity field that moves with the \
                     charge)",
                    )
                    .clicked()
            {
                (game.radiation_only, game.neglected_only) = (true, false);
            }
            if quasi_static
                && ui
                    .selectable_label(game.neglected_only, "left out by the model")
                    .on_hover_text(
                        "The full retarded field minus the fields the quasi-static beam \
                         interaction uses (present states continued with constant \
                         acceleration): what the dynamics leaves out (PHYSICS.md §3.3)",
                    )
                    .clicked()
            {
                (game.radiation_only, game.neglected_only) = (false, true);
            }
        });
        // A part is best seen on a linear scale: switch when the choice changes.
        let now = (game.radiation_only, game.neglected_only);
        if now != before {
            game.field_linear = now != (false, false);
        }
    }
    // The energy flow needs a finite speed of light (the magnetic energy of a finite B is
    // infinite at c = ∞ in these units).
    if level.physics.c.is_none() && game.field_quantity == FieldQuantity::S {
        game.field_quantity = FieldQuantity::E;
    }
    let flow = game.field_quantity == FieldQuantity::S;
    ui.horizontal_wrapped(|ui| {
        ui.label("Colour:")
            .on_hover_text("Which quantity colours the map");
        ui.selectable_value(&mut game.field_quantity, FieldQuantity::Bz, "B_z")
            .on_hover_text("B, perpendicular to the plane (all of B in the plane); signed");
        ui.selectable_value(&mut game.field_quantity, FieldQuantity::E, "|E|")
            .on_hover_text("The magnitude of E");
        if level.physics.c.is_some() {
            ui.selectable_value(&mut game.field_quantity, FieldQuantity::S, "S")
                .on_hover_text(
                    "The flow of field energy, the Poynting vector S = (c²/4π) E × B \
                     (Jackson §6.7): its magnitude colours the map, arrows and tracers show \
                     where it goes",
                );
        }
        ui.checkbox(
            &mut game.show_field_arrows,
            if flow { "S arrows" } else { "E arrows" },
        )
        .on_hover_text(if flow {
            "Arrows of the energy flow S: direction exact, length on the colour scale"
        } else {
            "Arrows of E: direction exact, length on the colour scale"
        });
    });
    if flow {
        flow_controls(ui, game, radiation);
    }
    ui.horizontal(|ui| {
        ui.label("Scale:").on_hover_text("How field values map to colour");
        ui.selectable_value(&mut game.field_linear, false, "log")
            .on_hover_text("Logarithmic: shows fields over many decades at once, but makes small ones look big");
        ui.selectable_value(&mut game.field_linear, true, "linear")
            .on_hover_text("Linear: brightness proportional to the field; honest sizes, small parts faint");
    });
    if game.field_linear {
        ui.add(egui::Slider::new(&mut game.field_gain_decades, 0.0..=8.0).text("gain (decades)"))
            .on_hover_text(
                "Amplifies the colours: full colour at 10^−gain of the full field's scale. \
                 The bar below states the values",
            );
    } else {
        ui.add(
            egui::Slider::new(&mut game.field_range_decades, 1.0..=14.0).text("range (decades)"),
        )
        .on_hover_text("How many decades below the full field's scale are still visible");
    }
    colour_bar(ui, game, radiation);
    // How large the part shown is (for the energy flow, each part has its own scale).
    if game.map == Some(MapMode::ParticleField)
        && (game.radiation_only || game.neglected_only)
        && !flow
    {
        let part = match game.field_quantity {
            FieldQuantity::Bz => radiation.part_b,
            FieldQuantity::E | FieldQuantity::S => radiation.part_e,
        };
        ui.label(
            egui::RichText::new(format!(
                "Largest values of this part: {} of the full field's scale.",
                percent(part)
            ))
            .color(egui::Color32::from_rgb(255, 210, 120)),
        );
    }
}

/// Controls of the energy-flow view: the part shown, tracers, the average over a period,
/// and what the view means (with its caveats).
pub(super) fn flow_controls(
    ui: &mut egui::Ui,
    game: &mut Game,
    radiation: &crate::radiation::RadiationView,
) {
    use crate::radiation::FlowPart;
    ui.horizontal_wrapped(|ui| {
        ui.label("Flow:")
            .on_hover_text("Which part of the energy flow: the view's field is the moving charges' field plus the rest");
        for (part, name, hint) in [
            (
                FlowPart::Total,
                "total",
                "All of it: the flow of the whole field shown",
            ),
            (
                FlowPart::Own,
                "own",
                "The moving charges' field alone: the energy moving with them, and their \
                 radiation leaving at c",
            ),
            (
                FlowPart::Exchange,
                "exchange",
                "The cross terms of the charges' field with the rest, (c²/4π)(E₁ × B₂ + \
                 E₂ × B₁): they carry the work the rest does on the charges, so they show \
                 where a particle's energy comes from",
            ),
            (
                FlowPart::External,
                "external",
                "The rest alone: the static sources (in the total view), antennas and waves",
            ),
        ] {
            ui.selectable_value(&mut game.flow_part, part, name)
                .on_hover_text(hint);
        }
    });
    ui.horizontal_wrapped(|ui| {
        ui.checkbox(&mut game.flow_tracers, "tracers")
            .on_hover_text(
                "Dots drifting with the energy's velocity S/u (u the energy density): never \
             faster than c, exactly c in a radiation field. The exchange terms have no \
             velocity of their own: there the dots move along S, at c at full colour",
            );
        if radiation.oscillates() {
            ui.checkbox(&mut game.flow_average, "average over a period")
                .on_hover_text(
                    "The antennas' and waves' flow averaged over their period (for several \
                     frequencies, each exactly), the charges held still: the energy that \
                     sloshes back and forth near an antenna drops out and the net outflow \
                     remains. Off: the flow at the moment shown",
                );
        }
    });
    explain(
        ui,
        "S = (c²/4π) E × B: where the field's energy flows",
        "S = (c²/4π) E × B, the flow of field energy (Jackson §6.7). Only its flux \
             through closed surfaces is unique: where the energy goes in between is \
             Poynting's convention. The map is a slice of a 3D flow; energy also leaves \
             the plane. Charges beside magnets keep a steady circulation (field momentum) \
             that carries energy nowhere. Near a particle part of the work done on it comes \
             from the exchange energy stored around it (a third at low speed), the rest \
             flows in.",
    );
}

/// The beam's energy budget at the animation time (the active disturbance's flight), in
/// units of the beam's launch kinetic energy: where the energy went.
pub(super) fn beam_energy(ui: &mut egui::Ui, game: &Game, level: &level::Level) {
    let d = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    let Some(p) = game.beams.get(d).and_then(crate::BeamView::shown) else {
        return;
    };
    let (Some(first), Some(last)) = (p.energy.first(), p.energy.last()) else {
        return;
    };
    let t = if game.animate { game.anim_time } else { last.t };
    // Linear interpolation between the recorded samples.
    let k = p
        .energy
        .partition_point(|e| e.t <= t)
        .clamp(1, p.energy.len().max(2) - 1);
    let e = if p.energy.len() < 2 || t >= last.t {
        *last
    } else {
        let (a, b) = (p.energy[k - 1], p.energy[k]);
        let f = if b.t > a.t {
            ((t - a.t) / (b.t - a.t)).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let mix = |x: f64, y: f64| x + (y - x) * f;
        physics::beam::EnergySample {
            t,
            kinetic: mix(a.kinetic, b.kinetic),
            potential: mix(a.potential, b.potential),
            interaction: mix(a.interaction, b.interaction),
            absorbed: mix(a.absorbed, b.absorbed),
            radiated: mix(a.radiated, b.radiated),
        }
    };
    let t0 = first.kinetic.max(1e-300);
    ui.label(egui::RichText::new("Energy of the whole beam").strong())
        .on_hover_text(
            "Where the beam's energy went, in units of its kinetic energy at launch T₀ \
         (all particles together)",
        );
    energy_bar(
        ui,
        "kinetic",
        e.kinetic / t0,
        egui::Color32::from_rgb(240, 200, 60),
    );
    energy_bar(
        ui,
        "potential",
        (e.potential - first.potential) / t0,
        egui::Color32::from_rgb(200, 90, 230),
    );
    if level.interacts() {
        energy_bar(
            ui,
            "mutual",
            (e.interaction - first.interaction) / t0,
            egui::Color32::from_rgb(230, 130, 90),
        );
    }
    energy_bar(
        ui,
        "absorbed",
        e.absorbed / t0,
        egui::Color32::from_rgb(160, 160, 160),
    );
    let rr = level.physics.radiation_reaction;
    if level.physics.c.is_some() {
        energy_bar(
            ui,
            "radiated",
            e.radiated / t0,
            egui::Color32::from_rgb(90, 200, 255),
        );
    }
    let mut total = e.kinetic
        + (e.potential - first.potential)
        + (e.interaction - first.interaction)
        + e.absorbed
        - first.kinetic;
    if rr {
        total += e.radiated;
    }
    energy_bar(
        ui,
        "total − T₀",
        total / t0,
        egui::Color32::from_rgb(120, 220, 120),
    );
    let note = match (level.physics.c.is_some(), rr, level.interacts()) {
        (false, _, _) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector. The total is conserved."
        }
        (true, true, _) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: carried away by the field (radiation reaction included, so \
             it is part of the total). Mutual: the Coulomb part of the particles' \
             interaction only (their magnetic and radiation field energy is not counted)."
        }
        (true, false, true) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: an estimate of what the particles would radiate, neglected \
             in the dynamics (not in the total). Mutual: the Coulomb part of the \
             interaction only."
        }
        (true, false, false) => {
            "Changes since launch. Absorbed: given to the bodies and the \
             detector; radiated: an estimate, neglected in the dynamics (not in the total)."
        }
    };
    ui.label(egui::RichText::new(note).small());
    ui.label(format!("t = {:.2}", e.t));
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]
/// The colour bar of the field views, with tick marks: the value of full colour, and
/// where tenths (linear) or decades (logarithmic) of it fall.
pub(super) fn colour_bar(
    ui: &mut egui::Ui,
    game: &Game,
    radiation: &crate::radiation::RadiationView,
) {
    use crate::radiation::{FieldQuantity, colour_value};
    let signed = game.field_quantity == FieldQuantity::Bz;
    let sat = match game.field_quantity {
        FieldQuantity::Bz => radiation.b_sat(),
        FieldQuantity::E => radiation.e_sat(),
        FieldQuantity::S => radiation.s_sat(game.flow_part),
    };
    let name = match game.field_quantity {
        FieldQuantity::Bz => "±B_z",
        FieldQuantity::E => "|E|",
        FieldQuantity::S => "|S|",
    };
    let linear = game.field_linear;
    let scale = if linear {
        10f64.powf(game.field_gain_decades)
    } else {
        10f64.powf(game.field_range_decades)
    };
    if sat <= 0.0 {
        ui.label(
            egui::RichText::new(format!(
                "{} = 0 everywhere: nothing to colour.",
                name.trim_start_matches('±')
            ))
            .small(),
        );
        return;
    }
    let width = ui.available_width().min(260.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 26.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let bar = egui::Rect::from_min_size(rect.min, egui::vec2(width, 10.0));
    let colour = |s: f64| {
        let a = (s.abs().powf(0.8) * 0.9) as f32;
        let base = if game.field_quantity == FieldQuantity::S {
            egui::Color32::from_rgb(140, 255, 158)
        } else if !signed {
            egui::Color32::from_rgb(255, 235, 140)
        } else if s >= 0.0 {
            egui::Color32::from_rgb(255, 150, 40)
        } else {
            egui::Color32::from_rgb(40, 190, 255)
        };
        egui::Color32::from_rgb(20, 22, 28).lerp_to_gamma(base, a)
    };
    let n = 64;
    for k in 0..n {
        let f = (f64::from(k) + 0.5) / f64::from(n);
        let s = if signed { 2.0 * f - 1.0 } else { f };
        let x0 = bar.left() + bar.width() * k as f32 / n as f32;
        let x1 = bar.left() + bar.width() * (k + 1) as f32 / n as f32;
        painter.rect_filled(
            egui::Rect::from_min_max(egui::pos2(x0, bar.top()), egui::pos2(x1, bar.bottom())),
            0.0,
            colour(s),
        );
    }
    // Ticks at values: position where the colour value of that field falls.
    let at = |v: f64| {
        let s = colour_value(v, sat, linear, scale);
        let f = if signed { 0.5 * (s + 1.0) } else { s };
        bar.left() + bar.width() * f as f32
    };
    let full = if linear { sat / scale } else { sat };
    let mut ticks: Vec<(f64, bool)> = Vec::new();
    if linear {
        for k in 1..=4 {
            ticks.push((full * f64::from(k) / 4.0, k == 4));
        }
    } else {
        let decades = game.field_range_decades.floor() as i32;
        for k in 0..=decades {
            ticks.push((sat * 10f64.powi(-k), k == 0));
        }
    }
    let text = ui.visuals().text_color();
    for (v, label) in ticks {
        for sign in if signed { vec![1.0, -1.0] } else { vec![1.0] } {
            let x = at(sign * v);
            painter.line_segment(
                [egui::pos2(x, bar.top()), egui::pos2(x, bar.bottom() + 3.0)],
                egui::Stroke::new(1.0, text),
            );
            if label && sign > 0.0 {
                painter.text(
                    egui::pos2(x.min(bar.right() - 30.0), bar.bottom() + 3.0),
                    egui::Align2::LEFT_TOP,
                    format!("{v:.1e}"),
                    egui::FontId::proportional(10.0),
                    text,
                );
            }
        }
    }
    let x0 = at(0.0);
    painter.line_segment(
        [
            egui::pos2(x0, bar.top()),
            egui::pos2(x0, bar.bottom() + 3.0),
        ],
        egui::Stroke::new(1.0, text),
    );
    ui.label(
        egui::RichText::new(if linear {
            format!(
                "Full colour at {} = {full:.2e} (the full field's scale ×10^−{:.1}); ticks every quarter.",
                name,
                game.field_gain_decades
            )
        } else {
            format!(
                "Full colour at {} = {sat:.2e}; ticks every decade below it.",
                name
            )
        })
        .small(),
    );
}

/// A radiation goal: what is required, what the flight radiates into it, and (with a
/// band) the spectrum over the goal's directions with the band marked. `system`: other
/// particles fly too, and the receiver sees all of them.
pub(super) fn radiation_goal(
    ui: &mut egui::Ui,
    goal: &level::RadiationGoal,
    radiation: Option<f64>,
    spectrum: &[(f64, f64)],
    spectrum_top: f64,
    system: bool,
) {
    let [axis, half] = goal.direction;
    let [lo, hi] = goal.energy;
    let band = goal.band.map_or("all frequencies".to_string(), |[a, b]| {
        format!("ω {a:.3}–{b:.3}")
    });
    // A steady receiver measures a mean power (per unit time), the others an energy.
    let unit = if goal.window.is_some() {
        "per sr per unit time"
    } else {
        "per sr"
    };
    let need = if lo > 0.0 {
        format!("{} – {} {unit}", fmt_si(lo), fmt_si(hi))
    } else {
        format!("below {} {unit}", fmt_si(hi))
    };
    if let Some([t1, t2]) = goal.window {
        explain(
            ui,
            &format!(
                "Steady receiver: into {axis:.0}° ± {half:.0}°, {band}, from t = {t1:.0} to {t2:.0}: {need}"
            ),
            &format!(
                "The receiver (the band outside the arena) stands far away and collects the \
                 light the particle scatters while a wave drives it: the mean power it \
                 receives over its window, after the start's transient has died away (its \
                 time is the arrival time less the light time from the centre of the arena). \
                 The particle need not arrive anywhere; it must keep scattering, bound, until \
                 the window has passed.{}",
                if system {
                    " It sees every particle: their fields add."
                } else {
                    ""
                }
            ),
        );
    } else {
        explain(
            ui,
            &format!("Radiation goal: into {axis:.0}° ± {half:.0}°, {band}: {need}"),
            &format!(
                "The receiver (the band outside the arena) stands far away: it collects the \
                 particle's radiation, not the particle. The particle itself must still end in \
                 its detector. The radiation leaves at the speed of light whenever the \
                 particle is accelerated, mostly in the direction it is heading, and reaches \
                 the receiver after the flight: switch the map to 'particle field' to watch it \
                 go.{}",
                if system {
                    " It sees every particle: their fields add, in phase they reinforce, in \
                     antiphase they cancel."
                } else {
                    ""
                }
            ),
        );
    }
    if let Some(e) = radiation {
        let ok = e >= lo && e <= hi;
        let color = if ok {
            egui::Color32::from_rgb(90, 240, 110)
        } else {
            egui::Color32::from_rgb(255, 170, 80)
        };
        let text = if goal.window.is_some() {
            format!("Received: {} {unit} (mean power)", fmt_si(e))
        } else {
            format!("Radiated into it: {} {unit}", fmt_si(e))
        };
        ui.colored_label(color, egui::RichText::new(text).small());
    }
    if spectrum.is_empty() {
        return;
    }
    // The spectrum d²I/dωdΩ from 0 to twice the band's top (the band shaded), or without a
    // band to where the flight's spectrum has fallen off.
    let (w, h) = (ui.available_width().min(300.0), 70.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, h), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_gray(25));
    let omega_max = spectrum_top.max(1e-300);
    #[allow(clippy::cast_possible_truncation)]
    let x_of = |omega: f64| rect.left() + (omega / omega_max) as f32 * rect.width();
    if let Some([a, b]) = goal.band {
        painter.rect_filled(
            egui::Rect::from_x_y_ranges(x_of(a)..=x_of(b), rect.y_range()),
            0.0,
            egui::Color32::from_rgba_unmultiplied(90, 240, 110, 40),
        );
    }
    let peak = spectrum.iter().map(|s| s.1).fold(0.0, f64::max);
    if peak > 0.0 {
        #[allow(clippy::cast_possible_truncation)]
        let points: Vec<egui::Pos2> = spectrum
            .iter()
            .map(|&(omega, v)| {
                egui::pos2(
                    x_of(omega),
                    rect.bottom() - 2.0 - (v / peak) as f32 * (h - 6.0),
                )
            })
            .collect();
        painter.add(egui::Shape::line(
            points,
            egui::Stroke::new(1.5, egui::Color32::from_rgb(240, 220, 120)),
        ));
    }
    ui.label(
        egui::RichText::new(format!(
            "Spectrum into these directions, ω from 0 to {omega_max:.3}{}; peak {} {unit} per unit ω{}",
            if goal.band.is_some() { " (band shaded)" } else { " (all of it counts)" },
            fmt_si(peak),
            if goal.window.is_some() { ", over the window" } else { "" }
        ))
        .small(),
    );
}
