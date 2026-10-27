//! A placement's results: the energy bars, the tube and beam levels' verdicts,
//! the exact verification against the quick preview.

use super::*;

/// A horizontal bar for an energy relative to T₀, on a scale of −2…2 T₀.
pub(super) fn energy_bar(ui: &mut egui::Ui, label: &str, value: f64, color: egui::Color32) {
    ui.horizontal(|ui| {
        ui.add_sized([70.0, 16.0], egui::Label::new(label));
        let (rect, _) = ui.allocate_exact_size(egui::vec2(170.0, 14.0), egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 2.0, egui::Color32::from_gray(35));
        let zero = rect.center().x;
        #[allow(clippy::cast_possible_truncation)]
        let x = zero + (value.clamp(-2.0, 2.0) as f32) * rect.width() / 4.0;
        let bar = egui::Rect::from_x_y_ranges(zero.min(x)..=zero.max(x), rect.y_range());
        painter.rect_filled(bar, 2.0, color);
        painter.line_segment(
            [
                egui::pos2(zero, rect.top()),
                egui::pos2(zero, rect.bottom()),
            ],
            egui::Stroke::new(1.0, egui::Color32::GRAY),
        );
        // Bounded width: a flight that falls onto a charge reaches 1e5 T₀, and a long
        // number widened the panel (and with it moved the board under the cursor).
        if value != 0.0 && !(1e-3..1e3).contains(&value.abs()) {
            ui.label(format!("{value:+.2e}"));
        } else {
            ui.label(format!("{value:+.3}"));
        }
    });
}

pub(super) fn tube_result(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    use level::tube::TubeStatus;
    let Some(spec) = &level.tube else {
        return;
    };
    ui.horizontal(|ui| {
        ui.checkbox(&mut game.animate, "Animate (A)")
            .on_hover_text("Play the preview: the particles leave the electrodes and fly");
        ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"))
            .on_hover_text("Playback speed: 1 plays 1 time unit per second");
    });
    if let Some(v) = &game.tube_preview
        && let Some(last) = v.run.frames.last()
    {
        let end = level.tube_display_end();
        ui.label(
            egui::RichText::new(if v.done {
                format!("t = {:.1} of {end:.0}", game.anim_time.min(last.t))
            } else {
                format!(
                    "t = {:.1} of {end:.0} (computed to {:.0}; playback waits for it)",
                    game.anim_time.min(last.t),
                    last.t
                )
            })
            .small(),
        );
    }
    ui.label(egui::RichText::new("Result").strong());
    if game.solved() {
        ui.label(
            egui::RichText::new("✔ SOLVED: every goal is met (verified)")
                .color(egui::Color32::from_rgb(90, 240, 110))
                .size(18.0),
        );
    }
    let colour = |s: TubeStatus| match s {
        TubeStatus::Met => (egui::Color32::from_rgb(90, 240, 110), "✔"),
        TubeStatus::Missed => (egui::Color32::from_rgb(255, 120, 90), "✖"),
        TubeStatus::Uncertain => (egui::Color32::from_rgb(255, 210, 90), "? too close to tell"),
    };
    let current = game.sent_revision;
    let view = game.tube_preview.as_ref().filter(|v| v.revision == current);
    let verdict = game
        .tube_verdict
        .as_ref()
        .filter(|(r, _)| *r == current)
        .map(|(_, v)| v);
    for (k, g) in spec.goals.iter().enumerate() {
        let whose = match g.carriers {
            // The cathode's carriers by their sign (electrons in every level so far).
            level::tube::Carriers::Cathode if spec.charge_per_mass < 0.0 => "electron current",
            level::tube::Carriers::Cathode => "current",
            level::tube::Carriers::Ions => "ion current",
            // A contact's lead: what it takes in less what it gives off.
            level::tube::Carriers::Net
                if spec
                    .metal_of(g.electrode)
                    .map_or(spec.charge_per_mass, |m| spec.metals[m].charge_per_mass)
                    < 0.0 =>
            {
                "net electron current"
            }
            level::tube::Carriers::Net => "net current",
        };
        ui.label(format!(
            "Goal: the {whose} into {} between {} and {}, averaged over t = {} to {}",
            level.electrode_ref(g.electrode),
            fmt_num(g.min),
            fmt_num(g.max),
            fmt_num(g.start),
            fmt_num(g.end)
        ));
        if let Some(v) = view {
            let run = &v.run;
            let current = run.currents.get(k).copied().flatten();
            ui.label(match current {
                Some(c) => format!("Preview: {}", fmt_num(c)),
                None => "Preview: computing…".into(),
            });
            ui.label(egui::RichText::new("Charge collected by the goal electrode:").small());
            if let Some(series) = run.collected.get(k) {
                circuit_plot(ui, series, &[g.start, g.end], "Q");
            }
            let quiet = g.carriers == level::tube::Carriers::Cathode
                && current.is_some()
                && run
                    .frames
                    .iter()
                    .take_while(|f| f.t <= g.end)
                    .all(|f| !f.kind.contains(&level::tube::KIND_NEGATIVE));
            if quiet {
                ui.label(
                    egui::RichText::new(
                        "No electrons leave the cathode: the field at its face pushes them \
                         back. A diode conducts one way only.",
                    )
                    .color(egui::Color32::from_rgb(255, 210, 120)),
                );
            }
        }
        match verdict.and_then(|v| v.currents.get(k).copied()) {
            Some(c) => {
                let (col, mark) = colour(c.status);
                ui.colored_label(
                    col,
                    format!(
                        "Verified: {} ± {} {mark}",
                        fmt_num(c.fine),
                        fmt_num(c.error)
                    ),
                );
            }
            None => {
                ui.label("Verified: computing…");
            }
        }
    }
    for (k, g) in spec.potential_goals.iter().enumerate() {
        let what = match g.measure {
            level::tube::Measure::Mean => "mean potential".to_string(),
            level::tube::Measure::Swing => "swing (half its peak-to-peak)".to_string(),
            level::tube::Measure::Amplitude => {
                format!("oscillation amplitude at ω = {}", fmt_num(g.omega))
            }
            level::tube::Measure::Phase => format!(
                "oscillation phase (degrees) at ω = {} relative to the drive",
                fmt_num(g.omega)
            ),
        };
        let against = g.reference.map_or(String::new(), |r| {
            format!(" against {}", level.electrode_ref(r))
        });
        ui.label(format!(
            "Goal: the {what} of {}{against} between {} and {}, over t = {} to {}",
            level.electrode_ref(g.electrode),
            fmt_num(g.min),
            fmt_num(g.max),
            fmt_num(g.start),
            fmt_num(g.end)
        ));
        if let Some(v) = view {
            ui.label(match v.run.potential.get(k).copied().flatten() {
                Some(p) => format!("Preview: {}", fmt_num(p)),
                None => "Preview: computing…".into(),
            });
        }
        match verdict.and_then(|v| v.potential.get(k).copied()) {
            Some(c) => {
                let (col, mark) = colour(c.status);
                ui.colored_label(
                    col,
                    format!(
                        "Verified: {} ± {} {mark}",
                        fmt_num(c.fine),
                        fmt_num(c.error)
                    ),
                );
            }
            None => {
                ui.label("Verified: computing…");
            }
        }
    }
    for (k, g) in spec.heat_goals.iter().enumerate() {
        ui.label(format!(
            "Goal: the heat of metal {} between {} and {} per time, over t = {} to {}",
            g.metal + 1,
            fmt_num(g.min),
            fmt_num(g.max),
            fmt_num(g.start),
            fmt_num(g.end)
        ))
        .on_hover_text(
            "The mean power the field gives the metal's carriers (per unit length along z): \
             what its lattice takes as heat, the Joule heat, in a steady state",
        );
        if let Some(v) = view {
            ui.label(match v.run.heat.get(k).copied().flatten() {
                Some(p) => format!("Preview: {}", fmt_num(p)),
                None => "Preview: computing…".into(),
            });
        }
        match verdict.and_then(|v| v.heat.get(k).copied()) {
            Some(c) => {
                let (col, mark) = colour(c.status);
                ui.colored_label(
                    col,
                    format!(
                        "Verified: {} ± {} {mark}",
                        fmt_num(c.fine),
                        fmt_num(c.error)
                    ),
                );
            }
            None => {
                ui.label("Verified: computing…");
            }
        }
    }
    // The radiated output: the field's energy flow out of the board (`em.rs`) against the
    // load's absorption.
    if let (Some(o), Some(v)) = (spec.output, view) {
        let load = spec
            .circuits
            .get(o.circuit)
            .and_then(|c| c.resistance)
            .unwrap_or(0.0);
        let t = if game.animate {
            game.anim_time
        } else {
            v.run.frames.last().map_or(0.0, |f| f.t)
        };
        let recent: Vec<f64> = v
            .run
            .frames
            .iter()
            .filter(|f| f.t <= t && f.t >= t - 2.0)
            .filter_map(|f| f.currents.get(o.circuit).map(|i| f64::from(*i).powi(2)))
            .collect();
        #[allow(clippy::cast_precision_loss)]
        let absorbed = load * recent.iter().sum::<f64>() / recent.len().max(1) as f64;
        let field_map = matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow));
        let text = match game.em_snapshot.as_ref().filter(|_| field_map) {
            Some(s) => format!(
                "Microwave output: {} per unit length leaves the board (the field's energy \
                 flow through its edge, over the last 2 time units); the load absorbs {} \
                 (R times the mean of I squared).",
                fmt_num(s.radiated),
                fmt_num(absorbed)
            ),
            None => format!(
                "Microwave output: the load absorbs {} (R times the mean of I squared). To \
                 see the waves the antenna radiates, choose the magnetic or energy flow map \
                 (the mouse wheel zooms, a drag or Ctrl + arrows pans, Z shows the whole \
                 board).",
                fmt_num(absorbed)
            ),
        };
        ui.label(egui::RichText::new(text).small());
    }
    // (The output's picture needed the game mutably.)
    let view = game.tube_preview.as_ref().filter(|v| v.revision == current);
    // The electrodes in circuits: their potentials over time (the goals' windows marked).
    let potentials: Vec<(usize, Vec<(f64, f64)>)> = view
        .map(|v| {
            spec.circuits
                .iter()
                // An electrode on a constant ideal source: a flat line, not plotted.
                .filter(|c| {
                    !(c.resistance == Some(0.0)
                        && c.inductance == 0.0
                        && matches!(c.source, level::Source::Dc { .. }))
                })
                .map(|c| {
                    // About 600 points over the run.
                    let every = (v.run.frames.len() / 600).max(1);
                    let samples = v
                        .run
                        .frames
                        .iter()
                        .step_by(every)
                        .filter_map(|f| f.potentials.get(c.electrode).map(|p| (f.t, f64::from(*p))))
                        .collect();
                    (c.electrode, samples)
                })
                .collect()
        })
        .unwrap_or_default();
    if !potentials.is_empty() {
        let spec = spec.clone();
        potential_plots(ui, game, &spec, &potentials);
    }
    let view = game.tube_preview.as_ref().filter(|v| v.revision == current);
    let verdict = game
        .tube_verdict
        .as_ref()
        .filter(|(r, _)| *r == current)
        .map(|(_, v)| v);
    for (j, p) in spec.projectiles.iter().enumerate() {
        if p.detector.is_none() {
            continue;
        }
        let preview = view.and_then(|v| v.run.arrivals.get(j).copied().flatten());
        let text = match verdict.and_then(|v| v.arrivals.iter().find(|a| a.0 == j)) {
            Some((_, _, fine, s)) => {
                let (col, mark) = colour(*s);
                ui.colored_label(
                    col,
                    format!(
                        "Projectile {}: {} {mark}",
                        j + 1,
                        match fine {
                            Some(t) => format!("reaches its detector at t = {t:.1} (verified)"),
                            None => "does not reach its detector (verified)".into(),
                        }
                    ),
                );
                continue;
            }
            None => match preview {
                Some(t) => format!(
                    "Projectile {}: reaches its detector at t = {t:.1}; verifying…",
                    j + 1
                ),
                None => format!("Projectile {}: computing…", j + 1),
            },
        };
        ui.label(text);
    }
    if verdict.is_some_and(|v| v.capped) {
        ui.colored_label(
            egui::Color32::from_rgb(255, 170, 80),
            format!(
                "Too many particles: a run reached {} and emission paused, so the verdict is uncertain.",
                physics::tube::MAX_PARTICLES
            ),
        );
    }
    if let Some(v) = view {
        let shown = game
            .tube_frame()
            .and_then(|k| v.run.frames.get(k))
            .map_or(0, |f| f.x.len());
        ui.label(egui::RichText::new(format!("{shown} particles in flight now")).small());
        // The Joule heat: the energy the metals' lattices took over the last time unit.
        if !spec.metals.is_empty()
            && let Some(k) = game.tube_frame()
            && let Some(f1) = v.run.frames.get(k)
        {
            // Over the last 10 time units: each collision's heat fluctuates by about kT,
            // over one unit the rate was uncertain by a third.
            let k0 = v.run.frames.partition_point(|f| f.t < f1.t - 10.0);
            if let Some(f0) = v.run.frames.get(k0).filter(|f0| f1.t > f0.t) {
                explain(
                    ui,
                    &format!(
                        "Joule heat: {} per time",
                        fmt_num((f1.heat - f0.heat) / (f1.t - f0.t))
                    ),
                    "The heat the metal's lattice takes, per unit length along z per time \
                     unit, over the last 10 time units (each collision's share fluctuates by \
                     about kT). The lattice takes it at the collisions: about nine tenths of \
                     the supplies' I V; the carriers bring the rest into the contact they \
                     reach, the energy the field gave them since their last collision.",
                );
            }
        }
    }
    explain(
        ui,
        if spec.metals.is_empty() {
            "A z-invariant tube: what is drawn"
        } else {
            "A metal (the Drude model): what is drawn"
        },
        if spec.metals.is_empty() {
            "A z-invariant tube: every electrode is a long prism, the electrons (blue), ions \
             (orange) and projectiles (large) are long lines of charge, and currents are per \
             unit length. Emitters release as much as their field allows: the particles' own \
             charge limits the current."
        } else {
            "A z-invariant metal (the Drude model): its carriers (blue; a few traced) fly at \
             the thermal speed and bounce off the ion cores (dots), a zigzag; the field adds \
             a slow drift along it, the current. The chevrons across the metal move with its \
             carriers as a whole (their mean displacement): the drift, a few tenths of a \
             unit per time unit beside the thermal speed, faint when it is only the thermal \
             motion's noise. Each collision hands the lattice the energy the field gave: the \
             Joule heat. Contacts take in the carriers that reach them and give off as many \
             as keep the metal next to them neutral (a traced carrier a contact takes is \
             replaced by another, so they stay as many). A classical model of a metal (its real \
             electrons are a quantum gas): see the physics notes."
        },
    );
}

/// Beam levels: verified transmission of every beam shot in every flight against its
/// requirement, and the diagnostics of the selected flight.
pub(super) fn beam_result(ui: &mut egui::Ui, game: &mut Game, level: &level::Level) {
    ui.horizontal(|ui| {
        ui.checkbox(&mut game.animate, "Animate (A)").on_hover_text(
            "Play the flights in time: the particles move along their paths and the \
             time-dependent maps (waves, particle fields) follow",
        );
        ui.add(egui::Slider::new(&mut game.playback_speed, 0.05..=4.0).text("speed"))
            .on_hover_text("Playback speed: 1 plays 4 time units per second");
    });
    ui.label(egui::RichText::new("Result").strong());
    if game.solved() {
        ui.label(
            egui::RichText::new("✔ SOLVED: every beam delivers its share (verified)")
                .color(egui::Color32::from_rgb(90, 240, 110))
                .size(18.0),
        );
    }
    let flights = game.beams.len().max(1);
    for d in 0..flights {
        if flights > 1 {
            ui.label(egui::RichText::new(format!("Disturbance {}", d + 1)).small());
        }
        for (s, shot) in level.shots.iter().enumerate() {
            let need = shot.beam.map_or(1.0, |b| b.transmission);
            let text = match game.beam_transmission(d, s) {
                None => format!("Shot {}: computing…", s + 1),
                Some((ok, n)) => {
                    #[allow(clippy::cast_precision_loss)]
                    let enough = ok as f64 >= need * n as f64 - 1e-9;
                    format!(
                        "Shot {}: {ok} of {n} arrive, verified (need {:.0} %) {}",
                        s + 1,
                        need * 100.0,
                        if enough { "✔" } else { "✖" }
                    )
                }
            };
            ui.colored_label(shot_color32(level, s), text);
        }
        // Goal particles: the level's free particles with a detector.
        let n = level.shots.len();
        for (k, f) in level.free_particles.iter().enumerate() {
            if f.detector.is_none() {
                continue;
            }
            let text = match game.beam_transmission(d, n + k) {
                None => format!("Particle {}: computing…", k + 1),
                Some((1, _)) => format!("Particle {}: arrives, verified ✔", k + 1),
                Some(_) => format!("Particle {}: does not arrive (verified) ✖", k + 1),
            };
            let c = crate::draw::free_goal_color(k).to_srgba();
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let color = egui::Color32::from_rgb(
                (c.red * 255.0) as u8,
                (c.green * 255.0) as u8,
                (c.blue * 255.0) as u8,
            );
            ui.colored_label(color, text);
        }
    }
    let d = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    let view = game.beams.get(d);
    if let Some(p) = view.and_then(crate::BeamView::shown) {
        // Lost: particles with a goal (the shots' and the level's free particles with a
        // detector) that do not arrive; the player's free charges have none.
        let n_shots = level.shots.len();
        let has_goal = |s: usize| {
            s < n_shots
                || (s != usize::MAX
                    && level
                        .free_particles
                        .get(s - n_shots)
                        .is_some_and(|f| f.detector.is_some()))
        };
        let lost = p
            .outcomes
            .iter()
            .zip(&p.shots)
            .filter(|(o, s)| has_goal(**s) && **o != Outcome::Arrived)
            .count();
        // Which flight the views show: the quick preview, or the verdict's own flight
        // once it has arrived (exact at finite c).
        let verdict_shown = view.is_some_and(|v| v.verdict_flight.is_some());
        let energy = if p.energy_rel_error.is_nan() {
            // Finite c: the particles exchange energy with the field.
            let mut s = format!(
                "largest radiated share of a particle {:.1e} ({})",
                p.radiated_max,
                if level.physics.radiation_reaction {
                    "included"
                } else {
                    "neglected"
                }
            );
            if p.retardation_max > 0.0 {
                s += &format!(
                    "; approximation of the fields (quasi-static), estimated error {:.1e} of the interaction",
                    p.retardation_max
                );
            }
            s
        } else {
            format!("energy drift of the whole beam {:.1e}", p.energy_rel_error)
        };
        // Radiation goals: every charge's radiation together.
        for r in &p.radiation {
            if let Some(goal) = level
                .shots
                .get(r.shot)
                .and_then(|s| s.detector.acceptance)
                .and_then(|a| a.radiation)
            {
                radiation_goal(ui, &goal, Some(r.value), &r.spectrum, r.spectrum_top, true);
            }
        }
        ui.label(
            egui::RichText::new(format!(
                "{} particles{}; {lost} lost{}; {energy}.",
                p.paths.len(),
                match (level.interacts(), level.physics.c) {
                    (false, _) => ", not interacting",
                    (true, None) => ", interacting",
                    (true, Some(_)) if p.retarded => ", interacting (exact retarded fields)",
                    (true, Some(_)) => ", interacting (quasi-static fields)",
                },
                if verdict_shown { "" } else { " in the preview" },
            ))
            .small(),
        );
        let shown_text = match view {
            Some(v) if verdict_shown => v.preview.as_ref().map(|quick| {
                let (worst, differ) = quick_vs_verdict(quick, p);
                let mut s = String::from(if p.retarded {
                    "Shown in every view: the verdict's own flight (exact retarded fields). "
                } else {
                    "Shown in every view: the verdict's own flight (metal at its verification \
                     resolution). "
                });
                // The sandbox's exact option flies the preview retarded too.
                s += if p.retarded && !quick.retarded {
                    "The quick preview (quasi-static)"
                } else {
                    "The preview (looser tolerance)"
                };
                s += &format!(" ended within {worst:.1e} cells of it");
                if quick.retardation_max > 0.0 {
                    s += &format!(
                        " (its estimated error: {:.1e} of the interaction)",
                        quick.retardation_max
                    );
                }
                if differ > 0 {
                    s += &format!("; {differ} particles ended differently");
                }
                s + "."
            }),
            Some(v) if v.verified.is_none() && exact_verdict(level) => Some(
                "Shown: the quick preview (quasi-static fields). The exact flight, the \
                 verdict's, replaces it in every view when it is ready."
                    .to_string(),
            ),
            _ => None,
        };
        if let Some(text) = shown_text {
            ui.label(egui::RichText::new(text).small());
        }
    }
    beam_energy(ui, game, level);
    ui.separator();
}

/// Whether the verdict of a beam level flies the exact retarded interaction while the
/// preview is quasi-static (`Level::verification_beam_scenarios`).
pub(super) fn exact_verdict(level: &level::Level) -> bool {
    level.has_beams()
        && level.interacts()
        && level.physics.c.is_some()
        && !level.physics.beam_retarded
}

/// How far the quick preview's end points are from the verdict flight's (particles with
/// the same outcome), and how many particles ended differently.
pub(super) fn quick_vs_verdict(
    quick: &crate::worker::BeamPreview,
    verdict: &crate::worker::BeamPreview,
) -> (f64, usize) {
    let mut worst: f64 = 0.0;
    let mut differ = 0;
    for (i, (a, b)) in quick.paths.iter().zip(&verdict.paths).enumerate() {
        if quick.outcomes.get(i) != verdict.outcomes.get(i) {
            differ += 1;
        } else if let (Some(x), Some(y)) = (a.last(), b.last()) {
            worst = worst.max((x.1 - y.1).length());
        }
    }
    (worst, differ)
}

/// The drawing colour of a shot, for egui.
pub(super) fn shot_color32(level: &level::Level, shot: usize) -> egui::Color32 {
    let c = crate::draw::shot_color(level, shot).to_srgba();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let b = |v: f32| (v.clamp(0.0, 1.0) * 255.0) as u8;
    egui::Color32::from_rgb(b(c.red), b(c.green), b(c.blue))
}
