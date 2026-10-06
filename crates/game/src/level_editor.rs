//! Property editors for every part of a level (sandbox).
//!
//! **Completeness invariant:** everything a level file can contain must be editable here
//! (or by a sandbox tool), so any shipped level can be rebuilt by hand. Every level type
//! is destructured *exhaustively* (no `..`): adding a field to the level format fails to
//! compile until it gets a control here, or is explicitly marked as not user-editable
//! with a reason.

use bevy_egui::egui;
use level::beam::{BeamSpec, Distribution};
use level::{
    Cloud, Coil, Conductor, ConductorBias, Detector, DetectorAcceptance, Disturbance, Drive,
    Electrode, Element, ElementKind, FreeParticle, Grid, Launch, Level, Limits, Node, ParticleSpec,
    RadiationGoal, Region2, Shot, Source, Switch, TolerancesSpec, Wave, WorldPhysics,
};

use crate::ui::{fmt_si, parse_si};

// Ranges the editors allow; `check_editable` uses the same constants, so every shipped
// level is checked to lie within reach of the editors.
const GRID_CELLS: std::ops::RangeInclusive<u32> = 4..=200;
const SUBDIVISION: std::ops::RangeInclusive<u32> = 1..=8;
const MAX_POSITIVE: f64 = 1e12;
const MIN_POSITIVE: f64 = 1e-12;
const MAX_RADIUS: f64 = 5.0;
/// Largest beam (particles; the interaction costs grow as N²).
const MAX_BEAM: u32 = 200;
const MAX_COUNT: u32 = 20;
const TOLERANCE: std::ops::RangeInclusive<f64> = 1e-18..=1.0;
/// Width of list text fields, so that grid rows fit the panel.
const LIST_WIDTH: f32 = 100.0;

/// Text buffers for list fields (edited as comma-separated text).
#[derive(Default)]
pub struct EditTexts {
    pub magnitudes: String,
    pub magnet_strengths: String,
    pub antenna_amplitudes: String,
    pub antenna_omegas: String,
    pub plate_voltages: String,
    pub supply_voltages: String,
    pub free_charges: String,
    pub free_speeds: String,
}

pub fn list_to_text(m: &[f64]) -> String {
    m.iter()
        .map(|v| format!("{v:e}"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// A list of signed values (potentials; 0 allowed), without duplicates.
fn parse_signed_list(text: &str) -> Vec<f64> {
    let mut out: Vec<f64> = text
        .split(',')
        .filter_map(|t| parse_si(t.trim()))
        .filter(|v| v.is_finite())
        .collect();
    out.sort_by(f64::total_cmp);
    out.dedup();
    out
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
fn non_negative(ui: &mut egui::Ui, v: &mut f64, speed: f64) -> bool {
    ui.add(
        egui::DragValue::new(v)
            .speed(speed)
            .range(0.0..=MAX_POSITIVE)
            .custom_formatter(|v, _| fmt_si(v))
            .custom_parser(parse_si),
    )
    .has_focus()
}

/// The circuit driving an electrode or a circular coil (PHYSICS.md §2.10): a source
/// through a resistance, with a series inductance or a switch (electrodes) or a series
/// capacitance (coils).
fn drive_editor(
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

/// A drive within the editors' ranges.
fn drive_editable(d: &Drive) -> bool {
    let Drive {
        source,
        resistance,
        inductance,
        capacitance,
        switch,
    } = d;
    let non_negative = |v: &f64| (0.0..=MAX_POSITIVE).contains(v);
    let source_ok = match source {
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
    };
    source_ok
        && (MIN_POSITIVE..=MAX_POSITIVE).contains(resistance)
        && non_negative(inductance)
        && non_negative(capacitance)
        && switch.is_none_or(|Switch { closed: _, at }| non_negative(&at))
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
        antenna_radius,
        rf_omega,
        radiation_reaction,
        beam_interaction,
        beam_retarded,
        instant_drain,
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
    // At most two values per row: a row wider than the panel makes egui draw a stray
    // full-height line.
    focus |= row(ui, "Radii: charge, magnet", |ui| {
        positive(ui, charge_radius, 0.01, MAX_RADIUS)
            | positive(ui, magnet_radius, 0.01, MAX_RADIUS)
    });
    focus |= row(ui, "Radii: wire, antenna", |ui| {
        positive(ui, wire_radius, 0.01, MAX_RADIUS) | positive(ui, antenna_radius, 0.01, MAX_RADIUS)
    });
    focus |= row(ui, "Antenna RF ω", |ui| {
        ui.add(
            egui::DragValue::new(rf_omega)
                .speed(0.01)
                .range(0.0..=MAX_POSITIVE),
        )
        .has_focus()
    });
    row(ui, "Radiation reaction", |ui| {
        ui.checkbox(radiation_reaction, "included")
            .on_hover_text("Landau–Lifshitz force: the particle loses the energy it radiates");
        false
    });
    row(ui, "Beam interaction", |ui| {
        ui.checkbox(beam_interaction, "particles interact").on_hover_text(
            "Coulomb for c = ∞; for finite c quasi-static: the Liénard–Wiechert fields of each particle's present state continued back with constant acceleration (exact in the velocities, first order in the accelerations)",
        );
        false
    });
    row(ui, "Beam fields", |ui| {
        ui.checkbox(beam_retarded, "exact retarded (slow)").on_hover_text(
            "Finite c: the exact Liénard–Wiechert fields at the retarded times instead of the quasi-static ones (steps shorter than the light time between particles)",
        );
        false
    });
    row(ui, "Beam detectors", |ui| {
        ui.checkbox(instant_drain, "instant drain").on_hover_text(
            "Off (default): a detector is the mouth of a deep grounded cup; an absorbed particle flies on into it while the cup screens its charge (it fades over the entry time). On: its charge vanishes at once",
        );
        false
    });
    focus |= row(ui, "Time limit", |ui| {
        positive(ui, t_max, 1.0, MAX_POSITIVE)
    });
    focus |= row(ui, "Tolerance prev., ver.", |ui| {
        sci(ui, preview) | sci(ui, verify)
    });
    focus
}

fn edit_shot(ui: &mut egui::Ui, s: &mut Shot, grid: &Grid) -> bool {
    let Shot {
        particle,
        launch,
        detector,
        beam,
    } = s;
    let ParticleSpec {
        charge,
        mass,
        radius,
        moment,
    } = particle;
    let Launch {
        node: launch_node,
        direction,
        kinetic_energy,
        time,
    } = launch;
    let Detector {
        min,
        max,
        acceptance,
    } = detector;
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
    focus |= row(ui, "Magnetic moment m_z", |ui| si(ui, moment, 1e-8));
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
    focus |= row(ui, "Launch time (lab)", |ui| {
        ui.add(egui::DragValue::new(time).speed(0.05)).has_focus()
    });
    focus |= row(ui, "Detector corner 1", |ui| node(ui, min, grid));
    focus |= row(ui, "Detector corner 2", |ui| node(ui, max, grid));
    // Optional conditions on the arriving particle.
    focus |= edit_acceptance(ui, acceptance, true);
    focus |= edit_beam(ui, beam);
    focus
}

/// A shot fired as a beam (`level::beam`): count, spreads, distribution, required
/// transmission and seed.
/// Optional conditions on a particle entering a detector or gate: direction and energy,
/// and (shot detectors only, `radiation_allowed`) the flight's radiation.
fn edit_acceptance(
    ui: &mut egui::Ui,
    acceptance: &mut Option<DetectorAcceptance>,
    radiation_allowed: bool,
) -> bool {
    let mut focus = false;
    let mut acc = acceptance.unwrap_or_default();
    let DetectorAcceptance {
        direction,
        kinetic,
        radiation,
    } = &mut acc;
    focus |= row(ui, "Accept direction", |ui| {
        let mut on = direction.is_some();
        ui.checkbox(&mut on, "");
        match (on, direction.is_some()) {
            (true, false) => *direction = Some([0.0, 10.0]),
            (false, true) => *direction = None,
            _ => {}
        }
        let mut f = false;
        if let Some([axis, half]) = direction {
            f |= ui
                .add(egui::DragValue::new(axis).speed(1.0).suffix("°"))
                .on_hover_text("Axis, degrees from +x")
                .has_focus();
            ui.label("±");
            f |= ui
                .add(
                    egui::DragValue::new(half)
                        .speed(0.5)
                        .range(0.0..=180.0)
                        .suffix("°"),
                )
                .has_focus();
        }
        f
    });
    focus |= row(ui, "Accept energy", |ui| {
        let mut on = kinetic.is_some();
        ui.checkbox(&mut on, "");
        match (on, kinetic.is_some()) {
            (true, false) => *kinetic = Some([0.0, 1.0]),
            (false, true) => *kinetic = None,
            _ => {}
        }
        let mut f = false;
        if let Some([lo, hi]) = kinetic {
            f |= si(ui, lo, 0.01);
            ui.label("–");
            f |= si(ui, hi, 0.01);
        }
        f
    });
    if radiation_allowed {
        focus |= edit_radiation(ui, radiation);
    }
    *acceptance = (acc != DetectorAcceptance::default()).then_some(acc);
    focus
}

/// A radiation goal: arc of directions, optional frequency band, energy per steradian;
/// or, with a receiving window, a steady receiver's power per steradian.
fn edit_radiation(ui: &mut egui::Ui, radiation: &mut Option<RadiationGoal>) -> bool {
    let mut focus = row(ui, "Radiation goal", |ui| {
        let mut on = radiation.is_some();
        ui.checkbox(&mut on, "")
            .on_hover_text("Energy per steradian radiated into an arc of directions (far zone)");
        match (on, radiation.is_some()) {
            (true, false) => {
                *radiation = Some(RadiationGoal {
                    direction: [0.0, 10.0],
                    band: None,
                    energy: [0.0, 1.0],
                    abrupt_stop: false,
                    window: None,
                });
            }
            (false, true) => *radiation = None,
            _ => {}
        }
        false
    });
    let Some(RadiationGoal {
        direction: [axis, half],
        band,
        energy: [lo, hi],
        abrupt_stop,
        window,
    }) = radiation
    else {
        return focus;
    };
    focus |= row(ui, "  directions", |ui| {
        let f = ui
            .add(egui::DragValue::new(axis).speed(1.0).suffix("°"))
            .on_hover_text("Axis, degrees from +x")
            .has_focus();
        ui.label("±");
        f | ui
            .add(
                egui::DragValue::new(half)
                    .speed(0.5)
                    .range(0.0..=180.0)
                    .suffix("°"),
            )
            .has_focus()
    });
    focus |= row(ui, "  band ω", |ui| {
        let mut on = band.is_some();
        ui.checkbox(&mut on, "");
        match (on, band.is_some()) {
            (true, false) => *band = Some([0.5, 1.5]),
            (false, true) => *band = None,
            _ => {}
        }
        let mut f = false;
        if let Some([a, b]) = band {
            f |= si(ui, a, 0.01);
            ui.label("–");
            f |= si(ui, b, 0.01);
        }
        f
    });
    focus |= row(ui, "  steady from t", |ui| {
        let mut on = window.is_some();
        ui.checkbox(&mut on, "").on_hover_text(
            "A steady receiver: the mean power per steradian it receives over this window of \
             its time, instead of the flight's whole energy (the particle need not arrive)",
        );
        match (on, window.is_some()) {
            (true, false) => *window = Some([100.0, 200.0]),
            (false, true) => *window = None,
            _ => {}
        }
        let mut f = false;
        if let Some([a, b]) = window {
            f |= si(ui, a, 0.01);
            ui.label("–");
            f |= si(ui, b, 0.01);
        }
        f
    });
    let label = if window.is_some() {
        "  power / sr"
    } else {
        "  energy / sr"
    };
    focus |= row(ui, label, |ui| {
        si(ui, lo, 0.01) | {
            ui.label("–");
            si(ui, hi, 0.01)
        }
    });
    row(ui, "  abrupt stop", |ui| {
        ui.checkbox(abrupt_stop, "")
            .on_hover_text(
                "The detector is a target that stops the particle at once; the stop's radiation counts (needs a band)",
            );
        false
    });
    focus
}

/// Whether an acceptance's values are in the editor's range (a radiation goal only where
/// `radiation_allowed`: shot detectors).
fn acceptance_in_range(a: &DetectorAcceptance, radiation_allowed: bool) -> bool {
    let DetectorAcceptance {
        direction,
        kinetic,
        radiation,
    } = a;
    let window = |[lo, hi]: [f64; 2]| lo.is_finite() && hi.is_finite() && lo < hi;
    direction.is_none_or(|[x, h]| x.is_finite() && (0.0..=180.0).contains(&h))
        && kinetic.is_none_or(window)
        && radiation.is_none_or(|r| {
            let RadiationGoal {
                direction: [x, h],
                band,
                energy,
                abrupt_stop,
                window: receiving,
            } = r;
            radiation_allowed
                && (!abrupt_stop || band.is_some())
                && !(abrupt_stop && receiving.is_some())
                && receiving.is_none_or(|w| window(w) && w[0] >= 0.0)
                && x.is_finite()
                && (0.0..=180.0).contains(&h)
                && band.is_none_or(|b| window(b) && b[0] >= 0.0)
                && window(energy)
                && energy[0] >= 0.0
        })
}

fn edit_free_particles(ui: &mut egui::Ui, list: &mut Vec<FreeParticle>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, f) in list.iter_mut().enumerate() {
        let FreeParticle {
            particle,
            node: at,
            velocity,
            detector,
        } = f;
        let ParticleSpec {
            charge,
            mass,
            radius,
            moment,
        } = particle;
        ui.push_id(("free", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Particle {}", i + 1));
                focus |= node(ui, at, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            focus |= row(ui, "  q, m, radius", |ui| {
                si(ui, charge, 1e-8) | positive(ui, mass, 0.01, MAX_POSITIVE) | {
                    ui.add(
                        egui::DragValue::new(radius)
                            .speed(0.01)
                            .range(0.0..=MAX_RADIUS),
                    )
                    .has_focus()
                }
            });
            focus |= row(ui, "  moment m_z", |ui| si(ui, moment, 1e-8));
            focus |= row(ui, "  velocity (x, y)", |ui| {
                ui.add(egui::DragValue::new(&mut velocity[0]).speed(0.01))
                    .has_focus()
                    | ui.add(egui::DragValue::new(&mut velocity[1]).speed(0.01))
                        .has_focus()
            });
            let mut goal = detector.is_some();
            row(ui, "  detector", |ui| {
                ui.checkbox(&mut goal, "must arrive");
                false
            });
            match (goal, detector.is_some()) {
                (true, false) => {
                    *detector = Some(Detector {
                        min: *at,
                        max: [at[0] + 1, at[1] + 1, 0],
                        acceptance: None,
                    });
                }
                (false, true) => *detector = None,
                _ => {}
            }
            if let Some(Detector {
                min,
                max,
                acceptance,
            }) = detector
            {
                focus |= row(ui, "  corner 1", |ui| node(ui, min, grid));
                focus |= row(ui, "  corner 2", |ui| node(ui, max, grid));
                focus |= edit_acceptance(ui, acceptance, false);
            }
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ free particle").clicked() {
        let m = grid.max_node();
        list.push(FreeParticle {
            particle: ParticleSpec {
                charge: 1e-6,
                mass: 1.0,
                radius: 0.3,
                moment: 0.0,
            },
            node: [m[0] / 2, m[1] / 2, 0],
            velocity: [0.0, 0.0],
            detector: None,
        });
    }
    focus
}

/// The instrument's gates (stages every flight must pass, in order): corners and
/// optional conditions.
fn edit_gates(ui: &mut egui::Ui, list: &mut Vec<Detector>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, gate) in list.iter_mut().enumerate() {
        let Detector {
            min,
            max,
            acceptance,
        } = gate;
        ui.push_id(("gate", i), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Gate {}", i + 1));
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
            });
            focus |= row(ui, "  corner 1", |ui| node(ui, min, grid));
            focus |= row(ui, "  corner 2", |ui| node(ui, max, grid));
            focus |= edit_acceptance(ui, acceptance, false);
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
    if list.len() < MAX_COUNT as usize && ui.small_button("+ gate").clicked() {
        let m = grid.max_node();
        list.push(Detector {
            min: [m[0] / 2, m[1] / 2 - 2, 0],
            max: [m[0] / 2 + 1, m[1] / 2 + 2, 0],
            acceptance: None,
        });
    }
    focus
}

fn edit_beam(ui: &mut egui::Ui, beam: &mut Option<BeamSpec>) -> bool {
    let mut on = beam.is_some();
    row(ui, "Beam", |ui| {
        ui.checkbox(&mut on, "fire as a beam")
            .on_hover_text("Many particles with spreads; the goal is a verified transmission");
        false
    });
    match (on, beam.is_some()) {
        (true, false) => *beam = Some(BeamSpec::default()),
        (false, true) => *beam = None,
        _ => {}
    }
    let Some(b) = beam else {
        return false;
    };
    let BeamSpec {
        count,
        energy_spread,
        angle_spread_deg,
        width,
        length,
        distribution,
        transmission,
        seed,
    } = b;
    let mut focus = false;
    focus |= row(ui, "  particles, transmission", |ui| {
        ui.add(egui::DragValue::new(count).range(1..=MAX_BEAM))
            .has_focus()
            | ui.add(
                egui::DragValue::new(transmission)
                    .speed(0.01)
                    .range(0.01..=1.0),
            )
            .on_hover_text("Fraction that must arrive, verified")
            .has_focus()
    });
    focus |= row(ui, "  spread: energy, angle", |ui| {
        ui.add(
            egui::DragValue::new(energy_spread)
                .speed(0.001)
                .range(0.0..=0.3),
        )
        .on_hover_text("Relative to T₀ (Gaussian: σ; uniform: half-width)")
        .has_focus()
            | ui.add(
                egui::DragValue::new(angle_spread_deg)
                    .speed(0.1)
                    .range(0.0..=45.0)
                    .suffix("°"),
            )
            .has_focus()
    });
    focus |= row(ui, "  spread: width, length", |ui| {
        ui.add(egui::DragValue::new(width).speed(0.01).range(0.0..=5.0))
            .on_hover_text("Across the direction, cells")
            .has_focus()
            | ui.add(egui::DragValue::new(length).speed(0.01).range(0.0..=5.0))
                .on_hover_text("Along the direction (bunch length), cells")
                .has_focus()
    });
    row(ui, "  distribution", |ui| {
        ui.selectable_value(distribution, Distribution::Gaussian, "Gaussian ±3σ");
        ui.selectable_value(distribution, Distribution::Uniform, "uniform");
        false
    });
    focus |= row(ui, "  sample seed", |ui| {
        ui.add(egui::DragValue::new(seed)).has_focus()
    });
    focus
}

fn kind_combo(ui: &mut egui::Ui, id: usize, kind: &mut ElementKind) {
    let text = |k: ElementKind| match k {
        ElementKind::Charge => "charge",
        ElementKind::Magnet => "magnet μ",
        ElementKind::Antenna => "antenna p₀",
        ElementKind::Plate => "plate V",
        ElementKind::Supply => "supply V",
        ElementKind::Free => "free q",
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

fn edit_elements(ui: &mut egui::Ui, elements: &mut Vec<Element>, grid: &Grid) -> bool {
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
                    rate,
                    drive,
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
                } => {
                    ui.horizontal(|ui| {
                        ui.label("polygon, κ");
                        focus |= si(ui, kappa, 100.0);
                        ui.label("ramp dκ/dt");
                        focus |= si(ui, rate, 1.0);
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

fn edit_conductors(ui: &mut egui::Ui, list: &mut Vec<Conductor>, grid: &Grid) -> bool {
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

fn edit_clouds(ui: &mut egui::Ui, list: &mut Vec<Cloud>, grid: &Grid) -> bool {
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

/// Bias selector shared by metal spheres and electrodes.
fn bias_editor(ui: &mut egui::Ui, id: (&str, usize), bias: &mut ConductorBias) -> bool {
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

fn edit_electrodes(ui: &mut egui::Ui, list: &mut Vec<Electrode>, grid: &Grid) -> bool {
    let mut focus = false;
    let mut remove = None;
    for (i, e) in list.iter_mut().enumerate() {
        let Electrode {
            center,
            length,
            thickness,
            height,
            angle_deg,
            bias,
            tunable,
            drive,
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
            bias: ConductorBias::Grounded,
            tunable: false,
            drive: None,
        });
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
        max_antennas,
        antenna_amplitudes,
        antenna_omegas,
        continuous,
        region,
        max_plates,
        plate_voltages,
        plate,
        supply_voltages,
        max_free,
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
        conductors,
        clouds,
        free_particles,
        electrodes,
        gates,
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
    egui::CollapsingHeader::new(format!("Electrodes ({})", electrodes.len()))
        .show(ui, |ui| focus |= edit_electrodes(ui, electrodes, &g));
    egui::CollapsingHeader::new(format!("Metal spheres ({})", conductors.len()))
        .show(ui, |ui| focus |= edit_conductors(ui, conductors, &g));
    egui::CollapsingHeader::new(format!("Charge clouds ({})", clouds.len()))
        .show(ui, |ui| focus |= edit_clouds(ui, clouds, &g));
    egui::CollapsingHeader::new(format!("Free particles ({})", free_particles.len()))
        .show(ui, |ui| {
            focus |= edit_free_particles(ui, free_particles, &g)
        })
        .header_response
        .on_hover_text(
            "Particles that move and interact with the shots (targets, partners); with a \
             detector they must arrive too",
        );
    egui::CollapsingHeader::new(format!("Disturbances ({})", disturbances.len()))
        .show(ui, |ui| focus |= edit_disturbances(ui, disturbances));
    egui::CollapsingHeader::new(format!("Gates ({})", gates.len()))
        .show(ui, |ui| focus |= edit_gates(ui, gates, &g))
        .header_response
        .on_hover_text("Stages every flight must pass, in order, before its detector counts");
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
        conductors,
        clouds,
        free_particles,
        electrodes,
        gates,
    } = level;
    if gates.len() > MAX_COUNT as usize {
        return Err("too many gates".into());
    }
    for g in gates {
        let Detector {
            min,
            max,
            acceptance,
        } = g;
        let acc_ok = acceptance.is_none_or(|a| acceptance_in_range(&a, false));
        if !grid.contains(*min) || !grid.contains(*max) || !acc_ok {
            return Err("gate outside the editor's range".into());
        }
    }
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
        antenna_radius,
        rf_omega,
        radiation_reaction: _,
        // Checkboxes.
        beam_interaction: _,
        beam_retarded: _,
        instant_drain: _,
        t_max,
        tolerances: TolerancesSpec { preview, verify },
    } = physics;
    if c.is_some_and(|c| !positive(c))
        || !radius(*charge_radius)
        || !radius(*magnet_radius)
        || !radius(*wire_radius)
        || !radius(*antenna_radius)
        || !(0.0..=MAX_POSITIVE).contains(rf_omega)
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
                    moment,
                },
            launch:
                Launch {
                    node,
                    direction,
                    kinetic_energy,
                    time,
                },
            detector:
                Detector {
                    min,
                    max,
                    acceptance,
                },
            beam,
        } = s;
        let beam_ok = beam.is_none_or(|b| {
            let BeamSpec {
                count,
                energy_spread,
                angle_spread_deg,
                width,
                length,
                distribution: _,
                transmission,
                seed: _,
            } = b;
            (1..=MAX_BEAM).contains(&count)
                && (0.0..=0.3).contains(&energy_spread)
                && (0.0..=45.0).contains(&angle_spread_deg)
                && (0.0..=5.0).contains(&width)
                && (0.0..=5.0).contains(&length)
                && (0.01..=1.0).contains(&transmission)
        });
        if !beam_ok {
            return Err(format!("shot {}: beam outside the editor's range", i + 1));
        }
        if !charge.is_finite()
            || !moment.is_finite()
            || !positive(*mass)
            || !(0.0..=MAX_RADIUS).contains(r)
            || !on_grid(node)
            || !direction.iter().all(|d| d.is_finite())
            || direction[2] != 0.0
            || !positive(*kinetic_energy)
            || !time.is_finite()
            || !on_grid(min)
            || !on_grid(max)
            || !acceptance.is_none_or(|a| acceptance_in_range(&a, true))
        {
            return Err(format!(
                "shot {} has a value outside the editor's range",
                i + 1
            ));
        }
    }
    for (e, reference) in elements
        .iter()
        .map(|e| (e, false))
        .chain(reference_solution.iter().map(|e| (e, true)))
    {
        let Element {
            node,
            kind,
            value,
            angle_deg,
            omega,
            speed,
        } = e;
        match kind {
            ElementKind::Charge | ElementKind::Magnet | ElementKind::Antenna => {}
            // Player-only kinds (the level's own electrodes are `electrodes`, its dynamic
            // particles `free_particles`).
            ElementKind::Plate | ElementKind::Supply | ElementKind::Free if reference => {}
            ElementKind::Plate | ElementKind::Supply | ElementKind::Free => {
                return fail("plates, power supplies and free charges are player elements");
            }
        }
        if !on_grid(node)
            || !value.is_finite()
            || !angle_deg.is_finite()
            || omega.is_some_and(|w| !(0.0..=MAX_POSITIVE).contains(&w))
            || speed.is_some_and(|v| !v.is_finite() || v < 0.0)
        {
            return fail("element outside the grid or with an invalid value");
        }
    }
    for c in coils {
        match c {
            Coil::Circle {
                center,
                radius: r,
                kappa,
                rate,
                drive,
            } => {
                if !on_grid(center)
                    || !positive(*r)
                    || !kappa.is_finite()
                    || !rate.is_finite()
                    || !drive.as_ref().is_none_or(drive_editable)
                {
                    return fail("circular coil outside the editor's range");
                }
            }
            Coil::Polygon {
                vertices,
                kappa,
                rate,
            } => {
                if vertices.len() < 3
                    || !vertices.iter().all(on_grid)
                    || !kappa.is_finite()
                    || !rate.is_finite()
                {
                    return fail(
                        "polygon coil not reproducible (needs 3 or more vertices on the grid)",
                    );
                }
            }
        }
    }
    if electrodes.len() > MAX_COUNT as usize {
        return fail("too many electrodes");
    }
    for e in electrodes {
        let Electrode {
            center,
            length,
            thickness,
            height,
            angle_deg,
            bias,
            // A checkbox.
            tunable: _,
            drive,
        } = e;
        let value_ok = match bias {
            ConductorBias::Grounded => true,
            ConductorBias::Charge(v) | ConductorBias::Potential(v) => v.is_finite(),
        };
        let size_ok = [length, thickness, height]
            .iter()
            .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(*v));
        if !on_grid(center)
            || !size_ok
            || !angle_deg.is_finite()
            || !value_ok
            || !drive.as_ref().is_none_or(drive_editable)
        {
            return fail("electrode outside the editor's range");
        }
    }
    if conductors.len() > MAX_COUNT as usize {
        return fail("too many metal spheres");
    }
    for c in conductors {
        let Conductor {
            center,
            radius,
            bias,
        } = c;
        let value_ok = match bias {
            ConductorBias::Grounded => true,
            ConductorBias::Charge(v) | ConductorBias::Potential(v) => v.is_finite(),
        };
        if !on_grid(center) || !(MIN_POSITIVE..=MAX_POSITIVE).contains(radius) || !value_ok {
            return fail("metal sphere outside the editor's range");
        }
    }
    if free_particles.len() > MAX_COUNT as usize {
        return fail("too many free particles");
    }
    for f in free_particles {
        let FreeParticle {
            particle:
                ParticleSpec {
                    charge,
                    mass,
                    radius: r,
                    moment,
                },
            node,
            velocity,
            detector,
        } = f;
        let detector_ok = detector.is_none_or(|d| {
            let Detector {
                min,
                max,
                acceptance,
            } = d;
            on_grid(&min)
                && on_grid(&max)
                && acceptance.is_none_or(|a| acceptance_in_range(&a, false))
        });
        if !on_grid(node)
            || !charge.is_finite()
            || !(MIN_POSITIVE..=MAX_POSITIVE).contains(mass)
            || !(0.0..=MAX_RADIUS).contains(r)
            || !moment.is_finite()
            || !velocity.iter().all(|v| v.is_finite())
            || !detector_ok
        {
            return fail("free particle outside the editor's range");
        }
    }
    if clouds.len() > MAX_COUNT as usize {
        return fail("too many charge clouds");
    }
    for c in clouds {
        let Cloud {
            center,
            radius,
            charge,
        } = c;
        if !on_grid(center)
            || !(MIN_POSITIVE..=MAX_POSITIVE).contains(radius)
            || !charge.is_finite()
        {
            return fail("charge cloud outside the editor's range");
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
        max_antennas,
        antenna_amplitudes,
        antenna_omegas,
        // Any value is editable (a checkbox).
        continuous: _,
        region,
        max_plates,
        plate_voltages,
        plate,
        supply_voltages,
        max_free,
        free_charges,
        free_speeds,
        free_mass,
        free_radius,
    } = limits;
    if *max_free > MAX_COUNT
        || !free_charges.iter().all(|v| v.is_finite())
        || !free_speeds.iter().all(|v| v.is_finite() && *v >= 0.0)
        || !(MIN_POSITIVE..=MAX_POSITIVE).contains(free_mass)
        || !(0.0..=MAX_RADIUS).contains(free_radius)
    {
        return fail("free-charge limits outside the editor's range");
    }
    let plate_ok = [plate.length, plate.thickness, plate.height]
        .iter()
        .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(v));
    if !plate_ok
        || !plate_voltages
            .iter()
            .chain(supply_voltages)
            .all(|v| v.is_finite())
    {
        return fail("plate limits outside the editor's range");
    }
    if *max_charges > MAX_COUNT
        || *max_magnets > MAX_COUNT
        || *max_antennas > MAX_COUNT
        || *max_plates > MAX_COUNT
        || !magnitudes
            .iter()
            .chain(magnet_strengths)
            .chain(antenna_amplitudes)
            .chain(antenna_omegas)
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
                if !level::is_level_file(&p) {
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
