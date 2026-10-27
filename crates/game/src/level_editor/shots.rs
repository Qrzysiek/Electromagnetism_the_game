//! The grid, the physics, the shots with their detectors' acceptance and radiation
//! goals, the free particles, the gates and the beams.

use super::*;

pub(super) fn edit_grid(ui: &mut egui::Ui, g: &mut Grid, refine_by: &mut Option<u32>) -> bool {
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

pub(super) fn edit_physics(ui: &mut egui::Ui, p: &mut WorldPhysics) -> bool {
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

pub(super) fn edit_shot(ui: &mut egui::Ui, s: &mut Shot, grid: &Grid) -> bool {
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
pub(super) fn edit_acceptance(
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
pub(super) fn edit_radiation(ui: &mut egui::Ui, radiation: &mut Option<RadiationGoal>) -> bool {
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
pub(super) fn acceptance_in_range(a: &DetectorAcceptance, radiation_allowed: bool) -> bool {
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

pub(super) fn edit_free_particles(
    ui: &mut egui::Ui,
    list: &mut Vec<FreeParticle>,
    grid: &Grid,
) -> bool {
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
pub(super) fn edit_gates(ui: &mut egui::Ui, list: &mut Vec<Detector>, grid: &Grid) -> bool {
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

pub(super) fn edit_beam(ui: &mut egui::Ui, beam: &mut Option<BeamSpec>) -> bool {
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
