//! A tube level's specification: emitters, metals, circuits, goals, resolution.

use super::*;

#[allow(clippy::too_many_lines)]
pub(super) fn edit_tube(
    ui: &mut egui::Ui,
    tube: &mut Option<TubeSpec>,
    electrodes: usize,
    grid: &Grid,
) -> bool {
    let mut focus = false;
    let mut on = tube.is_some();
    if ui
        .checkbox(&mut on, "Tube level (z-invariant)")
        .on_hover_text("Needs two electrodes: a cathode and the goal's electrode")
        .changed()
    {
        *tube = on.then(|| TubeSpec {
            cathode: 0,
            charge_per_mass: -1.0,
            emit_toward_deg: None,
            saturation: None,
            ions: None,
            projectiles: Vec::new(),
            b_z: 0.0,
            solenoid: None,
            solenoid_drive: None,
            particle_rate: TL_RATE,
            particles: 1.0,
            emitter_mesh: 1.0,
            step: DEFAULT_STEP,
            goals: vec![CurrentGoal {
                electrode: 1,
                min: 0.0,
                max: 1.0,
                start: 10.0,
                end: 20.0,
                carriers: Carriers::Cathode,
            }],
            circuits: Vec::new(),
            potential_goals: Vec::new(),
            heat_goals: Vec::new(),
            output: None,
            metals: Vec::new(),
        });
    }
    let Some(TubeSpec {
        cathode,
        charge_per_mass,
        emit_toward_deg,
        saturation,
        ions,
        projectiles,
        b_z,
        solenoid,
        solenoid_drive,
        particle_rate,
        particles,
        emitter_mesh,
        step,
        goals,
        circuits,
        potential_goals,
        heat_goals,
        output,
        metals,
    }) = tube
    else {
        return focus;
    };
    let last = electrodes.max(1) - 1;
    // Electrodes are numbered from 1 in the editor, as in their section.
    let index = |ui: &mut egui::Ui, label: &str, i: &mut usize| {
        ui.horizontal(|ui| {
            ui.label(label);
            let mut k = *i + 1;
            let r = ui.add(egui::DragValue::new(&mut k).range(1..=last + 1));
            *i = (k - 1).min(last);
            r.has_focus()
        })
        .inner
    };
    let face = |ui: &mut egui::Ui, a: &mut Option<f64>, what: &str| {
        let mut focus = false;
        ui.horizontal(|ui| {
            let mut one_face = a.is_some();
            ui.checkbox(&mut one_face, "Coated face toward")
                .on_hover_text(format!(
                    "Only the {what}'s faces facing this direction (within 60°) emit"
                ));
            if one_face {
                let v = a.get_or_insert(0.0);
                focus |= ui
                    .add(
                        egui::DragValue::new(v)
                            .speed(1.0)
                            .range(-360.0..=360.0)
                            .suffix("°"),
                    )
                    .has_focus();
            } else {
                *a = None;
            }
        });
        focus
    };
    let ratio = |ui: &mut egui::Ui, label: &str, v: &mut f64, positive: bool| {
        ui.horizontal(|ui| {
            ui.label(label);
            let lo = if positive {
                MIN_POSITIVE
            } else {
                -MAX_POSITIVE
            };
            let r = ui
                .add(egui::DragValue::new(v).speed(0.01).range(lo..=MAX_POSITIVE))
                .has_focus();
            if *v == 0.0 {
                *v = if positive { 1.0 } else { -1.0 };
            }
            r
        })
        .inner
    };
    let free = |ui: &mut egui::Ui, v: &mut f64, speed: f64| {
        ui.add(
            egui::DragValue::new(v)
                .speed(speed)
                .range(-MAX_POSITIVE..=MAX_POSITIVE),
        )
        .has_focus()
    };
    focus |= index(ui, "Cathode: electrode", cathode);
    focus |= ratio(ui, "Carriers' q/m", charge_per_mass, false);
    focus |= face(ui, emit_toward_deg, "cathode");
    ui.horizontal(|ui| {
        let mut limited = saturation.is_some();
        ui.checkbox(&mut limited, "Temperature-limited, at most")
            .on_hover_text(
                "The cathode's emission per unit length of its surface cannot exceed this \
                 (Richardson's limit at its temperature), whatever the field",
            );
        if limited {
            let j = saturation.get_or_insert(0.1);
            focus |= ui
                .add(
                    egui::DragValue::new(j)
                        .speed(0.001)
                        .range(0.0..=MAX_POSITIVE),
                )
                .on_hover_text("0: a cold electrode, which emits nothing")
                .has_focus();
        } else {
            *saturation = None;
        }
    });
    if saturation.is_some() {
        ui.horizontal(|ui| {
            ui.label("  macroparticles per time").on_hover_text(
                "At the saturation, in the preview (the verification twice as many): a beam \
                 modulated at a high frequency needs many per period",
            );
            focus |= ui
                .add(
                    egui::DragValue::new(particle_rate)
                        .speed(1.0)
                        .range(1.0..=MAX_POSITIVE),
                )
                .has_focus();
        });
    }
    ui.horizontal(|ui| {
        ui.label("Macroparticles ×").on_hover_text(
            "Every emitter's macroparticles times this: a magnetized sheath needs many more \
             than a diode's flow (each preview costs about as many times more)",
        );
        focus |= ui
            .add(
                egui::DragValue::new(particles)
                    .speed(0.1)
                    .range(0.1..=100.0),
            )
            .has_focus();
        ui.label("emitter mesh ×").on_hover_text(
            "The emitters' segments this many times finer: the particles' release distance \
             (a first-order error) shrinks with them",
        );
        focus |= ui
            .add(
                egui::DragValue::new(emitter_mesh)
                    .speed(0.1)
                    .range(1.0..=16.0),
            )
            .has_focus();
    });
    // The radiated output (a magnetron's): one circuit's load.
    let mut has_output = output.is_some();
    ui.checkbox(&mut has_output, "Output (radiated)")
        .on_hover_text(
            "One circuit's load stands for an antenna along the tube's axis: what it absorbs is \
         radiated, drawn as a wave at the wave speed c (display only)",
        );
    if has_output {
        let o = output.get_or_insert(TubeOutput {
            circuit: 0,
            c: 50.0,
        });
        ui.horizontal(|ui| {
            ui.label("  circuit");
            let mut k = o.circuit + 1;
            focus |= ui
                .add(egui::DragValue::new(&mut k).range(1..=circuits.len().max(1)))
                .has_focus();
            o.circuit = k - 1;
            ui.label("wave speed c");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut o.c)
                        .speed(1.0)
                        .range(MIN_POSITIVE..=MAX_POSITIVE),
                )
                .has_focus();
        });
    } else {
        *output = None;
    }
    // Ion source.
    let mut has_ions = ions.is_some();
    ui.checkbox(&mut has_ions, "Ion source (bipolar)")
        .on_hover_text("An electrode emitting positive ions under the space-charge limit");
    if has_ions {
        let src = ions.get_or_insert(IonSource {
            electrode: last.min(1),
            charge_per_mass: 0.01,
            emit_toward_deg: None,
        });
        focus |= index(ui, "  Ions from electrode", &mut src.electrode);
        focus |= ratio(ui, "  Ions' q/m", &mut src.charge_per_mass, true);
        focus |= face(ui, &mut src.emit_toward_deg, "ion source");
    } else {
        *ions = None;
    }
    ui.horizontal(|ui| {
        ui.label("B along z");
        focus |= free(ui, b_z, 0.01);
        ui.label("step");
        focus |= positive(ui, step, 0.001, MAX_POSITIVE);
    });
    let mut has_solenoid = solenoid.is_some();
    ui.checkbox(&mut has_solenoid, "Player's solenoid (supply sets B)")
        .on_hover_text(
            "A solenoid around the tube the player powers: a supply on its terminal node \
             sets B along z, over the range of the limits' solenoid fields",
        );
    if has_solenoid {
        let m = grid.max_node();
        let n = solenoid.get_or_insert([m[0], m[1], 0]);
        ui.horizontal(|ui| {
            ui.label("  terminal");
            focus |= node(ui, n, grid);
        });
    } else {
        *solenoid = None;
    }
    let mut changing = solenoid_drive.is_some();
    ui.checkbox(&mut changing, "Field changes in time")
        .on_hover_text(
            "B along z follows this source's shape, its level set by the solenoid's supply \
             (else by B along z): a DC level, a sine's amplitude or a pulse's high. Its \
             change induces an electric field circulating about the board's centre \
             (Faraday's law): eddy currents in the metals.",
        );
    if changing {
        let s = solenoid_drive.get_or_insert(Source::Sine {
            offset: 0.0,
            amplitude: 1.0,
            omega: 0.3,
            phase: 0.0,
        });
        focus |= edit_source(ui, ("solenoid_drive", 0), s);
    } else {
        *solenoid_drive = None;
    }
    // Metals (the Drude model): resistors, Hall probes.
    let mut remove_metal = None;
    for (k, m) in metals.iter_mut().enumerate() {
        ui.push_id(("metal", k), |ui| {
            let MetalSpec {
                corners,
                contacts,
                density,
                charge_per_mass,
                thermal_speed,
                mean_free_path,
                core_radius,
                carriers_per_area,
                seed,
                label,
            } = m;
            ui.horizontal(|ui| {
                ui.label(format!("Metal {}", k + 1)).on_hover_text(
                    "A conductor's interior: carriers over the ions' background, scattering \
                     off ion cores placed at random (PHYSICS.md §2.11, the Drude metal)",
                );
                if ui.small_button("×").clicked() {
                    remove_metal = Some(k);
                }
            });
            ui.horizontal(|ui| {
                ui.label("  corners");
                focus |= node(ui, &mut corners[0], grid);
                focus |= node(ui, &mut corners[1], grid);
            });
            let mut remove_contact = None;
            ui.horizontal(|ui| {
                ui.label("  contacts").on_hover_text(
                    "Box electrodes along its sides (each face within a side): they take in \
                     the carriers that reach them and give off carriers into the metal",
                );
                for (i, c) in contacts.iter_mut().enumerate() {
                    let mut e = *c + 1;
                    focus |= ui
                        .add(egui::DragValue::new(&mut e).range(1..=last + 1))
                        .has_focus();
                    *c = (e - 1).min(last);
                    if ui.small_button("×").clicked() {
                        remove_contact = Some(i);
                    }
                }
                if contacts.len() < MAX_COUNT as usize && ui.small_button("+").clicked() {
                    contacts.push(0);
                }
            });
            if let Some(i) = remove_contact {
                contacts.remove(i);
            }
            ui.horizontal(|ui| {
                ui.label("  charge density")
                    .on_hover_text("The carriers' charge per area (the ions' is opposite)");
                focus |= positive(ui, density, 0.01, MAX_POSITIVE);
            });
            focus |= ratio(ui, "  carriers' q/m", charge_per_mass, false);
            ui.horizontal(|ui| {
                ui.label("  thermal speed")
                    .on_hover_text("√(kT/m) of the carriers and the lattice");
                focus |= positive(ui, thermal_speed, 0.01, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                ui.label("  mean free path").on_hover_text(
                    "ℓ = 1/(2 a n): with the cores' radius it sets their number per area \
                     (a/2ℓ at most 0.005: dilute)",
                );
                focus |= positive(ui, mean_free_path, 0.01, MAX_POSITIVE);
                ui.label("core radius");
                focus |= positive(ui, core_radius, 0.0001, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                ui.label("  carriers per area").on_hover_text(
                    "Macroparticles per unit area in the preview (the verification's twice as \
                     many)",
                );
                focus |= positive(ui, carriers_per_area, 1.0, MAX_POSITIVE);
                ui.label("seed");
                focus |= ui.add(egui::DragValue::new(seed)).has_focus();
            });
            focus |= label_editor(ui, label, "channel, load, Hall plate…");
        });
    }
    if let Some(k) = remove_metal {
        metals.remove(k);
    }
    if metals.len() < MAX_COUNT as usize && ui.button("+ Metal").clicked() {
        let m = grid.max_node();
        metals.push(MetalSpec {
            corners: [[m[0] / 4, m[1] / 3, 0], [3 * m[0] / 4, 2 * m[1] / 3, 0]],
            contacts: Vec::new(),
            density: 1.0,
            charge_per_mass: -1.0,
            thermal_speed: 1.0,
            mean_free_path: 0.5,
            core_radius: 0.002,
            carriers_per_area: 100.0,
            seed: 0,
            label: String::new(),
        });
    }
    // Current goals: each counts the cathode's carriers, the ions, or a lead's net.
    let mut remove_goal = None;
    for (k, g) in goals.iter_mut().enumerate() {
        ui.push_id(("current goal", k), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Current goal {}", k + 1));
                if ui.small_button("×").clicked() {
                    remove_goal = Some(k);
                }
                egui::ComboBox::from_id_salt("carriers")
                    .selected_text(match g.carriers {
                        Carriers::Cathode => "cathode's carriers",
                        Carriers::Ions => "ions",
                        Carriers::Net => "net (a lead's)",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut g.carriers,
                            Carriers::Cathode,
                            "cathode's carriers",
                        );
                        ui.selectable_value(&mut g.carriers, Carriers::Ions, "ions")
                            .on_hover_text("The ion source's ions (needs an ion source)");
                        ui.selectable_value(&mut g.carriers, Carriers::Net, "net (a lead's)")
                            .on_hover_text(
                                "Every carrier the electrode takes in less what it gives off: \
                                 the current in a metal's contact's lead",
                            );
                    });
            });
            focus |= index(ui, "  Current into electrode", &mut g.electrode);
            ui.horizontal(|ui| {
                ui.label("  from");
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut g.min)
                            .speed(0.01)
                            .range(-MAX_POSITIVE..=MAX_POSITIVE),
                    )
                    .has_focus();
                ui.label("to");
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut g.max)
                            .speed(0.01)
                            .range(-MAX_POSITIVE..=MAX_POSITIVE),
                    )
                    .has_focus();
            });
            ui.horizontal(|ui| {
                ui.label("  averaged over t");
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut g.start)
                            .speed(0.1)
                            .range(0.0..=MAX_POSITIVE),
                    )
                    .has_focus();
                ui.label("to");
                focus |= positive(ui, &mut g.end, 0.1, MAX_POSITIVE);
            });
        });
    }
    if let Some(k) = remove_goal {
        goals.remove(k);
    }
    if goals.len() < MAX_COUNT as usize && ui.button("Add a current goal").clicked() {
        goals.push(CurrentGoal {
            electrode: last.min(1),
            min: 0.0,
            max: 1.0,
            start: 10.0,
            end: 20.0,
            carriers: Carriers::Cathode,
        });
    }
    // Projectiles.
    let mut remove = None;
    for (j, p) in projectiles.iter_mut().enumerate() {
        let ProjectileSpec {
            node: start,
            direction_deg,
            speed,
            launch,
            charge,
            mass,
            detector,
        } = p;
        ui.push_id(("projectile", j), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Projectile {}", j + 1));
                focus |= node(ui, start, grid);
                if ui.small_button("×").clicked() {
                    remove = Some(j);
                }
            });
            ui.horizontal(|ui| {
                ui.label("q");
                focus |= free(ui, charge, 0.01);
                ui.label("m");
                focus |= positive(ui, mass, 0.01, MAX_POSITIVE);
            });
            ui.horizontal(|ui| {
                focus |= ui
                    .add(egui::DragValue::new(direction_deg).speed(1.0).suffix("°"))
                    .has_focus();
                ui.label("speed");
                focus |= ui
                    .add(
                        egui::DragValue::new(speed)
                            .speed(0.01)
                            .range(0.0..=MAX_POSITIVE),
                    )
                    .has_focus();
                ui.label("at t");
                focus |= ui
                    .add(
                        egui::DragValue::new(launch)
                            .speed(0.1)
                            .range(0.0..=MAX_POSITIVE),
                    )
                    .has_focus();
            });
            let mut has_detector = detector.is_some();
            ui.checkbox(&mut has_detector, "Detector (goal)");
            if has_detector {
                let m = grid.max_node();
                let d = detector.get_or_insert([[m[0] - 2, 0, 0], [m[0], 2, 0]]);
                focus |= node(ui, &mut d[0], grid);
                focus |= node(ui, &mut d[1], grid);
            } else {
                *detector = None;
            }
        });
    }
    if let Some(j) = remove {
        projectiles.remove(j);
    }
    if projectiles.len() < MAX_COUNT as usize && ui.small_button("+ projectile").clicked() {
        let m = grid.max_node();
        projectiles.push(ProjectileSpec {
            node: [0, m[1] / 2, 0],
            direction_deg: 0.0,
            speed: 1.0,
            launch: 0.0,
            charge: 0.1,
            mass: 1.0,
            detector: None,
        });
    }
    // Circuits of electrodes.
    let mut remove = None;
    for (k, c) in circuits.iter_mut().enumerate() {
        let TubeCircuit {
            electrode,
            source,
            resistance,
            capacitance,
            supply_offset,
            linked,
            inductance,
            return_to,
        } = c;
        ui.push_id(("tube_circuit", k), |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Circuit {}", k + 1));
                if ui.small_button("×").clicked() {
                    remove = Some(k);
                }
            });
            focus |= index(ui, "  on electrode", electrode);
            focus |= edit_source(ui, ("tube_circuit_source", k), source);
            ui.horizontal(|ui| {
                let mut floating = resistance.is_none();
                ui.checkbox(&mut floating, "floating")
                    .on_hover_text("No source: the electrode keeps its charge and its capacitor's");
                if floating {
                    *resistance = None;
                } else {
                    let r = resistance.get_or_insert(1.0);
                    ui.label("R")
                        .on_hover_text("0: an ideal source the electrode follows");
                    focus |= ui
                        .add(
                            egui::DragValue::new(r)
                                .speed(0.01)
                                .range(0.0..=MAX_POSITIVE),
                        )
                        .has_focus();
                }
                ui.label("C to ground");
                focus |= ui
                    .add(
                        egui::DragValue::new(capacitance)
                            .speed(0.01)
                            .range(0.0..=MAX_POSITIVE),
                    )
                    .has_focus();
            });
            // A coil in the branch, and where the branch returns (not for a floating
            // electrode, whose branch carries nothing).
            if resistance.is_some() {
                ui.horizontal(|ui| {
                    ui.label("L").on_hover_text(
                        "A coil in series with R (a tank's): 0 none. Integrated so that an LC \
                         tank keeps its energy",
                    );
                    focus |= ui
                        .add(
                            egui::DragValue::new(inductance)
                                .speed(0.01)
                                .range(0.0..=MAX_POSITIVE),
                        )
                        .has_focus();
                    let mut back = return_to.is_some();
                    ui.checkbox(&mut back, "returns to another electrode")
                        .on_hover_text(
                            "The branch ends on another electrode instead of ground (a tank \
                             between two electrodes)",
                        );
                    if !back {
                        *return_to = None;
                    } else if return_to.is_none() {
                        *return_to = Some(if *electrode == 0 { 1.min(last) } else { 0 });
                    }
                });
                if let Some(r) = return_to {
                    focus |= index(ui, "  returns to electrode", r);
                }
            } else {
                *inductance = 0.0;
                *return_to = None;
            }
            if matches!(source, Source::Sine { .. }) {
                ui.checkbox(supply_offset, "the supply sets the offset (a bias)");
            } else {
                *supply_offset = false;
            }
            // Further electrodes wired to this one (a grid's wires on an ideal source, a
            // magnetron's strapped vanes on a tank), 1-based.
            ui.horizontal(|ui| {
                    ui.label("linked electrodes").on_hover_text(
                        "Wired to this electrode: on an ideal source they follow it (a grid's                          wires); on a tank or a resistance they are one node with it (a                          magnetron's strapped vanes)",
                    );
                    let mut text = linked
                        .iter()
                        .map(|e| (e + 1).to_string())
                        .collect::<Vec<_>>()
                        .join(", ");
                    let r = ui.add(egui::TextEdit::singleline(&mut text).desired_width(LIST_WIDTH));
                    if r.changed() {
                        *linked = text
                            .split(',')
                            .filter_map(|t| t.trim().parse::<usize>().ok())
                            .filter(|&e| (1..=last + 1).contains(&e))
                            .map(|e| e - 1)
                            .collect();
                    }
                    focus |= r.has_focus();
                });
        });
    }
    if let Some(k) = remove {
        circuits.remove(k);
    }
    if circuits.len() < MAX_COUNT as usize && ui.small_button("+ circuit").clicked() {
        circuits.push(TubeCircuit {
            electrode: last.min(1),
            source: Source::Dc { value: 0.0 },
            resistance: Some(1.0),
            capacitance: 1.0,
            supply_offset: false,
            linked: Vec::new(),
            inductance: 0.0,
            return_to: None,
        });
    }
    // Potential goals.
    let mut remove = None;
    for (k, g) in potential_goals.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("Potential goal {}", k + 1));
            if ui.small_button("×").clicked() {
                remove = Some(k);
            }
        });
        focus |= index(ui, "  electrode", &mut g.electrode);
        let mut relative = g.reference.is_some();
        ui.checkbox(&mut relative, "  against another electrode")
            .on_hover_text("A voltmeter across the two (a Hall probe's voltage)");
        if relative {
            let r = g.reference.get_or_insert(0);
            focus |= index(ui, "    reference", r);
        } else {
            g.reference = None;
        }
        ui.horizontal(|ui| {
            ui.selectable_value(&mut g.measure, Measure::Mean, "mean");
            ui.selectable_value(&mut g.measure, Measure::Swing, "swing")
                .on_hover_text("Half the peak-to-peak swing over the window");
            ui.selectable_value(&mut g.measure, Measure::Amplitude, "amplitude")
                .on_hover_text("The oscillation's amplitude at ω (a drift removed)");
            ui.selectable_value(&mut g.measure, Measure::Phase, "phase")
                .on_hover_text(
                    "The oscillation's phase at ω relative to sin(ωt), in degrees (modulo 360)",
                );
        });
        ui.horizontal(|ui| {
            if matches!(g.measure, Measure::Amplitude | Measure::Phase) {
                ui.label("  ω");
                focus |= ui
                    .add(
                        egui::DragValue::new(&mut g.omega)
                            .speed(0.01)
                            .range(0.0..=MAX_POSITIVE),
                    )
                    .has_focus();
            }
            ui.label("  from");
            focus |= free(ui, &mut g.min, 0.01);
            ui.label("to");
            focus |= free(ui, &mut g.max, 0.01);
        });
        ui.horizontal(|ui| {
            ui.label("  over t");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut g.start)
                        .speed(0.1)
                        .range(0.0..=MAX_POSITIVE),
                )
                .has_focus();
            ui.label("to");
            focus |= positive(ui, &mut g.end, 0.1, MAX_POSITIVE);
        });
    }
    if let Some(k) = remove {
        potential_goals.remove(k);
    }
    if potential_goals.len() < MAX_COUNT as usize && ui.small_button("+ potential goal").clicked() {
        potential_goals.push(PotentialGoal {
            electrode: last,
            reference: None,
            measure: Measure::Mean,
            omega: 0.0,
            min: 0.0,
            max: 1.0,
            start: 20.0,
            end: 40.0,
        });
    }
    // Heat goals: a metal's mean power over a window.
    let mut remove = None;
    let n_metals = metals.len();
    for (k, g) in heat_goals.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("Heat goal {}", k + 1)).on_hover_text(
                "A metal's heat: the mean power the field gives its carriers over the window \
                 (the Joule heat)",
            );
            if ui.small_button("×").clicked() {
                remove = Some(k);
            }
        });
        ui.horizontal(|ui| {
            ui.label("  metal");
            let mut m = g.metal + 1;
            focus |= ui
                .add(egui::DragValue::new(&mut m).range(1..=n_metals.max(1)))
                .has_focus();
            g.metal = m - 1;
        });
        ui.horizontal(|ui| {
            ui.label("  from");
            focus |= free(ui, &mut g.min, 0.01);
            ui.label("to");
            focus |= free(ui, &mut g.max, 0.01);
        });
        ui.horizontal(|ui| {
            ui.label("  over t");
            focus |= ui
                .add(
                    egui::DragValue::new(&mut g.start)
                        .speed(0.1)
                        .range(0.0..=MAX_POSITIVE),
                )
                .has_focus();
            ui.label("to");
            focus |= positive(ui, &mut g.end, 0.1, MAX_POSITIVE);
        });
    }
    if let Some(k) = remove {
        heat_goals.remove(k);
    }
    if n_metals > 0
        && heat_goals.len() < MAX_COUNT as usize
        && ui.small_button("+ heat goal").clicked()
    {
        heat_goals.push(HeatGoal {
            metal: 0,
            min: 0.0,
            max: 1.0,
            start: 20.0,
            end: 40.0,
        });
    }
    focus
}
