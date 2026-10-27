//! The field maps and the field lines, recomputed as the placement or the view
//! changes.

use super::*;

/// Updates the GPU field map for the active shot and setup.
pub(super) fn update_map(
    mut game: ResMut<Game>,
    quad: Res<PotentialQuad>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut buffers: ResMut<Assets<bevy::render::storage::ShaderBuffer>>,
    mut transforms: Query<&mut Transform>,
) {
    // A map the level does not offer (after a level switch or an edit) gives way to the
    // potential map: it stayed shown with no choice in the panel to switch it off (a
    // particle field's radiation part over a tube level).
    if !map_available(&game.editor.level, game.map) {
        game.map = Some(MapMode::Potential);
    }
    // The energy unit depends on the flight (particles launched at rest): redraw when it
    // arrives.
    let has_flight = if game.editor.level.is_tube() {
        // Tube levels: the map follows the preview's frames.
        game.tube_preview
            .as_ref()
            .is_some_and(|v| v.revision == game.sent_revision)
    } else {
        game.flights
            .get(game.active_flight())
            .is_some_and(|f| f.preview_revision == game.sent_revision && f.preview.is_some())
    };
    let disturbance = game
        .active_disturbance
        .min(game.editor.level.flights_per_shot() - 1);
    let key = (
        game.sent_revision,
        game.active_shot,
        disturbance,
        game.show_all_shots,
        game.map,
        has_flight,
        // Tube levels follow the animation time (interpolated pictures between the frames)
        // and the frames streamed up to it: a time held ahead of the run's stream (a paused
        // animation, a capture's EM_TIME) is drawn again as the stream reaches it (it kept
        // the map of the first frames that had arrived).
        if game.editor.level.is_tube() {
            let streamed = game.tube_preview.as_ref().map_or(0, |v| {
                let f = &v.run.frames;
                (f.partition_point(|f| f.t <= game.anim_time) + 1).min(f.len())
            });
            let picture = if game.animate {
                game.anim_time.to_bits()
            } else {
                u64::try_from(game.tube_frame().unwrap_or(0)).unwrap_or(u64::MAX)
            };
            picture ^ u64::try_from(streamed).unwrap_or(0).rotate_left(29)
        } else {
            0
        },
        (
            game.tube_gain_decades.to_bits() ^ game.map_tau.to_bits().rotate_left(7),
            game.tube_range_decades.to_bits()
                ^ (game.tube_particles as u64)
                ^ game
                    .em_snapshot
                    .as_ref()
                    .map_or(0, |s| s.t.to_bits().rotate_left(17)),
        ),
    );
    if key == game.map_key {
        return;
    }
    game.map_key = key;
    let Some(
        mode @ (MapMode::Potential
        | MapMode::Magnetic
        | MapMode::Electric
        | MapMode::EnergyFlow
        | MapMode::Current
        | MapMode::Heat),
    ) = game.map
    else {
        return;
    };
    // A tube level's maps averaged over time (the player's meter), before the level is
    // borrowed below.
    let (averaged, em_averaged) = if game.editor.level.is_tube() && game.map_tau > 0.0 {
        let (tau, revision) = (game.map_tau, game.sent_revision);
        let picture = game.tube_display_frame();
        let snap = game.em_snapshot.clone();
        let g = &mut *game;
        (
            picture.map(|p| g.map_average.frame(revision, &p, tau)),
            snap.map(|s| g.map_average.snapshot(&s, tau)),
        )
    } else {
        (None, None)
    };
    let level = &game.editor.level;
    let bounds = level.bounds();
    if let Ok(mut t) = transforms.get_mut(quad.entity) {
        let size = bounds.max - bounds.min;
        let center = (bounds.max + bounds.min) * 0.5;
        #[allow(clippy::cast_possible_truncation)]
        {
            t.scale = Vec3::new(size.x as f32, size.y as f32, 1.0);
            t.translation = Vec3::new(center.x as f32, center.y as f32, -10.0);
        }
    }
    // Tube levels: the frame shown (live while animating); nothing (a blank map) until
    // the preview arrives, and a level without shots shows a blank map (not the
    // previous level's).
    if level.is_tube() || level.shots.is_empty() {
        let mut new_units = None;
        let mut new_metal_map = None;
        let shown = game.tube_frame();
        let picture = game.tube_display_frame();
        let field = match (&level.tube, &game.tube_preview, shown, &picture) {
            (Some(spec), Some(v), Some(k), Some(picture)) if v.revision == game.sent_revision => 'arm: {
                let run = &v.run;
                // Many particles merged for the map (the shader and the arrows sum them
                // all); averaged over time, the meter's cells.
                let binned = match &averaged {
                    Some(a) => std::borrow::Cow::Borrowed(a),
                    None => potential::map_frame(picture),
                };
                let picture = averaged.as_ref().unwrap_or(picture);
                // The metals' maps (`level::tube::TubeRun::metal_rates`): the current
                // density or the field's power per area inside them, by cell, over the span
                // before the time shown.
                if matches!(mode, MapMode::Current | MapMode::Heat) {
                    let (span, t) = (metal_map_span(game.map_tau), picture.t);
                    let Some(maps) = run.metal_maps(t, span) else {
                        break 'arm None;
                    };
                    let rates = faded(&maps);
                    // The colour units, the largest so far on this level (the scale holds
                    // still as the run plays, and setups compare on it).
                    let now = (
                        metal_percentile(&rates, &|r, c| r.j[c][0].hypot(r.j[c][1])),
                        metal_percentile(&rates, &|r, c| r.power[c].abs()),
                    );
                    let units = match game
                        .metal_map
                        .as_ref()
                        .filter(|m| m.level == game.level_index)
                    {
                        Some(m) => (m.units.0.max(now.0), m.units.1.max(now.1)),
                        None => now,
                    };
                    let heat = mode == MapMode::Heat;
                    let scale = if heat { units.1 } else { units.0 }.max(1e-300);
                    let (snap, arrows) = metal_map_grid(
                        run,
                        &rates,
                        (bounds.min, bounds.max),
                        heat,
                        (v.revision, t),
                    );
                    let (params, items) = potential::tube_field_params(
                        run,
                        &binned,
                        picture,
                        spec.cathode,
                        &snap,
                        heat,
                        scale,
                        game.tube_range_decades,
                        game.tube_particles == TubeParticles::Density,
                    );
                    new_metal_map = Some(MetalMap {
                        revision: v.revision,
                        level: game.level_index,
                        units,
                        maps,
                    });
                    break 'arm Some((
                        TubeField {
                            revision: v.revision,
                            frame: k,
                            scale,
                            scale_final: false,
                            arrows,
                        },
                        params,
                        items,
                    ));
                }
                // The electromagnetic field's maps (`em.rs`).
                if matches!(mode, MapMode::Magnetic | MapMode::EnergyFlow) {
                    let Some(snap) = em_averaged.as_ref().or(game.em_snapshot.as_ref()) else {
                        break 'arm None;
                    };
                    // The colour unit: the 99th percentile over the board; for a tube that
                    // radiates (an output), over the board outside its rings: the radiation.
                    // (Over the whole board the tube's inside, a few per cent of it with the
                    // strongest fields, set the unit, and the radiation leaving the tube lay
                    // decades below: black.) The inside then saturates.
                    let max = level.grid.position(level.grid.max_node());
                    let rings: Vec<(physics::DVec3, f64)> = level
                        .electrodes
                        .iter()
                        .filter(|e| e.shape == level::ElectrodeShape::Ring)
                        .map(|e| {
                            (
                                level.grid.position(e.center),
                                0.5 * (e.length + e.thickness),
                            )
                        })
                        .collect();
                    let radiating = spec.output.is_some();
                    #[allow(clippy::cast_precision_loss)]
                    let counted = |n: usize| {
                        let x = snap.origin.0 + ((n % snap.nx) as f64 + 0.5) * snap.spacing;
                        let y = snap.origin.1 + ((n / snap.nx) as f64 + 0.5) * snap.spacing;
                        !radiating
                            || ((0.0..=max.x).contains(&x)
                                && (0.0..=max.y).contains(&y)
                                && rings.iter().all(|(c, r)| (x - c.x).hypot(y - c.y) > *r))
                    };
                    let peak = |f: &dyn Fn(&[f32; 3]) -> f64| {
                        let mut v: Vec<f64> = snap
                            .samples
                            .iter()
                            .enumerate()
                            .filter(|(n, _)| counted(*n))
                            .map(|(_, s)| f(s))
                            .collect();
                        v.sort_by(f64::total_cmp);
                        #[allow(
                            clippy::cast_possible_truncation,
                            clippy::cast_sign_loss,
                            clippy::cast_precision_loss
                        )]
                        v.get((0.99 * v.len() as f64) as usize)
                            .copied()
                            .unwrap_or(0.0)
                    };
                    let b = peak(&|s| f64::from(s[0].abs()));
                    let sm = peak(&|s| f64::from(s[1].hypot(s[2])));
                    let units = if game.em_units.0 == v.revision {
                        (v.revision, game.em_units.1.max(b), game.em_units.2.max(sm))
                    } else {
                        (v.revision, b, sm)
                    };
                    new_units = Some(units);
                    let magnetic = mode == MapMode::Magnetic;
                    let scale = if magnetic { units.1 } else { units.2 }.max(1e-300);
                    let arrows = if magnetic {
                        Vec::new()
                    } else {
                        tube_flow_arrows(level, snap)
                    };
                    let (params, items) = potential::tube_field_params(
                        run,
                        &binned,
                        picture,
                        spec.cathode,
                        snap,
                        magnetic,
                        scale,
                        game.tube_range_decades,
                        game.tube_particles == TubeParticles::Density,
                    );
                    break 'arm Some((
                        TubeField {
                            revision: v.revision,
                            frame: k,
                            scale,
                            scale_final: false,
                            arrows,
                        },
                        params,
                        items,
                    ));
                }
                let old = game
                    .tube_field
                    .as_ref()
                    .filter(|f| f.revision == v.revision);
                // The electric map's field, once per frame (summing every segment and
                // particle at every arrow: on a large board twice per frame was a
                // slideshow); its colour unit is the field at the verification window's end
                // once the run has passed it, the frame's own until then.
                let (now, arrows) = if mode == MapMode::Electric {
                    let (s, a) = tube_field(level, run, &binned);
                    (Some(s), a)
                } else {
                    (None, Vec::new())
                };
                let (scale, scale_final) = match old {
                    Some(f) if f.scale_final => (f.scale, true),
                    _ => match run
                        .frames
                        .last()
                        .filter(|f| f.t >= level.tube_verify_end())
                        .and_then(|_| run.frame_at(level.tube_verify_end()))
                    {
                        Some(f) if mode == MapMode::Electric => {
                            (tube_field(level, run, &potential::map_frame(f)).0, true)
                        }
                        _ => (now.or(old.map(|f| f.scale)).unwrap_or(1.0), false),
                    },
                };
                let (params, items) = potential::tube_params(
                    run,
                    &binned,
                    picture,
                    spec.cathode,
                    spec.charge_per_mass.signum(),
                    mode == MapMode::Electric,
                    scale,
                    10f64.powf(game.tube_gain_decades),
                    game.tube_range_decades,
                    game.tube_particles == TubeParticles::Density,
                );
                Some((
                    TubeField {
                        revision: v.revision,
                        frame: k,
                        scale,
                        scale_final,
                        arrows,
                    },
                    params,
                    items,
                ))
            }
            _ => None,
        };
        if let Some(u) = new_units {
            game.em_units = u;
        }
        if new_metal_map.is_some() {
            game.metal_map = new_metal_map;
        }
        let (params, items) = match field {
            Some((f, params, items)) => {
                game.tube_field = Some(f);
                (params, items)
            }
            None => (potential::PotentialParams::default(), Vec::new()),
        };
        if let Some(mut m) = materials.get_mut(&quad.material) {
            m.params = params;
            if let Some(mut b) = buffers.get_mut(&m.items) {
                *b = potential::buffer(items);
            }
        }
        return;
    }
    // The active flight's field: its shot under its disturbance.
    let scenario = level.display_scenario(game.active_shot, disturbance, &game.editor.placement);
    // The dark region: forbidden for every particle shown (the whole beam of a beam
    // shot, every shot in the collective view).
    let mut limits = vec![Vec::new(); level.shots.len()];
    for l in level.launches() {
        limits[l.shot].push(l);
    }
    if !game.show_all_shots {
        limits = vec![std::mem::take(&mut limits[game.active_shot])];
    }
    // Time-dependent fields (a ramped coil's induced field, antennas, waves) do work:
    // energy conservation forbids nothing.
    if !physics::field::FieldSolver::is_static(&scenario.field) {
        limits.clear();
    }
    let unit = game.energy_unit(game.active_shot).0;
    let (params, items) = potential::params(
        &scenario,
        &limits,
        (
            level.physics.charge_radius,
            level.physics.magnet_radius,
            level.physics.antenna_radius,
        ),
        mode,
        unit,
    );
    if let Some(mut m) = materials.get_mut(&quad.material) {
        m.params = params;
        if let Some(mut b) = buffers.get_mut(&m.items) {
            *b = potential::buffer(items);
        }
    }
}

/// Recomputes the field lines (the electric field does not depend on the shot, but on the
/// disturbance) when they are shown and the setup, spacing or disturbance changed.
pub(super) fn update_field_lines(mut game: ResMut<Game>) {
    use std::sync::atomic::Ordering;
    // Results of a finished job for the current key.
    let disturbance = game
        .active_disturbance
        .min(game.editor.level.flights_per_shot() - 1);
    let key = (
        game.sent_revision,
        game.field_line_spacing.to_bits(),
        disturbance,
    );
    let arrived = game.field_lines_job.as_ref().and_then(|(k, rx)| {
        (*k == key)
            .then(|| rx.lock().ok().and_then(|r| r.try_recv().ok()))
            .flatten()
    });
    if let Some(lines) = arrived {
        game.field_lines = lines;
        game.field_lines_job = None;
    }
    if game.shot_count() == 0 {
        // No shots (a tube level): the 3D slice's field lines do not apply; drop any left
        // from the previous level.
        game.field_lines.clear();
        return;
    }
    if !game.show_field_lines || key == game.field_lines_key {
        return;
    }
    game.field_lines_key = key;
    // Computed on a thread of its own (near metal it takes seconds); a newer key stops it.
    // The electrostatic field: the static part (no oscillating sources), without a
    // ramped coil's induced field (it has no potential; its lines close on themselves).
    let mut scenario = game
        .editor
        .level
        .display_scenario(0, disturbance, &game.editor.placement);
    scenario.field = potential::static_part(&scenario.field);
    for l in &mut scenario.field.loops {
        l.rate = 0.0;
    }
    for p in &mut scenario.field.polygons {
        p.rate = 0.0;
    }
    let radius = game.editor.level.physics.charge_radius;
    let spacing = game.field_line_spacing;
    let token = key.0 ^ key.1.rotate_left(17) ^ (key.2 as u64).rotate_left(41);
    game.field_lines_wanted.store(token, Ordering::Release);
    let wanted = game.field_lines_wanted.clone();
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let keep_going = || wanted.load(Ordering::Acquire) == token;
        if let Some(lines) =
            visuals::field_lines_cancellable(&scenario, radius, spacing, keep_going)
        {
            let drawn: Vec<DrawnFieldLine> = lines
                .into_iter()
                .map(|l| {
                    (
                        l.points.into_iter().map(draw::to_vec2).collect(),
                        l.arrows
                            .into_iter()
                            .map(|(p, d)| (draw::to_vec2(p), draw::to_vec2(d)))
                            .collect(),
                    )
                })
                .collect();
            let _ = tx.send(drawn);
        }
    });
    game.field_lines_job = Some((key, std::sync::Mutex::new(rx)));
}
