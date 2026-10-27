//! A tube or metal level on the board: the flow's arrows, the electromagnetic
//! field, the drift, the traced carriers, the Hall forces, the animation.

use super::*;

/// The energy flow map's arrows: the Poynting vector of the field's snapshot at the
/// arrow grid's points (the nearest sample).
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(super) fn tube_flow_arrows(level: &Level, snap: &em::EmSnapshot) -> TubeArrows {
    let m = level.grid.max_node();
    let mut arrows = Vec::new();
    if snap.nx == 0 || snap.ny == 0 {
        return arrows;
    }
    // About 24 arrows across the board.
    let spacing = usize::try_from(m[0].max(m[1]) / 24)
        .unwrap_or(0)
        .max(TUBE_ARROW_SPACING);
    for y in (0..=m[1]).step_by(spacing) {
        for x in (0..=m[0]).step_by(spacing) {
            let p = level.grid.position([x, y, 0]);
            let i = (((p.x - snap.origin.0) / snap.spacing) as usize).min(snap.nx - 1);
            let j = (((p.y - snap.origin.1) / snap.spacing) as usize).min(snap.ny - 1);
            let s = snap.samples[i + snap.nx * j];
            arrows.push((draw::to_vec2(p), Vec2::new(s[1], s[2])));
        }
    }
    arrows
}

/// Tube levels: the charge each lead carries along the preview run, reckoned up to the time
/// shown (`level::tube_em::LeadFlow`; a new run starts it over): the beads on the leads in
/// every view but the field maps, which take the field's own (`draw::draw_leads`).
fn update_lead_flow(game: &mut Game) {
    let Game {
        tube_preview,
        lead_flow,
        editor,
        animate,
        anim_time,
        ..
    } = game;
    let Some(v) = tube_preview.as_ref().filter(|v| v.run.frames.len() >= 2) else {
        *lead_flow = None;
        return;
    };
    if lead_flow.as_ref().is_none_or(|(r, _)| *r != v.revision) {
        *lead_flow = Some((
            v.revision,
            level::tube_em::LeadFlow::new(&editor.level, &v.run),
        ));
    }
    let t = if *animate {
        *anim_time
    } else {
        v.run.frames.last().map_or(0.0, |f| f.t)
    };
    if let Some((_, flow)) = lead_flow.as_mut() {
        flow.advance(&v.run, t);
    }
}

/// Tube levels' field maps: keeps the field's thread on the preview and the time shown,
/// and takes its snapshots.
pub(super) fn update_em(mut game: ResMut<Game>) {
    let game = &mut *game;
    update_lead_flow(game);
    update_swing(game);
    update_drift(game);
    update_traced(game);
    update_hall(game);
    let field_map = game.editor.level.is_tube()
        && matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow));
    // Once started for this run, the field follows the playback under the other maps too
    // (switching back shows it at once, rather than catching up from where it was left).
    let following =
        game.editor.level.is_tube() && game.em.is_some() && game.em_run == Some(game.sent_revision);
    if !field_map && !following {
        return;
    }
    let Some(v) = game
        .tube_preview
        .as_ref()
        .filter(|v| v.revision == game.sent_revision)
    else {
        return;
    };
    let w = game.em.get_or_insert_with(em::EmWorker::start);
    if game.em_run != Some(v.revision) {
        w.send(em::EmRequest::Reset {
            revision: v.revision,
            level: std::sync::Arc::new(game.editor.level.clone()),
            run: Box::new(v.run.clone()),
        });
        game.em_run = Some(v.revision);
        game.em_asked = -1.0;
        game.em_snapshot = None;
    }
    let t = if game.animate {
        game.anim_time
    } else {
        v.run.frames.last().map_or(0.0, |f| f.t)
    };
    if (t - game.em_asked).abs() > 1e-9 {
        w.send(em::EmRequest::Want {
            revision: v.revision,
            t,
        });
        game.em_asked = t;
    }
    if let Some(s) = w.latest()
        && Some(s.revision) == game.em_run
    {
        game.em_snapshot = Some(s);
    }
}

/// Whether a point lies in a metal of the run.
pub(super) fn in_metals(metals: &[level::tube::MetalView], x: [f32; 2]) -> Option<usize> {
    metals.iter().position(|m| {
        f64::from(x[0]) >= m.lo.x
            && f64::from(x[0]) <= m.hi.x
            && f64::from(x[1]) >= m.lo.y
            && f64::from(x[1]) <= m.hi.y
    })
}

/// Tube levels with metals: each metal's carriers' mean displacement, frame by frame from
/// the run's start (`Drift`, the sum of `level::tube::metal_steps`), extended as frames
/// arrive.
pub(super) fn update_drift(game: &mut Game) {
    let Some(v) = game
        .tube_preview
        .as_ref()
        .filter(|v| v.revision == game.sent_revision && !v.run.metals.is_empty())
    else {
        game.drift = None;
        return;
    };
    let (metals, frames) = (&v.run.metals, &v.run.frames);
    let d = game.drift.get_or_insert_with(Drift::default);
    if d.revision != v.revision || d.x.len() != metals.len() {
        *d = Drift {
            revision: v.revision,
            x: vec![Vec::new(); metals.len()],
        };
    }
    let done = d.x.first().map_or(0, Vec::len);
    for k in done..frames.len().min(done + DRIFT_FRAMES_PER_UPDATE) {
        if k == 0 {
            for x in &mut d.x {
                x.push([0.0; 2]);
            }
            continue;
        }
        let steps = level::tube::metal_steps(metals, &frames[k - 1], &frames[k]);
        for (x, step) in d.x.iter_mut().zip(steps) {
            let last = x.last().copied().unwrap_or([0.0; 2]);
            x.push([last[0] + step[0], last[1] + step[1]]);
        }
    }
}

/// Tube levels with metals: the carriers traced (`Traced`), `draw::METAL_SAMPLE` of them
/// spread evenly over the metals' carriers when the run or the playback starts, each one a
/// contact takes replaced by another carrier in the metals, drawn in a fixed order that has
/// nothing to do with place (a hash of identity): the sample stays spread over the metals.
/// (Replaced by the newest carrier, just given off next to a contact and soon taken again,
/// the traces gathered at the contacts: a fifth of them within half a unit.)
pub(super) fn update_traced(game: &mut Game) {
    let Some(revision) = game
        .tube_preview
        .as_ref()
        .filter(|v| v.revision == game.sent_revision && !v.run.metals.is_empty())
        .map(|v| v.revision)
    else {
        game.traced = None;
        return;
    };
    let Some(frame) = game.tube_display_frame() else {
        return;
    };
    let Some((metals, every)) = game
        .tube_preview
        .as_ref()
        .map(|v| (v.run.metals.clone(), v.run.fine_every))
    else {
        return;
    };
    // The carriers in the metals now, by identity: of those recorded at every step where
    // frames are steps apart (`level::tube::TubeFrame::paths`: their arcs show).
    let mut present: Vec<u32> = frame
        .id
        .iter()
        .zip(&frame.x)
        .zip(&frame.kind)
        .filter(|((id, x), kind)| {
            **kind != level::tube::KIND_PROJECTILE
                && (every == 0 || id.is_multiple_of(every))
                && in_metals(&metals, **x).is_some()
        })
        .map(|((id, _), _)| *id)
        .collect();
    present.sort_unstable();
    let restart = game
        .traced
        .as_ref()
        .is_none_or(|t| t.revision != revision || frame.t < t.t - 1e-9);
    if restart {
        let every = (present.len() / draw::METAL_SAMPLE).max(1);
        game.traced = Some(Traced {
            revision,
            t: frame.t,
            ids: present.iter().step_by(every).copied().collect(),
        });
        return;
    }
    let Some(traced) = game.traced.as_mut() else {
        return;
    };
    let alive: std::collections::HashSet<u32> = present.iter().copied().collect();
    let kept: std::collections::HashSet<u32> = traced
        .ids
        .iter()
        .copied()
        .filter(|id| alive.contains(id))
        .collect();
    if kept.len() < traced.ids.len() {
        let mut candidates: Vec<u32> = present
            .iter()
            .copied()
            .filter(|id| !kept.contains(id))
            .collect();
        candidates.sort_unstable_by_key(|id| scatter(*id));
        let mut candidates = candidates.into_iter();
        for id in &mut traced.ids {
            if !alive.contains(id)
                && let Some(n) = candidates.next()
            {
                *id = n;
            }
        }
    }
    traced.t = frame.t;
}

/// Tube levels with metals in a field: the forces on each strip's average carrier at the
/// time shown (`level::tube::TubeRun::hall_forces`, over the last `draw::DRIFT_SPAN`), for
/// the board's diagram and its labels. Only in the pieces of metals with contacts, a current
/// through them: in an isolated metal (eddy currents) the middle's mean drift is noise.
pub(super) fn update_hall(game: &mut Game) {
    let hall = (|| {
        let v = game
            .tube_preview
            .as_ref()
            .filter(|v| v.run.has_field() && !v.run.metals.is_empty())?;
        let level = &game.editor.level;
        let spec = level.tube.as_ref()?;
        // A piece (`run.metals`) belongs to the metal its centre lies in.
        let fed = |m: &level::tube::MetalView| {
            let c = (m.lo + m.hi) * 0.5;
            spec.metals.iter().any(|s| {
                let (lo, hi) = level.metal_rect(s);
                !s.contacts.is_empty() && c.x > lo.x && c.x < hi.x && c.y > lo.y && c.y < hi.y
            })
        };
        let frame = game.tube_display_frame()?;
        let before = v.run.frame_between(frame.t - draw::DRIFT_SPAN)?;
        if frame.t - before.t < 0.5 * draw::DRIFT_SPAN {
            return None;
        }
        Some(
            (0..v.run.metals.len())
                .filter(|&k| fed(&v.run.metals[k]))
                .filter_map(|k| v.run.hall_forces(&frame, &before, k))
                .collect(),
        )
    })();
    game.hall = hall.unwrap_or_default();
}

/// An order of identities unrelated to them (Fibonacci hashing: the golden ratio's 64-bit
/// multiple, its high bits).
pub(super) fn scatter(id: u32) -> u64 {
    u64::from(id)
        .wrapping_mul(0x9e37_79b9_7f4a_7c15)
        .rotate_left(32)
}

/// Tube levels: each electrode's potential's mean and swing over the preview run's frames
/// (`Swing`), recomputed when frames arrive.
#[allow(clippy::cast_precision_loss)]
pub(super) fn update_swing(game: &mut Game) {
    let Some(v) = game
        .tube_preview
        .as_ref()
        .filter(|v| v.revision == game.sent_revision)
    else {
        game.swing = None;
        return;
    };
    let frames = &v.run.frames;
    if game
        .swing
        .as_ref()
        .is_some_and(|s| s.revision == v.revision && s.frames == frames.len())
    {
        return;
    }
    let n = frames.iter().map(|f| f.potentials.len()).max().unwrap_or(0);
    let mut mean = vec![0.0; n];
    for f in frames {
        for (m, p) in mean.iter_mut().zip(&f.potentials) {
            *m += f64::from(*p);
        }
    }
    for m in &mut mean {
        *m /= frames.len().max(1) as f64;
    }
    let mut own = vec![0.0_f64; n];
    for f in frames {
        for ((o, m), p) in own.iter_mut().zip(&mean).zip(&f.potentials) {
            *o = o.max((f64::from(*p) - m).abs());
        }
    }
    let largest = own.iter().fold(0.0_f64, |a, b| a.max(*b));
    game.swing = Some(Swing {
        revision: v.revision,
        frames: frames.len(),
        mean,
        own,
        largest,
    });
}

/// The scale of a tube's electric field map, and its arrows: E at every
/// `TUBE_ARROW_SPACING`-th node. The scale is the 99th percentile of |E| over the nodes
/// (not the maximum: the singular corners would wash out the colours; the 90th
/// saturated most of the board).
pub(super) fn tube_field(
    level: &Level,
    run: &level::tube::TubeRun,
    frame: &level::tube::TubeFrame,
) -> (f64, TubeArrows) {
    let m = level.grid.max_node();
    // At most about `TUBE_ARROWS` across (a large board's every second node was thousands
    // of field sums per frame).
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let spacing = TUBE_ARROW_SPACING.max(m[0].max(m[1]).unsigned_abs() as usize / TUBE_ARROWS);
    let mut arrows = Vec::new();
    let mut magnitudes = Vec::new();
    for y in (0..=m[1]).step_by(spacing) {
        for x in (0..=m[0]).step_by(spacing) {
            let p = level.grid.position([x, y, 0]);
            let e = run.e_at(frame, p);
            magnitudes.push(e.length());
            arrows.push((draw::to_vec2(p), draw::to_vec2(e)));
        }
    }
    magnitudes.sort_by(f64::total_cmp);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let scale = magnitudes
        .get((0.99 * magnitudes.len() as f64) as usize)
        .copied()
        .unwrap_or(1.0)
        .max(1e-30);
    (scale, arrows)
}

pub(super) fn animate(time: Res<Time>, mut game: ResMut<Game>) {
    if !game.animate {
        return;
    }
    // Longest flight among the shown trajectories.
    let per_shot = game.editor.level.flights_per_shot();
    let end = game
        .flights
        .iter()
        .enumerate()
        .filter(|(i, _)| game.show_all_shots || i / per_shot == game.active_shot)
        .filter_map(|(_, s)| s.preview.as_ref().map(|p| p.flight_time))
        .fold(0.0, f64::max);
    // Beams: the last particle's end.
    let end = game
        .beams
        .iter()
        .filter_map(BeamView::shown)
        .flat_map(|p| {
            p.paths
                .iter()
                .filter_map(|path| path.last().map(|(t, _)| *t))
        })
        .fold(end, f64::max);
    // Tube levels: the preview run's end so far.
    let end = game
        .tube_preview
        .as_ref()
        .and_then(|v| v.run.frames.last().map(|f| f.t))
        .map_or(end, |t| end.max(t));
    // While it streams in, playback waits at its last frame (no loop yet).
    let streaming = game.tube_preview.as_ref().is_some_and(|v| !v.done);
    if end <= 0.0 {
        return;
    }
    // With a radiation goal, play on until the last radiation has crossed the arena, so
    // that it is seen leaving towards the (far) receiver.
    let level = &game.editor.level;
    let end = match level.physics.c {
        Some(c) if level.has_radiation_goal() => {
            let b = level.bounds();
            end + (b.max - b.min).length() / c
        }
        _ => end,
    };
    // Playback in internal time units per second, with a pause at the end. Tube levels
    // play at 1 unit per second (a run of 20 takes 20 s at speed 1; at 4 units per second
    // it was over before the flow could be watched).
    let tube = level.is_tube();
    let rate = if tube { 1.0 } else { 4.0 };
    let before = game.anim_time;
    game.anim_time += time.delta_secs_f64() * game.playback_speed * rate;
    // The field maps' playback waits for the field (`update_em`): it holds the time shown
    // until the field catches up, never sets it back (a field left behind while another
    // map was shown would pull the time back to it, and every earlier time asked for
    // restarts the field from t = 0: the playback stuck near the start). A capture at a
    // set time (`EM_TIME`) waits for it instead (`dev_capture`).
    // Once the field has reached the run's end it can go no further (a closed box's
    // smoothed field shows a time its smoothing's delay before that): the playback runs
    // on to the end's pause and loops.
    if tube
        && std::env::var_os("EM_TIME").is_none()
        && matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow))
        && let Some(s) = &game.em_snapshot
        && (streaming || !s.at_end)
    {
        let ahead = 0.25 * game.playback_speed.max(0.25);
        game.anim_time = game.anim_time.min((s.t + ahead).max(before));
    }
    if streaming {
        game.anim_time = game.anim_time.min(end);
    } else if game.anim_time > end + 2.0 * game.playback_speed.max(0.25) {
        game.anim_time = 0.0;
    }
}
