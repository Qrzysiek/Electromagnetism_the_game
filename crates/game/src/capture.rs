//! The capture mode (`EM_CAPTURE` and its companions, CLAUDE.md): a level opened,
//! set up as asked, its window saved at a frame, and the game closed.

use super::*;

/// Developer capture for testing without input: with `EM_CAPTURE=<file.png>` the game
/// opens level `EM_LEVEL` (1-based) or the level file `EM_LEVEL_FILE`, enters the sandbox if `EM_SANDBOX=1`, shows the
/// flight of shot `EM_SHOT` under disturbance `EM_DISTURBANCE` (1-based), selects the map
/// `EM_MAP` (potential, magnetic, waves, particle, off), places `EM_PLACE` (JSON list of
/// elements, or "reference"; in units of 1/`EM_GRID` cell with the grid refined 1-4×), holds
/// the animation at `EM_TIME`, lets the
/// physics settle, saves a screenshot of its own window and exits. With `EM_WAIT=1` the
/// screenshot waits until the verdicts are in. No clicks or keys are sent to the desktop.
pub(super) fn dev_capture(
    mut commands: Commands,
    mut game: ResMut<Game>,
    mut window: Single<&mut Window, With<PrimaryWindow>>,
    mut frame: Local<u32>,
    mut shot_at: Local<Option<u32>>,
    mut exit: MessageWriter<AppExit>,
) {
    let Ok(path) = std::env::var("EM_CAPTURE") else {
        return;
    };
    *frame += 1;
    // Frame of the screenshot (EM_FRAME, default 120); the exit follows 60 frames later.
    let shot = std::env::var("EM_FRAME")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(120);
    let wait = std::env::var("EM_WAIT").is_ok_and(|v| v == "1");
    // Tube levels: also the whole streamed display run (a refused setup has none: it is
    // settled at once; waiting for its stream timed out every refused case of the
    // pairwise check).
    let streamed = !game.editor.level.is_tube()
        || game.invalid_setup().is_some()
        || game.tube_preview.as_ref().is_some_and(|v| v.done);
    // A field map: its field at the time shown.
    let field = !game.editor.level.is_tube()
        || !matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow))
        || game
            .em_snapshot
            .as_ref()
            .is_some_and(|s| (s.t - game.anim_time).abs() < 0.05);
    let settled = game.progress() == Progress::Done && streamed && field;
    if shot_at.is_none() && *frame >= shot && (!wait || settled) {
        *shot_at = Some(*frame);
        commands
            .spawn(bevy::render::view::screenshot::Screenshot::primary_window())
            .observe(bevy::render::view::screenshot::save_to_disk(path.clone()));
    }
    if shot_at.is_some_and(|f| *frame == f + 60) {
        exit.write(AppExit::Success);
    }
    if let Some(t) = std::env::var("EM_TIME")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
    {
        game.anim_time = t;
    }
    // A zoom on the board's centre (tube levels: 1 the whole board).
    if let Some(z) = std::env::var("EM_ZOOM")
        .ok()
        .and_then(|v| v.parse::<f32>().ok())
    {
        game.zoom = z.clamp(0.1, 1.0);
    }
    match *frame {
        3 => {
            if let Some(i) = std::env::var("EM_LEVEL")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&i| i >= 1 && i <= game.levels.len())
            {
                game.select_level(i - 1);
            }
            // Any level file (the pairwise element check, scripts/pairs.py, writes them).
            if let Some(level) = std::env::var("EM_LEVEL_FILE")
                .ok()
                .and_then(|p| std::fs::read_to_string(p).ok())
                .and_then(|s| level::Level::from_json(&s).ok())
            {
                let i = game.level_index;
                game.levels[i] = level;
                game.load_level(i);
            }
            if std::env::var("EM_SANDBOX").is_ok_and(|v| v == "1") {
                sandbox::enter(&mut game);
            }
            // The flight shown (1-based): its shot and its disturbance.
            if let Some(i) = std::env::var("EM_SHOT")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&i| i >= 1)
            {
                game.active_shot = i - 1;
            }
            if let Some(i) = std::env::var("EM_DISTURBANCE")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&i| i >= 1)
            {
                game.active_disturbance = i - 1;
            }
            game.radiation_only = std::env::var("EM_RAD_ONLY").is_ok_and(|v| v == "1");
            game.neglected_only = std::env::var("EM_NEGLECTED").is_ok_and(|v| v == "1");
            if std::env::var("EM_LINES").is_ok_and(|v| v == "1") {
                game.show_field_lines = true;
            }
            if std::env::var("EM_HARDCORE").is_ok_and(|v| v == "1") {
                game.editor.set_continuous(true);
            }
            match std::env::var("EM_QUANTITY").as_deref() {
                Ok("E") => game.field_quantity = radiation::FieldQuantity::E,
                Ok("S") => game.field_quantity = radiation::FieldQuantity::S,
                _ => {}
            }
            if let Ok(p) = std::env::var("EM_PART") {
                game.flow_part = match p.as_str() {
                    "own" => radiation::FlowPart::Own,
                    "exchange" => radiation::FlowPart::Exchange,
                    "external" => radiation::FlowPart::External,
                    _ => radiation::FlowPart::Total,
                };
            }
            if std::env::var("EM_TRACERS").is_ok_and(|v| v == "1") {
                game.flow_tracers = true;
            }
            // A supply held in focus (its number in the panel, 1-based): its electrodes
            // highlighted, as when its slider is hovered.
            if let Some(k) = std::env::var("EM_FOCUS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                && let Some((n, _)) = game.editor.level.supply_targets().get(k.wrapping_sub(1))
            {
                let n = *n;
                game.focus_supply(n, 3600.0);
            }
            if std::env::var("EM_AVERAGE").is_ok_and(|v| v == "0") {
                game.flow_average = false;
            }
            if let Some(tau) = std::env::var("EM_TAU").ok().and_then(|v| v.parse().ok()) {
                game.map_tau = tau;
            }
            if let Some(d) = std::env::var("EM_RANGE").ok().and_then(|v| v.parse().ok()) {
                game.field_range_decades = d;
            }
            if std::env::var("EM_LINEAR").is_ok_and(|v| v == "1") {
                game.field_linear = true;
            }
            if let Some(g) = std::env::var("EM_GAIN").ok().and_then(|v| v.parse().ok()) {
                game.field_gain_decades = g;
                game.tube_gain_decades = g;
            }
            if let Some(d) = std::env::var("EM_RANGE").ok().and_then(|v| v.parse().ok()) {
                game.tube_range_decades = d;
            }
            if let Ok(m) = std::env::var("EM_MAP") {
                game.map = map_named(&m);
            }
            // The grid refinement (1-4, as the Grid buttons): EM_PLACE's nodes are then in
            // units of 1/n cell.
            if let Some(n) = std::env::var("EM_GRID")
                .ok()
                .and_then(|v| v.parse::<u32>().ok())
                .filter(|n| (1..=4).contains(n))
            {
                game.editor.set_refinement(n);
            }
            if let Ok(p) = std::env::var("EM_PLACE") {
                // Player elements as JSON, e.g. the level's reference solution: "reference".
                let placement = if p == "reference" {
                    game.editor.level.reference_solution.clone()
                } else {
                    serde_json::from_str(&p).unwrap_or_default()
                };
                game.editor.set_placement(placement);
            }
        }
        10 if std::env::var("EM_SPHERE_TEST").is_ok() => {
            game.editor.edit_level(|l| {
                l.conductors.push(level::Conductor {
                    center: [20, 14, 0],
                    radius: 2.0,
                    bias: level::ConductorBias::Grounded,
                });
            });
        }
        60 if std::env::var("EM_SPHERE_TEST").is_ok() => {
            game.editor.edit_level(|l| {
                l.conductors.clear();
            });
        }
        _ => {}
    }
    // EM_MAP_AT ("frame:map,..."): the map switched at those frames, as a player would (a
    // regression check: a field map shown again after another kept playing).
    if let Ok(seq) = std::env::var("EM_MAP_AT") {
        for item in seq.split(',') {
            if let Some((f, m)) = item.split_once(':')
                && f.trim().parse::<u32>() == Ok(*frame)
            {
                game.map = map_named(m.trim());
            }
        }
    }
    // EM_LEVEL_AT ("frame:level,..."): the level (1-based) selected at those frames, as
    // the level list does (a regression check: a view of the previous level once stayed
    // over the next one).
    if let Ok(seq) = std::env::var("EM_LEVEL_AT") {
        for item in seq.split(',') {
            if let Some((f, l)) = item.split_once(':')
                && f.trim().parse::<u32>() == Ok(*frame)
                && let Ok(i) = l.trim().parse::<usize>()
                && (1..=game.levels.len()).contains(&i)
            {
                game.select_level(i - 1);
            }
        }
    }
    // EM_MINIMIZE=1: the window minimized at frame 20 and restored at frame 40 (a
    // regression check: minimizing once crashed the game).
    if std::env::var("EM_MINIMIZE").is_ok_and(|v| v == "1") {
        match *frame {
            20 => window.set_minimized(true),
            40 => window.set_minimized(false),
            _ => {}
        }
    }
}
