//! The mouse and the keyboard on the board: placing, dragging, panning, zooming,
//! the views' keys.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn input(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    mouse: Res<ButtonInput<MouseButton>>,
    mut wheel: MessageReader<MouseWheel>,
    egui_input: Res<EguiWantsInput>,
    window: Single<&Window, With<PrimaryWindow>>,
    camera: Single<(&Camera, &GlobalTransform)>,
    mut game: ResMut<Game>,
) {
    let game = &mut *game;
    // Game shortcuts apply unless a text field has focus. (egui's own "wants keyboard"
    // flag is also set by a focused button, which would block shortcuts after a click.)
    if !game.text_focus {
        let step = if keys.pressed(KeyCode::ShiftLeft) || keys.pressed(KeyCode::ShiftRight) {
            5
        } else {
            1
        };
        // Ctrl + arrows pan a zoomed tube level (without a mouse: a laptop's touchpad).
        let ctrl = keys.pressed(KeyCode::ControlLeft) || keys.pressed(KeyCode::ControlRight);
        if ctrl && game.editor.level.is_tube() && game.zoom < 1.0 {
            let dir = Vec2::new(
                f32::from(u8::from(keys.pressed(KeyCode::ArrowRight)))
                    - f32::from(u8::from(keys.pressed(KeyCode::ArrowLeft))),
                f32::from(u8::from(keys.pressed(KeyCode::ArrowUp)))
                    - f32::from(u8::from(keys.pressed(KeyCode::ArrowDown))),
            );
            if dir != Vec2::ZERO {
                let b = game.editor.level.bounds();
                #[allow(clippy::cast_possible_truncation)]
                let span = (b.max.x - b.min.x).max(b.max.y - b.min.y) as f32;
                // Most of the view across per second.
                let shift = dir * (0.8 * game.zoom * span * time.delta_secs());
                pan_by(game, shift);
            }
        }
        let e = &mut game.editor;
        if !ctrl {
            if keys.just_pressed(KeyCode::ArrowLeft) {
                e.move_cursor(-step, 0);
            }
            if keys.just_pressed(KeyCode::ArrowRight) {
                e.move_cursor(step, 0);
            }
            if keys.just_pressed(KeyCode::ArrowUp) {
                e.move_cursor(0, step);
            }
            if keys.just_pressed(KeyCode::ArrowDown) {
                e.move_cursor(0, -step);
            }
        }
        if keys.just_pressed(KeyCode::KeyG) {
            if e.grabbed.is_some() {
                e.drop_grabbed();
            } else {
                e.grab();
            }
        }
        if keys.just_pressed(KeyCode::Escape) {
            e.cancel_grab();
        }
        if keys.just_pressed(KeyCode::Space) || keys.just_pressed(KeyCode::Enter) {
            if e.grabbed.is_some() {
                e.drop_grabbed();
            } else {
                let _ = e.place();
            }
        }
        if keys.just_pressed(KeyCode::Delete)
            || keys.just_pressed(KeyCode::Backspace)
            || keys.just_pressed(KeyCode::KeyX)
        {
            e.remove();
        }
        if keys.just_pressed(KeyCode::KeyS) {
            e.flip_sign();
        }
        if keys.just_pressed(KeyCode::KeyE) || keys.just_pressed(KeyCode::Tab) {
            e.cycle_magnitude(1);
        }
        if keys.just_pressed(KeyCode::KeyQ) {
            e.cycle_magnitude(-1);
        }
        if keys.just_pressed(KeyCode::KeyW) {
            e.cycle_omega(if step > 1 { -1 } else { 1 });
        }
        if keys.just_pressed(KeyCode::KeyR) {
            e.rotate(if step > 1 { -1 } else { 1 });
        }
        if keys.just_pressed(KeyCode::KeyM) {
            e.toggle_kind();
        }
        if keys.just_pressed(KeyCode::KeyC) {
            e.clear();
        }
        for (k, f) in [
            (KeyCode::Digit1, 1),
            (KeyCode::Digit2, 2),
            (KeyCode::Digit3, 3),
            (KeyCode::Digit4, 4),
        ] {
            if keys.just_pressed(k) {
                e.set_refinement(f);
            }
        }
        if keys.just_pressed(KeyCode::BracketRight) {
            game.cycle_shot(1);
        }
        if keys.just_pressed(KeyCode::BracketLeft) {
            game.cycle_shot(-1);
        }
        if keys.just_pressed(KeyCode::KeyH) {
            game.show_all_shots = !game.show_all_shots;
        }
        if keys.just_pressed(KeyCode::KeyF) {
            game.show_field_lines = !game.show_field_lines;
        }
        if keys.just_pressed(KeyCode::KeyZ) || keys.just_pressed(KeyCode::Home) {
            game.zoom = 1.0;
            game.zoom_center = None;
        }
        if keys.just_pressed(KeyCode::KeyV) {
            game.map = next_map(&game.editor.level, game.map);
        }
        if keys.just_pressed(KeyCode::KeyA) {
            game.animate = !game.animate;
        }
        if keys.just_pressed(KeyCode::KeyN) || keys.just_pressed(KeyCode::PageDown) {
            game.next_level(1);
        }
        if keys.just_pressed(KeyCode::KeyP) || keys.just_pressed(KeyCode::PageUp) {
            game.next_level(-1);
        }
    }
    // The board probe while I is held (a capture's EM_PROBE=1 holds it).
    game.probe = (!game.text_focus && keys.pressed(KeyCode::KeyI))
        || std::env::var("EM_PROBE").is_ok_and(|v| v == "1");
    if game.probe
        && !game.editor.level.is_tube()
        && game
            .probe_field
            .as_ref()
            .is_none_or(|(revision, _)| *revision != game.sent_revision)
    {
        let field = game.editor.level.field(&game.editor.placement).0;
        game.probe_field = Some((game.sent_revision, field));
    }

    if egui_input.wants_pointer_input() || egui_input.is_pointer_over_area() {
        wheel.clear();
        return;
    }
    let (cam, cam_transform) = *camera;
    // Captures inspect a point of the board without a pointer: EM_HOVER=x,y.
    game.hover_world = std::env::var("EM_HOVER").ok().and_then(|v| {
        let (x, y) = v.split_once(',')?;
        Some(Vec2::new(x.trim().parse().ok()?, y.trim().parse().ok()?))
    });
    if game.hover_world.is_some() {
        return;
    }
    let Some(cursor) = window.cursor_position() else {
        return;
    };
    // Ignore the pointer over the side panel.
    if game
        .panel_left_px
        .is_some_and(|left| cursor.x * window.scale_factor() >= left)
    {
        wheel.clear();
        return;
    }
    let Ok(world) = cam.viewport_to_world_2d(cam_transform, cursor) else {
        return;
    };
    if !game.mouse_drag && game.arrow_drag.is_none() && game.pan_anchor.is_none() {
        game.hover_world = Some(world);
    }
    // A zoomed tube level pans by dragging: the middle button anywhere, the left one on
    // empty board (where a press would place nothing); the point grabbed stays under the
    // pointer.
    let zoomed_tube = game.editor.level.is_tube() && game.zoom < 1.0;
    if zoomed_tube && mouse.just_pressed(MouseButton::Middle) {
        game.pan_anchor = Some(world);
    }
    if let Some(anchor) = game.pan_anchor {
        if zoomed_tube && (mouse.pressed(MouseButton::Middle) || mouse.pressed(MouseButton::Left)) {
            pan_by(game, anchor - world);
            return;
        }
        game.pan_anchor = None;
    }
    // The node under the pointer, if the pointer is on the grid.
    let node = nearest_node(game, world);
    let on_grid = game.editor.level.grid.contains(node);
    // A tunable electrode under the pointer highlights its supply (and its slider's badge).
    if let Some(n) = game.editor.tunable_at(node) {
        game.focus_supply(n, 0.15);
    }
    // Follow the mouse only when it moves, so keyboard control is not overridden.
    if on_grid
        && game
            .last_cursor_world
            .is_none_or(|last| last.distance(world) > 1e-3)
    {
        game.editor.set_cursor(node);
        game.last_cursor_world = Some(world);
    }
    let (pressed, released, right) = (
        mouse.just_pressed(MouseButton::Left),
        mouse.just_released(MouseButton::Left),
        mouse.just_pressed(MouseButton::Right),
    );
    if on_grid
        && game.sandbox.active
        && (pressed || released || right)
        && sandbox::pointer_at(game, world, pressed, released, right)
    {
        return;
    }
    // A free charge's velocity: press on the handle at the tip of its arrow and drag.
    let here = [f64::from(world.x), f64::from(world.y)];
    if pressed {
        game.arrow_drag = editor::arrow_handle_at(
            &game.editor.level,
            &game.editor.placement,
            here,
            0.35,
            game.arrow_measure,
        );
    }
    if let Some(i) = game.arrow_drag {
        if mouse.pressed(MouseButton::Left) {
            if let Some(e) = game.editor.placement.get(i).copied() {
                let p = game.editor.level.grid.position(e.node);
                let (angle, speed) = editor::dragged_velocity(
                    &game.editor.level,
                    [p.x, p.y],
                    here,
                    game.arrow_measure,
                );
                game.editor.set_velocity(i, angle, speed);
            }
        } else {
            game.arrow_drag = None;
        }
        return;
    }
    // A cut being drawn: it ends where the button is released.
    if let Some(a) = game.cut_start {
        if !mouse.pressed(MouseButton::Left) {
            game.cut_start = None;
            if on_grid && node != a {
                let _ = game.editor.place_cut(a, node);
            } else {
                game.editor.message =
                    Some("Drag from one point of a metal to another to draw a cut.".into());
            }
        }
        return;
    }
    // With the cut tool, a press away from the player's elements starts a cut.
    if on_grid
        && pressed
        && game.editor.kind == level::ElementKind::Cut
        && game.editor.left(level::ElementKind::Cut) > 0
        && (!game.sandbox.active || game.sandbox.tool == sandbox::Tool::PlayerElement)
    {
        game.editor.set_cursor(node);
        if game.editor.element_at_cursor().is_none() {
            game.cut_start = Some(node);
            return;
        }
    }
    // Left press on one of the player's elements picks it up (drag), elsewhere places;
    // where nothing can be placed on a zoomed tube level, it starts a pan.
    if on_grid && pressed {
        game.editor.set_cursor(node);
        if game.editor.grab() {
            game.mouse_drag = true;
        } else if game.editor.place().is_err() && zoomed_tube {
            game.pan_anchor = Some(world);
        }
    } else if !on_grid && pressed && zoomed_tube {
        game.pan_anchor = Some(world);
    }
    if game.mouse_drag && !mouse.pressed(MouseButton::Left) {
        game.mouse_drag = false;
        game.editor.drop_grabbed();
    }
    if on_grid && right {
        game.editor.set_cursor(node);
        game.editor.remove();
    }
    // Tube levels (nothing to place by magnitude): the wheel zooms about the pointer.
    if game.editor.level.is_tube() {
        for ev in wheel.read() {
            if ev.y == 0.0 {
                continue;
            }
            let old = game.zoom;
            let new = (old * 1.2_f32.powf(-ev.y.signum())).clamp(0.1, 1.0);
            let b = game.editor.level.bounds();
            #[allow(clippy::cast_possible_truncation)]
            let fit = Vec2::new(
                (0.5 * (b.max.x + b.min.x)) as f32,
                (0.5 * (b.max.y + b.min.y)) as f32,
            );
            let centre = game.zoom_center.unwrap_or(fit);
            game.zoom = new;
            game.zoom_center = (new < 1.0).then(|| world + (centre - world) * (new / old));
        }
        return;
    }
    for ev in wheel.read() {
        if ev.y > 0.0 {
            game.editor.cycle_magnitude(1);
        } else if ev.y < 0.0 {
            game.editor.cycle_magnitude(-1);
        }
    }
}
