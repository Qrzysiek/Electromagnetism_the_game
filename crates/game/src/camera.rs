//! The camera: panning and fitting the board.

use super::*;

/// Moves a zoomed level's view by `shift` (board units), its centre kept on the board.
pub(super) fn pan_by(game: &mut Game, shift: Vec2) {
    let b = game.editor.level.bounds();
    #[allow(clippy::cast_possible_truncation)]
    let (lo, hi) = (
        Vec2::new(b.min.x as f32, b.min.y as f32),
        Vec2::new(b.max.x as f32, b.max.y as f32),
    );
    let centre = game.zoom_center.unwrap_or((lo + hi) * 0.5);
    game.zoom_center = Some((centre + shift).clamp(lo, hi));
}

/// Fits the level into the part of the window left of the side panel. The camera covers
/// the whole window (egui renders through it), so the fit is done with the projection
/// size and a horizontal offset rather than a viewport.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub(super) fn fit_camera(
    game: Res<Game>,
    window: Single<&Window, With<PrimaryWindow>>,
    camera: Single<(&mut Projection, &mut Transform), With<Camera2d>>,
) {
    // A minimized window has no area (0 × 0 on Windows): nothing to fit, the camera stays
    // as it was for the restored window. (Its width once made the clamp below
    // `clamp(0, −49)`, which panics: the game crashed on minimizing.)
    if window.physical_width() == 0 || window.physical_height() == 0 {
        return;
    }
    let (mut proj, mut transform) = camera.into_inner();
    let b = game.editor.level.bounds();
    let size = b.max - b.min;
    let center = (b.max + b.min) * 0.5;
    let win_w = window.physical_width() as f32;
    let win_h = window.physical_height() as f32;
    let panel_px = (win_w
        - game
            .panel_left_px
            .unwrap_or(win_w - PANEL_WIDTH * window.scale_factor()))
    .clamp(0.0, (win_w - 50.0).max(0.0));
    let avail_w = win_w - panel_px;
    // Room for the receivers of radiation goals, drawn just outside the edge.
    let margin = if game.editor.level.has_radiation_goal() {
        3.5
    } else {
        1.0
    };
    let (w, h) = (size.x as f32 + margin, size.y as f32 + margin);
    // World units per physical pixel (tube levels: zoomed by the wheel).
    let s = (w / avail_w).max(h / win_h) * game.zoom;
    #[allow(clippy::cast_possible_truncation)]
    let center = match game.zoom_center {
        Some(c) => physics::DVec3::new(
            f64::from(c.x).clamp(b.min.x, b.max.x),
            f64::from(c.y).clamp(b.min.y, b.max.y),
            0.0,
        ),
        None => center,
    };
    if let Projection::Orthographic(o) = &mut *proj {
        o.scaling_mode = ScalingMode::Fixed {
            width: win_w * s,
            height: win_h * s,
        };
    }
    transform.translation = Vec3::new(center.x as f32 + 0.5 * panel_px * s, center.y as f32, 0.0);
}
