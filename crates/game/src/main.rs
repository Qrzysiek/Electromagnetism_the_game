//! Electromagnetism – the game. Stage 1: 2D levels (a slice of the 3D world).

mod editor;
mod potential;
mod ui;
mod visuals;
mod worker;

use std::path::PathBuf;

use bevy::camera::ScalingMode;
use bevy::input::mouse::MouseWheel;
use bevy::prelude::*;
use bevy::window::PrimaryWindow;
use bevy_egui::input::EguiWantsInput;
use bevy_egui::{EguiPlugin, EguiPrimaryContextPass};
use level::Level;
use physics::DVec3;
use physics::trajectory::Outcome;
use physics::verify::Status;

use editor::Editor;
use worker::{Preview, Request, Response, Worker};

/// Width of the side panel in logical pixels.
pub const PANEL_WIDTH: f32 = 340.0;

/// A field line for drawing: polyline and arrowheads (position, unit direction of E).
pub type DrawnFieldLine = (Vec<Vec2>, Vec<(Vec2, Vec2)>);

#[derive(Resource)]
pub struct Game {
    pub levels: Vec<Level>,
    pub level_index: usize,
    pub editor: Editor,
    pub preview: Option<Preview>,
    /// Verification verdict for the current revision.
    pub verdict: Option<(Status, Outcome)>,
    pub sent_revision: u64,
    pub visuals_revision: u64,
    pub field_lines: Vec<DrawnFieldLine>,
    /// Distance between neighbouring field lines, in cells.
    pub field_line_spacing: f64,
    pub field_line_opacity: f32,
    /// (setup revision, density) the current field lines were computed for.
    pub field_lines_key: (u64, u64),
    pub show_potential: bool,
    pub show_field_lines: bool,
    pub animate: bool,
    pub playback_speed: f64,
    pub anim_time: f64,
    pub last_cursor_world: Option<Vec2>,
    /// Left edge of the side panel in physical pixels (reported by the UI).
    pub panel_left_px: Option<f32>,
}

impl Game {
    fn load_level(&mut self, index: usize) {
        self.level_index = index;
        self.editor = Editor::new(self.levels[index].clone());
        self.preview = None;
        self.verdict = None;
        self.anim_time = 0.0;
    }

    pub fn next_level(&mut self, step: isize) {
        let n = self.levels.len() as isize;
        let i = (self.level_index as isize + step).rem_euclid(n) as usize;
        self.load_level(i);
    }

    pub fn select_level(&mut self, index: usize) {
        if index != self.level_index {
            self.load_level(index);
        }
    }

    /// Solved = verified arrival for the current setup.
    pub fn solved(&self) -> bool {
        matches!(self.verdict, Some((Status::Verified, Outcome::Arrived)))
    }
}

/// Gizmo group for field lines: translucent and without joints (joints overlap the
/// segments and would double-blend).
#[derive(Default, Reflect, GizmoConfigGroup)]
struct FieldLineGizmos;

#[derive(Resource)]
struct PhysicsWorker(Worker);

#[derive(Resource)]
struct PotentialQuad {
    material: Handle<potential::PotentialMaterial>,
    entity: Entity,
}

fn levels_dir() -> PathBuf {
    let candidates = [
        std::env::current_dir().ok().map(|d| d.join("levels")),
        std::env::current_exe()
            .ok()
            .and_then(|e| e.parent().map(|p| p.join("levels"))),
        Some(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels")),
    ];
    candidates
        .into_iter()
        .flatten()
        .find(|p| p.is_dir())
        .expect("levels directory not found")
}

fn load_levels() -> Vec<Level> {
    let dir = levels_dir();
    let mut paths: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("read levels directory")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension().is_some_and(|e| e == "json") && !p.ends_with("golden_hashes.json")
        })
        .collect();
    paths.sort();
    paths
        .iter()
        .filter_map(|p| {
            let text = std::fs::read_to_string(p).ok()?;
            match Level::from_json(&text) {
                Ok(l) => Some(l),
                Err(e) => {
                    eprintln!("skipping {}: {e}", p.display());
                    None
                }
            }
        })
        .collect()
}

fn main() {
    let levels = load_levels();
    assert!(!levels.is_empty(), "no levels found");
    let editor = Editor::new(levels[0].clone());
    App::new()
        .insert_resource(ClearColor(Color::srgb(0.06, 0.06, 0.08)))
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Electromagnetism – the game".into(),
                resolution: (1400, 860).into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin::default())
        .add_plugins(potential::PotentialPlugin)
        .init_gizmo_group::<FieldLineGizmos>()
        .insert_resource(Game {
            levels,
            level_index: 0,
            editor,
            preview: None,
            verdict: None,
            sent_revision: 0,
            visuals_revision: 0,
            field_lines: Vec::new(),
            field_line_spacing: 1.5,
            field_line_opacity: 0.2,
            field_lines_key: (0, 0),
            show_potential: true,
            show_field_lines: false,
            animate: true,
            playback_speed: 1.0,
            anim_time: 0.0,
            last_cursor_world: None,
            panel_left_px: None,
        })
        .insert_resource(PhysicsWorker(Worker::spawn()))
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                input,
                sync_physics,
                update_field_lines,
                poll_physics,
                animate,
                fit_camera,
                draw,
            )
                .chain(),
        )
        .add_systems(EguiPrimaryContextPass, ui::panel)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut gizmo_store: ResMut<GizmoConfigStore>,
) {
    commands.spawn((
        Camera2d,
        Projection::Orthographic(OrthographicProjection::default_2d()),
        Msaa::Sample8,
    ));
    // Smooth polylines: round joints avoid notches between segments.
    let (config, _) = gizmo_store.config_mut::<DefaultGizmoConfigGroup>();
    config.line.width = 2.0;
    config.line.joints = GizmoLineJoint::Round(4);
    let (config, _) = gizmo_store.config_mut::<FieldLineGizmos>();
    config.line.width = 1.6;
    config.line.joints = GizmoLineJoint::None;
    // Potential map: a unit quad scaled to the world bounds, shaded on the GPU.
    let material = materials.add(potential::PotentialMaterial {
        params: potential::PotentialParams {
            charges: [Vec4::ZERO; potential::MAX_CHARGES],
            count: 0,
            u_a: 0.0,
        },
    });
    let entity = commands
        .spawn((
            Mesh2d(meshes.add(Rectangle::new(1.0, 1.0))),
            MeshMaterial2d(material.clone()),
            Transform::from_xyz(0.0, 0.0, -10.0),
        ))
        .id();
    commands.insert_resource(PotentialQuad { material, entity });
}

/// World position (cell units) of the grid node nearest to `p`.
fn nearest_node(game: &Game, p: Vec2) -> [i64; 3] {
    let s = f64::from(game.editor.subdivision());
    #[allow(clippy::cast_possible_truncation)]
    let r = |v: f32| (f64::from(v) * s).round() as i64;
    [r(p.x), r(p.y), 0]
}

fn input(
    keys: Res<ButtonInput<KeyCode>>,
    mouse: Res<ButtonInput<MouseButton>>,
    mut wheel: MessageReader<MouseWheel>,
    egui_input: Res<EguiWantsInput>,
    window: Single<&Window, With<PrimaryWindow>>,
    camera: Single<(&Camera, &GlobalTransform)>,
    mut game: ResMut<Game>,
) {
    let game = &mut *game;
    // The panel has no text fields, so game shortcuts always apply (egui keeps keyboard
    // focus on the last clicked button, which would otherwise block them).
    {
        let step = if keys.pressed(KeyCode::ShiftLeft) || keys.pressed(KeyCode::ShiftRight) {
            5
        } else {
            1
        };
        let e = &mut game.editor;
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
        if keys.just_pressed(KeyCode::Space) || keys.just_pressed(KeyCode::Enter) {
            let _ = e.place();
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
        if keys.just_pressed(KeyCode::KeyF) {
            game.show_field_lines = !game.show_field_lines;
        }
        if keys.just_pressed(KeyCode::KeyV) {
            game.show_potential = !game.show_potential;
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

    if egui_input.wants_pointer_input() || egui_input.is_pointer_over_area() {
        wheel.clear();
        return;
    }
    let (cam, cam_transform) = *camera;
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
    // The node under the pointer, if the pointer is on the grid.
    let node = nearest_node(game, world);
    let on_grid = game.editor.level.grid.contains(node);
    // Follow the mouse only when it moves, so keyboard control is not overridden.
    if on_grid
        && game
            .last_cursor_world
            .is_none_or(|last| last.distance(world) > 1e-3)
    {
        game.editor.set_cursor(node);
        game.last_cursor_world = Some(world);
    }
    if on_grid && mouse.just_pressed(MouseButton::Left) {
        game.editor.set_cursor(node);
        let _ = game.editor.place();
    }
    if on_grid && mouse.just_pressed(MouseButton::Right) {
        game.editor.set_cursor(node);
        game.editor.remove();
    }
    for ev in wheel.read() {
        if ev.y > 0.0 {
            game.editor.cycle_magnitude(1);
        } else if ev.y < 0.0 {
            game.editor.cycle_magnitude(-1);
        }
    }
}

/// Sends changed setups to the physics thread and recomputes the visuals.
fn sync_physics(
    mut game: ResMut<Game>,
    worker: Res<PhysicsWorker>,
    quad: Res<PotentialQuad>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut transforms: Query<&mut Transform>,
) {
    let revision = game.editor.revision + ((game.level_index as u64) << 48);
    if revision == game.sent_revision {
        return;
    }
    game.sent_revision = revision;
    game.verdict = None;
    let level = game.editor.level.clone();
    let scenario = level.scenario(&game.editor.placement);
    worker.0.submit(Request {
        revision,
        scenario: scenario.clone(),
        tolerances: level.tolerances(),
    });

    if game.visuals_revision != revision {
        game.visuals_revision = revision;
        let charges: Vec<(DVec3, f64)> = level
            .level_charges
            .iter()
            .chain(&game.editor.placement)
            .map(|c| (level.grid.position(c.node), c.charge))
            .collect();
        if let Some(mut m) = materials.get_mut(&quad.material) {
            m.params = potential::params(&scenario, &charges);
        }
        let bounds = scenario.bounds.expect("bounds");
        if let Ok(mut t) = transforms.get_mut(quad.entity) {
            let size = bounds.max - bounds.min;
            let center = (bounds.max + bounds.min) * 0.5;
            #[allow(clippy::cast_possible_truncation)]
            {
                t.scale = Vec3::new(size.x as f32, size.y as f32, 1.0);
                t.translation = Vec3::new(center.x as f32, center.y as f32, -10.0);
            }
        }
    }
}

/// Recomputes the field lines when they are shown and the setup or density changed.
fn update_field_lines(mut game: ResMut<Game>) {
    let key = (game.sent_revision, game.field_line_spacing.to_bits());
    if !game.show_field_lines || key == game.field_lines_key {
        return;
    }
    game.field_lines_key = key;
    let scenario = game.editor.level.scenario(&game.editor.placement);
    game.field_lines = visuals::field_lines(&scenario, game.field_line_spacing)
        .into_iter()
        .map(|l| {
            (
                l.points.into_iter().map(to_vec2).collect(),
                l.arrows
                    .into_iter()
                    .map(|(p, d)| (to_vec2(p), to_vec2(d)))
                    .collect(),
            )
        })
        .collect();
}

fn poll_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    for r in worker.0.poll() {
        match r {
            Response::Preview(p) if p.revision == game.sent_revision => {
                game.preview = Some(p);
            }
            Response::Verified {
                revision,
                status,
                outcome,
            } if revision == game.sent_revision => {
                game.verdict = Some((status, outcome));
            }
            _ => {}
        }
    }
}

fn animate(time: Res<Time>, mut game: ResMut<Game>) {
    let Some(end) = game.preview.as_ref().map(|p| p.flight_time) else {
        return;
    };
    if !game.animate {
        return;
    }
    // Playback in internal time units per second, with a pause at the end.
    game.anim_time += time.delta_secs_f64() * game.playback_speed * 4.0;
    if game.anim_time > end + 2.0 * game.playback_speed.max(0.25) {
        game.anim_time = 0.0;
    }
}

/// Fits the level into the part of the window left of the side panel. The camera covers
/// the whole window (egui renders through it), so the fit is done with the projection
/// size and a horizontal offset rather than a viewport.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
fn fit_camera(
    game: Res<Game>,
    window: Single<&Window, With<PrimaryWindow>>,
    camera: Single<(&mut Projection, &mut Transform), With<Camera2d>>,
) {
    let (mut proj, mut transform) = camera.into_inner();
    let b = game.editor.level.bounds();
    let size = b.max - b.min;
    let center = (b.max + b.min) * 0.5;
    let win_w = window.physical_width().max(1) as f32;
    let win_h = window.physical_height().max(1) as f32;
    let panel_px = (win_w
        - game
            .panel_left_px
            .unwrap_or(win_w - PANEL_WIDTH * window.scale_factor()))
    .clamp(0.0, win_w - 50.0);
    let avail_w = win_w - panel_px;
    let (w, h) = (size.x as f32 + 1.0, size.y as f32 + 1.0);
    // World units per physical pixel.
    let s = (w / avail_w).max(h / win_h);
    if let Projection::Orthographic(o) = &mut *proj {
        o.scaling_mode = ScalingMode::Fixed {
            width: win_w * s,
            height: win_h * s,
        };
    }
    transform.translation = Vec3::new(center.x as f32 + 0.5 * panel_px * s, center.y as f32, 0.0);
}

#[allow(clippy::cast_possible_truncation)]
fn to_vec2(v: DVec3) -> Vec2 {
    Vec2::new(v.x as f32, v.y as f32)
}

fn draw(
    game: Res<Game>,
    mut gizmos: Gizmos,
    mut line_gizmos: Gizmos<FieldLineGizmos>,
    quad: Res<PotentialQuad>,
    mut vis: Query<&mut Visibility>,
) {
    if let Ok(mut v) = vis.get_mut(quad.entity) {
        *v = if game.show_potential {
            Visibility::Visible
        } else {
            Visibility::Hidden
        };
    }
    let level = &game.editor.level;
    let grid = level.grid;
    let m = grid.max_node();
    let top_right = to_vec2(grid.position(m));

    // Grid: cell lines, plus fainter refined nodes.
    let cell = Color::srgba(1.0, 1.0, 1.0, 0.035);
    for i in 0..=grid.nx {
        #[allow(clippy::cast_precision_loss)]
        let x = i as f32;
        gizmos.line_2d(Vec2::new(x, 0.0), Vec2::new(x, top_right.y), cell);
    }
    for j in 0..=grid.ny {
        #[allow(clippy::cast_precision_loss)]
        let y = j as f32;
        gizmos.line_2d(Vec2::new(0.0, y), Vec2::new(top_right.x, y), cell);
    }
    if grid.subdivision > 1 {
        let fine = Color::srgba(1.0, 1.0, 1.0, 0.015);
        for i in 0..=m[0] {
            let x = to_vec2(grid.position([i, 0, 0])).x;
            gizmos.line_2d(Vec2::new(x, 0.0), Vec2::new(x, top_right.y), fine);
        }
        for j in 0..=m[1] {
            let y = to_vec2(grid.position([0, j, 0])).y;
            gizmos.line_2d(Vec2::new(0.0, y), Vec2::new(top_right.x, y), fine);
        }
    }

    // World bounds.
    let b = level.bounds();
    let (bmin, bmax) = (to_vec2(b.min), to_vec2(b.max));
    gizmos.rect_2d(
        (bmin + bmax) * 0.5,
        bmax - bmin,
        Color::srgba(1.0, 0.3, 0.3, 0.5),
    );

    // Detector.
    let d0 = to_vec2(grid.position(level.detector.min));
    let d1 = to_vec2(grid.position(level.detector.max));
    let det_color = if game.solved() {
        Color::srgb(0.3, 1.0, 0.4)
    } else {
        Color::srgb(0.2, 0.8, 0.3)
    };
    gizmos.rect_2d((d0 + d1) * 0.5, (d1 - d0).abs(), det_color);
    gizmos.rect_2d(
        (d0 + d1) * 0.5,
        (d1 - d0).abs() - Vec2::splat(0.15),
        det_color.with_alpha(0.5),
    );

    // Field lines.
    if game.show_field_lines {
        let color = Color::srgba(0.95, 0.95, 0.75, game.field_line_opacity);
        let (head, spread) = (0.3, 0.45_f32);
        for (points, arrows) in &game.field_lines {
            line_gizmos.linestrip_2d(points.iter().copied(), color);
            // Arrowheads pointing along E.
            for &(p, d) in arrows {
                let back = -d * head;
                let tip = p + d * (0.5 * head);
                line_gizmos.line_2d(tip, tip + Vec2::from_angle(spread).rotate(back), color);
                line_gizmos.line_2d(tip, tip + Vec2::from_angle(-spread).rotate(back), color);
            }
        }
    }

    // Launch point and direction.
    let a = to_vec2(grid.position(level.launch.node));
    let dir = level.launch.direction;
    #[allow(clippy::cast_possible_truncation)]
    let dir = Vec2::new(dir[0] as f32, dir[1] as f32).normalize_or_zero();
    gizmos.circle_2d(a, 0.25, Color::srgb(1.0, 1.0, 1.0));
    gizmos.arrow_2d(a, a + dir * 1.5, Color::srgb(1.0, 1.0, 1.0));

    // Charges.
    #[allow(clippy::cast_possible_truncation)]
    let radius = level.physics.charge_radius as f32;
    let q_max = level.limits.magnitudes.iter().copied().fold(1.0, f64::max);
    let draw_charge = |gizmos: &mut Gizmos, c: &level::Charge, player: bool| {
        let p = to_vec2(grid.position(c.node));
        let color = if c.charge > 0.0 {
            Color::srgb(1.0, 0.35, 0.3)
        } else {
            Color::srgb(0.35, 0.6, 1.0)
        };
        // Filled disc from concentric circles; the halo grows with |Q|.
        for k in 1..=6 {
            #[allow(clippy::cast_precision_loss)]
            gizmos.circle_2d(p, radius * k as f32 / 6.0, color);
        }
        #[allow(clippy::cast_possible_truncation)]
        let halo = radius * (1.0 + 0.8 * (c.charge.abs() / q_max) as f32);
        gizmos.circle_2d(p, halo, color.with_alpha(0.4));
        let s = radius * 0.6;
        let ink = Color::srgb(0.05, 0.05, 0.05);
        gizmos.line_2d(p - Vec2::X * s, p + Vec2::X * s, ink);
        if c.charge > 0.0 {
            gizmos.line_2d(p - Vec2::Y * s, p + Vec2::Y * s, ink);
        }
        let ring = if player {
            Color::srgb(1.0, 1.0, 1.0)
        } else {
            Color::srgb(0.5, 0.5, 0.5)
        };
        gizmos.circle_2d(p, radius * 1.15, ring);
    };
    for c in &level.level_charges {
        draw_charge(&mut gizmos, c, false);
    }
    for c in &game.editor.placement {
        draw_charge(&mut gizmos, c, true);
    }

    // Cursor.
    let cur = to_vec2(grid.position(game.editor.cursor));
    let cursor_color = if game.editor.selected_charge() > 0.0 {
        Color::srgb(1.0, 0.6, 0.5)
    } else {
        Color::srgb(0.6, 0.8, 1.0)
    };
    gizmos.rect_2d(cur, Vec2::splat(radius * 2.8), cursor_color);

    // Trajectory.
    if let Some(p) = &game.preview {
        let color = match game.verdict {
            None => Color::srgb(0.95, 0.95, 0.95),
            Some((Status::Verified, Outcome::Arrived)) => Color::srgb(0.3, 1.0, 0.4),
            Some((Status::Verified, _)) => Color::srgb(1.0, 0.55, 0.2),
            Some(_) => Color::srgb(1.0, 0.9, 0.2),
        };
        gizmos.linestrip_2d(p.path.iter().map(|q| to_vec2(q.x)), color);
        let end = to_vec2(p.path.last().expect("path has points").x);
        if !matches!(p.outcome, Outcome::Arrived) {
            let s = 0.3;
            gizmos.line_2d(end - Vec2::splat(s), end + Vec2::splat(s), color);
            gizmos.line_2d(end + Vec2::new(-s, s), end + Vec2::new(s, -s), color);
        }
        if game.animate
            && let Some(pt) = ui::point_at(p, game.anim_time)
        {
            let x = to_vec2(pt.x);
            gizmos.circle_2d(x, 0.18, Color::srgb(1.0, 1.0, 0.6));
            gizmos.circle_2d(x, 0.1, Color::srgb(1.0, 1.0, 0.6));
            let f = to_vec2(pt.force);
            if f.length() > 1e-6 {
                // Force arrow, length ∝ log(1 + |F|), direction exact.
                let len = (1.0 + f.length()).ln() * 2.0;
                gizmos.arrow_2d(x, x + f.normalize() * len, Color::srgb(1.0, 0.8, 0.2));
            }
        }
    }
}
