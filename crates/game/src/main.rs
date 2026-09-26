//! Electromagnetism – the game. Stage 1: 2D levels (a slice of the 3D world).

mod draw;
mod editor;
mod potential;
mod sandbox;
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
use physics::trajectory::Outcome;
use physics::verify::Status;

use editor::Editor;
use potential::MapMode;
use worker::{Preview, Request, Response, Worker};

/// Width of the side panel in logical pixels.
pub const PANEL_WIDTH: f32 = 340.0;

/// A field line for drawing: polyline and arrowheads (position, unit direction of E).
pub type DrawnFieldLine = (Vec<Vec2>, Vec<(Vec2, Vec2)>);

/// Physics results of one shot for the current setup.
#[derive(Clone, Debug, Default)]
pub struct ShotView {
    /// Latest preview (may belong to an older revision until the new one arrives).
    pub preview: Option<Preview>,
    /// Verification verdict for the current revision.
    pub verdict: Option<(Status, Outcome)>,
}

#[derive(Resource)]
pub struct Game {
    pub levels: Vec<Level>,
    /// File of each level (same order as `levels`).
    pub level_paths: Vec<PathBuf>,
    pub level_index: usize,
    pub editor: Editor,
    pub shots: Vec<ShotView>,
    pub active_shot: usize,
    /// Show every shot's trajectory, not only the active one.
    pub show_all_shots: bool,
    pub sent_revision: u64,
    /// (revision, active shot, mode) the field map was computed for.
    pub map_key: (u64, usize, Option<MapMode>),
    pub field_lines: Vec<DrawnFieldLine>,
    /// Distance between neighbouring field lines, in cells.
    pub field_line_spacing: f64,
    pub field_line_opacity: f32,
    /// (setup revision, spacing) the current field lines were computed for.
    pub field_lines_key: (u64, u64),
    /// Field map shown under the scene (`None`: off).
    pub map: Option<MapMode>,
    pub show_field_lines: bool,
    pub animate: bool,
    pub playback_speed: f64,
    pub anim_time: f64,
    pub last_cursor_world: Option<Vec2>,
    pub sandbox: sandbox::Sandbox,
    /// Set by the UI while a text field has focus (game shortcuts are suspended).
    pub text_focus: bool,
    /// Left edge of the side panel in physical pixels (reported by the UI).
    pub panel_left_px: Option<f32>,
}

impl Game {
    fn load_level(&mut self, index: usize) {
        self.level_index = index;
        self.editor = Editor::new(self.levels[index].clone());
        self.shots.clear();
        self.active_shot = 0;
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
            if self.sandbox.active {
                sandbox::enter(self);
            }
        }
    }

    /// Re-reads the level files (after saving) and selects `select` if given. The level
    /// being edited stays loaded.
    pub fn reload_levels(&mut self, select: Option<&std::path::Path>) {
        let (levels, paths) = load_levels();
        if levels.is_empty() {
            return;
        }
        self.levels = levels;
        self.level_paths = paths;
        self.level_index = select
            .and_then(|s| self.level_paths.iter().position(|p| p == s))
            .unwrap_or(0);
    }

    pub fn shot_count(&self) -> usize {
        self.editor.level.shots.len()
    }

    pub fn select_shot(&mut self, shot: usize) {
        if shot < self.shot_count() {
            self.active_shot = shot;
            self.anim_time = 0.0;
        }
    }

    pub fn cycle_shot(&mut self, step: isize) {
        let n = self.shot_count().max(1) as isize;
        self.select_shot((self.active_shot as isize + step).rem_euclid(n) as usize);
    }

    /// Verdict of a shot (`None` while computing).
    pub fn verdict(&self, shot: usize) -> Option<(Status, Outcome)> {
        self.shots.get(shot).and_then(|s| s.verdict)
    }

    /// Solved = every shot arrives, verified.
    pub fn solved(&self) -> bool {
        let n = self.shot_count();
        n > 0
            && (0..n).all(|i| matches!(self.verdict(i), Some((Status::Verified, Outcome::Arrived))))
    }
}

/// Gizmo group for field lines: translucent and without joints (joints overlap the
/// segments and would double-blend).
#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct FieldLineGizmos;

#[derive(Resource)]
struct PhysicsWorker(Worker);

#[derive(Resource)]
pub struct PotentialQuad {
    material: Handle<potential::PotentialMaterial>,
    pub entity: Entity,
}

pub fn levels_dir() -> PathBuf {
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

/// Shipped levels from `levels/`, then custom ones from `levels/custom/`.
fn load_levels() -> (Vec<Level>, Vec<PathBuf>) {
    let dir = levels_dir();
    let list = |d: &std::path::Path| -> Vec<PathBuf> {
        let mut v: Vec<PathBuf> = std::fs::read_dir(d)
            .map(|rd| {
                rd.filter_map(|e| e.ok().map(|e| e.path()))
                    .filter(|p| {
                        p.extension().is_some_and(|e| e == "json")
                            && !p.ends_with("golden_hashes.json")
                    })
                    .collect()
            })
            .unwrap_or_default();
        v.sort();
        v
    };
    let mut levels = Vec::new();
    let mut paths = Vec::new();
    for p in list(&dir).into_iter().chain(list(&dir.join("custom"))) {
        let Ok(text) = std::fs::read_to_string(&p) else {
            continue;
        };
        match Level::from_json(&text) {
            Ok(l) if !l.shots.is_empty() => {
                levels.push(l);
                paths.push(p);
            }
            Ok(_) => eprintln!("skipping {}: no shots", p.display()),
            Err(e) => eprintln!("skipping {}: {e}", p.display()),
        }
    }
    (levels, paths)
}

fn main() {
    let (levels, level_paths) = load_levels();
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
            level_paths,
            level_index: 0,
            editor,
            shots: Vec::new(),
            active_shot: 0,
            show_all_shots: true,
            sent_revision: 0,
            map_key: (0, 0, None),
            field_lines: Vec::new(),
            field_line_spacing: 1.5,
            field_line_opacity: 0.2,
            field_lines_key: (0, 0),
            map: Some(MapMode::Potential),
            show_field_lines: false,
            animate: true,
            playback_speed: 1.0,
            anim_time: 0.0,
            last_cursor_world: None,
            sandbox: sandbox::Sandbox::default(),
            text_focus: false,
            panel_left_px: None,
        })
        .insert_resource(PhysicsWorker(Worker::spawn()))
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                input,
                sync_physics,
                update_map,
                update_field_lines,
                poll_physics,
                animate,
                fit_camera,
                draw::draw,
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
    // Field map: a unit quad scaled to the world bounds, shaded on the GPU.
    let material = materials.add(potential::PotentialMaterial {
        params: potential::PotentialParams::default(),
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
    // Game shortcuts apply unless a text field has focus. (egui's own "wants keyboard"
    // flag is also set by a focused button, which would block shortcuts after a click.)
    if !game.text_focus {
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
        if keys.just_pressed(KeyCode::KeyV) {
            game.map = match game.map {
                Some(MapMode::Potential) => Some(MapMode::Magnetic),
                Some(MapMode::Magnetic) => None,
                None => Some(MapMode::Potential),
            };
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
    let (pressed, released, right) = (
        mouse.just_pressed(MouseButton::Left),
        mouse.just_released(MouseButton::Left),
        mouse.just_pressed(MouseButton::Right),
    );
    if on_grid
        && game.sandbox.active
        && (pressed || released || right)
        && sandbox::pointer(game, node, world, pressed, released, right)
    {
        return;
    }
    if on_grid && pressed {
        game.editor.set_cursor(node);
        let _ = game.editor.place();
    }
    if on_grid && right {
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

/// Current setup revision (changes with every edit and level switch).
fn revision(game: &Game) -> u64 {
    game.editor.revision + ((game.level_index as u64) << 48)
}

/// Sends changed setups to the physics thread.
fn sync_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    let revision = revision(&game);
    if revision == game.sent_revision {
        return;
    }
    game.sent_revision = revision;
    let level = &game.editor.level;
    let n = level.shots.len();
    let scenarios = level.scenarios(&game.editor.placement);
    let tolerances = level.tolerances();
    // Keep old previews on screen until new ones arrive; verdicts are recomputed.
    game.shots.resize(n, ShotView::default());
    for s in &mut game.shots {
        s.verdict = None;
    }
    if game.active_shot >= n {
        game.active_shot = 0;
    }
    worker.0.submit(Request {
        revision,
        scenarios,
        tolerances,
    });
}

/// Updates the GPU field map for the active shot and setup.
fn update_map(
    mut game: ResMut<Game>,
    quad: Res<PotentialQuad>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut transforms: Query<&mut Transform>,
) {
    let key = (game.sent_revision, game.active_shot, game.map);
    if key == game.map_key {
        return;
    }
    game.map_key = key;
    let Some(mode) = game.map else {
        return;
    };
    let level = &game.editor.level;
    if level.shots.is_empty() {
        return;
    }
    let scenario = level.scenario(game.active_shot, &game.editor.placement);
    if let Some(mut m) = materials.get_mut(&quad.material) {
        m.params = potential::params(
            &scenario,
            level.physics.charge_radius,
            level.physics.magnet_radius,
            mode,
        );
    }
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
}

/// Recomputes the field lines (the electric field does not depend on the shot) when they
/// are shown and the setup or spacing changed.
fn update_field_lines(mut game: ResMut<Game>) {
    let key = (game.sent_revision, game.field_line_spacing.to_bits());
    if !game.show_field_lines || key == game.field_lines_key || game.shot_count() == 0 {
        return;
    }
    game.field_lines_key = key;
    let scenario = game.editor.level.scenario(0, &game.editor.placement);
    game.field_lines = visuals::field_lines(&scenario, game.field_line_spacing)
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
}

fn poll_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    sandbox::poll(&mut game);
    let current = game.sent_revision;
    for r in worker.0.poll() {
        match r {
            Response::Preview {
                revision,
                shot,
                preview,
            } if revision == current => {
                if let Some(s) = game.shots.get_mut(shot) {
                    s.preview = Some(preview);
                }
            }
            Response::Verified {
                revision,
                shot,
                status,
                outcome,
            } if revision == current => {
                if let Some(s) = game.shots.get_mut(shot) {
                    s.verdict = Some((status, outcome));
                }
            }
            _ => {}
        }
    }
}

fn animate(time: Res<Time>, mut game: ResMut<Game>) {
    if !game.animate {
        return;
    }
    // Longest flight among the shown trajectories.
    let end = game
        .shots
        .iter()
        .enumerate()
        .filter(|(i, _)| game.show_all_shots || *i == game.active_shot)
        .filter_map(|(_, s)| s.preview.as_ref().map(|p| p.flight_time))
        .fold(0.0, f64::max);
    if end <= 0.0 {
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
