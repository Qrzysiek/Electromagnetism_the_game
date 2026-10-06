//! Electromagnetism – the game. Stage 1: 2D levels (a slice of the 3D world).

#![allow(clippy::disallowed_methods)] // pictures only: the flights come from `physics`

mod curriculum;
mod draw;
mod editor;
mod level_editor;
mod math;
mod potential;
mod radiation;
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
/// Field lines computed off the render thread: the setup key, and where they arrive.
pub type FieldLinesJob = (
    (u64, u64, usize),
    std::sync::Mutex<std::sync::mpsc::Receiver<Vec<DrawnFieldLine>>>,
);

/// A beam flight (beam levels, one per disturbance): the preview of every particle and
/// their verdicts.
#[derive(Clone, Debug, Default)]
pub struct BeamView {
    pub preview: Option<worker::BeamPreview>,
    /// Setup revision the preview belongs to (older while a new one is computed).
    pub preview_revision: u64,
    pub verified: Option<Vec<(Status, Outcome)>>,
    /// The verdict's own flight, when its model differs from the preview's (the exact
    /// retarded interaction at finite c, metal at verification resolution).
    pub verdict_flight: Option<worker::BeamPreview>,
}

impl BeamView {
    /// The flight every view shows (paths, animation, energy budget, field views, counts):
    /// the verdict's own flight once it has arrived, the preview until then.
    pub fn shown(&self) -> Option<&worker::BeamPreview> {
        self.verdict_flight.as_ref().or(self.preview.as_ref())
    }
}

/// State of the computation for the current setup.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Progress {
    /// New flights are being computed; the paths shown belong to the previous setup.
    Flying,
    /// Flights shown are current; their verification is running.
    Verifying,
    Done,
}

/// Physics results of one flight (a shot under one disturbance) for the current setup.
#[derive(Clone, Debug, Default)]
pub struct FlightView {
    /// Latest preview (may belong to an older revision until the new one arrives).
    pub preview: Option<Preview>,
    /// Setup revision the preview belongs to.
    pub preview_revision: u64,
    /// Verification verdict for the current revision.
    pub verdict: Option<(Status, Outcome)>,
}

#[derive(Resource)]
pub struct Game {
    pub levels: Vec<Level>,
    /// File of each level (same order as `levels`).
    pub level_paths: Vec<PathBuf>,
    /// Arc and tier of each level (same order; `None` for custom levels).
    pub places: Vec<Option<curriculum::Place>>,
    pub level_index: usize,
    pub editor: Editor,
    /// One view per flight, shot-major (`Level::flight_of`).
    pub flights: Vec<FlightView>,
    /// Beam levels: one view per flight (disturbance).
    pub beams: Vec<BeamView>,
    pub active_shot: usize,
    /// Disturbance whose flight the details panel shows.
    pub active_disturbance: usize,
    /// Show every shot's trajectory, not only the active one.
    pub show_all_shots: bool,
    pub sent_revision: u64,
    /// Measured computational cost (sandbox resource meters) and the setup revision it
    /// belongs to; kept from the previous setup until the new one is measured.
    pub cost: Option<(u64, level::cost::Cost)>,
    /// Levels with circuits: the circuit view (plots, model notes) and the setup revision
    /// it belongs to.
    pub circuit: Option<(u64, worker::CircuitView)>,
    /// (revision, active shot, active disturbance, all shots shown, mode, flight shown)
    /// the field map was computed for.
    pub map_key: (u64, usize, usize, bool, Option<MapMode>, bool),
    pub field_lines: Vec<DrawnFieldLine>,
    /// Distance between neighbouring field lines, in cells.
    pub field_line_spacing: f64,
    pub field_line_opacity: f32,
    /// (setup revision, spacing) the current field lines were computed for.
    pub field_lines_key: (u64, u64, usize),
    /// Field lines being computed off the render thread (a few seconds near metal): the
    /// setup key they belong to, and where they arrive.
    pub field_lines_job: Option<FieldLinesJob>,
    /// The key of the field lines wanted now (older jobs stop).
    pub field_lines_wanted: std::sync::Arc<std::sync::atomic::AtomicU64>,
    /// Field map shown under the scene (`None`: off).
    pub map: Option<MapMode>,
    pub show_field_lines: bool,
    /// E arrows on the time-dependent field views.
    pub show_field_arrows: bool,
    /// Particle field view: only the radiation (acceleration) part of the field.
    pub radiation_only: bool,
    /// Beam field views: only the part the quasi-static interaction leaves out.
    pub neglected_only: bool,
    /// Hardcore sliders: logarithmic scale for magnitudes and for frequencies.
    pub log_magnitude: bool,
    pub log_omega: bool,
    /// Colour quantity of the field views.
    pub field_quantity: radiation::FieldQuantity,
    /// Energy-flow view: the part shown, tracers, and the average over a period of the
    /// antennas and waves.
    pub flow_part: radiation::FlowPart,
    pub flow_tracers: bool,
    pub flow_average: bool,
    /// Dynamic range of the field views, in decades.
    pub field_range_decades: f64,
    /// Linear colour scale (else logarithmic over `field_range_decades`).
    pub field_linear: bool,
    /// Gain of the linear scale, in decades (full colour at 10^−gain of the full field's
    /// scale).
    pub field_gain_decades: f64,
    pub animate: bool,
    pub playback_speed: f64,
    pub anim_time: f64,
    pub last_cursor_world: Option<Vec2>,
    pub sandbox: sandbox::Sandbox,
    /// A player element is being dragged with the mouse.
    pub mouse_drag: bool,
    /// The free charge whose velocity arrow is being dragged.
    pub arrow_drag: Option<usize>,
    /// What velocity arrows measure (speed, or rapidity in relativistic levels).
    pub arrow_measure: editor::ArrowMeasure,
    /// Set by the UI while a text field has focus (game shortcuts are suspended).
    pub text_focus: bool,
    /// Left edge of the side panel in physical pixels (reported by the UI).
    pub panel_left_px: Option<f32>,
}

impl Game {
    pub fn new(levels: Vec<Level>, level_paths: Vec<PathBuf>) -> Self {
        let editor = Editor::new(levels[0].clone());
        let arrow_measure = editor::default_measure(&editor.level);
        Self {
            levels,
            places: curriculum::places(&level_paths),
            level_paths,
            level_index: 0,
            editor,
            flights: Vec::new(),
            beams: Vec::new(),
            active_shot: 0,
            active_disturbance: 0,
            show_all_shots: true,
            sent_revision: 0,
            cost: None,
            circuit: None,
            map_key: (0, 0, 0, false, None, false),
            field_lines: Vec::new(),
            field_line_spacing: 1.5,
            field_line_opacity: 0.2,
            field_lines_key: (0, 0, 0),
            field_lines_job: None,
            field_lines_wanted: std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0)),
            map: Some(MapMode::Potential),
            show_field_lines: false,
            show_field_arrows: true,
            radiation_only: false,
            neglected_only: false,
            log_magnitude: true,
            log_omega: false,
            field_quantity: radiation::FieldQuantity::Bz,
            flow_part: radiation::FlowPart::Total,
            flow_tracers: false,
            flow_average: true,
            field_range_decades: 2.5,
            field_linear: false,
            field_gain_decades: 0.0,
            animate: true,
            playback_speed: 1.0,
            anim_time: 0.0,
            last_cursor_world: None,
            sandbox: sandbox::Sandbox::default(),
            mouse_drag: false,
            arrow_drag: None,
            arrow_measure,
            text_focus: false,
            panel_left_px: None,
        }
    }

    fn load_level(&mut self, index: usize) {
        self.level_index = index;
        self.editor = Editor::new(self.levels[index].clone());
        self.arrow_measure = editor::default_measure(&self.editor.level);
        self.arrow_drag = None;
        self.flights.clear();
        self.beams.clear();
        self.active_shot = 0;
        self.active_disturbance = 0;
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
        self.places = curriculum::places(&paths);
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

    /// Index of the flight of `shot` under disturbance `d`.
    pub fn flight_index(&self, shot: usize, d: usize) -> usize {
        shot * self.editor.level.flights_per_shot() + d
    }

    /// The flight shown in the details panel.
    /// What the physics thread is still doing for the current setup: new flights (the
    /// paths shown are from the previous setup), or the verification.
    pub fn progress(&self) -> Progress {
        let current = self.sent_revision;
        let (stale, verifying) = if self.editor.level.has_beams() {
            (
                self.beams.iter().any(|b| b.preview_revision != current),
                self.beams.iter().any(|b| b.verified.is_none()),
            )
        } else {
            (
                self.flights.iter().any(|f| f.preview_revision != current),
                self.flights.iter().any(|f| f.verdict.is_none()),
            )
        };
        if stale {
            Progress::Flying
        } else if verifying {
            Progress::Verifying
        } else {
            Progress::Done
        }
    }

    /// Energy unit of the displays for `shot` (energy bars, force arrow, potential map):
    /// its launch energy T₀, or, for a particle launched (almost) at rest, the largest
    /// kinetic energy it reaches in the flight shown (T₀ below 1 % of that). The flag says
    /// which.
    pub fn energy_unit(&self, shot: usize) -> (f64, bool) {
        let level = &self.editor.level;
        let t0 = level
            .shots
            .get(shot)
            .map_or(0.0, |s| s.launch.kinetic_energy);
        let d = self.active_disturbance.min(level.flights_per_shot() - 1);
        let peak = self
            .flights
            .get(self.flight_index(shot, d))
            .and_then(|f| f.preview.as_ref())
            .map_or(0.0, |p| {
                p.path.iter().map(|q| q.kinetic).fold(0.0, f64::max)
            });
        if t0 < 0.01 * peak {
            (peak, false)
        } else {
            (t0.max(1e-300), true)
        }
    }

    pub fn active_flight(&self) -> usize {
        let d = self
            .active_disturbance
            .min(self.editor.level.flights_per_shot() - 1);
        self.flight_index(self.active_shot, d)
    }

    pub fn flight_verdict(&self, flight: usize) -> Option<(Status, Outcome)> {
        self.flights.get(flight).and_then(|s| s.verdict)
    }

    /// Verdict of a shot over all its flights (`None` while computing): the first flight
    /// that does not arrive verified, or verified arrival if all do.
    pub fn verdict(&self, shot: usize) -> Option<(Status, Outcome)> {
        let n = self.editor.level.flights_per_shot();
        let mut all = Vec::with_capacity(n);
        for d in 0..n {
            all.push(self.flight_verdict(self.flight_index(shot, d))?);
        }
        all.iter()
            .copied()
            .find(|v| !matches!(v, (Status::Verified, Outcome::Arrived)))
            .or(all.first().copied())
    }

    /// Verdict of the level's disturbance `d` over all shots.
    pub fn disturbance_verdict(&self, d: usize) -> Option<(Status, Outcome)> {
        let mut all = Vec::new();
        for shot in 0..self.shot_count() {
            all.push(self.flight_verdict(self.flight_index(shot, d))?);
        }
        all.iter()
            .copied()
            .find(|v| !matches!(v, (Status::Verified, Outcome::Arrived)))
            .or(all.first().copied())
    }

    /// Solved = every shot arrives, verified (beams: every beam shot reaches its
    /// verified transmission in every flight).
    pub fn solved(&self) -> bool {
        let level = &self.editor.level;
        if level.has_beams() {
            return !self.beams.is_empty()
                && (0..self.beams.len()).all(|d| {
                    level.shots.iter().enumerate().all(|(s, shot)| {
                        let need = shot.beam.map_or(1.0, |b| b.transmission);
                        #[allow(clippy::cast_precision_loss)]
                        let ok = self
                            .beam_transmission(d, s)
                            .is_some_and(|(ok, n)| ok as f64 >= need * n as f64 - 1e-9);
                        ok
                    }) && self.goals_arrived(d)
                });
        }
        let n = self.shot_count();
        n > 0
            && (0..n).all(|i| matches!(self.verdict(i), Some((Status::Verified, Outcome::Arrived))))
    }

    /// Whether every goal particle (a level's free particle with a detector) arrives,
    /// verified, in flight `d`. Goal particle `k` is the beam's "shot" `shots + k`.
    pub fn goals_arrived(&self, d: usize) -> bool {
        let level = &self.editor.level;
        let n = level.shots.len();
        level
            .free_particles
            .iter()
            .enumerate()
            .filter(|(_, f)| f.detector.is_some())
            .all(|(k, _)| {
                self.beam_transmission(d, n + k)
                    .is_some_and(|(ok, _)| ok == 1)
            })
    }

    /// Verified arrivals and particle count of beam shot `shot` in flight `d`; `None`
    /// while it is being computed.
    pub fn beam_transmission(&self, d: usize, shot: usize) -> Option<(usize, usize)> {
        let v = self.beams.get(d)?;
        let (p, r) = (v.shown()?, v.verified.as_ref()?);
        if r.len() != p.shots.len() {
            return None;
        }
        let mine: Vec<usize> = (0..p.shots.len()).filter(|&i| p.shots[i] == shot).collect();
        let ok = mine
            .iter()
            .filter(|&&i| matches!(r[i], (Status::Verified, Outcome::Arrived)))
            .count();
        Some((ok, mine.len()))
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
                    .filter(|p| level::is_level_file(p))
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
    // A developer capture (`EM_CAPTURE`) opens its window off-screen, unfocused and
    // without a taskbar entry, so it never covers or takes focus from a running game.
    let capture = std::env::var_os("EM_CAPTURE").is_some();
    App::new()
        .insert_resource(ClearColor(Color::srgb(0.06, 0.06, 0.08)))
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Electromagnetism – the game".into(),
                resolution: (1400, 860).into(),
                position: if capture {
                    WindowPosition::At(IVec2::new(-20000, -20000))
                } else {
                    WindowPosition::Automatic
                },
                focused: !capture,
                skip_taskbar: capture,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin::default())
        .add_plugins(potential::PotentialPlugin)
        .add_plugins(radiation::FieldViewPlugin)
        .init_gizmo_group::<FieldLineGizmos>()
        .insert_resource(Game::new(levels, level_paths))
        .insert_resource(PhysicsWorker(Worker::spawn()))
        .add_systems(Startup, (setup, radiation::setup))
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
        .add_systems(Update, dev_capture)
        .add_systems(Update, radiation::update.after(animate).before(draw::draw))
        .run();
}

/// Next map in the V-key cycle, skipping views the level has nothing to show in.
pub fn next_map(level: &Level, map: Option<MapMode>) -> Option<MapMode> {
    let order = [
        Some(MapMode::Potential),
        Some(MapMode::Magnetic),
        Some(MapMode::Waves),
        Some(MapMode::ParticleField),
        Some(MapMode::Total),
        None,
    ];
    let available = |m: Option<MapMode>| match m {
        Some(MapMode::Waves) => radiation::waves_available(level),
        Some(MapMode::ParticleField) => radiation::particle_field_available(level),
        _ => true,
    };
    let i = order.iter().position(|&m| m == map).unwrap_or(0);
    (1..=order.len())
        .map(|k| order[(i + k) % order.len()])
        .find(|&m| available(m))
        .unwrap_or(None)
}

/// Developer capture for testing without input: with `EM_CAPTURE=<file.png>` the game
/// opens level `EM_LEVEL` (1-based), enters the sandbox if `EM_SANDBOX=1`, shows the
/// flight of shot `EM_SHOT` under disturbance `EM_DISTURBANCE` (1-based), selects the map
/// `EM_MAP` (potential, magnetic, waves, particle, off), places `EM_PLACE` (JSON list of
/// elements, or "reference"; in units of 1/`EM_GRID` cell with the grid refined 1-4×), holds
/// the animation at `EM_TIME`, lets the
/// physics settle, saves a screenshot of its own window and exits. With `EM_WAIT=1` the
/// screenshot waits until the verdicts are in. No clicks or keys are sent to the desktop.
fn dev_capture(
    mut commands: Commands,
    mut game: ResMut<Game>,
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
    if shot_at.is_none() && *frame >= shot && (!wait || game.progress() == Progress::Done) {
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
    match *frame {
        3 => {
            if let Some(i) = std::env::var("EM_LEVEL")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&i| i >= 1 && i <= game.levels.len())
            {
                game.select_level(i - 1);
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
            if std::env::var("EM_AVERAGE").is_ok_and(|v| v == "0") {
                game.flow_average = false;
            }
            if let Some(d) = std::env::var("EM_RANGE").ok().and_then(|v| v.parse().ok()) {
                game.field_range_decades = d;
            }
            if std::env::var("EM_LINEAR").is_ok_and(|v| v == "1") {
                game.field_linear = true;
            }
            if let Some(g) = std::env::var("EM_GAIN").ok().and_then(|v| v.parse().ok()) {
                game.field_gain_decades = g;
            }
            if let Ok(m) = std::env::var("EM_MAP") {
                game.map = match m.as_str() {
                    "magnetic" => Some(MapMode::Magnetic),
                    "waves" => Some(MapMode::Waves),
                    "particle" => Some(MapMode::ParticleField),
                    "total" => Some(MapMode::Total),
                    "off" => None,
                    _ => Some(MapMode::Potential),
                };
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
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut buffers: ResMut<Assets<bevy::render::storage::ShaderBuffer>>,
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
        items: buffers.add(potential::buffer(Vec::new())),
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
pub fn nearest_node(game: &Game, p: Vec2) -> [i64; 3] {
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
    // Left press on one of the player's elements picks it up (drag), elsewhere places.
    if on_grid && pressed {
        game.editor.set_cursor(node);
        if game.editor.grab() {
            game.mouse_drag = true;
        } else {
            let _ = game.editor.place();
        }
    }
    if game.mouse_drag && !mouse.pressed(MouseButton::Left) {
        game.mouse_drag = false;
        game.editor.drop_grabbed();
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
    let flights = level.flight_count();
    let per_shot = level.flights_per_shot();
    let request_level = level.clone();
    // Keep old previews on screen until new ones arrive; verdicts are recomputed.
    game.flights.resize(flights, FlightView::default());
    for s in &mut game.flights {
        s.verdict = None;
    }
    if request_level.has_beams() {
        game.beams.resize(per_shot, BeamView::default());
        for b in &mut game.beams {
            b.verified = None;
            // The old setup's best flight stays on screen (faded) until the new preview
            // arrives.
            if let Some(f) = b.verdict_flight.take() {
                b.preview = Some(f);
            }
        }
    } else {
        game.beams.clear();
    }
    if game.active_shot >= n {
        game.active_shot = 0;
    }
    if game.active_disturbance >= per_shot {
        game.active_disturbance = 0;
    }
    let placement = game.editor.placement.clone();
    worker.0.submit(Request {
        revision,
        level: request_level,
        placement,
    });
}

/// Updates the GPU field map for the active shot and setup.
fn update_map(
    mut game: ResMut<Game>,
    quad: Res<PotentialQuad>,
    mut materials: ResMut<Assets<potential::PotentialMaterial>>,
    mut buffers: ResMut<Assets<bevy::render::storage::ShaderBuffer>>,
    mut transforms: Query<&mut Transform>,
) {
    // The energy unit depends on the flight (particles launched at rest): redraw when it
    // arrives.
    let has_flight = game
        .flights
        .get(game.active_flight())
        .is_some_and(|f| f.preview_revision == game.sent_revision && f.preview.is_some());
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
    );
    if key == game.map_key {
        return;
    }
    game.map_key = key;
    let Some(mode @ (MapMode::Potential | MapMode::Magnetic)) = game.map else {
        return;
    };
    let level = &game.editor.level;
    if level.shots.is_empty() {
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

/// Recomputes the field lines (the electric field does not depend on the shot, but on the
/// disturbance) when they are shown and the setup, spacing or disturbance changed.
fn update_field_lines(mut game: ResMut<Game>) {
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
    if !game.show_field_lines || key == game.field_lines_key || game.shot_count() == 0 {
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

fn poll_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    sandbox::poll(&mut game);
    let current = game.sent_revision;
    for r in worker.0.poll() {
        match r {
            Response::Preview {
                revision,
                flight: shot,
                preview,
            } if revision == current => {
                if let Some(s) = game.flights.get_mut(shot) {
                    s.preview = Some(preview);
                    s.preview_revision = revision;
                }
            }
            Response::Verified {
                revision,
                flight: shot,
                status,
                outcome,
            } if revision == current => {
                if let Some(s) = game.flights.get_mut(shot) {
                    s.verdict = Some((status, outcome));
                }
            }
            Response::BeamPreview {
                revision,
                flight,
                preview,
            } if revision == current => {
                if let Some(b) = game.beams.get_mut(flight) {
                    b.preview = Some(preview);
                    b.preview_revision = revision;
                    b.verdict_flight = None;
                }
            }
            Response::BeamVerified {
                revision,
                flight,
                results,
                run,
            } if revision == current => {
                if let Some(b) = game.beams.get_mut(flight) {
                    b.verified = Some(results);
                    b.verdict_flight = run;
                }
            }
            Response::Cost { revision, cost } if revision == current => {
                game.cost = Some((revision, cost));
            }
            Response::Circuit { revision, view } if revision == current => {
                game.circuit = Some((revision, view));
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
    // Room for the receivers of radiation goals, drawn just outside the edge.
    let margin = if game.editor.level.has_radiation_goal() {
        3.5
    } else {
        1.0
    };
    let (w, h) = (size.x as f32 + margin, size.y as f32 + margin);
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
