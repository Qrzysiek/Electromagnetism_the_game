//! Electromagnetism – the game. Stage 1: 2D levels (a slice of the 3D world).

#![allow(clippy::disallowed_methods)] // pictures only: the flights come from `physics`

mod curriculum;
mod draw;
mod editor;
mod em;
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

/// Tube levels: the preview run as it streams in, for one setup revision.
pub struct TubeView {
    pub revision: u64,
    pub run: level::tube::TubeRun,
    /// The whole display run has arrived.
    pub done: bool,
}

/// Tube levels: the electric field map's colour scale and the arrows of the frame shown.
/// The scale is fixed from the frame at the goal window's end once it has arrived (so
/// that the colours do not flicker from frame to frame); before, from the frame shown.
pub struct TubeField {
    pub revision: u64,
    pub frame: usize,
    pub scale: f64,
    pub scale_final: bool,
    pub arrows: TubeArrows,
}

/// What a field map is computed for: revision, active shot, active disturbance, all
/// shots shown, mode, flight shown, the tube frame shown, and the tube maps' gain and range
/// (as bits).
/// What the field map was drawn for: setup revision, shot, disturbance, all shots, map,
/// flight arrived, a tube's picture (its frame, or while animating the time's bits), and a
/// tube map's gain and range.
pub type MapKey = (
    u64,
    usize,
    usize,
    bool,
    Option<MapMode>,
    bool,
    u64,
    (u64, u64),
);

/// How a tube level's particles are drawn: as the density of their charge (a glow in the
/// map, the default), as dots, or not at all (to see the map alone).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TubeParticles {
    Density,
    Dots,
    Hidden,
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
    /// A setup the physics cannot compute (`Level::setup_issues`) and its revision.
    pub setup_issues: Option<(u64, Vec<String>)>,
    /// Tube levels: the latest preview run and verdict, with their revisions (the
    /// preview is kept while the next one is computed).
    pub tube_preview: Option<TubeView>,
    pub tube_verdict: Option<(u64, level::tube::TubeVerdict)>,
    /// Tube levels: the electric field map's scale and arrows `(position, E)` of the
    /// preview's final state (computed when it arrives), with its revision.
    pub tube_field: Option<TubeField>,
    /// Tube levels' field maps: the field's thread (started on first use), the run it was
    /// given, the last time asked, the latest snapshot, and the run's colour units so far
    /// (its revision, the largest |B_z| and |S|).
    pub em: Option<em::EmWorker>,
    pub em_run: Option<u64>,
    pub em_asked: f64,
    pub em_snapshot: Option<em::EmSnapshot>,
    pub em_units: (u64, f64, f64),
    /// What the field map was computed for (`MapKey`).
    pub map_key: MapKey,
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
    /// Tube maps: the potential map's gain (decades: full colour at 10^−gain of the
    /// largest potential difference) and the electric map's range (decades below the top).
    pub tube_gain_decades: f64,
    /// The supply whose electrodes are highlighted (its centre node) and until when: its
    /// slider hovered or moved, or one of its electrodes hovered on the map.
    pub supply_focus: Option<(level::Node, std::time::Instant)>,
    pub tube_range_decades: f64,
    /// Tube levels: how the particles are drawn.
    pub tube_particles: TubeParticles,
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
    /// Tube levels: the view's zoom (1: the whole board, smaller: closer) and its centre
    /// when zoomed (world coordinates).
    pub zoom: f32,
    pub zoom_center: Option<Vec2>,
    /// A zoomed tube level dragged by the pointer: the board point held under it.
    pub pan_anchor: Option<Vec2>,
}

impl Game {
    /// Highlights a supply's electrodes for `seconds` from now.
    pub fn focus_supply(&mut self, centre: level::Node, seconds: f64) {
        let until = std::time::Instant::now() + std::time::Duration::from_secs_f64(seconds);
        if self
            .supply_focus
            .is_none_or(|(n, t)| n != centre || t < until)
        {
            self.supply_focus = Some((centre, until));
        }
    }

    /// The supply highlighted now, if any.
    pub fn focused_supply(&self) -> Option<level::Node> {
        self.supply_focus
            .filter(|(_, t)| std::time::Instant::now() < *t)
            .map(|(n, _)| n)
    }

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
            setup_issues: None,
            tube_preview: None,
            tube_verdict: None,
            tube_field: None,
            em: None,
            em_run: None,
            em_asked: -1.0,
            em_snapshot: None,
            em_units: (0, 0.0, 0.0),
            map_key: (0, 0, 0, false, None, false, 0, (0, 0)),
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
            tube_gain_decades: 0.0,
            supply_focus: None,
            tube_particles: TubeParticles::Density,
            tube_range_decades: 1.5,
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
            zoom: 1.0,
            zoom_center: None,
            pan_anchor: None,
        }
    }

    fn load_level(&mut self, index: usize) {
        self.level_index = index;
        // The new editor continues the old one's revisions: a fresh one starts again at 1,
        // and a level loaded again at the same index (capture mode's EM_LEVEL_FILE, which
        // replaces the level shown) had the revision already sent, so it was never
        // computed.
        let next = self.editor.revision + 1;
        self.zoom = 1.0;
        self.zoom_center = None;
        self.editor = Editor::new(self.levels[index].clone());
        self.editor.revision = self.editor.revision.max(next);
        self.arrow_measure = editor::default_measure(&self.editor.level);
        self.arrow_drag = None;
        self.flights.clear();
        self.beams.clear();
        self.tube_preview = None;
        self.tube_verdict = None;
        self.tube_field = None;
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
    /// The current setup's issues if the physics cannot compute it.
    pub fn invalid_setup(&self) -> Option<&[String]> {
        self.setup_issues
            .as_ref()
            .filter(|(r, _)| *r == self.sent_revision)
            .map(|(_, i)| i.as_slice())
    }

    /// Tube levels: the index of the frame shown, at the animation time (the last one
    /// when not animating).
    pub fn tube_frame(&self) -> Option<usize> {
        let v = self.tube_preview.as_ref()?;
        let n = v.run.frames.len();
        if n == 0 {
            return None;
        }
        Some(if self.animate {
            v.run
                .frames
                .partition_point(|f| f.t <= self.anim_time)
                .saturating_sub(1)
        } else {
            n - 1
        })
    }

    /// Tube levels: the picture shown, interpolated between the stored frames at the
    /// animation time (smooth at any playback speed); the last frame when not animating.
    pub fn tube_display_frame(&self) -> Option<level::tube::TubeFrame> {
        let v = self.tube_preview.as_ref()?;
        if self.animate {
            v.run.frame_between(self.anim_time)
        } else {
            v.run.frames.last().cloned()
        }
    }

    pub fn progress(&self) -> Progress {
        let current = self.sent_revision;
        if self.invalid_setup().is_some() {
            return Progress::Done;
        }
        let (stale, verifying) = if self.editor.level.is_tube() {
            (
                self.tube_preview
                    .as_ref()
                    .is_none_or(|v| v.revision != current),
                self.tube_verdict
                    .as_ref()
                    .is_none_or(|(r, _)| *r != current),
            )
        } else if self.editor.level.has_beams() {
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
        if level.is_tube() {
            return self.tube_verdict.as_ref().is_some_and(|(r, v)| {
                *r == self.sent_revision && v.status == level::tube::TubeStatus::Met
            });
        }
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
        v.sort_by_key(|p| level::level_order_key(p));
        v
    };
    let mut levels = Vec::new();
    let mut paths = Vec::new();
    for p in list(&dir).into_iter().chain(list(&dir.join("custom"))) {
        let Ok(text) = std::fs::read_to_string(&p) else {
            continue;
        };
        match Level::from_json(&text) {
            // Playable: shots to deliver, or a tube level (its goal is a current).
            Ok(l) if !l.shots.is_empty() || l.is_tube() => {
                levels.push(l);
                paths.push(p);
            }
            Ok(_) => eprintln!("skipping {}: no shots and no tube", p.display()),
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
                update_em,
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
        Some(MapMode::Electric),
        Some(MapMode::EnergyFlow),
        None,
    ];
    // Tube levels: their potential and electric field only (the other views are the 3D
    // slice's); the electric map only there.
    let available = |m: Option<MapMode>| match m {
        _ if level.is_tube() => matches!(
            m,
            Some(MapMode::Potential | MapMode::Electric | MapMode::Magnetic | MapMode::EnergyFlow)
                | None
        ),
        Some(MapMode::Waves) => radiation::waves_available(level),
        Some(MapMode::ParticleField) => radiation::particle_field_available(level),
        Some(MapMode::Electric | MapMode::EnergyFlow) => false,
        _ => true,
    };
    let i = order.iter().position(|&m| m == map).unwrap_or(0);
    (1..=order.len())
        .map(|k| order[(i + k) % order.len()])
        .find(|&m| available(m))
        .unwrap_or(None)
}

/// Developer capture for testing without input: with `EM_CAPTURE=<file.png>` the game
/// opens level `EM_LEVEL` (1-based) or the level file `EM_LEVEL_FILE`, enters the sandbox if `EM_SANDBOX=1`, shows the
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
}

/// A map by its capture name (`EM_MAP`).
fn map_named(m: &str) -> Option<MapMode> {
    match m {
        "magnetic" => Some(MapMode::Magnetic),
        "waves" => Some(MapMode::Waves),
        "particle" => Some(MapMode::ParticleField),
        "total" => Some(MapMode::Total),
        "electric" => Some(MapMode::Electric),
        "energy" => Some(MapMode::EnergyFlow),
        "off" => None,
        _ => Some(MapMode::Potential),
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

#[allow(clippy::too_many_arguments)]
fn input(
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
        // Tube levels follow the animation time (interpolated pictures between the frames).
        if game.editor.level.is_tube() {
            if game.animate {
                game.anim_time.to_bits()
            } else {
                u64::try_from(game.tube_frame().unwrap_or(0)).unwrap_or(u64::MAX)
            }
        } else {
            0
        },
        (
            game.tube_gain_decades.to_bits(),
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
        mode @ (MapMode::Potential | MapMode::Magnetic | MapMode::Electric | MapMode::EnergyFlow),
    ) = game.map
    else {
        return;
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
        let shown = game.tube_frame();
        let picture = game.tube_display_frame();
        let field = match (&level.tube, &game.tube_preview, shown, &picture) {
            (Some(spec), Some(v), Some(k), Some(picture)) if v.revision == game.sent_revision => 'arm: {
                let run = &v.run;
                // Many particles merged for the map (the shader and the arrows sum them all).
                let binned = potential::map_frame(picture);
                // The electromagnetic field's maps (`em.rs`).
                if matches!(mode, MapMode::Magnetic | MapMode::EnergyFlow) {
                    let Some(snap) = game.em_snapshot.as_ref() else {
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
            Response::Invalid { revision, issues } if revision == current => {
                game.setup_issues = Some((revision, issues));
            }
            Response::TubeChunk {
                revision,
                chunk,
                done,
            } if revision == current => {
                // The field's thread follows the same run.
                if game.em_run == Some(revision)
                    && let Some(w) = &game.em
                {
                    w.send(em::EmRequest::Frames {
                        revision,
                        chunk: chunk.clone(),
                    });
                }
                let chunk = std::sync::Arc::try_unwrap(chunk).unwrap_or_else(|a| (*a).clone());
                match &mut game.tube_preview {
                    Some(v) if v.revision == revision => {
                        v.run.extend(chunk);
                        v.done = done;
                    }
                    _ => {
                        game.tube_preview = Some(TubeView {
                            revision,
                            run: chunk,
                            done,
                        });
                    }
                }
            }
            Response::TubeVerified { revision, verdict } if revision == current => {
                game.tube_verdict = Some((revision, verdict));
            }
            _ => {}
        }
    }
}

/// Arrow spacing of the tube's electric field map (cells).
const TUBE_ARROW_SPACING: usize = 2;

/// The most arrows across a tube's electric map (larger boards space them wider).
const TUBE_ARROWS: usize = 32;

/// The energy flow map's arrows: the Poynting vector of the field's snapshot at the
/// arrow grid's points (the nearest sample).
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
fn tube_flow_arrows(level: &Level, snap: &em::EmSnapshot) -> TubeArrows {
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

/// Tube levels' field maps: keeps the field's thread on the preview and the time shown,
/// and takes its snapshots.
fn update_em(mut game: ResMut<Game>) {
    let game = &mut *game;
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

/// Arrows `(position, E)` of a tube's electric field map.
pub type TubeArrows = Vec<(Vec2, Vec2)>;

/// The scale of a tube's electric field map, and its arrows: E at every
/// `TUBE_ARROW_SPACING`-th node. The scale is the 99th percentile of |E| over the nodes
/// (not the maximum: the singular corners would wash out the colours; the 90th
/// saturated most of the board).
fn tube_field(
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
    if tube
        && std::env::var_os("EM_TIME").is_none()
        && matches!(game.map, Some(MapMode::Magnetic | MapMode::EnergyFlow))
        && let Some(s) = &game.em_snapshot
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

/// Moves a zoomed level's view by `shift` (board units), its centre kept on the board.
fn pan_by(game: &mut Game, shift: Vec2) {
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
