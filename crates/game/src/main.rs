//! Electromagnetism – the game. Stage 1: 2D levels (a slice of the 3D world).

#![allow(clippy::disallowed_methods)] // pictures only: the flights come from `physics`

mod average;
mod curriculum;
mod draw;
mod editor;
mod em;
mod inspect;
mod legend;
mod level_docs;
mod level_editor;
mod math;
mod potential;
mod probe;
mod radiation;
mod sandbox;
mod scales;
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

mod camera;
mod capture;
mod keys;
mod maps;
mod metal_maps;
mod physics_sync;
mod tube_views;
use camera::{fit_camera, pan_by};
use capture::dev_capture;
use keys::input;
use maps::{update_field_lines, update_map};
pub use metal_maps::metal_map_span;
use metal_maps::{faded, metal_map_grid, metal_percentile};
use physics_sync::{poll_physics, sync_physics};
use tube_views::{animate, tube_field, tube_flow_arrows, update_em};

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

/// The metal maps (`level::tube::TubeRun::metal_maps`) as last drawn: the setup and the
/// level, the colour units (|j|'s and the power's, the largest so far on this level: a cut
/// core shows dimmer than the whole one it was) and the maps, their values with their
/// standard errors.
pub struct MetalMap {
    pub revision: u64,
    pub level: usize,
    pub units: (f64, f64),
    pub maps: level::tube::MetalMaps,
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
    /// The maps' time average (a meter's time constant; 0: the instantaneous field) and
    /// its state (`average.rs`).
    pub map_tau: f64,
    pub map_average: average::MapAverage,
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
    /// A cut being drawn: the node it started from (the left button held).
    pub cut_start: Option<level::Node>,
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
    /// The board point under the pointer (None over the panel, off the window, or while
    /// something is dragged): what the tooltip inspects (`inspect.rs`).
    pub hover_world: Option<Vec2>,
    /// The board probe is on (I held, or `EM_PROBE=1` in a capture): the tooltip shows
    /// everything known at the pointer (`probe.rs`).
    pub probe: bool,
    /// Particle levels: the static field the probe samples, with the setup revision it was
    /// built for (built when the probe is first used after a change).
    pub probe_field: Option<(u64, physics::field::LevelField)>,
    /// Tube levels: how each electrode's potential swings in the preview run (`Swing`),
    /// for the glow of the oscillating ones.
    pub swing: Option<Swing>,
    /// Tube levels with metals: the carriers traced with their paths (`Traced`) and each
    /// metal's carriers' mean displacement along it (`Drift`, the ticks across it).
    pub traced: Option<Traced>,
    pub drift: Option<Drift>,
    /// The metal maps' values at the time shown, for the tooltip (`MetalMap`).
    pub metal_map: Option<MetalMap>,
    /// Tube levels with metals in a field: each strip's Hall diagram at the time shown
    /// (`update_hall`).
    pub hall: Vec<level::tube::HallForces>,
    /// How the panel groups the electrodes' potential plots, and the electrodes whose lines
    /// it hides (the level's; cleared with it).
    pub plot_layout: ui::PlotLayout,
    pub plot_hidden: Vec<usize>,
    /// An electrode outlined on the board in a plot line's colour (its chip hovered),
    /// until then.
    pub plot_focus: Option<(usize, [u8; 3], std::time::Instant)>,
}

/// The carriers traced in a run's metals, by identity, constant in number: one a contact
/// absorbs is replaced by the newest carrier in the metals (just given off by a contact).
/// Started afresh when the run or the playback starts again.
#[derive(Clone, Debug, Default)]
pub struct Traced {
    pub revision: u64,
    pub t: f64,
    pub ids: Vec<u32>,
}

/// Each metal's carriers' mean displacement (x, y) from the run's start, at every frame so
/// far (the run's revision and the frames counted).
#[derive(Clone, Debug, Default)]
pub struct Drift {
    pub revision: u64,
    pub x: Vec<Vec<[f64; 2]>>,
}

/// Tube levels: each electrode's mean potential over the run's frames so far and its
/// largest departure from it, and the largest of all (the run's revision and the frames
/// counted).
#[derive(Clone, Debug, Default)]
pub struct Swing {
    pub revision: u64,
    pub frames: usize,
    pub mean: Vec<f64>,
    pub own: Vec<f64>,
    pub largest: f64,
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

    /// Outlines electrode `i` in a plot line's colour for `seconds` from now.
    pub fn focus_electrode(&mut self, i: usize, rgb: [u8; 3], seconds: f64) {
        let until = std::time::Instant::now() + std::time::Duration::from_secs_f64(seconds);
        self.plot_focus = Some((i, rgb, until));
    }

    /// A tube level's uniform field B along z at the time shown, with its largest size over
    /// the run (a swinging field's peak): the run's, else the player's solenoid supply or the
    /// level's own field. None without a field.
    pub fn tube_field_now(&self) -> Option<(f64, f64)> {
        let level = &self.editor.level;
        let spec = level.tube.as_ref()?;
        if let Some(v) = self.tube_preview.as_ref().filter(|v| v.run.has_field()) {
            let t = if self.animate {
                self.anim_time
            } else {
                v.run.frames.last().map_or(0.0, |f| f.t)
            };
            let peak = v.run.b_wave.map_or(v.run.b_z.abs(), |w| w.peak());
            return Some((v.run.b_at(t), peak));
        }
        let b = spec
            .solenoid
            .and_then(|n| {
                self.editor
                    .placement
                    .iter()
                    .find(|e| e.kind == level::ElementKind::Supply && e.node == n)
            })
            .map_or(spec.b_z, |e| e.value);
        (b != 0.0 && spec.solenoid_drive.is_none()).then_some((b, b.abs()))
    }

    /// The electrode a plot's chip outlines now, with its colour.
    pub fn focused_electrode(&self) -> Option<(usize, [u8; 3])> {
        self.plot_focus
            .filter(|(_, _, t)| std::time::Instant::now() < *t)
            .map(|(i, c, _)| (i, c))
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
            map_tau: 0.0,
            map_average: average::MapAverage::default(),
            animate: true,
            playback_speed: 1.0,
            anim_time: 0.0,
            last_cursor_world: None,
            sandbox: sandbox::Sandbox::default(),
            mouse_drag: false,
            arrow_drag: None,
            cut_start: None,
            arrow_measure,
            text_focus: false,
            panel_left_px: None,
            zoom: 1.0,
            zoom_center: None,
            pan_anchor: None,
            hover_world: None,
            probe: false,
            probe_field: None,
            swing: None,
            traced: None,
            drift: None,
            metal_map: None,
            hall: Vec::new(),
            plot_layout: ui::PlotLayout::default(),
            plot_hidden: Vec::new(),
            plot_focus: None,
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
        // The previous level's electromagnetic field (a tube's).
        self.em_run = None;
        self.em_snapshot = None;
        self.active_shot = 0;
        self.active_disturbance = 0;
        self.anim_time = 0.0;
        self.plot_hidden.clear();
        self.plot_focus = None;
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

/// Whether a level offers a map (its panel's choices). Tube levels: their potential,
/// magnetic, electric and energy-flow maps (the other views are the 3D slice's); those
/// two only there.
pub fn map_available(level: &Level, map: Option<MapMode>) -> bool {
    match map {
        // The metal maps where there is metal.
        Some(MapMode::Current | MapMode::Heat) if level.is_tube() => {
            level.tube.as_ref().is_some_and(|t| !t.metals.is_empty())
        }
        _ if level.is_tube() => matches!(
            map,
            Some(MapMode::Potential | MapMode::Electric | MapMode::Magnetic | MapMode::EnergyFlow)
                | None
        ),
        Some(MapMode::Waves) => radiation::waves_available(level),
        Some(MapMode::ParticleField) => radiation::particle_field_available(level),
        Some(MapMode::Electric | MapMode::EnergyFlow | MapMode::Current | MapMode::Heat) => false,
        _ => true,
    }
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
        Some(MapMode::Current),
        Some(MapMode::Heat),
        None,
    ];
    let i = order.iter().position(|&m| m == map).unwrap_or(0);
    (1..=order.len())
        .map(|k| order[(i + k) % order.len()])
        .find(|&m| map_available(level, m))
        .unwrap_or(None)
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
        "current" => Some(MapMode::Current),
        "heat" => Some(MapMode::Heat),
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

/// Arrow spacing of the tube's electric field map (cells).
const TUBE_ARROW_SPACING: usize = 2;

/// The most arrows across a tube's electric map (larger boards space them wider).
const TUBE_ARROWS: usize = 32;

/// The frames `update_drift` takes on per display frame (a run's frames hold at most
/// `level::tube::FRAME_BUDGET` carriers in all: a run that arrives at once is caught up
/// with over a few display frames, not in one).
const DRIFT_FRAMES_PER_UPDATE: usize = 20;

/// Arrows `(position, E)` of a tube's electric field map.
pub type TubeArrows = Vec<(Vec2, Vec2)>;
