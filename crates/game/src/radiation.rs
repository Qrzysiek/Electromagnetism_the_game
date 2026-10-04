//! Time-dependent field views (visual only): the fields of antennas and plane waves
//! ("Waves"), the Liénard–Wiechert field of the particle itself, or of every particle of
//! a beam ("Particle field"), and everything at once ("Total": static sources, antennas,
//! waves, disturbances and the particles), at the animation time. For a beam at finite c
//! with the quasi-static interaction, the view can show only what that model leaves out
//! of the dynamics: the full retarded field minus the quasi-static fields of the
//! particles' present states (PHYSICS.md §3.3).
//!
//! The map is drawn on the GPU (`radiation.wgsl`), every pixel at every frame: the
//! analytic antenna and wave fields, and the retarded fields of the charges from their
//! world-line samples. The static part of the total field (charges, magnets, coils,
//! metal, uniform stray fields) does not change in time: it is computed once per setup on
//! the CPU in f64 with the tested physics code and handed to the shader on a grid. The
//! colour scales and the E arrows are computed on the CPU (`sample`, `charges`).

use bevy::asset::{AssetPath, embedded_asset, embedded_path};
use bevy::prelude::*;
use bevy::render::render_resource::{AsBindGroup, ShaderType};
use bevy::render::storage::ShaderBuffer;
use bevy::shader::ShaderRef;
use bevy::sprite_render::{AlphaMode2d, Material2d, Material2dPlugin};
use level::Level;
use physics::DVec3;
use physics::dynamics::Kinematics;
use physics::external::External;
use physics::field::{FieldSolver, LevelField};
use physics::geometry::Shape;
use physics::lienard::{self, SampledWorldline};
use rayon::prelude::*;

use crate::Game;
use crate::potential::MapMode;
use crate::worker::Preview;

/// Grid points per cell of the static part of the total field.
const STATIC_PER_CELL: f64 = 6.0;
/// Spacing of the E arrows, in cells.
const ARROW_SPACING: f64 = 1.5;
/// Quantity used for the colour of the field views.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldQuantity {
    /// B perpendicular to the plane (all of B in the plane), signed.
    Bz,
    /// Magnitude of E.
    E,
    /// The flow of field energy, the Poynting vector `S = (c²/4π) E × B` (Jackson §6.7):
    /// in the plane `(c²/4π) B_z (E_y, −E_x)`. Its magnitude colours the map; arrows and
    /// tracers show its direction. Finite c only.
    S,
}

/// Which part of the energy flow the map shows. The view's field is split into the moving
/// charges' field and the rest (static sources, antennas, waves: what the view shows
/// besides the charges); the flow into each part's own terms and the exchange terms
/// (`physics::poynting`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum FlowPart {
    #[default]
    Total,
    /// The charges' field alone.
    Own,
    /// The cross terms of the charges' field with the rest: they carry the work the rest
    /// does on the charges.
    Exchange,
    /// The rest alone.
    External,
}

impl FlowPart {
    pub const ALL: [FlowPart; 4] = [
        FlowPart::Total,
        FlowPart::Own,
        FlowPart::Exchange,
        FlowPart::External,
    ];

    fn index(self) -> usize {
        match self {
            FlowPart::Total => 0,
            FlowPart::Own => 1,
            FlowPart::Exchange => 2,
            FlowPart::External => 3,
        }
    }
}

/// An energy tracer: a dot drifting with the energy's velocity `S/u`.
#[derive(Clone, Copy, Debug)]
struct Tracer {
    x: DVec3,
    /// Seconds (real time) since it appeared.
    age: f32,
}

/// Number of energy tracers, and how long each lives (s, real time).
const TRACERS: usize = 400;
const TRACER_LIFE: f32 = 4.0;

/// A moving charge of the field views: its world line, its charge, and when and where it
/// was absorbed (as in the beam dynamics: its field then disappears where the light cone
/// of that moment passes, or, in a screening cup, fades as it flies on).
struct Source {
    line: SampledWorldline,
    /// The samples `(t, x, v, a)` of the world line (for the GPU).
    samples: Vec<(f64, DVec3, DVec3, DVec3)>,
    charge: f64,
    mass: f64,
    /// Magnetic moment along z (its dipole field B_z = −m/r³ in the plane).
    moment: f64,
    end: Option<(f64, DVec3)>,
    /// Whether its charge stays where it was absorbed (a body; `Fate::Stop`) rather than
    /// being drained (the detector).
    stays: bool,
    /// Its fade in a screening cup (`Fate::Cup`): the world line goes on uniformly into
    /// the cup and the charge seen from outside fades, at the retarded time.
    fade: Option<physics::beam::Fade>,
}

#[derive(ShaderType, Debug, Clone, Copy, Default)]
pub struct FieldParams {
    /// min.x, min.y, size.x, size.y of the arena.
    pub area: Vec4,
    /// c (0 for infinite), B saturation, E saturation, dynamic range.
    pub scales: Vec4,
    /// Static grid width, height, 1 if present; number of charges.
    pub grid: UVec4,
    /// Antennas, waves, flags (1 radiation only, 2 left out by the model, 4 colour |E|,
    /// 8 linear, 16 colour |S|, 32 colour of the charges' field alone), energy flow: part
    /// (bits 0–1: total, own, exchange, external), 4 averaged over a period, number of
    /// frequency groups (bits 8–15).
    pub counts: UVec4,
    /// Energy flow: saturation of the part shown, 0, 0, 0.
    pub flow: Vec4,
}

#[derive(Asset, TypePath, AsBindGroup, Debug, Clone)]
pub struct FieldMaterial {
    #[uniform(0)]
    pub params: FieldParams,
    #[storage(1, read_only)]
    pub samples: Handle<ShaderBuffer>,
    #[storage(2, read_only)]
    pub items: Handle<ShaderBuffer>,
    #[storage(3, read_only)]
    pub statics: Handle<ShaderBuffer>,
}

impl Material2d for FieldMaterial {
    fn fragment_shader() -> ShaderRef {
        ShaderRef::Path(
            AssetPath::from_path_buf(embedded_path!("radiation.wgsl")).with_source("embedded"),
        )
    }

    fn alpha_mode(&self) -> AlphaMode2d {
        AlphaMode2d::Blend
    }
}

pub struct FieldViewPlugin;

impl Plugin for FieldViewPlugin {
    fn build(&self, app: &mut App) {
        embedded_asset!(app, "radiation.wgsl");
        app.add_plugins(Material2dPlugin::<FieldMaterial>::default());
    }
}

/// (revision, flight, mode, flight time of the preview, radiation only, neglected only,
/// the verdict's own flight shown, energy flow shown, flow averaged): what the view's
/// scales and world lines were computed for.
type ViewKey = (u64, usize, MapMode, u64, bool, bool, bool, bool, bool);

#[derive(Resource)]
pub struct RadiationView {
    material: Handle<FieldMaterial>,
    entity: Entity,
    key: Option<ViewKey>,
    radiation_only: bool,
    /// Show only the part of the beam's field the quasi-static interaction leaves out.
    neglected_only: bool,
    /// The moving charges shown: one particle, or every particle of a beam.
    sources: Vec<Source>,
    b_sat: f64,
    e_sat: f64,
    /// Saturation (99th percentile) of each part of the energy flow |S| (`FlowPart`
    /// order), for the averaging chosen.
    s_sat: [f64; 4],
    /// Largest (99th percentile) B_z and |E| of the part shown (radiation part, or left
    /// out by the model) relative to the full field's scale; 1 for the full field.
    pub part_b: f64,
    pub part_e: f64,
    /// The static part of the view's field, and its oscillating sources by frequency
    /// (antennas and plane waves of the same ω together): for the energy flow averaged
    /// over a period.
    still: LevelField,
    groups: Vec<(f64, LevelField)>,
    /// The energy flow averaged over a period of the oscillating sources.
    average: bool,
    /// Time and style (quantity, range, arrows, flow part) of the current arrows.
    time: f64,
    style: (FieldQuantity, u64, bool, FlowPart),
    tracers: Vec<Tracer>,
    /// State of the tracers' random placement (SplitMix64).
    tracer_seed: u64,
    /// The level's speed of light (for the tracers of the exchange flow).
    c: f64,
    /// Obstacles of the shown flight (the field is not drawn inside sources).
    obstacles: Vec<Shape>,
    /// The flight's whole field (what the particles feel: for the quasi-static
    /// interaction's continued pasts, `source_fields`).
    source_field: LevelField,
    pub arrows: Vec<(Vec2, Vec2)>,
}

impl RadiationView {
    /// Value of full colour for B_z (the full field's 99th percentile).
    pub fn b_sat(&self) -> f64 {
        self.b_sat
    }

    /// Value of full colour for |E|.
    pub fn e_sat(&self) -> f64 {
        self.e_sat
    }

    /// Value of full colour for |S| of a part of the energy flow.
    pub fn s_sat(&self, part: FlowPart) -> f64 {
        self.s_sat[part.index()]
    }

    /// Whether the view's field has oscillating sources (antennas, plane waves), whose
    /// energy flow can be averaged over a period.
    pub fn oscillates(&self) -> bool {
        !self.groups.is_empty()
    }
}

/// A storage buffer must not be empty.
fn buffer(v: Vec<[f32; 4]>) -> ShaderBuffer {
    ShaderBuffer::from(if v.is_empty() { vec![[0.0f32; 4]] } else { v })
}

pub fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<FieldMaterial>>,
    mut buffers: ResMut<Assets<ShaderBuffer>>,
) {
    let material = materials.add(FieldMaterial {
        params: FieldParams::default(),
        samples: buffers.add(buffer(Vec::new())),
        items: buffers.add(buffer(Vec::new())),
        statics: buffers.add(buffer(Vec::new())),
    });
    let entity = commands
        .spawn((
            Mesh2d(meshes.add(Rectangle::new(1.0, 1.0))),
            MeshMaterial2d(material.clone()),
            Transform::from_xyz(0.0, 0.0, -9.0),
            Visibility::Hidden,
        ))
        .id();
    commands.insert_resource(RadiationView {
        material,
        entity,
        key: None,
        radiation_only: false,
        neglected_only: false,
        sources: Vec::new(),
        b_sat: 1.0,
        e_sat: 1.0,
        s_sat: [1.0; 4],
        part_b: 1.0,
        part_e: 1.0,
        still: LevelField::default(),
        groups: Vec::new(),
        average: true,
        time: f64::NAN,
        style: (FieldQuantity::Bz, 0, false, FlowPart::Total),
        tracers: Vec::new(),
        tracer_seed: 0x5EED,
        c: f64::INFINITY,
        obstacles: Vec::new(),
        source_field: LevelField::default(),
        arrows: Vec::new(),
    });
}

/// Whether the level has time-dependent sources to show.
pub fn waves_available(level: &Level) -> bool {
    level.c().is_finite()
        && (level.physics.rf_omega != 0.0
            || level
                .disturbances
                .iter()
                .any(|d| d.waves.iter().any(|w| w.omega != 0.0)))
}

/// The particle's own field needs a finite speed of light (retardation).
pub fn particle_field_available(level: &Level) -> bool {
    // Retarded for finite c, instantaneous (exact) for c = ∞; charges and moments.
    level
        .shots
        .iter()
        .any(|s| s.particle.charge != 0.0 || s.particle.moment != 0.0)
}

/// Only the time-dependent sources of a level field: antennas and waves.
fn time_dependent(field: &LevelField) -> LevelField {
    LevelField {
        antennas: field.antennas.clone(),
        external: field
            .external
            .iter()
            .copied()
            .filter(|e| matches!(e, External::Wave(w) if w.omega != 0.0))
            .collect(),
        time_offset: 0.0,
        ..LevelField::default()
    }
}

/// World line of a computed flight: positions, velocities and accelerations (from the
/// Lorentz force) at the preview's path points.
#[allow(clippy::type_complexity)]
fn worldline(
    preview: &Preview,
    mass: f64,
    c: f64,
) -> Option<(SampledWorldline, Vec<(f64, DVec3, DVec3, DVec3)>)> {
    let kin = Kinematics::new(mass, c);
    let mut samples: Vec<(f64, DVec3, DVec3, DVec3)> = Vec::new();
    for p in &preview.path {
        if samples.last().is_some_and(|s| p.t <= s.0) {
            continue;
        }
        let gamma = kin.gamma(p.p);
        let v = kin.velocity(p.p);
        let a = (p.force - v * (v.dot(p.force) / (c * c))) / (gamma * mass);
        samples.push((p.t, p.x, v, a));
    }
    (!samples.is_empty()).then(|| (SampledWorldline::new(&samples), samples))
}

/// A beam particle's world line with its motion before launch prepended: its launch
/// acceleration continued back as the beam dynamics does (`physics::beam::tapered`:
/// constant, then fading so that the velocity changes by at most 0.1c; PHYSICS.md §3.3)
/// over the light time across the arena, so that no switch-on shell appears.
fn with_past(
    line: &[(f64, DVec3, DVec3, DVec3)],
    bounds: &physics::geometry::Aabb,
    c: f64,
) -> Vec<(f64, DVec3, DVec3, DVec3)> {
    let Some(&(t0, x0, v0, a0)) = line.first() else {
        return Vec::new();
    };
    let span = (bounds.max - bounds.min).length() / c * 1.2;
    // Beyond 3.7 of the taper's time scale the acceleration has faded (sech² < 1e-8):
    // uniform motion, the world line's own continuation before its first sample.
    let span = span.min(3.7 * 0.1 * c / a0.length().max(1e-300));
    if !(span > 0.0 && span.is_finite()) {
        return line.to_vec();
    }
    let mut out: Vec<_> = (1..=64)
        .rev()
        .map(|k| {
            let dt = -span * f64::from(k) / 64.0;
            let (x, v, a) = physics::beam::tapered(x0, v0, a0, dt, c);
            (t0 + dt, x, v, a)
        })
        .collect();
    out.extend_from_slice(line);
    out
}

/// 99th percentile of a list of magnitudes (ignoring non-finite values).
/// 0 when the quantity vanishes everywhere (e.g. no magnetic field at all): nothing to
/// colour. (It was 1e-300: in the shader's f32 that is 0, and 0/0 painted the whole map.)
/// Where the field is confined to a small area the 99th percentile is 0: then the largest
/// value is used.
fn percentile99(mut v: Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    let i = (v.len() - 1) * 99 / 100;
    let p = if v[i] > 0.0 { v[i] } else { v[v.len() - 1] };
    if p > 1e-30 { p } else { 0.0 }
}

/// Signed compression into [−1, 1] on an asinh scale: values from `sat / range` to `sat`
/// are distinguishable (`range` = 10^decades).
fn compress(v: f64, sat: f64, range: f64) -> f64 {
    let r = sat / range;
    ((v / r).asinh() / (sat / r).asinh()).clamp(-1.0, 1.0)
}

/// The colour scale: logarithmic over a range, or linear with a gain (`scale` is the
/// range or the gain). Signed, in [−1, 1].
pub fn colour_value(v: f64, sat: f64, linear: bool, scale: f64) -> f64 {
    if sat <= 0.0 {
        // The quantity vanishes everywhere.
        return 0.0;
    }
    if linear {
        (v * scale / sat).clamp(-1.0, 1.0)
    } else {
        compress(v, sat, scale)
    }
}

#[allow(
    clippy::too_many_lines,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub fn update(
    game: Res<Game>,
    time: Res<Time>,
    mut view: ResMut<RadiationView>,
    mut materials: ResMut<Assets<FieldMaterial>>,
    mut buffers: ResMut<Assets<ShaderBuffer>>,
    mut quads: Query<(&mut Visibility, &mut Transform)>,
) {
    let mode = match game.map {
        Some(m @ (MapMode::Waves | MapMode::ParticleField | MapMode::Total)) => m,
        _ => {
            if let Ok((mut v, _)) = quads.get_mut(view.entity) {
                *v = Visibility::Hidden;
            }
            view.arrows.clear();
            view.tracers.clear();
            return;
        }
    };
    let level = &game.editor.level;
    if level.shots.is_empty() || (!level.c().is_finite() && mode == MapMode::Waves) {
        return;
    }
    let c = level.c();
    let flight = game.active_flight();
    let (shot, _) = level.flight_of(flight);
    let preview = game.flights.get(flight).and_then(|f| f.preview.as_ref());
    // A beam: every particle's world line, from the active disturbance's flight (the
    // verdict's own flight once it has arrived: exact at finite c).
    let d_beam = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    let beam = level
        .has_beams()
        .then(|| game.beams.get(d_beam).and_then(crate::BeamView::shown));
    let verdict_shown = game
        .beams
        .get(d_beam)
        .is_some_and(|b| b.verdict_flight.is_some());
    let t_final = match beam {
        Some(b) => b.map_or(0.0, |b| {
            b.worldlines
                .iter()
                .filter_map(|l| l.last().map(|s| s.0))
                .fold(0.0, f64::max)
        }),
        None => preview.map_or(0.0, |p| p.flight_time),
    };
    let neglected_only = game.neglected_only && beam.is_some();
    let key = (
        game.sent_revision,
        flight,
        mode,
        t_final.to_bits(),
        game.radiation_only,
        neglected_only,
        verdict_shown,
        game.field_quantity == FieldQuantity::S,
        game.flow_average,
    );
    let bounds = level.bounds();
    let size = bounds.max - bounds.min;
    let field = {
        let mut scn = level.scenarios_at(
            &game.editor.placement,
            physics::conductor::Resolution::Display,
        );
        let scn = scn.swap_remove(flight);
        view.obstacles = scn.obstacles;
        view.source_field = scn.field.clone();
        match mode {
            // The full field; its time offset is the shot's launch time, so it is
            // sampled at the flight time.
            MapMode::Total => scn.field,
            _ => time_dependent(&scn.field),
        }
    };
    let launch_time = level.shots[shot].launch.time;
    let Some((samples_h, items_h, statics_h)) = materials
        .get(&view.material)
        .map(|m| (m.samples.clone(), m.items.clone(), m.statics.clone()))
    else {
        return;
    };

    // Once per setup, flight and mode: world lines, colour scales, the static field.
    if view.key != Some(key) {
        view.key = Some(key);
        view.time = f64::NAN;
        view.radiation_only = game.radiation_only;
        view.sources = match beam {
            Some(b) => {
                // Species by position in the flight (dynamic particles follow the shots).
                let species = level.beam_species(&game.editor.placement);
                b.map_or_else(Vec::new, |b| {
                    b.worldlines
                        .iter()
                        .zip(&b.ends)
                        .zip(&species)
                        .zip(&b.outcomes)
                        .zip(&b.fades)
                        .filter(|((((l, _), _), _), _)| !l.is_empty())
                        .map(|((((l, &end), sp), &o), &fade)| {
                            let samples = with_past(l, &bounds, c);
                            Source {
                                line: SampledWorldline::new(&samples),
                                samples,
                                charge: sp.charge,
                                mass: sp.mass,
                                moment: sp.moment,
                                end,
                                stays: matches!(o, physics::trajectory::Outcome::Collided(_)),
                                fade,
                            }
                        })
                        .collect()
                })
            }
            None => preview
                .and_then(|p| {
                    use physics::trajectory::Outcome;
                    let (line, samples) = worldline(p, level.shots[shot].particle.mass, c)?;
                    let last = p.path.last()?;
                    // Absorbed by a body (the charge stays) or by the detector (drained);
                    // otherwise it flies on.
                    let (end, stays) = match p.outcome {
                        Outcome::Collided(_) => (Some((last.t, last.x)), true),
                        Outcome::Arrived | Outcome::Rejected | Outcome::SkippedGate(_) => {
                            (Some((last.t, last.x)), false)
                        }
                        _ => (None, false),
                    };
                    Some(Source {
                        line,
                        samples,
                        charge: level.shots[shot].particle.charge,
                        mass: level.shots[shot].particle.mass,
                        moment: level.shots[shot].particle.moment,
                        end,
                        stays,
                        fade: None,
                    })
                })
                .into_iter()
                .collect(),
        };
        // Scales from the full field (the baseline): the parts (radiation, left out by the
        // model) are shown on the same colour scale, so that their true size is seen.
        view.neglected_only = false;
        view.radiation_only = false;
        let (t0, t1) = match mode {
            MapMode::Waves => {
                let w = field
                    .antennas
                    .iter()
                    .map(|a| a.omega)
                    .chain(field.external.iter().filter_map(|e| match e {
                        External::Wave(w) => Some(w.omega),
                        External::Uniform { .. } => None,
                    }))
                    .fold(0.0, f64::max);
                let period = if w > 0.0 {
                    2.0 * std::f64::consts::PI / w
                } else {
                    1.0
                };
                (launch_time, launch_time + period)
            }
            _ => (0.0, if t_final > 0.0 { t_final } else { 1.0 }),
        };
        // Sampled over 8 times on the cell grid, in parallel.
        let points: Vec<(f64, DVec3)> = (0..8)
            .flat_map(|k| {
                let t = t0 + (t1 - t0) * (f64::from(k) + 0.5) / 8.0;
                let (nx, ny) = (size.x.floor() as u32, size.y.floor() as u32);
                (0..=ny).flat_map(move |j| {
                    (0..=nx).map(move |i| {
                        (
                            t,
                            DVec3::new(
                                bounds.min.x + f64::from(i),
                                bounds.min.y + f64::from(j),
                                0.0,
                            ),
                        )
                    })
                })
            })
            .collect();
        let v = &*view;
        let fr = &field;
        let (bs, es): (Vec<f64>, Vec<f64>) = points
            .par_iter()
            .filter_map(|&(t, x)| sample(v, fr, mode, x, t, c).map(|(e, b)| (b.abs(), e.length())))
            .unzip();
        view.b_sat = percentile99(bs);
        view.e_sat = percentile99(es);
        view.neglected_only = neglected_only;
        view.radiation_only = game.radiation_only;
        // The size of the part shown, relative to the baseline.
        if view.neglected_only || view.radiation_only {
            let v = &*view;
            let (bs, es): (Vec<f64>, Vec<f64>) = points
                .par_iter()
                .filter_map(|&(t, x)| {
                    sample(v, fr, mode, x, t, c).map(|(e, b)| (b.abs(), e.length()))
                })
                .unzip();
            let ratio = |p: f64, sat: f64| if sat > 0.0 { p / sat } else { 0.0 };
            view.part_b = ratio(percentile99(bs), view.b_sat);
            view.part_e = ratio(percentile99(es), view.e_sat);
        } else {
            view.part_b = 1.0;
            view.part_e = 1.0;
        }
        // The energy flow: the static part of the rest (with antennas that do not
        // oscillate) and its oscillating sources by frequency, for the average over a
        // period; the scale of each part, for the flow shown (finite c only).
        view.still = static_part(&field);
        view.still.antennas = field
            .antennas
            .iter()
            .filter(|a| a.omega == 0.0)
            .copied()
            .collect();
        view.groups = frequency_groups(&field);
        view.average = game.flow_average;
        view.c = c;
        view.tracers.clear();
        if game.field_quantity == FieldQuantity::S && c.is_finite() {
            let v = &*view;
            let per_point: Vec<[f64; 4]> = points
                .par_iter()
                .filter_map(|&(t, x)| {
                    flows_at(v, fr, mode, x, t, c).map(|f| f.map(|(s, _)| s.length()))
                })
                .collect();
            for part in FlowPart::ALL {
                view.s_sat[part.index()] =
                    percentile99(per_point.iter().map(|p| p[part.index()]).collect());
            }
        }
        // The static part of the total field, once (it does not change in time).
        // The static part of the total field on a grid; for the particle field only the
        // mask of the bodies (inside metal and other bodies nothing is drawn).
        let grid = if matches!(mode, MapMode::Total | MapMode::ParticleField) {
            let (w, h) = (
                (size.x * STATIC_PER_CELL).ceil().max(2.0) as u32,
                (size.y * STATIC_PER_CELL).ceil().max(2.0) as u32,
            );
            let still = static_part(&field);
            let obstacles = &view.obstacles;
            let data: Vec<[f32; 4]> = (0..h)
                .into_par_iter()
                .flat_map_iter(|j| {
                    let y = bounds.min.y + (f64::from(j) + 0.5) / f64::from(h) * size.y;
                    let still = &still;
                    (0..w).map(move |i| {
                        let x = bounds.min.x + (f64::from(i) + 0.5) / f64::from(w) * size.x;
                        let x = DVec3::new(x, y, 0.0);
                        if obstacles.iter().any(|o| o.signed_distance(x) < 0.0) {
                            return [0.0, 0.0, 0.0, 1.0];
                        }
                        if mode != MapMode::Total {
                            return [0.0; 4];
                        }
                        let f = still.sample(x, 0.0);
                        [f.e.x as f32, f.e.y as f32, f.b.z as f32, 0.0]
                    })
                })
                .collect();
            if let Some(mut b) = buffers.get_mut(&statics_h) {
                *b = buffer(data);
            }
            UVec4::new(w, h, 1, 0)
        } else {
            UVec4::ZERO
        };
        if let Some(mut m) = materials.get_mut(&view.material) {
            m.params.grid = grid;
        }
        if std::env::var("EM_CAPTURE").is_ok() {
            eprintln!(
                "radiation scales: B {:.3e}, E {:.3e}",
                view.b_sat, view.e_sat
            );
        }
    }

    // Time shown: lab time for waves, flight time for the particles' own fields.
    let t_anim = if game.animate {
        game.anim_time
    } else {
        t_final
    };
    let t = match mode {
        MapMode::Waves => launch_time + t_anim,
        _ => t_anim,
    };
    let linear = game.field_linear;
    let range = if linear {
        10f64.powf(game.field_gain_decades)
    } else {
        10f64.powf(game.field_range_decades)
    };
    let quantity = game.field_quantity;

    // GPU inputs for this frame: times relative to t, phases at t.
    let started = std::time::Instant::now();
    let mut samples: Vec<[f32; 4]> = Vec::new();
    let mut items: Vec<[f32; 4]> = Vec::new();
    for (i, s) in view.sources.iter().enumerate() {
        items.push([
            (samples.len() / 2) as f32,
            s.samples.len() as f32,
            s.charge as f32,
            if s.end.is_some() { 1.0 } else { 0.0 },
        ]);
        let (te, xe) = s.end.unwrap_or((f64::INFINITY, DVec3::ZERO));
        items.push([
            (te - t).min(1e30) as f32,
            xe.x as f32,
            xe.y as f32,
            if s.stays { 1.0 } else { 0.0 },
        ]);
        items.push([
            s.moment as f32,
            s.fade.map_or(0.0, |f| f.rate as f32),
            0.0,
            0.0,
        ]);
        // The quasi-static interaction's continued past (`beam::continued_fields`): its
        // motion in the fields it feels now, `(U₀, ω²)`, `(Λ U₀, valid)`, `(Λ² U₀, 0)`.
        let motion = (view.neglected_only && s.moment == 0.0)
            .then(|| source_fields(&view, i, t, c))
            .flatten()
            .and_then(|(r, v, _, e, b)| {
                physics::field_motion::FieldMotion::new(s.charge, s.mass, c, 0.0, r, v, e, b)
            });
        match motion.map(|m| m.parameters()) {
            Some((u0, v1, v2, w2)) => {
                items.push([u0[0] as f32, u0[1] as f32, u0[2] as f32, w2 as f32]);
                items.push([v1[0] as f32, v1[1] as f32, v1[2] as f32, 1.0]);
                items.push([v2[0] as f32, v2[1] as f32, v2[2] as f32, 0.0]);
            }
            None => items.extend([[0.0; 4]; 3]),
        }
        for &(ts, x, v, a) in &s.samples {
            samples.push([(ts - t) as f32, x.x as f32, x.y as f32, v.x as f32]);
            samples.push([v.y as f32, a.x as f32, a.y as f32, 0.0]);
        }
    }
    // Phases reduced in f64 (f32 would lose them over long times).
    let tau = std::f64::consts::TAU;
    let t_field = t + field.time_offset;
    for a in &field.antennas {
        let ph = (a.omega * t_field + a.phase).rem_euclid(tau);
        items.push([
            a.position.x as f32,
            a.position.y as f32,
            a.amplitude.x as f32,
            a.amplitude.y as f32,
        ]);
        // Its frequency group, for the flow averaged over a period (static: 255).
        let group = if a.omega == 0.0 {
            STATIC_GROUP
        } else {
            group_index(&view.groups, a.omega)
        };
        items.push([a.omega as f32, ph as f32, a.radius as f32, group as f32]);
    }
    let mut n_waves = 0u32;
    for e in &field.external {
        if let External::Wave(w) = e
            && w.omega != 0.0
        {
            let ph = (w.omega * t_field + w.phase).rem_euclid(tau);
            items.push([
                w.direction.x as f32,
                w.direction.y as f32,
                w.polarization.x as f32,
                w.polarization.y as f32,
            ]);
            items.push([
                w.amplitude as f32,
                w.omega as f32,
                ph as f32,
                group_index(&view.groups, w.omega) as f32,
            ]);
            n_waves += 1;
        }
    }
    if let Some(mut b) = buffers.get_mut(&samples_h) {
        *b = buffer(samples);
    }
    if let Some(mut b) = buffers.get_mut(&items_h) {
        *b = buffer(items);
    }
    let flow = quantity == FieldQuantity::S && c.is_finite();
    // In the particle-field view the colour (B_z or |E|) is the moving charges' field
    // alone, as its arrows and scale are; the antennas and waves are uploaded only as the
    // rest of the energy flow's exchange and external parts. (Until 2026-10 they were
    // added to the colour too, so a level's light wave filled the view, saturated on the
    // particle's far smaller scale: level 84.)
    let charges_only = mode == MapMode::ParticleField;
    let flags = u32::from(view.radiation_only)
        | (u32::from(view.neglected_only) << 1)
        | (u32::from(quantity == FieldQuantity::E) << 2)
        | (u32::from(linear) << 3)
        | (u32::from(flow) << 4)
        | (u32::from(charges_only) << 5);
    let part = game.flow_part;
    let averaged = view.average && !view.groups.is_empty();
    let flow_bits = part.index() as u32
        | (u32::from(averaged) << 2)
        | ((view.groups.len().min(254) as u32) << 8);
    if let Some(mut m) = materials.get_mut(&view.material) {
        m.params.area = Vec4::new(
            bounds.min.x as f32,
            bounds.min.y as f32,
            size.x as f32,
            size.y as f32,
        );
        m.params.scales = Vec4::new(
            if c.is_finite() { c as f32 } else { 0.0 },
            view.b_sat as f32,
            view.e_sat as f32,
            range as f32,
        );
        m.params.grid.w = view.sources.len() as u32;
        m.params.counts = UVec4::new(field.antennas.len() as u32, n_waves, flags, flow_bits);
        // The flow's saturation, with the constant c²/4π (the GPU forms E × B_z).
        m.params.flow = Vec4::new(
            (view.s_sat[part.index()] * 4.0 * std::f64::consts::PI / (c * c)) as f32,
            0.0,
            0.0,
            0.0,
        );
    }
    if let Ok((mut vis, mut tr)) = quads.get_mut(view.entity) {
        *vis = Visibility::Visible;
        let center = (bounds.max + bounds.min) * 0.5;
        tr.translation = Vec3::new(center.x as f32, center.y as f32, -9.0);
        tr.scale = Vec3::new(size.x as f32, size.y as f32, 1.0);
    }

    // Energy tracers: every frame, with the field at the time shown, at the animation's
    // pace (4 time units per second at speed 1, as `animate`).
    let real_dt = time.delta_secs();
    advance_tracers(
        &mut view,
        &field,
        mode,
        &bounds,
        t,
        part,
        flow && game.flow_tracers,
        f64::from(real_dt) * 4.0 * game.playback_speed,
        real_dt,
        linear,
        range,
    );

    // E arrows, or S arrows for the energy flow (CPU): only when the time or the style
    // changed.
    let style = (
        game.field_quantity,
        game.field_range_decades.to_bits()
            ^ game.field_gain_decades.to_bits().rotate_left(1)
            ^ u64::from(linear),
        game.show_field_arrows,
        part,
    );
    if t.to_bits() == view.time.to_bits() && style == view.style {
        return;
    }
    view.time = t;
    view.style = style;
    let mut arrows = Vec::new();
    if game.show_field_arrows {
        let mut points = Vec::new();
        let mut y = bounds.min.y + 0.5 * ARROW_SPACING;
        while y < bounds.max.y {
            let mut x = bounds.min.x + 0.5 * ARROW_SPACING;
            while x < bounds.max.x {
                points.push(DVec3::new(x, y, 0.0));
                x += ARROW_SPACING;
            }
            y += ARROW_SPACING;
        }
        let v = &*view;
        let fr = &field;
        arrows = points
            .par_iter()
            .filter_map(|&x| {
                let (a, sat) = if flow {
                    (
                        flow_at(v, fr, mode, x, t, c, part)?.0,
                        v.s_sat[part.index()],
                    )
                } else {
                    (sample(v, fr, mode, x, t, c)?.0, v.e_sat)
                };
                let len =
                    (colour_value(a.length(), sat, linear, range) * 0.9 * ARROW_SPACING) as f32;
                let d = Vec2::new(a.x as f32, a.y as f32).normalize_or_zero();
                (len > 0.05).then_some((Vec2::new(x.x as f32, x.y as f32), d * len))
            })
            .collect();
    }
    view.arrows = arrows;
    if std::env::var("EM_CAPTURE").is_ok() {
        eprintln!(
            "radiation view: {:.2} ms",
            started.elapsed().as_secs_f64() * 1e3
        );
    }
}

/// The static part of a level field: without antennas and oscillating waves.
fn static_part(field: &LevelField) -> LevelField {
    let mut f = field.clone();
    f.antennas.clear();
    f.external
        .retain(|e| !matches!(e, External::Wave(w) if w.omega != 0.0));
    f
}

/// Frequency group index of the GPU's antennas that do not oscillate (ω = 0): static.
const STATIC_GROUP: u32 = 255;

/// The oscillating sources of a level field by frequency: antennas and plane waves of the
/// same ω together, each group a field of its own (for the flow averaged over a period:
/// the cross terms of different frequencies and of the static part average out).
fn frequency_groups(field: &LevelField) -> Vec<(f64, LevelField)> {
    let mut groups: Vec<(f64, LevelField)> = Vec::new();
    let group = |groups: &mut Vec<(f64, LevelField)>, w: f64| -> usize {
        if let Some(k) = groups.iter().position(|g| g.0.to_bits() == w.to_bits()) {
            k
        } else {
            groups.push((
                w,
                LevelField {
                    time_offset: field.time_offset,
                    ..LevelField::default()
                },
            ));
            groups.len() - 1
        }
    };
    for a in field.antennas.iter().filter(|a| a.omega != 0.0) {
        let k = group(&mut groups, a.omega);
        groups[k].1.antennas.push(*a);
    }
    for e in &field.external {
        if let External::Wave(w) = e
            && w.omega != 0.0
        {
            let k = group(&mut groups, w.omega);
            groups[k].1.external.push(*e);
        }
    }
    groups
}

/// The index of an oscillating source's frequency group (as `frequency_groups` orders
/// them), for the GPU.
fn group_index(groups: &[(f64, LevelField)], omega: f64) -> u32 {
    #[allow(clippy::cast_possible_truncation)]
    let k = groups
        .iter()
        .position(|g| g.0.to_bits() == omega.to_bits())
        .unwrap_or(0) as u32;
    k
}

/// The view's field at a point, split into the moving charges' field and the rest (what
/// the view shows besides them: static sources in the total view, antennas and waves).
#[derive(Clone, Copy, Debug)]
struct Split {
    e_rest: DVec3,
    b_rest: DVec3,
    e_charges: DVec3,
    b_charges: DVec3,
}

/// `Split` at `x`, time `t`; `None` inside a source (as `sample`).
fn sample_split(
    view: &RadiationView,
    field: &LevelField,
    mode: MapMode,
    x: DVec3,
    t: f64,
    c: f64,
) -> Option<Split> {
    if field
        .antennas
        .iter()
        .any(|a| (x - a.position).length() < a.radius)
    {
        return None;
    }
    let f = field.sample(x, t);
    let (e, bz) = match mode {
        MapMode::Waves => (DVec3::ZERO, 0.0),
        _ => {
            // Inside bodies nothing is drawn (the GPU's mask of the static grid).
            if view.obstacles.iter().any(|o| o.signed_distance(x) < 0.0) {
                return None;
            }
            charges(view, x, t, c)?
        }
    };
    Some(Split {
        e_rest: f.e,
        b_rest: DVec3::Z * f.b.z,
        e_charges: e,
        b_charges: DVec3::Z * bz,
    })
}

/// The energy flow `S` and energy density `u` of every part (`physics::poynting`, in
/// `FlowPart` order), at `x` and time `t`: all, the charges' own terms, the exchange terms,
/// the rest's own terms. Averaged over a period of the oscillating sources (the charges
/// held at time `t`): the rest's flow and energy are its static part's plus, for each
/// frequency, `½ (S(t) + S(t + T/4))` (exact for one frequency); the cross terms of
/// different frequencies, of the oscillating and the static fields, and of the charges
/// with the oscillating fields average out.
fn flows_at(
    view: &RadiationView,
    field: &LevelField,
    mode: MapMode,
    x: DVec3,
    t: f64,
    c: f64,
) -> Option<[(DVec3, f64); 4]> {
    use physics::poynting::{energy_density, exchange, poynting};
    let s = sample_split(view, field, mode, x, t, c)?;
    let own = (
        poynting(s.e_charges, s.b_charges, c),
        energy_density(s.e_charges, s.b_charges, c),
    );
    let (rest, (ux, sx)) = if view.average && !view.groups.is_empty() {
        let st = view.still.sample(x, t);
        let (es, bs) = (st.e, DVec3::Z * st.b.z);
        let mut rest = (poynting(es, bs, c), energy_density(es, bs, c));
        for (w, g) in &view.groups {
            for dt in [0.0, std::f64::consts::FRAC_PI_2 / w] {
                let f = g.sample(x, t + dt);
                let (e, b) = (f.e, DVec3::Z * f.b.z);
                rest.0 += poynting(e, b, c) * 0.5;
                rest.1 += energy_density(e, b, c) * 0.5;
            }
        }
        (rest, exchange(s.e_charges, s.b_charges, es, bs, c))
    } else {
        (
            (
                poynting(s.e_rest, s.b_rest, c),
                energy_density(s.e_rest, s.b_rest, c),
            ),
            exchange(s.e_charges, s.b_charges, s.e_rest, s.b_rest, c),
        )
    };
    Some([
        (own.0 + sx + rest.0, own.1 + ux + rest.1),
        own,
        (sx, ux),
        rest,
    ])
}

/// `flows_at` for one part.
fn flow_at(
    view: &RadiationView,
    field: &LevelField,
    mode: MapMode,
    x: DVec3,
    t: f64,
    c: f64,
    part: FlowPart,
) -> Option<(DVec3, f64)> {
    flows_at(view, field, mode, x, t, c).map(|f| f[part.index()])
}

/// `(E, B_z)` of the shown sources at `x`, time `t`; `None` inside a source.
fn sample(
    view: &RadiationView,
    field: &LevelField,
    mode: MapMode,
    x: DVec3,
    t: f64,
    c: f64,
) -> Option<(DVec3, f64)> {
    match mode {
        MapMode::Waves => {
            if field
                .antennas
                .iter()
                .any(|a| (x - a.position).length() < a.radius)
            {
                return None;
            }
            let f = field.sample(x, t);
            Some((f.e, f.b.z))
        }
        MapMode::Total => {
            if view.obstacles.iter().any(|o| o.signed_distance(x) < 0.0) {
                return None;
            }
            let f = field.sample(x, t);
            let (e, bz) = charges(view, x, t, c)?;
            Some((f.e + e, f.b.z + bz))
        }
        _ => charges(view, x, t, c),
    }
}

/// `(E, B_z)` of the moving charges at `x`, time `t` (`None` right at one): retarded
/// (Liénard–Wiechert) for finite c, Coulomb for c = ∞. Options: only the radiation
/// (acceleration) part; or only what the quasi-static beam interaction leaves out (the
/// retarded field minus the quasi-static fields of the present states).
fn charges(view: &RadiationView, x: DVec3, t: f64, c: f64) -> Option<(DVec3, f64)> {
    use physics::lienard::Worldline;
    let (mut e, mut bz) = (DVec3::ZERO, 0.0);
    for (i, s) in view.sources.iter().enumerate() {
        let w = &s.line;
        // A magnetic moment's dipole field, B_z = −m/r³ in the plane: from the position
        // now for c = ∞ (exact), from the retarded position otherwise (the moving
        // dipole's velocity and radiation terms, O(v/c), are left out); where it
        // stopped once absorbed there, none once drained. Not a radiation field.
        if s.moment != 0.0 && !view.radiation_only {
            let at = if !c.is_finite() {
                match s.end {
                    Some((te, xe)) if t >= te => s.stays.then_some(xe),
                    _ => Some(w.state(t).0),
                }
            } else {
                match s.end {
                    Some((te, xe)) if c * (t - te) >= (x - xe).length() => s.stays.then_some(xe),
                    _ => Some(w.state(lienard::retarded_time(w, c, x, t)).0),
                }
            };
            if let Some(r0) = at {
                let r = (x - r0).length();
                if r <= 0.15 {
                    return None;
                }
                bz -= s.moment / (r * r * r);
            }
        }
        if s.charge == 0.0 {
            continue;
        }
        if !c.is_finite() {
            if let Some((te, xe)) = s.end
                && t >= te
            {
                // Absorbed: its charge at rest where it stopped, fading as it flies on
                // into a screening cup, or drained.
                if s.stays {
                    let d = x - xe;
                    if d.length() <= 0.15 {
                        return None;
                    }
                    e += d * (s.charge / d.length().powi(3));
                } else if let Some(fd) = s.fade {
                    // Screened inside the cup: not masked like a free charge.
                    let d = x - fd.position(t);
                    e += d * (s.charge * fd.factor(t) / d.length().max(1e-6).powi(3));
                }
                continue;
            }
            let d = x - w.state(t).0;
            if d.length() <= 0.15 {
                return None;
            }
            e += d * (s.charge / d.length().powi(3));
            continue;
        }
        // Once the light cone of its absorption has passed (as in the dynamics): a charge
        // at rest where it stopped, or nothing if it was drained. One flying on into a
        // screening cup keeps its world line, its charge fading at the retarded time.
        let absorbed = s
            .end
            .filter(|&(te, xe)| s.fade.is_none() && c * (t - te) >= (x - xe).length());
        if let Some((_, xe)) = absorbed {
            if s.stays && !view.radiation_only {
                let d = x - xe;
                if d.length() <= 0.15 {
                    return None;
                }
                e += d * (s.charge / d.length().powi(3));
            }
        } else {
            let f = lienard::fields(w, s.charge, c, x, t);
            let r = x - w.state(f.retarded_time).0;
            // Seen inside the cup (screened): not masked like a free charge.
            let in_cup = s.fade.is_some_and(|fd| f.retarded_time > fd.t_off);
            if r.length() <= 0.15 && !in_cup {
                return None;
            }
            let k = s.fade.map_or(1.0, |fd| fd.factor(f.retarded_time));
            if view.radiation_only {
                // B = n × E / c holds for each part separately.
                e += f.e_radiation * k;
                bz += r.normalize().cross(f.e_radiation).z / c * k;
            } else {
                e += f.e() * k;
                bz += f.b.z * k;
            }
        }
        // What the quasi-static interaction uses instead: the fields of the present
        // state continued back with constant acceleration, while the particle flies.
        if view.neglected_only
            && let Some((te, xe)) = s.end
            && t >= te
        {
            // The dynamics has the absorbed charge at rest at once, fading in its cup as
            // it moves on uniformly, or drained.
            if s.stays {
                let d = x - xe;
                e -= d * (s.charge / d.length().powi(3));
            } else if let Some(fd) = s.fade {
                let (eq, bq) = fd.quasi_static_fields(s.charge, c, x, t);
                e -= eq;
                bz -= bq.z;
            }
        } else if view.neglected_only {
            let (r, v, a) = w.state(t);
            if (x - r).length() <= 0.15 {
                return None;
            }
            let (eq, bq) = match (s.moment == 0.0)
                .then(|| source_fields(view, i, t, c))
                .flatten()
            {
                Some((_, _, _, es, bs)) => {
                    physics::beam::continued_fields(s.charge, s.mass, c, x, r, v, a, (es, bs))
                }
                None => physics::beam::accelerated_fields(s.charge, c, x, 0.0, r, v, a),
            };
            e -= eq;
            bz -= bq.z;
        }
    }
    Some((e, bz))
}

/// Source `i` at time `t` while it flies: position, velocity, acceleration, and the fields
/// it feels as the quasi-static interaction takes them (its continued past moves in them):
/// the flight's own field and the other particles' (those flying by their uniform motion,
/// those stopped by a body at rest where they stopped, those in a screening cup fading).
#[allow(clippy::type_complexity)]
fn source_fields(
    view: &RadiationView,
    i: usize,
    t: f64,
    c: f64,
) -> Option<(DVec3, DVec3, DVec3, DVec3, DVec3)> {
    use physics::lienard::Worldline;
    let s = &view.sources[i];
    if s.end.is_some_and(|(te, _)| t >= te) {
        return None;
    }
    let (r, v, a) = s.line.state(t);
    let f = view.source_field.sample(r, t);
    let (mut e, mut b) = (f.e, f.b);
    for (j, o) in view.sources.iter().enumerate() {
        if j == i || o.charge == 0.0 {
            continue;
        }
        match o.end {
            Some((te, xe)) if t >= te => {
                if o.stays {
                    let d = r - xe;
                    e += d * (o.charge / d.length().max(1e-6).powi(3));
                } else if let Some(fd) = o.fade {
                    let (ef, bf) = fd.quasi_static_fields(o.charge, c, r, t);
                    e += ef;
                    b += bf;
                }
            }
            _ => {
                let (ro, vo, _) = o.line.state(t);
                let (eh, bh) = physics::beam::heaviside_fields(o.charge, c, r, ro, vo);
                e += eh;
                b += bh;
            }
        }
    }
    Some((r, v, a, e, b))
}

/// The velocity a tracer moves with: the energy's velocity `S/u` of the part shown (at most
/// c; exactly c in a radiation field). The exchange terms have no energy velocity of their
/// own (their density can vanish or be negative): there, along S at a display speed, c at
/// full colour. `None` inside a body or where the flow is below the colour scale.
#[allow(clippy::too_many_arguments)]
fn tracer_velocity(
    view: &RadiationView,
    field: &LevelField,
    mode: MapMode,
    x: DVec3,
    t: f64,
    part: FlowPart,
    linear: bool,
    range: f64,
) -> Option<DVec3> {
    let c = view.c;
    let (s, u) = flow_at(view, field, mode, x, t, c, part)?;
    let shown = colour_value(s.length(), view.s_sat[part.index()], linear, range);
    if shown < 0.05 {
        return None;
    }
    Some(match part {
        FlowPart::Exchange => s.normalize_or_zero() * (c * shown),
        _ if u > 0.0 => s / u,
        _ => DVec3::ZERO,
    })
}

/// Moves the energy tracers (midpoint rule) with `tracer_velocity` over `dt` (time units)
/// at time `t`. A tracer is placed anew, at random where the flow shows, when it has lived
/// `TRACER_LIFE` seconds, left the arena, met a body or reached a place where the flow is
/// below the colour scale. `on`: tracers shown (else removed).
#[allow(clippy::too_many_arguments)]
fn advance_tracers(
    view: &mut RadiationView,
    field: &LevelField,
    mode: MapMode,
    bounds: &physics::geometry::Aabb,
    t: f64,
    part: FlowPart,
    on: bool,
    dt: f64,
    real_dt: f32,
    linear: bool,
    range: f64,
) {
    if !on || !view.c.is_finite() {
        view.tracers.clear();
        return;
    }
    let mut tracers = std::mem::take(&mut view.tracers);
    let mut seed = view.tracer_seed;
    let mut random = || {
        // SplitMix64.
        seed = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = seed;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        #[allow(clippy::cast_precision_loss)]
        let r = (z >> 11) as f64 / (1u64 << 53) as f64;
        r
    };
    let size = bounds.max - bounds.min;
    let place = |random: &mut dyn FnMut() -> f64| {
        bounds.min + DVec3::new(random() * size.x, random() * size.y, 0.0)
    };
    // New tracers start at random ages, so that they do not all renew at once.
    while tracers.len() < TRACERS {
        #[allow(clippy::cast_possible_truncation)]
        let age = (random() * f64::from(TRACER_LIFE)) as f32;
        tracers.push(Tracer {
            x: place(&mut random),
            age,
        });
    }
    let v = &*view;
    let inside = |x: DVec3| {
        x.x >= bounds.min.x && x.x <= bounds.max.x && x.y >= bounds.min.y && x.y <= bounds.max.y
    };
    // Move every tracer (in parallel); those to renew come back as `None`.
    let moved: Vec<Option<Tracer>> = tracers
        .par_iter()
        .map(|tr| {
            if tr.age + real_dt > TRACER_LIFE {
                return None;
            }
            let v0 = tracer_velocity(v, field, mode, tr.x, t, part, linear, range)?;
            let mid = tr.x + v0 * (0.5 * dt);
            let v1 = tracer_velocity(v, field, mode, mid, t + 0.5 * dt, part, linear, range)?;
            let x = tr.x + v1 * dt;
            inside(x).then_some(Tracer {
                x,
                age: tr.age + real_dt,
            })
        })
        .collect();
    let mut out = Vec::with_capacity(TRACERS);
    for m in moved {
        match m {
            Some(tr) => out.push(tr),
            None => {
                // Anew where the flow shows (a few tries), else anywhere.
                let mut x = place(&mut random);
                for _ in 0..8 {
                    if tracer_velocity(v, field, mode, x, t, part, linear, range).is_some() {
                        break;
                    }
                    x = place(&mut random);
                }
                out.push(Tracer { x, age: 0.0 });
            }
        }
    }
    view.tracers = out;
    view.tracer_seed = seed;
}

/// Draws the energy tracers (called from `draw::draw`): dots fading in and out.
pub fn draw_tracers(gizmos: &mut Gizmos, view: &RadiationView) {
    for tr in &view.tracers {
        let a = (tr.age / 0.5)
            .min((TRACER_LIFE - tr.age) / 0.5)
            .clamp(0.0, 1.0)
            * 0.85;
        #[allow(clippy::cast_possible_truncation)]
        let p = Vec2::new(tr.x.x as f32, tr.x.y as f32);
        gizmos.circle_2d(p, 0.07, Color::srgba(0.85, 1.0, 0.85, a));
    }
}

/// Draws the E arrows (called from `draw::draw`).
pub fn draw_arrows(gizmos: &mut Gizmos, view: &RadiationView) {
    let color = Color::srgba(1.0, 1.0, 1.0, 0.7);
    for &(p, d) in &view.arrows {
        let tip = p + d * 0.5;
        let tail = p - d * 0.5;
        gizmos.line_2d(tail, tip, color);
        let back = -d.normalize_or_zero() * (d.length() * 0.35).min(0.3);
        gizmos.line_2d(tip, tip + Vec2::from_angle(0.5).rotate(back), color);
        gizmos.line_2d(tip, tip + Vec2::from_angle(-0.5).rotate(back), color);
    }
}
