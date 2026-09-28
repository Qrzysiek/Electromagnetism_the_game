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
}

/// A moving charge of the field views: its world line, its charge, and when and where it
/// was absorbed (its field then disappears where the light cone of that moment passes,
/// as in the beam dynamics).
struct Source {
    line: SampledWorldline,
    /// The samples `(t, x, v, a)` of the world line (for the GPU).
    samples: Vec<(f64, DVec3, DVec3, DVec3)>,
    charge: f64,
    /// Magnetic moment along z (its dipole field B_z = −m/r³ in the plane).
    moment: f64,
    end: Option<(f64, DVec3)>,
    /// Whether its charge stays where it was absorbed (a body; `Fate::Stop`) rather than
    /// being drained (the detector).
    stays: bool,
}

#[derive(ShaderType, Debug, Clone, Copy, Default)]
pub struct FieldParams {
    /// min.x, min.y, size.x, size.y of the arena.
    pub area: Vec4,
    /// c (0 for infinite), B saturation, E saturation, dynamic range.
    pub scales: Vec4,
    /// Static grid width, height, 1 if present; number of charges.
    pub grid: UVec4,
    /// Antennas, waves, flags (1 radiation only, 2 left out by the model, 4 colour |E|).
    pub counts: UVec4,
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

#[derive(Resource)]
pub struct RadiationView {
    material: Handle<FieldMaterial>,
    entity: Entity,
    /// (revision, flight, mode, flight time of the preview, radiation only, neglected
    /// only) the scales and world lines were computed for.
    key: Option<(u64, usize, MapMode, u64, bool, bool)>,
    radiation_only: bool,
    /// Show only the part of the beam's field the quasi-static interaction leaves out.
    neglected_only: bool,
    /// The moving charges shown: one particle, or every particle of a beam.
    sources: Vec<Source>,
    b_sat: f64,
    e_sat: f64,
    /// Largest (99th percentile) B_z and |E| of the part shown (radiation part, or left
    /// out by the model) relative to the full field's scale; 1 for the full field.
    pub part_b: f64,
    pub part_e: f64,
    /// Time and style (quantity, range, arrows) of the current arrows.
    time: f64,
    style: (FieldQuantity, u64, bool),
    /// Obstacles of the shown flight (the field is not drawn inside sources).
    obstacles: Vec<Shape>,
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
        part_b: 1.0,
        part_e: 1.0,
        time: f64::NAN,
        style: (FieldQuantity::Bz, 0, false),
        obstacles: Vec::new(),
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
/// acceleration continued back (as the exact retarded beam dynamics assumes, PHYSICS.md
/// §3.3) over the light time across the arena, so that no switch-on shell appears.
fn with_past(
    line: &[(f64, DVec3, DVec3, DVec3)],
    bounds: &physics::geometry::Aabb,
    c: f64,
) -> Vec<(f64, DVec3, DVec3, DVec3)> {
    let Some(&(t0, x0, v0, a0)) = line.first() else {
        return Vec::new();
    };
    let span = (bounds.max - bounds.min).length() / c * 1.2;
    // Constant acceleration while its velocity change stays below 0.1c.
    let span = span.min(0.1 * c / a0.length().max(1e-300));
    let mut out: Vec<_> = (1..=32)
        .rev()
        .map(|k| {
            let dt = -span * f64::from(k) / 32.0;
            (
                t0 + dt,
                x0 + v0 * dt + a0 * (0.5 * dt * dt),
                v0 + a0 * dt,
                a0,
            )
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
    // A beam: every particle's world line, from the active disturbance's flight.
    let beam = level.has_beams().then(|| {
        let d = game
            .active_disturbance
            .min(game.beams.len().saturating_sub(1));
        game.beams.get(d).and_then(|b| b.preview.as_ref())
    });
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
                let shots = level.beam_shots();
                b.map_or_else(Vec::new, |b| {
                    b.worldlines
                        .iter()
                        .zip(&b.ends)
                        .zip(&shots)
                        .zip(&b.outcomes)
                        .filter(|(((l, _), _), _)| !l.is_empty())
                        .map(|(((l, &end), &s), &o)| {
                            let samples = with_past(l, &bounds, c);
                            Source {
                                line: SampledWorldline::new(&samples),
                                samples,
                                charge: level.shots[s].particle.charge,
                                moment: level.shots[s].particle.moment,
                                end,
                                stays: matches!(o, physics::trajectory::Outcome::Collided(_)),
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
                        moment: level.shots[shot].particle.moment,
                        end,
                        stays,
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
    for s in &view.sources {
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
        items.push([s.moment as f32, 0.0, 0.0, 0.0]);
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
        items.push([a.omega as f32, ph as f32, a.radius as f32, 0.0]);
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
            items.push([w.amplitude as f32, w.omega as f32, ph as f32, 0.0]);
            n_waves += 1;
        }
    }
    if let Some(mut b) = buffers.get_mut(&samples_h) {
        *b = buffer(samples);
    }
    if let Some(mut b) = buffers.get_mut(&items_h) {
        *b = buffer(items);
    }
    let flags = u32::from(view.radiation_only)
        | (u32::from(view.neglected_only) << 1)
        | (u32::from(quantity == FieldQuantity::E) << 2)
        | (u32::from(linear) << 3);
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
        m.params.counts = UVec4::new(field.antennas.len() as u32, n_waves, flags, 0);
    }
    if let Ok((mut vis, mut tr)) = quads.get_mut(view.entity) {
        *vis = Visibility::Visible;
        let center = (bounds.max + bounds.min) * 0.5;
        tr.translation = Vec3::new(center.x as f32, center.y as f32, -9.0);
        tr.scale = Vec3::new(size.x as f32, size.y as f32, 1.0);
    }

    // E arrows (CPU): only when the time or the style changed.
    let style = (
        game.field_quantity,
        game.field_range_decades.to_bits()
            ^ game.field_gain_decades.to_bits().rotate_left(1)
            ^ u64::from(linear),
        game.show_field_arrows,
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
                let (e, _) = sample(v, fr, mode, x, t, c)?;
                let len =
                    (colour_value(e.length(), v.e_sat, linear, range) * 0.9 * ARROW_SPACING) as f32;
                let d = Vec2::new(e.x as f32, e.y as f32).normalize_or_zero();
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
    for s in &view.sources {
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
                // Absorbed: its charge at rest where it stopped, or drained.
                if s.stays {
                    let d = x - xe;
                    if d.length() <= 0.15 {
                        return None;
                    }
                    e += d * (s.charge / d.length().powi(3));
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
        // at rest where it stopped, or nothing if it was drained.
        let absorbed = s.end.filter(|&(te, xe)| c * (t - te) >= (x - xe).length());
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
            if r.length() <= 0.15 {
                return None;
            }
            if view.radiation_only {
                // B = n × E / c holds for each part separately.
                e += f.e_radiation;
                bz += r.normalize().cross(f.e_radiation).z / c;
            } else {
                e += f.e();
                bz += f.b.z;
            }
        }
        // What the quasi-static interaction uses instead: the fields of the present
        // state continued back with constant acceleration, while the particle flies.
        if view.neglected_only
            && let Some((te, xe)) = s.end
            && t >= te
        {
            // The dynamics has the absorbed charge at rest at once (or drained).
            if s.stays {
                let d = x - xe;
                e -= d * (s.charge / d.length().powi(3));
            }
        } else if view.neglected_only {
            let (r, v, a) = w.state(t);
            if (x - r).length() <= 0.15 {
                return None;
            }
            let (eq, bq) = physics::beam::accelerated_fields(s.charge, c, x, 0.0, r, v, a);
            e -= eq;
            bz -= bq.z;
        }
    }
    Some((e, bz))
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
