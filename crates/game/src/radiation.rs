//! Time-dependent field views (visual only): the fields of antennas and plane waves
//! ("Waves"), the Liénard–Wiechert field of the particle itself, or of every particle of
//! a beam ("Particle field"), and everything at once ("Total": static sources, antennas,
//! waves, disturbances and the particles), at the animation time. For a beam at finite c
//! with the quasi-static interaction, the view can show only what that model leaves out
//! of the dynamics: the full retarded field minus the quasi-static fields of the
//! particles' present states (PHYSICS.md §3.3). Computed in f64 on the CPU with the tested physics code
//! (`physics::antenna`, `physics::external`, `physics::lienard`), shown as a texture of
//! B_z (in the plane B is exactly perpendicular to it) plus optional E arrows.

use bevy::asset::RenderAssetUsages;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
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

/// Texture pixels per cell.
const PX_PER_CELL: f64 = 5.0;
/// Colour texels per field sample (smooth interpolation of the sampled values).
const UPSAMPLE: u32 = 4;
/// Largest change of the compressed value between neighbouring samples that is still
/// interpolated; steeper cells are evaluated exactly at every colour texel.
const REFINE_STEP: f32 = 0.3;
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
    charge: f64,
    end: Option<(f64, DVec3)>,
}

#[derive(Resource)]
pub struct RadiationView {
    image: Handle<Image>,
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
    /// Time and style (quantity, range) of the current texture.
    time: f64,
    style: (FieldQuantity, u64, bool),
    /// Obstacles of the shown flight (the field is not drawn inside sources).
    obstacles: Vec<Shape>,
    pub arrows: Vec<(Vec2, Vec2)>,
}

pub fn setup(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    let image = images.add(blank(1, 1));
    let entity = commands
        .spawn((
            Sprite {
                image: image.clone(),
                custom_size: Some(Vec2::ONE),
                ..default()
            },
            Transform::from_xyz(0.0, 0.0, -9.0),
            Visibility::Hidden,
        ))
        .id();
    commands.insert_resource(RadiationView {
        image,
        entity,
        key: None,
        radiation_only: false,
        neglected_only: false,
        sources: Vec::new(),
        b_sat: 1.0,
        e_sat: 1.0,
        time: f64::NAN,
        style: (FieldQuantity::Bz, 0, false),
        obstacles: Vec::new(),
        arrows: Vec::new(),
    });
}

fn blank(w: u32, h: u32) -> Image {
    Image::new_fill(
        Extent3d {
            width: w,
            height: h,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        &[0, 0, 0, 0],
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::default(),
    )
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
    level.c().is_finite()
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
fn worldline(preview: &Preview, mass: f64, c: f64) -> Option<SampledWorldline> {
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
    (!samples.is_empty()).then(|| SampledWorldline::new(&samples))
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
fn percentile99(mut v: Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() {
        return 1.0;
    }
    v.sort_by(f64::total_cmp);
    let i = (v.len() - 1) * 99 / 100;
    v[i].max(1e-300)
}

/// Signed compression into [−1, 1] on an asinh scale: values from `sat / range` to `sat`
/// are distinguishable (`range` = 10^decades).
fn compress(v: f64, sat: f64, range: f64) -> f64 {
    let r = sat / range;
    ((v / r).asinh() / (sat / r).asinh()).clamp(-1.0, 1.0)
}

#[allow(
    clippy::too_many_arguments,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub fn update(
    game: Res<Game>,
    mut view: ResMut<RadiationView>,
    mut images: ResMut<Assets<Image>>,
    mut sprites: Query<(&mut Visibility, &mut Transform, &mut Sprite)>,
) {
    let mode = match game.map {
        Some(m @ (MapMode::Waves | MapMode::ParticleField | MapMode::Total)) => m,
        _ => {
            if let Ok((mut v, _, _)) = sprites.get_mut(view.entity) {
                *v = Visibility::Hidden;
            }
            view.arrows.clear();
            return;
        }
    };
    let level = &game.editor.level;
    if level.shots.is_empty() || (!level.c().is_finite() && mode != MapMode::Total) {
        return;
    }
    let c = level.c();
    let flight = game.active_flight();
    let (shot, disturbance) = level.flight_of(flight);
    let preview = game.flights.get(flight).and_then(|f| f.preview.as_ref());

    // Scales and world line: once per setup, flight and mode.
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
    let _ = disturbance;
    if view.key != Some(key) {
        view.key = Some(key);
        view.time = f64::NAN;
        view.radiation_only = game.radiation_only;
        // Scales from the full field: the part left out by the model is shown on the same
        // colour scale, so that its true size is seen (set after the scales).
        view.neglected_only = false;
        view.sources = match beam {
            Some(b) => {
                let shots = level.beam_shots();
                b.map_or_else(Vec::new, |b| {
                    b.worldlines
                        .iter()
                        .zip(&b.ends)
                        .zip(&shots)
                        .filter(|((l, _), _)| !l.is_empty())
                        .map(|((l, &end), &s)| Source {
                            line: SampledWorldline::new(&with_past(l, &bounds, c)),
                            charge: level.shots[s].particle.charge,
                            end,
                        })
                        .collect()
                })
            }
            None => preview
                .and_then(|p| worldline(p, level.shots[shot].particle.mass, c))
                .map(|line| Source {
                    line,
                    charge: level.shots[shot].particle.charge,
                    end: None,
                })
                .into_iter()
                .collect(),
        };
        // Sample over one RF period (waves) or over the flight (particle).
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
        let mut bs = Vec::new();
        let mut es = Vec::new();
        for k in 0..8 {
            let t = t0 + (t1 - t0) * (f64::from(k) + 0.5) / 8.0;
            let mut y = bounds.min.y;
            while y <= bounds.max.y {
                let mut x = bounds.min.x;
                while x <= bounds.max.x {
                    if let Some((e, b)) = sample(&view, &field, mode, DVec3::new(x, y, 0.0), t, c) {
                        bs.push(b.abs());
                        es.push(e.length());
                    }
                    x += 1.0;
                }
                y += 1.0;
            }
        }
        view.b_sat = percentile99(bs);
        view.e_sat = percentile99(es);
        view.neglected_only = neglected_only;
        if std::env::var("EM_CAPTURE").is_ok() {
            eprintln!(
                "radiation scales: B {:.3e}, E {:.3e}",
                view.b_sat, view.e_sat
            );
        }
    }

    // Time shown: lab time for waves, flight time for the particle's own field.
    let t_anim = if game.animate {
        game.anim_time
    } else {
        t_final
    };
    let t = match mode {
        MapMode::Waves => launch_time + t_anim,
        _ => t_anim,
    };
    let style = (
        game.field_quantity,
        game.field_range_decades.to_bits(),
        game.show_field_arrows,
    );
    if t.to_bits() == view.time.to_bits() && style == view.style {
        return;
    }
    view.time = t;
    view.style = style;
    let range = 10f64.powf(game.field_range_decades);
    let quantity = game.field_quantity;

    // Texture.
    let started = std::time::Instant::now();
    let size = bounds.max - bounds.min;
    let (w, h) = (
        (size.x * PX_PER_CELL).ceil().max(1.0) as u32,
        (size.y * PX_PER_CELL).ceil().max(1.0) as u32,
    );
    let (b_sat, e_sat) = (view.b_sat, view.e_sat);
    let v = &*view;
    let fr = &field;
    // The field is evaluated at the texel centres of a (w × h) grid. The compressed,
    // signed value is then interpolated bilinearly onto an UPSAMPLE× finer texture
    // before colouring, so zero lines (sign changes) and saturation edges come out
    // smooth instead of breaking into texel-sized dots. Visual only: no new physics.
    let values: Vec<f32> = (0..h)
        .into_par_iter()
        .flat_map_iter(|j| {
            let y = bounds.max.y - (f64::from(j) + 0.5) / f64::from(h) * size.y;
            (0..w).map(move |i| {
                let x = bounds.min.x + (f64::from(i) + 0.5) / f64::from(w) * size.x;
                let f = sample(v, fr, mode, DVec3::new(x, y, 0.0), t, c);
                (match quantity {
                    FieldQuantity::Bz => f.map_or(0.0, |(_, b)| compress(b, b_sat, range)),
                    FieldQuantity::E => f.map_or(0.0, |(e, _)| compress(e.length(), e_sat, range)),
                }) as f32
            })
        })
        .collect();
    let (wu, hu) = (w * UPSAMPLE, h * UPSAMPLE);
    let at = |i: i64, j: i64| {
        let i = i.clamp(0, i64::from(w) - 1) as usize;
        let j = j.clamp(0, i64::from(h) - 1) as usize;
        values[j * w as usize + i]
    };
    // Exact value at a fine pixel (for cells where interpolation is not good enough).
    let exact = |ii: u32, jj: u32| -> f32 {
        let x = bounds.min.x + (f64::from(ii) + 0.5) / f64::from(wu) * size.x;
        let y = bounds.max.y - (f64::from(jj) + 0.5) / f64::from(hu) * size.y;
        let f = sample(v, fr, mode, DVec3::new(x, y, 0.0), t, c);
        (match quantity {
            FieldQuantity::Bz => f.map_or(0.0, |(_, b)| compress(b, b_sat, range)),
            FieldQuantity::E => f.map_or(0.0, |(e, _)| compress(e.length(), e_sat, range)),
        }) as f32
    };
    let pixels: Vec<u8> = (0..hu)
        .into_par_iter()
        .flat_map_iter(|jj| {
            // Position in texel units, texel centres at integers.
            let fy = (f64::from(jj) + 0.5) / f64::from(UPSAMPLE) - 0.5;
            let (j0, ty) = (fy.floor() as i64, (fy - fy.floor()) as f32);
            (0..wu).flat_map(move |ii| {
                let fx = (f64::from(ii) + 0.5) / f64::from(UPSAMPLE) - 0.5;
                let (i0, tx) = (fx.floor() as i64, (fx - fx.floor()) as f32);
                let corners = [
                    at(i0, j0),
                    at(i0 + 1, j0),
                    at(i0, j0 + 1),
                    at(i0 + 1, j0 + 1),
                ];
                let lo = corners.iter().copied().fold(f32::INFINITY, f32::min);
                let hi = corners.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                // Adaptive refinement: across a sign change or a steep step the value is
                // computed exactly at this pixel; elsewhere it is interpolated.
                let s = if (lo < 0.0 && hi > 0.0) || hi - lo > REFINE_STEP {
                    exact(ii, jj)
                } else {
                    let top = corners[0] * (1.0 - tx) + corners[1] * tx;
                    let bottom = corners[2] * (1.0 - tx) + corners[3] * tx;
                    top * (1.0 - ty) + bottom * ty
                };
                let a = (s.abs().powf(0.8) * 230.0) as u8;
                match quantity {
                    FieldQuantity::Bz if s >= 0.0 => [255, 150, 40, a],
                    FieldQuantity::Bz => [40, 190, 255, a],
                    FieldQuantity::E => [255, 235, 140, a],
                }
            })
        })
        .collect();
    let (w, h) = (wu, hu);
    let mut img = blank(w, h);
    img.data = Some(pixels);
    if let Some(mut slot) = images.get_mut(&view.image) {
        *slot = img;
    }
    if let Ok((mut vis, mut tr, mut sprite)) = sprites.get_mut(view.entity) {
        *vis = Visibility::Visible;
        let center = (bounds.max + bounds.min) * 0.5;
        tr.translation = Vec3::new(center.x as f32, center.y as f32, -9.0);
        sprite.custom_size = Some(Vec2::new(size.x as f32, size.y as f32));
    }

    // E arrows.
    let mut arrows = Vec::new();
    if game.show_field_arrows {
        let mut y = bounds.min.y + 0.5 * ARROW_SPACING;
        while y < bounds.max.y {
            let mut x = bounds.min.x + 0.5 * ARROW_SPACING;
            while x < bounds.max.x {
                if let Some((e, _)) = sample(&view, &field, mode, DVec3::new(x, y, 0.0), t, c) {
                    let len =
                        (compress(e.length(), view.e_sat, range) * 0.9 * ARROW_SPACING) as f32;
                    let d = Vec2::new(e.x as f32, e.y as f32).normalize_or_zero();
                    if len > 0.05 {
                        arrows.push((Vec2::new(x as f32, y as f32), d * len));
                    }
                }
                x += ARROW_SPACING;
            }
            y += ARROW_SPACING;
        }
    }
    view.arrows = arrows;
    if std::env::var("EM_CAPTURE").is_ok() {
        eprintln!(
            "radiation view: {:.2} ms",
            started.elapsed().as_secs_f64() * 1e3
        );
    }
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
        if !c.is_finite() {
            if s.end.is_some_and(|(te, _)| t >= te) {
                continue;
            }
            let d = x - w.state(t).0;
            if d.length() <= 0.15 {
                return None;
            }
            e += d * (s.charge / d.length().powi(3));
            continue;
        }
        // Gone once the light cone of its absorption has passed (as in the dynamics).
        let gone = s
            .end
            .is_some_and(|(te, xe)| c * (t - te) >= (x - xe).length());
        if !gone {
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
        if view.neglected_only && s.end.is_none_or(|(te, _)| t < te) {
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
