//! GPU field maps (potential or magnetic field): a quad covering the world bounds with a
//! fragment shader that evaluates the field per pixel (visual only).

use bevy::asset::{AssetPath, embedded_asset, embedded_path};
use bevy::prelude::*;
use bevy::render::render_resource::{AsBindGroup, ShaderType};
use bevy::render::storage::ShaderBuffer;
use bevy::shader::ShaderRef;
use bevy::sprite_render::{Material2d, Material2dPlugin};
use level::beam::Launch;
use physics::external::External;
use physics::field::{FieldSolver, LevelField};
use physics::trajectory::Scenario;

/// Most groups of particles with a dark region; must match `potential.wgsl`.
pub const MAX_LIMITS: usize = 64;

/// Which field the map shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MapMode {
    Potential,
    Magnetic,
    /// Fields of antennas and plane waves at the animation time (`radiation.rs`).
    Waves,
    /// The particle's own Liénard–Wiechert field at the animation time (`radiation.rs`).
    ParticleField,
    /// Everything at once: static sources, antennas, waves, disturbances and the
    /// particle's own field, at the animation time (`radiation.rs`).
    Total,
    /// Tube levels: |E| of the preview's final state, with arrows.
    Electric,
    /// Tube levels: the energy flow (Poynting vector) of the run's electromagnetic field,
    /// with arrows (`em.rs`; `Magnetic` shows its B_z).
    EnergyFlow,
    /// Tube levels with metals: the current density inside them, |j| with arrows, over a
    /// span of time (`level::tube::TubeRun::metal_rates`).
    Current,
    /// Tube levels with metals: the power per area the field gives their carriers (the heat
    /// goal's), signed.
    Heat,
}

/// Uniform parameters of the maps. The sources themselves are in the storage buffer
/// `PotentialMaterial::items`, one category after the other (`Items`).
#[derive(ShaderType, Debug, Clone, Copy)]
pub struct PotentialParams {
    /// Energy limits, one per group of particles (see `params`).
    pub limits: [Vec4; MAX_LIMITS],
    /// Numbers of solid charges, charge clouds, induced charges and magnets.
    pub counts: UVec4,
    /// Numbers of circular coils, straight coil segments, static antennas and electrode
    /// panels.
    pub counts2: UVec4,
    /// The uniform stray fields: `E_x`, `E_y` (potential `−E·(x − origin)`),
    /// `B_z / b_ref`, 0.
    pub uniform: Vec4,
    /// `origin` of the uniform field's potential (x, y), antenna body radius, coil wire
    /// radius.
    pub origin: Vec4,
    pub limit_count: u32,
    pub phi_weight: f32,
    pub u_a: f32,
    pub mode: u32,
    pub moment_weight: f32,
    /// Tube levels: numbers of charged segments and of line charges.
    pub line_segments: u32,
    pub line_charges: u32,
    /// Electric field map: |E| at the top of its colour scale, and the decades below it
    /// still coloured.
    pub e_ref: f32,
    pub e_decades: f32,
    /// Tube levels: the particles drawn as the density of their charge, a glow (z > 0):
    /// x = 1/(2 s²) for a Gaussian of width s, y and w = 1/(|q| at full glow) of each
    /// species, z = 1 + the number of line charges that glow (the particles outside
    /// metals: the metals' carriers and backgrounds after them do not).
    pub glow: Vec4,
    /// Tube levels' field maps (mode 3): the display grid's origin (x, y), its spacing
    /// and the colour scale's unit.
    pub grid: Vec4,
    /// Its samples across and up, the quantity (0: B_z, 1: |S|) and the decades coloured
    /// below the unit (as bits).
    pub grid_n: UVec4,
}

impl Default for PotentialParams {
    fn default() -> Self {
        Self {
            limits: [Vec4::ZERO; MAX_LIMITS],
            counts: UVec4::ZERO,
            counts2: UVec4::ZERO,
            uniform: Vec4::ZERO,
            origin: Vec4::new(0.0, 0.0, 0.3, 0.1),
            limit_count: 0,
            phi_weight: 0.0,
            u_a: 0.0,
            mode: 0,
            moment_weight: 0.0,
            line_segments: 0,
            line_charges: 0,
            e_ref: 1.0,
            e_decades: 1.5,
            glow: Vec4::ZERO,
            grid: Vec4::ZERO,
            grid_n: UVec4::ZERO,
        }
    }
}

/// Width of the particles' density glow (cells): a beam of macroparticles a fraction of a
/// cell apart reads as one stream.
const GLOW_WIDTH: f32 = 0.3;
/// Overlapping macroparticles at full glow, at least: a species whose crowded places hold
/// more glows fully only there (its 99th percentile of the charge per glow-sized cell), so
/// that each species shows its own gradients and a thin, fast stream beside a dense, slow
/// cloud of the other species still shows.
const GLOW_FULL: f32 = 4.0;
/// The map draws at most about this many particles one by one; above it, each species'
/// particles merged per cell of `MAP_BIN` cells (doubled until under it), at their
/// charge's centroid: the map softens the space charge over `TUBE_MAP_SOFTENING` anyway,
/// and the shader sums every item at every pixel.
const MAP_PARTICLES: usize = 1500;
const MAP_BIN: f64 = 0.25;

/// A frame for the map (its potential, field arrows and glow): the particles merged per
/// species and cell when they are many (`MAP_PARTICLES`), projectiles always one by one;
/// a frame with few particles as it is.
#[allow(clippy::cast_possible_truncation)]
pub fn map_frame(frame: &level::tube::TubeFrame) -> std::borrow::Cow<'_, level::tube::TubeFrame> {
    if frame.x.len() <= MAP_PARTICLES {
        return std::borrow::Cow::Borrowed(frame);
    }
    let mut bin = MAP_BIN;
    loop {
        let mut cells: std::collections::BTreeMap<(i64, i64, u8), [f64; 3]> =
            std::collections::BTreeMap::new();
        let mut out = level::tube::TubeFrame {
            self_b: frame.self_b.clone(),
            t: frame.t,
            x: Vec::new(),
            q: Vec::new(),
            kind: Vec::new(),
            id: Vec::new(),
            sigma: frame.sigma.clone(),
            potentials: frame.potentials.clone(),
            currents: frame.currents.clone(),
            heat: frame.heat,
            paths: Vec::new(),
            cells: Vec::new(),
        };
        for ((x, q), kind) in frame.x.iter().zip(&frame.q).zip(&frame.kind) {
            if *kind == level::tube::KIND_PROJECTILE {
                out.x.push(*x);
                out.q.push(*q);
                out.kind.push(*kind);
                continue;
            }
            let key = (
                (f64::from(x[0]) / bin).floor() as i64,
                (f64::from(x[1]) / bin).floor() as i64,
                *kind,
            );
            let c = cells.entry(key).or_insert([0.0; 3]);
            let q = f64::from(*q);
            c[0] += q;
            c[1] += q * f64::from(x[0]);
            c[2] += q * f64::from(x[1]);
        }
        if cells.len() <= MAP_PARTICLES || bin > 4.0 {
            for ((_, _, kind), c) in cells {
                if c[0] != 0.0 {
                    out.x.push([(c[1] / c[0]) as f32, (c[2] / c[0]) as f32]);
                    out.q.push(c[0] as f32);
                    out.kind.push(kind);
                }
            }
            return std::borrow::Cow::Owned(out);
        }
        bin *= 2.0;
    }
}

/// Each species' charge at full glow (negative, positive): see `GLOW_FULL`.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn glow_full(frame: &level::tube::TubeFrame) -> (f32, f32) {
    let mut out = [0.0_f32; 2];
    for (s, kind) in [level::tube::KIND_NEGATIVE, level::tube::KIND_POSITIVE]
        .iter()
        .enumerate()
    {
        let qs: Vec<f32> = frame
            .q
            .iter()
            .zip(&frame.kind)
            .filter(|(_, k)| *k == kind)
            .map(|(q, _)| q.abs())
            .collect();
        if qs.is_empty() {
            continue;
        }
        let mean = qs.iter().sum::<f32>() / qs.len() as f32;
        let mut cells: std::collections::BTreeMap<(i64, i64), f32> =
            std::collections::BTreeMap::new();
        let w = f64::from(GLOW_WIDTH);
        for ((x, q), k) in frame.x.iter().zip(&frame.q).zip(&frame.kind) {
            if k == kind {
                let key = (
                    (f64::from(x[0]) / w).floor() as i64,
                    (f64::from(x[1]) / w).floor() as i64,
                );
                *cells.entry(key).or_insert(0.0) += q.abs();
            }
        }
        let mut per: Vec<f32> = cells.into_values().collect();
        per.sort_by(f32::total_cmp);
        let p99 = per[(per.len() * 99 / 100).min(per.len() - 1)];
        out[s] = p99.max(mean * GLOW_FULL);
    }
    (out[0], out[1])
}

#[derive(Asset, TypePath, AsBindGroup, Debug, Clone)]
pub struct PotentialMaterial {
    #[uniform(0)]
    pub params: PotentialParams,
    /// The sources (`Items`), as `vec4<f32>`.
    #[storage(1, read_only)]
    pub items: Handle<ShaderBuffer>,
}

impl Material2d for PotentialMaterial {
    fn fragment_shader() -> ShaderRef {
        ShaderRef::Path(
            AssetPath::from_path_buf(embedded_path!("potential.wgsl")).with_source("embedded"),
        )
    }
}

pub struct PotentialPlugin;

impl Plugin for PotentialPlugin {
    fn build(&self, app: &mut App) {
        embedded_asset!(app, "potential.wgsl");
        app.add_plugins(Material2dPlugin::<PotentialMaterial>::default());
    }
}

/// A storage buffer must not be empty.
pub fn buffer(v: Vec<[f32; 4]>) -> ShaderBuffer {
    ShaderBuffer::from(if v.is_empty() { vec![[0.0f32; 4]] } else { v })
}

fn count(n: usize) -> u32 {
    u32::try_from(n).expect("fits")
}

/// Field for which the shot's particle circles with this gyroradius (cells): the unit of
/// the magnetic map.
pub const GYRO_REFERENCE: f64 = 5.0;

/// The static part of a level field, which the maps show: without the oscillating
/// antennas and plane waves (a static antenna, ω = 0, is an electrostatic dipole; a
/// static wave term a uniform field). Ramped coils stay, at t = 0.
pub fn static_part(field: &LevelField) -> LevelField {
    let mut f = field.clone();
    f.antennas.retain(|a| a.omega == 0.0);
    f.external.retain(External::is_static);
    f
}

/// The sources of the map, in the order of the storage buffer:
/// - solid charges `(x, y, Q, radius)`, in the plane;
/// - charge clouds `(x, y, Q, R)` (the uniform sphere's potential inside);
/// - charges induced on metal spheres `(x, y, Q, z)` (images and equivalent charges);
/// - magnets `(x, y, μ / b_ref, radius)`;
/// - circular coils `(x, y, radius, κ / b_ref)`;
/// - straight coil segments, two each: `(a_x, a_y, b_x, b_y)`, `(κ / b_ref, 0, 0, 0)`;
/// - static antennas `(x, y, p_x, p_y)` (potential `n·p / r²`);
/// - electrode panels (upper halves, display mesh), three each: `(a, σ)`, `(b, size)`,
///   `(c, area)` with the vertices' (x, y, z): the shader integrates the panels near a
///   pixel exactly and the others by the three-point rule, as `Electrodes::sample` does
///   for pictures (test E6).
struct Items {
    out: Vec<[f32; 4]>,
}

/// Softening of the particles' line charges on the tube map (cells): the preview's
/// segment length, the resolution of its space charge.
const TUBE_MAP_SOFTENING: f64 = level::tube::PREVIEW_SEGMENT;

/// Shader parameters and sources for a tube level's map (z-invariant): the potential of
/// the frame shown (the electrodes' surface charge, applied and induced, and the
/// particles' space charge), as the carriers' energy relative to the cathode's potential
/// then, in units of `|q|` times the largest potential difference to it over the run
/// (`TubeRun::potential_range`: the electrodes' own potentials and their sources' over a
/// period), divided by `gain` (contours every quarter). The electric map colours `decades`
/// below `e_ref`.
#[allow(clippy::cast_possible_truncation, clippy::too_many_arguments)]
pub fn tube_params(
    run: &level::tube::TubeRun,
    frame: &level::tube::TubeFrame,
    raw: &level::tube::TubeFrame,
    cathode: usize,
    sign: f64,
    electric: bool,
    e_ref: f64,
    gain: f64,
    decades: f64,
    glow: bool,
) -> (PotentialParams, Vec<[f32; 4]>) {
    let vc = frame
        .potentials
        .get(cathode)
        .map(|v| f64::from(*v))
        .or_else(|| run.potentials.get(cathode).copied())
        .unwrap_or(0.0);
    // Until the scale was measured at the AC sources' extremes it came from the potentials
    // at t = 0, all 0 with AC anodes: an infinite gain, every colour saturated, a sharp
    // wall between blue and red.
    let (lo, hi) = run.potential_range;
    let scale = (hi - vc).abs().max((lo - vc).abs()).max(1e-9) / gain.max(1e-9);
    let mut items = Vec::new();
    for ((a, b), s) in run.segments.iter().zip(&frame.sigma) {
        items.push(v4(a.x, a.y, b.x, b.y));
        items.push(v4(f64::from(*s), 0.0, 0.0, 0.0));
    }
    let eps2 = TUBE_MAP_SOFTENING * TUBE_MAP_SOFTENING;
    // The particles: first those that glow, then a metal's carriers (its density is
    // uniform: they would only veil the potential).
    let in_metal = |x: &[f32; 2]| {
        run.metals.iter().any(|m| {
            (m.lo.x..=m.hi.x).contains(&f64::from(x[0]))
                && (m.lo.y..=m.hi.y).contains(&f64::from(x[1]))
        })
    };
    let mut glowing = 0;
    for metal in [false, true] {
        #[allow(clippy::cast_possible_truncation)]
        for (x, q) in frame.x.iter().zip(&frame.q) {
            if in_metal(x) == metal {
                items.push([x[0], x[1], *q, eps2 as f32]);
                glowing += usize::from(!metal);
            }
        }
    }
    // The metals' backgrounds (the ions' charge): fixed line charges on the cells the
    // carriers are binned on (`map_frame`), so that a neutral metal reads neutral.
    let mut charges = frame.x.len();
    for m in &run.metals {
        #[allow(clippy::cast_possible_truncation)]
        let (i0, i1, j0, j1) = (
            (m.lo.x / MAP_BIN).floor() as i64,
            (m.hi.x / MAP_BIN).ceil() as i64,
            (m.lo.y / MAP_BIN).floor() as i64,
            (m.hi.y / MAP_BIN).ceil() as i64,
        );
        for j in j0..j1 {
            for i in i0..i1 {
                #[allow(clippy::cast_precision_loss)]
                let (x0, y0) = (i as f64 * MAP_BIN, j as f64 * MAP_BIN);
                let (ax, bx) = (x0.max(m.lo.x), (x0 + MAP_BIN).min(m.hi.x));
                let (ay, by) = (y0.max(m.lo.y), (y0 + MAP_BIN).min(m.hi.y));
                if bx > ax && by > ay {
                    items.push(v4(
                        0.5 * (ax + bx),
                        0.5 * (ay + by),
                        m.background * (bx - ax) * (by - ay),
                        eps2,
                    ));
                    charges += 1;
                }
            }
        }
    }
    // The glow's scale per species (see `GLOW_FULL`), from the particles one by one.
    let (full_neg, full_pos) = glow_full(raw);
    let inv = |f: f32| if f > 0.0 { 1.0 / f } else { 0.0 };
    let glow = if glow && (full_neg > 0.0 || full_pos > 0.0) {
        #[allow(clippy::cast_precision_loss)]
        Vec4::new(
            1.0 / (2.0 * GLOW_WIDTH * GLOW_WIDTH),
            inv(full_neg),
            1.0 + glowing as f32,
            inv(full_pos),
        )
    } else {
        Vec4::ZERO
    };
    let out = PotentialParams {
        glow,
        phi_weight: (sign / scale) as f32,
        u_a: (sign * vc / scale) as f32,
        line_segments: count(run.segments.len()),
        line_charges: count(charges),
        mode: u32::from(electric) * 2,
        e_ref: e_ref as f32,
        e_decades: decades as f32,
        ..PotentialParams::default()
    };
    (out, items)
}

#[allow(clippy::cast_possible_truncation)]
fn v4(a: f64, b: f64, c: f64, d: f64) -> [f32; 4] {
    [a as f32, b as f32, c as f32, d as f32]
}

/// Shader parameters and sources for one shot. The shader sums the potential `Φ` of the
/// static sources and the field `B_z / b_ref` of the magnetic ones (pre-scaled by
/// `1/b_ref`, with `b_ref = p / (|q| r_ref)` for the momentum `p` of the energy unit `T₀`,
/// or for a neutral particle with a magnetic moment `m`, `b_ref = T₀ / |m|`); the colours
/// show the shot's
/// `U/T₀ = (q Φ − m B_z)/T₀`. `Φ` is the potential of the flight's own field (its static
/// part): point charges, charge clouds, metal, electrodes, static antennas, and the
/// uniform stray field of the active disturbance, `−E·x` up to a constant (the shader's
/// zero is the launch point, which keeps its f32 values small).
///
/// The dark region is where every particle in `limits` is forbidden by energy
/// conservation: each group (one shot, or one beam) contributes the highest total energy
/// `E = T + U(x₀)` of its particles, and a point is dark only when `U > E` there for
/// every group. More than `MAX_LIMITS` groups: no dark region (never too large).
#[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
pub fn params(
    scn: &Scenario<LevelField>,
    limits: &[Vec<Launch>],
    radii: (f64, f64, f64),
    mode: MapMode,
    energy_unit: f64,
) -> (PotentialParams, Vec<[f32; 4]>) {
    let (charge_radius, magnet_radius, antenna_radius) = radii;
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    // Colours in units of the shot's energy unit (its T₀, or for a particle launched at
    // rest the kinetic energy it reaches, `Game::energy_unit`).
    let t0 = kin.kinetic_energy(scn.p0).max(energy_unit).max(1e-300);
    let (q, m) = (scn.particle.charge, scn.particle.moment);
    // The momentum of the energy unit (a particle launched at rest: of the largest
    // kinetic energy it reaches; with |p₀| it once gave b_ref = 1e-300 and an infinite map).
    let p_unit = kin
        .momentum_from_kinetic_energy(t0, physics::DVec3::X)
        .length();
    let b_ref = if q != 0.0 {
        (p_unit / (q.abs() * GYRO_REFERENCE)).max(1e-300)
    } else if m != 0.0 {
        t0 / m.abs()
    } else {
        1.0
    };
    let f = static_part(&scn.field);
    // The uniform stray fields (static disturbance terms); their potential is −E·x on the
    // CPU, −E·(x − origin) in the shader: CPU potentials are shifted by E·origin to match.
    let (mut e_u, mut b_u) = (physics::DVec3::ZERO, 0.0);
    for ext in &f.external {
        let s = ext.sample(physics::DVec3::ZERO, 0.0);
        e_u += s.e;
        b_u += s.b.z;
    }
    let origin = scn.x0;
    let shift = e_u.dot(origin);
    let energy_at = |x: physics::DVec3, qx: f64, mx: f64| {
        let s = f.sample(x, 0.0);
        qx * (s.phi + shift) - mx * s.b.z
    };
    let mut out = PotentialParams {
        u_a: (energy_at(scn.x0, q, m) / t0) as f32,
        phi_weight: (q / t0) as f32,
        moment_weight: (-m * b_ref / t0) as f32,
        mode: match mode {
            MapMode::Potential => 0,
            MapMode::Magnetic | MapMode::Waves | MapMode::ParticleField | MapMode::Total => 1,
            MapMode::Electric | MapMode::EnergyFlow | MapMode::Current | MapMode::Heat => 2,
        },
        uniform: Vec4::new(e_u.x as f32, e_u.y as f32, (b_u / b_ref) as f32, 0.0),
        origin: Vec4::new(
            origin.x as f32,
            origin.y as f32,
            antenna_radius as f32,
            f.loops
                .first()
                .map(|l| l.wire_radius)
                .or_else(|| f.polygons.first().map(|p| p.wire_radius))
                .unwrap_or(0.1) as f32,
        ),
        ..PotentialParams::default()
    };
    let mut items = Items { out: Vec::new() };
    let mut push = |v: [f32; 4]| items.out.push(v);
    // Point charges (solid), then clouds.
    let (mut solid, mut clouds) = (0, 0);
    for (pos, qc, radius) in f.coulomb.sources() {
        if radius == 0.0 {
            push(v4(pos.x, pos.y, qc, charge_radius));
            solid += 1;
        }
    }
    for (pos, qc, radius) in f.coulomb.sources() {
        if radius > 0.0 {
            push(v4(pos.x, pos.y, qc, radius));
            clouds += 1;
        }
    }
    // Charges induced on metal spheres (display resolution), with their z.
    let mut induced = 0;
    for (pos, qc) in f.conductors.induced.charges() {
        push(v4(pos.x, pos.y, qc, pos.z));
        induced += 1;
    }
    for d in &f.dipoles {
        push(v4(
            d.position.x,
            d.position.y,
            d.moment.z / b_ref,
            magnet_radius,
        ));
    }
    for l in &f.loops {
        push(v4(l.center.x, l.center.y, l.radius, l.kappa / b_ref));
    }
    let mut segments = 0;
    for poly in &f.polygons {
        let v = &poly.vertices;
        for i in 0..v.len() {
            let (a, b) = (v[i], v[(i + 1) % v.len()]);
            push(v4(a.x, a.y, b.x, b.y));
            push(v4(poly.kappa / b_ref, 0.0, 0.0, 0.0));
            segments += 1;
        }
    }
    for a in &f.antennas {
        let (p, _, _) = a.moment(0.0);
        push(v4(a.position.x, a.position.y, p.x, p.y));
    }
    let mut panels = 0;
    for (t, s) in f.electrodes.panels() {
        push(v4(t.a.x, t.a.y, t.a.z, s));
        push(v4(t.b.x, t.b.y, t.b.z, t.size()));
        push(v4(t.c.x, t.c.y, t.c.z, t.area()));
        panels += 1;
    }
    out.counts = UVec4::new(
        count(solid),
        count(clouds),
        count(induced),
        count(f.dipoles.len()),
    );
    out.counts2 = UVec4::new(
        count(f.loops.len()),
        count(segments),
        count(f.antennas.len()),
        count(panels),
    );
    if limits.len() <= MAX_LIMITS {
        for (i, group) in limits.iter().enumerate() {
            // Highest total energy of the group, in units of its highest kinetic energy.
            let (mut e, mut t) = (f64::NEG_INFINITY, 1e-300_f64);
            let (mut qg, mut mg) = (0.0, 0.0);
            for l in group {
                let kin = physics::dynamics::Kinematics::new(l.particle.mass, scn.c);
                let tk = kin.kinetic_energy(l.p0);
                (qg, mg) = (l.particle.charge, l.particle.moment);
                e = e.max(tk + energy_at(l.x0, qg, mg));
                t = t.max(tk);
            }
            out.limits[i] = Vec4::new(
                (qg / t) as f32,
                (-mg * b_ref / t) as f32,
                (e / t) as f32,
                0.0,
            );
        }
        out.limit_count = count(limits.len());
    }
    (out, items.out)
}

/// Shader parameters and sources for a tube level's field map (mode 3): the snapshot of
/// its electromagnetic field (`em.rs`), B_z (`magnetic`) or the Poynting vector's
/// magnitude, coloured on a log scale over `decades` below `scale`; the particles' glow
/// over it as on the potential map (`tube_params`).
#[allow(clippy::cast_possible_truncation, clippy::too_many_arguments)]
pub fn tube_field_params(
    run: &level::tube::TubeRun,
    frame: &level::tube::TubeFrame,
    raw: &level::tube::TubeFrame,
    cathode: usize,
    snapshot: &crate::em::EmSnapshot,
    magnetic: bool,
    scale: f64,
    decades: f64,
    glow: bool,
) -> (PotentialParams, Vec<[f32; 4]>) {
    let (mut params, mut items) = tube_params(
        run, frame, raw, cathode, -1.0, false, 1.0, 1.0, decades, glow,
    );
    items.extend(snapshot.samples.iter().map(|s| [s[0], s[1], s[2], 0.0]));
    params.mode = 3;
    params.grid = Vec4::new(
        snapshot.origin.0 as f32,
        snapshot.origin.1 as f32,
        snapshot.spacing as f32,
        scale as f32,
    );
    params.grid_n = UVec4::new(
        u32::try_from(snapshot.nx).unwrap_or(0),
        u32::try_from(snapshot.ny).unwrap_or(0),
        u32::from(!magnetic),
        (decades as f32).to_bits(),
    );
    (params, items)
}
