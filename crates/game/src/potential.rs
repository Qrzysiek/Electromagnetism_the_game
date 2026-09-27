//! GPU field maps (potential or magnetic field): a quad covering the world bounds with a
//! fragment shader that evaluates the field per pixel (visual only).

use bevy::asset::{AssetPath, embedded_asset, embedded_path};
use bevy::prelude::*;
use bevy::render::render_resource::{AsBindGroup, ShaderType};
use bevy::shader::ShaderRef;
use bevy::sprite_render::{Material2d, Material2dPlugin};
use level::beam::Launch;
use physics::field::{FieldSolver, LevelField};
use physics::trajectory::Scenario;

/// Array sizes; must match `potential.wgsl`.
pub const MAX_CHARGES: usize = 1024;
pub const MAX_MAGNETS: usize = 64;
pub const MAX_LOOPS: usize = 16;
pub const MAX_SEGMENTS: usize = 64;
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
}

#[derive(ShaderType, Debug, Clone, Copy)]
pub struct PotentialParams {
    pub charges: [Vec4; MAX_CHARGES],
    pub magnets: [Vec4; MAX_MAGNETS],
    pub loops: [Vec4; MAX_LOOPS],
    pub segments: [Vec4; MAX_SEGMENTS],
    pub segment_kappa: [Vec4; MAX_SEGMENTS],
    pub limits: [Vec4; MAX_LIMITS],
    pub counts: UVec4,
    pub limit_count: u32,
    pub phi_weight: f32,
    pub u_a: f32,
    pub solid: u32,
    pub mode: u32,
    pub wire: f32,
    pub moment_weight: f32,
}

impl Default for PotentialParams {
    fn default() -> Self {
        Self {
            charges: [Vec4::ZERO; MAX_CHARGES],
            magnets: [Vec4::ZERO; MAX_MAGNETS],
            loops: [Vec4::ZERO; MAX_LOOPS],
            segments: [Vec4::ZERO; MAX_SEGMENTS],
            segment_kappa: [Vec4::ZERO; MAX_SEGMENTS],
            limits: [Vec4::ZERO; MAX_LIMITS],
            counts: UVec4::ZERO,
            limit_count: 0,
            phi_weight: 0.0,
            u_a: 0.0,
            solid: 0,
            mode: 0,
            wire: 0.1,
            moment_weight: 0.0,
        }
    }
}

#[derive(Asset, TypePath, AsBindGroup, Debug, Clone)]
pub struct PotentialMaterial {
    #[uniform(0)]
    pub params: PotentialParams,
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

fn count(n: usize) -> u32 {
    u32::try_from(n).expect("fits")
}

/// Field for which the shot's particle circles with this gyroradius (cells): the unit of
/// the magnetic map.
pub const GYRO_REFERENCE: f64 = 5.0;

/// Shader parameters for one shot. The shader sums `Φ = Σ Q/r` over the charges and the
/// field `B_z / b_ref` of the magnetic sources (pre-scaled by `1/b_ref`, with
/// `b_ref = |p₀| / (|q| r_ref)`, or for a neutral particle with a magnetic moment `m`,
/// `b_ref = T₀ / |m|`); the colours show the shot's `U/T₀ = (q Φ − m B_z)/T₀`.
///
/// The dark region is where every particle in `limits` is forbidden by energy
/// conservation: each group (one shot, or one beam) contributes the highest total energy
/// `E = T + U(x₀)` of its particles, and a point is dark only when `U > E` there for
/// every group. More than `MAX_LIMITS` groups: no dark region (never too large).
#[allow(clippy::cast_possible_truncation)]
pub fn params(
    scn: &Scenario<LevelField>,
    limits: &[Vec<Launch>],
    charge_radius: f64,
    magnet_radius: f64,
    mode: MapMode,
) -> PotentialParams {
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let t0 = kin.kinetic_energy(scn.p0).max(1e-300);
    let (q, m) = (scn.particle.charge, scn.particle.moment);
    let b_ref = if q != 0.0 {
        (scn.p0.length() / (q.abs() * GYRO_REFERENCE)).max(1e-300)
    } else if m != 0.0 {
        t0 / m.abs()
    } else {
        1.0
    };
    let f = &scn.field;
    let at_a = f.sample(scn.x0, 0.0);
    let mut out = PotentialParams {
        u_a: ((q * at_a.phi - m * at_a.b.z) / t0) as f32,
        phi_weight: (q / t0) as f32,
        moment_weight: (-m * b_ref / t0) as f32,
        mode: match mode {
            MapMode::Potential => 0,
            MapMode::Magnetic | MapMode::Waves | MapMode::ParticleField | MapMode::Total => 1,
        },
        wire: f
            .loops
            .first()
            .map(|l| l.wire_radius)
            .or_else(|| f.polygons.first().map(|p| p.wire_radius))
            .unwrap_or(0.1) as f32,
        ..PotentialParams::default()
    };
    // Fixed charges (solid, in the plane), then the charges induced on metal (display
    // resolution): sphere images and equivalent charges, and electrode panels as point
    // charges at their centroids and mirror images (with their z).
    let fixed: Vec<(physics::DVec3, f64)> = f.coulomb.charges().collect();
    let mut induced: Vec<(physics::DVec3, f64)> = f.conductors.induced.charges().collect();
    for (t, s) in f.electrodes.panels() {
        let c = t.centroid();
        let q_panel = s * t.area();
        induced.push((c, q_panel));
        induced.push((physics::DVec3::new(c.x, c.y, -c.z), q_panel));
    }
    let mut n = 0;
    for (pos, qc) in fixed.iter().take(MAX_CHARGES) {
        out.charges[n] = Vec4::new(pos.x as f32, pos.y as f32, *qc as f32, charge_radius as f32);
        n += 1;
    }
    out.solid = count(n);
    for (pos, qc) in induced.iter().take(MAX_CHARGES - n) {
        out.charges[n] = Vec4::new(pos.x as f32, pos.y as f32, *qc as f32, pos.z as f32);
        n += 1;
    }
    out.counts.x = count(n);
    for (i, d) in f.dipoles.iter().take(MAX_MAGNETS).enumerate() {
        out.magnets[i] = Vec4::new(
            d.position.x as f32,
            d.position.y as f32,
            (d.moment.z / b_ref) as f32,
            magnet_radius as f32,
        );
        out.counts.y = count(i + 1);
    }
    for (i, l) in f.loops.iter().take(MAX_LOOPS).enumerate() {
        out.loops[i] = Vec4::new(
            l.center.x as f32,
            l.center.y as f32,
            l.radius as f32,
            (l.kappa / b_ref) as f32,
        );
        out.counts.z = count(i + 1);
    }
    let mut n = 0;
    for poly in &f.polygons {
        let v = &poly.vertices;
        for i in 0..v.len() {
            if n == MAX_SEGMENTS {
                break;
            }
            let (a, b) = (v[i], v[(i + 1) % v.len()]);
            out.segments[n] = Vec4::new(a.x as f32, a.y as f32, b.x as f32, b.y as f32);
            out.segment_kappa[n] = Vec4::new((poly.kappa / b_ref) as f32, 0.0, 0.0, 0.0);
            n += 1;
        }
    }
    out.counts.w = count(n);
    if limits.len() <= MAX_LIMITS {
        for (i, group) in limits.iter().enumerate() {
            // Highest total energy of the group, in units of its highest kinetic energy.
            let (mut e, mut t) = (f64::NEG_INFINITY, 1e-300_f64);
            let (mut qg, mut mg) = (0.0, 0.0);
            for l in group {
                let kin = physics::dynamics::Kinematics::new(l.particle.mass, scn.c);
                let tk = kin.kinetic_energy(l.p0);
                let at = f.sample(l.x0, 0.0);
                (qg, mg) = (l.particle.charge, l.particle.moment);
                e = e.max(tk + qg * at.phi - mg * at.b.z);
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
    out
}
