//! GPU field maps (potential or magnetic field): a quad covering the world bounds with a
//! fragment shader that evaluates the field per pixel (visual only).

use bevy::asset::{AssetPath, embedded_asset, embedded_path};
use bevy::prelude::*;
use bevy::render::render_resource::{AsBindGroup, ShaderType};
use bevy::shader::ShaderRef;
use bevy::sprite_render::{Material2d, Material2dPlugin};
use physics::field::{FieldSolver, LevelField};
use physics::trajectory::Scenario;

/// Array sizes; must match `potential.wgsl`.
pub const MAX_CHARGES: usize = 1024;
pub const MAX_MAGNETS: usize = 64;
pub const MAX_LOOPS: usize = 16;
pub const MAX_SEGMENTS: usize = 64;

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
    pub counts: UVec4,
    pub u_a: f32,
    pub mode: u32,
    pub wire: f32,
}

impl Default for PotentialParams {
    fn default() -> Self {
        Self {
            charges: [Vec4::ZERO; MAX_CHARGES],
            magnets: [Vec4::ZERO; MAX_MAGNETS],
            loops: [Vec4::ZERO; MAX_LOOPS],
            segments: [Vec4::ZERO; MAX_SEGMENTS],
            segment_kappa: [Vec4::ZERO; MAX_SEGMENTS],
            counts: UVec4::ZERO,
            u_a: 0.0,
            mode: 0,
            wire: 0.1,
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

/// Shader parameters for one shot. Charges are pre-scaled to `w = qQ/T₀` (the shader
/// computes `U/T₀`), magnetic sources to `1/b_ref` with `b_ref = |p₀| / (|q| r_ref)`.
#[allow(clippy::cast_possible_truncation)]
pub fn params(
    scn: &Scenario<LevelField>,
    charge_radius: f64,
    magnet_radius: f64,
    mode: MapMode,
) -> PotentialParams {
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let t0 = kin.kinetic_energy(scn.p0).max(1e-300);
    let q = scn.particle.charge;
    let b_ref = (scn.p0.length() / (q.abs() * GYRO_REFERENCE)).max(1e-300);
    let f = &scn.field;
    let mut out = PotentialParams {
        u_a: (q * f.sample(scn.x0, 0.0).phi / t0) as f32,
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
    // Fixed charges and the charges they induce on metal spheres (display resolution).
    let charges = f.coulomb.charges().chain(f.conductors.induced.charges());
    for (i, (pos, qc)) in charges.take(MAX_CHARGES).enumerate() {
        out.charges[i] = Vec4::new(
            pos.x as f32,
            pos.y as f32,
            (q * qc / t0) as f32,
            charge_radius as f32,
        );
        out.counts.x = count(i + 1);
    }
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
    out
}
