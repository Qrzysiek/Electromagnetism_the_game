//! GPU potential map: a quad covering the world bounds with a fragment shader that
//! evaluates the potential per pixel (visual only).

use bevy::asset::{AssetPath, embedded_asset, embedded_path};
use bevy::prelude::*;
use bevy::render::render_resource::{AsBindGroup, ShaderType};
use bevy::shader::ShaderRef;
use bevy::sprite_render::{Material2d, Material2dPlugin};
use physics::DVec3;
use physics::field::{Coulomb, FieldSolver};
use physics::trajectory::Scenario;

/// Must match the array size in `potential.wgsl`.
pub const MAX_CHARGES: usize = 256;

#[derive(ShaderType, Debug, Clone, Copy)]
pub struct PotentialParams {
    pub charges: [Vec4; MAX_CHARGES],
    pub count: u32,
    pub u_a: f32,
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

/// Shader parameters for a scenario: charges pre-scaled to `w = qQ/T₀` so the shader
/// computes `U/T₀` directly. `charges` are (position, Q) in field summation order.
#[allow(clippy::cast_possible_truncation)]
pub fn params(scn: &Scenario<Coulomb>, charges: &[(DVec3, f64)]) -> PotentialParams {
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let t0 = kin.kinetic_energy(scn.p0).max(1e-300);
    let q = scn.particle.charge;
    let mut out = PotentialParams {
        charges: [Vec4::ZERO; MAX_CHARGES],
        count: 0,
        u_a: (q * scn.field.sample(scn.x0, 0.0).phi / t0) as f32,
    };
    for (i, (&(pos, qc), sphere)) in charges
        .iter()
        .zip(&scn.obstacles)
        .take(MAX_CHARGES)
        .enumerate()
    {
        out.charges[i] = Vec4::new(
            pos.x as f32,
            pos.y as f32,
            (q * qc / t0) as f32,
            sphere.radius as f32,
        );
        out.count = u32::try_from(i + 1).expect("fits");
    }
    out
}
