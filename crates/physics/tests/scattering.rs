//! Validation test MS1: multiple Coulomb scattering in the plane (PHYSICS.md §9). Run with
//! `cargo test --release -p physics --test scattering -- --nocapture`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use std::f64::consts::PI;

use common::Rng;
use physics::DVec3;
use physics::dynamics::Particle;
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{Outcome, RunSettings, Scenario, run};
use rayon::prelude::*;

/// MS1, Jackson §13.6–13.7 in the plane. A particle crossing a strip of point charges that
/// lie in its own plane (areal density n, random signs, thickness L) is deflected by
/// `θ = K/b` (K = 2qQ/pv) by each, with the impact parameters b uniform on the line:
/// single deflections have a `1/θ²` tail, and their sum is a Cauchy distribution of
/// half-width `Γ = π n L K` (its characteristic function is `exp(−π n L K |k|)`), growing
/// linearly with the thickness. A foil in 3D (Jackson's multiple scattering) weights the
/// impact parameters by `b db` instead: a Gaussian of width ∝ √L (times a logarithm). The
/// slice cannot hold a foil. Here 16 random strips of each thickness (L = 2, 4; n = 2,
/// |Q| = 1e-3, radius 4e-3, height 20), each crossed by 250 particles (q = m = v = 1,
/// Newtonian) from 20 cells before to 20 after; particles that hit a charge count as the
/// large deflections they are. Required (set before measuring; the medians' statistical
/// error is about 4 %): the median |θ| within 15 % of Γ, and the ratio of the medians of
/// L = 4 and L = 2 within 1.8–2.2 (Cauchy 2; a Gaussian gives 1.41).
#[test]
fn ms1_multiple_scattering_in_the_plane_is_cauchy() {
    let (q, m, v) = (1.0, 1.0, 1.0);
    let (big_q, density, radius, height) = (1e-3, 2.0, 4e-3, 20.0);
    let k = 2.0 * q * big_q / (m * v * v);
    let mut medians = Vec::new();
    for thickness in [2.0, 4.0] {
        let gamma = PI * density * thickness * k;
        let mut angles: Vec<f64> = Vec::new();
        let mut lost = 0usize;
        for seed in 0..16u64 {
            let mut rng = Rng::new(1000 + seed);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let count = (density * thickness * height).round() as usize;
            let charges: Vec<FixedCharge> = (0..count)
                .map(|_| FixedCharge {
                    position: DVec3::new(
                        rng.range(0.0, thickness),
                        rng.range(-0.5 * height, 0.5 * height),
                        0.0,
                    ),
                    charge: if rng.next_f64() < 0.5 { big_q } else { -big_q },
                    radius,
                })
                .collect();
            let field = Coulomb::new(&charges);
            let obstacles: Vec<Shape> = charges
                .iter()
                .map(|c| {
                    Shape::Sphere(Sphere {
                        center: c.position,
                        radius,
                    })
                })
                .collect();
            let flights: Vec<Option<f64>> = (0..250)
                .into_par_iter()
                .map(|i| {
                    let y = -5.0 + 10.0 * (f64::from(i) + 0.5) / 250.0;
                    let tr = run(
                        &Scenario {
                            field: &field,
                            obstacles: obstacles.clone(),
                            particle: Particle {
                                charge: q,
                                mass: m,
                                radius: 0.0,
                                moment: 0.0,
                            },
                            c: f64::INFINITY,
                            x0: DVec3::new(-20.0, y, 0.0),
                            p0: DVec3::new(m * v, 0.0, 0.0),
                            detector: None,
                            bounds: None,
                            t_max: (thickness + 40.0) / v,
                            radiation_reaction: false,
                            acceptance: None,
                            gates: Vec::new(),
                        },
                        &RunSettings::with_tolerance(1e-10),
                    );
                    matches!(tr.outcome, Outcome::Timeout)
                        .then(|| tr.end.p.y.atan2(tr.end.p.x).abs())
                })
                .collect();
            for f in flights {
                match f {
                    Some(a) => angles.push(a),
                    None => lost += 1,
                }
            }
        }
        // A particle that hits a charge passed it closer than its radius: a deflection
        // beyond K/radius = 0.5 rad, far in the tail. Counted there (dropping them would
        // bias the median low by 5–9 %).
        angles.sort_by(f64::total_cmp);
        let median = angles[(angles.len() + lost) / 2];
        // The Gaussian's width with the same variance would not describe the core: the
        // variance is set by the closest encounters (∝ 1/b_min).
        println!(
            "MS1 L = {thickness}: median |θ| = {median:.5} ({} flights, {lost} absorbed), \
             Cauchy half-width πnLK = {gamma:.5}: {:.1e}",
            angles.len(),
            median / gamma - 1.0
        );
        assert!((median / gamma - 1.0).abs() < 0.15, "MS1 L = {thickness}");
        medians.push(median);
    }
    let ratio = medians[1] / medians[0];
    println!("MS1 median ratio for twice the thickness: {ratio:.3} (Cauchy 2, Gaussian 1.414)");
    assert!((1.8..=2.2).contains(&ratio), "MS1 ratio {ratio:.3}");
}
