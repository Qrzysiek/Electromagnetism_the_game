//! Frame-budget benchmark (SPEC §10, Stage 1): one trajectory through a level with 50
//! fixed charges, at preview and verification tolerances.

use criterion::{Criterion, criterion_group, criterion_main};
use physics::DVec3;
use physics::dynamics::Particle;
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Aabb, Region, Shape, Sphere};
use physics::trajectory::{RunSettings, Scenario, run};

/// 50 charges on integer nodes of a 40 × 30 grid (2D slice), deterministic layout.
fn level_50() -> Scenario<Coulomb> {
    let mut state = 0x1234_5678_u64;
    let mut next = move |n: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % n
    };
    let mut charges = Vec::new();
    while charges.len() < 50 {
        #[allow(clippy::cast_precision_loss)]
        let pos = DVec3::new(next(34) as f64 + 3.0, next(30) as f64, 0.0);
        if (pos - DVec3::new(0.0, 15.0, 0.0)).length() < 4.0
            || charges.iter().any(|c: &FixedCharge| c.position == pos)
        {
            continue;
        }
        let sign = if next(2) == 0 { -1.0 } else { 1.0 };
        #[allow(clippy::cast_precision_loss)]
        let magnitude = (next(3) + 1) as f64 * 0.05;
        charges.push(FixedCharge {
            position: pos,
            charge: sign * magnitude,
            radius: 0.25,
        });
    }
    Scenario {
        obstacles: charges
            .iter()
            .map(|c| {
                Shape::Sphere(Sphere {
                    center: c.position,
                    radius: c.radius,
                })
            })
            .collect(),
        field: Coulomb::new(&charges),
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
        },
        c: 100.0,
        x0: DVec3::new(0.0, 15.0, 0.0),
        p0: DVec3::new(1.0, 0.1, 0.0),
        detector: Some(Region::Box(Aabb {
            min: DVec3::new(38.0, 12.0, -1.0),
            max: DVec3::new(40.0, 18.0, 1.0),
        })),
        bounds: Some(Aabb {
            min: DVec3::new(-1.0, -1.0, -1.0),
            max: DVec3::new(41.0, 31.0, 1.0),
        }),
        t_max: 400.0,
        radiation_reaction: false,
    }
}

fn bench(c: &mut Criterion) {
    let scn = level_50();
    for tol in [1e-8, 1e-11] {
        let rs = RunSettings::with_tolerance(tol);
        let tr = run(&scn, &rs);
        println!(
            "tol {tol:.0e}: outcome {:?} at t = {:.2}, {} steps, {} field evaluations",
            tr.outcome, tr.end.t, tr.stats.n_accept, tr.stats.n_fcn
        );
        c.bench_function(&format!("trajectory_50_charges_tol_{tol:.0e}"), |b| {
            b.iter(|| run(&scn, &rs))
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
