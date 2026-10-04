//! Property tests T7 (speed limit) and T8 (plane of symmetry), PHYSICS.md §9.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

mod common;

use common::cube;
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Shape, Sphere};
use physics::trajectory::{RunSettings, Scenario, run};
use proptest::prelude::*;

fn scenario(charges: Vec<FixedCharge>, c: f64, x0: DVec3, p0: DVec3) -> Scenario<Coulomb> {
    Scenario {
        field: Coulomb::new(&charges),
        obstacles: charges
            .iter()
            .map(|q| {
                Shape::Sphere(Sphere {
                    center: q.position,
                    radius: q.radius,
                })
            })
            .collect(),
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0,
        p0,
        detector: None,
        bounds: Some(cube(50.0)),
        t_max: 20.0,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    }
}

fn charge_in(z_range: impl Strategy<Value = f64>) -> impl Strategy<Value = FixedCharge> {
    (
        -10.0..10.0f64,
        -10.0..10.0f64,
        z_range,
        prop_oneof![-100.0..-1.0f64, 1.0..100.0f64],
    )
        .prop_map(|(x, y, z, q)| FixedCharge {
            position: DVec3::new(x, y, z),
            charge: q,
            radius: 0.3,
        })
}

/// Momentum magnitude, log-uniform between 1e-3 and 1e7 (γ up to 1e7 for m = c = 1).
fn momentum_magnitude() -> impl Strategy<Value = f64> {
    (-3.0..7.0f64).prop_map(|e| 10f64.powf(e))
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 64,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// T7: the speed never exceeds c. The state is (x, p), so |v| = |p|c/sqrt(m²c² + p²)
    /// cannot exceed c; strictly below c it holds up to γ ≈ 1e7, beyond which the
    /// floating-point quotient can round to exactly c.
    #[test]
    fn t7_speed_never_exceeds_c(
        charges in prop::collection::vec(charge_in(-3.0..=3.0f64), 1..10),
        c in 1.0..10.0f64,
        p in momentum_magnitude(),
        dir in (-1.0..1.0f64, -1.0..1.0f64, -1.0..1.0f64),
    ) {
        let d = DVec3::new(dir.0, dir.1, dir.2);
        prop_assume!(d.length() > 0.1);
        let scn = scenario(charges, c, DVec3::new(-20.0, 0.5, 0.2), d.normalize() * p);
        let tr = run(&scn, &RunSettings::with_tolerance(1e-9));
        let kin = Kinematics::new(1.0, c);
        for s in &tr.samples {
            let v = kin.velocity(s.p).length();
            prop_assert!(s.p.is_finite() && s.x.is_finite());
            prop_assert!(v <= c, "v = {v}, c = {c}");
            if kin.gamma(s.p) < 1e7 {
                prop_assert!(v < c, "v = {v}, c = {c}, γ = {}", kin.gamma(s.p));
            }
        }
    }

    /// T8: with all charges in the plane z = 0 and the particle starting in it with
    /// in-plane momentum, z and p_z stay exactly zero.
    #[test]
    fn t8_symmetry_plane_is_preserved_exactly(
        charges in prop::collection::vec(charge_in(Just(0.0)), 1..10),
        c in prop_oneof![Just(f64::INFINITY), 1.0..10.0f64],
        p in momentum_magnitude(),
        angle in 0.0..std::f64::consts::TAU,
    ) {
        let p0 = DVec3::new(angle.cos(), angle.sin(), 0.0) * p;
        let scn = scenario(charges, c, DVec3::new(-20.0, 0.5, 0.0), p0);
        let tr = run(&scn, &RunSettings::with_tolerance(1e-9));
        for s in tr.samples.iter().chain([&tr.end]) {
            // ±0 are both exactly zero.
            prop_assert_eq!(s.x.z.abs().to_bits(), 0, "z = {:e}", s.x.z);
            prop_assert_eq!(s.p.z.abs().to_bits(), 0, "p_z = {:e}", s.p.z);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 32,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// T8 with magnets (M5): dipoles perpendicular to the plane and in-plane coils give a
    /// field exactly perpendicular to the plane there, so z and p_z stay exactly zero.
    #[test]
    fn t8_symmetry_plane_is_preserved_with_magnets(
        charges in prop::collection::vec(charge_in(Just(0.0)), 0..5),
        moments in prop::collection::vec((-10.0..10.0f64, -10.0..10.0f64, -50.0..50.0f64), 1..4),
        kappa in -3.0..3.0f64,
        c in prop_oneof![Just(f64::INFINITY), 1.0..10.0f64],
        p in momentum_magnitude(),
        angle in 0.0..std::f64::consts::TAU,
    ) {
        use physics::field::LevelField;
        use physics::magnetic::{CircularLoop, MagneticDipole};
        let dipoles: Vec<MagneticDipole> = moments
            .iter()
            .map(|&(x, y, m)| MagneticDipole {
                position: DVec3::new(x, y, 0.0),
                moment: DVec3::new(0.0, 0.0, m),
                radius: 0.3,
            })
            .collect();
        let mut obstacles: Vec<Shape> = charges
            .iter()
            .map(|q| Shape::Sphere(Sphere { center: q.position, radius: q.radius }))
            .collect();
        obstacles.extend(dipoles.iter().map(|d| Shape::Sphere(Sphere { center: d.position, radius: d.radius })));
        let field = LevelField {
            coulomb: Coulomb::new(&charges),
            dipoles,
            loops: vec![CircularLoop {
                center: DVec3::new(0.0, 0.0, 0.0),
                normal: DVec3::Z,
                radius: 40.0,
                kappa,
                wire_radius: 0.1,
                rate: 0.0,
            }],
            polygons: vec![],
            ..LevelField::default()
        };
        let scn = Scenario {
            field,
            obstacles,
            particle: Particle { charge: 1.0, mass: 1.0, radius: 0.0, moment: 0.0 },
            c,
            x0: DVec3::new(-20.0, 0.5, 0.0),
            p0: DVec3::new(angle.cos(), angle.sin(), 0.0) * p,
            detector: None,
            bounds: Some(cube(50.0)),
            t_max: 20.0,
            radiation_reaction: false,
            acceptance: None,
            gates: Vec::new(),
        };
        let tr = run(&scn, &RunSettings::with_tolerance(1e-9));
        for s in tr.samples.iter().chain([&tr.end]) {
            prop_assert_eq!(s.x.z.abs().to_bits(), 0, "z = {:e}", s.x.z);
            prop_assert_eq!(s.p.z.abs().to_bits(), 0, "p_z = {:e}", s.p.z);
        }
    }
}
