//! Validation tests P1–P3 for the field energy and its flow (PHYSICS.md §10, §9). Run with
//! `cargo test --release -p physics --test poynting -- --nocapture`.

mod common;

use common::sphere_quadrature;
use physics::DVec3;
use physics::lienard::{Worldline, fields, fields_from};
use physics::poynting::{energy_density, energy_velocity, exchange, poynting};

/// A charge circling at radius 1 (ω = 1.2) while oscillating across the plane
/// (`z = 0.3 sin 2.1t`): β up to 0.68 at c = 2, as in test L4.
struct Bobbing;

impl Worldline for Bobbing {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        let (s, c) = (1.2 * t).sin_cos();
        let (sz, cz) = (2.1 * t).sin_cos();
        (
            DVec3::new(c, s, 0.3 * sz),
            DVec3::new(-1.2 * s, 1.2 * c, 0.63 * cz),
            DVec3::new(-1.44 * c, -1.44 * s, -1.323 * sz),
        )
    }
}

struct Uniform {
    x0: DVec3,
    v: DVec3,
}

impl Worldline for Uniform {
    fn state(&self, t: f64) -> (DVec3, DVec3, DVec3) {
        (self.x0 + self.v * t, self.v, DVec3::ZERO)
    }
}

// --- P1: Poynting's theorem ------------------------------------------------------------------

/// P1, Poynting's theorem (Jackson §6.7): away from charges `∂u/∂t + ∇·S = 0`, for the
/// retarded field of a charge in relativistic three-dimensional motion (`Bobbing`), and for
/// its exchange terms with a static field (a fixed charge's Coulomb field, a uniform E and a
/// uniform B, each a vacuum solution away from its sources). Central differences over 1e-4
/// in space and time, at points 1.6 to 50 from the charge, off the plane. Required: the
/// residual below 1e-6 of the sum of the terms' magnitudes (set before measuring: the
/// differences' own error is about 1e-8, as in test L3).
#[test]
fn p1_poynting_theorem() {
    let (q, c) = (1.0, 2.0);
    let (qf, xf) = (0.7, DVec3::new(2.5, -1.5, 0.4));
    let (e0, b0) = (DVec3::new(0.05, 0.02, -0.03), DVec3::new(0.1, -0.2, 0.3));
    // (u, S) of the charge's own field and of its exchange terms with the static field.
    let parts = |x: DVec3, t: f64| {
        let f = fields(&Bobbing, q, c, x, t);
        let (e, b) = (f.e(), f.b);
        let d = x - xf;
        let (ee, be) = (e0 + d * (qf / d.length().powi(3)), b0);
        [
            (energy_density(e, b, c), poynting(e, b, c)),
            exchange(e, b, ee, be, c),
        ]
    };
    let h = 1e-4;
    let axes = [DVec3::X, DVec3::Y, DVec3::Z];
    let mut worst = [0.0_f64; 2];
    for (x, t) in [
        (DVec3::new(3.0, 0.5, 0.2), 1.1),
        (DVec3::new(-6.0, 4.0, -2.0), 1.1),
        (DVec3::new(30.0, -25.0, 10.0), 1.1),
        (DVec3::new(0.5, -1.5, 1.0), 2.3),
        (DVec3::new(0.0, 0.0, 50.0), 2.3),
    ] {
        let (later, earlier) = (parts(x, t + h), parts(x, t - h));
        let plus: Vec<_> = axes.iter().map(|&a| parts(x + a * h, t)).collect();
        let minus: Vec<_> = axes.iter().map(|&a| parts(x - a * h, t)).collect();
        for (k, name) in ["own", "exchange"].iter().enumerate() {
            let dudt = (later[k].0 - earlier[k].0) / (2.0 * h);
            let terms: Vec<f64> = (0..3)
                .map(|i| (plus[i][k].1[i] - minus[i][k].1[i]) / (2.0 * h))
                .collect();
            let div: f64 = terms.iter().sum();
            let scale = dudt.abs() + terms.iter().map(|v| v.abs()).sum::<f64>();
            let res = (dudt + div).abs() / scale;
            println!("P1 at {x}, t = {t}: {name} residual {res:.2e} (∂u/∂t = {dudt:.3e})");
            worst[k] = worst[k].max(res);
        }
    }
    assert!(worst[0] < 1e-6 && worst[1] < 1e-6, "P1 {worst:?}");
}

// --- P2: radiated power ----------------------------------------------------------------------

/// P2, the radiated power (Jackson 14.38): the radiation part's flow through a sphere of
/// radius R around the retarded position, weighted by `κ = 1 − n·β` (`dt/dt′`: the energy
/// crossing the sphere was emitted over a shorter time), is Liénard's power at the retarded
/// time, `(2q²/3c³) γ⁶ (a² − |v × a|²/c²)`, for any R. Three states (β = 0.3, 0.6, 0.9; the
/// acceleration across, oblique to and nearly along the motion), R = 1 and 25; the sphere
/// quadrature has its polar axis along v (48 × 96 nodes). Required: 1e-10 (set before
/// measuring). Also: the energy velocity `S/u` of the radiation field is c, outwards.
#[test]
fn p2_radiated_power() {
    let (q, c) = (1.0, 2.0);
    let r = DVec3::new(0.3, -0.2, 0.1);
    let states = [
        (DVec3::new(0.6, 0.0, 0.0), DVec3::new(0.0, 0.7, 0.0)),
        (
            DVec3::new(1.2 / 2f64.sqrt(), 1.2 / 2f64.sqrt(), 0.0),
            DVec3::new(0.4, -0.2, 0.3),
        ),
        (DVec3::new(0.0, 1.8, 0.0), DVec3::new(0.1, 0.5, 0.05)),
    ];
    let mut worst: f64 = 0.0;
    let mut worst_speed: f64 = 0.0;
    for (v, a) in states {
        let gamma = 1.0 / (1.0 - v.length_squared() / (c * c)).sqrt();
        let lienard = 2.0 * q * q / (3.0 * c.powi(3))
            * gamma.powi(6)
            * (a.length_squared() - v.cross(a).length_squared() / (c * c));
        for radius in [1.0, 25.0] {
            let mut flow = 0.0;
            for (n, w) in sphere_quadrature(v, 48, 96) {
                let f = fields_from(q, c, r + n * radius, 0.0, r, v, a);
                let b = n.cross(f.e_radiation) / c;
                let kappa = 1.0 - n.dot(v) / c;
                flow += w * kappa * radius * radius * poynting(f.e_radiation, b, c).dot(n);
                if f.e_radiation.length() > 0.0 {
                    let speed = energy_velocity(f.e_radiation, b, c);
                    worst_speed = worst_speed.max((speed - n * c).length() / c);
                }
            }
            let err = flow / lienard - 1.0;
            println!(
                "P2 β = {:.1}, R = {radius}: flow {flow:.12e}, Liénard {lienard:.12e}: {err:.1e}",
                v.length() / c
            );
            worst = worst.max(err.abs());
        }
    }
    println!("P2 energy velocity of the radiation field: c n to {worst_speed:.1e}");
    assert!(worst < 1e-10, "P2 {worst:.3e}");
    assert!(worst_speed < 1e-12, "P2 energy velocity {worst_speed:.3e}");
}

// --- P3: where a particle's energy comes from ---------------------------------------------

/// P3, the work on a particle in the exchange terms: a charge in uniform motion (β = v/c)
/// through a uniform external E (and a uniform B, which does no work) receives the power
/// `q v·E`. The exchange flow into a sphere centred on its present position brings only
/// `q v·E [1 − (1 − β²)(artanh β − β)/β³]`: 2/3 of it at low speed, all of it as β → 1. The
/// rest comes from the exchange energy stored inside the sphere, which falls as the charge
/// moves off the centre at the rate `(1/4π) E·∮ E_p (v·n) dA = q v·E (1 − β²)(artanh β − β)/β³`
/// (the Heaviside field; the E_p × B term has no flux through the sphere, since that field
/// points away from the present position). β = 0.1, 0.5, 0.9; radii 0.5 and 5; the field
/// along, across and oblique to the motion. Required: 1e-10 of `q |v| |E|` (set before
/// measuring).
#[test]
fn p3_work_on_a_particle() {
    let (q, c) = (1.0, 2.0);
    let b_ext = DVec3::new(0.2, -0.1, 0.5);
    let dir = DVec3::new(0.6, 0.8, 0.0);
    let mut worst: f64 = 0.0;
    for beta in [0.1_f64, 0.5, 0.9] {
        let v = dir * (beta * c);
        let line = Uniform {
            x0: DVec3::new(0.3, -0.2, 0.1),
            v,
        };
        let t = 0.7;
        let p = line.state(t).0;
        let stored = (1.0 - beta * beta) * (beta.atanh() - beta) / beta.powi(3);
        for e_ext in [
            dir * 0.4,
            DVec3::new(-0.32, 0.24, 0.3),
            DVec3::new(0.3, 0.2, -0.25),
        ] {
            let expected = q * v.dot(e_ext) * (1.0 - stored);
            for radius in [0.5, 5.0] {
                let mut inflow = 0.0;
                for (n, w) in sphere_quadrature(v, 48, 96) {
                    let x = p + n * radius;
                    let f = fields(&line, q, c, x, t);
                    let (_, s) = exchange(f.e(), f.b, e_ext, b_ext, c);
                    inflow -= w * radius * radius * s.dot(n);
                }
                let err = (inflow - expected).abs() / (q * v.length() * e_ext.length());
                println!(
                    "P3 β = {beta}, E = {e_ext}, radius {radius}: inflow {inflow:.12e}, \
                     expected {expected:.12e} (work {:.6e}): {err:.1e}",
                    q * v.dot(e_ext)
                );
                worst = worst.max(err);
            }
        }
    }
    assert!(worst < 1e-10, "P3 {worst:.3e}");
}
