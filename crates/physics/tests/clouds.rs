//! Charge clouds: uniformly charged spheres particles can fly through (Thomson's atom,
//! PHYSICS.md §2.1). Run with `cargo test --release -p physics --test clouds -- --nocapture`.

mod common;

use std::f64::consts::PI;

use common::cube;
use physics::DVec3;
use physics::dynamics::{Kinematics, Particle};
use physics::field::{ChargeCloud, Coulomb, FieldSolver, FixedCharge};
use physics::trajectory::{RunSettings, Scenario, run};

const TOL: f64 = 1e-12;

/// C1: the cloud's field is minus the gradient of its potential, obeys Gauss's law
/// inside (div E = 4πρ with ρ = 3Q/(4πR³), i.e. div E = 3Q/R³), and is continuous at the
/// surface; outside it equals a point charge's. A point charge beside it is unchanged.
#[test]
fn c1_cloud_field_potential_gauss_and_continuity() {
    let cloud = ChargeCloud {
        position: DVec3::new(1.0, -2.0, 0.0),
        charge: 2.5,
        radius: 3.0,
    };
    let point = FixedCharge {
        position: DVec3::new(9.0, 4.0, 0.0),
        charge: -1.2,
        radius: 0.3,
    };
    let f = Coulomb::with_clouds(&[point], &[cloud]);
    let h = 1e-5;
    let mut worst_grad: f64 = 0.0;
    let mut worst_div: f64 = 0.0;
    for x in [
        DVec3::new(1.5, -1.0, 0.3),
        DVec3::new(-0.7, -3.1, -1.0),
        DVec3::new(6.0, 1.0, 0.0),
    ] {
        let s = f.sample(x, 0.0);
        let phi = |y: DVec3| f.sample(y, 0.0).phi;
        let e_of = |y: DVec3| f.sample(y, 0.0).e;
        let grad = DVec3::new(
            phi(x + DVec3::X * h) - phi(x - DVec3::X * h),
            phi(x + DVec3::Y * h) - phi(x - DVec3::Y * h),
            phi(x + DVec3::Z * h) - phi(x - DVec3::Z * h),
        ) / (2.0 * h);
        worst_grad = worst_grad.max((s.e + grad).length() / s.e.length());
        let div = (e_of(x + DVec3::X * h).x - e_of(x - DVec3::X * h).x + e_of(x + DVec3::Y * h).y
            - e_of(x - DVec3::Y * h).y
            + e_of(x + DVec3::Z * h).z
            - e_of(x - DVec3::Z * h).z)
            / (2.0 * h);
        let inside = (x - cloud.position).length() < cloud.radius;
        let expected = if inside {
            3.0 * cloud.charge / cloud.radius.powi(3)
        } else {
            0.0
        };
        worst_div =
            worst_div.max((div - expected).abs() / (3.0 * cloud.charge / cloud.radius.powi(3)));
    }
    // Continuity across the surface (field and potential), and the point-charge value
    // outside.
    let dir = DVec3::new(0.6, 0.8, 0.0);
    let at = |r: f64| f.sample(cloud.position + dir * r, 0.0);
    let (a, b) = (
        at(cloud.radius * (1.0 - 1e-12)),
        at(cloud.radius * (1.0 + 1e-12)),
    );
    let jump = (a.e - b.e).length() / b.e.length() + (a.phi - b.phi).abs() / b.phi.abs();
    let only_cloud = Coulomb::with_clouds(&[], &[cloud]);
    let far = cloud.position + dir * 7.0;
    let pc = Coulomb::new(&[FixedCharge {
        position: cloud.position,
        charge: cloud.charge,
        radius: 0.0,
    }]);
    let outside = (only_cloud.sample(far, 0.0).e - pc.sample(far, 0.0).e).length();
    println!(
        "C1: |E + ∇φ|/|E| ≤ {worst_grad:.1e}; Gauss |div E − 4πρ| / (3Q/R³) ≤ {worst_div:.1e}; \
         surface jump {jump:.1e}; outside vs point charge {outside:.1e} (4π = {:.4})",
        4.0 * PI
    );
    assert!(worst_grad < 1e-8 && worst_div < 1e-5 && jump < 1e-10 && outside == 0.0);
}

/// C2, Jackson Problem 16.1 (and §16.7): a charge bound harmonically (here inside a
/// charge cloud, ω₀² = |qQ|/(mR³)) radiates, and its oscillation energy decays as
/// `e^{−Γt}` with `Γ = ω₀² τ`, `τ = 2q²/(3mc³)`. Landau–Lifshitz agrees with the
/// Abraham–Lorentz result to O(ω₀τ) (here 0.01): the fitted decay rate over five decay
/// times must agree within 2 %. And exactly: the end state after the five decay times
/// (~600 oscillations) against an independent integration of the full relativistic
/// Landau–Lifshitz equation written from Landau & Lifshitz §76
/// (`scripts/wolfram/c2_radiating_oscillator.wls`, Wolfram Engine 14.2, NDSolve at 30
/// digits; the same to 17 digits with a tighter goal): position and momentum within
/// 1e-6 of the end amplitudes (set before measuring).
#[test]
fn c2_radiating_oscillator_decays_at_gamma() {
    let (q, m, big_q, r, c) = (-1.0_f64, 1.0_f64, 1.0_f64, 4.0_f64, 2.0_f64);
    let cloud = ChargeCloud {
        position: DVec3::ZERO,
        charge: big_q,
        radius: r,
    };
    let omega0 = (q.abs() * big_q / (m * r.powi(3))).sqrt();
    let tau = 2.0 * q * q / (3.0 * m * c.powi(3));
    let gamma = omega0 * omega0 * tau;
    let field = Coulomb::with_clouds(&[], &[cloud]);
    // An in-plane ellipse of amplitude ~1 (well inside the cloud).
    let v0 = 0.6 * omega0;
    let scn = Scenario {
        field,
        obstacles: vec![],
        particle: Particle {
            charge: q,
            mass: m,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0: DVec3::new(1.0, 0.0, 0.0),
        p0: Kinematics::new(m, c).momentum_from_kinetic_energy(0.5 * m * v0 * v0, DVec3::Y),
        detector: None,
        bounds: Some(cube(100.0)),
        t_max: 5.0 / gamma,
        radiation_reaction: true,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(TOL));
    // End state (t = 5/Γ = 3840) from the Wolfram script: x, y, px, py.
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    const END: [f64; 4] = [
        -0.059818079583771145,
        0.032384524846739647,
        -0.0069765044047088572,
        -0.0046503164591507392,
    ];
    let x_ref = DVec3::new(END[0], END[1], 0.0);
    let p_ref = DVec3::new(END[2], END[3], 0.0);
    let ex = (tr.end.x - x_ref).length() / x_ref.length();
    let ep = (tr.end.p - p_ref).length() / p_ref.length();
    println!(
        "C2 end state at t = {:.1}: position {ex:.1e}, momentum {ep:.1e} from the independent LL integration",
        tr.end.t
    );
    assert!((tr.end.t - 3840.0).abs() < 1e-9 && ex < 1e-6 && ep < 1e-6);
    let kin = Kinematics::new(m, c);
    let phi0 = scn.field.sample(DVec3::ZERO, 0.0).phi;
    // ln(oscillation energy) against t: least-squares slope.
    let pts: Vec<(f64, f64)> = tr
        .samples
        .iter()
        .map(|s| {
            let w = kin.kinetic_energy(s.p) + q * (scn.field.sample(s.x, 0.0).phi - phi0);
            (s.t, w.ln())
        })
        .collect();
    #[allow(clippy::cast_precision_loss)]
    let n = pts.len() as f64;
    let (mt, my) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let slope = pts.iter().map(|&(t, y)| (t - mt) * (y - my)).sum::<f64>()
        / pts.iter().map(|&(t, _)| (t - mt).powi(2)).sum::<f64>();
    println!(
        "C2: ω₀ = {omega0:.4}, ω₀τ = {:.4}: fitted decay rate {:.6e}, Jackson Γ = ω₀²τ = {gamma:.6e}, \
         ratio {:.4} ({} samples, max |F_RR|/|F_L| = {:.1e})",
        omega0 * tau,
        -slope,
        -slope / gamma,
        pts.len(),
        tr.reaction_ratio_max
    );
    assert!((-slope / gamma - 1.0).abs() < 0.02);
}
