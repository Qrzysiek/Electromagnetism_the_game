//! Validation tests Z10–Z13 of circuits driving electrodes and coils (`drive.rs`,
//! PHYSICS.md §2.10): the one-way coupling, against closed forms. Common requirement, set
//! before measuring: within 1e-8 of the scale (the circuit's accuracy, tests Z4–Z9, with a
//! margin for the flights' tolerance); the canonical angular momentum, an exact invariant,
//! within 1e-10 as for a linear ramp (M10).

#![allow(clippy::disallowed_methods)] // references

mod common;

use std::sync::Arc;

use common::cube;
use physics::DVec3;
use physics::bem::{Bias, BoxElectrode, Electrodes};
use physics::circuit::Waveform;
use physics::drive::{Chain, Drives};
use physics::dynamics::Particle;
use physics::field::{FieldSolver, LevelField};
use physics::magnetic::CircularLoop;
use physics::trajectory::{RunSettings, Scenario, run};

/// Preview panels (cells).
const PANEL: f64 = 0.5;

fn plate(y: f64, bias: Bias) -> BoxElectrode {
    BoxElectrode {
        center: DVec3::new(0.0, y, 0.0),
        angle: 0.0,
        half_length: 4.0,
        half_thickness: 0.25,
        half_height: 4.0,
        bias,
    }
}

fn rc(wave: Waveform, r: f64) -> Chain {
    Chain {
        wave,
        r,
        l: 0.0,
        c: 0.0,
        switch: None,
    }
}

/// Z10: a plate charged through a resistor (V₀ = 2, R = 0.5), facing a grounded one
/// across a gap of 1.5: its potential `V₀(1 − e^{−t/RC})` with C its coefficient of
/// capacitance from the boundary elements (the other plate held at 0); and the field of
/// the driven plates at three points and three times against the static solution with
/// the plate at that potential (linearity).
#[test]
fn z10_plate_charged_through_a_resistor() {
    let (v0, r) = (2.0, 0.5);
    let boxes = vec![
        plate(1.0, Bias::Potential(0.0)),
        plate(-1.0, Bias::Potential(0.0)),
    ];
    let electrodes = Electrodes::with_panel_size(boxes.clone(), &[], PANEL);
    let c = electrodes.capacitance().expect("electrodes")[0][0];
    let tau = r * c;
    let drives = Drives::build(
        &electrodes,
        &[(0, rc(Waveform::Dc(v0), r))],
        &[],
        &[],
        f64::INFINITY,
        10.0 * tau,
    )
    .expect("solves");
    let node = drives.electrodes[0].1;
    let want = |t: f64| v0 * (1.0 - (-t / tau).exp());
    let worst_v = (0..=2000)
        .map(|k| 10.0 * tau * f64::from(k) / 2000.0)
        .map(|t| (drives.solution.potential(node, t) - want(t)).abs() / v0)
        .fold(0.0, f64::max);
    let field = LevelField {
        electrodes: electrodes.clone(),
        drives: Some(Arc::new(drives)),
        ..LevelField::default()
    };
    let mut worst_e: f64 = 0.0;
    for t in [0.3 * tau, tau, 3.0 * tau] {
        let mut b = boxes.clone();
        b[0].bias = Bias::Potential(want(t));
        let reference = LevelField {
            electrodes: Electrodes::with_panel_size(b, &[], PANEL),
            ..LevelField::default()
        };
        for x in [
            DVec3::ZERO,
            DVec3::new(1.5, 0.4, 0.0),
            DVec3::new(6.0, 2.0, 0.0),
        ] {
            let d = (field.sample(x, t).e - reference.sample(x, t).e).length();
            worst_e = worst_e.max(d / (v0 / 1.5));
        }
    }
    println!(
        "Z10 RC plate (C = {c:.6}, τ = {tau:.6}): the potential within {worst_v:.1e} of V₀, \
         the field within {worst_e:.1e} of V₀/d"
    );
    assert!(
        worst_v < 1e-8 && worst_e < 1e-8,
        "Z10 {worst_v:.3e} {worst_e:.3e}"
    );
}

/// Z11: a coil (radius 5, wire radius 0.1, c = 4) driven through R = 5 by a trapezoid
/// pulse (0 → 10 in 1 from t = 0.5, held 4, back in 1), its inductance `L = N/c²` from its
/// self-flux (Z3). (a) Its current against the superposed ramp responses of the RL
/// circuit, `(k/R)[(t − t₀) − τ(1 − e^{−(t − t₀)/τ})]`, τ = L/R. (b) A charge at rest
/// 2 from the axis: its canonical angular momentum `x p_y − y p_x + q(x A_y − y A_x)` about
/// the axis is conserved exactly in any axially symmetric field, so it stays 0 while the
/// induced field `−κ̇ A` accelerates the charge and κB bends it; the four corners are
/// breakpoints, where the flight restarts (its samples end one ulp before each).
#[test]
fn z11_coil_driven_by_a_pulse() {
    let c_light = 4.0;
    let coil = CircularLoop {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        radius: 5.0,
        kappa: 0.0,
        wire_radius: 0.1,
        rate: 0.0,
    };
    let wave = Waveform::Pulse {
        low: 0.0,
        high: 10.0,
        delay: 0.5,
        rise: 1.0,
        width: 4.0,
        fall: 1.0,
        period: 0.0,
    };
    let r = 5.0;
    let t_end = 20.0;
    let drives = Drives::build(
        &Electrodes::default(),
        &[],
        &[coil],
        &[(0, rc(wave, r))],
        c_light,
        t_end,
    )
    .expect("solves");
    let l = physics::inductance::ring_self(&coil) / (c_light * c_light);
    let tau = l / r;
    let ramp = |t: f64, t0: f64, k: f64| {
        if t <= t0 {
            0.0
        } else {
            k / r * ((t - t0) - tau * (1.0 - (-(t - t0) / tau).exp()))
        }
    };
    let want =
        |t: f64| ramp(t, 0.5, 10.0) - ramp(t, 1.5, 10.0) - ramp(t, 5.5, 10.0) + ramp(t, 6.5, 10.0);
    let scale = 10.0 / r;
    let worst_i = (0..=4000)
        .map(|k| t_end * f64::from(k) / 4000.0)
        .map(|t| {
            let (kappa, _) = drives.loop_strength(0, t).expect("driven");
            (kappa * c_light * c_light - want(t)).abs() / scale
        })
        .fold(0.0, f64::max);
    let breaks = drives.breakpoints().to_vec();
    let drives = Arc::new(drives);
    let field = LevelField {
        loops: vec![coil],
        drives: Some(drives.clone()),
        ..LevelField::default()
    };
    let scn = Scenario {
        field,
        obstacles: Vec::new(),
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c: c_light,
        x0: DVec3::new(2.0, 0.0, 0.0),
        p0: DVec3::ZERO,
        detector: None,
        bounds: Some(cube(20.0)),
        t_max: t_end,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
    let (mut worst_l, mut l_scale, mut rho_max): (f64, f64, f64) = (0.0, 0.0, 0.0);
    for s in &tr.samples {
        let (kappa, _) = drives.loop_strength(0, s.t).expect("driven");
        let a = coil.unit_vector_potential(s.x) * kappa;
        let mechanical = s.x.x * s.p.y - s.x.y * s.p.x;
        let canonical = mechanical + (s.x.x * a.y - s.x.y * a.x);
        worst_l = worst_l.max(canonical.abs());
        l_scale = l_scale.max(mechanical.abs());
        rho_max = rho_max.max(s.x.truncate().length());
    }
    let rel = worst_l / l_scale;
    // The flight's samples end one ulp below each breakpoint.
    let restarts = breaks
        .iter()
        .filter(|&&b| tr.samples.iter().any(|s| s.t.to_bits() + 1 == b.to_bits()))
        .count();
    println!(
        "Z11 coil (L = {l:.6}, τ = {tau:.6}): the current within {worst_i:.1e} of V₀/R; a \
         charge from rest: {:?} at t = {}, ρ ≤ {rho_max:.3}, |mechanical L_z| ≤ \
         {l_scale:.4e}, canonical within {rel:.1e} of it; restarts at {restarts} of {} \
         breakpoints ({} steps)",
        tr.outcome,
        tr.end.t,
        breaks.len(),
        tr.stats.n_accept
    );
    assert!(worst_i < 1e-8, "Z11 current {worst_i:.3e}");
    assert!(rel < 1e-10, "Z11 canonical momentum {rel:.3e}");
    assert_eq!(restarts, breaks.len(), "Z11 restarts");
}

/// Z12: a floating plate (charge 0) between a driven plate (DC V₀ = 1 through R = 1) and a
/// grounded one. With C the Maxwell matrix, charge conservation on the floating plate
/// makes `V_F = −(C_FA/C_FF) V_A`, and the driven plate charges through the series
/// capacitance `C_AA − C_AF C_FA/C_FF`: both against the closed forms; and the field at a
/// time against the static solution with the floating plate free (charge 0).
#[test]
fn z12_floating_plate_follows() {
    let (v0, r) = (1.0, 1.0);
    let boxes = vec![
        BoxElectrode {
            center: DVec3::new(0.0, 2.0, 0.0),
            ..plate(0.0, Bias::Potential(0.0))
        },
        plate(0.0, Bias::Charge(0.0)),
        BoxElectrode {
            center: DVec3::new(0.0, -2.0, 0.0),
            ..plate(0.0, Bias::Potential(0.0))
        },
    ];
    let electrodes = Electrodes::with_panel_size(boxes.clone(), &[], PANEL);
    let cm = electrodes.capacitance().expect("electrodes").to_vec();
    let c_series = cm[0][0] - cm[0][1] * cm[1][0] / cm[1][1];
    let tau = r * c_series;
    let drives = Drives::build(
        &electrodes,
        &[(0, rc(Waveform::Dc(v0), r))],
        &[],
        &[],
        f64::INFINITY,
        10.0 * tau,
    )
    .expect("solves");
    let node = |e: usize| {
        drives
            .electrodes
            .iter()
            .find(|&&(k, _)| k == e)
            .expect("in the circuit")
            .1
    };
    let (na, nf) = (node(0), node(1));
    let va = |t: f64| v0 * (1.0 - (-t / tau).exp());
    let mut worst: f64 = 0.0;
    for k in 0..=2000 {
        let t = 10.0 * tau * f64::from(k) / 2000.0;
        let a = (drives.solution.potential(na, t) - va(t)).abs();
        let f = (drives.solution.potential(nf, t) + cm[1][0] / cm[1][1] * va(t)).abs();
        worst = worst.max(a.max(f) / v0);
    }
    let t = 0.7 * tau;
    let field = LevelField {
        electrodes: electrodes.clone(),
        drives: Some(Arc::new(drives)),
        ..LevelField::default()
    };
    let mut b = boxes.clone();
    b[0].bias = Bias::Potential(va(t));
    let reference = LevelField {
        electrodes: Electrodes::with_panel_size(b, &[], PANEL),
        ..LevelField::default()
    };
    let worst_e = [DVec3::new(0.5, 1.0, 0.0), DVec3::new(2.0, -1.0, 0.0)]
        .iter()
        .map(|&x| (field.sample(x, t).e - reference.sample(x, t).e).length() / (v0 / 1.5))
        .fold(0.0, f64::max);
    println!(
        "Z12 floating plate: V_F/V_A = {:.6}, series C = {c_series:.6}; the potentials \
         within {worst:.1e} of V₀, the field within {worst_e:.1e}",
        -cm[1][0] / cm[1][1]
    );
    assert!(
        worst < 1e-8 && worst_e < 1e-8,
        "Z12 {worst:.3e} {worst_e:.3e}"
    );
}

/// Z13: the impulse of a driven plate's field. A heavy charge (m = 1e12, c = ∞, so it
/// barely moves) at rest between the plates of Z10, launched at lab time 0.5 (the field's
/// time offset) and flown for 3: its momentum is `q E₁(x₀) ∫ V dt` over [0.5, 3.5], with
/// E₁ the field per unit potential of the driven plate and `∫ V dt =
/// V₀[(t₂ − t₁) + τ(e^{−t₂/τ} − e^{−t₁/τ})]`.
#[test]
fn z13_impulse_of_a_charging_plate() {
    let (v0, r) = (2.0, 0.5);
    let boxes = vec![
        plate(1.0, Bias::Potential(0.0)),
        plate(-1.0, Bias::Potential(0.0)),
    ];
    let electrodes = Electrodes::with_panel_size(boxes.clone(), &[], PANEL);
    let tau = r * electrodes.capacitance().expect("electrodes")[0][0];
    let (t1, t2) = (0.5, 3.5);
    let drives = Drives::build(
        &electrodes,
        &[(0, rc(Waveform::Dc(v0), r))],
        &[],
        &[],
        f64::INFINITY,
        t2,
    )
    .expect("solves");
    let mut unit = boxes.clone();
    unit[0].bias = Bias::Potential(1.0);
    let e1 = LevelField {
        electrodes: Electrodes::with_panel_size(unit, &[], PANEL),
        ..LevelField::default()
    }
    .sample(DVec3::ZERO, 0.0)
    .e;
    let want = e1 * (v0 * ((t2 - t1) + tau * ((-t2 / tau).exp() - (-t1 / tau).exp())));
    let scn = Scenario {
        field: LevelField {
            electrodes,
            drives: Some(Arc::new(drives)),
            time_offset: t1,
            ..LevelField::default()
        },
        obstacles: Vec::new(),
        particle: Particle {
            charge: 1.0,
            mass: 1e12,
            radius: 0.0,
            moment: 0.0,
        },
        c: f64::INFINITY,
        x0: DVec3::ZERO,
        p0: DVec3::ZERO,
        detector: None,
        bounds: Some(cube(20.0)),
        t_max: t2 - t1,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    let tr = run(&scn, &RunSettings::with_tolerance(1e-12));
    let rel = (tr.end.p - want).length() / want.length();
    println!(
        "Z13 impulse: p = {:.12e} (closed form {:.12e}), {rel:.1e}; {:?} at t = {}",
        tr.end.p.y, want.y, tr.outcome, tr.end.t
    );
    assert!(rel < 1e-8, "Z13 {rel:.3e}");
}
