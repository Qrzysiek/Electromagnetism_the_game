//! Validation of the conductors' interiors (`physics::lattice`, PHYSICS.md §2.11, the Drude
//! metals) against the Boltzmann equation of the model (derived in PHYSICS.md: specular
//! hard-disc collisions at the rate `κ v`, `κ = 2 a n_ion`, `⟨cos θ⟩ = −1/3`, the speed
//! redrawn from the collision-weighted thermal law). D1 the conductivity and the collision
//! rate, D3 the Joule heat, D4 the conductivity in a magnetic field (the Hall effect), D6
//! the Einstein relation: carriers in a box periodic both ways, stepped by the lattice's
//! own collision step (`Lattice::advance`, as `Tube` steps them). The self-consistent
//! resistor (Ohm's law with the carriers' own field and the contacts) is a measurement of
//! a few minutes, run with `--ignored`.

use physics::DVec3;
use physics::bem::Bias;
use physics::circuit::Waveform;
use physics::lattice::Lattice;
use physics::tube::{Circuit, Emitter, Tube};
use physics::zinv::{Electrode, Electrodes, Section};
use rayon::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};

/// The cores' radius and the mean free path of the tests: `a² n_ion = a/(2ℓ) = 0.00125`
/// (dilute: the medium's own correction, `≈ 1.1 a² n ln(a² n)` measured, is −0.7 %
/// here and −2.5 % at 0.005; PHYSICS.md).
const A: f64 = 0.00125;
const MEAN_FREE_PATH: f64 = 0.5;
/// `⟨cos θ⟩` of a specular bounce off a disc at a uniform impact parameter.
const COS_THETA: f64 = -1.0 / 3.0;

/// What a box run measured (per carrier, thermal speed 1).
struct BoxRun {
    /// `κ = 2 a n_ion` of the cores placed.
    kappa: f64,
    /// The mean velocity over the window, from the displacements.
    drift: DVec3,
    /// Collisions per carrier per time over the window.
    rate: f64,
    /// Over the whole run, per mass: the heat the lattice took, the field's work, the
    /// kinetic energy's change.
    heat: f64,
    work: f64,
    kinetic: f64,
    /// The mean square displacement along x at the window's middle and end.
    msd: (f64, f64),
}

/// `carriers` carriers (thermal speed 1) in a box `side × side`, periodic both ways, over
/// random cores (`A`, `MEAN_FREE_PATH`), under the acceleration `acc` along x and the
/// rotation `omega` (`q B/m`: the velocity turns at `−omega`), for `t1` with the
/// window from `t0`.
#[allow(clippy::cast_precision_loss)]
fn run_box(carriers: usize, dt: f64, acc: f64, omega: f64, t0: f64, t1: f64) -> BoxRun {
    let side = 8.0;
    let (lo, hi) = (DVec3::ZERO, DVec3::new(side, side, 0.0));
    let n_ion = 1.0 / (2.0 * A * MEAN_FREE_PATH);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(side * side * n_ion) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 5);
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, 11);
    let kappa = 2.0 * A * lattice.cores.len() as f64 / (side * side);
    // Position, stored (half step) velocity, displacement over the window.
    let kick = DVec3::new(acc * dt, 0.0, 0.0);
    let mut state: Vec<(DVec3, DVec3, DVec3)> = lattice
        .fill()
        .into_iter()
        .map(|(x, v)| (x, v - kick * 0.5, DVec3::ZERO))
        .collect();
    let kinetic0: f64 = state
        .iter()
        .map(|s| 0.5 * (s.1 + kick * 0.5).length_squared())
        .sum();
    let (sn, cs) = libm::sincos(-omega * dt);
    let (mut flux, mut hits, mut heat, mut work) = (DVec3::ZERO, 0u64, 0.0, 0.0);
    let mut msd_mid = 0.0;
    let mut t = 0.0;
    let mut step = 0u64;
    let mut window = 0.0;
    while t < t1 - 0.5 * dt {
        step += 1;
        let inside = t >= t0 - 0.5 * dt;
        let moves: Vec<(DVec3, f64, u32)> = state
            .par_iter_mut()
            .enumerate()
            .map(|(id, (x, v, d))| {
                // Boris: half the kick, the rotation, the other half.
                let v0 = *v;
                let m = v0 + kick * 0.5;
                let r = DVec3::new(m.x * cs - m.y * sn, m.x * sn + m.y * cs, 0.0);
                let v1 = r + kick * 0.5;
                let a = *x;
                let id = u32::try_from(id).unwrap_or(u32::MAX);
                let (b, v_end, h, k) = match lattice.advance(a, v0, v1, dt, id, step) {
                    Some(m) => (m.x, m.v, m.heat_per_mass, m.collisions),
                    None => (a + v1 * dt, v1, 0.0, 0),
                };
                let dx = b - a;
                if inside {
                    *d += dx;
                }
                *x = DVec3::new(b.x.rem_euclid(side), b.y.rem_euclid(side), 0.0);
                *v = v_end;
                (dx, h, k)
            })
            .collect();
        for (dx, h, k) in &moves {
            if inside {
                flux += *dx;
                hits += u64::from(*k);
            }
            heat += h;
            work += acc * dx.x;
        }
        if inside {
            window += dt;
        }
        t += dt;
        if msd_mid == 0.0 && window >= 0.5 * (t1 - t0) {
            msd_mid = state.iter().map(|s| s.2.x * s.2.x).sum::<f64>() / carriers as f64;
        }
    }
    let n = carriers as f64;
    let kinetic1: f64 = state
        .iter()
        .map(|s| 0.5 * (s.1 + kick * 0.5).length_squared())
        .sum();
    BoxRun {
        kappa,
        drift: flux / (window * n),
        rate: hits as f64 / (window * n),
        heat: heat / n,
        work: work / n,
        kinetic: (kinetic1 - kinetic0) / n,
        msd: (
            msd_mid,
            state.iter().map(|s| s.2.x * s.2.x).sum::<f64>() / n,
        ),
    }
}

/// The Boltzmann equation's mean velocity per unit acceleration (along x) of carriers of
/// thermal speed 1 colliding at `κ v`, rotating at `−ω` (PHYSICS.md §2.11):
/// `(v_x + i v_y)/a = (π/κ) P/(1 − ⟨cos θ⟩ Q)`, `P = κ ∫ v³ f₀/(κ v + i ω) dv` with
/// `f₀ = e^{−v²/2}/2π`, `Q = κ ∫ v g/(κ v + i ω) dv` with the drawn speeds' law
/// `g = √(2/π) v² e^{−v²/2}`. At `ω = 0`: `(3/8) √(π/2)/κ`. Simpson's rule on [0, 14].
fn mobility(kappa: f64, omega: f64) -> (f64, f64) {
    let (n, top) = (28_000, 14.0);
    let h = top / f64::from(n);
    let (mut p, mut q) = ((0.0, 0.0), (0.0, 0.0));
    for k in 0..=n {
        let v = h * f64::from(k);
        let w = if k == 0 || k == n {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        // κ/(κv + iω) = κ (κv − iω)/(κ²v² + ω²).
        let d = kappa * kappa * v * v + omega * omega;
        if d == 0.0 {
            continue;
        }
        let (re, im) = (kappa * kappa * v / d, -kappa * omega / d);
        let e = libm::exp(-0.5 * v * v);
        let fp = v * v * v * e / (2.0 * PI);
        let fq = v * libm::sqrt(2.0 / PI) * v * v * e;
        p = (p.0 + w * fp * re, p.1 + w * fp * im);
        q = (q.0 + w * fq * re, q.1 + w * fq * im);
    }
    let (p, q) = (
        (p.0 * h / 3.0, p.1 * h / 3.0),
        (q.0 * h / 3.0, q.1 * h / 3.0),
    );
    // (π/κ) p / (1 − c q).
    let den = (1.0 - COS_THETA * q.0, -COS_THETA * q.1);
    let dd = den.0 * den.0 + den.1 * den.1;
    let ratio = (
        (p.0 * den.0 + p.1 * den.1) / dd,
        (p.1 * den.0 - p.0 * den.1) / dd,
    );
    (PI / kappa * ratio.0, PI / kappa * ratio.1)
}

/// D1, the conductivity: 100 000 carriers over cores with `a² n_ion = 0.005` and
/// `ℓ = 0.5` (κ = 2), an acceleration 0.2 (the field gives `qEℓ/kT = 0.1` over a free
/// path: the linear response), over t = 5–40 at dt 0.02. The drift (from the
/// displacements) against the Boltzmann equation's `(3/8)√(π/2) a/κ`, that is
/// `σ = (3/16) √(π/2) ρ_c |q/m| / (a n_ion σ_v)`: required 5 % (the medium's corrections
/// at this density, the field's nonlinearity, each ~2 %; the statistics, 0.8 %). The
/// collision rate against `κ ⟨v⟩`, `⟨v⟩ = √(π/2)`: required 1 % (the free paths of
/// randomly placed discs are exactly exponential at this rate; statistics 0.1 %).
#[test]
fn d1_conductivity_and_collision_rate() {
    let acc = 0.2;
    let r = run_box(200_000, 0.02, acc, 0.0, 5.0, 40.0);
    let theory = mobility(r.kappa, 0.0).0;
    let closed = 3.0 / 8.0 * libm::sqrt(FRAC_PI_2) / r.kappa;
    let rate = r.kappa * libm::sqrt(FRAC_PI_2);
    println!(
        "D1: mobility {:.5} against {theory:.5} (closed form {closed:.5}): {:+.2} %; collision rate {:.4} against {rate:.4}: {:+.2} %",
        r.drift.x / acc,
        100.0 * (r.drift.x / acc / theory - 1.0),
        r.rate,
        100.0 * (r.rate / rate - 1.0)
    );
    assert!(
        (theory / closed - 1.0).abs() < 1e-6,
        "quadrature {theory} vs {closed}"
    );
    assert!((r.drift.x / acc / theory - 1.0).abs() < 0.05);
    assert!((r.rate / rate - 1.0).abs() < 0.01);
}

/// D3, the Joule heat: the heat the lattice takes at the collisions equals the work the
/// field does on the carriers, less their kinetic energy's change (energy is conserved
/// collision by collision: the path between collisions is the step's parabola). Required
/// `1e-9` of the work (the rounding of ~10⁸ terms). Measured before the bounce was put
/// on the parabola (on the chord): 3.7e-4, the predicted `dt² κ⟨v⟩/(12 μ)`.
#[test]
fn d3_joule_heat_is_the_fields_work() {
    let r = run_box(20_000, 0.02, 0.2, 0.0, 5.0, 40.0);
    let balance = r.work - r.heat - r.kinetic;
    println!(
        "D3: work {:.6}, heat {:.6}, kinetic change {:.6}: balance {:.3e} of the work",
        r.work,
        r.heat,
        r.kinetic,
        balance / r.work
    );
    assert!((balance / r.work).abs() < 1e-9);
}

/// D4, the Hall effect: the same carriers rotating at `ω = 1.25` (`ω τ ≈ 0.5`,
/// `τ = 1/(κ ⟨v⟩)`): both components of the drift against the Boltzmann equation at this
/// field (`mobility`), required 5 % of its magnitude (as D1). In a strip the Hall field
/// cancels the transverse current: `R_H = Im(1/σ̃)/B`, `(4/π)/(n q)` in a weak field (the
/// Hall factor 4/π of a collision rate ∝ v), checked on the quadrature to 1e-4.
#[test]
fn d4_conductivity_in_a_magnetic_field() {
    let (acc, omega) = (0.2, 1.25);
    let r = run_box(200_000, 0.02, acc, omega, 5.0, 40.0);
    let (tx, ty) = mobility(r.kappa, omega);
    let (mx, my) = (r.drift.x / acc, r.drift.y / acc);
    let size = libm::sqrt(tx * tx + ty * ty);
    let error = libm::sqrt((mx - tx) * (mx - tx) + (my - ty) * (my - ty)) / size;
    // The Hall factor in a weak field: R_H n q = −Im(1/μ̃)/ω with μ̃ the mobility above
    // (σ̃ = n q² μ̃/m).
    let w = 1e-3;
    let (wx, wy) = mobility(r.kappa, w);
    let hall = -(wy / (wx * wx + wy * wy)) / w;
    println!(
        "D4: drift/a ({mx:.5}, {my:.5}) against ({tx:.5}, {ty:.5}): {:.2} % of its size; weak-field Hall factor {hall:.5} (4/π = {:.5})",
        100.0 * error,
        4.0 / PI
    );
    assert!(error < 0.05);
    assert!((hall / (4.0 / PI) - 1.0).abs() < 1e-4, "Hall factor {hall}");
}

/// D6, the Einstein relation: without a field, the carriers spread with
/// `⟨Δx²⟩ = 2 D t` (the slope between the window's middle and end, past the ballistic
/// start), `D = μ kT/m` with D1's `μ`: required 5 % (as D1).
#[test]
fn d6_einstein_relation() {
    let r = run_box(200_000, 0.02, 0.0, 0.0, 5.0, 45.0);
    let d = (r.msd.1 - r.msd.0) / (2.0 * 20.0);
    let theory = mobility(r.kappa, 0.0).0;
    println!(
        "D6: diffusion {d:.5} against mu kT/m {theory:.5}: {:+.2} %",
        100.0 * (d / theory - 1.0)
    );
    assert!((d / theory - 1.0).abs() < 0.05);
}

/// A conductor `l × w` between two contacts (boxes one cell thick at its ends, the right
/// at `v`, the left grounded, both emitting carriers into it where the field pulls them
/// out), random cores (`A`, `MEAN_FREE_PATH`), the background 1 per area, carriers of
/// `q/m = −1` and thermal speed 1.
fn resistor(l: f64, w: f64, v: f64, carriers: usize, dt: f64) -> Tube {
    let contact = |x: f64, bias: Bias| Electrode {
        section: Section::Rect {
            center: DVec3::new(x, 0.5 * w, 0.0),
            angle: 0.0,
            half_length: 0.5,
            half_thickness: 0.5 * w,
        },
        bias,
    };
    let electrodes = Electrodes::new(
        vec![
            contact(-0.5, Bias::Grounded),
            contact(l + 0.5, Bias::Potential(v)),
        ],
        &[],
        0.25,
    );
    let (lo, hi) = (DVec3::ZERO, DVec3::new(l, w, 0.0));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(l * w / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 7);
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, 11);
    let weight = lattice.carrier_charge().abs();
    let emitter = |electrode: usize, toward: f64| Emitter {
        electrode,
        charge_per_mass: -1.0,
        weight,
        emit_toward: Some(DVec3::new(toward, 0.0, 0.0)),
        saturation: None,
    };
    let source = |electrode: usize, v: f64| Circuit {
        electrode,
        source: Waveform::Dc(v),
        resistance: 0.0,
        capacitance: 0.0,
        inductance: 0.0,
        to: None,
    };
    Tube {
        electrodes,
        emitters: vec![emitter(0, 1.0), emitter(1, -1.0)],
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: Some((
            DVec3::new(-2.0, -2.0, 0.0),
            DVec3::new(l + 2.0, w + 2.0, 0.0),
        )),
        b_z: 0.0,
        circuits: vec![source(0, 0.0), source(1, v)],
        lattices: vec![lattice],
    }
}

/// D1 in a resistor, Ohm's law with the carriers' own field: a conductor 16 × 4 between
/// two contacts at 3.2 V (a bulk field ~0.2, as D1), 16 000 carriers, dt 0.02, over
/// t = 8–136 in blocks of 4. The current density from the carriers' moves (the charge
/// carried; a carrier's stored velocity after a collision is not its mean over the step),
/// the field averaged over the middle half's area; `σ = j/E` against the Boltzmann
/// equation's, required 5 % (as D1; the statistics, from the blocks' scatter, ~1 %).
#[test]
#[ignore = "a measurement (about 100 s in release): run with --ignored (PHYSICS.md §2.11)"]
fn d1_ohms_law_in_a_resistor() {
    let (l, w, n, dt) = (16.0, 4.0, 16_000, 0.02);
    let tube = resistor(l, w, 3.2, n, dt);
    #[allow(clippy::cast_precision_loss)]
    let kappa = 2.0 * A * tube.lattices[0].cores.len() as f64 / (l * w);
    let theory = mobility(kappa, 0.0).0;
    let mut s = tube.state();
    let (x0, x1) = (0.25 * l, 0.75 * l);
    let (mut bj, mut be, mut bn, mut ben) = (0.0, 0.0, 0.0, 0.0);
    let mut blocks: Vec<(f64, f64)> = Vec::new();
    let mut step = 0u64;
    while s.t < 136.0 {
        let before: std::collections::HashMap<u32, f64> =
            s.id.iter().zip(&s.x).map(|(id, x)| (*id, x.x)).collect();
        tube.step(&mut s);
        step += 1;
        if s.t < 8.0 {
            continue;
        }
        bj +=
            s.id.iter()
                .zip(&s.x)
                .zip(&s.charge)
                .filter_map(|((id, x), q)| before.get(id).map(|x0| q * (x.x - x0)))
                .sum::<f64>()
                / (dt * l * w);
        bn += 1.0;
        if step.is_multiple_of(10) {
            be += (0..33)
                .flat_map(|i| (0..9).map(move |k| (i, k)))
                .map(|(i, k)| {
                    let x = x0 + (x1 - x0) * f64::from(i) / 32.0;
                    let y = w * (f64::from(k) + 0.5) / 9.0;
                    tube.field_at(&s, DVec3::new(x, y, 0.0)).x
                })
                .sum::<f64>()
                / (33.0 * 9.0);
            ben += 1.0;
        }
        #[allow(clippy::cast_precision_loss)]
        if s.t - 8.0 >= 4.0 * (blocks.len() + 1) as f64 {
            blocks.push((bj / bn, be / ben));
            (bj, be, bn, ben) = (0.0, 0.0, 0.0, 0.0);
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let m = blocks.len() as f64;
    let j = blocks.iter().map(|b| b.0).sum::<f64>() / m;
    let e = blocks.iter().map(|b| b.1).sum::<f64>() / m;
    let scatter = |f: &dyn Fn(&(f64, f64)) -> f64, mean: f64| {
        libm::sqrt(blocks.iter().map(|b| (f(b) - mean).powi(2)).sum::<f64>() / ((m - 1.0) * m))
            / mean.abs()
    };
    let (sj, se) = (scatter(&|b| b.0, j), scatter(&|b| b.1, e));
    // ρ_c |q/m| = 1: σ = j/E against the mobility.
    let sigma = j / e;
    println!(
        "D1 resistor: j {j:.5e}, E {e:.5e}: sigma {sigma:.5} against {theory:.5}: {:+.2} % (statistics {:.2} %)",
        100.0 * (sigma / theory - 1.0),
        100.0 * libm::sqrt(sj * sj + se * se)
    );
    assert!((sigma / theory - 1.0).abs() < 0.05);
}
