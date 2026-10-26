//! Validation of the conductors' interiors (`physics::lattice`, PHYSICS.md §2.11, the Drude
//! metals) against the Boltzmann equation of the model (derived in PHYSICS.md: specular
//! hard-disc collisions at the rate `κ v`, `κ = 2 a n_ion`, `⟨cos θ⟩ = −1/3`, the speed
//! redrawn from the collision-weighted thermal law). D1 the conductivity and the collision
//! rate, D2 the AC conductivity, D3 the Joule heat, D4 the conductivity in a magnetic field
//! (the Hall effect), D6 the Einstein relation: carriers in a box periodic both ways,
//! stepped by the lattice's own collision step (`Lattice::advance`, as `Tube` steps them).
//! D5 screening: a metal in `Tube`, its carriers' own field screening two plates' within the
//! Debye length. The self-consistent resistor (Ohm's law with the carriers' own field and
//! the contacts), D7, its thermal noise against Nyquist's, and D8, eddy currents in a
//! changing field, are measurements of minutes, run with `--ignored`.

use physics::DVec3;
use physics::bem::Bias;
use physics::circuit::Waveform;
use physics::lattice::Lattice;
use physics::tube::{Circuit, Emitter, OhmicContact, SolenoidWave, Tube};
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
    /// Under an oscillating acceleration: the drift's parts in phase with it and a
    /// quarter period behind (its projections on cos ωt and sin ωt over the window).
    ac: (f64, f64),
}

/// `carriers` carriers (thermal speed 1) in a box `side × side`, periodic both ways, over
/// random cores (`A`, `MEAN_FREE_PATH`), under the acceleration `acc` along x (times
/// `cos(ac t)` if `ac` > 0) and the rotation `omega` (`q B/m`: the velocity turns at
/// `−omega`), for `t1` with the window from `t0`.
#[allow(clippy::cast_precision_loss, clippy::too_many_arguments)]
fn run_box(carriers: usize, dt: f64, acc: f64, omega: f64, ac: f64, t0: f64, t1: f64) -> BoxRun {
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
    let (mut in_phase, mut quadrature) = (0.0, 0.0);
    let mut msd_mid = 0.0;
    let mut t = 0.0;
    let mut step = 0u64;
    let mut window = 0.0;
    while t < t1 - 0.5 * dt {
        step += 1;
        let inside = t >= t0 - 0.5 * dt;
        // The step's kick (an oscillating drive at the step's middle).
        let (drive_sin, drive_cos) = libm::sincos(ac * (t + 0.5 * dt));
        let kick = if ac > 0.0 { kick * drive_cos } else { kick };
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
                in_phase += dx.x * drive_cos;
                quadrature += dx.x * drive_sin;
            }
            heat += h;
            work += kick.x / dt * dx.x;
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
        ac: (
            2.0 * in_phase / (window * n),
            2.0 * quadrature / (window * n),
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
    let r = run_box(200_000, 0.02, acc, 0.0, 0.0, 5.0, 40.0);
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

/// D2, the AC conductivity: the same carriers under an acceleration `0.2 cos ωt`, ω = 2.5
/// (about the collision rate: `ω τ ≈ 1`), over 14 periods from t = 5. The drift's parts in
/// phase with the drive and a quarter period behind it against the Boltzmann equation's
/// at that frequency: with time going as `e^{−iωt}` its collision term `κ v` becomes
/// `κ v − iω`, so the mobility is `mobility` at `−ω` (Drude's `σ₀/(1 − iωτ)` for a rate
/// independent of the speed). Required 5 % of its magnitude (as D1).
#[test]
fn d2_ac_conductivity() {
    let (acc, w) = (0.2, 2.5);
    let periods = 14.0 * std::f64::consts::TAU / w;
    let r = run_box(200_000, 0.02, acc, 0.0, w, 5.0, 5.0 + periods);
    let (tr, ti) = mobility(r.kappa, -w);
    let (mr, mi) = (r.ac.0 / acc, r.ac.1 / acc);
    let size = libm::sqrt(tr * tr + ti * ti);
    let error = libm::sqrt((mr - tr) * (mr - tr) + (mi - ti) * (mi - ti)) / size;
    println!(
        "D2: drift/a in phase {mr:.5}, behind {mi:.5} against {tr:.5}, {ti:.5}: {:.2} % of its size (DC {:.5})",
        100.0 * error,
        mobility(r.kappa, 0.0).0
    );
    assert!(error < 0.05);
}

/// D3, the Joule heat: the heat the lattice takes at the collisions equals the work the
/// field does on the carriers, less their kinetic energy's change (energy is conserved
/// collision by collision: the path between collisions is the step's parabola). Required
/// `1e-9` of the work (the rounding of ~10⁸ terms). Measured before the bounce was put
/// on the parabola (on the chord): 3.7e-4, the predicted `dt² κ⟨v⟩/(12 μ)`.
#[test]
fn d3_joule_heat_is_the_fields_work() {
    let r = run_box(20_000, 0.02, 0.2, 0.0, 0.0, 5.0, 40.0);
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
    let r = run_box(200_000, 0.02, acc, omega, 0.0, 5.0, 40.0);
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

/// D5, screening: a metal 16 × 4 (background 1, `q/m = −1`, thermal speed 1, the cores of
/// the tests: collisions keep it at the lattice's temperature), no contacts (its carriers
/// conserved), between a plate 1 above it at +0.5 and one 1 below at −0.5. The plates' field
/// is screened within the Debye length `λ_D = σ_v/√(4πρ|q/m|)` = 0.282: the carriers'
/// density departs from the background by `∝ e^{−d/λ_D}` at a distance d from either wall
/// (the linearized Poisson–Boltzmann equation; `eφ/kT` ≈ 0.14 across the layer here). 32 000
/// carriers (40 per Debye area: mean field), over t = 5–30; the decay length fitted to the
/// departure's logarithm (bands weighted by their inverse variance), both walls together
/// (the bottom's sign turned), over d < 0.6. The carriers' softening 0.025 (`h/λ_D` = 0.09,
/// as a level's verification: converged, the same at 0.0125; at 0.05, about the
/// macroparticles' spacing, the screening came out 8–13 % longer). Required 10 % (the
/// finite plasma parameter and the nonlinearity, of order 1 %; the statistics ~4 %).
#[test]
fn d5_screening() {
    let (l, w, v) = (16.0, 4.0, 0.5);
    let plate = |y: f64, bias: Bias| Electrode {
        section: Section::Rect {
            center: DVec3::new(0.5 * l, y, 0.0),
            angle: 0.0,
            half_length: 0.5 * l,
            half_thickness: 0.25,
        },
        bias,
    };
    let electrodes = Electrodes::new(
        vec![
            plate(w + 1.25, Bias::Potential(v)),
            plate(-1.25, Bias::Potential(-v)),
        ],
        &[],
        0.25,
    );
    let (lo, hi) = (DVec3::ZERO, DVec3::new(l, w, 0.0));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(l * w / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 9);
    let carriers = 32_000;
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, 13);
    let tube = Tube {
        arena: Some((
            DVec3::new(-2.0, -3.0, 0.0),
            DVec3::new(l + 2.0, w + 3.0, 0.0),
        )),
        lattices: vec![lattice],
        ..Tube::new(electrodes, 0.025, 0.02)
    };
    let mut s = tube.state();
    let band = 0.05;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let bands = libm::round(w / band) as usize;
    let mut density = vec![0.0; bands];
    let mut samples = 0.0;
    let mut next = 5.0;
    while s.t < 30.0 {
        tube.step(&mut s);
        if s.t < next {
            continue;
        }
        next += 0.2;
        samples += 1.0;
        // The middle half of the length (away from the strip's ends).
        for (x, q) in s.x.iter().zip(&s.charge) {
            if x.x > 0.25 * l && x.x < 0.75 * l && x.y >= 0.0 && x.y < w {
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let b = ((x.y / band) as usize).min(bands - 1);
                density[b] += q.abs();
            }
        }
    }
    let area = 0.5 * l * band;
    let n: Vec<f64> = density.iter().map(|d| d / (samples * area)).collect();
    // The departure from the background near each wall (the top gains carriers, the
    // bottom loses them), and the decay length fitted to its logarithm by least squares,
    // each band weighted by the inverse of its logarithm's variance (∝ departure², the
    // density's noise being about the same in every band; unweighted, the noisy tail set
    // the slope: +14 %).
    let lambda_d = 1.0 / libm::sqrt(4.0 * std::f64::consts::PI);
    let fit_bands = 12;
    let (mut sw, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut profile = Vec::new();
    for k in 0..fit_bands {
        let d = (f64::from(u32::try_from(k).unwrap_or(0)) + 0.5) * band;
        let excess = 0.5 * ((n[bands - 1 - k] - 1.0) - (n[k] - 1.0));
        profile.push(excess);
        let y = libm::log(excess.max(1e-9));
        let wt = excess * excess;
        sw += wt;
        sx += wt * d;
        sy += wt * y;
        sxx += wt * d * d;
        sxy += wt * d * y;
    }
    let slope = (sw * sxy - sx * sy) / (sw * sxx - sx * sx);
    let lambda = -1.0 / slope;
    println!(
        "D5: decay length {lambda:.4} against the Debye length {lambda_d:.4} ({:+.2} %); departure by band: {}",
        100.0 * (lambda / lambda_d - 1.0),
        profile
            .iter()
            .map(|e| format!("{e:.4}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    assert!((lambda / lambda_d - 1.0).abs() < 0.10);
}

/// D6, the Einstein relation: without a field, the carriers spread with
/// `⟨Δx²⟩ = 2 D t` (the slope between the window's middle and end, past the ballistic
/// start), `D = μ kT/m` with D1's `μ`: required 5 % (as D1).
#[test]
fn d6_einstein_relation() {
    let r = run_box(200_000, 0.02, 0.0, 0.0, 0.0, 5.0, 45.0);
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
/// `q/m = −1` and thermal speed 1; `seed` places the cores and the carriers.
fn resistor(l: f64, w: f64, v: f64, carriers: usize, dt: f64, seed: u64) -> Tube {
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
    let cores = Lattice::random_cores(lo, hi, count, A, seed);
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, seed + 4);
    let weight = lattice.carrier_charge().abs();
    let emitter = |electrode: usize, toward: f64| Emitter {
        electrode,
        charge_per_mass: -1.0,
        weight,
        emit_toward: Some(DVec3::new(toward, 0.0, 0.0)),
        saturation: None,
        thermal_speed: 1.0,
        ohmic: Some(OhmicContact {
            density: 1.0,
            depth: 0.5 * MEAN_FREE_PATH,
        }),
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
        emitters: vec![emitter(0, 1.0), emitter(1, -1.0)],
        arena: Some((
            DVec3::new(-2.0, -2.0, 0.0),
            DVec3::new(l + 2.0, w + 2.0, 0.0),
        )),
        circuits: vec![source(0, 0.0), source(1, v)],
        lattices: vec![lattice],
        ..Tube::new(electrodes, 0.05, dt)
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
    let tube = resistor(l, w, 3.2, n, dt, 7);
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

/// The charge the carriers carry along a conductor `l` long, `Σ q Δx / l` from their moves,
/// summed from `t0` and sampled every step until `t1`.
fn carried(tube: &Tube, l: f64, t0: f64, t1: f64) -> Vec<f64> {
    let mut s = tube.state();
    let mut q = 0.0;
    let mut out = Vec::new();
    while s.t < t1 {
        let before: std::collections::HashMap<u32, f64> =
            s.id.iter().zip(&s.x).map(|(id, x)| (*id, x.x)).collect();
        tube.step(&mut s);
        if s.t < t0 {
            continue;
        }
        q +=
            s.id.iter()
                .zip(&s.x)
                .zip(&s.charge)
                .filter_map(|((id, x), c)| before.get(id).map(|x0| c * (x.x - x0)))
                .sum::<f64>()
                / l;
        out.push(q);
    }
    out
}

/// `⟨(Q(t + lag) − Q(t))²⟩` over a series sampled every `dt`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
fn spread(q: &[f64], lag: f64, dt: f64) -> f64 {
    let n = libm::round(lag / dt) as usize;
    let m = q.len() - n;
    (0..m).map(|k| (q[k + n] - q[k]).powi(2)).sum::<f64>() / m as f64
}

/// D7, Johnson–Nyquist (fluctuation–dissipation in the model): the resistor of D1 (16 × 4,
/// 2000 carriers, dt 0.02) at 1.6 V (`qEℓ/kT` = 0.05: the carriers' warming, which lowers
/// σ by about 1.3 % at 0.1, is 0.3 % here), 16 seeds of the cores and carriers, each over
/// t = 8–1608. The charge carried along it, `Q(t)`, gives the resistance `R = V/⟨I⟩` and,
/// about its mean drift, the spread `⟨(ΔQ − ⟨I⟩T)²⟩`, which grows as `2 kT T/R` once the
/// current's correlations have decayed (by T ≈ 0.5 here; the slope between T = 1 and 4,
/// which drops the short-time offset). `kT` is a macroparticle's: `|q_c| σ_v²/|q/m|`
/// (the noise of a line of charge, which a z-invariant world's real carriers would not
/// have: PHYSICS.md). Required 5 % (as D1: a transport coefficient; the statistics, from
/// the seeds' scatter, ~1.7 %).
#[test]
#[ignore = "a measurement (about 20 min in release): run with --ignored (PHYSICS.md §2.11)"]
#[allow(clippy::cast_precision_loss)]
fn d7_johnson_nyquist() {
    let (l, w, n, dt, v) = (16.0, 4.0, 2000, 0.02, 1.6);
    let (t0, t1) = (8.0, 1608.0);
    let started = std::time::Instant::now();
    let runs: Vec<(f64, f64)> = (0..16_u64)
        .into_par_iter()
        .map(|k| {
            let tube = resistor(l, w, v, n, dt, 101 + 17 * k);
            let kt = tube.lattices[0].carrier_charge().abs();
            let q = carried(&tube, l, t0, t1);
            let current = q.last().copied().unwrap_or(0.0) / (t1 - t0);
            let r = v / current.abs();
            // About the mean drift.
            let detrended: Vec<f64> = q
                .iter()
                .enumerate()
                .map(|(i, x)| x - current * dt * (i + 1) as f64)
                .collect();
            let slope = (spread(&detrended, 4.0, dt) - spread(&detrended, 1.0, dt)) / 3.0;
            (r, slope * r / (2.0 * kt))
        })
        .collect();
    let m = runs.len() as f64;
    let mean = runs.iter().map(|r| r.1).sum::<f64>() / m;
    let error =
        libm::sqrt(runs.iter().map(|r| (r.1 - mean).powi(2)).sum::<f64>() / (m * (m - 1.0)));
    let r = runs.iter().map(|r| r.0).sum::<f64>() / m;
    println!(
        "D7: the charge's spread against 2kT/R: {mean:.4} +- {error:.4} over {} seeds (R \
         {r:.3} on average; each: {}); {:.0} s",
        runs.len(),
        runs.iter()
            .map(|r| format!("{:.3}", r.1))
            .collect::<Vec<_>>()
            .join(" "),
        started.elapsed().as_secs_f64()
    );
    assert!((mean - 1.0).abs() < 0.05, "D7 {mean} +- {error}");
}

/// `∫ u dA` of the torsion function of a rectangle `w` × `b` (`∇²u = −1`, `u = 0` on its
/// edge; Saint-Venant): `(w³ b/12) [1 − (192 w/π⁵ b) Σ_{n odd} tanh(nπb/2w)/n⁵]`.
fn torsion_area(w: f64, b: f64) -> f64 {
    let sum: f64 = (0..200)
        .map(|k| f64::from(2 * k + 1))
        .map(|n| libm::tanh(n * PI * b / (2.0 * w)) / libm::pow(n, 5.0))
        .sum();
    w * w * w * b / 12.0 * (1.0 - 192.0 * w / (libm::pow(PI, 5.0) * b) * sum)
}

/// `∫ u dy` along the middle line of the torsion function of a rectangle `w` (across) ×
/// `b` (along the line): `w² b/8 − (8 w³/π⁴) Σ_{n odd} (−1)^((n−1)/2) tanh(nπb/2w)/n⁴`.
fn torsion_midline(w: f64, b: f64) -> f64 {
    let sum: f64 = (0..200)
        .map(|k| (if k % 2 == 0 { 1.0 } else { -1.0 }, f64::from(2 * k + 1)))
        .map(|(sign, n)| sign * libm::tanh(n * PI * b / (2.0 * w)) / libm::pow(n, 4.0))
        .sum();
    w * w * b / 8.0 - 8.0 * w * w * w / libm::pow(PI, 4.0) * sum
}

/// An isolated metal, the rectangle `lo`..`hi` (no contacts), in a field along z that
/// changes in a triangle between `±b0` with period `period`, about `axis`: its eddy currents.
/// Per ramp (but `skip` after each corner), the mean of the power the induced field gives
/// the carriers, `Σ q E·Δx/dt`, and of the current through the rectangle's left half,
/// `Σ q Δy/dt` (signed by the ramp: positive for a rising field): one block per ramp.
/// An isolated metal of the tests' (`carriers` in the rectangle `lo`..`hi`, no contacts) in
/// a field along z that changes in a triangle between `±b0` with period `period`, about
/// `axis` (D8, D9).
fn plate(
    lo: DVec3,
    hi: DVec3,
    axis: DVec3,
    carriers: usize,
    dt: f64,
    b0: f64,
    period: f64,
) -> Tube {
    let (w, b) = (hi.x - lo.x, hi.y - lo.y);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(w * b / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 31);
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, 47);
    // One small electrode away from the metal (the electrostatics needs one), floating with
    // no charge: grounded in this open z-invariant world it could draw charge from the log
    // gauge's "infinity" (r = 1), and a floating metal beside it runs away (OP-1).
    // Far from the metal: beside it, even floating, its induced dipole bent the strip's field.
    let corner = lo - DVec3::new(40.0, 40.0, 0.0);
    let electrodes = Electrodes::new(
        vec![Electrode {
            section: Section::Rect {
                center: corner,
                angle: 0.0,
                half_length: 0.25,
                half_thickness: 0.25,
            },
            bias: Bias::Charge(0.0),
        }],
        &[],
        0.25,
    );
    let wave = SolenoidWave {
        waveform: Waveform::Pulse {
            low: -b0,
            high: b0,
            delay: 0.0,
            rise: 0.5 * period,
            width: 0.0,
            fall: 0.5 * period,
            period,
        },
        axis,
    };
    Tube {
        arena: Some((
            lo - DVec3::new(2.0, 2.0, 0.0),
            hi + DVec3::new(2.0, 2.0, 0.0),
        )),
        b_wave: Some(wave),
        lattices: vec![lattice],
        ..Tube::new(electrodes, 0.05, dt)
    }
}

#[allow(clippy::too_many_arguments)]
fn eddy(
    lo: DVec3,
    hi: DVec3,
    axis: DVec3,
    carriers: usize,
    dt: f64,
    b0: f64,
    period: f64,
    end: f64,
    skip: f64,
) -> Vec<(f64, f64)> {
    let tube = plate(lo, hi, axis, carriers, dt, b0, period);
    let wave = tube.b_wave.expect("a swinging field");
    let mut s = tube.state();
    let (mut power, mut current, mut n) = (0.0, 0.0, 0.0);
    let mut blocks = Vec::new();
    let mut ramp = 0_i64;
    let middle = 0.5 * (lo.x + hi.x);
    while s.t < end {
        let t0 = s.t;
        let before: std::collections::HashMap<u32, DVec3> =
            s.id.iter().zip(&s.x).map(|(id, x)| (*id, *x)).collect();
        tube.step(&mut s);
        // A new ramp: close the last one's block.
        #[allow(clippy::cast_possible_truncation)]
        let this = libm::floor(t0 / (0.5 * period)) as i64;
        if this != ramp {
            if n > 0.0 {
                blocks.push((power / n, current / n));
            }
            (power, current, n, ramp) = (0.0, 0.0, 0.0, this);
        }
        let phase = t0 - period * libm::floor(t0 / period);
        let rising = phase < 0.5 * period;
        let since = if rising { phase } else { phase - 0.5 * period };
        if since < skip || 0.5 * period - since < dt {
            continue;
        }
        let sign = if rising { 1.0 } else { -1.0 };
        let (mut p, mut j) = (0.0, 0.0);
        for ((id, x), q) in s.id.iter().zip(&s.x).zip(&s.charge) {
            let Some(x0) = before.get(id) else {
                continue;
            };
            let d = *x - *x0;
            p += q * wave.induced(t0, *x0).dot(d);
            if x0.x < middle {
                j += q * d.y;
            }
        }
        power += p / dt;
        current += sign * j / dt;
        n += 1.0;
    }
    if n > 0.0 {
        blocks.push((power / n, current / n));
    }
    blocks
}

/// D8, eddy currents: an isolated metal of the tests' (background 1, `q/m = −1`, thermal
/// speed 1, `ℓ = 0.5`, 20 000 carriers) centred on the solenoid's axis, in a field along z
/// ramping between ±0.25 in a triangle of period 20 (`|Ḃ| = 0.05`; `ω_c τ ≤ 0.1`), over
/// t = 0–200 at dt 0.02, 2 time units skipped after each corner (the currents settle within a
/// few `τ`): a square 16 × 16 and a strip 8 × 32. Ohm's law in the induced field
/// `E = −½ Ḃ ẑ × (x − x₀)`: the current's stream function `ψ = −σ Ḃ u`, `u` the rectangle's
/// torsion function (`∇²u = −1`, 0 on its edge), so the power `σ Ḃ² ∫ u dA` and the current
/// through the left half `σ Ḃ ∫ u dy` along the middle line (clockwise for a rising field:
/// up on the left, Lenz's law); σ from the Boltzmann equation (D1). Required 5 % each (as D1:
/// the medium's correction ~2 %, the statistics ~1–2.5 %). At twice the rate (`Ḃ = 0.1`, the
/// field 0.4 at the square's edges, `qEℓ/kT ≈ 0.2`) the carriers warm and the square came
/// out 5 % low: the linear response needs the weaker field. The electrode the electrostatics
/// needs floats: a grounded one in this open world made the strip run away (OP-1).
#[test]
#[ignore = "a measurement (about 6 minutes in release): run with --ignored (PHYSICS.md §2.11)"]
fn d8_eddy_currents() {
    let (n, dt, b0, period) = (20_000, 0.02, 0.25, 20.0);
    // The cores' density sets the collision rate: kappa = 2 a n_ion = 1/l.
    let sigma = mobility(1.0 / MEAN_FREE_PATH, 0.0).0;
    let bdot = 4.0 * b0 / period;
    let mut results = Vec::new();
    for (name, w, b) in [("square 16 x 16", 16.0, 16.0), ("strip 8 x 32", 8.0, 32.0)] {
        let axis = DVec3::new(0.5 * w, 0.5 * b, 0.0);
        let blocks = eddy(
            DVec3::ZERO,
            DVec3::new(w, b, 0.0),
            axis,
            n,
            dt,
            b0,
            period,
            200.0,
            2.0,
        );
        #[allow(clippy::cast_precision_loss)]
        let m = blocks.len() as f64;
        let mean = |f: &dyn Fn(&(f64, f64)) -> f64| blocks.iter().map(f).sum::<f64>() / m;
        let error = |f: &dyn Fn(&(f64, f64)) -> f64, mu: f64| {
            libm::sqrt(blocks.iter().map(|x| (f(x) - mu).powi(2)).sum::<f64>() / ((m - 1.0) * m))
        };
        let (power, current) = (mean(&|x| x.0), mean(&|x| x.1));
        let (sp, sc) = (error(&|x| x.0, power), error(&|x| x.1, current));
        let p_theory = sigma * bdot * bdot * torsion_area(w, b);
        let j_theory = sigma * bdot * torsion_midline(w, b);
        println!(
            "D8, {name}, {m} ramps: power {power:.5} ± {sp:.5} against {p_theory:.5} \
             ({:+.2} ± {:.2} %), current through the left half {current:.5} ± {sc:.5} \
             against {j_theory:.5} ({:+.2} ± {:.2} %)",
            100.0 * (power / p_theory - 1.0),
            100.0 * sp / p_theory,
            100.0 * (current / j_theory - 1.0),
            100.0 * sc / j_theory
        );
        results.push((name, power / p_theory - 1.0, current / j_theory - 1.0));
    }
    for (name, p, c) in results {
        assert!(p.abs() < 0.05, "{name}: power {p}");
        assert!(c.abs() < 0.05, "{name}: current {c}");
    }
}

/// D9, the display's accounts by cell (`physics::tube::CellAccounts`): an isolated plate 6 × 4
/// of the tests' metal (4000 carriers) in D8's swinging field, 400 steps of 0.02: the cells'
/// works add up to the lattice's (`field_work`, `induced_work`), and their charge carried to
/// every carrier's own `Σ q Δx` over the steps it started inside, here booked independently,
/// each to 1e-9 of the sums' scale (the sum of the cells' sizes: rounding over 1.6 million
/// terms is ~1e-10 of it; a carrier missed or booked twice would be ~1e-4).
#[test]
fn d9_cell_accounts_add_up() {
    let (lo, hi) = (DVec3::ZERO, DVec3::new(6.0, 4.0, 0.0));
    let tube = plate(lo, hi, DVec3::new(3.0, 2.0, 0.0), 4000, 0.02, 0.25, 20.0);
    let lattice = &tube.lattices[0];
    let mut s = tube.state();
    let mut flux = [0.0_f64; 2];
    for _ in 0..400 {
        let before: std::collections::HashMap<u32, (DVec3, f64)> =
            s.id.iter()
                .zip(&s.x)
                .zip(&s.charge)
                .filter(|((_, x), _)| lattice.contains(**x))
                .map(|((id, x), q)| (*id, (*x, *q)))
                .collect();
        tube.step(&mut s);
        for (id, x) in s.id.iter().zip(&s.x) {
            if let Some((x0, q)) = before.get(id) {
                flux[0] += q * (x.x - x0.x);
                flux[1] += q * (x.y - x0.y);
            }
        }
    }
    let cells = &s.cells[0];
    let (nx, ny) = lattice.account_size();
    assert_eq!((nx, ny, cells.work.len()), (6, 4, 24));
    let check = |what: &str, cells: &[f64], total: f64| {
        let sum: f64 = cells.iter().sum();
        let scale: f64 = cells.iter().map(|c| c.abs()).sum::<f64>().max(1e-300);
        println!(
            "D9 {what}: cells {sum:.12e}, total {total:.12e}, off by {:.1e} of their scale {scale:.3e}",
            (sum - total).abs() / scale
        );
        assert!((sum - total).abs() < 1e-9 * scale, "{what}: {sum} {total}");
    };
    check("field's work", &cells.work, s.field_work[0]);
    check("induced work", &cells.induced, s.induced_work[0]);
    let fx: Vec<f64> = cells.flux.iter().map(|f| f[0]).collect();
    let fy: Vec<f64> = cells.flux.iter().map(|f| f[1]).collect();
    check("charge carried along x", &fx, flux[0]);
    check("charge carried along y", &fy, flux[1]);
}

/// The carriers' kinetic energy (of their stored velocities).
fn kinetic(s: &physics::tube::TubeState) -> f64 {
    s.v.iter()
        .zip(&s.charge)
        .zip(&s.qm)
        .map(|((v, q), qm)| 0.5 * (q / qm).abs() * v.length_squared())
        .sum()
}

/// What the carriers' energy account leaves over: their kinetic energy's change less the
/// field's work on them, plus the lattice's take and the contacts' net take, since `start`
/// (kinetic energy, work, heat, contacts' take).
fn unbooked(s: &physics::tube::TubeState, start: (f64, f64, f64, f64)) -> f64 {
    let work: f64 = s.field_work.iter().sum();
    let contacts: f64 = s.contact_energy.iter().sum();
    (kinetic(s) - start.0) - ((work - start.1) - (s.joule - start.2) - (contacts - start.3))
}

/// The account's terms now (`unbooked`'s start).
fn account(s: &physics::tube::TubeState) -> (f64, f64, f64, f64) {
    (
        kinetic(s),
        s.field_work.iter().sum(),
        s.joule,
        s.contact_energy.iter().sum(),
    )
}

/// D10, the carriers' energy account (PHYSICS.md §9, Symmetries: energy). The kinetic energy
/// of the carriers' stored velocities changes by the field's work on them (the Boris push
/// keeps its own account exactly: the magnetic rotation keeps |v|, and `½m(v₁² − v₀²) =
/// (q dt/2) E·(v₀ + v₁)`), less the lattice's take in collisions, less what the contacts take
/// in, plus what they give off (`TubeState::contact_energy`). (a) Without cores (no
/// collisions; D1's resistor at 3.2, 4000 carriers, 500 steps of 0.02): exactly, to 1e-10
/// of the largest term. (b) With them: a collision's heat is taken at the velocity of its
/// moment in the step, the stored velocity is the step's middle's (the leapfrog's), and the
/// two differ by the part of the step's kick `g` between them: each collision books
/// `m (½ − φ) g·(u′ − v)` too much or too little (zero mean; `g` mostly the carriers' own
/// grainy field) and `m g² (½ − φ)²` too little: a random walk over the collisions, ∝ dt.
/// Over 10 time units in 40 blocks, the blocks' rms residual at dt = 0.02 against 0.005 must
/// fall by at least 2.5 (∝ dt: 4; the rms of 40 blocks is known to ~11 %), and the run's
/// residual at 0.005 stay below 5e-3 of the field's work (estimated ~1e-3: ~1e5 collisions
/// of ~5e-5 each).
#[test]
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn d10_carriers_energy_account() {
    // (a) No cores.
    let mut tube = resistor(16.0, 4.0, 3.2, 4000, 0.02, 7);
    tube.lattices[0] = Lattice::new(
        DVec3::ZERO,
        DVec3::new(16.0, 4.0, 0.0),
        1.0,
        Vec::new(),
        A,
        -1.0,
        1.0,
        4000,
        11,
    );
    let mut s = tube.state();
    let start = account(&s);
    for _ in 0..500 {
        tube.step(&mut s);
    }
    let end = account(&s);
    let residual = unbooked(&s, start);
    let scale = [end.0, end.1 - start.1, end.3 - start.3]
        .into_iter()
        .fold(0.0_f64, |m, x| m.max(x.abs()));
    println!(
        "D10 (a) without cores: kinetic {:.6} -> {:.6}, the field's work {:.6}, the contacts' \
         net take {:.6}, no collisions ({}): unbooked {:.2e} of {scale:.3}",
        start.0,
        end.0,
        end.1 - start.1,
        end.3 - start.3,
        s.collisions,
        residual.abs() / scale
    );
    assert_eq!(s.collisions, 0);
    assert!(residual.abs() < 1e-10 * scale, "(a) {residual}");
    // (b) With cores, at two steps.
    let blocks = |dt: f64| {
        let tube = resistor(16.0, 4.0, 3.2, 4000, dt, 7);
        let mut s = tube.state();
        let per = (0.25 / dt).round() as usize;
        let first = account(&s);
        let mut out = Vec::new();
        for _ in 0..40 {
            let start = account(&s);
            for _ in 0..per {
                tube.step(&mut s);
            }
            out.push(unbooked(&s, start));
        }
        let total = unbooked(&s, first);
        let work = account(&s).1 - first.1;
        let rms = (out.iter().map(|r| r * r).sum::<f64>() / out.len() as f64).sqrt();
        let mean = out.iter().sum::<f64>() / out.len() as f64;
        (rms, mean, total, work, s.collisions)
    };
    let (coarse, fine) = (blocks(0.02), blocks(0.005));
    for (dt, b) in [(0.02, coarse), (0.005, fine)] {
        println!(
            "D10 (b) dt {dt}: blocks' rms residual {:.3e}, mean {:+.2e}; the run's {:+.4} \
             against the field's work {:.3} ({:+.2e}); {} collisions",
            b.0,
            b.1,
            b.2,
            b.3,
            b.2 / b.3,
            b.4
        );
    }
    println!(
        "D10 (b) the rms falls {:.2} times for a step 4 times shorter",
        coarse.0 / fine.0
    );
    assert!(coarse.0 / fine.0 > 2.5, "(b) {} {}", coarse.0, fine.0);
    assert!((fine.2 / fine.3).abs() < 5e-3, "(b) {} {}", fine.2, fine.3);
}

/// The circuit's energy through D1's resistor (`Σ V (given off − taken in)` over its contacts)
/// against the field's work on the carriers `W` (D10) and the contacts' potential steps `S`
/// (`TubeState::potential_steps`: a carrier is given off and taken in where the potential
/// differs from its contact's, in the layer its contact keeps): `circuit = W + S − H`, the
/// remainder `H` the step's own (the leapfrog's work near the contacts; the exact motion
/// conserves the energy). Over t = 5–25 in 8 blocks: at V = 0, no current, `H` must fall at
/// least 3 times for a step 4 times shorter (first order: 4); at V = 3.2, with the current,
/// it must be V = 0's within 3 standard errors (the blocks' scatter): it does not depend on
/// the current, so the steps carry the circuit's share exactly.
#[test]
#[allow(clippy::cast_precision_loss)]
fn d11_circuit_energy_account() {
    let run = |v: f64, dt: f64| {
        let tube = resistor(16.0, 4.0, v, 4000, dt, 7);
        let mut s = tube.state();
        s.book_potential_steps = true;
        let terms = |s: &physics::tube::TubeState| {
            let circuit: f64 = (0..2)
                .map(|e| {
                    let given: f64 = tube
                        .emitters
                        .iter()
                        .zip(&s.released)
                        .filter(|(em, _)| em.electrode == e)
                        .map(|(_, q)| q)
                        .sum();
                    s.potentials[e] * (given - s.collected[e].iter().sum::<f64>())
                })
                .sum();
            let work: f64 = s.field_work.iter().sum();
            let steps: f64 = s.potential_steps.iter().sum();
            (circuit, work, steps)
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let per = libm::round(2.5 / dt) as usize;
        for _ in 0..2 * per {
            tube.step(&mut s);
        }
        let first = terms(&s);
        let mut blocks = Vec::new();
        for _ in 0..8 {
            let a = terms(&s);
            for _ in 0..per {
                tube.step(&mut s);
            }
            let b = terms(&s);
            blocks.push(((b.1 - a.1) + (b.2 - a.2) - (b.0 - a.0)) / 2.5);
        }
        let end = terms(&s);
        let rest = blocks.iter().sum::<f64>() / 8.0;
        let spread = (blocks.iter().map(|r| (r - rest).powi(2)).sum::<f64>() / 7.0).sqrt();
        let rate = |k: fn(&(f64, f64, f64)) -> f64| (k(&end) - k(&first)) / 20.0;
        (
            rest,
            spread / 8.0_f64.sqrt(),
            rate(|t| t.0),
            rate(|t| t.1),
            rate(|t| t.2),
        )
    };
    let cases: Vec<(f64, f64)> = vec![(0.0, 0.02), (0.0, 0.005), (3.2, 0.02)];
    let out: Vec<_> = cases.par_iter().map(|&(v, dt)| run(v, dt)).collect();
    for ((v, dt), (rest, err, circuit, work, steps)) in cases.iter().zip(&out) {
        println!(
            "D11 V {v} dt {dt}: per time the circuit gives {circuit:.4}, the field works {work:.4}, \
             the steps {steps:+.4}: the remainder {rest:+.4} ± {err:.4}"
        );
    }
    let (coarse, fine, current) = (out[0], out[1], out[2]);
    println!(
        "D11 the remainder falls {:.2} times for a step 4 times shorter; with the current it is \
         {:+.4} off V = 0's ({:.1} standard errors)",
        coarse.0 / fine.0,
        current.0 - coarse.0,
        (current.0 - coarse.0).abs() / libm::hypot(coarse.1, current.1)
    );
    assert!(coarse.0 / fine.0 > 3.0, "{coarse:?} {fine:?}");
    assert!(
        (current.0 - coarse.0).abs() < 3.0 * libm::hypot(coarse.1, current.1),
        "{current:?} {coarse:?}"
    );
}
