//! Validation of the charges' own magnetic field in the tube engine (`physics::darwin2d`,
//! `Tube::darwin`; PHYSICS.md §2.12, docs/DARWIN.md). V-D5: a metal ring (a one-turn
//! solenoid along z) in a ramping external field: its current rises at `L/R` to the final
//! current it reaches at once without the self-field, `L` from its magnetic energy and `R`
//! from its Joule heat, and its own field in the hole is `μ₀` times its current (Ampère: a
//! solenoid's field is uniform inside and none outside).

use physics::DVec3;
use physics::circuit::Waveform;
use physics::darwin2d::Darwin2d;
use physics::lattice::{Lattice, Material};
use physics::tube::{SolenoidWave, Tube, TubeState};
use std::f64::consts::PI;

mod common;
use common::{Complex64, ci0, csqrt, measured, mobility, significant};

/// The tests' metal (as `drude.rs`): cores of radius `A` for a mean free path of 0.5.
const A: f64 = 0.00125;
const MEAN_FREE_PATH: f64 = 0.5;

/// A run's sample: (t, the ring's current, the magnetic energy, the external field's work).
type Sample = (f64, f64, f64, f64, f64);

/// The ring: 12 × 12 with a square hole leaving a wall `wall` wide, background charge
/// density 1, carriers of q/m −1 and thermal speed 1, in a field `b0 sin ωt` about its
/// centre; the self-field at the speed of light `c` (None: none).
fn ring(
    c: Option<f64>,
    carriers: usize,
    wall: f64,
    b0: f64,
    omega: f64,
    dt: f64,
    seed: u64,
) -> Tube {
    let (lo, hi) = (DVec3::ZERO, DVec3::new(12.0, 12.0, 0.0));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(144.0 / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, seed ^ 31);
    let lattice = Lattice::new(lo, hi, common::metal(A), cores, carriers, seed).with_hole(
        DVec3::new(wall, wall, 0.0),
        DVec3::new(12.0 - wall, 12.0 - wall, 0.0),
    );
    // One small floating electrode far away (the electrostatics needs one; OP-1).
    let electrodes = common::small_floating_electrode(lo - DVec3::new(40.0, 40.0, 0.0));
    Tube {
        arena: Some((
            lo - DVec3::new(2.0, 2.0, 0.0),
            hi + DVec3::new(2.0, 2.0, 0.0),
        )),
        b_wave: Some(SolenoidWave {
            waveform: Waveform::Sine {
                offset: 0.0,
                amplitude: b0,
                omega,
                phase: 0.0,
            },
            axis: DVec3::new(6.0, 6.0, 0.0),
        }),
        lattices: vec![lattice],
        darwin: c.map(|c| Darwin2d::new(127, 127, 0.25, (-10.0, -10.0), c)),
        ..Tube::new(electrodes, 0.05, dt)
    }
}

/// One run to `end`: samples every `every` of (t, the ring's current counter-clockwise,
/// the magnetic energy, the external field's work on the carriers so far). The current is
/// read from the self-field in the hole, `B/μ₀` (Ampère: a solenoid's field inside is `μ₀`
/// times its current per length): every carrier round the ring in one number. (Counted from
/// one arm's carriers, ~1500 drifting at 0.05 against a thermal speed of 1, it scattered by
/// ~25 % a sample, and their velocities after a step lie `dt/τ`, ~4 %, below the drift the
/// field is solved from, the collisions' share.)
fn run(tube: &Tube, end: f64, every: f64) -> Vec<Sample> {
    let mut s = tube.state();
    let mut out = Vec::new();
    let mut next = 0.0;
    while s.t < end {
        tube.step(&mut s);
        if s.t >= next {
            let d = s.darwin.as_ref().expect("the field");
            let b = d.b_z();
            out.push((
                s.t,
                d.b_at(&b, 6.0, 6.0) / d.mu0,
                d.magnetic_energy(),
                s.induced_work[0],
                s.joule,
            ));
            next += every;
        }
    }
    out
}

/// V-D5, the ring's impedance. A thin ring (its wall 1, two mean free paths, against a skin
/// depth `1/√(ωμ₀σ)` ≈ 2: a lumped inductor, the current nearly uniform across the wall) in a
/// field `0.125 sin(0.8 t)` (an EMF of amplitude ~12 around it): its current, read from the
/// field in the hole, lags the EMF by `φ`,
/// `tan φ = ωL_tot/R`, `L_tot` the carriers' inertia (a kinetic inductance, the run without a
/// self-field: c = 10⁴) plus, with the self-field (c = 3), the geometric inductance L. So
/// `(tan φ_self − tan φ_free)/ω = L/R`, against L from the field's energy (`⟨W⟩` less the
/// undriven twin's thermal `⟨W⟩`, `= L ⟨I²⟩/2`) and R from the external field's work on the
/// carriers (`P = R ⟨I²⟩`) in the same run. Lock-in at ω over 8 periods after one, two seeds.
/// (On a ring with a wall of 3, about the skin depth at ω = 0.4, R came out 13–17 % higher
/// with the self-field than without, whichever the kick's order: the current crowds across
/// the wall, and R by the external field's work, which grows as r across it, follows; a
/// lumped R does not exist there.)
/// And the lattice's heat against that work, the Joule heat with the self-field: its ratio to
/// the work within 5 % of the one without (the coupling's target; set 2026-10-27 before
/// measuring it; the free run holds the plain model's own gap, if any).
/// Required (set before measuring from the noise: each phase ~0.03 rad, L's energy ~5 %):
/// L/R within 15 %; the currents' magnitudes in the ratio their phases give,
/// `|I_free/I_self| = cos φ_free/cos φ_self` (the same EMF and R), within 5 %; R the same with
/// and without the self-field within 15 % (first 5 %, set without estimating R's noise: the
/// thermal current's work against the drive and the driven amplitude's noise put each run's R
/// at ±6–9 % over 8 periods, the two seeds' ratio at ±7 %). (Until 2026-10-26 a ramp's rise time `∫(I∞ − I)dt/I∞`: its
/// thermal noise, ±1 on τ ≈ 2 a seed, left nine seeds undecided.)
#[test]
#[ignore = "four minutes: six runs of 2444 carriers over 9 periods"]
fn ring_impedance_is_r_plus_i_omega_l() {
    let (carriers, wall, dt, b0, omega) = (2444, 1.0, 0.02, 0.125, 0.8);
    let c = 3.0;
    let period = 2.0 * PI / omega;
    let (start, end) = (period, 9.0 * period);
    // The current's phasor `C_I = (2/T)∫ I e^{−iωt} dt` against the EMF's, `∝ −Ḃ` (Lenz: a
    // resistor carries `I ∝ −cos ωt`): `z = C_I/(−b0 ω)`, `arg z = −arg Z`; the mean power and
    // energy over the window; and the driven current's `⟨I²⟩ = |C_I|²/2` (the raw mean square
    // also holds the ring's thermal current, `kT/L_tot` by equipartition: 40 % of the driven
    // one without the self-field, whose `L_tot` is only the carriers' inertia, 7 % with it).
    let lock_in = |v: &[Sample]| {
        let w: Vec<&Sample> = v.iter().filter(|r| r.0 >= start).collect();
        #[allow(clippy::cast_precision_loss)]
        let n = w.len() as f64;
        let (mut re, mut im, mut energy) = (0.0, 0.0, 0.0);
        for r in &w {
            re += r.1 * libm::cos(omega * r.0);
            im -= r.1 * libm::sin(omega * r.0);
            energy += r.2;
        }
        let phasor = Complex64::new(2.0 * re / n, 2.0 * im / n);
        let z = phasor / Complex64::new(-b0 * omega, 0.0);
        let (first, last) = (w[0], w[w.len() - 1]);
        let power = (last.3 - first.3) / (last.0 - first.0);
        let heat = (last.4 - first.4) / (last.0 - first.0);
        let square = 0.5 * (phasor.re * phasor.re + phasor.im * phasor.im);
        (z, power, energy / n, square, heat)
    };
    let mut rows = Vec::new();
    for seed in [47_u64, 53] {
        let free = run(
            &ring(Some(1e4), carriers, wall, b0, omega, dt, seed),
            end,
            0.1,
        );
        let selfish = run(
            &ring(Some(c), carriers, wall, b0, omega, dt, seed),
            end,
            0.1,
        );
        let quiet = run(
            &ring(Some(c), carriers, wall, 0.0, omega, dt, seed),
            end,
            0.1,
        );
        let (zf, pf, _, i2f, hf) = lock_in(&free);
        let (zs, ps, ws, i2s, hs) = lock_in(&selfish);
        let (_, _, wq, _, _) = lock_in(&quiet);
        let (tf, ts) = (-zf.im / zf.re, -zs.im / zs.re);
        let (rf, rs) = (pf / i2f, ps / i2s);
        let l = 2.0 * (ws - wq) / i2s;
        let measured = (ts - tf) / omega;
        let magnitudes = libm::hypot(zf.re, zf.im) / libm::hypot(zs.re, zs.im);
        let cosines = libm::sqrt((1.0 + ts * ts) / (1.0 + tf * tf));
        println!(
            "V-D5 seed {seed}: tan φ free {tf:.4}, self {ts:.4}: L/R {measured:.4}; L {l:.2} (μ₀ × {:.1}), R {rs:.3} (free {rf:.3}): L/R {:.4}; |I_free/I_self| {magnitudes:.4} vs {cosines:.4}",
            l * c * c / (4.0 * PI),
            l / rs
        );
        println!(
            "V-D5 seed {seed}: heat/work {:.4} with the self-field ({hs:.3} of {ps:.3}), {:.4} without ({hf:.3} of {pf:.3})",
            hs / ps,
            hf / pf
        );
        rows.push((
            measured,
            l / rs,
            rs / rf,
            magnitudes / cosines,
            (hs / ps) / (hf / pf),
        ));
    }
    #[allow(clippy::cast_precision_loss)]
    let n = rows.len() as f64;
    let mean = |f: &dyn Fn(&Sample) -> f64| rows.iter().map(f).sum::<f64>() / n;
    let (measured, predicted) = (mean(&|r| r.0), mean(&|r| r.1));
    let (r_ratio, m_ratio) = (mean(&|r| r.2), mean(&|r| r.3));
    let heat_ratio = mean(&|r| r.4);
    println!(
        "V-D5: L/R from the phases {measured:.4}, from the energy and power {predicted:.4} ({:+.1} %); R with/without the self-field {r_ratio:.4}; magnitudes/phases {m_ratio:.4}; heat/work with/without {heat_ratio:.4}",
        100.0 * (measured / predicted - 1.0)
    );
    common::measured("V-D5.lr_phases", measured, 3);
    common::measured("V-D5.lr_energy", predicted, 3);
    common::measured("V-D5.lr_percent", 100.0 * (measured / predicted - 1.0), 1);
    common::measured("V-D5.r_ratio", r_ratio, 4);
    common::measured("V-D5.magnitudes", m_ratio, 4);
    common::measured("V-D5.heat_ratio", heat_ratio, 4);
    assert!(
        (measured / predicted - 1.0).abs() < 0.15,
        "L/R {measured} vs {predicted}"
    );
    assert!((r_ratio - 1.0).abs() < 0.15, "R ratio {r_ratio}");
    assert!(
        (m_ratio - 1.0).abs() < 0.05,
        "magnitudes against phases {m_ratio}"
    );
    assert!(
        (heat_ratio - 1.0).abs() < 0.05,
        "heat/work with/without the self-field {heat_ratio}"
    );
}

/// The disc: radius 6 about the origin, the tests' metal (`Lattice::with_disc`), in a field
/// `B₀ sin ωt` about its centre, the self-field at the speed of light `c`.
fn disc(c: f64, carriers: usize, b0: f64, omega: f64, dt: f64, seed: u64) -> Tube {
    let radius = 6.0;
    let half = DVec3::new(radius * 1.001, radius * 1.001, 0.0);
    let (lo, hi) = (-half, half);
    // The cores' density over the bounding square: those outside the outline go.
    let square = (hi.x - lo.x) * (hi.y - lo.y);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(square / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, seed ^ 31);
    let lattice = Lattice::new(lo, hi, common::metal(A), cores, carriers, seed)
        .with_disc(DVec3::ZERO, radius);
    let electrodes = common::small_floating_electrode(DVec3::new(-40.0, -40.0, 0.0));
    Tube {
        arena: Some((
            lo - DVec3::new(2.0, 2.0, 0.0),
            hi + DVec3::new(2.0, 2.0, 0.0),
        )),
        b_wave: Some(SolenoidWave {
            waveform: Waveform::Sine {
                offset: 0.0,
                amplitude: b0,
                omega,
                phase: 0.0,
            },
            axis: DVec3::ZERO,
        }),
        lattices: vec![lattice],
        // Nodes every 0.25 over ±16, one on the centre: the grid keeps the disc's symmetry
        // under quarter turns and reflections.
        darwin: Some(Darwin2d::new(127, 127, 0.25, (-16.0, -16.0), c)),
        ..Tube::new(electrodes, 0.05, dt)
    }
}

/// V-D6, the skin effect. A disc (radius R = 6) in a field `B₀ sin ωt` about its centre (B₀
/// 0.2, ω 1): its eddy currents run round it (no charge builds anywhere: the induced field is
/// along them) and keep the field out of its inside as a cylinder's,
/// `B(r) = B(R) I₀(kr)/I₀(kR)`, `k² = iωμ₀σ(ω)` (`e^{iωt}`), with the model's AC conductivity
/// σ(ω) (the Boltzmann equation's, D2; ωτ ≈ 0.27, the mean free path 0.5 against the skin
/// depth `1/|k|` ≈ 3: local), and leave it alone outside (their own field 0 beyond the rim).
/// c = 5 gives `|k|R` ≈ 2: the centre's field 0.67 of the rim's, 46° behind. The total field,
/// the solenoid's plus the self-field, at r = 0, 1.5, 3, 4.5 (each the mean of four points a
/// quarter turn apart) and 8, Fourier at ω over 20 periods after one, against the solenoid's
/// over the same samples: within 0.05 (complex, of the solenoid's; set before measuring:
/// σ(ω) is the model's to 2–3 % (D1, D2, the medium's density), which moves the centre's
/// value by ~0.03). The four points' spread: the symmetry under quarter turns.
#[test]
#[ignore = "minutes: 12 000 carriers over 20 periods"]
fn disc_skin_effect() {
    let (c, periods) = (5.0, 20.0);
    let (omega, b0, dt, carriers, radius) = (1.0, 0.2, 0.02, 12_000, 6.0);
    let mu0 = 4.0 * PI / (c * c);
    let period = 2.0 * PI / omega;
    let (start, end) = (period, (1.0 + periods) * period);
    let radii = [0.0, 1.5, 3.0, 4.5, 8.0];
    let turns = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];
    let tube = disc(c, carriers, b0, omega, dt, 47);
    let mut s = tube.state();
    // (t, the self-field at each radius and quarter turn).
    let mut samples: Vec<(f64, Vec<f64>)> = Vec::new();
    while s.t < end {
        tube.step(&mut s);
        if s.t >= start {
            let d = s.darwin.as_ref().expect("the field");
            let b = d.b_z();
            let row = radii
                .iter()
                .flat_map(|r| turns.iter().map(move |(cx, cy)| (r * cx, r * cy)))
                .map(|(x, y)| d.b_at(&b, x, y))
                .collect();
            samples.push((s.t, row));
        }
    }
    // C = (2/T) ∫ B e^{−iωt} dt, the solenoid's added to the self-field.
    #[allow(clippy::cast_precision_loss)]
    let n = samples.len() as f64;
    let amplitude = |k: Option<usize>| -> Complex64 {
        let mut c = Complex64::new(0.0, 0.0);
        for (t, row) in &samples {
            let b = b0 * libm::sin(omega * t) + k.map_or(0.0, |k| row[k]);
            c.re += b * libm::cos(omega * t);
            c.im -= b * libm::sin(omega * t);
        }
        Complex64::new(2.0 * c.re / n, 2.0 * c.im / n)
    };
    let ext = amplitude(None);
    let sigma = mobility(1.0 / MEAN_FREE_PATH, omega);
    let k = csqrt(Complex64::new(0.0, omega * mu0) * sigma);
    let rim = ci0(k * radius);
    let (mut worst, mut spread): (f64, f64) = (0.0, 0.0);
    for (i, r) in radii.iter().enumerate() {
        let at: Vec<Complex64> = (0..4).map(|q| amplitude(Some(4 * i + q)) / ext).collect();
        let mean = Complex64::new(
            at.iter().map(|z| z.re).sum::<f64>() / 4.0,
            at.iter().map(|z| z.im).sum::<f64>() / 4.0,
        );
        let want = if *r < radius {
            ci0(k * *r) / rim
        } else {
            Complex64::new(1.0, 0.0)
        };
        let err = libm::hypot(mean.re - want.re, mean.im - want.im);
        let apart = at
            .iter()
            .map(|z| libm::hypot(z.re - mean.re, z.im - mean.im))
            .fold(0.0, f64::max);
        println!(
            "V-D6 r {r}: measured ({:.3}, {:.3}) |{:.3}|, cylinder ({:.3}, {:.3}) |{:.3}|: {err:.3} (quarter turns within {apart:.3})",
            mean.re,
            mean.im,
            libm::hypot(mean.re, mean.im),
            want.re,
            want.im,
            libm::hypot(want.re, want.im)
        );
        worst = worst.max(err);
        spread = spread.max(apart);
    }
    println!(
        "V-D6: σ(ω) ({:.4}, {:.4}), k ({:.3}, {:.3}), |k|R {:.2}; worst {worst:.3}, quarter turns {spread:.3}",
        sigma.re,
        sigma.im,
        k.re,
        k.im,
        libm::hypot(k.re, k.im) * radius
    );
    measured("V-D6.worst", worst, 2);
    measured("V-D6.spread", spread, 2);
    assert!(worst < 0.05, "worst {worst}");
}

/// The measured values' formatting (DEC-22): significant digits, fixed-point in the middle
/// range, scientific beyond, the exponent taken after rounding.
#[test]
fn measured_values_keep_their_significant_digits() {
    assert_eq!(significant(0.027_3, 2), "0.027");
    assert_eq!(significant(1.071_3, 4), "1.071");
    assert_eq!(significant(0.845_6, 3), "0.846");
    assert_eq!(significant(-0.31, 1), "-0.3");
    assert_eq!(significant(1.2e-8, 2), "1.2e-8");
    assert_eq!(significant(9.996, 3), "10.0");
    assert_eq!(significant(123_456.0, 2), "1.2e5");
}

/// The symmetry tests' small metal: a disc of radius 3 about `center`, its bounding square
/// `half` each way, 1500 carriers over a background (charge per area) with charge-to-mass
/// `qm` and thermal speed `thermal`, in a field `b_static + b_wave sin(ωt)` about its centre,
/// a floating electrode at `electrode`, the self-field at the speed of light `c` (None:
/// without it) on 63 × 63 nodes 0.25 apart, centred on the disc.
#[derive(Clone)]
struct Little {
    center: DVec3,
    half: f64,
    cores: Vec<DVec3>,
    walls: bool,
    background: f64,
    qm: f64,
    thermal: f64,
    b_static: f64,
    b_wave: f64,
    omega: f64,
    electrode: DVec3,
    c: Option<f64>,
    dt: f64,
    carriers: usize,
}

impl Little {
    fn new(cores: bool) -> Self {
        let center = DVec3::ZERO;
        let half = 3.003;
        let h = DVec3::new(half, half, 0.0);
        let square = 4.0 * half * half;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = libm::round(square / (2.0 * A * MEAN_FREE_PATH)) as usize;
        Self {
            center,
            half,
            cores: if cores {
                Lattice::random_cores(center - h, center + h, count, A, 5)
            } else {
                Vec::new()
            },
            walls: true,
            background: 1.0,
            qm: -1.0,
            thermal: 1.0,
            b_static: 0.3,
            b_wave: 0.4,
            omega: 2.0,
            electrode: DVec3::new(-40.0, -40.0, 0.0),
            c: Some(3.0),
            dt: 0.02,
            carriers: 1500,
        }
    }

    fn tube(&self) -> Tube {
        let h = DVec3::new(self.half, self.half, 0.0);
        let (lo, hi) = (self.center - h, self.center + h);
        let mut lattice = Lattice::new(
            lo,
            hi,
            Material {
                background: self.background,
                core_radius: A,
                charge_per_mass: self.qm,
                thermal_speed: self.thermal,
            },
            self.cores.clone(),
            self.carriers,
            11,
        )
        .with_disc(self.center, 3.0);
        if !self.walls {
            lattice.walls.clear();
        }
        let electrodes = common::small_floating_electrode(self.electrode);
        Tube {
            arena: Some((
                lo - DVec3::new(2.0, 2.0, 0.0),
                hi + DVec3::new(2.0, 2.0, 0.0),
            )),
            b_wave: Some(SolenoidWave {
                waveform: Waveform::Sine {
                    offset: self.b_static,
                    amplitude: self.b_wave,
                    omega: self.omega,
                    phase: 0.0,
                },
                axis: self.center,
            }),
            lattices: vec![lattice],
            darwin: self.c.map(|c| {
                Darwin2d::new(63, 63, 0.25, (self.center.x - 8.0, self.center.y - 8.0), c)
            }),
            ..Tube::new(electrodes, 0.05, self.dt)
        }
    }
}

/// A run's start with its carriers moved and turned by `pos` and `vel` from another's (the
/// same carriers, in the same order: their collisions draw the same numbers).
fn mapped(
    tube: &Tube,
    from: &TubeState,
    pos: &dyn Fn(DVec3) -> DVec3,
    vel: &dyn Fn(DVec3) -> DVec3,
) -> TubeState {
    let mut s = tube.state();
    assert_eq!(s.x.len(), from.x.len());
    assert_eq!(s.id, from.id);
    for k in 0..s.x.len() {
        s.x[k] = pos(from.x[k]);
        s.v[k] = vel(from.v[k]);
        if !s.u.is_empty() {
            s.u[k] = vel(from.u[k]);
        }
    }
    s
}

/// The self-field's B_z at points about the disc's centre.
fn field_at(s: &TubeState, points: &[DVec3]) -> Vec<f64> {
    let d = s.darwin.as_ref().expect("the field");
    let b = d.b_z();
    points.iter().map(|p| d.b_at(&b, p.x, p.y)).collect()
}

/// Points inside and about the disc: a ring of 24 at radii 1, 2.5 and 4.
#[allow(clippy::cast_precision_loss)]
fn probe_points(center: DVec3) -> Vec<DVec3> {
    let mut out = Vec::new();
    for r in [1.0, 2.5, 4.0] {
        for k in 0..24 {
            let a = 2.0 * PI * (k as f64 + 0.37) / 24.0;
            out.push(center + DVec3::new(libm::cos(a), libm::sin(a), 0.0) * r);
        }
    }
    out
}

/// Y21, charge conjugation of the Darwin coupling: every charge reversed (carriers positive
/// over a negative background) and the external field reversed; with collisions. The
/// carriers' paths the same and the self-field reversed, bit for bit (every product in the
/// push, the deposit, the solve and the interpolation changes sign exactly).
#[test]
fn y21_darwin_charge_conjugation() {
    let p = Little::new(true);
    let mut q = p.clone();
    (q.background, q.qm, q.b_static, q.b_wave) = (-p.background, -p.qm, -p.b_static, -p.b_wave);
    let (a, b) = (p.tube(), q.tube());
    let (mut sa, mut sb) = (a.state(), b.state());
    for _ in 0..100 {
        a.step(&mut sa);
        b.step(&mut sb);
    }
    assert!(sa.collisions > 0);
    assert_eq!(sa.x, sb.x);
    assert_eq!(sa.v, sb.v);
    let points = probe_points(p.center);
    let (fa, fb) = (field_at(&sa, &points), field_at(&sb, &points));
    assert!(
        fa.iter()
            .zip(&fb)
            .all(|(x, y)| x.to_bits() == (-*y).to_bits())
    );
    println!(
        "Y21: {} collisions; paths identical, self-field reversed",
        sa.collisions
    );
}

/// Y22, the coupling's scaling: speeds and c doubled, the step halved, the charges (the
/// background's and so the carriers') ×4 at the same q/m, the external field and its
/// frequency doubled: the same paths in half the time, velocities and fields doubled, bit for
/// bit (all factors powers of two; the solve's tolerance relative); with collisions.
#[test]
fn y22_darwin_scaling() {
    let p = Little::new(true);
    let mut q = p.clone();
    q.background = 4.0 * p.background;
    q.thermal = 2.0 * p.thermal;
    q.b_static = 2.0 * p.b_static;
    q.b_wave = 2.0 * p.b_wave;
    q.omega = 2.0 * p.omega;
    q.c = p.c.map(|c| 2.0 * c);
    q.dt = 0.5 * p.dt;
    let (a, b) = (p.tube(), q.tube());
    let (mut sa, mut sb) = (a.state(), b.state());
    for _ in 0..100 {
        a.step(&mut sa);
        b.step(&mut sb);
    }
    assert!(sa.collisions > 0);
    assert_eq!(sa.collisions, sb.collisions);
    assert_eq!(sa.x, sb.x);
    assert!(sa.v.iter().zip(&sb.v).all(|(x, y)| *x * 2.0 == *y));
    let points = probe_points(p.center);
    let (fa, fb) = (field_at(&sa, &points), field_at(&sb, &points));
    assert!(
        fa.iter()
            .zip(&fb)
            .all(|(x, y)| (2.0 * *x).to_bits() == y.to_bits())
    );
    assert_eq!(sa.t.to_bits(), (2.0 * sb.t).to_bits());
    println!(
        "Y22: {} collisions; the same paths in half the time",
        sa.collisions
    );
}

/// The largest distance between a run's carriers mapped and another's, and between their
/// self-fields (the second's at the mapped points, `sign` the field's parity), of the first's
/// largest field.
fn mismatch(
    sa: &TubeState,
    sb: &TubeState,
    pos: &dyn Fn(DVec3) -> DVec3,
    sign: f64,
    center: DVec3,
) -> (f64, f64) {
    let dx =
        sa.x.iter()
            .zip(&sb.x)
            .map(|(a, b)| (pos(*a) - *b).length())
            .fold(0.0, f64::max);
    let points = probe_points(center);
    let moved: Vec<DVec3> = points.iter().map(|p| pos(*p)).collect();
    let (fa, fb) = (field_at(sa, &points), field_at(sb, &moved));
    let top = fa.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let db = fa
        .iter()
        .zip(&fb)
        .map(|(a, b)| (sign * a - b).abs())
        .fold(0.0, f64::max);
    (dx, db / top)
}

/// A map of the plane (positions or velocities).
type Map<'a> = &'a dyn Fn(DVec3) -> DVec3;

/// Y23, the coupling's quarter turn, mirror and translation: the collisionless disc (its walls,
/// a static and a swinging field, the floating electrode) turned a quarter about its centre,
/// mirrored in x (B_z a pseudoscalar: the external field reversed, the self-field must come
/// out reversed), and moved by whole cells (1, −2) with its grid; 100 steps (2 time units).
/// Against the original run mapped: positions within 1e-8, self-fields within 1e-8 of the
/// largest (the grid shares the disc's symmetries; round-off from the sums' orders and the
/// outline's sines, amplified by the carriers' encounters).
#[test]
fn y23_darwin_rotation_mirror_translation() {
    let p = Little::new(false);
    let a = p.tube();
    let s0 = a.state();
    let mut sa = s0.clone();
    for _ in 0..100 {
        a.step(&mut sa);
    }
    let turn = |x: DVec3| DVec3::new(-x.y, x.x, 0.0);
    let mirror = |x: DVec3| DVec3::new(-x.x, x.y, 0.0);
    let shift = DVec3::new(1.0, -2.0, 0.0);
    let moved = |x: DVec3| x + shift;
    let same = |x: DVec3| x;
    let mut worst: (f64, f64) = (0.0, 0.0);
    let cases: [(&str, Map, Map, f64); 3] = [
        ("quarter turn", &turn, &turn, 1.0),
        ("mirror", &mirror, &mirror, -1.0),
        ("translation", &moved, &same, 1.0),
    ];
    for (name, pos, vel, sign) in cases {
        let mut q = p.clone();
        q.center = pos(p.center);
        q.electrode = pos(p.electrode);
        q.b_static = sign * p.b_static;
        q.b_wave = sign * p.b_wave;
        let b = q.tube();
        let mut sb = mapped(&b, &s0, pos, vel);
        for _ in 0..100 {
            b.step(&mut sb);
        }
        let (dx, db) = mismatch(&sa, &sb, pos, sign, p.center);
        println!("Y23 {name}: positions within {dx:.2e}, self-field within {db:.2e}");
        worst = (worst.0.max(dx), worst.1.max(db));
    }
    assert!(worst.0 < 1e-8 && worst.1 < 1e-8, "{worst:?}");
}

/// The collisionless disc without walls (the background's well holds its carriers; its
/// square 4.5 each way) in a static field `b`, run `n` steps forward, then reversed (each
/// velocity, canonical velocity and vector potential reversed, the field reversed; leapfrog
/// positions and velocities half a step apart, so the turn is taken after one more step from
/// the positions before it) and run `n` back: how far each carrier ends from its start, at
/// most; and how far the self-field moved them over the forward run (against c = ∞).
fn reversal(c: Option<f64>, dt: f64, n: usize) -> (f64, f64) {
    let mut p = Little::new(false);
    (p.walls, p.half, p.b_wave, p.c, p.dt) = (false, 4.5, 0.0, c, dt);
    let forward = p.tube();
    let mut s = forward.state();
    let start = s.x.clone();
    let mut free = p.clone();
    free.c = None;
    let free_tube = free.tube();
    let mut sf = free_tube.state();
    for _ in 0..n {
        forward.step(&mut s);
        free_tube.step(&mut sf);
    }
    let effect =
        s.x.iter()
            .zip(&sf.x)
            .map(|(a, b)| (*a - *b).length())
            .fold(0.0, f64::max);
    let turn = s.x.clone();
    forward.step(&mut s);
    let mut q = p.clone();
    q.b_static = -p.b_static;
    let back = q.tube();
    s.x = turn;
    for v in &mut s.v {
        *v = -*v;
    }
    for u in &mut s.u {
        *u = -*u;
    }
    if let Some(a) = s.self_a.as_mut() {
        a.x.iter_mut().chain(a.y.iter_mut()).for_each(|z| *z = -*z);
    }
    if let Some(d) = s.darwin.as_mut() {
        d.chi.iter_mut().for_each(|z| *z = -*z);
    }
    for _ in 0..n {
        back.step(&mut s);
    }
    let back_err =
        s.x.iter()
            .zip(&start)
            .map(|(a, b)| (*a - *b).length())
            .fold(0.0, f64::max);
    (back_err, effect)
}

/// Y24, time reversal: the collisionless disc in a static field (0.3), 0.2 time units forward
/// and back (inside the time the carriers' encounters amplify any difference by orders of
/// magnitude: over 2 units round-off grows 10⁵-fold). Without the self-field the leapfrog is
/// reversible: back at the start within 1e-8 (over 2 units). With it (c = 3) the splitting is
/// first order (the canonical kick at the step's start position, the drift with the velocity
/// there, `v = u − (q/m)A(x)` depending on x): the return error falling at least 1.6 times when
/// the step is halved, and within the coupling's 5 % target (V-D6) of the self-field's effect
/// on the paths over that time at dt = 0.01 (at 0.02 it is 7.8 %: PHYSICS.md §2.12).
#[test]
fn y24_darwin_time_reversal() {
    let (plain, _) = reversal(None, 0.02, 100);
    let (coarse, effect) = reversal(Some(3.0), 0.02, 10);
    let (fine, effect_fine) = reversal(Some(3.0), 0.01, 20);
    println!(
        "Y24: without the self-field back within {plain:.2e}; with it {coarse:.3e} (dt 0.02) and {fine:.3e} (0.01), the self-field's effect {effect:.3e} and {effect_fine:.3e}: ratios {:.4} and {:.4}, falling {:.2}×",
        coarse / effect,
        fine / effect_fine,
        coarse / fine
    );
    assert!(plain < 1e-8, "plain {plain}");
    assert!(fine < 0.05 * effect_fine, "{fine} against {effect_fine}");
    assert!(coarse / fine >= 1.6, "{coarse} then {fine}");
}

/// The electrostatic energy of the carriers among themselves (the engine's rods: `−q q' ln r²`
/// beyond the softening, its parabola inside) and with the backgrounds, less constants.
fn electrostatic(tube: &Tube, x: &[DVec3], charge: &[f64]) -> f64 {
    let h2 = tube.softening * tube.softening;
    let mut e = 0.0;
    for i in 0..x.len() {
        for j in (i + 1)..x.len() {
            let (dx, dy) = (x[i].x - x[j].x, x[i].y - x[j].y);
            let r2 = dx * dx + dy * dy;
            let g = if r2 >= h2 {
                -libm::log(r2)
            } else {
                -(r2 / h2 - 1.0 + libm::log(h2))
            };
            e += charge[i] * charge[j] * g;
        }
        for l in &tube.lattices {
            e += charge[i] * l.background_at(x[i]).0;
        }
    }
    e
}

/// A closed system holding a current: the collisionless disc of `carriers` held by its
/// background alone (no walls: a bounce changes v at once, which the field can only follow by
/// losing energy, below), in a static field 0.3 (no work), its carriers turned at 0.3 about the
/// centre at the start. Its total energy `E_kin + E_es + U` at whole steps (the kinetic energy
/// the mean of the half steps' about each, U the mean of the solves' about it) from t = 0.3 to
/// 3.1: the drift per time unit, the mean U and the kinetic energy.
fn energy_drift(c: Option<f64>, carriers: usize, dt: f64, walls: bool) -> (f64, f64, f64) {
    let mut p = Little::new(false);
    (p.c, p.dt, p.b_wave, p.carriers) = (c, dt, 0.0, carriers);
    if !walls {
        (p.walls, p.half) = (false, 4.5);
    }
    let tube = p.tube();
    let mut s = tube.state();
    for k in 0..s.x.len() {
        let w = DVec3::new(-s.x[k].y, s.x[k].x, 0.0) * 0.3;
        s.v[k] += w;
        if !s.u.is_empty() {
            s.u[k] += w;
        }
    }
    let kinetic = |v: &[DVec3], s: &TubeState| {
        v.iter()
            .zip(&s.charge)
            .zip(&s.qm)
            .map(|((v, q), qm)| 0.5 * (q / qm).abs() * v.length_squared())
            .sum::<f64>()
    };
    let magnetic = |s: &TubeState| s.darwin.as_ref().map_or(0.0, |d| d.magnetic_energy());
    let total = |s: &mut TubeState| {
        let (x, v, u) = (s.x.clone(), s.v.clone(), magnetic(s));
        tube.step(s);
        let um = 0.5 * (u + magnetic(s));
        let ek = 0.5 * (kinetic(&v, s) + kinetic(&s.v, s));
        (ek + electrostatic(&tube, &x, &s.charge) + um, um, ek)
    };
    while s.t < 0.3 - 0.5 * dt {
        tube.step(&mut s);
    }
    let (t0, (e0, _, _)) = (s.t, total(&mut s));
    while s.t < 3.1 - 0.5 * dt {
        tube.step(&mut s);
    }
    let (t1, (e1, u1, k1)) = (s.t, total(&mut s));
    ((e1 - e0) / (t1 - t0), u1, k1)
}

/// Y25, the coupling's energy (the Darwin limit's: `E_kin + E_es + ½∫B²/μ₀` constant in a
/// closed system), on `energy_drift`'s spun disc over 2.8 time units. Without the self-field
/// the leapfrog's: within 1e-3 of the kinetic energy per time unit. With it (c = 3, U ≈ 10)
/// the drift per time unit within the coupling's 5 % target (V-D6) of U at 1500 carriers
/// (dt 0.02), the walled disc (its carriers bouncing off the rim) too, and falling at least
/// 1.6 times when both the step halves and the carriers quadruple: the drift is the
/// discretisation's, and goes with it. (Until 2026-10-27 the criterion was 1.6 times with the
/// carriers alone: the impulses then lost each carrier's own inductance, an error `∝ 1/N` that
/// dominated, 3.5 % and 2.1 %, and the walled disc lost 14.6 %; with a uniform medium's
/// self-inductance added (DEC-24) the drift fitted `α dt + β/N`, α ≈ 80 % of U per time per
/// unit dt, β ≈ 980 %, the walled disc 4.5 %. With a velocity set at once entering the next
/// solve as a fixed current (DEC-25) the carriers' count no longer matters: the step's first
/// order is all of it, and the spin is the carriers' current from the start, U ≈ 10 against
/// 1.4 when it went into `u` and was screened.)
#[test]
fn y25_darwin_energy() {
    let (plain, _, k_plain) = energy_drift(None, 1500, 0.02, false);
    let (coarse, u, k) = energy_drift(Some(3.0), 1500, 0.02, false);
    let (fine, u_fine, _) = energy_drift(Some(3.0), 6000, 0.02, false);
    let (half_dt, u_half, _) = energy_drift(Some(3.0), 1500, 0.01, false);
    let (refined, u_refined, _) = energy_drift(Some(3.0), 6000, 0.01, false);
    let (walled, u_walled, _) = energy_drift(Some(3.0), 1500, 0.02, true);
    let pct = |d: f64, u: f64| 100.0 * d.abs() / u;
    println!(
        "Y25: without the self-field {plain:.2e} per time ({:.1e} of the kinetic {k_plain:.2}); with it, % of U per time: dt 0.02 {:.2} at 1500 carriers (U = {u:.3}, kinetic {k:.2}), {:.2} at 6000; dt 0.01 {:.2} and {:.2}: refined {:.2}×; walled {:.2}",
        plain.abs() / k_plain,
        pct(coarse, u),
        pct(fine, u_fine),
        pct(half_dt, u_half),
        pct(refined, u_refined),
        coarse / refined,
        pct(walled, u_walled)
    );
    measured("Y25.coarse_percent", pct(coarse, u), 2);
    measured("Y25.fine_percent", pct(fine, u_fine), 2);
    measured("Y25.half_dt_percent", pct(half_dt, u_half), 2);
    measured("Y25.refined_percent", pct(refined, u_refined), 2);
    measured("Y25.walled_percent", pct(walled, u_walled), 2);
    assert!(plain.abs() < 1e-3 * k_plain, "plain {plain}");
    assert!(coarse.abs() < 0.05 * u, "{coarse} of {u}");
    assert!(
        walled.abs() < 0.05 * u_walled,
        "walled {walled} of {u_walled}"
    );
    assert!(coarse / refined >= 1.6, "{coarse} then {refined}");
}

/// The energy budget with collisions: `energy_drift`'s disc at rest (thermal carriers, its
/// walls) with its lattice's cores, `E_kin + E_es + ½∫B²/μ₀` plus the heat the lattice booked,
/// per time unit over t = 0.3–6.1, and the kinetic energy.
fn collision_budget(c: Option<f64>) -> (f64, f64) {
    let dt = 0.02;
    let mut p = Little::new(true);
    (p.c, p.dt, p.b_wave, p.carriers) = (c, dt, 0.0, 1500);
    let tube = p.tube();
    let mut s = tube.state();
    let kinetic = |v: &[DVec3], s: &TubeState| {
        v.iter()
            .zip(&s.charge)
            .zip(&s.qm)
            .map(|((v, q), qm)| 0.5 * (q / qm).abs() * v.length_squared())
            .sum::<f64>()
    };
    let magnetic = |s: &TubeState| s.darwin.as_ref().map_or(0.0, |d| d.magnetic_energy());
    let total = |s: &mut TubeState| {
        let (x, v, u, j) = (s.x.clone(), s.v.clone(), magnetic(s), s.joule);
        tube.step(s);
        let um = 0.5 * (u + magnetic(s));
        let ek = 0.5 * (kinetic(&v, s) + kinetic(&s.v, s));
        (
            ek + electrostatic(&tube, &x, &s.charge) + um + 0.5 * (j + s.joule),
            ek,
        )
    };
    while s.t < 0.3 - 0.5 * dt {
        tube.step(&mut s);
    }
    let (t0, (e0, _)) = (s.t, total(&mut s));
    while s.t < 6.1 - 0.5 * dt {
        tube.step(&mut s);
    }
    let (t1, (e1, k1)) = (s.t, total(&mut s));
    ((e1 - e0) / (t1 - t0), k1)
}

/// Y26, the coupling's energy with collisions: `collision_budget`'s disc, its carriers
/// colliding with the lattice's cores, the total and the lattice's booked heat together
/// constant within 1e-3 of the kinetic energy per time unit, Y25's bound for the leapfrog,
/// without the self-field and with it (c = 3, dt 0.02). Set 2026-10-27 after the measurement
/// that found the impulses' self-inductance 13 % short (DEC-24, superseded by DEC-25): the
/// disc then lost 5.1e-3 of its kinetic energy per time unit, the lattice making it up as
/// heat it booked negative.
#[test]
fn y26_darwin_energy_with_collisions() {
    let (plain, k_plain) = collision_budget(None);
    let (coupled, k) = collision_budget(Some(3.0));
    println!(
        "Y26: E + heat drifts by {:.2e} of the kinetic energy {k_plain:.2} per time without the self-field, {:.2e} of {k:.2} with it",
        plain / k_plain,
        coupled / k
    );
    measured("Y26.plain", (plain / k_plain).abs(), 2);
    measured("Y26.coupled", (coupled / k).abs(), 2);
    assert!(plain.abs() < 1e-3 * k_plain, "plain {plain} of {k_plain}");
    assert!(coupled.abs() < 1e-3 * k, "coupled {coupled} of {k}");
}

/// V-D7's tube: a closed secondary 20 × 20 and, with `primary`, a closed primary 14 × 14 inside
/// it (3–17), both rings with a wall of 1 (thin against the skin depth ≈ 2 at ω = 0.8, as
/// V-D5's ring), of the tests' metal at V-D5's 55.5 carriers per area, in a field
/// `b0 sin ωt` about their common centre (10, 10); the self-field at c = 3.
fn transformer(primary: bool, b0: f64, omega: f64, dt: f64, seed: u64) -> Tube {
    let ring = |lo: DVec3, side: f64, carriers: usize, seed: u64| {
        let hi = lo + DVec3::new(side, side, 0.0);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = libm::round(side * side / (2.0 * A * MEAN_FREE_PATH)) as usize;
        let cores = Lattice::random_cores(lo, hi, count, A, seed ^ 31);
        Lattice::new(lo, hi, common::metal(A), cores, carriers, seed).with_hole(
            lo + DVec3::new(1.0, 1.0, 0.0),
            hi - DVec3::new(1.0, 1.0, 0.0),
        )
    };
    let mut lattices = vec![ring(DVec3::ZERO, 20.0, 4218, seed)];
    if primary {
        lattices.push(ring(DVec3::new(3.0, 3.0, 0.0), 14.0, 2886, seed + 1));
    }
    // One small floating electrode far away (the electrostatics needs one; OP-1).
    let electrodes = common::small_floating_electrode(DVec3::new(-40.0, -40.0, 0.0));
    Tube {
        arena: Some((DVec3::new(-2.0, -2.0, 0.0), DVec3::new(22.0, 22.0, 0.0))),
        b_wave: Some(SolenoidWave {
            waveform: Waveform::Sine {
                offset: 0.0,
                amplitude: b0,
                omega,
                phase: 0.0,
            },
            axis: DVec3::new(10.0, 10.0, 0.0),
        }),
        lattices,
        // Nodes 0.25 apart over −6…25.5.
        darwin: Some(Darwin2d::new(127, 127, 0.25, (-6.0, -6.0), 3.0)),
        ..Tube::new(electrodes, 0.05, dt)
    }
}

/// One V-D7 run to `end`, sampled every `every`: (t, the self-field in the gap between the
/// rings, the mean of four points midway along its sides, and the self-field's whole flux
/// `∫B dA`).
fn transformer_run(tube: &Tube, end: f64, every: f64) -> Vec<(f64, f64, f64)> {
    let gap = [(18.0, 10.0), (2.0, 10.0), (10.0, 18.0), (10.0, 2.0)];
    let mut s = tube.state();
    let mut out = Vec::new();
    let mut next = 0.0;
    while s.t < end {
        tube.step(&mut s);
        if s.t >= next {
            let d = s.darwin.as_ref().expect("the field");
            let b = d.b_z();
            let g = gap.iter().map(|&(x, y)| d.b_at(&b, x, y)).sum::<f64>() / 4.0;
            out.push((s.t, g, b.iter().sum::<f64>() * d.h * d.h));
            next += every;
        }
    }
    out
}

/// V-D7, the transformer. Two closed thin rings, the primary inside the secondary, in a field
/// `0.125 sin(0.8 t)`; runs with and without the primary (the same secondary, seed by seed).
/// In z-invariance the primary's field lies wholly inside it, so all of its flux `Φ₁` threads
/// the secondary, whose circuit law is `Z₂ I₂ = −iω D S₂ − iω Φ₁` (`D = −i b0` the drive's
/// phasor, `S₂` the area its current links) and `Z₂ I₂′ = −iω D S₂` without the primary.
/// Without assuming the impedance, inductances or resistance: `G − G′ = Φ₁ G′/(D S₂)`, `G = μ₀
/// I₂` the field in the gap between the rings (inside every filament of the secondary's
/// current), `Φ₁` the self-field's whole flux less the secondary's own `G S₂`. `S₂` is measured,
/// not assumed: the secondary's own flux per gap field without the primary, `F′/G′` (the area
/// its current links, weighted by how the current spreads across the wall: the same weight for
/// the drive's flux and its own). (First taken as a current uniform across the wall,
/// `(a² + ab + b²)/3` = 361: the prediction moved by its own size between the wall's sides,
/// 324 and 400.) Lock-in at ω over 8 periods after one, two seeds. Required (set before
/// measuring, from V-D5's noise: each phasor ~3 %, the coupling ~30 % of G, so one seed's
/// difference to ~14 %): within 20 % of the measured `G − G′`, the seeds' mean.
#[test]
#[ignore = "eight minutes: four runs of 4200 and 7100 carriers over 9 periods"]
fn transformer_secondary_sees_the_primary_flux() {
    let (b0, omega, dt) = (0.125, 0.8, 0.02);
    let period = 2.0 * PI / omega;
    let drive = Complex64::new(0.0, -b0);
    let lock_in = |rows: &[(f64, f64, f64)]| {
        let mut g = Complex64::new(0.0, 0.0);
        let mut f = Complex64::new(0.0, 0.0);
        let mut n = 0.0;
        for &(t, gap, flux) in rows.iter().filter(|r| r.0 >= period) {
            let e = Complex64::new(libm::cos(omega * t), -libm::sin(omega * t));
            g += e * gap;
            f += e * flux;
            n += 1.0;
        }
        (g * (2.0 / n), f * (2.0 / n))
    };
    let mut measured_sum = Complex64::new(0.0, 0.0);
    let mut predicted_sum = Complex64::new(0.0, 0.0);
    for seed in [61_u64, 67] {
        let end = 9.0 * period;
        let alone = transformer_run(&transformer(false, b0, omega, dt, seed), end, 0.1);
        let both = transformer_run(&transformer(true, b0, omega, dt, seed), end, 0.1);
        let (g_alone, flux_alone) = lock_in(&alone);
        let (g_both, flux_both) = lock_in(&both);
        let s2 = flux_alone / g_alone;
        let phi1 = flux_both - g_both * s2;
        let measured = g_both - g_alone;
        let predicted = phi1 * g_alone / (drive * s2);
        println!(
            "V-D7 seed {seed}: S₂ {:.1}∠{:.1}°; G′ {:.4}∠{:.1}°, G {:.4}∠{:.1}°, Φ₁ {:.2}∠{:.1}°; G − G′ {:.4}∠{:.1}° against Φ₁G′/(DS₂) {:.4}∠{:.1}°",
            s2.norm(),
            s2.arg().to_degrees(),
            g_alone.norm(),
            g_alone.arg().to_degrees(),
            g_both.norm(),
            g_both.arg().to_degrees(),
            phi1.norm(),
            phi1.arg().to_degrees(),
            measured.norm(),
            measured.arg().to_degrees(),
            predicted.norm(),
            predicted.arg().to_degrees()
        );
        measured_sum += measured;
        predicted_sum += predicted;
    }
    let err = (measured_sum - predicted_sum).norm() / measured_sum.norm();
    println!(
        "V-D7: G − G′ against Φ₁G′/(DS₂), the seeds' mean: {:.1} % apart",
        100.0 * err
    );
    measured("V-D7.percent", 100.0 * err, 2);
    assert!(err < 0.20, "{err}");
}

/// Y27's tube: a lead circuit round a metal ring, in a grounded envelope (−2.5…22.5 square:
/// the electrostatics' reference, as every level's; without one, the open world's free-space
/// gauge gives two circuit nodes' common charge a negative capacitance, OP-1, and the circuit
/// ran away). Two small boxes 12 apart at y = 18.5, each with its lead from the ground
/// terminal (10, −2.5) on the envelope, along it and up (a loop 12 × 20.5, closed across the
/// boxes' gap by their capacitance); the left one driven by `30 sin(0.3 t)` (an ideal source
/// behind its lead), the right one grounded. Inside the loop a closed ring 8 × 8 (wall 1) of
/// the tests' metal at V-D5's 55.5 carriers per area; the self-field at c = 3.
fn lead_circuit(dt: f64) -> Tube {
    use physics::bem::Bias;
    use physics::tube::Circuit;
    use physics::zinv::{Electrode, Electrodes, Section};
    let (lo, hi) = (DVec3::new(6.0, 6.0, 0.0), DVec3::new(14.0, 14.0, 0.0));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(64.0 / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 71);
    let ring = Lattice::new(lo, hi, common::metal(A), cores, 1554, 73).with_hole(
        lo + DVec3::new(1.0, 1.0, 0.0),
        hi - DVec3::new(1.0, 1.0, 0.0),
    );
    let rect = |cx: f64| Section::Rect {
        center: DVec3::new(cx, 18.5, 0.0),
        angle: 0.0,
        half_length: 1.0,
        half_thickness: 0.5,
    };
    let electrodes = Electrodes::new(
        vec![
            Electrode {
                section: rect(4.0),
                bias: Bias::Grounded,
            },
            Electrode {
                section: rect(16.0),
                bias: Bias::Grounded,
            },
            Electrode {
                section: Section::Frame {
                    center: DVec3::new(10.0, 10.0, 0.0),
                    half_width: 12.5,
                    half_height: 12.5,
                    segment: 0.5,
                },
                bias: Bias::Grounded,
            },
        ],
        &[],
        0.25,
    );
    let g = DVec3::new(10.0, -2.5, 0.0);
    let lead = |x: f64| Some(vec![g, DVec3::new(x, -2.5, 0.0), DVec3::new(x, 18.0, 0.0)]);
    Tube {
        arena: Some((DVec3::new(-2.5, -2.5, 0.0), DVec3::new(22.5, 22.5, 0.0))),
        circuits: vec![Circuit {
            electrode: 0,
            source: Waveform::Sine {
                offset: 0.0,
                amplitude: 30.0,
                omega: 0.3,
                phase: 0.0,
            },
            resistance: 0.0,
            capacitance: 0.0,
            inductance: 0.0,
            to: None,
        }],
        lattices: vec![ring],
        darwin: Some(Darwin2d::new(127, 127, 0.25, (-6.0, -6.0), 3.0)),
        leads: vec![lead(4.0), lead(16.0)],
        ..Tube::new(electrodes, 0.05, dt)
    }
}

/// `lead_circuit`'s energy budget from t = 0.3 to `end`: the drift per time unit of
/// `E_kin + E_es + ½VᵀCV + ½∫B²/μ₀ + heat − W`, `W` the source's work as the circuit books it
/// (`½ (S₀ + S₁) h (I₀ + I₁)` a step), and the mean field energy; the electrodes' energy from
/// their potentials and capacitances (the ring, neutral and far, induces little on them).
fn lead_budget(dt: f64, end: f64) -> (f64, f64) {
    let tube = lead_circuit(dt);
    let source = tube.circuits[0].source;
    let mut s = tube.state();
    let cap = tube.electrodes.capacitance.clone();
    let kinetic = |v: &[DVec3], s: &TubeState| {
        v.iter()
            .zip(&s.charge)
            .zip(&s.qm)
            .map(|((v, q), qm)| 0.5 * (q / qm).abs() * v.length_squared())
            .sum::<f64>()
    };
    let magnetic = |s: &TubeState| s.darwin.as_ref().map_or(0.0, |d| d.magnetic_energy());
    let electrodes = |s: &TubeState| {
        let v = &s.potentials;
        let mut w = 0.0;
        for i in 0..v.len() {
            for j in 0..v.len() {
                w += 0.5 * v[i] * cap[i][j] * v[j];
            }
        }
        w
    };
    let current = |s: &TubeState| s.leads.as_ref().expect("the leads").current[0];
    let mut work = 0.0;
    let step = |s: &mut TubeState, work: &mut f64| {
        let (t0, i0) = (s.t, current(s));
        tube.step(s);
        let i1 = current(s);
        *work += 0.5 * (source.value(t0) + source.value(s.t)) * 0.5 * dt * (i0 + i1);
    };
    // The total at a whole step: the kinetic energy the mean of the half steps' about it, the
    // field's the mean of the solves' about it.
    let total = |s: &mut TubeState, work: &mut f64| {
        let (x, v, u, j, we) = (
            s.x.clone(),
            s.v.clone(),
            magnetic(s),
            s.joule,
            electrodes(s),
        );
        let w0 = *work;
        step(s, work);
        let um = 0.5 * (u + magnetic(s));
        let ek = 0.5 * (kinetic(&v, s) + kinetic(&s.v, s));
        let e = ek + electrostatic(&tube, &x, &s.charge) + we + um + 0.5 * (j + s.joule)
            - 0.5 * (w0 + *work);
        (e, um)
    };
    while s.t < 0.3 - 0.5 * dt {
        step(&mut s, &mut work);
    }
    let (t0, (e0, _)) = (s.t, total(&mut s, &mut work));
    let (mut field, mut n) = (0.0, 0.0);
    while s.t < end - 0.5 * dt {
        step(&mut s, &mut work);
        field += magnetic(&s);
        n += 1.0;
    }
    let (t1, (e1, _)) = (s.t, total(&mut s, &mut work));
    ((e1 - e0) / (t1 - t0), field / n)
}

/// Y27, the leads' circuits' energy (DEC-27): a lead circuit driven round a metal ring
/// (`lead_circuit`), its supply's work against the ring's heat and every energy, the field's
/// included, over 1.3 periods of the drive: the drift per time unit within 5 % of the mean
/// field energy, Y25's bound, at dt 0.01 and 0.005 (set before measuring). The circuit takes
/// the leads' own inductances implicitly and the ring's flux along them a step behind
/// (`physics::tube` `leads.rs`): the lag's error is what this bounds.
#[test]
#[ignore = "four minutes: a ring of 1554 carriers with a lead circuit, 27 time units at two steps"]
fn y27_lead_circuit_energy() {
    let end = 27.6;
    let (coarse, u) = lead_budget(0.01, end);
    let (fine, u_fine) = lead_budget(0.005, end);
    println!(
        "Y27: the budget drifts by {coarse:.3e} per time at dt 0.01 ({:.3} % of the mean field energy {u:.2}), {fine:.3e} at 0.005 ({:.3} % of {u_fine:.2})",
        100.0 * coarse.abs() / u,
        100.0 * fine.abs() / u_fine
    );
    measured("Y27.coarse_percent", 100.0 * coarse.abs() / u, 2);
    measured("Y27.fine_percent", 100.0 * fine.abs() / u_fine, 2);
    assert!(coarse.abs() < 0.05 * u, "{coarse} of {u}");
    assert!(fine.abs() < 0.05 * u_fine, "{fine} of {u_fine}");
}
