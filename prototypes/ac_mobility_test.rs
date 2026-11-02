/// `mobility`'s complex number, minimal: (re, im).
#[derive(Clone, Copy, Debug)]
struct C(f64, f64);

impl C {
    fn mul(self, o: C) -> C {
        C(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn div(self, o: C) -> C {
        let d = o.0 * o.0 + o.1 * o.1;
        C(
            (self.0 * o.0 + self.1 * o.1) / d,
            (self.1 * o.0 - self.0 * o.1) / d,
        )
    }
    fn abs(self) -> f64 {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }
}

/// The Boltzmann equation's mobility of the Drude metal's carriers, as D2 takes it
/// (`physics/tests/common` `mobility`, PHYSICS.md §2.11): carriers of thermal speed 1
/// colliding at `κ v` with hard discs (`⟨cos θ⟩ = −1/3`), their speeds drawn anew from the
/// flux's law at each collision, driven `∝ e^{iωt}`: `(π/κ) P/(1 − ⟨cos θ⟩ Q)`. A thermal
/// speed `s` is `κ s` for `κ` (the time scales with `1/(κ s)`).
fn mobility(kappa: f64, omega: f64) -> C {
    use std::f64::consts::PI;
    let cos_theta = -1.0 / 3.0;
    let (n, top) = (28_000, 14.0);
    let h = top / f64::from(n);
    let (mut p, mut q) = (C(0.0, 0.0), C(0.0, 0.0));
    for k in 0..=n {
        let v = h * f64::from(k);
        let w = if k == 0 || k == n {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let d = kappa * kappa * v * v + omega * omega;
        if d == 0.0 {
            continue;
        }
        let (re, im) = (kappa * kappa * v / d, -kappa * omega / d);
        let e = (-0.5 * v * v).exp();
        let fp = v * v * v * e / (2.0 * PI);
        let fq = v * (2.0 / PI).sqrt() * v * v * e;
        p = C(p.0 + w * fp * re, p.1 + w * fp * im);
        q = C(q.0 + w * fq * re, q.1 + w * fq * im);
    }
    let (p, q) = (C(p.0 * h / 3.0, p.1 * h / 3.0), C(q.0 * h / 3.0, q.1 * h / 3.0));
    let den = C(1.0 - cos_theta * q.0, -cos_theta * q.1);
    p.div(den).mul(C(PI / kappa, 0.0))
}

/// The least-squares fit of `a + Re(Y e^{iωt})` to samples, each the mean of the signal
/// over `[t0, t1]` (`t0 = t1`: a value at an instant): `Y`.
fn fit(samples: &[(f64, f64, f64)], omega: f64) -> C {
    // Basis: 1, cos, sin, each averaged over the sample's span.
    let basis = |t0: f64, t1: f64| -> [f64; 3] {
        if t1 - t0 < 1e-12 {
            [1.0, (omega * t0).cos(), (omega * t0).sin()]
        } else {
            let span = omega * (t1 - t0);
            [
                1.0,
                ((omega * t1).sin() - (omega * t0).sin()) / span,
                ((omega * t0).cos() - (omega * t1).cos()) / span,
            ]
        }
    };
    let mut m = [[0.0; 3]; 3];
    let mut r = [0.0; 3];
    for &(t0, t1, y) in samples {
        let b = basis(t0, t1);
        for i in 0..3 {
            r[i] += b[i] * y;
            for j in 0..3 {
                m[i][j] += b[i] * b[j];
            }
        }
    }
    // Gaussian elimination (3 × 3, well conditioned over whole periods).
    for i in 0..3 {
        for k in i + 1..3 {
            let f = m[k][i] / m[i][i];
            for j in 0..3 {
                m[k][j] -= f * m[i][j];
            }
            r[k] -= f * r[i];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let s: f64 = (i + 1..3).map(|j| m[i][j] * x[j]).sum();
        x[i] = (r[i] - s) / m[i][i];
    }
    // c cos + s sin = Re(Y e^{iωt}) with Y = c − i s.
    C(x[1], -x[2])
}

/// The AC conductivity through a level (DRUDE.md Phase 2 item 2): Ohm's law's strip (16 × 4,
/// q/m −1, thermal speed s = 2, mean free path 1, here 200 macroparticles per area: 12 800)
/// driven by a sine whose frequency its supply sets, amplitude 10 (`a τ/s ≈ 0.13`: the
/// field's heating, second order, a per cent or two). Its cores ten times smaller than the
/// levels' (radius 0.001, `a²n = 0.0005`: the medium's own correction to the Boltzmann
/// equation, `1.1 a²n ln(a²n)`, −0.4 % against −2.9 % at the levels' 0.005; PHYSICS.md
/// §2.11). Over the strip's middle (the central three fifths along it, a quarter from its
/// sides: `TubeRun::hall_forces`' region): the carriers' mean drift along it from frame to
/// frame, and the mean field along it, the potential's drop across the region (seven lines
/// across it, `TubeRun::field_at`: smooth where the field at points is grainy), each fitted
/// with a sine at ω from t = 5 over whole periods, give the mobility `μ(ω) = v / ((q/m) E)`,
/// averaged over three seeds of the carriers and cores. Against the Boltzmann equation's
/// (D2's `mobility`, at `κ s = 2` for `κ = 1/ℓ = 1`), at ω = 0.8 and 2.4 (`ω/κs` = 0.4 and
/// 1.2: the drift 12° and 30° behind the field). Required: within 5 % of its size (the
/// target for the model's transport coefficients, as D1 and D2).
#[test]
#[ignore = "about 4 min"]
fn ac_mobility_through_the_level() {
    let base = {
        let mut level = level::shipped("ohms_law");
        let spec = level.tube.as_mut().expect("tube");
        spec.metals[0].carriers_per_area = 200.0;
        spec.metals[0].core_radius = 0.001;
        // The right contact on an ideal sine its supply sets.
        spec.circuits = vec![TubeCircuit {
            electrode: 1,
            source: Source::Sine {
                offset: 0.0,
                amplitude: 0.0,
                omega: 1.0,
                phase: 0.0,
            },
            resistance: Some(0.0),
            capacitance: 0.0,
            supply_offset: false,
            linked: Vec::new(),
            inductance: 0.0,
            return_to: None,
        }];
        spec.goals.clear();
        level.limits.supply_omegas = vec![0.1, 10.0];
        level.limits.supply_voltages = vec![-40.0, 40.0];
        level
    };
    let kappa_s = 2.0;
    let q_over_m = -1.0;
    let phase = |c: C| c.1.atan2(c.0).to_degrees();
    let mut worst: f64 = 0.0;
    for (omega, periods) in [(0.8, 16.0), (2.4, 24.0)] {
        let supply = Element {
            omega: Some(omega),
            ..Element::supply(base.electrodes[1].center, 10.0)
        };
        base.check_placement(&[supply]).expect("placement");
        let (t0, t1) = (5.0, 5.0 + periods * std::f64::consts::TAU / omega);
        let mut sum = C(0.0, 0.0);
        let seeds = [1_u64, 2, 3];
        for seed in seeds {
            let mut level = base.clone();
            level.tube.as_mut().expect("tube").metals[0].seed = seed;
            let run = level.tube_run(&[supply], 1, t1 + 0.5);
            let (centre, axis, normal, length, width) = run.metals[0].strip();
            let half = 0.3 * length;
            let lines: Vec<f64> = (-3..=3).map(|k| 0.5 * f64::from(k)).collect();
            assert!(lines.iter().all(|y| y.abs() < 0.5 * width - 0.25));
            let mut drift = Vec::new();
            let mut field = Vec::new();
            for pair in run.frames.windows(2) {
                let (before, now) = (&pair[0], &pair[1]);
                if before.t < t0 || now.t > t1 {
                    continue;
                }
                let then: std::collections::HashMap<u32, [f32; 2]> =
                    before.id.iter().zip(&before.x).map(|(i, x)| (*i, *x)).collect();
                let (mut sum, mut n) = (0.0, 0.0);
                for (x, id) in now.x.iter().zip(&now.id) {
                    let p = physics::DVec3::new(f64::from(x[0]), f64::from(x[1]), 0.0);
                    let d = p - centre;
                    if d.dot(axis).abs() >= half || d.dot(normal).abs() >= 0.5 * width - 0.25 {
                        continue;
                    }
                    let Some(x0) = then.get(id) else { continue };
                    let p0 = physics::DVec3::new(f64::from(x0[0]), f64::from(x0[1]), 0.0);
                    sum += (p - p0).dot(axis);
                    n += 1.0;
                }
                drift.push((before.t, now.t, sum / n / (now.t - before.t)));
                // The mean field along the axis over the region: the potential's drop.
                let drop: f64 = lines
                    .iter()
                    .map(|&y| {
                        let at = |u: f64| run.field_at(now, centre + axis * u + normal * y).0;
                        at(-half) - at(half)
                    })
                    .sum::<f64>()
                    / lines.len() as f64;
                field.push((now.t, now.t, drop / (2.0 * half)));
            }
            let v = fit(&drift, omega);
            let e = fit(&field, omega);
            let mu = v.div(C(q_over_m * e.0, q_over_m * e.1));
            println!(
                "  ω = {omega}, seed {seed}: μ {:.4} ∠{:.1}° ({} drift samples, field \
                 amplitude {:.4} ∠{:.1}°, drift {:.4})",
                mu.abs(),
                phase(mu),
                drift.len(),
                e.abs(),
                phase(e),
                v.abs()
            );
            sum = C(sum.0 + mu.0, sum.1 + mu.1);
        }
        let n = seeds.len() as f64;
        let measured = C(sum.0 / n, sum.1 / n);
        let theory = mobility(kappa_s, omega);
        let error = C(measured.0 - theory.0, measured.1 - theory.1).abs() / theory.abs();
        println!(
            "ω = {omega} (ω/κs = {:.2}): μ {:.4} ∠{:.1}° against {:.4} ∠{:.1}°: {:.2} % of its \
             size",
            omega / kappa_s,
            measured.abs(),
            phase(measured),
            theory.abs(),
            phase(theory),
            100.0 * error,
        );
        worst = worst.max(error);
    }
    assert!(worst < 0.05, "worst {:.2} %", 100.0 * worst);
}

#[test]
#[ignore = "diagnostic"]
fn zz_frame_timing() {
    let mut level = level::shipped("ohms_law");
    let spec = level.tube.as_mut().expect("tube");
    spec.circuits = vec![TubeCircuit {
        electrode: 1,
        source: Source::Sine { offset: 0.0, amplitude: 0.0, omega: 1.0, phase: 0.0 },
        resistance: Some(0.0),
        capacitance: 0.0,
        supply_offset: false,
        linked: Vec::new(),
        inductance: 0.0,
        return_to: None,
    }];
    spec.goals.clear();
    level.limits.supply_omegas = vec![0.1, 10.0];
    level.limits.supply_voltages = vec![-40.0, 40.0];
    let omega = 2.4;
    let supply = Element { omega: Some(omega), ..Element::supply(level.electrodes[1].center, 10.0) };
    let run = level.tube_run(&[supply], 1, 6.0);
    let step = level.tube.as_ref().unwrap().step;
    for f in run.frames.iter().take(8) {
        let v = f64::from(f.potentials[1]);
        let expect = 10.0 * (omega * f.t).sin();
        // The time at which the drive had this value (nearest).
        let shift = ((v / 10.0).asin() / omega) - f.t;
        println!("t {:.4} V {:.5} expected {:.5} (diff {:.2e}; frames {} steps apart?) shift {:.4} (step {step})", f.t, v, expect, v - expect, (f.t / step).round(), shift);
    }
}
