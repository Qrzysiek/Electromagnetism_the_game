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
#[ignore = "minutes: four runs of 4200 and 7100 carriers over 9 periods"]
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
