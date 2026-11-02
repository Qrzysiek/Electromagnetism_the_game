"""Replaces the leak exploration with a power-over-the-ramp profile."""
p = "crates/physics/tests/drude.rs"
s = open(p, encoding="utf-8").read()
i = s.index("/// TEMP exploration: carriers crossing between touching strips.")
s = s[:i] + '''/// TEMP exploration: the induced field's power over each ramp, and each strip's net current
/// along it.
#[test]
#[ignore = "exploration"]
fn eddy_profile_explore() {
    let (a, carriers, dt, period) = (16.0, 20_000_usize, 0.02, 20.0);
    for strips in [1_u32, 4] {
        let w = a / f64::from(strips);
        let lattices: Vec<Lattice> = (0..strips)
            .map(|k| {
                let lo = DVec3::new(w * f64::from(k), 0.0, 0.0);
                let hi = lo + DVec3::new(w, a, 0.0);
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let count = libm::round(w * a / (2.0 * A * MEAN_FREE_PATH)) as usize;
                let cores = Lattice::random_cores(lo, hi, count, A, 31 + u64::from(k));
                Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers / strips as usize, 47 + u64::from(k))
            })
            .collect();
        let electrodes = Electrodes::new(
            vec![Electrode {
                section: Section::Rect {
                    center: DVec3::new(-1.5, -1.5, 0.0),
                    angle: 0.0,
                    half_length: 0.25,
                    half_thickness: 0.25,
                },
                bias: Bias::Grounded,
            }],
            &[],
            0.25,
        );
        let wave = SolenoidWave {
            waveform: Waveform::Pulse {
                low: -0.5,
                high: 0.5,
                delay: 0.0,
                rise: 0.5 * period,
                width: 0.0,
                fall: 0.5 * period,
                period,
            },
            axis: DVec3::new(0.5 * a, 0.5 * a, 0.0),
        };
        let tube = Tube {
            electrodes,
            emitters: Vec::new(),
            projectiles: Vec::new(),
            softening: 0.05,
            dt,
            arena: Some((DVec3::new(-2.0, -2.0, 0.0), DVec3::new(a + 2.0, a + 2.0, 0.0))),
            b_z: 0.0,
            b_wave: Some(wave),
            circuits: Vec::new(),
            lattices,
        };
        let mut st = tube.state();
        let bins = 20;
        let mut power = vec![0.0; bins];
        let mut along = vec![vec![0.0; bins]; strips as usize];
        let mut count = vec![0.0; bins];
        while st.t < 60.0 {
            let t0 = st.t;
            let before: std::collections::HashMap<u32, DVec3> =
                st.id.iter().zip(&st.x).map(|(id, x)| (*id, *x)).collect();
            tube.step(&mut st);
            let phase = t0 - period * libm::floor(t0 / period);
            let rising = phase < 0.5 * period;
            let since = if rising { phase } else { phase - 0.5 * period };
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let b = ((since / (0.5 * period) * bins as f64) as usize).min(bins - 1);
            let sign = if rising { 1.0 } else { -1.0 };
            for ((id, x), q) in st.id.iter().zip(&st.x).zip(&st.charge) {
                let Some(x0) = before.get(id) else { continue };
                let d = *x - *x0;
                power[b] += q * wave.induced(t0, *x0).dot(d) / dt;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let k = ((x0.x / w) as usize).min(strips as usize - 1);
                along[k][b] += sign * q * d.y / dt;
            }
            count[b] += 1.0;
        }
        println!("{strips} strip(s): power over the ramp (bins of 0.5):");
        let row: Vec<String> = (0..bins).map(|b| format!("{:.2}", power[b] / count[b])).collect();
        println!("  {}", row.join(" "));
        for (k, a) in along.iter().enumerate() {
            let row: Vec<String> = (0..bins).map(|b| format!("{:.2}", a[b] / count[b])).collect();
            println!("  strip {k} net current along y: {}", row.join(" "));
        }
    }
}
'''
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
