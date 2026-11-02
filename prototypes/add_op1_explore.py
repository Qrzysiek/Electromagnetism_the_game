"""Adds a temporary OP-1 exploration: the particles' strip with the electrode anywhere."""
p = "crates/physics/tests/drude.rs"
s = open(p, encoding="utf-8").read()
s += '''
/// TEMP exploration (OP-1): a centred rectangle's left-half current and power over the ramp
/// (EM_W, EM_B, EM_PERIOD, EM_B0, EM_END; the grounded electrode at EM_EX, EM_EY).
#[test]
#[ignore = "exploration"]
fn op1_explore() {
    let get = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
    let (w, b) = (get("EM_W", 8.0), get("EM_B", 32.0));
    let (period, b0, end) = (get("EM_PERIOD", 20.0), get("EM_B0", 0.25), get("EM_END", 60.0));
    let (ex, ey) = (get("EM_EX", -1.5), get("EM_EY", -1.5));
    let (carriers, dt) = (20_000_usize, 0.02);
    let (lo, hi) = (DVec3::ZERO, DVec3::new(w, b, 0.0));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = libm::round(w * b / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, 31);
    let lattice = Lattice::new(lo, hi, 1.0, cores, A, -1.0, 1.0, carriers, 47);
    let electrodes = Electrodes::new(
        vec![Electrode {
            section: Section::Rect {
                center: DVec3::new(ex, ey, 0.0),
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
            low: -b0,
            high: b0,
            delay: 0.0,
            rise: 0.5 * period,
            width: 0.0,
            fall: 0.5 * period,
            period,
        },
        axis: DVec3::new(0.5 * w, 0.5 * b, 0.0),
    };
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: Some((lo - DVec3::new(2.0, 2.0, 0.0), hi + DVec3::new(2.0, 2.0, 0.0))),
        b_z: 0.0,
        b_wave: Some(wave),
        circuits: Vec::new(),
        lattices: vec![lattice],
    };
    let mut st = tube.state();
    let n0 = st.x.len();
    let bins = 20;
    let (mut left, mut power, mut count_b) = (vec![0.0; bins], vec![0.0; bins], vec![0.0; bins]);
    while st.t < end {
        let t0 = st.t;
        let before: std::collections::HashMap<u32, DVec3> =
            st.id.iter().zip(&st.x).map(|(id, x)| (*id, *x)).collect();
        tube.step(&mut st);
        let phase = t0 - period * libm::floor(t0 / period);
        let rising = phase < 0.5 * period;
        let since = if rising { phase } else { phase - 0.5 * period };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = ((since / (0.5 * period) * bins as f64) as usize).min(bins - 1);
        let sign = if rising { 1.0 } else { -1.0 };
        for ((id, x), q) in st.id.iter().zip(&st.x).zip(&st.charge) {
            let Some(x0) = before.get(id) else { continue };
            let d = *x - *x0;
            power[k] += q * wave.induced(t0, *x0).dot(d) / dt;
            if x0.x < 0.5 * w {
                left[k] += sign * q * d.y / dt;
            }
        }
        count_b[k] += 1.0;
    }
    let row = |v: &[f64]| (0..bins).map(|k| format!("{:.2}", v[k] / count_b[k])).collect::<Vec<_>>().join(" ");
    println!(
        "particles {w} x {b}, electrode at ({ex}, {ey}), ramps of {}, b0 {b0}: carriers {n0} -> {}, lost {:.3}",
        0.5 * period,
        st.x.len(),
        st.lost
    );
    println!("  left-half current: {}", row(&left));
    println!("  power: {}", row(&power));
}
'''
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
