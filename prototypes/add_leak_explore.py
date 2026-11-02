"""Adds a temporary test: do carriers cross between touching strips?"""
p = "crates/physics/tests/drude.rs"
s = open(p, encoding="utf-8").read()
s += '''
/// TEMP exploration: carriers crossing between touching strips.
#[test]
#[ignore = "exploration"]
fn strips_leak_explore() {
    let (a, strips, carriers, dt) = (16.0, 4_u32, 20_000_usize, 0.02);
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
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: Some((DVec3::new(-2.0, -2.0, 0.0), DVec3::new(a + 2.0, a + 2.0, 0.0))),
        b_z: 0.0,
        b_wave: None,
        circuits: Vec::new(),
        lattices,
    };
    let mut s = tube.state();
    let strip = |x: DVec3| libm::floor(x.x / w);
    let start: std::collections::HashMap<u32, f64> =
        s.id.iter().zip(&s.x).map(|(id, x)| (*id, strip(*x))).collect();
    while s.t < 10.0 {
        tube.step(&mut s);
    }
    let moved = s
        .id
        .iter()
        .zip(&s.x)
        .filter(|(id, x)| start.get(id).is_some_and(|k| *k != strip(**x)))
        .count();
    println!("{moved} of {} carriers in another strip after t = 10", s.x.len());
}
'''
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
