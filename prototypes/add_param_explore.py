"""Adds a temporary parameterized D8 exploration (EM_STRIPS, EM_B0, EM_PERIOD, EM_SKIP, EM_END)."""
p = "crates/physics/tests/drude.rs"
s = open(p, encoding="utf-8").read()
s += '''
/// TEMP exploration: D8 with its parameters from the environment.
#[test]
#[ignore = "exploration"]
fn eddy_param_explore() {
    let get = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let strips = get("EM_STRIPS", 1.0) as u32;
    let (b0, period, skip, end) = (get("EM_B0", 0.5), get("EM_PERIOD", 20.0), get("EM_SKIP", 2.0), get("EM_END", 60.0));
    let (a, n, dt) = (16.0, 20_000, 0.02);
    let sigma = mobility(1.0 / MEAN_FREE_PATH, 0.0).0;
    let bdot = 4.0 * b0 / period;
    let blocks = eddy(a, strips, n, dt, b0, period, end, skip);
    #[allow(clippy::cast_precision_loss)]
    let m = blocks.len() as f64;
    let power = blocks.iter().map(|b| b.0).sum::<f64>() / m;
    let current = blocks.iter().map(|b| b.1).sum::<f64>() / m;
    let sp = libm::sqrt(blocks.iter().map(|b| (b.0 - power).powi(2)).sum::<f64>() / ((m - 1.0) * m));
    let w = a / f64::from(strips);
    let p_theory = f64::from(strips) * sigma * bdot * bdot * torsion_area(w, a);
    let j_theory = sigma * bdot * torsion_midline(w, a);
    println!(
        "strips {strips}, b0 {b0}, period {period}, skip {skip}: {m} ramps: power {power:.5} ± {sp:.5} vs {p_theory:.5} ({:+.2} %), current {current:.5} vs {j_theory:.5} ({:+.2} %)",
        100.0 * (power / p_theory - 1.0),
        100.0 * (current / j_theory - 1.0)
    );
}
'''
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
