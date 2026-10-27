use super::plots::ticks;
use super::*;

/// Axis ticks: round steps (1, 2 or 5 times a power of ten) inside the range, 0 exact.
#[test]
fn axis_ticks_are_round() {
    assert_eq!(ticks(0.0, 10.0, 5), vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    assert_eq!(ticks(-0.3, 1.2, 4), vec![0.0, 0.5, 1.0]);
    assert_eq!(ticks(-1.0, 1.0, 4), vec![-1.0, -0.5, 0.0, 0.5, 1.0]);
    assert_eq!(ticks(0.0, 200.0, 5), vec![0.0, 50.0, 100.0, 150.0, 200.0]);
    assert_eq!(ticks(3.0, 3.0, 4), vec![3.0]);
}

#[test]
fn si_round_trip() {
    for v in [2e6, 1.5e3, 1e-6, -4e6, 0.5, 12.0] {
        let back = parse_si(&fmt_si(v)).unwrap();
        assert!(
            (back - v).abs() <= 1e-12 * v.abs(),
            "{v} -> {} -> {back}",
            fmt_si(v)
        );
    }
    assert_eq!(parse_si("3u").map(f64::to_bits), Some(3e-6f64.to_bits()));
    assert_eq!(parse_si("1e6").map(f64::to_bits), Some(1e6f64.to_bits()));
}
