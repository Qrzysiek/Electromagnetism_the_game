//! Location of the first zero crossing of an event function on one integration step
//! (PHYSICS.md §6).
//!
//! Event functions are signed distances, which are 1-Lipschitz in position, so
//! `|dg/dt| ≤ v_max` along the trajectory. On an interval `[a, b]` with `g(a), g(b) > 0`
//! this gives the lower bound `min g ≥ (g(a) + g(b) − v_max (b − a)) / 2`. If the bound is
//! positive the interval is certified crossing-free; otherwise it is bisected, left half
//! first, so the *earliest* crossing is found. A grazing pass that dips below zero between
//! two samples therefore cannot be missed, as long as `v_max` bounds the speed.
//!
//! The crossing time is resolved to adjacent floating-point numbers: the result is the
//! first representable time at which the interpolated `g` is `≤ 0`.

/// Outcome of the search on one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossingSearch {
    /// First time with `g ≤ 0`, if any.
    pub time: Option<f64>,
    /// Smallest value of `g` among all points evaluated (including the endpoints).
    pub min_sampled: f64,
}

/// Finds the first `t ∈ (ta, tb]` with `g(t) ≤ 0`, given `ga = g(ta) > 0`, `gb = g(tb)`,
/// and `|dg/dt| ≤ v_max` on the interval.
pub fn first_crossing(
    g: &mut impl FnMut(f64) -> f64,
    ta: f64,
    tb: f64,
    ga: f64,
    gb: f64,
    v_max: f64,
) -> CrossingSearch {
    debug_assert!(ga > 0.0);
    let mut min_sampled = ga.min(gb);
    let time = search(g, ta, tb, ga, gb, v_max, &mut min_sampled);
    CrossingSearch { time, min_sampled }
}

fn search(
    g: &mut impl FnMut(f64) -> f64,
    ta: f64,
    tb: f64,
    ga: f64,
    gb: f64,
    v_max: f64,
    min_sampled: &mut f64,
) -> Option<f64> {
    if gb > 0.0 && ga + gb - v_max * (tb - ta) > 0.0 {
        return None;
    }
    let tm = 0.5 * (ta + tb);
    if tm <= ta || tm >= tb {
        // Interval cannot be subdivided further.
        return (gb <= 0.0).then_some(tb);
    }
    let gm = g(tm);
    *min_sampled = min_sampled.min(gm);
    if gm <= 0.0 {
        return search(g, ta, tm, ga, gm, v_max, min_sampled).or(Some(tm));
    }
    search(g, ta, tm, ga, gm, v_max, min_sampled)
        .or_else(|| search(g, tm, tb, gm, gb, v_max, min_sampled))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_simple_root_to_float_resolution() {
        // g(t) = 1 - t, root at t = 1.
        let r = first_crossing(&mut |t| 1.0 - t, 0.0, 2.0, 1.0, -1.0, 1.0);
        let t = r.time.unwrap();
        assert!(1.0 - t <= 0.0 && 1.0 - f64::from_bits(t.to_bits() - 1) > 0.0);
    }

    #[test]
    fn finds_dip_between_positive_endpoints() {
        // |g'| ≤ 2 on [0, 1]; the dip below zero lies between the endpoints.
        let g = |t: f64| (t - 0.5).powi(2) - 1e-6;
        let r = first_crossing(&mut { g }, 0.0, 1.0, g(0.0), g(1.0), 2.0);
        let t = r.time.expect("dip missed");
        assert!((t - (0.5 - 1e-3)).abs() < 1e-12);
    }

    #[test]
    fn certifies_near_miss() {
        let g = |t: f64| (t - 0.5).powi(2) + 1e-6;
        let r = first_crossing(&mut { g }, 0.0, 1.0, g(0.0), g(1.0), 2.0);
        assert_eq!(r.time, None);
        assert!(r.min_sampled < 1e-5);
    }

    #[test]
    fn returns_earliest_of_several_crossings() {
        let g = |t: f64| (10.0 * t).cos() + 0.5;
        let r = first_crossing(&mut { g }, 0.0, 3.0, g(0.0), g(3.0), 10.0);
        let expected = (2.0 * std::f64::consts::PI / 3.0) / 10.0;
        assert!((r.time.unwrap() - expected).abs() < 1e-13);
    }
}
