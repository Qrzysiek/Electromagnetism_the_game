//! Numbers and states as the panel writes them: SI prefixes, potentials, outcomes,
//! boundaries, energies, disturbances.

use super::*;

/// Compact number with an SI-style suffix (2e6 → "2M", 1e-6 → "1µ").
pub fn fmt_si(v: f64) -> String {
    let a = v.abs();
    let (scale, suffix) = [
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "k"),
        (1.0, ""),
        (1e-3, "m"),
        (1e-6, "µ"),
        (1e-9, "n"),
    ]
    .into_iter()
    .find(|(s, _)| a >= *s)
    .unwrap_or((1.0, ""));
    let m = v / scale;
    let text = format!("{m:.3}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    format!("{text}{suffix}")
}

/// Parses a number with an optional SI-style suffix ("2M", "1.5k", "3µ" or "3u", "1e6").
pub fn parse_si(text: &str) -> Option<f64> {
    let t = text.trim();
    let (num, scale) = match t.chars().last()? {
        'G' => (&t[..t.len() - 1], 1e9),
        'M' => (&t[..t.len() - 1], 1e6),
        'k' => (&t[..t.len() - 1], 1e3),
        'm' => (&t[..t.len() - 1], 1e-3),
        'u' => (&t[..t.len() - 1], 1e-6),
        'µ' => (&t[..t.len() - 'µ'.len_utf8()], 1e-6),
        'n' => (&t[..t.len() - 1], 1e-9),
        _ => (t, 1.0),
    };
    num.trim().parse::<f64>().ok().map(|v| v * scale)
}

/// Shortens `s` to at most `n` characters, with an ellipsis.
pub(super) fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        s.chars().take(n - 1).collect::<String>() + "…"
    }
}

/// Decimals that resolve a direction tolerance of `half` degrees, for it and the arrival
/// angle: one from 1° up, one more per decade below (0.05° shows as 0.050°).
pub(super) fn angle_decimals(half: f64) -> usize {
    // At most 6: finer tolerances are not met by any flight anyway.
    let (mut digits, mut step) = (1, 1.0);
    while half > 0.0 && half < step && digits < 6 {
        step /= 10.0;
        digits += 1;
    }
    digits
}

pub(super) fn outcome_text(game: &Game, o: Outcome) -> String {
    match o {
        Outcome::Arrived => "reached the detector".into(),
        Outcome::Rejected => {
            "entered the detector outside its acceptance (direction or energy)".into()
        }
        Outcome::Collided(i) => format!("hit {}", obstacle_name(game, i)),
        Outcome::LeftBounds => "left the map".into(),
        Outcome::Timeout => "ran out of time".into(),
        Outcome::Failed(e) => format!("integration failed ({e})"),
        Outcome::SkippedGate(i) => format!("entered the detector without passing gate {}", i + 1),
    }
}

pub(super) fn boundary_text(b: Boundary) -> &'static str {
    match b {
        Boundary::Obstacle(_) => "a charge",
        Boundary::Bounds => "the map edge",
        Boundary::Detector => "the detector edge",
        Boundary::Acceptance => "the detector's direction/energy window",
        Boundary::TimeLimit => "the time limit",
        Boundary::Gate(_) => "a gate's edge",
        Boundary::GateAcceptance(_) => "a gate's direction/energy window",
    }
}

/// Tube levels: the goal electrode's current, preview and verified, against the goal,
/// and its collected charge over time (the window's ends marked).
/// A tube level's numbers (currents, their errors, times): plain decimals with about
/// three significant digits (SI prefixes, as `fmt_si` writes them, read oddly here:
/// "850m" for 0.85).
pub(crate) fn fmt_num(v: f64) -> String {
    let a = v.abs();
    if a != 0.0 && !(1e-3..1e4).contains(&a) {
        format!("{v:.2e}")
    } else if a >= 100.0 {
        format!("{v:.0}")
    } else if a >= 10.0 {
        format!("{v:.1}")
    } else if a >= 1.0 {
        format!("{v:.2}")
    } else {
        format!("{v:.3}")
    }
}

/// A launch energy: plain for ordinary values, with SI prefixes for the large ones of
/// scaled beams.
pub(super) fn energy_text(t: f64) -> String {
    if t.abs() >= 1e4 {
        fmt_si(t)
    } else {
        format!("{t:.3}")
    }
}

/// A potential for display: "grounded" or its value.
pub fn fmt_potential(v: f64) -> String {
    if v == 0.0 {
        "grounded".into()
    } else {
        format!("V = {}", fmt_si(v))
    }
}

/// One-line description of a disturbance.
pub fn disturbance_text(d: &level::Disturbance) -> String {
    let mut parts = Vec::new();
    if d.e != [0.0, 0.0] {
        parts.push(format!(
            "stray E = ({}, {})",
            fmt_si(d.e[0]),
            fmt_si(d.e[1])
        ));
    }
    if d.bz != 0.0 {
        parts.push(format!("stray B_z = {}", fmt_si(d.bz)));
    }
    for w in &d.waves {
        parts.push(if w.omega == 0.0 {
            format!(
                "static field {} at {:.0}°",
                fmt_si(w.amplitude),
                w.direction_deg + 90.0
            )
        } else {
            format!(
                "wave E₀ = {}, ω = {}, towards {:.0}°, phase {:.0}°",
                fmt_si(w.amplitude),
                fmt_si(w.omega),
                w.direction_deg,
                w.phase_deg
            )
        });
    }
    if parts.is_empty() {
        "undisturbed".into()
    } else {
        parts.join("; ")
    }
}

/// A fraction as a percentage with sensible digits.
pub(super) fn percent(f: f64) -> String {
    let p = f * 100.0;
    if p >= 10.0 {
        format!("{p:.0} %")
    } else if p >= 0.1 {
        format!("{p:.1} %")
    } else {
        format!("{p:.1e} %")
    }
}
