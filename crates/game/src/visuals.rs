//! Visual-only computations: potential map and field lines. Nothing here feeds back into
//! gameplay.

use physics::DVec3;
use physics::field::{Coulomb, FieldSolver};
use physics::trajectory::Scenario;

pub struct PotentialImage {
    pub width: u32,
    pub height: u32,
    pub rgba: Vec<u8>,
}

/// Potential map over the world bounds, in the plane z = 0, from the particle's point of
/// view: colour shows `U = q(φ − φ(A))` relative to the launch kinetic energy `T₀`
/// (red: uphill, blue: downhill), with one-pixel contours every `T₀/4`. Regions with
/// `U > T₀` are classically forbidden (energy conservation, exact also relativistically)
/// and are drawn dark; their boundary, the turning line `U = T₀`, is drawn bright.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
pub fn potential_image(scn: &Scenario<Coulomb>, px_per_cell: u32) -> PotentialImage {
    let b = scn.bounds.expect("level scenarios have bounds");
    let size = b.max - b.min;
    let width = (size.x * f64::from(px_per_cell)).round() as u32;
    let height = (size.y * f64::from(px_per_cell)).round() as u32;
    let kin = physics::dynamics::Kinematics::new(scn.particle.mass, scn.c);
    let t0 = kin.kinetic_energy(scn.p0).max(1e-12);
    let q = scn.particle.charge;
    let phi_a = scn.field.sample(scn.x0, 0.0).phi;

    // U / T0 at pixel centres (NaN inside charges). Row 0 is the top of the image.
    let (w, h) = (width as usize, height as usize);
    let mut u = vec![0.0f64; w * h];
    for j in 0..h {
        let y = b.max.y - (j as f64 + 0.5) / f64::from(px_per_cell);
        for i in 0..w {
            let x = b.min.x + (i as f64 + 0.5) / f64::from(px_per_cell);
            let p = DVec3::new(x, y, 0.0);
            let inside = scn.obstacles.iter().any(|s| s.signed_distance(p) < 0.0);
            u[j * w + i] = if inside {
                f64::NAN
            } else {
                q * (scn.field.sample(p, 0.0).phi - phi_a) / t0
            };
        }
    }

    // Contour band index (integer-valued).
    #[allow(clippy::cast_possible_truncation)]
    let band = |v: f64| (v * 4.0).floor() as i64;
    let mut rgba = Vec::with_capacity(w * h * 4);
    for j in 0..h {
        for i in 0..w {
            let v = u[j * w + i];
            // A contour passes between this pixel and its right or lower neighbour.
            let mut contour = None;
            for (ni, nj) in [(i + 1, j), (i, j + 1)] {
                if ni < w && nj < h {
                    let nv = u[nj * w + ni];
                    if v.is_finite() && nv.is_finite() && band(v) != band(nv) {
                        let turning = (v > 1.0) != (nv > 1.0);
                        // Inside the forbidden region only the turning line is drawn.
                        if turning || v <= 1.0 {
                            contour = Some(contour.unwrap_or(false) || turning);
                        }
                    }
                }
            }
            rgba.extend_from_slice(&colour(v, contour));
        }
    }
    PotentialImage {
        width,
        height,
        rgba,
    }
}

/// `contour`: `None` for no line, `Some(true)` for the turning line `U = T₀`.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn colour(u: f64, contour: Option<bool>) -> [u8; 4] {
    if !u.is_finite() {
        return [20, 20, 24, 255];
    }
    let s = (u / 1.5).tanh();
    let base = [0.10, 0.11, 0.14];
    let mut c = if s >= 0.0 {
        [base[0] + 0.55 * s, base[1] + 0.10 * s, base[2] + 0.02 * s]
    } else {
        [
            base[0] + 0.02 * -s,
            base[1] + 0.18 * -s,
            base[2] + 0.55 * -s,
        ]
    };
    if u > 1.0 {
        // Forbidden: darken strongly.
        c = [c[0] * 0.3, c[1] * 0.3, c[2] * 0.3];
    }
    match contour {
        Some(true) => c = [0.95, 0.95, 0.85],
        Some(false) => c = [c[0] + 0.10, c[1] + 0.10, c[2] + 0.10],
        None => {}
    }
    let to_u8 = |v: f64| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    [to_u8(c[0]), to_u8(c[1]), to_u8(c[2]), 255]
}

/// A traced field line with arrowheads along it.
pub struct FieldLine {
    pub points: Vec<DVec3>,
    /// Arrowhead positions and unit directions of **E** there.
    pub arrows: Vec<(DVec3, DVec3)>,
}

/// Arc length between arrowheads, in cells.
const ARROW_SPACING: f64 = 3.0;
/// Largest direction change allowed per step (radians): keeps curves smooth.
const MAX_TURN: f64 = 0.02;

/// Field lines in the plane z = 0, starting on a small circle around each charge, with a
/// number of lines proportional to |Q| (`lines_for_max` for the largest charge). Traced
/// with RK4 on the unit field direction, with a step that shrinks where the line bends.
/// Lines from negative charges that end on a positive charge are dropped, because the
/// positive charges already draw them. In the 2D slice of a 3D field, line density does
/// not represent field strength (SPEC §4): lines show direction only.
pub fn field_lines(
    scn: &Scenario<Coulomb>,
    charges: &[(DVec3, f64)],
    lines_for_max: u32,
) -> Vec<FieldLine> {
    let q_max = charges.iter().map(|c| c.1.abs()).fold(0.0, f64::max);
    if q_max == 0.0 {
        return Vec::new();
    }
    let b = scn.bounds.expect("bounds");
    let dir_at = |x: DVec3, sign: f64| -> Option<DVec3> {
        let e = scn.field.sample(x, 0.0).e;
        (e.length() > 0.0 && e.is_finite()).then(|| e.normalize() * sign)
    };
    let mut lines = Vec::new();
    for (k, &(pos, q)) in charges.iter().enumerate() {
        let r0 = scn.obstacles[k].radius * 1.2;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = ((f64::from(lines_for_max) * q.abs() / q_max).round() as u32).max(2);
        // Trace away from the charge: along E from positive, against E from negative.
        let sign = q.signum();
        for i in 0..n {
            let a = std::f64::consts::TAU * (f64::from(i) + 0.5) / f64::from(n);
            let mut x = pos + DVec3::new(a.cos(), a.sin(), 0.0) * r0;
            let mut points = vec![x];
            let mut ends_on: Option<usize> = None;
            let mut ds: f64 = 0.02;
            let mut length = 0.0;
            'trace: for _ in 0..20_000 {
                let Some(d0) = dir_at(x, sign) else { break };
                // RK4 step, retried with a smaller step while the line turns too fast.
                let (x_new, d_new) = loop {
                    let k1 = d0;
                    let Some(k2) = dir_at(x + k1 * (0.5 * ds), sign) else {
                        break 'trace;
                    };
                    let Some(k3) = dir_at(x + k2 * (0.5 * ds), sign) else {
                        break 'trace;
                    };
                    let Some(k4) = dir_at(x + k3 * ds, sign) else {
                        break 'trace;
                    };
                    let xn = x + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (ds / 6.0);
                    let Some(dn) = dir_at(xn, sign) else {
                        break 'trace;
                    };
                    let turn = d0.angle_between(dn);
                    if turn > MAX_TURN && ds > 1e-3 {
                        ds *= 0.5;
                        continue;
                    }
                    if turn < 0.3 * MAX_TURN {
                        ds = (ds * 1.5).min(0.5);
                    }
                    break (xn, dn);
                };
                let _ = d_new;
                length += (x_new - x).length();
                x = x_new;
                points.push(x);
                if let Some(j) = scn
                    .obstacles
                    .iter()
                    .enumerate()
                    .position(|(j, s)| j != k && s.signed_distance(x) < 0.0)
                {
                    ends_on = Some(j);
                    break;
                }
                let outside = x.x < b.min.x || x.y < b.min.y || x.x > b.max.x || x.y > b.max.y;
                if outside || length > 500.0 {
                    break;
                }
            }
            if q < 0.0 && ends_on.is_some_and(|j| charges[j].1 > 0.0) {
                continue;
            }
            let arrows = arrows_along(&points, |x| dir_at(x, 1.0));
            lines.push(FieldLine { points, arrows });
        }
    }
    lines
}

/// Arrowheads every `ARROW_SPACING` of arc length (the first at half the spacing), each
/// pointing along E at its position.
fn arrows_along(points: &[DVec3], e_dir: impl Fn(DVec3) -> Option<DVec3>) -> Vec<(DVec3, DVec3)> {
    let mut out = Vec::new();
    let mut next = 0.5 * ARROW_SPACING;
    let mut s = 0.0;
    for w in points.windows(2) {
        let seg = (w[1] - w[0]).length();
        while seg > 0.0 && s + seg >= next {
            let p = w[0].lerp(w[1], (next - s) / seg);
            if let Some(d) = e_dir(p) {
                out.push((p, d));
            }
            next += ARROW_SPACING;
        }
        s += seg;
    }
    out
}
