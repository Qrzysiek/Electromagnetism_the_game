//! Far-zone radiation of a computed flight: the energy radiated per unit solid angle in a
//! direction, in total and per unit frequency (PHYSICS.md §3.4; Jackson §14.2, §14.5).
//!
//! Units `k = 1`, force `q(E + v×B)` (as `lienard.rs`), where Jackson's Gaussian formulas
//! hold unchanged. With `β = v/c`, `β̇ = a/c`, `κ = 1 − n·β`:
//!
//! ```text
//! dW/dΩ      = (q²/4πc) ∫ |n × ((n − β) × β̇)|² / κ⁵ dt            (Liénard, 14.38)
//! d²I/dω dΩ  = (q²/4π²c) |∫ n × ((n − β) × β̇) / κ² e^{iω(t − n·r/c)} dt|²   (14.65)
//! ```
//!
//! and `∫₀^∞ d²I/dω dΩ dω = dW/dΩ` (Parseval). The integrals run over the flight: the
//! particle moves uniformly before its launch and after its end (the acceleration form
//! then needs no boundary terms, and no radiation is attributed to the launch or to
//! the detector absorbing it). In the plane `z = 0` the directions `n` are in-plane.
//!
//! The time integral of 14.65 is done piecewise between dense samples of the flight, with
//! the amplitude and the phase `ω(t − n·r/c)` linear on each piece and integrated exactly
//! (Filon's idea): the pieces are short enough that the phase advances by at most
//! `PHASE_STEP` at the top frequency, so the result does not depend on the phase being
//! resolved by the samples themselves.

use glam::DVec3;
use std::f64::consts::PI;

/// Largest phase advance per piece at the highest frequency of a spectrum (radians). The
/// pieces are integrated exactly for a linear phase, so 1 rad is enough (test S3's high
/// band: 8.5e-5 with 1 rad, 2.1e-4 with 0.5, 7.8e-5 with 0.25); half the cost of 0.5.
pub const PHASE_STEP: f64 = 1.0;

/// A sample of the flight: time, position, velocity, acceleration.
pub type Emission = (f64, DVec3, DVec3, DVec3);

/// A radiation goal: the energy per steradian radiated into an arc of in-plane directions
/// (averaged over the arc), in all frequencies or in a band, must lie in a window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiationWindow {
    /// Unit axis of the arc, in the plane.
    pub axis: DVec3,
    /// Half-width of the arc, radians (0: the axis alone).
    pub half_angle: f64,
    /// Frequency band `[ω_min, ω_max]`; `None`: all frequencies.
    pub band: Option<(f64, f64)>,
    /// Allowed energy per steradian `[min, max]`.
    pub energy: (f64, f64),
}

impl RadiationWindow {
    /// Directions sampling the arc, at most 2° apart (at most 91), symmetric about the
    /// axis, an odd number (for Simpson's rule).
    pub fn directions(&self) -> Vec<DVec3> {
        let base = self.axis.y.atan2(self.axis.x);
        if self.half_angle <= 0.0 {
            return vec![DVec3::new(base.cos(), base.sin(), 0.0)];
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = ((self.half_angle / 2f64.to_radians()).ceil() as usize).clamp(1, 45);
        (0..=2 * k)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let a = base + self.half_angle * (i as f64 / k as f64 - 1.0);
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect()
    }

    /// Signed relative margin of a measured energy per steradian: the smaller of
    /// `(E − min)/min` (when `min > 0`) and `(max − E)/max`. Negative: outside.
    pub fn margin(&self, e: f64) -> f64 {
        let (lo, hi) = self.energy;
        let mut m = (hi - e) / hi.abs().max(1e-300);
        if lo > 0.0 {
            m = m.min((e - lo) / lo);
        }
        m
    }

    /// Highest frequency the measurement needs resolved (for the sampling of the flight):
    /// the band's top, or 0 without a band.
    pub fn top_frequency(&self) -> f64 {
        self.band.map_or(0.0, |(_, hi)| hi)
    }

    /// Largest rate `d(t − n·r/c)/dt = 1 − n·v/c` over the arc's directions at velocity
    /// `v` (the smallest `n·v` is at an end of the arc, or on the axis for a point).
    pub fn phase_rate(&self, v: DVec3, c: f64) -> f64 {
        let base = self.axis.y.atan2(self.axis.x);
        [-self.half_angle, self.half_angle]
            .iter()
            .map(|d| 1.0 - DVec3::new((base + d).cos(), (base + d).sin(), 0.0).dot(v) / c)
            .fold(0.0, f64::max)
    }

    /// Measured energy per steradian, averaged over the arc (Simpson's rule in angle).
    pub fn measure(&self, samples: &[Emission], q: f64, c: f64) -> f64 {
        let dirs = self.directions();
        let values: Vec<f64> = dirs
            .iter()
            .map(|&n| match self.band {
                None => lienard_energy(samples, q, c, n),
                Some((lo, hi)) => band_energy(samples, q, c, n, lo, hi),
            })
            .collect();
        if values.len() == 1 {
            return values[0];
        }
        let m = values.len() - 1;
        let sum: f64 = values
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let w = if i == 0 || i == m {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * v
            })
            .sum();
        #[allow(clippy::cast_precision_loss)]
        let mean = sum / (3.0 * m as f64);
        mean
    }
}

/// `n × ((n − β) × β̇)` for velocity `v` and acceleration `a`.
fn radiation_vector(n: DVec3, v: DVec3, a: DVec3, c: f64) -> DVec3 {
    n.cross((n - v / c).cross(a / c))
}

/// Energy radiated per steradian in direction `n` over the samples (Liénard; Simpson's
/// rule on consecutive pairs of intervals, in its form for unequal widths, and the
/// trapezoidal rule on a last odd interval).
pub fn lienard_energy(samples: &[Emission], q: f64, c: f64, n: DVec3) -> f64 {
    if !c.is_finite() || samples.len() < 2 {
        return 0.0;
    }
    let f: Vec<f64> = samples
        .iter()
        .map(|&(_, _, v, a)| {
            let kappa = 1.0 - n.dot(v) / c;
            radiation_vector(n, v, a, c).length_squared() / kappa.powi(5)
        })
        .collect();
    let t: Vec<f64> = samples.iter().map(|s| s.0).collect();
    let mut sum = 0.0;
    let mut i = 0;
    while i + 2 < t.len() {
        let (h0, h1) = (t[i + 1] - t[i], t[i + 2] - t[i + 1]);
        if h0 > 0.0 && h1 > 0.0 {
            sum += (h0 + h1) / 6.0
                * ((2.0 - h1 / h0) * f[i]
                    + (h0 + h1) * (h0 + h1) / (h0 * h1) * f[i + 1]
                    + (2.0 - h0 / h1) * f[i + 2]);
        }
        i += 2;
    }
    if i + 1 < t.len() {
        sum += 0.5 * (f[i] + f[i + 1]) * (t[i + 1] - t[i]);
    }
    q * q / (4.0 * PI * c) * sum
}

/// Minimal complex number for the radiation integral.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    fn cis(phase: f64) -> Self {
        let (s, c) = phase.sin_cos();
        Self { re: c, im: s }
    }

    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
}

/// `φ₀(D) = ∫₀¹ e^{iDs} ds` and `φ₁(D) = ∫₀¹ s e^{iDs} ds`, given `e^{iD}`.
fn filon_weights(d: f64, e: Cx) -> (Cx, Cx) {
    if d.abs() < 0.05 {
        let d2 = d * d;
        let phi0 = Cx {
            re: 1.0 - d2 / 6.0 + d2 * d2 / 120.0,
            im: d / 2.0 - d * d2 / 24.0,
        };
        let phi1 = Cx {
            re: 0.5 - d2 / 8.0 + d2 * d2 / 144.0,
            im: d / 3.0 - d * d2 / 30.0,
        };
        return (phi0, phi1);
    }
    // φ₀ = (e − 1)/(iD); φ₁ = e/(iD) + (e − 1)/D².
    let phi0 = Cx {
        re: e.im / d,
        im: (1.0 - e.re) / d,
    };
    let phi1 = Cx {
        re: e.im / d + (e.re - 1.0) / (d * d),
        im: -e.re / d + e.im / (d * d),
    };
    (phi0, phi1)
}

/// `d²I/dω dΩ` in direction `n` at the frequencies `ω₀ + k δω` (`k < count`).
pub fn spectrum(
    samples: &[Emission],
    q: f64,
    c: f64,
    n: DVec3,
    omega0: f64,
    d_omega: f64,
    count: usize,
) -> Vec<f64> {
    if !c.is_finite() || samples.len() < 2 || count == 0 {
        return vec![0.0; count];
    }
    // Amplitude f = n × ((n − β) × β̇)/κ² and retarded phase time τ = t − n·r/c.
    let pts: Vec<(f64, f64, DVec3)> = samples
        .iter()
        .map(|&(t, x, v, a)| {
            let kappa = 1.0 - n.dot(v) / c;
            (
                t,
                t - n.dot(x) / c,
                radiation_vector(n, v, a, c) / (kappa * kappa),
            )
        })
        .collect();
    let mut re = vec![DVec3::ZERO; count];
    let mut im = vec![DVec3::ZERO; count];
    for w in pts.windows(2) {
        let ((ta, tau_a, fa), (tb, tau_b, fb)) = (w[0], w[1]);
        let h = tb - ta;
        let dtau = tau_b - tau_a;
        let df = fb - fa;
        // e^{iωτ_a} and e^{iωΔτ} along the frequency grid, by recurrence.
        let mut ea = Cx::cis(omega0 * tau_a);
        let step_a = Cx::cis(d_omega * tau_a);
        let mut ed = Cx::cis(omega0 * dtau);
        let step_d = Cx::cis(d_omega * dtau);
        #[allow(clippy::cast_precision_loss)]
        for k in 0..count {
            let d = (omega0 + k as f64 * d_omega) * dtau;
            let (p0, p1) = filon_weights(d, ed);
            // h e^{iωτ_a} [f_a φ₀ + (f_b − f_a) φ₁]
            let w0 = ea.mul(p0);
            let w1 = ea.mul(p1);
            re[k] += (fa * w0.re + df * w1.re) * h;
            im[k] += (fa * w0.im + df * w1.im) * h;
            ea = ea.mul(step_a);
            ed = ed.mul(step_d);
        }
    }
    let scale = q * q / (4.0 * PI * PI * c);
    re.iter()
        .zip(&im)
        .map(|(r, i)| scale * (r.length_squared() + i.length_squared()))
        .collect()
}

/// Frequency spacing that resolves the spectrum of a flight whose retarded phase time
/// spans `span`: a line from a finite flight is `2π/span` wide; 4 points across it (the
/// trapezoidal rule on a smooth, band-limited `|A(ω)|²`; tests S1–S3 hold as with 8).
pub fn resolving_step(span: f64) -> f64 {
    PI / (2.0 * span.max(1e-9))
}

/// Energy per steradian in direction `n` within the band `[lo, hi]`: the spectrum on a
/// grid that resolves the flight's lines (`resolving_step`, at most 4097 points),
/// integrated by the trapezoidal rule.
pub fn band_energy(samples: &[Emission], q: f64, c: f64, n: DVec3, lo: f64, hi: f64) -> f64 {
    if !c.is_finite() || samples.len() < 2 || hi <= lo {
        return 0.0;
    }
    let tau = |s: &Emission| s.0 - n.dot(s.1) / c;
    let span = (tau(&samples[samples.len() - 1]) - tau(&samples[0])).abs();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let intervals = (((hi - lo) / resolving_step(span)).ceil() as usize).clamp(8, 4096);
    #[allow(clippy::cast_precision_loss)]
    let d_omega = (hi - lo) / intervals as f64;
    let s = spectrum(samples, q, c, n, lo, d_omega, intervals + 1);
    let inner: f64 = s[1..intervals].iter().sum();
    (inner + 0.5 * (s[0] + s[intervals])) * d_omega
}

/// A frequency range for displaying a flight's spectrum without a band: three times the
/// largest critical frequency `(3/2) γ³ c |a⊥| / v²` along it (Jackson 14.85, with the
/// bend radius `v²/|a⊥|`), where the spectrum of a bend has fallen off.
pub fn display_range(samples: &[Emission], c: f64) -> f64 {
    samples
        .iter()
        .map(|&(_, _, v, a)| {
            let v2 = v.length_squared();
            if v2 <= 0.0 || !c.is_finite() {
                return 0.0;
            }
            let gamma = 1.0 / (1.0 - v2 / (c * c)).max(1e-300).sqrt();
            let a_perp = (a - v * (a.dot(v) / v2)).length();
            4.5 * gamma.powi(3) * c * a_perp / v2
        })
        .fold(0.0, f64::max)
}

/// Spectrum for display: `d²I/dω dΩ` averaged over (up to 9 of) the window's directions,
/// from 0 to `omega_max`, computed on a grid that resolves the flight's lines and averaged
/// into `bins` bins (so a narrow line keeps its energy). Returns the bins' centres and
/// mean values.
pub fn arc_spectrum(
    window: &RadiationWindow,
    samples: &[Emission],
    q: f64,
    c: f64,
    omega_max: f64,
    bins: usize,
) -> Vec<(f64, f64)> {
    if !c.is_finite() || samples.len() < 2 || bins == 0 || omega_max <= 0.0 {
        return Vec::new();
    }
    let dirs = window.directions();
    let stride = dirs.len().div_ceil(9);
    let chosen: Vec<DVec3> = dirs.iter().step_by(stride).copied().collect();
    let tau_span = |n: DVec3| {
        let tau = |s: &Emission| s.0 - n.dot(s.1) / c;
        (tau(&samples[samples.len() - 1]) - tau(&samples[0])).abs()
    };
    let span = chosen.iter().map(|&n| tau_span(n)).fold(0.0, f64::max);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let per_bin = ((omega_max / bins as f64 / resolving_step(span)).ceil() as usize).clamp(1, 64);
    let count = bins * per_bin;
    #[allow(clippy::cast_precision_loss)]
    let d_omega = omega_max / count as f64;
    let mut total = vec![0.0; count];
    for &n in &chosen {
        // Grid points at the fine cells' centres.
        let s = spectrum(samples, q, c, n, 0.5 * d_omega, d_omega, count);
        for (t, v) in total.iter_mut().zip(s) {
            *t += v;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let norm = (chosen.len() * per_bin) as f64;
    (0..bins)
        .map(|b| {
            let sum: f64 = total[b * per_bin..(b + 1) * per_bin].iter().sum();
            #[allow(clippy::cast_precision_loss)]
            let centre = (b as f64 + 0.5) * omega_max / bins as f64;
            (centre, sum / norm)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filon_weights_match_their_series() {
        for d in [0.049_999, 0.050_001, 0.3, -0.7] {
            let e = Cx::cis(d);
            let (a0, a1) = filon_weights(d, e);
            // Direct quadrature of ∫₀¹ e^{iDs} and ∫₀¹ s e^{iDs}.
            let m = 20_000;
            let (mut b0, mut b1) = (Cx::default(), Cx::default());
            for j in 0..m {
                let s = (f64::from(j) + 0.5) / f64::from(m);
                let z = Cx::cis(d * s);
                b0.re += z.re / f64::from(m);
                b0.im += z.im / f64::from(m);
                b1.re += s * z.re / f64::from(m);
                b1.im += s * z.im / f64::from(m);
            }
            for (x, y) in [(a0, b0), (a1, b1)] {
                assert!(
                    (x.re - y.re).abs() < 1e-8 && (x.im - y.im).abs() < 1e-8,
                    "{d}"
                );
            }
        }
    }

    #[test]
    fn arc_directions_are_symmetric_and_fine() {
        let w = RadiationWindow {
            axis: DVec3::Y,
            half_angle: 10f64.to_radians(),
            band: None,
            energy: (0.0, 1.0),
        };
        let d = w.directions();
        assert_eq!(d.len(), 11);
        assert!((d[5] - DVec3::Y).length() < 1e-12);
        assert!((d[0].x + d[10].x).abs() < 1e-12);
    }
}
