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
//! The integral of 14.65 is taken in the phase time `τ = t − n·r/c` (`dτ = κ dt`), where
//! the phase `ωτ` is exactly linear: with the amplitude `n × ((n − β) × β̇)/κ³` quadratic
//! in τ over each pair of pieces between dense samples of the flight, each pair is
//! integrated exactly against `e^{iωτ}` (Filon's rule), so the samples need to resolve the
//! motion but not the phase, however high the frequency.
//!
//! Several particles (a *system*: `system_*`): their far fields add. The field received
//! at the time `T` after the light time `R/c` from the origin comes from each particle at
//! its own phase time `τ = T − R/c`, so the amplitudes add at equal `τ` (and the spectral
//! amplitudes at equal ω) before squaring (Jackson Pr. 14.23):
//!
//! ```text
//! dW/dΩ      = (1/4πc) ∫ |Σⱼ qⱼ gⱼ(τ)|² dτ,   gⱼ = n × ((n − βⱼ) × β̇ⱼ)/κⱼ³
//! d²I/dω dΩ  = (1/4π²c) |Σⱼ qⱼ ∫ gⱼ e^{iωτ} dτ|²
//! ```

use glam::DVec3;
use std::f64::consts::PI;

/// A sample of the flight: time, position, velocity, acceleration.
pub type Emission = (f64, DVec3, DVec3, DVec3);

/// A particle of a system: its charge and the samples of its flight.
pub type Source<'a> = (f64, &'a [Emission]);

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
    /// The detector stops the particle abruptly (a target): the stop's radiation counts
    /// (`spectrum`). Its spectrum is flat to infinite frequency, so this needs a band.
    pub abrupt_stop: bool,
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

    /// The band's top, or 0 without a band.
    pub fn top_frequency(&self) -> f64 {
        self.band.map_or(0.0, |(_, hi)| hi)
    }

    /// Measured energy per steradian, averaged over the arc (Simpson's rule in angle).
    pub fn measure(&self, samples: &[Emission], q: f64, c: f64) -> f64 {
        let dirs = self.directions();
        let values: Vec<f64> = dirs
            .iter()
            .map(|&n| match self.band {
                None => lienard_energy(samples, q, c, n),
                Some((lo, hi)) => band_energy(samples, q, c, n, lo, hi, self.abrupt_stop),
            })
            .collect();
        arc_mean(&values)
    }

    /// The same measure for a system of particles, whose far fields add (coherently). An
    /// abrupt stop is not supported here (it is ignored).
    pub fn measure_system(&self, sources: &[Source<'_>], c: f64) -> f64 {
        let dirs = self.directions();
        let values: Vec<f64> = dirs
            .iter()
            .map(|&n| match self.band {
                None => system_lienard_energy(sources, c, n),
                Some((lo, hi)) => system_band_energy(sources, c, n, lo, hi),
            })
            .collect();
        arc_mean(&values)
    }
}

/// The mean over the arc of values at `directions()`: Simpson's rule in angle.
fn arc_mean(values: &[f64]) -> f64 {
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

/// The moments `φₖ(D) = ∫₀¹ uᵏ e^{iDu} du`, k = 0, 1, 2, given `e^{iD}`: by their series
/// `Σⱼ (iD)ʲ/(j! (j + k + 1))` for |D| < 0.25 (ten terms, error < 1e-14), otherwise in
/// closed form, `φ₀ = (e − 1)/(iD)`, `φₖ = e/(iD) − k φₖ₋₁/(iD)`.
fn filon_moments(d: f64, e: Cx) -> [Cx; 3] {
    if d.abs() < 0.25 {
        let mut out = [Cx::default(); 3];
        for (k, phi) in out.iter_mut().enumerate() {
            // term_j = (iD)^j / j!
            let mut term = Cx { re: 1.0, im: 0.0 };
            for j in 0..10 {
                #[allow(clippy::cast_precision_loss)]
                let w = 1.0 / (j + k + 1) as f64;
                phi.re += term.re * w;
                phi.im += term.im * w;
                #[allow(clippy::cast_precision_loss)]
                let next = term.mul(Cx {
                    re: 0.0,
                    im: d / (j + 1) as f64,
                });
                term = next;
            }
        }
        return out;
    }
    // x/(iD) = (x.im/D, −x.re/D)
    let over_id = |x: Cx| Cx {
        re: x.im / d,
        im: -x.re / d,
    };
    let phi0 = over_id(Cx {
        re: e.re - 1.0,
        im: e.im,
    });
    let e_id = over_id(e);
    let p0 = over_id(phi0);
    let phi1 = Cx {
        re: e_id.re - p0.re,
        im: e_id.im - p0.im,
    };
    let p1 = over_id(phi1);
    let phi2 = Cx {
        re: e_id.re - 2.0 * p1.re,
        im: e_id.im - 2.0 * p1.im,
    };
    [phi0, phi1, phi2]
}

/// `d²I/dω dΩ` in direction `n` at the frequencies `ω₀ + k δω` (`k < count`). With
/// `abrupt_stop` the particle stops instantly at the last sample: the integrand is the
/// derivative of `F = n × (n × β)/κ` (Jackson, before 14.66), so the jump of `F` to 0
/// adds `−F_end e^{iω(t − n·r/c)}` at the stop to the integral (the sudden-stop
/// radiation of Jackson §15.2, flat in frequency). Valid for frequencies far below the
/// inverse of the real stopping time.
#[allow(clippy::too_many_arguments)]
pub fn spectrum(
    samples: &[Emission],
    q: f64,
    c: f64,
    n: DVec3,
    omega0: f64,
    d_omega: f64,
    count: usize,
    abrupt_stop: bool,
) -> Vec<f64> {
    if !c.is_finite() || samples.len() < 2 || count == 0 {
        return vec![0.0; count];
    }
    let (re, im) = amplitude(samples, c, n, omega0, d_omega, count, abrupt_stop);
    let scale = q * q / (4.0 * PI * PI * c);
    re.iter()
        .zip(&im)
        .map(|(r, i)| scale * (r.length_squared() + i.length_squared()))
        .collect()
}

/// `d²I/dω dΩ` of a system in direction `n` at the frequencies `ω₀ + k δω`: the particles'
/// amplitudes `qⱼ ∫ gⱼ e^{iωτ} dτ` added before squaring.
pub fn system_spectrum(
    sources: &[Source<'_>],
    c: f64,
    n: DVec3,
    omega0: f64,
    d_omega: f64,
    count: usize,
) -> Vec<f64> {
    let mut re = vec![DVec3::ZERO; count];
    let mut im = vec![DVec3::ZERO; count];
    if !c.is_finite() || count == 0 {
        return vec![0.0; count];
    }
    for &(q, samples) in sources {
        if q == 0.0 || samples.len() < 2 {
            continue;
        }
        let (r, i) = amplitude(samples, c, n, omega0, d_omega, count, false);
        for k in 0..count {
            re[k] += r[k] * q;
            im[k] += i[k] * q;
        }
    }
    let scale = 1.0 / (4.0 * PI * PI * c);
    re.iter()
        .zip(&im)
        .map(|(r, i)| scale * (r.length_squared() + i.length_squared()))
        .collect()
}

/// The phase times and amplitudes `(τ, g)` of the samples in direction `n`:
/// `τ = t − n·r/c`, `g = n × ((n − β) × β̇)/κ³`.
fn phase_points(samples: &[Emission], c: f64, n: DVec3) -> Vec<(f64, DVec3)> {
    samples
        .iter()
        .map(|&(t, x, v, a)| {
            let kappa = 1.0 - n.dot(v) / c;
            (
                t - n.dot(x) / c,
                radiation_vector(n, v, a, c) / (kappa * kappa * kappa),
            )
        })
        .collect()
}

/// A piece of the amplitude's model in the phase time (`pieces`).
type Piece = (f64, f64, [DVec3; 3]);

/// The pieces of the amplitude's model in the phase time: `(τ₀, H, [c₀, c₁, c₂])` with
/// `g = c₀ + c₁ u + c₂ u²`, `u = (τ − τ₀)/H`, quadratic over each pair of pieces between
/// samples (Lagrange basis on the pair's unequal widths), linear on a last odd piece (and
/// on a degenerate pair).
fn pieces(pts: &[(f64, DVec3)]) -> Vec<Piece> {
    let mut out = Vec::new();
    let mut i = 0;
    while i + 2 < pts.len() {
        let ((t0, g0), (t1, g1), (t2, g2)) = (pts[i], pts[i + 1], pts[i + 2]);
        let h = t2 - t0;
        let u = (t1 - t0) / h;
        if h > 0.0 && u > 0.0 && u < 1.0 {
            // Lagrange basis on u = 0, u₁, 1, as powers of u.
            let c0 = g0;
            let c1 = -g0 * ((u + 1.0) / u) - g1 / (u * (u - 1.0)) - g2 * (u / (1.0 - u));
            let c2 = g0 / u + g1 / (u * (u - 1.0)) + g2 / (1.0 - u);
            out.push((t0, h, [c0, c1, c2]));
        } else {
            out.push((t0, t1 - t0, [g0, g1 - g0, DVec3::ZERO]));
            out.push((t1, t2 - t1, [g1, g2 - g1, DVec3::ZERO]));
        }
        i += 2;
    }
    if i + 1 < pts.len() {
        let ((t0, g0), (t1, g1)) = (pts[i], pts[i + 1]);
        out.push((t0, t1 - t0, [g0, g1 - g0, DVec3::ZERO]));
    }
    out
}

/// The spectral amplitude `∫ g e^{iωτ} dτ` (real and imaginary parts) at the frequencies
/// `ω₀ + k δω` (`spectrum`).
fn amplitude(
    samples: &[Emission],
    c: f64,
    n: DVec3,
    omega0: f64,
    d_omega: f64,
    count: usize,
    abrupt_stop: bool,
) -> (Vec<DVec3>, Vec<DVec3>) {
    // In the phase time τ = t − n·r/c (dτ = κ dt, κ > 0) the phase is exactly linear:
    // ∫ f e^{iωτ} dt = ∫ g e^{iωτ} dτ with g = n × ((n − β) × β̇)/κ³. g is taken quadratic
    // in τ over each pair of pieces (Filon's rule, on the pieces' unequal widths) and
    // integrated exactly against e^{iωτ}; a last odd piece is taken linear.
    let pts = phase_points(samples, c, n);
    let mut re = vec![DVec3::ZERO; count];
    let mut im = vec![DVec3::ZERO; count];
    // Adds H e^{iωτ₀} Σₖ cₖ φₖ(ωH) over the frequency grid.
    let mut add = |tau0: f64, h: f64, coeffs: [DVec3; 3], order: usize| {
        if h <= 0.0 {
            return;
        }
        // e^{iωτ₀} and e^{iωH} along the frequency grid, by recurrence.
        let mut ea = Cx::cis(omega0 * tau0);
        let step_a = Cx::cis(d_omega * tau0);
        let mut ed = Cx::cis(omega0 * h);
        let step_d = Cx::cis(d_omega * h);
        #[allow(clippy::cast_precision_loss)]
        for k in 0..count {
            let d = (omega0 + k as f64 * d_omega) * h;
            let phi = filon_moments(d, ed);
            for (c_k, p) in coeffs.iter().zip(phi).take(order + 1) {
                let w = ea.mul(p);
                re[k] += *c_k * (w.re * h);
                im[k] += *c_k * (w.im * h);
            }
            ea = ea.mul(step_a);
            ed = ed.mul(step_d);
        }
    };
    for (t0, h, coeffs) in pieces(&pts) {
        let order = if coeffs[2] == DVec3::ZERO { 1 } else { 2 };
        add(t0, h, coeffs, order);
    }
    if abrupt_stop {
        let &(t, x, v, _) = samples.last().expect("at least two samples");
        let kappa = 1.0 - n.dot(v) / c;
        let f_end = n.cross(n.cross(v / c)) / kappa;
        let tau = t - n.dot(x) / c;
        let mut e = Cx::cis(omega0 * tau);
        let step = Cx::cis(d_omega * tau);
        for k in 0..count {
            re[k] -= f_end * e.re;
            im[k] -= f_end * e.im;
            e = e.mul(step);
        }
    }
    (re, im)
}

/// Energy radiated per steradian in direction `n` by a system: `(1/4πc) ∫ |Σⱼ qⱼ gⱼ|² dτ`
/// over the phase time, each `gⱼ` in its piecewise quadratic model (`pieces`, as in the
/// spectrum; 0 outside its flight, where it moves uniformly). On every interval between
/// the pieces' ends of all particles the integrand is a polynomial of degree 4:
/// three-point Gauss–Legendre integrates it exactly.
pub fn system_lienard_energy(sources: &[Source<'_>], c: f64, n: DVec3) -> f64 {
    if !c.is_finite() {
        return 0.0;
    }
    let models: Vec<(f64, Vec<Piece>)> = sources
        .iter()
        .filter(|(q, s)| *q != 0.0 && s.len() >= 2)
        .map(|&(q, s)| (q, pieces(&phase_points(s, c, n))))
        .filter(|(_, p)| !p.is_empty())
        .collect();
    let mut cuts: Vec<f64> = models
        .iter()
        .flat_map(|(_, p)| p.iter().flat_map(|&(t0, h, _)| [t0, t0 + h]))
        .collect();
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    // Per particle, the index of the piece at the current interval (they only advance).
    let mut at = vec![0usize; models.len()];
    let r = 0.5 * (0.6f64).sqrt();
    let gauss = [
        (0.5 - r, 5.0 / 18.0),
        (0.5, 8.0 / 18.0),
        (0.5 + r, 5.0 / 18.0),
    ];
    let mut sum = 0.0;
    for w in cuts.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b <= a {
            continue;
        }
        let mid = 0.5 * (a + b);
        // The piece of each particle containing the interval, if any.
        let active: Vec<(f64, f64, f64, [DVec3; 3])> = models
            .iter()
            .zip(at.iter_mut())
            .filter_map(|((q, p), k)| {
                while *k < p.len() && p[*k].0 + p[*k].1 <= mid {
                    *k += 1;
                }
                let (t0, h, co) = *p.get(*k)?;
                (t0 <= mid).then_some((*q, t0, h, co))
            })
            .collect();
        if active.is_empty() {
            continue;
        }
        for (u, weight) in gauss {
            let tau = a + u * (b - a);
            let total: DVec3 = active
                .iter()
                .map(|&(q, t0, h, co)| {
                    let s = (tau - t0) / h;
                    (co[0] + co[1] * s + co[2] * (s * s)) * q
                })
                .sum();
            sum += weight * total.length_squared() * (b - a);
        }
    }
    sum / (4.0 * PI * c)
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
pub fn band_energy(
    samples: &[Emission],
    q: f64,
    c: f64,
    n: DVec3,
    lo: f64,
    hi: f64,
    abrupt_stop: bool,
) -> f64 {
    if !c.is_finite() || samples.len() < 2 || hi <= lo {
        return 0.0;
    }
    let tau = |s: &Emission| s.0 - n.dot(s.1) / c;
    let span = (tau(&samples[samples.len() - 1]) - tau(&samples[0])).abs();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let intervals = (((hi - lo) / resolving_step(span)).ceil() as usize).clamp(8, 4096);
    #[allow(clippy::cast_precision_loss)]
    let d_omega = (hi - lo) / intervals as f64;
    let s = spectrum(samples, q, c, n, lo, d_omega, intervals + 1, abrupt_stop);
    let inner: f64 = s[1..intervals].iter().sum();
    (inner + 0.5 * (s[0] + s[intervals])) * d_omega
}

/// The phase-time span of a system's flight in direction `n` (for the frequency grid).
fn system_span(sources: &[Source<'_>], c: f64, n: DVec3) -> f64 {
    let tau = |s: &Emission| s.0 - n.dot(s.1) / c;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (_, s) in sources {
        if let (Some(first), Some(last)) = (s.first(), s.last()) {
            lo = lo.min(tau(first));
            hi = hi.max(tau(last));
        }
    }
    if hi > lo { hi - lo } else { 0.0 }
}

/// Energy per steradian of a system in direction `n` within the band `[lo, hi]`
/// (`band_energy` for the system's spectrum).
pub fn system_band_energy(sources: &[Source<'_>], c: f64, n: DVec3, lo: f64, hi: f64) -> f64 {
    if !c.is_finite() || hi <= lo {
        return 0.0;
    }
    let span = system_span(sources, c, n);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let intervals = (((hi - lo) / resolving_step(span)).ceil() as usize).clamp(8, 4096);
    #[allow(clippy::cast_precision_loss)]
    let d_omega = (hi - lo) / intervals as f64;
    let s = system_spectrum(sources, c, n, lo, d_omega, intervals + 1);
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
    let tau_span = |n: DVec3| {
        let tau = |s: &Emission| s.0 - n.dot(s.1) / c;
        (tau(&samples[samples.len() - 1]) - tau(&samples[0])).abs()
    };
    arc_spectrum_of(window, tau_span, omega_max, bins, |n, w0, dw, count| {
        spectrum(samples, q, c, n, w0, dw, count, window.abrupt_stop)
    })
}

/// `arc_spectrum` for a system of particles (their far fields add).
pub fn system_arc_spectrum(
    window: &RadiationWindow,
    sources: &[Source<'_>],
    c: f64,
    omega_max: f64,
    bins: usize,
) -> Vec<(f64, f64)> {
    if !c.is_finite() || bins == 0 || omega_max <= 0.0 {
        return Vec::new();
    }
    arc_spectrum_of(
        window,
        |n| system_span(sources, c, n),
        omega_max,
        bins,
        |n, w0, dw, count| system_spectrum(sources, c, n, w0, dw, count),
    )
}

/// The display spectrum from a spectrum function of direction and frequency grid.
fn arc_spectrum_of(
    window: &RadiationWindow,
    tau_span: impl Fn(DVec3) -> f64,
    omega_max: f64,
    bins: usize,
    spectrum_at: impl Fn(DVec3, f64, f64, usize) -> Vec<f64>,
) -> Vec<(f64, f64)> {
    let dirs = window.directions();
    let stride = dirs.len().div_ceil(9);
    let chosen: Vec<DVec3> = dirs.iter().step_by(stride).copied().collect();
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
        let s = spectrum_at(n, 0.5 * d_omega, d_omega, count);
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
    fn filon_moments_match_quadrature() {
        for d in [0.249_999, 0.250_001, 0.03, 1.7, -3.2, 40.0] {
            let got = filon_moments(d, Cx::cis(d));
            // Midpoint quadrature of ∫₀¹ uᵏ e^{iDu} du (error ~ D²/(24 m²)).
            let m = 200_000;
            let mut want = [Cx::default(); 3];
            for j in 0..m {
                let u = (f64::from(j) + 0.5) / f64::from(m);
                let z = Cx::cis(d * u);
                for (k, w) in want.iter_mut().enumerate() {
                    let uk = u.powi(i32::try_from(k).unwrap());
                    w.re += uk * z.re / f64::from(m);
                    w.im += uk * z.im / f64::from(m);
                }
            }
            for k in 0..3 {
                let (a, b) = (got[k], want[k]);
                assert!(
                    (a.re - b.re).abs() < 1e-9 && (a.im - b.im).abs() < 1e-9,
                    "D = {d}, k = {k}: {a:?} vs {b:?}"
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
            abrupt_stop: false,
        };
        let d = w.directions();
        assert_eq!(d.len(), 11);
        assert!((d[5] - DVec3::Y).length() < 1e-12);
        assert!((d[0].x + d[10].x).abs() < 1e-12);
    }
}
