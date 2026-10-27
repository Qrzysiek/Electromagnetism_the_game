//! Shared helpers for the validation tests.

#![allow(dead_code)]
#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)

use physics::DVec3;
use physics::dynamics::Particle;
use physics::field::{Coulomb, FixedCharge};
use physics::geometry::{Aabb, Shape, Sphere};
use physics::trajectory::{Scenario, StepView};

pub const UNIT_PARTICLE: Particle = Particle {
    charge: 1.0,
    mass: 1.0,
    radius: 0.0,
    moment: 0.0,
};

/// One fixed charge `q` at the origin (radius `r`), unit test particle.
pub fn central_charge(
    q: f64,
    r: f64,
    c: f64,
    x0: DVec3,
    p0: DVec3,
    t_max: f64,
) -> Scenario<Coulomb> {
    let charge = FixedCharge {
        position: DVec3::ZERO,
        charge: q,
        radius: r,
    };
    Scenario {
        field: Coulomb::new(&[charge]),
        obstacles: vec![Shape::Sphere(Sphere {
            center: DVec3::ZERO,
            radius: r,
        })],
        particle: UNIT_PARTICLE,
        c,
        x0,
        p0,
        detector: None,
        bounds: None,
        t_max,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    }
}

pub fn cube(half: f64) -> Aabb {
    Aabb {
        min: DVec3::splat(-half),
        max: DVec3::splat(half),
    }
}

/// Angle between two vectors in [0, π], accurate for small and large angles.
pub fn angle_between(a: DVec3, b: DVec3) -> f64 {
    a.cross(b).length().atan2(a.dot(b))
}

/// Gauss–Legendre nodes and weights on [−1, 1] (Newton iteration on P_n).
pub fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            // P_n(x) and P_n'(x) by the three-term recurrence.
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                #[allow(clippy::cast_precision_loss)]
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            dp = nf * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

/// Quadrature over the unit sphere, polar axis along `axis`: Gauss–Legendre in cos θ
/// (`n_mu` nodes) times the trapezoid rule in φ (`n_phi`). Directions and weights (their
/// sum is 4π); exponentially convergent for integrands analytic on the sphere.
pub fn sphere_quadrature(axis: DVec3, n_mu: usize, n_phi: usize) -> Vec<(DVec3, f64)> {
    let e3 = axis.normalize();
    let e1 = e3.any_orthonormal_vector();
    let e2 = e3.cross(e1);
    #[allow(clippy::cast_precision_loss)]
    let dphi = std::f64::consts::TAU / n_phi as f64;
    let mut out = Vec::with_capacity(n_mu * n_phi);
    for (mu, w) in gauss_legendre(n_mu) {
        let s = (1.0 - mu * mu).sqrt();
        for k in 0..n_phi {
            #[allow(clippy::cast_precision_loss)]
            let (sp, cp) = (k as f64 * dphi).sin_cos();
            out.push((e1 * (s * cp) + e2 * (s * sp) + e3 * mu, w * dphi));
        }
    }
    out
}

/// Deterministic pseudo-random numbers (SplitMix64) for reproducible configurations.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        #[allow(clippy::cast_precision_loss)]
        let r = (z >> 11) as f64 / (1u64 << 53) as f64;
        r
    }

    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// Analytic relativistic Coulomb orbit (PHYSICS.md §9): `u(φ) = A cos(Γφ) + B₀`, with
/// `φ = 0` at the periapsis, for potential energy `κ/r`.
#[derive(Clone, Copy, Debug)]
pub struct CoulombOrbit {
    pub gamma: f64,
    pub a: f64,
    pub b0: f64,
}

impl CoulombOrbit {
    /// From the conserved quantities of the initial state (`c` may be infinite).
    pub fn from_state(kappa: f64, m: f64, c: f64, x: DVec3, p: DVec3) -> Self {
        let l = x.cross(p).length();
        if c.is_infinite() {
            // Newtonian limit: Γ = 1, B₀ = −κm/L², A = sqrt(2mE/L² + B₀²).
            let e = p.length_squared() / (2.0 * m) + kappa / x.length();
            let b0 = -kappa * m / (l * l);
            return Self {
                gamma: 1.0,
                a: (2.0 * m * e / (l * l) + b0 * b0).sqrt(),
                b0,
            };
        }
        let gamma_lorentz = (1.0 + p.length_squared() / (m * m * c * c)).sqrt();
        let w = gamma_lorentz * m * c * c + kappa / x.length();
        let lc2 = l * l * c * c;
        let g2 = 1.0 - kappa * kappa / lc2;
        let b0 = -kappa * w / (lc2 * g2);
        // W² − m²c⁴ = (W − mc²)(W + mc²), written to avoid cancellation.
        let mc2 = m * c * c;
        let a = ((w - mc2) * (w + mc2) / (lc2 * g2) + b0 * b0).sqrt();
        Self {
            gamma: g2.sqrt(),
            a,
            b0,
        }
    }

    /// Scattering deflection angle `|π − (2/Γ) arccos(−B₀/A)|`.
    pub fn deflection(&self) -> f64 {
        (std::f64::consts::PI - 2.0 / self.gamma * (-self.b0 / self.a).acos()).abs()
    }

    /// Periapsis distance `1/(A + B₀)`.
    pub fn r_min(&self) -> f64 {
        1.0 / (self.a + self.b0)
    }

    /// Radius at polar angle `φ` measured from the periapsis.
    pub fn r(&self, phi: f64) -> f64 {
        1.0 / (self.a * (self.gamma * phi).cos() + self.b0)
    }
}

/// Times at which the radial velocity `x·p` changes sign from negative to positive
/// (periapsis passages) within one step, located by bisection on the dense output.
pub fn periapsis_times<F: physics::field::FieldSolver>(step: &StepView<'_, F>) -> Vec<f64> {
    let s = |t: f64| {
        let (x, p) = step.state(t);
        x.dot(p)
    };
    let (ta, tb) = (step.t_start(), step.t_end());
    let n = 16;
    let mut out = Vec::new();
    let mut t_prev = ta;
    let mut s_prev = s(ta);
    for i in 1..=n {
        let t = ta + (tb - ta) * f64::from(i) / f64::from(n);
        let st = s(t);
        if s_prev < 0.0 && st >= 0.0 {
            let (mut lo, mut hi) = (t_prev, t);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    break;
                }
                if s(mid) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            out.push(hi);
        }
        t_prev = t;
        s_prev = st;
    }
    out
}

/// A measured value that a document quotes (DEC-22): written to `target/measured/<id>` to
/// `digits` significant digits, from where `python scripts/measured.py update` puts it
/// between the document's markers `<!--m:<id>-->…<!--/m-->`.
pub fn measured(id: &str, value: f64, digits: usize) {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/measured");
    if std::fs::create_dir_all(&dir).is_ok() {
        let _ = std::fs::write(dir.join(id), significant(value, digits));
    }
}

/// `value` to `digits` significant digits: fixed-point from 1e-3 to 1e4 (`0.027`, `1.071`),
/// scientific beyond (`1.2e-8`).
pub fn significant(value: f64, digits: usize) -> String {
    let sci = format!("{:.*e}", digits.saturating_sub(1), value);
    // The exponent after rounding (9.99 to two digits is 1.0e1).
    let exp: i32 = sci
        .split('e')
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap_or(0);
    if (-3..4).contains(&exp) {
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            clippy::cast_sign_loss
        )]
        let decimals = (digits as i32 - 1 - exp).max(0) as usize;
        format!("{value:.decimals$}")
    } else {
        sci
    }
}

// The Drude metal's and the self-field's tests (drude.rs, darwin.rs): their metal, the
// small floating electrode that anchors their electrostatics, the Boltzmann equation's AC
// mobility (their reference for σ(ω)), and complex helpers whose square root and Bessel
// series stay on `libm` (num-complex's own square root goes through the platform's
// library).

pub use num_complex::Complex64;
use physics::bem::Bias;
use physics::lattice::Material;
use physics::zinv::{Electrode, Electrodes, Section};
use std::f64::consts::PI;

/// The mean collision angle's cosine (specular hard discs).
pub const COS_THETA: f64 = -1.0 / 3.0;

/// The tests' metal: background 1 (the carriers' charge density −1), carriers of `q/m` −1
/// and thermal speed 1, cores of radius `core_radius`.
#[must_use]
pub fn metal(core_radius: f64) -> Material {
    Material {
        background: 1.0,
        core_radius,
        charge_per_mass: -1.0,
        thermal_speed: 1.0,
    }
}

/// A small floating square electrode (side 0.5, no net charge) at `center`: the reference the
/// z-invariant electrostatics needs (docs/OPEN_PROBLEMS.md, OP-1), far from the metal.
#[must_use]
pub fn small_floating_electrode(center: DVec3) -> Electrodes {
    Electrodes::new(
        vec![Electrode {
            section: Section::Rect {
                center,
                angle: 0.0,
                half_length: 0.25,
                half_thickness: 0.25,
            },
            bias: Bias::Charge(0.0),
        }],
        &[],
        0.25,
    )
}

/// The Boltzmann equation's mean velocity per unit acceleration (along x) of carriers of
/// thermal speed 1 colliding at `κ v`, rotating at `−ω` (PHYSICS.md §2.11; equally, the
/// response to a drive `∝ e^{iωt}`, as V-D6 uses it): `(v_x + i v_y)/a = (π/κ) P/(1 − ⟨cos θ⟩ Q)`,
/// `P = κ ∫ v³ f₀/(κ v + i ω) dv` with `f₀ = e^{−v²/2}/2π`, `Q = κ ∫ v g/(κ v + i ω) dv` with
/// the drawn speeds' law `g = √(2/π) v² e^{−v²/2}`. At `ω = 0`: `(3/8) √(π/2)/κ`. Simpson's
/// rule on [0, 14]. With charge density 1 and |q/m| 1 it is the conductivity σ(ω).
#[must_use]
pub fn mobility(kappa: f64, omega: f64) -> Complex64 {
    let (n, top) = (28_000, 14.0);
    let h = top / f64::from(n);
    let (mut p, mut q) = (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0));
    for k in 0..=n {
        let v = h * f64::from(k);
        let w = if k == 0 || k == n {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        // κ/(κv + iω) = κ (κv − iω)/(κ²v² + ω²).
        let d = kappa * kappa * v * v + omega * omega;
        if d == 0.0 {
            continue;
        }
        let (re, im) = (kappa * kappa * v / d, -kappa * omega / d);
        let e = libm::exp(-0.5 * v * v);
        let fp = v * v * v * e / (2.0 * PI);
        let fq = v * libm::sqrt(2.0 / PI) * v * v * e;
        p += Complex64::new(w * fp * re, w * fp * im);
        q += Complex64::new(w * fq * re, w * fq * im);
    }
    let (p, q) = (p * h / 3.0, q * h / 3.0);
    // (π/κ) p / (1 − c q).
    let den = Complex64::new(1.0 - COS_THETA * q.re, -COS_THETA * q.im);
    p / den * (PI / kappa)
}

/// The principal square root, in `libm`.
#[must_use]
pub fn csqrt(a: Complex64) -> Complex64 {
    let r = libm::hypot(a.re, a.im);
    let re = libm::sqrt(0.5 * (r + a.re));
    let im = libm::sqrt(0.5 * (r - a.re)).copysign(a.im);
    Complex64::new(re, im)
}

/// The modified Bessel function `I₀(z) = Σ (z²/4)^m/(m!)²` of a complex `z` (|z| ≲ 5).
#[must_use]
pub fn ci0(z: Complex64) -> Complex64 {
    let q = z * z * 0.25;
    let (mut sum, mut term) = (Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0));
    for m in 1..40 {
        let m = f64::from(m);
        term = term * q / (m * m);
        sum += term;
    }
    sum
}
