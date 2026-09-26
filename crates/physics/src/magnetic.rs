//! Static magnetic field sources (PHYSICS.md §2.2).
//!
//! Source strengths are given directly in field units, i.e. already multiplied by
//! `μ₀/4π`: a dipole by `μ = μ₀ m / 4π` and a coil by `κ = μ₀ I / 4π`. In internal units
//! (`k = 1/(4πε₀) = 1`) `μ₀/4π = 1/c²`, so this keeps magnets meaningful in Newtonian
//! levels (`c = ∞`) while remaining exactly the SI relations. The unit of B is
//! `M₀ / (Q₀ T₀)`.

use glam::DVec3;

/// Uniformly magnetized sphere. Outside the sphere its field is exactly that of a point
/// dipole: `B = μ (3 (m̂·r̂) r̂ − m̂) / r³`, with `moment` = `μ m̂`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagneticDipole {
    pub position: DVec3,
    pub moment: DVec3,
    pub radius: f64,
}

impl MagneticDipole {
    pub fn field(&self, x: DVec3) -> DVec3 {
        let r = x - self.position;
        let r2 = r.length_squared();
        let inv_r = 1.0 / r2.sqrt();
        let inv_r3 = inv_r / r2;
        let r_hat = r * inv_r;
        (r_hat * (3.0 * self.moment.dot(r_hat)) - self.moment) * inv_r3
    }
}

/// Circular coil (thin wire of radius `wire_radius`) of radius `radius` around `center`
/// in the plane orthogonal to `normal`, carrying current counter-clockwise about
/// `normal`, with strength `kappa = μ₀ I / 4π`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularLoop {
    pub center: DVec3,
    pub normal: DVec3,
    pub radius: f64,
    pub kappa: f64,
    pub wire_radius: f64,
}

/// Complete elliptic integrals `(K(m), E(m))` of parameter `m = k²`, given `m` and
/// `m1 = 1 − m` (passed separately to avoid cancellation near `m = 1`). Computed with the
/// arithmetic–geometric mean, which converges quadratically to full precision
/// (DLMF 19.8.1 and 19.8.6).
pub fn elliptic_ke(m: f64, m1: f64) -> (f64, f64) {
    let mut a = 1.0f64;
    let mut b = m1.sqrt();
    let mut c = m.sqrt();
    let mut sum = 0.5 * c * c; // 2^(n-1) c_n^2 for n = 0
    let mut pow2 = 0.5;
    for _ in 0..40 {
        if (c.abs()) <= f64::EPSILON * a {
            break;
        }
        let an = 0.5 * (a + b);
        let bn = (a * b).sqrt();
        c = 0.5 * (a - b);
        a = an;
        b = bn;
        pow2 *= 2.0;
        sum += pow2 * c * c;
    }
    let k = std::f64::consts::PI / (2.0 * a);
    (k, k * (1.0 - sum))
}

/// `g(m) = (1 − m/2) E(m) − (1 − m) K(m)`, which is `O(m²)`: evaluated directly it
/// cancels catastrophically near the loop axis, so for small `m` it is summed as a power
/// series. With `c_n = (C(2n, n) / 4ⁿ)²`, `K = π/2 Σ c_n mⁿ` and
/// `E = π/2 Σ c_n mⁿ / (1 − 2n)` (DLMF 19.5.1, 19.5.2); the coefficients of `m⁰` and `m¹`
/// in `g` vanish identically.
fn loop_bracket(m: f64, m1: f64) -> f64 {
    if m < 0.05 {
        let mut c = 1.0f64; // c_0
        let mut e_prev = 1.0; // E coefficient of m^(n-1)
        let mut k_prev = 1.0; // K coefficient of m^(n-1)
        let mut sum = 0.0;
        let mut mp = 1.0;
        for n in 1..=16u32 {
            let nf = f64::from(n);
            // c_n = c_(n-1) ((2n − 1) / (2n))²
            c *= ((2.0 * nf - 1.0) / (2.0 * nf)).powi(2);
            let e_n = c / (1.0 - 2.0 * nf);
            let k_n = c;
            mp *= m;
            if n >= 2 {
                sum += mp * (e_n - 0.5 * e_prev - k_n + k_prev);
            }
            e_prev = e_n;
            k_prev = k_n;
        }
        std::f64::consts::FRAC_PI_2 * sum
    } else {
        let (k, e) = elliptic_ke(m, m1);
        (1.0 - 0.5 * m) * e - m1 * k
    }
}

impl CircularLoop {
    /// Exact field of a circular current loop (e.g. Simpson et al., NASA/TM-2001-209961):
    /// with `α² = (a−ρ)² + z²`, `β² = (a+ρ)² + z²`, `k² = 1 − α²/β²`,
    /// `B_ρ = 2κ z / (α² β ρ) [(a² + ρ² + z²) E − α² K]`,
    /// `B_z = 2κ / (α² β) [(a² − ρ² − z²) E + α² K]`.
    pub fn field(&self, x: DVec3) -> DVec3 {
        let n = self.normal;
        let d = x - self.center;
        let z = d.dot(n);
        let radial = d - n * z;
        let rho = radial.length();
        let a = self.radius;
        let alpha2 = (a - rho) * (a - rho) + z * z;
        let beta2 = (a + rho) * (a + rho) + z * z;
        let beta = beta2.sqrt();
        let m1 = alpha2 / beta2;
        let (k, e) = elliptic_ke(4.0 * a * rho / beta2, m1);
        let bz =
            2.0 * self.kappa / (alpha2 * beta) * ((a * a - rho * rho - z * z) * e + alpha2 * k);
        // B_ρ vanishes on the axis and (exactly, since z = 0) in the plane of the loop.
        // With a² + ρ² + z² = β²(1 − m/2) and α² = β² m1, the bracket
        // [(a² + ρ² + z²) E − α² K] equals β² g(m).
        let b_rho = if rho > 0.0 && z != 0.0 {
            let m = 4.0 * a * rho / beta2;
            2.0 * self.kappa * z * beta / (alpha2 * rho) * loop_bracket(m, m1)
        } else {
            0.0
        };
        let rho_hat = if rho > 0.0 { radial / rho } else { DVec3::ZERO };
        n * bz + rho_hat * b_rho
    }
}

/// Closed polygonal coil through `vertices` (the last connects back to the first),
/// current flowing in vertex order, strength `kappa = μ₀ I / 4π`.
#[derive(Clone, Debug, PartialEq)]
pub struct PolygonCoil {
    pub vertices: Vec<DVec3>,
    pub kappa: f64,
    pub wire_radius: f64,
}

/// Exact field of a straight segment from `a` to `b` (current from `a` to `b`) at the
/// origin-relative vectors `ra = a − x`, `rb = b − x`:
/// `B = κ (ra × rb)(|ra| + |rb|) / (|ra| |rb| (|ra| |rb| + ra·rb))`.
fn segment_field(ra: DVec3, rb: DVec3, kappa: f64) -> DVec3 {
    let la = ra.length();
    let lb = rb.length();
    let denom = la * lb * (la * lb + ra.dot(rb));
    ra.cross(rb) * (kappa * (la + lb) / denom)
}

impl PolygonCoil {
    pub fn field(&self, x: DVec3) -> DVec3 {
        let n = self.vertices.len();
        let mut b = DVec3::ZERO;
        for i in 0..n {
            let a = self.vertices[i];
            let c = self.vertices[(i + 1) % n];
            b += segment_field(a - x, c - x, self.kappa);
        }
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loop_bracket_series_matches_direct_formula() {
        // Leading term (π/2)(3/16) m², and continuity where the method switches.
        let m = 1e-6;
        let lead = std::f64::consts::FRAC_PI_2 * 3.0 / 16.0 * m * m;
        assert!((loop_bracket(m, 1.0 - m) / lead - 1.0).abs() < 1e-5);
        for m in [0.02, 0.049_999, 0.05] {
            let (k, e) = elliptic_ke(m, 1.0 - m);
            let direct = (1.0 - 0.5 * m) * e - (1.0 - m) * k;
            let series = if m < 0.05 {
                loop_bracket(m, 1.0 - m)
            } else {
                direct
            };
            // The direct formula cancels (g ≈ 1e-3 while K, E ≈ 1.6), so it is only good
            // to a few ulp of K and E in absolute terms.
            assert!(
                (series - direct).abs() <= 1e-14,
                "m = {m}: {series} vs {direct}"
            );
        }
    }

    #[test]
    fn elliptic_integrals_match_reference_values() {
        // K(0) = E(0) = π/2; K(1/2) and E(1/2) from DLMF / A&S Table 17.1.
        let (k0, e0) = elliptic_ke(0.0, 1.0);
        assert!((k0 - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
        assert!((e0 - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
        let (k, e) = elliptic_ke(0.5, 0.5);
        assert!((k - 1.854_074_677_301_371_9).abs() < 1e-15, "K = {k}");
        assert!((e - 1.350_643_881_047_675_5).abs() < 1e-15, "E = {e}");
        // Legendre's relation at m = 1/2: 2 K E − K² = π/2.
        assert!((2.0 * k * e - k * k - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
    }
}
