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

    /// Exact gradient of `B_z`, from `B_z = 3 (M·r) r_z / r⁵ − M_z / r³`:
    /// `∇B_z = 3 (M r_z + (M·r) ẑ) / r⁵ − 15 (M·r) r_z r / r⁷ + 3 M_z r / r⁵`.
    pub fn grad_bz(&self, x: DVec3) -> DVec3 {
        let r = x - self.position;
        let r2 = r.length_squared();
        let inv_r5 = 1.0 / (r2 * r2 * r2.sqrt());
        let m = self.moment;
        let mr = m.dot(r);
        (m * r.z + DVec3::Z * mr) * (3.0 * inv_r5) - r * (15.0 * mr * r.z * inv_r5 / r2)
            + r * (3.0 * m.z * inv_r5)
    }
}

/// Circular coil (thin wire of radius `wire_radius`) of radius `radius` around `center`
/// in the plane orthogonal to `normal`, carrying current counter-clockwise about
/// `normal`, with strength `kappa = μ₀ I / 4π` at t = 0 and `kappa + rate·t` later (a
/// ramped current; its induced field is `−rate · A_unit`, quasi-static, PHYSICS.md §2.2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularLoop {
    pub center: DVec3,
    pub normal: DVec3,
    pub radius: f64,
    pub kappa: f64,
    pub wire_radius: f64,
    /// Ramp rate dκ/dt (0: a steady current).
    pub rate: f64,
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

/// `h(m) = (1 − m/2) K(m) − E(m)`, the bracket of the loop's vector potential, which is
/// `O(m²)` (`π m²/32` to leading order): summed as a power series for small `m`, where the
/// direct form cancels. The coefficient of `mⁿ` is `π/2 [c_n 2n/(2n − 1) − c_(n−1)/2]`
/// (from the series of K and E, see `loop_bracket`); it vanishes for n = 1.
fn potential_bracket(m: f64, m1: f64) -> f64 {
    if m < 0.05 {
        let mut c_prev = 1.0f64; // c_0
        let mut sum = 0.0;
        let mut mp = 1.0;
        for n in 1..=16u32 {
            let nf = f64::from(n);
            let c = c_prev * ((2.0 * nf - 1.0) / (2.0 * nf)).powi(2);
            mp *= m;
            if n >= 2 {
                sum += mp * (c * 2.0 * nf / (2.0 * nf - 1.0) - 0.5 * c_prev);
            }
            c_prev = c;
        }
        std::f64::consts::FRAC_PI_2 * sum
    } else {
        let (k, e) = elliptic_ke(m, m1);
        (1.0 - 0.5 * m) * k - e
    }
}

impl CircularLoop {
    /// Strength at time `t`.
    pub fn kappa_at(&self, t: f64) -> f64 {
        self.kappa + self.rate * t
    }

    /// Field at time `t` (for a steady current exactly `field`).
    pub fn field_at(&self, x: DVec3, t: f64) -> DVec3 {
        if self.rate == 0.0 {
            return self.field(x);
        }
        Self {
            kappa: self.kappa_at(t),
            ..*self
        }
        .field(x)
    }

    /// Vector potential per unit strength (κ = 1), azimuthal: with `m = k² = 4aρ/β²`,
    /// `A_φ = (4/k) √(a/ρ) h(m) = 2β h(m)/ρ`, `h = (1 − m/2) K − E` (e.g. Jackson §5.5,
    /// eq. 5.37, with μ₀I/4π = κ). Its curl is `field` (test M8).
    pub fn unit_vector_potential(&self, x: DVec3) -> DVec3 {
        let n = self.normal;
        let d = x - self.center;
        let z = d.dot(n);
        let radial = d - n * z;
        let rho = radial.length();
        if rho == 0.0 {
            return DVec3::ZERO;
        }
        let a = self.radius;
        let alpha2 = (a - rho) * (a - rho) + z * z;
        let beta2 = (a + rho) * (a + rho) + z * z;
        let m = 4.0 * a * rho / beta2;
        let a_phi = 2.0 * beta2.sqrt() * potential_bracket(m, alpha2 / beta2) / rho;
        n.cross(radial / rho) * a_phi
    }

    /// Induced electric field of the ramp, `−∂A/∂t = −rate · A_unit` (quasi-static).
    pub fn induced_e(&self, x: DVec3) -> DVec3 {
        if self.rate == 0.0 {
            DVec3::ZERO
        } else {
            self.unit_vector_potential(x) * (-self.rate)
        }
    }

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

    /// Gradient of the component of B along `normal`, at a point in the loop's plane
    /// (outside the wire). There B is curl-free, so `∂B_z/∂ρ = ∂B_ρ/∂z`; with
    /// `B_ρ = 2κ z β g(m) / (α² ρ)` and α, β, m even in z, at z = 0 this is
    /// `2κ β g(m) / (α² ρ)`, along ρ̂ (zero on the axis). `g` is the series-evaluated
    /// bracket of `field`, so the result is accurate near the axis too.
    pub fn grad_bz_in_plane(&self, x: DVec3) -> DVec3 {
        let n = self.normal;
        let d = x - self.center;
        let radial = d - n * d.dot(n);
        let rho = radial.length();
        if rho == 0.0 {
            return DVec3::ZERO;
        }
        let a = self.radius;
        let alpha2 = (a - rho) * (a - rho);
        let beta2 = (a + rho) * (a + rho);
        let beta = beta2.sqrt();
        let m = 4.0 * a * rho / beta2;
        let dbz_drho = 2.0 * self.kappa * beta / (alpha2 * rho) * loop_bracket(m, alpha2 / beta2);
        radial / rho * dbz_drho
    }
}

/// Closed polygonal coil through `vertices` (the last connects back to the first),
/// current flowing in vertex order, strength `kappa = μ₀ I / 4π` at t = 0 and
/// `kappa + rate·t` later (see `CircularLoop`).
#[derive(Clone, Debug, PartialEq)]
pub struct PolygonCoil {
    pub vertices: Vec<DVec3>,
    pub kappa: f64,
    pub wire_radius: f64,
    /// Ramp rate dκ/dt (0: a steady current).
    pub rate: f64,
}

/// Vector potential of a straight segment per unit strength, `∫ dl / |x − l|` along it:
/// `û ln((|rb| + rb·û) / (|ra| + ra·û))` with `ra = a − x`, `rb = b − x`. Behind the
/// segment's start both terms of that ratio are tiny; the equal form
/// `(|ra| − ra·û) / (|rb| − rb·û)` (both products are |r⊥|²) is used there.
fn segment_potential(ra: DVec3, rb: DVec3) -> DVec3 {
    let d = rb - ra;
    let u = d / d.length();
    let (la, lb) = (ra.length(), rb.length());
    let (pa, pb) = (ra.dot(u), rb.dot(u));
    let ratio = if pa + pb >= 0.0 {
        (lb + pb) / (la + pa)
    } else {
        (la - pa) / (lb - pb)
    };
    u * ratio.ln()
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
    /// Strength at time `t`.
    pub fn kappa_at(&self, t: f64) -> f64 {
        self.kappa + self.rate * t
    }

    /// Field at time `t` (for a steady current exactly `field`).
    pub fn field_at(&self, x: DVec3, t: f64) -> DVec3 {
        if self.rate == 0.0 {
            return self.field(x);
        }
        let n = self.vertices.len();
        let mut b = DVec3::ZERO;
        let k = self.kappa_at(t);
        for i in 0..n {
            b += segment_field(self.vertices[i] - x, self.vertices[(i + 1) % n] - x, k);
        }
        b
    }

    /// Vector potential per unit strength (κ = 1): the sum of the segments'.
    pub fn unit_vector_potential(&self, x: DVec3) -> DVec3 {
        let n = self.vertices.len();
        let mut a = DVec3::ZERO;
        for i in 0..n {
            a += segment_potential(self.vertices[i] - x, self.vertices[(i + 1) % n] - x);
        }
        a
    }

    /// Induced electric field of the ramp, `−rate · A_unit` (quasi-static).
    pub fn induced_e(&self, x: DVec3) -> DVec3 {
        if self.rate == 0.0 {
            DVec3::ZERO
        } else {
            self.unit_vector_potential(x) * (-self.rate)
        }
    }

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

    /// Gradient of `B_z` at a point in the coil's plane (z = const, outside the wires).
    /// The closed coil's field is curl-free there, so `∇B_z = ∂B_in-plane/∂z`. For one
    /// segment at height h above the point, `ra × rb = u × v − h (u − v) × ẑ`, and the
    /// scalar factor of `segment_field` is even in h, so `∂B/∂z = −F (a − b) × ẑ` with
    /// `F = κ (|u| + |v|) / (|u| |v| (|u| |v| + u·v))`. (Per segment this is not the
    /// gradient of that segment's B_z, whose field alone is not curl-free; summed over
    /// the closed polygon it is.)
    pub fn grad_bz_in_plane(&self, x: DVec3) -> DVec3 {
        let n = self.vertices.len();
        let mut g = DVec3::ZERO;
        for i in 0..n {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % n];
            let (u, v) = (a - x, b - x);
            let (lu, lv) = (u.length(), v.length());
            let f = self.kappa * (lu + lv) / (lu * lv * (lu * lv + u.dot(v)));
            g -= (a - b).cross(DVec3::Z) * f;
        }
        g
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
