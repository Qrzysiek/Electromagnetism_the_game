//! Dimensionless internal units and their conversion to SI (PHYSICS.md §1).
//!
//! A world is fixed by a reference length `L₀`, charge `Q₀` and mass `M₀`. Internally the
//! Coulomb constant `k = 1/(4πε₀)` is 1, which fixes the time unit
//! `T₀ = sqrt(M₀ L₀³ / (k Q₀²))`. The speed of light becomes the dimensionless number
//! `c = c_SI T₀ / L₀`.

/// Speed of light in vacuum, m/s (exact).
pub const C_SI: f64 = 299_792_458.0;
/// Coulomb constant `1/(4πε₀)`, N·m²/C² (CODATA 2018).
pub const K_SI: f64 = 8.987_551_792_3e9;
/// Elementary charge, C (exact).
pub const E_CHARGE_SI: f64 = 1.602_176_634e-19;
/// Electron mass, kg (CODATA 2018).
pub const M_ELECTRON_SI: f64 = 9.109_383_701_5e-31;
/// Proton mass, kg (CODATA 2018).
pub const M_PROTON_SI: f64 = 1.672_621_923_69e-27;

/// SI reference values defining a world's unit system.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorldScale {
    /// `L₀`, metres (grid cell size).
    pub length_m: f64,
    /// `Q₀`, coulombs.
    pub charge_c: f64,
    /// `M₀`, kilograms.
    pub mass_kg: f64,
}

impl WorldScale {
    /// `T₀`, seconds.
    pub fn time_s(&self) -> f64 {
        (self.mass_kg * self.length_m.powi(3) / (K_SI * self.charge_c * self.charge_c)).sqrt()
    }

    /// `E₀ = k Q₀² / L₀`, joules.
    pub fn energy_j(&self) -> f64 {
        K_SI * self.charge_c * self.charge_c / self.length_m
    }

    /// `V₀ = L₀ / T₀`, m/s.
    pub fn velocity_m_s(&self) -> f64 {
        self.length_m / self.time_s()
    }

    /// `p₀ = M₀ V₀`, kg·m/s.
    pub fn momentum_si(&self) -> f64 {
        self.mass_kg * self.velocity_m_s()
    }

    /// The speed of light in internal units.
    pub fn c(&self) -> f64 {
        C_SI / self.velocity_m_s()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_units_are_consistent() {
        let s = WorldScale {
            length_m: 1e-6,
            charge_c: 1e3 * E_CHARGE_SI,
            mass_kg: M_PROTON_SI,
        };
        // E₀ = M₀ V₀² follows from k = 1 in internal units.
        let e_from_kinematics = s.mass_kg * s.velocity_m_s().powi(2);
        assert!((s.energy_j() / e_from_kinematics - 1.0).abs() < 1e-14);
        assert!((s.c() * s.velocity_m_s() / C_SI - 1.0).abs() < 1e-15);
    }
}
