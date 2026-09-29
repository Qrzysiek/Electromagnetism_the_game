//! Field energy and its flow (Jackson §6.7), in the game's units (Coulomb's constant 1,
//! μ₀/4π = 1/c²): the energy density `u = (E² + c²B²)/8π` and the Poynting vector
//! `S = (c²/4π) E × B`, which satisfy Poynting's theorem `∂u/∂t + ∇·S = −J·E`. For the
//! energy-flow view (PHYSICS.md §10) and its tests (P1–P3). Finite c only: at c = ∞ the
//! magnetic energy of a finite B is infinite in these units.
//!
//! A field made of two parts (a particle's and the rest) has, besides each part's own terms,
//! the exchange terms `u₁₂ = (E₁·E₂ + c² B₁·B₂)/4π` and `S₁₂ = (c²/4π)(E₁ × B₂ + E₂ × B₁)`,
//! with `∂u₁₂/∂t + ∇·S₁₂ = −(J₁·E₂ + J₂·E₁)`: they carry the work that each part does on
//! the other's charges.

use std::f64::consts::PI;

use glam::DVec3;

/// Energy density `(E² + c²B²)/8π`.
pub fn energy_density(e: DVec3, b: DVec3, c: f64) -> f64 {
    (e.length_squared() + c * c * b.length_squared()) / (8.0 * PI)
}

/// Poynting vector `(c²/4π) E × B`: the flow of field energy.
pub fn poynting(e: DVec3, b: DVec3, c: f64) -> DVec3 {
    e.cross(b) * (c * c / (4.0 * PI))
}

/// Exchange terms of two fields `(E₁, B₁)` and `(E₂, B₂)`: the energy density
/// `(E₁·E₂ + c² B₁·B₂)/4π` and its flow `(c²/4π)(E₁ × B₂ + E₂ × B₁)`.
pub fn exchange(e1: DVec3, b1: DVec3, e2: DVec3, b2: DVec3, c: f64) -> (f64, DVec3) {
    (
        (e1.dot(e2) + c * c * b1.dot(b2)) / (4.0 * PI),
        (e1.cross(b2) + e2.cross(b1)) * (c * c / (4.0 * PI)),
    )
}

/// The velocity of the field energy, `S/u`. Never faster than light (`|S| ≤ c u`, since
/// `2|E| c|B| ≤ E² + c²B²`), and exactly c in a radiation field (`|E| = c|B|`, `E ⊥ B`).
/// Zero where there is no field.
pub fn energy_velocity(e: DVec3, b: DVec3, c: f64) -> DVec3 {
    let u = energy_density(e, b, c);
    if u > 0.0 {
        poynting(e, b, c) / u
    } else {
        DVec3::ZERO
    }
}
