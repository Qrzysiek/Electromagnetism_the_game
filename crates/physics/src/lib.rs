//! Physics core: field sources, relativistic particle dynamics, event detection and
//! outcome verification. Pure library with no graphics or I/O dependencies.
//!
//! The implemented physics is documented in `PHYSICS.md` at the repository root.

pub mod antenna;
pub mod dynamics;
pub mod events;
pub mod external;
pub mod field;
pub mod geometry;
pub mod integrator;
pub mod magnetic;
pub mod trajectory;
pub mod units;
pub mod verify;

pub use glam::DVec3;
