//! Physics core: field sources, relativistic particle dynamics, event detection and
//! outcome verification. Pure library with no graphics or I/O dependencies.
//!
//! The implemented physics is documented in `PHYSICS.md` at the repository root.

pub mod antenna;
pub mod beam;
pub mod bem;
pub mod cancel;
pub mod circuit;
pub mod conductor;
pub mod drive;
pub mod dynamics;
pub mod events;
pub mod external;
pub mod field;
pub mod field_motion;
pub mod geometry;
pub mod inductance;
pub mod integrator;
pub mod lienard;
pub mod magnetic;
pub mod panel;
pub mod poynting;
pub mod spectrum;
pub mod trajectory;
pub mod tube;
pub mod units;
pub mod verify;
pub mod zinv;

pub use glam::DVec3;
