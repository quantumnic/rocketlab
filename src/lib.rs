#![allow(
    clippy::manual_strip,
    clippy::manual_pattern_char_comparison,
    clippy::collapsible_else_if,
    clippy::if_same_then_else,
    clippy::type_complexity,
    clippy::assign_op_pattern,
    clippy::option_if_let_else,
    clippy::too_many_arguments,
    clippy::new_without_default,
    clippy::wrong_self_convention,
    clippy::should_implement_trait,
    clippy::result_unit_err,
    clippy::manual_unwrap_or,
    dead_code
)]

pub mod atmosphere;
pub mod attitude;
pub mod constants;
pub mod constellation;
pub mod cr3bp;

pub mod delta_v_budget;
pub mod engines;
pub mod flyby;
pub mod heat_shield;
pub mod heating;
pub mod interplanetary;
pub mod kepler;
pub mod lambert;
pub mod landing;
pub mod launch_sites;
pub mod low_thrust;
pub mod mission;
pub mod nbody;
pub mod nozzle;
pub mod orbit_determination;
pub mod orbit_lifetime;
pub mod orbits;
pub mod parachute;
pub mod perturbations;
pub mod plotting;
pub mod propulsion;
pub mod reentry;
pub mod rendezvous;
pub mod sgp4;
pub mod soi;
pub mod trajectory;

#[cfg(feature = "wasm")]
pub mod wasm;
