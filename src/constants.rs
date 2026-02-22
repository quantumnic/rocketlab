//! Physical and astronomical constants used throughout rocketlab.
//!
//! Sources: IAU 2012, WGS-84, NASA Fact Sheets

/// Gravitational parameter of Earth (km³/s²)
pub const MU_EARTH: f64 = 398_600.441_8;

/// Gravitational parameter of the Sun (km³/s²)
pub const MU_SUN: f64 = 1.327_124_400_18e11;

/// Gravitational parameter of Mars (km³/s²)
pub const MU_MARS: f64 = 42_828.375_214;

/// Gravitational parameter of the Moon (km³/s²)
pub const MU_MOON: f64 = 4_902.800_066;

/// Mean radius of Earth (km)
pub const R_EARTH: f64 = 6_371.0;

/// Equatorial radius of Earth (km) — WGS-84
pub const R_EARTH_EQUATORIAL: f64 = 6_378.137;

/// Mean radius of Mars (km)
pub const R_MARS: f64 = 3_389.5;

/// Standard gravitational acceleration at Earth's surface (m/s²)
pub const G0: f64 = 9.80665;

/// Speed of light in vacuum (km/s)
pub const C_LIGHT: f64 = 299_792.458;

/// Earth's J2 oblateness coefficient — WGS-84
pub const J2_EARTH: f64 = 1.082_63e-3;

/// Earth's angular rotation rate (rad/s)
pub const OMEGA_EARTH: f64 = 7.2921159e-5;

/// Boltzmann constant (J/K)
pub const K_BOLTZMANN: f64 = 1.380649e-23;

/// Universal gas constant (J/(mol·K))
pub const R_GAS: f64 = 8.314_462_618;

/// Standard atmospheric pressure at sea level (Pa)
pub const P0_ATM: f64 = 101_325.0;

/// Standard temperature at sea level (K)
pub const T0_ATM: f64 = 288.15;

/// Mean molecular weight of air (kg/mol)
pub const M_AIR: f64 = 0.0289644;

/// Astronomical unit (km)
pub const AU_KM: f64 = 1.495_978_707e8;

/// Pi constant (for clarity in formulas)
pub const PI: f64 = std::f64::consts::PI;

/// Two pi
pub const TWO_PI: f64 = 2.0 * PI;
