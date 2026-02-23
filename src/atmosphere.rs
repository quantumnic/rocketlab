//! US Standard Atmosphere 1976 model.
//!
//! Computes temperature, pressure, and density as a function of geometric altitude.
//! Valid from 0 to 86 km (mesopause). Uses the seven-layer piecewise model with
//! linear temperature lapse rates.
//!
//! Reference: NASA-TM-X-74335, "U.S. Standard Atmosphere, 1976"

use crate::constants::{G0, M_AIR, R_EARTH, R_GAS};

/// Atmospheric layer boundaries and lapse rates for US Standard Atmosphere 1976.
/// Each tuple: (base_altitude_m, base_temperature_K, lapse_rate_K_per_m)
const LAYERS: [(f64, f64, f64); 7] = [
    (0.0, 288.15, -0.0065),      // Troposphere
    (11_000.0, 216.65, 0.0),     // Tropopause
    (20_000.0, 216.65, 0.001),   // Stratosphere lower
    (32_000.0, 228.65, 0.0028),  // Stratosphere upper
    (47_000.0, 270.65, 0.0),     // Stratopause
    (51_000.0, 270.65, -0.0028), // Mesosphere lower
    (71_000.0, 214.65, -0.002),  // Mesosphere upper → 86 km
];

/// Base pressures at each layer boundary (Pa), precomputed from the 1976 standard.
const BASE_PRESSURES: [f64; 7] = [
    101_325.0,
    22_632.064_0,
    5_474.889_4,
    868.018_7,
    110.906_2,
    66.938_9,
    3.956_4,
];

/// Result of an atmospheric query.
#[derive(Debug, Clone, Copy)]
pub struct AtmosphereResult {
    /// Geometric altitude (m)
    pub altitude: f64,
    /// Temperature (K)
    pub temperature: f64,
    /// Pressure (Pa)
    pub pressure: f64,
    /// Density (kg/m³)
    pub density: f64,
    /// Speed of sound (m/s)
    pub speed_of_sound: f64,
    /// Dynamic viscosity (Pa·s) — Sutherland's law
    pub viscosity: f64,
}

/// Atmospheric model selector
#[derive(Debug, Clone, Copy)]
pub enum AtmosphereModel {
    /// US Standard Atmosphere 1976 (0-86 km)
    UsStandard1976,
    /// Mars atmosphere model (placeholder)
    Mars,
}

impl AtmosphereModel {
    /// Create US Standard 1976 atmosphere model
    pub fn us_standard_1976() -> Self {
        AtmosphereModel::UsStandard1976
    }

    /// Create Mars atmosphere model (placeholder)
    pub fn mars() -> Self {
        AtmosphereModel::Mars
    }

    /// Get atmospheric properties at given altitude
    pub fn query(&self, altitude_m: f64) -> Option<AtmosphereResult> {
        match self {
            AtmosphereModel::UsStandard1976 => us_standard_atmosphere(altitude_m),
            AtmosphereModel::Mars => {
                // Placeholder: use simplified Mars atmosphere
                // Real implementation would need proper Mars model
                if altitude_m > 80_000.0 { return None; }
                Some(AtmosphereResult {
                    altitude: altitude_m,
                    temperature: 210.0 - altitude_m * 0.001, // Rough approximation
                    pressure: 610.0 * (-altitude_m / 10_400.0).exp(), // Scale height ~10.4 km
                    density: 0.02 * (-altitude_m / 10_400.0).exp(),
                    speed_of_sound: 240.0, // Rough approximation
                    viscosity: 1.3e-5,     // Rough approximation
                })
            }
        }
    }

    /// Get atmospheric density at given altitude
    pub fn density(&self, altitude_m: f64) -> Result<f64, String> {
        match self.query(altitude_m) {
            Some(result) => Ok(result.density),
            None => Err("Altitude out of range".to_string()),
        }
    }

    /// Get atmospheric density, pressure, and temperature at given altitude
    pub fn density_pressure_temperature(&self, altitude_m: f64) -> Result<(f64, f64, f64), String> {
        match self.query(altitude_m) {
            Some(result) => Ok((result.density, result.pressure, result.temperature)),
            None => Err("Altitude out of range".to_string()),
        }
    }
}

/// Ratio of specific heats for air
const GAMMA_AIR: f64 = 1.4;

/// Sutherland reference viscosity (Pa·s)
const MU_REF: f64 = 1.716e-5;
/// Sutherland reference temperature (K)
const T_REF: f64 = 273.15;
/// Sutherland constant (K)
const S_SUTH: f64 = 110.4;

/// Query the US Standard Atmosphere 1976 at a given geometric altitude.
///
/// # Arguments
/// * `altitude_m` - Geometric altitude in meters (0 to 86,000)
///
/// # Returns
/// `AtmosphereResult` with temperature, pressure, density, speed of sound, viscosity.
///
/// Returns `None` if altitude is out of range.
pub fn us_standard_atmosphere(altitude_m: f64) -> Option<AtmosphereResult> {
    if !(0.0..=86_000.0).contains(&altitude_m) {
        return None;
    }

    // Convert geometric to geopotential altitude
    let h = geopotential_altitude(altitude_m);

    // Find the layer
    let layer_idx = find_layer(h);
    let (h_b, t_b, lapse) = LAYERS[layer_idx];
    let p_b = BASE_PRESSURES[layer_idx];

    let dh = h - h_b;
    let temperature = t_b + lapse * dh;

    let pressure = if lapse.abs() < 1e-10 {
        // Isothermal layer: P = P_b * exp(-g₀·M·Δh / (R·T_b))
        p_b * (-G0 * M_AIR * dh / (R_GAS * t_b)).exp()
    } else {
        // Gradient layer: P = P_b * (T/T_b)^(-g₀·M / (R·λ))
        let exponent = -G0 * M_AIR / (R_GAS * lapse);
        p_b * (temperature / t_b).powf(exponent)
    };

    let density = pressure * M_AIR / (R_GAS * temperature);
    let speed_of_sound = (GAMMA_AIR * R_GAS * temperature / M_AIR).sqrt();
    let viscosity = sutherland_viscosity(temperature);

    Some(AtmosphereResult {
        altitude: altitude_m,
        temperature,
        pressure,
        density,
        speed_of_sound,
        viscosity,
    })
}

/// Convert geometric altitude to geopotential altitude.
///
/// h_gp = R_E · h_geo / (R_E + h_geo)
fn geopotential_altitude(geometric_m: f64) -> f64 {
    let r_e = R_EARTH * 1000.0; // km to m
    r_e * geometric_m / (r_e + geometric_m)
}

/// Find which atmospheric layer a geopotential altitude belongs to.
fn find_layer(h_geopotential: f64) -> usize {
    for i in (0..LAYERS.len()).rev() {
        if h_geopotential >= LAYERS[i].0 {
            return i;
        }
    }
    0
}

/// Sutherland's law for dynamic viscosity.
fn sutherland_viscosity(t: f64) -> f64 {
    MU_REF * (t / T_REF).powf(1.5) * (T_REF + S_SUTH) / (t + S_SUTH)
}

/// Mach number from velocity and altitude.
pub fn mach_number(velocity_m_s: f64, altitude_m: f64) -> Option<f64> {
    us_standard_atmosphere(altitude_m).map(|atm| velocity_m_s / atm.speed_of_sound)
}

/// Dynamic pressure q = 0.5 · ρ · v²  (Pa)
pub fn dynamic_pressure(velocity_m_s: f64, density: f64) -> f64 {
    0.5 * density * velocity_m_s * velocity_m_s
}

/// Calculate Max-Q for a trajectory given as (time, altitude_m, velocity_m_s) points.
///
/// Returns (time_at_max_q, altitude_at_max_q, max_q_pa).
pub fn find_max_q(trajectory: &[(f64, f64, f64)]) -> Option<(f64, f64, f64)> {
    let mut max_q = 0.0_f64;
    let mut result = None;

    for &(t, alt, vel) in trajectory {
        if let Some(atm) = us_standard_atmosphere(alt) {
            let q = dynamic_pressure(vel, atm.density);
            if q > max_q {
                max_q = q;
                result = Some((t, alt, q));
            }
        }
    }

    result
}

/// Drag force magnitude (N).
///
/// F_drag = 0.5 · ρ · v² · C_d · A
pub fn drag_force(velocity_m_s: f64, density: f64, cd: f64, area_m2: f64) -> f64 {
    dynamic_pressure(velocity_m_s, density) * cd * area_m2
}

/// Sutton-Graves stagnation-point heating rate (W/m²).
///
/// q̇ = k · sqrt(ρ/r_n) · v³
///
/// where k ≈ 1.7415e-4 for Earth (air), r_n is nose radius (m).
pub fn sutton_graves_heating(velocity_m_s: f64, density: f64, nose_radius_m: f64) -> f64 {
    let k = 1.7415e-4;
    k * (density / nose_radius_m).sqrt() * velocity_m_s.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sea_level() {
        let atm = us_standard_atmosphere(0.0).unwrap();
        assert!(
            (atm.temperature - 288.15).abs() < 0.1,
            "T: {}",
            atm.temperature
        );
        assert!(
            (atm.pressure - 101_325.0).abs() < 1.0,
            "P: {}",
            atm.pressure
        );
        assert!((atm.density - 1.225).abs() < 0.01, "ρ: {}", atm.density);
        assert!(
            (atm.speed_of_sound - 340.3).abs() < 0.5,
            "a: {}",
            atm.speed_of_sound
        );
    }

    #[test]
    fn test_tropopause() {
        // At 11 km: T ≈ 216.65 K, P ≈ 22632 Pa
        let atm = us_standard_atmosphere(11_000.0).unwrap();
        assert!(
            (atm.temperature - 216.65).abs() < 1.0,
            "T: {}",
            atm.temperature
        );
        assert!(
            (atm.pressure - 22_632.0).abs() < 100.0,
            "P: {}",
            atm.pressure
        );
    }

    #[test]
    fn test_stratosphere() {
        // At 20 km: T ≈ 216.65 K, P ≈ 5474.9 Pa, ρ ≈ 0.0880 kg/m³
        let atm = us_standard_atmosphere(20_000.0).unwrap();
        assert!(
            (atm.temperature - 216.65).abs() < 1.0,
            "T: {}",
            atm.temperature
        );
        assert!(
            (atm.pressure - 5_474.9).abs() < 100.0,
            "P: {}",
            atm.pressure
        );
        assert!((atm.density - 0.0880).abs() < 0.005, "ρ: {}", atm.density);
    }

    #[test]
    fn test_high_altitude() {
        // At 50 km: T ≈ 270.65 K, P ≈ 79.78 Pa
        let atm = us_standard_atmosphere(50_000.0).unwrap();
        assert!(
            (atm.temperature - 270.65).abs() < 2.0,
            "T: {}",
            atm.temperature
        );
        assert!((atm.pressure - 79.78).abs() < 5.0, "P: {}", atm.pressure);
    }

    #[test]
    fn test_out_of_range() {
        assert!(us_standard_atmosphere(-1.0).is_none());
        assert!(us_standard_atmosphere(87_000.0).is_none());
    }

    #[test]
    fn test_dynamic_pressure() {
        // At sea level, 100 m/s: q = 0.5 * 1.225 * 10000 = 6125 Pa
        let q = dynamic_pressure(100.0, 1.225);
        assert!((q - 6125.0).abs() < 1.0);
    }

    #[test]
    fn test_mach_sea_level() {
        let m = mach_number(340.3, 0.0).unwrap();
        assert!((m - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_max_q() {
        let traj = vec![
            (0.0, 0.0, 0.0),
            (30.0, 5000.0, 300.0),
            (60.0, 15000.0, 500.0),
            (90.0, 40000.0, 1500.0),
            (120.0, 80000.0, 3000.0),
        ];
        let (t, alt, q) = find_max_q(&traj).unwrap();
        // Max-Q should be at the point with best balance of density and velocity
        assert!(q > 0.0);
        assert!(t > 0.0);
        assert!(alt > 0.0);
    }
}
