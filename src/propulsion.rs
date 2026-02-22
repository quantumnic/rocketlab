//! Propulsion calculations: Tsiolkovsky equation, specific impulse, thrust.
//!
//! Reference: Sutton & Biblarz, "Rocket Propulsion Elements"

use crate::constants::G0;

/// Tsiolkovsky rocket equation: Δv = Isp · g₀ · ln(m₀/mf)
///
/// # Arguments
/// * `isp` - Specific impulse (seconds)
/// * `mass_ratio` - m₀/mf (initial mass / final mass)
///
/// # Returns
/// Delta-v (m/s)
pub fn tsiolkovsky_delta_v(isp: f64, mass_ratio: f64) -> f64 {
    isp * G0 * mass_ratio.ln()
}

/// Inverse Tsiolkovsky: given Δv and Isp, find required mass ratio.
///
/// mass_ratio = exp(Δv / (Isp · g₀))
pub fn required_mass_ratio(delta_v: f64, isp: f64) -> f64 {
    (delta_v / (isp * G0)).exp()
}

/// Effective exhaust velocity from specific impulse.
///
/// v_e = Isp · g₀ (m/s)
pub fn exhaust_velocity(isp: f64) -> f64 {
    isp * G0
}

/// Specific impulse from exhaust velocity.
///
/// Isp = v_e / g₀ (seconds)
pub fn isp_from_exhaust_velocity(v_e: f64) -> f64 {
    v_e / G0
}

/// Thrust from mass flow rate and exhaust velocity.
///
/// F = ṁ · v_e (N)
pub fn thrust(mass_flow_rate: f64, v_e: f64) -> f64 {
    mass_flow_rate * v_e
}

/// Mass flow rate from thrust and specific impulse.
///
/// ṁ = F / (Isp · g₀) (kg/s)
pub fn mass_flow_rate(thrust_n: f64, isp: f64) -> f64 {
    thrust_n / (isp * G0)
}

/// Propellant mass required for a given delta-v.
///
/// m_prop = m_payload · (exp(Δv/(Isp·g₀)) - 1)
pub fn propellant_mass(delta_v: f64, isp: f64, payload_mass: f64) -> f64 {
    payload_mass * (required_mass_ratio(delta_v, isp) - 1.0)
}

/// Burn time for a given delta-v at constant thrust.
///
/// t = (m₀ · v_e / F) · (1 - exp(-Δv/v_e))
pub fn burn_time(delta_v: f64, isp: f64, initial_mass: f64, thrust_n: f64) -> f64 {
    let v_e = exhaust_velocity(isp);
    (initial_mass * v_e / thrust_n) * (1.0 - (-delta_v / v_e).exp())
}

/// Common propellant combinations with typical Isp values (sea level / vacuum)
#[derive(Debug, Clone)]
pub struct PropellantCombo {
    pub name: &'static str,
    pub oxidizer: &'static str,
    pub fuel: &'static str,
    pub isp_sea_level: f64,
    pub isp_vacuum: f64,
    pub mixture_ratio: f64,
}

/// Database of common propellant combinations.
pub fn propellant_database() -> Vec<PropellantCombo> {
    vec![
        PropellantCombo {
            name: "LOX/RP-1",
            oxidizer: "LOX",
            fuel: "RP-1",
            isp_sea_level: 263.0,
            isp_vacuum: 311.0,
            mixture_ratio: 2.56,
        },
        PropellantCombo {
            name: "LOX/LH2",
            oxidizer: "LOX",
            fuel: "LH2",
            isp_sea_level: 366.0,
            isp_vacuum: 451.0,
            mixture_ratio: 5.0,
        },
        PropellantCombo {
            name: "N2O4/UDMH",
            oxidizer: "N2O4",
            fuel: "UDMH",
            isp_sea_level: 254.0,
            isp_vacuum: 311.0,
            mixture_ratio: 2.6,
        },
        PropellantCombo {
            name: "LOX/CH4",
            oxidizer: "LOX",
            fuel: "CH4 (Methane)",
            isp_sea_level: 299.0,
            isp_vacuum: 363.0,
            mixture_ratio: 3.6,
        },
        PropellantCombo {
            name: "Solid (APCP)",
            oxidizer: "AP (integral)",
            fuel: "HTPB (integral)",
            isp_sea_level: 242.0,
            isp_vacuum: 268.0,
            mixture_ratio: 0.0, // N/A for solid
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tsiolkovsky() {
        // Saturn V first stage: Isp=263s, mass ratio≈3.5 → Δv≈3232 m/s
        let dv = tsiolkovsky_delta_v(263.0, 3.5);
        assert!((dv - 3232.0).abs() < 50.0, "dv: {dv}");
    }

    #[test]
    fn test_mass_ratio_roundtrip() {
        let isp = 300.0;
        let mr = 5.0;
        let dv = tsiolkovsky_delta_v(isp, mr);
        let mr2 = required_mass_ratio(dv, isp);
        assert!((mr - mr2).abs() < 1e-10, "mr: {mr} vs {mr2}");
    }

    #[test]
    fn test_exhaust_velocity() {
        // Isp 300s → v_e ≈ 2942 m/s
        let v_e = exhaust_velocity(300.0);
        assert!((v_e - 2941.995).abs() < 0.1, "v_e: {v_e}");
    }

    #[test]
    fn test_thrust_calculation() {
        // F = mdot * v_e
        let v_e = exhaust_velocity(300.0);
        let mdot = 100.0; // kg/s
        let f = thrust(mdot, v_e);
        assert!((f - 294_199.5).abs() < 1.0, "F: {f}");
    }

    #[test]
    fn test_propellant_mass() {
        // 1000 m/s dv, 300s Isp, 1000 kg payload
        let m_prop = propellant_mass(1000.0, 300.0, 1000.0);
        let expected = 1000.0 * (required_mass_ratio(1000.0, 300.0) - 1.0);
        assert!((m_prop - expected).abs() < 0.01);
    }

    #[test]
    fn test_propellant_database() {
        let db = propellant_database();
        assert!(db.len() >= 4);
        assert!(db.iter().any(|p| p.name == "LOX/RP-1"));
    }
}
