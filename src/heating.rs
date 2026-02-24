//! Aerodynamic heating models for atmospheric entry.
//!
//! Stagnation-point convective and radiative heating correlations:
//! - **Sutton-Graves** convective heating (already in reentry, expanded here)
//! - **Fay-Riddell** stagnation-point convective heating
//! - **Tauber-Sutton** radiative heating for Earth entry
//! - Integrated heat load computation
//! - Thermal protection system (TPS) material database
//!
//! References:
//! - Fay, J.A. & Riddell, F.R. "Theory of Stagnation Point Heat Transfer in
//!   Dissociated Air", Journal of the Aeronautical Sciences, 25(2), 1958
//! - Sutton, K. & Graves, R.A. "A General Stagnation-Point Convective-Heating
//!   Equation for Arbitrary Gas Mixtures", NASA TR R-376, 1971
//! - Tauber, M.E. & Sutton, K. "Stagnation-Point Radiative Heating Relations
//!   for Earth and Mars Entries", J. Spacecraft, 28(1), 1991

/// Sutton-Graves convective heating correlation.
///
/// q_conv = k * sqrt(rho / r_n) * V³
///
/// where k = 1.7415e-4 (W·s³/(kg^0.5·m)) for Earth air.
///
/// # Arguments
/// * `rho` — freestream density (kg/m³)
/// * `v` — velocity (m/s)
/// * `r_n` — nose radius (m)
///
/// # Returns
/// Convective heat flux (W/m²)
pub fn sutton_graves_convective(rho: f64, v: f64, r_n: f64) -> f64 {
    let k = 1.7415e-4; // Earth air constant
    k * (rho / r_n).sqrt() * v.powi(3)
}

/// Fay-Riddell stagnation-point convective heating.
///
/// More accurate than Sutton-Graves; accounts for wall enthalpy and
/// Lewis number effects. Simplified form:
///
/// q_conv = 0.94 * (rho_s * mu_s)^0.4 * (rho_inf * mu_inf)^0.1
///          * sqrt(du_e/dx) * (h_s - h_w)
///
/// In practice, the correlation is often simplified to:
/// q_FR = C * sqrt(rho_inf / r_n) * V^n * f(h_w/h_s)
///
/// We use the engineering form from Anderson, "Hypersonic and High-Temperature
/// Gas Dynamics", 2nd ed:
///
/// q_FR = 1.83e-4 * sqrt(rho / r_n) * V^3 * (1 - h_w/h_0)
///
/// where h_0 = V²/2 + h_inf is total enthalpy, h_w is wall enthalpy.
///
/// # Arguments
/// * `rho` — freestream density (kg/m³)
/// * `v` — velocity (m/s)
/// * `r_n` — nose radius (m)
/// * `t_wall` — wall temperature (K)
///
/// # Returns
/// Convective heat flux (W/m²)
pub fn fay_riddell_convective(rho: f64, v: f64, r_n: f64, t_wall: f64) -> f64 {
    let k = 1.83e-4;
    let cp_air = 1005.0; // J/(kg·K) at moderate temperatures
    let h_w = cp_air * t_wall;
    let h_0 = v * v / 2.0 + cp_air * 300.0; // total enthalpy (assume T_inf ~ 300K for simplicity)
    let wall_factor = (1.0 - h_w / h_0).max(0.0);
    k * (rho / r_n).sqrt() * v.powi(3) * wall_factor
}

/// Tauber-Sutton radiative heating for Earth entry.
///
/// q_rad = C * r_n^a * rho^b * f(V)
///
/// For Earth entry at velocities > 9 km/s:
/// q_rad = 4.736e4 * r_n^0.5 * rho^1.22 * V^(8.5)  (simplified form, CGS-derived)
///
/// More precisely (Tauber & Sutton 1991):
/// q_rad = C * r_n^a * rho^b * V^d
///
/// where C, a, b, d depend on the velocity regime.
///
/// Valid for V > 7.5 km/s, rho > 1e-5 kg/m³
///
/// # Arguments
/// * `rho` — freestream density (kg/m³)
/// * `v` — velocity (m/s)
/// * `r_n` — nose radius (m)
///
/// # Returns
/// Radiative heat flux (W/m²)
pub fn tauber_sutton_radiative(rho: f64, v: f64, r_n: f64) -> f64 {
    let v_kms = v / 1000.0;

    if v_kms < 7.5 || rho < 1e-6 {
        // Below threshold, radiative heating negligible
        return 0.0;
    }

    // Tauber-Sutton correlation coefficients for Earth entry
    // From Table 1 in Tauber & Sutton (1991)
    let (c, a, b, d) = if v_kms < 10.0 {
        // 7.5 – 10 km/s regime
        (1.072e6, 0.5, 1.22, 6.0)
    } else if v_kms < 12.0 {
        // 10 – 12 km/s
        (2.0e5, 0.5, 1.30, 7.0)
    } else {
        // > 12 km/s (superorbital)
        (4.0e4, 0.5, 1.35, 8.5)
    };

    c * r_n.powf(a) * rho.powf(b) * v_kms.powf(d)
}

/// Total stagnation-point heating (convective + radiative).
pub fn total_stagnation_heating(rho: f64, v: f64, r_n: f64, t_wall: f64) -> HeatingResult {
    let q_conv = fay_riddell_convective(rho, v, r_n, t_wall);
    let q_rad = tauber_sutton_radiative(rho, v, r_n);
    HeatingResult {
        q_convective: q_conv,
        q_radiative: q_rad,
        q_total: q_conv + q_rad,
    }
}

/// Result of a heating calculation.
#[derive(Debug, Clone, Copy)]
pub struct HeatingResult {
    /// Convective heat flux (W/m²)
    pub q_convective: f64,
    /// Radiative heat flux (W/m²)
    pub q_radiative: f64,
    /// Total heat flux (W/m²)
    pub q_total: f64,
}

/// TPS material properties.
#[derive(Debug, Clone)]
pub struct TPSMaterial {
    pub name: &'static str,
    /// Density (kg/m³)
    pub density: f64,
    /// Thermal conductivity (W/(m·K))
    pub conductivity: f64,
    /// Specific heat (J/(kg·K))
    pub specific_heat: f64,
    /// Maximum service temperature (K)
    pub max_temperature: f64,
    /// Ablation energy (J/kg) — for ablative materials
    pub ablation_energy: Option<f64>,
    /// Emissivity (dimensionless)
    pub emissivity: f64,
}

/// Common TPS materials database.
pub fn tps_materials() -> Vec<TPSMaterial> {
    vec![
        TPSMaterial {
            name: "PICA (Phenolic Impregnated Carbon Ablator)",
            density: 265.0,
            conductivity: 0.21,
            specific_heat: 1600.0,
            max_temperature: 3300.0,
            ablation_energy: Some(12.0e6),
            emissivity: 0.9,
        },
        TPSMaterial {
            name: "SLA-561V (Super Lightweight Ablator)",
            density: 256.0,
            conductivity: 0.12,
            specific_heat: 1250.0,
            max_temperature: 2500.0,
            ablation_energy: Some(8.0e6),
            emissivity: 0.85,
        },
        TPSMaterial {
            name: "Avcoat (Apollo heat shield)",
            density: 544.0,
            conductivity: 0.40,
            specific_heat: 1260.0,
            max_temperature: 3500.0,
            ablation_energy: Some(15.0e6),
            emissivity: 0.9,
        },
        TPSMaterial {
            name: "RCC (Reinforced Carbon-Carbon)",
            density: 1600.0,
            conductivity: 40.0,
            specific_heat: 710.0,
            max_temperature: 1920.0,
            ablation_energy: None, // reusable
            emissivity: 0.85,
        },
        TPSMaterial {
            name: "HRSI tiles (Space Shuttle)",
            density: 352.0,
            conductivity: 0.07,
            specific_heat: 630.0,
            max_temperature: 1530.0,
            ablation_energy: None, // reusable
            emissivity: 0.86,
        },
        TPSMaterial {
            name: "SPAM (Silicone Phenolic Ablative Material)",
            density: 480.0,
            conductivity: 0.30,
            specific_heat: 1100.0,
            max_temperature: 3000.0,
            ablation_energy: Some(10.0e6),
            emissivity: 0.88,
        },
    ]
}

/// Estimate ablative TPS thickness required.
///
/// Simple 1-D energy balance: Q_total = rho_tps * h_abl * thickness
///
/// # Arguments
/// * `total_heat_load` — integrated heat load (J/m²)
/// * `material` — TPS material properties
/// * `safety_factor` — typically 1.5–2.0
///
/// # Returns
/// Required TPS thickness (m)
pub fn ablative_tps_thickness(
    total_heat_load: f64,
    material: &TPSMaterial,
    safety_factor: f64,
) -> f64 {
    match material.ablation_energy {
        Some(h_abl) => safety_factor * total_heat_load / (material.density * h_abl),
        None => {
            // For reusable TPS, use thermal soak model (simplified)
            // thickness ~ Q / (rho * cp * delta_T)
            let delta_t = material.max_temperature - 400.0; // assume structure limit ~400K
            safety_factor * total_heat_load / (material.density * material.specific_heat * delta_t)
        }
    }
}

/// Equilibrium wall temperature from radiative cooling.
///
/// q_in * (1 - alpha) = epsilon * sigma * T_w^4
///
/// Solving: T_w = (q_in / (epsilon * sigma))^(1/4)
///
/// # Arguments
/// * `q_in` — incident heat flux (W/m²)
/// * `emissivity` — surface emissivity
///
/// # Returns
/// Equilibrium wall temperature (K)
pub fn equilibrium_wall_temperature(q_in: f64, emissivity: f64) -> f64 {
    use crate::constants::SIGMA_SB;
    (q_in / (emissivity * SIGMA_SB)).powf(0.25)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sutton_graves_sea_level_low_speed() {
        // At very low speed, heating should be small
        let q = sutton_graves_convective(1.225, 100.0, 1.0);
        // 1.7415e-4 * sqrt(1.225) * 100^3 = 1.7415e-4 * 1.107 * 1e6 ≈ 192.7 W/m²
        assert!(
            (q - 192.8).abs() < 5.0,
            "Low speed heating = {q} W/m², expected ~193"
        );
    }

    #[test]
    fn test_sutton_graves_reentry_conditions() {
        // Typical LEO reentry: V=7500 m/s, alt~60km (rho≈3e-4), r_n=1m
        let q = sutton_graves_convective(3e-4, 7500.0, 1.0);
        // ~1.2 MW/m² — reasonable for LEO reentry
        assert!(
            q > 0.5e6 && q < 5e6,
            "LEO reentry heating = {q} W/m², expected ~1-2 MW/m²"
        );
    }

    #[test]
    fn test_fay_riddell_with_wall_temp() {
        // Hot wall reduces heating
        let q_cold = fay_riddell_convective(3e-4, 7500.0, 1.0, 300.0);
        let q_hot = fay_riddell_convective(3e-4, 7500.0, 1.0, 2000.0);
        assert!(
            q_hot < q_cold,
            "Hot wall should reduce heating: {q_hot} vs {q_cold}"
        );
    }

    #[test]
    fn test_tauber_sutton_below_threshold() {
        // V < 7.5 km/s → no radiative heating
        let q = tauber_sutton_radiative(1e-4, 6000.0, 1.0);
        assert!(q == 0.0, "No radiative heating below 7.5 km/s");
    }

    #[test]
    fn test_tauber_sutton_superorbital() {
        // Superorbital entry: V=12 km/s, moderate density
        let q = tauber_sutton_radiative(1e-4, 12000.0, 1.0);
        assert!(q > 0.0, "Should have radiative heating at 12 km/s");
    }

    #[test]
    fn test_nose_radius_effect() {
        // Larger nose radius reduces convective heating
        let q_small = sutton_graves_convective(3e-4, 7500.0, 0.5);
        let q_large = sutton_graves_convective(3e-4, 7500.0, 2.0);
        assert!(
            q_small > q_large,
            "Smaller nose radius = more heating: {q_small} vs {q_large}"
        );
    }

    #[test]
    fn test_equilibrium_wall_temp() {
        // At 1 MW/m² with emissivity 0.9, T_w should be ~2000-3000 K
        let t_w = equilibrium_wall_temperature(1.0e6, 0.9);
        assert!(
            t_w > 1800.0 && t_w < 3500.0,
            "Equilibrium wall temp = {t_w} K, expected ~2000-3000 K"
        );
    }

    #[test]
    fn test_tps_thickness_pica() {
        let materials = tps_materials();
        let pica = &materials[0]; // PICA

        // Apollo-like entry: ~50 MJ/m² total heat load
        let thickness = ablative_tps_thickness(50.0e6, pica, 1.5);
        // Should be a few cm
        assert!(
            thickness > 0.01 && thickness < 0.10,
            "PICA thickness = {thickness} m for 50 MJ/m², expected 2-8 cm"
        );
    }

    #[test]
    fn test_tps_materials_count() {
        let mats = tps_materials();
        assert!(mats.len() >= 5, "Should have at least 5 TPS materials");
    }

    #[test]
    fn test_total_heating() {
        let result = total_stagnation_heating(3e-4, 7500.0, 1.0, 1500.0);
        assert!(result.q_convective > 0.0);
        assert!(result.q_total >= result.q_convective);
    }
}
