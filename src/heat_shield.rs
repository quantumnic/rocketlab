//! Ablative heat shield sizing
//!
//! Models thermal protection system (TPS) materials, ablation recession rates,
//! thermal response, and heat shield mass estimation.
//!
//! Reference: Tauber & Sutton (1991), NASA TN-D series,
//! Wright et al. "Afterbody Aeroheating Flight Data" (2006)

use std::f64::consts::PI;

/// TPS material properties
#[derive(Debug, Clone)]
pub struct TpsMaterial {
    pub name: &'static str,
    /// Density (kg/m³)
    pub density: f64,
    /// Effective heat of ablation (J/kg) — energy absorbed per kg of material ablated
    pub heat_of_ablation: f64,
    /// Thermal conductivity (W/m·K)
    pub conductivity: f64,
    /// Specific heat capacity (J/kg·K)
    pub specific_heat: f64,
    /// Maximum service temperature (K)
    pub max_temperature: f64,
    /// Emissivity
    pub emissivity: f64,
}

/// Common TPS materials database
pub fn tps_materials() -> Vec<TpsMaterial> {
    vec![
        TpsMaterial {
            name: "PICA (Phenolic Impregnated Carbon Ablator)",
            density: 270.0,
            heat_of_ablation: 25.0e6,
            conductivity: 0.21,
            specific_heat: 1260.0,
            max_temperature: 3300.0,
            emissivity: 0.85,
        },
        TpsMaterial {
            name: "PICA-X (SpaceX variant)",
            density: 280.0,
            heat_of_ablation: 28.0e6,
            conductivity: 0.20,
            specific_heat: 1300.0,
            max_temperature: 3400.0,
            emissivity: 0.85,
        },
        TpsMaterial {
            name: "SLA-561V (Mars missions)",
            density: 256.0,
            heat_of_ablation: 12.0e6,
            conductivity: 0.12,
            specific_heat: 1047.0,
            max_temperature: 2200.0,
            emissivity: 0.80,
        },
        TpsMaterial {
            name: "Avcoat (Apollo/Orion)",
            density: 513.0,
            heat_of_ablation: 15.0e6,
            conductivity: 0.40,
            specific_heat: 1260.0,
            max_temperature: 3000.0,
            emissivity: 0.90,
        },
        TpsMaterial {
            name: "Carbon-Carbon (Shuttle leading edge)",
            density: 1600.0,
            heat_of_ablation: 50.0e6,
            conductivity: 40.0,
            specific_heat: 710.0,
            max_temperature: 1920.0,
            emissivity: 0.80,
        },
        TpsMaterial {
            name: "SIRCA (Silicone Impregnated Reusable Ceramic)",
            density: 192.0,
            heat_of_ablation: 8.0e6,
            conductivity: 0.08,
            specific_heat: 1046.0,
            max_temperature: 1800.0,
            emissivity: 0.85,
        },
        TpsMaterial {
            name: "Cork (sounding rockets)",
            density: 500.0,
            heat_of_ablation: 3.0e6,
            conductivity: 0.065,
            specific_heat: 1800.0,
            max_temperature: 900.0,
            emissivity: 0.70,
        },
    ]
}

/// Heat shield sizing result
#[derive(Debug, Clone)]
pub struct HeatShieldResult {
    /// Required ablator thickness (m)
    pub ablator_thickness: f64,
    /// Recession depth (m)
    pub recession: f64,
    /// Thermal soak depth (m) — additional margin for in-depth conduction
    pub soak_depth: f64,
    /// Safety margin thickness (m)
    pub margin: f64,
    /// Total shield mass (kg)
    pub total_mass: f64,
    /// Shield area (m²)
    pub area: f64,
    /// Material used
    pub material_name: &'static str,
    /// Peak surface temperature (K)
    pub peak_surface_temp: f64,
}

/// Compute ablation recession for a given total heat load.
///
/// recession = Q_total / (ρ × h_ablation)
///
/// where Q_total is the integrated heat flux (J/m²)
pub fn compute_recession(total_heat_load: f64, material: &TpsMaterial) -> f64 {
    total_heat_load / (material.density * material.heat_of_ablation)
}

/// Compute thermal soak depth (how far heat penetrates during entry).
///
/// Uses semi-infinite solid approximation:
/// δ = 2 × √(α × t)
///
/// where α = k/(ρ×cp) is thermal diffusivity
pub fn thermal_soak_depth(material: &TpsMaterial, duration: f64) -> f64 {
    let alpha = material.conductivity / (material.density * material.specific_heat);
    2.0 * (alpha * duration).sqrt()
}

/// Radiation equilibrium temperature for a given heat flux.
///
/// T = (q / (ε × σ))^(1/4)
///
/// σ = 5.670374419e-8 W/m²·K⁴ (Stefan-Boltzmann)
pub fn radiation_equilibrium_temp(heat_flux: f64, emissivity: f64) -> f64 {
    let sigma = 5.670374419e-8;
    (heat_flux / (emissivity * sigma)).powf(0.25)
}

/// Size an ablative heat shield.
///
/// # Arguments
/// * `total_heat_load` - Integrated heat load over entry (J/m²)
/// * `peak_heat_flux` - Peak heat flux (W/m²)
/// * `entry_duration` - Duration of significant heating (s)
/// * `shield_diameter` - Heat shield diameter (m)
/// * `material` - TPS material
/// * `safety_factor` - Thickness safety factor (typically 1.3-1.5)
pub fn size_heat_shield(
    total_heat_load: f64,
    peak_heat_flux: f64,
    entry_duration: f64,
    shield_diameter: f64,
    material: &TpsMaterial,
    safety_factor: f64,
) -> HeatShieldResult {
    // Recession from ablation
    let recession = compute_recession(total_heat_load, material);

    // Thermal soak depth
    let soak = thermal_soak_depth(material, entry_duration);

    // Minimum thickness = recession + soak depth
    let min_thickness = recession + soak;
    let margin = min_thickness * (safety_factor - 1.0);
    let total_thickness = min_thickness + margin;

    // Shield area (spherical cap approximation — ~hemisphere front face)
    let area = PI / 4.0 * shield_diameter * shield_diameter;

    // Mass
    let mass = material.density * area * total_thickness;

    // Peak surface temperature
    let peak_temp = radiation_equilibrium_temp(peak_heat_flux, material.emissivity);

    HeatShieldResult {
        ablator_thickness: total_thickness,
        recession,
        soak_depth: soak,
        margin,
        total_mass: mass,
        area,
        material_name: material.name,
        peak_surface_temp: peak_temp,
    }
}

/// Estimate total heat load for a ballistic entry.
///
/// Uses Allen-Eggers approximation:
/// Q_total ≈ ½ × ρ_surface × v_entry³ × β_eff / (β_ballistic × sin(γ))
///
/// Simplified: Q ≈ ½ × m × v² / (Cd × A) for the kinetic energy
/// dissipated per unit area.
///
/// # Arguments
/// * `v_entry` - Entry velocity (m/s)
/// * `mass` - Vehicle mass (kg)
/// * `cd` - Drag coefficient
/// * `area` - Reference area (m²)
/// * `gamma` - Entry flight path angle (rad, negative for descent)
pub fn ballistic_heat_load(v_entry: f64, mass: f64, cd: f64, area: f64, gamma: f64) -> f64 {
    let beta = mass / (cd * area); // Ballistic coefficient
                                   // Fraction of kinetic energy going to heating (empirical: ~0.5 for blunt bodies)
    let heat_fraction = 0.5;
    // Total KE per unit area / convective heating fraction
    heat_fraction * beta * v_entry * v_entry * gamma.abs().sin()
}

/// Estimate peak heat flux for a ballistic entry (Sutton-Graves).
///
/// q_peak = k × sqrt(ρ) × V³
/// Simplified: q_peak ∝ V³ × sqrt(ρ_peak)
///
/// ρ_peak occurs at altitude where max deceleration happens.
///
/// # Arguments
/// * `v_entry` - Entry velocity (m/s)  
/// * `nose_radius` - Nose radius (m)
/// * `density_at_peak` - Atmospheric density at peak heating (kg/m³)
pub fn peak_heat_flux_stagnation(v_entry: f64, nose_radius: f64, density_at_peak: f64) -> f64 {
    // Sutton-Graves constant for Earth (air): k ≈ 1.7415e-4
    let k_sg = 1.7415e-4;
    k_sg * (density_at_peak / nose_radius).sqrt() * v_entry.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tps_database() {
        let db = tps_materials();
        assert!(db.len() >= 6);
        for mat in &db {
            assert!(mat.density > 0.0);
            assert!(mat.heat_of_ablation > 0.0);
            assert!(mat.emissivity > 0.0 && mat.emissivity <= 1.0);
        }
    }

    #[test]
    fn test_recession() {
        let pica = &tps_materials()[0]; // PICA
                                        // 100 MJ/m² heat load
        let rec = compute_recession(100.0e6, pica);
        // recession = 100e6 / (270 * 25e6) = 0.0148 m ≈ 1.5 cm
        assert!((rec - 0.0148).abs() < 0.002, "recession = {rec} m");
    }

    #[test]
    fn test_thermal_soak() {
        let pica = &tps_materials()[0];
        let soak = thermal_soak_depth(pica, 300.0);
        // α = 0.21 / (270 * 1260) ≈ 6.17e-7 m²/s
        // δ = 2 * sqrt(6.17e-7 * 300) ≈ 0.027 m
        assert!(soak > 0.01 && soak < 0.05, "soak = {soak} m");
    }

    #[test]
    fn test_radiation_equilibrium() {
        // At 1 MW/m², typical LEO return
        let t = radiation_equilibrium_temp(1.0e6, 0.85);
        // T = (1e6 / (0.85 * 5.67e-8))^0.25 ≈ 2070 K
        assert!(t > 1900.0 && t < 2200.0, "T = {t} K");
    }

    #[test]
    fn test_size_heat_shield_leo_return() {
        let pica = &tps_materials()[0];
        // LEO return: ~100 MJ/m² total, ~1 MW/m² peak, ~300s heating, 4m diameter
        let result = size_heat_shield(100.0e6, 1.0e6, 300.0, 4.0, pica, 1.4);

        assert!(result.recession > 0.01, "recession = {}", result.recession);
        assert!(
            result.ablator_thickness > 0.03 && result.ablator_thickness < 0.15,
            "thickness = {} m",
            result.ablator_thickness
        );
        assert!(
            result.total_mass > 5.0 && result.total_mass < 200.0,
            "mass = {} kg",
            result.total_mass
        );
        assert!(result.peak_surface_temp > 1500.0);
    }

    #[test]
    fn test_size_heat_shield_lunar_return() {
        // Lunar return: ~300 MJ/m² total, ~5 MW/m² peak, ~200s, 5m diameter
        let avcoat = &tps_materials()[3]; // Avcoat
        let result = size_heat_shield(300.0e6, 5.0e6, 200.0, 5.0, avcoat, 1.5);

        // Avcoat should be thicker for lunar return
        assert!(
            result.ablator_thickness > 0.03,
            "thickness = {} m",
            result.ablator_thickness
        );
        assert!(result.total_mass > 50.0, "mass = {} kg", result.total_mass);
    }

    #[test]
    fn test_ballistic_heat_load() {
        // LEO return: v = 7800 m/s, 5000 kg, Cd = 1.2, 12.5 m² area, γ = -3°
        let q = ballistic_heat_load(7800.0, 5000.0, 1.2, 12.5, (-3.0_f64).to_radians());
        // Should be on the order of 100 MJ/m²
        assert!(q > 10.0e6 && q < 1000.0e6, "Q = {} MJ/m²", q / 1e6);
    }

    #[test]
    fn test_peak_heat_flux_stagnation() {
        // LEO return: v = 7800 m/s, nose R = 2m, ρ at 50km ≈ 1e-3 kg/m³
        let q_peak = peak_heat_flux_stagnation(7800.0, 2.0, 1.0e-3);
        // Should be ~100s of kW/m² to low MW/m²
        assert!(
            q_peak > 100_000.0 && q_peak < 10_000_000.0,
            "q_peak = {} W/m²",
            q_peak
        );
    }

    #[test]
    fn test_mars_entry_heat_shield() {
        // Mars entry: ~60 MJ/m², 2 MW/m² peak, 100s, 4.5m
        let sla = &tps_materials()[2]; // SLA-561V
        let result = size_heat_shield(60.0e6, 2.0e6, 100.0, 4.5, sla, 1.3);
        assert!(result.recession > 0.01);
        assert!(result.total_mass > 10.0);
    }
}
