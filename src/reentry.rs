//! Re-entry physics: ballistic entry, skip entry, heating, thermal protection
//!
//! Comprehensive models for spacecraft atmospheric entry including:
//! - Allen-Eggers ballistic entry equations  
//! - Skip entry trajectory analysis
//! - Sutton-Graves stagnation point heating
//! - Convective and radiative heat transfer
//! - Thermal Protection System (TPS) sizing
//! - Entry corridor analysis
//!
//! References:
//! - Allen, H.J. & Eggers, A.J. "A study of the motion and aerodynamic heating
//!   of missiles entering the earth's atmosphere at high supersonic speeds" (1958)
//! - Sutton & Graves "A General Stagnation-Point Convective Heating Equation" (1971)  
//! - Tauber & Sutton "Stagnation-Point Radiative Heating Relations for Earth Entry" (1991)
//! - Howe et al. "Hypervelocity Atmospheric Flight: Real Gas Flow Fields" (1989)

use crate::atmosphere::AtmosphereModel;
use crate::constants::{G0, R_EARTH};
use nalgebra::Vector3;

/// Atmospheric entry vehicle configuration
#[derive(Debug, Clone)]
pub struct EntryVehicle {
    /// Vehicle mass (kg)
    pub mass: f64,
    /// Reference area for drag calculation (m²)  
    pub area: f64,
    /// Drag coefficient (dimensionless)
    pub cd: f64,
    /// Ballistic coefficient (kg/m²) - mass/(cd * area)
    pub ballistic_coefficient: f64,
    /// Nose radius for heating calculations (m)
    pub nose_radius: f64,
    /// Vehicle shape factor for heating (sphere=1.0, blunt cone=0.94, etc.)
    pub shape_factor: f64,
    /// Heat capacity of vehicle structure (J/kg/K)
    pub heat_capacity: f64,
    /// Vehicle dimensions
    pub dimensions: VehicleDimensions,
    /// Thermal protection system properties
    pub tps: Option<TPSProperties>,
}

/// Vehicle geometric dimensions
#[derive(Debug, Clone)]
pub struct VehicleDimensions {
    /// Overall length (m)
    pub length: f64,
    /// Maximum diameter (m)  
    pub diameter: f64,
    /// Heat shield diameter (m)
    pub heat_shield_diameter: f64,
    /// Wetted area for heating (m²)
    pub wetted_area: f64,
}

/// Thermal Protection System properties
#[derive(Debug, Clone)]
pub struct TPSProperties {
    /// TPS material type
    pub material: TPSMaterial,
    /// Thermal conductivity (W/m/K)
    pub thermal_conductivity: f64,
    /// Density (kg/m³)
    pub density: f64,
    /// Specific heat capacity (J/kg/K)
    pub specific_heat: f64,
    /// Emissivity (dimensionless)
    pub emissivity: f64,
    /// Ablation temperature (K) - for ablative TPS
    pub ablation_temp: Option<f64>,
    /// Heat of ablation (J/kg) - for ablative TPS  
    pub heat_of_ablation: Option<f64>,
}

/// TPS material types
#[derive(Debug, Clone, PartialEq)]
pub enum TPSMaterial {
    /// Space Shuttle-type tiles (silica)
    ReusableTile,
    /// PICA-X (SpaceX Dragon)
    PicaX,
    /// Avcoat (Apollo, Orion)
    Avcoat,
    /// Carbon-Carbon (Space Shuttle nose/wing leading edges)
    CarbonCarbon,
    /// SIRCA (Stardust)
    Sirca,
    /// Custom material
    Custom(String),
}

/// Atmospheric entry state
#[derive(Debug, Clone)]
pub struct EntryState {
    /// Time since entry interface (s)
    pub time: f64,
    /// Position vector (m) - planetocentric
    pub position: Vector3<f64>,
    /// Velocity vector (m/s) - inertial
    pub velocity: Vector3<f64>,
    /// Altitude above surface (m)
    pub altitude: f64,
    /// Velocity magnitude (m/s)
    pub speed: f64,
    /// Flight path angle (radians) - negative is descending
    pub gamma: f64,
    /// Atmospheric density (kg/m³)
    pub density: f64,
    /// Atmospheric pressure (Pa)
    pub pressure: f64,
    /// Atmospheric temperature (K)
    pub temperature: f64,
    /// Dynamic pressure (Pa) - 0.5 * rho * V²
    pub dynamic_pressure: f64,
    /// Deceleration magnitude (g's)
    pub deceleration_g: f64,
}

/// Heat transfer calculations at entry state
#[derive(Debug, Clone)]
pub struct EntryHeating {
    /// Stagnation point convective heat flux (W/m²)
    pub q_conv_stagnation: f64,
    /// Stagnation point radiative heat flux (W/m²)  
    pub q_rad_stagnation: f64,
    /// Total stagnation point heat flux (W/m²)
    pub q_total_stagnation: f64,
    /// Peak heat flux over vehicle surface (W/m²)
    pub q_peak: f64,
    /// Stagnation point temperature (K)
    pub temp_stagnation: f64,
    /// Heat load integrated over time (J/m²)
    pub heat_load: f64,
    /// TPS surface temperature (K)
    pub tps_surface_temp: f64,
    /// TPS back wall temperature (K)
    pub tps_backwall_temp: f64,
}

/// Entry trajectory propagation result
#[derive(Debug, Clone)]
pub struct EntryTrajectory {
    /// Time history (s)
    pub time: Vec<f64>,
    /// State history
    pub states: Vec<EntryState>,
    /// Heating history
    pub heating: Vec<EntryHeating>,
    /// Maximum deceleration experienced (g's)
    pub max_deceleration: f64,
    /// Maximum heat flux (W/m²)
    pub max_heat_flux: f64,
    /// Maximum heat load (J/m²)
    pub max_heat_load: f64,
    /// Total entry duration (s)
    pub entry_duration: f64,
    /// Final velocity at landing/parachute deployment (m/s)
    pub final_velocity: f64,
}

/// Entry analysis parameters
#[derive(Debug, Clone)]
pub struct EntryAnalysis {
    /// Entry interface altitude (m) - typically 120-125 km
    pub entry_altitude: f64,
    /// Entry velocity magnitude (m/s)
    pub entry_velocity: f64,
    /// Entry flight path angle (radians) - negative for descent
    pub entry_gamma: f64,
    /// Entry azimuth angle (radians)
    pub entry_azimuth: f64,
    /// Target altitude for analysis end (m) - e.g., parachute deployment
    pub target_altitude: f64,
    /// Integration time step (s)
    pub time_step: f64,
    /// Maximum integration time (s)
    pub max_time: f64,
}

/// Entry corridor bounds for skip entry analysis
#[derive(Debug, Clone)]
pub struct EntryCorridor {
    /// Minimum entry angle for successful entry (radians)
    pub gamma_min: f64,
    /// Maximum entry angle to avoid excessive heating (radians)  
    pub gamma_max: f64,
    /// Entry corridor width (radians)
    pub corridor_width: f64,
    /// Skip-out velocity threshold (m/s)
    pub skip_out_velocity: f64,
}

/// Main entry physics simulator
pub struct EntrySimulator {
    /// Vehicle configuration
    pub vehicle: EntryVehicle,
    /// Atmospheric model
    pub atmosphere: AtmosphereModel,
    /// Planet radius (m) - defaults to Earth
    pub planet_radius: f64,
    /// Planet gravitational parameter (m³/s²)
    pub mu: f64,
}

impl EntryVehicle {
    /// Create Apollo Command Module configuration
    pub fn apollo_cm() -> Self {
        Self {
            mass: 5800.0, // CM mass at entry
            area: 12.02,  // Heat shield area (m²)
            cd: 1.3,      // Blunt body drag coefficient
            ballistic_coefficient: 5800.0 / (1.3 * 12.02),
            nose_radius: 4.7 / 2.0, // Heat shield radius
            shape_factor: 0.94,     // Blunt cone shape factor
            heat_capacity: 900.0,   // Aluminum/steel structure
            dimensions: VehicleDimensions {
                length: 3.2,
                diameter: 3.9,
                heat_shield_diameter: 3.9,
                wetted_area: 35.0, // Total wetted area estimate
            },
            tps: Some(TPSProperties::avcoat()),
        }
    }

    /// Create SpaceX Dragon configuration
    pub fn dragon_v2() -> Self {
        Self {
            mass: 6400.0, // Dragon 2 mass estimate
            area: 12.56,  // Heat shield area (π * r²)
            cd: 1.2,
            ballistic_coefficient: 6400.0 / (1.2 * 12.56),
            nose_radius: 2.0,
            shape_factor: 0.94,
            heat_capacity: 900.0,
            dimensions: VehicleDimensions {
                length: 8.1,
                diameter: 4.0,
                heat_shield_diameter: 4.0,
                wetted_area: 40.0,
            },
            tps: Some(TPSProperties::pica_x()),
        }
    }

    /// Create Orion MPCV configuration  
    pub fn orion() -> Self {
        Self {
            mass: 8500.0, // Orion CM mass
            area: 19.63,  // Heat shield area
            cd: 1.2,
            ballistic_coefficient: 8500.0 / (1.2 * 19.63),
            nose_radius: 2.5,
            shape_factor: 0.94,
            heat_capacity: 900.0,
            dimensions: VehicleDimensions {
                length: 3.3,
                diameter: 5.0,
                heat_shield_diameter: 5.0,
                wetted_area: 50.0,
            },
            tps: Some(TPSProperties::avcoat()),
        }
    }

    /// Create Mars entry vehicle (Mars Science Laboratory style)
    pub fn mars_entry_vehicle() -> Self {
        Self {
            mass: 2400.0, // MSL entry mass
            area: 15.9,   // Heat shield area
            cd: 1.24,
            ballistic_coefficient: 2400.0 / (1.24 * 15.9),
            nose_radius: 2.25,
            shape_factor: 0.94,
            heat_capacity: 900.0,
            dimensions: VehicleDimensions {
                length: 3.0,
                diameter: 4.5,
                heat_shield_diameter: 4.5,
                wetted_area: 42.0,
            },
            tps: Some(TPSProperties::pica_x()),
        }
    }
}

impl TPSProperties {
    /// Avcoat (Apollo/Orion TPS)
    pub fn avcoat() -> Self {
        Self {
            material: TPSMaterial::Avcoat,
            thermal_conductivity: 0.11,
            density: 520.0,
            specific_heat: 1050.0,
            emissivity: 0.85,
            ablation_temp: Some(600.0 + 273.15), // 600°C
            heat_of_ablation: Some(2.1e6),
        }
    }

    /// PICA-X (SpaceX Dragon)
    pub fn pica_x() -> Self {
        Self {
            material: TPSMaterial::PicaX,
            thermal_conductivity: 0.576,
            density: 270.0,
            specific_heat: 1200.0,
            emissivity: 0.89,
            ablation_temp: Some(1850.0), // Very high temperature resistance
            heat_of_ablation: Some(8.0e6),
        }
    }

    /// Space Shuttle tiles
    pub fn shuttle_tiles() -> Self {
        Self {
            material: TPSMaterial::ReusableTile,
            thermal_conductivity: 0.12,
            density: 350.0,
            specific_heat: 1050.0,
            emissivity: 0.85,
            ablation_temp: None, // Reusable, doesn't ablate
            heat_of_ablation: None,
        }
    }
}

impl EntrySimulator {
    /// Create Earth entry simulator
    pub fn earth(vehicle: EntryVehicle) -> Self {
        Self {
            vehicle,
            atmosphere: AtmosphereModel::us_standard_1976(),
            planet_radius: R_EARTH * 1000.0, // Convert to meters
            mu: 3.986004418e14,              // Earth gravitational parameter
        }
    }

    /// Create Mars entry simulator  
    pub fn mars(vehicle: EntryVehicle) -> Self {
        Self {
            vehicle,
            atmosphere: AtmosphereModel::mars(), // Would need Mars atmosphere implementation
            planet_radius: 3389.5e3,
            mu: 4.2828e13,
        }
    }

    /// Propagate entry trajectory using Allen-Eggers equations
    pub fn propagate_ballistic_entry(
        &self,
        analysis: &EntryAnalysis,
    ) -> Result<EntryTrajectory, String> {
        let mut time = Vec::new();
        let mut states = Vec::new();
        let mut heating = Vec::new();

        // Initial conditions
        let mut t = 0.0;
        let mut alt = analysis.entry_altitude;
        let mut v = analysis.entry_velocity;
        let mut gamma = analysis.entry_gamma;

        let mut max_deceleration: f64 = 0.0;
        let mut max_heat_flux: f64 = 0.0;
        let mut max_heat_load: f64 = 0.0;

        while t < analysis.max_time && alt > analysis.target_altitude && v > 100.0 {
            // Get atmospheric properties
            let (rho, p, temp) = match self.atmosphere.density_pressure_temperature(alt) {
                Ok(values) => values,
                Err(_) => break, // Outside atmosphere model range
            };

            // Calculate current state
            let _speed = v;
            let dynamic_pressure = 0.5 * rho * v * v;
            let drag = self.vehicle.cd * self.vehicle.area * dynamic_pressure;
            let deceleration = drag / self.vehicle.mass;
            let deceleration_g = deceleration / G0;

            // Position and velocity vectors (simplified 2D)
            let position = Vector3::new(0.0, 0.0, self.planet_radius + alt);
            let velocity = Vector3::new(v * gamma.cos(), 0.0, v * gamma.sin());

            let state = EntryState {
                time: t,
                position,
                velocity,
                altitude: alt,
                speed: v,
                gamma,
                density: rho,
                pressure: p,
                temperature: temp,
                dynamic_pressure,
                deceleration_g,
            };

            // Calculate heating
            let heating_calc = self.calculate_heating(&state);

            // Update trajectory using Allen-Eggers equations
            let dt = analysis.time_step;
            let radius = self.planet_radius + alt;
            let g = self.mu / (radius * radius);

            // Drag deceleration along velocity vector
            let dvdt_drag = -drag / self.vehicle.mass;

            // Gravity component along velocity vector
            let dvdt_gravity = -g * gamma.sin();

            // Gravity component perpendicular to velocity
            let dgammadt = -(g / v) * gamma.cos() + (v / radius) * gamma.cos();

            // Total velocity change
            let dv = (dvdt_drag + dvdt_gravity) * dt;
            let dgamma = dgammadt * dt;
            let dalt = v * gamma.sin() * dt;

            // Update state
            v += dv;
            gamma += dgamma;
            alt += dalt;
            t += dt;

            // Track maximums
            max_deceleration = max_deceleration.max(deceleration_g);
            max_heat_flux = max_heat_flux.max(heating_calc.q_total_stagnation);
            max_heat_load = max_heat_load.max(heating_calc.heat_load);

            // Check for unrealistic conditions before storing
            if deceleration_g > 50.0 || heating_calc.q_total_stagnation > 1e7 {
                break; // Vehicle would not survive
            }

            // Store results
            time.push(t);
            states.push(state);
            heating.push(heating_calc);
        }

        Ok(EntryTrajectory {
            time,
            states,
            heating,
            max_deceleration,
            max_heat_flux,
            max_heat_load,
            entry_duration: t,
            final_velocity: v,
        })
    }

    /// Calculate heating using Sutton-Graves and radiative heating models
    fn calculate_heating(&self, state: &EntryState) -> EntryHeating {
        let v = state.speed;
        let rho = state.density;
        let temp = state.temperature;
        let rn = self.vehicle.nose_radius;

        // Sutton-Graves convective heating at stagnation point
        // q_conv = 1.7415e-4 * sqrt(rho/rn) * V^3.15 (W/m²)
        let q_conv_stagnation = 1.7415e-4 * (rho / rn).sqrt() * v.powf(3.15);

        // Radiative heating (Tauber-Sutton correlation)
        // Significant at velocities > 10 km/s
        let q_rad_stagnation = if v > 8000.0 {
            // Simplified radiative heating model
            let c_rad = 1.59e-13; // Radiative heating coefficient
            c_rad * rho.sqrt() * v.powf(7.5) / rn
        } else {
            0.0
        };

        let q_total_stagnation = q_conv_stagnation + q_rad_stagnation;

        // Peak heating over vehicle surface (typically ~0.7x stagnation point)
        let q_peak = q_total_stagnation * 0.7;

        // Stagnation point temperature (simplified)
        let sigma = 5.67e-8; // Stefan-Boltzmann constant
        let emissivity = self
            .vehicle
            .tps
            .as_ref()
            .map(|tps| tps.emissivity)
            .unwrap_or(0.85);

        // Energy balance: q_in = q_radiation_out
        // q = ε * σ * T^4, so T = (q / (ε * σ))^0.25
        let temp_stagnation = (q_total_stagnation / (emissivity * sigma)).powf(0.25);

        // Heat load integration (simplified - would need proper integration)
        let heat_load = q_total_stagnation * state.time; // Rough approximation

        // TPS temperatures (simplified 1D conduction)
        let (tps_surface_temp, tps_backwall_temp) = if let Some(tps) = &self.vehicle.tps {
            self.calculate_tps_temperatures(q_total_stagnation, tps)
        } else {
            (temp_stagnation, temp)
        };

        EntryHeating {
            q_conv_stagnation,
            q_rad_stagnation,
            q_total_stagnation,
            q_peak,
            temp_stagnation,
            heat_load,
            tps_surface_temp,
            tps_backwall_temp,
        }
    }

    /// Simplified TPS temperature calculation
    fn calculate_tps_temperatures(&self, heat_flux: f64, tps: &TPSProperties) -> (f64, f64) {
        // Very simplified 1D steady-state conduction
        // In reality, this would be a complex transient heat transfer problem

        let thickness = 0.05; // Assume 5cm TPS thickness
        let k = tps.thermal_conductivity;
        let sigma = 5.67e-8;
        let emissivity = tps.emissivity;

        // Surface temperature from energy balance
        // q_in = q_conduction = q_radiation
        // This is a simplification - real analysis needs coupled equations
        let t_surface = (heat_flux / (emissivity * sigma)).powf(0.25);

        // Back wall temperature (very simplified)
        let temp_drop = heat_flux * thickness / k;
        let t_backwall = t_surface - temp_drop;

        (t_surface, t_backwall.max(300.0)) // Minimum room temperature
    }

    /// Analyze entry corridor for skip entry
    pub fn analyze_entry_corridor(
        &self,
        _entry_velocity: f64,
        entry_altitude: f64,
    ) -> EntryCorridor {
        // Shallow entry angle limit (skip-out boundary)
        let gamma_min = -0.5_f64.to_radians(); // Very shallow

        // Steep entry angle limit (excessive heating/deceleration)
        let gamma_max = -15.0_f64.to_radians(); // Steep but survivable

        // Skip-out velocity (vehicle escapes atmosphere)
        let escape_velocity = (2.0 * self.mu / (self.planet_radius + entry_altitude)).sqrt();
        let skip_out_velocity = 0.9 * escape_velocity; // 90% of escape velocity

        EntryCorridor {
            gamma_min,
            gamma_max,
            corridor_width: (gamma_min - gamma_max).abs(),
            skip_out_velocity,
        }
    }
}

/// Pre-defined entry analysis configurations
pub mod scenarios {
    use super::*;

    /// Apollo 11 lunar return entry
    pub fn apollo_11_entry() -> EntryAnalysis {
        EntryAnalysis {
            entry_altitude: 120_000.0, // 120 km entry interface
            entry_velocity: 11_100.0,  // Lunar return velocity
            entry_gamma: -6.5_f64.to_radians(),
            entry_azimuth: 0.0,
            target_altitude: 10_000.0, // Parachute deployment
            time_step: 1.0,
            max_time: 800.0, // ~13 minutes
        }
    }

    /// Dragon ISS return entry
    pub fn dragon_iss_return() -> EntryAnalysis {
        EntryAnalysis {
            entry_altitude: 125_000.0,
            entry_velocity: 7_800.0, // LEO deorbit velocity
            entry_gamma: -1.2_f64.to_radians(),
            entry_azimuth: 0.0,
            target_altitude: 15_000.0, // Drogue chute deployment
            time_step: 0.5,
            max_time: 600.0,
        }
    }

    /// Mars entry (MSL-style)
    pub fn mars_entry() -> EntryAnalysis {
        EntryAnalysis {
            entry_altitude: 125_000.0,           // Mars atmospheric entry
            entry_velocity: 5_500.0,             // Mars approach velocity
            entry_gamma: -15.5_f64.to_radians(), // Steep Mars entry
            entry_azimuth: 0.0,
            target_altitude: 10_000.0, // Parachute deployment
            time_step: 0.1,
            max_time: 360.0, // 6 minutes typical Mars entry
        }
    }

    /// High-speed asteroid sample return
    pub fn asteroid_sample_return() -> EntryAnalysis {
        EntryAnalysis {
            entry_altitude: 125_000.0,
            entry_velocity: 12_800.0, // High-speed interplanetary return
            entry_gamma: -8.0_f64.to_radians(),
            entry_azimuth: 0.0,
            target_altitude: 5_000.0,
            time_step: 0.5,
            max_time: 900.0, // Longer entry due to high speed
        }
    }
}

#[cfg(test)]
mod tests {
    use super::scenarios::*;
    use super::*;

    #[test]
    fn test_apollo_cm_configuration() {
        let apollo = EntryVehicle::apollo_cm();

        assert!(apollo.mass > 5000.0);
        assert!(apollo.ballistic_coefficient > 0.0);
        assert_eq!(apollo.tps.unwrap().material, TPSMaterial::Avcoat);
    }

    #[test]
    fn test_dragon_configuration() {
        let dragon = EntryVehicle::dragon_v2();

        assert!(dragon.area > 10.0);
        assert_eq!(dragon.tps.unwrap().material, TPSMaterial::PicaX);
        assert!(dragon.nose_radius > 1.0);
    }

    #[test]
    fn test_entry_simulator_creation() {
        let vehicle = EntryVehicle::apollo_cm();
        let sim = EntrySimulator::earth(vehicle);

        assert!((sim.planet_radius - R_EARTH * 1000.0).abs() < 1000.0);
        assert!(sim.mu > 3.9e14 && sim.mu < 4.0e14);
    }

    #[test]
    fn test_sutton_graves_heating() {
        let vehicle = EntryVehicle::apollo_cm();
        let sim = EntrySimulator::earth(vehicle);

        // High-speed entry state
        let state = EntryState {
            time: 100.0,
            position: Vector3::new(0.0, 0.0, R_EARTH * 1000.0 + 80000.0),
            velocity: Vector3::new(10000.0, 0.0, -1000.0),
            altitude: 80000.0,
            speed: 10050.0, // ~10 km/s
            gamma: -6.0_f64.to_radians(),
            density: 1e-5, // Thin atmosphere at 80km
            pressure: 1.0,
            temperature: 200.0,
            dynamic_pressure: 0.5e-5 * 10050.0 * 10050.0,
            deceleration_g: 5.0,
        };

        let heating = sim.calculate_heating(&state);

        assert!(heating.q_conv_stagnation > 0.0);
        assert!(heating.q_total_stagnation > heating.q_conv_stagnation);
        assert!(heating.temp_stagnation > 1000.0); // High temperature
    }

    #[test]
    fn test_entry_corridor_analysis() {
        let vehicle = EntryVehicle::apollo_cm();
        let sim = EntrySimulator::earth(vehicle);

        let corridor = sim.analyze_entry_corridor(11100.0, 120000.0);

        assert!(corridor.gamma_min < 0.0); // Negative angle for descent
        assert!(corridor.gamma_max < corridor.gamma_min); // More negative = steeper
        assert!(corridor.corridor_width > 0.0);
        assert!(corridor.skip_out_velocity > 8000.0);
    }

    #[test]
    fn test_apollo_entry_scenario() {
        let analysis = apollo_11_entry();

        assert!((analysis.entry_velocity - 11100.0).abs() < 100.0);
        assert!(analysis.entry_gamma < 0.0); // Descending
        assert!(analysis.entry_altitude > 100000.0);
        assert!(analysis.target_altitude > 5000.0);
    }

    #[test]
    fn test_mars_entry_scenario() {
        let analysis = mars_entry();

        assert!(analysis.entry_velocity < 7000.0); // Lower than Earth speeds
        assert!(analysis.entry_gamma < -10.0_f64.to_radians()); // Steep Mars entry
        assert!(analysis.time_step < 1.0); // Need fine resolution for Mars
    }

    #[test]
    fn test_tps_properties() {
        let avcoat = TPSProperties::avcoat();
        let pica_x = TPSProperties::pica_x();

        assert!(avcoat.density > 400.0);
        assert!(pica_x.thermal_conductivity > 0.4);
        assert!(pica_x.ablation_temp.unwrap() > avcoat.ablation_temp.unwrap());
    }

    #[test]
    fn test_ballistic_coefficient_calculation() {
        let vehicle = EntryVehicle::apollo_cm();
        let expected_bc = vehicle.mass / (vehicle.cd * vehicle.area);

        assert!((vehicle.ballistic_coefficient - expected_bc).abs() < 1.0);
        assert!(vehicle.ballistic_coefficient > 300.0); // Reasonable value
    }

    #[test]
    fn test_heating_temperature_calculation() {
        let vehicle = EntryVehicle::dragon_v2();
        let sim = EntrySimulator::earth(vehicle);

        let tps = TPSProperties::pica_x();
        let (t_surface, t_backwall) = sim.calculate_tps_temperatures(1e6, &tps); // 1 MW/m²

        assert!(t_surface > t_backwall);
        assert!(t_surface > 1000.0); // High surface temperature
        assert!(t_backwall >= 300.0); // At or above room temperature
    }

    #[test]
    fn test_radiative_heating_threshold() {
        let vehicle = EntryVehicle::apollo_cm();
        let sim = EntrySimulator::earth(vehicle);

        // Low speed state (no radiative heating)
        let low_speed_state = EntryState {
            time: 100.0,
            position: Vector3::zeros(),
            velocity: Vector3::new(5000.0, 0.0, 0.0),
            altitude: 50000.0,
            speed: 5000.0,
            gamma: -6.0_f64.to_radians(),
            density: 1e-3,
            pressure: 100.0,
            temperature: 250.0,
            dynamic_pressure: 12.5,
            deceleration_g: 2.0,
        };

        let heating_low = sim.calculate_heating(&low_speed_state);

        // High speed state (with radiative heating)
        let mut high_speed_state = low_speed_state.clone();
        high_speed_state.speed = 12000.0;
        high_speed_state.velocity = Vector3::new(12000.0, 0.0, 0.0);

        let heating_high = sim.calculate_heating(&high_speed_state);

        assert!((heating_low.q_rad_stagnation - 0.0).abs() < 1e6); // Minimal radiative
        assert!(heating_high.q_rad_stagnation > 1e6); // Significant radiative
        assert!(heating_high.q_total_stagnation > heating_low.q_total_stagnation);
    }
}
