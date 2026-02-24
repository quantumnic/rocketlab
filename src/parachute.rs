//! Parachute descent systems
//!
//! Models drogue and main parachute staging, opening shock loads,
//! terminal velocity computation, and descent timeline analysis.
//!
//! Reference: Knacke "Parachute Recovery Systems Design Manual" (1991),
//! Ewing, Bixby & Knacke "Recovery Systems Design Guide" (1978)

use std::f64::consts::PI;

use crate::atmosphere;

/// Parachute type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParachuteType {
    /// Drogue: stabilization & initial deceleration (Cd ≈ 0.5-0.65)
    Drogue,
    /// Conical ribbon (Cd ≈ 0.5-0.55)
    ConicalRibbon,
    /// Ringsail (Cd ≈ 0.75-0.85)
    Ringsail,
    /// Flat circular (Cd ≈ 0.75-0.80)
    FlatCircular,
    /// Cross (Cd ≈ 0.6-0.75)
    Cross,
    /// Hemisflo (Cd ≈ 0.3-0.4, supersonic)
    Hemisflo,
    /// Disk-Gap-Band (Cd ≈ 0.5-0.6, used on Viking/MSL)
    DiskGapBand,
}

impl ParachuteType {
    /// Nominal drag coefficient
    pub fn cd_nominal(&self) -> f64 {
        match self {
            Self::Drogue => 0.55,
            Self::ConicalRibbon => 0.52,
            Self::Ringsail => 0.80,
            Self::FlatCircular => 0.78,
            Self::Cross => 0.67,
            Self::Hemisflo => 0.35,
            Self::DiskGapBand => 0.55,
        }
    }

    /// Opening load factor Ck (ratio of opening shock to steady-state drag)
    pub fn opening_load_factor(&self) -> f64 {
        match self {
            Self::Drogue => 1.0,
            Self::ConicalRibbon => 1.05,
            Self::Ringsail => 1.15,
            Self::FlatCircular => 1.80,
            Self::Cross => 1.20,
            Self::Hemisflo => 1.0,
            Self::DiskGapBand => 1.10,
        }
    }

    /// Typical inflation time in seconds at 1 atm (scales with dynamic pressure)
    pub fn fill_time_nominal(&self) -> f64 {
        match self {
            Self::Drogue => 0.5,
            Self::ConicalRibbon => 1.0,
            Self::Ringsail => 3.0,
            Self::FlatCircular => 2.5,
            Self::Cross => 2.0,
            Self::Hemisflo => 0.3,
            Self::DiskGapBand => 1.5,
        }
    }
}

/// A single parachute stage
#[derive(Debug, Clone)]
pub struct ParachuteStage {
    /// Parachute type
    pub chute_type: ParachuteType,
    /// Number of canopies
    pub num_canopies: u32,
    /// Nominal diameter of each canopy (m)
    pub diameter: f64,
    /// Deployment altitude (m)
    pub deploy_altitude: f64,
    /// Optional: maximum deployment Mach number
    pub max_deploy_mach: Option<f64>,
    /// Optional: maximum deployment dynamic pressure (Pa)
    pub max_deploy_q: Option<f64>,
}

impl ParachuteStage {
    /// Reference area of all canopies combined (m²)
    pub fn reference_area(&self) -> f64 {
        self.num_canopies as f64 * PI / 4.0 * self.diameter * self.diameter
    }

    /// Effective drag area CdS (m²)
    pub fn cds(&self) -> f64 {
        self.chute_type.cd_nominal() * self.reference_area()
    }

    /// Compute terminal velocity at a given altitude (m/s)
    pub fn terminal_velocity(&self, mass: f64, altitude: f64) -> f64 {
        let atm =
            atmosphere::us_standard_atmosphere(altitude).unwrap_or(atmosphere::AtmosphereResult {
                altitude,
                temperature: 0.0,
                pressure: 0.0,
                density: 0.0,
                speed_of_sound: 0.0,
                viscosity: 0.0,
            });
        if atm.density <= 0.0 {
            return f64::INFINITY;
        }
        (2.0 * mass * crate::constants::G0 / (atm.density * self.cds())).sqrt()
    }

    /// Opening shock force (N)
    ///
    /// F_shock = Ck × q × CdS
    /// where q = ½ρv²
    pub fn opening_shock(&self, velocity: f64, altitude: f64) -> f64 {
        let atm =
            atmosphere::us_standard_atmosphere(altitude).unwrap_or(atmosphere::AtmosphereResult {
                altitude,
                temperature: 0.0,
                pressure: 0.0,
                density: 0.0,
                speed_of_sound: 0.0,
                viscosity: 0.0,
            });
        let q = 0.5 * atm.density * velocity * velocity;
        self.chute_type.opening_load_factor() * q * self.cds()
    }

    /// Opening shock in g-loads for a given payload mass
    pub fn opening_shock_g(&self, velocity: f64, altitude: f64, mass: f64) -> f64 {
        self.opening_shock(velocity, altitude) / (mass * crate::constants::G0)
    }
}

/// Multi-stage parachute system
#[derive(Debug, Clone)]
pub struct ParachuteSystem {
    /// Stages in deployment order (first = drogue, last = main)
    pub stages: Vec<ParachuteStage>,
    /// Payload mass including parachute system (kg)
    pub payload_mass: f64,
}

/// Result of a descent simulation
#[derive(Debug, Clone)]
pub struct DescentResult {
    /// Timeline: (time_s, altitude_m, velocity_m_s, deceleration_g, active_stage)
    pub timeline: Vec<(f64, f64, f64, f64, usize)>,
    /// Peak deceleration (g)
    pub peak_decel_g: f64,
    /// Opening shock per stage (g)
    pub opening_shocks_g: Vec<f64>,
    /// Terminal velocity at landing (m/s)
    pub landing_velocity: f64,
    /// Total descent time (s)
    pub descent_time: f64,
}

impl ParachuteSystem {
    /// Simulate descent from initial conditions.
    ///
    /// # Arguments
    /// * `h0` - Initial altitude (m)
    /// * `v0` - Initial velocity magnitude (m/s, downward positive)
    /// * `dt` - Time step (s)
    pub fn simulate_descent(&self, h0: f64, v0: f64, dt: f64) -> DescentResult {
        let mut h = h0;
        let mut v = v0;
        let mut t = 0.0;
        let mut timeline = Vec::new();
        let mut peak_g = 0.0;
        let mut opening_shocks = Vec::new();
        let mut current_stage: Option<usize> = None;
        let g = crate::constants::G0;

        // Sort stages by deployment altitude (highest first)
        let mut stage_indices: Vec<usize> = (0..self.stages.len()).collect();
        stage_indices.sort_by(|a, b| {
            self.stages[*b]
                .deploy_altitude
                .partial_cmp(&self.stages[*a].deploy_altitude)
                .unwrap()
        });

        let mut stage_deployed = vec![false; self.stages.len()];

        while h > 0.0 {
            // Check for stage deployments
            for &idx in &stage_indices {
                if !stage_deployed[idx] && h <= self.stages[idx].deploy_altitude {
                    stage_deployed[idx] = true;
                    current_stage = Some(idx);
                    let shock_g = self.stages[idx].opening_shock_g(v, h, self.payload_mass);
                    opening_shocks.push(shock_g);
                    if shock_g > peak_g {
                        peak_g = shock_g;
                    }
                }
            }

            // Compute drag from all deployed chutes
            let atm =
                atmosphere::us_standard_atmosphere(h).unwrap_or(atmosphere::AtmosphereResult {
                    altitude: h,
                    temperature: 0.0,
                    pressure: 0.0,
                    density: 0.0,
                    speed_of_sound: 0.0,
                    viscosity: 0.0,
                });
            let mut total_cds = 0.0;
            for (idx, stage) in self.stages.iter().enumerate() {
                if stage_deployed[idx] {
                    total_cds += stage.cds();
                }
            }

            let drag_accel = if self.payload_mass > 0.0 && atm.density > 0.0 {
                0.5 * atm.density * v * v * total_cds / self.payload_mass
            } else {
                0.0
            };

            let decel_g = drag_accel / g;
            if decel_g > peak_g {
                peak_g = decel_g;
            }

            let active = current_stage.unwrap_or(0);
            timeline.push((t, h, v, decel_g, active));

            // Update velocity: gravity - drag
            let a = g - drag_accel;
            v += a * dt;
            if v < 0.0 {
                v = 0.0;
            }
            h -= v * dt;
            t += dt;

            if t > 100_000.0 {
                break; // Safety cutoff
            }
        }

        let landing_v = if let Some(last) = timeline.last() {
            last.2
        } else {
            v0
        };

        DescentResult {
            timeline,
            peak_decel_g: peak_g,
            opening_shocks_g: opening_shocks,
            landing_velocity: landing_v,
            descent_time: t,
        }
    }
}

/// Size a parachute for a target landing velocity.
///
/// Returns required canopy diameter (m).
///
/// # Arguments
/// * `mass` - Payload mass (kg)
/// * `v_landing` - Target landing velocity (m/s)
/// * `chute_type` - Parachute type
/// * `num_canopies` - Number of canopies
/// * `altitude` - Landing altitude (m), default 0
pub fn size_for_landing_velocity(
    mass: f64,
    v_landing: f64,
    chute_type: ParachuteType,
    num_canopies: u32,
    altitude: f64,
) -> f64 {
    let atm =
        atmosphere::us_standard_atmosphere(altitude).unwrap_or(atmosphere::AtmosphereResult {
            altitude,
            temperature: 0.0,
            pressure: 0.0,
            density: 0.0,
            speed_of_sound: 0.0,
            viscosity: 0.0,
        });
    let g = crate::constants::G0;
    let cd = chute_type.cd_nominal();

    // v_t = sqrt(2mg / (ρ CdS))  →  S = 2mg / (ρ Cd v²)
    let total_area = 2.0 * mass * g / (atm.density * cd * v_landing * v_landing);
    let area_per_canopy = total_area / num_canopies as f64;
    (4.0 * area_per_canopy / PI).sqrt()
}

/// Apollo-style parachute system (3 main + 2 drogue)
pub fn apollo_parachute_system(capsule_mass: f64) -> ParachuteSystem {
    ParachuteSystem {
        stages: vec![
            ParachuteStage {
                chute_type: ParachuteType::ConicalRibbon,
                num_canopies: 2,
                diameter: 5.2,
                deploy_altitude: 7600.0,
                max_deploy_mach: Some(0.7),
                max_deploy_q: Some(10000.0),
            },
            ParachuteStage {
                chute_type: ParachuteType::Ringsail,
                num_canopies: 3,
                diameter: 25.4,
                deploy_altitude: 3000.0,
                max_deploy_mach: None,
                max_deploy_q: None,
            },
        ],
        payload_mass: capsule_mass,
    }
}

/// Crew Dragon parachute system (4 main + 2 drogue)
pub fn dragon_parachute_system(capsule_mass: f64) -> ParachuteSystem {
    ParachuteSystem {
        stages: vec![
            ParachuteStage {
                chute_type: ParachuteType::ConicalRibbon,
                num_canopies: 2,
                diameter: 3.7,
                deploy_altitude: 5500.0,
                max_deploy_mach: Some(0.6),
                max_deploy_q: None,
            },
            ParachuteStage {
                chute_type: ParachuteType::Ringsail,
                num_canopies: 4,
                diameter: 35.0,
                deploy_altitude: 1800.0,
                max_deploy_mach: None,
                max_deploy_q: None,
            },
        ],
        payload_mass: capsule_mass,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_terminal_velocity() {
        // A 100 kg mass on a 10m diameter flat circular chute at sea level
        let stage = ParachuteStage {
            chute_type: ParachuteType::FlatCircular,
            num_canopies: 1,
            diameter: 10.0,
            deploy_altitude: 5000.0,
            max_deploy_mach: None,
            max_deploy_q: None,
        };

        let vt = stage.terminal_velocity(100.0, 0.0);
        // v_t = sqrt(2 * 100 * 9.81 / (1.225 * 0.78 * π/4 * 100)) ≈ 5.7 m/s
        assert!(vt > 4.0 && vt < 8.0, "v_t = {vt} m/s");
    }

    #[test]
    fn test_opening_shock() {
        let stage = ParachuteStage {
            chute_type: ParachuteType::Ringsail,
            num_canopies: 3,
            diameter: 25.0,
            deploy_altitude: 3000.0,
            max_deploy_mach: None,
            max_deploy_q: None,
        };

        let shock_g = stage.opening_shock_g(80.0, 3000.0, 5000.0);
        // Should be several g's
        assert!(
            shock_g > 1.0 && shock_g < 150.0,
            "Opening shock = {shock_g} g"
        );
    }

    #[test]
    fn test_size_for_landing() {
        // Size a chute for 100 kg payload at 7 m/s landing
        let d = size_for_landing_velocity(100.0, 7.0, ParachuteType::Ringsail, 1, 0.0);
        // Should be a few meters
        assert!(d > 2.0 && d < 15.0, "diameter = {d} m");

        // Verify by computing terminal velocity with sized chute
        let stage = ParachuteStage {
            chute_type: ParachuteType::Ringsail,
            num_canopies: 1,
            diameter: d,
            deploy_altitude: 1000.0,
            max_deploy_mach: None,
            max_deploy_q: None,
        };
        let vt = stage.terminal_velocity(100.0, 0.0);
        assert!((vt - 7.0).abs() < 0.1, "v_t = {vt}, expected 7.0");
    }

    #[test]
    fn test_apollo_system() {
        let system = apollo_parachute_system(5500.0);
        assert_eq!(system.stages.len(), 2);
        assert_eq!(system.stages[0].num_canopies, 2); // 2 drogues
        assert_eq!(system.stages[1].num_canopies, 3); // 3 mains

        // Terminal velocity with mains should be ~8-9 m/s (Apollo was ~8.5 m/s)
        let vt = system.stages[1].terminal_velocity(5500.0, 0.0);
        assert!(vt > 5.0 && vt < 15.0, "Apollo v_t = {vt} m/s");
    }

    #[test]
    fn test_descent_simulation() {
        let system = ParachuteSystem {
            stages: vec![
                ParachuteStage {
                    chute_type: ParachuteType::Drogue,
                    num_canopies: 1,
                    diameter: 3.0,
                    deploy_altitude: 5000.0,
                    max_deploy_mach: None,
                    max_deploy_q: None,
                },
                ParachuteStage {
                    chute_type: ParachuteType::Ringsail,
                    num_canopies: 1,
                    diameter: 15.0,
                    deploy_altitude: 2000.0,
                    max_deploy_mach: None,
                    max_deploy_q: None,
                },
            ],
            payload_mass: 500.0,
        };

        let result = system.simulate_descent(10000.0, 100.0, 0.5);

        assert!(result.descent_time > 0.0);
        assert!(result.landing_velocity < 100.0, "Should decelerate");
        assert!(result.peak_decel_g > 0.0);
        assert!(!result.opening_shocks_g.is_empty());
    }

    #[test]
    fn test_dragon_system() {
        let system = dragon_parachute_system(12000.0);
        assert_eq!(system.stages[1].num_canopies, 4); // 4 mains
    }
}
