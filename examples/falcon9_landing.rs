//! Falcon 9 boostback and landing simulation
//!
//! Simulates a complete Falcon 9 first stage mission including:
//! - Ascent and stage separation  
//! - Boostback burn for return trajectory
//! - Entry burn for aerodynamic control
//! - Landing burn with supersonic retropropulsion

use rocketlab::{
    landing::{PDGProblem, PoweredDescentGuidance, ConstraintViolations, PDGSolution},
    trajectory::{TrajectoryAnalysis, TrajectoryOptimizer, IntegrationMethod},
    atmosphere::{us_standard_atmosphere, dynamic_pressure, mach_number},
    engines::EngineDatabase,
    propulsion::{rocket_equation, propellant_mass},
    reentry::{ReentrySimulator, ReentryVehicle, HeatingModel, EntryProfile},
    constants::{MU_EARTH, R_EARTH, G0},
};
use nalgebra::Vector3;
use chrono::{DateTime, Utc, Duration, TimeZone};

/// Falcon 9 first stage specifications
struct Falcon9FirstStage {
    pub mass_dry: f64,          // kg
    pub mass_fuel: f64,         // kg (LOX/RP-1)
    pub num_engines: u8,        // Number of Merlin 1D engines
    pub engine_thrust_sl: f64,  // N (sea level)
    pub engine_thrust_vac: f64, // N (vacuum)
    pub engine_isp_sl: f64,     // s (sea level)
    pub engine_isp_vac: f64,    // s (vacuum)
    pub diameter: f64,          // m
    pub length: f64,            // m
    pub grid_fins: bool,        // Grid fins for aerodynamic control
    pub legs_deployed: bool,    // Landing legs status
}

/// Falcon 9 mission profile phases
#[derive(Debug, Clone)]
enum FlightPhase {
    Liftoff,
    MaxQ,
    MECO,           // Main Engine Cut Off
    StageSequence,  // Stage separation sequence 
    Boostback,      // Boostback burn
    CoastPhase,     // Ballistic coast
    EntryBurn,      // Entry burn
    Aerodynamic,    // Aerodynamic descent with grid fins
    LandingBurn,    // Final landing burn
    TouchdownTargeted,
}

/// Landing site options
#[derive(Debug, Clone)]
enum LandingSite {
    ASDS {          // Autonomous Spaceport Drone Ship
        name: String,
        position: Vector3<f64>, // Position in ECEF coordinates
        distance_km: f64,       // Distance from launch site
    },
    RTLS {          // Return to Launch Site
        name: String,
        position: Vector3<f64>,
    },
}

/// Falcon 9 mission simulator
struct Falcon9Mission {
    stage1: Falcon9FirstStage,
    launch_time: DateTime<Utc>,
    current_phase: FlightPhase,
    landing_site: LandingSite,
    trajectory_data: Vec<(f64, Vector3<f64>, Vector3<f64>)>, // time, position, velocity
}

impl Falcon9Mission {
    fn new_rtls_mission() -> Self {
        // Falcon 9 Block 5 specifications (conservative estimates)
        let stage1 = Falcon9FirstStage {
            mass_dry: 25_600.0,       // kg (with landing hardware)
            mass_fuel: 395_700.0,     // kg (LOX: ~276,000 kg, RP-1: ~119,700 kg)
            num_engines: 9,
            engine_thrust_sl: 845_000.0,   // N per Merlin 1D (total ~7.6 MN)
            engine_thrust_vac: 934_000.0,  // N per engine in vacuum
            engine_isp_sl: 282.0,     // seconds
            engine_isp_vac: 311.0,    // seconds
            diameter: 3.7,            // meters
            length: 47.0,             // meters (first stage)
            grid_fins: true,
            legs_deployed: false,
        };

        let landing_site = LandingSite::RTLS {
            name: "LZ-1 (Landing Zone 1)".to_string(),
            position: Vector3::new(0.0, 0.0, 0.0), // Simplified coordinates
        };

        // Typical launch time
        let launch_time = Utc.ymd(2024, 1, 15).and_hms(14, 30, 0);

        Self {
            stage1,
            launch_time,
            current_phase: FlightPhase::Liftoff,
            landing_site,
            trajectory_data: Vec::new(),
        }
    }

    fn new_asds_mission() -> Self {
        let mut mission = Self::new_rtls_mission();
        mission.landing_site = LandingSite::ASDS {
            name: "Of Course I Still Love You".to_string(),
            position: Vector3::new(600_000.0, 0.0, 0.0), // ~600 km downrange
            distance_km: 600.0,
        };
        mission
    }

    /// Simulate ascent to MECO
    fn simulate_ascent(&mut self) -> Vector3<f64> {
        println!("🚀 Falcon 9 Liftoff - {}", self.launch_time);
        println!("9 Merlin 1D engines ignited");
        
        // Simplified ascent profile
        let meco_time = 162.0; // seconds (typical for Falcon 9)
        let meco_altitude = 80.0; // km
        let meco_velocity = 2.8; // km/s
        
        println!("📈 Ascent profile:");
        println!("   Max-Q at T+{:.0}s", 70.0);
        println!("   MECO at T+{:.0}s", meco_time);
        println!("   MECO altitude: {:.1} km", meco_altitude);
        println!("   MECO velocity: {:.1} km/s", meco_velocity);
        
        // Calculate fuel consumption during ascent
        let ascent_fuel_fraction = 0.85; // ~85% of fuel used for ascent
        let fuel_consumed = self.stage1.mass_fuel * ascent_fuel_fraction;
        self.stage1.mass_fuel -= fuel_consumed;
        
        println!("   Fuel consumed: {:.0} kg", fuel_consumed);
        println!("   Fuel remaining: {:.0} kg", self.stage1.mass_fuel);
        
        self.current_phase = FlightPhase::MECO;
        
        // MECO state vector (simplified)
        let position = Vector3::new(R_EARTH + meco_altitude, 0.0, 0.0);
        let velocity = Vector3::new(0.5 * meco_velocity, meco_velocity * 0.866, 0.0); // Angled trajectory
        
        self.trajectory_data.push((meco_time, position, velocity));
        velocity
    }

    /// Simulate boostback burn
    fn simulate_boostback(&mut self, meco_velocity: Vector3<f64>) -> Vector3<f64> {
        println!("\n🔄 Boostback Burn");
        
        match &self.landing_site {
            LandingSite::RTLS { name, .. } => {
                println!("Target: {} (Return to Launch Site)", name);
                
                // RTLS requires significant velocity reversal
                let boostback_delta_v = 1.8; // km/s (typical for RTLS)
                let burn_duration = 40.0;     // seconds
                
                // Calculate required fuel
                let current_mass = self.stage1.mass_dry + self.stage1.mass_fuel;
                let fuel_needed = propellant_mass(current_mass, boostback_delta_v, self.stage1.engine_isp_vac);
                
                if fuel_needed > self.stage1.mass_fuel {
                    panic!("❌ Insufficient fuel for boostback burn! Need {:.0} kg, have {:.0} kg", 
                           fuel_needed, self.stage1.mass_fuel);
                }
                
                self.stage1.mass_fuel -= fuel_needed;
                
                println!("✅ Boostback burn complete:");
                println!("   Delta-V: {:.1} km/s", boostback_delta_v);
                println!("   Burn duration: {:.0}s", burn_duration);
                println!("   Fuel consumed: {:.0} kg", fuel_needed);
                println!("   Fuel remaining: {:.0} kg", self.stage1.mass_fuel);
                
                // New velocity after boostback (simplified - mostly reverses horizontal component)
                Vector3::new(-0.8 * meco_velocity.x, -0.5 * meco_velocity.y, meco_velocity.z)
            }
            
            LandingSite::ASDS { name, distance_km, .. } => {
                println!("Target: {} ({:.0} km downrange)", name, distance_km);
                
                // ASDS requires less delta-V (no full reversal)
                let boostback_delta_v = 0.6; // km/s
                let burn_duration = 20.0;     // seconds
                
                let current_mass = self.stage1.mass_dry + self.stage1.mass_fuel;
                let fuel_needed = propellant_mass(current_mass, boostback_delta_v, self.stage1.engine_isp_vac);
                
                self.stage1.mass_fuel -= fuel_needed;
                
                println!("✅ Boostback burn complete:");
                println!("   Delta-V: {:.1} km/s", boostback_delta_v);
                println!("   Fuel consumed: {:.0} kg", fuel_needed);
                
                // Adjust trajectory for downrange landing
                Vector3::new(0.2 * meco_velocity.x, meco_velocity.y, meco_velocity.z)
            }
        }
    }

    /// Simulate atmospheric entry with grid fins
    fn simulate_entry(&mut self) -> Vector3<f64> {
        println!("\n🌊 Atmospheric Entry");
        println!("Grid fins deployed for aerodynamic control");
        
        // Entry interface at ~70 km
        let entry_altitude = 70.0; // km
        let entry_velocity = 2.5;  // km/s (typical re-entry speed)
        
        // Create reentry vehicle model for Falcon 9
        let vehicle = ReentryVehicle {
            mass: self.stage1.mass_dry + self.stage1.mass_fuel,
            reference_area: std::f64::consts::PI * (self.stage1.diameter / 2.0).powi(2), // m²
            drag_coefficient: 0.6, // With grid fins extended
            ballistic_coefficient: 0.0, // Will be calculated
            lift_to_drag_ratio: 0.1, // Grid fins provide some lift
            nose_radius: 1.0,      // m
        };
        
        // Entry profile
        let profile = EntryProfile {
            entry_velocity,
            entry_altitude,
            entry_gamma: -45.0,    // degrees (steep entry)
            target_altitude: 15.0, // km (where landing burn starts)
            max_time: 300.0,       // seconds
        };
        
        println!("Entry conditions:");
        println!("   Entry altitude: {:.1} km", entry_altitude);
        println!("   Entry velocity: {:.1} km/s", entry_velocity);
        println!("   Ballistic coefficient: {:.0} kg/m²", vehicle.mass / vehicle.reference_area);
        
        // Simplified atmospheric descent simulation
        let mut altitude = entry_altitude;
        let mut velocity = entry_velocity;
        let mut time = 0.0;
        
        while altitude > 15.0 && time < 300.0 {
            // Get atmospheric conditions
            if let Some(atm) = us_standard_atmosphere(altitude * 1000.0) {
                // Calculate aerodynamic forces
                let q = dynamic_pressure(velocity * 1000.0, atm.density); // Pa
                let mach = mach_number(velocity * 1000.0, altitude * 1000.0).unwrap_or(0.0);
                
                // Simplified dynamics (Euler integration)
                let drag_accel = q * vehicle.reference_area * vehicle.drag_coefficient / vehicle.mass;
                let gravity = MU_EARTH / (R_EARTH + altitude).powi(2);
                
                velocity -= (drag_accel / 1000.0 + gravity / 1000.0) * 1.0; // 1-second time step
                altitude -= velocity * 1.0;
                time += 1.0;
                
                // Print key milestones
                if time as i32 % 30 == 0 {
                    println!("   T+{:.0}s: Alt {:.1} km, Vel {:.1} km/s, Mach {:.1}, Q {:.0} Pa", 
                             time, altitude, velocity, mach, q);
                }
            }
        }
        
        self.current_phase = FlightPhase::Aerodynamic;
        println!("✅ Atmospheric descent complete");
        println!("   Final altitude: {:.1} km", altitude);
        println!("   Final velocity: {:.1} km/s", velocity);
        
        Vector3::new(0.0, 0.0, -velocity) // Mostly vertical descent
    }

    /// Simulate landing burn using supersonic retropropulsion  
    fn simulate_landing_burn(&mut self, entry_velocity: Vector3<f64>) -> Result<PDGSolution, String> {
        println!("\n🔥 Landing Burn - Supersonic Retropropulsion");
        println!("Single Merlin 1D engine restart");
        
        // Landing burn typically starts around 15 km altitude
        let initial_altitude = 15.0; // km
        let initial_position = Vector3::new(0.0, 0.0, R_EARTH + initial_altitude);
        
        // Landing constraints for Falcon 9
        let problem = PDGProblem {
            // Initial state
            r0: initial_position,
            v0: entry_velocity,
            m0: self.stage1.mass_dry + self.stage1.mass_fuel,
            
            // Final state (soft touchdown)
            rf: Vector3::new(0.0, 0.0, R_EARTH + 0.001), // Just above surface
            vf: Vector3::new(0.0, 0.0, -0.01), // Gentle touchdown
            
            // Engine constraints (single Merlin 1D)
            thrust_max: self.stage1.engine_thrust_sl, // Single engine
            thrust_min: 0.4 * self.stage1.engine_thrust_sl, // 40% throttle minimum
            isp: self.stage1.engine_isp_sl, // Sea level performance
            
            // Landing constraints
            glideslope_angle: 0.0, // Vertical landing
            max_tilt_angle: 5.0_f64.to_radians(), // Nearly vertical
            
            // Environmental
            gravity: Vector3::new(0.0, 0.0, -G0), // Earth gravity
            
            // Flight time (landing burn is fast!)
            tf_guess: 30.0, // seconds
            
            // Solver parameters
            max_iterations: 100,
            tolerance: 1e-6,
        };
        
        println!("Landing burn parameters:");
        println!("   Initial altitude: {:.1} km", initial_altitude);
        println!("   Entry velocity: {:.2} km/s", entry_velocity.magnitude());
        println!("   Vehicle mass: {:.0} kg", problem.m0);
        println!("   Engine thrust: {:.0} kN", problem.thrust_max / 1000.0);
        println!("   Min throttle: {:.0}%", 40.0);
        
        // Solve the landing guidance problem
        let mut pdg = PoweredDescentGuidance::new(problem);
        let solution = pdg.solve()?;
        
        // Validate the solution
        let violations = ConstraintViolations::validate(&pdg.problem, &solution);
        
        let final_time = solution.time[solution.time.len() - 1];
        let fuel_used = pdg.problem.m0 - solution.mass[solution.mass.len() - 1];
        let final_velocity = solution.velocity[solution.velocity.len() - 1].magnitude();
        
        println!("✅ Landing burn solution:");
        println!("   Burn duration: {:.1} s", final_time);
        println!("   Fuel consumed: {:.0} kg", fuel_used);
        println!("   Final velocity: {:.2} m/s", final_velocity);
        
        // Check if we have enough fuel
        if fuel_used > self.stage1.mass_fuel {
            println!("❌ Insufficient fuel! Need {:.0} kg, have {:.0} kg", fuel_used, self.stage1.mass_fuel);
            return Err("Fuel shortage during landing burn".to_string());
        }
        
        self.stage1.mass_fuel -= fuel_used;
        
        if violations.has_violations() {
            println!("⚠️  Landing constraints violated:");
            if !violations.thrust_max_violations.is_empty() {
                println!("   - Engine thrust exceeded at {} points", violations.thrust_max_violations.len());
            }
            if !violations.thrust_min_violations.is_empty() {
                println!("   - Minimum throttle violated at {} points", violations.thrust_min_violations.len());
            }
        } else {
            println!("✅ All landing constraints satisfied");
        }
        
        self.current_phase = FlightPhase::TouchdownTargeted;
        Ok(solution)
    }

    /// Print mission summary
    fn print_mission_summary(&self) {
        println!("\n🎯 Falcon 9 Mission Summary");
        println!("=========================");
        println!("Launch Time: {}", self.launch_time);
        println!("Mission Type: {:?}", self.landing_site);
        println!("Current Phase: {:?}", self.current_phase);
        
        println!("\nFirst Stage Status:");
        println!("Dry mass: {:.0} kg", self.stage1.mass_dry);
        println!("Fuel remaining: {:.0} kg", self.stage1.mass_fuel);
        
        let fuel_fraction_used = 1.0 - (self.stage1.mass_fuel / 395_700.0);
        println!("Fuel utilization: {:.1}%", fuel_fraction_used * 100.0);
        
        println!("\n🚀 SpaceX Reusability Facts:");
        println!("• First successful landing: December 2015 (Falcon 9 Flight 20)");
        println!("• Landing accuracy: Typically within 10 meters of target");
        println!("• Reuse record: Some boosters have flown 15+ times");
        println!("• Cost savings: ~30% reduction in launch costs");
        println!("• Landing platforms: Land-based LZ-1, LZ-2 and drone ships ASDS");
        
        match self.current_phase {
            FlightPhase::TouchdownTargeted => {
                println!("\n🎉 The Falcon has landed! Another successful recovery.");
                println!("Stage ready for refurbishment and reuse.");
            }
            _ => {
                println!("\n⏳ Mission in progress...");
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Falcon 9 First Stage Recovery Simulation");
    println!("========================================\n");
    
    // Choose mission type
    println!("Select mission type:");
    println!("1. RTLS (Return to Launch Site) - more challenging");
    println!("2. ASDS (Autonomous Spaceport Drone Ship) - easier");
    
    // For this example, simulate an RTLS mission
    println!("\nSimulating RTLS mission...\n");
    
    let mut falcon9 = Falcon9Mission::new_rtls_mission();
    
    // Phase 1: Ascent to MECO
    let meco_velocity = falcon9.simulate_ascent();
    
    // Phase 2: Boostback burn  
    let post_boostback_velocity = falcon9.simulate_boostback(meco_velocity);
    
    // Phase 3: Atmospheric entry
    let entry_velocity = falcon9.simulate_entry();
    
    // Phase 4: Landing burn
    match falcon9.simulate_landing_burn(entry_velocity) {
        Ok(_solution) => {
            println!("🎯 Successful touchdown at Landing Zone 1!");
        }
        Err(e) => {
            println!("❌ Landing failed: {}", e);
            println!("Booster lost in the Atlantic. Back to the drawing board!");
        }
    }
    
    // Mission summary
    falcon9.print_mission_summary();
    
    Ok(())
}