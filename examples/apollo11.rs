//! Apollo 11 trajectory recreation
//! 
//! Simulates the complete Apollo 11 mission from launch to lunar landing,
//! including trans-lunar injection, lunar orbit insertion, and landing descent.

use rocketlab::{
    mission::{hohmann_transfer, gravity_assist},
    orbits::{OrbitalElements, StateVector},
    lambert::solve_lambert,
    trajectory::TrajectoryAnalysis,
    kepler::solve_kepler,
    constants::{MU_EARTH, MU_MOON, R_EARTH, R_MOON},
    atmosphere::us_standard_atmosphere,
    landing::{PDGProblem, PoweredDescentGuidance, ConstraintViolations, PDGSolution},
    propulsion::{rocket_equation, propellant_mass},
};
use nalgebra::Vector3;
use chrono::{DateTime, Utc, TimeZone};

/// Apollo Command/Service Module specifications
struct CSMSpecs {
    pub mass_dry: f64,      // kg
    pub mass_fuel: f64,     // kg
    pub isp_sps: f64,       // seconds (Service Propulsion System)
    pub thrust_sps: f64,    // Newtons
}

/// Apollo Lunar Module specifications  
struct LMSpecs {
    pub mass_dry: f64,      // kg (descent + ascent stages)
    pub mass_fuel_descent: f64,  // kg
    pub mass_fuel_ascent: f64,   // kg
    pub isp_descent: f64,   // seconds
    pub isp_ascent: f64,    // seconds
    pub thrust_descent: f64, // Newtons
    pub thrust_ascent: f64,  // Newtons
}

/// Apollo mission trajectory phases
#[derive(Debug, Clone)]
enum MissionPhase {
    Launch,
    EarthParkingOrbit,
    TransLunarInjection,
    TransLunarCoast,
    LunarOrbitInsertion,
    LunarOrbit,
    LunarDescentInitiation,
    PoweredDescent,
    LunarSurface,
    LunarAscentInitiation,
    PoweredAscent,
    LunarOrbitRendezvous,
    TransEarthInjection,
    TransEarthCoast,
    EarthReentry,
}

/// Apollo 11 mission simulator
struct Apollo11Mission {
    csm: CSMSpecs,
    lm: LMSpecs,
    launch_date: DateTime<Utc>,
    current_phase: MissionPhase,
}

impl Apollo11Mission {
    fn new() -> Self {
        // Historical Apollo specifications
        let csm = CSMSpecs {
            mass_dry: 11_700.0,     // Command + Service Module dry mass
            mass_fuel: 18_400.0,    // SPS fuel
            isp_sps: 314.0,         // SPS specific impulse
            thrust_sps: 97_400.0,   // SPS thrust
        };

        let lm = LMSpecs {
            mass_dry: 4_700.0,      // LM dry mass (both stages)
            mass_fuel_descent: 8_200.0,  // Descent stage fuel
            mass_fuel_ascent: 2_350.0,   // Ascent stage fuel  
            isp_descent: 311.0,     // Descent engine Isp
            isp_ascent: 311.0,      // Ascent engine Isp
            thrust_descent: 45_040.0,    // Descent engine max thrust
            thrust_ascent: 15_600.0,     // Ascent engine thrust
        };

        // Apollo 11 launched July 16, 1969, 13:32 UTC
        let launch_date = Utc.ymd(1969, 7, 16).and_hms(13, 32, 0);

        Self {
            csm,
            lm,
            launch_date,
            current_phase: MissionPhase::Launch,
        }
    }

    /// Simulate launch to Earth parking orbit
    fn simulate_launch(&mut self) -> StateVector {
        println!("🚀 Apollo 11 Launch - {}", self.launch_date);
        println!("Saturn V ignition at Kennedy Space Center LC-39A");
        
        // Saturn V performance to LEO
        let leo_altitude = 185.0; // km (actual Apollo 11 parking orbit)
        let leo_radius = R_EARTH + leo_altitude;
        
        // Circular parking orbit
        let velocity = (MU_EARTH / leo_radius).sqrt();
        let position = Vector3::new(leo_radius, 0.0, 0.0);
        let velocity_vec = Vector3::new(0.0, velocity, 0.0);
        
        println!("✅ Earth parking orbit achieved:");
        println!("   Altitude: {:.1} km", leo_altitude);
        println!("   Velocity: {:.1} km/s", velocity);
        
        self.current_phase = MissionPhase::EarthParkingOrbit;
        StateVector { position, velocity: velocity_vec }
    }

    /// Simulate Trans-Lunar Injection (TLI)
    fn simulate_tli(&mut self, parking_state: &StateVector) -> StateVector {
        println!("\n🌙 Trans-Lunar Injection (TLI)");
        println!("Third stage S-IVB restart for lunar trajectory");
        
        // TLI delta-V (historical value)
        let tli_delta_v = 3.150; // km/s
        
        // Add velocity in prograde direction
        let velocity_unit = parking_state.velocity.normalize();
        let new_velocity = parking_state.velocity + velocity_unit * tli_delta_v;
        
        println!("✅ TLI burn complete:");
        println!("   Delta-V: {:.3} km/s", tli_delta_v);
        println!("   Velocity after TLI: {:.3} km/s", new_velocity.magnitude());
        
        // Check if trajectory is hyperbolic (escaping Earth)
        let specific_energy = new_velocity.magnitude_squared() / 2.0 - MU_EARTH / parking_state.position.magnitude();
        if specific_energy > 0.0 {
            println!("   ✅ Trajectory is hyperbolic - escaping Earth gravity");
        }
        
        self.current_phase = MissionPhase::TransLunarInjection;
        StateVector { 
            position: parking_state.position, 
            velocity: new_velocity 
        }
    }

    /// Simulate trans-lunar coast and lunar approach
    fn simulate_lunar_approach(&mut self) -> StateVector {
        println!("\n🌙 Trans-Lunar Coast (3.25 days)");
        println!("Coasting through cislunar space...");
        
        // Simplified: place spacecraft in lunar sphere of influence
        // Real calculation would integrate the trajectory
        
        // Historical Apollo 11 lunar approach
        let lunar_distance = 384_400.0; // km (Earth-Moon distance)
        let approach_altitude = 100.0;   // km above lunar surface
        let lunar_orbit_radius = R_MOON + approach_altitude;
        
        // Position at lunar arrival
        let position = Vector3::new(lunar_orbit_radius, 0.0, 0.0);
        
        // Hyperbolic approach velocity (approximation)
        let v_infinity = 0.9; // km/s (excess velocity relative to Moon)
        let approach_velocity = (v_infinity.powi(2) + 2.0 * MU_MOON / lunar_orbit_radius).sqrt();
        let velocity = Vector3::new(0.0, approach_velocity, 0.0);
        
        println!("✅ Lunar approach:");
        println!("   Distance from Moon: {:.1} km", lunar_orbit_radius - R_MOON);
        println!("   Approach velocity: {:.3} km/s", approach_velocity);
        
        self.current_phase = MissionPhase::TransLunarCoast;
        StateVector { position, velocity }
    }

    /// Simulate Lunar Orbit Insertion (LOI)
    fn simulate_lunar_orbit_insertion(&mut self, approach_state: &StateVector) -> StateVector {
        println!("\n🛰️ Lunar Orbit Insertion (LOI)");
        println!("SPS burn to capture into lunar orbit");
        
        // Target circular orbit at 100 km altitude
        let orbit_radius = R_MOON + 100.0;
        let circular_velocity = (MU_MOON / orbit_radius).sqrt();
        
        // Calculate required delta-V (retrograde burn)
        let current_speed = approach_state.velocity.magnitude();
        let delta_v_required = current_speed - circular_velocity;
        
        // Check CSM fuel availability
        let fuel_mass_needed = propellant_mass(
            self.csm.mass_dry + self.csm.mass_fuel, 
            delta_v_required, 
            self.csm.isp_sps
        );
        
        println!("✅ LOI burn analysis:");
        println!("   Delta-V required: {:.3} km/s", delta_v_required);
        println!("   Fuel consumed: {:.0} kg", fuel_mass_needed);
        println!("   Fuel remaining: {:.0} kg", self.csm.mass_fuel - fuel_mass_needed);
        
        // Update fuel mass
        self.csm.mass_fuel -= fuel_mass_needed;
        
        // New orbital state
        let position = Vector3::new(orbit_radius, 0.0, 0.0);
        let velocity = Vector3::new(0.0, circular_velocity, 0.0);
        
        self.current_phase = MissionPhase::LunarOrbitInsertion;
        StateVector { position, velocity }
    }

    /// Simulate lunar landing using powered descent guidance
    fn simulate_lunar_landing(&mut self) -> Result<PDGSolution, String> {
        println!("\n🌕 Lunar Module Powered Descent");
        println!("Descent from 100 km to lunar surface");
        
        // Landing site: Sea of Tranquility (Apollo 11 actual site)
        let landing_site = Vector3::new(0.0, 0.0, 0.0); // Simplified: at lunar center reference
        
        // Initial conditions for powered descent
        let initial_altitude = 100.0; // km
        let initial_position = Vector3::new(0.0, 0.0, R_MOON + initial_altitude);
        let initial_velocity = Vector3::new(-1.5, 0.0, -0.1); // km/s (tangential + descent)
        
        // Landing constraints and vehicle parameters
        let problem = PDGProblem {
            // Initial state
            r0: initial_position,
            v0: initial_velocity,
            m0: self.lm.mass_dry + self.lm.mass_fuel_descent,
            
            // Final state (soft landing)
            rf: Vector3::new(0.0, 0.0, R_MOON + 0.01), // Just above surface
            vf: Vector3::new(0.0, 0.0, -0.01), // Nearly zero velocity
            
            // Vehicle constraints
            thrust_max: self.lm.thrust_descent,
            thrust_min: 0.1 * self.lm.thrust_descent, // 10% throttling
            isp: self.lm.isp_descent,
            
            // Landing constraints
            glideslope_angle: 20.0_f64.to_radians(), // degrees to radians
            max_tilt_angle: 30.0_f64.to_radians(),
            
            // Environmental
            gravity: Vector3::new(0.0, 0.0, -MU_MOON / (R_MOON * R_MOON)) * 1000.0, // m/s²
            
            // Flight time estimate
            tf_guess: 750.0, // seconds (12.5 minutes - Apollo typical)
            
            // Solver parameters
            max_iterations: 50,
            tolerance: 1e-6,
        };
        
        println!("Landing constraints:");
        println!("   Initial altitude: {:.1} km", initial_altitude);
        println!("   LM mass: {:.0} kg", problem.m0);
        println!("   Max thrust: {:.0} N", problem.thrust_max);
        println!("   Flight time estimate: {:.0} s", problem.tf_guess);
        
        // Solve the powered descent guidance problem
        let mut pdg = PoweredDescentGuidance::new(problem);
        let solution = pdg.solve()?;
        
        // Validate constraints
        let violations = ConstraintViolations::validate(&pdg.problem, &solution);
        
        println!("✅ Lunar landing solution:");
        println!("   Flight time: {:.1} s", solution.time[solution.time.len()-1]);
        println!("   Fuel consumed: {:.0} kg", pdg.problem.m0 - solution.mass[solution.mass.len()-1]);
        println!("   Final velocity: {:.3} km/s", solution.velocity[solution.velocity.len()-1].magnitude() / 1000.0);
        
        if violations.has_violations() {
            println!("⚠️  Constraint violations detected");
            if !violations.thrust_max_violations.is_empty() {
                println!("   - Thrust limit exceeded at {} points", violations.thrust_max_violations.len());
            }
            if !violations.glideslope_violations.is_empty() {
                println!("   - Glideslope violated at {} points", violations.glideslope_violations.len());
            }
        } else {
            println!("✅ All constraints satisfied");
        }
        
        self.current_phase = MissionPhase::PoweredDescent;
        Ok(solution)
    }

    /// Print mission summary
    fn print_mission_summary(&self) {
        println!("\n📊 Apollo 11 Mission Summary");
        println!("==========================================");
        println!("Launch Date: {}", self.launch_date);
        println!("Current Phase: {:?}", self.current_phase);
        println!("\nVehicle Status:");
        println!("CSM fuel remaining: {:.0} kg", self.csm.mass_fuel);
        println!("LM descent fuel: {:.0} kg", self.lm.mass_fuel_descent);
        println!("LM ascent fuel: {:.0} kg", self.lm.mass_fuel_ascent);
        
        // Historical Apollo 11 facts
        println!("\n🏆 Historical Apollo 11 Facts:");
        println!("• Launch: July 16, 1969, 13:32 UTC");
        println!("• Lunar landing: July 20, 1969, 20:17 UTC");
        println!("• Landing site: Sea of Tranquility (0°40'27\"N 23°28'23\"E)");
        println!("• Surface EVA time: 21 hours 36 minutes");
        println!("• Lunar surface samples: 21.5 kg");
        println!("• Total mission time: 8 days 3 hours 18 minutes");
        println!("• Splashdown: July 24, 1969, Pacific Ocean");
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Apollo 11 Mission Simulation");
    println!("============================\n");
    
    // Initialize mission
    let mut apollo11 = Apollo11Mission::new();
    
    // Phase 1: Launch to LEO
    let parking_orbit = apollo11.simulate_launch();
    
    // Phase 2: Trans-Lunar Injection
    let tli_state = apollo11.simulate_tli(&parking_orbit);
    
    // Phase 3: Lunar approach
    let lunar_approach = apollo11.simulate_lunar_approach();
    
    // Phase 4: Lunar orbit insertion
    let lunar_orbit = apollo11.simulate_lunar_orbit_insertion(&lunar_approach);
    
    // Phase 5: Lunar landing
    match apollo11.simulate_lunar_landing() {
        Ok(_landing_solution) => {
            println!("🎉 Eagle has landed! Houston, Tranquility Base here.");
        }
        Err(e) => {
            println!("❌ Landing guidance failed: {}", e);
            println!("Houston, we have a problem with the landing sequence.");
        }
    }
    
    // Mission summary
    apollo11.print_mission_summary();
    
    Ok(())
}