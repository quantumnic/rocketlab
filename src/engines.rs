//! Comprehensive rocket engine database with real specifications
//! 
//! Contains actual performance data for engines used by:
//! - SpaceX (Merlin, Raptor, SuperDraco)
//! - NASA (RS-25, RL-10, J-2)
//! - Russia (RD-180, RD-191, NK-33)  
//! - ESA (Vulcain 2, Vinci)
//! - Blue Origin (BE-3, BE-4)
//! - ULA (various heritage engines)
//!
//! Data sourced from official specifications, NASA databases, and verified technical literature.

use std::collections::HashMap;
use serde::{Deserialize, Serialize};

/// Complete rocket engine specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RocketEngine {
    /// Engine name
    pub name: String,
    /// Manufacturer/operator
    pub manufacturer: String,
    /// Engine family/series
    pub family: String,
    /// Propellant combination
    pub propellant: PropellantType,
    /// Engine cycle type
    pub cycle: EngineType,
    /// Thrust at sea level (kN), None if not applicable
    pub thrust_sl: Option<f64>,
    /// Thrust in vacuum (kN)
    pub thrust_vac: f64,
    /// Specific impulse at sea level (s), None if not applicable
    pub isp_sl: Option<f64>,
    /// Specific impulse in vacuum (s)
    pub isp_vac: f64,
    /// Dry mass (kg)
    pub mass: f64,
    /// Chamber pressure (bar)
    pub chamber_pressure: f64,
    /// Mixture ratio (oxidizer/fuel by mass)
    pub mixture_ratio: f64,
    /// Thrust-to-weight ratio (vacuum)
    pub thrust_to_weight: f64,
    /// Engine dimensions
    pub dimensions: EngineDimensions,
    /// First flight year
    pub first_flight: u32,
    /// Current operational status
    pub status: EngineStatus,
    /// Notable features
    pub features: Vec<String>,
    /// Engine variants (different configurations)
    pub variants: Vec<String>,
    /// Applications (rockets that use this engine)
    pub applications: Vec<String>,
}

/// Engine dimensions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineDimensions {
    /// Length (m)
    pub length: f64,
    /// Diameter (m)  
    pub diameter: f64,
    /// Nozzle exit diameter (m)
    pub exit_diameter: f64,
    /// Expansion ratio
    pub expansion_ratio: f64,
}

/// Propellant combination
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum PropellantType {
    /// Liquid Oxygen + Rocket Propellant 1 (refined kerosene)
    LoxRp1,
    /// Liquid Oxygen + Liquid Methane
    LoxCh4,
    /// Liquid Oxygen + Liquid Hydrogen
    LoxLh2,
    /// Nitrogen Tetroxide + Unsymmetrical Dimethylhydrazine
    NtoUdmh,
    /// Nitrogen Tetroxide + Monomethylhydrazine
    NtoMmh,
    /// Nitrogen Tetroxide + Aerozine 50
    NtoAz50,
    /// Solid propellant
    Solid,
    /// Hypergolic (generic)
    Hypergolic,
    /// Monopropellant (hydrazine)
    Monoprop,
}

/// Engine cycle type
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum EngineType {
    /// Gas generator cycle
    GasGenerator,
    /// Staged combustion cycle (fuel rich)
    StagedCombustionFuelRich,
    /// Staged combustion cycle (oxygen rich)
    StagedCombustionOxRich,
    /// Full flow staged combustion
    FullFlow,
    /// Expander cycle
    Expander,
    /// Pressure fed
    PressureFed,
    /// Electric propulsion
    Electric,
    /// Solid motor
    Solid,
    /// Hybrid
    Hybrid,
}

/// Engine operational status
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum EngineStatus {
    /// Currently operational
    Operational,
    /// In development
    Development,
    /// Retired/legacy
    Retired,
    /// Under test
    Testing,
    /// Cancelled
    Cancelled,
}

/// Comprehensive engine database
pub struct EngineDatabase {
    engines: HashMap<String, RocketEngine>,
    by_manufacturer: HashMap<String, Vec<String>>,
    by_propellant: HashMap<PropellantType, Vec<String>>,
    by_thrust_class: HashMap<ThrustClass, Vec<String>>,
}

/// Thrust classification
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ThrustClass {
    /// < 10 kN (attitude control, small satellites)
    Micro,
    /// 10-100 kN (upper stages, small boosters)  
    Small,
    /// 100-1000 kN (medium boosters, main engines)
    Medium,
    /// 1000-10000 kN (large boosters, heavy lift)
    Large,
    /// > 10000 kN (super heavy lift)
    SuperHeavy,
}

impl EngineDatabase {
    /// Create new database with all engines loaded
    pub fn new() -> Self {
        let mut database = Self {
            engines: HashMap::new(),
            by_manufacturer: HashMap::new(),
            by_propellant: HashMap::new(),
            by_thrust_class: HashMap::new(),
        };
        
        database.load_all_engines();
        database.build_indices();
        database
    }
    
    /// Load all engine specifications
    fn load_all_engines(&mut self) {
        // SpaceX Engines
        self.add_engine(merlin_1d_sea_level());
        self.add_engine(merlin_1d_vacuum());
        self.add_engine(merlin_1d_plus());
        self.add_engine(raptor_1());
        self.add_engine(raptor_2());
        self.add_engine(raptor_3());
        self.add_engine(super_draco());
        self.add_engine(draco());
        
        // NASA Engines
        self.add_engine(rs25());
        self.add_engine(rl10b2());
        self.add_engine(j2_engine());
        self.add_engine(rs68());
        
        // Russian Engines
        self.add_engine(rd180());
        self.add_engine(rd191());
        self.add_engine(nk33());
        self.add_engine(rd170());
        self.add_engine(rd253());
        
        // European (ESA) Engines
        self.add_engine(vulcain_2());
        self.add_engine(vinci());
        self.add_engine(aestus());
        
        // Blue Origin
        self.add_engine(be3());
        self.add_engine(be4());
        
        // Other Notable Engines
        self.add_engine(f1_engine());
        self.add_engine(h1_engine());
        self.add_engine(ssme_rs25());
    }
    
    fn add_engine(&mut self, engine: RocketEngine) {
        self.engines.insert(engine.name.clone(), engine);
    }
    
    /// Build search indices
    fn build_indices(&mut self) {
        for (name, engine) in &self.engines {
            // By manufacturer
            self.by_manufacturer
                .entry(engine.manufacturer.clone())
                .or_default()
                .push(name.clone());
            
            // By propellant
            self.by_propellant
                .entry(engine.propellant.clone())
                .or_default()
                .push(name.clone());
            
            // By thrust class
            let thrust_class = classify_thrust(engine.thrust_vac);
            self.by_thrust_class
                .entry(thrust_class)
                .or_default()
                .push(name.clone());
        }
    }
    
    /// Get engine by name
    pub fn get_engine(&self, name: &str) -> Option<&RocketEngine> {
        self.engines.get(name)
    }
    
    /// Get all engines by manufacturer
    pub fn get_by_manufacturer(&self, manufacturer: &str) -> Vec<&RocketEngine> {
        self.by_manufacturer
            .get(manufacturer)
            .unwrap_or(&vec![])
            .iter()
            .filter_map(|name| self.engines.get(name))
            .collect()
    }
    
    /// Get all engines by propellant type
    pub fn get_by_propellant(&self, propellant: PropellantType) -> Vec<&RocketEngine> {
        self.by_propellant
            .get(&propellant)
            .unwrap_or(&vec![])
            .iter()
            .filter_map(|name| self.engines.get(name))
            .collect()
    }
    
    /// Get engines in thrust range (kN)
    pub fn get_by_thrust_range(&self, min_thrust: f64, max_thrust: f64) -> Vec<&RocketEngine> {
        self.engines
            .values()
            .filter(|engine| {
                engine.thrust_vac >= min_thrust && engine.thrust_vac <= max_thrust
            })
            .collect()
    }
    
    /// Get all engine names
    pub fn list_engines(&self) -> Vec<&String> {
        self.engines.keys().collect()
    }
    
    /// Engine comparison metrics
    pub fn compare_engines(&self, engine1: &str, engine2: &str) -> Option<EngineComparison> {
        let e1 = self.engines.get(engine1)?;
        let e2 = self.engines.get(engine2)?;
        
        Some(EngineComparison {
            engine1: e1.clone(),
            engine2: e2.clone(),
            thrust_ratio: e1.thrust_vac / e2.thrust_vac,
            isp_ratio: e1.isp_vac / e2.isp_vac,
            twr_ratio: e1.thrust_to_weight / e2.thrust_to_weight,
            mass_ratio: e1.mass / e2.mass,
        })
    }
}

#[derive(Debug, Clone)]
pub struct EngineComparison {
    pub engine1: RocketEngine,
    pub engine2: RocketEngine,
    pub thrust_ratio: f64,
    pub isp_ratio: f64,
    pub twr_ratio: f64,
    pub mass_ratio: f64,
}

fn classify_thrust(thrust: f64) -> ThrustClass {
    if thrust < 10.0 {
        ThrustClass::Micro
    } else if thrust < 100.0 {
        ThrustClass::Small
    } else if thrust < 1000.0 {
        ThrustClass::Medium
    } else if thrust < 10000.0 {
        ThrustClass::Large
    } else {
        ThrustClass::SuperHeavy
    }
}

// SpaceX Engines

fn merlin_1d_sea_level() -> RocketEngine {
    RocketEngine {
        name: "Merlin 1D".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Merlin".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(845.0),
        thrust_vac: 934.0,
        isp_sl: Some(282.0),
        isp_vac: 311.0,
        mass: 470.0,
        chamber_pressure: 97.0,
        mixture_ratio: 2.36,
        thrust_to_weight: 203.0,
        dimensions: EngineDimensions {
            length: 2.92,
            diameter: 0.92,
            exit_diameter: 0.92,
            expansion_ratio: 16.0,
        },
        first_flight: 2013,
        status: EngineStatus::Operational,
        features: vec![
            "Throttleable 70%-100%".to_string(),
            "Restart capable".to_string(),
            "Grid fins for landing".to_string(),
            "Gas generator cycle".to_string(),
        ],
        variants: vec![
            "Merlin 1D".to_string(),
            "Merlin 1D+".to_string(),
            "Merlin 1D Block 5".to_string(),
        ],
        applications: vec![
            "Falcon 9".to_string(),
            "Falcon Heavy".to_string(),
        ],
    }
}

fn merlin_1d_vacuum() -> RocketEngine {
    RocketEngine {
        name: "Merlin 1D Vacuum".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Merlin".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::GasGenerator,
        thrust_sl: None,
        thrust_vac: 981.0,
        isp_sl: None,
        isp_vac: 348.0,
        mass: 490.0,
        chamber_pressure: 97.0,
        mixture_ratio: 2.36,
        thrust_to_weight: 204.0,
        dimensions: EngineDimensions {
            length: 4.0,
            diameter: 0.92,
            exit_diameter: 2.4,
            expansion_ratio: 117.0,
        },
        first_flight: 2013,
        status: EngineStatus::Operational,
        features: vec![
            "Extendable nozzle".to_string(),
            "Multiple restart capable".to_string(),
            "Deep throttling".to_string(),
        ],
        variants: vec!["MVac".to_string(), "MVac+".to_string()],
        applications: vec![
            "Falcon 9 second stage".to_string(),
            "Falcon Heavy center core".to_string(),
        ],
    }
}

fn merlin_1d_plus() -> RocketEngine {
    RocketEngine {
        name: "Merlin 1D+".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Merlin".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(854.0),
        thrust_vac: 943.0,
        isp_sl: Some(285.0),
        isp_vac: 313.0,
        mass: 470.0,
        chamber_pressure: 108.0,
        mixture_ratio: 2.36,
        thrust_to_weight: 205.0,
        dimensions: EngineDimensions {
            length: 2.92,
            diameter: 0.92,
            exit_diameter: 0.92,
            expansion_ratio: 16.0,
        },
        first_flight: 2018,
        status: EngineStatus::Operational,
        features: vec![
            "Improved performance".to_string(),
            "Enhanced reliability".to_string(),
            "Block 5 configuration".to_string(),
        ],
        variants: vec!["M1D+ v1.1".to_string()],
        applications: vec![
            "Falcon 9 Block 5".to_string(),
            "Falcon Heavy".to_string(),
        ],
    }
}

fn raptor_1() -> RocketEngine {
    RocketEngine {
        name: "Raptor 1".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Raptor".to_string(),
        propellant: PropellantType::LoxCh4,
        cycle: EngineType::FullFlow,
        thrust_sl: Some(1800.0),
        thrust_vac: 2000.0,
        isp_sl: Some(327.0),
        isp_vac: 363.0,
        mass: 1600.0,
        chamber_pressure: 300.0,
        mixture_ratio: 3.6,
        thrust_to_weight: 127.0,
        dimensions: EngineDimensions {
            length: 3.1,
            diameter: 1.3,
            exit_diameter: 1.3,
            expansion_ratio: 40.0,
        },
        first_flight: 2019,
        status: EngineStatus::Testing,
        features: vec![
            "Full-flow staged combustion".to_string(),
            "Deep throttling 40%-100%".to_string(),
            "Rapid reusability".to_string(),
            "Autogenous pressurization".to_string(),
        ],
        variants: vec!["Raptor".to_string(), "Raptor Vacuum".to_string()],
        applications: vec!["Starship".to_string(), "Super Heavy".to_string()],
    }
}

fn raptor_2() -> RocketEngine {
    RocketEngine {
        name: "Raptor 2".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Raptor".to_string(),
        propellant: PropellantType::LoxCh4,
        cycle: EngineType::FullFlow,
        thrust_sl: Some(2300.0),
        thrust_vac: 2560.0,
        isp_sl: Some(327.0),
        isp_vac: 363.0,
        mass: 1500.0,
        chamber_pressure: 350.0,
        mixture_ratio: 3.6,
        thrust_to_weight: 174.0,
        dimensions: EngineDimensions {
            length: 3.0,
            diameter: 1.3,
            exit_diameter: 1.3,
            expansion_ratio: 40.0,
        },
        first_flight: 2022,
        status: EngineStatus::Operational,
        features: vec![
            "Simplified design".to_string(),
            "Reduced part count".to_string(),
            "Improved manufacturing".to_string(),
            "Higher chamber pressure".to_string(),
        ],
        variants: vec![
            "Raptor 2".to_string(),
            "Raptor 2 Vacuum".to_string(),
        ],
        applications: vec!["Starship".to_string(), "Super Heavy".to_string()],
    }
}

fn raptor_3() -> RocketEngine {
    RocketEngine {
        name: "Raptor 3".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Raptor".to_string(),
        propellant: PropellantType::LoxCh4,
        cycle: EngineType::FullFlow,
        thrust_sl: Some(2600.0),
        thrust_vac: 2890.0,
        isp_sl: Some(330.0),
        isp_vac: 366.0,
        mass: 1300.0,
        chamber_pressure: 400.0,
        mixture_ratio: 3.6,
        thrust_to_weight: 226.0,
        dimensions: EngineDimensions {
            length: 2.8,
            diameter: 1.3,
            exit_diameter: 1.3,
            expansion_ratio: 45.0,
        },
        first_flight: 2024,
        status: EngineStatus::Development,
        features: vec![
            "Ultra-high chamber pressure".to_string(),
            "Optimized combustion".to_string(),
            "Reduced mass".to_string(),
            "Enhanced performance".to_string(),
        ],
        variants: vec![
            "Raptor 3".to_string(),
            "Raptor 3 Vacuum".to_string(),
        ],
        applications: vec!["Starship V2".to_string(), "Super Heavy V2".to_string()],
    }
}

fn super_draco() -> RocketEngine {
    RocketEngine {
        name: "SuperDraco".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Draco".to_string(),
        propellant: PropellantType::NtoMmh,
        cycle: EngineType::PressureFed,
        thrust_sl: Some(73.0),
        thrust_vac: 73.0,
        isp_sl: Some(240.0),
        isp_vac: 235.0,
        mass: 150.0,
        chamber_pressure: 57.0,
        mixture_ratio: 1.85,
        thrust_to_weight: 49.6,
        dimensions: EngineDimensions {
            length: 0.9,
            diameter: 0.3,
            exit_diameter: 0.3,
            expansion_ratio: 2.0,
        },
        first_flight: 2015,
        status: EngineStatus::Operational,
        features: vec![
            "Hypergolic propellants".to_string(),
            "3D printed combustor".to_string(),
            "Deep throttling 20%-100%".to_string(),
            "Emergency abort capability".to_string(),
        ],
        variants: vec!["SuperDraco".to_string()],
        applications: vec!["Dragon 2".to_string(), "Crew Dragon".to_string()],
    }
}

fn draco() -> RocketEngine {
    RocketEngine {
        name: "Draco".to_string(),
        manufacturer: "SpaceX".to_string(),
        family: "Draco".to_string(),
        propellant: PropellantType::NtoMmh,
        cycle: EngineType::PressureFed,
        thrust_sl: Some(0.4),
        thrust_vac: 0.4,
        isp_sl: Some(300.0),
        isp_vac: 300.0,
        mass: 2.0,
        chamber_pressure: 15.0,
        mixture_ratio: 1.85,
        thrust_to_weight: 20.4,
        dimensions: EngineDimensions {
            length: 0.3,
            diameter: 0.05,
            exit_diameter: 0.05,
            expansion_ratio: 1.0,
        },
        first_flight: 2012,
        status: EngineStatus::Operational,
        features: vec![
            "Attitude control thruster".to_string(),
            "Hypergolic ignition".to_string(),
            "High reliability".to_string(),
        ],
        variants: vec!["Draco".to_string()],
        applications: vec![
            "Dragon 1".to_string(),
            "Dragon 2".to_string(),
            "Falcon 9 second stage RCS".to_string(),
        ],
    }
}

// NASA Engines

fn rs25() -> RocketEngine {
    RocketEngine {
        name: "RS-25".to_string(),
        manufacturer: "Aerojet Rocketdyne".to_string(),
        family: "SSME".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::StagedCombustionFuelRich,
        thrust_sl: Some(1860.0),
        thrust_vac: 2279.0,
        isp_sl: Some(366.0),
        isp_vac: 452.3,
        mass: 3527.0,
        chamber_pressure: 206.8,
        mixture_ratio: 6.0,
        thrust_to_weight: 65.9,
        dimensions: EngineDimensions {
            length: 4.3,
            diameter: 2.4,
            exit_diameter: 2.4,
            expansion_ratio: 69.0,
        },
        first_flight: 1981,
        status: EngineStatus::Operational,
        features: vec![
            "Staged combustion cycle".to_string(),
            "Throttleable 67%-109%".to_string(),
            "High performance hydrogen engine".to_string(),
            "Space Shuttle heritage".to_string(),
        ],
        variants: vec![
            "SSME Block I".to_string(),
            "SSME Block II".to_string(),
            "RS-25D".to_string(),
        ],
        applications: vec![
            "Space Shuttle".to_string(),
            "SLS".to_string(),
        ],
    }
}

fn rl10b2() -> RocketEngine {
    RocketEngine {
        name: "RL10B-2".to_string(),
        manufacturer: "Aerojet Rocketdyne".to_string(),
        family: "RL10".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::Expander,
        thrust_sl: None,
        thrust_vac: 110.1,
        isp_sl: None,
        isp_vac: 462.4,
        mass: 277.0,
        chamber_pressure: 43.0,
        mixture_ratio: 5.88,
        thrust_to_weight: 40.5,
        dimensions: EngineDimensions {
            length: 4.14,
            diameter: 2.13,
            exit_diameter: 2.13,
            expansion_ratio: 285.0,
        },
        first_flight: 1998,
        status: EngineStatus::Operational,
        features: vec![
            "Expander cycle".to_string(),
            "Extendable nozzle".to_string(),
            "Multiple restart capability".to_string(),
            "High specific impulse".to_string(),
        ],
        variants: vec![
            "RL10A-4-2".to_string(),
            "RL10B-2".to_string(),
            "RL10C-1".to_string(),
        ],
        applications: vec![
            "Delta IV upper stage".to_string(),
            "Centaur upper stage".to_string(),
        ],
    }
}

fn j2_engine() -> RocketEngine {
    RocketEngine {
        name: "J-2".to_string(),
        manufacturer: "Rocketdyne".to_string(),
        family: "J-2".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::GasGenerator,
        thrust_sl: None,
        thrust_vac: 1033.1,
        isp_sl: None,
        isp_vac: 421.0,
        mass: 1788.0,
        chamber_pressure: 51.7,
        mixture_ratio: 5.5,
        thrust_to_weight: 58.9,
        dimensions: EngineDimensions {
            length: 3.37,
            diameter: 2.0,
            exit_diameter: 2.0,
            expansion_ratio: 27.5,
        },
        first_flight: 1966,
        status: EngineStatus::Retired,
        features: vec![
            "Apollo program heritage".to_string(),
            "Restart capability".to_string(),
            "Gas generator cycle".to_string(),
        ],
        variants: vec!["J-2".to_string(), "J-2S".to_string()],
        applications: vec!["Saturn V S-II".to_string(), "Saturn V S-IVB".to_string()],
    }
}

fn rs68() -> RocketEngine {
    RocketEngine {
        name: "RS-68A".to_string(),
        manufacturer: "Aerojet Rocketdyne".to_string(),
        family: "RS-68".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(3312.0),
        thrust_vac: 3560.0,
        isp_sl: Some(362.0),
        isp_vac: 409.0,
        mass: 6696.0,
        chamber_pressure: 91.0,
        mixture_ratio: 6.0,
        thrust_to_weight: 54.2,
        dimensions: EngineDimensions {
            length: 5.2,
            diameter: 2.43,
            exit_diameter: 2.43,
            expansion_ratio: 21.5,
        },
        first_flight: 2002,
        status: EngineStatus::Retired,
        features: vec![
            "Simplified gas generator".to_string(),
            "Low cost design".to_string(),
            "High thrust".to_string(),
        ],
        variants: vec!["RS-68".to_string(), "RS-68A".to_string()],
        applications: vec!["Delta IV Heavy".to_string()],
    }
}

// Russian Engines

fn rd180() -> RocketEngine {
    RocketEngine {
        name: "RD-180".to_string(),
        manufacturer: "NPO Energomash".to_string(),
        family: "RD-170".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::StagedCombustionOxRich,
        thrust_sl: Some(3830.0),
        thrust_vac: 4152.0,
        isp_sl: Some(311.3),
        isp_vac: 338.0,
        mass: 5307.0,
        chamber_pressure: 266.8,
        mixture_ratio: 2.72,
        thrust_to_weight: 79.8,
        dimensions: EngineDimensions {
            length: 3.56,
            diameter: 3.15,
            exit_diameter: 1.56,
            expansion_ratio: 36.87,
        },
        first_flight: 2000,
        status: EngineStatus::Operational,
        features: vec![
            "Oxygen-rich staged combustion".to_string(),
            "Twin chamber design".to_string(),
            "High chamber pressure".to_string(),
            "Excellent reliability record".to_string(),
        ],
        variants: vec!["RD-180".to_string()],
        applications: vec!["Atlas V".to_string()],
    }
}

fn rd191() -> RocketEngine {
    RocketEngine {
        name: "RD-191".to_string(),
        manufacturer: "NPO Energomash".to_string(),
        family: "RD-170".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::StagedCombustionOxRich,
        thrust_sl: Some(2085.0),
        thrust_vac: 2260.0,
        isp_sl: Some(311.0),
        isp_vac: 337.0,
        mass: 2290.0,
        chamber_pressure: 267.0,
        mixture_ratio: 2.6,
        thrust_to_weight: 92.7,
        dimensions: EngineDimensions {
            length: 3.78,
            diameter: 1.8,
            exit_diameter: 1.3,
            expansion_ratio: 32.0,
        },
        first_flight: 2013,
        status: EngineStatus::Operational,
        features: vec![
            "Single chamber RD-180 derivative".to_string(),
            "Throttling capability".to_string(),
            "High reliability".to_string(),
        ],
        variants: vec!["RD-191".to_string(), "RD-181".to_string()],
        applications: vec!["Angara".to_string(), "Antares".to_string()],
    }
}

fn nk33() -> RocketEngine {
    RocketEngine {
        name: "NK-33".to_string(),
        manufacturer: "Kuznetsov Design Bureau".to_string(),
        family: "NK".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::StagedCombustionOxRich,
        thrust_sl: Some(1512.0),
        thrust_vac: 1667.0,
        isp_sl: Some(297.0),
        isp_vac: 331.0,
        mass: 1222.0,
        chamber_pressure: 145.0,
        mixture_ratio: 2.8,
        thrust_to_weight: 126.3,
        dimensions: EngineDimensions {
            length: 3.7,
            diameter: 2.0,
            exit_diameter: 1.3,
            expansion_ratio: 28.9,
        },
        first_flight: 2013,
        status: EngineStatus::Retired,
        features: vec![
            "Soviet N1 rocket heritage".to_string(),
            "High thrust-to-weight ratio".to_string(),
            "Staged combustion cycle".to_string(),
        ],
        variants: vec!["NK-33".to_string(), "AJ26".to_string()],
        applications: vec!["Antares".to_string(), "Soyuz 2-1v".to_string()],
    }
}

fn rd170() -> RocketEngine {
    RocketEngine {
        name: "RD-170".to_string(),
        manufacturer: "NPO Energomash".to_string(),
        family: "RD-170".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::StagedCombustionOxRich,
        thrust_sl: Some(7904.0),
        thrust_vac: 8083.0,
        isp_sl: Some(309.0),
        isp_vac: 337.0,
        mass: 9750.0,
        chamber_pressure: 245.0,
        mixture_ratio: 2.63,
        thrust_to_weight: 82.4,
        dimensions: EngineDimensions {
            length: 4.0,
            diameter: 4.15,
            exit_diameter: 2.4,
            expansion_ratio: 36.4,
        },
        first_flight: 1987,
        status: EngineStatus::Retired,
        features: vec![
            "Most powerful rocket engine ever flown".to_string(),
            "Four-chamber design".to_string(),
            "Energia rocket".to_string(),
        ],
        variants: vec!["RD-170".to_string(), "RD-171".to_string()],
        applications: vec!["Energia".to_string(), "Zenit".to_string()],
    }
}

fn rd253() -> RocketEngine {
    RocketEngine {
        name: "RD-253".to_string(),
        manufacturer: "Valentin Glushko".to_string(),
        family: "RD-250".to_string(),
        propellant: PropellantType::NtoUdmh,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(1635.0),
        thrust_vac: 1747.0,
        isp_sl: Some(285.0),
        isp_vac: 316.0,
        mass: 1280.0,
        chamber_pressure: 83.0,
        mixture_ratio: 2.6,
        thrust_to_weight: 130.3,
        dimensions: EngineDimensions {
            length: 2.6,
            diameter: 1.5,
            exit_diameter: 1.4,
            expansion_ratio: 15.9,
        },
        first_flight: 1965,
        status: EngineStatus::Operational,
        features: vec![
            "Hypergolic propellants".to_string(),
            "Storable propellants".to_string(),
            "High reliability".to_string(),
        ],
        variants: vec!["RD-253".to_string(), "RD-275".to_string()],
        applications: vec!["Proton".to_string()],
    }
}

// ESA Engines

fn vulcain_2() -> RocketEngine {
    RocketEngine {
        name: "Vulcain 2".to_string(),
        manufacturer: "ArianeGroup".to_string(),
        family: "Vulcain".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(1340.0),
        thrust_vac: 1390.0,
        isp_sl: Some(318.0),
        isp_vac: 434.0,
        mass: 1800.0,
        chamber_pressure: 116.0,
        mixture_ratio: 6.1,
        thrust_to_weight: 78.7,
        dimensions: EngineDimensions {
            length: 3.0,
            diameter: 1.7,
            exit_diameter: 2.1,
            expansion_ratio: 58.5,
        },
        first_flight: 2005,
        status: EngineStatus::Operational,
        features: vec![
            "Cryogenic main engine".to_string(),
            "Improved Vulcain".to_string(),
            "European technology".to_string(),
        ],
        variants: vec!["Vulcain 2".to_string(), "Vulcain 2.1".to_string()],
        applications: vec!["Ariane 5".to_string()],
    }
}

fn vinci() -> RocketEngine {
    RocketEngine {
        name: "Vinci".to_string(),
        manufacturer: "ArianeGroup".to_string(),
        family: "Vinci".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::Expander,
        thrust_sl: None,
        thrust_vac: 180.0,
        isp_sl: None,
        isp_vac: 464.0,
        mass: 550.0,
        chamber_pressure: 61.0,
        mixture_ratio: 5.8,
        thrust_to_weight: 33.4,
        dimensions: EngineDimensions {
            length: 4.2,
            diameter: 2.2,
            exit_diameter: 2.2,
            expansion_ratio: 240.0,
        },
        first_flight: 2023,
        status: EngineStatus::Operational,
        features: vec![
            "Multiple restart capability".to_string(),
            "High specific impulse".to_string(),
            "Expander cycle".to_string(),
        ],
        variants: vec!["Vinci".to_string()],
        applications: vec!["Ariane 6 upper stage".to_string()],
    }
}

fn aestus() -> RocketEngine {
    RocketEngine {
        name: "Aestus".to_string(),
        manufacturer: "ArianeGroup".to_string(),
        family: "Aestus".to_string(),
        propellant: PropellantType::NtoMmh,
        cycle: EngineType::PressureFed,
        thrust_sl: None,
        thrust_vac: 29.5,
        isp_sl: None,
        isp_vac: 324.0,
        mass: 111.0,
        chamber_pressure: 11.4,
        mixture_ratio: 2.04,
        thrust_to_weight: 27.1,
        dimensions: EngineDimensions {
            length: 1.8,
            diameter: 0.6,
            exit_diameter: 1.34,
            expansion_ratio: 84.0,
        },
        first_flight: 1999,
        status: EngineStatus::Operational,
        features: vec![
            "Storable propellants".to_string(),
            "Multiple restart".to_string(),
            "High expansion ratio".to_string(),
        ],
        variants: vec!["Aestus".to_string(), "Aestus II".to_string()],
        applications: vec!["Ariane 5 ESV".to_string()],
    }
}

// Blue Origin Engines

fn be3() -> RocketEngine {
    RocketEngine {
        name: "BE-3".to_string(),
        manufacturer: "Blue Origin".to_string(),
        family: "BE".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(490.0),
        thrust_vac: 550.0,
        isp_sl: Some(365.0),
        isp_vac: 410.0,
        mass: 750.0,
        chamber_pressure: 46.0,
        mixture_ratio: 6.0,
        thrust_to_weight: 74.7,
        dimensions: EngineDimensions {
            length: 2.7,
            diameter: 1.4,
            exit_diameter: 1.4,
            expansion_ratio: 40.0,
        },
        first_flight: 2015,
        status: EngineStatus::Operational,
        features: vec![
            "Throttleable 18%-100%".to_string(),
            "Gimbaling".to_string(),
            "Restart capable".to_string(),
        ],
        variants: vec!["BE-3".to_string(), "BE-3U".to_string()],
        applications: vec!["New Shepard".to_string()],
    }
}

fn be4() -> RocketEngine {
    RocketEngine {
        name: "BE-4".to_string(),
        manufacturer: "Blue Origin".to_string(),
        family: "BE".to_string(),
        propellant: PropellantType::LoxCh4,
        cycle: EngineType::StagedCombustionOxRich,
        thrust_sl: Some(2400.0),
        thrust_vac: 2670.0,
        isp_sl: Some(310.0),
        isp_vac: 343.0,
        mass: 2700.0,
        chamber_pressure: 134.0,
        mixture_ratio: 3.3,
        thrust_to_weight: 90.6,
        dimensions: EngineDimensions {
            length: 3.7,
            diameter: 1.7,
            exit_diameter: 1.7,
            expansion_ratio: 27.0,
        },
        first_flight: 2024,
        status: EngineStatus::Operational,
        features: vec![
            "Oxygen-rich staged combustion".to_string(),
            "Methane fuel".to_string(),
            "Throttling capability".to_string(),
        ],
        variants: vec!["BE-4".to_string()],
        applications: vec!["New Glenn".to_string(), "Vulcan Centaur".to_string()],
    }
}

// Historic Engines

fn f1_engine() -> RocketEngine {
    RocketEngine {
        name: "F-1".to_string(),
        manufacturer: "Rocketdyne".to_string(),
        family: "F-1".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(6770.0),
        thrust_vac: 7740.5,
        isp_sl: Some(263.0),
        isp_vac: 304.0,
        mass: 8391.0,
        chamber_pressure: 70.0,
        mixture_ratio: 2.27,
        thrust_to_weight: 94.1,
        dimensions: EngineDimensions {
            length: 5.79,
            diameter: 3.76,
            exit_diameter: 3.76,
            expansion_ratio: 16.0,
        },
        first_flight: 1967,
        status: EngineStatus::Retired,
        features: vec![
            "Most powerful single-chamber engine".to_string(),
            "Apollo Saturn V".to_string(),
            "Gas generator cycle".to_string(),
            "Iconic moon rocket engine".to_string(),
        ],
        variants: vec!["F-1".to_string()],
        applications: vec!["Saturn V S-IC".to_string()],
    }
}

fn h1_engine() -> RocketEngine {
    RocketEngine {
        name: "H-1".to_string(),
        manufacturer: "Rocketdyne".to_string(),
        family: "H-1".to_string(),
        propellant: PropellantType::LoxRp1,
        cycle: EngineType::GasGenerator,
        thrust_sl: Some(947.0),
        thrust_vac: 1030.0,
        isp_sl: Some(289.0),
        isp_vac: 315.0,
        mass: 635.0,
        chamber_pressure: 57.6,
        mixture_ratio: 2.23,
        thrust_to_weight: 152.0,
        dimensions: EngineDimensions {
            length: 2.29,
            diameter: 1.1,
            exit_diameter: 1.1,
            expansion_ratio: 8.0,
        },
        first_flight: 1961,
        status: EngineStatus::Retired,
        features: vec![
            "Saturn I booster".to_string(),
            "Clustered configuration".to_string(),
            "Apollo program heritage".to_string(),
        ],
        variants: vec!["H-1".to_string()],
        applications: vec!["Saturn I".to_string(), "Saturn IB".to_string()],
    }
}

fn ssme_rs25() -> RocketEngine {
    RocketEngine {
        name: "SSME (RS-25)".to_string(),
        manufacturer: "Rocketdyne".to_string(),
        family: "SSME".to_string(),
        propellant: PropellantType::LoxLh2,
        cycle: EngineType::StagedCombustionFuelRich,
        thrust_sl: Some(1860.0),
        thrust_vac: 2278.0,
        isp_sl: Some(366.0),
        isp_vac: 452.3,
        mass: 3177.0,
        chamber_pressure: 206.8,
        mixture_ratio: 6.0,
        thrust_to_weight: 73.1,
        dimensions: EngineDimensions {
            length: 4.3,
            diameter: 2.4,
            exit_diameter: 2.4,
            expansion_ratio: 69.0,
        },
        first_flight: 1981,
        status: EngineStatus::Retired,
        features: vec![
            "Space Shuttle Main Engine".to_string(),
            "Reusable".to_string(),
            "High performance".to_string(),
            "Throttleable".to_string(),
        ],
        variants: vec!["SSME Phase I".to_string(), "SSME Phase II".to_string()],
        applications: vec!["Space Shuttle".to_string()],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_engine_database_creation() {
        let db = EngineDatabase::new();
        assert!(db.engines.len() > 20, "Should have many engines loaded");
    }
    
    #[test]
    fn test_get_spacex_engines() {
        let db = EngineDatabase::new();
        let spacex_engines = db.get_by_manufacturer("SpaceX");
        
        assert!(spacex_engines.len() >= 6, "Should have multiple SpaceX engines");
        
        let merlin = spacex_engines.iter()
            .find(|e| e.name.contains("Merlin"))
            .expect("Should have Merlin engine");
        
        assert_eq!(merlin.propellant, PropellantType::LoxRp1);
        assert_eq!(merlin.cycle, EngineType::GasGenerator);
    }
    
    #[test]
    fn test_raptor_engine_specs() {
        let db = EngineDatabase::new();
        let raptor2 = db.get_engine("Raptor 2").expect("Should have Raptor 2");
        
        assert_eq!(raptor2.propellant, PropellantType::LoxCh4);
        assert_eq!(raptor2.cycle, EngineType::FullFlow);
        assert!(raptor2.chamber_pressure > 300.0, "Raptor should have very high chamber pressure");
        assert!(raptor2.thrust_vac > 2000.0, "Raptor 2 should have >2000 kN thrust");
    }
    
    #[test]
    fn test_engine_comparison() {
        let db = EngineDatabase::new();
        let comparison = db.compare_engines("Raptor 2", "Merlin 1D")
            .expect("Should be able to compare engines");
        
        assert!(comparison.thrust_ratio > 2.0, "Raptor should have much more thrust");
        assert!(comparison.isp_ratio > 1.1, "Raptor should have better ISP");
    }
    
    #[test]
    fn test_propellant_filtering() {
        let db = EngineDatabase::new();
        let methalox_engines = db.get_by_propellant(PropellantType::LoxCh4);
        
        assert!(methalox_engines.len() >= 3, "Should have multiple methalox engines");
        
        let has_raptor = methalox_engines.iter()
            .any(|e| e.name.contains("Raptor"));
        assert!(has_raptor, "Should include Raptor engines");
        
        let has_be4 = methalox_engines.iter()
            .any(|e| e.name.contains("BE-4"));
        assert!(has_be4, "Should include BE-4 engine");
    }
    
    #[test]
    fn test_thrust_classification() {
        assert_eq!(classify_thrust(5.0), ThrustClass::Micro);
        assert_eq!(classify_thrust(50.0), ThrustClass::Small);
        assert_eq!(classify_thrust(500.0), ThrustClass::Medium);
        assert_eq!(classify_thrust(5000.0), ThrustClass::Large);
        assert_eq!(classify_thrust(15000.0), ThrustClass::SuperHeavy);
    }
    
    #[test]
    fn test_thrust_range_query() {
        let db = EngineDatabase::new();
        let high_thrust = db.get_by_thrust_range(2000.0, 5000.0);
        
        assert!(high_thrust.len() > 0, "Should find high-thrust engines");
        
        let has_raptor = high_thrust.iter()
            .any(|e| e.name.contains("Raptor"));
        assert!(has_raptor, "Should include Raptor engines in high thrust range");
    }
    
    #[test]
    fn test_rs25_specifications() {
        let db = EngineDatabase::new();
        let rs25 = db.get_engine("RS-25").expect("Should have RS-25");
        
        assert_eq!(rs25.propellant, PropellantType::LoxLh2);
        assert_eq!(rs25.cycle, EngineType::StagedCombustionFuelRich);
        assert!(rs25.isp_vac > 450.0, "RS-25 should have very high ISP");
        assert!(rs25.chamber_pressure > 200.0, "RS-25 should have high chamber pressure");
    }
    
    #[test]
    fn test_f1_engine_legacy() {
        let db = EngineDatabase::new();
        let f1 = db.get_engine("F-1").expect("Should have F-1 engine");
        
        assert_eq!(f1.status, EngineStatus::Retired);
        assert!(f1.thrust_sl.unwrap() > 6000.0, "F-1 should have massive thrust");
        assert!(f1.applications.contains(&"Saturn V S-IC".to_string()));
    }
}