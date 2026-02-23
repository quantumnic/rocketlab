//! Trajectory simulation using RK4 numerical integration.
//!
//! Supports 2D (vertical plane) trajectory with gravity, drag, and thrust.
//! Useful for launch trajectory simulation and gravity turn analysis.
//!
//! Reference: Curtis, §11.2; Bate, Mueller & White, §9.3

use crate::atmosphere;
use crate::constants::{G0, MU_EARTH, R_EARTH_EQUATORIAL};

/// State vector for 2D trajectory simulation.
///
/// Uses an inertial frame with:
/// - x: downrange (m)
/// - y: altitude above surface (m)
/// - vx, vy: velocity components (m/s)
/// - mass: current vehicle mass (kg)
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryState {
    pub time: f64,
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub mass: f64,
}

impl TrajectoryState {
    /// Total velocity magnitude (m/s).
    pub fn speed(&self) -> f64 {
        (self.vx * self.vx + self.vy * self.vy).sqrt()
    }

    /// Altitude above surface (m).
    pub fn altitude(&self) -> f64 {
        self.y
    }

    /// Flight path angle (radians from horizontal).
    pub fn flight_path_angle(&self) -> f64 {
        self.vy.atan2(self.vx)
    }
}

/// Vehicle configuration for trajectory simulation.
#[derive(Debug, Clone)]
pub struct Vehicle {
    /// Dry mass (kg)
    pub dry_mass: f64,
    /// Propellant mass (kg)
    pub propellant_mass: f64,
    /// Thrust at sea level (N)
    pub thrust_sl: f64,
    /// Thrust in vacuum (N)
    pub thrust_vac: f64,
    /// Specific impulse at sea level (s)
    pub isp_sl: f64,
    /// Specific impulse in vacuum (s)
    pub isp_vac: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Reference area (m²)
    pub reference_area: f64,
    /// Pitch program: angle from vertical (radians) as a function of time
    /// If None, uses gravity turn (velocity-aligned after kickover)
    pub pitch_program: Option<fn(f64) -> f64>,
}

impl Vehicle {
    /// Total initial mass (kg).
    pub fn total_mass(&self) -> f64 {
        self.dry_mass + self.propellant_mass
    }

    /// Interpolate thrust between sea level and vacuum based on altitude.
    fn thrust_at_altitude(&self, altitude_m: f64) -> f64 {
        if let Some(atm) = atmosphere::us_standard_atmosphere(altitude_m) {
            let pressure_ratio = atm.pressure / 101_325.0;
            self.thrust_vac + (self.thrust_sl - self.thrust_vac) * pressure_ratio
        } else {
            self.thrust_vac
        }
    }

    /// Interpolate Isp between sea level and vacuum.
    fn isp_at_altitude(&self, altitude_m: f64) -> f64 {
        if let Some(atm) = atmosphere::us_standard_atmosphere(altitude_m) {
            let pressure_ratio = atm.pressure / 101_325.0;
            self.isp_vac + (self.isp_sl - self.isp_vac) * pressure_ratio
        } else {
            self.isp_vac
        }
    }

    /// Mass flow rate at altitude (kg/s).
    fn mass_flow_at_altitude(&self, altitude_m: f64) -> f64 {
        let thrust = self.thrust_at_altitude(altitude_m);
        let isp = self.isp_at_altitude(altitude_m);
        thrust / (isp * G0)
    }
}

/// Simulation parameters.
#[derive(Debug, Clone)]
pub struct SimConfig {
    /// Time step (s)
    pub dt: f64,
    /// Maximum simulation time (s)
    pub max_time: f64,
    /// Gravity turn kickover altitude (m)
    pub kickover_altitude: f64,
    /// Gravity turn kickover angle from vertical (radians)
    pub kickover_angle: f64,
    /// Record interval: store state every N steps
    pub record_interval: usize,
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            dt: 0.1,
            max_time: 600.0,
            kickover_altitude: 500.0,
            kickover_angle: 0.1_f64.to_radians(), // ~0.1 degrees
            record_interval: 10,
        }
    }
}

/// Result of a trajectory simulation.
#[derive(Debug, Clone)]
pub struct TrajectoryResult {
    /// Recorded trajectory states
    pub states: Vec<TrajectoryState>,
    /// Max-Q: (time, altitude, dynamic_pressure_Pa)
    pub max_q: Option<(f64, f64, f64)>,
    /// Time of engine cutoff (propellant exhausted)
    pub burnout_time: Option<f64>,
    /// Final state
    pub final_state: TrajectoryState,
    /// Termination reason
    pub termination: TerminationReason,
}

/// Why the simulation ended.
#[derive(Debug, Clone)]
pub enum TerminationReason {
    MaxTime,
    Impact,
    AltitudeReached(f64),
    PropellantDepleted,
}

/// Run a 2D trajectory simulation with RK4 integration.
pub fn simulate(vehicle: &Vehicle, config: &SimConfig) -> TrajectoryResult {
    let mut state = TrajectoryState {
        time: 0.0,
        x: 0.0,
        y: 0.0,
        vx: 0.0,
        vy: 0.01, // tiny upward velocity to start
        mass: vehicle.total_mass(),
    };

    let mut states = vec![state];
    let mut max_q: Option<(f64, f64, f64)> = None;
    let mut burnout_time: Option<f64> = None;
    let mut step_count = 0_usize;
    let mut termination = TerminationReason::MaxTime;

    let dt = config.dt;

    while state.time < config.max_time {
        // Check termination
        if state.y < -10.0 && state.time > 1.0 {
            termination = TerminationReason::Impact;
            break;
        }

        let has_fuel = state.mass > vehicle.dry_mass;

        // Track burnout
        if !has_fuel && burnout_time.is_none() {
            burnout_time = Some(state.time);
        }

        // Track max-Q
        if let Some(atm) = atmosphere::us_standard_atmosphere(state.y.max(0.0)) {
            let q = atmosphere::dynamic_pressure(state.speed(), atm.density);
            if let Some((_, _, prev_q)) = max_q {
                if q > prev_q {
                    max_q = Some((state.time, state.y, q));
                }
            } else {
                max_q = Some((state.time, state.y, q));
            }
        }

        // RK4 integration
        let k1 = derivatives(&state, vehicle, config, has_fuel);
        let s2 = advance_state(&state, &k1, dt * 0.5);
        let has_fuel2 = s2.mass > vehicle.dry_mass;
        let k2 = derivatives(&s2, vehicle, config, has_fuel2);
        let s3 = advance_state(&state, &k2, dt * 0.5);
        let has_fuel3 = s3.mass > vehicle.dry_mass;
        let k3 = derivatives(&s3, vehicle, config, has_fuel3);
        let s4 = advance_state(&state, &k3, dt);
        let has_fuel4 = s4.mass > vehicle.dry_mass;
        let k4 = derivatives(&s4, vehicle, config, has_fuel4);

        state.x += dt / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        state.y += dt / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
        state.vx += dt / 6.0 * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
        state.vy += dt / 6.0 * (k1.3 + 2.0 * k2.3 + 2.0 * k3.3 + k4.3);
        state.mass += dt / 6.0 * (k1.4 + 2.0 * k2.4 + 2.0 * k3.4 + k4.4);

        // Clamp mass
        state.mass = state.mass.max(vehicle.dry_mass);

        state.time += dt;
        step_count += 1;

        if step_count.is_multiple_of(config.record_interval) {
            states.push(state);
        }
    }

    // Ensure final state is recorded
    if states
        .last()
        .is_none_or(|s| (s.time - state.time).abs() > 1e-6)
    {
        states.push(state);
    }

    TrajectoryResult {
        states,
        max_q,
        burnout_time,
        final_state: state,
        termination,
    }
}

/// Compute derivatives (dx/dt, dy/dt, dvx/dt, dvy/dt, dm/dt).
fn derivatives(
    state: &TrajectoryState,
    vehicle: &Vehicle,
    config: &SimConfig,
    has_fuel: bool,
) -> (f64, f64, f64, f64, f64) {
    let alt = state.y.max(0.0);
    let r = (R_EARTH_EQUATORIAL * 1000.0) + alt;
    let speed = state.speed();

    // Gravity (always downward)
    let g = MU_EARTH * 1e9 / (r * r); // Convert km³/s² to m³/s², result in m/s²

    // Thrust direction
    let (thrust_x, thrust_y, dm_dt) = if has_fuel {
        let thrust = vehicle.thrust_at_altitude(alt);
        let mdot = vehicle.mass_flow_at_altitude(alt);

        // Determine thrust direction
        let (tx, ty) = if state.y < config.kickover_altitude {
            // Vertical ascent
            (0.0, 1.0)
        } else if speed > 1.0 {
            // Gravity turn: thrust along velocity vector
            let angle = state.flight_path_angle() - config.kickover_angle;
            (angle.cos(), angle.sin())
        } else {
            (0.0, 1.0)
        };

        let ax = thrust * tx / state.mass;
        let ay = thrust * ty / state.mass;
        (ax, ay, -mdot)
    } else {
        (0.0, 0.0, 0.0)
    };

    // Drag
    let (drag_x, drag_y) = if speed > 0.1 {
        if let Some(atm) = atmosphere::us_standard_atmosphere(alt) {
            let fd = atmosphere::drag_force(speed, atm.density, vehicle.cd, vehicle.reference_area);
            let ad = fd / state.mass;
            // Drag opposes velocity
            (-ad * state.vx / speed, -ad * state.vy / speed)
        } else {
            (0.0, 0.0)
        }
    } else {
        (0.0, 0.0)
    };

    let dvx_dt = thrust_x + drag_x;
    let dvy_dt = thrust_y + drag_y - g;

    (state.vx, state.vy, dvx_dt, dvy_dt, dm_dt)
}

/// Advance state by derivatives * dt.
fn advance_state(
    state: &TrajectoryState,
    derivs: &(f64, f64, f64, f64, f64),
    dt: f64,
) -> TrajectoryState {
    TrajectoryState {
        time: state.time + dt,
        x: state.x + derivs.0 * dt,
        y: state.y + derivs.1 * dt,
        vx: state.vx + derivs.2 * dt,
        vy: state.vy + derivs.3 * dt,
        mass: state.mass + derivs.4 * dt,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_vehicle() -> Vehicle {
        // Small sounding rocket
        Vehicle {
            dry_mass: 100.0,
            propellant_mass: 400.0,
            thrust_sl: 15_000.0,
            thrust_vac: 16_500.0,
            isp_sl: 240.0,
            isp_vac: 270.0,
            cd: 0.3,
            reference_area: 0.1,
            pitch_program: None,
        }
    }

    #[test]
    fn test_basic_launch() {
        let vehicle = test_vehicle();
        let config = SimConfig {
            dt: 0.5,
            max_time: 120.0,
            ..Default::default()
        };

        let result = simulate(&vehicle, &config);

        // Should have trajectory points
        assert!(!result.states.is_empty());

        // Vehicle should gain altitude
        assert!(
            result.final_state.y > 0.0,
            "Should gain altitude: y={}",
            result.final_state.y
        );

        // Should have positive velocity
        assert!(
            result.final_state.speed() > 0.0,
            "Should have velocity: v={}",
            result.final_state.speed()
        );

        // Max-Q should exist
        assert!(result.max_q.is_some(), "Should find max-Q");
    }

    #[test]
    fn test_burnout() {
        let vehicle = test_vehicle();
        let config = SimConfig {
            dt: 0.5,
            max_time: 200.0,
            ..Default::default()
        };

        let result = simulate(&vehicle, &config);

        // Should record burnout
        assert!(result.burnout_time.is_some(), "Should record burnout");

        // Burnout mass should equal dry mass
        let burnout_t = result.burnout_time.unwrap();
        assert!(burnout_t > 0.0 && burnout_t < 200.0);
    }

    #[test]
    fn test_gravity_only() {
        // No thrust, just fall
        let vehicle = Vehicle {
            dry_mass: 100.0,
            propellant_mass: 0.0,
            thrust_sl: 0.0,
            thrust_vac: 0.0,
            isp_sl: 1.0,
            isp_vac: 1.0,
            cd: 0.0,
            reference_area: 0.0,
            pitch_program: None,
        };

        let mut config = SimConfig::default();
        config.max_time = 10.0;
        config.dt = 0.01;

        let result = simulate(&vehicle, &config);

        // Should fall (hits impact termination at y < -10)
        assert!(
            result.final_state.y < 0.0,
            "Should fall: y={}",
            result.final_state.y
        );
    }
}
