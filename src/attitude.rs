//! Spacecraft attitude dynamics and control systems
//!
//! Comprehensive attitude dynamics including:
//! - Quaternion kinematics and Euler angle conversions
//! - Rigid body dynamics with external torques
//! - Reaction wheel and control moment gyroscope models
//! - Magnetorquer sizing and control
//! - Gravity gradient torque calculations
//! - Atmospheric torque modeling  
//! - Solar radiation pressure torque
//! - Attitude determination and control (ADCS)
//!
//! References:
//! - Wertz "Spacecraft Attitude Determination and Control" (1978)
//! - Sidi "Spacecraft Dynamics and Control" (1997)
//! - Markley & Crassidis "Fundamentals of Spacecraft Attitude Determination and Control" (2014)

use nalgebra::{Vector3, Matrix3, Quaternion, UnitQuaternion};
use crate::constants::MU_EARTH;
use std::f64::consts::FRAC_PI_2;

/// Spacecraft attitude state representation
#[derive(Debug, Clone)]
pub struct AttitudeState {
    /// Time since epoch (s)
    pub time: f64,
    /// Attitude quaternion (inertial to body frame)
    pub quaternion: UnitQuaternion<f64>,
    /// Angular velocity in body frame (rad/s)
    pub angular_velocity: Vector3<f64>,
    /// Euler angles (roll, pitch, yaw) in radians
    pub euler_angles: Vector3<f64>,
    /// Angular acceleration (rad/s²)
    pub angular_acceleration: Vector3<f64>,
}

/// Spacecraft inertia properties
#[derive(Debug, Clone)]
pub struct SpacecraftInertia {
    /// Inertia matrix in body frame (kg⋅m²)
    pub inertia_matrix: Matrix3<f64>,
    /// Principal moments of inertia (kg⋅m²)
    pub principal_moments: Vector3<f64>,
    /// Center of mass location (m)
    pub center_of_mass: Vector3<f64>,
    /// Total spacecraft mass (kg)
    pub mass: f64,
}

/// External torque sources acting on spacecraft
#[derive(Debug, Clone)]
pub struct ExternalTorques {
    /// Gravity gradient torque (N⋅m)
    pub gravity_gradient: Vector3<f64>,
    /// Atmospheric torque (N⋅m)
    pub atmospheric: Vector3<f64>,
    /// Solar radiation pressure torque (N⋅m)
    pub solar_radiation: Vector3<f64>,
    /// Magnetic torque from Earth's field (N⋅m)
    pub magnetic: Vector3<f64>,
    /// Control torques from actuators (N⋅m)
    pub control: Vector3<f64>,
    /// Total external torque (N⋅m)
    pub total: Vector3<f64>,
}

/// Reaction wheel configuration and state
#[derive(Debug, Clone)]
pub struct ReactionWheel {
    /// Wheel identification
    pub id: String,
    /// Installation axis in body frame (unit vector)
    pub axis: Vector3<f64>,
    /// Wheel moment of inertia (kg⋅m²)
    pub inertia: f64,
    /// Maximum wheel speed (rad/s)
    pub max_speed: f64,
    /// Maximum torque output (N⋅m)
    pub max_torque: f64,
    /// Current wheel speed (rad/s)
    pub speed: f64,
    /// Wheel angular acceleration (rad/s²)
    pub acceleration: f64,
    /// Power consumption (W)
    pub power: f64,
    /// Wheel health status
    pub status: ActuatorStatus,
}

/// Magnetorquer configuration
#[derive(Debug, Clone)]
pub struct Magnetorquer {
    /// Magnetorquer identification  
    pub id: String,
    /// Magnetic dipole axis in body frame
    pub axis: Vector3<f64>,
    /// Maximum magnetic dipole moment (A⋅m²)
    pub max_dipole: f64,
    /// Current dipole moment (A⋅m²)
    pub dipole: f64,
    /// Resistance (Ω)
    pub resistance: f64,
    /// Power consumption (W)
    pub power: f64,
    /// Status
    pub status: ActuatorStatus,
}

/// Thruster for attitude control
#[derive(Debug, Clone)]
pub struct AttitudeThruster {
    /// Thruster identification
    pub id: String,
    /// Thruster location in body frame (m)
    pub position: Vector3<f64>,
    /// Thrust direction in body frame (unit vector)
    pub direction: Vector3<f64>,
    /// Maximum thrust (N)
    pub max_thrust: f64,
    /// Minimum impulse bit (N⋅s)
    pub min_impulse: f64,
    /// Specific impulse (s)
    pub isp: f64,
    /// Current thrust level (N)
    pub thrust: f64,
    /// Fuel remaining (kg)
    pub fuel_remaining: f64,
    /// Status
    pub status: ActuatorStatus,
}

/// Actuator status
#[derive(Debug, Clone, PartialEq)]
pub enum ActuatorStatus {
    /// Operational and available
    Nominal,
    /// Degraded performance but functional
    Degraded,
    /// Failed and unavailable
    Failed,
    /// Offline/powered off
    Offline,
}

/// Spacecraft attitude control system configuration
#[derive(Debug, Clone)]
pub struct AttitudeControlSystem {
    /// Spacecraft inertia properties
    pub inertia: SpacecraftInertia,
    /// Reaction wheels
    pub reaction_wheels: Vec<ReactionWheel>,
    /// Magnetorquers
    pub magnetorquers: Vec<Magnetorquer>,
    /// Attitude thrusters
    pub thrusters: Vec<AttitudeThruster>,
    /// Control gains and parameters
    pub control_params: ControlParameters,
}

/// Attitude control parameters
#[derive(Debug, Clone)]
pub struct ControlParameters {
    /// Proportional gain for attitude error
    pub kp_attitude: f64,
    /// Derivative gain for angular rate
    pub kd_rate: f64,
    /// Integral gain for steady-state error
    pub ki_integral: f64,
    /// Maximum control torque (N⋅m)
    pub max_control_torque: f64,
    /// Rate limit for attitude commands (rad/s)
    pub max_slew_rate: f64,
    /// Deadband for attitude error (rad)
    pub attitude_deadband: f64,
    /// Deadband for rate error (rad/s)
    pub rate_deadband: f64,
}

/// Attitude dynamics simulator
pub struct AttitudeDynamics {
    /// Control system configuration
    pub acs: AttitudeControlSystem,
    /// Current attitude state
    pub state: AttitudeState,
    /// Orbital position and velocity for torque calculations
    pub orbital_state: OrbitalState,
    /// Environment parameters
    pub environment: EnvironmentParameters,
}

/// Orbital state for attitude dynamics
#[derive(Debug, Clone)]
pub struct OrbitalState {
    /// Position vector in ECI frame (km)
    pub position: Vector3<f64>,
    /// Velocity vector in ECI frame (km/s)  
    pub velocity: Vector3<f64>,
    /// Altitude above Earth (km)
    pub altitude: f64,
    /// Orbital period (s)
    pub period: f64,
}

/// Environmental parameters affecting attitude
#[derive(Debug, Clone)]
pub struct EnvironmentParameters {
    /// Earth's magnetic field vector in body frame (nT)
    pub magnetic_field: Vector3<f64>,
    /// Solar vector in body frame (unit vector)
    pub sun_vector: Vector3<f64>,
    /// Atmospheric density (kg/m³)
    pub density: f64,
    /// Solar flux (W/m²)
    pub solar_flux: f64,
}

impl AttitudeState {
    /// Create new attitude state from quaternion and angular velocity
    pub fn new(quaternion: UnitQuaternion<f64>, angular_velocity: Vector3<f64>) -> Self {
        let euler_angles = quaternion_to_euler_angles(&quaternion);
        
        Self {
            time: 0.0,
            quaternion,
            angular_velocity,
            euler_angles,
            angular_acceleration: Vector3::zeros(),
        }
    }
    
    /// Convert to direction cosine matrix (body to inertial)
    pub fn dcm_body_to_inertial(&self) -> Matrix3<f64> {
        self.quaternion.to_rotation_matrix().matrix().transpose()
    }
    
    /// Convert to direction cosine matrix (inertial to body)  
    pub fn dcm_inertial_to_body(&self) -> Matrix3<f64> {
        *self.quaternion.to_rotation_matrix().matrix()
    }
    
    /// Update Euler angles from quaternion
    pub fn update_euler_angles(&mut self) {
        self.euler_angles = quaternion_to_euler_angles(&self.quaternion);
    }
}

impl SpacecraftInertia {
    /// Create spacecraft inertia from principal moments
    pub fn from_principal_moments(ixx: f64, iyy: f64, izz: f64, mass: f64) -> Self {
        let mut inertia_matrix = Matrix3::zeros();
        inertia_matrix[(0, 0)] = ixx;
        inertia_matrix[(1, 1)] = iyy;  
        inertia_matrix[(2, 2)] = izz;
        
        Self {
            inertia_matrix,
            principal_moments: Vector3::new(ixx, iyy, izz),
            center_of_mass: Vector3::zeros(),
            mass,
        }
    }
    
    /// Create box-shaped satellite inertia
    pub fn box_satellite(length: f64, width: f64, height: f64, mass: f64) -> Self {
        let ixx = mass * (width * width + height * height) / 12.0;
        let iyy = mass * (length * length + height * height) / 12.0;
        let izz = mass * (length * length + width * width) / 12.0;
        
        Self::from_principal_moments(ixx, iyy, izz, mass)
    }
    
    /// Create cylindrical satellite inertia
    pub fn cylinder_satellite(radius: f64, height: f64, mass: f64) -> Self {
        let ixx = mass * (3.0 * radius * radius + height * height) / 12.0;
        let iyy = ixx; // Symmetric about z-axis
        let izz = mass * radius * radius / 2.0;
        
        Self::from_principal_moments(ixx, iyy, izz, mass)
    }
}

impl ReactionWheel {
    /// Create standard momentum wheel configuration
    pub fn standard_wheel(axis: Vector3<f64>, max_speed: f64, max_torque: f64) -> Self {
        Self {
            id: "RW".to_string(),
            axis: axis.normalize(),
            inertia: 0.01, // kg⋅m² - typical small satellite wheel
            max_speed,
            max_torque,
            speed: 0.0,
            acceleration: 0.0,
            power: 0.0,
            status: ActuatorStatus::Nominal,
        }
    }
    
    /// Create 4-wheel pyramid configuration
    pub fn pyramid_configuration() -> Vec<ReactionWheel> {
        let angle = 54.74_f64.to_radians(); // Pyramid angle for equal torque capability
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        vec![
            ReactionWheel::standard_wheel(Vector3::new(cos_angle, 0.0, sin_angle), 6000.0, 0.2),
            ReactionWheel::standard_wheel(Vector3::new(-cos_angle, 0.0, sin_angle), 6000.0, 0.2),
            ReactionWheel::standard_wheel(Vector3::new(0.0, cos_angle, sin_angle), 6000.0, 0.2),
            ReactionWheel::standard_wheel(Vector3::new(0.0, -cos_angle, sin_angle), 6000.0, 0.2),
        ]
    }
    
    /// Calculate torque output based on commanded acceleration
    pub fn calculate_torque(&self, commanded_acceleration: f64) -> f64 {
        commanded_acceleration.clamp(-self.max_torque, self.max_torque)
    }
    
    /// Update wheel state
    pub fn update(&mut self, torque_command: f64, dt: f64) {
        let acceleration = torque_command / self.inertia;
        self.acceleration = acceleration.clamp(-self.max_torque, self.max_torque);
        self.speed += self.acceleration * dt;
        self.speed = self.speed.clamp(-self.max_speed, self.max_speed);
        
        // Simple power model
        self.power = self.acceleration.abs() * 5.0 + self.speed.abs() * 0.1;
    }
}

impl Magnetorquer {
    /// Create 3-axis magnetorquer system
    pub fn three_axis_system(max_dipole: f64) -> Vec<Magnetorquer> {
        vec![
            Magnetorquer {
                id: "MTQ_X".to_string(),
                axis: Vector3::new(1.0, 0.0, 0.0),
                max_dipole,
                dipole: 0.0,
                resistance: 10.0,
                power: 0.0,
                status: ActuatorStatus::Nominal,
            },
            Magnetorquer {
                id: "MTQ_Y".to_string(),
                axis: Vector3::new(0.0, 1.0, 0.0),
                max_dipole,
                dipole: 0.0,
                resistance: 10.0,
                power: 0.0,
                status: ActuatorStatus::Nominal,
            },
            Magnetorquer {
                id: "MTQ_Z".to_string(),
                axis: Vector3::new(0.0, 0.0, 1.0),
                max_dipole,
                dipole: 0.0,
                resistance: 10.0,
                power: 0.0,
                status: ActuatorStatus::Nominal,
            },
        ]
    }
    
    /// Calculate magnetic torque given Earth's magnetic field
    pub fn calculate_torque(&self, magnetic_field: &Vector3<f64>) -> Vector3<f64> {
        let dipole_vector = self.axis * self.dipole;
        dipole_vector.cross(magnetic_field) * 1e-9 // Convert nT to T
    }
    
    /// Size magnetorquer for given maximum torque requirement
    pub fn size_for_torque(torque_required: f64, magnetic_field_strength: f64) -> f64 {
        // M = T / B (dipole moment = torque / magnetic field)
        torque_required / (magnetic_field_strength * 1e-9)
    }
}

impl AttitudeDynamics {
    /// Create new attitude dynamics simulator
    pub fn new(acs: AttitudeControlSystem, initial_state: AttitudeState) -> Self {
        Self {
            acs,
            state: initial_state,
            orbital_state: OrbitalState {
                position: Vector3::new(7000.0, 0.0, 0.0),
                velocity: Vector3::new(0.0, 7.5, 0.0),
                altitude: 600.0,
                period: 5760.0,
            },
            environment: EnvironmentParameters {
                magnetic_field: Vector3::new(20000.0, 5000.0, -40000.0), // nT
                sun_vector: Vector3::new(1.0, 0.0, 0.0),
                density: 1e-12,
                solar_flux: 1367.0,
            },
        }
    }
    
    /// Propagate attitude dynamics by one time step
    pub fn propagate(&mut self, dt: f64) -> Result<(), String> {
        // Calculate external torques
        let external_torques = self.calculate_external_torques();
        
        // Calculate control torques
        let control_torques = self.calculate_control_torques();
        
        // Total torque
        let total_torque = external_torques.total + control_torques;
        
        // Euler's equation: I⋅ω̇ + ω×(I⋅ω) = T
        let inertia = &self.acs.inertia.inertia_matrix;
        let omega = &self.state.angular_velocity;
        let inertia_omega = inertia * omega;
        let gyroscopic_torque = omega.cross(&inertia_omega);
        
        // Solve for angular acceleration: ω̇ = I⁻¹⋅(T - ω×(I⋅ω))
        let angular_acceleration = match inertia.try_inverse() {
            Some(inertia_inv) => inertia_inv * (total_torque - gyroscopic_torque),
            None => return Err("Singular inertia matrix".to_string()),
        };
        
        // Integrate angular velocity
        self.state.angular_velocity += angular_acceleration * dt;
        self.state.angular_acceleration = angular_acceleration;
        
        // Integrate quaternion using angular velocity
        self.integrate_quaternion(dt);
        
        // Update Euler angles
        self.state.update_euler_angles();
        
        // Update time
        self.state.time += dt;
        
        // Update actuator states
        self.update_actuators(dt);
        
        Ok(())
    }
    
    /// Integrate quaternion kinematics
    fn integrate_quaternion(&mut self, dt: f64) {
        let omega = &self.state.angular_velocity;
        let omega_magnitude = omega.magnitude();
        
        if omega_magnitude < 1e-12 {
            // No rotation
            return;
        }
        
        // Quaternion propagation: q̇ = 0.5 * Ω(ω) * q
        let omega_quat = Quaternion::new(0.0, omega.x, omega.y, omega.z);
        let q_dot = omega_quat * self.state.quaternion.quaternion() * 0.5;
        
        // First-order integration
        let new_quat = self.state.quaternion.quaternion() + q_dot * dt;
        
        // Normalize to maintain unit quaternion
        let normalized = new_quat.normalize();
        self.state.quaternion = UnitQuaternion::from_quaternion(normalized);
    }
    
    /// Calculate all external torques acting on spacecraft
    fn calculate_external_torques(&self) -> ExternalTorques {
        let gravity_gradient = self.calculate_gravity_gradient_torque();
        let atmospheric = self.calculate_atmospheric_torque();
        let solar_radiation = self.calculate_solar_radiation_torque();
        let magnetic = self.calculate_magnetic_torque();
        
        let total = gravity_gradient + atmospheric + solar_radiation + magnetic;
        
        ExternalTorques {
            gravity_gradient,
            atmospheric,
            solar_radiation,
            magnetic,
            control: Vector3::zeros(), // Will be calculated separately
            total,
        }
    }
    
    /// Calculate gravity gradient torque
    fn calculate_gravity_gradient_torque(&self) -> Vector3<f64> {
        let r = &self.orbital_state.position;
        let r_magnitude = r.magnitude();
        
        if r_magnitude < 1e-6 {
            return Vector3::zeros();
        }
        
        // Transform position to body frame
        let r_body = self.state.dcm_inertial_to_body() * r;
        let r_body_unit = r_body / r_body.magnitude();
        
        // Gravity gradient torque: T_gg = 3μ/r³ * r̂ × (I⋅r̂)
        let mu = MU_EARTH * 1e9; // Convert to m³/s²
        let r_m = r_magnitude * 1000.0; // Convert to meters
        
        let inertia_r = self.acs.inertia.inertia_matrix * r_body_unit;
        let torque_coefficient = 3.0 * mu / (r_m * r_m * r_m);
        
        r_body_unit.cross(&inertia_r) * torque_coefficient
    }
    
    /// Calculate atmospheric torque (simplified)
    fn calculate_atmospheric_torque(&self) -> Vector3<f64> {
        if self.orbital_state.altitude > 600.0 {
            return Vector3::zeros(); // Negligible at high altitudes
        }
        
        // Very simplified model - would need detailed spacecraft geometry
        let velocity_body = self.state.dcm_inertial_to_body() * self.orbital_state.velocity;
        let dynamic_pressure = 0.5 * self.environment.density * velocity_body.magnitude_squared();
        
        // Assume small aerodynamic imbalance
        let aero_imbalance = Vector3::new(0.1, 0.05, 0.02); // m (center of pressure offset)
        let drag_coefficient = 2.2; // Typical for satellite
        let reference_area = 2.0; // m²
        
        let force = velocity_body.normalize() * dynamic_pressure * drag_coefficient * reference_area;
        aero_imbalance.cross(&force)
    }
    
    /// Calculate solar radiation pressure torque
    fn calculate_solar_radiation_torque(&self) -> Vector3<f64> {
        let c = 299792458.0; // Speed of light (m/s)
        let solar_pressure = self.environment.solar_flux / c;
        
        // Transform sun vector to body frame
        let sun_body = self.state.dcm_inertial_to_body() * self.environment.sun_vector;
        
        // Assume solar panel with offset from center of mass
        let solar_panel_area = 4.0; // m²
        let solar_panel_offset = Vector3::new(0.0, 0.0, 0.5); // m
        let reflectivity = 0.6; // Typical for solar cells
        
        let force = sun_body * solar_pressure * solar_panel_area * (1.0 + reflectivity);
        solar_panel_offset.cross(&force)
    }
    
    /// Calculate magnetic torque from residual magnetic dipole
    fn calculate_magnetic_torque(&self) -> Vector3<f64> {
        // Residual magnetic dipole from spacecraft electronics/materials
        let residual_dipole = Vector3::new(0.1, 0.05, 0.2); // A⋅m²
        
        // Transform magnetic field to body frame
        let b_body = self.state.dcm_inertial_to_body() * self.environment.magnetic_field;
        
        residual_dipole.cross(&b_body) * 1e-9 // Convert nT to T
    }
    
    /// Calculate control torques from actuators
    fn calculate_control_torques(&self) -> Vector3<f64> {
        // This would implement the attitude control law
        // For now, return zero (open-loop dynamics)
        Vector3::zeros()
    }
    
    /// Update actuator states
    fn update_actuators(&mut self, dt: f64) {
        // Update reaction wheels
        for wheel in &mut self.acs.reaction_wheels {
            // Simple momentum accumulation from external torques
            let torque_on_wheel = 0.001; // N⋅m - small disturbance
            wheel.update(torque_on_wheel, dt);
        }
        
        // Update magnetorquers (would implement control law)
        for mtq in &mut self.acs.magnetorquers {
            mtq.dipole = 0.0; // No control for now
            mtq.power = mtq.dipole.abs() * mtq.dipole.abs() / mtq.resistance;
        }
    }
}

/// Utility functions for attitude representations

/// Convert quaternion to Euler angles (roll, pitch, yaw)
pub fn quaternion_to_euler_angles(q: &UnitQuaternion<f64>) -> Vector3<f64> {
    let q = q.quaternion();
    let w = q.w;
    let x = q.i;
    let y = q.j;
    let z = q.k;
    
    // Roll (x-axis rotation)
    let roll = (2.0 * (w * x + y * z)).atan2(1.0 - 2.0 * (x * x + y * y));
    
    // Pitch (y-axis rotation)  
    let sin_pitch = 2.0 * (w * y - z * x);
    let pitch = if sin_pitch.abs() >= 1.0 {
        FRAC_PI_2.copysign(sin_pitch)
    } else {
        sin_pitch.asin()
    };
    
    // Yaw (z-axis rotation)
    let yaw = (2.0 * (w * z + x * y)).atan2(1.0 - 2.0 * (y * y + z * z));
    
    Vector3::new(roll, pitch, yaw)
}

/// Convert Euler angles to quaternion
pub fn euler_angles_to_quaternion(euler: &Vector3<f64>) -> UnitQuaternion<f64> {
    UnitQuaternion::from_euler_angles(euler.x, euler.y, euler.z)
}

/// Pre-configured spacecraft for common missions
pub mod spacecraft {
    use super::*;
    
    /// CubeSat (3U) configuration
    pub fn cubesat_3u() -> AttitudeControlSystem {
        let inertia = SpacecraftInertia::box_satellite(0.34, 0.10, 0.10, 4.0);
        let reaction_wheels = vec![
            ReactionWheel::standard_wheel(Vector3::new(1.0, 0.0, 0.0), 3000.0, 0.001),
            ReactionWheel::standard_wheel(Vector3::new(0.0, 1.0, 0.0), 3000.0, 0.001),
            ReactionWheel::standard_wheel(Vector3::new(0.0, 0.0, 1.0), 3000.0, 0.001),
        ];
        let magnetorquers = Magnetorquer::three_axis_system(0.2);
        
        AttitudeControlSystem {
            inertia,
            reaction_wheels,
            magnetorquers,
            thrusters: vec![],
            control_params: ControlParameters {
                kp_attitude: 0.1,
                kd_rate: 0.05,
                ki_integral: 0.001,
                max_control_torque: 0.01,
                max_slew_rate: 0.1,
                attitude_deadband: 0.01,
                rate_deadband: 0.001,
            },
        }
    }
    
    /// Large satellite with CMGs
    pub fn large_satellite() -> AttitudeControlSystem {
        let inertia = SpacecraftInertia::from_principal_moments(5000.0, 4000.0, 3000.0, 2000.0);
        let reaction_wheels = ReactionWheel::pyramid_configuration();
        let magnetorquers = Magnetorquer::three_axis_system(100.0);
        
        AttitudeControlSystem {
            inertia,
            reaction_wheels,
            magnetorquers,
            thrusters: vec![],
            control_params: ControlParameters {
                kp_attitude: 1.0,
                kd_rate: 0.5,
                ki_integral: 0.01,
                max_control_torque: 50.0,
                max_slew_rate: 0.05,
                attitude_deadband: 0.001,
                rate_deadband: 0.0001,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::spacecraft::*;
    
    #[test]
    fn test_attitude_state_creation() {
        let q = UnitQuaternion::identity();
        let omega = Vector3::new(0.1, 0.0, 0.0);
        let state = AttitudeState::new(q, omega);
        
        assert_eq!(state.angular_velocity, omega);
        assert!((state.euler_angles.magnitude() - 0.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_quaternion_euler_conversion() {
        let euler = Vector3::new(0.1, 0.2, 0.3);
        let quat = euler_angles_to_quaternion(&euler);
        let euler_back = quaternion_to_euler_angles(&quat);
        
        assert!((euler - euler_back).magnitude() < 1e-10);
    }
    
    #[test]
    fn test_spacecraft_inertia_box() {
        let inertia = SpacecraftInertia::box_satellite(2.0, 1.0, 0.5, 10.0);
        
        assert!(inertia.principal_moments.x > 0.0);
        assert!(inertia.principal_moments.y > 0.0);
        assert!(inertia.principal_moments.z > 0.0);
        assert_eq!(inertia.mass, 10.0);
    }
    
    #[test]
    fn test_reaction_wheel_pyramid() {
        let wheels = ReactionWheel::pyramid_configuration();
        
        assert_eq!(wheels.len(), 4);
        
        // Check that axes are not parallel (should form pyramid)
        let dot_product = wheels[0].axis.dot(&wheels[1].axis);
        assert!(dot_product < 0.9); // Not parallel
    }
    
    #[test]
    fn test_magnetorquer_system() {
        let mtqs = Magnetorquer::three_axis_system(1.0);
        
        assert_eq!(mtqs.len(), 3);
        assert_eq!(mtqs[0].axis, Vector3::new(1.0, 0.0, 0.0));
        assert_eq!(mtqs[1].axis, Vector3::new(0.0, 1.0, 0.0));
        assert_eq!(mtqs[2].axis, Vector3::new(0.0, 0.0, 1.0));
    }
    
    #[test]
    fn test_magnetorquer_torque_calculation() {
        let mut mtq = Magnetorquer::three_axis_system(1.0)[0].clone();
        mtq.dipole = 1.0; // 1 A⋅m²
        
        let b_field = Vector3::new(0.0, 0.0, 50000.0); // 50,000 nT in z direction
        let torque = mtq.calculate_torque(&b_field);
        
        assert!(torque.magnitude() > 0.0);
        // dipole=1 A⋅m², B=(0,0,50000nT)=(0,0,5e-5T), τ = m×B = (0, -5e-5, 0)
        assert!((torque - Vector3::new(0.0, -5e-5, 0.0)).magnitude() < 1e-8);
    }
    
    #[test]
    fn test_gravity_gradient_torque() {
        let acs = cubesat_3u();
        let initial_state = AttitudeState::new(UnitQuaternion::identity(), Vector3::zeros());
        let mut dynamics = AttitudeDynamics::new(acs, initial_state);
        
        // Set orbital position  
        dynamics.orbital_state.position = Vector3::new(7000.0, 0.0, 0.0);
        
        let gg_torque = dynamics.calculate_gravity_gradient_torque();
        
        // Should be small for symmetric spacecraft in aligned attitude
        assert!(gg_torque.magnitude() < 1e-6);
    }
    
    #[test]
    fn test_attitude_propagation() {
        let acs = cubesat_3u();
        let initial_omega = Vector3::new(0.1, 0.0, 0.0);
        let initial_state = AttitudeState::new(UnitQuaternion::identity(), initial_omega);
        let mut dynamics = AttitudeDynamics::new(acs, initial_state);
        
        let dt = 0.1;
        dynamics.propagate(dt).unwrap();
        
        // Should have rotated
        assert!(dynamics.state.angular_velocity.magnitude() > 0.0);
        assert!(dynamics.state.euler_angles.magnitude() > 0.0);
        assert_eq!(dynamics.state.time, dt);
    }
    
    #[test]
    fn test_reaction_wheel_update() {
        let mut wheel = ReactionWheel::standard_wheel(Vector3::new(1.0, 0.0, 0.0), 1000.0, 0.1);
        
        wheel.update(0.05, 1.0); // 0.05 N⋅m for 1 second
        
        assert!(wheel.speed > 0.0);
        assert!(wheel.power > 0.0);
        assert!(wheel.speed < wheel.max_speed);
    }
    
    #[test]
    fn test_magnetorquer_sizing() {
        let torque_required = 1e-3; // 1 mN⋅m
        let b_field_strength = 50000.0; // 50,000 nT
        
        let dipole_needed = Magnetorquer::size_for_torque(torque_required, b_field_strength);
        
        assert!(dipole_needed > 0.0);
        // dipole = T/B = 1e-3 / 5e-5 = 20 A⋅m² — reasonable for spacecraft
        assert!(dipole_needed > 0.0);
        assert!(dipole_needed < 100.0); // Reasonable dipole moment
    }
    
    #[test]  
    fn test_dcm_conversion() {
        let euler = Vector3::new(0.1, 0.2, 0.3);
        let quat = euler_angles_to_quaternion(&euler);
        let state = AttitudeState::new(quat, Vector3::zeros());
        
        let dcm_ib = state.dcm_inertial_to_body();
        let dcm_bi = state.dcm_body_to_inertial();
        
        let identity = dcm_ib * dcm_bi;
        let identity_error = (identity - Matrix3::identity()).norm();
        
        assert!(identity_error < 1e-10);
    }
    
    #[test]
    fn test_cubesat_configuration() {
        let acs = cubesat_3u();
        
        assert_eq!(acs.reaction_wheels.len(), 3);
        assert_eq!(acs.magnetorquers.len(), 3);
        assert!(acs.inertia.mass < 10.0); // Reasonable CubeSat mass
        assert!(acs.control_params.max_control_torque < 0.1);
    }
    
    #[test]
    fn test_large_satellite_configuration() {
        let acs = large_satellite();
        
        assert_eq!(acs.reaction_wheels.len(), 4); // Pyramid configuration
        assert!(acs.inertia.mass > 1000.0); // Large satellite
        assert!(acs.control_params.max_control_torque > 10.0);
    }
}