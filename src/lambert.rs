//! Lambert's problem solver: Given two position vectors and time of flight,
//! find the velocity vectors that connect them.
//! 
//! Implementation uses the Universal Variable method from Vallado.
//! This is THE algorithm for interplanetary transfer calculations.
//!
//! References:
//! - Vallado "Fundamentals of Astrodynamics and Applications", Chapter 5
//! - Bate, Mueller & White "Fundamentals of Astrodynamics", Chapter 5
//! - Battin "Mathematics and Methods of Astrodynamics", Chapter 7

use nalgebra::Vector3;
use crate::constants::PI;

/// Solution to Lambert's problem
#[derive(Debug, Clone)]
pub struct LambertSolution {
    /// Initial velocity vector (km/s)
    pub v1: Vector3<f64>,
    /// Final velocity vector (km/s)
    pub v2: Vector3<f64>,
    /// Transfer angle (radians)
    pub transfer_angle: f64,
    /// Orbital energy (km²/s²)
    pub energy: f64,
    /// Semi-major axis of transfer orbit (km)
    pub a: f64,
    /// Eccentricity of transfer orbit
    pub e: f64,
    /// Number of iterations to converge
    pub iterations: u32,
}

/// Lambert's problem solver using Universal Variables
/// 
/// Solves for the orbit that connects two position vectors in a given time.
/// This is the fundamental algorithm for:
/// - Interplanetary transfers (Earth → Mars)
/// - Orbital rendezvous (chase spacecraft → target)
/// - Mission planning (launch windows, pork chop plots)
/// 
/// # Arguments
/// * `r1` - Initial position vector (km)
/// * `r2` - Final position vector (km)
/// * `dt` - Time of flight (seconds)
/// * `mu` - Gravitational parameter (km³/s²)
/// * `prograde` - True for prograde (short way), false for retrograde
/// 
/// # Returns
/// Result containing Lambert solution or error message
/// 
/// # Example
/// ```
/// use nalgebra::Vector3;
/// use rocketlab::{lambert::solve_lambert, constants::MU_EARTH};
/// 
/// let r1 = Vector3::new(7000.0, 0.0, 0.0);  // LEO position
/// let r2 = Vector3::new(0.0, 8000.0, 0.0);  // Target position
/// let dt = 3600.0;  // 1 hour transfer
/// 
/// match solve_lambert(r1, r2, dt, MU_EARTH, true) {
///     Ok(solution) => println!("ΔV = {:.3} km/s", solution.v1.magnitude()),
///     Err(e) => println!("Lambert solver failed: {}", e),
/// }
/// ```
pub fn solve_lambert(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    dt: f64,
    mu: f64,
    prograde: bool,
) -> Result<LambertSolution, String> {
    const MAX_ITERATIONS: u32 = 50;
    const TOLERANCE: f64 = 1e-8;
    
    // Position magnitudes
    let r1_mag = r1.magnitude();
    let r2_mag = r2.magnitude();
    
    if r1_mag < 1e-6 || r2_mag < 1e-6 {
        return Err("Position vectors too small".to_string());
    }
    
    if dt <= 0.0 {
        return Err("Time of flight must be positive".to_string());
    }
    
    // Transfer angle
    let cos_dnu = r1.dot(&r2) / (r1_mag * r2_mag);
    let cos_dnu = cos_dnu.clamp(-1.0, 1.0); // Numerical safety
    let dnu = cos_dnu.acos();
    
    // Determine actual transfer angle based on direction
    let cross = r1.cross(&r2);
    let transfer_angle = if prograde {
        if cross.z >= 0.0 { dnu } else { 2.0 * PI - dnu }
    } else {
        if cross.z >= 0.0 { 2.0 * PI - dnu } else { dnu }
    };
    
    // Lambert geometry parameters
    let c = (r1 - r2).magnitude();
    let s = (r1_mag + r2_mag + c) / 2.0;
    let am = s / 2.0; // Minimum energy semi-major axis
    let lambda = ((r1_mag * r2_mag).sqrt() * (1.0 + cos_dnu)) / s;
    
    // Initial guess for universal anomaly
    let mut x = if transfer_angle <= PI {
        // Short way transfer - parabolic initial guess
        0.0
    } else {
        // Long way transfer - hyperbolic initial guess  
        (2.0 * transfer_angle / dt.sqrt()).ln()
    };
    
    let mut iterations = 0;
    let mut dt_calculated;
    
    // Newton-Raphson iteration on universal anomaly
    loop {
        iterations += 1;
        if iterations > MAX_ITERATIONS {
            return Err("Lambert solver failed to converge".to_string());
        }
        
        let x2 = x * x;
        let x3 = x2 * x;
        
        // Stumpff functions
        let (c2, c3) = stumpff_functions(x2);
        
        // Current semi-major axis
        let a = am / (1.0 - x2 * c2 / s);
        
        if a <= 0.0 && iterations > 5 {
            return Err("Negative semi-major axis - no solution".to_string());
        }
        
        // Universal anomaly functions
        let sqrt_mu_a = (mu * a.abs()).sqrt();
        let dt_universal = if a > 0.0 {
            // Elliptical case
            (r1_mag * r2_mag * (x3 * c3 + x) / sqrt_mu_a) / lambda
        } else {
            // Hyperbolic case
            (r1_mag * r2_mag * (x3 * c3 - x) / sqrt_mu_a) / lambda
        };
        
        dt_calculated = s * (x3 * c3 + x) / sqrt_mu_a + dt_universal;
        
        let dt_error = dt_calculated - dt;
        
        if dt_error.abs() < TOLERANCE {
            // Converged! Calculate velocities
            let f = 1.0 - a * (1.0 - cos_dnu) / r1_mag;
            let g = r1_mag * r2_mag * transfer_angle.sin() / sqrt_mu_a;
            let fdot = sqrt_mu_a * (transfer_angle.cos() - 1.0 + (1.0 - cos_dnu) * a / r1_mag) / 
                       (r1_mag * r2_mag * transfer_angle.sin());
            let gdot = 1.0 - a * (1.0 - cos_dnu) / r2_mag;
            
            // Velocity vectors
            let v1 = (r2 - f * r1) / g;
            let v2 = fdot * r1 + gdot * v1;
            
            // Orbital parameters
            let energy = -mu / (2.0 * a);
            let h = r1.cross(&v1);
            let e_vec = (v1.cross(&h) / mu) - (r1 / r1_mag);
            let e = e_vec.magnitude();
            
            return Ok(LambertSolution {
                v1,
                v2,
                transfer_angle,
                energy,
                a,
                e,
                iterations,
            });
        }
        
        // Newton-Raphson update
        let dt_derivative = dt_derivative_dx(x, a, r1_mag, r2_mag, s, lambda, mu);
        if dt_derivative.abs() < 1e-12 {
            return Err("Derivative too small - cannot converge".to_string());
        }
        
        x -= dt_error / dt_derivative;
        
        // Prevent divergence
        x = x.clamp(-4.0 * PI, 4.0 * PI);
    }
}

/// Stumpff functions C₂(z) and C₃(z) used in universal variable formulation
fn stumpff_functions(z: f64) -> (f64, f64) {
    if z > 1e-6 {
        // Elliptical case: z > 0
        let sqrt_z = z.sqrt();
        let c2 = (1.0 - sqrt_z.cos()) / z;
        let c3 = (sqrt_z - sqrt_z.sin()) / (z * sqrt_z);
        (c2, c3)
    } else if z < -1e-6 {
        // Hyperbolic case: z < 0
        let sqrt_minus_z = (-z).sqrt();
        let c2 = (sqrt_minus_z.cosh() - 1.0) / (-z);
        let c3 = (sqrt_minus_z.sinh() - sqrt_minus_z) / ((-z) * sqrt_minus_z);
        (c2, c3)
    } else {
        // Parabolic case: z ≈ 0, use series expansion
        let c2 = 1.0 / 2.0 - z / 24.0 + z * z / 720.0;
        let c3 = 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
        (c2, c3)
    }
}

/// Derivative of time of flight with respect to universal anomaly x
fn dt_derivative_dx(x: f64, a: f64, r1: f64, r2: f64, s: f64, lambda: f64, mu: f64) -> f64 {
    let x2 = x * x;
    let (c2, c3) = stumpff_functions(x2);
    let sqrt_mu_a = (mu * a.abs()).sqrt();
    
    // Derivative of Stumpff functions
    let (dc2, dc3) = if x2 > 1e-6 {
        let sqrt_z = x2.sqrt();
        let dc2 = (sqrt_z.sin() - sqrt_z) / (x2 * sqrt_z);
        let dc3 = (sqrt_z.cos() - 1.0 + sqrt_z * sqrt_z / 2.0) / (x2 * x2);
        (dc2, dc3)
    } else if x2 < -1e-6 {
        let sqrt_minus_z = (-x2).sqrt();
        let dc2 = (sqrt_minus_z - sqrt_minus_z.sinh()) / (x2 * sqrt_minus_z);
        let dc3 = (sqrt_minus_z.cosh() - 1.0 - x2 / 2.0) / (x2 * x2);
        (dc2, dc3)
    } else {
        let dc2 = -1.0 / 12.0 + x2 / 360.0;
        let dc3 = -1.0 / 60.0 + x2 / 2520.0;
        (dc2, dc3)
    };
    
    // Complex derivative calculation (simplified)
    let term1 = s * (3.0 * c3 * x2 + 1.0) / sqrt_mu_a;
    let term2 = r1 * r2 * (3.0 * dc3 * x2 + c3) / (sqrt_mu_a * lambda);
    
    term1 + term2
}

/// Calculate ΔV required for Lambert transfer
/// 
/// Returns the total velocity change needed at departure and arrival
pub fn lambert_delta_v(
    r1: Vector3<f64>,
    v1_initial: Vector3<f64>,
    r2: Vector3<f64>, 
    v2_target: Vector3<f64>,
    dt: f64,
    mu: f64,
    prograde: bool,
) -> Result<(f64, f64, f64), String> {
    let solution = solve_lambert(r1, r2, dt, mu, prograde)?;
    
    let dv1 = (solution.v1 - v1_initial).magnitude();
    let dv2 = (v2_target - solution.v2).magnitude();
    let dv_total = dv1 + dv2;
    
    Ok((dv1, dv2, dv_total))
}

/// Multi-revolution Lambert solver
/// 
/// Solves Lambert's problem for transfers that complete multiple orbits
/// before reaching the target. Used for minimum-energy transfers.
pub fn solve_lambert_multirev(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    dt: f64,
    mu: f64,
    revolutions: i32,
) -> Result<Vec<LambertSolution>, String> {
    let mut solutions = Vec::new();
    
    // Try both prograde and retrograde for each revolution count
    for rev in 0..=revolutions.abs() {
        for &prograde in &[true, false] {
            // Adjust time for multiple revolutions
            let period_est = 2.0 * PI * ((r1.magnitude() + r2.magnitude()) / 4.0).powf(1.5) / mu.sqrt();
            let dt_adjusted = dt + rev as f64 * period_est;
            
            if let Ok(solution) = solve_lambert(r1, r2, dt_adjusted, mu, prograde) {
                solutions.push(solution);
            }
        }
    }
    
    if solutions.is_empty() {
        Err("No multi-revolution solutions found".to_string())
    } else {
        Ok(solutions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MU_EARTH, MU_SUN};
    use nalgebra::Vector3;
    
    #[test]
    fn test_lambert_circular_orbit() {
        // Test case: quarter orbit transfer in circular orbit
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, 7000.0, 0.0);
        let dt = 0.25 * 2.0 * PI * (7000.0_f64.powf(3.0) / MU_EARTH).sqrt(); // Quarter period
        
        let solution = solve_lambert(r1, r2, dt, MU_EARTH, true).unwrap();
        
        // For circular orbit, semi-major axis should be radius
        assert!((solution.a - 7000.0).abs() < 10.0, "Semi-major axis mismatch");
        
        // Low eccentricity expected for circular orbit transfer
        assert!(solution.e < 0.1, "Eccentricity too high for circular transfer");
    }
    
    #[test]
    fn test_lambert_hohmann_transfer() {
        // Hohmann transfer from 7000 km to 10000 km circular orbits
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(-10000.0, 0.0, 0.0);
        
        // Hohmann transfer time
        let a_transfer: f64 = (7000.0 + 10000.0) / 2.0;
        let dt = PI * (a_transfer.powf(3.0) / MU_EARTH).sqrt();
        
        let solution = solve_lambert(r1, r2, dt, MU_EARTH, true).unwrap();
        
        // Check transfer orbit semi-major axis
        assert!((solution.a - a_transfer).abs() < 50.0, "Hohmann transfer a mismatch");
        
        // Check transfer orbit eccentricity
        let expected_e = (10000.0 - 7000.0) / (10000.0 + 7000.0);
        assert!((solution.e - expected_e).abs() < 0.01, "Hohmann transfer e mismatch");
    }
    
    #[test]
    fn test_lambert_vallado_example() {
        // Vallado Example 5-2: Earth to Mars transfer
        // Simplified positions (actual would use ephemeris)
        let r1 = Vector3::new(149597870.7, 0.0, 0.0); // Earth orbit (1 AU)
        let r2 = Vector3::new(-227939366.0, 0.0, 0.0); // Mars orbit (1.52 AU, opposite side)
        
        let dt = 207.0 * 24.0 * 3600.0; // 207 days transfer
        
        let solution = solve_lambert(r1, r2, dt, MU_SUN, true).unwrap();
        
        // Transfer should be elliptical with reasonable semi-major axis
        assert!(solution.a > 149597870.7, "Transfer orbit too small");
        assert!(solution.a < 300000000.0, "Transfer orbit too large");
        assert!(solution.e > 0.1, "Transfer should be elliptical");
    }
    
    #[test]
    fn test_lambert_edge_cases() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(8000.0, 0.0, 0.0);
        
        // Zero time should fail
        assert!(solve_lambert(r1, r2, 0.0, MU_EARTH, true).is_err());
        
        // Negative time should fail
        assert!(solve_lambert(r1, r2, -100.0, MU_EARTH, true).is_err());
        
        // Zero position should fail
        let r_zero = Vector3::new(0.0, 0.0, 0.0);
        assert!(solve_lambert(r_zero, r2, 3600.0, MU_EARTH, true).is_err());
    }
    
    #[test]
    fn test_stumpff_functions() {
        // Test Stumpff functions for various cases
        
        // Parabolic case (z ≈ 0)
        let (c2, c3) = stumpff_functions(1e-10);
        assert!((c2 - 0.5).abs() < 1e-6, "Stumpff C₂ parabolic case");
        assert!((c3 - 1.0/6.0).abs() < 1e-6, "Stumpff C₃ parabolic case");
        
        // Elliptical case (z > 0)
        let z = 0.5;
        let (c2, c3) = stumpff_functions(z);
        let expected_c2 = (1.0 - z.sqrt().cos()) / z;
        assert!((c2 - expected_c2).abs() < 1e-10, "Stumpff C₂ elliptical case");
        
        // Hyperbolic case (z < 0)
        let z = -0.5;
        let (c2, c3) = stumpff_functions(z);
        let expected_c2 = ((-z).sqrt().cosh() - 1.0) / (-z);
        assert!((c2 - expected_c2).abs() < 1e-10, "Stumpff C₂ hyperbolic case");
    }
    
    #[test]
    fn test_delta_v_calculation() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, 7000.0, 0.0);
        let v1_initial = Vector3::new(0.0, 7.546, 0.0); // Circular velocity
        let v2_target = Vector3::new(-7.546, 0.0, 0.0); // Circular velocity
        let dt = 1800.0; // 30 minutes
        
        let (dv1, dv2, dv_total) = lambert_delta_v(
            r1, v1_initial, r2, v2_target, dt, MU_EARTH, true
        ).unwrap();
        
        assert!(dv1 > 0.0, "Initial ΔV should be positive");
        assert!(dv2 > 0.0, "Final ΔV should be positive");
        assert!((dv_total - dv1 - dv2).abs() < 1e-10, "Total ΔV calculation");
    }
}