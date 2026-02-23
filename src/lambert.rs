//! Lambert's problem solver: Given two position vectors and time of flight,
//! find the velocity vectors that connect them.
//! 
//! Implementation uses the f and g series with universal variable formulation
//! from Bate, Mueller & White and Curtis.
//!
//! References:
//! - Curtis "Orbital Mechanics for Engineering Students", Algorithm 5.2
//! - Bate, Mueller & White "Fundamentals of Astrodynamics", Chapter 5
//! - Vallado "Fundamentals of Astrodynamics and Applications", Chapter 5

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

/// Stumpff functions C(z) and S(z)
fn stumpff_c(z: f64) -> f64 {
    if z > 1e-6 {
        (1.0 - z.sqrt().cos()) / z
    } else if z < -1e-6 {
        ((-z).sqrt().cosh() - 1.0) / (-z)
    } else {
        1.0 / 2.0 - z / 24.0 + z * z / 720.0
    }
}

fn stumpff_s(z: f64) -> f64 {
    if z > 1e-6 {
        let sz = z.sqrt();
        (sz - sz.sin()) / (z * sz)
    } else if z < -1e-6 {
        let sz = (-z).sqrt();
        (sz.sinh() - sz) / ((-z) * sz)
    } else {
        1.0 / 6.0 - z / 120.0 + z * z / 5040.0
    }
}

/// Lambert's problem solver using universal variable method
/// 
/// Algorithm from Curtis "Orbital Mechanics for Engineering Students" Algorithm 5.2
/// Uses Newton-Raphson iteration on the universal anomaly z.
///
/// # Arguments
/// * `r1` - Initial position vector (km)
/// * `r2` - Final position vector (km)  
/// * `dt` - Time of flight (seconds)
/// * `mu` - Gravitational parameter (km³/s²)
/// * `prograde` - True for prograde (short way), false for retrograde
pub fn solve_lambert(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    dt: f64,
    mu: f64,
    prograde: bool,
) -> Result<LambertSolution, String> {
    const MAX_ITERATIONS: u32 = 5000;
    const TOLERANCE: f64 = 1e-8;
    
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
    let cos_dnu = cos_dnu.clamp(-1.0, 1.0);
    
    let cross = r1.cross(&r2);
    let transfer_angle = {
        let dnu = cos_dnu.acos();
        if prograde {
            if cross.z >= 0.0 { dnu } else { 2.0 * PI - dnu }
        } else {
            if cross.z >= 0.0 { 2.0 * PI - dnu } else { dnu }
        }
    };
    
    // A parameter (Curtis Eq 5.35)
    let sin_dnu = transfer_angle.sin();
    let cos_dnu_eff = transfer_angle.cos();
    if (1.0 - cos_dnu_eff).abs() < 1e-12 {
        // Transfer angle ≈ 0° or 360°: degenerate
        return Err("Transfer angle too close to 0 or 2π".to_string());
    }
    if sin_dnu.abs() < 1e-10 {
        // 180° transfer: orbit plane is undetermined for coplanar vectors
        return Err("180° transfer: orbit plane undetermined, add out-of-plane component".to_string());
    }
    let a_param = sin_dnu * (r1_mag * r2_mag / (1.0 - cos_dnu_eff)).sqrt();
    
    // Function F(z) = ratio * S(z)^(3/2) + A * sqrt(y(z)) - sqrt(mu) * dt
    // y(z) = r1 + r2 + A * (z*S(z) - 1) / sqrt(C(z))
    
    let y = |z: f64| -> f64 {
        let c = stumpff_c(z);
        let s = stumpff_s(z);
        r1_mag + r2_mag + a_param * (z * s - 1.0) / c.sqrt()
    };
    
    // Find z lower bound where y > 0
    // For z < 0 (hyperbolic), y can go negative
    let mut z_low = -4.0 * PI * PI; // Start very hyperbolic
    while y(z_low) < 0.0 {
        z_low += 0.1;
        if z_low > 50.0 {
            return Err("Cannot find valid starting z".to_string());
        }
    }
    
    // F(z) and F'(z) for Newton-Raphson
    let _f_and_fp = |z: f64| -> (f64, f64) {
        let c = stumpff_c(z);
        let s = stumpff_s(z);
        let y_val = r1_mag + r2_mag + a_param * (z * s - 1.0) / c.sqrt();
        
        if y_val < 0.0 {
            return (f64::MAX, 1.0);
        }
        
        let x = (y_val / c).sqrt();
        let f_val = x.powi(3) * s + a_param * y_val.sqrt() - dt * mu.sqrt();
        
        // Derivative
        let f_prime = if z.abs() > 1e-6 {
            let dy_dz = if z.abs() > 1e-6 {
                let _term1 = a_param / c.sqrt();
                let _dc_dz = if z > 1e-6 {
                    (1.0 / (2.0 * z)) * (1.0 - z * s - 2.0 * c)
                } else if z < -1e-6 {
                    (1.0 / (2.0 * z)) * (1.0 - z * s - 2.0 * c)
                } else {
                    -1.0 / 12.0
                };
                // Simplified: use finite difference for robustness
                let eps = 1e-6;
                let y_plus = r1_mag + r2_mag + a_param * ((z + eps) * stumpff_s(z + eps) - 1.0) / stumpff_c(z + eps).sqrt();
                let y_minus = r1_mag + r2_mag + a_param * ((z - eps) * stumpff_s(z - eps) - 1.0) / stumpff_c(z - eps).sqrt();
                (y_plus - y_minus) / (2.0 * eps)
            } else {
                0.0
            };
            
            x.powi(3) * (s - 3.0 * s * dy_dz / (2.0 * y_val)) + 
            a_param / 8.0 * (3.0 * s * y_val.sqrt() / c + a_param / x)
        } else {
            // z ≈ 0: use series
            (2.0 / 40.0) * y_val.powf(1.5) + a_param / 8.0 * (y_val.sqrt() + a_param * (1.0 / (2.0 * y_val.sqrt())))
        };
        
        (f_val, f_prime)
    };
    
    // Use bisection to find the root of F(z) = 0
    // F(z) = x³S(z) + A√y - √μ·Δt  where x = √(y/C)
    let f_eval = |z: f64| -> f64 {
        let c = stumpff_c(z);
        let s = stumpff_s(z);
        let y_val = r1_mag + r2_mag + a_param * (z * s - 1.0) / c.sqrt();
        if y_val < 0.0 { return f64::MAX; }
        let x = (y_val / c).sqrt();
        x.powi(3) * s + a_param * y_val.sqrt() - dt * mu.sqrt()
    };
    
    let mut z_lo = z_low;
    let mut z_hi = 4.0 * PI * PI * 400.0;
    
    // Ensure bracket: find z_hi where F > 0
    while f_eval(z_hi) < 0.0 && z_hi < 1e8 {
        z_hi *= 2.0;
    }
    // Ensure z_lo where F < 0
    while f_eval(z_lo) > 0.0 && z_lo > -1e8 {
        z_lo -= 10.0;
    }
    
    let mut z = (z_lo + z_hi) / 2.0;
    let mut iterations = 0u32;
    
    loop {
        iterations += 1;
        if iterations > MAX_ITERATIONS {
            return Err("Lambert solver failed to converge".to_string());
        }
        
        let f_val = f_eval(z);
        
        if f_val.abs() < TOLERANCE || (z_hi - z_lo).abs() < 1e-14 {
            break;
        }
        
        if f_val < 0.0 {
            z_lo = z;
        } else {
            z_hi = z;
        }
        
        z = (z_lo + z_hi) / 2.0;
    }
    
    // Compute f, g, gdot from converged z
    let c = stumpff_c(z);
    let s = stumpff_s(z);
    let y_val = r1_mag + r2_mag + a_param * (z * s - 1.0) / c.sqrt();
    
    let f = 1.0 - y_val / r1_mag;
    let g = a_param * (y_val / mu).sqrt();
    let gdot = 1.0 - y_val / r2_mag;
    
    // Velocity vectors (Curtis Eq 5.46)
    let v1 = (r2 - f * r1) / g;
    let v2 = (gdot * r2 - r1) / g;
    
    // Orbital elements of transfer
    let energy = v1.magnitude_squared() / 2.0 - mu / r1_mag;
    let a = -mu / (2.0 * energy);
    let h = r1.cross(&v1);
    let e_vec = (v1.cross(&h) / mu) - (r1 / r1_mag);
    let e = e_vec.magnitude();
    
    Ok(LambertSolution {
        v1,
        v2,
        transfer_angle,
        energy,
        a,
        e,
        iterations,
    })
}

/// Calculate ΔV required for Lambert transfer
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
pub fn solve_lambert_multirev(
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    dt: f64,
    mu: f64,
    revolutions: i32,
) -> Result<Vec<LambertSolution>, String> {
    let mut solutions = Vec::new();
    
    for rev in 0..=revolutions.abs() {
        for &prograde in &[true, false] {
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
    fn test_stumpff_functions() {
        // Parabolic case (z ≈ 0)
        let c = stumpff_c(1e-10);
        let s = stumpff_s(1e-10);
        assert!((c - 0.5).abs() < 1e-6, "Stumpff C parabolic case");
        assert!((s - 1.0/6.0).abs() < 1e-6, "Stumpff S parabolic case");
        
        // Elliptical case (z > 0)
        let z = 0.5;
        let c = stumpff_c(z);
        let expected_c = (1.0 - z.sqrt().cos()) / z;
        assert!((c - expected_c).abs() < 1e-10);
        
        // Hyperbolic case (z < 0)
        let z = -0.5;
        let c = stumpff_c(z);
        let expected_c = ((-z).sqrt().cosh() - 1.0) / (-z);
        assert!((c - expected_c).abs() < 1e-10);
    }
    
    #[test]
    fn test_lambert_circular_orbit() {
        // Quarter orbit transfer in circular orbit
        let r: f64 = 7000.0;
        let r1 = Vector3::new(r, 0.0, 0.0);
        let r2 = Vector3::new(0.0, r, 0.0);
        let dt = 0.25 * 2.0 * PI * (r.powi(3) / MU_EARTH).sqrt();
        
        let solution = solve_lambert(r1, r2, dt, MU_EARTH, true).unwrap();
        
        // For circular orbit transfer, semi-major axis ≈ radius
        assert!((solution.a - r).abs() / r < 0.02, 
            "Semi-major axis mismatch: got {:.1}, expected ~{:.1}", solution.a, r);
        assert!(solution.e < 0.05, "Eccentricity too high: {}", solution.e);
    }
    
    #[test]
    fn test_lambert_hohmann_transfer() {
        // Hohmann transfer from 7000 km to 10000 km circular orbits
        let r1_r: f64 = 7000.0;
        let r2_r: f64 = 10000.0;
        let r1 = Vector3::new(r1_r, 0.0, 0.0);
        // Slightly off 180° to avoid degenerate coplanar case
        let r2 = Vector3::new(-r2_r, 10.0, 0.0); // Near-Hohmann transfer
        
        let a_transfer = (r1_r + r2_r) / 2.0;
        let dt = PI * (a_transfer.powi(3) / MU_EARTH).sqrt();
        
        let solution = solve_lambert(r1, r2, dt, MU_EARTH, true).unwrap();
        
        assert!((solution.a - a_transfer).abs() < 100.0, 
            "Hohmann a: got {:.1}, expected {:.1}", solution.a, a_transfer);
        
        let expected_e = (r2_r - r1_r) / (r1_r + r2_r);
        assert!((solution.e - expected_e).abs() < 0.02, 
            "Hohmann e: got {:.4}, expected {:.4}", solution.e, expected_e);
    }
    
    #[test]
    fn test_lambert_vallado_example() {
        // Earth to Mars Hohmann-like transfer
        let r1 = Vector3::new(149597870.7, 0.0, 0.0); // 1 AU
        let r2 = Vector3::new(-227939366.0, 1000.0, 0.0); // ~1.52 AU, near-opposite
        let dt = 207.0 * 24.0 * 3600.0; // 207 days
        
        let solution = solve_lambert(r1, r2, dt, MU_SUN, true).unwrap();
        
        assert!(solution.a > 149597870.7, "Transfer orbit too small: {}", solution.a);
        assert!(solution.a < 300000000.0, "Transfer orbit too large: {}", solution.a);
        assert!(solution.e > 0.1, "Transfer should be elliptical: e={}", solution.e);
    }
    
    #[test]
    fn test_lambert_edge_cases() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(8000.0, 0.0, 0.0);
        
        assert!(solve_lambert(r1, r2, 0.0, MU_EARTH, true).is_err());
        assert!(solve_lambert(r1, r2, -100.0, MU_EARTH, true).is_err());
        
        let r_zero = Vector3::new(0.0, 0.0, 0.0);
        assert!(solve_lambert(r_zero, r2, 3600.0, MU_EARTH, true).is_err());
    }
    
    #[test]
    fn test_delta_v_calculation() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(0.0, 7000.0, 0.0);
        let v_circ = (MU_EARTH / 7000.0).sqrt();
        let v1_initial = Vector3::new(0.0, v_circ, 0.0);
        let v2_target = Vector3::new(-v_circ, 0.0, 0.0);
        let dt = 0.25 * 2.0 * PI * (7000.0_f64.powi(3) / MU_EARTH).sqrt();
        
        let (dv1, dv2, dv_total) = lambert_delta_v(
            r1, v1_initial, r2, v2_target, dt, MU_EARTH, true
        ).unwrap();
        
        assert!(dv1 >= 0.0, "Initial ΔV should be non-negative");
        assert!(dv2 >= 0.0, "Final ΔV should be non-negative");
        assert!((dv_total - dv1 - dv2).abs() < 1e-10);
        // For quarter-orbit at same radius, ΔV should be very small
        assert!(dv_total < 1.0, "ΔV too large for same-orbit transfer: {}", dv_total);
    }
}
