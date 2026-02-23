//! Mission design and trajectory planning.
//!
//! This module provides tools for designing spacecraft missions including:
//! - Hohmann and bi-elliptic transfers
//! - Gravity assist maneuvers (patched conics)
//! - Pork chop plots for trajectory optimization
//! - Launch window analysis

use nalgebra::Vector3;
use crate::constants::{MU_EARTH, MU_SUN, AU_KM, TWO_PI};
use crate::lambert::solve_lambert;
use chrono::{DateTime, Utc, Duration};
use std::f64::consts::PI;

/// Result of a Hohmann transfer calculation
#[derive(Debug, Clone)]
pub struct HohmannTransfer {
    /// Transfer semi-major axis
    pub transfer_sma: f64,
    /// Transfer eccentricity
    pub transfer_ecc: f64,
    /// Delta-V at initial orbit
    pub delta_v1: f64,
    /// Delta-V at final orbit
    pub delta_v2: f64,
    /// Total delta-V required
    pub total_delta_v: f64,
    /// Transfer time (seconds)
    pub transfer_time: f64,
    /// Phase angle at launch
    pub phase_angle: f64,
}

/// Result of a bi-elliptic transfer calculation
#[derive(Debug, Clone)]
pub struct BiEllipticTransfer {
    /// First transfer semi-major axis
    pub transfer1_sma: f64,
    /// Second transfer semi-major axis
    pub transfer2_sma: f64,
    /// Aphelion radius of transfer orbits
    pub aphelion_radius: f64,
    /// Delta-V at initial orbit
    pub delta_v1: f64,
    /// Delta-V at aphelion
    pub delta_v2: f64,
    /// Delta-V at final orbit
    pub delta_v3: f64,
    /// Total delta-V required
    pub total_delta_v: f64,
    /// Total transfer time (seconds)
    pub total_time: f64,
}

/// Gravity assist result
#[derive(Debug, Clone)]
pub struct GravityAssist {
    /// Incoming velocity relative to planet
    pub v_inf_in: f64,
    /// Outgoing velocity relative to planet
    pub v_inf_out: f64,
    /// Turn angle achieved
    pub turn_angle: f64,
    /// Flyby altitude
    pub flyby_altitude: f64,
    /// Delta-V equivalent provided by gravity assist
    pub delta_v_equivalent: f64,
    /// Minimum safe flyby radius
    pub periapsis: f64,
}

/// Pork chop plot point
#[derive(Debug, Clone)]
pub struct PorkChopPoint {
    /// Launch date
    pub launch_date: DateTime<Utc>,
    /// Arrival date
    pub arrival_date: DateTime<Utc>,
    /// Flight time (days)
    pub flight_time_days: f64,
    /// Departure delta-V (km/s)
    pub c3_departure: f64,
    /// Arrival delta-V (km/s)
    pub c3_arrival: f64,
    /// Total characteristic energy
    pub total_c3: f64,
    /// Departure declination angle (rad)
    pub declination: f64,
}

/// Launch window analysis result
#[derive(Debug, Clone)]
pub struct LaunchWindow {
    /// Window opening date
    pub window_start: DateTime<Utc>,
    /// Window closing date
    pub window_end: DateTime<Utc>,
    /// Optimal launch date in window
    pub optimal_date: DateTime<Utc>,
    /// Minimum delta-V in window (km/s)
    pub min_delta_v: f64,
    /// Window duration (days)
    pub duration_days: f64,
    /// Launch azimuth range (degrees)
    pub azimuth_range: (f64, f64),
}

/// Calculate Hohmann transfer between two circular orbits
pub fn hohmann_transfer(r1: f64, r2: f64, mu: f64) -> HohmannTransfer {
    // Transfer orbit semi-major axis
    let a_transfer = (r1 + r2) / 2.0;
    
    // Transfer orbit eccentricity
    let e_transfer = (r2 - r1) / (r1 + r2);
    
    // Orbital velocities
    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();
    
    // Transfer orbit velocities
    let v_transfer_p = (mu * (2.0 / r1 - 1.0 / a_transfer)).sqrt(); // At perigee
    let v_transfer_a = (mu * (2.0 / r2 - 1.0 / a_transfer)).sqrt(); // At apogee
    
    // Delta-Vs
    let delta_v1 = v_transfer_p - v1;
    let delta_v2 = v2 - v_transfer_a;
    let total_delta_v = delta_v1.abs() + delta_v2.abs();
    
    // Transfer time (half period)
    let transfer_time = PI * (a_transfer.powi(3) / mu).sqrt();
    
    // Phase angle calculation
    let n1 = (mu / r1.powi(3)).sqrt(); // Mean motion of inner orbit
    let n2 = (mu / r2.powi(3)).sqrt(); // Mean motion of outer orbit
    let phase_angle = PI - (n2 - n1) * transfer_time;
    
    HohmannTransfer {
        transfer_sma: a_transfer,
        transfer_ecc: e_transfer,
        delta_v1,
        delta_v2,
        total_delta_v,
        transfer_time,
        phase_angle,
    }
}

/// Calculate bi-elliptic transfer between two circular orbits
pub fn bi_elliptic_transfer(r1: f64, r2: f64, r_aphelion: f64, mu: f64) -> BiEllipticTransfer {
    // First transfer orbit (r1 to r_aphelion)
    let a1 = (r1 + r_aphelion) / 2.0;
    let v1_initial = (mu / r1).sqrt();
    let v1_transfer = (mu * (2.0 / r1 - 1.0 / a1)).sqrt();
    let delta_v1 = v1_transfer - v1_initial;
    
    // Second transfer orbit (r_aphelion to r2)
    let a2 = (r_aphelion + r2) / 2.0;
    let v2_transfer_out = (mu * (2.0 / r_aphelion - 1.0 / a1)).sqrt();
    let v2_transfer_in = (mu * (2.0 / r_aphelion - 1.0 / a2)).sqrt();
    let delta_v2 = v2_transfer_in - v2_transfer_out;
    
    // Final orbit insertion
    let v2_final = (mu / r2).sqrt();
    let v2_arrival = (mu * (2.0 / r2 - 1.0 / a2)).sqrt();
    let delta_v3 = v2_final - v2_arrival;
    
    let total_delta_v = delta_v1.abs() + delta_v2.abs() + delta_v3.abs();
    
    // Transfer times
    let time1 = PI * (a1.powi(3) / mu).sqrt();
    let time2 = PI * (a2.powi(3) / mu).sqrt();
    let total_time = time1 + time2;
    
    BiEllipticTransfer {
        transfer1_sma: a1,
        transfer2_sma: a2,
        aphelion_radius: r_aphelion,
        delta_v1,
        delta_v2,
        delta_v3,
        total_delta_v,
        total_time,
    }
}

/// Calculate gravity assist parameters
pub fn gravity_assist(
    v_inf_in: f64,
    planet_mu: f64,
    flyby_altitude: f64,
    planet_radius: f64,
) -> GravityAssist {
    let periapsis = planet_radius + flyby_altitude;
    
    // Calculate turn angle using vis-viva and angular momentum conservation
    let e_hyperbola = 1.0 + (periapsis * v_inf_in.powi(2)) / planet_mu;
    let turn_angle = 2.0 * (1.0 / e_hyperbola).asin();
    
    // For elastic encounter, |v_inf_out| = |v_inf_in|
    let v_inf_out = v_inf_in;
    
    // Delta-V equivalent (what it would cost with propulsion)
    let delta_v_equivalent = 2.0 * v_inf_in * (turn_angle / 2.0).sin();
    
    GravityAssist {
        v_inf_in,
        v_inf_out,
        turn_angle,
        flyby_altitude,
        delta_v_equivalent,
        periapsis,
    }
}

/// Generate pork chop plot data for interplanetary transfer
pub fn generate_pork_chop_data(
    departure_planet_pos: fn(f64) -> Vector3<f64>,
    arrival_planet_pos: fn(f64) -> Vector3<f64>,
    start_date: DateTime<Utc>,
    end_date: DateTime<Utc>,
    flight_time_min_days: f64,
    flight_time_max_days: f64,
    step_days: f64,
) -> Vec<PorkChopPoint> {
    let mut points = Vec::new();
    let start_mjd = date_to_mjd(start_date);
    let end_mjd = date_to_mjd(end_date);
    
    let mut departure_mjd = start_mjd;
    while departure_mjd <= end_mjd {
        let mut flight_time = flight_time_min_days;
        
        while flight_time <= flight_time_max_days {
            let arrival_mjd = departure_mjd + flight_time;
            
            // Get planet positions
            let r1 = departure_planet_pos(departure_mjd);
            let r2 = arrival_planet_pos(arrival_mjd);
            
            // Solve Lambert problem
            let tof = flight_time * 86400.0; // Convert days to seconds
            if let Ok(solution) = solve_lambert(r1, r2, tof, MU_SUN, true) {
                // Calculate characteristic energies (C3)
                let v1: Vector3<f64> = solution.v1;
                let v2: Vector3<f64> = solution.v2;
                let c3_departure = v1.magnitude_squared() - 2.0 * MU_SUN / r1.magnitude();
                let c3_arrival = v2.magnitude_squared() - 2.0 * MU_SUN / r2.magnitude();
                
                // Calculate declination (simplified)
                let declination = (v1.z / v1.magnitude()).asin();
                
                points.push(PorkChopPoint {
                    launch_date: mjd_to_date(departure_mjd),
                    arrival_date: mjd_to_date(arrival_mjd),
                    flight_time_days: flight_time,
                    c3_departure: c3_departure / 1e6, // Convert to km²/s²
                    c3_arrival: c3_arrival / 1e6,
                    total_c3: (c3_departure + c3_arrival) / 1e6,
                    declination,
                });
            }
            
            flight_time += step_days;
        }
        
        departure_mjd += step_days;
    }
    
    points
}

/// Find launch windows with delta-V constraints
pub fn find_launch_windows(
    pork_chop_data: &[PorkChopPoint],
    max_delta_v: f64,
    min_window_days: f64,
) -> Vec<LaunchWindow> {
    let mut windows = Vec::new();
    let mut in_window = false;
    let mut window_start = None;
    let mut window_points = Vec::new();
    
    for point in pork_chop_data {
        let delta_v = (point.c3_departure + point.c3_arrival).sqrt();
        
        if delta_v <= max_delta_v {
            if !in_window {
                // Start new window
                in_window = true;
                window_start = Some(point.launch_date);
                window_points.clear();
            }
            window_points.push(point);
        } else if in_window {
            // End current window
            in_window = false;
            
            if let Some(start) = window_start {
                let duration = (point.launch_date - start).num_days() as f64;
                
                if duration >= min_window_days && !window_points.is_empty() {
                    // Find optimal point in window
                    let optimal_point = window_points
                        .iter()
                        .min_by(|a, b| {
                            (a.c3_departure + a.c3_arrival)
                                .partial_cmp(&(b.c3_departure + b.c3_arrival))
                                .unwrap()
                        })
                        .unwrap();
                    
                    let min_delta_v = (optimal_point.c3_departure + optimal_point.c3_arrival).sqrt();
                    
                    // Calculate azimuth range (simplified)
                    let azimuth_range = (0.0, 360.0); // Placeholder - real calculation would be more complex
                    
                    windows.push(LaunchWindow {
                        window_start: start,
                        window_end: window_points.last().unwrap().launch_date,
                        optimal_date: optimal_point.launch_date,
                        min_delta_v,
                        duration_days: duration,
                        azimuth_range,
                    });
                }
            }
        }
    }
    
    windows
}

/// Calculate synodic period between two planets
pub fn synodic_period(planet1_period: f64, planet2_period: f64) -> f64 {
    1.0 / ((1.0 / planet1_period) - (1.0 / planet2_period)).abs()
}

/// Calculate type-I vs type-II transfer comparison
pub fn compare_transfer_types(r1: f64, r2: f64, flight_time: f64, mu: f64) -> (f64, f64) {
    // For now, return placeholder values
    // Real implementation would solve Lambert problem for both short and long way
    let type1_dv = 5.0; // km/s
    let type2_dv = 7.0; // km/s
    (type1_dv, type2_dv)
}

/// Convert DateTime to Modified Julian Date
fn date_to_mjd(date: DateTime<Utc>) -> f64 {
    // Simplified conversion - real implementation would be more precise
    let epoch = DateTime::parse_from_rfc3339("2000-01-01T00:00:00Z").unwrap();
    let days_since_j2000 = (date.signed_duration_since(epoch)).num_days() as f64;
    51544.5 + days_since_j2000 // MJD of J2000.0 is 51544.5
}

/// Convert Modified Julian Date to DateTime
fn mjd_to_date(mjd: f64) -> DateTime<Utc> {
    let j2000_mjd = 51544.5;
    let days_since_j2000 = mjd - j2000_mjd;
    let epoch = DateTime::parse_from_rfc3339("2000-01-01T00:00:00Z").unwrap();
    epoch.with_timezone(&Utc) + Duration::days(days_since_j2000 as i64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hohmann_transfer_earth_mars() {
        // Earth to Mars transfer (approximate)
        let r_earth = 1.0 * AU_KM * 1000.0; // Convert to meters
        let r_mars = 1.52 * AU_KM * 1000.0;
        
        let transfer = hohmann_transfer(r_earth, r_mars, MU_SUN);
        
        // Check transfer semi-major axis
        assert_relative_eq!(transfer.transfer_sma, 1.26 * AU_KM * 1000.0, epsilon = 0.01);
        
        // Transfer time should be about 259 days
        let transfer_days = transfer.transfer_time / 86400.0;
        assert_relative_eq!(transfer_days, 259.0, epsilon = 10.0);
        
        // Total delta-V should be reasonable (3-6 km/s)
        assert!(transfer.total_delta_v > 3000.0);
        assert!(transfer.total_delta_v < 6000.0);
    }

    #[test]
    fn test_bi_elliptic_transfer() {
        let r1 = 6678.0e3; // LEO
        let r2 = 42164.0e3; // GEO
        let r_aphelion = 80000.0e3; // High aphelion
        
        let transfer = bi_elliptic_transfer(r1, r2, r_aphelion, MU_EARTH);
        
        // Bi-elliptic should have 3 delta-V components
        assert!(transfer.delta_v1 > 0.0);
        assert!(transfer.delta_v2 != 0.0);
        assert!(transfer.delta_v3 < 0.0);
        
        // Total transfer time should be longer than Hohmann
        let hohmann = hohmann_transfer(r1, r2, MU_EARTH);
        assert!(transfer.total_time > hohmann.transfer_time);
    }

    #[test]
    fn test_gravity_assist() {
        // Earth gravity assist parameters
        let v_inf = 5000.0; // 5 km/s hyperbolic excess velocity
        let flyby_alt = 200.0e3; // 200 km altitude
        let earth_radius = 6371.0e3;
        
        let assist = gravity_assist(v_inf, MU_EARTH, flyby_alt, earth_radius);
        
        // Turn angle should be reasonable
        assert!(assist.turn_angle > 0.0);
        assert!(assist.turn_angle < PI);
        
        // Delta-V equivalent should be positive
        assert!(assist.delta_v_equivalent > 0.0);
        
        // V-infinity should be conserved in magnitude
        assert_relative_eq!(assist.v_inf_out, v_inf, epsilon = 1e-6);
    }

    #[test]
    fn test_synodic_period() {
        // Earth-Mars synodic period
        let earth_period = 365.25; // days
        let mars_period = 686.98; // days
        
        let synodic = synodic_period(earth_period, mars_period);
        
        // Should be about 779 days (26 months)
        assert_relative_eq!(synodic, 779.0, epsilon = 10.0);
    }

    #[test]
    fn test_date_conversions() {
        let test_date = DateTime::parse_from_rfc3339("2020-01-01T00:00:00Z")
            .unwrap()
            .with_timezone(&Utc);
        
        let mjd = date_to_mjd(test_date);
        let converted_back = mjd_to_date(mjd);
        
        // Should round-trip approximately
        let diff_seconds = (test_date - converted_back).num_seconds().abs();
        assert!(diff_seconds < 86400); // Within a day
    }
}