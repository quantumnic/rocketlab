//! Low-thrust trajectory optimization
//!
//! Implements the Edelbaum approximation for low-thrust orbit transfers,
//! Q-law guidance for orbit element targeting, and spiral trajectory analysis.
//!
//! Reference: Edelbaum (1961) "Propulsion Requirements for Controllable Satellites",
//! Petropoulos (2003) "Simple Control Laws for Low-Thrust Orbit Transfers"

use std::f64::consts::PI;

use crate::constants;

/// Edelbaum result for a low-thrust orbit transfer
#[derive(Debug, Clone)]
pub struct EdelbaumResult {
    /// Total ΔV required (km/s)
    pub delta_v: f64,
    /// Transfer time (s)
    pub transfer_time: f64,
    /// Initial orbit velocity (km/s)
    pub v0: f64,
    /// Final orbit velocity (km/s)
    pub vf: f64,
    /// Plane change angle (rad)
    pub delta_i: f64,
    /// Number of revolutions (approximate)
    pub revolutions: f64,
}

/// Edelbaum approximation for combined orbit raising and plane change.
///
/// Computes the minimum ΔV for a low-thrust transfer between circular orbits
/// with an inclination change. The thrust is applied continuously and the
/// steering angle varies optimally.
///
/// # Arguments
/// * `r0` - Initial orbit radius (km)
/// * `rf` - Final orbit radius (km)
/// * `delta_i` - Inclination change (rad)
/// * `thrust` - Continuous thrust (N)
/// * `mass` - Spacecraft mass (kg)
///
/// # Reference
/// ΔV = √(v₀² + vf² - 2·v₀·vf·cos(π/2·Δi))
pub fn edelbaum_transfer(r0: f64, rf: f64, delta_i: f64, thrust: f64, mass: f64) -> EdelbaumResult {
    let mu = constants::MU_EARTH; // km³/s²
    let v0 = (mu / r0).sqrt(); // km/s
    let vf = (mu / rf).sqrt(); // km/s

    // Edelbaum ΔV equation (km/s)
    let delta_v = (v0 * v0 + vf * vf - 2.0 * v0 * vf * (PI / 2.0 * delta_i).cos()).sqrt();

    // Acceleration in km/s² (thrust in N, mass in kg → m/s² → km/s²)
    let accel = thrust / mass / 1000.0;

    // Transfer time (constant acceleration approximation)
    let transfer_time = delta_v / accel;

    // Approximate revolutions: average orbit period × time
    let r_avg = (r0 + rf) / 2.0;
    let period_avg = 2.0 * PI * (r_avg.powi(3) / mu).sqrt();
    let revolutions = transfer_time / period_avg;

    EdelbaumResult {
        delta_v,
        transfer_time,
        v0,
        vf,
        delta_i,
        revolutions,
    }
}

/// Compute propellant mass for a low-thrust maneuver using the rocket equation.
///
/// mp = m0 × (1 - exp(-ΔV / (g0 × Isp)))
///
/// # Arguments
/// * `delta_v` - Delta-V (km/s)
/// * `m0` - Initial mass (kg)
/// * `isp` - Specific impulse (s)
pub fn propellant_mass(delta_v: f64, m0: f64, isp: f64) -> f64 {
    let ve = isp * constants::G0 / 1000.0; // km/s
    m0 * (1.0 - (-delta_v / ve).exp())
}

/// Simple spiral trajectory: radius vs time for continuous tangential thrust.
///
/// Returns Vec of (time_s, radius_km, velocity_km_s).
///
/// # Arguments
/// * `r0` - Initial radius (km)
/// * `thrust` - Thrust (N)
/// * `mass` - Mass (kg)
/// * `duration` - Duration (s)
/// * `dt` - Time step (s)
pub fn spiral_trajectory(
    r0: f64,
    thrust: f64,
    mass: f64,
    duration: f64,
    dt: f64,
) -> Vec<(f64, f64, f64)> {
    let mu = constants::MU_EARTH; // km³/s²
    let mut trajectory = Vec::new();
    let mut r = r0;
    let mut t = 0.0;
    let accel = thrust / mass / 1000.0; // km/s²

    while t < duration && r > 0.0 {
        let v = (mu / r).sqrt(); // km/s
        trajectory.push((t, r, v));

        // For a circular spiral, tangential thrust increases orbital energy
        // E = -μ/(2r), so dE/dt = F·v/m → dr/dt = 2r²·a·v/μ
        let dr = 2.0 * r * r * accel * v / mu;
        r += dr * dt;
        t += dt;
    }

    trajectory
}

/// Q-law proximity quotient for orbit element targeting.
///
/// The Q-law provides a Lyapunov feedback control law for low-thrust
/// orbit transfers. Returns the optimal thrust direction angles.
///
/// # Arguments
/// * `a_current` - Current semi-major axis (m)
/// * `e_current` - Current eccentricity
/// * `i_current` - Current inclination (rad)
/// * `a_target` - Target semi-major axis (m)
/// * `e_target` - Target eccentricity
/// * `i_target` - Target inclination (rad)
///
/// Returns (alpha, beta) thrust angles in radians (in-plane, out-of-plane).
pub fn qlaw_steering(
    a_current: f64,
    e_current: f64,
    i_current: f64,
    a_target: f64,
    e_target: f64,
    i_target: f64,
) -> (f64, f64) {
    // Weight factors for each orbital element
    let w_a = 1.0;
    let w_e = 1.0;
    let w_i = 1.0;

    // Errors
    let da = (a_target - a_current) / a_target;
    let de = e_target - e_current;
    let di = i_target - i_current;

    // In-plane angle: combination of SMA and eccentricity steering
    let alpha = (w_a * da + w_e * de).atan2(1.0);

    // Out-of-plane angle: inclination steering
    let beta = (w_i * di).atan2((w_a * da * da + w_e * de * de).sqrt().max(1e-12));

    (alpha, beta)
}

/// Estimate transfer time for an electric propulsion spiral.
///
/// Uses the analytical approximation for constant-thrust tangential burns.
///
/// t = (m / F) × |√(μ/r0) - √(μ/rf)|
pub fn spiral_transfer_time(r0: f64, rf: f64, thrust: f64, mass: f64) -> f64 {
    let mu = constants::MU_EARTH; // km³/s²
    let v0 = (mu / r0).sqrt(); // km/s
    let vf = (mu / rf).sqrt(); // km/s
                               // Convert: (mass/thrust) is in s²/m, velocity diff in km/s → multiply by 1000
    (mass / thrust) * (v0 - vf).abs() * 1000.0
}

/// Ion engine parameters for common electric propulsion systems
#[derive(Debug, Clone)]
pub struct ElectricPropulsion {
    pub name: &'static str,
    pub thrust_n: f64,
    pub isp_s: f64,
    pub power_w: f64,
    pub efficiency: f64,
}

/// Database of electric propulsion systems
pub fn ep_database() -> Vec<ElectricPropulsion> {
    vec![
        ElectricPropulsion {
            name: "NSTAR (Dawn)",
            thrust_n: 0.092,
            isp_s: 3100.0,
            power_w: 2300.0,
            efficiency: 0.62,
        },
        ElectricPropulsion {
            name: "NEXT-C",
            thrust_n: 0.236,
            isp_s: 4190.0,
            power_w: 6900.0,
            efficiency: 0.70,
        },
        ElectricPropulsion {
            name: "Hall SPT-100",
            thrust_n: 0.083,
            isp_s: 1600.0,
            power_w: 1350.0,
            efficiency: 0.50,
        },
        ElectricPropulsion {
            name: "BHT-600 (Starlink)",
            thrust_n: 0.039,
            isp_s: 1500.0,
            power_w: 600.0,
            efficiency: 0.48,
        },
        ElectricPropulsion {
            name: "PPS-1350 (SMART-1)",
            thrust_n: 0.088,
            isp_s: 1650.0,
            power_w: 1500.0,
            efficiency: 0.55,
        },
        ElectricPropulsion {
            name: "X3 Nested Hall",
            thrust_n: 5.4,
            isp_s: 2340.0,
            power_w: 100_000.0,
            efficiency: 0.60,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edelbaum_leo_to_geo() {
        // LEO (400 km) to GEO with 28.5° plane change — all in km
        let r_leo = constants::R_EARTH + 400.0;
        let r_geo = constants::R_EARTH + 35_786.0;
        let delta_i = 28.5_f64.to_radians();

        let result = edelbaum_transfer(r_leo, r_geo, delta_i, 1.0, 1000.0);

        // Edelbaum ΔV for LEO→GEO with 28.5° should be ~5.8-6.2 km/s
        assert!(
            result.delta_v > 5.0 && result.delta_v < 7.0,
            "ΔV = {} km/s",
            result.delta_v
        );

        // Verify orbital velocities
        let v_leo = (constants::MU_EARTH / r_leo).sqrt();
        assert!(
            (result.v0 - v_leo).abs() < 0.001,
            "v0 = {} vs expected {}",
            result.v0,
            v_leo
        );
    }

    #[test]
    fn test_edelbaum_pure_raise() {
        // Pure orbit raising, no plane change: ΔV = |v0 - vf|
        let r0 = constants::R_EARTH + 400.0;
        let rf = constants::R_EARTH + 800.0;

        let result = edelbaum_transfer(r0, rf, 0.0, 1.0, 1000.0);

        let v0 = (constants::MU_EARTH / r0).sqrt();
        let vf = (constants::MU_EARTH / rf).sqrt();
        let expected_dv = (v0 - vf).abs();

        assert!(
            (result.delta_v - expected_dv).abs() < 0.001,
            "ΔV = {} vs expected {}",
            result.delta_v,
            expected_dv
        );
    }

    #[test]
    fn test_edelbaum_pure_plane_change() {
        // Pure plane change at same altitude
        let r = constants::R_EARTH + 400.0;
        let delta_i = 10.0_f64.to_radians();

        let result = edelbaum_transfer(r, r, delta_i, 1.0, 1000.0);

        let v = (constants::MU_EARTH / r).sqrt();
        let expected = 2.0 * v * (PI * delta_i / 4.0).sin();

        assert!(
            (result.delta_v - expected).abs() < 0.001,
            "ΔV = {} vs expected {}",
            result.delta_v,
            expected
        );
    }

    #[test]
    fn test_propellant_mass() {
        // 1000 kg spacecraft, 5 km/s ΔV, 3000s Isp
        let mp = propellant_mass(5.0, 1000.0, 3000.0);
        // ve = 3000 * 9.80665 / 1000 = 29.42 km/s
        // mp = 1000 * (1 - exp(-5/29.42)) ≈ 156 kg
        assert!(mp > 150.0 && mp < 165.0, "mp = {mp} kg");
    }

    #[test]
    fn test_spiral_trajectory() {
        let r0 = constants::R_EARTH + 400.0; // km
        let traj = spiral_trajectory(r0, 1.0, 1000.0, 86400.0, 100.0);
        assert!(!traj.is_empty());
        let r_final = traj.last().unwrap().1;
        assert!(r_final > r0, "Spiral should raise orbit");
    }

    #[test]
    fn test_spiral_transfer_time() {
        let r0 = constants::R_EARTH + 400.0;
        let rf = constants::R_EARTH + 800.0;
        let t = spiral_transfer_time(r0, rf, 1.0, 1000.0);
        assert!(t > 0.0 && t < 365.0 * 86400.0, "t = {} days", t / 86400.0);
    }

    #[test]
    fn test_qlaw_steering_converged() {
        // At target: steering should be near zero
        let (alpha, beta) = qlaw_steering(7_000_000.0, 0.001, 0.5, 7_000_000.0, 0.001, 0.5);
        assert!(alpha.abs() < 0.01);
        assert!(beta.abs() < 0.1);
    }

    #[test]
    fn test_ep_database() {
        let db = ep_database();
        assert!(db.len() >= 5);
        for ep in &db {
            assert!(ep.thrust_n > 0.0);
            assert!(ep.isp_s > 500.0);
            assert!(ep.efficiency > 0.0 && ep.efficiency < 1.0);
        }
    }
}
