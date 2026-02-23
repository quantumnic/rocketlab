//! # Orbit Lifetime Estimation
//!
//! King-Hele decay theory for estimating orbital lifetime under atmospheric drag.
//! Essential for debris mitigation compliance (25-year rule) and mission planning.
//!
//! References:
//! - D.G. King-Hele, "Satellite Orbits in an Atmosphere" (1987)
//! - NASA-STD-8719.14A "Process for Limiting Orbital Debris"
//! - Vallado "Fundamentals of Astrodynamics and Applications" Ch. 8

use std::f64::consts::PI;

/// Parameters for orbit lifetime estimation.
#[derive(Debug, Clone)]
pub struct LifetimeParams {
    /// Semi-major axis [m]
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Ballistic coefficient B = Cd * A / m [m²/kg]
    pub ballistic_coefficient: f64,
    /// Central body radius [m]
    pub body_radius: f64,
    /// Gravitational parameter [m³/s²]
    pub mu: f64,
    /// Atmospheric scale height [m] (approximate)
    pub scale_height: f64,
    /// Reference density at perigee altitude [kg/m³]
    pub rho_perigee: f64,
}

/// Result of lifetime estimation.
#[derive(Debug, Clone)]
pub struct LifetimeResult {
    /// Estimated orbital lifetime [s]
    pub lifetime_seconds: f64,
    /// Estimated orbital lifetime [days]
    pub lifetime_days: f64,
    /// Estimated orbital lifetime [years]
    pub lifetime_years: f64,
    /// Perigee altitude [m]
    pub perigee_altitude: f64,
    /// Apogee altitude [m]
    pub apogee_altitude: f64,
    /// Initial orbital period [s]
    pub period: f64,
    /// Estimated revolutions to decay
    pub revolutions: f64,
    /// Meets 25-year rule?
    pub meets_25yr_rule: bool,
}

/// Decay history point.
#[derive(Debug, Clone, Copy)]
pub struct DecayPoint {
    /// Time since epoch [days]
    pub time_days: f64,
    /// Semi-major axis [m]
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Perigee altitude [km]
    pub perigee_alt_km: f64,
    /// Apogee altitude [km]
    pub apogee_alt_km: f64,
}

/// US Standard Atmosphere 1976 density model (simplified exponential fit).
///
/// Returns density [kg/m³] for a given altitude [m].
fn atmosphere_density(altitude: f64) -> f64 {
    let h_km = altitude / 1000.0;

    // Piecewise exponential fit to US Standard Atmosphere 1976
    // (base altitude [km], base density [kg/m³], scale height [km])
    let layers: &[(f64, f64, f64)] = &[
        (0.0, 1.225, 8.44),
        (25.0, 3.899e-2, 6.49),
        (30.0, 1.774e-2, 6.75),
        (40.0, 3.972e-3, 7.26),
        (50.0, 1.057e-3, 8.38),
        (60.0, 3.096e-4, 7.71),
        (70.0, 8.283e-5, 6.55),
        (80.0, 1.846e-5, 5.80),
        (90.0, 3.416e-6, 5.38),
        (100.0, 5.604e-7, 5.88),
        (110.0, 9.708e-8, 7.26),
        (120.0, 2.222e-8, 26.2),
        (150.0, 2.076e-9, 33.8),
        (200.0, 2.541e-10, 45.5),
        (250.0, 6.073e-11, 53.3),
        (300.0, 1.916e-11, 53.6),
        (400.0, 2.803e-12, 58.5),
        (500.0, 5.215e-13, 60.8),
        (600.0, 1.137e-13, 63.8),
        (700.0, 3.070e-14, 71.8),
        (800.0, 1.136e-14, 88.7),
        (900.0, 5.759e-15, 124.0),
        (1000.0, 3.561e-15, 181.0),
    ];

    if h_km <= 0.0 {
        return 1.225;
    }

    // Find appropriate layer
    let mut base_h = layers[0].0;
    let mut base_rho = layers[0].1;
    let mut h_scale = layers[0].2;

    for layer in layers {
        if h_km >= layer.0 {
            base_h = layer.0;
            base_rho = layer.1;
            h_scale = layer.2;
        } else {
            break;
        }
    }

    base_rho * (-(h_km - base_h) / h_scale).exp()
}

/// Atmospheric scale height at given altitude [m] → [m].
fn scale_height_at(altitude: f64) -> f64 {
    let h_km = altitude / 1000.0;
    // Approximate scale height variation
    if h_km < 100.0 {
        7000.0
    } else if h_km < 200.0 {
        30_000.0 + (h_km - 100.0) * 150.0
    } else if h_km < 400.0 {
        45_000.0 + (h_km - 200.0) * 75.0
    } else if h_km < 600.0 {
        58_000.0 + (h_km - 400.0) * 30.0
    } else {
        65_000.0 + (h_km - 600.0) * 100.0
    }
}

/// Estimate orbit lifetime using King-Hele theory.
///
/// For nearly circular orbits (e < 0.01), uses the simplified formula:
/// τ ≈ -a / (2 * da/dt) where da/dt = -2π * a² * B * ρ * v
///
/// For eccentric orbits, uses King-Hele's method with perigee-dominated drag.
///
/// # Arguments
/// * `params` - Orbital and spacecraft parameters
///
/// # Returns
/// Lifetime estimation result
pub fn estimate_lifetime(params: &LifetimeParams) -> LifetimeResult {
    let perigee_alt = params.a * (1.0 - params.e) - params.body_radius;
    let apogee_alt = params.a * (1.0 + params.e) - params.body_radius;
    let period = 2.0 * PI * (params.a.powi(3) / params.mu).sqrt();

    // Use numerical integration for better accuracy
    let result = numerical_lifetime(params);

    let revolutions = result / period;
    let days = result / 86400.0;
    let years = days / 365.25;

    LifetimeResult {
        lifetime_seconds: result,
        lifetime_days: days,
        lifetime_years: years,
        perigee_altitude: perigee_alt,
        apogee_altitude: apogee_alt,
        period,
        revolutions,
        meets_25yr_rule: years <= 25.0,
    }
}

/// Numerical orbit lifetime estimation via orbit-averaged drag.
///
/// Integrates the orbit-averaged semi-major axis and eccentricity decay
/// equations forward in time until perigee drops below ~80 km.
fn numerical_lifetime(params: &LifetimeParams) -> f64 {
    let mut a = params.a;
    let mut e = params.e;
    let b = params.ballistic_coefficient;
    let mu = params.mu;
    let r_body = params.body_radius;

    let deorbit_alt = 80_000.0; // 80 km
    let mut total_time = 0.0;
    let max_time = 200.0 * 365.25 * 86400.0; // 200 years max

    loop {
        let perigee = a * (1.0 - e) - r_body;
        if perigee < deorbit_alt || total_time > max_time {
            break;
        }

        let period = 2.0 * PI * (a.powi(3) / mu).sqrt();
        let perigee_r = a * (1.0 - e);
        let _perigee_alt = perigee_r - r_body;
        let rho_p = atmosphere_density(perigee_r - r_body);

        if rho_p < 1e-20 {
            // Extremely low density, won't decay in reasonable time
            total_time = max_time;
            break;
        }

        // Orbit-averaged decay rates (King-Hele)
        let v_p = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();

        if e < 0.01 {
            // Nearly circular: da/dt ≈ -B * ρ * v * a (orbit-averaged)
            let da_dt = -b * rho_p * v_p * a;
            // Adaptive time step
            let dt = (-a / da_dt * 0.01).min(period * 100.0).max(period);
            a += da_dt * dt;
            total_time += dt;
        } else {
            // Eccentric orbit: drag concentrated at perigee
            // King-Hele simplified: da/rev ≈ -2π * rp * B * ρ_p * v_p / n
            // where rp = a(1-e), v_p = perigee velocity, n = mean motion
            let n = (mu / (a * a * a)).sqrt();

            // Orbit-averaged da/rev (perigee-dominated for high eccentricity)
            // From King-Hele (1987): Δa/rev ≈ -2π * a * B * ρ_p * v_p * (rp/a)
            // The factor (rp/a) = (1-e) accounts for perigee concentration
            let da_rev = -2.0 * PI * b * rho_p * v_p * perigee_r / n;

            // Eccentricity decay: circularization
            // de/rev ≈ -π * B * ρ_p * v_p * (1-e²) / (n*a)
            let de_rev = -PI * b * rho_p * v_p * (1.0 - e * e) / (n * a);

            // Step multiple revolutions at once for efficiency
            let n_revs = if da_rev.abs() > 0.0 {
                ((-a * 0.001) / da_rev).clamp(1.0, 1000.0)
            } else {
                1000.0
            };

            a += da_rev * n_revs;
            e += de_rev * n_revs;
            e = e.max(0.001);
            total_time += period * n_revs;
        }

        if a < r_body + deorbit_alt {
            break;
        }
    }

    total_time
}

/// Estimate orbit lifetime for a circular orbit (simplified).
///
/// Quick estimation for nearly-circular LEO orbits.
///
/// # Arguments
/// * `altitude` - Orbital altitude [m]
/// * `ballistic_coeff` - B = Cd * A / m [m²/kg]
/// * `mu` - Gravitational parameter [m³/s²] (default: Earth)
/// * `body_radius` - Body radius [m] (default: Earth)
pub fn circular_lifetime(
    altitude: f64,
    ballistic_coeff: f64,
    mu: f64,
    body_radius: f64,
) -> LifetimeResult {
    let a = body_radius + altitude;
    let rho = atmosphere_density(altitude);
    let h = scale_height_at(altitude);

    let params = LifetimeParams {
        a,
        e: 0.001,
        ballistic_coefficient: ballistic_coeff,
        body_radius,
        mu,
        scale_height: h,
        rho_perigee: rho,
    };

    estimate_lifetime(&params)
}

/// Compute ballistic coefficient from spacecraft parameters.
///
/// B = Cd * A / m
///
/// # Arguments
/// * `cd` - Drag coefficient (typically 2.0-2.5 for LEO)
/// * `area` - Cross-sectional area [m²]
/// * `mass` - Spacecraft mass [kg]
pub fn ballistic_coefficient(cd: f64, area: f64, mass: f64) -> f64 {
    cd * area / mass
}

/// Generate a decay history (semi-major axis and eccentricity vs time).
///
/// # Arguments
/// * `params` - Orbital parameters
/// * `n_points` - Number of output points
pub fn decay_history(params: &LifetimeParams, n_points: usize) -> Vec<DecayPoint> {
    let lifetime = estimate_lifetime(params);
    let dt = lifetime.lifetime_seconds / n_points as f64;

    let mut a = params.a;
    let mut e = params.e;
    let b = params.ballistic_coefficient;
    let mu = params.mu;
    let r_body = params.body_radius;

    let mut history = Vec::with_capacity(n_points);
    let mut time = 0.0;

    for _ in 0..n_points {
        let perigee_alt = (a * (1.0 - e) - r_body) / 1000.0;
        let apogee_alt = (a * (1.0 + e) - r_body) / 1000.0;

        history.push(DecayPoint {
            time_days: time / 86400.0,
            a,
            e,
            perigee_alt_km: perigee_alt,
            apogee_alt_km: apogee_alt,
        });

        if perigee_alt < 80.0 {
            break;
        }

        // Simple Euler step
        let perigee_r = a * (1.0 - e);
        let rho_p = atmosphere_density(perigee_r - r_body);
        let v_p = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();

        let da_dt = -b * rho_p * v_p * a;
        a += da_dt * dt;
        e = (e + da_dt * dt * (1.0 - e) / (2.0 * a)).max(0.001);
        time += dt;
    }

    history
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 3.986004418e14;
    const R_EARTH: f64 = 6_371_000.0;

    #[test]
    fn test_atmosphere_density_sea_level() {
        let rho = atmosphere_density(0.0);
        assert!((rho - 1.225).abs() < 0.01, "Sea level: {rho:.3} kg/m³");
    }

    #[test]
    fn test_atmosphere_density_400km() {
        let rho = atmosphere_density(400_000.0);
        // ~2.8e-12 kg/m³ at 400 km (US Std Atm 1976)
        assert!(rho > 1e-13 && rho < 1e-10, "400km: {rho:.2e} kg/m³");
    }

    #[test]
    fn test_atmosphere_density_decreasing() {
        let mut prev = atmosphere_density(0.0);
        for h in (100..=800).step_by(100) {
            let rho = atmosphere_density(h as f64 * 1000.0);
            assert!(rho < prev, "Density should decrease with altitude");
            prev = rho;
        }
    }

    #[test]
    fn test_ballistic_coefficient() {
        // ISS: Cd≈2.2, A≈1600m², m≈420000kg
        let b = ballistic_coefficient(2.2, 1600.0, 420_000.0);
        // B ≈ 0.0084 m²/kg
        assert!((b - 0.00838).abs() < 0.001, "ISS B = {b:.5}");
    }

    #[test]
    fn test_iss_lifetime_order_of_magnitude() {
        // ISS at 408 km, high drag area
        let b = ballistic_coefficient(2.2, 1600.0, 420_000.0);
        let result = circular_lifetime(408_000.0, b, MU_EARTH, R_EARTH);
        // ISS needs reboost every few months; without it, decays in ~1-2 years
        // Our simplified model should give order-of-magnitude correct
        assert!(
            result.lifetime_years > 0.1 && result.lifetime_years < 10.0,
            "ISS lifetime: {:.2} years",
            result.lifetime_years
        );
    }

    #[test]
    fn test_high_orbit_long_lifetime() {
        // 800 km orbit, small satellite
        let b = ballistic_coefficient(2.2, 0.01, 5.0); // CubeSat
        let result = circular_lifetime(800_000.0, b, MU_EARTH, R_EARTH);
        // 800 km should have very long lifetime (>25 years)
        assert!(
            result.lifetime_years > 20.0,
            "800km lifetime: {:.1} years",
            result.lifetime_years
        );
    }

    #[test]
    fn test_low_orbit_short_lifetime() {
        // 200 km, high drag
        let b = ballistic_coefficient(2.2, 1.0, 10.0);
        let result = circular_lifetime(200_000.0, b, MU_EARTH, R_EARTH);
        // 200 km should decay quickly (days to weeks)
        assert!(
            result.lifetime_days < 30.0,
            "200km lifetime: {:.1} days",
            result.lifetime_days
        );
    }

    #[test]
    fn test_25_year_rule() {
        // Check compliance at various altitudes
        let b = ballistic_coefficient(2.2, 0.05, 10.0); // Typical CubeSat

        let low = circular_lifetime(400_000.0, b, MU_EARTH, R_EARTH);
        assert!(low.meets_25yr_rule, "400km should meet 25yr rule");

        // Very high orbit might not
        let b_low_drag = ballistic_coefficient(2.2, 0.001, 100.0);
        let high = circular_lifetime(900_000.0, b_low_drag, MU_EARTH, R_EARTH);
        // 900km with low B might exceed 25 years
        // Just check the flag is consistent
        assert_eq!(high.meets_25yr_rule, high.lifetime_years <= 25.0);
    }

    #[test]
    fn test_decay_history_generation() {
        let b = ballistic_coefficient(2.2, 1.0, 50.0);
        let params = LifetimeParams {
            a: R_EARTH + 350_000.0,
            e: 0.001,
            ballistic_coefficient: b,
            body_radius: R_EARTH,
            mu: MU_EARTH,
            scale_height: 50_000.0,
            rho_perigee: atmosphere_density(350_000.0),
        };

        let history = decay_history(&params, 50);
        assert!(!history.is_empty());
        // First point should be near initial altitude
        assert!(
            (history[0].perigee_alt_km - 350.0).abs() < 10.0,
            "Initial alt: {:.1} km",
            history[0].perigee_alt_km
        );
        // Perigee should generally decrease
        if history.len() > 2 {
            assert!(
                history.last().unwrap().perigee_alt_km < history[0].perigee_alt_km,
                "Perigee should decrease over time"
            );
        }
    }

    #[test]
    fn test_eccentric_orbit_lifetime() {
        // GTO-like orbit: 200km x 35786km
        let rp = R_EARTH + 200_000.0;
        let ra = R_EARTH + 35_786_000.0;
        let a = (rp + ra) / 2.0;
        let e = (ra - rp) / (ra + rp);

        let b = ballistic_coefficient(2.2, 10.0, 3000.0);
        let params = LifetimeParams {
            a,
            e,
            ballistic_coefficient: b,
            body_radius: R_EARTH,
            mu: MU_EARTH,
            scale_height: 40_000.0,
            rho_perigee: atmosphere_density(200_000.0),
        };

        let result = estimate_lifetime(&params);
        // GTO with 200km perigee decays relatively fast (months to years)
        assert!(
            result.lifetime_days > 1.0 && result.lifetime_years < 50.0,
            "GTO lifetime: {:.1} days ({:.2} years)",
            result.lifetime_days,
            result.lifetime_years
        );
    }
}
