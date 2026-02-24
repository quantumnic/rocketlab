//! Sphere of Influence calculations.
//!
//! Two classical models for determining the radius within which a secondary
//! body's gravity dominates over the primary:
//!
//! - **Laplace SOI**: r_SOI = a * (m_body / m_primary)^(2/5)
//! - **Hill sphere**: r_Hill = a * (m_body / (3 * m_primary))^(1/3)
//!
//! References:
//! - Battin, R.H. "An Introduction to the Mathematics and Methods of
//!   Astrodynamics", AIAA, 1999, §8.3
//! - Bate, Mueller & White, "Fundamentals of Astrodynamics", §8.4
//! - Roy, A.E. "Orbital Motion", 4th ed., §7.6

use crate::constants::*;

/// Laplace sphere of influence radius (km).
///
/// r_SOI = a * (mu_body / mu_primary)^(2/5)
///
/// where `a` is the semi-major axis of the body's orbit around the primary,
/// `mu_body` is the gravitational parameter of the orbiting body, and
/// `mu_primary` is the gravitational parameter of the central body.
pub fn laplace_soi(a: f64, mu_body: f64, mu_primary: f64) -> f64 {
    a * (mu_body / mu_primary).powf(2.0 / 5.0)
}

/// Hill sphere radius (km).
///
/// r_Hill = a * (mu_body / (3 * mu_primary))^(1/3)
pub fn hill_radius(a: f64, mu_body: f64, mu_primary: f64) -> f64 {
    a * (mu_body / (3.0 * mu_primary)).powf(1.0 / 3.0)
}

/// Pre-computed SOI data for solar system bodies relative to the Sun.
#[derive(Debug, Clone, Copy)]
pub struct BodySOI {
    pub name: &'static str,
    pub laplace_soi_km: f64,
    pub hill_radius_km: f64,
}

/// Get sphere of influence data for common solar system bodies.
pub fn solar_system_soi() -> Vec<BodySOI> {
    let bodies: Vec<(&str, f64, f64)> = vec![
        ("Earth", A_EARTH, MU_EARTH),
        ("Mars", A_MARS, MU_MARS),
        ("Venus", A_VENUS, MU_VENUS),
        ("Jupiter", A_JUPITER, MU_JUPITER),
        ("Moon (around Earth)", 384_400.0, MU_MOON),
    ];

    bodies
        .iter()
        .map(|(name, a, mu_body)| {
            // Moon uses Earth as primary; all others use Sun
            let mu_primary = if *name == "Moon (around Earth)" {
                MU_EARTH
            } else {
                MU_SUN
            };
            BodySOI {
                name,
                laplace_soi_km: laplace_soi(*a, *mu_body, mu_primary),
                hill_radius_km: hill_radius(*a, *mu_body, mu_primary),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_earth_soi_laplace() {
        // Earth's Laplace SOI ~ 924,000 km (Vallado: ~929,000 km)
        let soi = laplace_soi(A_EARTH, MU_EARTH, MU_SUN);
        assert!(
            (soi - 924_000.0).abs() < 30_000.0,
            "Earth Laplace SOI = {soi} km, expected ~924,000 km"
        );
    }

    #[test]
    fn test_earth_hill_radius() {
        // Earth's Hill sphere ~ 1,500,000 km
        let rh = hill_radius(A_EARTH, MU_EARTH, MU_SUN);
        assert!(
            (rh - 1_500_000.0).abs() < 100_000.0,
            "Earth Hill radius = {rh} km, expected ~1,500,000 km"
        );
    }

    #[test]
    fn test_moon_soi() {
        // Moon's SOI around Earth ~ 66,100 km (Vallado)
        let soi = laplace_soi(384_400.0, MU_MOON, MU_EARTH);
        assert!(
            (soi - 66_100.0).abs() < 3_000.0,
            "Moon Laplace SOI = {soi} km, expected ~66,100 km"
        );
    }

    #[test]
    fn test_jupiter_soi() {
        // Jupiter's SOI ~ 48,200,000 km
        let soi = laplace_soi(A_JUPITER, MU_JUPITER, MU_SUN);
        assert!(
            (soi - 48_200_000.0).abs() < 2_000_000.0,
            "Jupiter Laplace SOI = {soi} km, expected ~48,200,000 km"
        );
    }

    #[test]
    fn test_mars_soi() {
        // Mars SOI ~ 577,000 km
        let soi = laplace_soi(A_MARS, MU_MARS, MU_SUN);
        assert!(
            (soi - 577_000.0).abs() < 30_000.0,
            "Mars Laplace SOI = {soi} km, expected ~577,000 km"
        );
    }

    #[test]
    fn test_solar_system_soi_count() {
        let bodies = solar_system_soi();
        assert_eq!(bodies.len(), 5);
    }

    #[test]
    fn test_hill_always_larger_than_laplace() {
        // Hill radius > Laplace SOI for any body
        let soi = laplace_soi(A_EARTH, MU_EARTH, MU_SUN);
        let hill = hill_radius(A_EARTH, MU_EARTH, MU_SUN);
        assert!(hill > soi, "Hill ({hill}) should be > Laplace SOI ({soi})");
    }
}
