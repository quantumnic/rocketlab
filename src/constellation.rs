//! Constellation design: Walker delta and star patterns.
//!
//! The Walker notation T/P/F describes a constellation where:
//! - T = total number of satellites
//! - P = number of equally spaced orbital planes
//! - F = phasing parameter (0 ≤ F < P)
//!
//! References:
//! - Walker, J.G., "Satellite Constellations", 1984
//! - Wertz, "Space Mission Engineering: The New SMAD", Ch. 7
//! - Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed.

use crate::constants::*;

/// A satellite's position within a constellation.
#[derive(Debug, Clone)]
pub struct ConstellationSatellite {
    /// Plane index (0-based)
    pub plane: usize,
    /// Satellite index within the plane (0-based)
    pub slot: usize,
    /// RAAN of the orbital plane (radians)
    pub raan: f64,
    /// Mean anomaly (radians)
    pub mean_anomaly: f64,
    /// Semi-major axis (km)
    pub semi_major_axis: f64,
    /// Inclination (radians)
    pub inclination: f64,
}

/// Walker constellation configuration.
#[derive(Debug, Clone)]
pub struct WalkerConstellation {
    /// Total number of satellites
    pub total: usize,
    /// Number of orbital planes
    pub planes: usize,
    /// Phasing parameter (0 ≤ F < P)
    pub phasing: usize,
    /// Semi-major axis (km)
    pub semi_major_axis: f64,
    /// Inclination (radians)
    pub inclination: f64,
    /// Whether this is a delta (0-180° RAAN spread) or star (0-360°) pattern
    pub is_star: bool,
}

impl WalkerConstellation {
    /// Create a Walker Delta constellation (RAAN spread over 180°).
    ///
    /// # Arguments
    /// * `t` — total satellites
    /// * `p` — number of planes
    /// * `f` — phasing parameter
    /// * `sma` — semi-major axis (km)
    /// * `inc_deg` — inclination (degrees)
    pub fn delta(t: usize, p: usize, f: usize, sma: f64, inc_deg: f64) -> Self {
        assert!(
            t.is_multiple_of(p),
            "Total satellites must be divisible by number of planes"
        );
        assert!(
            f < p,
            "Phasing parameter must be less than number of planes"
        );
        Self {
            total: t,
            planes: p,
            phasing: f,
            semi_major_axis: sma,
            inclination: inc_deg.to_radians(),
            is_star: false,
        }
    }

    /// Create a Walker Star constellation (RAAN spread over 360°).
    ///
    /// # Arguments
    /// * `t` — total satellites
    /// * `p` — number of planes
    /// * `f` — phasing parameter
    /// * `sma` — semi-major axis (km)
    /// * `inc_deg` — inclination (degrees)
    pub fn star(t: usize, p: usize, f: usize, sma: f64, inc_deg: f64) -> Self {
        assert!(
            t.is_multiple_of(p),
            "Total satellites must be divisible by number of planes"
        );
        assert!(
            f < p,
            "Phasing parameter must be less than number of planes"
        );
        Self {
            total: t,
            planes: p,
            phasing: f,
            semi_major_axis: sma,
            inclination: inc_deg.to_radians(),
            is_star: true,
        }
    }

    /// Number of satellites per plane.
    pub fn sats_per_plane(&self) -> usize {
        self.total / self.planes
    }

    /// Generate all satellite positions in the constellation.
    pub fn generate(&self) -> Vec<ConstellationSatellite> {
        let s = self.sats_per_plane();
        let raan_spread = if self.is_star { 2.0 * PI } else { PI };
        let delta_raan = raan_spread / self.planes as f64;
        let delta_ma = TWO_PI / s as f64;
        let phase_offset = TWO_PI * self.phasing as f64 / self.total as f64;

        let mut sats = Vec::with_capacity(self.total);

        for plane in 0..self.planes {
            let raan = plane as f64 * delta_raan;
            for slot in 0..s {
                let ma = slot as f64 * delta_ma + plane as f64 * phase_offset;
                let ma_normalized = ma % TWO_PI;

                sats.push(ConstellationSatellite {
                    plane,
                    slot,
                    raan,
                    mean_anomaly: ma_normalized,
                    semi_major_axis: self.semi_major_axis,
                    inclination: self.inclination,
                });
            }
        }

        sats
    }

    /// Orbital period of one satellite (seconds).
    pub fn period(&self) -> f64 {
        TWO_PI * (self.semi_major_axis.powi(3) / MU_EARTH).sqrt()
    }

    /// Altitude above Earth's equatorial surface (km).
    pub fn altitude(&self) -> f64 {
        self.semi_major_axis - R_EARTH_EQUATORIAL
    }

    /// Ground track repeat period approximation (seconds).
    ///
    /// For the constellation to provide repeating coverage.
    pub fn coverage_repeat_period(&self) -> f64 {
        self.period() * self.sats_per_plane() as f64
    }
}

/// Well-known constellation configurations.
pub mod presets {
    use super::*;

    /// GPS constellation: 24/6/1 Walker Star at 20,180 km alt, 55° inc
    pub fn gps() -> WalkerConstellation {
        WalkerConstellation::star(24, 6, 1, R_EARTH_EQUATORIAL + 20_180.0, 55.0)
    }

    /// Galileo constellation: 24/3/1 Walker Delta at 23,222 km alt, 56° inc
    pub fn galileo() -> WalkerConstellation {
        WalkerConstellation::delta(24, 3, 1, R_EARTH_EQUATORIAL + 23_222.0, 56.0)
    }

    /// Iridium constellation: 66/6/2 Walker Star at 780 km alt, 86.4° inc
    pub fn iridium() -> WalkerConstellation {
        WalkerConstellation::star(66, 6, 2, R_EARTH_EQUATORIAL + 780.0, 86.4)
    }

    /// Starlink shell 1: 72 planes × 22 sats = 1584 at 550 km, 53° inc
    pub fn starlink_shell1() -> WalkerConstellation {
        WalkerConstellation::star(1584, 72, 1, R_EARTH_EQUATORIAL + 550.0, 53.0)
    }

    /// OneWeb: 36/12/1 at 1200 km, 87.9° inc (Phase 1 partial)
    pub fn oneweb() -> WalkerConstellation {
        WalkerConstellation::delta(36, 12, 1, R_EARTH_EQUATORIAL + 1200.0, 87.9)
    }
}

/// Minimum elevation angle coverage analysis.
///
/// Estimates the maximum gap in coverage at a given latitude for
/// a Walker constellation, based on simple geometric analysis.
///
/// # Arguments
/// * `constellation` — the Walker constellation
/// * `min_elevation_deg` — minimum elevation angle for a satellite to be "visible" (degrees)
/// * `latitude_deg` — ground observer latitude (degrees)
///
/// # Returns
/// Approximate fraction of time with at least one satellite visible (0.0 to 1.0)
pub fn coverage_fraction(
    constellation: &WalkerConstellation,
    min_elevation_deg: f64,
    _latitude_deg: f64,
) -> f64 {
    let h = constellation.altitude();
    let el = min_elevation_deg.to_radians();

    // Earth central angle subtended by one satellite's coverage circle
    // ρ = arcsin(R_E * cos(el) / (R_E + h)) (nadir angle)
    // λ = 90° - el - ρ (Earth central angle)
    let sin_rho = R_EARTH_EQUATORIAL * el.cos() / (R_EARTH_EQUATORIAL + h);
    if sin_rho > 1.0 {
        return 0.0;
    }
    let rho = sin_rho.asin();
    let lambda = (PI / 2.0 - el - rho).max(0.0);

    // Coverage area of one satellite as fraction of hemisphere
    let single_coverage = (1.0 - lambda.cos()) / 2.0;

    // Total coverage (approximate, ignoring overlap for low density)
    let total = single_coverage * constellation.total as f64;

    // Clamp to 1.0 for dense constellations
    total.min(1.0)
}

/// Compute the minimum number of satellites needed for continuous single
/// coverage of a latitude band at a given altitude and elevation angle.
///
/// Uses the streets-of-coverage approach.
///
/// # Arguments
/// * `altitude_km` — orbital altitude (km)
/// * `min_elevation_deg` — minimum elevation angle (degrees)
/// * `inclination_deg` — orbital inclination (degrees)
///
/// # Returns
/// Approximate minimum number of satellites for single global coverage
pub fn min_sats_for_coverage(
    altitude_km: f64,
    min_elevation_deg: f64,
    _inclination_deg: f64,
) -> usize {
    let h = altitude_km;
    let el = min_elevation_deg.to_radians();

    let sin_rho = R_EARTH_EQUATORIAL * el.cos() / (R_EARTH_EQUATORIAL + h);
    if sin_rho > 1.0 {
        return usize::MAX; // Altitude too low for this elevation angle
    }
    let rho = sin_rho.asin();
    let lambda = (PI / 2.0 - el - rho).max(0.001);

    // Earth's surface area / coverage area of one satellite
    let n = (4.0 * PI) / (2.0 * PI * (1.0 - lambda.cos()));

    (n.ceil()) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_walker_delta_generation() {
        // Simple 6/3/1 constellation
        let c = WalkerConstellation::delta(6, 3, 1, R_EARTH_EQUATORIAL + 500.0, 45.0);
        let sats = c.generate();
        assert_eq!(sats.len(), 6);
        assert_eq!(c.sats_per_plane(), 2);

        // Check planes are equally spaced (delta: 180° spread)
        let raans: Vec<f64> = (0..3)
            .map(|p| sats.iter().find(|s| s.plane == p).unwrap().raan)
            .collect();
        let delta_raan = raans[1] - raans[0];
        assert!(
            (delta_raan - PI / 3.0).abs() < 1e-10,
            "RAAN spacing should be 60° for 3 planes in delta pattern"
        );
    }

    #[test]
    fn test_walker_star_generation() {
        // 8/4/0 star constellation
        let c = WalkerConstellation::star(8, 4, 0, R_EARTH_EQUATORIAL + 1000.0, 90.0);
        let sats = c.generate();
        assert_eq!(sats.len(), 8);

        // Check planes are equally spaced (star: 360° spread)
        let raans: Vec<f64> = (0..4)
            .map(|p| sats.iter().find(|s| s.plane == p).unwrap().raan)
            .collect();
        let delta_raan = raans[1] - raans[0];
        assert!(
            (delta_raan - PI / 2.0).abs() < 1e-10,
            "RAAN spacing should be 90° for 4 planes in star pattern"
        );
    }

    #[test]
    fn test_gps_constellation() {
        let gps = presets::gps();
        assert_eq!(gps.total, 24);
        assert_eq!(gps.planes, 6);

        let alt = gps.altitude();
        assert!(
            (alt - 20_180.0).abs() < 1.0,
            "GPS altitude should be ~20,180 km"
        );

        let period = gps.period();
        let period_hr = period / 3600.0;
        // GPS period ≈ 11.97 hours (half sidereal day)
        assert!(
            period_hr > 11.5 && period_hr < 12.5,
            "GPS period {period_hr} hr should be ~12 hr"
        );
    }

    #[test]
    fn test_iridium_constellation() {
        let iridium = presets::iridium();
        assert_eq!(iridium.total, 66);
        assert_eq!(iridium.planes, 6);
        assert_eq!(iridium.sats_per_plane(), 11);

        let sats = iridium.generate();
        assert_eq!(sats.len(), 66);
    }

    #[test]
    fn test_starlink_shell1() {
        let sl = presets::starlink_shell1();
        assert_eq!(sl.total, 1584);
        let alt = sl.altitude();
        assert!((alt - 550.0).abs() < 1.0);
    }

    #[test]
    fn test_coverage_high_constellation() {
        // GPS should provide very high coverage fraction
        let gps = presets::gps();
        let frac = coverage_fraction(&gps, 10.0, 45.0);
        assert!(
            frac > 0.9,
            "GPS should provide >90% coverage at 10° elevation"
        );
    }

    #[test]
    fn test_min_sats_leo() {
        // At 780 km (Iridium altitude), min elevation 8.2°
        let n = min_sats_for_coverage(780.0, 8.2, 86.4);
        // Iridium uses 66 sats; minimum should be in the ballpark
        assert!(
            n > 20 && n < 200,
            "Minimum sats at 780 km should be reasonable, got {n}"
        );
    }

    #[test]
    fn test_period_altitude_consistency() {
        let c = WalkerConstellation::delta(6, 3, 0, R_EARTH_EQUATORIAL + 400.0, 51.6);
        let period = c.period();
        // ISS-like orbit: ~92.65 min
        let period_min = period / 60.0;
        assert!(
            (period_min - 92.65).abs() < 1.0,
            "Period {period_min} min should be ~92.65 for 400 km orbit"
        );
    }
}
