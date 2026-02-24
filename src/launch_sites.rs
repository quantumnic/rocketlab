//! Launch site database and latitude effects on delta-v.
//!
//! Includes all major operational launch sites with geodetic coordinates,
//! azimuth constraints, and delta-v calculations accounting for Earth's rotation.
//!
//! References:
//! - Wertz, "Space Mission Engineering: The New SMAD", Ch. 18
//! - Curtis, "Orbital Mechanics for Engineering Students", Ch. 5

use crate::constants::*;

/// A launch site with its properties and constraints.
#[derive(Debug, Clone)]
pub struct LaunchSite {
    /// Site name
    pub name: &'static str,
    /// Country/operator
    pub country: &'static str,
    /// Geodetic latitude (degrees, positive north)
    pub latitude_deg: f64,
    /// Geodetic longitude (degrees, positive east)
    pub longitude_deg: f64,
    /// Minimum allowed launch azimuth (degrees from north, clockwise)
    pub min_azimuth_deg: f64,
    /// Maximum allowed launch azimuth (degrees from north, clockwise)
    pub max_azimuth_deg: f64,
    /// Minimum achievable inclination (degrees) — limited by latitude
    pub min_inclination_deg: f64,
}

impl LaunchSite {
    /// Earth's rotational velocity boost at this site (km/s).
    ///
    /// This is the eastward velocity component from Earth's rotation at the
    /// launch site's latitude: v = ω_Earth × R_eq × cos(lat)
    pub fn rotation_velocity(&self) -> f64 {
        let lat_rad = self.latitude_deg.to_radians();
        OMEGA_EARTH * R_EARTH_EQUATORIAL * lat_rad.cos()
    }

    /// Delta-v penalty/benefit for launching to a given inclination.
    ///
    /// Returns the velocity component from Earth's rotation that can be used
    /// for the given launch azimuth.
    ///
    /// # Arguments
    /// * `inclination_deg` — target orbit inclination (degrees)
    ///
    /// # Returns
    /// Rotational velocity benefit in km/s (positive = benefit from rotation)
    pub fn rotation_benefit(&self, inclination_deg: f64) -> f64 {
        let lat = self.latitude_deg.to_radians();
        let inc = inclination_deg.to_radians();

        // Launch azimuth for direct ascent to target inclination
        let sin_az = inc.cos() / lat.cos();
        if sin_az.abs() > 1.0 {
            // Can't reach this inclination directly from this latitude
            return 0.0;
        }
        let az = sin_az.asin();

        // Velocity boost = V_rot * sin(azimuth) (eastward component projected onto velocity)
        self.rotation_velocity() * az.sin()
    }

    /// Minimum delta-v to LEO (approximate, km/s) from this site.
    ///
    /// Uses empirical formula accounting for gravity losses (~1.5 km/s),
    /// drag losses (~0.1-0.3 km/s), and rotation benefit.
    ///
    /// # Arguments
    /// * `altitude_km` — target circular orbit altitude (km)
    /// * `inclination_deg` — target inclination (degrees)
    pub fn delta_v_to_orbit(&self, altitude_km: f64, inclination_deg: f64) -> f64 {
        let r = R_EARTH_EQUATORIAL + altitude_km;
        let v_orbit = (MU_EARTH / r).sqrt();

        // Gravity loss: typically 1.0-1.5 km/s depending on T/W
        let gravity_loss = 1.3;
        // Drag loss: typically 0.1-0.3 km/s
        let drag_loss = 0.2;
        // Steering loss
        let steering_loss = 0.1;

        let v_rot = self.rotation_benefit(inclination_deg);

        v_orbit + gravity_loss + drag_loss + steering_loss - v_rot
    }

    /// Achievable inclination range from this site.
    ///
    /// Returns (min_inclination, max_inclination) in degrees.
    /// The minimum is limited by latitude (direct ascent), and retrograde
    /// orbits (180° - lat) are also possible.
    pub fn inclination_range(&self) -> (f64, f64) {
        let lat = self.latitude_deg.abs();
        (lat, 180.0 - lat)
    }
}

/// Kennedy Space Center / Cape Canaveral, Florida
pub fn kennedy() -> LaunchSite {
    LaunchSite {
        name: "Kennedy Space Center / Cape Canaveral",
        country: "USA",
        latitude_deg: 28.524,
        longitude_deg: -80.651,
        min_azimuth_deg: 35.0,
        max_azimuth_deg: 120.0,
        min_inclination_deg: 28.5,
    }
}

/// Vandenberg Space Force Base, California
pub fn vandenberg() -> LaunchSite {
    LaunchSite {
        name: "Vandenberg Space Force Base",
        country: "USA",
        latitude_deg: 34.632,
        longitude_deg: -120.611,
        min_azimuth_deg: 147.0,
        max_azimuth_deg: 245.0,
        min_inclination_deg: 56.0, // Polar and SSO only
    }
}

/// Guiana Space Centre (CSG), Kourou, French Guiana
pub fn kourou() -> LaunchSite {
    LaunchSite {
        name: "Guiana Space Centre (CSG)",
        country: "France/ESA",
        latitude_deg: 5.232,
        longitude_deg: -52.775,
        min_azimuth_deg: -10.5,
        max_azimuth_deg: 93.5,
        min_inclination_deg: 5.2,
    }
}

/// Baikonur Cosmodrome, Kazakhstan
pub fn baikonur() -> LaunchSite {
    LaunchSite {
        name: "Baikonur Cosmodrome",
        country: "Kazakhstan/Russia",
        latitude_deg: 45.965,
        longitude_deg: 63.305,
        min_azimuth_deg: 45.0,
        max_azimuth_deg: 99.0,
        min_inclination_deg: 46.0,
    }
}

/// Tanegashima Space Center, Japan
pub fn tanegashima() -> LaunchSite {
    LaunchSite {
        name: "Tanegashima Space Center",
        country: "Japan",
        latitude_deg: 30.400,
        longitude_deg: 130.969,
        min_azimuth_deg: 70.0,
        max_azimuth_deg: 114.0,
        min_inclination_deg: 30.0,
    }
}

/// Satish Dhawan Space Centre (SHAR), Sriharikota, India
pub fn sriharikota() -> LaunchSite {
    LaunchSite {
        name: "Satish Dhawan Space Centre (SHAR)",
        country: "India",
        latitude_deg: 13.720,
        longitude_deg: 80.230,
        min_azimuth_deg: 92.0,
        max_azimuth_deg: 140.0,
        min_inclination_deg: 13.7,
    }
}

/// Plesetsk Cosmodrome, Russia
pub fn plesetsk() -> LaunchSite {
    LaunchSite {
        name: "Plesetsk Cosmodrome",
        country: "Russia",
        latitude_deg: 62.927,
        longitude_deg: 40.577,
        min_azimuth_deg: 0.0,
        max_azimuth_deg: 110.0,
        min_inclination_deg: 62.8,
    }
}

/// Jiuquan Satellite Launch Center, China
pub fn jiuquan() -> LaunchSite {
    LaunchSite {
        name: "Jiuquan Satellite Launch Center",
        country: "China",
        latitude_deg: 40.958,
        longitude_deg: 100.291,
        min_azimuth_deg: 60.0,
        max_azimuth_deg: 160.0,
        min_inclination_deg: 40.0,
    }
}

/// Wenchang Space Launch Site, Hainan, China
pub fn wenchang() -> LaunchSite {
    LaunchSite {
        name: "Wenchang Space Launch Site",
        country: "China",
        latitude_deg: 19.614,
        longitude_deg: 110.951,
        min_azimuth_deg: 50.0,
        max_azimuth_deg: 175.0,
        min_inclination_deg: 19.6,
    }
}

/// Rocket Lab Launch Complex 1, Mahia Peninsula, New Zealand
pub fn mahia() -> LaunchSite {
    LaunchSite {
        name: "Rocket Lab LC-1, Mahia Peninsula",
        country: "New Zealand",
        latitude_deg: -39.262,
        longitude_deg: 177.864,
        min_azimuth_deg: 0.0,
        max_azimuth_deg: 360.0,
        min_inclination_deg: 39.0,
    }
}

/// SpaceX Starbase, Boca Chica, Texas
pub fn starbase() -> LaunchSite {
    LaunchSite {
        name: "SpaceX Starbase, Boca Chica",
        country: "USA",
        latitude_deg: 25.997,
        longitude_deg: -97.157,
        min_azimuth_deg: 93.0,
        max_azimuth_deg: 113.0,
        min_inclination_deg: 26.0,
    }
}

/// Get all launch sites in the database.
pub fn all_sites() -> Vec<LaunchSite> {
    vec![
        kennedy(),
        vandenberg(),
        kourou(),
        baikonur(),
        tanegashima(),
        sriharikota(),
        plesetsk(),
        jiuquan(),
        wenchang(),
        mahia(),
        starbase(),
    ]
}

/// Compute launch azimuth for a direct ascent to target inclination.
///
/// # Arguments
/// * `latitude_deg` — launch site latitude (degrees)
/// * `inclination_deg` — target inclination (degrees)
///
/// # Returns
/// Launch azimuth in degrees from north (ascending node), or None if unreachable.
///
/// # Reference
/// Curtis (2014), Eq. 5.22
pub fn launch_azimuth(latitude_deg: f64, inclination_deg: f64) -> Option<f64> {
    let lat = latitude_deg.to_radians();
    let inc = inclination_deg.to_radians();

    let sin_az = inc.cos() / lat.cos();
    if sin_az.abs() > 1.0 {
        return None; // Inclination not reachable from this latitude
    }

    Some(sin_az.asin().to_degrees())
}

/// Compare delta-v cost from multiple sites to a given orbit.
///
/// Returns a sorted vector of (site_name, delta_v) from cheapest to most expensive.
pub fn compare_sites(altitude_km: f64, inclination_deg: f64) -> Vec<(&'static str, f64)> {
    let mut results: Vec<(&'static str, f64)> = all_sites()
        .iter()
        .filter(|s| {
            let (min_i, max_i) = s.inclination_range();
            inclination_deg >= min_i && inclination_deg <= max_i
        })
        .map(|s| (s.name, s.delta_v_to_orbit(altitude_km, inclination_deg)))
        .collect();

    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rotation_velocity() {
        // Equatorial site: v_rot ≈ 0.465 km/s
        let kourou = kourou();
        let v = kourou.rotation_velocity();
        assert!(
            (v - 0.465).abs() < 0.01,
            "Kourou rotation velocity {v} should be ~0.465 km/s"
        );

        // Higher latitude = less benefit
        let baikonur = baikonur();
        assert!(baikonur.rotation_velocity() < kourou.rotation_velocity());
    }

    #[test]
    fn test_kennedy_leo_delta_v() {
        // Kennedy to 200 km LEO at 28.5°: typically ~9.3-9.5 km/s total
        let ksc = kennedy();
        let dv = ksc.delta_v_to_orbit(200.0, 28.5);
        assert!(
            dv > 9.0 && dv < 10.0,
            "KSC to LEO delta-v {dv} should be ~9.3-9.5 km/s"
        );
    }

    #[test]
    fn test_kourou_advantage() {
        // Kourou should have the lowest delta-v for equatorial orbits
        let results = compare_sites(200.0, 5.5);
        assert!(!results.is_empty());
        assert_eq!(
            results[0].0, "Guiana Space Centre (CSG)",
            "Kourou should be cheapest for equatorial orbit"
        );
    }

    #[test]
    fn test_launch_azimuth_iss() {
        // KSC to ISS (51.6°): azimuth ≈ 44° (northeast)
        let az = launch_azimuth(28.524, 51.6).unwrap();
        assert!(
            az > 35.0 && az < 55.0,
            "KSC to ISS azimuth {az} should be ~44°"
        );
    }

    #[test]
    fn test_inclination_range() {
        let ksc = kennedy();
        let (min_i, max_i) = ksc.inclination_range();
        assert!((min_i - 28.524).abs() < 0.01);
        assert!((max_i - 151.476).abs() < 0.01);
    }

    #[test]
    fn test_polar_orbit_from_vandenberg() {
        let vafb = vandenberg();
        let dv = vafb.delta_v_to_orbit(500.0, 90.0);
        // No rotation benefit for polar orbit
        assert!(
            dv > 9.0 && dv < 10.5,
            "Vandenberg to polar orbit delta-v {dv} should be ~9.5 km/s"
        );
    }

    #[test]
    fn test_unreachable_inclination() {
        // Baikonur can't do equatorial orbits directly
        let az = launch_azimuth(45.965, 10.0);
        assert!(
            az.is_none(),
            "10° inclination should be unreachable from Baikonur"
        );
    }
}
