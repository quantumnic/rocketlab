//! Interplanetary trajectory design.
//!
//! Patched-conic method for interplanetary transfer computation:
//! - Synodic period and launch window calculation
//! - Heliocentric transfer orbit (Hohmann approximation)
//! - Departure and arrival hyperbolic excess velocities
//! - C3 (characteristic energy) computation
//! - Phase angle requirements
//!
//! References:
//! - Curtis, "Orbital Mechanics for Engineering Students", Ch. 8
//! - Bate, Mueller & White, §8.6–8.8
//! - Vallado, "Fundamentals of Astrodynamics and Applications", Ch. 12

use crate::constants::*;

/// Result of an interplanetary transfer computation.
#[derive(Debug, Clone)]
pub struct InterplanetaryTransfer {
    /// Name of departure body
    pub departure_body: &'static str,
    /// Name of arrival body
    pub arrival_body: &'static str,
    /// Heliocentric transfer semi-major axis (km)
    pub a_transfer: f64,
    /// Transfer orbit eccentricity
    pub e_transfer: f64,
    /// Time of flight (seconds)
    pub tof: f64,
    /// Time of flight (days)
    pub tof_days: f64,
    /// Departure v-infinity (km/s)
    pub v_inf_departure: f64,
    /// Arrival v-infinity (km/s)
    pub v_inf_arrival: f64,
    /// Departure C3 (km²/s²)
    pub c3_departure: f64,
    /// Arrival C3 (km²/s²)
    pub c3_arrival: f64,
    /// Required phase angle at departure (radians)
    pub phase_angle: f64,
    /// Total delta-v if departing from circular parking orbit and
    /// capturing into circular orbit at destination (km/s)
    pub total_delta_v: f64,
}

/// Planet data for interplanetary calculations.
#[derive(Debug, Clone, Copy)]
pub struct PlanetData {
    pub name: &'static str,
    /// Semi-major axis of orbit around Sun (km)
    pub a_sun: f64,
    /// Gravitational parameter (km³/s²)
    pub mu: f64,
    /// Mean body radius (km)
    pub radius: f64,
    /// Typical parking orbit altitude (km)
    pub parking_alt: f64,
}

/// Standard planet data.
pub const EARTH_DATA: PlanetData = PlanetData {
    name: "Earth",
    a_sun: A_EARTH,
    mu: MU_EARTH,
    radius: R_EARTH,
    parking_alt: 200.0,
};

pub const MARS_DATA: PlanetData = PlanetData {
    name: "Mars",
    a_sun: A_MARS,
    mu: MU_MARS,
    radius: R_MARS,
    parking_alt: 300.0,
};

pub const VENUS_DATA: PlanetData = PlanetData {
    name: "Venus",
    a_sun: A_VENUS,
    mu: MU_VENUS,
    radius: R_VENUS,
    parking_alt: 300.0,
};

pub const JUPITER_DATA: PlanetData = PlanetData {
    name: "Jupiter",
    a_sun: A_JUPITER,
    mu: MU_JUPITER,
    radius: R_JUPITER,
    parking_alt: 1000.0,
};

/// Compute synodic period between two planets (seconds).
///
/// T_syn = 1 / |1/T1 - 1/T2|
///
/// where T1, T2 are orbital periods around the Sun.
pub fn synodic_period(a1: f64, a2: f64) -> f64 {
    let t1 = TWO_PI * (a1.powi(3) / MU_SUN).sqrt();
    let t2 = TWO_PI * (a2.powi(3) / MU_SUN).sqrt();
    let n1 = TWO_PI / t1;
    let n2 = TWO_PI / t2;
    TWO_PI / (n1 - n2).abs()
}

/// Compute synodic period in days.
pub fn synodic_period_days(a1: f64, a2: f64) -> f64 {
    synodic_period(a1, a2) / 86400.0
}

/// Hohmann-like interplanetary transfer using patched conics.
///
/// Assumes coplanar circular orbits for both planets (first approximation).
/// Computes departure and arrival v-infinity, C3, phase angle, and total delta-v
/// including escape from parking orbit and capture at destination.
pub fn hohmann_interplanetary(
    departure: &PlanetData,
    arrival: &PlanetData,
) -> InterplanetaryTransfer {
    let r1 = departure.a_sun;
    let r2 = arrival.a_sun;

    // Transfer orbit parameters
    let a_t = (r1 + r2) / 2.0;
    let e_t = (r2 - r1).abs() / (r1 + r2);

    // Time of flight (half the transfer orbit period)
    let tof = PI * (a_t.powi(3) / MU_SUN).sqrt();

    // Heliocentric velocities of planets (circular)
    let v_planet_dep = (MU_SUN / r1).sqrt();
    let v_planet_arr = (MU_SUN / r2).sqrt();

    // Transfer orbit velocities at departure and arrival (vis-viva)
    let v_dep_helio = (MU_SUN * (2.0 / r1 - 1.0 / a_t)).sqrt();
    let v_arr_helio = (MU_SUN * (2.0 / r2 - 1.0 / a_t)).sqrt();

    // V-infinity at departure and arrival
    let v_inf_dep = (v_dep_helio - v_planet_dep).abs();
    let v_inf_arr = (v_arr_helio - v_planet_arr).abs();

    // C3 = v_inf²
    let c3_dep = v_inf_dep * v_inf_dep;
    let c3_arr = v_inf_arr * v_inf_arr;

    // Phase angle at departure
    // α = π - n_target * TOF, where n_target is target's mean motion
    let n_arr = (MU_SUN / arrival.a_sun.powi(3)).sqrt();
    let phase_angle = PI - n_arr * tof;
    // Normalize to [0, 2π)
    let phase_angle = ((phase_angle % TWO_PI) + TWO_PI) % TWO_PI;

    // Delta-v from parking orbit (departure)
    let r_park_dep = departure.radius + departure.parking_alt;
    let v_park_dep = (departure.mu / r_park_dep).sqrt();
    let v_escape_dep = (v_inf_dep * v_inf_dep + 2.0 * departure.mu / r_park_dep).sqrt();
    let dv_dep = v_escape_dep - v_park_dep;

    // Delta-v for capture into parking orbit (arrival)
    let r_park_arr = arrival.radius + arrival.parking_alt;
    let v_park_arr = (arrival.mu / r_park_arr).sqrt();
    let v_arrival_hyp = (v_inf_arr * v_inf_arr + 2.0 * arrival.mu / r_park_arr).sqrt();
    let dv_arr = v_arrival_hyp - v_park_arr;

    InterplanetaryTransfer {
        departure_body: departure.name,
        arrival_body: arrival.name,
        a_transfer: a_t,
        e_transfer: e_t,
        tof,
        tof_days: tof / 86400.0,
        v_inf_departure: v_inf_dep,
        v_inf_arrival: v_inf_arr,
        c3_departure: c3_dep,
        c3_arrival: c3_arr,
        phase_angle,
        total_delta_v: dv_dep + dv_arr,
    }
}

/// Compute departure delta-v from a parking orbit given C3.
///
/// Δv = sqrt(v_inf² + 2*μ/r) - sqrt(μ/r)
pub fn departure_dv_from_c3(c3: f64, mu_body: f64, r_parking: f64) -> f64 {
    let v_park = (mu_body / r_parking).sqrt();
    let v_hyp = (c3 + 2.0 * mu_body / r_parking).sqrt();
    v_hyp - v_park
}

/// Compute capture delta-v into a parking orbit given v-infinity.
pub fn capture_dv(v_inf: f64, mu_body: f64, r_parking: f64) -> f64 {
    let v_park = (mu_body / r_parking).sqrt();
    let v_hyp = (v_inf * v_inf + 2.0 * mu_body / r_parking).sqrt();
    v_hyp - v_park
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_earth_mars_synodic_period() {
        // Earth-Mars synodic period ≈ 780 days (2.135 years)
        let t_syn = synodic_period_days(A_EARTH, A_MARS);
        assert!(
            (t_syn - 780.0).abs() < 10.0,
            "Earth-Mars synodic period = {t_syn} days, expected ~780"
        );
    }

    #[test]
    fn test_earth_venus_synodic_period() {
        // Earth-Venus synodic period ≈ 584 days
        let t_syn = synodic_period_days(A_EARTH, A_VENUS);
        assert!(
            (t_syn - 584.0).abs() < 10.0,
            "Earth-Venus synodic period = {t_syn} days, expected ~584"
        );
    }

    #[test]
    fn test_earth_mars_hohmann() {
        let transfer = hohmann_interplanetary(&EARTH_DATA, &MARS_DATA);

        // TOF should be ~259 days (8.5 months)
        assert!(
            (transfer.tof_days - 259.0).abs() < 10.0,
            "Earth-Mars TOF = {} days, expected ~259",
            transfer.tof_days
        );

        // Departure v_inf ~ 2.94 km/s
        assert!(
            (transfer.v_inf_departure - 2.94).abs() < 0.2,
            "Departure v_inf = {} km/s, expected ~2.94",
            transfer.v_inf_departure
        );

        // C3 ~ 8.65 km²/s²
        assert!(
            (transfer.c3_departure - 8.65).abs() < 1.0,
            "Departure C3 = {} km²/s², expected ~8.65",
            transfer.c3_departure
        );

        // Total delta-v ~ 5.6 km/s (departure + capture)
        assert!(
            transfer.total_delta_v > 4.5 && transfer.total_delta_v < 7.0,
            "Total delta-v = {} km/s, expected ~5.6",
            transfer.total_delta_v
        );
    }

    #[test]
    fn test_earth_venus_hohmann() {
        let transfer = hohmann_interplanetary(&EARTH_DATA, &VENUS_DATA);

        // TOF ~ 146 days
        assert!(
            (transfer.tof_days - 146.0).abs() < 10.0,
            "Earth-Venus TOF = {} days, expected ~146",
            transfer.tof_days
        );
    }

    #[test]
    fn test_earth_jupiter_hohmann() {
        let transfer = hohmann_interplanetary(&EARTH_DATA, &JUPITER_DATA);

        // TOF ~ 997 days (2.73 years)
        assert!(
            (transfer.tof_days - 997.0).abs() < 20.0,
            "Earth-Jupiter TOF = {} days, expected ~997",
            transfer.tof_days
        );
    }

    #[test]
    fn test_phase_angle_mars() {
        let transfer = hohmann_interplanetary(&EARTH_DATA, &MARS_DATA);
        // Phase angle for Mars ~ 44° (0.77 rad)
        let phase_deg = transfer.phase_angle.to_degrees();
        assert!(
            (phase_deg - 44.0).abs() < 5.0,
            "Mars phase angle = {phase_deg}°, expected ~44°"
        );
    }

    #[test]
    fn test_departure_dv_from_c3() {
        // LEO 200 km, C3 = 10 km²/s²
        let r_park = R_EARTH + 200.0;
        let dv = departure_dv_from_c3(10.0, MU_EARTH, r_park);
        // Should be ~3.6 km/s
        assert!(
            dv > 3.0 && dv < 4.5,
            "Departure dv from C3=10 should be ~3.6 km/s, got {dv}"
        );
    }

    #[test]
    fn test_transfer_eccentricity() {
        let transfer = hohmann_interplanetary(&EARTH_DATA, &MARS_DATA);
        // e = (r2-r1)/(r1+r2) ≈ 0.208
        assert!(
            transfer.e_transfer > 0.15 && transfer.e_transfer < 0.25,
            "Transfer eccentricity = {}, expected ~0.208",
            transfer.e_transfer
        );
    }
}
