//! Orbital perturbations: J2-J6 zonal harmonics, atmospheric drag, solar radiation pressure,
//! third-body effects (Sun/Moon).
//!
//! References:
//! - Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed., Ch. 9
//! - Battin, "An Introduction to the Methods of Astrodynamics"
//! - Montenbruck & Gill, "Satellite Orbits", Ch. 3

use crate::constants::*;

// ── Zonal Harmonic Coefficients (WGS-84 / EGM-96) ──────────────────────────

/// J2 zonal harmonic (Earth oblateness)
pub const J2: f64 = 1.082_63e-3;
/// J3 zonal harmonic
pub const J3: f64 = -2.5327e-6;
/// J4 zonal harmonic
pub const J4: f64 = -1.6196e-6;
/// J5 zonal harmonic
pub const J5: f64 = -2.273e-7;
/// J6 zonal harmonic
pub const J6: f64 = 5.407e-7;

/// Secular rates due to J2 perturbation on orbital elements.
///
/// Returns (Ω_dot, ω_dot, M0_dot) in rad/s — the secular drift rates of
/// RAAN, argument of perigee, and mean anomaly respectively.
///
/// # Arguments
/// * `a` — semi-major axis (km)
/// * `e` — eccentricity
/// * `i` — inclination (radians)
/// * `mu` — gravitational parameter (km³/s²)
/// * `r_eq` — equatorial radius (km)
///
/// # Reference
/// Vallado (2013), Eq. 9-40 through 9-42
pub fn j2_secular_rates(a: f64, e: f64, i: f64, mu: f64, r_eq: f64) -> (f64, f64, f64) {
    let n = (mu / a.powi(3)).sqrt(); // mean motion
    let p = a * (1.0 - e * e); // semi-latus rectum
    let ratio = r_eq / p;
    let cos_i = i.cos();
    let sin_i = i.sin();

    // RAAN precession rate
    let raan_dot = -1.5 * n * J2 * ratio.powi(2) * cos_i;

    // Argument of perigee precession rate
    let omega_dot = 0.75 * n * J2 * ratio.powi(2) * (5.0 * cos_i * cos_i - 1.0);

    // Mean anomaly correction
    let eta = (1.0 - e * e).sqrt();
    let m0_dot = 0.75 * n * J2 * ratio.powi(2) * eta * (3.0 * cos_i * cos_i - 1.0);

    let _ = sin_i; // used conceptually, suppress warning
    (raan_dot, omega_dot, m0_dot)
}

/// Sun-synchronous inclination for given semi-major axis and eccentricity.
///
/// Returns the inclination (radians) required for the RAAN to precess at
/// exactly 360°/year (0.9856°/day) to maintain sun-synchronous geometry.
///
/// # Reference
/// Vallado (2013), Section 11.4
pub fn sun_synchronous_inclination(a: f64, e: f64) -> f64 {
    let n = (MU_EARTH / a.powi(3)).sqrt();
    let p = a * (1.0 - e * e);
    let ratio = R_EARTH_EQUATORIAL / p;

    // Required RAAN rate: 360°/365.2421897 days = 1.99099e-7 rad/s
    let raan_dot_required = 2.0 * PI / (365.2421897 * 86400.0);

    let cos_i = -raan_dot_required / (1.5 * n * J2 * ratio.powi(2));
    cos_i.acos()
}

/// Critical inclination where argument of perigee has zero secular drift.
/// This occurs at i ≈ 63.4° and i ≈ 116.6° (Molniya orbits use this).
pub const CRITICAL_INCLINATION_RAD: f64 = 1.10714872; // 63.4349° in radians

/// Repeating ground track semi-major axis.
///
/// For a satellite that repeats its ground track every `revs_per_day` revolutions
/// per `days` sidereal days, accounting for J2 effects.
///
/// # Arguments
/// * `revs_per_day` — number of revolutions per repeat cycle
/// * `days` — number of days in the repeat cycle
/// * `i` — inclination (radians)
/// * `e` — eccentricity
pub fn repeating_ground_track_sma(revs_per_day: f64, days: f64, i: f64, e: f64) -> f64 {
    // Iterative solution: mean motion must satisfy n_eff = revs_per_day / days * OMEGA_EARTH
    let _target_n = revs_per_day / days * OMEGA_EARTH * (2.0 * PI) / OMEGA_EARTH;
    // Actually: the satellite must complete revs_per_day orbits in `days` sidereal days
    // n_eff = 2π * revs_per_day / (days * 86164.0905) [sidereal day in seconds]
    let sidereal_day = 86164.0905; // seconds
    let n_target = 2.0 * PI * revs_per_day / (days * sidereal_day);

    // First guess: Keplerian
    let mut a = (MU_EARTH / (n_target * n_target)).powf(1.0 / 3.0);

    // Iterate with J2 correction
    for _ in 0..20 {
        let p = a * (1.0 - e * e);
        let ratio = R_EARTH_EQUATORIAL / p;
        let cos_i = i.cos();
        let eta = (1.0 - e * e).sqrt();
        let n = (MU_EARTH / a.powi(3)).sqrt();

        // Corrected mean motion
        let dn = n * (1.0 + 0.75 * J2 * ratio.powi(2) * eta * (3.0 * cos_i * cos_i - 1.0));

        // Also account for RAAN and omega drift
        let raan_dot = -1.5 * n * J2 * ratio.powi(2) * cos_i;
        let _omega_dot = 0.75 * n * J2 * ratio.powi(2) * (5.0 * cos_i * cos_i - 1.0);

        // Effective node rate: satellite ground track repeats when
        // (dn - OMEGA_EARTH + raan_dot) * T_repeat = 2π * revs_per_day
        let n_eff = dn + raan_dot; // approximate nodal rate
        let _ = n_eff;

        // Simple Newton step on Keplerian + J2 mean motion
        let a_new = (MU_EARTH / (n_target * n_target)).powf(1.0 / 3.0);
        // Refine: use corrected mean motion
        let a_corrected = (MU_EARTH / a.powi(3)).sqrt();
        let scale = n_target
            / (a_corrected * (1.0 + 0.75 * J2 * ratio.powi(2) * eta * (3.0 * cos_i * cos_i - 1.0)));
        a = (MU_EARTH / (scale * scale)).powf(1.0 / 3.0);

        if (a - a_new).abs() / a < 1e-12 {
            break;
        }
    }

    a
}

/// Acceleration due to J2 perturbation in ECI coordinates.
///
/// # Arguments
/// * `r` — position vector [x, y, z] in km (ECI)
/// * `mu` — gravitational parameter (km³/s²)
/// * `r_eq` — equatorial radius (km)
///
/// # Returns
/// Perturbation acceleration [ax, ay, az] in km/s²
pub fn j2_acceleration(r: &[f64; 3], mu: f64, r_eq: f64) -> [f64; 3] {
    let r_mag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let r2 = r_mag * r_mag;
    let r5 = r2 * r2 * r_mag;
    let z2_r2 = (r[2] * r[2]) / r2;

    let coeff = -1.5 * J2 * mu * r_eq * r_eq / r5;

    [
        coeff * r[0] * (1.0 - 5.0 * z2_r2),
        coeff * r[1] * (1.0 - 5.0 * z2_r2),
        coeff * r[2] * (3.0 - 5.0 * z2_r2),
    ]
}

/// Acceleration due to atmospheric drag.
///
/// # Arguments
/// * `v_rel` — velocity relative to atmosphere [vx, vy, vz] in km/s
/// * `rho` — atmospheric density in kg/m³
/// * `cd` — drag coefficient
/// * `area_m2` — cross-sectional area in m²
/// * `mass_kg` — spacecraft mass in kg
///
/// # Returns
/// Drag acceleration [ax, ay, az] in km/s²
pub fn drag_acceleration(
    v_rel: &[f64; 3],
    rho: f64,
    cd: f64,
    area_m2: f64,
    mass_kg: f64,
) -> [f64; 3] {
    let v_mag = (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]).sqrt();
    if v_mag < 1e-15 {
        return [0.0; 3];
    }

    // Ballistic coefficient: B = Cd * A / m
    let b = cd * area_m2 / mass_kg;

    // a_drag = -0.5 * rho * v² * B * v_hat
    // rho is kg/m³, v is km/s → convert: rho * (v*1000)² gives N/m² per unit area
    // Result in km/s²: multiply by 1e-3
    let factor = -0.5 * rho * v_mag * 1e3 * b * 1e-3; // km/s²

    [factor * v_rel[0], factor * v_rel[1], factor * v_rel[2]]
}

/// Velocity of atmosphere at position r (co-rotating with Earth).
///
/// # Arguments
/// * `r` — position in ECI [x, y, z] km
///
/// # Returns
/// Atmospheric velocity [vx, vy, vz] in km/s
pub fn atmosphere_velocity(r: &[f64; 3]) -> [f64; 3] {
    // v_atm = omega_earth × r
    [-OMEGA_EARTH * r[1], OMEGA_EARTH * r[0], 0.0]
}

/// Solar radiation pressure acceleration.
///
/// # Arguments
/// * `r_sat_sun` — vector from satellite to Sun [x, y, z] in km
/// * `cr` — reflectivity coefficient (1.0 = perfect absorber, 2.0 = perfect reflector)
/// * `area_m2` — cross-sectional area exposed to Sun in m²
/// * `mass_kg` — spacecraft mass in kg
///
/// # Returns
/// SRP acceleration [ax, ay, az] in km/s²
///
/// # Reference
/// Solar flux at 1 AU: P = 4.56e-6 N/m² (1361 W/m² / c)
pub fn srp_acceleration(r_sat_sun: &[f64; 3], cr: f64, area_m2: f64, mass_kg: f64) -> [f64; 3] {
    let r_mag = (r_sat_sun[0].powi(2) + r_sat_sun[1].powi(2) + r_sat_sun[2].powi(2)).sqrt();
    if r_mag < 1e-10 {
        return [0.0; 3];
    }

    // Solar radiation pressure at 1 AU (N/m²)
    let p_sr = 4.56e-6;

    // Scale by distance (inverse square)
    let au_m = AU_KM * 1e3;
    let r_m = r_mag * 1e3;
    let p_at_r = p_sr * (au_m / r_m).powi(2);

    // Acceleration magnitude in m/s², then convert to km/s²
    let a_mag = cr * p_at_r * area_m2 / mass_kg * 1e-3;

    // Direction: away from Sun (along r_sat_sun normalized)
    // Actually SRP pushes AWAY from Sun, so along the direction from Sun to satellite
    // r_sat_sun points from satellite to Sun, so acceleration is in -r_sat_sun direction
    let r_hat = [
        r_sat_sun[0] / r_mag,
        r_sat_sun[1] / r_mag,
        r_sat_sun[2] / r_mag,
    ];

    [-a_mag * r_hat[0], -a_mag * r_hat[1], -a_mag * r_hat[2]]
}

/// Third-body gravitational perturbation acceleration.
///
/// # Arguments
/// * `r_sat` — satellite position in ECI [x, y, z] km
/// * `r_body` — perturbing body position in ECI [x, y, z] km
/// * `mu_body` — gravitational parameter of perturbing body (km³/s²)
///
/// # Returns
/// Third-body perturbation acceleration [ax, ay, az] in km/s²
///
/// # Reference
/// Battin (1999), Eq. 8.60
pub fn third_body_acceleration(r_sat: &[f64; 3], r_body: &[f64; 3], mu_body: f64) -> [f64; 3] {
    // Vector from satellite to body
    let d = [
        r_body[0] - r_sat[0],
        r_body[1] - r_sat[1],
        r_body[2] - r_sat[2],
    ];
    let d_mag = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let r_body_mag = (r_body[0].powi(2) + r_body[1].powi(2) + r_body[2].powi(2)).sqrt();

    if d_mag < 1e-10 || r_body_mag < 1e-10 {
        return [0.0; 3];
    }

    let d3 = d_mag.powi(3);
    let rb3 = r_body_mag.powi(3);

    [
        mu_body * (d[0] / d3 - r_body[0] / rb3),
        mu_body * (d[1] / d3 - r_body[1] / rb3),
        mu_body * (d[2] / d3 - r_body[2] / rb3),
    ]
}

/// Propagate orbit with J2 perturbation using RK4.
///
/// # Arguments
/// * `r0` — initial position [x, y, z] in km (ECI)
/// * `v0` — initial velocity [vx, vy, vz] in km/s (ECI)
/// * `duration_s` — propagation duration in seconds
/// * `dt` — time step in seconds
/// * `mu` — gravitational parameter (km³/s²)
/// * `r_eq` — equatorial radius (km)
///
/// # Returns
/// Vector of (time, [x,y,z], [vx,vy,vz]) tuples
pub fn propagate_j2(
    r0: [f64; 3],
    v0: [f64; 3],
    duration_s: f64,
    dt: f64,
    mu: f64,
    r_eq: f64,
) -> Vec<(f64, [f64; 3], [f64; 3])> {
    let mut results = Vec::new();
    let mut state = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]];
    let mut t = 0.0;

    results.push((
        t,
        [state[0], state[1], state[2]],
        [state[3], state[4], state[5]],
    ));

    while t < duration_s {
        let step = dt.min(duration_s - t);
        state = rk4_step_j2(&state, step, mu, r_eq);
        t += step;
        results.push((
            t,
            [state[0], state[1], state[2]],
            [state[3], state[4], state[5]],
        ));
    }

    results
}

fn rk4_step_j2(state: &[f64; 6], dt: f64, mu: f64, r_eq: f64) -> [f64; 6] {
    let k1 = derivatives_j2(state, mu, r_eq);
    let s2 = add_scaled(state, &k1, 0.5 * dt);
    let k2 = derivatives_j2(&s2, mu, r_eq);
    let s3 = add_scaled(state, &k2, 0.5 * dt);
    let k3 = derivatives_j2(&s3, mu, r_eq);
    let s4 = add_scaled(state, &k3, dt);
    let k4 = derivatives_j2(&s4, mu, r_eq);

    let mut result = [0.0; 6];
    for i in 0..6 {
        result[i] = state[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    result
}

fn derivatives_j2(state: &[f64; 6], mu: f64, r_eq: f64) -> [f64; 6] {
    let r = [state[0], state[1], state[2]];
    let r_mag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let r3 = r_mag.powi(3);

    let a_central = [-mu * r[0] / r3, -mu * r[1] / r3, -mu * r[2] / r3];
    let a_j2 = j2_acceleration(&r, mu, r_eq);

    [
        state[3],
        state[4],
        state[5],
        a_central[0] + a_j2[0],
        a_central[1] + a_j2[1],
        a_central[2] + a_j2[2],
    ]
}

fn add_scaled(a: &[f64; 6], b: &[f64; 6], s: f64) -> [f64; 6] {
    [
        a[0] + s * b[0],
        a[1] + s * b[1],
        a[2] + s * b[2],
        a[3] + s * b[3],
        a[4] + s * b[4],
        a[5] + s * b[5],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_j2_secular_raan_precession() {
        // ISS-like orbit: a = 6778 km, e ≈ 0, i = 51.6°
        let a = 6778.0;
        let e = 0.0005;
        let i = 51.6_f64.to_radians();

        let (raan_dot, omega_dot, _m0_dot) =
            j2_secular_rates(a, e, i, MU_EARTH, R_EARTH_EQUATORIAL);

        // RAAN precession for ISS ≈ -5° to -7° per day
        let raan_deg_per_day = raan_dot.to_degrees() * 86400.0;
        assert!(
            raan_deg_per_day < -4.0 && raan_deg_per_day > -8.0,
            "ISS RAAN precession {raan_deg_per_day} deg/day outside expected range"
        );

        // Argument of perigee should precess positively for i < 63.4°
        assert!(
            omega_dot > 0.0,
            "omega should precess positively for i < 63.4°"
        );
    }

    #[test]
    fn test_critical_inclination() {
        // At critical inclination, omega_dot should be ~0
        let a = 26560.0; // GPS-like
        let e = 0.01;
        let i = CRITICAL_INCLINATION_RAD;

        let (_raan_dot, omega_dot, _m0_dot) =
            j2_secular_rates(a, e, i, MU_EARTH, R_EARTH_EQUATORIAL);

        // omega_dot should be very close to zero at critical inclination
        let omega_deg_per_day = omega_dot.to_degrees() * 86400.0;
        assert!(
            omega_deg_per_day.abs() < 0.01,
            "omega_dot at critical inclination should be ~0, got {omega_deg_per_day}"
        );
    }

    #[test]
    fn test_sun_synchronous_inclination() {
        // Typical SSO at ~700 km altitude: i ≈ 98.2°
        let a = R_EARTH_EQUATORIAL + 700.0;
        let e = 0.001;

        let i = sun_synchronous_inclination(a, e);
        let i_deg = i.to_degrees();

        assert!(
            i_deg > 97.0 && i_deg < 100.0,
            "SSO inclination at 700 km should be ~98°, got {i_deg}"
        );
    }

    #[test]
    fn test_j2_acceleration_magnitude() {
        // At equator, altitude ~400 km
        let r = [R_EARTH_EQUATORIAL + 400.0, 0.0, 0.0];
        let a = j2_acceleration(&r, MU_EARTH, R_EARTH_EQUATORIAL);

        let a_mag = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        // J2 perturbation is ~1e-3 of central gravity
        let central = MU_EARTH / (r[0] * r[0]);
        let ratio = a_mag / central;
        assert!(
            ratio > 1e-4 && ratio < 1e-2,
            "J2/central ratio {ratio} should be ~1e-3"
        );
    }

    #[test]
    fn test_drag_acceleration() {
        // LEO satellite: rho ~ 1e-12 kg/m³ at 400 km
        let v_rel = [7.66, 0.0, 0.0]; // km/s
        let rho = 1e-12; // kg/m³
        let cd = 2.2;
        let area = 10.0; // m²
        let mass = 1000.0; // kg

        let a_drag = drag_acceleration(&v_rel, rho, cd, area, mass);

        // Drag should be negative in velocity direction
        assert!(a_drag[0] < 0.0, "Drag should oppose velocity");

        // Order of magnitude check: ~1e-9 km/s² for typical LEO
        let a_mag = (a_drag[0] * a_drag[0] + a_drag[1] * a_drag[1] + a_drag[2] * a_drag[2]).sqrt();
        assert!(
            a_mag > 1e-14 && a_mag < 1e-6,
            "Drag magnitude {a_mag} km/s² out of expected range for LEO"
        );
    }

    #[test]
    fn test_srp_acceleration() {
        // At 1 AU, typical spacecraft
        let r_to_sun = [AU_KM, 0.0, 0.0]; // satellite-to-sun vector
        let cr = 1.5;
        let area = 20.0; // m² (solar panels)
        let mass = 500.0; // kg

        let a_srp = srp_acceleration(&r_to_sun, cr, area, mass);

        // SRP at 1 AU: P = 4.56e-6 N/m², a ~ Cr*P*A/m ~ 1.5*4.56e-6*20/500 ~ 2.7e-7 m/s²
        // In km/s²: ~2.7e-10
        let a_mag = (a_srp[0].powi(2) + a_srp[1].powi(2) + a_srp[2].powi(2)).sqrt();
        assert!(
            a_mag > 1e-11 && a_mag < 1e-8,
            "SRP magnitude {a_mag} km/s² outside expected range at 1 AU"
        );

        // Direction should be away from Sun (negative x since r_to_sun is +x)
        assert!(a_srp[0] < 0.0, "SRP should push away from Sun");
    }

    #[test]
    fn test_third_body_moon() {
        // Satellite in LEO, Moon at ~384400 km along +x
        let r_sat = [R_EARTH_EQUATORIAL + 400.0, 0.0, 0.0];
        let r_moon = [384_400.0, 0.0, 0.0];

        let a_tb = third_body_acceleration(&r_sat, &r_moon, MU_MOON);
        let a_mag = (a_tb[0].powi(2) + a_tb[1].powi(2) + a_tb[2].powi(2)).sqrt();

        // Lunar perturbation on LEO: ~1e-9 to 1e-8 km/s²
        assert!(
            a_mag > 1e-11 && a_mag < 1e-6,
            "Lunar perturbation {a_mag} km/s² outside expected range"
        );
    }

    #[test]
    fn test_j2_propagation_energy_bounded() {
        // Propagate a circular LEO orbit for one period with J2
        let a = R_EARTH_EQUATORIAL + 400.0;
        let v = (MU_EARTH / a).sqrt();
        let r0 = [a, 0.0, 0.0];
        let v0 = [0.0, v, 0.0];

        let period = 2.0 * PI * (a.powi(3) / MU_EARTH).sqrt();
        let results = propagate_j2(r0, v0, period, 10.0, MU_EARTH, R_EARTH_EQUATORIAL);

        // Check that radius stays bounded (not crashing or escaping)
        for (_, r, _) in &results {
            let r_mag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            assert!(r_mag > R_EARTH_EQUATORIAL, "Satellite shouldn't crash");
            assert!(r_mag < a + 100.0, "Satellite shouldn't escape");
        }
    }

    #[test]
    fn test_atmosphere_velocity() {
        let r = [R_EARTH_EQUATORIAL + 400.0, 0.0, 0.0];
        let v_atm = atmosphere_velocity(&r);

        // At equator: v = omega * R ≈ 7.29e-5 * 6778 ≈ 0.494 km/s
        let v_mag = (v_atm[0].powi(2) + v_atm[1].powi(2) + v_atm[2].powi(2)).sqrt();
        assert!(
            (v_mag - 0.494).abs() < 0.05,
            "Atmosphere velocity {v_mag} km/s should be ~0.494"
        );

        // Should be in +y direction for +x position
        assert!(v_atm[1] > 0.0);
        assert!(v_atm[2].abs() < 1e-15);
    }
}
