//! Gravity assist / flyby mechanics.
//!
//! Models for planetary flyby trajectory design including:
//! - Hyperbolic excess velocity (v-infinity) computation
//! - Turn angle and periapsis calculations
//! - B-plane targeting (impact parameter)
//! - Post-flyby velocity via the crank angle
//! - Delta-v equivalent from gravity assist
//!
//! References:
//! - Bate, Mueller & White, "Fundamentals of Astrodynamics", Ch. 8
//! - Vallado, "Fundamentals of Astrodynamics and Applications", §12.6
//! - Strange & Longuski, "Graphical Method for Gravity-Assist Trajectory Design" (2002)

use nalgebra::Vector3;

/// Flyby geometry and result.
#[derive(Debug, Clone)]
pub struct FlybyResult {
    /// Turn angle (radians)
    pub turn_angle: f64,
    /// Periapsis radius (km)
    pub r_periapsis: f64,
    /// Semi-major axis of hyperbola (km, negative)
    pub a_hyp: f64,
    /// Eccentricity of hyperbola
    pub e_hyp: f64,
    /// B-plane impact parameter (km)
    pub b_parameter: f64,
    /// Incoming v-infinity magnitude (km/s)
    pub v_inf_in: f64,
    /// Outgoing v-infinity magnitude (km/s) — same as incoming for unpowered flyby
    pub v_inf_out: f64,
    /// Heliocentric velocity after flyby (km/s)
    pub v_helio_out: Vector3<f64>,
    /// Delta-v equivalent gained (km/s)
    pub delta_v_equivalent: f64,
}

/// Compute the turn angle for a hyperbolic flyby.
///
/// δ = 2 * arcsin(1/e)
///
/// where e is the eccentricity of the hyperbolic trajectory.
pub fn turn_angle_from_eccentricity(e: f64) -> f64 {
    assert!(e > 1.0, "Flyby requires hyperbolic orbit (e > 1)");
    2.0 * (1.0 / e).asin()
}

/// Compute hyperbolic eccentricity from periapsis and v-infinity.
///
/// e = 1 + r_p * v_inf² / mu
pub fn hyperbolic_eccentricity(r_periapsis: f64, v_inf: f64, mu_body: f64) -> f64 {
    1.0 + r_periapsis * v_inf * v_inf / mu_body
}

/// B-plane impact parameter (km).
///
/// b = (r_p / e) * sqrt(e² - 1) = a * sqrt(e² - 1)
///
/// where a = -mu / v_inf² (semi-major axis of hyperbola, negative).
pub fn b_plane_parameter(r_periapsis: f64, v_inf: f64, mu_body: f64) -> f64 {
    let e = hyperbolic_eccentricity(r_periapsis, v_inf, mu_body);
    let a = mu_body / (v_inf * v_inf); // magnitude (positive here)
    a * (e * e - 1.0).sqrt()
}

/// Periapsis radius from B-plane parameter and v-infinity.
///
/// Given b and v_inf, find r_p by solving:
/// b² = r_p² + 2 * r_p * mu / v_inf²
///
/// This is a quadratic in r_p.
pub fn periapsis_from_b_parameter(b: f64, v_inf: f64, mu_body: f64) -> f64 {
    let a = mu_body / (v_inf * v_inf);
    // r_p² + 2*a*r_p - b² = 0
    // r_p = -a + sqrt(a² + b²)
    -a + (a * a + b * b).sqrt()
}

/// Compute a full unpowered gravity assist.
///
/// # Arguments
/// * `v_body` — heliocentric velocity of the flyby body (km/s)
/// * `v_sc_in` — heliocentric velocity of spacecraft before flyby (km/s)
/// * `r_periapsis` — closest approach distance from body center (km)
/// * `mu_body` — gravitational parameter of flyby body (km³/s²)
/// * `rotation_plane_normal` — unit normal defining the flyby plane
///   (the v-infinity vector is rotated by the turn angle in this plane)
///
/// # Returns
/// Complete flyby result including outgoing heliocentric velocity.
pub fn unpowered_flyby(
    v_body: &Vector3<f64>,
    v_sc_in: &Vector3<f64>,
    r_periapsis: f64,
    mu_body: f64,
    rotation_plane_normal: &Vector3<f64>,
) -> FlybyResult {
    // Incoming v-infinity (body-relative)
    let v_inf_vec = v_sc_in - v_body;
    let v_inf = v_inf_vec.norm();

    // Hyperbolic parameters
    let e = hyperbolic_eccentricity(r_periapsis, v_inf, mu_body);
    let delta = turn_angle_from_eccentricity(e);
    let a_hyp = -mu_body / (v_inf * v_inf);
    let b = b_plane_parameter(r_periapsis, v_inf, mu_body);

    // Rotate v-infinity vector by turn angle around the given normal
    let v_inf_unit = v_inf_vec.normalize();
    let n = rotation_plane_normal.normalize();

    // Rodrigues' rotation formula
    let v_inf_out_unit = v_inf_unit * delta.cos()
        + n.cross(&v_inf_unit) * delta.sin()
        + n * n.dot(&v_inf_unit) * (1.0 - delta.cos());

    let v_inf_out_vec = v_inf_out_unit * v_inf;

    // Heliocentric velocity after flyby
    let v_helio_out = v_body + v_inf_out_vec;

    // Delta-v equivalent
    let delta_v_equivalent = (v_helio_out - v_sc_in).norm();

    FlybyResult {
        turn_angle: delta,
        r_periapsis,
        a_hyp,
        e_hyp: e,
        b_parameter: b,
        v_inf_in: v_inf,
        v_inf_out: v_inf,
        v_helio_out,
        delta_v_equivalent,
    }
}

/// Maximum turn angle achievable at a given v-infinity and minimum periapsis.
///
/// This determines the "steering capability" of a flyby body.
pub fn max_turn_angle(v_inf: f64, r_periapsis_min: f64, mu_body: f64) -> f64 {
    let e = hyperbolic_eccentricity(r_periapsis_min, v_inf, mu_body);
    turn_angle_from_eccentricity(e)
}

/// Delta-v equivalent of an unpowered flyby (maximum, for a given v_inf).
///
/// Δv = 2 * v_inf * sin(δ/2)
pub fn max_delta_v_flyby(v_inf: f64, r_periapsis_min: f64, mu_body: f64) -> f64 {
    let delta = max_turn_angle(v_inf, r_periapsis_min, mu_body);
    2.0 * v_inf * (delta / 2.0).sin()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::*;

    #[test]
    fn test_turn_angle() {
        // e = 2.0 → δ = 2*arcsin(0.5) = 60° = π/3
        let delta = turn_angle_from_eccentricity(2.0);
        assert!(
            (delta - std::f64::consts::PI / 3.0).abs() < 1e-10,
            "Turn angle for e=2 should be 60°, got {} deg",
            delta.to_degrees()
        );
    }

    #[test]
    fn test_hyperbolic_eccentricity() {
        // For a flyby of Jupiter at 5 R_J with v_inf = 10 km/s
        let r_p = 5.0 * R_JUPITER;
        let v_inf = 10.0;
        let e = hyperbolic_eccentricity(r_p, v_inf, MU_JUPITER);
        assert!(e > 1.0, "Must be hyperbolic");
        // e = 1 + r_p * v_inf^2 / mu = 1 + 349555 * 100 / 1.267e8 ≈ 1.276
        assert!((e - 1.276).abs() < 0.01, "Expected e ≈ 1.276, got {e}");
    }

    #[test]
    fn test_b_plane_parameter() {
        let r_p = 5.0 * R_JUPITER;
        let v_inf = 10.0;
        let b = b_plane_parameter(r_p, v_inf, MU_JUPITER);
        // b should be positive and larger than r_p
        assert!(b > r_p, "B-parameter should exceed periapsis radius");
    }

    #[test]
    fn test_periapsis_roundtrip() {
        let r_p_orig = 200_000.0; // 200,000 km
        let v_inf = 8.0;
        let b = b_plane_parameter(r_p_orig, v_inf, MU_JUPITER);
        let r_p_calc = periapsis_from_b_parameter(b, v_inf, MU_JUPITER);
        assert!(
            (r_p_calc - r_p_orig).abs() < 1e-6,
            "Roundtrip failed: {r_p_calc} vs {r_p_orig}"
        );
    }

    #[test]
    fn test_unpowered_flyby_conserves_v_inf() {
        // v-infinity magnitude must be conserved in unpowered flyby
        let v_body = Vector3::new(13.07, 0.0, 0.0); // Jupiter helio velocity ~13 km/s
        let v_sc_in = Vector3::new(8.0, 6.0, 0.0);
        let r_p = 5.0 * R_JUPITER;
        let normal = Vector3::new(0.0, 0.0, 1.0);

        let result = unpowered_flyby(&v_body, &v_sc_in, r_p, MU_JUPITER, &normal);
        assert!(
            (result.v_inf_in - result.v_inf_out).abs() < 1e-10,
            "v_inf must be conserved"
        );
    }

    #[test]
    fn test_flyby_changes_heliocentric_velocity() {
        let v_body = Vector3::new(13.07, 0.0, 0.0);
        let v_sc_in = Vector3::new(8.0, 6.0, 0.0);
        let r_p = 5.0 * R_JUPITER;
        let normal = Vector3::new(0.0, 0.0, 1.0);

        let result = unpowered_flyby(&v_body, &v_sc_in, r_p, MU_JUPITER, &normal);
        let v_in_mag = v_sc_in.norm();
        let v_out_mag = result.v_helio_out.norm();
        // Heliocentric speed should change (that's the point of a flyby)
        assert!(
            (v_out_mag - v_in_mag).abs() > 0.01,
            "Flyby should change heliocentric speed"
        );
        assert!(
            result.delta_v_equivalent > 0.0,
            "Should have non-zero delta-v equivalent"
        );
    }

    #[test]
    fn test_max_delta_v_jupiter_flyby() {
        // Voyager-like Jupiter flyby: v_inf ≈ 10 km/s, r_p = 5 R_J
        let dv = max_delta_v_flyby(10.0, 5.0 * R_JUPITER, MU_JUPITER);
        // Should be significant, several km/s
        assert!(
            dv > 3.0 && dv < 20.0,
            "Jupiter flyby dv should be 3-20 km/s, got {dv}"
        );
    }
}
