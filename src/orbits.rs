//! Orbital mechanics: elements ↔ state vectors, vis-viva, orbital parameters.
//!
//! Reference: Curtis, §4.5–4.7; Bate, Mueller & White, §2.5

use crate::constants::{PI, TWO_PI};
use nalgebra::Vector3;

/// Classical orbital elements (Keplerian)
#[derive(Debug, Clone, Copy)]
pub struct OrbitalElements {
    /// Semi-major axis (km)
    pub a: f64,
    /// Eccentricity (dimensionless)
    pub e: f64,
    /// Inclination (radians)
    pub i: f64,
    /// Right ascension of ascending node (radians)
    pub raan: f64,
    /// Argument of periapsis (radians)
    pub omega: f64,
    /// True anomaly (radians)
    pub nu: f64,
}

/// Position and velocity state vector
#[derive(Debug, Clone, Copy)]
pub struct StateVector {
    /// Position vector (km)
    pub r: Vector3<f64>,
    /// Velocity vector (km/s)
    pub v: Vector3<f64>,
}

impl OrbitalElements {
    /// Create elements with angles in degrees (convenience).
    pub fn from_degrees(
        a: f64,
        e: f64,
        i_deg: f64,
        raan_deg: f64,
        omega_deg: f64,
        nu_deg: f64,
    ) -> Self {
        let to_rad = PI / 180.0;
        Self {
            a,
            e,
            i: i_deg * to_rad,
            raan: raan_deg * to_rad,
            omega: omega_deg * to_rad,
            nu: nu_deg * to_rad,
        }
    }

    /// Convert orbital elements to state vector (position + velocity).
    ///
    /// Uses the standard rotation from perifocal frame to ECI.
    /// Reference: Curtis, Algorithm 4.5
    pub fn to_state_vector(self, mu: f64) -> StateVector {
        let p = self.a * (1.0 - self.e * self.e); // semi-latus rectum
        let r_mag = p / (1.0 + self.e * self.nu.cos());

        // Position and velocity in perifocal frame (PQW)
        let r_pqw = Vector3::new(r_mag * self.nu.cos(), r_mag * self.nu.sin(), 0.0);

        let v_factor = (mu / p).sqrt();
        let v_pqw = Vector3::new(
            -v_factor * self.nu.sin(),
            v_factor * (self.e + self.nu.cos()),
            0.0,
        );

        // Rotation from perifocal to ECI
        let (cos_o, sin_o) = (self.raan.cos(), self.raan.sin());
        let (cos_w, sin_w) = (self.omega.cos(), self.omega.sin());
        let (cos_i, sin_i) = (self.i.cos(), self.i.sin());

        let rot = nalgebra::Matrix3::new(
            cos_o * cos_w - sin_o * sin_w * cos_i,
            -cos_o * sin_w - sin_o * cos_w * cos_i,
            sin_o * sin_i,
            sin_o * cos_w + cos_o * sin_w * cos_i,
            -sin_o * sin_w + cos_o * cos_w * cos_i,
            -cos_o * sin_i,
            sin_w * sin_i,
            cos_w * sin_i,
            cos_i,
        );

        StateVector {
            r: rot * r_pqw,
            v: rot * v_pqw,
        }
    }

    /// Orbital period (seconds).
    pub fn period(&self, mu: f64) -> f64 {
        TWO_PI * (self.a.powi(3) / mu).sqrt()
    }

    /// Specific orbital energy (km²/s²).
    pub fn energy(&self, mu: f64) -> f64 {
        -mu / (2.0 * self.a)
    }

    /// Velocity at a given radius using vis-viva equation (km/s).
    pub fn velocity_at_radius(&self, r: f64, mu: f64) -> f64 {
        vis_viva(r, self.a, mu)
    }

    /// Periapsis radius (km).
    pub fn periapsis(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Apoapsis radius (km).
    pub fn apoapsis(&self) -> f64 {
        self.a * (1.0 + self.e)
    }
}

impl StateVector {
    /// Convert state vector back to orbital elements.
    ///
    /// Reference: Curtis, Algorithm 4.2
    pub fn to_orbital_elements(self, mu: f64) -> OrbitalElements {
        let r = &self.r;
        let v = &self.v;
        let r_mag = r.norm();
        let v_mag = v.norm();

        // Specific angular momentum
        let h = r.cross(v);
        let h_mag = h.norm();

        // Node vector
        let k = Vector3::new(0.0, 0.0, 1.0);
        let n = k.cross(&h);
        let n_mag = n.norm();

        // Eccentricity vector
        let e_vec = ((v_mag * v_mag - mu / r_mag) * r - r.dot(v) * v) / mu;
        let e = e_vec.norm();

        // Semi-major axis
        let energy = v_mag * v_mag / 2.0 - mu / r_mag;
        let a = -mu / (2.0 * energy);

        // Inclination
        let i = (h[2] / h_mag).acos();

        // RAAN
        let raan = if n_mag > 1e-10 {
            let val = (n[0] / n_mag).acos();
            if n[1] >= 0.0 {
                val
            } else {
                TWO_PI - val
            }
        } else {
            0.0
        };

        // Argument of periapsis
        let omega = if n_mag > 1e-10 && e > 1e-10 {
            let val = (n.dot(&e_vec) / (n_mag * e)).acos();
            if e_vec[2] >= 0.0 {
                val
            } else {
                TWO_PI - val
            }
        } else {
            0.0
        };

        // True anomaly
        let nu = if e > 1e-10 {
            let val = (e_vec.dot(r) / (e * r_mag)).clamp(-1.0, 1.0).acos();
            if r.dot(v) >= 0.0 {
                val
            } else {
                TWO_PI - val
            }
        } else {
            0.0
        };

        OrbitalElements {
            a,
            e,
            i,
            raan,
            omega,
            nu,
        }
    }
}

/// Vis-viva equation: v = sqrt(μ(2/r - 1/a))
///
/// # Arguments
/// * `r` - Current radius (km)
/// * `a` - Semi-major axis (km)
/// * `mu` - Gravitational parameter (km³/s²)
///
/// # Returns
/// Orbital velocity (km/s)
pub fn vis_viva(r: f64, a: f64, mu: f64) -> f64 {
    (mu * (2.0 / r - 1.0 / a)).sqrt()
}

/// Escape velocity at a given radius (km/s).
///
/// v_esc = sqrt(2μ/r)
pub fn escape_velocity(r: f64, mu: f64) -> f64 {
    (2.0 * mu / r).sqrt()
}

/// Circular orbit velocity at radius r (km/s).
pub fn circular_velocity(r: f64, mu: f64) -> f64 {
    (mu / r).sqrt()
}

/// Orbital period for semi-major axis a (seconds).
pub fn orbital_period(a: f64, mu: f64) -> f64 {
    TWO_PI * (a.powi(3) / mu).sqrt()
}

/// Hohmann transfer delta-v values.
///
/// Returns (delta_v1, delta_v2, total_delta_v, transfer_time) in (km/s, km/s, km/s, seconds).
///
/// Reference: Bate, Mueller & White, §6.2
pub fn hohmann_transfer(r1: f64, r2: f64, mu: f64) -> (f64, f64, f64, f64) {
    let a_transfer = (r1 + r2) / 2.0;

    let v_c1 = circular_velocity(r1, mu);
    let v_c2 = circular_velocity(r2, mu);

    let v_t1 = vis_viva(r1, a_transfer, mu);
    let v_t2 = vis_viva(r2, a_transfer, mu);

    let dv1 = (v_t1 - v_c1).abs();
    let dv2 = (v_c2 - v_t2).abs();

    let transfer_time = PI * (a_transfer.powi(3) / mu).sqrt();

    (dv1, dv2, dv1 + dv2, transfer_time)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MU_EARTH, R_EARTH_EQUATORIAL};

    #[test]
    fn test_escape_velocity_earth() {
        // Earth escape velocity from surface ≈ 11.186 km/s
        let v_esc = escape_velocity(R_EARTH_EQUATORIAL, MU_EARTH);
        assert!((v_esc - 11.186).abs() < 0.01, "Escape velocity: {v_esc}");
    }

    #[test]
    fn test_iss_orbital_period() {
        // ISS altitude ≈ 408 km, period ≈ 92.68 min
        let r = R_EARTH_EQUATORIAL + 408.0;
        let period_min = orbital_period(r, MU_EARTH) / 60.0;
        assert!(
            (period_min - 92.68).abs() < 0.5,
            "ISS period: {period_min} min"
        );
    }

    #[test]
    fn test_hohmann_leo_to_geo() {
        // LEO (200 km) to GEO (35786 km)
        let r_leo = R_EARTH_EQUATORIAL + 200.0;
        let r_geo = R_EARTH_EQUATORIAL + 35_786.0;
        let (dv1, dv2, total, _time) = hohmann_transfer(r_leo, r_geo, MU_EARTH);

        // Textbook values: dv1 ≈ 2.46 km/s, dv2 ≈ 1.48 km/s, total ≈ 3.94 km/s
        assert!((dv1 - 2.46).abs() < 0.02, "dv1: {dv1}");
        assert!((dv2 - 1.48).abs() < 0.02, "dv2: {dv2}");
        assert!((total - 3.94).abs() < 0.03, "total: {total}");
    }

    #[test]
    fn test_circular_velocity_leo() {
        let r = R_EARTH_EQUATORIAL + 200.0;
        let v = circular_velocity(r, MU_EARTH);
        // LEO circular velocity ≈ 7.78 km/s
        assert!((v - 7.78).abs() < 0.02, "v_circ: {v}");
    }

    #[test]
    fn test_elements_state_roundtrip() {
        let elements = OrbitalElements::from_degrees(
            6778.0, // ISS-like
            0.001, 51.6, 30.0, 45.0, 60.0,
        );

        let state = elements.to_state_vector(MU_EARTH);
        let elements2 = state.to_orbital_elements(MU_EARTH);

        assert!(
            (elements.a - elements2.a).abs() < 0.01,
            "a: {} vs {}",
            elements.a,
            elements2.a
        );
        assert!(
            (elements.e - elements2.e).abs() < 1e-6,
            "e: {} vs {}",
            elements.e,
            elements2.e
        );
        assert!(
            (elements.i - elements2.i).abs() < 1e-6,
            "i: {} vs {}",
            elements.i,
            elements2.i
        );
    }

    #[test]
    fn test_vis_viva_circular() {
        let r = R_EARTH_EQUATORIAL + 400.0;
        let v = vis_viva(r, r, MU_EARTH); // circular: a = r
        let v_circ = circular_velocity(r, MU_EARTH);
        assert!((v - v_circ).abs() < 1e-10);
    }
}
