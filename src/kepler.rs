//! Kepler's equation solver using Newton-Raphson iteration.
//!
//! Solves M = E - e·sin(E) for eccentric anomaly E given mean anomaly M
//! and eccentricity e.
//!
//! Reference: Bate, Mueller & White, §2.4

use crate::constants::TWO_PI;

/// Maximum iterations for Newton-Raphson convergence
const MAX_ITER: usize = 50;

/// Convergence tolerance (radians)
const TOLERANCE: f64 = 1e-12;

/// Solve Kepler's equation for eccentric anomaly.
///
/// # Arguments
/// * `mean_anomaly` - Mean anomaly M (radians)
/// * `eccentricity` - Orbital eccentricity e (0 ≤ e < 1 for elliptical)
///
/// # Returns
/// Eccentric anomaly E (radians)
///
/// # Panics
/// Panics if eccentricity is negative or ≥ 1 (hyperbolic not yet supported here).
pub fn solve_kepler(mean_anomaly: f64, eccentricity: f64) -> Result<f64, KeplerError> {
    if eccentricity < 0.0 {
        return Err(KeplerError::InvalidEccentricity(eccentricity));
    }

    if eccentricity >= 1.0 {
        return solve_kepler_hyperbolic(mean_anomaly, eccentricity);
    }

    // Normalize M to [0, 2π)
    let m = mean_anomaly.rem_euclid(TWO_PI);
    let e = eccentricity;

    // Initial guess: E₀ = M + e·sin(M) for low eccentricity
    // For high eccentricity, use E₀ = π
    let mut ea = if e < 0.8 {
        m + e * m.sin()
    } else {
        std::f64::consts::PI
    };

    for _ in 0..MAX_ITER {
        let f = ea - e * ea.sin() - m;
        let fp = 1.0 - e * ea.cos();

        let delta = f / fp;
        ea -= delta;

        if delta.abs() < TOLERANCE {
            return Ok(ea);
        }
    }

    Err(KeplerError::NoConvergence)
}

/// Solve the hyperbolic Kepler's equation: M = e·sinh(H) - H
///
/// # Arguments
/// * `mean_anomaly` - Hyperbolic mean anomaly M
/// * `eccentricity` - Eccentricity e > 1
pub fn solve_kepler_hyperbolic(mean_anomaly: f64, eccentricity: f64) -> Result<f64, KeplerError> {
    if eccentricity <= 1.0 {
        return Err(KeplerError::InvalidEccentricity(eccentricity));
    }

    let m = mean_anomaly;
    let e = eccentricity;

    // Initial guess
    let mut h = m;

    for _ in 0..MAX_ITER {
        let f = e * h.sinh() - h - m;
        let fp = e * h.cosh() - 1.0;

        let delta = f / fp;
        h -= delta;

        if delta.abs() < TOLERANCE {
            return Ok(h);
        }
    }

    Err(KeplerError::NoConvergence)
}

/// Convert eccentric anomaly to true anomaly.
///
/// # Arguments
/// * `eccentric_anomaly` - E (radians)
/// * `eccentricity` - e
///
/// # Returns
/// True anomaly ν (radians)
pub fn eccentric_to_true_anomaly(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let ea = eccentric_anomaly;

    let y = (1.0 + e).sqrt() * (ea / 2.0).sin();
    let x = (1.0 - e).sqrt() * (ea / 2.0).cos();

    2.0 * y.atan2(x)
}

/// Convert true anomaly to eccentric anomaly.
pub fn true_to_eccentric_anomaly(true_anomaly: f64, eccentricity: f64) -> f64 {
    let nu = true_anomaly;
    let e = eccentricity;

    let y = (1.0 - e).sqrt() * (nu / 2.0).sin();
    let x = (1.0 + e).sqrt() * (nu / 2.0).cos();

    2.0 * y.atan2(x)
}

/// Convert eccentric anomaly to mean anomaly.
pub fn eccentric_to_mean_anomaly(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    eccentric_anomaly - eccentricity * eccentric_anomaly.sin()
}

/// Errors from Kepler equation solving
#[derive(Debug, Clone)]
pub enum KeplerError {
    /// Newton-Raphson did not converge
    NoConvergence,
    /// Invalid eccentricity value
    InvalidEccentricity(f64),
}

impl std::fmt::Display for KeplerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KeplerError::NoConvergence => write!(f, "Kepler equation solver did not converge"),
            KeplerError::InvalidEccentricity(e) => write!(f, "Invalid eccentricity: {e}"),
        }
    }
}

impl std::error::Error for KeplerError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular_orbit() {
        // For e=0, E = M
        let m = 1.5;
        let e = solve_kepler(m, 0.0).unwrap();
        assert!((e - m).abs() < 1e-10);
    }

    #[test]
    fn test_low_eccentricity() {
        // ISS-like orbit, e ≈ 0.0007
        let m = 1.0;
        let ecc = 0.0007;
        let ea = solve_kepler(m, ecc).unwrap();

        // Verify: M = E - e*sin(E)
        let m_check = ea - ecc * ea.sin();
        assert!((m_check - m).abs() < 1e-10);
    }

    #[test]
    fn test_high_eccentricity() {
        // Molniya-like orbit, e ≈ 0.74
        let m = 2.0;
        let ecc = 0.74;
        let ea = solve_kepler(m, ecc).unwrap();

        let m_check = ea - ecc * ea.sin();
        assert!((m_check - m).abs() < 1e-10);
    }

    #[test]
    fn test_anomaly_roundtrip() {
        let ecc = 0.3;
        let nu = 1.2;
        let ea = true_to_eccentric_anomaly(nu, ecc);
        let nu2 = eccentric_to_true_anomaly(ea, ecc);
        assert!((nu - nu2).abs() < 1e-10);
    }

    #[test]
    fn test_mean_anomaly_roundtrip() {
        let ecc = 0.5;
        let m_orig = 1.8;
        let ea = solve_kepler(m_orig, ecc).unwrap();
        let m_back = eccentric_to_mean_anomaly(ea, ecc);
        assert!((m_orig - m_back).abs() < 1e-10);
    }

    #[test]
    fn test_hyperbolic() {
        let ecc = 1.5;
        let m = 2.0;
        let h = solve_kepler(m, ecc).unwrap();

        // Verify: M = e*sinh(H) - H
        let m_check = ecc * h.sinh() - h;
        assert!((m_check - m).abs() < 1e-10);
    }
}
