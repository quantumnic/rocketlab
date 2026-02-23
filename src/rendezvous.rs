//! # Clohessy-Wiltshire (Hill) Equations for Relative Motion
//!
//! Linearized equations of motion for a deputy spacecraft relative to a chief
//! in a circular reference orbit. Fundamental to spacecraft rendezvous and
//! proximity operations (Bate, Mueller & White Ch. 7; Vallado Ch. 6).
//!
//! The Hill frame (LVLH) has:
//! - x: radial (outward)
//! - y: along-track (in velocity direction)
//! - z: cross-track (normal to orbit plane)

use std::f64::consts::PI;

/// State in the Hill (LVLH) frame: position [m] and velocity [m/s].
#[derive(Debug, Clone, Copy)]
pub struct HillState {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

/// Result of a CW maneuver computation.
#[derive(Debug, Clone)]
pub struct CWManeuver {
    /// Delta-v at departure [m/s]
    pub dv1: [f64; 3],
    /// Delta-v at arrival [m/s]
    pub dv2: [f64; 3],
    /// Total delta-v magnitude [m/s]
    pub total_dv: f64,
    /// Transfer time [s]
    pub transfer_time: f64,
}

impl HillState {
    pub fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self {
            x,
            y,
            z,
            vx,
            vy,
            vz,
        }
    }

    /// Position magnitude [m].
    pub fn range(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Velocity magnitude [m/s].
    pub fn speed(&self) -> f64 {
        (self.vx * self.vx + self.vy * self.vy + self.vz * self.vz).sqrt()
    }
}

/// Propagate a Hill-frame state using the Clohessy-Wiltshire equations.
///
/// # Arguments
/// * `state` - Initial relative state in Hill frame
/// * `n` - Mean motion of the reference orbit [rad/s] (sqrt(mu/a^3))
/// * `dt` - Propagation time [s]
///
/// # Returns
/// The propagated Hill state.
///
/// # Reference
/// Clohessy & Wiltshire (1960), "Terminal Guidance System for Satellite Rendezvous"
pub fn cw_propagate(state: &HillState, n: f64, dt: f64) -> HillState {
    let nt = n * dt;
    let s = nt.sin();
    let c = nt.cos();

    let x0 = state.x;
    let y0 = state.y;
    let z0 = state.z;
    let vx0 = state.vx;
    let vy0 = state.vy;
    let vz0 = state.vz;

    // CW state transition matrix (Vallado Eq. 6-64 through 6-69)
    let x = (4.0 - 3.0 * c) * x0 + s / n * vx0 + 2.0 * (1.0 - c) / n * vy0;
    let y = 6.0 * (s - nt) * x0 + y0 - 2.0 * (1.0 - c) / n * vx0 + (4.0 * s - 3.0 * nt) / n * vy0;
    let z = z0 * c + vz0 * s / n;

    let vx = 3.0 * n * s * x0 + c * vx0 + 2.0 * s * vy0;
    let vy = -6.0 * n * (1.0 - c) * x0 - 2.0 * s * vx0 + (4.0 * c - 3.0) * vy0;
    let vz = -z0 * n * s + vz0 * c;

    HillState {
        x,
        y,
        z,
        vx,
        vy,
        vz,
    }
}

/// Compute the mean motion from orbital altitude and gravitational parameter.
///
/// # Arguments
/// * `altitude` - Orbital altitude [m]
/// * `body_radius` - Central body radius [m]
/// * `mu` - Gravitational parameter [m³/s²]
pub fn mean_motion(altitude: f64, body_radius: f64, mu: f64) -> f64 {
    let a = body_radius + altitude;
    (mu / (a * a * a)).sqrt()
}

/// Compute a two-impulse CW rendezvous maneuver.
///
/// Given an initial relative state and desired final relative state,
/// compute the two delta-v impulses needed for the transfer.
///
/// # Arguments
/// * `initial` - Initial relative state
/// * `target` - Desired final relative state (usually origin for docking)
/// * `n` - Mean motion [rad/s]
/// * `dt` - Transfer time [s]
pub fn cw_two_impulse(initial: &HillState, target: &HillState, n: f64, dt: f64) -> CWManeuver {
    let nt = n * dt;
    let s = nt.sin();
    let c = nt.cos();

    // State transition sub-matrices
    // Phi_rr (position to position)
    let phi_rr = [
        [4.0 - 3.0 * c, 0.0, 0.0],
        [6.0 * (s - nt), 1.0, 0.0],
        [0.0, 0.0, c],
    ];
    // Phi_rv (velocity to position)
    let phi_rv = [
        [s / n, 2.0 * (1.0 - c) / n, 0.0],
        [-2.0 * (1.0 - c) / n, (4.0 * s - 3.0 * nt) / n, 0.0],
        [0.0, 0.0, s / n],
    ];
    // Phi_vr (position to velocity)
    let phi_vr = [
        [3.0 * n * s, 0.0, 0.0],
        [-6.0 * n * (1.0 - c), 0.0, 0.0],
        [0.0, 0.0, -n * s],
    ];
    // Phi_vv (velocity to velocity)
    let phi_vv = [
        [c, 2.0 * s, 0.0],
        [-2.0 * s, 4.0 * c - 3.0, 0.0],
        [0.0, 0.0, c],
    ];

    // Invert Phi_rv (3x3 block-diagonal-ish)
    // det of xy block
    let det_xy = phi_rv[0][0] * phi_rv[1][1] - phi_rv[0][1] * phi_rv[1][0];
    let inv_rv = if det_xy.abs() > 1e-15 && phi_rv[2][2].abs() > 1e-15 {
        [
            [phi_rv[1][1] / det_xy, -phi_rv[0][1] / det_xy, 0.0],
            [-phi_rv[1][0] / det_xy, phi_rv[0][0] / det_xy, 0.0],
            [0.0, 0.0, 1.0 / phi_rv[2][2]],
        ]
    } else {
        // Degenerate case (half-orbit etc.)
        return CWManeuver {
            dv1: [0.0; 3],
            dv2: [0.0; 3],
            total_dv: f64::INFINITY,
            transfer_time: dt,
        };
    };

    let r0 = [initial.x, initial.y, initial.z];
    let rf = [target.x, target.y, target.z];

    // Required transfer velocity: v_transfer = Phi_rv^-1 * (rf - Phi_rr * r0)
    let mut rhs = [0.0; 3];
    for i in 0..3 {
        let mut phi_rr_r0 = 0.0;
        for j in 0..3 {
            phi_rr_r0 += phi_rr[i][j] * r0[j];
        }
        rhs[i] = rf[i] - phi_rr_r0;
    }
    let mut v_transfer_0 = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            v_transfer_0[i] += inv_rv[i][j] * rhs[j];
        }
    }

    // Final velocity from the transfer
    let mut v_transfer_f = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            v_transfer_f[i] += phi_vr[i][j] * r0[j] + phi_vv[i][j] * v_transfer_0[j];
        }
    }

    let v0 = [initial.vx, initial.vy, initial.vz];
    let vf = [target.vx, target.vy, target.vz];

    let dv1 = [
        v_transfer_0[0] - v0[0],
        v_transfer_0[1] - v0[1],
        v_transfer_0[2] - v0[2],
    ];
    let dv2 = [
        vf[0] - v_transfer_f[0],
        vf[1] - v_transfer_f[1],
        vf[2] - v_transfer_f[2],
    ];

    let mag1 = (dv1[0] * dv1[0] + dv1[1] * dv1[1] + dv1[2] * dv1[2]).sqrt();
    let mag2 = (dv2[0] * dv2[0] + dv2[1] * dv2[1] + dv2[2] * dv2[2]).sqrt();

    CWManeuver {
        dv1,
        dv2,
        total_dv: mag1 + mag2,
        transfer_time: dt,
    }
}

/// Compute a V-bar approach trajectory (along the velocity vector).
///
/// Standard ISS approach: deputy approaches from behind along the V-bar.
///
/// # Arguments
/// * `y_start` - Initial along-track offset [m] (positive = behind)
/// * `n` - Mean motion [rad/s]
/// * `dt` - Transfer time [s]
pub fn vbar_approach(y_start: f64, n: f64, dt: f64) -> CWManeuver {
    let initial = HillState::new(0.0, y_start, 0.0, 0.0, 0.0, 0.0);
    let target = HillState::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    cw_two_impulse(&initial, &target, n, dt)
}

/// Compute an R-bar approach trajectory (along the radial vector).
///
/// Deputy approaches from below along the R-bar.
///
/// # Arguments
/// * `x_start` - Initial radial offset [m] (negative = below)
/// * `n` - Mean motion [rad/s]
/// * `dt` - Transfer time [s]
pub fn rbar_approach(x_start: f64, n: f64, dt: f64) -> CWManeuver {
    let initial = HillState::new(x_start, 0.0, 0.0, 0.0, 0.0, 0.0);
    let target = HillState::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    cw_two_impulse(&initial, &target, n, dt)
}

/// Compute the natural drift rate for a relative orbit.
///
/// For a deputy with radial offset x0, the along-track drift rate is:
/// dy/dt = -3/2 * n * x0 (secular term from CW equations)
///
/// This is critical for formation flying and station-keeping.
pub fn along_track_drift_rate(x_offset: f64, n: f64) -> f64 {
    -1.5 * n * x_offset
}

/// Compute a football (2:1 ellipse) relative orbit.
///
/// The natural CW solution produces a 2:1 ellipse in the x-y plane
/// when initialized with the proper velocity conditions.
///
/// # Arguments
/// * `amplitude` - Semi-major axis of the relative ellipse [m]
/// * `n` - Mean motion [rad/s]
///
/// # Returns
/// Initial Hill state that produces a centered 2:1 ellipse.
pub fn football_orbit(amplitude: f64, n: f64) -> HillState {
    // For a centered 2:1 ellipse (no drift):
    // x(t) = A * sin(nt)
    // y(t) = -2A * cos(nt)
    // Requires: vy0 = -2*n*x0 to cancel drift, x0=0, y0=-2A, vx0=A*n
    HillState::new(0.0, -2.0 * amplitude, 0.0, amplitude * n, 0.0, 0.0)
}

/// Generate a CW trajectory (series of states at regular intervals).
///
/// # Arguments
/// * `state` - Initial Hill state
/// * `n` - Mean motion [rad/s]
/// * `total_time` - Total propagation time [s]
/// * `steps` - Number of output steps
pub fn cw_trajectory(
    state: &HillState,
    n: f64,
    total_time: f64,
    steps: usize,
) -> Vec<(f64, HillState)> {
    let dt = total_time / steps as f64;
    (0..=steps)
        .map(|i| {
            let t = dt * i as f64;
            (t, cw_propagate(state, n, t))
        })
        .collect()
}

/// Station-keeping delta-v per orbit for a deputy with radial offset.
///
/// From CW equations, the along-track drift from a radial offset x0 over
/// one orbital period is: Δy = -3π * x0
/// The delta-v to cancel this drift each orbit is approximately:
/// Δv ≈ 3π * n * |x0| / 2 (split into two maneuvers)
pub fn stationkeeping_dv_per_orbit(x_offset: f64, n: f64) -> f64 {
    3.0 * PI * n * x_offset.abs() / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 3.986004418e14; // m³/s²
    const R_EARTH: f64 = 6.371e6; // m
    const ISS_ALT: f64 = 408_000.0; // m

    fn iss_mean_motion() -> f64 {
        mean_motion(ISS_ALT, R_EARTH, MU_EARTH)
    }

    #[test]
    fn test_mean_motion_iss() {
        let n = iss_mean_motion();
        let period = 2.0 * PI / n;
        // ISS period ~92.65 min = 5559 s
        assert!((period - 5559.0).abs() < 30.0, "ISS period: {period:.1} s");
    }

    #[test]
    fn test_cw_propagate_origin_stays_origin() {
        let n = iss_mean_motion();
        let state = HillState::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let propagated = cw_propagate(&state, n, 3600.0);
        assert!(propagated.range() < 1e-10);
    }

    #[test]
    fn test_cw_no_drift_condition() {
        // If vy0 = -2*n*x0 and vx0 = 0, no along-track drift
        let n = iss_mean_motion();
        let x0 = 100.0; // 100m radial offset
        let state = HillState::new(x0, 0.0, 0.0, 0.0, -2.0 * n * x0, 0.0);
        // After one full orbit, y should return to ~0 (periodic, no secular drift)
        let period = 2.0 * PI / n;
        let propagated = cw_propagate(&state, n, period);
        assert!(
            propagated.y.abs() < 0.1,
            "y after one orbit: {:.3} m",
            propagated.y
        );
    }

    #[test]
    fn test_cw_football_orbit_2_to_1() {
        let n = iss_mean_motion();
        let amp = 200.0; // 200m amplitude
        let state = football_orbit(amp, n);
        let period = 2.0 * PI / n;
        // At t=T/4, x should be at max (~amplitude), y~0
        let quarter = cw_propagate(&state, n, period / 4.0);
        assert!(
            (quarter.x.abs() - amp).abs() < 1.0,
            "x at T/4: {:.1}",
            quarter.x
        );
        // After full orbit, should return to initial state
        let full = cw_propagate(&state, n, period);
        assert!(
            (full.x - state.x).abs() < 0.1,
            "x drift: {:.3}",
            full.x - state.x
        );
        assert!(
            (full.y - state.y).abs() < 0.1,
            "y drift: {:.3}",
            full.y - state.y
        );
    }

    #[test]
    fn test_cw_two_impulse_rendezvous() {
        let n = iss_mean_motion();
        let period = 2.0 * PI / n;
        // Deputy 1km behind, rendezvous in half orbit
        let initial = HillState::new(0.0, -1000.0, 0.0, 0.0, 0.0, 0.0);
        let target = HillState::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let maneuver = cw_two_impulse(&initial, &target, n, period / 2.0);
        // Should have finite, reasonable delta-v
        assert!(
            maneuver.total_dv > 0.0 && maneuver.total_dv < 10.0,
            "Rendezvous dv: {:.3} m/s",
            maneuver.total_dv
        );
    }

    #[test]
    fn test_vbar_approach() {
        let n = iss_mean_motion();
        let period = 2.0 * PI / n;
        let maneuver = vbar_approach(500.0, n, period * 0.75);
        assert!(maneuver.total_dv > 0.0 && maneuver.total_dv < 5.0);
    }

    #[test]
    fn test_along_track_drift() {
        let n = iss_mean_motion();
        let drift = along_track_drift_rate(100.0, n);
        // Positive x (above) → negative drift (falls behind)
        assert!(drift < 0.0);
        // Drift magnitude: ~0.17 m/s for 100m offset at ISS
        assert!(drift.abs() < 1.0);
    }

    #[test]
    fn test_cross_track_oscillation() {
        let n = iss_mean_motion();
        let state = HillState::new(0.0, 0.0, 500.0, 0.0, 0.0, 0.0);
        let period = 2.0 * PI / n;
        // Cross-track is simple harmonic: z(t) = z0*cos(nt)
        let half = cw_propagate(&state, n, period / 2.0);
        assert!(
            (half.z + 500.0).abs() < 1.0,
            "z at half period: {:.1}",
            half.z
        );
        let full = cw_propagate(&state, n, period);
        assert!(
            (full.z - 500.0).abs() < 1.0,
            "z at full period: {:.1}",
            full.z
        );
    }

    #[test]
    fn test_trajectory_generation() {
        let n = iss_mean_motion();
        let state = football_orbit(100.0, n);
        let period = 2.0 * PI / n;
        let traj = cw_trajectory(&state, n, period, 100);
        assert_eq!(traj.len(), 101);
        assert!((traj[0].0 - 0.0).abs() < 1e-10);
        assert!((traj[100].0 - period).abs() < 1.0);
    }
}
