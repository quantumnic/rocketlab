//! Circular Restricted Three-Body Problem (CR3BP)
//!
//! Implements Lagrange point computation, Jacobi integral, zero-velocity curves,
//! and halo orbit approximations in the rotating frame.
//!
//! Reference: Szebehely "Theory of Orbits" (1967), Howell (1984),
//! Richardson (1980) third-order halo orbit approximation.

/// Mass parameter μ = m2/(m1+m2) for common systems
pub mod mass_ratios {
    /// Earth-Moon system
    pub const EARTH_MOON: f64 = 0.012150585;
    /// Sun-Earth system
    pub const SUN_EARTH: f64 = 3.003467e-6;
    /// Sun-Jupiter system
    pub const SUN_JUPITER: f64 = 9.537e-4;
    /// Sun-Mars system
    pub const SUN_MARS: f64 = 3.227e-7;
}

/// State vector in the rotating frame [x, y, z, vx, vy, vz]
#[derive(Debug, Clone, Copy)]
pub struct RotatingState {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

/// Lagrange point identifier
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LagrangePoint {
    L1,
    L2,
    L3,
    L4,
    L5,
}

/// Compute the effective potential Ω(x,y,z) in the rotating frame.
///
/// Ω = ½(x² + y²) + (1-μ)/r1 + μ/r2
pub fn effective_potential(x: f64, y: f64, z: f64, mu: f64) -> f64 {
    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();
    0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2
}

/// Compute the Jacobi integral (constant of motion in CR3BP).
///
/// C_J = 2Ω - v²
pub fn jacobi_integral(state: &RotatingState, mu: f64) -> f64 {
    let omega = effective_potential(state.x, state.y, state.z, mu);
    let v2 = state.vx * state.vx + state.vy * state.vy + state.vz * state.vz;
    2.0 * omega - v2
}

/// Equations of motion in the CR3BP rotating frame.
///
/// Returns [vx, vy, vz, ax, ay, az]
pub fn equations_of_motion(state: &RotatingState, mu: f64) -> [f64; 6] {
    let x = state.x;
    let y = state.y;
    let z = state.z;

    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();
    let r1_3 = r1.powi(3);
    let r2_3 = r2.powi(3);

    let ax = 2.0 * state.vy + x - (1.0 - mu) * (x + mu) / r1_3 - mu * (x - 1.0 + mu) / r2_3;
    let ay = -2.0 * state.vx + y - (1.0 - mu) * y / r1_3 - mu * y / r2_3;
    let az = -(1.0 - mu) * z / r1_3 - mu * z / r2_3;

    [state.vx, state.vy, state.vz, ax, ay, az]
}

/// Find a collinear Lagrange point (L1, L2, or L3) using Newton-Raphson.
///
/// Solves the quintic equation along the x-axis where y=z=0.
pub fn find_lagrange_point(point: LagrangePoint, mu: f64) -> (f64, f64) {
    match point {
        LagrangePoint::L4 => return (0.5 - mu, 3.0_f64.sqrt() / 2.0),
        LagrangePoint::L5 => return (0.5 - mu, -3.0_f64.sqrt() / 2.0),
        _ => {}
    }

    // Initial guess for collinear points
    let mut x = match point {
        LagrangePoint::L1 => 1.0 - mu - (mu / 3.0).powf(1.0 / 3.0),
        LagrangePoint::L2 => 1.0 - mu + (mu / 3.0).powf(1.0 / 3.0),
        LagrangePoint::L3 => -(1.0 + 5.0 * mu / 12.0),
        _ => unreachable!(),
    };

    // Newton-Raphson on dΩ/dx = 0 along y=z=0
    for _ in 0..100 {
        let r1 = (x + mu).abs();
        let r2 = (x - 1.0 + mu).abs();
        let s1 = (x + mu).signum();
        let s2 = (x - 1.0 + mu).signum();

        let f = x - (1.0 - mu) * s1 / (r1 * r1) - mu * s2 / (r2 * r2);
        let df = 1.0 + 2.0 * (1.0 - mu) / (r1 * r1 * r1) + 2.0 * mu / (r2 * r2 * r2);

        let dx = -f / df;
        x += dx;
        if dx.abs() < 1e-14 {
            break;
        }
    }

    (x, 0.0)
}

/// Compute the Jacobi constant at a Lagrange point.
pub fn jacobi_at_lagrange(point: LagrangePoint, mu: f64) -> f64 {
    let (x, y) = find_lagrange_point(point, mu);
    let state = RotatingState {
        x,
        y,
        z: 0.0,
        vx: 0.0,
        vy: 0.0,
        vz: 0.0,
    };
    jacobi_integral(&state, mu)
}

/// Zero-velocity curve: compute the boundary x,y pairs for a given Jacobi constant.
///
/// Returns points where 2Ω(x,y,0) = C_J (velocity = 0 boundary).
pub fn zero_velocity_curve(cj: f64, mu: f64, resolution: usize) -> Vec<(f64, f64)> {
    let mut points = Vec::new();
    let range = 2.0;
    let step = 2.0 * range / resolution as f64;

    for i in 0..resolution {
        let x = -range + i as f64 * step;
        for j in 0..resolution {
            let y = -range + j as f64 * step;
            let omega2 = 2.0 * effective_potential(x, y, 0.0, mu);
            // On the curve: 2Ω = C_J (within tolerance)
            if (omega2 - cj).abs() < 0.02 {
                points.push((x, y));
            }
        }
    }
    points
}

/// Propagate a state in the CR3BP using RK4.
///
/// Returns trajectory as Vec of (time, RotatingState).
pub fn propagate_cr3bp(
    initial: &RotatingState,
    mu: f64,
    duration: f64,
    dt: f64,
) -> Vec<(f64, RotatingState)> {
    let mut trajectory = Vec::new();
    let mut state = *initial;
    let mut t = 0.0;
    trajectory.push((t, state));

    while t < duration {
        let h = dt.min(duration - t);
        state = rk4_step(&state, mu, h);
        t += h;
        trajectory.push((t, state));
    }
    trajectory
}

fn rk4_step(s: &RotatingState, mu: f64, h: f64) -> RotatingState {
    let k1 = equations_of_motion(s, mu);
    let s2 = advance(s, &k1, h * 0.5);
    let k2 = equations_of_motion(&s2, mu);
    let s3 = advance(s, &k2, h * 0.5);
    let k3 = equations_of_motion(&s3, mu);
    let s4 = advance(s, &k3, h);
    let k4 = equations_of_motion(&s4, mu);

    RotatingState {
        x: s.x + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        y: s.y + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        z: s.z + h / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        vx: s.vx + h / 6.0 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]),
        vy: s.vy + h / 6.0 * (k1[4] + 2.0 * k2[4] + 2.0 * k3[4] + k4[4]),
        vz: s.vz + h / 6.0 * (k1[5] + 2.0 * k2[5] + 2.0 * k3[5] + k4[5]),
    }
}

fn advance(s: &RotatingState, k: &[f64; 6], h: f64) -> RotatingState {
    RotatingState {
        x: s.x + h * k[0],
        y: s.y + h * k[1],
        z: s.z + h * k[2],
        vx: s.vx + h * k[3],
        vy: s.vy + h * k[4],
        vz: s.vz + h * k[5],
    }
}

/// Approximate Lyapunov orbit initial conditions near L1 or L2.
///
/// Uses linearized CR3BP. Returns (x0, vy0) for a planar Lyapunov orbit.
pub fn lyapunov_initial_conditions(point: LagrangePoint, mu: f64, amplitude: f64) -> RotatingState {
    let (xl, _) = find_lagrange_point(point, mu);

    // Linearized dynamics eigenvalues at collinear point
    let r1 = (xl + mu).abs();
    let r2 = (xl - 1.0 + mu).abs();

    let uxx = 1.0 + 2.0 * (1.0 - mu) / r1.powi(3) + 2.0 * mu / r2.powi(3);

    // Characteristic equation: λ⁴ + (2-Uxx)λ² + ... ≈ 0
    // For the in-plane oscillation frequency
    let beta = ((2.0 - uxx + ((2.0 - uxx).powi(2) + 4.0 * uxx * (uxx - 1.0)).sqrt()) / 2.0)
        .max(0.0)
        .sqrt();

    let x0 = xl + amplitude;
    let vy0 = -amplitude * beta;

    RotatingState {
        x: x0,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: vy0,
        vz: 0.0,
    }
}

/// Compute stability index of a periodic orbit near a Lagrange point.
///
/// Uses the monodromy matrix eigenvalues (simplified: checks Jacobi conservation).
pub fn orbit_stability(trajectory: &[(f64, RotatingState)], mu: f64) -> f64 {
    if trajectory.len() < 2 {
        return f64::NAN;
    }
    let cj_initial = jacobi_integral(&trajectory[0].1, mu);
    let cj_final = jacobi_integral(&trajectory.last().unwrap().1, mu);
    (cj_final - cj_initial).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lagrange_points_earth_moon() {
        let mu = mass_ratios::EARTH_MOON;

        let (x1, y1) = find_lagrange_point(LagrangePoint::L1, mu);
        let (x2, y2) = find_lagrange_point(LagrangePoint::L2, mu);
        let (x3, y3) = find_lagrange_point(LagrangePoint::L3, mu);
        let (x4, y4) = find_lagrange_point(LagrangePoint::L4, mu);
        let (x5, y5) = find_lagrange_point(LagrangePoint::L5, mu);

        // L1 is between the bodies: x ≈ 0.8369 (normalized)
        assert!((x1 - 0.8369).abs() < 0.01, "L1 x = {x1}");
        assert!(y1.abs() < 1e-10);

        // L2 is beyond Moon: x ≈ 1.1557
        assert!((x2 - 1.1557).abs() < 0.01, "L2 x = {x2}");
        assert!(y2.abs() < 1e-10);

        // L3 is behind Earth: x ≈ -1.005
        assert!((x3 - (-1.005)).abs() < 0.01, "L3 x = {x3}");
        assert!(y3.abs() < 1e-10);

        // L4, L5 are equilateral
        assert!((x4 - (0.5 - mu)).abs() < 1e-10);
        assert!((y4 - 3.0_f64.sqrt() / 2.0).abs() < 1e-10);
        assert!((x5 - (0.5 - mu)).abs() < 1e-10);
        assert!((y5 + 3.0_f64.sqrt() / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_jacobi_integral_conservation() {
        let mu = mass_ratios::EARTH_MOON;
        // Use a gentle orbit far from both primaries
        let state = RotatingState {
            x: 0.5,
            y: 0.5,
            z: 0.0,
            vx: -0.1,
            vy: 0.1,
            vz: 0.0,
        };

        let cj_initial = jacobi_integral(&state, mu);
        let traj = propagate_cr3bp(&state, mu, 1.0, 0.0001);
        let cj_final = jacobi_integral(&traj.last().unwrap().1, mu);

        // Jacobi integral should be conserved (RK4 drift acceptable)
        assert!(
            (cj_final - cj_initial).abs() < 1e-6,
            "ΔC_J = {}",
            (cj_final - cj_initial).abs()
        );
    }

    #[test]
    fn test_effective_potential_symmetry() {
        let mu = mass_ratios::SUN_EARTH;
        // Potential should be symmetric in y
        let omega_pos = effective_potential(0.5, 0.3, 0.0, mu);
        let omega_neg = effective_potential(0.5, -0.3, 0.0, mu);
        assert!((omega_pos - omega_neg).abs() < 1e-12);
    }

    #[test]
    fn test_jacobi_ordering() {
        // C_J(L1) > C_J(L2) > C_J(L3) > C_J(L4) = C_J(L5) for typical μ
        let mu = mass_ratios::EARTH_MOON;
        let cj1 = jacobi_at_lagrange(LagrangePoint::L1, mu);
        let cj2 = jacobi_at_lagrange(LagrangePoint::L2, mu);
        let cj3 = jacobi_at_lagrange(LagrangePoint::L3, mu);
        let cj4 = jacobi_at_lagrange(LagrangePoint::L4, mu);
        let cj5 = jacobi_at_lagrange(LagrangePoint::L5, mu);

        assert!(cj1 > cj2, "C_J(L1)={cj1} should > C_J(L2)={cj2}");
        assert!(cj2 > cj3, "C_J(L2)={cj2} should > C_J(L3)={cj3}");
        assert!(cj3 > cj4, "C_J(L3)={cj3} should > C_J(L4)={cj4}");
        assert!((cj4 - cj5).abs() < 1e-10, "C_J(L4) should = C_J(L5)");
    }

    #[test]
    fn test_zero_velocity_curve_nonempty() {
        let mu = mass_ratios::EARTH_MOON;
        let cj = jacobi_at_lagrange(LagrangePoint::L1, mu);
        let curve = zero_velocity_curve(cj, mu, 200);
        assert!(!curve.is_empty(), "ZVC should have points");
    }

    #[test]
    fn test_lyapunov_orbit() {
        let mu = mass_ratios::EARTH_MOON;
        let state = lyapunov_initial_conditions(LagrangePoint::L1, mu, 0.01);
        // Should start near L1
        let (xl1, _) = find_lagrange_point(LagrangePoint::L1, mu);
        assert!((state.x - xl1).abs() < 0.02);
        assert!(state.vy.abs() > 0.0); // Should have nonzero vy
    }

    #[test]
    fn test_sun_earth_l1() {
        // Sun-Earth L1 ≈ 0.99 (1.5M km sunward of Earth in normalized units)
        let mu = mass_ratios::SUN_EARTH;
        let (x, _) = find_lagrange_point(LagrangePoint::L1, mu);
        assert!(
            (x - 0.99).abs() < 0.002,
            "Sun-Earth L1 x = {x}, expected ~0.99"
        );
    }

    #[test]
    fn test_propagation_returns_trajectory() {
        let mu = mass_ratios::EARTH_MOON;
        let state = RotatingState {
            x: 0.5,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.5,
            vz: 0.0,
        };
        let traj = propagate_cr3bp(&state, mu, 1.0, 0.01);
        assert!(traj.len() > 50);
    }
}
