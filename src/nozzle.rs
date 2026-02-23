//! # Rocket Nozzle Design
//!
//! Design and analysis of rocket nozzles: De Laval (converging-diverging),
//! conical, bell (Rao optimum), and aerospike/plug nozzles.
//!
//! References:
//! - Sutton & Biblarz "Rocket Propulsion Elements" Ch. 3-4
//! - G.V.R. Rao "Exhaust Nozzle Contour for Optimum Thrust" (1958)
//! - NASA SP-8120 "Liquid Rocket Engine Nozzles"

use std::f64::consts::PI;

/// Nozzle geometry and performance parameters.
#[derive(Debug, Clone)]
pub struct NozzleDesign {
    /// Throat radius [m]
    pub throat_radius: f64,
    /// Exit radius [m]
    pub exit_radius: f64,
    /// Expansion ratio (Ae/At)
    pub expansion_ratio: f64,
    /// Nozzle half-angle at exit [rad] (conical) or initial/final angles (bell)
    pub half_angle: f64,
    /// Nozzle length [m]
    pub length: f64,
    /// Thrust coefficient (CF)
    pub thrust_coefficient: f64,
    /// Divergence loss factor (lambda)
    pub divergence_factor: f64,
    /// Nozzle type
    pub nozzle_type: NozzleType,
}

/// Types of rocket nozzles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NozzleType {
    Conical,
    Bell,      // Rao optimum / parabolic bell
    Aerospike, // Plug/aerospike (truncated)
}

impl std::fmt::Display for NozzleType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NozzleType::Conical => write!(f, "Conical"),
            NozzleType::Bell => write!(f, "Bell (Rao)"),
            NozzleType::Aerospike => write!(f, "Aerospike"),
        }
    }
}

/// Isentropic flow relations for a calorically perfect gas.
#[derive(Debug, Clone, Copy)]
pub struct IsentropicFlow {
    /// Ratio of specific heats (gamma)
    pub gamma: f64,
}

impl IsentropicFlow {
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }

    /// Temperature ratio T/T0 as a function of Mach number.
    pub fn temperature_ratio(&self, mach: f64) -> f64 {
        1.0 / (1.0 + (self.gamma - 1.0) / 2.0 * mach * mach)
    }

    /// Pressure ratio P/P0 as a function of Mach number.
    pub fn pressure_ratio(&self, mach: f64) -> f64 {
        self.temperature_ratio(mach)
            .powf(self.gamma / (self.gamma - 1.0))
    }

    /// Density ratio rho/rho0 as a function of Mach number.
    pub fn density_ratio(&self, mach: f64) -> f64 {
        self.temperature_ratio(mach).powf(1.0 / (self.gamma - 1.0))
    }

    /// Area ratio A/A* as a function of Mach number (supersonic branch).
    ///
    /// A/A* = (1/M) * [(2/(γ+1)) * (1 + (γ-1)/2 * M²)]^((γ+1)/(2(γ-1)))
    pub fn area_ratio(&self, mach: f64) -> f64 {
        let g = self.gamma;
        let term = 2.0 / (g + 1.0) * (1.0 + (g - 1.0) / 2.0 * mach * mach);
        (1.0 / mach) * term.powf((g + 1.0) / (2.0 * (g - 1.0)))
    }

    /// Solve for Mach number given area ratio (supersonic solution).
    ///
    /// Uses Newton-Raphson iteration.
    pub fn mach_from_area_ratio(&self, area_ratio: f64, supersonic: bool) -> f64 {
        // Initial guess
        let mut m = if supersonic {
            1.0 + (area_ratio - 1.0).sqrt()
        } else {
            0.5
        };

        for _ in 0..100 {
            let ar = self.area_ratio(m);
            let f = ar - area_ratio;

            // Derivative dA*/dM (numerical)
            let dm = 1e-8;
            let dar = (self.area_ratio(m + dm) - ar) / dm;
            if dar.abs() < 1e-20 {
                break;
            }
            let correction = f / dar;
            m -= correction;
            if m < 0.01 {
                m = 0.01;
            }
            if correction.abs() < 1e-10 {
                break;
            }
        }
        m
    }

    /// Thrust coefficient CF for a nozzle exhausting to ambient pressure.
    ///
    /// CF = sqrt(2γ²/(γ-1) * (2/(γ+1))^((γ+1)/(γ-1)) * [1-(Pe/Pc)^((γ-1)/γ)])
    ///      + (Pe/Pc - Pa/Pc) * Ae/At
    ///
    /// # Arguments
    /// * `pe_pc` - Exit pressure / chamber pressure ratio
    /// * `pa_pc` - Ambient pressure / chamber pressure ratio
    /// * `area_ratio` - Nozzle expansion ratio Ae/At
    pub fn thrust_coefficient(&self, pe_pc: f64, pa_pc: f64, area_ratio: f64) -> f64 {
        let g = self.gamma;
        let term1 = 2.0 * g * g / (g - 1.0);
        let term2 = (2.0 / (g + 1.0)).powf((g + 1.0) / (g - 1.0));
        let term3 = 1.0 - pe_pc.powf((g - 1.0) / g);
        let momentum = (term1 * term2 * term3).sqrt();
        let pressure = (pe_pc - pa_pc) * area_ratio;
        momentum + pressure
    }

    /// Characteristic velocity c* [m/s].
    ///
    /// c* = sqrt(γ R T_c) / (γ * sqrt((2/(γ+1))^((γ+1)/(γ-1))))
    ///
    /// # Arguments
    /// * `tc` - Chamber temperature [K]
    /// * `molar_mass` - Propellant molar mass [g/mol]
    pub fn characteristic_velocity(&self, tc: f64, molar_mass: f64) -> f64 {
        let g = self.gamma;
        let r_specific = 8314.46 / molar_mass; // J/(kg·K), molar_mass in g/mol ≡ kg/kmol
        let num = (g * r_specific * tc).sqrt();
        let den = g * (2.0 / (g + 1.0)).powf((g + 1.0) / (2.0 * (g - 1.0)));
        num / den
    }
}

/// Design a conical nozzle.
///
/// The simplest supersonic nozzle. Divergence losses are accounted for
/// by the lambda factor: λ = (1 + cos α) / 2
///
/// # Arguments
/// * `throat_radius` - Throat radius [m]
/// * `expansion_ratio` - Ae/At
/// * `half_angle_deg` - Cone half-angle [degrees] (typically 15°)
/// * `gamma` - Ratio of specific heats
/// * `pe_pc` - Exit/chamber pressure ratio
/// * `pa_pc` - Ambient/chamber pressure ratio
pub fn design_conical(
    throat_radius: f64,
    expansion_ratio: f64,
    half_angle_deg: f64,
    gamma: f64,
    pe_pc: f64,
    pa_pc: f64,
) -> NozzleDesign {
    let alpha = half_angle_deg.to_radians();
    let throat_area = PI * throat_radius * throat_radius;
    let exit_area = throat_area * expansion_ratio;
    let exit_radius = (exit_area / PI).sqrt();

    // Nozzle length
    let length = (exit_radius - throat_radius) / alpha.tan();

    // Divergence loss factor
    let divergence_factor = (1.0 + alpha.cos()) / 2.0;

    let flow = IsentropicFlow::new(gamma);
    let cf_ideal = flow.thrust_coefficient(pe_pc, pa_pc, expansion_ratio);
    let cf = cf_ideal * divergence_factor;

    NozzleDesign {
        throat_radius,
        exit_radius,
        expansion_ratio,
        half_angle: alpha,
        length,
        thrust_coefficient: cf,
        divergence_factor,
        nozzle_type: NozzleType::Conical,
    }
}

/// Design a bell (Rao optimum) nozzle.
///
/// The bell nozzle uses a parabolic contour optimized by G.V.R. Rao (1958).
/// It achieves near-ideal performance with shorter length than conical.
/// Typical bell nozzles are 80% the length of an equivalent 15° conical nozzle.
///
/// # Arguments
/// * `throat_radius` - Throat radius [m]
/// * `expansion_ratio` - Ae/At
/// * `percent_bell` - Fraction of equivalent 15° cone length (typically 0.80)
/// * `gamma` - Ratio of specific heats
/// * `pe_pc` - Exit/chamber pressure ratio
/// * `pa_pc` - Ambient/chamber pressure ratio
pub fn design_bell(
    throat_radius: f64,
    expansion_ratio: f64,
    percent_bell: f64,
    gamma: f64,
    pe_pc: f64,
    pa_pc: f64,
) -> NozzleDesign {
    let throat_area = PI * throat_radius * throat_radius;
    let exit_area = throat_area * expansion_ratio;
    let exit_radius = (exit_area / PI).sqrt();

    // Reference: 15° conical nozzle length
    let conical_length = (exit_radius - throat_radius) / (15.0_f64.to_radians().tan());
    let length = conical_length * percent_bell;

    // Bell nozzle divergence factor is higher than conical
    // Rao optimum approaches λ ≈ 0.985-0.995 for typical designs
    // Approximate from percent_bell and expansion ratio
    let divergence_factor = 0.985 + 0.01 * percent_bell.min(1.0);

    let flow = IsentropicFlow::new(gamma);
    let cf_ideal = flow.thrust_coefficient(pe_pc, pa_pc, expansion_ratio);
    let cf = cf_ideal * divergence_factor;

    // Exit angle for Rao bell (approximate empirical correlation)
    // For 80% bell: θe ≈ 7-8° for ε=20-40
    let exit_angle = (8.0 - 0.05 * (expansion_ratio - 20.0).max(0.0)).max(3.0);

    NozzleDesign {
        throat_radius,
        exit_radius,
        expansion_ratio,
        half_angle: exit_angle.to_radians(),
        length,
        thrust_coefficient: cf,
        divergence_factor,
        nozzle_type: NozzleType::Bell,
    }
}

/// Compute nozzle throat area from required thrust.
///
/// At = F / (CF * Pc)
///
/// # Arguments
/// * `thrust` - Required thrust [N]
/// * `cf` - Thrust coefficient
/// * `chamber_pressure` - Chamber pressure [Pa]
///
/// # Returns
/// Throat area [m²]
pub fn throat_area_from_thrust(thrust: f64, cf: f64, chamber_pressure: f64) -> f64 {
    thrust / (cf * chamber_pressure)
}

/// Compute exit Mach number from expansion ratio.
pub fn exit_mach(expansion_ratio: f64, gamma: f64) -> f64 {
    let flow = IsentropicFlow::new(gamma);
    flow.mach_from_area_ratio(expansion_ratio, true)
}

/// Compute optimal expansion ratio for a given altitude.
///
/// The optimal expansion ratio is when exit pressure equals ambient pressure
/// (Pe = Pa), eliminating the pressure thrust term.
///
/// # Arguments
/// * `chamber_pressure` - Chamber pressure [Pa]
/// * `ambient_pressure` - Ambient pressure [Pa]
/// * `gamma` - Ratio of specific heats
pub fn optimal_expansion_ratio(chamber_pressure: f64, ambient_pressure: f64, gamma: f64) -> f64 {
    let pe_pc = ambient_pressure / chamber_pressure;
    let flow = IsentropicFlow::new(gamma);

    // Find Mach number where P/P0 = pe_pc
    // P/P0 = (1 + (γ-1)/2 * M²)^(-γ/(γ-1))
    let g = gamma;
    let m_exit = ((pe_pc.powf(-(g - 1.0) / g) - 1.0) * 2.0 / (g - 1.0)).sqrt();

    flow.area_ratio(m_exit)
}

/// Compute mass flow rate through a choked nozzle.
///
/// ṁ = Pc * At * γ * sqrt(2/(γ+1))^((γ+1)/(γ-1)) / sqrt(γ * R * Tc)
///
/// Simplified: ṁ = Pc * At / c*
///
/// # Arguments
/// * `chamber_pressure` - Chamber pressure [Pa]
/// * `throat_area` - Throat area [m²]
/// * `c_star` - Characteristic velocity [m/s]
pub fn mass_flow_rate(chamber_pressure: f64, throat_area: f64, c_star: f64) -> f64 {
    chamber_pressure * throat_area / c_star
}

/// Specific impulse from thrust coefficient and characteristic velocity.
///
/// Isp = CF * c* / g0
///
/// # Arguments
/// * `cf` - Thrust coefficient
/// * `c_star` - Characteristic velocity [m/s]
pub fn specific_impulse(cf: f64, c_star: f64) -> f64 {
    cf * c_star / 9.80665
}

/// Generate nozzle contour points for a bell nozzle.
///
/// Approximate Rao parabolic contour using the method of characteristics
/// simplification: upstream circular arc + downstream parabola.
///
/// # Arguments
/// * `design` - Nozzle design parameters
/// * `n_points` - Number of contour points
///
/// # Returns
/// Vec of (x, r) coordinates [m] from throat to exit
pub fn bell_contour(design: &NozzleDesign, n_points: usize) -> Vec<(f64, f64)> {
    let rt = design.throat_radius;
    let re = design.exit_radius;
    let l = design.length;

    // Upstream circular arc (from Sutton & Biblarz)
    // The throat has a circular arc with radius ~ 1.5 * rt (upstream) and 0.382 * rt (downstream)
    let r_curve_dn = 0.382 * rt; // Downstream throat radius of curvature

    // Initial expansion angle (Rao bell)
    let theta_n = if design.expansion_ratio < 10.0 {
        30.0_f64.to_radians()
    } else if design.expansion_ratio < 30.0 {
        25.0_f64.to_radians()
    } else {
        20.0_f64.to_radians()
    };

    let theta_e = design.half_angle; // Exit angle

    // Inflection point (end of circular arc)
    let x_n = r_curve_dn * theta_n.sin();
    let r_n = rt + r_curve_dn * (1.0 - theta_n.cos());

    let mut points = Vec::with_capacity(n_points);

    // Circular arc section (throat to inflection)
    let n_arc = n_points / 4;
    for i in 0..=n_arc {
        let theta = theta_n * (i as f64 / n_arc as f64);
        let x = r_curve_dn * theta.sin();
        let r = rt + r_curve_dn * (1.0 - theta.cos());
        points.push((x, r));
    }

    // Parabolic section (inflection to exit)
    // Parametric parabola matching slope at inflection and exit
    let n_parabola = n_points - n_arc - 1;
    for i in 1..=n_parabola {
        let t = i as f64 / n_parabola as f64;
        // Quadratic Bezier: P = (1-t)²P_n + 2t(1-t)P_q + t²P_e
        // Control point from tangent intersection
        let x_e = l;

        // Tangent intersection point
        let x_q = (re - r_n + x_n * theta_n.tan() - x_e * theta_e.tan())
            / (theta_n.tan() - theta_e.tan());
        let r_q = r_n + (x_q - x_n) * theta_n.tan();

        let x = (1.0 - t).powi(2) * x_n + 2.0 * t * (1.0 - t) * x_q + t * t * x_e;
        let r = (1.0 - t).powi(2) * r_n + 2.0 * t * (1.0 - t) * r_q + t * t * re;
        points.push((x, r));
    }

    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isentropic_relations() {
        let flow = IsentropicFlow::new(1.4);

        // At Mach 1: T/T0 = 2/(γ+1) = 0.8333
        let t_ratio = flow.temperature_ratio(1.0);
        assert!(
            (t_ratio - 0.8333).abs() < 0.001,
            "T/T0 at M=1: {t_ratio:.4}"
        );

        // At Mach 1: P/P0 = (2/(γ+1))^(γ/(γ-1)) = 0.5283
        let p_ratio = flow.pressure_ratio(1.0);
        assert!(
            (p_ratio - 0.5283).abs() < 0.001,
            "P/P0 at M=1: {p_ratio:.4}"
        );

        // At Mach 1: A/A* = 1.0
        let ar = flow.area_ratio(1.0);
        assert!((ar - 1.0).abs() < 0.001, "A/A* at M=1: {ar:.4}");
    }

    #[test]
    fn test_area_ratio_mach_2() {
        let flow = IsentropicFlow::new(1.4);
        // At Mach 2: A/A* = 1.6875 (standard gas dynamics tables)
        let ar = flow.area_ratio(2.0);
        assert!((ar - 1.6875).abs() < 0.001, "A/A* at M=2: {ar:.4}");
    }

    #[test]
    fn test_mach_from_area_ratio() {
        let flow = IsentropicFlow::new(1.4);
        let m = flow.mach_from_area_ratio(1.6875, true);
        assert!((m - 2.0).abs() < 0.01, "M for A/A*=1.6875: {m:.4}");
    }

    #[test]
    fn test_thrust_coefficient_vacuum() {
        let flow = IsentropicFlow::new(1.2);
        // Vacuum (Pa=0), typical LOX/LH2 conditions
        // For γ=1.2, ε=40: CF_vac ≈ 1.88-1.92 (Sutton Table 3-2)
        let pe_pc = flow.pressure_ratio(flow.mach_from_area_ratio(40.0, true));
        let cf = flow.thrust_coefficient(pe_pc, 0.0, 40.0);
        assert!(cf > 1.80 && cf < 2.00, "CF_vac for ε=40, γ=1.2: {cf:.3}");
    }

    #[test]
    fn test_thrust_coefficient_sea_level() {
        let flow = IsentropicFlow::new(1.2);
        // Sea level, Pc=7MPa, Pa=101325Pa
        let pa_pc = 101325.0 / 7e6;
        let epsilon = 12.0;
        let me = flow.mach_from_area_ratio(epsilon, true);
        let pe_pc = flow.pressure_ratio(me);
        let cf = flow.thrust_coefficient(pe_pc, pa_pc, epsilon);
        // Sea level CF typically 1.3-1.6 for reasonable expansion ratios
        assert!(cf > 1.2 && cf < 1.8, "CF_SL for ε=12: {cf:.3}");
    }

    #[test]
    fn test_conical_nozzle_divergence() {
        // 15° half-angle: λ = (1 + cos 15°)/2 = 0.9830
        let design = design_conical(0.1, 20.0, 15.0, 1.2, 0.01, 0.0);
        assert!(
            (design.divergence_factor - 0.9830).abs() < 0.001,
            "λ = {:.4}",
            design.divergence_factor
        );
        assert_eq!(design.nozzle_type, NozzleType::Conical);
    }

    #[test]
    fn test_bell_shorter_than_conical() {
        let conical = design_conical(0.1, 25.0, 15.0, 1.2, 0.01, 0.0);
        let bell = design_bell(0.1, 25.0, 0.8, 1.2, 0.01, 0.0);
        assert!(
            bell.length < conical.length,
            "Bell length {:.3}m should be < conical {:.3}m",
            bell.length,
            conical.length
        );
        // Bell should be ~80% of 15° conical
        let ratio = bell.length / conical.length;
        assert!((ratio - 0.8).abs() < 0.01, "Length ratio: {ratio:.3}");
    }

    #[test]
    fn test_bell_better_divergence() {
        let conical = design_conical(0.1, 25.0, 15.0, 1.2, 0.01, 0.0);
        let bell = design_bell(0.1, 25.0, 0.8, 1.2, 0.01, 0.0);
        assert!(
            bell.divergence_factor > conical.divergence_factor,
            "Bell λ={:.4} should be > conical λ={:.4}",
            bell.divergence_factor,
            conical.divergence_factor
        );
    }

    #[test]
    fn test_optimal_expansion_ratio() {
        // Sea level: Pc=7MPa, Pa=101325Pa, γ=1.2
        let eps = optimal_expansion_ratio(7e6, 101325.0, 1.2);
        // For these conditions, optimal ε is typically 8-15
        assert!(eps > 5.0 && eps < 25.0, "Optimal ε = {eps:.1}");
    }

    #[test]
    fn test_throat_area_sizing() {
        // Size throat for 1 MN thrust, CF=1.5, Pc=10 MPa
        let at = throat_area_from_thrust(1e6, 1.5, 10e6);
        let rt = (at / PI).sqrt();
        // At = 1e6 / (1.5 * 1e7) = 0.0667 m²
        assert!((at - 0.0667).abs() < 0.001, "At = {at:.4} m²");
        assert!(rt > 0.1 && rt < 0.2, "rt = {rt:.3} m");
    }

    #[test]
    fn test_characteristic_velocity() {
        let flow = IsentropicFlow::new(1.2);
        // LOX/LH2: Tc ≈ 3500K, M ≈ 12 g/mol → c* ≈ 2300-2400 m/s
        let c_star = flow.characteristic_velocity(3500.0, 12.0);
        assert!(c_star > 2000.0 && c_star < 2800.0, "c* = {c_star:.0} m/s");
    }

    #[test]
    fn test_specific_impulse() {
        // RS-25: CF ≈ 1.9 (vacuum), c* ≈ 2330 m/s → Isp ≈ 452 s
        let isp = specific_impulse(1.9, 2330.0);
        assert!((isp - 451.0).abs() < 15.0, "Isp = {isp:.1} s");
    }

    #[test]
    fn test_bell_contour_generation() {
        let design = design_bell(0.1, 25.0, 0.8, 1.2, 0.01, 0.0);
        let contour = bell_contour(&design, 100);
        assert!(contour.len() > 50);
        // First point near throat
        assert!((contour[0].1 - design.throat_radius).abs() < 0.01);
        // Last point near exit radius
        let last = contour.last().unwrap();
        assert!(
            (last.1 - design.exit_radius).abs() / design.exit_radius < 0.02,
            "Exit radius: {:.4} vs {:.4}",
            last.1,
            design.exit_radius
        );
        // Monotonically increasing radius
        for w in contour.windows(2) {
            assert!(w[1].1 >= w[0].1 - 1e-10, "Radius not monotonic");
        }
    }

    #[test]
    fn test_mass_flow() {
        // Pc=7MPa, At=0.05m², c*=2300 m/s
        let mdot = mass_flow_rate(7e6, 0.05, 2300.0);
        // mdot = 7e6 * 0.05 / 2300 ≈ 152 kg/s
        assert!((mdot - 152.2).abs() < 1.0, "mdot = {mdot:.1} kg/s");
    }
}
