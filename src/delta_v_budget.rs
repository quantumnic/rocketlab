//! Delta-v budget builder for mission planning.
//!
//! Computes modular delta-v budgets from launch to final orbit,
//! including launch losses, plane changes, circularization, and
//! mission-specific maneuvers.
//!
//! References:
//! - Wertz, "Space Mission Engineering: The New SMAD", Ch. 7, 9, 18
//! - Curtis, "Orbital Mechanics for Engineering Students", Ch. 6

use crate::constants::*;
use std::fmt;

/// A single delta-v maneuver in a mission budget.
#[derive(Debug, Clone)]
pub struct Maneuver {
    /// Name of the maneuver
    pub name: String,
    /// Delta-v in km/s
    pub delta_v: f64,
    /// Description/notes
    pub notes: String,
}

/// Complete delta-v budget for a mission.
#[derive(Debug, Clone)]
pub struct DeltaVBudget {
    /// Mission name
    pub mission_name: String,
    /// Ordered list of maneuvers
    pub maneuvers: Vec<Maneuver>,
    /// Margin factor (e.g., 1.05 for 5% margin)
    pub margin_factor: f64,
}

impl DeltaVBudget {
    /// Create a new empty budget.
    pub fn new(name: &str) -> Self {
        Self {
            mission_name: name.to_string(),
            maneuvers: Vec::new(),
            margin_factor: 1.0,
        }
    }

    /// Set a percentage margin on the total delta-v.
    pub fn with_margin(mut self, percent: f64) -> Self {
        self.margin_factor = 1.0 + percent / 100.0;
        self
    }

    /// Add a maneuver to the budget.
    pub fn add(&mut self, name: &str, delta_v: f64, notes: &str) -> &mut Self {
        self.maneuvers.push(Maneuver {
            name: name.to_string(),
            delta_v,
            notes: notes.to_string(),
        });
        self
    }

    /// Total delta-v without margin (km/s).
    pub fn total(&self) -> f64 {
        self.maneuvers.iter().map(|m| m.delta_v).sum()
    }

    /// Total delta-v with margin (km/s).
    pub fn total_with_margin(&self) -> f64 {
        self.total() * self.margin_factor
    }

    /// Compute payload mass fraction given Isp (s) and structural coefficient.
    ///
    /// Uses the rocket equation applied to the total budget.
    ///
    /// # Arguments
    /// * `isp` — specific impulse (seconds)
    /// * `structural_coeff` — structural coefficient ε = m_struct / (m_struct + m_prop)
    /// * `num_stages` — number of stages (equal delta-v split assumed)
    ///
    /// # Returns
    /// Payload mass fraction (m_payload / m_initial)
    pub fn payload_fraction(&self, isp: f64, structural_coeff: f64, num_stages: usize) -> f64 {
        let v_e = isp * G0 * 1e-3; // km/s
        let dv_per_stage = self.total_with_margin() / num_stages as f64;
        let mass_ratio_per_stage = (dv_per_stage / v_e).exp();

        // For each stage: lambda_i = (1 - eps * R) / R
        // where R = mass ratio, eps = structural coefficient
        // Payload fraction = product of (1 - eps * R_i) / R_i for all stages
        // More precisely: for each stage, payload ratio = (1/R - eps) / (1 - eps)
        let mut payload_frac = 1.0;
        for _ in 0..num_stages {
            let r = mass_ratio_per_stage;
            let stage_payload = (1.0 / r - structural_coeff) / (1.0 - structural_coeff);
            if stage_payload <= 0.0 {
                return 0.0; // Impossible mission
            }
            payload_frac *= stage_payload;
        }

        payload_frac
    }
}

impl fmt::Display for DeltaVBudget {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "╔══════════════════════════════════════════════════════╗"
        )?;
        writeln!(f, "║  Delta-V Budget: {:<36} ║", self.mission_name)?;
        writeln!(
            f,
            "╠══════════════════════════════════════════════════════╣"
        )?;
        writeln!(
            f,
            "║  {:<30} {:>8}  {:<10} ║",
            "Maneuver", "ΔV (km/s)", "Notes"
        )?;
        writeln!(
            f,
            "╠══════════════════════════════════════════════════════╣"
        )?;

        for m in &self.maneuvers {
            let notes_short = if m.notes.len() > 10 {
                &m.notes[..10]
            } else {
                &m.notes
            };
            writeln!(
                f,
                "║  {:<30} {:>8.3}  {:<10} ║",
                m.name, m.delta_v, notes_short
            )?;
        }

        writeln!(
            f,
            "╠══════════════════════════════════════════════════════╣"
        )?;
        writeln!(f, "║  {:<30} {:>8.3}  {:<10} ║", "TOTAL", self.total(), "")?;
        if (self.margin_factor - 1.0).abs() > 1e-6 {
            writeln!(
                f,
                "║  {:<30} {:>8.3}  {:<10} ║",
                format!("TOTAL + {:.0}% margin", (self.margin_factor - 1.0) * 100.0),
                self.total_with_margin(),
                ""
            )?;
        }
        writeln!(
            f,
            "╚══════════════════════════════════════════════════════╝"
        )?;
        Ok(())
    }
}

// ── Pre-built mission budgets ────────────────────────────────────────────────

/// LEO to GEO mission budget (Hohmann transfer).
pub fn leo_to_geo() -> DeltaVBudget {
    let r_leo = R_EARTH_EQUATORIAL + 200.0;
    let r_geo = R_EARTH_EQUATORIAL + 35_786.0;

    let v_leo = (MU_EARTH / r_leo).sqrt();
    let v_geo = (MU_EARTH / r_geo).sqrt();

    let a_transfer = (r_leo + r_geo) / 2.0;
    let v_transfer_peri = (MU_EARTH * (2.0 / r_leo - 1.0 / a_transfer)).sqrt();
    let v_transfer_apo = (MU_EARTH * (2.0 / r_geo - 1.0 / a_transfer)).sqrt();

    let dv1 = v_transfer_peri - v_leo;
    let dv2 = v_geo - v_transfer_apo;

    let mut budget = DeltaVBudget::new("LEO to GEO (Hohmann)");
    budget.add("Launch to LEO (200 km)", 9.4, "incl. losses");
    budget.add("GTO injection (perigee burn)", dv1, "Hohmann");
    budget.add("GEO circularization (apogee burn)", dv2, "Hohmann");
    budget.add("Station-keeping (15 yr)", 0.05 * 15.0, "~50 m/s/yr");
    budget.margin_factor = 1.05;
    budget
}

/// LEO to Moon (TLI + LOI) budget.
pub fn leo_to_lunar_orbit() -> DeltaVBudget {
    let mut budget = DeltaVBudget::new("LEO to Lunar Orbit");
    budget.add("Launch to LEO (200 km)", 9.4, "incl. losses");
    budget.add("Trans-Lunar Injection (TLI)", 3.13, "from 200 km LEO");
    budget.add("Lunar Orbit Insertion (LOI)", 0.82, "100 km circ.");
    budget.add("Descent to surface", 1.72, "powered desc.");
    budget.margin_factor = 1.05;
    budget
}

/// LEO to Mars (TMI + MOI + EDL).
pub fn leo_to_mars() -> DeltaVBudget {
    let mut budget = DeltaVBudget::new("LEO to Mars Surface");
    budget.add("Launch to LEO (200 km)", 9.4, "incl. losses");
    budget.add("Trans-Mars Injection (TMI)", 3.6, "Hohmann min.");
    budget.add("Mars Orbit Insertion (MOI)", 0.9, "250 km circ.");
    budget.add("Deorbit burn", 0.05, "entry traj.");
    budget.add("EDL (aero + powered)", 0.6, "terminal desc.");
    budget.margin_factor = 1.10;
    budget
}

/// ISS resupply mission from Kennedy.
pub fn iss_resupply() -> DeltaVBudget {
    let mut budget = DeltaVBudget::new("ISS Resupply (KSC)");
    budget.add("Launch to 200 km parking orbit", 9.4, "28.5° inc");
    budget.add("Plane change (28.5° → 51.6°)", 1.8, "combined w/ circ");
    budget.add("Phasing orbit raise", 0.05, "rendezvous");
    budget.add("Terminal approach", 0.01, "proximity ops");
    budget.add("Deorbit burn", 0.12, "return");
    budget.margin_factor = 1.05;
    budget
}

/// Sun-synchronous orbit from Vandenberg.
pub fn sso_from_vandenberg() -> DeltaVBudget {
    let mut budget = DeltaVBudget::new("SSO from Vandenberg");
    budget.add("Launch to 500 km SSO", 9.6, "~97.4° inc");
    budget.add("Orbit maintenance (5 yr)", 0.025, "drag comp.");
    budget.add("Deorbit", 0.15, "25-yr rule");
    budget.margin_factor = 1.05;
    budget
}

/// Simple plane change delta-v.
///
/// # Arguments
/// * `v` — orbital velocity at the maneuver point (km/s)
/// * `delta_i_deg` — inclination change (degrees)
pub fn plane_change_dv(v: f64, delta_i_deg: f64) -> f64 {
    2.0 * v * (delta_i_deg.to_radians() / 2.0).sin()
}

/// Combined plane change with orbit raise (more efficient).
///
/// # Arguments
/// * `v1` — velocity in initial orbit (km/s)
/// * `v2` — velocity in final orbit (km/s)
/// * `delta_i_deg` — inclination change (degrees)
pub fn combined_plane_change_dv(v1: f64, v2: f64, delta_i_deg: f64) -> f64 {
    let di = delta_i_deg.to_radians();
    (v1 * v1 + v2 * v2 - 2.0 * v1 * v2 * di.cos()).sqrt()
}

/// Deorbit delta-v from a circular orbit to target perigee.
///
/// # Arguments
/// * `r_orbit` — current circular orbit radius (km)
/// * `r_perigee` — target perigee radius (km), typically R_Earth + ~50 km
/// * `mu` — gravitational parameter (km³/s²)
pub fn deorbit_dv(r_orbit: f64, r_perigee: f64, mu: f64) -> f64 {
    let v_circ = (mu / r_orbit).sqrt();
    let a_transfer = (r_orbit + r_perigee) / 2.0;
    let v_transfer = (mu * (2.0 / r_orbit - 1.0 / a_transfer)).sqrt();
    (v_circ - v_transfer).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leo_to_geo_budget() {
        let budget = leo_to_geo();

        // GTO injection ≈ 2.457 km/s (standard Hohmann value)
        let gto_dv = budget.maneuvers[1].delta_v;
        assert!(
            (gto_dv - 2.457).abs() < 0.02,
            "GTO injection {gto_dv} should be ~2.457 km/s"
        );

        // GEO circularization ≈ 1.478 km/s
        let geo_dv = budget.maneuvers[2].delta_v;
        assert!(
            (geo_dv - 1.478).abs() < 0.02,
            "GEO circ {geo_dv} should be ~1.478 km/s"
        );

        // Total without launch should be ~4.7 km/s
        let orbit_dv: f64 = budget.maneuvers[1..].iter().map(|m| m.delta_v).sum();
        assert!(
            orbit_dv > 4.0 && orbit_dv < 5.5,
            "LEO to GEO orbital dv {orbit_dv} should be ~4.7 km/s"
        );
    }

    #[test]
    fn test_plane_change() {
        // 90° plane change at LEO speed (~7.7 km/s): dv = 2*v*sin(45°) = v*sqrt(2) ≈ 10.9 km/s
        let v = 7.7;
        let dv = plane_change_dv(v, 90.0);
        assert!(
            (dv - v * 2.0_f64.sqrt()).abs() < 0.01,
            "90° plane change dv {dv} should be {}",
            v * 2.0_f64.sqrt()
        );
    }

    #[test]
    fn test_combined_plane_change_efficiency() {
        // Combined maneuver should be cheaper than separate
        let v1 = 7.7; // LEO
        let v2 = 3.07; // GEO

        // Separate: Hohmann + plane change at GEO
        let hohmann_dv = (v1 as f64 - v2 as f64).abs(); // simplified
        let separate_dv = hohmann_dv + plane_change_dv(v2, 28.5);

        // Combined
        let combined_dv = combined_plane_change_dv(v1, v2, 28.5);

        // Combined should be less than or comparable to separate
        // (In reality combined at apogee is always cheaper for plane changes)
        assert!(
            combined_dv < separate_dv + 1.0,
            "Combined {combined_dv} should be efficient vs separate {separate_dv}"
        );
    }

    #[test]
    fn test_deorbit_dv() {
        // From ISS orbit (408 km) to entry (50 km perigee)
        let r_orbit = R_EARTH_EQUATORIAL + 408.0;
        let r_perigee = R_EARTH_EQUATORIAL + 50.0;
        let dv = deorbit_dv(r_orbit, r_perigee, MU_EARTH);

        // Deorbit from ISS: ~100-120 m/s
        assert!(
            dv > 0.08 && dv < 0.15,
            "ISS deorbit dv {dv} km/s should be ~0.1 km/s"
        );
    }

    #[test]
    fn test_payload_fraction() {
        let budget = leo_to_geo();
        // LOX/LH2, Isp=450s, eps=0.1, 2 stages
        let frac = budget.payload_fraction(450.0, 0.10, 2);
        assert!(
            frac > 0.0 && frac < 0.15,
            "Payload fraction {frac} should be reasonable for 2-stage LOX/LH2"
        );
    }

    #[test]
    fn test_budget_display() {
        let budget = leo_to_geo();
        let display = format!("{budget}");
        assert!(display.contains("LEO to GEO"));
        assert!(display.contains("TOTAL"));
    }

    #[test]
    fn test_lunar_budget() {
        let budget = leo_to_lunar_orbit();
        let total = budget.total();
        // Total should be ~15+ km/s (launch + TLI + LOI + descent)
        assert!(
            total > 14.0 && total < 17.0,
            "Lunar mission total {total} km/s should be ~15 km/s"
        );
    }

    #[test]
    fn test_mars_budget() {
        let budget = leo_to_mars();
        let total = budget.total();
        // Total should be ~14.5+ km/s
        assert!(
            total > 13.0 && total < 16.0,
            "Mars mission total {total} km/s should be ~14.5 km/s"
        );
    }
}
