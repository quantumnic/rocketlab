//! Apollo 11 Mission Analysis Example
//! 
//! Demonstrates orbital mechanics calculations for a lunar mission:
//! - Trans-Lunar Injection (TLI) 
//! - Lunar Orbit Insertion (LOI)
//! - Mission timeline and delta-V budget

use rocketlab::constants::*;
use rocketlab::orbits::*;

fn main() {
    println!("═══════════════════════════════════════════════════════════");
    println!("  🚀 APOLLO 11 MISSION ANALYSIS — July 16-24, 1969");
    println!("═══════════════════════════════════════════════════════════\n");

    // --- Phase 1: Earth Parking Orbit ---
    let parking_alt = 185.0; // km
    let r_parking = R_EARTH + parking_alt;
    let v_circ = circular_velocity(r_parking, MU_EARTH);
    let period_min = orbital_period(r_parking, MU_EARTH) / 60.0;
    
    println!("Phase 1: Earth Parking Orbit");
    println!("  Altitude:  {:.0} km", parking_alt);
    println!("  Velocity:  {:.3} km/s", v_circ);
    println!("  Period:    {:.1} min\n", period_min);

    // --- Phase 2: Trans-Lunar Injection ---
    let r_moon = 384400.0; // km (average Earth-Moon distance)
    let transfer = hohmann_transfer(r_parking, r_moon, MU_EARTH);
    
    println!("Phase 2: Trans-Lunar Injection (TLI)");
    println!("  ΔV₁ (TLI burn):    {:.3} km/s", transfer.0);
    println!("  Transfer time:      {:.1} hours", 
        std::f64::consts::PI * ((r_parking + r_moon).powi(3) / (8.0 * MU_EARTH)).sqrt() / 3600.0);
    
    // --- Phase 3: Lunar Orbit ---
    let lunar_orbit_alt: f64 = 111.0; // km
    let r_lunar_orbit: f64 = 1737.4 + lunar_orbit_alt; // Moon radius + alt
    let mu_moon: f64 = 4902.8; // km³/s²
    let v_lunar = (mu_moon / r_lunar_orbit).sqrt();
    
    println!("\nPhase 3: Lunar Orbit");
    println!("  Altitude:  {:.0} km above Moon", lunar_orbit_alt);
    println!("  Velocity:  {:.3} km/s", v_lunar);
    println!("  Period:    {:.1} min\n", 2.0 * std::f64::consts::PI * (r_lunar_orbit.powi(3) / mu_moon).sqrt() / 60.0);

    // --- Delta-V Budget ---
    let escape_v = escape_velocity(r_parking, MU_EARTH);
    println!("Delta-V Budget:");
    println!("  Earth escape velocity: {:.3} km/s", escape_v);
    println!("  TLI ΔV:               {:.3} km/s", transfer.0);
    println!("  LOI ΔV (estimate):     ~{:.3} km/s", 0.820);
    println!("  TEI ΔV (estimate):     ~{:.3} km/s", 0.740);
    
    let total_dv = transfer.0 + 0.820 + 0.740;
    println!("  ──────────────────────────────────");
    println!("  Total mission ΔV:      ~{:.3} km/s", total_dv);
    
    // --- Propulsion ---
    println!("\nS-IVB Third Stage (TLI burn):");
    let isp = 421.0; // s (J-2 engine)
    let g0 = 9.80665e-3; // km/s²
    let mass_ratio = (transfer.0 / (isp * g0)).exp();
    println!("  Engine:     J-2 (Rocketdyne)");
    println!("  Isp:        {:.0} s", isp);
    println!("  Mass ratio: {:.3}", mass_ratio);
    
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  \"That's one small step for man, one giant leap for mankind.\"");
    println!("═══════════════════════════════════════════════════════════");
}
