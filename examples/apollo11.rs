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
    let r_parking = R_EARTH_EQUATORIAL + parking_alt;
    let v_circ = circular_velocity(r_parking, MU_EARTH);
    let period_min = orbital_period(r_parking, MU_EARTH) / 60.0;
    
    println!("Phase 1: Earth Parking Orbit");
    println!("  Altitude:  {:.0} km", parking_alt);
    println!("  Velocity:  {:.3} km/s", v_circ);
    println!("  Period:    {:.1} min\n", period_min);

    // --- Phase 2: Trans-Lunar Injection ---
    // Use Moon's sphere of influence radius instead of orbital radius
    let r_moon_orbit = 384400.0; // km (average Earth-Moon distance)
    let transfer = hohmann_transfer(r_parking, r_moon_orbit, MU_EARTH);
    
    println!("Phase 2: Trans-Lunar Injection (TLI)");
    println!("  ΔV₁ (TLI burn):    {:.3} km/s", transfer.0);
    println!("  Transfer time:      {:.1} hours", transfer.3 / 3600.0);
    
    // --- Phase 3: Lunar Orbit ---
    let lunar_orbit_alt = 111.0; // km
    let r_moon = 1737.4; // Moon radius in km
    let r_lunar_orbit = r_moon + lunar_orbit_alt;
    let v_lunar = circular_velocity(r_lunar_orbit, MU_MOON);
    let lunar_period = orbital_period(r_lunar_orbit, MU_MOON);
    
    println!("\nPhase 3: Lunar Orbit");
    println!("  Altitude:  {:.0} km above Moon", lunar_orbit_alt);
    println!("  Velocity:  {:.3} km/s", v_lunar);
    println!("  Period:    {:.1} min\n", lunar_period / 60.0);

    // --- Delta-V Budget ---
    let escape_v = escape_velocity(r_parking, MU_EARTH);
    let loi_estimate = 0.820; // Lunar Orbit Insertion estimate
    let tei_estimate = 0.740; // Trans-Earth Injection estimate
    
    println!("Delta-V Budget:");
    println!("  Earth escape velocity: {:.3} km/s", escape_v);
    println!("  TLI ΔV:               {:.3} km/s", transfer.0);
    println!("  LOI ΔV (estimate):     ~{:.3} km/s", loi_estimate);
    println!("  TEI ΔV (estimate):     ~{:.3} km/s", tei_estimate);
    
    let total_dv = transfer.0 + loi_estimate + tei_estimate;
    println!("  ──────────────────────────────────");
    println!("  Total mission ΔV:      ~{:.3} km/s", total_dv);
    
    // --- Propulsion Analysis ---
    println!("\nS-IVB Third Stage (TLI burn):");
    let isp = 421.0; // s (J-2 engine vacuum ISP)
    let ve = isp * G0 / 1000.0; // Convert to km/s
    let mass_ratio = (transfer.0 / ve).exp();
    println!("  Engine:     J-2 (Rocketdyne)");
    println!("  Isp:        {:.0} s", isp);
    println!("  Mass ratio: {:.3}", mass_ratio);
    println!("  Exhaust vel: {:.3} km/s", ve);
    
    // --- Mission Timeline ---
    println!("\nMission Timeline:");
    let tli_duration = transfer.3; // Transfer time in seconds
    println!("  Earth to Moon:     {:.1} days", tli_duration / 86400.0);
    println!("  Lunar orbit stay:  ~{:.1} days", 6.5); // Historical
    println!("  Moon to Earth:     ~{:.1} days", 2.5); // Historical
    println!("  Total duration:    ~{:.1} days", 8.0);
    
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  \"That's one small step for man, one giant leap for mankind.\"");
    println!("               — Neil Armstrong, July 20, 1969");
    println!("═══════════════════════════════════════════════════════════");
}