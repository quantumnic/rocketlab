//! Falcon 9 First Stage Landing Analysis
//!
//! Demonstrates propulsion and trajectory calculations for
//! SpaceX-style booster recovery using actual engine specifications.

use rocketlab::engines::EngineDatabase;
use rocketlab::constants::*;
use rocketlab::propulsion::*;

fn main() {
    println!("═══════════════════════════════════════════════════════════");
    println!("  🚀 FALCON 9 FIRST STAGE RECOVERY ANALYSIS");
    println!("═══════════════════════════════════════════════════════════\n");

    // Get engine specifications from database
    let db = EngineDatabase::new();
    let merlin = db.get_engine("Merlin 1D+")
        .or_else(|| db.get_engine("Merlin 1D"))
        .expect("Could not find Merlin engine in database");
    
    println!("Engine: {}", merlin.name);
    println!("  Manufacturer: {}", merlin.manufacturer);
    println!("  Thrust (SL):  {:.0} kN", merlin.thrust_sl.unwrap_or(0.0));
    println!("  Thrust (vac): {:.0} kN", merlin.thrust_vac);
    println!("  Isp (SL):     {:.0} s", merlin.isp_sl.unwrap_or(0.0));
    println!("  Isp (vac):    {:.0} s", merlin.isp_vac);
    println!("  Propellant:   {:?}", merlin.propellant);

    // First stage parameters based on public SpaceX data
    let m_dry = 22200.0; // kg (first stage dry mass)
    let m_prop_landing = 25000.0; // kg (fuel reserved for landing burns)
    let m_total_landing = m_dry + m_prop_landing;
    
    println!("\nFirst Stage Landing Parameters:");
    println!("  Dry mass:           {:.0} kg", m_dry);
    println!("  Landing fuel:       {:.0} kg", m_prop_landing);
    println!("  Landing mass:       {:.0} kg", m_total_landing);
    
    // Engine performance
    let thrust_sl = merlin.thrust_sl.unwrap_or(845.0) * 1000.0; // Convert kN to N
    let isp_sl = merlin.isp_sl.unwrap_or(282.0);
    let ve = exhaust_velocity(isp_sl);
    
    println!("\nEngine Performance:");
    println!("  Single engine thrust: {:.0} kN", thrust_sl / 1000.0);
    println!("  Exhaust velocity:     {:.0} m/s", ve);
    println!("  Single engine flow:   {:.1} kg/s", mass_flow_rate(thrust_sl, isp_sl));
    
    // Available delta-V for all landing maneuvers
    let dv_available = tsiolkovsky_delta_v(isp_sl, m_total_landing / m_dry);
    println!("  Total available ΔV:   {:.0} m/s", dv_available);
    
    // Landing burn sequence analysis with realistic values
    println!("\n═══ LANDING BURN SEQUENCE ANALYSIS ═══");
    
    // Boostback burn: 3 engines, moderate throttle
    let n_engines_bb = 3;
    let throttle_bb = 0.7; // 70% throttle
    let thrust_bb = thrust_sl * n_engines_bb as f64 * throttle_bb;
    let burn_time_bb = 30.0; // seconds
    let fuel_used_bb = mass_flow_rate(thrust_bb, isp_sl) * burn_time_bb;
    let dv_bb = tsiolkovsky_delta_v(isp_sl, m_total_landing / (m_total_landing - fuel_used_bb));
    
    println!("\n1. Boostback Burn:");
    println!("   Engines:        {} Merlin 1D+ @ {:.0}% throttle", n_engines_bb, throttle_bb * 100.0);
    println!("   Total thrust:   {:.0} kN", thrust_bb / 1000.0);
    println!("   Burn duration:  {:.0} s", burn_time_bb);
    println!("   Fuel consumed:  {:.0} kg", fuel_used_bb);
    println!("   Delta-V:        {:.0} m/s", dv_bb);
    
    // Entry burn: 3 engines, short duration, moderate throttle
    let m_after_bb = m_total_landing - fuel_used_bb;
    let throttle_entry = 0.5; // 50% throttle
    let thrust_entry = thrust_sl * n_engines_bb as f64 * throttle_entry;
    let burn_time_entry = 15.0; // seconds
    let fuel_used_entry = mass_flow_rate(thrust_entry, isp_sl) * burn_time_entry;
    let dv_entry = tsiolkovsky_delta_v(isp_sl, m_after_bb / (m_after_bb - fuel_used_entry));
    
    println!("\n2. Entry Burn:");
    println!("   Engines:        {} Merlin 1D+ @ {:.0}% throttle", n_engines_bb, throttle_entry * 100.0);
    println!("   Total thrust:   {:.0} kN", thrust_entry / 1000.0);
    println!("   Burn duration:  {:.0} s", burn_time_entry);
    println!("   Fuel consumed:  {:.0} kg", fuel_used_entry);
    println!("   Delta-V:        {:.0} m/s", dv_entry);
    
    // Landing burn: 1 engine, variable throttle for precision
    let n_engines_land = 1;
    let throttle_land = 0.8; // 80% throttle
    let thrust_land = thrust_sl * throttle_land;
    let m_before_landing = m_after_bb - fuel_used_entry;
    let burn_time_land = 25.0; // seconds
    let fuel_used_land = mass_flow_rate(thrust_land, isp_sl) * burn_time_land;
    let dv_land = tsiolkovsky_delta_v(isp_sl, m_before_landing / (m_before_landing - fuel_used_land));
    
    println!("\n3. Landing Burn:");
    println!("   Engines:        {} Merlin 1D+ @ {:.0}% throttle", n_engines_land, throttle_land * 100.0);
    println!("   Total thrust:   {:.0} kN", thrust_land / 1000.0);
    println!("   Burn duration:  {:.0} s", burn_time_land);
    println!("   Fuel consumed:  {:.0} kg", fuel_used_land);
    println!("   Delta-V:        {:.0} m/s", dv_land);
    
    // Summary
    let total_fuel_used = fuel_used_bb + fuel_used_entry + fuel_used_land;
    let total_dv_used = dv_bb + dv_entry + dv_land;
    let final_mass = m_total_landing - total_fuel_used;
    
    println!("\n═══ MISSION SUMMARY ═══");
    println!("Total fuel consumed:    {:.0} kg ({:.1}% of landing fuel)", 
             total_fuel_used, (total_fuel_used / m_prop_landing) * 100.0);
    println!("Total ΔV used:          {:.0} m/s", total_dv_used);
    println!("Remaining fuel:         {:.0} kg", m_prop_landing - total_fuel_used);
    println!("Final landing mass:     {:.0} kg", final_mass);
    
    if total_fuel_used <= m_prop_landing {
        println!("Fuel margin:            {:.1}%", 
                 ((m_prop_landing - total_fuel_used) / m_prop_landing) * 100.0);
    } else {
        println!("⚠️  FUEL SHORTAGE:      {:.0} kg over budget", total_fuel_used - m_prop_landing);
    }
    
    // Performance metrics
    println!("\n═══ PERFORMANCE METRICS ═══");
    let final_twr = thrust_land / (final_mass * G0); // Landing TWR
    let max_deceleration = thrust_bb / m_before_landing; // Max deceleration capability
    
    println!("Landing thrust-to-weight: {:.2}", final_twr);
    println!("Hover capability:       {}", if final_twr > 1.0 { "YES" } else { "NO" });
    println!("Max deceleration:       {:.1} m/s² ({:.1} g)", max_deceleration, max_deceleration / G0);
    
    // Gravitational losses during burns
    let gravity_loss_bb = G0 * burn_time_bb;
    let gravity_loss_entry = G0 * burn_time_entry;
    let gravity_loss_land = G0 * burn_time_land;
    let total_gravity_loss = gravity_loss_bb + gravity_loss_entry + gravity_loss_land;
    
    println!("\nGravitational Losses:");
    println!("  Boostback:      {:.0} m/s", gravity_loss_bb);
    println!("  Entry:          {:.0} m/s", gravity_loss_entry);
    println!("  Landing:        {:.0} m/s", gravity_loss_land);
    println!("  Total:          {:.0} m/s", total_gravity_loss);
    
    // Mission feasibility assessment
    println!("\n═══ MISSION ASSESSMENT ═══");
    if total_fuel_used <= m_prop_landing && final_twr > 0.7 {
        println!("✅ MISSION FEASIBLE");
        println!("   ΔV margin: {:.0} m/s ({:.1}%)", 
                 dv_available - total_dv_used,
                 ((dv_available - total_dv_used) / dv_available) * 100.0);
    } else if total_fuel_used > m_prop_landing {
        println!("❌ INSUFFICIENT FUEL");
        println!("   Need {:.0} kg more fuel", total_fuel_used - m_prop_landing);
    } else {
        println!("⚠️  LOW THRUST-TO-WEIGHT");
        println!("   Landing TWR: {:.2} (recommend > 1.0)", final_twr);
    }
    
    // Real mission comparison
    println!("\n═══ HISTORICAL COMPARISON ═══");
    println!("First successful Falcon 9 landing:");
    println!("  Mission: Orbcomm OG2-2");
    println!("  Date: December 21, 2015");
    println!("  Landing: Land-based at Cape Canaveral");
    println!("  Estimated fuel used: ~{:.0} kg ({:.1}%)",
             total_fuel_used, (total_fuel_used / 133000.0) * 100.0); // 133t first stage prop mass

    println!("\n═══════════════════════════════════════════════════════════");
    println!("  \"The Falcon has landed.\"");
    println!("               — SpaceX, December 21, 2015");
    println!("═══════════════════════════════════════════════════════════");
}