//! Falcon 9 First Stage Landing Analysis
//!
//! Demonstrates propulsion and trajectory calculations for
//! SpaceX-style booster recovery.

use rocketlab::engines::EngineDatabase;
use rocketlab::constants::*;

fn main() {
    println!("═══════════════════════════════════════════════════════════");
    println!("  🚀 FALCON 9 FIRST STAGE RECOVERY ANALYSIS");
    println!("═══════════════════════════════════════════════════════════\n");

    // Merlin 1D engine specs
    let db = EngineDatabase::new();
    let merlin = db.get_engine("Merlin 1D+").or_else(|| db.get_engine("Merlin 1D")).unwrap();
    
    let thrust_sl = merlin.thrust_sl.unwrap_or(merlin.thrust_vac * 0.9);
    let isp_sl_val = merlin.isp_sl.unwrap_or(merlin.isp_vac * 0.85);
    
    println!("Engine: {}", merlin.name);
    println!("  Thrust (SL):  {:.0} kN", thrust_sl);
    println!("  Thrust (vac): {:.0} kN", merlin.thrust_vac);
    println!("  Isp (SL):     {:.0} s", isp_sl_val);
    println!("  Isp (vac):    {:.0} s", merlin.isp_vac);
    println!("  Propellant:   {:?}", merlin.propellant);

    // Stage parameters
    let m_dry: f64 = 22_200.0; // kg (first stage dry mass)
    let m_prop_landing: f64 = 5_000.0; // kg (reserved for landing)
    let m_total_landing = m_dry + m_prop_landing;
    
    println!("\nFirst Stage Landing Parameters:");
    println!("  Dry mass:           {:.0} kg", m_dry);
    println!("  Landing fuel:       {:.0} kg", m_prop_landing);
    println!("  Total mass:         {:.0} kg", m_total_landing);
    
    // Tsiolkovsky rocket equation for landing ΔV budget
    let g0 = G0;
    let isp_sl = isp_sl_val;
    let ve = isp_sl * g0; // Effective exhaust velocity (m/s)
    let dv_available = ve * (m_total_landing / m_dry).ln();
    
    println!("\nLanding ΔV Budget:");
    println!("  Exhaust velocity:   {:.1} m/s", ve);
    println!("  Available ΔV:       {:.1} m/s", dv_available);
    
    // Landing burn phases
    println!("\nLanding Burn Phases:");
    
    // Boostback burn (3 engines, ~60s)
    let n_engines_bb = 3;
    let thrust_bb = thrust_sl * n_engines_bb as f64 * 1000.0; // N
    let mdot_bb = thrust_bb / ve;
    let burn_time_bb = 60.0; // s
    let dv_bb = ve * ((m_total_landing) / (m_total_landing - mdot_bb * burn_time_bb)).ln();
    println!("  Boostback:  {:.0} m/s ΔV ({} engines, {:.0}s)", dv_bb, n_engines_bb, burn_time_bb);
    
    // Entry burn (3 engines, ~20s)
    let burn_time_entry = 20.0;
    let m_after_bb = m_total_landing - mdot_bb * burn_time_bb;
    let dv_entry = ve * (m_after_bb / (m_after_bb - mdot_bb * burn_time_entry)).ln();
    println!("  Entry burn: {:.0} m/s ΔV ({} engines, {:.0}s)", dv_entry, n_engines_bb, burn_time_entry);
    
    // Landing burn (1 engine, ~30s)
    let thrust_lb = thrust_sl * 1000.0; // N, 1 engine
    let mdot_lb = thrust_lb / ve;
    let burn_time_lb = 30.0;
    let m_before_lb = m_after_bb - mdot_bb * burn_time_entry;
    let dv_lb = ve * (m_before_lb / (m_before_lb - mdot_lb * burn_time_lb)).ln();
    println!("  Landing:    {:.0} m/s ΔV (1 engine, {:.0}s)", dv_lb, burn_time_lb);
    
    println!("\n  Total ΔV used: ~{:.0} m/s", dv_bb + dv_entry + dv_lb);

    println!("\n═══════════════════════════════════════════════════════════");
    println!("  \"The Falcon has landed.\"");
    println!("═══════════════════════════════════════════════════════════");
}
