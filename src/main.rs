use clap::{Parser, Subcommand};
use rocketlab::{constants, kepler, orbits, propulsion};

#[derive(Parser)]
#[command(name = "rocketlab")]
#[command(about = "Rocket trajectory simulator and aerospace toolkit")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate delta-v using Tsiolkovsky equation
    Deltav {
        /// Specific impulse (seconds)
        #[arg(long)]
        isp: f64,
        /// Mass ratio (m0/mf)
        #[arg(long, name = "mass-ratio")]
        mass_ratio: f64,
    },
    /// Display orbital elements and derived quantities
    Orbit {
        /// Orbital elements: 'a,e,i,raan,omega,nu' (km, deg)
        #[arg(long)]
        elements: String,
    },
    /// Compute Hohmann transfer between circular orbits
    Transfer {
        /// Initial orbit altitude (km)
        #[arg(long)]
        from: f64,
        /// Target orbit altitude (km)
        #[arg(long)]
        to: f64,
    },
    /// Solve Kepler's equation
    Kepler {
        /// Mean anomaly (degrees)
        #[arg(long, name = "mean-anomaly")]
        mean_anomaly: f64,
        /// Eccentricity
        #[arg(long)]
        eccentricity: f64,
    },
    /// List propellant combinations
    Propellants,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Deltav { isp, mass_ratio } => {
            let dv = propulsion::tsiolkovsky_delta_v(isp, mass_ratio);
            let v_e = propulsion::exhaust_velocity(isp);
            println!("━━━ Delta-V Calculation ━━━");
            println!("  Isp:             {isp:.1} s");
            println!("  Mass ratio:      {mass_ratio:.3}");
            println!("  Exhaust vel:     {v_e:.1} m/s");
            println!("  Delta-v:         {dv:.1} m/s ({:.3} km/s)", dv / 1000.0);
        }
        Commands::Orbit { elements } => {
            let parts: Vec<f64> = elements
                .split(',')
                .map(|s| s.trim().parse::<f64>().expect("Invalid number"))
                .collect();
            if parts.len() != 6 {
                eprintln!("Expected 6 elements: a,e,i,raan,omega,nu");
                std::process::exit(1);
            }
            let oe = orbits::OrbitalElements::from_degrees(
                parts[0], parts[1], parts[2], parts[3], parts[4], parts[5],
            );
            let mu = constants::MU_EARTH;
            let state = oe.to_state_vector(mu);
            let period_min = oe.period(mu) / 60.0;
            let energy = oe.energy(mu);

            println!("━━━ Orbital Elements ━━━");
            println!("  Semi-major axis: {:.3} km", oe.a);
            println!("  Eccentricity:    {:.6}", oe.e);
            println!("  Inclination:     {:.4}°", oe.i.to_degrees());
            println!(
                "  Periapsis:       {:.3} km (alt {:.3} km)",
                oe.periapsis(),
                oe.periapsis() - constants::R_EARTH_EQUATORIAL
            );
            println!(
                "  Apoapsis:        {:.3} km (alt {:.3} km)",
                oe.apoapsis(),
                oe.apoapsis() - constants::R_EARTH_EQUATORIAL
            );
            println!("  Period:          {period_min:.2} min");
            println!("  Energy:          {energy:.4} km²/s²");
            println!();
            println!("━━━ State Vector (ECI) ━━━");
            println!(
                "  r = [{:.4}, {:.4}, {:.4}] km",
                state.r[0], state.r[1], state.r[2]
            );
            println!(
                "  v = [{:.6}, {:.6}, {:.6}] km/s",
                state.v[0], state.v[1], state.v[2]
            );
            println!("  |r| = {:.4} km", state.r.norm());
            println!("  |v| = {:.6} km/s", state.v.norm());
        }
        Commands::Transfer { from, to } => {
            let mu = constants::MU_EARTH;
            let r1 = constants::R_EARTH_EQUATORIAL + from;
            let r2 = constants::R_EARTH_EQUATORIAL + to;
            let (dv1, dv2, total, time) = orbits::hohmann_transfer(r1, r2, mu);

            println!("━━━ Hohmann Transfer ━━━");
            println!("  From:       {from:.0} km altitude (r = {r1:.1} km)");
            println!("  To:         {to:.0} km altitude (r = {r2:.1} km)");
            println!("  Δv₁:        {dv1:.4} km/s");
            println!("  Δv₂:        {dv2:.4} km/s");
            println!("  Total Δv:   {total:.4} km/s");
            println!(
                "  Transfer:   {:.1} min ({:.2} hours)",
                time / 60.0,
                time / 3600.0
            );
        }
        Commands::Kepler {
            mean_anomaly,
            eccentricity,
        } => {
            let m_rad = mean_anomaly.to_radians();
            match kepler::solve_kepler(m_rad, eccentricity) {
                Ok(ea) => {
                    let nu = kepler::eccentric_to_true_anomaly(ea, eccentricity);
                    println!("━━━ Kepler's Equation ━━━");
                    println!("  Mean anomaly:      {mean_anomaly:.6}° ({m_rad:.6} rad)");
                    println!("  Eccentricity:      {eccentricity:.6}");
                    println!("  Eccentric anomaly: {:.6}° ({ea:.6} rad)", ea.to_degrees());
                    println!("  True anomaly:      {:.6}° ({nu:.6} rad)", nu.to_degrees());
                }
                Err(e) => {
                    eprintln!("Error: {e}");
                    std::process::exit(1);
                }
            }
        }
        Commands::Propellants => {
            let db = propulsion::propellant_database();
            println!("━━━ Propellant Combinations ━━━");
            println!(
                "{:<15} {:<8} {:<8} {:<8} {:<10}",
                "Name", "Isp(SL)", "Isp(Vac)", "O/F", "Fuel"
            );
            println!("{}", "─".repeat(55));
            for p in &db {
                println!(
                    "{:<15} {:<8.0} {:<8.0} {:<8.1} {}",
                    p.name, p.isp_sea_level, p.isp_vacuum, p.mixture_ratio, p.fuel
                );
            }
        }
    }
}
