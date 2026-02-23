use clap::{Parser, Subcommand};
use rocketlab::{atmosphere, constants, kepler, orbits, plotting, propulsion, trajectory};

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
        /// Intermediate altitude for bi-elliptic (km). If set, uses bi-elliptic transfer.
        #[arg(long, name = "via")]
        via: Option<f64>,
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
    /// Query US Standard Atmosphere 1976
    Atmosphere {
        /// Altitude in meters
        #[arg(long)]
        altitude: f64,
    },
    /// Simulate a sounding rocket trajectory
    Simulate {
        /// Dry mass (kg)
        #[arg(long, default_value = "100")]
        dry_mass: f64,
        /// Propellant mass (kg)
        #[arg(long, default_value = "400")]
        propellant: f64,
        /// Thrust at sea level (N)
        #[arg(long, default_value = "15000")]
        thrust: f64,
        /// Specific impulse at sea level (s)
        #[arg(long, default_value = "240")]
        isp: f64,
        /// Simulation duration (s)
        #[arg(long, default_value = "300")]
        duration: f64,
        /// Show ASCII plot
        #[arg(long)]
        plot: bool,
    },
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
        Commands::Transfer { from, to, via } => {
            let mu = constants::MU_EARTH;
            let r1 = constants::R_EARTH_EQUATORIAL + from;
            let r2 = constants::R_EARTH_EQUATORIAL + to;

            if let Some(via_alt) = via {
                let r_int = constants::R_EARTH_EQUATORIAL + via_alt;
                let (dv1, dv2, dv3, total, time) = orbits::bi_elliptic_transfer(r1, r2, r_int, mu);

                println!("━━━ Bi-Elliptic Transfer ━━━");
                println!("  From:       {from:.0} km altitude (r = {r1:.1} km)");
                println!("  Via:        {via_alt:.0} km altitude (r = {r_int:.1} km)");
                println!("  To:         {to:.0} km altitude (r = {r2:.1} km)");
                println!("  Δv₁:        {dv1:.4} km/s");
                println!("  Δv₂:        {dv2:.4} km/s");
                println!("  Δv₃:        {dv3:.4} km/s");
                println!("  Total Δv:   {total:.4} km/s");
                println!(
                    "  Transfer:   {:.1} min ({:.2} hours)",
                    time / 60.0,
                    time / 3600.0
                );

                // Compare with Hohmann
                let (_, _, dv_hoh, time_hoh) = orbits::hohmann_transfer(r1, r2, mu);
                println!();
                println!("  ── Comparison with Hohmann ──");
                println!("  Hohmann Δv: {dv_hoh:.4} km/s");
                println!(
                    "  Savings:    {:.4} km/s ({:.1}%)",
                    dv_hoh - total,
                    (dv_hoh - total) / dv_hoh * 100.0
                );
                println!(
                    "  Hohmann t:  {:.1} min ({:.2} hours)",
                    time_hoh / 60.0,
                    time_hoh / 3600.0
                );
            } else {
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
        Commands::Atmosphere { altitude } => match atmosphere::us_standard_atmosphere(altitude) {
            Some(atm) => {
                println!("━━━ US Standard Atmosphere 1976 ━━━");
                println!(
                    "  Altitude:      {:.1} m ({:.2} km)",
                    altitude,
                    altitude / 1000.0
                );
                println!(
                    "  Temperature:   {:.2} K ({:.2} °C)",
                    atm.temperature,
                    atm.temperature - 273.15
                );
                println!(
                    "  Pressure:      {:.4} Pa ({:.4} atm)",
                    atm.pressure,
                    atm.pressure / 101_325.0
                );
                println!("  Density:       {:.6} kg/m³", atm.density);
                println!("  Speed of sound:{:.2} m/s", atm.speed_of_sound);
                println!("  Viscosity:     {:.4e} Pa·s", atm.viscosity);
            }
            None => {
                eprintln!("Altitude out of range (0–86,000 m)");
                std::process::exit(1);
            }
        },
        Commands::Simulate {
            dry_mass,
            propellant,
            thrust,
            isp,
            duration,
            plot,
        } => {
            let vehicle = trajectory::Vehicle {
                dry_mass,
                propellant_mass: propellant,
                thrust_sl: thrust,
                thrust_vac: thrust * 1.1,
                isp_sl: isp,
                isp_vac: isp * 1.12,
                cd: 0.3,
                reference_area: 0.1,
                pitch_program: None,
            };

            let config = trajectory::SimConfig {
                dt: 0.1,
                max_time: duration,
                record_interval: 50,
                ..Default::default()
            };

            let result = trajectory::simulate(&vehicle, &config);

            println!("━━━ Trajectory Simulation ━━━");
            println!(
                "  Vehicle:     {:.0} kg dry + {:.0} kg prop",
                dry_mass, propellant
            );
            println!("  Thrust:      {:.0} N (SL), Isp: {:.0} s", thrust, isp);
            println!("  Duration:    {:.0} s", duration);
            println!();

            if let Some(bt) = result.burnout_time {
                println!("  Burnout:     {bt:.1} s");
            }

            if let Some((t, alt, q)) = result.max_q {
                println!("  Max-Q:       {q:.0} Pa at t={t:.1}s, alt={:.0}m", alt);
            }

            let fs = &result.final_state;
            println!();
            println!("  ── Final State ──");
            println!("  Altitude:    {:.0} m ({:.2} km)", fs.y, fs.y / 1000.0);
            println!("  Downrange:   {:.0} m ({:.2} km)", fs.x, fs.x / 1000.0);
            println!(
                "  Speed:       {:.1} m/s ({:.3} km/s)",
                fs.speed(),
                fs.speed() / 1000.0
            );
            println!("  FPA:         {:.2}°", fs.flight_path_angle().to_degrees());
            println!("  Mass:        {:.1} kg", fs.mass);

            if plot {
                let alt_data: Vec<(f64, f64)> = result
                    .states
                    .iter()
                    .map(|s| (s.time, s.y / 1000.0))
                    .collect();
                let vel_data: Vec<(f64, f64)> =
                    result.states.iter().map(|s| (s.time, s.speed())).collect();

                println!();
                let alt_plot = plotting::AsciiPlot::new("Altitude vs Time", "Time (s)", "Alt (km)")
                    .plot(&alt_data);
                println!("{alt_plot}");

                let vel_plot =
                    plotting::AsciiPlot::new("Velocity vs Time", "Time (s)", "Vel (m/s)")
                        .plot(&vel_data);
                println!("{vel_plot}");
            }
        }
    }
}
