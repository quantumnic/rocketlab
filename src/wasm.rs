//! WebAssembly bindings for rocketlab
//!
//! Exposes core simulation functions to JavaScript via wasm-bindgen.
//! Build with: `wasm-pack build --target web --features wasm`

use wasm_bindgen::prelude::*;
use serde::Serialize;

use crate::atmosphere;
use crate::constants::*;
use crate::engines::EngineDatabase;
use crate::kepler;
use crate::orbits;
use crate::propulsion;

// ---------- Data transfer structs ----------

#[derive(Serialize)]
struct TrajectoryPoint {
    t: f64,
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    alt: f64,
    vel: f64,
    mass: f64,
    q: f64,
}

#[derive(Serialize)]
struct OrbitPoint {
    x: f64,
    y: f64,
    r: f64,
    theta: f64,
    v: f64,
}

#[derive(Serialize)]
struct HohmannResult {
    dv1: f64,
    dv2: f64,
    total_dv: f64,
    transfer_time: f64,
    a_transfer: f64,
}

#[derive(Serialize)]
struct ReentryPoint {
    t: f64,
    alt: f64,
    vel: f64,
    decel_g: f64,
    heat_flux: f64,
    temp: f64,
}

#[derive(Serialize)]
struct LandingPoint {
    t: f64,
    x: f64,
    z: f64,
    vx: f64,
    vz: f64,
    mass: f64,
    thrust: f64,
}

#[derive(Serialize)]
struct LambertResult {
    v1x: f64,
    v1y: f64,
    v1z: f64,
    v2x: f64,
    v2y: f64,
    v2z: f64,
    a: f64,
    e: f64,
    transfer_angle: f64,
}

#[derive(Serialize)]
struct EngineInfo {
    name: String,
    manufacturer: String,
    thrust_sl: Option<f64>,
    thrust_vac: f64,
    isp_sl: Option<f64>,
    isp_vac: f64,
    mass: f64,
    propellant: String,
    cycle: String,
    chamber_pressure: f64,
    thrust_to_weight: f64,
    first_flight: u32,
    status: String,
}

// ---------- Exported functions ----------

/// Simulate a 2-D launch trajectory using RK4 integration.
#[wasm_bindgen]
pub fn simulate_trajectory(
    thrust: f64,
    total_mass: f64,
    dry_mass: f64,
    drag_coeff: f64,
    angle_deg: f64,
    dt: f64,
    steps: usize,
) -> JsValue {
    let vehicle = crate::trajectory::Vehicle {
        dry_mass,
        propellant_mass: total_mass - dry_mass,
        thrust_sl: thrust,
        thrust_vac: thrust * 1.1,
        isp_sl: 282.0,
        isp_vac: 311.0,
        cd: drag_coeff,
        reference_area: 10.0,
        pitch_program: None,
    };

    let config = crate::trajectory::SimConfig {
        dt,
        max_time: dt * steps as f64,
        kickover_altitude: 500.0,
        kickover_angle: angle_deg.to_radians(),
        record_interval: 1,
    };

    let result = crate::trajectory::simulate(&vehicle, &config);

    let points: Vec<TrajectoryPoint> = result
        .states
        .iter()
        .map(|s| {
            let q = atmosphere::us_standard_atmosphere(s.y.max(0.0))
                .map(|a| atmosphere::dynamic_pressure(s.speed(), a.density))
                .unwrap_or(0.0);
            TrajectoryPoint {
                t: s.time,
                x: s.x,
                y: s.y,
                vx: s.vx,
                vy: s.vy,
                alt: s.altitude(),
                vel: s.speed(),
                mass: s.mass,
                q,
            }
        })
        .collect();

    serde_wasm_bindgen::to_value(&points).unwrap_or(JsValue::NULL)
}

/// Compute orbital positions for given Keplerian elements.
#[wasm_bindgen]
pub fn compute_orbit(semi_major: f64, eccentricity: f64, inclination_deg: f64, steps: usize) -> JsValue {
    let mu = MU_EARTH;
    let e = eccentricity.clamp(0.0, 0.999);
    let _inc = inclination_deg.to_radians();
    let p = semi_major * (1.0 - e * e);

    let points: Vec<OrbitPoint> = (0..=steps)
        .map(|i| {
            let theta = 2.0 * PI * i as f64 / steps as f64;
            let r = p / (1.0 + e * theta.cos());
            let v = orbits::vis_viva(r, semi_major, mu);
            OrbitPoint {
                x: r * theta.cos(),
                y: r * theta.sin(),
                r,
                theta,
                v,
            }
        })
        .collect();

    serde_wasm_bindgen::to_value(&points).unwrap_or(JsValue::NULL)
}

/// Hohmann transfer between two circular orbits (radii in km).
#[wasm_bindgen]
pub fn hohmann_transfer(r1: f64, r2: f64) -> JsValue {
    let (dv1, dv2, total, time) = orbits::hohmann_transfer(r1, r2, MU_EARTH);
    let result = HohmannResult {
        dv1,
        dv2,
        total_dv: total,
        transfer_time: time,
        a_transfer: (r1 + r2) / 2.0,
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Simplified atmospheric re-entry simulation.
#[wasm_bindgen]
pub fn reentry_sim(
    entry_velocity: f64,
    entry_angle_deg: f64,
    mass: f64,
    cd: f64,
    nose_radius: f64,
    steps: usize,
) -> JsValue {
    let dt = 1.0;
    let gamma = entry_angle_deg.to_radians();
    let area = PI * nose_radius * nose_radius;
    let mut alt = 120_000.0;
    let mut v = entry_velocity;
    let mut g = gamma;
    let mut t = 0.0;

    let mut points = Vec::with_capacity(steps);

    for _ in 0..steps {
        if alt < 0.0 || v < 50.0 {
            break;
        }
        let r = R_EARTH_EQUATORIAL * 1000.0 + alt;
        let grav = MU_EARTH * 1e9 / (r * r);

        let (rho, _pres, _temp) = atmosphere::us_standard_atmosphere(alt.max(0.0).min(85999.0))
            .map(|a| (a.density, a.pressure, a.temperature))
            .unwrap_or((0.0, 0.0, 200.0));

        let drag = 0.5 * rho * v * v * cd * area;
        let decel = drag / mass;
        let decel_g = decel / G0;
        let heat_flux = 1.7415e-4 * (rho / nose_radius).sqrt() * v.powi(3);
        let sigma = 5.67e-8;
        let temp_stag = (heat_flux / (0.85 * sigma)).powf(0.25);

        points.push(ReentryPoint {
            t,
            alt,
            vel: v,
            decel_g,
            heat_flux,
            temp: temp_stag,
        });

        let dv_drag = -decel;
        let dv_grav = -grav * g.sin();
        let dg = -(grav / v) * g.cos() + (v / r) * g.cos();

        v += (dv_drag + dv_grav) * dt;
        g += dg * dt;
        alt += v * g.sin() * dt;
        t += dt;
    }

    serde_wasm_bindgen::to_value(&points).unwrap_or(JsValue::NULL)
}

/// Falcon 9 style powered descent / landing simulation.
#[wasm_bindgen]
pub fn falcon9_landing(
    altitude: f64,
    velocity: f64,
    fuel_mass: f64,
    thrust: f64,
    steps: usize,
) -> JsValue {
    let dry_mass = 22200.0;
    let isp = 282.0;
    let ve = isp * G0;
    let dt = 0.1;
    let mut z = altitude;
    let mut vz = -velocity.abs();
    let mut m = dry_mass + fuel_mass;
    let mut t = 0.0;
    let x = 0.0;
    let mut vx = 0.0_f64;

    let mut points = Vec::with_capacity(steps);

    for _ in 0..steps {
        if z <= 0.0 {
            points.push(LandingPoint { t, x, z: 0.0, vx, vz, mass: m, thrust: 0.0 });
            break;
        }
        let has_fuel = m > dry_mass;
        let thr = if has_fuel && z < altitude * 0.6 { thrust.min(m * G0 * 3.0) } else { 0.0 };
        let mdot = if thr > 0.0 { thr / ve } else { 0.0 };

        points.push(LandingPoint { t, x, z, vx, vz, mass: m, thrust: thr });

        let az = if thr > 0.0 { thr / m } else { 0.0 } - G0;
        vz += az * dt;
        z += vz * dt;
        m -= mdot * dt;
        if m < dry_mass { m = dry_mass; }
        t += dt;
    }

    serde_wasm_bindgen::to_value(&points).unwrap_or(JsValue::NULL)
}

/// Solve Lambert's problem (positions in km, time of flight in seconds).
#[wasm_bindgen]
pub fn lambert_solve(
    r1x: f64, r1y: f64, r1z: f64,
    r2x: f64, r2y: f64, r2z: f64,
    tof: f64,
) -> JsValue {
    let r1 = nalgebra::Vector3::new(r1x, r1y, r1z);
    let r2 = nalgebra::Vector3::new(r2x, r2y, r2z);
    match crate::lambert::solve_lambert(r1, r2, tof, MU_SUN, true) {
        Ok(sol) => {
            let result = LambertResult {
                v1x: sol.v1.x, v1y: sol.v1.y, v1z: sol.v1.z,
                v2x: sol.v2.x, v2y: sol.v2.y, v2z: sol.v2.z,
                a: sol.a,
                e: sol.e,
                transfer_angle: sol.transfer_angle,
            };
            serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
        }
        Err(_) => JsValue::NULL,
    }
}

/// Return the full engine database as JSON.
#[wasm_bindgen]
pub fn engine_database() -> JsValue {
    let db = EngineDatabase::new();
    let engines: Vec<EngineInfo> = db
        .list_engines()
        .into_iter()
        .filter_map(|name| {
            let e = db.get_engine(name)?;
            Some(EngineInfo {
                name: e.name.clone(),
                manufacturer: e.manufacturer.clone(),
                thrust_sl: e.thrust_sl,
                thrust_vac: e.thrust_vac,
                isp_sl: e.isp_sl,
                isp_vac: e.isp_vac,
                mass: e.mass,
                propellant: format!("{:?}", e.propellant),
                cycle: format!("{:?}", e.cycle),
                chamber_pressure: e.chamber_pressure,
                thrust_to_weight: e.thrust_to_weight,
                first_flight: e.first_flight,
                status: format!("{:?}", e.status),
            })
        })
        .collect();
    serde_wasm_bindgen::to_value(&engines).unwrap_or(JsValue::NULL)
}

/// Propagate a TLE with (simplified) SGP4 and return position/velocity.
#[wasm_bindgen]
pub fn propagate_sgp4(tle_line0: &str, tle_line1: &str, tle_line2: &str, minutes: f64) -> JsValue {
    let tle_str = format!("{}\n{}\n{}", tle_line0.trim(), tle_line1.trim(), tle_line2.trim());
    match crate::sgp4::TLE::from_str(&tle_str) {
        Ok(tle) => match crate::sgp4::SGP4::new(tle) {
            Ok(sgp4) => match sgp4.propagate_from_epoch(minutes) {
                Ok(state) => {
                    #[derive(Serialize)]
                    struct Sgp4Result {
                        x: f64, y: f64, z: f64,
                        vx: f64, vy: f64, vz: f64,
                        altitude: f64,
                    }
                    let r = &state.position;
                    let v = &state.velocity;
                    let alt = r.magnitude() - R_EARTH;
                    let res = Sgp4Result {
                        x: r.x, y: r.y, z: r.z,
                        vx: v.x, vy: v.y, vz: v.z,
                        altitude: alt,
                    };
                    serde_wasm_bindgen::to_value(&res).unwrap_or(JsValue::NULL)
                }
                Err(_) => JsValue::NULL,
            },
            Err(_) => JsValue::NULL,
        },
        Err(_) => JsValue::NULL,
    }
}
