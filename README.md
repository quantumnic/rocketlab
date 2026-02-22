# 🚀 rocketlab

Rocket trajectory simulator and aerospace engineering toolkit — orbital mechanics, propulsion, atmospheric models, pure Rust.

## Features

- **Orbital Mechanics**: Kepler's equation, orbital elements ↔ state vectors, vis-viva, Hohmann transfers
- **Propulsion**: Tsiolkovsky equation, specific impulse, thrust calculations
- **CLI**: Interactive command-line interface for quick calculations

## Usage

```bash
# Delta-v calculation
rocketlab deltav --isp 300 --mass-ratio 8

# Orbital elements to state vectors
rocketlab orbit --elements '6778,0.001,51.6,0,0,0'
```

## Build

```bash
cargo build --release
```

## References

- Bate, Mueller & White — *Fundamentals of Astrodynamics*
- Curtis — *Orbital Mechanics for Engineering Students*
