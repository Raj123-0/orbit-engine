# orbit-engine

> High-precision orbital propagators, Lambert targeting, and transfer optimization for astrodynamics

[![CI](https://github.com/Raj123-0/orbit-engine/actions/workflows/ci.yml/badge.svg)](https://github.com/Raj123-0/orbit-engine/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

## Mathematical Foundations

A two-body Keplerian orbit is governed by Newton's law of universal gravitation:
d^2r/dt^2 = -mu * r / |r|^3

Kepler's equation relates mean anomaly M and eccentric anomaly E:
M = E - e * sin(E)

This package provides high-accuracy conversions between classical orbital elements (a, e, i, Omega, omega, nu) and Cartesian state vectors (r, v), along with fast Newton-Raphson Kepler equation solvers.

## Installation

```bash
git clone https://github.com/Raj123-0/orbit-engine.git
cd orbit-engine
pip install -e .
```

## Quickstart

```python
from orbitengine import OrbitalElements, KeplerOrbit, propagate_orbit

orbit = OrbitalElements(
    semi_major_axis=7000.0,  # km
    eccentricity=0.001,
    inclination=0.9,         # radians
    raan=0.0,
    arg_periapsis=0.0,
    true_anomaly=0.0
)

print(f"Orbital period: {orbit.period() / 60:.2f} minutes")

# Propagate forward by 30 minutes
future_orbit = propagate_orbit(orbit, delta_time=1800.0)
cart_r, cart_v = KeplerOrbit(future_orbit).to_cartesian()
print(f"Position at t+30m: {cart_r}")
```

## Running Tests

```bash
pytest -v tests/
```

## License
MIT License. Created autonomously by [Auto'd](https://github.com/Raj123-0/autod).
