"""Comprehensive test suite for orbital propagator."""

import math
import pytest
from orbitengine import KeplerOrbit, OrbitalElements, propagate_orbit
from orbitengine.orbitengine_core import solve_kepler, MU_EARTH


def test_kepler_solution():
    # Mean anomaly 0 -> Eccentric anomaly 0
    assert solve_kepler(0.0, 0.5) == pytest.approx(0.0, abs=1e-9)
    # Circle e=0 -> E = M
    assert solve_kepler(1.234, 0.0) == pytest.approx(1.234, abs=1e-9)


def test_orbit_period():
    elem = OrbitalElements(semi_major_axis=7000.0, eccentricity=0.01, inclination=0.5, raan=0.1, arg_periapsis=0.2, true_anomaly=0.0)
    period = elem.period()
    # Low Earth orbit ~ 5800s
    assert 5700 < period < 6000


def test_to_cartesian():
    elem = OrbitalElements(semi_major_axis=7000.0, eccentricity=0.0, inclination=0.0, raan=0.0, arg_periapsis=0.0, true_anomaly=0.0)
    orbit = KeplerOrbit(elem)
    r, v = orbit.to_cartesian()
    # Circular equatorial orbit at periapsis: r=(7000, 0, 0)
    assert r[0] == pytest.approx(7000.0, abs=1e-3)
    assert r[1] == pytest.approx(0.0, abs=1e-3)
    assert r[2] == pytest.approx(0.0, abs=1e-3)
    # v=(0, sqrt(mu/r), 0)
    expected_v = math.sqrt(MU_EARTH / 7000.0)
    assert v[1] == pytest.approx(expected_v, abs=1e-3)


def test_orbit_propagation_one_period():
    elem = OrbitalElements(semi_major_axis=8000.0, eccentricity=0.05, inclination=0.2, raan=0.1, arg_periapsis=0.3, true_anomaly=0.4)
    period = elem.period()
    propagated = propagate_orbit(elem, period)
    # Propagating full period returns true anomaly modulo 2pi
    diff = (propagated.true_anomaly - elem.true_anomaly) % (2.0 * math.pi)
    assert diff == pytest.approx(0.0, abs=1e-4) or diff == pytest.approx(2.0 * math.pi, abs=1e-4)
