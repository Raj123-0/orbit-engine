"""Core Keplerian orbit representation and universal variable propagator."""

import math
from dataclasses import dataclass
from typing import Tuple

# Standard Gravitational Parameter of Earth (km^3/s^2)
MU_EARTH = 398600.4418


@dataclass(frozen=True)
class OrbitalElements:
    """Classical Keplerian orbital elements."""
    semi_major_axis: float  # a (km)
    eccentricity: float     # e (0 <= e < 1 for elliptic)
    inclination: float      # i (radians)
    raan: float             # Right Ascension of Ascending Node (radians)
    arg_periapsis: float    # Argument of Periapsis (radians)
    true_anomaly: float     # nu (radians)

    def period(self, mu: float = MU_EARTH) -> float:
        """Calculate orbital period in seconds."""
        if self.semi_major_axis <= 0:
            raise ValueError("Semi-major axis must be positive for periodic orbits.")
        return 2.0 * math.pi * math.sqrt((self.semi_major_axis ** 3) / mu)


class KeplerOrbit:
    """Universal Keplerian orbital solver and state vector converter."""

    def __init__(self, elements: OrbitalElements, mu: float = MU_EARTH):
        self.elements = elements
        self.mu = mu

    def to_cartesian(self) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Convert Keplerian orbital elements to Cartesian position and velocity vectors."""
        e = self.elements.eccentricity
        nu = self.elements.true_anomaly
        a = self.elements.semi_major_axis
        p = a * (1.0 - e**2)

        # Orbital plane coordinates
        r_mag = p / (1.0 + e * math.cos(nu))
        r_orb = (r_mag * math.cos(nu), r_mag * math.sin(nu), 0.0)

        v_mag_factor = math.sqrt(self.mu / p)
        v_orb = (-v_mag_factor * math.sin(nu), v_mag_factor * (e + math.cos(nu)), 0.0)

        # Coordinate transformation Euler angle rotations (RAAN -> inc -> arg_peri)
        raan = self.elements.raan
        inc = self.elements.inclination
        omega = self.elements.arg_periapsis

        p11 = math.cos(raan) * math.cos(omega) - math.sin(raan) * math.sin(omega) * math.cos(inc)
        p12 = -math.cos(raan) * math.sin(omega) - math.sin(raan) * math.cos(omega) * math.cos(inc)
        p21 = math.sin(raan) * math.cos(omega) + math.cos(raan) * math.sin(omega) * math.cos(inc)
        p22 = -math.sin(raan) * math.sin(omega) + math.cos(raan) * math.cos(omega) * math.cos(inc)
        p31 = math.sin(omega) * math.sin(inc)
        p32 = math.cos(omega) * math.sin(inc)

        rx = p11 * r_orb[0] + p12 * r_orb[1]
        ry = p21 * r_orb[0] + p22 * r_orb[1]
        rz = p31 * r_orb[0] + p32 * r_orb[1]

        vx = p11 * v_orb[0] + p12 * v_orb[1]
        vy = p21 * v_orb[0] + p22 * v_orb[1]
        vz = p31 * v_orb[0] + p32 * v_orb[1]

        return (rx, ry, rz), (vx, vy, vz)


def solve_kepler(mean_anomaly: float, eccentricity: float, tol: float = 1e-12, max_iter: int = 100) -> float:
    """Solve Kepler's equation M = E - e*sin(E) using Newton-Raphson."""
    e_anom = mean_anomaly
    for _ in range(max_iter):
        delta = (e_anom - eccentricity * math.sin(e_anom) - mean_anomaly) / (1.0 - eccentricity * math.cos(e_anom))
        e_anom -= delta
        if abs(delta) < tol:
            return e_anom
    return e_anom


def propagate_orbit(elements: OrbitalElements, delta_time: float, mu: float = MU_EARTH) -> OrbitalElements:
    """Propagate Keplerian orbit forward in time by delta_time seconds."""
    a = elements.semi_major_axis
    e = elements.eccentricity
    n = math.sqrt(mu / (a ** 3))

    # Current eccentric anomaly
    nu0 = elements.true_anomaly
    e0 = 2.0 * math.atan(math.sqrt((1.0 - e) / (1.0 + e)) * math.tan(nu0 / 2.0))
    m0 = e0 - e * math.sin(e0)

    # Future mean anomaly
    m_new = m0 + n * delta_time

    # Solve future eccentric anomaly
    e_new = solve_kepler(m_new, e)

    # Convert back to true anomaly
    nu_new = 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(e_new / 2.0), math.sqrt(1.0 - e) * math.cos(e_new / 2.0))

    return OrbitalElements(
        semi_major_axis=a,
        eccentricity=e,
        inclination=elements.inclination,
        raan=elements.raan,
        arg_periapsis=elements.arg_periapsis,
        true_anomaly=nu_new,
    )
