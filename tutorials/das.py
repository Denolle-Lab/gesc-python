"""Ideal straight horizontal cable response; no cable coupling or excitation model."""

import numpy as np


def horizontal_p_sv_response(wavenumber_rad_m, angle_rad, gauge_length_m):
    """Relative axial-strain transfer per unit seabed longitudinal displacement.

    For ux=i U exp(i k x-i omega t), strain amplitude is -k U cos²(theta)
    times the centered boxcar gauge response sinc(k L cos(theta)/2).
    Returns the real coefficient multiplying U. np.sinc uses a pi-normalized input.
    """
    k, angle = np.asarray(wavenumber_rad_m), np.asarray(angle_rad)
    if (
        not np.all(np.isfinite(k))
        or np.any(k <= 0)
        or not np.all(np.isfinite(angle))
        or not np.isfinite(gauge_length_m)
        or gauge_length_m <= 0
    ):
        raise ValueError("Use finite positive wavenumber/gauge length and finite angle.")
    direction = np.cos(angle)
    return -k * direction**2 * np.sinc(k * gauge_length_m * direction / (2 * np.pi))


def along_cable_velocity(phase_velocity_m_s, angle_rad):
    c, angle = np.asarray(phase_velocity_m_s), np.asarray(angle_rad)
    if not np.all(np.isfinite(c)) or np.any(c <= 0) or not np.all(np.isfinite(angle)):
        raise ValueError("Use finite positive phase velocity and finite angle.")
    direction = np.abs(np.cos(angle))
    if np.any(direction < 1e-8):
        raise ValueError("Broadside propagation has no resolved along-cable phase speed.")
    return c / direction
