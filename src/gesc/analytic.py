"""Independent homogeneous-half-space secular equations for verification.

Finite pressure-release water over an elastic half-space: acoustic load is
rho_w omega² tan(q h)/q. The equation below removes tan poles using sin/cos.
Only the lowest trapped root is returned; it is not a general mode tracker.
"""

import numpy as np
from scipy.optimize import brentq

from .models import AcousticWaterColumn, IsotropicMaterial


def rayleigh_velocity(solid: IsotropicMaterial):
    ratio = (solid.vs_m_s / solid.vp_m_s) ** 2

    def polynomial(y):
        return y**3 - 8 * y**2 + (24 - 16 * ratio) * y - 16 * (1 - ratio)

    return solid.vs_m_s * np.sqrt(brentq(polynomial, 0, 1))


def water_secular(c: float, frequency: float, solid: IsotropicMaterial, water: AcousticWaterColumn):
    if not 0 < c < solid.vs_m_s or not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("Use positive frequency and a trapped speed 0 < c < Vs_halfspace.")
    x = c / solid.vs_m_s
    alpha = np.sqrt(1 - (c / solid.vp_m_s) ** 2)
    beta = np.sqrt(1 - x**2)
    ray = (2 - x**2) ** 2 - 4 * alpha * beta
    kh = 2 * np.pi * frequency * water.depth_m / c
    r = (c / water.sound_speed_m_s) ** 2 - 1
    coefficient = water.density_kg_m3 / solid.density_kg_m3 * x**4 * alpha
    if r < 0:
        s = np.sqrt(-r)
        return (ray + coefficient * np.tanh(kh * s) / s) / x**2
    s = np.sqrt(r)
    # sin(kh*s)/s tends to kh, including at c=c_w.
    return (ray * np.cos(kh * s) + coefficient * kh * np.sinc(kh * s / np.pi)) / x**2


def scholte_velocity(frequency: float, solid: IsotropicMaterial, water: AcousticWaterColumn):
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("Frequency must be finite and positive.")
    grid = np.linspace(
        0.02 * min(solid.vs_m_s, water.sound_speed_m_s), solid.vs_m_s * (1 - 1e-9), 4096
    )
    values = [water_secular(c, frequency, solid, water) for c in grid]
    for a, b, fa, fb in zip(grid, grid[1:], values, values[1:], strict=False):
        if fa * fb < 0:
            return brentq(lambda c: water_secular(c, frequency, solid, water), a, b, xtol=1e-9)
    raise ValueError("No bracketed trapped root; refine the analytic root search for this model.")
