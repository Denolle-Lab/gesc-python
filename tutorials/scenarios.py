"""Explicit synthetic geology and acquisition choices for the worked notebooks."""

import numpy as np

from gesc import AcousticWaterColumn, EarthModel, IsotropicMaterial, ModeRequest, SolidRegion
from gesc.observations import PhaseObservation
from gesc.sensitivity import PhaseOperator


def material(vp, vs, density):
    return IsotropicMaterial(vp_m_s=vp, vs_m_s=vs, density_kg_m3=density)


def layered_model(depths, materials, provenance):
    regions, top = [], 0
    for depth, value in zip(depths, materials[:-1], strict=True):
        regions.append(SolidRegion(top_m=top, bottom_m=depth, top=value, bottom=value))
        top = depth
    return EarthModel(
        regions=tuple(regions),
        halfspace=materials[-1],
        provenance=provenance,
        density_strategy="explicit",
    )


def critical_zone():
    model = layered_model(
        (10, 50),
        (material(800, 250, 1700), material(1800, 650, 2100), material(3000, 1500, 2400)),
        "Synthetic weathered critical zone; chosen teaching values, not a Shale Hills model",
    )
    request = ModeRequest(
        frequencies_hz=(10.3, 12, 14, 16, 18),
        families=("rayleigh",),
        max_modes=1,
        bottom_depth_m=800,
        points_per_region=36,
    )
    return model, request


def volcano():
    model = layered_model(
        (300, 1500),
        (material(2000, 900, 2000), material(3600, 1800, 2400), material(5000, 2800, 2700)),
        "Synthetic basaltic edifice; chosen teaching values, not a Piton de la Fournaise model",
    )
    request = ModeRequest(
        frequencies_hz=(0.5, 0.65, 0.85, 1.1, 1.45, 1.8),
        families=("rayleigh",),
        max_modes=1,
        bottom_depth_m=16000,
        points_per_region=36,
    )
    return model, request


def coastal():
    model = EarthModel(
        halfspace=material(4000, 2000, 2500),
        water=AcousticWaterColumn(
            depth_m=500,
            sound_speed_m_s=1500,
            density_kg_m3=1000,
            surface_boundary="pressure_release",
        ),
        density_strategy="explicit",
        provenance="Synthetic coastal volcanic-rock halfspace and finite water column",
    )
    request = ModeRequest(
        frequencies_hz=(0.5, 0.75, 1),
        families=("rayleigh",),
        max_modes=1,
        bottom_depth_m=12000,
        points_per_region=48,
    )
    return model, request


def phase_observation(request, sigma=3e-4):
    return PhaseObservation(
        quantity="fundamental_phase_fractional",
        family="rayleigh",
        frequencies_hz=request.frequencies_hz,
        standard_deviation=(sigma,) * len(request.frequencies_hz),
        units="fraction",
        sign="positive_is_faster",
        window="direct_fundamental",
        reference_id="synthetic-reference-model",
    )


def correlated_covariance(size, sigma=3e-4, correlation=0.6):
    if size < 1 or not np.isfinite(sigma) or sigma <= 0 or not -1 < correlation < 1:
        raise ValueError("Use positive size/sigma and |correlation| < 1.")
    indexes = np.arange(size)
    return sigma**2 * correlation ** np.abs(indexes[:, None] - indexes[None, :])


def select_band(operator, indexes):
    indexes = np.asarray(indexes, dtype=int)
    if len(indexes) == 0 or np.any(np.diff(indexes) <= 0):
        raise ValueError("Select a nonempty, strictly increasing band.")
    payload = operator.model_dump()
    for key in ("g", "h", "numerical_phase_error_estimate"):
        payload[key] = [payload[key][i] for i in indexes]
    for key in ("frequencies_hz", "standard_deviation"):
        payload["observation"][key] = [payload["observation"][key][i] for i in indexes]
    return PhaseOperator.model_validate(payload)


def cook_inlet():
    model = layered_model(
        (100,),
        (material(1700, 350, 1900), material(2400, 800, 2200)),
        "Synthetic KKFL-S-inspired sediment; 80 m segment-average water depth "
        "from Ni et al. 2024; geology assumed",
    )
    model = EarthModel.model_validate(
        {
            **model.model_dump(),
            "water": {
                "depth_m": 80,
                "sound_speed_m_s": 1480,
                "density_kg_m3": 1025,
                "surface_boundary": "pressure_release",
            },
        }
    )
    request = ModeRequest(
        frequencies_hz=(1.5, 1.8, 2.1, 2.5, 2.8),
        families=("rayleigh",),
        max_modes=1,
        bottom_depth_m=4000,
        points_per_region=36,
    )
    return model, request
