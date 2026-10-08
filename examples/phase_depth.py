"""Explicit fundamental phase change -> two depth coefficients, synthetic only."""

import numpy as np

from gesc import EarthModel, IsotropicMaterial, ModeRequest, SolidRegion, solve_modes
from gesc.inversion import LinearInverseRequest, infer_time_lapse
from gesc.observations import PhaseObservation
from gesc.sensitivity import linearize_phase, perturb_model

upper = IsotropicMaterial(vp_m_s=2200, vs_m_s=1000, density_kg_m3=2200)
lower = IsotropicMaterial(vp_m_s=4000, vs_m_s=2000, density_kg_m3=2500)
model = EarthModel(
    regions=(SolidRegion(top_m=0, bottom_m=500, top=upper, bottom=upper),),
    halfspace=lower,
    provenance="synthetic two-coefficient monitoring example",
    density_strategy="explicit",
)
request = ModeRequest(
    frequencies_hz=(0.5, 0.75, 1.0),
    families=("rayleigh",),
    max_modes=1,
    bottom_depth_m=12000,
    points_per_region=36,
)
observation = PhaseObservation(
    quantity="fundamental_phase_fractional",
    family="rayleigh",
    frequencies_hz=request.frequencies_hz,
    standard_deviation=(0.001,) * 3,
    units="fraction",
    sign="positive_is_faster",
    window="direct_fundamental",
    reference_id="synthetic-background",
)
operator = linearize_phase(model, request, observation)
truth = np.array([0.0005, -0.0003])
baseline = solve_modes(model, request).fundamental()
changed = solve_modes(perturb_model(model, truth), request).fundamental()
data = tuple(
    b.phase_velocity_m_s / a.phase_velocity_m_s - 1 for a, b in zip(baseline, changed, strict=True)
)
report = infer_time_lapse(
    operator,
    LinearInverseRequest(
        data_fractional=data,
        data_covariance=tuple(map(tuple, np.eye(3) * 1e-6)),
        prior_mean=(0.0, 0.0),
        prior_covariance=((1e-4, 0.0), (0.0, 1e-4)),
    ),
)
print("True regional delta log Vs:", truth)
print("Linear prediction:", operator.predict(truth))
print("Re-solved fractional phase changes:", data)
print("Posterior mean:", report.mean_log_vs)
print("95% intervals:", tuple(zip(report.lower_95_log_vs, report.upper_95_log_vs, strict=True)))
print("Depth supports:", [(p.top_m, p.bottom_m) for p in report.depth_parameters])
print("Prior-dominated flags:", report.prior_dominated)
