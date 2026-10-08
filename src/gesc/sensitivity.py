"""Fundamental phase-change linearization on physical region coefficients.

Finite differences are the first transparent reference, not an adjoint claim.
Each coefficient is delta log(Vs) throughout one existing region or half-space,
with Vp/density held fixed. Water nuisance coefficients are fractional log changes.
"""

import hashlib
from typing import Literal

import numpy as np
from pydantic import Field, model_validator

from .backends import solve_modes
from .backends.collocation import Candidate, overlap
from .models import Contract, EarthModel, IsotropicMaterial, Positive, SolidRegion
from .modes import ModeRequest
from .observations import PhaseObservation


class DepthParameter(Contract):
    name: str
    top_m: float = Field(ge=0)
    bottom_m: Positive | None
    region_index: int = Field(ge=0)


class PhaseOperator(Contract):
    model_fingerprint: str
    observation: PhaseObservation
    parameters: tuple[DepthParameter, ...]
    g: tuple[tuple[float, ...], ...]
    nuisance_names: tuple[str, ...]
    h: tuple[tuple[float, ...], ...]
    derivative_step: Positive
    step_halving_max_difference: float = Field(ge=0)
    numerical_phase_error_estimate: tuple[float, ...]
    parameterization: Literal["delta_log_vs_fixed_vp_density"] = "delta_log_vs_fixed_vp_density"

    @model_validator(mode="after")
    def shape(self):
        n = len(self.observation.frequencies_hz)
        if len(self.g) != n or any(len(r) != len(self.parameters) for r in self.g):
            raise ValueError("G must have shape (frequency, depth parameter).")
        if len(self.h) != n or any(len(r) != len(self.nuisance_names) for r in self.h):
            raise ValueError("H must have shape (frequency, nuisance parameter).")
        if len(self.numerical_phase_error_estimate) != n:
            raise ValueError("One numerical error diagnostic required per frequency.")
        values = [*np.ravel(self.g), *np.ravel(self.h), *self.numerical_phase_error_estimate]
        if not np.all(np.isfinite(values)) or any(
            x < 0 for x in self.numerical_phase_error_estimate
        ):
            raise ValueError("Operator entries/errors must be finite; errors nonnegative.")
        return self

    def predict(self, coefficients, nuisance=None):
        a = np.asarray(coefficients, dtype=float)
        b = (
            np.zeros(len(self.nuisance_names))
            if nuisance is None
            else np.asarray(nuisance, dtype=float)
        )
        if a.shape != (len(self.parameters),) or b.shape != (len(self.nuisance_names),):
            raise ValueError("Coefficient lengths must match the operator's parameter contracts.")
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise ValueError("Perturbation coefficients must be finite.")
        return np.array(self.g) @ a + np.array(self.h) @ b


def perturb_model(model: EarthModel, coefficients, water_log_changes=None):
    a = np.asarray(coefficients, dtype=float)
    if a.shape != (len(model.regions) + 1,) or not np.all(np.isfinite(a)):
        raise ValueError(
            "One finite log-Vs coefficient per region plus explicit half-space required."
        )

    def change(material, value):
        return IsotropicMaterial.model_validate(
            {
                **material.model_dump(),
                "vs_m_s": material.vs_m_s * np.exp(value),
            }
        )

    regions = tuple(
        SolidRegion.model_validate(
            {
                **r.model_dump(),
                "top": change(r.top, a[i]),
                "bottom": change(r.bottom, a[i]),
            }
        )
        for i, r in enumerate(model.regions)
    )
    water = model.water
    if water_log_changes is not None:
        b = np.asarray(water_log_changes, dtype=float)
        if water is None or b.shape != (3,) or not np.all(np.isfinite(b)):
            raise ValueError(
                "Water perturbations require a column and three finite log coefficients."
            )
        fields = ("depth_m", "sound_speed_m_s", "density_kg_m3")
        water = type(water).model_validate(
            {
                **water.model_dump(),
                **{key: getattr(water, key) * np.exp(v) for key, v in zip(fields, b, strict=False)},
            }
        )
    return EarthModel.model_validate(
        {
            **model.model_dump(),
            "regions": regions,
            "halfspace": change(model.halfspace, a[-1]),
            "water": water,
        }
    )


def _compatible_modes(reference, perturbed, threshold):
    def candidate(mode):
        return Candidate(
            mode.wavenumber_rad_m,
            mode.phase_velocity_m_s,
            *(
                np.array(getattr(mode, key))
                for key in (
                    "depth_m",
                    "horizontal",
                    "vertical",
                    "transverse",
                    "fluid_depth_m",
                    "fluid_pressure_pa",
                    "fluid_vertical_m",
                    "fluid_horizontal_m",
                    "shear_traction_pa",
                    "normal_traction_pa",
                )
            ),
            mode.qc.generalized_residual,
        )

    if any(
        overlap(candidate(a), candidate(b)) < threshold
        for a, b in zip(reference, perturbed, strict=True)
    ):
        raise ValueError(
            "Perturbation changed branch identity; reduce step or constrain background."
        )


def linearize_phase(
    model: EarthModel,
    request: ModeRequest,
    observation: PhaseObservation,
    *,
    step: float = 1e-3,
    step_tolerance: float = 1e-3,
) -> PhaseOperator:
    if not np.isfinite(step) or not 0 < step <= 0.01:
        raise ValueError("Choose finite 0 < finite-difference step <= 0.01 in log units.")
    if not np.isfinite(step_tolerance) or step_tolerance <= 0:
        raise ValueError("Step comparison tolerance must be finite and positive.")
    if (
        request.frequencies_hz != observation.frequencies_hz
        or observation.family not in request.families
    ):
        raise ValueError("Request family/frequencies must match the phase observation exactly.")
    base = solve_modes(model, request).fundamental(observation.family)
    # Refinement deltas are empirical diagnostics, not a rigorous mathematical bound.
    error = np.array(
        [max(x.qc.refinement_relative_error, x.qc.bottom_relative_error) for x in base]
    )
    if np.any(error > 0.1 * np.array(observation.standard_deviation)):
        raise ValueError("Numerical phase error exceeds 0.1 observation sigma; refine mesh/bottom.")
    count = len(model.regions) + 1
    nuisance_names = (
        ("log_water_depth", "log_water_sound_speed", "log_water_density") if model.water else ()
    )
    columns, smaller_columns = [], []
    for j in range(count + len(nuisance_names)):
        derivatives = []
        for delta in (step, step / 2):
            speeds = []
            for sign in (1, -1):
                a, b = np.zeros(count), np.zeros(3)
                if j < count:
                    a[j] = sign * delta
                else:
                    b[j - count] = sign * delta
                changed = perturb_model(model, a, b if model.water else None)
                modes = solve_modes(changed, request).fundamental(observation.family)
                _compatible_modes(base, modes, request.overlap_min)
                changed_error = np.array(
                    [max(m.qc.refinement_relative_error, m.qc.bottom_relative_error) for m in modes]
                )
                if np.any(changed_error > 0.1 * np.array(observation.standard_deviation)):
                    raise ValueError(
                        "Perturbed numerical error exceeds observation budget; refine."
                    )
                speeds.append(np.log([m.phase_velocity_m_s for m in modes]))
            derivatives.append((speeds[0] - speeds[1]) / (2 * delta))
        columns.append(derivatives[0])
        smaller_columns.append(derivatives[1])
    matrix, fine = np.array(columns).T, np.array(smaller_columns).T
    discrepancy = float(np.max(np.abs(matrix - fine)))
    if discrepancy > step_tolerance:
        raise ValueError("Finite-difference step check failed; improve resolution or change step.")
    parameters = tuple(
        DepthParameter(
            name=f"log_vs_region_{i}",
            region_index=i,
            top_m=r.top_m,
            bottom_m=r.bottom_m,
        )
        for i, r in enumerate(model.regions)
    ) + (
        DepthParameter(
            name="log_vs_halfspace",
            region_index=count - 1,
            top_m=model.halfspace_top_m,
            bottom_m=None,
        ),
    )
    return PhaseOperator(
        model_fingerprint=hashlib.sha256(model.model_dump_json().encode()).hexdigest(),
        observation=observation,
        parameters=parameters,
        g=tuple(map(tuple, fine[:, :count])),
        h=tuple(map(tuple, fine[:, count:])),
        nuisance_names=nuisance_names,
        derivative_step=step / 2,
        step_halving_max_difference=discrepancy,
        numerical_phase_error_estimate=tuple(error),
    )
