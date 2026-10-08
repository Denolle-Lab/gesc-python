"""Small Gaussian depth-basis inverse with correlated data and water nuisance priors."""

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .models import Contract
from .sensitivity import DepthParameter, PhaseOperator


class LinearInverseRequest(Contract):
    data_fractional: tuple[float, ...]
    data_covariance: tuple[tuple[float, ...], ...]
    prior_mean: tuple[float, ...]
    prior_covariance: tuple[tuple[float, ...], ...]
    nuisance_prior_covariance: tuple[tuple[float, ...], ...] | None = None


class PosteriorReport(Contract):
    parameter_names: tuple[str, ...]
    depth_parameters: tuple[DepthParameter, ...]
    mean_log_vs: tuple[float, ...]
    covariance_log_vs: tuple[tuple[float, ...], ...]
    lower_95_log_vs: tuple[float, ...]
    upper_95_log_vs: tuple[float, ...]
    predicted_fractional: tuple[float, ...]
    residual_fractional: tuple[float, ...]
    averaging_matrix: tuple[tuple[float, ...], ...]
    whitened_singular_values: tuple[float, ...]
    prior_dominated: tuple[bool, ...]
    # Nuisance-conditioned posterior prediction includes the water conditional mean.
    nuisance_mean: tuple[float, ...]


def _covariance(value, n, name):
    array = np.asarray(value, dtype=float)
    if array.shape != (n, n) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite {n} x {n} covariance.")
    if not np.allclose(array, array.T, rtol=1e-10, atol=0):
        raise ValueError(f"{name} must be symmetric.")
    try:
        cho_factor(array, lower=True)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be positive definite.") from error
    return array


def infer_time_lapse(operator: PhaseOperator, request: LinearInverseRequest) -> PosteriorReport:
    g, h = np.array(operator.g), np.array(operator.h)
    n, p = g.shape
    d, mean = np.array(request.data_fractional), np.array(request.prior_mean)
    if d.shape != (n,) or mean.shape != (p,) or not np.all(np.isfinite([*d, *mean])):
        raise ValueError("Data/prior mean must match the operator dimensions and be finite.")
    cd = _covariance(request.data_covariance, n, "Data covariance")
    ca = _covariance(request.prior_covariance, p, "Prior covariance")
    if not np.allclose(
        np.sqrt(np.diag(cd)), operator.observation.standard_deviation, rtol=1e-8, atol=0
    ):
        raise ValueError(
            "Data covariance diagonal must match the observation uncertainty contract."
        )
    cb = None
    ce = cd.copy()
    if h.shape[1]:
        if request.nuisance_prior_covariance is None:
            raise ValueError("Marine operators require an explicit water nuisance covariance.")
        cb = _covariance(request.nuisance_prior_covariance, h.shape[1], "Nuisance covariance")
        ce += h @ cb @ h.T
    elif request.nuisance_prior_covariance is not None:
        raise ValueError("This operator has no water nuisance parameters.")
    # Conditioning form avoids explicit inverse matrices, including prior inversion.
    total = g @ ca @ g.T + ce
    factor = cho_factor(total, lower=True)
    gain = cho_solve(factor, g @ ca).T
    post_mean = mean + gain @ (d - g @ mean)
    covariance = ca - gain @ g @ ca
    covariance = (covariance + covariance.T) / 2
    diagonal = np.diag(covariance)
    if np.any(diagonal < -1e-12 * np.max(np.diag(ca))):
        raise ValueError("Posterior covariance lost positivity; inspect conditioning/scales.")
    sd = np.sqrt(np.maximum(diagonal, 0))
    water_mean = np.zeros(h.shape[1])
    if cb is not None:
        water_mean = cb @ h.T @ cho_solve(factor, d - g @ mean)
    prediction = g @ post_mean + h @ water_mean
    whitened = np.linalg.solve(np.linalg.cholesky(ce), g) @ np.linalg.cholesky(ca)
    return PosteriorReport(
        parameter_names=tuple(p.name for p in operator.parameters),
        depth_parameters=operator.parameters,
        mean_log_vs=tuple(post_mean),
        covariance_log_vs=tuple(map(tuple, covariance)),
        lower_95_log_vs=tuple(post_mean - 1.95996398454 * sd),
        upper_95_log_vs=tuple(post_mean + 1.95996398454 * sd),
        predicted_fractional=tuple(prediction),
        residual_fractional=tuple(d - prediction),
        averaging_matrix=tuple(map(tuple, gain @ g)),
        whitened_singular_values=tuple(np.linalg.svd(whitened, compute_uv=False)),
        prior_dominated=tuple(diagonal / np.diag(ca) > 0.8),
        nuisance_mean=tuple(water_mean),
    )
