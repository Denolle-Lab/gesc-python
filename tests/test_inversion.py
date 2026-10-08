import numpy as np
import pytest

from gesc.inversion import LinearInverseRequest, infer_time_lapse
from gesc.observations import PhaseObservation
from gesc.sensitivity import DepthParameter, PhaseOperator


def operator(g, h=None):
    n, p = np.shape(g)
    h = np.zeros((n, 0)) if h is None else np.array(h)
    return PhaseOperator(
        model_fingerprint="synthetic-linear-fixture",
        observation=PhaseObservation(
            quantity="fundamental_phase_fractional",
            family="rayleigh",
            frequencies_hz=tuple(range(1, n + 1)),
            standard_deviation=(0.01,) * n,
            units="fraction",
            sign="positive_is_faster",
            window="direct_fundamental",
            reference_id="synthetic",
        ),
        parameters=tuple(
            DepthParameter(name=f"q{i}", top_m=i * 100, bottom_m=(i + 1) * 100, region_index=i)
            for i in range(p)
        ),
        g=tuple(map(tuple, g)),
        h=tuple(map(tuple, h)),
        nuisance_names=tuple(f"water{i}" for i in range(h.shape[1])),
        derivative_step=1e-3,
        step_halving_max_difference=0,
        numerical_phase_error_estimate=(0.0,) * n,
    )


def test_scalar_gaussian_conditioning_and_prior_dominated():
    op = operator([[1.0, 0.0]])
    report = infer_time_lapse(
        op,
        LinearInverseRequest(
            data_fractional=(0.03,),
            data_covariance=((0.0001,),),
            prior_mean=(0.0, 0.0),
            prior_covariance=((0.0004, 0.0), (0.0, 0.0004)),
        ),
    )
    assert report.mean_log_vs == pytest.approx((0.024, 0.0))
    assert report.covariance_log_vs[0][0] == pytest.approx(0.00008)
    assert report.prior_dominated == (False, True)
    assert report.averaging_matrix[0][0] == pytest.approx(0.8)


def test_correlated_errors_and_water_marginalization():
    g, h = np.array([[1.0, 0.2], [0.3, 1.0]]), np.array([[1.0], [1.0]])
    ca = np.diag([0.0004, 0.0009])
    cd, cb = np.array([[0.0001, 0.00003], [0.00003, 0.0001]]), np.array([[0.0002]])
    d = np.array([0.01, -0.02])
    report = infer_time_lapse(
        operator(g, h),
        LinearInverseRequest(
            data_fractional=tuple(d),
            data_covariance=tuple(map(tuple, cd)),
            prior_mean=(0.0, 0.0),
            prior_covariance=tuple(map(tuple, ca)),
            nuisance_prior_covariance=tuple(map(tuple, cb)),
        ),
    )
    # Independent joint block conditioning, then extract the shear marginal.
    design = np.hstack([g, h])
    prior = np.diag([0.0004, 0.0009, 0.0002])
    covariance = np.linalg.solve(
        np.linalg.solve(prior, np.eye(3)) + design.T @ np.linalg.solve(cd, design), np.eye(3)
    )
    mean = covariance @ design.T @ np.linalg.solve(cd, d)
    np.testing.assert_allclose(report.mean_log_vs, mean[:2], atol=1e-12)
    np.testing.assert_allclose(report.covariance_log_vs, covariance[:2, :2], atol=1e-12)
    np.testing.assert_allclose(report.nuisance_mean, mean[2:], atol=1e-12)


def test_gaussian_interval_coverage():
    rng = np.random.default_rng(12)
    op = operator([[1.0]])
    inside = 0
    for _ in range(300):
        truth = rng.normal(0, 0.02)
        datum = truth + rng.normal(0, 0.01)
        report = infer_time_lapse(
            op,
            LinearInverseRequest(
                data_fractional=(datum,),
                data_covariance=((0.0001,),),
                prior_mean=(0.0,),
                prior_covariance=((0.0004,),),
            ),
        )
        inside += report.lower_95_log_vs[0] <= truth <= report.upper_95_log_vs[0]
    assert 0.9 < inside / 300 < 0.99


def test_invalid_covariance_and_missing_water_prior():
    op = operator([[1.0]], [[1.0]])
    request = LinearInverseRequest(
        data_fractional=(0.0,),
        data_covariance=((0.0001,),),
        prior_mean=(0.0,),
        prior_covariance=((0.0004,),),
    )
    with pytest.raises(ValueError, match="water"):
        infer_time_lapse(op, request)
    invalid = request.model_copy(update={"prior_covariance": ((-1.0,),)})
    with pytest.raises(ValueError, match="positive definite"):
        infer_time_lapse(op, invalid)
