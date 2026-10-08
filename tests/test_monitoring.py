import numpy as np
import pytest

from gesc import ModeRequest, solve_modes
from gesc.analytic import scholte_velocity
from gesc.observations import PhaseObservation
from gesc.sensitivity import PhaseOperator, linearize_phase, perturb_model


def observation(mode_request):
    return PhaseObservation(
        quantity="fundamental_phase_fractional",
        family="rayleigh",
        frequencies_hz=mode_request.frequencies_hz,
        standard_deviation=(0.001,) * len(mode_request.frequencies_hz),
        units="fraction",
        sign="positive_is_faster",
        window="direct_fundamental",
        reference_id="fixture",
    )


def test_depth_operator_predicts_resolved_perturbation(layered, mode_request):
    req = ModeRequest.model_validate(
        {**mode_request.model_dump(), "frequencies_hz": (0.5, 0.75, 1)}
    )
    op = linearize_phase(layered, req, observation(req))
    assert [(x.top_m, x.bottom_m) for x in op.parameters] == [(0, 500), (500, None)]
    a = np.array([0.0005, -0.0003])
    baseline = solve_modes(layered, req).fundamental()
    perturbed = solve_modes(perturb_model(layered, a), req).fundamental()
    exact = np.array(
        [
            b.phase_velocity_m_s / a.phase_velocity_m_s - 1
            for a, b in zip(baseline, perturbed, strict=True)
        ]
    )
    np.testing.assert_allclose(op.predict(a), exact, atol=1e-6)
    assert op.g[0][0] < op.g[-1][0]  # Higher f is more sensitive to the shallow unit here.
    assert PhaseOperator.model_validate_json(op.model_dump_json()) == op


def test_water_nuisance_derivatives_against_analytic(marine, mode_request):
    op = linearize_phase(marine, mode_request, observation(mode_request))
    assert len(op.nuisance_names) == 3
    for j in range(4):
        values = []
        for sign in (1, -1):
            a, b = np.zeros(1), np.zeros(3)
            (a if j == 0 else b)[0 if j == 0 else j - 1] = sign * 1e-4
            model = perturb_model(marine, a, b)
            values.append(np.log(scholte_velocity(0.5, model.halfspace, model.water)))
        exact = (values[0] - values[1]) / 2e-4
        actual = op.g[0][0] if j == 0 else op.h[0][j - 1]
        assert actual == pytest.approx(exact, abs=3e-4)


def test_numerical_budget_blocks_tiny_measurements(dry, mode_request):
    obs = PhaseObservation.model_validate(
        {**observation(mode_request).model_dump(), "standard_deviation": (1e-12,)}
    )
    with pytest.raises(ValueError, match="sigma"):
        linearize_phase(dry, mode_request, obs)


def test_observable_and_coefficient_mismatch(dry, mode_request):
    obs = PhaseObservation.model_validate(
        {**observation(mode_request).model_dump(), "frequencies_hz": (1,)}
    )
    with pytest.raises(ValueError, match="exactly"):
        linearize_phase(dry, mode_request, obs)
    with pytest.raises(ValueError, match="coefficient"):
        perturb_model(dry, [0, 0])
