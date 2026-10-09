"""Field-scale examples must match the pinned CPS-derived fundamental branches."""

import json
from pathlib import Path

import numpy as np
import pytest

from gesc import solve_modes
from tutorials.scenarios import cook_inlet, critical_zone, volcano


@pytest.mark.parametrize("factory", (critical_zone, volcano, cook_inlet))
def test_tutorial_fundamental_matches_cps_derived_fixture(factory):
    payload = json.loads((Path(__file__).parent / "data/tutorial_phase_reference.json").read_text())
    fixture = next(case for case in payload["cases"] if case["name"] == factory.__name__)
    model, request = factory()
    assert model.model_dump(mode="json") == fixture["model"]
    assert list(request.frequencies_hz) == fixture["frequencies_hz"]
    modes = solve_modes(model, request).fundamental()
    # The root search increment is not its refined root precision. Require 5 mm/s,
    # below 0.1 sigma even for the slowest tutorial phase speed.
    np.testing.assert_allclose(
        [m.phase_velocity_m_s for m in modes], fixture["phase_velocity_m_s"], rtol=0, atol=0.005
    )
