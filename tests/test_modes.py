import numpy as np
import pytest
from scipy.optimize import brentq

from gesc import EarthModel, ModeRequest, ModeSet, SolidRegion, solve_modes
from gesc.analytic import rayleigh_velocity, scholte_velocity
from gesc.backends.assembly import chebyshev


def test_chebyshev_derivative_sign():
    z, d = chebyshev(20, 0, 3)
    np.testing.assert_allclose(d @ z**4, 4 * z**3, atol=1e-10)


def test_dry_rayleigh_analytic_and_missing_love(dry, mode_request):
    req = ModeRequest.model_validate(
        {**mode_request.model_dump(), "families": ("love", "rayleigh")}
    )
    result = solve_modes(dry, req)
    assert result.rows[0].status == "no_trapped_modes"
    mode = result.fundamental()[0]
    assert mode.phase_velocity_m_s == pytest.approx(rayleigh_velocity(dry.halfspace), rel=1e-6)
    assert mode.qc.boundary_residual < 1e-8
    assert max(abs(np.array(mode.normal_traction_pa))) > 0
    scale = dry.halfspace.density_kg_m3 * dry.halfspace.vs_m_s**2 / mode_request.bottom_depth_m
    assert abs(mode.normal_traction_pa[0]) / scale < 1e-7
    assert abs(mode.shear_traction_pa[0]) / scale < 1e-7
    assert ModeSet.model_validate_json(result.model_dump_json()) == result


@pytest.mark.parametrize("depth", [50, 500, 1500])
def test_finite_water_independent_determinant(marine, mode_request, depth):
    water = type(marine.water).model_validate({**marine.water.model_dump(), "depth_m": depth})
    model = EarthModel.model_validate({**marine.model_dump(), "water": water})
    req = ModeRequest.model_validate({**mode_request.model_dump(), "points_per_region": 48})
    mode = solve_modes(model, req).fundamental()[0]
    exact = scholte_velocity(0.5, model.halfspace, model.water)
    assert mode.phase_velocity_m_s == pytest.approx(exact, rel=2e-6)
    scale = max(abs(np.array(mode.fluid_pressure_pa)))
    assert abs(mode.fluid_pressure_pa[0]) / scale < 1e-7
    assert abs(mode.normal_traction_pa[0] + mode.fluid_pressure_pa[-1]) / scale < 1e-6
    assert mode.fluid_vertical_m[-1] == pytest.approx(mode.vertical[0], abs=1e-6)
    assert abs(mode.shear_traction_pa[0]) / scale < 1e-6


def test_analytic_dry_and_deep_water_limits(solid, water):
    thin = type(water).model_validate({**water.model_dump(), "depth_m": 1e-5})
    assert scholte_velocity(0.5, solid, thin) == pytest.approx(rayleigh_velocity(solid), rel=1e-7)
    deep = type(water).model_validate({**water.model_dump(), "depth_m": 100000})
    deeper = type(water).model_validate({**water.model_dump(), "depth_m": 200000})
    a, b = scholte_velocity(0.5, solid, deep), scholte_velocity(0.5, solid, deeper)
    assert a < min(solid.vs_m_s, water.sound_speed_m_s)
    assert a == pytest.approx(b, rel=1e-10)


def test_love_layer_analytic_and_interface(layered, mode_request):
    req = ModeRequest.model_validate({**mode_request.model_dump(), "families": ("love",)})
    mode = solve_modes(layered, req).fundamental("love")[0]
    upper, lower = layered.regions[0].top, layered.halfspace
    omega, height = np.pi, 500

    def secular(c):
        a = omega * np.sqrt(1 / upper.vs_m_s**2 - 1 / c**2)
        b = omega * np.sqrt(1 / c**2 - 1 / lower.vs_m_s**2)
        mu1 = upper.density_kg_m3 * upper.vs_m_s**2
        mu2 = lower.density_kg_m3 * lower.vs_m_s**2
        return mu1 * a * np.sin(a * height) - mu2 * b * np.cos(a * height)

    exact = brentq(secular, upper.vs_m_s * (1 + 1e-8), lower.vs_m_s * (1 - 1e-8))
    assert mode.phase_velocity_m_s == pytest.approx(exact, rel=2e-6)
    n = req.points_per_region
    assert mode.transverse[n - 1] == pytest.approx(mode.transverse[n], abs=1e-7)
    assert mode.shear_traction_pa[n - 1] == pytest.approx(mode.shear_traction_pa[n], rel=1e-6)


def test_gradient_and_interface_convergence(layered, mode_request):
    original = layered.regions[0]
    region = SolidRegion.model_validate(
        {
            **original.model_dump(),
            "bottom": {**original.bottom.model_dump(), "vs_m_s": 1400},
        }
    )
    model = EarthModel.model_validate({**layered.model_dump(), "regions": (region,)})
    mode = solve_modes(model, mode_request).fundamental()[0]
    assert mode.qc.status == "passed"
    n = mode_request.points_per_region
    assert mode.horizontal[n - 1] == pytest.approx(mode.horizontal[n], abs=1e-7)
    assert mode.vertical[n - 1] == pytest.approx(mode.vertical[n], abs=1e-7)
    assert mode.normal_traction_pa[n - 1] == pytest.approx(mode.normal_traction_pa[n], rel=1e-5)


def test_unconverged_bottom_cannot_enter_monitoring(dry, mode_request):
    shallow = ModeRequest.model_validate({**mode_request.model_dump(), "bottom_depth_m": 1000})
    result = solve_modes(dry, shallow)
    with pytest.raises(ValueError, match="QC"):
        result.fundamental()


def test_invalid_bottom_is_hard_stop(layered, mode_request):
    invalid = ModeRequest.model_validate({**mode_request.model_dump(), "bottom_depth_m": 400})
    with pytest.raises(ValueError, match="below"):
        solve_modes(layered, invalid)


def test_branch_continuation_over_frequencies(dry, mode_request):
    req = ModeRequest.model_validate(
        {**mode_request.model_dump(), "frequencies_hz": (0.5, 0.75, 1)}
    )
    modes = solve_modes(dry, req).fundamental()
    assert all(mode.branch == 0 for mode in modes)
    assert all(mode.qc.tracking_overlap > 0.8 for mode in modes)
