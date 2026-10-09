"""Conservation, reversibility, and likelihood equivalence of teaching tomography."""

import numpy as np
import pytest

from tutorials.tomography import gaussian_condition, path_lengths


def test_ray_conserves_length_and_reverses_at_interfaces():
    edges = (0.0, 12000.0, 24000.0)
    for start, end in [
        ((0, 0), (24000, 24000)),
        ((0, 9000), (24000, 17000)),
        ((12000, 0), (12000, 24000)),
    ]:
        lengths = path_lengths(start, end, edges, edges)
        np.testing.assert_allclose(lengths.sum(), np.linalg.norm(np.subtract(end, start)))
        np.testing.assert_allclose(lengths, path_lengths(end, start, edges, edges))
        assert np.all(lengths >= 0)
    np.testing.assert_allclose(
        path_lengths((0, 6000), (24000, 6000), edges, edges), (12000, 12000, 0, 0)
    )


def test_invalid_ray_geometry_stops():
    with pytest.raises(ValueError, match="inside"):
        path_lengths((-1, 0), (2, 2), (0, 1, 2), (0, 1, 2))
    with pytest.raises(ValueError, match="distinct"):
        path_lengths((0, 0), (0, 0), (0, 1), (0, 1))


def test_two_stage_full_covariance_matches_joint_likelihood():
    # Overdetermined rays, a correlated map covariance, and a different depth design.
    ray = np.array([[1, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
    depth = np.array([[0.8, 0.2], [0.4, 0.6]])
    data = np.array([0.1, -0.15, -0.03, 0.06])
    cd = np.eye(4) * 0.01
    ca = np.eye(2) * 0.04
    joint = gaussian_condition(ray @ depth, data, cd, ca)
    left_inverse = np.linalg.solve(ray.T @ ray, ray.T)
    mapped = left_inverse @ data
    mapped_covariance = left_inverse @ cd @ left_inverse.T
    staged = gaussian_condition(depth, mapped, mapped_covariance, ca)
    for a, b in zip(joint, staged, strict=True):
        np.testing.assert_allclose(a, b, atol=1e-12)
    diagonal_only = gaussian_condition(depth, mapped, np.diag(np.diag(mapped_covariance)), ca)
    assert not np.allclose(joint[1], diagonal_only[1])


def test_das_directionality_gauge_null_and_apparent_velocity():
    from tutorials.das import along_cable_velocity, horizontal_p_sv_response

    k = 2 * np.pi / 100
    np.testing.assert_allclose(horizontal_p_sv_response(k, 0, 100), 0, atol=1e-15)
    np.testing.assert_allclose(horizontal_p_sv_response(k, np.pi / 2, 20), 0, atol=1e-15)
    np.testing.assert_allclose(horizontal_p_sv_response(k, 0, 1e-5), -k, rtol=1e-12)
    np.testing.assert_allclose(along_cable_velocity(400, np.pi / 3), 800)
    with pytest.raises(ValueError, match="Broadside"):
        along_cable_velocity(400, np.pi / 2)
