"""Teaching-only straight-ray path integrals and linear Gaussian conditioning.

This is a four-cell acquisition experiment, not a production tomography backend.
Coordinates/lengths are metres. Rays are straight and the background is laterally uniform.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve


def path_lengths(start, end, x_edges, y_edges):
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    x_edges, y_edges = np.asarray(x_edges, dtype=float), np.asarray(y_edges, dtype=float)
    if (
        start.shape != (2,)
        or end.shape != (2,)
        or not np.all(np.isfinite([start, end]))
        or not np.all(np.isfinite(x_edges))
        or not np.all(np.isfinite(y_edges))
        or len(x_edges) < 2
        or len(y_edges) < 2
        or np.any(np.diff(x_edges) <= 0)
        or np.any(np.diff(y_edges) <= 0)
    ):
        raise ValueError("Use finite 2-D endpoints and strictly increasing cell edges.")
    for point in (start, end):
        if not (x_edges[0] <= point[0] <= x_edges[-1] and y_edges[0] <= point[1] <= y_edges[-1]):
            raise ValueError("Both stations must lie inside the cell domain.")
    delta = end - start
    distance = np.linalg.norm(delta)
    if distance == 0:
        raise ValueError("Stations must be distinct.")
    intersections = [0.0, 1.0]
    for axis, edges in enumerate((x_edges, y_edges)):
        if delta[axis] != 0:
            values = (edges - start[axis]) / delta[axis]
            intersections.extend(values[(values > 0) & (values < 1)])
    intersections = np.unique(intersections)
    lengths = np.zeros((len(y_edges) - 1, len(x_edges) - 1))
    for left, right in zip(intersections[:-1], intersections[1:], strict=True):
        middle = start + (left + right) / 2 * delta
        col = np.clip(np.searchsorted(x_edges, middle[0], side="right") - 1, 0, len(x_edges) - 2)
        row = np.clip(np.searchsorted(y_edges, middle[1], side="right") - 1, 0, len(y_edges) - 2)
        lengths[row, col] += (right - left) * distance
    return lengths.ravel()


def gaussian_condition(design, data, covariance, prior_covariance):
    design, data = np.asarray(design), np.asarray(data)
    covariance, prior = np.asarray(covariance), np.asarray(prior_covariance)
    n, p = design.shape
    if data.shape != (n,) or covariance.shape != (n, n) or prior.shape != (p, p):
        raise ValueError("Data, covariance, and prior must match the design dimensions.")
    if not all(np.all(np.isfinite(value)) for value in (design, data, covariance, prior)):
        raise ValueError("All inversion entries must be finite.")
    for value in (covariance, prior):
        if not np.allclose(value, value.T, rtol=1e-10, atol=0):
            raise ValueError("Covariances must be symmetric positive definite.")
        cho_factor(value, lower=True)
    factor = cho_factor(design @ prior @ design.T + covariance, lower=True)
    gain = cho_solve(factor, design @ prior).T
    mean = gain @ data
    posterior = prior - gain @ design @ prior
    return mean, (posterior + posterior.T) / 2, gain @ design
