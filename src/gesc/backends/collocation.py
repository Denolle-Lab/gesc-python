"""Small dense reference backend. Every returned mode includes refinement QC."""

import hashlib
from dataclasses import dataclass

import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import eig
from scipy.optimize import linear_sum_assignment

from ..models import EarthModel
from ..modes import FrequencyModes, Mode, ModeQC, ModeRequest, ModeSet
from .assembly import Mesh, discretize, love_pencil, rayleigh_pencil


@dataclass
class Candidate:
    k: float
    c: float
    z: np.ndarray
    u: np.ndarray
    w: np.ndarray
    t: np.ndarray
    wz: np.ndarray
    p: np.ndarray
    ww: np.ndarray
    wu: np.ndarray
    shear: np.ndarray
    normal: np.ndarray
    residual: float
    boundary_residual: float = 0.0


def _spectrum(model: EarthModel, mesh: Mesh, frequency: float, family: str):
    omega = 2 * np.pi * frequency * mesh.length / mesh.speed
    n = len(mesh.blocks[0].z)
    z = np.concatenate([b.z for b in mesh.blocks]) * mesh.length
    if family == "love":
        q0, q2 = love_pencil(mesh, frequency)
        row_scale = np.maximum(
            np.max(np.abs(q0), axis=1) + omega**2 * np.max(np.abs(q2), axis=1), 1e-300
        )
        roots, vectors = eig(q0 / row_scale[:, None], q2 / row_scale[:, None])
        q1 = np.zeros_like(q0)
        # For Love A u = k² B u; retain the polynomial convention below.
        q2 = -q2
        roots = np.sqrt(roots.astype(complex))
    else:
        q0, q1, q2 = rayleigh_pencil(model, mesh, frequency)
        size = len(q0)
        zero, identity = np.zeros_like(q0), np.eye(size)
        # Left row equilibration leaves the polynomial and right eigenvectors
        # unchanged, while avoiding shallow-water derivative rows dominating QZ.
        row_scale = np.maximum(
            np.max(np.abs(q0), axis=1)
            + omega * np.max(np.abs(q1), axis=1)
            + omega**2 * np.max(np.abs(q2), axis=1),
            1e-300,
        )
        a = np.block([[zero, identity], [-q0 / row_scale[:, None], -q1 / row_scale[:, None]]])
        b = np.block([[identity, zero], [zero, q2 / row_scale[:, None]]])
        roots, full_vectors = eig(a, b, check_finite=True)
        vectors = full_vectors[:size]
    candidates = []
    for j, root in enumerate(roots):
        if not np.isfinite(root) or root.real <= 0 or abs(root.imag) > 1e-7 * root.real:
            continue
        k = float(root.real)
        c = omega / k * mesh.speed
        # A locked bottom also creates body/box modes; retain trapped candidates only.
        if c >= model.halfspace.vs_m_s * (1 - 1e-7):
            continue
        v = vectors[:, j]
        solid_size = len(z) if family == "love" else 2 * len(z)
        anchor = np.argmax(np.abs(v[:solid_size]))
        if abs(v[anchor]) == 0:
            continue
        v = v * np.exp(-1j * np.angle(v[anchor]))
        if np.linalg.norm(v.imag) > 1e-6 * np.linalg.norm(v.real):
            continue
        v = v.real
        q = q0 + k * q1 + k**2 * q2
        denominator = (
            np.linalg.norm(q0, axis=1)
            + k * np.linalg.norm(q1, axis=1)
            + k**2 * np.linalg.norm(q2, axis=1)
        ) * np.linalg.norm(v)
        scaled = np.abs(q @ v) / np.maximum(denominator, 1e-300)
        residual = float(np.max(scaled))
        boundary_rows = [0, len(q0) - 1]
        if family == "rayleigh":
            boundary_rows += [n, 2 * len(z) - n - 1, 2 * len(z) - 1]
            if model.water:
                boundary_rows += [2 * len(z)]
            for i in range(len(mesh.blocks) - 1):
                boundary_rows += [
                    2 * n * i + n - 1,
                    2 * n * i + 2 * n - 1,
                    2 * n * (i + 1),
                    2 * n * (i + 1) + n,
                ]
        else:
            for i in range(len(mesh.blocks) - 1):
                boundary_rows += [(i + 1) * n - 1, (i + 1) * n]
        boundary_residual = float(np.max(scaled[boundary_rows]))
        # Reject assembly/eigen artifacts before branch matching.
        if residual > 1e-6:
            continue
        zero = np.zeros_like(z)
        wz, p, ww, wu = (np.array([]) for _ in range(4))
        if family == "love":
            u, w, t = zero, zero, v
        else:
            u = np.concatenate([v[2 * n * i : 2 * n * i + n] for i in range(len(mesh.blocks))])
            w = np.concatenate(
                [v[2 * n * i + n : 2 * n * i + 2 * n] for i in range(len(mesh.blocks))]
            )
            t = zero
            if model.water:
                p = v[solid_size:]
                rw = model.water.density_kg_m3 / mesh.density
                ww = mesh.water_d @ p / (rw * omega**2)
                wu = k * p / (rw * omega**2)
                wz = mesh.water_z * mesh.length
        # Physical displacement gauge: maximum solid vector displacement = 1 m.
        scale = np.max(np.sqrt(u**2 + w**2 + t**2))
        u, w, t = u / scale, w / scale, t / scale
        if model.water and family == "rayleigh":
            p = p / scale * mesh.density * mesh.speed**2 / mesh.length
            ww = ww / scale
            wu = wu / scale
        shear, normal = [], []
        traction_scale = mesh.density * mesh.speed**2 / mesh.length
        for i, block in enumerate(mesh.blocks):
            s = slice(i * n, (i + 1) * n)
            if family == "love":
                shear.extend(traction_scale * block.mu * (block.d @ t[s]))
                normal.extend(np.zeros(n))
            else:
                shear.extend(traction_scale * block.mu * (block.d @ u[s] + k * w[s]))
                normal.extend(
                    traction_scale
                    * ((block.lam + 2 * block.mu) * (block.d @ w[s]) - k * block.lam * u[s])
                )
        candidates.append(
            Candidate(
                k / mesh.length,
                c,
                z,
                u,
                w,
                t,
                wz,
                p,
                ww,
                wu,
                np.array(shear),
                np.array(normal),
                residual,
                boundary_residual,
            )
        )
    return sorted(candidates, key=lambda mode: mode.c)


def overlap(first: Candidate, second: Candidate):
    """Solid displacement overlap only; no fluid-energy normalization claim."""
    end = min(first.z[-1], second.z[-1])
    # Cluster near the surface to resolve decaying modes even with a deep bottom.
    z = end * np.linspace(0, 1, 600) ** 2

    def fields(candidate):
        unique, indexes = np.unique(candidate.z, return_index=True)
        return np.array(
            [np.interp(z, unique, x[indexes]) for x in (candidate.u, candidate.w, candidate.t)]
        )

    a, b = fields(first), fields(second)
    numerator = trapezoid(np.sum(a * b, axis=0), z)
    denominator = np.sqrt(trapezoid(np.sum(a * a, axis=0), z) * trapezoid(np.sum(b * b, axis=0), z))
    return min(1.0, float(abs(numerator) / denominator)) if denominator > 0 else 0.0


def _match(reference, candidates):
    if not reference or not candidates:
        return {}
    score = np.array([[overlap(a, b) for b in candidates] for a in reference])
    rows, columns = linear_sum_assignment(-score)
    return {i: (candidates[j], score[i, j]) for i, j in zip(rows, columns, strict=False)}


def solve_modes(model: EarthModel, request: ModeRequest) -> ModeSet:
    """Solve dry/marine isotropic trapped candidates with independent family counts.

    Initial branches are ordered by phase speed; later frequencies use overlap assignment.
    This is a reference for simple isolated branches, not a certified mode census near
    degeneracy/cutoff. A failed or absent branch cannot enter the monitoring operator.
    """
    meshes = (
        discretize(model, request.bottom_depth_m, request.points_per_region),
        discretize(model, request.bottom_depth_m, request.points_per_region + 8),
        discretize(model, request.bottom_depth_m * 1.5, request.points_per_region + 16),
    )
    rows = []
    for family in request.families:
        previous = []
        for frequency in request.frequencies_hz:
            spectra = [_spectrum(model, mesh, frequency, family) for mesh in meshes]
            base = spectra[0][: request.max_modes]
            tracking = {i: 1.0 for i in range(len(base))}
            fundamental_rank_changed = False
            if previous:
                matches = _match(previous, spectra[0])
                # Keep branch slots, rather than relabeling an overtone as fundamental.
                if len(matches) < len(previous):
                    rows.append(
                        FrequencyModes(
                            frequency_hz=frequency,
                            family=family,
                            modes=(),
                            status="incomplete",
                            message="Tracked branch absent; shorten band or refine resolution.",
                        )
                    )
                    continue
                base = [matches[i][0] for i in range(len(previous))]
                tracking = {i: matches[i][1] for i in range(len(previous))}
                # A high overlap can jump to an overtone across an avoided crossing.
                # Branch zero denotes the lowest-speed candidate in this reference.
                # Require both continuation and rank; do not certify the overtone
                # as fundamental when frequency steps skip a rapid shape change.
                fundamental_rank_changed = base[0] is not spectra[0][0]
            fine, deep = (_match(base, spectrum) for spectrum in spectra[1:])
            modes = []
            for i, candidate in enumerate(base):
                refined, of = fine.get(i, (None, 0))
                deeper, od = deep.get(i, (None, 0))
                refinement = abs(refined.c / candidate.c - 1) if refined else 1.0
                bottom_error = abs(deeper.c / candidate.c - 1) if deeper else 1.0
                score = min(tracking[i], of, od)
                ambiguous = score < request.overlap_min or (i == 0 and fundamental_rank_changed)
                # Competing base eigenfunctions must be distinguishable by overlap.
                if any(
                    overlap(candidate, other) > 0.98
                    for other in spectra[0]
                    if other is not candidate
                ):
                    ambiguous = True
                passed = (
                    candidate.residual <= request.residual_tolerance
                    and max(refinement, bottom_error) <= request.convergence_tolerance
                )
                qc = ModeQC(
                    generalized_residual=candidate.residual,
                    boundary_residual=candidate.boundary_residual,
                    refinement_relative_error=refinement,
                    bottom_relative_error=bottom_error,
                    tracking_overlap=score,
                    status="ambiguous" if ambiguous else "passed" if passed else "failed",
                )
                modes.append(
                    Mode(
                        frequency_hz=frequency,
                        family=family,
                        branch=i,
                        phase_velocity_m_s=candidate.c,
                        wavenumber_rad_m=candidate.k,
                        depth_m=tuple(candidate.z),
                        horizontal=tuple(candidate.u),
                        vertical=tuple(candidate.w),
                        transverse=tuple(candidate.t),
                        shear_traction_pa=tuple(candidate.shear),
                        normal_traction_pa=tuple(candidate.normal),
                        fluid_depth_m=tuple(candidate.wz),
                        fluid_pressure_pa=tuple(candidate.p),
                        fluid_vertical_m=tuple(candidate.ww),
                        qc=qc,
                        fluid_horizontal_m=tuple(candidate.wu),
                    )
                )
            previous = base
            status = "passed"
            if not modes:
                status = "no_trapped_modes"
            elif any(m.qc.status == "ambiguous" for m in modes):
                status = "ambiguous"
            elif any(m.qc.status == "failed" for m in modes):
                status = "failed"
            elif len(modes) < request.max_modes:
                status = "incomplete"
            rows.append(
                FrequencyModes(
                    frequency_hz=frequency,
                    family=family,
                    modes=tuple(modes),
                    status=status,
                    message=(
                        "Refine points/bottom or reduce band if QC fails. Counts are independent; "
                        "requested-mode completeness is not a proof of the entire spectrum."
                    ),
                )
            )
    fingerprint = hashlib.sha256(model.model_dump_json().encode()).hexdigest()
    return ModeSet(model_fingerprint=fingerprint, request=request, rows=tuple(rows))
