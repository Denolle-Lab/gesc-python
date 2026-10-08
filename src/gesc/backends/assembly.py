"""Dimensionless strong-form collocation with displacement and acoustic pressure.

Convention: exp(i k x - i omega t), ux=i U, uz=W. Q0+k Q1+k² Q2=0.
Finite bottom is locked; no claim of exact radiation or body-wave completeness.
"""

from dataclasses import dataclass

import numpy as np

from ..models import EarthModel


def chebyshev(n: int, top: float, bottom: float):
    x = np.cos(np.pi * np.arange(n) / (n - 1))
    c = np.ones(n)
    c[[0, -1]] = 2
    c *= (-1.0) ** np.arange(n)
    difference = x[:, None] - x[None, :]
    d = (c[:, None] / c[None, :]) / (difference + np.eye(n))
    d -= np.diag(d.sum(axis=1))
    # z increases as x decreases, so the derivative includes a minus sign.
    return top + (1 - x) * (bottom - top) / 2, -2 * d / (bottom - top)


@dataclass
class Block:
    z: np.ndarray
    d: np.ndarray
    rho: np.ndarray
    mu: np.ndarray
    lam: np.ndarray


@dataclass
class Mesh:
    blocks: list[Block]
    length: float
    speed: float
    density: float
    water_z: np.ndarray | None
    water_d: np.ndarray | None


def discretize(model: EarthModel, bottom: float, n: int):
    if bottom <= model.halfspace_top_m:
        raise ValueError("Place the numerical bottom below the explicit half-space interface.")
    length = bottom
    speed = model.halfspace.vs_m_s
    density = model.halfspace.density_kg_m3
    blocks = []
    for region in (*model.regions, None):
        top = model.halfspace_top_m if region is None else region.top_m
        end = bottom if region is None else region.bottom_m
        z, d = chebyshev(n, top / length, end / length)
        if region is None:
            vp, vs, rho = (
                np.full(n, getattr(model.halfspace, key))
                for key in ("vp_m_s", "vs_m_s", "density_kg_m3")
            )
        else:
            # Protect only endpoint roundoff; never extrapolate a physical profile.
            vp, vs, rho = region.sample(np.clip(z * length, top, end))
        rho = rho / density
        mu = rho * (vs / speed) ** 2
        lam = rho * (vp / speed) ** 2 - 2 * mu
        blocks.append(Block(z, d, rho, mu, lam))
    wz = wd = None
    if model.water:
        wz, wd = chebyshev(n, -model.water.depth_m / length, 0)
    return Mesh(blocks, length, speed, density, wz, wd)


def rayleigh_pencil(model: EarthModel, mesh: Mesh, frequency: float):
    n = len(mesh.blocks[0].z)
    solid_size = 2 * n * len(mesh.blocks)
    size = solid_size + (n if model.water else 0)
    q0, q1, q2 = (np.zeros((size, size)) for _ in range(3))
    omega = 2 * np.pi * frequency * mesh.length / mesh.speed
    for i, block in enumerate(mesh.blocks):
        u = slice(2 * n * i, 2 * n * i + n)
        w = slice(u.stop, u.stop + n)
        d = block.d
        mu, lam, rho = (np.diag(x) for x in (block.mu, block.lam, block.rho))
        c = lam + 2 * mu
        q0[u, u] = d @ mu @ d + omega**2 * rho
        q0[w, w] = d @ c @ d + omega**2 * rho
        q1[u, w] = lam @ d + d @ mu
        q1[w, u] = -mu @ d - d @ lam
        q2[u, u], q2[w, w] = -c, -mu

    def clear(row):
        for q in (q0, q1, q2):
            q[row, :] = 0

    def traction(row, i, endpoint, normal, sign=1):
        block = mesh.blocks[i]
        start = 2 * n * i
        u, w = slice(start, start + n), slice(start + n, start + 2 * n)
        if normal:
            q0[row, w] += sign * (block.lam[endpoint] + 2 * block.mu[endpoint]) * block.d[endpoint]
            q1[row, start + endpoint] -= sign * block.lam[endpoint]
        else:
            q0[row, u] += sign * block.mu[endpoint] * block.d[endpoint]
            q1[row, start + n + endpoint] += sign * block.mu[endpoint]

    # Solid surface: shear free; dry normal free or normal stress + p = 0.
    clear(0)
    clear(n)
    traction(0, 0, 0, False)
    traction(n, 0, 0, True)
    if model.water:
        p = slice(solid_size, size)
        rw = model.water.density_kg_m3 / mesh.density
        cw = model.water.sound_speed_m_s / mesh.speed
        q0[p, p] = mesh.water_d @ mesh.water_d + omega**2 / cw**2 * np.eye(n)
        q2[p, p] = -np.eye(n)
        q0[n, size - 1] = 1
        # Pressure-release sea surface; displacement continuity at seabed.
        clear(solid_size)
        q0[solid_size, solid_size] = 1
        clear(size - 1)
        q0[size - 1, n] = 1
        q0[size - 1, p] = -mesh.water_d[-1] / (rw * omega**2)
    # Duplicate interface nodes carry one-sided material values.
    for i in range(len(mesh.blocks) - 1):
        start = 2 * n * i
        next_start = start + 2 * n
        for row, offset in ((start + n - 1, 0), (start + 2 * n - 1, n)):
            clear(row)
            q0[row, start + offset + n - 1] = 1
            q0[row, next_start + offset] = -1
        for row, normal in ((next_start, False), (next_start + n, True)):
            clear(row)
            traction(row, i, n - 1, normal)
            traction(row, i + 1, 0, normal, -1)
    start = solid_size - 2 * n
    for row in (start + n - 1, solid_size - 1):
        clear(row)
        q0[row, row] = 1
    if not all(np.all(np.isfinite(q)) for q in (q0, q1, q2)):
        raise ValueError("Nonfinite Rayleigh assembly; inspect material scales and mesh.")
    return q0, q1, q2


def love_pencil(mesh: Mesh, frequency: float):
    n = len(mesh.blocks[0].z)
    size = n * len(mesh.blocks)
    a, b = np.zeros((size, size)), np.zeros((size, size))
    omega = 2 * np.pi * frequency * mesh.length / mesh.speed
    for i, block in enumerate(mesh.blocks):
        s = slice(i * n, (i + 1) * n)
        a[s, s] = block.d @ np.diag(block.mu) @ block.d + np.diag(block.rho) * omega**2
        b[s, s] = np.diag(block.mu)
    a[0] = b[0] = 0
    a[0, :n] = mesh.blocks[0].mu[0] * mesh.blocks[0].d[0]
    a[-1] = b[-1] = 0
    a[-1, -1] = 1
    for i in range(len(mesh.blocks) - 1):
        left, right = (i + 1) * n - 1, (i + 1) * n
        a[left] = b[left] = a[right] = b[right] = 0
        a[left, left], a[left, right] = 1, -1
        a[right, i * n : (i + 1) * n] = mesh.blocks[i].mu[-1] * mesh.blocks[i].d[-1]
        a[right, (i + 1) * n : (i + 2) * n] = -mesh.blocks[i + 1].mu[0] * mesh.blocks[i + 1].d[0]
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Nonfinite Love assembly; inspect material scales and mesh.")
    return a, b
