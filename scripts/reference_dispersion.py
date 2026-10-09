"""Regenerate tiny CPS-derived disba phase fixtures; optional, offline once installed."""

import json
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from tutorials.scenarios import cook_inlet, critical_zone, volcano  # noqa: E402


def main():
    from disba import PhaseDispersion

    cases = []
    for factory in (critical_zone, volcano, cook_inlet):
        model, request = factory()
        materials = [r.top for r in model.regions] + [model.halfspace]
        thickness = [r.bottom_m - r.top_m for r in model.regions] + [0]
        vp = [m.vp_m_s for m in materials]
        vs = [m.vs_m_s for m in materials]
        density = [m.density_kg_m3 for m in materials]
        if model.water:
            thickness.insert(0, model.water.depth_m)
            vp.insert(0, model.water.sound_speed_m_s)
            vs.insert(0, 0)
            density.insert(0, model.water.density_kg_m3)
        solver = PhaseDispersion(
            np.array(thickness) / 1000,
            np.array(vp) / 1000,
            np.array(vs) / 1000,
            np.array(density) / 1000,
            algorithm="dunkin",
            dc=0.00005,
        )
        periods = 1 / np.array(request.frequencies_hz)
        order = np.argsort(periods)
        curve = solver(periods[order], mode=0, wave="rayleigh")
        if len(curve.period) != len(periods):
            raise RuntimeError("CPS-derived reference did not find every requested fundamental.")
        np.testing.assert_allclose(curve.period, periods[order])
        velocities = np.empty_like(periods)
        velocities[order] = curve.velocity * 1000
        cases.append(
            {
                "name": factory.__name__,
                "model": model.model_dump(),
                "frequencies_hz": list(request.frequencies_hz),
                "phase_velocity_m_s": velocities.tolist(),
            }
        )
    payload = {
        "reference": "disba, CPS-derived Dunkin phase-dispersion algorithm",
        "disba_version": version("disba"),
        "numba_version": version("numba"),
        "numpy_version": version("numpy"),
        "velocity_search_increment_km_s": 0.00005,
        "scope": "Isotropic fundamental phase speed only; no standalone CPS timing claim",
        "cases": cases,
    }
    path = ROOT / "tests/data/tutorial_phase_reference.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {path}: {sum(len(c['frequencies_hz']) for c in cases)} phase values")


if __name__ == "__main__":
    main()
