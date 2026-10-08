"""Tiny finite-water demonstration; no waveform or dataset generation."""

from gesc import AcousticWaterColumn, EarthModel, IsotropicMaterial, ModeRequest, solve_modes
from gesc.analytic import scholte_velocity

solid = IsotropicMaterial(vp_m_s=4000, vs_m_s=2000, density_kg_m3=2500)
water = AcousticWaterColumn(
    depth_m=500,
    sound_speed_m_s=1500,
    density_kg_m3=1000,
    surface_boundary="pressure_release",
)
model = EarthModel(
    halfspace=solid,
    water=water,
    provenance="homogeneous marine verification example",
    density_strategy="explicit",
)
request = ModeRequest(
    frequencies_hz=(0.5, 0.75, 1.0),
    families=("love", "rayleigh"),
    max_modes=1,
    bottom_depth_m=12000,
    points_per_region=48,
)
for row in solve_modes(model, request).rows:
    print(row.frequency_hz, row.family, row.status)
    for mode in row.modes:
        exact = scholte_velocity(mode.frequency_hz, solid, water)
        print(f"  c={mode.phase_velocity_m_s:.6f} m/s, analytic={exact:.6f}, QC={mode.qc.status}")
