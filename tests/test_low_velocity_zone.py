from gesc import EarthModel, IsotropicMaterial, ModeRequest, SolidRegion, solve_modes


def test_buried_low_velocity_zone_keeps_independent_family_counts(solid):
    cap = IsotropicMaterial(vp_m_s=3000, vs_m_s=1500, density_kg_m3=2300)
    low = IsotropicMaterial(vp_m_s=2000, vs_m_s=800, density_kg_m3=2200)
    model = EarthModel(
        regions=(
            SolidRegion(top_m=0, bottom_m=200, top=cap, bottom=cap),
            SolidRegion(top_m=200, bottom_m=600, top=low, bottom=low),
        ),
        halfspace=solid,
        provenance="buried LVZ verification fixture",
        density_strategy="explicit",
    )
    request = ModeRequest(
        frequencies_hz=(0.75,),
        families=("love", "rayleigh"),
        max_modes=4,
        bottom_depth_m=12000,
        points_per_region=32,
    )
    love, rayleigh = solve_modes(model, request).rows
    assert len(love.modes) == 1
    assert len(rayleigh.modes) == 2
    assert love.status == rayleigh.status == "incomplete"
    assert all(m.qc.status == "passed" for m in (*love.modes, *rayleigh.modes))
