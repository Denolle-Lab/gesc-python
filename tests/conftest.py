import pytest

from gesc import AcousticWaterColumn, EarthModel, IsotropicMaterial, ModeRequest, SolidRegion


@pytest.fixture
def solid():
    return IsotropicMaterial(vp_m_s=4000, vs_m_s=2000, density_kg_m3=2500)


@pytest.fixture
def water():
    return AcousticWaterColumn(
        depth_m=500,
        sound_speed_m_s=1500,
        density_kg_m3=1000,
        surface_boundary="pressure_release",
    )


@pytest.fixture
def dry(solid):
    return EarthModel(
        halfspace=solid, provenance="homogeneous fixture", density_strategy="explicit"
    )


@pytest.fixture
def marine(dry, water):
    return EarthModel.model_validate({**dry.model_dump(), "water": water})


@pytest.fixture
def layered(solid):
    upper = IsotropicMaterial(vp_m_s=2200, vs_m_s=1000, density_kg_m3=2200)
    return EarthModel(
        regions=(SolidRegion(top_m=0, bottom_m=500, top=upper, bottom=upper),),
        halfspace=solid,
        provenance="one-layer fixture",
        density_strategy="explicit",
    )


@pytest.fixture
def mode_request():
    return ModeRequest(
        frequencies_hz=(0.5,),
        families=("rayleigh",),
        max_modes=1,
        bottom_depth_m=12000,
        points_per_region=36,
    )
