import numpy as np
import pytest
from pydantic import ValidationError

from gesc import EarthModel, IsotropicMaterial, ModeRequest, SolidRegion
from gesc.observations import PhaseObservation


@pytest.mark.parametrize("value", [0, -1, np.inf, np.nan])
def test_invalid_solid_is_not_water(value):
    with pytest.raises(ValidationError):
        IsotropicMaterial(vp_m_s=4000, vs_m_s=value, density_kg_m3=2500)


def test_physical_admissibility():
    with pytest.raises(ValidationError, match="bulk modulus"):
        IsotropicMaterial(vp_m_s=1000, vs_m_s=1000, density_kg_m3=2500)


def test_required_provenance_and_density(dry):
    data = dry.model_dump()
    del data["density_strategy"]
    with pytest.raises(ValidationError):
        EarthModel.model_validate(data)
    with pytest.raises(ValidationError, match="provenance"):
        EarthModel.model_validate({**dry.model_dump(), "provenance": " "})


def test_no_implicit_region_gaps_or_extrapolation(solid, dry):
    region = SolidRegion(top_m=10, bottom_m=100, top=solid, bottom=solid)
    with pytest.raises(ValidationError, match="contiguous"):
        EarthModel.model_validate({**dry.model_dump(), "regions": (region,)})
    with pytest.raises(ValueError, match="within"):
        region.sample([9])


def test_one_sided_interface_values(layered):
    assert layered.regions[0].sample([500])[1][0] == 1000
    assert layered.halfspace.vs_m_s == 2000


def test_request_requires_monotonic_frequencies(mode_request):
    with pytest.raises(ValidationError, match="increasing"):
        ModeRequest.model_validate({**mode_request.model_dump(), "frequencies_hz": (1, 0.5)})


def test_reject_undefined_dvv_interpretation():
    data = dict(
        quantity="fundamental_phase_fractional",
        family="rayleigh",
        frequencies_hz=(1,),
        standard_deviation=(0.001,),
        units="fraction",
        sign="positive_is_faster",
        window="direct_fundamental",
        reference_id="reference",
    )
    for key, value in (("window", "coda"), ("quantity", "stretching"), ("units", "percent")):
        with pytest.raises(ValidationError):
            PhaseObservation(**{**data, key: value})
