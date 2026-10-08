"""Explicit SI physical models; no mesh coordinates in the geology contract."""

from typing import Annotated, Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

Positive = Annotated[float, Field(gt=0, allow_inf_nan=False)]
Nonnegative = Annotated[float, Field(ge=0, allow_inf_nan=False)]


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, allow_inf_nan=False)


class IsotropicMaterial(Contract):
    vp_m_s: Positive
    vs_m_s: Positive
    density_kg_m3: Positive

    @model_validator(mode="after")
    def admissible(self):
        if self.vp_m_s**2 <= 4 / 3 * self.vs_m_s**2:
            raise ValueError("Positive bulk modulus required: set Vp² > 4 Vs² / 3.")
        return self


class SolidRegion(Contract):
    """A linear profile of Vp, Vs, and density with one-sided interface values."""

    top_m: Nonnegative
    bottom_m: Positive
    top: IsotropicMaterial
    bottom: IsotropicMaterial

    @model_validator(mode="after")
    def ordered(self):
        if self.bottom_m <= self.top_m:
            raise ValueError("Region bottom must be deeper than its top.")
        return self

    def sample(self, depth_m):
        z = np.asarray(depth_m, dtype=float)
        if not np.all(np.isfinite(z)) or np.any((z < self.top_m) | (z > self.bottom_m)):
            raise ValueError("Sample within this region; select an interface side explicitly.")
        t = (z - self.top_m) / (self.bottom_m - self.top_m)
        return tuple(
            getattr(self.top, key) * (1 - t) + getattr(self.bottom, key) * t
            for key in ("vp_m_s", "vs_m_s", "density_kg_m3")
        )


class AcousticWaterColumn(Contract):
    depth_m: Positive
    sound_speed_m_s: Positive
    density_kg_m3: Positive
    surface_boundary: Literal["pressure_release"]


class EarthModel(Contract):
    regions: tuple[SolidRegion, ...] = ()
    halfspace: IsotropicMaterial
    water: AcousticWaterColumn | None = None
    provenance: str = Field(min_length=1)
    density_strategy: Literal["explicit"]
    units: Literal["SI"] = "SI"
    symmetry: Literal["isotropic"] = "isotropic"

    @model_validator(mode="after")
    def contiguous(self):
        if not self.provenance.strip():
            raise ValueError("Provide nonblank model provenance.")
        expected = 0.0
        for region in self.regions:
            if region.top_m != expected:
                raise ValueError("Solid regions must be contiguous, ordered, and start at z=0.")
            expected = region.bottom_m
        return self

    @property
    def halfspace_top_m(self):
        return self.regions[-1].bottom_m if self.regions else 0.0
