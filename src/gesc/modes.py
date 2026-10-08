"""Requests, QC, and ragged modal outputs. Normalization is not excitation."""

from typing import Literal

from pydantic import Field, model_validator

from .models import Contract, Positive


class ModeRequest(Contract):
    frequencies_hz: tuple[Positive, ...] = Field(min_length=1)
    families: tuple[Literal["love", "rayleigh"], ...]
    max_modes: int = Field(ge=1, le=16)
    bottom_depth_m: Positive
    points_per_region: int = Field(ge=12, le=160)
    residual_tolerance: Positive = 1e-8
    convergence_tolerance: Positive = 1e-3
    overlap_min: float = Field(default=0.8, gt=0, lt=1)

    @model_validator(mode="after")
    def ordered(self):
        if any(b <= a for a, b in zip(self.frequencies_hz, self.frequencies_hz[1:], strict=False)):
            raise ValueError("Frequencies must be unique and strictly increasing.")
        if not self.families or len(set(self.families)) != len(self.families):
            raise ValueError("Request at least one unique family.")
        return self


class ModeQC(Contract):
    generalized_residual: float = Field(ge=0)
    boundary_residual: float = Field(ge=0)
    refinement_relative_error: float = Field(ge=0)
    bottom_relative_error: float = Field(ge=0)
    tracking_overlap: float = Field(ge=0, le=1.0000001)
    status: Literal["passed", "failed", "ambiguous"]


class Mode(Contract):
    frequency_hz: Positive
    family: Literal["love", "rayleigh"]
    branch: int = Field(ge=0)
    phase_velocity_m_s: Positive
    wavenumber_rad_m: Positive
    depth_m: tuple[float, ...]
    horizontal: tuple[float, ...]
    vertical: tuple[float, ...]
    transverse: tuple[float, ...]
    shear_traction_pa: tuple[float, ...]
    normal_traction_pa: tuple[float, ...]
    fluid_depth_m: tuple[float, ...] = ()
    fluid_pressure_pa: tuple[float, ...] = ()
    fluid_vertical_m: tuple[float, ...] = ()
    fluid_horizontal_m: tuple[float, ...] = ()
    normalization: Literal["solid_displacement_max_1m"] = "solid_displacement_max_1m"
    qc: ModeQC


class FrequencyModes(Contract):
    frequency_hz: Positive
    family: Literal["love", "rayleigh"]
    modes: tuple[Mode, ...]
    status: Literal["passed", "failed", "no_trapped_modes", "incomplete", "ambiguous"]
    message: str
    # Resolution/bottom checks do not prove a mathematically exhaustive spectrum.
    completeness: Literal["requested_modes_only"] = "requested_modes_only"


class ModeSet(Contract):
    backend: Literal["collocation_reference"] = "collocation_reference"
    model_fingerprint: str
    request: ModeRequest
    rows: tuple[FrequencyModes, ...]

    def fundamental(self, family: Literal["love", "rayleigh"] = "rayleigh"):
        rows = [r for r in self.rows if r.family == family]
        if not rows or any(not r.modes or r.modes[0].qc.status != "passed" for r in rows):
            raise ValueError("Fundamental branch failed QC or is absent; refine mesh/bottom/band.")
        return tuple(r.modes[0] for r in rows)
