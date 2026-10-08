"""This release accepts fundamental phase changes only; other estimators are explicit gaps."""

from typing import Literal

from pydantic import Field, model_validator

from .models import Contract, Positive


class PhaseObservation(Contract):
    quantity: Literal["fundamental_phase_fractional"]
    family: Literal["rayleigh", "love"]
    frequencies_hz: tuple[Positive, ...] = Field(min_length=1)
    standard_deviation: tuple[Positive, ...] = Field(min_length=1)
    units: Literal["fraction"]
    sign: Literal["positive_is_faster"]
    window: Literal["direct_fundamental"]
    reference_id: str = Field(min_length=1)

    @model_validator(mode="after")
    def shape(self):
        if len(self.frequencies_hz) != len(self.standard_deviation):
            raise ValueError("Provide one fractional uncertainty per frequency.")
        if any(b <= a for a, b in zip(self.frequencies_hz, self.frequencies_hz[1:], strict=False)):
            raise ValueError("Observation frequencies must be increasing and unique.")
        if not self.reference_id.strip():
            raise ValueError("Provide a nonblank reference identifier.")
        return self
