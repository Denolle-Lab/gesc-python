# GESC project conventions

- Keep physical models, numerical discretization, and observations separate.
- Python owns deterministic validation, units, QC, and inference. Use Pydantic v2 contracts.
- Public API uses SI units, depth positive downward, and seabed/dry surface at z=0.
- Water is acoustic; never obtain it by evaluating inverse shear modulus at zero shear speed.
- Reject invalid assembly or unvalidated branches. Never replace nonfinite coefficients with zero.
- Do not claim waveform, anisotropy, coda sensitivity, or CPS speedup from dispersion tests.
- Tests must run offline. Optional external benchmarks require pinned manifests and explicit scope.
- Preserve a transparent numerical reference before optimizing or generating large datasets.
