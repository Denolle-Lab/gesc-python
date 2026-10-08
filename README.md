# gesc-python

A transparent Python reference for surface-wave modes and fundamental
phase-change depth inference, developed by the Denolle lab.

**Status: first validation slice.** This is an isotropic, elastic, flat 1-D
collocation reference with finite locked-bottom convergence checks. It includes
acoustic water from the outset. It does not yet synthesize physical waveforms,
implement VTI/general anisotropy, or establish a speed advantage over CPS.

## Implemented

- Explicit SI, piecewise-linear solid profiles and one-sided interfaces;
  independently specified homogeneous half-space, density, and provenance.
- Love displacement and Rayleigh displacement/acoustic-pressure collocation.
- Finite, homogeneous water with pressure-release surface and fluid-solid coupling.
- Solid displacement/traction and fluid pressure/displacement mode fields.
- Independent family counts, overlap continuation, residuals, mesh refinement,
  and deeper-bottom comparisons. Ambiguous/unconverged modes cannot enter monitoring.
- Independent dry Rayleigh and finite-water secular equations for verification.
- Finite-difference fundamental phase sensitivities to region log-Vs and water
  depth/sound speed/density; step-halving and numerical measurement-budget checks.
- A small Gaussian depth-basis inverse with correlated errors, explicit water
  nuisance covariance, intervals, averaging matrix, and prior-dominance diagnostics.

The first depth basis is one constant log-Vs change per physical solid region,
plus the whole homogeneous half-space. Vp and solid density are held fixed.
This is not a continuous pointwise sensitivity-kernel implementation.

## Install and run

```sh
python -m pip install -e '.[dev]'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q
python examples/marine_modes.py
python examples/phase_depth.py
```

Runtime dependencies are NumPy, SciPy, and Pydantic v2. Tests run offline after
dependencies are installed. Plotting, HDF5, CPS, SpecSWD, MATLAB, and ML frameworks
are not required. The import name is `gesc`; the distribution name is `gesc-python`.

## Conventions and limits

Public inputs use meters, meters/second, kilograms/cubic meter, seconds, and hertz.
Depth increases downward. The dry surface or seabed is z=0; water occupies
[-water_depth, 0]. Wave fields use exp(i k x - i omega t), with returned real
horizontal U corresponding to physical ux=i U and returned vertical W to uz=W.
The same convention applies in water. Love transverse displacement is real;
returned shear traction has the corresponding implicit i factor for Rayleigh.

Modes use a maximum solid vector-displacement gauge of 1 m. Pressure is in Pa
and traction in Pa under that gauge. These are relative eigenfunction fields,
not excited physical amplitudes. Energy/residue normalization is a later gate.

The initial spectrum is ordered by phase speed; continuation uses solid
displacement overlap and one-to-one assignment. The initial fundamental label is
validated on simple isolated dry/layered/marine fixtures, not guaranteed for every
pathological profile. This is not a complete mode census or a leaky-mode solver.
Use closely spaced frequencies near crossings/cutoffs. For multimode work,
inspect row completeness and every mode's QC; an absent mode is never padded.

Refinement and deeper-bottom deltas are empirical numerical-error estimates,
not rigorous bounds. The bottom extension changes both depth and point count;
it is a conservative additional comparison, not an isolated bottom-error proof.
Tiny monitoring changes require a tighter, validated numerical budget.

`PhaseObservation` accepts only direct fundamental fractional phase changes,
positive for faster propagation. Group delay, stretching/MWCS, coda windows,
percent values, attenuation, and azimuthal anisotropy are rejected or outside this
release. Water nuisance priors are required for marine inference. Posterior
coverage under the linear simulator does not validate field-data adequacy.

## Documentation

- [Numerical formulation](docs/NUMERICS.md)
- [Pinned source audit and reuse decisions](docs/SOURCE_AUDIT.md)
- [Accuracy-matched CPS benchmark contract](docs/BENCHMARK_PLAN.md)
- [Requirements and phased milestones](docs/ROADMAP.md)
- [First-slice validation record](docs/VALIDATION.md)

## Scientific foundation

Denolle, Dunham, and Beroza (2012), *Solving the Surface-Wave Eigenproblem with
Chebyshev Spectral Collocation*, BSSA 102(3), 1214-1223,
[doi:10.1785/0120110183](https://doi.org/10.1785/0120110183).
This implementation rederives a displacement quadratic pencil; it is not a
literal translation of the paper's modified-stress linear pencil or the MATLAB
code. No existing solver source or copyrighted paper PDF is vendored here.
