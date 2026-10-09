# Requirements, provisional decisions, and milestones

## Recorded user requirements

- 1-D piecewise-smooth models with gradients, interfaces, and low-velocity zones.
- Physically excited multimode time-domain surface waves and paired model/data generation.
- Aim to outperform CPS at matched accuracy, without assuming a speedup.
- Probabilistic wavefield inversion for anisotropic shear structure.
- Lean Python API; robust efficient numerical core may use another language.
- Scholte waves with water as a core capability.
- Fundamental-mode frequency-dependent dv/v forward prediction and depth inference.

Approval covers the first numerical slice and a new public Denolle-Lab repository.
It does not settle the field dv/v estimator, anisotropy symmetry, deployment band,
hardware, or dataset scale. This release's direct phase observable is an explicit
supported contract, not an assumption about the user's unprovided measurements.

## Provisional architecture decisions

Share physical models/modes across synthesis, observation operators, sensitivities,
and inversion. Keep a transparent NumPy/SciPy reference and optional production
adapters. Pydantic v2 owns public validation contracts. No LLM owns numerical QC.
VTI is the recommended next anisotropy scope, not a user-selected symmetry.
SpecSWD is a reuse candidate pending license, build, independent numerical,
mode/residue, and waveform interface validation. No external backend is mandatory.
Plotting, HDF5, seismic formats, and ML remain optional dependencies.

## Milestones and gates

| Phase | Deliverable | Gate |
| --- | --- | --- |
| 1 (current) | Contracts, dry/marine reference modes, small fundamental phase depth operator/inverse | Offline analytical, boundary, convergence, derivative, and Gaussian conditioning tests |
| 2 | Broader reference/branch verification, group derivatives, baseline uncertainty | LVZ/cutoff/degeneracy tests, independent comparator cases, uncertainty budget |
| 3 | Force/moment/acoustic excitation and physical multimode sensor waveforms | Residue/energy normalization, reciprocity, FFT, amplitude/phase, paper/CPS fixtures |
| 4 | VTI and production backend evaluation | Isotropic limit, complete stiffness validation, matched accuracy, measured workloads |
| 5 | Estimator-specific monitoring and full probabilistic wavefield inference | Actual bands/windows/covariance, nonlinear checks, posterior prediction/coverage |
| 6 | Small paired dataset pilot; optional learning acceleration | Background-model splits, structural holdouts, calibration, learning-curve justification |
| 7 | General anisotropy, attenuation, source illumination, other physics | Independent gates for each extension |

The current depth basis is regional log-Vs, not independent depth pixels. A
kernel-per-depth API, flexible interface perturbations, and continuous basis
projection follow after derivative verification. Correlated bands and explicit
water nuisance priors are already supported in the small linear inverse.

## Decisions required before field acceptance

1. What is measured dv/v: phase, group/envelope, stretching/MWCS, direct arrival or coda?
2. Baseline model ensemble, reference definition, covariance, units, and target precision.
3. Band/depth scale, water-depth range, sources, geometry, and sensors.
4. VTI radial anisotropy first or immediate tilted/azimuthal anisotropy.
5. Hardware and forward workload for the CPS performance objective.

Time-lapse priors and regularization cannot create missing vertical resolution.
Do not combine waveform/curve/image likelihoods as independent evidence when
derived from the same observations. Ambient correlations require an illumination
operator; they are not automatically point-source Green functions.

### Completed tutorial slice

Five worked synthetic seismology notebooks and a static GitHub Pages build cover
critical-zone monitoring, volcanic elastic-change hypotheses, water–solid
tradeoffs, covariance in time-lapse tomography, and KKFL-S-informed offshore DAS.
The final notebook includes a limited absolute-dispersion Vs MAP experiment.
The tutorials do not turn the core into a production lateral tomography,
structural inversion, coda sensitivity, cable-mechanics, or waveform package.

Field integration can next accept locally constrained background models and
measured dispersion/covariance, including channel geometry and bathymetry for
KKFL-S. Absolute/baseline model uncertainty, general structural inversion,
branch continuation near degeneracy, and estimator-specific coda sensitivity
remain separate milestones.
