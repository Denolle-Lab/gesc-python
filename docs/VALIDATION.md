# First-slice validation record

Local verification on 2026-10-08: Python 3.12.13, NumPy 2.5.3, SciPy 1.18.1,
Pydantic 2.14.0, one OpenBLAS/OMP thread. `python -m pytest -q`: **30 passed**
in 40.65 seconds. Ruff lint/format and Git whitespace checks passed. Editable
package build/install and both examples passed. CI verifies the same suite on
Linux across Python 3.11, 3.12, and 3.13; check the workflow for remote results.

Coverage includes required model/provenance/density fields, nonfinite/invalid
materials, one-sided interfaces, Chebyshev derivative sign, dry Rayleigh cubic,
Love layer secular equation, finite-water determinant at 50/500/1500 m,
analytical dry/deep-water limits, acoustic-solid displacement/traction continuity,
gradient refinement, a buried LVZ, independent family counts, bottom rejection,
frequency continuation, step-halving, analytical water derivatives, and rejection
of numerical error larger than the stated measurement budget.

The monitoring tests compare the linear prediction with freshly solved
perturbed models. In the two-depth-coefficient example, the imposed log-Vs
changes are 0.0005 and -0.0003. Linear fractional phase predictions at
0.5/0.75/1 Hz are approximately -2.1708e-4, -5.3628e-5, and 6.2618e-4;
re-solved values are -2.1709e-4, -5.3935e-5, and 6.2596e-4.

Inference tests cover scalar Gaussian conditioning, correlated errors, agreement
with independent joint water/shear conditioning, and nominal interval coverage
under a small known linear Gaussian simulator. They do not demonstrate
field-data adequacy or continuous-depth resolving power.

The finite-water example at 500 m gives c=1754.356339 m/s versus an independent
secular root of 1754.356333 m/s at 0.5 Hz. Its 0.75 and 1 Hz examples also pass
QC. These are dispersion/eigenfunction checks, not physical waveform validation.

Shallow-water derivative scaling initially failed the residual gate. Left row
equilibration fixed the case without relaxing the acceptance threshold. All
residual checks continue to use the original unscaled equations.

Not executed: CPS/SpecSWD comparison, physical source/residue waveform synthesis,
VTI/general anisotropy, attenuation, measured coda/stretching/MWCS inversion,
large data generation, or performance benchmarking. See the benchmark plan and
roadmap for those acceptance gates.

## Worked seismology tutorials

The notebook build executes five synthetic experiments in fresh kernels, with
exceptions and scientific assertions blocking publication. It checks finite-
water analytical dispersion, perturbed-model linear prediction, numerical
error relative to observation uncertainty, posterior information after band
selection, full-covariance joint/two-stage tomography equivalence, and an
absolute-dispersion seabed MAP experiment. The Cook Inlet notebook separates
published acquisition facts from synthetic geology and picks.

The tiny `tests/data/tutorial_phase_reference.json` fixture contains fundamental
Rayleigh/Scholte phase speeds for the three layered tutorial models. It was
computed with disba 0.7.0's CPS-derived Dunkin implementation, using a 0.00005
km/s velocity search increment. The model and frequencies are stored with the
values. Tests require agreement within 0.005 m/s and passing GESC mode QC; they
run offline without disba or Numba. This is a phase-only cross-check against a
CPS-derived implementation, not a standalone CPS build/timing benchmark or
waveform validation.

Regenerate the reference intentionally with Python 3.12:

```bash
python -m pip install -r requirements/reference.txt
python scripts/reference_dispersion.py
```

Geometry tests check conservation/reversal of cell path lengths, joint/staged
Gaussian equivalence, DAS broadside and gauge nulls, the short-gauge limit, and
apparent-velocity projection. A regression test prevents a maximum-overlap
continuation from publishing an overtone as a QC-passing fundamental.

`requirements/tutorials.txt` pins the optional notebook environment separately
from core dependency constraints. The build writes package versions, source
hashes, figures, assertion metrics, and execution times to `validation.json`.
The GitHub Actions workflow executes the notebooks and tests before deploying
GitHub Pages. Scientific calculations do not use network access; installation
and deployment do.
