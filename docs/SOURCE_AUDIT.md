# Source audit, 2026-10-08

The user approved the first implementation slice and a new public
Denolle-Lab/gesc-python repository. Handoff recommendations were not treated as
proof of capability. Static audits below are distinct from numerical tests of
this new package. No external backend was built or timed.

| Source | Pinned revision | Reuse decision |
| --- | --- | --- |
| [mdenolle/gesc](https://github.com/mdenolle/gesc) | 92bf1655e5e52a9b204fbb274b1a0853a6cdac3f | Mathematical/history reference; no source vendored |
| [Denolle-Lab/PyGesc](https://github.com/Denolle-Lab/PyGesc) | 5acbfa2550af6c98c6cc71cbc815c1d21d96b236 | Unfinished translation; no validated backend reused |
| [Denolle-Lab/surface_wave_sensitivity](https://github.com/Denolle-Lab/surface_wave_sensitivity) | ec500e58650abf30a9f86b4ff1b265d7154ab88a | disba-based example/model ideas, not a monitoring inverse |
| [nqdu/SpecSWD](https://github.com/nqdu/SpecSWD) | e5dc8ce6138c29384a12ec180c74eace6a8773c9 | Production candidate after independent validation/license check |
| [keurfonluu/disba](https://github.com/keurfonluu/disba) | 1a9e26d49dbaadead84cdf427d6f51b3f5c3707b | Optional isotropic diagnostic; CPS-derived algorithm |

## Paper

Denolle et al. (2012), DOI 10.1785/0120110183, was read from the supplied PDF.
SHA256: 1d7aebf9b2d34da54ddd056b7bddd8c1642fcc5d8f0ab5a0734f322d77bf07ab.
It matches the audited GESC repository's paper copy. Printed pp. 1215-1217
describe singular A u=k B u, modified stresses, interfaces, locked bottom, and
MATLAB/LAPACK QZ. Pp. 1218-1221 give analytical/layered/CPS tests and far-field
waveforms. P. 1223 explicitly excludes implemented anelasticity. These examples
do not establish Scholte/general-anisotropy support or a measured CPS speedup.
The paper is not redistributed with this software.

## MATLAB GESC

Paths below are relative to gesc/sub at the pinned revision.

- `EigRW.m:36,80` and `EigLW.m:21,56`: full eigendecomposition before requested-mode retention.
- `EigRW.m:8-11`, `EigLW.m:11-14`: nonfinite matrix entries replaced by zero.
- `EigRW.m:102-106`, `EigLW.m:69-71`: normalization by a surface component.
- `gesc.m:186-193`: skips failures and couples family count through a minimum.
- `gesc.m:119-127,181`: fundamental-only adaptation and point-count early return.
- `gesc.m:143-145`: half-space coefficients taken from profile maxima.
- `RayInnerMat.m:10`, `LoveInnerMat.m:16`: inverse mu; no valid zero-shear acoustic block.
- `RayBC.m:53-57`: water branch makes no active change; water edits elsewhere are commented.
- `get_integrals_sw.m:34-47`: dry energy integrals and phase/group quantities exist.
- `get_sw_gf.m:132`: multimode Rayleigh entry uses Love source eigenfunction `uyi`.
- `get_sw_gf.m:186`, `get_sw_dc_farfield.m:244`: multimode inverse FFT selects one frequency.
- `get_sw_gf_farfield.m:89-96`: structure variable C overwritten; S2 is undefined.
- Source/moment/modal Green routines exist beyond the paper's far-field tests;
  their existence is not validation of all branches, amplitudes, or geometries.

No explicit license file or automated numerical test suite was found in that tree.
The new backend rederives the equations and does not copy these routines.

## Existing lab Python work

PyGesc has an explicit MIT license. `gesc.py:395-408` calls missing EigRW,
EigLW, LoveInnerMat, LoveBC, and get_integrals_sw routines. The driver also retains
MATLAB identifiers, dict/attribute mismatches, and undefined bigC construction.
`gesc.py:430,435-437` uses elementwise multiplication where the MATLAB expression
requires matrix products. It cannot serve as a validated end-to-end backend.

The sensitivity repository has MIT licensing and two SensKernels notebooks.
They invoke disba PhaseDispersion, EigenFunction, and PhaseSensitivity for
isotropic layered examples. Useful examples exist, but no estimator-aware
probabilistic dv/v-to-depth API was found in those inspected notebooks.
Geospatial datasets and notebook graphics are not copied into this repository.

## SpecSWD and disba

SpecSWD includes fluid/coupling assembly (`src/vti/assemble.cpp`), fluid and solid
displacement conversion (`src/vti/transform.cpp`), anisotropic assembly, and
phase/group kernel code. Tests include analytical dry Rayleigh roots and
finite-difference acoustic/elastic kernel comparisons in `tests/unit_tests/rayl_op.cpp`.
Tests were inspected, not run. `src/shared/schur.cpp:113-127` uses Euclidean
eigenvector normalization; physical excitation residues need separate validation.
No complete physical waveform-synthesis API was identified in the inspected Python package.

Its build requires C++17, Fortran, Eigen, LAPACK/LAPACKE, HDF5, and pybind11.
Python metadata requires NumPy/SciPy/h5py/Matplotlib. Metadata declares MIT,
but no explicit license file was found and separately attributed code is present.
Adoption requires upstream clarification and file-level provenance checks.
GoogleTest is fetched from unpinned main by CMake; reproducible evaluation must pin it.

disba has an explicit BSD-3-Clause license; dispersion, eigenfunction, ellipticity,
and sensitivity APIs and water-layer tests exist. Its runtime includes NumPy,
Numba, Matplotlib, and importlib_metadata. It is not an anisotropic waveform solver.

## CPS

The [official site](https://rbherrmann.github.io/ComputerProgramsSeismology/) lists
the April 25, 2026 distribution and an explicit
[MIT license](https://rbherrmann.github.io/ComputerProgramsSeismology/mitlicense.html).
The [program inventory](https://rbherrmann.github.io/ComputerProgramsSeismology/cps_programs.html)
documents isotropic/VTI modal synthesis and generalized-reflection eigenfunctions.
[Ocean examples](https://rbherrmann.github.io/ComputerProgramsSeismology/TUTORIAL/OCEAN1tut/DOITSW.html)
demonstrate water-layer dispersion. The source archive has not been audited/built
here, and pressure/source normalization parity remains a benchmark gate.
