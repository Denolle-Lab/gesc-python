# Accuracy-matched CPS benchmark contract

Status: specification only. No CPS timings or speedup measurements exist yet.

Before execution freeze model coefficients, units, provenance, frequency band,
source/receiver geometry, modes, time sampling, filters, attenuation, boundaries,
hardware, compiler/BLAS, threads, comparator archive SHA256, versions, and I/O policy.

## Comparators

- Isotropic: sprep96, sdisp96, slegn96/sregn96, spulse96.
- VTI: tprep96, tdisp96, tlegn96/tregn96, tpulse96.
- Test rlegn96/rregn96 for low-velocity-zone eigenfunctions; document selection.
- Use independent scalar dry/marine equations before cross-solver comparisons.
- Verify marine sensor/source conventions before timing equivalent synthesis.
- Full wavenumber-integration body/near-field responses are separate physical
  deliverables, not an interchangeable timing target for a trapped modal sum.

## Verification matrix

1. Homogeneous dry half-space; locked homogeneous layer; layer over half-space.
2. Paper Table 1, first four Love/Rayleigh eigenfunctions at periods 1 and 3 s.
3. Paper far-field moment-source geometry: source depth 5 km, distance 300 km,
   azimuth 50 degrees, 4096 samples over 1000 s. Obtain complete source pulse
   and source tensor conventions before labeling this an exact reproduction.
4. Gradients, low-velocity zones, thin layers, and near-cutoff higher modes.
5. Finite-water secular benchmark, dry/deep-water limits, and water/sediment contrasts.
6. VTI isotropic limit and joint Love/Rayleigh trade-offs after that backend exists.

Initial acceptance proposals: phase/group relative error 1e-4, normalized
eigenfunction error 1e-3, waveform relative L2 error 1%, plus absolute amplitude
and spectral phase checks. Tighten as acquisition/measurement precision requires.
Use r_max |delta k| to control long-range phase error. For monitoring, numerical
predicted phase-change error should be below 0.1 measurement sigma; derivative
and direct perturbed-model tests are required. Correlation is not amplitude validation.

For smooth profiles refine CPS's constant-layer approximation until it meets
the same physical target. Refine collocation/SEM depth resolution, bottom treatment,
frequency interpolation, and modal truncation separately. Match physical source,
sensor, filtering, attenuation, and modal content. Report all failures/cutoffs.

## Timing matrix

Sweep 64/256/1024 frequencies and 1/4/8 requested modes where supported.
Report dispersion-only, eigenfunctions/excitation, physical synthesis, repeated
synthesis with cached modes, derivative computation, and repeated monitoring
inverse costs separately. Later report cost per effective posterior sample.

Use fixed hardware, one controlled BLAS thread first, then explicit parallel
scaling. Include cold/warm median and variability, peak memory, and failure/mode
completeness rates. Report both kernel-only and end-to-end I/O costs. Amortize
modal caches over source/wavelet/acquisition realizations for both competitors
when supported. A lower-accuracy or failed result does not count as a speedup.

No large model sweep is authorized by this first slice; start with roughly
20-30 verification cases and a small workload pilot after the waveform gate.
