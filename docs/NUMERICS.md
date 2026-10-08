# Numerical reference contract

## Solid and water equations

Use exp(i k x - i omega t), ux=i U, uz=W, mu=rho Vs², and
C=lambda+2mu=rho Vp². A prime is a depth derivative. The isotropic solid equations are

```text
(mu U')' + rho omega² U - C k² U + k[lambda W' + (mu W)'] = 0
(C W')' + rho omega² W - mu k² W - k[mu U' + (lambda U)'] = 0
```

Tangential traction is i mu(U'+kW); normal traction is C W'-k lambda U.
Solid interfaces enforce continuity of U, W, and both tractions, preserving
one-sided material values at duplicate collocation nodes.

Water obeys p''+(omega²/cw²-k²)p=0, with homogeneous cw and rho_w.
Its displacements are ux=i k p/(rho_w omega²), uz=p'/(rho_w omega²).
At the seabed impose W=p'/(rho_w omega²), normal stress=-p, and shear=0.
There is no tangential displacement constraint across an inviscid interface.
At z=-h impose p=0. Surface-gravity dynamics are excluded.

Love waves satisfy (mu T')'+rho omega² T=k² mu T, with free shear at the
solid top, including a water-loaded top. Inviscid water has no transverse field.

## Discretization and eigenproblem

Each finite region and the truncated homogeneous half-space has its own
Chebyshev-Lobatto grid. Affine scaling preserves the derivative sign for depth
increasing downward. Sampling never infers the half-space from profile maxima.
Profiles interpolate Vp, Vs, and density linearly within each physical region.

Length, speed, and density scales are the numerical bottom depth, half-space Vs,
and half-space density. Public results are converted back to SI. Rayleigh/water
assembly yields Q0+k Q1+k² Q2=0, in dimensionless variables. The companion pencil is

```text
[ 0   I ] [u ]     [ I  0  ] [u ]
[-Q0 -Q1] [ku] = k [ 0  Q2 ] [ku]
```

Constraint rows make the pencil singular. SciPy/LAPACK generalized eigensolves
handle it directly; B^-1 A is never formed. Infinite, non-real, and nonpositive
roots are filtered, as are candidates faster than half-space Vs. This trapped
filter alone is insufficient; residual, refinement, and bottom tests follow.
No nonfinite assembly entry is replaced by zero.
Left row equilibration is applied before QZ, preserving the polynomial roots and
right eigenvectors while reducing shallow-water derivative scale imbalance.
Residuals are evaluated on the original unscaled polynomial.

The row-scaled polynomial residual is |Q u| divided by
(||Q0_row||+|k| ||Q1_row||+k² ||Q2_row||) ||u||. Boundary QC uses only constraint
rows of that same residual. Physical traction/interface tests provide an
additional, independent check of signs and unit conversion.

## Branch and convergence QC

For every frequency/family, solve the base grid, n+8 grid at the same bottom,
and n+16 grid with a 1.5-times-deeper bottom. Match by displacement overlap;
report phase-speed differences. Track existing branches between frequencies
with one-to-one overlap assignment; do not close missing branch slots by sorting.
Near-identical competing eigenfunctions or low overlaps are marked ambiguous.

Overlap uses interpolated solid displacements, not fluid-energy orthogonality.
Initial ordering and degeneracy handling require broader validation before a
production mode census. No group velocity, excitation residue, or full continuum
is inferred from the current eigenvector normalization.
Only initially found branches are continued; new cut-on branches are not
automatically assigned identities. Request rows remain incomplete in that case.

## Independent marine secular equation

For a homogeneous half-space and finite homogeneous water, define
x=c/Vs, a=sqrt(1-c²/Vp²), b=sqrt(1-c²/Vs²), q²=(omega/c)²(c²/cw²-1).
The solid Rayleigh factor is R=(2-x²)²-4ab. Pressure-release water contributes
an impedance rho_w omega² tan(qh)/q, giving

```text
R + (rho_w/rho_s) x^4 a tan(qh)/sqrt(c²/cw²-1) = 0
```

The implementation uses sin/cos above cw to remove tangent poles, and tanh
below cw to avoid growing hyperbolic factors. It brackets the lowest trapped
root for simple verification models. Its finite grid can miss extremely close
roots; this is a limited verification tool, not a general production root finder.
The dry reference uses the separate cubic Rayleigh equation.

## Monitoring and uncertainty

Coefficients a_j are delta log Vs in an existing physical region or whole
half-space. Vp and density are fixed; interface positions are fixed. Water
nuisance coefficients are delta log h, delta log cw, and delta log rho_w.

G_ij = d log c0(fi)/d a_j is computed by centered re-solves at step h and h/2.
The returned derivative uses h/2; their discrepancy is recorded and checked.
Perturbed/reference fundamental branches must agree by overlap and pass QC.
The empirical phase error estimate must be below 0.1 observation sigma.
This threshold is a provisional numerical budget, not an achieved guarantee for
arbitrary backgrounds. Direct re-solves test the linearization range.

For data d=G a+H b+epsilon, Gaussian zero-mean water prior Cb yields
Ce=Cd+H Cb H^T. Conditioning uses Cholesky solves with G Ca G^T+Ce.
Data covariance must agree with the observation's marginal standard deviations;
correlations remain explicit. Reports include shear intervals, conditional water
mean, posterior predicted phase changes, averaging matrix, singular values of
Ce^-1/2 G Ca^1/2, and a prior-dominance flag when posterior/prior variance >0.8.
That flag threshold is descriptive and provisional. Baseline uncertainty,
nonlinear inversion, and temporal regularization are later work.
