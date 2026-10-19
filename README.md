# ELES.jl — exact layered Ekman solver

Exact solutions of the steady, horizontally homogeneous boundary-layer system

```
d/dz ( K(z) du/dz ) = M u,      M = [0 -a; b 0],   a, b > 0
```

on a partition of `[z_b, z_t]`, with the layers matched by continuity of the
velocity and of the stress `K du/dz`. Diagonalising `M` reduces each layer to
a scalar equation `(K φ')' = γ φ` with purely imaginary `γ`, which has a
closed-form fundamental pair when the layer viscosity is a polynomial of
degree at most two, the cubic with one double and one simple root, or an
exponential. An energy identity makes the assembled matching system
nonsingular for any finite partition.

This package accompanies **A differentiable exact solver for layered Ekman flows
with piecewise-continuous eddy viscosity**, a manuscript submitted to
*Applied Mathematics Letters* (AML), and reproduces its three examples.

## Install and test

```sh
git clone https://github.com/J2Y/ELES.jl.git
cd ELES.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The bundled environment was tested with Julia 1.12.6.

## Examples

The examples are Julia notebooks under `notebooks/`, one per manuscript
example, and they are stored with their output so they can be read without
being run:

| notebook | manuscript |
|---|---|
| `01_catalogue_and_verification.ipynb` | Example 1 — the catalogue, Abel's identity, a ten-layer mixed-degree profile against an independent Runge–Kutta reference |
| `02_operational_kprofile.ipynb` | Example 2 — `K = κ u∗ z (1 − z/h)²`, the profile of the Troen–Mahrt and K-profile schemes |
| `03_viscosity_inversion.ipynb` | Example 3 — inversion at interpolation degrees 0, 1, 2 with a matched parameter budget |
| `04_katabatic_slope_flow.ipynb` | the slope-flow remark in Section 2 |

Figure 1 uses the continuously differentiable, single-peaked profile specified
in `data/figure1_profile.csv`, loaded by `examples/forward_verification.jl`.
It retains four quadratic, four affine and two constant layers with 13 parameters.

To run them, register the kernel once and open Jupyter:

```sh
julia --project=. scripts/setup_notebook.jl
julia --project=. scripts/start_notebooks.jl
```

Each notebook activates this project in its first cell, so a supported Julia kernel
can be used. They can also be re-executed in place:

```sh
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=eles notebooks/*.ipynb
```

## Reproducing the paper tables and uncertainty analysis

```sh
julia --project=. scripts/reproduce_paper.jl
python3 -m pip install -r reproduction/requirements.txt
python3 reproduction/check_closed_forms.py
```

Results are written to `results/paper/`. The Julia command reproduces Table 2,
the fourteen one-/four-layer catalogue checks plus quadratic root cases, and
Table 3/Figure 2 with 95% local Gaussian intervals. Velocity errors in Table 2
are maxima over **201 uniformly spaced evaluation points**, not certified
continuous suprema. The independent mpmath check compares the exported Julia
velocities with 40-digit Bessel, Legendre, power-law and hypergeometric solutions.
Only that optional check needs Python; the solver and plotting remain Julia.

The default reevaluates the archived fit parameters in `data/`; use
`julia --project=. scripts/reproduce_paper.jl --refit` to optimize from the
symmetric initial profile again. Uncertainty uses a fixed Julia RNG seed;
Monte Carlo endpoints can differ slightly from the earlier NumPy-generated
figure. The likelihood approximation, noise scale and parameter budget agree.
Historical MATLAB `bvp4c` and SciPy comparisons are not rerun by these commands.
Changes and measured validation results are recorded in [PROGRESS.md](PROGRESS.md).

## Using the solver

Two entry points, for two different situations.

**Interpolating forms.** `LayerModel(interfaces; kind)` parameterises the
viscosity by its values at fixed nodes, which is what an inversion varies:

| `kind` | degree | parameters | notes |
|---|---|---|---|
| `:constant` | 0 | one per layer | may jump at an interface |
| `:linear` | 1 | one per knot | continuous |
| `:quadratic` | 2 | one per knot **and** one per layer midpoint | the general quadratic |
| `:exponential` | — | one per knot | interpolates `log K` |

`mixed_model(interfaces, kinds)` mixes them; adjacent nodes of a constant
layer share a parameter automatically, so the profile stays continuous.

```julia
using ELES
edges = collect(range(0, 1, length = 11))
model = LayerModel(edges; kind = :linear)
z = collect(range(0, 1, length = 101))
W = ekman_forward(model, collect(range(0.08, 0.2, length = 11)), z; f = 1.0)
```

**A profile already known analytically** is passed coefficient-wise to
`solve_poly_layers`; the node-based API has no cubic kind. `kcoefs[k]` lists `[k_0, k_1, …]` of `K(x) = Σ k_i x^i` on layer `k`,
in the midpoint-normalised coordinate `x ∈ [-1, 1]`. Any polynomial degree is
admissible there; see `notebooks/02_operational_kprofile.ipynb`.

Other exports: `field`, `flux`, `viscosity`, `interface_residuals`,
`coupled_forward` (slope flow, with `v` the temperature perturbation),
`fit_viscosity` (box-constrained inversion by Zygote and L-BFGS),
`series_tail`, `layer_radius`, `LAYER_KINDS`.

## Two conditions a layer must satisfy

`K ≥ κ > 0` on the layer, which the existence result needs, and a Taylor
series that converges on the layer, which the numerical evaluation needs.
**The second does not follow from the first.** The series converges for
`|x| < layer_radius(kcoef)`, the distance from the layer midpoint to the
nearest zero of `K` in the *complex* plane. A general quadratic with real
coefficients can carry a conjugate pair of zeros of modulus `sqrt(k0/k2)`,
which may be below one while `K` stays positive everywhere real:

```julia
julia> layer_radius([0.07, 0.0, 0.43])   # knots 0.5, 0.5 and midpoint 0.07
0.4035...                                # K ≥ 0.07 on the layer, yet < 1
```

Such a layer is refused with an explanatory error. For a radius only slightly
above one the series converges but slowly, and `series_tail(solution)` reports
the largest coefficient in the last `max(2, degree(K)+1)` retained orders,
normalized per fundamental solution. This catches sparse coefficient sequences
that a two-term tail check misses. It is not a remainder or velocity-error
bound: compare orders or partitions as well. A large diagnostic calls for
refinement. An affine layer always has radius
`(K_k + K_{k-1})/|K_k − K_{k-1}| > 1`; an exponential layer is entire.

## Differentiability

With the interfaces fixed, the solution is continuously differentiable in the
viscosity parameters. Each Taylor coefficient is a rational function of the
layer coefficients with `k_0 = K(midpoint) > 0`, so the dependence stays
analytic through the zero-slope limit at which the closed-form Bessel and
power-law bases degenerate; no branch switching is needed and reverse-mode
automatic differentiation applies directly. The admissibility check is kept
off the AD tape, so it still throws during a gradient evaluation without
interfering with valid derivatives. `fit_viscosity` rejects inadmissible
line-search trials rather than aborting; invalid initial profiles still throw.
Its `diagnostics` reports the rejection count and final Taylor tail. If
`refinement_required` is true, `converged` is false and the partition or order
needs refinement.

## Layout

```
src/ELES.jl                  the solver
examples/reference_shooting.jl  independent RK4 shooting references
notebooks/                   the four examples, with stored output
test/runtests.jl             the test suite
data/synthetic_velocity.csv  observations used by Example 3
scripts/                     Jupyter kernel registration and launch
```

`examples/reference_shooting.jl` integrates the first-order flux system rather
than the analytic layer bases, and restates the interpolation rules from their
definitions rather than calling the solver, so the tests and notebooks compare
against something independent.

## Licence and citation

MIT, see `LICENSE`. Citation metadata in `CITATION.cff`.
