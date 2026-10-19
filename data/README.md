# Synthetic velocity observations

`synthetic_velocity.csv` contains the manuscript's 161 observation heights,
true viscosity, clean complex-velocity components, and noisy components:

```
z,K_true,u_true,v_true,u_obs,v_obs
```

The true viscosity is `K(s)=0.08+1.6*s*(1-s)^2` on `[0,1]`.
Clean velocities were generated independently by SciPy `solve_bvp` in flux form,
with tolerance `1e-10`, gamma `i`, and endpoint values `-i` and `0`.
Independent NumPy normal real/imaginary samples (seed `20261007`) were rescaled
to exactly 1% relative complex-velocity L2 norm.

Only `z`, `u_obs`, and `v_obs` enter the inverse loss. `K_true` evaluates recovery;
clean velocities determine the known synthetic noise scale for the optional
working-likelihood uncertainty calculation. No experimental/LES data are used.

`fitted_parameters.csv` archives the eleven fitted viscosity parameters per
model used in the manuscript (constant/linear/quadratic). The columns are
`model,index,node,value`, in the solver's parameter order. It is a reproducibility
checkpoint, not input to the inverse loss. `scripts/reproduce_paper.jl` evaluates
these values by default; `--refit` recomputes them from the noisy velocities.

`figure1_profile.csv` specifies the manuscript's smooth forward-verification
profile. Each row provides a layer's endpoints, degree, three viscosity values
and endpoint slopes. The degree mix is `(2,1,2,1,2,1,1,2,0,0)`, with 13 nodal
parameters. Viscosity and its first derivative are continuous at every interface.
The idealized profile has a single maximum at `z=0.305`, remains positive, and
approaches a small constant background in the upper two layers.
`examples/forward_verification.jl` loads this same specification for notebook 01.
