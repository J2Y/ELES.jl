module ELES

using LinearAlgebra
using Zygote
using Optim
using ChainRulesCore: ignore_derivatives

export LayerModel, mixed_model, viscosity, solve_layered, field, flux,
       ekman_forward, coupled_forward, interface_residuals, fit_viscosity,
       series_tail, solve_poly_layers, layer_radius, LAYER_KINDS

"""
Admissible layer kinds, matching Table 1 of the manuscript: a layer viscosity
that is a polynomial of degree at most two, or an exponential.

* `:constant` --- one value per layer; may jump at an interface.
* `:linear` --- one value per knot; affine on each layer.
* `:quadratic` --- the **general** quadratic. One value per knot plus one per
  layer midpoint, giving three coefficients per layer, which is as far as the
  classical closed forms reach. Positivity is not automatic here and is
  checked layer by layer.
* `:exponential` --- one value per knot, interpolating `log K` linearly.

These are the *interpolating* forms: the parameters are viscosity values at
fixed nodes, which is what an inversion varies. A profile that is already
known analytically, such as the K-profile cubic of the manuscript's Example 2,
is instead passed coefficient-wise to [`solve_poly_layers`]. The node-based
`LayerModel` API has no cubic kind; the coefficient-based entry point accepts
polynomials of any degree, subject to positivity and Taylor convergence.
"""
const LAYER_KINDS = (:constant, :linear, :quadratic, :exponential)

"""
    LayerModel(interfaces; kind=:linear)

Fixed layers with positive independent viscosity parameters. `:constant` uses
one value per layer (jumps allowed); `:linear` and `:exponential` use one value
per knot; `:quadratic` uses one per knot and per layer midpoint (`2N+1` for
`N` layers). Use `mixed_model` for continuous mixtures of kinds.
"""
struct LayerModel
    interfaces::Vector{Float64}
    kinds::Vector{Symbol}
    left_index::Vector{Int}
    right_index::Vector{Int}
    mid_index::Vector{Int}      # 0 where the layer has no midpoint parameter
    nparameters::Int
end

function _interfaces(z)
    edges = Float64.(collect(z))
    length(edges) >= 2 || throw(ArgumentError("at least two interfaces are required"))
    all(isfinite, edges) && all(diff(edges) .> 0) ||
        throw(ArgumentError("interfaces must be finite and strictly increasing"))
    edges
end

function LayerModel(z; kind::Symbol=:linear)
    kind in LAYER_KINDS || throw(ArgumentError("unknown layer kind"))
    edges = _interfaces(z); n = length(edges)-1
    kind == :constant &&
        return LayerModel(edges,fill(kind,n),collect(1:n),collect(1:n),zeros(Int,n),n)
    mid = kind == :quadratic ? collect(n+2:2n+1) : zeros(Int,n)
    LayerModel(edges,fill(kind,n),collect(1:n),collect(2:n+1),mid,
               kind == :quadratic ? 2n+1 : n+1)
end

"""
    mixed_model(interfaces, kinds)

Continuous mixed profile. Adjacent nodes of each constant layer share a
parameter automatically. `left_index`/`right_index` document the mapping
from independent parameters to layer endpoints; no hard-coded layer count.
"""
function mixed_model(z, kinds)
    edges = _interfaces(z); types = Symbol.(kinds); n = length(edges)-1
    length(types) == n || throw(DimensionMismatch("one kind is required per layer"))
    all(t -> t in LAYER_KINDS, types) || throw(ArgumentError("unknown layer kind"))
    group = ones(Int,n+1)
    for k in 1:n
        group[k+1] = group[k]+(types[k] == :constant ? 0 : 1)
    end
    nknot = group[end]; mid = zeros(Int,n); extra = 0
    for k in 1:n
        if types[k] == :quadratic
            extra += 1; mid[k] = nknot+extra
        end
    end
    LayerModel(edges,types,group[1:end-1],group[2:end],mid,nknot+extra)
end

function _values(model,K)
    length(K) == model.nparameters || throw(DimensionMismatch("wrong number of viscosity values"))
    all(x -> isfinite(x) && x > 0,K) || throw(ArgumentError("viscosity values must be finite and positive"))
    nothing
end

function _layer(model,z)
    model.interfaces[1] <= z <= model.interfaces[end] ||
        throw(ArgumentError("evaluation height is outside the model"))
    clamp(searchsortedlast(model.interfaces,z),1,length(model.kinds))
end

# Fixed truncation of the exponential viscosity series. A fixed length keeps
# the construction free of data-dependent control flow, so reverse-mode AD
# traces one path; `_kcoef` rejects contrasts for which it is not enough.
const _EXP_TERMS = 24
const _FACT = [i == 0 ? 1.0 : prod(1.0:i) for i in 0:_EXP_TERMS]

"""
    _kcoef(model, K, k)

Coefficients `[k_0, k_1, ...]` of the layer viscosity in the midpoint-normalised
local coordinate `x` in `[-1,1]`, so `K(x) = sum_i k_i x^i`. Polynomial kinds
give an exact finite vector; `:exponential` uses degree 24 with a contrast
check on the retained viscosity coefficient.
"""
function _kcoef(model,K,k)
    kind = model.kinds[k]
    left = K[model.left_index[k]]; right = K[model.right_index[k]]
    kind == :constant && return [left]
    if kind == :exponential
        alpha = (log(left)+log(right))/2
        beta  = (log(right)-log(left))/2
        abs(beta)^_EXP_TERMS/_FACT[end] < 1e-14 ||
            throw(ArgumentError("exponential layer contrast is too large for the truncated series"))
        return [exp(alpha)*beta^i/_FACT[i+1] for i in 0:_EXP_TERMS]
    end
    kind == :linear && return [(left+right)/2,(right-left)/2]
    # General quadratic through the two knot values and the layer midpoint,
    # in the local coordinate x in [-1,1]: K(-1)=left, K(0)=mid, K(1)=right.
    mid = K[model.mid_index[k]]
    c = [mid,(right-left)/2,(right+left)/2-mid]
    ignore_derivatives() do
        _admissible(c,k)
    end
    c
end

"""
    layer_radius(kcoef)

Distance from the layer midpoint to the nearest zero of `K` in the **complex**
plane, measured in the local coordinate `x`. The Taylor recurrence converges
for `|x| < layer_radius`, so a layer needs `layer_radius > 1`.

Requiring only that `K > 0` on the real interval is *not* enough. A general
quadratic with real coefficients can have a complex-conjugate pair of zeros,
whose modulus is `sqrt(k0/k2)`; that pair limits the radius while leaving `K`
positive everywhere real. For instance knot values `0.5, 0.5` with midpoint
`0.07` keep `K >= 0.07` yet put the zeros at `0.40`, and the series then
diverges on the layer.
"""
function layer_radius(kcoef)
    degree = something(findlast(!iszero, kcoef), 1) - 1
    degree < 1 && return Inf
    leading = kcoef[degree+1]
    companion = zeros(ComplexF64, degree, degree)
    for i in 2:degree
        companion[i, i-1] = 1
    end
    for i in 1:degree
        companion[i, degree] = -kcoef[i]/leading
    end
    minimum(abs, eigvals(companion))
end

"""True when `sum kcoef[i+1] x^i > 0` for every `x` in `[-1,1]`."""
function _positive_on_layer(c)
    _Kat(c, -1.0) > 0 && _Kat(c, 0.0) > 0 && _Kat(c, 1.0) > 0 || return false
    # A real zero inside the layer is exactly a radius below one with a real
    # minimiser, so the radius test below covers the interior; sample anyway.
    all(x -> _Kat(c, x) > 0, range(-1, 1; length = 65))
end

# The check inspects the coefficients but contributes nothing to the
# derivative, and `eigvals` mutates, so it is kept off the AD tape with
# `ignore_derivatives`. The primal still runs, so an inadmissible layer still
# throws during a gradient evaluation.
"""Reject a layer on which the series cannot converge, or `K` changes sign."""
function _admissible(c, k)
    _positive_on_layer(c) ||
        throw(ArgumentError("layer $k: K is not positive on the layer"))
    radius = layer_radius(c)
    radius > 1 || throw(ArgumentError(
        "layer $k: the nearest zero of K in the complex plane lies at " *
        "|x| = $(round(radius; digits=4)) <= 1, so the Taylor series does not " *
        "converge on this layer; use a finer partition"))
    nothing
end

_Kat(kcoef,x) = begin           # Horner, differentiable
    acc = kcoef[end]
    for i in length(kcoef)-1:-1:1
        acc = acc*x+kcoef[i]
    end
    acc
end

function _local(model,z,k)
    lo = model.interfaces[k]; hi = model.interfaces[k+1]
    half = (hi-lo)/2
    return (z-(lo+hi)/2)/half,half
end

"""Evaluate K(z); values at a jump use the right layer, except the top endpoint."""
function viscosity(model::LayerModel,K,z::Real)
    _values(model,K); k = _layer(model,z); t,_ = _local(model,z,k)
    _Kat(_kcoef(model,K,k),t)
end
viscosity(model::LayerModel,K,z::AbstractVector) = map(x -> viscosity(model,K,x),z)

"""Viscosity actually used by a solution; works for `solve_poly_layers` too."""
function viscosity(sol,z::Real)
    k = _layer(sol.model,z); t,_ = _local(sol.model,z,k)
    _Kat(sol.kcoefs[k],t)
end
viscosity(sol,z::AbstractVector) = map(x -> viscosity(sol,x),z)

# Midpoint-normalised analytic fundamental solutions for an arbitrary
# power-series viscosity. Collecting the coefficient of x^n in
# (K phi')' = gamma*half^2*phi gives
#   (n+1) * sum_i k_i (n-i+2) a_{n-i+2} = gamma*half^2 a_n,
# which is the recurrence below; at degree two it reduces to the
# squared-affine/affine recurrence used previously. The finite series is a
# numerical evaluation of the exact layer ODE, not an exact-arithmetic claim.
function _fundamental(kcoef,gamma,half,order)
    d = length(kcoef)-1
    g = gamma*half^2
    k0 = kcoef[1]
    coeff = Zygote.Buffer(zeros(ComplexF64,order+1,2))
    for i in 1:order+1; coeff[i,1] = 0; coeff[i,2] = 0; end
    coeff[1,1] = 1; coeff[2,2] = 1
    for n in 0:order-2
        acc1 = zero(ComplexF64); acc2 = zero(ComplexF64)
        for i in 1:d
            j = n-i+3                      # holds a_{n-i+2}
            j >= 1 || continue
            w = kcoef[i+1]*(n-i+2)
            acc1 += w*coeff[j,1]; acc2 += w*coeff[j,2]
        end
        den = k0*(n+1)*(n+2)
        coeff[n+3,1] = (g*coeff[n+1,1]-(n+1)*acc1)/den
        coeff[n+3,2] = (g*coeff[n+1,2]-(n+1)*acc2)/den
    end
    copy(coeff)
end

function _basis(coeff,t,half)
    order = size(coeff,1)-1
    value = transpose(coeff)*(t .^ (0:order))
    deriv = transpose(coeff[2:end,:])*((1:order).*t.^(0:order-1))/half
    value[1],value[2],deriv[1],deriv[2]
end

struct LayerSolution
    model::LayerModel
    values::Vector{Float64}
    kcoefs::Vector{Vector{Float64}}
    gamma::ComplexF64
    bases::Vector{Matrix{ComplexF64}}
    coefficients::Vector{ComplexF64}
    order::Int
end

"""
    series_tail(sol)

Largest relative coefficient magnitude in the last `max(2, degree(K)+1)`
retained orders, normalized separately for each fundamental solution. The
degree-sized window catches sparse coefficient sequences that a two-term
check can miss. This is a diagnostic, not a bound on the omitted remainder
or the velocity error: compare orders or refine the partition as well.
"""
function series_tail(sol::LayerSolution)
    worst = 0.0
    for (basis,kcoef) in zip(sol.bases,sol.kcoefs)
        degree = something(findlast(!iszero,kcoef),1)-1
        first = max(1,size(basis,1)-max(2,degree+1)+1)
        for j in axes(basis,2)
            scale = maximum(abs,view(basis,:,j))
            scale > 0 || continue
            tail = maximum(abs,view(basis,first:size(basis,1),j))
            worst = max(worst,tail/scale)
        end
    end
    worst
end

"""
    solve_layered(model, values, gamma, Wb, Wt; order=160)

Solve (K W')'=gamma*W with endpoint values Wb,Wt. The matching system
preserves both W and K W' at interfaces, including jumps in K. The
homogeneous, fixed-interface solver is differentiable in positive values.
The implemented numerical evaluation uses an order-`order` Taylor series;
check `series_tail` and compare orders or partitions before using extreme viscosity
contrasts or layer widths.
"""
function solve_layered(model::LayerModel,K,gamma,Wb,Wt; order::Int=160)
    _values(model,K)
    _assemble(model,collect(K),map(k -> _kcoef(model,K,k),1:length(model.kinds)),
              gamma,Wb,Wt,order)
end

"""
    solve_poly_layers(interfaces, kcoefs, gamma, Wb, Wt; order=160)

General-degree entry point. `kcoefs[k]` lists `[k_0, k_1, ...]` of
`K(x) = sum_i k_i x^i` on layer `k`, in the midpoint-normalised coordinate
`x` in `[-1,1]`. Unlike `solve_layered`, the viscosity is given directly
rather than through knot values and a layer kind, so any polynomial (or
truncated power series) is admissible, including general quadratics and
cubics that are not perfect powers. `K > 0` on each closed layer is the
implementation checks positivity and rejects complex roots at radius at most
one; finite-order accuracy still needs an order or partition comparison.
"""
function solve_poly_layers(interfaces,kcoefs,gamma,Wb,Wt; order::Int=160)
    edges = _interfaces(interfaces); n = length(edges)-1
    length(kcoefs) == n || throw(DimensionMismatch("one coefficient vector is required per layer"))
    coefs = [collect(float.(c)) for c in kcoefs]
    for (k,c) in enumerate(coefs)
        isempty(c) && throw(ArgumentError("a layer needs at least a constant term"))
        all(isfinite,c) || throw(ArgumentError("viscosity coefficients must be finite"))
        ignore_derivatives() do
            _admissible(c,k)
        end
    end
    model = LayerModel(edges,fill(:given,n),collect(1:n),collect(1:n),zeros(Int,n),n)
    _assemble(model,[_Kat(c,0.0) for c in coefs],coefs,gamma,Wb,Wt,order)
end

function _assemble(model::LayerModel,values,kcoefs,gamma,Wb,Wt,order::Int)
    order >= 4 || throw(ArgumentError("series order must be at least four"))
    n = length(model.kinds)
    bases = map(1:n) do k
        half = (model.interfaces[k+1]-model.interfaces[k])/2
        _fundamental(kcoefs[k],gamma,half,order)
    end
    basis(k,z) = ((t,h) = _local(model,z,k); _basis(bases[k],t,h))
    layerK(k,z) = ((t,_) = _local(model,z,k); _Kat(kcoefs[k],t))
    M = Zygote.Buffer(zeros(ComplexF64,2n,2n))
    rhs = Zygote.Buffer(zeros(ComplexF64,2n))
    for i in 1:2n
        rhs[i] = 0
        for j in 1:2n; M[i,j] = 0; end
    end
    v1,v2,_,_ = basis(1,model.interfaces[1])
    M[1,1] = v1; M[1,2] = v2; rhs[1] = Wb
    for k in 1:n-1
        z = model.interfaces[k+1]
        l1,l2,dl1,dl2 = basis(k,z)
        r1,r2,dr1,dr2 = basis(k+1,z)
        KL,KR = layerK(k,z),layerK(k+1,z)
        M[2k,2k-1] = l1; M[2k,2k] = l2
        M[2k,2k+1] = -r1; M[2k,2k+2] = -r2
        M[2k+1,2k-1] = KL*dl1; M[2k+1,2k] = KL*dl2
        M[2k+1,2k+1] = -KR*dr1; M[2k+1,2k+2] = -KR*dr2
    end
    v1,v2,_,_ = basis(n,model.interfaces[end])
    M[2n,2n-1] = v1; M[2n,2n] = v2; rhs[2n] = Wt
    coefficients = copy(M) \ copy(rhs)
    LayerSolution(model,collect(values),[collect(float.(c)) for c in kcoefs],
                  ComplexF64(gamma),bases,coefficients,order)
end

function _state(sol,z,k)
    t,h = _local(sol.model,z,k)
    v1,v2,d1,d2 = _basis(sol.bases[k],t,h)
    c1,c2 = sol.coefficients[2k-1],sol.coefficients[2k]
    (value=c1*v1+c2*v2,flux=_Kat(sol.kcoefs[k],t)*(c1*d1+c2*d2))
end
field(sol::LayerSolution,z::Real) = _state(sol,z,_layer(sol.model,z)).value
field(sol::LayerSolution,z::AbstractVector) = map(x -> field(sol,x),z)
flux(sol::LayerSolution,z::Real) = _state(sol,z,_layer(sol.model,z)).flux
flux(sol::LayerSolution,z::AbstractVector) = map(x -> flux(sol,x),z)

"""One-sided maximum interface errors in W and K W'."""
function interface_residuals(sol::LayerSolution)
    errors = map(1:length(sol.model.kinds)-1) do k
        z = sol.model.interfaces[k+1]
        l = _state(sol,z,k); r = _state(sol,z,k+1)
        (abs(l.value-r.value),abs(l.flux-r.flux))
    end
    isempty(errors) && return (value=0.0,flux=0.0)
    (value=maximum(first.(errors)),flux=maximum(last.(errors)))
end

"""Ekman complex velocity defect W=u+i*v; add the geostrophic velocity separately."""
function ekman_forward(model,K,z; f=1.0,lower=-1im,upper=0im,order=160)
    field(solve_layered(model,K,im*f,lower,upper;order),z)
end

"""
    coupled_forward(model,K,z; a,b,lower=(0.,-5.),upper=(0.,0.),order=160)

Solve (K u')'=-a*v, (K v')'=b*u with a,b>0, using
W=u+i*sqrt(a/b)*v and gamma=i*sqrt(a*b). For katabatic flow v is the
potential-temperature perturbation. Equal momentum and heat diffusivities
(unit turbulent Prandtl number) are required for this reduction.
"""
function coupled_forward(model,K,z; a,b,lower=(0.0,-5.0),upper=(0.0,0.0),order=160)
    a > 0 && b > 0 || throw(ArgumentError("a and b must be positive"))
    r = sqrt(a/b)
    sol = solve_layered(model,K,im*sqrt(a*b),lower[1]+im*r*lower[2],upper[1]+im*r*upper[2];order)
    W = field(sol,z)
    (u=real.(W),v=imag.(W)/r,solution=sol)
end

"""
    fit_viscosity(model,z,observed; initial,lower=0.06,upper=0.50,...)

Box-constrained AD inversion from complex field observations, using log
viscosity parameters and Zygote/Fminbox(L-BFGS). No true K values enter
this loss. Bounds may be scalar or one value per independent parameter.
An inadmissible initial profile throws; inadmissible line-search trials
return an infinite loss and are rejected. The returned `diagnostics` reports
these rejections and the final Taylor tail. `converged` is false if that tail
exceeds `tail_tol`; refine the partition rather than disabling validation.
"""
function fit_viscosity(model,z,observed; initial,lower=0.06,upper=0.50,
                       gamma=1im,bottom=-1im,top=0im,order=160,
                       iterations=1000,outer_iterations=20,g_tol=1e-11,
                       tail_tol=1e-10)
    length(z) == length(observed) && !isempty(z) || throw(DimensionMismatch("data sizes differ or are empty"))
    _values(model,initial)
    lo = lower isa Real ? fill(Float64(lower),model.nparameters) : collect(lower)
    hi = upper isa Real ? fill(Float64(upper),model.nparameters) : collect(upper)
    length(lo) == length(hi) == model.nparameters || throw(DimensionMismatch("bound sizes differ"))
    all(0 .< lo .< initial .< hi) || throw(ArgumentError("initial values must lie strictly inside positive bounds"))
    isfinite(tail_tol) && tail_tol > 0 || throw(ArgumentError("tail_tol must be finite and positive"))
    # Validate the initial model once, outside AD and outside trial rejection.
    solve_layered(model,initial,gamma,bottom,top;order)
    function admissible(th)
        values = exp.(th)
        all(x -> isfinite(x) && x > 0,values) || return false
        try
            for k in eachindex(model.kinds)
                _kcoef(model,values,k)
            end
        catch err
            err isa ArgumentError || rethrow()
            return false
        end
        true
    end
    forward(th) = field(solve_layered(model,exp.(th),gamma,bottom,top;order),z)
    rawloss(th) = sum(abs2,forward(th).-observed)/length(z)
    rejected = Ref(0)
    function loss(th)
        if !admissible(th)
            rejected[] += 1
            return Inf
        end
        rawloss(th)
    end
    function gradient!(G,th)
        if !admissible(th)
            # Hager-Zhang contracts trials with nonfinite value/slope.
            fill!(G,NaN)
        else
            G .= Zygote.gradient(rawloss,th)[1]
        end
        G
    end
    result = Optim.optimize(loss,gradient!,log.(lo),log.(hi),log.(initial),
                            Optim.Fminbox(Optim.LBFGS()),
                            Optim.Options(;iterations,outer_iterations,g_tol))
    th = Optim.minimizer(result)
    finaltail = series_tail(solve_layered(model,exp.(th),gamma,bottom,top;order))
    (values=exp.(th),prediction=forward(th),loss=Optim.minimum(result),
     converged=Optim.converged(result) && finaltail <= tail_tol,optimizer=result,
     diagnostics=(admissibility_rejections=rejected[],series_tail=finaltail,
                  refinement_required=finaltail > tail_tol))
end

end # module
