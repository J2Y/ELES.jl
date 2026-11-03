using LinearAlgebra
"""
Independent fixed-step RK4 shooting reference for (K W')'=gamma*W.
Integrates the first-order flux system rather than analytic layer bases.
Every interface is an integration boundary, including discontinuous K.
"""
function reference_field(model,K,gamma,bottom,top,zgrid; steps_per_layer=600)
    P = Matrix{ComplexF64}(I,2,2)
    snapshots = Dict{Float64,Matrix{ComplexF64}}()
    snapshots[round(model.interfaces[1];digits=12)] = copy(P)
    for k in eachindex(model.kinds)
        lo,hi = model.interfaces[k:k+1]
        zl = zgrid[(lo .<= zgrid) .& (zgrid .<= hi)]
        mesh = sort(unique([lo; collect(range(lo,hi;length=steps_per_layer+1)); zl; hi]))
        left,right = K[model.left_index[k]],K[model.right_index[k]]
        middle = model.mid_index[k] == 0 ? nothing : K[model.mid_index[k]]
        # The interpolation rule is restated here from its definition rather
        # than taken from the solver, so the RK4 integration stays an
        # independent check of the analytic layer bases.
        function Klocal(z)
            s = (z-lo)/(hi-lo)
            kind = model.kinds[k]
            kind == :constant && return left
            kind == :exponential && return exp((1-s)*log(left)+s*log(right))
            kind == :linear && return (1-s)*left+s*right
            kind == :quadratic || throw(ArgumentError("reference does not know this layer kind"))
            # General quadratic through left, middle, right at s = 0, 1/2, 1.
            x = 2s-1
            middle + (right-left)/2*x + ((right+left)/2-middle)*x^2
        end
        rhs(z,Q) = [Q[2,1]/Klocal(z) Q[2,2]/Klocal(z); gamma*Q[1,1] gamma*Q[1,2]]
        for j in 1:length(mesh)-1
            z,h = mesh[j],mesh[j+1]-mesh[j]
            a = rhs(z,P); b = rhs(z+h/2,P+h*a/2)
            c = rhs(z+h/2,P+h*b/2); d = rhs(z+h,P+h*c)
            P = P+h*(a+2b+2c+d)/6
            snapshots[round(mesh[j+1];digits=12)] = copy(P)
        end
    end
    qbottom = (top-P[1,1]*bottom)/P[1,2]
    map(zgrid) do z
        Q = snapshots[round(z;digits=12)]
        Q[1,1]*bottom+Q[1,2]*qbottom
    end
end

"""
    reference_arbitrary(K, gamma, bottom, top, zgrid; steps=20000)

Fixed-step RK4 shooting reference for `(K W')' = gamma W` with an arbitrary
positive `K(z)` given as a function, for profiles that no `LayerModel`
describes. Like `reference_field` it integrates the first-order flux system
rather than the analytic layer bases, so it is an independent check.
"""
function reference_arbitrary(Kfun, gamma, bottom, top, zgrid; steps = 20000)
    mesh = sort(unique([collect(range(zgrid[1], zgrid[end]; length = steps + 1)); zgrid]))
    P = Matrix{ComplexF64}(I, 2, 2)
    snapshots = Dict{Float64,Matrix{ComplexF64}}()
    snapshots[round(mesh[1]; digits = 12)] = copy(P)
    rhs(z, Q) = [Q[2, 1]/Kfun(z) Q[2, 2]/Kfun(z); gamma*Q[1, 1] gamma*Q[1, 2]]
    for j in 1:length(mesh)-1
        z, h = mesh[j], mesh[j+1] - mesh[j]
        a = rhs(z, P);               b = rhs(z + h/2, P + h*a/2)
        c = rhs(z + h/2, P + h*b/2); d = rhs(z + h, P + h*c)
        P = P + h*(a + 2b + 2c + d)/6
        snapshots[round(mesh[j+1]; digits = 12)] = copy(P)
    end
    q = (top - P[1, 1]*bottom)/P[1, 2]
    map(z -> (Q = snapshots[round(z; digits = 12)]; Q[1, 1]*bottom + Q[1, 2]*q), zgrid)
end
