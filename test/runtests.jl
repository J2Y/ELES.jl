using Test, ELES, LinearAlgebra, Zygote
using ELES: _Kat
include(joinpath(@__DIR__, "..", "examples", "reference_shooting.jl"))

"""Parameter vector making every layer interpolate `f` at its own nodes."""
function nodal(model, f)
    values = zeros(model.nparameters)
    for k in eachindex(model.kinds)
        lo, hi = model.interfaces[k], model.interfaces[k+1]
        if model.kinds[k] == :constant
            values[model.left_index[k]] = f((lo + hi)/2)
        else
            values[model.left_index[k]] = f(lo)
            values[model.right_index[k]] = f(hi)
            model.mid_index[k] == 0 || (values[model.mid_index[k]] = f((lo + hi)/2))
        end
    end
    values
end

function central(loss, th; h = 1e-5)
    Iden = Matrix{Float64}(I, length(th), length(th))
    [(loss(th + h*Iden[:, j]) - loss(th - h*Iden[:, j]))/(2h) for j in eachindex(th)]
end

@testset "Closed-form constant Ekman solution and flat-slope AD" begin
    edges = collect(range(0., 1.; length = 11)); z = collect(range(0., 1.; length = 31))
    gamma = 1im; q = sqrt(gamma/.1)
    exact = -1im .* sinh.(q .* (1 .- z)) ./ sinh(q)
    for kind in LAYER_KINDS
        model = LayerModel(edges; kind)
        K = fill(.1, model.nparameters)
        @test maximum(abs.(ekman_forward(model, K, z) .- exact)) < 1e-12
        # A flat initial profile must still have valid parameter derivatives.
        loss(th) = sum(abs2, ekman_forward(model, exp.(th), z) .- exact .* .9)
        th = log.(K); ad = Zygote.gradient(loss, th)[1]
        @test norm(ad - central(loss, th))/norm(ad) < 1e-6
    end
end

@testset "Interface flux and independent RK4 reference" begin
    edges = collect(range(0., 1.; length = 11)); z = collect(range(0., 1.; length = 51))
    for kind in LAYER_KINDS
        model = LayerModel(edges; kind)
        K = nodal(model, z -> .08 + .12z + .05*sinpi(z))
        sol = solve_layered(model, K, 1im, -1im, 0im)
        jumps = interface_residuals(sol)
        @test max(jumps.value, jumps.flux) < 1e-12
        @test abs(field(sol, 0.) + 1im) < 1e-12
        @test abs(field(sol, 1.)) < 1e-12
        @test maximum(abs.(field(sol, z) .-
                           reference_field(model, K, 1im, -1im, 0im, z))) < 1e-8
        @test series_tail(sol) < 1e-14
    end
end

@testset "Abel's identity on every layer basis" begin
    edges = collect(range(0., 1.; length = 6))
    kinds = [:linear, :quadratic, :constant, :exponential, :quadratic]
    model = mixed_model(edges, kinds)
    K = nodal(model, z -> .07 + .3z*(1 - z) + .05z)
    sol = solve_layered(model, K, 1im, -1im, 0im)
    for k in eachindex(kinds)
        lo, hi = model.interfaces[k], model.interfaces[k+1]
        half = (hi - lo)/2
        coefficients = sol.bases[k]
        order = size(coefficients, 1) - 1
        products = map(range(lo + 1e-9, hi - 1e-9; length = 97)) do z
            t = (z - (lo + hi)/2)/half
            value = transpose(coefficients)*(t .^ (0:order))
            slope = transpose(coefficients[2:end, :])*((1:order) .* t .^ (0:order-1))/half
            wronskian = value[1]*slope[2] - value[2]*slope[1]
            (_Kat(sol.kcoefs[k], t)*wronskian, wronskian)
        end
        scaled = first.(products)
        # K times the Wronskian is constant, and the Wronskian never vanishes.
        @test maximum(abs.(scaled .- scaled[1]))/abs(scaled[1]) < 1e-12
        @test minimum(abs.(last.(products))) > 1e-6
    end
end

@testset "Mixed models and parameter sharing" begin
    types = [:quadratic, :linear, :constant, :exponential]
    mixed = mixed_model(range(0., 1.; length = 5), types)
    @test mixed.nparameters == 5                  # 4 knots + 1 quadratic midpoint
    @test mixed.left_index[3] == mixed.right_index[3]
    @test mixed.mid_index == [5, 0, 0, 0]
    z = collect(range(0., 1.; length = 41))
    sol = solve_layered(mixed, [.1, .15, .12, .08, .13], 1im, -1im, 0im)
    @test maximum(abs.(field(sol, z) .-
                       reference_field(mixed, sol.values, 1im, -1im, 0im, z))) < 1e-8
    # Two adjacent constant layers share one parameter, so K is flat across both.
    pair = mixed_model(range(0., 1.; length = 5),
                       [:linear, :constant, :constant, :linear])
    @test pair.nparameters == 3
    @test pair.left_index[2] == pair.right_index[2] == pair.left_index[3]
end

@testset "Constant layers keep flux continuous across a jump in K" begin
    edges = collect(range(0., 1.; length = 5))
    model = LayerModel(edges; kind = :constant)
    K = [.05, .25, .09, .3]                       # deliberately discontinuous
    sol = solve_layered(model, K, 1im, -1im, 0im)
    jumps = interface_residuals(sol)
    @test jumps.value < 1e-12 && jumps.flux < 1e-12
    for k in 2:4                                  # K itself does jump
        @test abs(viscosity(sol, edges[k]) - viscosity(sol, edges[k] - 1e-9)) > 1e-3
    end
end

@testset "General quadratic layers" begin
    edges = collect(range(0., 1.; length = 5)); z = collect(range(0., 1.; length = 41))
    model = LayerModel(edges; kind = :quadratic)
    @test model.nparameters == 9                  # 5 knots + 4 midpoints
    knots = [.10, .18, .26, .20, .12]
    # Midpoints on the chord reproduce the affine layer exactly.
    chord = [(knots[k] + knots[k+1])/2 for k in 1:4]
    q = solve_layered(model, [knots; chord], 1im, -1im, 0im)
    lin = solve_layered(LayerModel(edges; kind = :linear), knots, 1im, -1im, 0im)
    @test maximum(abs.(field(q, z) .- field(lin, z))) < 1e-12
    # Curved midpoints: check against the independent reference.
    curved = chord .+ [.03, -.02, .04, -.03]
    sol = solve_layered(model, [knots; curved], 1im, -1im, 0im)
    @test maximum(abs.(field(sol, z) .-
                       reference_field(model, sol.values, 1im, -1im, 0im, z))) < 1e-8
    @test max(interface_residuals(sol).value, interface_residuals(sol).flux) < 1e-12
    @test series_tail(sol) < 1e-14
    # Admissibility is set by the distance to the nearest zero of K in the
    # COMPLEX plane, not by positivity on the real interval. With left = 0.10
    # and right = 0.18 the zeros are a conjugate pair of modulus
    # sqrt(m/(0.14-m)) for midpoints m in (0.0029, 0.137), so the radius
    # exceeds one exactly when m > 0.07 -- even though every such K stays
    # positive on the layer.
    @test layer_radius([.06, .04, .08]) ≈ sqrt(.06/.08)
    @test layer_radius([.10, .04, .04]) ≈ sqrt(.10/.04)
    for m in (.001, .030, .060)                     # radius below one
        @test_throws ArgumentError solve_layered(model, [knots; [m; chord[2:4]]],
                                                 1im, -1im, 0im)
    end
    for m in (.100, .120)                           # radius above one
        @test solve_layered(model, [knots; [m; chord[2:4]]], 1im, -1im, 0im) isa Any
    end
    # Parameter derivatives survive the extra midpoint parameters.
    th = log.([knots; curved])
    loss(t) = sum(abs2, ekman_forward(model, exp.(t), z) .- field(sol, z) .* .9)
    ad = Zygote.gradient(loss, th)[1]
    @test norm(ad - central(loss, th))/norm(ad) < 1e-6
end

@testset "Interpolation orders: quadratic beats linear by one power of h" begin
    truth(z) = .025 + .72z*(1 - z) + .15z^3       # C^3, K_min = 0.025
    zf = collect(range(0., 1.; length = 1001))
    fine = LayerModel(collect(range(0., 1.; length = 129)); kind = :quadratic)
    exact = field(solve_layered(fine, nodal(fine, truth), 1im, -1im, 0im), zf)
    errors = Dict(:linear => Float64[], :quadratic => Float64[])
    for N in (8, 16, 32), kind in (:linear, :quadratic)
        model = LayerModel(collect(range(0., 1.; length = N + 1)); kind)
        sol = solve_layered(model, nodal(model, truth), 1im, -1im, 0im)
        push!(errors[kind], maximum(abs.(field(sol, zf) .- exact)))
    end
    @test all(>(3.2), errors[:linear][1:end-1] ./ errors[:linear][2:end])     # h^2
    @test all(>(6.0), errors[:quadratic][1:end-1] ./ errors[:quadratic][2:end])  # h^3
    @test all(errors[:quadratic] .< errors[:linear])
end

@testset "solve_poly_layers on the operational K-profile" begin
    # K = kappa u* z (1 - z/h)^2 : simple root at the surface, double root at h.
    kappa, ustar, h, f = .4, .3, 1000., 1e-4
    z0, ztop = .1, .99h
    k1 = kappa*ustar/h^2
    gc = [0., k1*h*h, -2k1*h, k1]
    K(z) = gc[1] + gc[2]z + gc[3]z^2 + gc[4]z^3
    build(edges) = [[sum(gc[j+1]*binomial(j, i)*((edges[k] + edges[k+1])/2)^(j-i)*
                         ((edges[k+1] - edges[k])/2)^i for j in i:3) for i in 0:3]
                    for k in 1:length(edges)-1]
    z = collect(range(z0, ztop; length = 121))
    reference = reference_arbitrary(K, 1im*f, -1im, 0im, z)
    # Layers chosen so both roots stay at |x| >= R in the local coordinate.
    function partition(R)
        edges = [z0]
        while edges[end] < ztop - 1e-9
            a = edges[end]
            push!(edges, min(a*(R + 1)/(R - 1), (2h + a*(R - 1))/(R + 1), ztop))
        end
        edges
    end
    respecting = partition(1.5)
    @test length(respecting) - 1 == 8
    sol = solve_poly_layers(respecting, build(respecting), 1im*f, -1im, 0im)
    @test series_tail(sol) < 1e-13
    @test max(interface_residuals(sol).value, interface_residuals(sol).flux) < 1e-12
    # Solver-side invariants: the answer must not depend on the truncation
    # order, nor on which admissible partition is used.
    high = solve_poly_layers(respecting, build(respecting), 1im*f, -1im, 0im; order = 240)
    @test maximum(abs.(field(sol, z) .- field(high, z))) < 1e-13
    coarser = partition(3.0)
    @test length(coarser) - 1 > length(respecting) - 1
    @test maximum(abs.(field(sol, z) .-
                       field(solve_poly_layers(coarser, build(coarser),
                                               1im*f, -1im, 0im), z))) < 1e-12
    # Over this 1496:1 contrast the fixed-step RK4 reference, not the solver,
    # is the accuracy bottleneck: refining it by four converges on the layer
    # solution at the expected fourth order.
    coarse_gap = maximum(abs.(field(sol, z) .- reference))
    fine_gap = maximum(abs.(field(sol, z) .-
                            reference_arbitrary(K, 1im*f, -1im, 0im, z; steps = 80000)))
    @test fine_gap < coarse_gap/50
    @test fine_gap < 1e-6
    # A uniform partition leaves the double root too close. The series tail,
    # which is what a user is told to check, flags it.
    uniform = collect(range(z0, ztop; length = 17))
    bad = solve_poly_layers(uniform, build(uniform), 1im*f, -1im, 0im)
    @test series_tail(bad) > 1e-6
    @test maximum(abs.(field(bad, z) .- reference)) > 1e-3
    # Coefficients that put a root inside a layer are rejected.
    @test_throws ArgumentError solve_poly_layers([0., 1.], [[-.1, 1., 0., 0.]],
                                                 1im, -1im, 0im)
end

@testset "solve_poly_layers reproduces the catalogue rows" begin
    z = collect(range(.1, .9; length = 61))
    cases = ["affine" => [.12, .30, 0., 0.],
             "double-root quadratic" => [.16, .40, .25, 0.],   # (0.5 z + 0.4)^2
             "general quadratic" => [.39, .80, .30, 0.],
             "general cubic" => [.10, .35, .20, .40]]
    for (name, coefficients) in cases
        Kfun(x) = sum(coefficients[i+1]*x^i for i in 0:3)
        @test minimum(Kfun.(z)) > 0
        for N in (1, 4)
            edges = collect(range(.1, .9; length = N + 1))
            local_coefficients = [[sum(coefficients[j+1]*binomial(j, i)*
                                       ((edges[k] + edges[k+1])/2)^(j-i)*
                                       ((edges[k+1] - edges[k])/2)^i for j in i:3)
                                   for i in 0:3] for k in 1:N]
            sol = solve_poly_layers(edges, local_coefficients, 1im, -1im, 0im)
            got = field(sol, z)
            @test maximum(abs.(got .- reference_arbitrary(Kfun, 1im, -1im, 0im, z))) < 1e-8
            # Raising the truncation order must not move the answer.
            high = solve_poly_layers(edges, local_coefficients, 1im, -1im, 0im; order = 240)
            @test maximum(abs.(got .- field(high, z))) < 1e-12
            @test series_tail(sol) < 1e-14
        end
    end
end

@testset "Constant Prandtl profile with dimensional variables" begin
    a = 9.81*sinpi(6/180)/280; b = .003*sinpi(6/180)
    z = collect(range(0., 100.; length = 51))
    model = LayerModel(range(0., 100.; length = 11))
    K = .1; th0 = -5.; delta = sqrt(2K)/(a*b)^(1/4)
    u(z) = th0*sqrt(a/b)*exp(-z/delta)*sin(z/delta)
    theta(z) = th0*exp(-z/delta)*cos(z/delta)
    got = coupled_forward(model, fill(K, 11), z; a, b,
                          lower = (0., th0), upper = (u(100.), theta(100.)))
    @test maximum(abs.(got.u .- u.(z))) < 1e-11
    @test maximum(abs.(got.v .- theta.(z))) < 1e-11
    @test_throws ArgumentError coupled_forward(model, fill(K, 11), z; a = -1., b)
end

@testset "Viscosity evaluation and the right-layer convention" begin
    edges = [0., .5, 1.]
    model = LayerModel(edges; kind = :constant)
    sol = solve_layered(model, [.2, .4], 1im, -1im, 0im)
    @test viscosity(sol, .25) ≈ .2
    @test viscosity(sol, .5) ≈ .4            # a jump takes the right layer
    @test viscosity(sol, 1.) ≈ .4            # except at the top endpoint
    @test viscosity(model, [.2, .4], .25) ≈ .2
    linear = LayerModel(edges; kind = :linear)
    @test viscosity(linear, [.1, .2, .3], .5) ≈ .2
end

@testset "Determinism and truncation order" begin
    edges = collect(range(0., 1.; length = 6))
    model = LayerModel(edges; kind = :quadratic)
    K = nodal(model, z -> .1 + .2z*(1 - z))
    z = collect(range(0., 1.; length = 33))
    first_run = field(solve_layered(model, K, 1im, -1im, 0im), z)
    @test first_run == field(solve_layered(model, K, 1im, -1im, 0im), z)
    @test maximum(abs.(first_run .-
                       field(solve_layered(model, K, 1im, -1im, 0im; order = 200), z))) < 1e-13
end

@testset "Input validation" begin
    @test_throws ArgumentError LayerModel([0., 0., 1.])
    @test_throws ArgumentError LayerModel([1., 0.])
    @test_throws ArgumentError LayerModel([0.])
    @test_throws ArgumentError LayerModel([0., 1.]; kind = :cubic)
    @test_throws ArgumentError mixed_model([0., .5, 1.], [:linear, :sextic])
    @test_throws DimensionMismatch mixed_model([0., .5, 1.], [:linear])
    m = LayerModel([0., 1.])
    @test_throws ArgumentError solve_layered(m, [-.1, .1], 1im, -1im, 0im)
    @test_throws ArgumentError solve_layered(m, [0., .1], 1im, -1im, 0im)
    @test_throws ArgumentError solve_layered(m, [NaN, .1], 1im, -1im, 0im)
    @test_throws DimensionMismatch solve_layered(m, [.1], 1im, -1im, 0im)
    @test_throws ArgumentError solve_layered(m, [.1, .1], 1im, -1im, 0im; order = 3)
    sol = solve_layered(m, [.1, .1], 1im, -1im, 0im)
    @test_throws ArgumentError field(sol, 2.)
    @test_throws ArgumentError field(sol, -1.)
    # An exponential layer spanning too large a contrast is refused rather than
    # silently truncated below the series tail.
    wide = LayerModel([0., 1.]; kind = :exponential)
    @test_throws ArgumentError solve_layered(wide, [1e-9, 1e9], 1im, -1im, 0im)
    @test_throws DimensionMismatch solve_poly_layers([0., 1.], [[.1], [.1]], 1im, -1im, 0im)
    @test_throws ArgumentError solve_poly_layers([0., 1.], [Float64[]], 1im, -1im, 0im)
    @test_throws DimensionMismatch fit_viscosity(m, [0., 1.], [1im]; initial = [.1, .1])
    @test_throws ArgumentError fit_viscosity(m, [0., 1.], [1im, 2im];
                                             initial = [.1, .1], lower = .2, upper = .5)
end

@testset "Sparse high-degree Taylor diagnostics" begin
    # The old last-two-coefficient diagnostic was 2.66e-18, although this
    # admissible positive profile has a 6.47e-4 one-layer truncation error.
    c=zeros(17);c[1]=1.;c[end]=.9
    z=collect(range(0.,1.;length=81));K(z)=1+.9*(2z-1)^16
    coarse=solve_poly_layers([0.,1.],[c],1im,-1im,0im)
    reference=reference_arbitrary(K,1im,-1im,0im,z)
    @test maximum(abs.(field(coarse,z).-reference))>1e-4
    @test series_tail(coarse)>1e-4
    edges=collect(range(0.,1.;length=9))
    coeffs=[[(i==0 ? 1. : 0.)+.9*binomial(16,i)*(edges[k]+edges[k+1]-1)^(16-i)*
              (edges[k+1]-edges[k])^i for i in 0:16] for k in 1:8]
    fine=solve_poly_layers(edges,coeffs,1im,-1im,0im)
    @test maximum(abs.(field(fine,z).-reference))<1e-11
    @test series_tail(fine)<1e-10
    # Coefficients beyond the retained order must not disappear from diagnosis.
    long=zeros(170);long[1]=1;long[end]=.1
    limited=solve_poly_layers([0.,1.],[long],1im,-1im,0im)
    @test series_tail(limited)>0.1
end

@testset "Inverse line search rejects inadmissible trials" begin
    model=LayerModel([0.,1.];kind=:quadratic)
    z=collect(range(0.,1.;length=41))
    K(z)=.5*(2z-1)*(z-1)+4*.07*z*(1-z)+.5*z*(2z-1)
    observed=reference_arbitrary(K,1im,-1im,0im,z)
    fit=fit_viscosity(model,z,observed;initial=[.1,.1,.1],iterations=100,outer_iterations=5)
    @test isfinite(fit.loss)
    @test all(isfinite,fit.prediction)
    @test fit.diagnostics.admissibility_rejections>0
    @test fit.diagnostics.refinement_required
    @test !fit.converged
    @test_throws ArgumentError fit_viscosity(model,z,observed;initial=[.49,.49,.07])
    # A representable target must still be recovered rather than rejected.
    truth=[.1,.1,.09];target=ekman_forward(model,truth,z)
    recovered=fit_viscosity(model,z,target;initial=[.1,.1,.1],iterations=150)
    @test recovered.loss<1e-12
    @test maximum(abs.(recovered.values.-truth))<1e-4
end

@testset "Paper uncertainty intervals and noise convention" begin
    include(joinpath(@__DIR__,"..","examples","paper_reproduction.jl"))
    model=LayerModel([0.,.5,1.];kind=:constant);z=collect(range(0.,1.;length=41))
    values=[.12,.18];noise=fill(.001+.002im,length(z))
    u=PaperReproduction.local_uncertainty(model,values,z,noise;draws=2000)
    K=viscosity(model,values,z)
    @test u.rank==2
    @test u.sigma≈norm(noise)/sqrt(2length(z))
    @test all(u.lower.<=K) && all(K.<=u.upper)
    @test u.fraction_admissible>0.99
    @test PaperReproduction.profile_matrix(model,[.5])*values≈[.18]
end
