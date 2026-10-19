module PaperReproduction
using ELES, LinearAlgebra, Random, DelimitedFiles, Printf, Plots
include("reference_shooting.jl")
export interpolation_table, catalogue_checks, local_uncertainty, inversion_figure,
       saved_inversion, fit_models, write_inversion

function nodal_values(model,f)
    p=zeros(model.nparameters)
    for k in eachindex(model.kinds)
        lo,hi=model.interfaces[k:k+1]
        if model.kinds[k]==:constant
            p[model.left_index[k]]=f((lo+hi)/2)
        else
            p[model.left_index[k]]=f(lo);p[model.right_index[k]]=f(hi)
            model.mid_index[k]==0 || (p[model.mid_index[k]]=f((lo+hi)/2))
        end
    end
    p
end
shift(c,m,h)=[sum(c[j+1]*binomial(j,i)*m^(j-i)*h^i for j in i:length(c)-1) for i in 0:length(c)-1]
function expand_factors(scale,roots)
    c=[scale]
    for r in roots
        b=zeros(length(c)+1)
        for i in eachindex(c);b[i]+=r*c[i];b[i+1]+=c[i];end
        c=b
    end
    c
end
# Exact extrema of a cubic interpolation error in the local coordinate.
function cubic_sup(c)
    points=[-1.,1.];a,b,d=3c[4],2c[3],c[2]
    if a==0
        b==0 || push!(points,-d/b)
    else
        disc=b*b-4a*d
        disc<0 || append!(points,[(-b+sqrt(disc))/(2a),(-b-sqrt(disc))/(2a)])
    end
    maximum(abs(sum(c[i+1]*x^i for i in 0:3)) for x in points if -1<=x<=1)
end

"""Table 2; velocity errors are sampled on 201 uniformly spaced heights."""
function interpolation_table(output)
    mkpath(output);K(z)=.025+.72z*(1-z)+.15z^3
    z=collect(range(0.,1.;length=201));reference=reference_arbitrary(K,1im,-1im,0im,z)
    records=NamedTuple[]
    open(joinpath(output,"table2_convergence.csv"),"w") do io
        println(io,"degree,layers,K_sup,ratio,bound,velocity_max_201_points")
        for (degree,kind) in enumerate((:constant,:linear,:quadratic))
            d=degree-1;previous=NaN
            for n in (8,16,32,64)
                model=LayerModel(range(0.,1.;length=n+1);kind)
                sol=solve_layered(model,nodal_values(model,K),1im,-1im,0im)
                err=maximum(eachindex(model.kinds)) do k
                    lo,hi=model.interfaces[k:k+1]
                    c=shift([.025,.72,-.72,.15],(lo+hi)/2,(hi-lo)/2)
                    c[1:length(sol.kcoefs[k])].-=sol.kcoefs[k]
                    cubic_sup(c)
                end
                bound=(.5,.125,1/(72sqrt(3)))[degree]*(.72,1.44,.9)[degree]/n^(d+1)
                velocity=maximum(abs.(field(sol,z).-reference));ratio=previous/err
                push!(records,(degree=d,layers=n,K_sup=err,ratio=ratio,bound=bound,velocity=velocity))
                @printf(io,"%d,%d,%.17g,%.17g,%.17g,%.17g\n",d,n,err,ratio,bound,velocity)
                previous=err
            end
        end
    end
    records
end

"""Seven profiles on one/four layers, plus the quadratic root configurations.
Exports current Julia velocities for an independent mpmath closed-form check.
"""
function catalogue_checks(output)
    mkpath(output);z=collect(range(0.,1.;length=81))
    polys=[("affine",[.1,.3]),
           ("double_root_quadratic",expand_factors(.25,[.8,.8])),
           ("general_quadratic",expand_factors(.6,[.5,1.3])),
           ("cubic_power",expand_factors(.4,[.6,.6,.6])),
           ("restricted_cubic",expand_factors(.5,[.6,.6,1.4])),
           ("quartic_power",expand_factors(.3,[.7,.7,.7,.7]))]
    cases=NamedTuple[(name=name,K=z->sum(c[i+1]*z^i for i in 0:length(c)-1),
            build=(m,h)->shift(c,m,h),layers=(1,4)) for (name,c) in polys]
    push!(cases,(name="exponential",K=z->.12exp(1.8z),
                 build=(m,h)->[.12exp(1.8m)*(1.8h)^i/prod(1.0:i) for i in 0:24],layers=(1,4)))
    for (name,c) in [("roots_straddle",expand_factors(-.6,[.2,-1.2])),
                     ("roots_above",expand_factors(.5,[-1.4,-2.2])),
                     ("roots_complex",[.2125,-.05,.05])]
        push!(cases,(name=name,K=z->sum(c[i+1]*z^i for i in 0:length(c)-1),
                     build=(m,h)->shift(c,m,h),layers=(1,)))
    end
    records=NamedTuple[]
    open(joinpath(output,"catalogue_metrics.csv"),"w") do metrics
        println(metrics,"case,layers,vs_RK4,interface_jump,series_tail")
        open(joinpath(output,"catalogue_velocities.csv"),"w") do values
            println(values,"case,layers,z,real,imag")
            for case in cases,n in case.layers
                edges=collect(range(0.,1.;length=n+1))
                coeffs=[case.build((edges[k]+edges[k+1])/2,(edges[k+1]-edges[k])/2) for k in 1:n]
                sol=solve_poly_layers(edges,coeffs,1im,-1im,0im)
                w=field(sol,z);ref=reference_arbitrary(case.K,1im,-1im,0im,z)
                err=maximum(abs.(w.-ref));j=interface_residuals(sol);jump=max(j.value,j.flux)
                @assert err<1e-10 && jump<1e-12
                tail=series_tail(sol)
                push!(records,(case=case.name,layers=n,error=err,jump=jump,tail=tail))
                @printf(metrics,"%s,%d,%.17g,%.17g,%.17g\n",case.name,n,err,jump,tail)
                for i in eachindex(z)
                    @printf(values,"%s,%d,%.17g,%.17g,%.17g\n",case.name,n,z[i],real(w[i]),imag(w[i]))
                end
            end
        end
    end
    records
end

const LAYERS=Dict(:constant=>11,:linear=>10,:quadratic=>5)
function saved_inversion(root)
    data=readdlm(joinpath(root,"data","synthetic_velocity.csv"),',',Float64;skipstart=1)
    params=readdlm(joinpath(root,"data","fitted_parameters.csv"),',',Any;skipstart=1)
    z=data[:,1];truth=data[:,2];observed=data[:,5].+im.*data[:,6];results=Dict()
    for kind in (:constant,:linear,:quadratic)
        rows=findall(i->params[i,1]==string(kind),axes(params,1))
        values=Float64.(params[rows,4]);model=LayerModel(range(0.,1.;length=LAYERS[kind]+1);kind)
        sol=solve_layered(model,values,1im,-1im,0im);K=viscosity(sol,z);w=field(sol,z)
        fit=(values=values,prediction=w,loss=sum(abs2,w.-observed)/length(z))
        results[kind]=(K,fit,norm(K.-truth)/norm(truth))
    end
    data,results
end
function fit_models(data)
    z=data[:,1];observed=data[:,5].+im.*data[:,6];results=Dict()
    for kind in (:constant,:linear,:quadratic)
        model=LayerModel(range(0.,1.;length=LAYERS[kind]+1);kind)
        initial=nodal_values(model,z->.085+.94z*(1-z))
        fit=fit_viscosity(model,z,observed;initial)
        @assert !fit.diagnostics.refinement_required
        K=viscosity(solve_layered(model,fit.values,1im,-1im,0im),z)
        results[kind]=(K,fit,norm(K.-data[:,2])/norm(data[:,2]))
    end
    results
end
function profile_matrix(model,z)
    A=zeros(length(z),model.nparameters)
    for (i,v) in enumerate(z)
        k=clamp(searchsortedlast(model.interfaces,v),1,length(model.kinds))
        lo,hi=model.interfaces[k:k+1];x=(2v-lo-hi)/(hi-lo)
        left,right=model.left_index[k],model.right_index[k]
        if model.kinds[k]==:constant;A[i,left]=1
        elseif model.kinds[k]==:linear;A[i,left]=(1-x)/2;A[i,right]=(1+x)/2
        elseif model.kinds[k]==:quadratic
            A[i,left]=x*(x-1)/2;A[i,right]=x*(x+1)/2;A[i,model.mid_index[k]]=1-x*x
        else;throw(ArgumentError("uncertainty helper supports degrees 0,1,2"))
        end
    end
    A
end
function percentile(v,p)
    s=sort(v);x=1+(length(s)-1)*p;i=floor(Int,x);j=min(i+1,length(s))
    s[i]+(x-i)*(s[j]-s[i])
end
"""Box-conditioned local Gauss-Newton intervals; not a nonlinear posterior.
`noise` supplies the known synthetic residual, never the inverse objective.
"""
function local_uncertainty(model,values,z,noise;draws=20000,rng=MersenneTwister(20261017))
    draws>=2 || throw(ArgumentError("at least two draws are required"))
    theta=log.(values);p=length(theta);Iden=Matrix{Float64}(I,p,p)
    stacked(th)=begin w=ekman_forward(model,exp.(th),z);[real.(w);imag.(w)] end
    jac(h)=hcat([(stacked(theta+h*Iden[:,j])-stacked(theta-h*Iden[:,j]))/(2h) for j in 1:p]...)
    J=jac(1e-5);stepcheck=norm(J-jac(5e-6))/norm(J)
    @assert rank(J)==p && stepcheck<1e-6
    sigma=norm(noise)/sqrt(2length(z));cov=sigma^2*inv(Symmetric(J'J))
    logs=theta.+cholesky(Symmetric(cov)).L*randn(rng,p,draws)
    in_box=vec(all((log(.06).<=logs).&(logs.<=log(.50));dims=1))
    accepted=Int[]
    for j in findall(in_box)
        sample=exp.(logs[:,j])
        try
            for k in eachindex(model.kinds);ELES._kcoef(model,sample,k);end
            push!(accepted,j)
        catch err
            err isa ArgumentError || rethrow()
        end
    end
    length(accepted)>=2 || throw(ArgumentError("too few admissible local Gaussian samples"))
    curves=profile_matrix(model,z)*exp.(logs[:,accepted])
    lower=[percentile(collect(row),.025) for row in eachrow(curves)]
    upper=[percentile(collect(row),.975) for row in eachrow(curves)]
    selected=argmin(abs.(z.-.5))
    (lower=lower,upper=upper,selected=curves[selected,:],covariance=cov,
     rank=rank(J),step_check=stepcheck,sigma=sigma,
     fraction_in_box=count(in_box)/draws,fraction_admissible=length(accepted)/draws)
end
function inversion_figure(data,results;draws=20000)
    z=data[:,1];truth=data[:,2];clean=data[:,3].+im.*data[:,4];obs=data[:,5].+im.*data[:,6]
    rng=MersenneTwister(20261017);uncertainty=Dict();kinds=(:constant,:linear,:quadratic)
    colors=(:purple,:red,:teal);styles=(:dashdot,:dash,:solid)
    a=plot(truth,z;label="true K",lc=:black,xlabel="K(z)",ylabel="z",title="(a) recovery; 95% local intervals")
    b=plot(;xlabel="K(0.5+)",ylabel="density",title="(b) local uncertainty")
    c=scatter(real.(obs),z;label="u data",ms=2,msw=0,mc=:gray,xlabel="velocity defect",ylabel="z",title="(c) observations and fits")
    scatter!(c,imag.(obs),z;label="v data",ms=2,msw=0,mc=:lightgray)
    for (i,kind) in enumerate(kinds)
        K,fit,EK=results[kind];model=LayerModel(range(0.,1.;length=LAYERS[kind]+1);kind)
        u=local_uncertainty(model,fit.values,z,obs.-clean;draws,rng);uncertainty[kind]=u
        label=@sprintf("degree %d (%0.2f%%)",i-1,100EK)
        band=Shape(vcat(u.lower,reverse(u.upper)),vcat(z,reverse(z)))
        plot!(a,band;color=colors[i],fillalpha=.12,linealpha=0,label="")
        plot!(a,K,z;lc=colors[i],ls=styles[i],label)
        histogram!(b,u.selected;bins=60,normalize=:pdf,seriestype=:stephist,lc=colors[i],label="degree $(i-1)")
        plot!(c,real.(fit.prediction),z;lc=colors[i],ls=styles[i],label="")
        plot!(c,imag.(fit.prediction),z;lc=colors[i],ls=styles[i],alpha=.6,label="")
    end
    vline!(b,[truth[argmin(abs.(z.-.5))]];lc=:black,label="true K")
    plot(a,b,c;layout=(1,3),size=(1150,360),legendfontsize=7,titlefontsize=9),uncertainty
end
function write_inversion(output,results,uncertainty,z)
    mkpath(output)
    open(joinpath(output,"table3_inversion.csv"),"w") do io
        println(io,"degree,layers,parameters,E_K,loss,jacobian_rank,FD_step_check,fraction_in_box,fraction_admissible")
        for (d,kind) in enumerate((:constant,:linear,:quadratic))
            K,fit,EK=results[kind];u=uncertainty[kind]
            @printf(io,"%d,%d,%d,%.17g,%.17g,%d,%.17g,%.17g,%.17g\n",d-1,LAYERS[kind],length(fit.values),EK,fit.loss,u.rank,u.step_check,u.fraction_in_box,u.fraction_admissible)
        end
    end
    open(joinpath(output,"local_intervals.csv"),"w") do io
        println(io,"model,z,K_fit,lower95,upper95")
        for kind in (:constant,:linear,:quadratic),i in eachindex(z)
            u=uncertainty[kind]
            @printf(io,"%s,%.17g,%.17g,%.17g,%.17g\n",string(kind),z[i],results[kind][1][i],u.lower[i],u.upper[i])
        end
    end
end
end
