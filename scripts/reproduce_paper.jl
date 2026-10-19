# From the package root: julia --project=. scripts/reproduce_paper.jl [--refit] [--output=DIR]
using Pkg
root=dirname(@__DIR__);Pkg.activate(root)
ENV["GKSwstype"]="100"
include(joinpath(root,"examples","paper_reproduction.jl"))
using .PaperReproduction, Plots
function main(args=ARGS)
output=joinpath(root,"results","paper")
for arg in args
    startswith(arg,"--output=") && (output=abspath(arg[10:end]))
    arg=="--refit" || startswith(arg,"--output=") || error("unknown argument: $arg")
end
interpolation_table(output)
catalogue_checks(output)
data,results=saved_inversion(root)
if "--refit" in args
    results=fit_models(data)
    println("Inversion optimized again from the documented symmetric initial profile.")
else
    println("Reevaluating archived fitted parameters; use --refit to optimize again.")
end
figure,intervals=inversion_figure(data,results)
savefig(figure,joinpath(output,"figure2_inversion.png"));savefig(figure,joinpath(output,"figure2_inversion.pdf"))
write_inversion(output,results,intervals,data[:,1])
println("Written: ",output)
println("Independent closed forms: python3 reproduction/check_closed_forms.py --results=",output)

end
main()
