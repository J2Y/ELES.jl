# Run once: julia --project=. scripts/setup_notebook.jl
using Pkg
root = dirname(@__DIR__)
Pkg.activate(root)
Pkg.instantiate()
using IJulia
IJulia.installkernel("ELES", "--project=$(root)";
                    specname="eles",displayname="Julia - ELES",
                    env=Dict("GKSwstype"=>"100"))
println("Select the Julia kernel named ELES, then Run All in notebooks/.")
println("To open Jupyter from Julia: using IJulia; IJulia.notebook(dir=joinpath(\"$(root)\",\"notebooks\"))")
