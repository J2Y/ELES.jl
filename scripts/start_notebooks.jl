# One-command setup and launch: julia --project=. scripts/start_notebooks.jl
include("setup_notebook.jl")
IJulia.notebook(dir=joinpath(root,"notebooks"))
