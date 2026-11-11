# One-time setup: creates a Julia environment in this folder and installs
# the packages the experiments use. Run from the repository root:
#     julia setup.jl
using Pkg
Pkg.activate(@__DIR__)
Pkg.add([
    "MAT",                  # read/write .mat files (PRblur data, restored images)
    "Plots",                # convergence plots
    "BenchmarkTools",       # timing
    "StochasticRounding",   # Float16sr and stochastic_round
    "RegularizationTools",  # regularization utilities
    "ToeplitzMatrices",     # Toeplitz matrix for the inverse heat problem
])
println("Done. Run experiments from Main/ with:  julia --project=.. kron_problem.jl")
