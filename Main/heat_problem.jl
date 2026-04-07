# call on other files for pre-defined functions
include("refinement_funcs.jl")
include("other_funcs.jl")
include("kronecker/main.jl")
# libraries 
using .KronMatrixModule
using LinearAlgebra
using MAT
using Plots
using BenchmarkTools
using StochasticRounding
using RegularizationTools
using Random
using ToeplitzMatrices

###

# define problem 
n = 256
A_ill, b_ill, x_ill = heat(n, 4.0) # original 
A, b, x = AUG(A_ill, b_ill, x_ill, 0.06) # block reg 

x0 = ones(2*n) # initial guess for regularized case
# print basic info 
println("Original condition num: ", cond(A_ill))
println("Regularize condition num: ", cond(A))

###

println()

## LU Factorization 

println("LU Factorization")

# RTN reg 
start_time = time()
x_rtn, error_rtn = IR_fact(A, b, x0; SR = false, fact = lu, x_true=x, imax=10)
elapsed = time() - start_time
println("RTN LU Run time: ", elapsed) # runtime check
println("Minimum relative error is $(minimum(error_rtn)) at iteration $(argmin(error_rtn))") # error check

# SR reg 
n_runs = 3
x_sr = Vector{AbstractArray}(undef, n_runs)
errors_sr = Vector{Vector{Float64}}(undef, n_runs)

for i in 1:n_runs
    start_time = time()
    x_sr[i], errors_sr[i] = IR_fact(A, b, x0; SR = true, fact = lu, x_true = x, imax = 10)
    elapsed = time() - start_time

    println("Run $i SR LU Run time: ", elapsed)
    println("Minimum relative error is $(minimum(errors_sr[i])) at iteration $(argmin(errors_sr[i]))")
end

###

println()

## GMRES

println("GMRES")

# RTN reg 
start_time = time()
x_rtn_g, error_rtn_g = IR_GMRES(A, b, x0; SR = false, x_true=x, imax=10, atol=1e-5, GMRES_loop = 15)
elapsed = time() - start_time
println("RTN GMRES Run time: ", elapsed) # runtime check
println("Minimum relative error is $(minimum(error_rtn_g)) at iteration $(argmin(error_rtn_g))") # error check

# SR reg 
n_runs = 3
x_sr_g = Vector{AbstractArray}(undef, n_runs)
errors_sr_g = Vector{Vector{Float64}}(undef, n_runs)

for i in 1:n_runs
    start_time = time()
    x_sr_g[i], errors_sr_g[i] = IR_GMRES(A, b, x0; SR = true, x_true=x, imax=10, atol=1e-5, GMRES_loop = 15)
    elapsed = time() - start_time

    println("Run $i SR GMRES Run time: ", elapsed)
    println("Minimum relative error is $(minimum(errors_sr_g[i])) at iteration $(argmin(errors_sr_g[i]))")
end

###

## plot

# LU IR
p1 = plot(
    error_rtn;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "LU_IR inverse heat matrix",
    lw = 2,
    marker = :circle,
    color = :orange,
    label = "RTN"
)

colors = [:navy, :royalblue, :skyblue]

for i in 1:length(errors_sr)
    plot!(p1, errors_sr[i];
        lw = 1,
        marker = :circle,
        color = colors[i],
        label = "SR$i"
    )
end

# GMRES IR
p2 = plot(
    error_rtn_g;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "GMRES_IR inverse heat matrix",
    lw = 2,
    marker = :circle,
    color = :orange,
    label = "RTN"
)

for i in 1:length(errors_sr_g)
    plot!(p2, errors_sr_g[i];
        lw = 1,
        marker = :circle,
        color = colors[i],
        label = "SR$i"
    )
end

p = plot(p1, p2; layout = (1,2), size = (1000,500), margin = 10Plots.mm)
savefig(p, "inverse_heat.png")
