include("kronecker/main.jl")
include("refinement_funcs.jl")
include("other_funcs.jl")
using .KronMatrixModule
using LinearAlgebra
using MAT
using Plots
using BenchmarkTools
using StochasticRounding
using RegularizationTools
using Random

# input 
data = matread("PRblurData_new.mat")
A1 = data["A1"]
A2 = data["A2"]
btrue = data["btrue"]
xtrue = data["xtrue"]
A = KronMatrix2(A1, A2)
x0 = ones(size(xtrue))

### 

# define noise level (or no noise) and reg values

# btrue = bnoise(btrue, level=1e-3) # small noise
# btrue = bnoise(btrue, level=1e-2) # medium noise 
# btrue = bnoise(btrue, level=1e-1) # large noise 
lam = 0.06 # type of reg

###

### Single precision RTN
println("Single-precision RTN")
println()
# Here contain algrothm for only RTN, single precision will be use for non-reg and reg

println("NO reg output for single precision")
println()

println("Float64")
start_time = time()
x_64, error_64 = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = false,
ul = Float64, u = Float64, ur = Float64, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for RTN_fl64: ", elapsed)
println("Minimum relative error is $(minimum(error_64)) at iteration $(argmin(error_64))")
println()

println("Float32")
start_time = time()
x_32, error_32 = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = false,
ul = Float32, u = Float32, ur = Float32, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for RTN_fl32: ", elapsed)
println("Minimum relative error is $(minimum(error_32)) at iteration $(argmin(error_32))")
println()

println("Float16")
start_time = time()
x_16, error_16 = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = false, 
ul = Float16, u = Float16, ur = Float16, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for RTN_fl16: ", elapsed)
println("Minimum relative error is $(minimum(error_16)) at iteration $(argmin(error_16))")
println()

###

println("REG output for single precision")
println()

println("Float64")
start_time = time()
x_64reg, error_64reg = IR_kronreg(A, btrue, x0; x_true=xtrue, SR = false, lambda = lam, 
ul = Float64, u = Float64, ur = Float64, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for reg_RTN_fl64: ", elapsed)
println("Minimum relative error is $(minimum(error_64reg)) at iteration $(argmin(error_64reg))")
println()

println("Float32")
start_time = time()
x_32reg, error_32reg = IR_kronreg(A, btrue, x0; x_true=xtrue, SR = false, lambda = lam, 
ul = Float32, u = Float32, ur = Float32, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for reg_RTN_fl32: ", elapsed)
println("Minimum relative error is $(minimum(error_32reg)) at iteration $(argmin(error_32reg))")
println()

println("Float16")
start_time = time()
x_16reg, error_16reg = IR_kronreg(A, btrue, x0; x_true=xtrue, SR = false, lambda = lam, 
ul = Float16, u = Float16, ur = Float16, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for reg_RTN_fl16: ", elapsed)
println("Minimum relative error is $(minimum(error_16reg)) at iteration $(argmin(error_16reg))")
println()

###

### Mix-precision no reg
println("Mix-precision NO reg")
println()

# IR >> fixed result
println("RTN")
start_time = time()
x, error = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = false, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for mix_RTN: ", elapsed)
println("Minimum relative error is $(minimum(error)) at iteration $(argmin(error))")
println()

###

# SR 
println("SR")
n_runs = 3
x_sr = Vector{AbstractArray}(undef, n_runs)
error_sr = Vector{Vector{Float64}}(undef, n_runs)

for i in 1:n_runs
    println(i)
    start_time = time()
    x_sr[i], error_sr[i] = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = true, imax = 20, GMRES_loop = 15)
    elapsed = time() - start_time
    println("Run $i mix_SR_$(i) runtime: ", elapsed)

    println("Minimum relative error is $(minimum(error_sr[i])) at iteration $(argmin(error_sr[i]))")

    println()
end

###

### Mix-precision reg
println("Mix-precision with reg")
println()

# Here are RTN and SR mix-precision with regularization codes. 

# IR >> fixed result
println("RTNreg")
start_time = time()
x_reg, error_reg = IR_kronreg(A, btrue, x0; x_true=xtrue, SR = false, lambda = lam, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for mix_RTN_reg: ", elapsed)
println("Minimum relative error is $(minimum(error_reg)) at iteration $(argmin(error_reg))")
println()

###

# SR 
println("SRreg")
n_runs = 3
x_sr_reg = Vector{AbstractArray}(undef, n_runs)
error_sr_reg = Vector{Vector{Float64}}(undef, n_runs)

for i in 1:3
    println(i)
    start_time = time()
    x_sr_reg[i], error_sr_reg[i] = IR_kronreg(A, btrue, x0; x_true=xtrue, SR = true, lambda = lam, imax = 20, GMRES_loop = 15)
    elapsed = time() - start_time
    println("Run $i mix_SR_reg_$(i) runtime: ", elapsed)

    println("Minimum relative error is $(minimum(error_sr_reg[i])) at iteration $(argmin(error_sr_reg[i]))")

    println()
end

###

# plot 

# single precision 

# No reg
p1 = plot(
    error_64;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "NO Regularization",
    lw = 3,
    marker = :circle,
    label = "Float64"
)

plot!(p1, error_32;
    lw = 2,
    marker = :circle,
    label = "Float32"
)

plot!(p1, error_16;
    lw = 1,
    marker = :circle,
    label = "Float16"
)

# Reg
p2 = plot(
    error_64reg;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "With Regularization",
    lw = 3,
    marker = :circle,
    label = "Float64"
)

plot!(p2, error_32reg;
    lw = 2,
    marker = :circle,
    label = "Float32"
)

plot!(p2, error_16reg;
    lw = 1,
    marker = :circle,
    label = "Float16"
)


ps = plot(
    p1, p2;
    layout = (1,2),
    size = (1200,450),
    margin = 10Plots.mm,
    plot_title = "Noise free single precision PRblur",
    plot_titlefontsize = 16
)
savefig(ps, "single_error.png")

# mix precison 

# No reg
p3 = plot(
    error;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "No Regularization",
    lw = 2,
    marker = :circle,
    color = :orange,
    label = "RTN"
)

# colors for SR runs
colors = [:navy, :royalblue, :skyblue]

for i in 1:length(error_sr)
    plot!(p3, error_sr[i];
        lw = 1,
        marker = :circle,
        color = colors[i],
        label = "SR$i"
    )
end


# Reg
p4 = plot(
    error_reg;
    yscale = :log10,
    xlabel = "Iteration",
    ylabel = "Relative Error",
    title = "With Regularization",
    lw = 2,
    marker = :circle,
    color = :orange,
    label = "RTN"
)

for i in 1:length(error_sr_reg)
    plot!(p4, error_sr_reg[i];
        lw = 1,
        marker = :circle,
        color = colors[i],
        label = "SR$i"
    )
end


pm = plot(
    p3, p4;
    layout = (1,2),
    size = (1200,450),
    margin = 10Plots.mm,
    plot_title = "Noise free mix precision PRblur",
    plot_titlefontsize = 16
)
savefig(pm, "mix_error.png")
