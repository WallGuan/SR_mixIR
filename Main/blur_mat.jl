#=
blur_mat.jl

Runs the PRblur deblurring problem (noise-free, and large noise with Tikhonov
regularization) with RTN and SR mixed-precision IR, then saves the restored
images as .mat files so they can be displayed in MATLAB.

Run from the Main/ folder:  julia blur_mat.jl
=#

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

Random.seed!(2026) # fix the noise draw so runs are repeatable (SR itself stays random)
btrue_large = bnoise(btrue, level=1e-1) # large noise: ||noise|| / ||b|| = 0.1

###

### No nosie 

# IR >> fixed result
start_time = time()
x_rtn, error = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = false, imax = 10, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for RTN no noise: ", elapsed)
println("Minimum relative error is $(minimum(error)) at iteration $(argmin(error))")

###

# SR 
start_time = time()
x_sr, error = IR_GMRES(A, btrue, x0; x_true=xtrue, SR = true, imax = 10, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for SR no noise: ", elapsed)
println("Minimum relative error is $(minimum(error)) at iteration $(argmin(error))")

### Large nosie 

# IR >> fixed result
start_time = time()
x_reg_rtn, error = IR_kronreg(A, btrue_large, x0; x_true=xtrue, SR = false, lambda = 0.1, imax = 20, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for RTN: ", elapsed)
println("Minimum relative error is $(minimum(error)) at iteration $(argmin(error))")

###

# SR 
start_time = time()
x_reg_sr, error = IR_kronreg(A, btrue_large, x0; x_true=xtrue, SR = true, lambda = 0.1, imax = 5, GMRES_loop = 15)
elapsed = time() - start_time
println("Run time for SR: ", elapsed)
println("Minimum relative error is $(minimum(error)) at iteration $(argmin(error))")

# save in .mat
# Note: the key "x_true" holds the RESTORED solution (name kept for the MATLAB
# display script); "b" is the blurred image that was actually solved for.
matwrite("no_RTN.mat", Dict(
    "A1" => A1,
    "A2" => A2,
    "b" => btrue,
    "x_true" => x_rtn
))

matwrite("no_SR.mat", Dict(
    "A1" => A1,
    "A2" => A2,
    "b" => btrue,
    "x_true" => x_sr
))

matwrite("large_RTN.mat", Dict(
    "A1" => A1,
    "A2" => A2,
    "b" => btrue_large, # was btrue: save the noisy data this case actually used
    "x_true" => x_reg_rtn
))

matwrite("large_SR.mat", Dict(
    "A1" => A1,
    "A2" => A2,
    "b" => btrue_large, # was btrue: save the noisy data this case actually used
    "x_true" => x_reg_sr
))