#=
refinement_funcs.jl

Mixed-precision iterative refinement (IR) solvers used in the thesis experiments.

Every solver follows the same three-precision pattern (Carson & Higham style IR):

    repeat for i = 1, ..., imax
        r = b - A*x          residual          -> computed in high precision  `ur`
        solve A*d = r        correction step   -> computed in low precision   `ul`
        x = x + d            update            -> computed in working precision `u`

Stochastic rounding (SR) is switched on with `SR = true`. In that case the
low-precision quantities (the matrix used in the correction solve and each
residual) are rounded to `Float16sr` from StochasticRounding.jl: a Float16
number whose arithmetic rounds up or down at random, with probability
proportional to the distance to each neighbour. When `SR = false`, standard
round-to-nearest (RTN) conversion to `ul` is used instead.

Note: when `SR = true` the low precision is always Float16 (via `Float16sr`);
the `ul` keyword only applies to the RTN runs.

Returned `error` vectors start with the relative error of the initial guess
(iteration 0), so `error[k]` is the error after k-1 refinement steps.
=#

"""
    IR_fact(A, b, x0; SR=true, fact=lu, x_true=nothing,
            ul=Float16, u=Float32, ur=Float64, imax=10, atol=1e-12)

Iterative refinement with a direct (factorization-based) correction solve.
`A` is factorized once in low precision with `fact` (default `lu`); the
factorization is reused for every correction step.

Returns `(x, error)` where `error` holds the relative errors
`norm(x - x_true) / norm(x_true)` (only filled if `x_true` is given).
"""
function IR_fact(A, b, x0; SR = true, fact = lu, x_true=nothing,
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12)

    x = x0; error = ur[]  # initialize x and relative errors
    # low-precision copy of A: stochastic rounding or round-to-nearest
    if SR; A_low = stochastic_round.(Float16sr, A); else; A_low = ul.(A); end
    FACT = fact(A_low) # factorize once in low precision, reuse every iteration

    # relative error of the initial guess (iteration 0) to standardize the graphs
    if x_true !== nothing
        push!(error, norm(x - x_true) / norm(x_true))
    end

    for i in 1:imax
        r = ur.(b) - ur.(A) * ur.(x) # residual in high precision
        if SR; r = stochastic_round.(Float16sr, r); else; r = ul.(r); end # round residual to low precision
        d = FACT \ r # correction solve in low precision
        x = u.(x) + u.(d) # update in working precision

        if x_true !== nothing # track relative error if the true solution is known
            push!(error, norm(x - x_true) / norm(x_true))
        end

        if abs(norm(r)) < atol # stop once the residual is below the absolute tolerance
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

"""
    IR_GMRES(A, b, x0; SR=false, x_true=nothing,
             ul=Float16, u=Float32, ur=Float64, imax=10, atol=1e-12, GMRES_loop=10)

Iterative refinement where each correction step is solved by `GMRES_loop`
(restart length) iterations of GMRES in low precision. Works with dense
matrices and with `KronMatrix2` operators (A = A1 ⊗ A2), in which case only
the two small Kronecker factors are rounded and the full matrix is never formed.
"""
function IR_GMRES(A, b, x0; SR = false, x_true=nothing,
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12, GMRES_loop = 10)

    x = x0; error = ur[]  # initialize x and relative errors
    # low-precision copy of A (for Kronecker operators, round the two factors only)
    if A isa KronMatrix2 && SR
        a_sr = stochastic_round.(Float16sr, A.a[1])
        b_sr = stochastic_round.(Float16sr, A.b[1])
        A_low = KronMatrix2(a_sr, b_sr)
    elseif SR
        A_low = stochastic_round.(Float16sr, A)
    else
        A_low = ul.(A)
    end

    # relative error of the initial guess (iteration 0) to standardize the graphs
    if x_true !== nothing
        push!(error, norm(x - x_true) / norm(x_true))
    end

    for i in 1:imax
        r = ur.(b) - ur.(A) * ur.(x) # residual in high precision
        if SR; r = stochastic_round.(Float16sr, r); else; r = ul.(r); end # round residual to low precision
        d = GMRES(A_low, r, m = GMRES_loop) # correction solve: GMRES in low precision
        x = u.(x) + u.(d) # update in working precision

        if x_true !== nothing # track relative error if the true solution is known
            push!(error, norm(x - x_true) / norm(x_true))
        end

        if abs(norm(r)) < atol # stop once the residual is below the absolute tolerance
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

"""
    IR_kronreg(A, b, x0; SR=false, x_true=nothing, lambda=0.06,
               ul=Float16, u=Float32, ur=Float64, imax=10, atol=1e-12, GMRES_loop=10)

Iterative refinement with Tikhonov regularization, written as the augmented
(block) system

    [ I      A        ] [ r ]   [ b ]
    [ A'  -λ² I       ] [ x ] = [ 0 ]

so the regularized least-squares solution is the second block of the unknown.
The augmented matrix is never formed: products with it are computed by
`AugMatVecMult`, which keeps the Kronecker structure of `A`.
Assumes A is square (as in the PRblur problem), so both blocks have length(x0).

Change from the thesis version: the regularization parameter is now kept in
two copies, `lambda_r` (high precision, used in the residual) and `lambda_low`
(low precision, used in GMRES). Previously the low-precision value overwrote
`lambda` after the first iteration, so later residuals used the rounded λ.
Re-running may therefore differ slightly from the logs in `logs_graphs/`.
"""
function IR_kronreg(A, b, x0; SR = false, x_true=nothing, lambda = 0.06,
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12, GMRES_loop = 10)

    x = [zero(x0); x0]; b = vcat(b, zero(b)) # augmented unknown [r; x] and right-hand side [b; 0]
    error = ur[]  # initialize relative errors
    # low-precision copy of A (for Kronecker operators, round the two factors only)
    if A isa KronMatrix2 && SR
        a_sr = stochastic_round.(Float16sr, A.a[1])
        b_sr = stochastic_round.(Float16sr, A.b[1])
        A_low = KronMatrix2(a_sr, b_sr)
    elseif SR
        A_low = stochastic_round.(Float16sr, A)
    else
        A_low = ul.(A)
    end

    # regularization parameter: exact copy for the residual, rounded copy for the solve
    lambda_r = ur(lambda)
    if SR; lambda_low = stochastic_round(Float16sr, lambda); else; lambda_low = ul(lambda); end

    # relative error of the initial guess (iteration 0); x is the second block
    if x_true !== nothing
        push!(error, norm(x[length(x0)+1:end] - x_true) / norm(x_true))
    end

    for i in 1:imax
        r = ur.(b) - AugMatVecMult(ur.(A), lambda_r, ur.(x)) # residual in high precision
        if SR; r = stochastic_round.(Float16sr, r); else; r = ul.(r); end # round residual to low precision
        d = GMRES_kronreg(A_low, r, lambda_low, m = GMRES_loop) # correction solve in low precision
        x = u.(x) + u.(d) # update in working precision

        if x_true !== nothing # track relative error of the x block
            push!(error, norm(x[length(x0)+1:end] - x_true) / norm(x_true))
        end

        if abs(norm(r)) < atol # stop once the residual is below the absolute tolerance
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

#=
Exploratory GPU version (not used in the thesis). The custom Float16sr type
is not CUDA-compatible, so the final experiments ran on CPU nodes instead.

function IR_LU_GPU(A, b, x0; x_true=nothing, ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-5)

    A = cu(A)
    b = cu(b)
    x = cu(x0)
    error = Float64[]  # store relative errors
    A_low = ul.(A)
    LU = lu(A_low)
    LU = cu(LU)


    for i in 1:imax
        # Step 1: compute residual in ur precision
        r = b - A * x # GPU

        # Step 2: solve d in ul precision
        r_low = convert(CuArray{ul}, r)
        d = LU \ r_low # GPU calculate \ in fl32 or fl64, here in fl32 instead of fl16 as intended

        # Step 3: update x in u precision
        x = convert(CuArray{u}, x) + convert(CuArray{u}, d) # GPU

        # Step 4: compute relative error (if true solution given)
        if x_true !== nothing
            push!(error, norm(x - x_true) / norm(x_true))
        end

        if abs(norm(r)) < atol
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

=#

"""
    GMRES(A, r; tol=1e-8, m=10)

Minimal GMRES (no restarts) for A*d = r, starting from d = 0. Builds an
Arnoldi basis of at most `m` vectors with modified Gram-Schmidt and solves the
small Hessenberg least-squares problem at each step. All arithmetic runs in
the element type of `r`, so passing a low-precision residual keeps the whole
solve in low precision. Stops early when the residual drops below `tol * norm(r)`.
"""
function GMRES(A, r; tol = 1e-8, m = 10)
    n = length(r)
    T = eltype(r)

    # Krylov basis V in R^{n×(m+1)}, Hessenberg H in R^{(m+1)×m}
    V = zeros(T, n, m + 1)
    H = zeros(T, m + 1, m)

    # Step 0: normalize initial residual to get v1
    beta = norm(r)
    V[:, 1] = r / beta

    # RHS for small LS: minimize ||beta * e1 - H * y||
    e1 = zeros(T, m + 1)
    e1[1] = beta

    d = zeros(T, n)  # will hold current correction

    for j in 1:m

        # 1) Arnoldi step: w = A * v_j, then orthogonalize
        w = A * V[:, j]
        for i in 1:j
            H[i, j] = dot(V[:, i], w)
            w -= H[i, j] * V[:, i]
        end
        H[j + 1, j] = norm(w)

        # happy breakdown: Krylov subspace closed
        if H[j + 1, j] == 0
            break
        end

        V[:, j + 1] = w / H[j + 1, j]

        # 2) Solve small least-squares problem:
        #    minimize ||beta*e1 - H(1:j+1, 1:j) * y||
        y = H[1:j+1, 1:j] \ e1[1:j+1]

        # 3) Form GMRES approximation d_j = V_j * y
        d = V[:, 1:j] * y

        # 4) Check residual: A*d ≈ r
        if norm(r - A * d) ≤ tol * beta
            break
        end
    end
    return d
end

###

"""
    GMRES_kronreg(A, r, lambda; tol=1e-6, m=20)

Same as `GMRES`, but for the augmented Tikhonov system: every product with
the block matrix is computed by `AugMatVecMult(A, lambda, v)`.
"""
function GMRES_kronreg(A, r, lambda; tol = 1e-6, m = 20)
    n = length(r)
    T = eltype(r)

    # Krylov basis V in R^{n×(m+1)}, Hessenberg H in R^{(m+1)×m}
    V = zeros(T, n, m + 1)

    H = zeros(T, m + 1, m)

    # Step 0: normalize initial residual to get v1
    beta = norm(r)
    V[:, 1] = r / beta

    # RHS for small LS: minimize ||beta * e1 - H * y||
    e1 = zeros(T, m + 1)
    e1[1] = beta

    d = zeros(T, n)  # will hold current correction

    for j in 1:m
        # 1) Arnoldi step: w = A * v_j, then orthogonalize
        w = AugMatVecMult(A, lambda, V[:, j])
        for i in 1:j
            H[i, j] = dot(V[:, i], w)
            w -= H[i, j] * V[:, i]
        end
        H[j + 1, j] = norm(w)

        # happy breakdown: Krylov subspace closed
        if H[j + 1, j] == 0
            break
        end

        V[:, j + 1] = w / H[j + 1, j]

        # 2) Solve small least-squares problem:
        #    minimize ||beta*e1 - H(1:j+1, 1:j) * y||
        y = H[1:j+1, 1:j] \ e1[1:j+1]

        # 3) Form GMRES approximation d_j = V_j * y
        d = V[:, 1:j] * y

        # 4) Check residual: A*d ≈ r
        if norm(r - AugMatVecMult(A, lambda, d)) ≤ tol * beta
            break
        end
    end

    return d
end
