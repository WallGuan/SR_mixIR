function IR_fact(A, b, x0; SR = true, fact = lu, x_true=nothing, 
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12)

    x = x0; error = ur[]  # initialize x and relative errors
    if SR; A_low = stochastic_round.(Float16sr, A); else; A_low = ul.(A); end # check if SR
    FACT = fact(A_low) # fact A

    # compute relative error from initial guess to standardize the graph
    if x_true !== nothing # compute error if x_true exist
        push!(error, norm(x - x_true) / norm(x_true))
    end

    for i in 1:imax
        r = ur.(b) - ur.(A) * ur.(x) # residual 
        if SR; r = stochastic_round.(Float16sr, r); else; r = ul.(r); end # check if SR
        d = FACT \ r # solve d
        x = u.(x) + u.(d) # update x

        if x_true !== nothing # compute error of x_true exist
            push!(error, norm(x - x_true) / norm(x_true))
        end
        
        if abs(norm(r)) < atol # count iterations
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

function IR_GMRES(A, b, x0; SR = false, x_true=nothing, 
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12, GMRES_loop = 10)

    x = x0; error = ur[]  # initialize x and relative errors
    if A isa KronMatrix2 && SR
        a_sr = stochastic_round.(Float16sr, A.a[1])
        b_sr = stochastic_round.(Float16sr, A.b[1])
        A_low = KronMatrix2(a_sr, b_sr)
    elseif SR 
        A_low = stochastic_round.(Float16sr, A) 
    else
        A_low = ul.(A)
    end

    # compute relative error from initial guess to standardize the graph
    if x_true !== nothing # compute error if x_true exist
        push!(error, norm(x - x_true) / norm(x_true))
    end

    for i in 1:imax
        r = ur.(b) - ur.(A) * ur.(x) # residual 
        if SR; r = stochastic_round.(Float16sr, r); else; r = ul.(r); end # check if SR
        d = GMRES(A_low, r, m = GMRES_loop) # GMRES
        x = u.(x) + u.(d) # update x

        if x_true !== nothing # compute error of x_true exist
            push!(error, norm(x - x_true) / norm(x_true))
        end
        
        if abs(norm(r)) < atol # count iterations
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

function IR_kronreg(A, b, x0; SR = false, x_true=nothing, lambda = 0.06, 
    ul = Float16, u = Float32, ur = Float64, imax=10, atol=1e-12, GMRES_loop = 10)

    x = [zero(x0); x0]; b = vcat(b, zero(b)) # initialize x and b in augmented size
    error = ur[]  # initialize relative errors
    if A isa KronMatrix2 && SR
        a_sr = stochastic_round.(Float16sr, A.a[1])
        b_sr = stochastic_round.(Float16sr, A.b[1])
        A_low = KronMatrix2(a_sr, b_sr)
    elseif SR 
        A_low = stochastic_round.(Float16sr, A) 
    else
        A_low = ul.(A)
    end

    # compute relative error from initial guess to standardize the graph
    if x_true !== nothing # compute error if x_true exist
        push!(error, norm(x[length(x0)+1:end] - x_true) / norm(x_true)) 
    end
    
    for i in 1:imax
        r = ur.(b) - AugMatVecMult(ur.(A), ur.(lambda), ur.(x)) # residual 
        if SR; r = stochastic_round.(Float16sr, r); lambda = stochastic_round.(Float16sr, lambda); 
        else; r = ul.(r); lambda = ul.(lambda); end # check if SR
        d = GMRES_kronreg(A_low, r, lambda, m = GMRES_loop) # 
        x = u.(x) + u.(d) # update x

        if x_true !== nothing # compute error of x_true exist
            push!(error, norm(x[length(x0)+1:end] - x_true) / norm(x_true)) 
        end
        
        if abs(norm(r)) < atol # count iterations
            println("Converged at iter $(i)")
            break
        end
    end

    return x, error
end

###

#= 

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
