#=
other_funcs.jl

Helper functions for the experiments: test-matrix generators, the Tikhonov
augmented system, noise, and plotting of error histories.
=#

"""
    ill_A(m, n, cond)

Random m×n test matrix with a prescribed condition number `cond`: random
orthogonal U and V, and singular values spaced evenly on a log scale between
σ_max and σ_max/cond.
"""
function ill_A(m, n, cond)
    U, _ = qr(randn(m, m))
    V, _ = qr(randn(n, n))

    k = min(m, n) # num 
    sigma_max = 10^(rand() * 4 - 2) # generate 
    sigma_min = sigma_max / cond
    s = exp.(LinRange(log(sigma_max), log(sigma_min), k)) # exponentially spaced
    S = zeros(m, n)
    for i in 1:k # rectangle S
        S[i, i] = s[i]
    end

    A = U * S * V'

    return A
end

"""
    error_analysis(error; logscale=true, saveplot=true, savepath="error_plot.png")

Print the minimum relative error and where it occurs, and optionally save a
convergence plot. Note `min_iter` is a 1-based index: index 1 is the initial
guess (iteration 0).
"""
function error_analysis(error; logscale=true, saveplot=true, savepath="error_plot.png") # fixed typo: default was `ture`
    min_error = minimum(error)
    min_iter = argmin(error)
    println("Minimum relative error = $(min_error) at iteration $(min_iter)")

    # Plot
    if saveplot
        p = plot(
            error;
            yscale = logscale ? :log10 : :identity,
            xlabel = "Iteration",
            ylabel = "Relative Error",
            title = logscale ? "Error Convergence (log scale)" : "Error Convergence (linear scale)",
            legend = false,
            lw = 2,
            marker = :circle,
        )

        # Mark the minimum point
        scatter!(p, [min_iter], [min_error], color = :red, label = "Min Error")

        # Save or display
        savefig(p, savepath)
    end

    return (min_error = min_error, min_iter = min_iter)
end

"""
    AUG(A, b, x, lamda)

Explicitly form the Tikhonov augmented system K = [I A; A' -λ²I] with
right-hand side [b; 0] and reference solution [0; x]. Used for the small
inverse heat problem; the large PRblur problem uses `AugMatVecMult` instead,
so K is never formed.
"""
function AUG(A, b, x, lamda)
    m, n = size(A)
    K = [ I(m)   A; A'   -lamda^2 * I(n) ]
    rhs = [ b; zeros(n) ]
    sol = [ zeros(n); x ]
    println("Original condition num: ", cond(A))
    println("Regularize condition num: ", cond(K))
    return K, rhs, sol
end

"""
    AugMatVecMult(A, alpha, v)

Matrix-free product z = K*v with K = [I A; A' -α²I], where v = [v1; v2].
Only products with A and A' are needed, so a `KronMatrix2` operator keeps
its Kronecker structure. The result has the element type of `v`.
"""
function AugMatVecMult(A, alpha, v)
    T = eltype(v)
    alpha2 = T(alpha^2)
    m, n = size(A)

    z = zeros(T, m + n)

    z[1:m]       = v[1:m] + A * v[m+1:m+n]
    z[m+1:m+n]   = A' * v[1:m] - alpha2 * v[m+1:m+n]

    return z
end

"""
    bnoise(b; level=1e-3)

Add white Gaussian noise scaled so that norm(noise) / norm(b) = `level`
(1e-3, 1e-2, 1e-1 = the small, medium and large noise cases).
"""
function bnoise(b; level=1e-3)
    r = randn(eltype(b), length(b))
    noise = (level * norm(b) / norm(r)) * r
    return b + noise
end

"""
    heat(n, kappa=1.0)

Inverse heat equation test problem (Volterra integral equation), following
`heat` from Hansen's Regularization Tools. Returns (A, b, x_true) with A a
lower-triangular Toeplitz matrix; κ controls how ill-conditioned it is.

Note: here the quadrature points are `range(h/2, stop=1, length=n)`, whereas
Regularization Tools uses `t = h/2 : h : 1` (spacing exactly h, ending at 1 - h/2).
The two grids differ slightly, so A is close to, but not identical to, the
MATLAB version. Kept as is so the thesis results stay reproducible.
"""
function heat(n::Int, kappa::Float64=1.0)
    h = 1 / n
    t = range(h/2, stop=1, length=n)
    c = h / (2 * kappa * sqrt(pi))
    d = 1 / (4 * kappa^2)

    k = c .* t.^(-1.5) .* exp.(-d ./ t)
    r = zeros(Float64, n)
    r[1] = k[1]
    A = Toeplitz(k, r) |> Matrix

    x = zeros(Float64, n)
    for i in 1:(n ÷ 2)
        ti = i * 20 / n
        if ti < 2
            x[i] = 0.75 * ti^2 / 4
        elseif ti < 3
            x[i] = 0.75 + (ti - 2) * (3 - ti)
        else
            x[i] = 0.75 * exp(-(ti - 3) * 2)
        end
    end

    b = A * x
    return A, b, x
end
