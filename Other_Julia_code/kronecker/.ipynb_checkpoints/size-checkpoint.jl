import Base: size

# Core size function
function size(K::KronMatrix2)
    A1, B1 = K.a[1], K.b[1]
    s = (size(A1, 1) * size(B1, 1), size(A1, 2) * size(B1, 2))

    for i in 2:length(K.a)
        Ai, Bi = K.a[i], K.b[i]
        sc = (size(Ai, 1) * size(Bi, 1), size(Ai, 2) * size(Bi, 2))
        if sc != s
            error("Inconsistent term sizes in KronMatrix2 object.")
        end
    end

    return s
end

# Dimension-based size
function size(K::KronMatrix2, dim::Int)
    s = size(K)
    if dim == 1 || dim == 2
        return s[dim]
    else
        error("Dimension must be 1 or 2.")
    end
end

# Simulate MATLAB-like dynamic outputs
function size_with_outputs(K::KronMatrix2, outputs::Int; dim::Union{Int,Nothing}=nothing)
    s = dim === nothing ? size(K) : (size(K, dim),)
    if outputs == 0
        println("\nans =\n")
        println(s)
    elseif outputs == 1
        return s
    elseif outputs == 2
        if dim !== nothing
            error("Too many output arguments.")
        else
            return s[1], s[2]
        end
    else
        error("Too many output arguments.")
    end
end