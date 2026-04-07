import Base: *

function *(K::KronMatrix2, M)
    if isa(M, AbstractArray) || isa(M, Number)
        return left_mtimes(K, M)
    elseif isa(M, KronMatrix2)
        if size(K.a[1], 2) == size(M.a[1], 1) && size(K.b[1], 2) == size(M.b[1], 1)
            Anew = K.a[1] * M.a[1]
            Bnew = K.b[1] * M.b[1]
            return KronMatrix2(Anew, Bnew)
        else
            error("Kron factors must be of compatible sizes for multiplication")
        end
    else
        error("Unsupported right-hand-side type for multiplication with KronMatrix2")
    end
end

function *(M, K::KronMatrix2)
    if isa(M, AbstractArray) || isa(M, Number)
        return right_mtimes(M, K)
    else
        error("Unsupported left-hand-side type for multiplication with KronMatrix2")
    end
end

function left_mtimes(K::KronMatrix2, x::AbstractArray)
    A = K.a[1]
    B = K.b[1]
    if isa(x, Number)
        y = KronMatrix2(A, x * B)
    elseif size(x, 1) == size(A, 2) * size(B, 2)
        T = eltype(A)
        y = zeros(T, size(A, 1) * size(B, 1), size(x, 2))
        for j in 1:size(x, 2)
            xj = reshape(x[:, j], size(B, 2), size(A, 2))
            y[:, j] = vec(B * xj * A')
        end
    else
        error("Dimension mismatch, incompatible array sizes")
    end
    return y
end


function right_mtimes(x::AbstractArray, K::KronMatrix2)
    A = K.a[1]
    B = K.b[1]
    y = transpose(left_mtimes(KronMatrix2(transpose(A), transpose(B)), transpose(x)))
    return y
end
