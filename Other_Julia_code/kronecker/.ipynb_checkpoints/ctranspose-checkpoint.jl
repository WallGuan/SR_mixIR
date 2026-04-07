import Base: adjoint

function adjoint(K::KronMatrix2)
    A_t = [Matrix(adjoint(A)) for A in K.a]
    B_t = [Matrix(adjoint(B)) for B in K.b]
    return KronMatrix2(A_t[1], B_t[1]) # ensure input is in matrix form
end
