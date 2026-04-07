import LinearAlgebra: lu

function lu(K::KronMatrix2)
    FA = lu(K.a[1])
    FB = lu(K.b[1])

    P = KronMatrix2(Matrix(FA.P), Matrix(FB.P))
    L = KronMatrix2(FA.L, FB.L)
    U = KronMatrix2(FA.U, FB.U)

    return P, L, U
end
