import LinearAlgebra: svd, Diagonal

function svd(K::KronMatrix2; full::Bool = true)
    Ua, Sa, Va = svd(K.a[1])
    Ub, Sb, Vb = svd(K.b[1])
    # change S from vector to matrix form
    Sa = Diagonal(Sa) 
    Sb = Diagonal(Sb)

    if full # reture all three
        U = KronMatrix2(Ua, Ub)
        S = KronMatrix2(Sa, Sb)
        V = KronMatrix2(Va, Vb)
        return U, S, V
    else # only s
        s = kron(Matrix(Sa), Matrix(Sb))
        s_vector = sort(vec(s); rev=true)
        return s_vector
    end
end
