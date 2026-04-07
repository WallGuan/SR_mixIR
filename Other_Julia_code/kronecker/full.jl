function full(K::KronMatrix2)
    A = kron(K.a[1], K.b[1])
    for i in 2:length(K.a)
        A += kron(K.a[i], K.b[i])
    end
    return A
end
