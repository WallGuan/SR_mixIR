import Base: \

function \(K::KronMatrix2, y::AbstractMatrix)
    A = K.a[1]
    B = K.b[1]

    ma, na = size(A)
    mb, nb = size(B)

    Y = reshape(y, mb, ma)
    Z = A \ transpose(Y)
    X = B \ transpose(Z)
    x = reshape(X, nb * na)
    return x
end
