import Base.Broadcast: broadcasted

function broadcasted(::Type{T}, K::KronMatrix2) where {T}
    aT = K.a[1]
    bT = K.b[1]
    return KronMatrix2(T, aT, bT)
end