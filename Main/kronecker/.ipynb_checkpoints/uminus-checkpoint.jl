import Base: -

function -(Kin::KronMatrix2)
    A_neg = -Kin.a[1]
    B_copy = Kin.b[1]
    return KronMatrix2(A_neg, B_copy)
end
