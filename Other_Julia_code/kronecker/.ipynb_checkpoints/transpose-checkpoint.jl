import Base: transpose

function transpose(Kin::KronMatrix2)
    A_t = [Matrix(transpose(A)) for A in Kin.a]
    B_t = [Matrix(transpose(B)) for B in Kin.b]
    return KronMatrix2(A_t[1], B_t[1]) # ensure input are in matrix form
end
