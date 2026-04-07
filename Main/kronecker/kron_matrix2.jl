mutable struct KronMatrix2{T}
    a::Vector{Matrix{T}}  # store matrices A{i}
    b::Vector{Matrix{T}}  # store matrices B{i}

    # No-arg constructor
    function KronMatrix2{T}() where {T}
        new{T}(Matrix{T}[], Matrix{T}[])
    end

    # Copy constructor
    function KronMatrix2{T}(G::KronMatrix2{T}) where {T}
        new{T}(copy(G.a), copy(G.b))
    end

    # Two-matrix constructor
    function KronMatrix2(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
        new{T}([A], [B])
    end

    # Define precision 
    function KronMatrix2(::Type{T}, A::AbstractMatrix, B::AbstractMatrix) where {T}
        new{T}([convert(Matrix{T}, A)], [convert(Matrix{T}, B)])
    end
    
end
