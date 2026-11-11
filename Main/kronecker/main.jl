#=
KronMatrixModule: Julia port of the MATLAB KronMatrix2 class
(see ../../KroneckerProductTools_matlab/).

A KronMatrix2 stores A = A1 ⊗ A2 by its two factors only. Products, solves,
transposes, LU and SVD work on the small factors, so an N² × N² operator costs
O(N²) memory instead of O(N⁴). Broadcasting a number type over it
(e.g. `Float16.(K)`) converts both factors (see precision.jl).
=#
module KronMatrixModule

import Base: adjoint, size, -, *, \
import Base.Broadcast: broadcasted
export KronMatrix2, adjoint, full, size, svd, lu, transpose, -, *, \

# Include the parts of the module
include("kron_matrix2.jl")
include("ctranspose.jl")
include("full.jl")
include("mdivide.jl")
include("mtimes.jl")
include("size.jl")
include("svd.jl")
include("transpose.jl")
include("uminus.jl")
include("lu.jl")
include("precision.jl")

end
