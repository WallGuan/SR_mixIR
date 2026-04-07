module KronMatrixModule

import Base: adjoint, size, -, *, \
export KronMatrix2, adjoint, full, size, svd, transpose, -, *, \

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

end
