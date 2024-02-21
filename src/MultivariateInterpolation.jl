module MultivariateInterpolation

export MultiIndexSet, spatial_dimension, polynomial_degree, lp_degree
export contains_exponent

export Grid

using LinearAlgebra

include("core/utils.jl")
include("core/multi_index.jl")
include("core/grid.jl")

end
