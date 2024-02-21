using MultivariateInterpolation
using Test
using SafeTestsets

@testset "MultivariateInterpolation.jl" begin
    @time @safetestset "multi index set" begin
        include("core/multi_index.jl")
    end
    @time @safetestset "grid" begin
        include("core/grid.jl")
    end
end
