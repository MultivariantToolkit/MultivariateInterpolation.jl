using Random
using MultivariateInterpolation
using FastChebInterp

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

ORDER = [1, rand(RNG, 2:10)]
DIMS = [1 3]
DEGS = [0 1 4]
LPS = [0.5 1.0 2.0 Inf]

PARAM_COMB = Iterators.product(DIMS, DEGS, LPS)

@testset "leja ordering" begin
    @testset "order: $n" for n in ORDER
        leja_points = MultivariateInterpolation.leja_ordered_chebychev_2nd_order(n)
        cheb_points = MultivariateInterpolation.chebychev_2nd_order(n)
        @test isapprox(sort(leja_points), sort(cheb_points))
    end
end

@testset "chebyshev nodes" begin
    @testset "order: $n" for n in ORDER
        test_nodes = MultivariateInterpolation.chebychev_2nd_order(n)
        groundtruth = chebpoints(n - 1, -1, 1)

        @test isapprox(test_nodes, groundtruth)
    end
end

@testset "$dim $deg $lp" for (dim, deg, lp) in PARAM_COMB
    mi = MultiIndexSet(dim, deg, lp)
    default_nodes = MultivariateInterpolation.leja_ordered_chebychev_2nd_order(deg + 1)

    @testset "constructor custom nodes" begin
        custom_nodes = rand(RNG, deg + 1)
        test_grid_custom_nodes = Grid(custom_nodes, mi)

        @test isapprox(test_grid_custom_nodes.multi_index_set, mi)
        @test isapprox(test_grid_custom_nodes.nodes, custom_nodes)
    end

    @testset "constructor default nodes" begin
        test_grid_default_nodes = Grid(mi)

        @test isapprox(test_grid_default_nodes.multi_index_set, mi)
        @test isapprox(test_grid_default_nodes.nodes, default_nodes)
    end

    @testset "constructor from degree" begin
        test_grid_from_degree = Grid(dim, deg, lp)

        @test isapprox(test_grid_from_degree.multi_index_set, mi)
        @test isapprox(test_grid_from_degree.nodes, default_nodes)
    end

    @testset "constructor fail" begin
        nodes_fail = rand(RNG, deg)
        @test_throws ArgumentError Grid(nodes_fail, mi)
    end

    @testset "delegations" begin
        test_grid = Grid(mi)
        @test spatial_dimension(test_grid) == dim
        @test polynomial_degree(test_grid) == deg
        @test lp_degree(test_grid) == lp
    end
end
