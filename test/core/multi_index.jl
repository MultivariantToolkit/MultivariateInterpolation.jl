using Random
using MultivariateInterpolation

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

DIMS = [1 3]
DIMS_FAIL = [0 -1]
DEGS = [0 1 4]
DEGS_FAIL = [-1]
LPS = [0.5 1.0 2.0 Inf]
LPS_FAIL = [0 -1]

PARAM_COMB = Iterators.product(DIMS, DEGS, LPS)
PARAM_COMB_FAIL = Iterators.product(DIMS_FAIL, DEGS_FAIL, LPS_FAIL)

@testset "$dim $deg $lp" for (dim, deg, lp) in PARAM_COMB
    @testset "from degree" begin
        multi_index_set = MultiIndexSet(dim, deg, lp)
        exponents = MultivariateInterpolation._generate_exponents_static(dim, deg, lp)
        groundtruth = MultiIndexSet(exponents)
        @test isapprox(multi_index_set, groundtruth)
    end
    @testset " failing input" begin
        @testset "$dim_fail $deg_fail $lp_fail" for (dim_fail, deg_fail, lp_fail) in
                                                    PARAM_COMB_FAIL
            @test_throws ArgumentError MultiIndexSet(dim_fail, deg, lp)
            @test_throws ArgumentError MultiIndexSet(dim_fail, deg_fail, lp)
            @test_throws ArgumentError MultiIndexSet(dim_fail, deg_fail, lp_fail)
            @test_throws ArgumentError MultiIndexSet(dim, deg_fail, lp)
            @test_throws ArgumentError MultiIndexSet(dim, deg_fail, lp_fail)
            @test_throws ArgumentError MultiIndexSet(dim, deg, lp_fail)
            @test_throws ArgumentError MultiIndexSet(dim_fail, deg, lp_fail)
        end
    end
end
