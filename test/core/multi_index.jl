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
        groundtruth = MultiIndexSet(exponents, lp)
        @test isapprox(multi_index_set, groundtruth, atol=ATOL, rtol=RTOL)
    end

    @testset "utility" begin
        multi_index_set = MultiIndexSet(dim, deg, lp)
        @test spatial_dimension(multi_index_set) == dim
        @test polynomial_degree(multi_index_set, lp) == deg
        @test lp_degree(multi_index_set) == lp
    end

    @testset "contains_exponent" begin
        multi_index_set = MultiIndexSet(dim, deg, lp)
        @test contains_exponent(multi_index_set, zeros(dim))

        @testset "dim: $i" for i in 1:dim
            only_highest_deg = zeros(dim)
            only_highest_deg[i] = deg
            @test contains_exponent(multi_index_set, only_highest_deg)
        end

        @testset "lp: $lp2" for lp2 in LPS
            mi2 = MultiIndexSet(dim, deg, lp2)
            if lp2 <= lp || dim == 1 || deg == 0 || ((deg == 1) && lp2 < Inf)
                @test all(contains_exponent(multi_index_set, mi2))
            else
                @test !all(contains_exponent(multi_index_set, mi2))
            end
        end

        if lp >= 1
            @testset "deg: $deg2" for deg2 in 1:maximum(DEGS)
                mi2 = MultiIndexSet(dim, deg2, lp)
                if deg2 <= deg
                    @test all(contains_exponent(multi_index_set, mi2))
                else
                    @test !all(contains_exponent(multi_index_set, mi2))
                end
            end
        end
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
