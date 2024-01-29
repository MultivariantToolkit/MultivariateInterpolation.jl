"""
Here we store all functionality regarding the multi-index set. 

The initial version is based on `minterpy==0.2.0a0`
"""

struct MultiIndexSet{MT,Tlp} <: AbstractArray{Int,2}
    data::MT
    lp::Tlp
    function MultiIndexSet(mat::MT, lp::Tlp) where {MT<:AbstractMatrix,Tlp<:Real}
        _assert_valid_lp(lp)
        return new{MT,Tlp}(mat, lp)
    end
end

"""
    lp_degree(m::MultiIndexSet)

TBW
"""
function lp_degree(m::MultiIndexSet)
    return m.lp
end

# array interface
Base.size(mi::MultiIndexSet) = size(mi.data)
Base.getindex(mi::MultiIndexSet, i::Int) = getindex(mi.data, i)
Base.getindex(mi::MultiIndexSet, i::Int, j::Int) = getindex(mi.data, i, j)

@inline function _assert_valid_dimension(dim)
    dim > 0 ||
        throw(ArgumentError("Dimension needs to be a positive integer! <$dim> given."))
    return nothing
end

@inline function _assert_valid_degree(deg)
    deg >= 0 || throw(
        ArgumentError(
            "Polynomial degree needs to be a non-negative integer! <$deg> given."
        ),
    )
    return nothing
end

@inline function _assert_valid_lp(lp)
    lp > zero(lp) ||
        throw(ArgumentError("The lp-degree needs to be a positive integer! <$lp> given."))
    return nothing
end

@inline function _assert_valid_multi_index_degree(dim, deg)
    _assert_valid_dimension(dim)
    _assert_valid_degree(deg)
    return nothing
end

function _get_candidates(dim, deg)
    temp = (0:deg for el in 1:dim)
    return Iterators.product(temp...)
end

function _flat(x)
    n = length(first(x))
    return reshape(collect(Iterators.flatten(x)), n, :)
end

function _generate_exponents_static(dim, deg, lp)
    candidates_iter = _get_candidates(dim, deg)
    if lp < Inf
        result_iter = Iterators.filter(n -> norm(n, lp) <= deg, candidates_iter)
    else
        result_iter = candidates_iter
    end
    return _flat(result_iter)
end

"""
    MultiIndexSet(dim::T, deg::T, lp::Tlp) where {T<:Integer,Tlp<:Real}

TBW
"""
function MultiIndexSet(dim::T, deg::T, lp::Tlp) where {T<:Integer,Tlp<:Real}
    _assert_valid_multi_index_degree(dim, deg)
    return MultiIndexSet(_generate_exponents_static(dim, deg, lp), lp)
end

### utility

"""
    spatial_dimension(m::MultiIndexSet)

TBW
"""
function spatial_dimension(m::MultiIndexSet)
    return size(m, 1)
end

"""
    polynomial_degree(multi_index_set,lp)

TBW
"""
function polynomial_degree(multi_index_set, lp)
    norms = norm.(eachcol(multi_index_set), lp)
    return ceil(maximum(norms))
end

"""
    contains_exponent(mi::MultiIndexSet, x::AbstractVector{T}) where {T<:Real}

TBW
"""
function contains_exponent(mi::MultiIndexSet, x::AbstractVector{T}) where {T<:Real}
    if length(x) != spatial_dimension(mi)
        return false
    end
    for col in eachcol(mi)
        if isapprox(col, x)
            return true
        end
    end
    return false
end

"""
    contains_exponent(mi::MultiIndexSet, x::AbstractMatrix{T}) where {T<:Real}

TBW
"""
function contains_exponent(mi::MultiIndexSet, x::AbstractMatrix{T}) where {T<:Real}
    res = zeros(Bool, size(x, 2))
    for i in 1:size(x, 2)
        if contains_exponent(mi, view(x, :, i))
            res[i] = true
        end
    end
    return res
end
