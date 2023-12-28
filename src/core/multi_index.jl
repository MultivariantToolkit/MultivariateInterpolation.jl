"""
Here we store all functionality regarding the multi-index set. 

The initial version is based on `minterpy==0.2.0a0`
"""

struct MultiIndexSet{MT} <: AbstractArray{Int,2}
    data::MT
    MultiIndexSet(mat::MT) where {MT<:AbstractMatrix} = new{MT}(mat)
end

Base.size(mi::MultiIndexSet) = size(mi.data)
Base.getindex(mi::MultiIndexSet, i::Int) = getindex(mi.data, i)
Base.getindex(mi::MultiIndexSet, i::Int, j::Int) = getindex(mi.data, i, j)

@inline function _assert_valid_dimension(dim)
    return dim > 0 ||
        throw(ArgumentError("Dimension needs to be a positive integer! <$dim> given."))
end

@inline function _assert_valid_degree(deg)
    return deg >= 0 || throw(
        ArgumentError(
            "Polynomial degree needs to be a non-negative integer! <$deg> given."
        ),
    )
end

@inline function _assert_valid_lp(lp)
    return lp > zero(lp) ||
        throw(ArgumentError("The lp-degree needs to be a positive integer! <$lp> given."))
end

@inline function _assert_valid_multi_index_degree(dim, deg, lp)
    _assert_valid_dimension(dim)
    _assert_valid_degree(deg)
    return _assert_valid_lp(lp)
end

function _get_candidates(dim, deg)
    temp = (0:deg for el in 1:dim)
    return Iterators.product(temp...)
end

function _flat(x)
    n = length(first(x))
    return reshape(collect(Iterators.flatten(x)), :, n)
end

function _generate_exponents_static(dim, deg, lp)
    canidates_iter = _get_candidates(dim, deg)
    if lp < Inf
        result_iter = Iterators.filter(n -> norm(n, lp) <= deg, canidates_iter)
    else
        result_iter = canidates_iter
    end
    return _flat(result_iter)
end

function MultiIndexSet(dim::T, deg::T, lp::Tlp) where {T<:Integer,Tlp<:Real}
    _assert_valid_multi_index_degree(dim, deg, lp)
    return MultiIndexSet(_generate_exponents_static(dim, deg, lp))
end
