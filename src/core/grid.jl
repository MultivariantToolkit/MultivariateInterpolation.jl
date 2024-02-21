#######
# Grid functionality
#######

struct Grid{GV<:AbstractVector}
    nodes::GV # 1D generating nodes == generating values
    multi_index_set::MultiIndexSet
    function Grid(
        nodes::N, mi::M
    ) where {N<:AbstractVector{T},M<:MultiIndexSet} where {T<:Real}
        length(nodes) == polynomial_degree(mi) + 1 || throw(
            ArgumentError(
                """The number of nodes must be equal to the  polynomial_degree of the multi-index set plus one:\n
                length(nodes) == $(length(nodes))\n
                polynomial_degree(multi_index_set) == $(polynomial_degree(mi))
                """,
            ),
        )
        return new{N}(nodes, mi)
    end
end

function Grid(mi::MultiIndexSet)
    deg = polynomial_degree(mi)
    return Grid(leja_ordered_chebychev_2nd_order(deg + 1), mi)
end

Grid(dim, deg, lp) = Grid(MultiIndexSet(dim, deg, lp))

# delegation to multi_index_set
lp_degree(g::Grid) = lp_degree(g.multi_index_set)
spatial_dimension(g::Grid) = spatial_dimension(g.multi_index_set)
polynomial_degree(g::Grid) = polynomial_degree(g.multi_index_set)

@inline function _generating_points_view(nodes, dim, i)
    if dim % 2 == 0
        return nodes[i]
    else
        return -nodes[i]
    end
end

unisolvent_nodes(grid::Grid) = unisolvent_nodes(grid.nodes, grid.multi_index_set.data)

"""
    unisolvent_nodes(values::AbstractVector{T},exps) where {T}

TBW
"""
function unisolvent_nodes(
    values::AbstractVector{T}, exps::AbstractMatrix{TE}
) where {T,TE<:Int}
    dim, num_nodes = size(exps)
    nodes = Matrix{T}(undef, dim, num_nodes)

    for i in 1:num_nodes
        for j in 1:dim
            nodes[j, i] = _generating_points_view(values, j, exps[j, i] + 1)
        end
    end
    return nodes
end

"""
    chebychev_2nd_order(order::Int)

TBW
"""
function chebychev_2nd_order(order::Int)
    # consider using @memorize
    order >= 0 || throw(
        ArgumentError(
            "The order of Chebychev nodes need to be a positive non-zero integer. <$order> given.",
        ),
    )
    if order == 1
        return [0.0]
    end
    return cos.((0:(order - 1)) * pi / (order - 1))
end

"""
    leja_ordered_chebychev_2nd_order(order::Int)

TBW
"""
function leja_ordered_chebychev_2nd_order(order::Int)
    return leja!(chebychev_2nd_order(order))
end
