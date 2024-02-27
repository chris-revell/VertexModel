#
#  TopologyChange.jl
#  VertexModel
#
#  Created by Christopher Revell on 01/02/2021.
#
#
# Function to recalculate derived matrices (Ā, Aᵀ etc.) for any change in vertex network topology (ie any change to A or B matrices).

module TopologyChange

# Julia packages
using LinearAlgebra
using UnPack
using SparseArrays
using FastBroadcast
using FromFile
using DrWatson

# Local modules
@from "SenseCheck.jl" using SenseCheck
@from "OrderAroundCell.jl" using OrderAroundCell

function topologyChange!(matrices)

    @unpack A,
        B,
        Aᵀ,
        Ā,
        Āᵀ,
        Bᵀ,
        B̄,
        B̄ᵀ,
        C,
        cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        boundaryVertices,
        boundaryEdges = matrices

    # Find adjacency matrices from incidence matrices
    @.. thread=false Ā .= abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    @.. thread=false B̄ .= abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)

    # C adjacency matrix. Rows => cells; Columns => vertices. C .= B̄*Ā.÷2 (NB Integer division)
    C .= B̄*Ā./2

    # Update transpose matrices
    Aᵀ .= sparse(transpose(A))
    Āᵀ .= abs.(Aᵀ)
    Bᵀ .= sparse(transpose(B))
    B̄ᵀ .= abs.(Bᵀ)

    dropzeros!(A)
    dropzeros!(B)
    dropzeros!(C)
    dropzeros!(Ā)
    dropzeros!(B̄)
    dropzeros!(Aᵀ)
    dropzeros!(Āᵀ)
    dropzeros!(Bᵀ)
    dropzeros!(B̄ᵀ)

    # Calculate additional topology data
    # Number of edges around each cell found by summing columns of B̄
    cellEdgeCount .= sum.(eachrow(B̄))  # FastBroadcast doesn't work for this line; not sure why

    # Find boundary vertices
    # Summing each column of B finds boundary edges (for all other edges, cell orientations on either side cancel);
    # multiplying by Aᵀ gives nonzero values only where a vertex (row) has nonzero values at columns (edges) corresponding to nonzero values in the list of boundary edges.
    # Note that the abs is needed in case the direction of boundary edges cancel at a vertex
    boundaryVertices .= Āᵀ*abs.(sum.(eachcol(B))).÷2  # FastBroadcast doesn't work for this line; not sure why

    # Find list of edges at system periphery
    boundaryEdges .= abs.([sum(x) for x in eachcol(B)])

    for i=1:length(cellVertexOrders)
        cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(matrices,i)
    end

    return nothing

end

export topologyChange!

end
