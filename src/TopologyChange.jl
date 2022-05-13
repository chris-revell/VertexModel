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

@views function topologyChange!(matrices)

    @unpack A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices = matrices

    # Find adjacency matrices from incidence matrices
    @.. thread=true Ā .= abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    dropzeros!(Ā)

    @.. thread=true B̄ .= abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)
    dropzeros!(B̄)

    # C .= B̄*Ā.÷2     # C adjacency matrix. Rows => cells; Columns => vertices (NB Integer division)
    mul!(C,B̄,Ā)
    @.. thread=true C .÷= 2
    dropzeros!(C)

    # Update transpose matrices
    Aᵀ .= sparse(transpose(A))
    dropzeros!(Aᵀ)
    Āᵀ .= abs.(Aᵀ)
    dropzeros!(Āᵀ)
    Bᵀ .= sparse(transpose(B))
    dropzeros!(Bᵀ)
    B̄ᵀ .= abs.(Bᵀ)
    dropzeros!(B̄ᵀ)

    # Calculate additional topology data
    # Number of edges around each cell found by summing columns of B̄
    cellEdgeCount .= sum.(eachrow(B̄))  # FastBroadcast doesn't work for this line; not sure why

    # Find boundary vertices
    # Summing each column of B finds boundary edges (for all other edges, cell orientations on either side cancel);
    # multiplying by Aᵀ gives nonzero values only where a vertex (row) has nonzero values at columns (edges) corresponding to nonzero values in the list of boundary edges.
    # Note that the abs is needed in case the direction of boundary edges cancel at a vertex
    boundaryVertices .= Āᵀ*abs.(sum.(eachcol(B))).÷2  # FastBroadcast doesn't work for this line; not sure why


    # Test for inconsistencies in the incidence matrices
    # test = B*A
    # dropzeros!(test)
    # length(findnz(test)[1]) > 0 ? throw() : nothing

    return nothing

end

export topologyChange!

end
