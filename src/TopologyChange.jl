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

@views function topologyChange!(matrices)

    @unpack A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices = matrices

    # Find adjacency matrices from incidence matrices
    Ā .= abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    dropzeros!(Ā)
    B̄ .= abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)
    dropzeros!(B̄)
    C .= B̄*Ā.÷2     # C adjacency matrix. Rows => cells; Columns => vertices (NB Integer division)
    dropzeros!(C)

    # Update transpose matrices
    Aᵀ .= transpose(A)
    dropzeros!(Aᵀ)
    Āᵀ .= abs.(Aᵀ)
    dropzeros!(Āᵀ)
    Bᵀ .= transpose(B)
    dropzeros!(Bᵀ)
    B̄ᵀ .= abs.(Bᵀ)
    dropzeros!(B̄ᵀ)

    # Calculate additional topology data
    cellEdgeCount    .= sum(B̄,dims=2)[:,1]           # Number of edges around each cell found by summing columns of B̄
    boundaryVertices .= Āᵀ*abs.(sum(Bᵀ,dims=2))[:,1] # Find the vertices at the boundary

    # Test for inconsistencies in the incidence matrices
    test = B*A
    dropzeros!(test)
    length(findnz(test)[1]) > 0 ? throw() : nothing

    return nothing

end

export topologyChange!

end
