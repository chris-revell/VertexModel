#
#  TopologyChange.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 01/02/2021.
#
#
# Function to recalculate derived matrices (Ā, Aᵀ etc.) for any change in vertex network topology (ie any change to A or B matrices).

module TopologyChange

# Julia packages
using LinearAlgebra

@inline @views function topologyChange!(A,Ā,Aᵀ,Āᵀ,B,B̄,Bᵀ,B̄ᵀ,C,cellEdgeCount,boundaryVertices,vertexEdges,nVerts)

    # Find adjacency matrices from incidence matrices
    Ā .= abs.(A)    # All -1 components converted to +1 (Adjacency matrix - vertices to edges)
    B̄ .= abs.(B)    # All -1 components converted to +1 (Adjacency matrix - cells to edges)
    C .= 0.5*B̄*Ā    # C adjacency matrix. Rows => cells; Columns => vertices

    # Establish transpose matrices
    Aᵀ .= transpose(A)
    Āᵀ .= abs.(Aᵀ)
    Bᵀ .= transpose(B)
    B̄ᵀ .= abs.(Bᵀ)

    # Calculate additional topology data
    cellEdgeCount    .= sum(B̄,dims=2)           # Number of edges around each cell found by summing columns of B̄
    boundaryVertices .= Āᵀ*abs.(sum(Bᵀ,dims=2)) # Find the vertices at the boundary

    # Use adjacency matrix Ā to find edges j that intersect vertex i
    for i=1:nVerts
        vertexEdges[i,:] .= findall(j->j!=0,Ā[:,i])
    end

    return nothing

end

export topologyChange!

end
