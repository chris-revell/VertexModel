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

@views function topologyChange!(A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices)

    # @unpack A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices = matrices

    # Find adjacency matrices from incidence matrices
    Ā .= abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    B̄ .= abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)
    C .= B̄*Ā.÷2     # C adjacency matrix. Rows => cells; Columns => vertices (NB Integer division)

    # Update transpose matrices
    Aᵀ .= transpose(A)
    Āᵀ .= abs.(Aᵀ)
    Bᵀ .= transpose(B)
    B̄ᵀ .= abs.(Bᵀ)

    # Calculate additional topology data
    cellEdgeCount    .= sum(B̄,dims=2)[:,1]           # Number of edges around each cell found by summing columns of B̄
    boundaryVertices .= Āᵀ*abs.(sum(Bᵀ,dims=2))[:,1] # Find the vertices at the boundary

    # Use adjacency matrix Ā to find edges j that intersect vertex i,
    # and arrange in polar angle order (ordering will not change
    # without topology change, even if polar angles themselves do change)
    # for i=1:nVerts
    #     # Currently vertexEdges skips boundary vertices, since these are not used in torque calculations, and can have 2 edges rather than 3
    #     if boundaryVertices[i] == 0
    #         vertexEdges[i,:] .= findall(j->j!=0,Ā[:,i])
    #         polarAngles = zeros(3)
    #         for k = 1:3
    #             # Multiply edge tangent vectors by corresponding incidence matrix component to ensure we consider all vectors point into the vertex
    #             polarAngles[k] = atan(A[vertexEdges[i,k],k]*edgeTangents[vertexEdges[i,k],1],A[vertexEdges[i,k],k]*edgeTangents[vertexEdges[i,k], 2])
    #         end
    #         # Find order of indices around vertex.
    #         orderAroundVertex = sortperm(polarAngles)
    #         vertexEdges[i,:] .= vertexEdges[i,orderAroundVertex]
    #     end
    # end

    return nothing

end

export topologyChange!

end
