#
#  CalculateForce.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 11/02/2021.
#
#
# Function to calculate force vectors on vertex k from cell i (Fᵢₖ) for all vertices.

module CalculateForce

# Julia packages
using LinearAlgebra

@inline @views function calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ)

    fill!(F,0.0)
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                F[k,i,:] .+= 0.5*cellPressures[i]*B[i,j]*Ā[j,k].*(ϵ*edgeTangents[j,:]) .+ cellTensions[i]*B̄[i,j]*A[j,k].*edgeTangents[j,:]./edgeLengths[j]
            end
        end
    end

end

export calculateForce!

end
