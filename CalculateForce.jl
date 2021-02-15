#
#  CalculateForce.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 11/02/2021.
#
#
#

module CalculateForce

# Julia packages
using LinearAlgebra

# Local modules

@inline @views function calculateForce!(F,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges)

    ϵᵢ = [0.0 1.0
          -1.0 0.0]

    fill!(F,0.0)
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                F[k,:] .+= 0.5*cellPressures[i]*B[i,j]*Ā[j,k].*(ϵᵢ*edgeTangents[j,:]) .+ cellTensions[i]*B̄[i,j]*A[j,k].*edgeTangents[j,:]./edgeLengths[j]
            end
        end
    end


end

export calculateForce!

end
