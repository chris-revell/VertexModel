#
#  CalculateForce.jl
#  VertexModel
#
#  Created by Christopher Revell on 11/02/2021.
#
#
# Function to calculate force vectors on vertex k from cell i (Fᵢₖ) for all vertices.

module CalculateForce

# Julia packages
using LinearAlgebra
using StaticArrays
using LoopVectorization

@inline function calculateForce!(F,externalF,A,Ā,B,B̄,cellPressures,cellTensions,edgeTangents,edgeLengths,nVerts,nCells,nEdges,ϵ,pressureExternal,boundaryVertices)

    fill!(F,SVector{2}(zeros(2)))
    fill!(externalF,SVector{2}(zeros(2)))

    # Internal forces
    # NB This iteration could be improved to better leverage sparse arrays
    @turbo for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                F[k] += 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j]) + cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
                externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
            end
        end
    end

end

export calculateForce!

end
