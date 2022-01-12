#
#  CalculateForce.jl
#  VertexModel
#
#  Created by Christopher Revell on 11/02/2021.
#
#
# Function to calculate force vector on vertex k from cell i (Fᵢₖ) for all vertices.

module CalculateForce

# Julia packages
using LinearAlgebra
using StaticArrays
using LoopVectorization
using UnPack
using Plots
using GeometryBasics

function calculateForce!(R,params,matrices)

    @unpack A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ,boundaryVertices = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal = params

    fill!(F,@SVector zeros(2))

    # Internal forces
    # NB This iteration could be improved to better leverage sparse arrays
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                F[k] += 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j]) + cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
                # if isnan(F[k][1]) || isnan(F[k][2])
                #     println("cellPressures[i], $(cellPressures[i])")
                #     println("B[i,j], $(B[i,j])")
                #     println("Ā[j,k], $(Ā[j,k])")
                #     println("edgeTangents[j], $(edgeTangents[j])")
                #     println("cellTensions[i], $(cellTensions[i])")
                #     println("B̄[i,j], $(B̄[i,j])")
                #     println("A[j,k], $(A[j,k])")
                #     println("edgeTangents[j], $(edgeTangents[j])")
                #     println("edgeLengths[j], $(edgeLengths[j])")
                # end
                #externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
                F[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
            end
        end
    end

end

export calculateForce!

end
