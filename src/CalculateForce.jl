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
using UnPack
using Plots
using GeometryBasics

function calculateForce!(R,params,matrices)

    @unpack A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,externalF,ϵ,boundaryVertices = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal = params

    fill!(F,@SVector zeros(2))
    fill!(externalF,@SVector zeros(2))

    # Internal forces
    # NB This iteration could probably be improved to better leverage sparse arrays
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                # Pressure term
                F[k,i] += 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])
                # Tension term
                F[k,i] += cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
                # External force
                externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
            end
        end
    end

end

export calculateForce!

end
