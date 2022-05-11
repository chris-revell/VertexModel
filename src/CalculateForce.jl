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
using SparseArrays

function calculateForce!(R,params,matrices)

    @unpack A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,externalF,ϵ,boundaryVertices = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal = params

    fill!(F,@SVector zeros(2))
    fill!(externalF,@SVector zeros(2))

    # Internal forces
    for k=1:nVerts
        for j in nzrange(A,k)
            for i in nzrange(B,rowvals(A)[j])
                # display("rowvals(B)[i]")
                # display(rowvals(B)[i])
                # display("rowvals(A)[j]")
                # display(rowvals(A)[j])
                # display("k")
                # display(k)
                F[k,rowvals(B)[i]] += 0.5*cellPressures[rowvals(B)[i]]*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k]*(ϵ*edgeTangents[rowvals(A)[j]])
                F[k,rowvals(B)[i]] += cellTensions[rowvals(B)[i]]*B̄[rowvals(B)[i],rowvals(A)[j]]*A[rowvals(A)[j],k]*edgeTangents[rowvals(A)[j]]/edgeLengths[rowvals(A)[j]]
                externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k]*(ϵ*edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
        # for i=1:nCells
        #     for j=1:nEdges
        #         # Pressure term
        #         F[k,i] += 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])
        #         # Tension term
        #         F[k,i] += cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
        #         # External force
        #         externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j])) # 0 unless boundaryVertices != 0
            end
        end
    end

end

export calculateForce!

end
