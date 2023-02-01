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
using .Threads

function calculateForce!(R,params,matrices)

    @unpack A,Aᵀ,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,externalF,ϵ,boundaryVertices,boundaryEdges = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal,peripheralTension = params

    fill!(F,@SVector zeros(2))
    fill!(externalF,@SVector zeros(2))

    # Internal forces
    for k=1:nVerts
        for j in nzrange(A,k)
            for i in nzrange(B,rowvals(A)[j])
                F[k,rowvals(B)[i]] += 0.5*cellPressures[rowvals(B)[i]]*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k].*(ϵ*edgeTangents[rowvals(A)[j]])
                F[k,rowvals(B)[i]] += cellTensions[rowvals(B)[i]]*B̄[rowvals(B)[i],rowvals(A)[j]]*A[rowvals(A)[j],k].*edgeTangents[rowvals(A)[j]]./edgeLengths[rowvals(A)[j]]
                externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k].*(ϵ*edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
            end
        end
    end

    peripheryLength = sum(boundaryEdges.*edgeLengths)    
    for j in 1:nEdges #findall(x->x!=0,boundaryEdges)
        for k in nzrange(Aᵀ,j)            
            externalF[rowvals(Aᵀ)[k]] -= boundaryEdges[j]*peripheralTension*(peripheryLength-sqrt(π*nCells))*Aᵀ[rowvals(Aᵀ)[k],j].*edgeTangents[j]./edgeLengths[j]
        end
    end

end

export calculateForce!

end
