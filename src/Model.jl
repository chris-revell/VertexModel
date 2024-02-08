#
#  Model.jl
#  VertexModel
#
#  Created by Christopher Revell on 15/02/2023.
#
#
# Function to calculate force vector on vertex k from cell i (Fᵢₖ) for all vertices and update vertex positions.

module Model

# Julia packages
using LinearAlgebra
using StaticArrays
using UnPack
using SparseArrays
using .Threads
using FromFile 
using DrWatson

# Local modules
@from "SpatialData.jl" using SpatialData
@from "Laplacians.jl" using Laplacians

function model!(du, u, p, t)

    R_i, params, matrices = p
    @unpack A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,externalF,ϵ,boundaryVertices,boundaryEdges,vertexAreas = matrices
    @unpack nVerts,nCells,nEdges,pressureExternal,peripheralTension, tMax = params

    spatialData!(u,params,matrices)

    fill!(F,@SVector zeros(2))
    fill!(externalF,@SVector zeros(2))

    peripheryLength = sum(boundaryEdges.*edgeLengths)

    #stretch monolayer, map R_x->(1 + \lambda)R_x, Ry->R-y/(1+\lambda)
    λ=0.2
    Λ=@SMatrix[
        1+λ 0.0
        0.0 (1+λ)
    ]

    # Λ=@SMatrix[
    #    1/(1+λ) 0.0
    #     0.0 (1+λ)
    # ]

    stretch=Λ-I(2)
    
    for k=1:nVerts
        for j in nzrange(A,k)
            for i in nzrange(B,rowvals(A)[j])
                # Force components from cell pressure perpendicular to edge tangents 
                F[k,rowvals(B)[i]] += 0.5*cellPressures[rowvals(B)[i]]*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k].*(ϵ*edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k,rowvals(B)[i]] += cellTensions[rowvals(B)[i]]*B̄[rowvals(B)[i],rowvals(A)[j]]*A[rowvals(A)[j],k].*edgeTangents[rowvals(A)[j]]./edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure 
                externalF[k] += boundaryVertices[k]*(0.5*pressureExternal*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k].*(ϵ*edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
            end
            # Force on vertex from peripheral tension
            externalF[k] -= boundaryEdges[rowvals(A)[j]]*peripheralTension*(peripheryLength-sqrt(π*nCells))*A[rowvals(A)[j],k].*edgeTangents[rowvals(A)[j]]./edgeLengths[rowvals(A)[j]]
        end
    end

    du .=((sum.(eachrow(matrices.F)).+externalF)./(100.0.*vertexAreas)) .+  ([stretch*x for x in R_i] )./(tMax)
    
end

function modeltest!(du, u, p, t)

    params, matrices = p
    @unpack vertexAreas, g_vec = matrices
    #@unpack nVerts,nCells,nEdges,pressureExternal,peripheralTension = params

    spatialData!(u,params,matrices)

    M=makeM(matrices)

    du .= -(M'*g_vec)./(100.0.*vertexAreas)
end

export model!, modeltest!

end
