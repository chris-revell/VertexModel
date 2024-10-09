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
@from "RotationMatrix.jl" using RotationMatrix

function model!(du, u, p, t)

    params, matrices = p
    @unpack A,
        B,
        Ā,
        B̄,
        cellTensions,
        cellPressures,
        cellPerpAxes,
        edgeLengths,
        edgeTangents,
        F,
        externalF,
        boundaryVertices,
        boundaryEdges,
        vertexAreas = matrices
    @unpack nVerts,
        nCells,
        nEdges,
        pressureExternal,
        peripheralTension,
        surfaceCentre,
        surfaceRadius = params

    spatialData!(u, params, matrices)

    fill!(F, @SVector zeros(3))
    dropzeros!(F)
    fill!(externalF, @SVector zeros(3))

    # peripheryLength = sum(boundaryEdges .* edgeLengths)

    for k = 1:nVerts
        for j in nzrange(A, k)
            for i in nzrange(B, rowvals(A)[j])
                # Force components from cell pressure perpendicular to edge tangents 
                F[k, rowvals(B)[i]] -= 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ(v=cellPerpAxes[rowvals(B)[i]]) * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] += cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure 
                # externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ(v=cellPerpAxes[rowvals(B)[i]]) * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
            end
            # Force on vertex from peripheral tension
            # externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
        end
        externalF[k] -= 100.0.*(norm(u[k] .- surfaceCentre)-surfaceRadius).*normalize(u[k].-surfaceCentre)
        du[k] = (sum(@view F[k, :]) .+ externalF[k]) ./ vertexAreas[k]
    end
    
    return du
end

export model!

end
