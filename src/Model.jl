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
@from "Stretch.jl" using Stretch

function model!(du, u, p, t)

    params, matrices = p
    @unpack A,
        B,
        Ā,
        B̄,
        cellTensions,
        cellPressures,
        edgeLengths,
        edgeTangents,
        F,
        externalF,
        ϵ,
        boundaryVertices,
        boundaryEdges,
        vertexAreas = matrices
    @unpack nVerts,
        nCells,
        nEdges,
        pressureExternal,
        peripheralTension,
        vertexWeighting,
        stretchType,
        tStretch,
        κ = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    dR = reinterpret(SVector{2,Float64}, du)

    spatialData!(R, params, matrices)

    fill!(F, @SVector zeros(2))
    dropzeros!(F)
    fill!(externalF, @SVector zeros(2))

    peripheryLength = sum(boundaryEdges .* edgeLengths)

    for k = 1:nVerts
        for j in nzrange(A, k)
            for i in nzrange(B, rowvals(A)[j])
                # Force components from cell pressure perpendicular to edge tangents 
                F[k, rowvals(B)[i]] += 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure 
                externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
            end
            # Force on vertex from peripheral tension
            externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
        end
        
        dR[k] = (sum(@view F[k, :]) .+ externalF[k])
    end
    #if weighting drag by vertex area divide by vertex areas
    vertexWeighting && (dR ./= vertexAreas)
    # dR accesses the same underlying data as du, so by altering dR we have already updated du appropriately

    if stretchType != "none"
        Rt, R_final=stretchCells(R,t, params, matrices)
        #κ=1 #spring const assuming that the vertices are anchored to the stretched membrane layer by springs

        if t<= tStretch
            #dR .+= dR .-κ.*(R .- Rt)
            dR .+= -κ.*(R .- Rt)

        else
            dR .+=-κ.*(R .- R_final)
        end
    end


    return du
end

export model!

end
