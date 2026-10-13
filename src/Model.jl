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
using Distributions

# Local modules
@from "SpatialData.jl" using SpatialData

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
        vertexAreas,
        Λs = matrices
    @unpack initialSystem, 
        boundaryType,
        nVerts,
        nCells,
        nEdges,
        pressureExternal,
        peripheralTension,
        vertexWeighting,
        energyModel = params

    # Reinterpret state vector as a vector of SVectors 
    R = reinterpret(SVector{2,Float64}, u)
    dR = reinterpret(SVector{2,Float64}, du)

    spatialData!(R, params, matrices)

    fill!(F, @SVector zeros(2)) # internal forces on vertices 
    dropzeros!(F)
    fill!(externalF, @SVector zeros(2))
    

    peripheryLength = sum(boundaryEdges .* edgeLengths)

    for k = 1:nVerts
        for j in nzrange(A, k) # iterate over the nonzero entries for vertex k 
            for i in nzrange(B, rowvals(A)[j]) # rowvals(A) gives the row indices of nonzero entries of A
                
                # Force components from cell pressure perpendicular to edge tangents - the area derivative wrt. vertex position of Energy from pressure
                F[k, rowvals(B)[i]] += 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure -- only applies to boundary vertices 
                if boundaryType == "free"
                    externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                end
                
                if energyModel == "quadratic2pops"
                    # We have the separate edge tension term in that case: 
                    # Factor of 1/2 because this is added for each cell that meets at j
                    F[k, rowvals(B)[i]] -=  0.5 * Λs[rowvals(A)[j]] * B̄[rowvals(B)[i], rowvals(A)[j]] *  A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                end 
                
            end
            # Force on vertex from peripheral tension -- only for boundary edges 
            if boundaryType == "free"
                externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            end
        end
        if boundaryType == "free"
            dR[k] = (sum(@view F[k, :]) .+ externalF[k])
        elseif boundaryType == "periodic"
            dR[k] = sum(@view F[k, :])
        end
        
        
    end
    
    vertexWeighting == 1 ? dR ./= vertexAreas : nothing 

    # dR accesses the same underlying data as du, so by altering dR we have already updated du appropriately
    return du
end

function g!(du, u, p, t)
    params, matrices = p
    @unpack noiseWeighting, β, nVerts = params
    @unpack vertexAreas = matrices
    dR = reinterpret(SVector{2,Float64}, du)
    # If weighting noise, scale with sqrt vertex area: 
    for k=1:nVerts
        noiseWeighting == 1 ? dR[k] = SVector(β * sqrt(vertexAreas[k]),β * sqrt(vertexAreas[k])) : dR[k] = SVector(β, β)
    end
end

function bilayerModel!(du, u, p, t)

    # In this case, the state vector u = [x₁,y₁,x₂,y₂,x₃,y₃,x₅,y₅,x₆,y₆,κ,λ]
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
        vertexAreas,
        Λs = matrices
    @unpack initialSystem, 
        boundaryType,
        nVerts,
        nCells,
        nEdges,
        pressureExternal,
        peripheralTension,
        vertexWeighting,
        energyModel = params

    
    OptVec = reinterpret(SVector{Float64}, u)
    dOptVec = reinterpret(SVector{Float64}, du)

    # Vertex vector: 
    R = SVector{2,Float64}[6]
    dR = SVector{2,Float64}[6]
    for k = 1:6
        if k==4 
            R[k] = [-OptVec[6],-OptVec[7]]
        else
            R[k] = [OptVec[2*k-1],OptVec[2*k]]
        end
    end
    for k=1:3
        R[k] = [OptVec[2*k-1],OptVec[2*k]]
    end
    R[4] = [-OptVec[5],-OptVec[6]]
    for k=5:6
        R[k] = [OptVec[2*(k-1)-1],OptVec[2*(k-1)]]
    end

    κ = OptVec[11]
    λ = OptVec[12]
    dκ = []
    dλ = []
    
    spatialDataBilayer!(R, κ, λ, params, matrices)

    fill!(F, @SVector zeros(2)) # internal forces on vertices 
    dropzeros!(F)
    fill!(externalF, @SVector zeros(2))

    peripheryLength = sum(boundaryEdges .* edgeLengths)

    for k = 1:nVerts
        for j in nzrange(A, k) # iterate over the nonzero entries for vertex k 
            for i in nzrange(B, rowvals(A)[j]) # rowvals(A) gives the row indices of nonzero entries of A
                
                # Force components from cell pressure perpendicular to edge tangents - the area derivative wrt. vertex position of Energy from pressure
                F[k, rowvals(B)[i]] += 0.5 * cellPressures[rowvals(B)[i]] * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure -- only applies to boundary vertices 
                if boundaryType == "free"
                    externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                end
                
                if energyModel == "quadratic2pops"
                    # We have the separate edge tension term in that case: 
                    # Factor of 1/2 because this is added for each cell that meets at j
                    F[k, rowvals(B)[i]] -=  0.5 * Λs[rowvals(A)[j]] * B̄[rowvals(B)[i], rowvals(A)[j]] *  A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                end 
                
            end
            # Force on vertex from peripheral tension -- only for boundary edges 
            if boundaryType == "free"
                externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
            end
        end
        dR[k] = (sum(@view F[k, :]) .+ externalF[k])
        
        
        
    end
    vertexWeighting == 1 ? dR ./= vertexAreas : nothing

    # Update values of dOptVec for vertices (entries 1-10)
    for k=1:6
        if k!=4
            push!(dOptVec, dR[k][1])
            push!(dOptVec, dR[k][2])
        end
    end


    dκ = dUdκ

    # Update values of dOptVec for κ and λ (entries 11 and 12)
    push!(dOptVec, dκ)
    push!(dOptVec, dλ)

end

export model!, g!, bilayerModel! 

end