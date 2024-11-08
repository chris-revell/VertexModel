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

    params, matrices = p
    @unpack A,
        B,
        Ā,
        B̄,
        C,
        cellAreas,
        cellPerimeters,
        cellHeights,
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
        Γa,
        ΓA,
        ΓL,
        L₀ = params

    spatialData!(u, params, matrices)

    fill!(F, @SVector zeros(2))
    dropzeros!(F)
    fill!(externalF, @SVector zeros(2))


    peripheryLength = sum(boundaryEdges .* edgeLengths)

    #R,H=u
    # dHdA=((-2*Γa +ΓA).*cellAreas.* cellPerimeters .+(Γa.*cellPerimeters.^2).*(1 .- 2.0.*Γa.*cellPerimeters).+ 
    #   (cellAreas.^2).* (-1 .+ 2.0.*Γa.*cellPerimeters))./(((cellAreas.^2) .+Γa.*cellPerimeters.^2).^2)

    # dHdL=((2.0.*Γa -ΓA).*cellAreas.^2 .- 4.0.*Γa.*cellAreas.^3 .+Γa.*(-2.0.*Γa .+ΓA).*cellPerimeters.^2 .+ 4.0.*Γa .* cellAreas.*cellPerimeters.*(-1 .+Γa.*cellPerimeters))
    # ./(2.0.*((cellAreas.^2) .+Γa.*cellPerimeters.^2).^2)

    for k = 1:nVerts
        for j in nzrange(A, k)
            for i in nzrange(B, rowvals(A)[j])
                # Force components from cell pressure perpendicular to edge tangents 
                F[k, rowvals(B)[i]] += 0.5 * (cellPressures[rowvals(B)[i]] )* B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                F[k, rowvals(B)[i]] -= cellTensions[rowvals(B)[i]] * B̄[rowvals(B)[i], rowvals(A)[j]] * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
                # Force on vertex from external pressure 
                externalF[k] += boundaryVertices[k] * (0.5 * pressureExternal * B[rowvals(B)[i], rowvals(A)[j]] * Ā[rowvals(A)[j], k] .* (ϵ * edgeTangents[rowvals(A)[j]])) # 0 unless boundaryVertices != 0
                

            end
            # Force on vertex from peripheral tension
            externalF[k] -= boundaryEdges[rowvals(A)[j]] * peripheralTension * (peripheryLength - sqrt(π * nCells)) * A[rowvals(A)[j], k] .* edgeTangents[rowvals(A)[j]] ./ edgeLengths[rowvals(A)[j]]
        end

        du[k] = (sum(@view F[k, :]) .+ externalF[k]) 

    end

   
    
    
    #dAdr=-1/2*B*Diagonal([(ϵ*T)' for T in edgeTangents])*Ā
    #@show (dAdr'*cellHeights)./vertexAreas

    #du.+= (dAdr'*cellHeights)./vertexAreas
    
    #du[nVerts+1]=
    #@show u[nVerts+1]
    #test=(C'*(cellAreas.*(cellAreas.*u[nVerts+1][1].-1.0).+cellPerimeters.*(0.5.*ΓA .+ Γa.*(2.0.*cellAreas.+u[nVerts+1][1].*cellPerimeters.-1.0))))
    #@show test
    du[nVerts+1]=SVector(-sum(C'*(cellAreas.*(cellAreas.*u[nVerts+1][1].-1.0).+cellPerimeters.*(0.5.*ΓA .+ Γa.*(2.0.*cellAreas.+u[nVerts+1][1].*cellPerimeters.-1.0))))./nVerts, 0.0)
    #@show du[nVerts+1]
    return du
end



function modeltest!(du, u, p, t)

    params, matrices = p
    @unpack vertexAreas, g_vec = matrices
    #@unpack nVerts,nCells,nEdges,pressureExternal,peripheralTension = params

    spatialData!(u,params,matrices)

    M=makeM(matrices)

    du .= -(M'*g_vec)./(vertexAreas)
end

export model!, modeltest!

end
