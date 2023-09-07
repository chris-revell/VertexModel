#
#  SpatialData.jl
#  VertexModel
#
#  Created by Christopher Revell on 09/02/2021.
#
#
# Function to calculate spatial data including tangents, lengths, midpoints, tensions etc from incidence and vertex position matrices.

module SpatialData

# Julia packages
using LinearAlgebra
using StaticArrays
using UnPack
using FastBroadcast
using SparseArrays

function spatialData!(R,params,matrices)

    @unpack A,B,Ā,B̄,Bᵀ,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints = matrices
    @unpack nCells,γ,L₀,A₀ = params

    # cellPositions  .= C*R./cellEdgeCount
    mul!(cellPositions,C,R)
    @.. thread=true cellPositions ./= cellEdgeCount

    # edgeTangents   .= A*R
    mul!(edgeTangents,A,R)

    @.. thread=true edgeLengths .= norm.(edgeTangents)

    # edgeMidpoints  .= 0.5.*Ā*R
    mul!(edgeMidpoints,Ā,R)
    @.. thread=true edgeMidpoints .*= 0.5

    # cellPerimeters .= B̄*edgeLengths
    mul!(cellPerimeters,B̄,edgeLengths)

    # Calculate cell boundary tensions
    @.. thread=true cellTensions   .= γ.*(L₀ .- cellPerimeters)

    # Calculate oriented cell areas
    fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
    for i=1:nCells
        for j in nzrange(Bᵀ,i)
            cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
        end
        cellAreas[i] = cellOrientedAreas[i][1,2]
    end

    # Calculate cell internal pressures
    @.. thread=true cellPressures  .= cellAreas .- A₀

    return nothing

end

export spatialData!

end
