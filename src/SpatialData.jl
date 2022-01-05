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
using LoopVectorization
using UnPack

function spatialData!(R,params,matrices)

    @unpack A,B,Ā,B̄,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints = matrices
    @unpack nCells,nEdges,γ,preferredPerimeter,preferredArea = params

    #cellPositions  .= C*R./cellEdgeCount    
    mul!(cellPositions,C,R)
    cellPositions .= cellPositions./cellEdgeCount
    #edgeTangents   .= A*R
    mul!(edgeTangents,A,R)
    edgeLengths    .= norm.(edgeTangents)
    #edgeMidpoints  .= 0.5.*Ā*R
    mul!(edgeMidpoints,0.5.*Ā,R)
    #cellPerimeters .= B̄*edgeLengths
    mul!(cellPerimeters,B̄,edgeLengths)
    cellTensions   .= γ.*(preferredPerimeter .- cellPerimeters)
    cellPressures  .= cellAreas .- preferredArea

    # Calculate oriented cell areas
    fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
    for i=1:nCells
        for j=1:nEdges
            cellOrientedAreas[i] += B[i,j]*edgeTangents[j]*edgeMidpoints[j]'
        end
        cellAreas[i] = cellOrientedAreas[i][1,2]
    end

    return nothing

end

export spatialData!

end
