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

function spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,γ,preferredPerimeter,preferredArea)

    cellPositions  .= C*R./cellEdgeCount
    edgeTangents   .= A*R
    edgeLengths    .= norm.(edgeTangents)
    edgeMidpoints  .= 0.5.*Ā*R
    cellPerimeters .= B̄*edgeLengths
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
