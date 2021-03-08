#
#  SpatialData.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 09/02/2021.
#
#
# Function to calculate spatial data including tangents, lengths, midpoints, tensions etc from incidence and vertex position matrices.

module SpatialData

# Julia packages
using LinearAlgebra

@inline @views function spatialData!(A,Ā,B,B̄,C,R,nCells,nEdges,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)

    cellPositions  .= C*R./cellEdgeCount
    edgeTangents   .= A*R
    edgeLengths    .= sqrt.(sum(edgeTangents.^2,dims=2))
    edgeMidpoints  .= 0.5.*Ā*R
    cellPerimeters .= B̄*edgeLengths
    cellTensions   .= gamma.*(preferredPerimeter .- cellPerimeters)
    cellPressures  .= cellAreas .- 1.0

    # Calculate oriented cell areas
    fill!(cellOrientedAreas,0.0)
    for i=1:nCells
        for j=1:nEdges
            cellOrientedAreas[i,:,:] .+= B[i,j].*edgeTangents[j,:]*edgeMidpoints[j,:]'
        end
    end
    cellAreas .= cellOrientedAreas[:,1,2]

    return nothing

end

export spatialData!

end
