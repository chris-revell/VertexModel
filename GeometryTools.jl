#
#  GeometryTools.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 09/02/2021.
#
#
#

module GeometryTools

# Julia packages
using LinearAlgebra

# Local modules


@inline @views function geometryTools!(A,Ā,B,B̄,C,R,nCells,nEdges,nVerts,cellPositions,cellEdgeCount,cellAreas,cellOrientedAreas,cellPerimeters,cellTensions,cellPressures,cellEffectivePressures,edgeLengths,edgeMidpoints,edgeTangents,gamma,preferredPerimeter)

    cellPositions .= C*R./cellEdgeCount
    edgeTangents  .= A*R
    edgeLengths   .= sqrt.(sum(edgeTangents.^2,dims=2))
    edgeMidpoints .= 0.5.*Ā*R
    cellPerimeters.= B̄*edgeLengths

    fill!(cellOrientedAreas,0.0)
    for i=1:nCells
        for j=1:nEdges
            cellOrientedAreas[i,:,:] .+= B[i,j].*edgeTangents[j,:]*edgeMidpoints[j,:]'
        end
    end
    cellAreas .= cellOrientedAreas[:,1,2]

    # Effective pressure of each cell
    #cellEffectivePressures .= cellAreas .- 1.0 .+ gamma.*cellPerimeters.*(cellPerimeters.-preferredPerimeter)./(2.0.*cellAreas)    
    cellTensions           .= gamma.*(preferredPerimeter .- cellPerimeters)
    cellPressures          .= cellAreas .- 1.0

    return nothing

end

export geometryTools!

end
