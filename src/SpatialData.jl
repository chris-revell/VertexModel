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
using FromFile
using LinearAlgebra
using StaticArrays
using UnPack
using FastBroadcast
using SparseArrays
using GeometryBasics

@from "OrderAroundCell.jl" using OrderAroundCell
# @from "AnalysisFunctions.jl" using AnalysisFunctions

function spatialData!(R,params,matrices)

    @unpack A,B,Ā,B̄,Bᵀ,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints = matrices
    @unpack nCells,γ,L₀,A₀ = params

    # cellPositions  .= C*R./cellEdgeCount
    mul!(cellPositions,C,R)
    @.. thread=false cellPositions ./= cellEdgeCount

    # edgeTangents   .= A*R
    mul!(edgeTangents,A,R)

    @.. thread=false edgeLengths .= norm.(edgeTangents)

    # edgeMidpoints  .= 0.5.*Ā*R
    mul!(edgeMidpoints,Ā,R)
    @.. thread=false edgeMidpoints .*= 0.5

    # cellPerimeters .= B̄*edgeLengths
    mul!(cellPerimeters,B̄,edgeLengths)

    # Calculate cell boundary tensions
    @.. thread=false cellTensions   .= γ*L₀.*log.(L₀./cellPerimeters) #γ.*(L₀ .- cellPerimeters)

    # Calculate oriented cell areas
    fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
    # for i=1:nCells
    #     for j in nzrange(Bᵀ,i)
    #         cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
    #     end
    #     cellAreas[i] = cellOrientedAreas[i][1,2]
    # end

    for i=1:nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices,i)
        cellAreas[i] = abs(area(Point2f.(R[orderedVertices])))
    end
    # cellAreas .= abs.(area.(makeCellPolygons(R,params,matrices)))

    # Calculate cell internal pressures
    @.. thread=false cellPressures  .= A₀.*log.(cellAreas./A₀) #cellAreas .- A₀

    return nothing

end

export spatialData!

end
