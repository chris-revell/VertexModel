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
@from "RotationMatrix.jl" using RotationMatrix
# @from "AnalysisFunctions.jl" using AnalysisFunctions

function spatialData!(R,params,matrices)

    @unpack A,
        B,
        Ā,
        B̄,
        Bᵀ,
        C,
        cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        cellPerimeters,
        cellOrientedAreas,
        # cellShapeTensor,
        cellAreas,
        cellTensions,
        cellPressures,
        cellPerpAxes,
        edgeLengths,
        edgeTangents,
        edgeMidpoints,
        edgeMidpointLinks,
        vertexAreas,
        μ,
        Γ = matrices
    @unpack nCells,
        nEdges,
        nVerts,
        γ,
        L₀,
        A₀ = params

    cellPositions  .= C*R./cellEdgeCount
    
    edgeTangents   .= A*R
    
    @.. thread=false edgeLengths .= norm.(edgeTangents)

    edgeMidpoints  .= 0.5.*Ā*R
    
    fill!(edgeMidpointLinks, SVector{3,Float64}(zeros(3)))
    dropzeros!(edgeMidpointLinks)
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1], nzC[2])
    for (i, k) in ikPairs
        for j in cellEdgeOrders[i]
            edgeMidpointLinks[i, k] = edgeMidpointLinks[i, k] .+ 0.5 .* B[i, j] .* edgeTangents[j] .* Ā[j, k]
        end
    end

    # Find vertex areas, with special consideration of peripheral vertices with 1 or 2 adjacent cells
    for k = 1:nVerts
        k_is = findall(x -> x != 0, @view C[:, k])
        if length(k_is) == 1
            k_js = findall(x -> x != 0, A[:, k])
            vertexAreas[k] = 0.5^3 * norm(edgeTangents[k_js[1]] × edgeTangents[k_js[2]])
        elseif length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x -> x != 0, B[k_is[1], :] .* A[:, k])
            vertexAreas[k] = 0.5^3 * norm(edgeTangents[edgesSharedBy_i1_And_k[1]] × edgeTangents[edgesSharedBy_i1_And_k[2]])
            edgesSharedBy_i2_And_k = findall(x -> x != 0, B[k_is[2], :] .* A[:, k])
            vertexAreas[k] += 0.5^3 * norm(edgeTangents[edgesSharedBy_i2_And_k[1]] × edgeTangents[edgesSharedBy_i2_And_k[2]])
        else
            vertexAreas[k] = 0.5 * norm(edgeMidpointLinks[k_is[1], k] × edgeMidpointLinks[k_is[2], k])
        end
    end

    cellPerimeters .= B̄ * edgeLengths

    # Clockwise ordering of edges and vertices around cell face used in B 
    # means that the cross product of adjacent edge tangents defines a perpendicular 
    # vector into the cell face, with edge ordering then following the right hand rule 
    # around this perpendicular vector. 
    for i=1:nCells 
        cellPerpAxes[i] = (B[i,cellEdgeOrders[i][1]].*edgeTangents[cellEdgeOrders[i][1]])×(B[i,cellEdgeOrders[i][2]].*edgeTangents[cellEdgeOrders[i][2]]) # Don't need to normalize() this vector now because that is done later in the calculation of the rotation matrix
    end

    ϵFlatten = zeros(3,3)
    # Find cell areas
    for i = 1:nCells
        # cellAreas[i] = 0.0
        # for k=1:length(cellVertexOrders)
        #     cellAreas[i] += abs(area(Point{3,Float64}.( [cellPositions[i], R[cellVertexOrders[i][k]], R[cellVertexOrders[i][k+1]] ] )))
        # end
        crossVec = matrices.cellPerpAxes[i]×[1,0,0]
        ϵFlatten .= ϵ(v=crossVec, θ=asin(norm(crossVec)/(norm(matrices.cellPerpAxes[i]))))
        rotatedPoints = [Point{2,Float64}((ϵFlatten*pt)[2:end]) for pt in R[cellVertexOrders[i][0:end]]]
        cellAreas[i] = abs(area(rotatedPoints))
        # cellAreas[i] = abs(area(Point{3,Float64}.(R[cellVertexOrders[i]])))
    end

    
    # Calculate cell boundary tensions
    # @.. thread = false cellTensions .= Γ .* L₀ .* log.(cellPerimeters ./ L₀)
    @.. thread = false cellTensions .= Γ .* L₀ .* log.(L₀./cellPerimeters)

    # Calculate cell internal pressures
    # @.. thread = false cellPressures .= A₀ .* μ .* log.(cellAreas ./ A₀)
    @.. thread = false cellPressures .= A₀ .* μ .* log.(A₀./cellAreas)

    return nothing

end

export spatialData!

end


# Calculate oriented cell areas
# fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
# for i=1:nCells
#     for j in nzrange(Bᵀ,i)
#         cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
#     end
#     cellAreas[i] = cellOrientedAreas[i][1,2]
# end
