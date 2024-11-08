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
using Statistics

@from "OrderAroundCell.jl" using OrderAroundCell
# @from "AnalysisFunctions.jl" using AnalysisFunctions

function spatialData!(RH,params,matrices)

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
        cellShapeTensor,
        cellAreas,
        cellTensions,
        cellPressures,
        cellHeights,
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
        Γa,
        ΓA,
        ΓL,
        L₀,
        A₀ = params

    R=RH[1:nVerts]
    
    
    cellPositions  .= C*R./cellEdgeCount
    
    edgeTangents   .= A*R
    
    @.. thread=false edgeLengths .= norm.(edgeTangents)

    edgeMidpoints  .= 0.5.*Ā*R
    
    fill!(edgeMidpointLinks, SVector{2,Float64}(zeros(2)))
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
            vertexAreas[k] = 0.5^3 * norm([edgeTangents[k_js[1]]..., 0.0] × [edgeTangents[k_js[2]]..., 0.0])
        elseif length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x -> x != 0, B[k_is[1], :] .* A[:, k])
            vertexAreas[k] = 0.5^3 * norm([edgeTangents[edgesSharedBy_i1_And_k[1]]..., 0.0] × [edgeTangents[edgesSharedBy_i1_And_k[2]]..., 0.0])
            edgesSharedBy_i2_And_k = findall(x -> x != 0, B[k_is[2], :] .* A[:, k])
            vertexAreas[k] += 0.5^3 * norm([edgeTangents[edgesSharedBy_i2_And_k[1]]..., 0.0] × [edgeTangents[edgesSharedBy_i2_And_k[2]]..., 0.0])
        else
            vertexAreas[k] = 0.5 * norm([edgeMidpointLinks[k_is[1], k]..., 0.0] × [edgeMidpointLinks[k_is[2], k]..., 0.0])
        end
    end

    cellPerimeters .= B̄ * edgeLengths

    for i = 1:nCells
        cellAreas[i] = abs(area(Point{2,Float64}.(R[cellVertexOrders[i]])))

        #cell_verts=findall(x-> x!=0, C[i,:])
        Rα=[R[cellVertexOrders[i]][y]-matrices.cellPositions[i] for y in 1:length(cellVertexOrders[i])]
        cellShapeTensor[i]=sum(Rα.*transpose.(Rα))/Float64.(cellEdgeCount[i])
    end
    #@show(RH[nVerts+1])
    if length(RH)>nVerts
        cellHeights .= RH[nVerts+1][1]
    else
        cellHeights .= 1.0
    end
    #@show H
    # Calculate cell boundary tensions
    @.. thread = false cellTensions .= (Γa.*(2.0 .*cellAreas .+ cellPerimeters.*cellHeights.-1.0) .+ (0.5*ΓA)).*cellHeights .+ (2*ΓL).*(cellPerimeters-L₀)

    # Calculate cell internal pressures
    @.. thread = false cellPressures .= cellHeights.*(cellAreas.*cellHeights-1) .+ (2*Γa).*(2.0 .*cellAreas .+ cellPerimeters.*cellHeights-1)

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
