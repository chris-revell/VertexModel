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

    @unpack A,B,Ā,B̄,Bᵀ,C,cellEdgeCount,cellVertexOrders,cellEdgeOrders,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,edgeMidpointLinks,vertexAreas,μ,Γ = matrices
    @unpack nCells,nEdges,nVerts,γ,L₀,A₀ = params

    # cellPositions  .= C*R./cellEdgeCount
    mul!(cellPositions,C,R)
    @.. thread=false cellPositions ./= cellEdgeCount

    # edgeTangents   .= A*R
    mul!(edgeTangents,A,R)

    @.. thread=false edgeLengths .= norm.(edgeTangents)

    # edgeMidpoints  .= 0.5.*Ā*R
    mul!(edgeMidpoints,Ā,R)
    @.. thread=false edgeMidpoints .*= 0.5

    fill!(edgeMidpointLinks, SVector{2,Float64}(zeros(2)))
    dropzeros!(edgeMidpointLinks)
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1],nzC[2])
    for (i,k) in ikPairs
        # for j in findnz(A[:,k])[1]
        for j in cellEdgeOrders[i]
            edgeMidpointLinks[i,k] = edgeMidpointLinks[i,k] .+ 0.5.*B[i,j].*edgeTangents[j].*Ā[j,k]
        end
    end
    
    # Find vertex areas, with special consideration of peripheral vertices with 1 or 2 adjacent cells
    for k=1:nVerts
        k_is = findall(x->x!=0, C[:,k])
        if length(k_is) == 1
            k_js = findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[k_js[1]]...,0.0]×[edgeTangents[k_js[2]]...,0.0])
        elseif length(k_is) == 2
            # edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            edgesSharedBy_i1_And_k = findall(x->x!=0,B[k_is[1],:].*A[:,k]) # This line does the same as the commented out line above but seems to allocate less memory
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
            # edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
            edgesSharedBy_i2_And_k = findall(x->x!=0,B[k_is[2],:].*A[:,k]) # This line does the same as the commented out line above but seems to allocate less memory
            vertexAreas[k] += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        else
            vertexAreas[k] = 0.5*norm([edgeMidpointLinks[k_is[1], k]...,0.0]×[edgeMidpointLinks[k_is[2],k]...,0.0])
        end
    end

    # cellPerimeters .= B̄*edgeLengths
    mul!(cellPerimeters,B̄,edgeLengths)

    # Calculate cell boundary tensions
    # @.. thread=false cellTensions   .= γ*L₀.*log.(L₀./cellPerimeters) #γ.*(L₀ .- cellPerimeters)
    cellTensions   .= Γ.*L₀.*log.(L₀./cellPerimeters) #γ.*(L₀ .- cellPerimeters)

    # Calculate oriented cell areas
    # fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
    # for i=1:nCells
    #     for j in nzrange(Bᵀ,i)
    #         cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
    #     end
    #     cellAreas[i] = cellOrientedAreas[i][1,2]
    # end

    for i=1:nCells
        cellAreas[i] = abs(area(Point{2,Float64}.(R[cellVertexOrders[i]])))
    end

    # Calculate cell internal pressures
    # @.. thread=false cellPressures  .= ones(nCells).*log.(cellAreas./A₀) #cellAreas .- A₀
    cellPressures  .= A₀.*μ.*log.(cellAreas./A₀) #cellAreas .- A₀

    return nothing

end

export spatialData!

end

#47, 83
