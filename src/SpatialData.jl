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

    @unpack A,B,Ā,B̄,Bᵀ,C,cellEdgeCount,cellPositions,cellPerimeters,cellOrientedAreas,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,edgeMidpointLinks,vertexAreas = matrices
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
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1],nzC[2])
    for (i,k) in ikPairs
        for j=1:nEdges
            edgeMidpointLinks[i,k] = edgeMidpointLinks[i,k] .+ 0.5.*B[i,j]*edgeTangents[j]*Ā[j,k]
        end
    end
    
    for k=1:nVerts
        k_is = findall(x->x!=0, C[:,k])
        if length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
            edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        else
            vertexAreas[k] = 0.5*norm([edgeMidpointLinks[k_is[1], k]...,0.0]×[edgeMidpointLinks[k_is[2],k]...,0.0])
        end
    end

    # cellPerimeters .= B̄*edgeLengths
    mul!(cellPerimeters,B̄,edgeLengths)

    # Calculate cell boundary tensions
    @.. thread=false cellTensions   .= γ.*(L₀ .- cellPerimeters)

    # Calculate oriented cell areas
    fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
    for i=1:nCells
        for j in nzrange(Bᵀ,i)
            cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
        end
        cellAreas[i] = cellOrientedAreas[i][1,2]
    end

    # Calculate cell internal pressures
    @.. thread=false cellPressures  .= cellAreas .- A₀

    return nothing

end

export spatialData!

end
