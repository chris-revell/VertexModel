#
#  SpatialData.jl
#  VertexModel
#
#  Function to calculate spatial data including tangents, lengths, midpoints, tensions etc from incidence and vertex position matrices.

module SpatialData

# Julia packages
using FromFile
using LinearAlgebra
using StaticArrays
using UnPack
using FastBroadcast
using SparseArrays
using GeometryBasics
using StaticArrays

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
        cellϵs,
        cellAreas,
        cellL₀s,
        cellA₀s,
        cellTensions,
        cellPressures,
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
        A₀,
        energyModel = params

    cellPositions  .= C*R./cellEdgeCount

    # for i=1:nCells
    #     cellA₀s[i] =  A₀ #*(1.0 + 0.1*norm(cellPositions[i][1:2]))
    #     cellL₀s[i] =  L₀ #*(1.0 + 0.3*params.surfaceRadius*acos(((cellPositions[i]-params.surfaceCentre)[3]/norm(cellPositions[i]-params.surfaceCentre))))
    #     # Γ[i] = γ#/(1.0 + 0.001*params.surfaceRadius*acos(((cellPositions[i]-params.surfaceCentre)[3]/norm(cellPositions[i]-params.surfaceCentre))))
    #     # μ[i] = 1.0#/(1.0 + 0.001*params.surfaceRadius*acos(((cellPositions[i]-params.surfaceCentre)[3]/norm(cellPositions[i]-params.surfaceCentre))))
    # end
    
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
            vertexAreas[k] = 0.125 * norm(edgeTangents[k_js[1]] × edgeTangents[k_js[2]]) # 0.125 = 0.5^3
        elseif length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x -> x != 0, B[k_is[1], :] .* A[:, k])
            vertexAreas[k] = 0.125 * norm(edgeTangents[edgesSharedBy_i1_And_k[1]] × edgeTangents[edgesSharedBy_i1_And_k[2]]) # 0.125 = 0.5^3
            edgesSharedBy_i2_And_k = findall(x -> x != 0, B[k_is[2], :] .* A[:, k])
            vertexAreas[k] += 0.125 * norm(edgeTangents[edgesSharedBy_i2_And_k[1]] × edgeTangents[edgesSharedBy_i2_And_k[2]]) # 0.125 = 0.5^3
        else
            vertexAreas[k] = 0.5 * norm(edgeMidpointLinks[k_is[1], k] × edgeMidpointLinks[k_is[2], k])
        end
    end

    cellPerimeters .= B̄ * edgeLengths
    
    # Find cell areas
    # crossVec = zeros(3)
    for i = 1:nCells
        # crossVec .= matrices.cellPerpAxes[i]×[1,0,0]
        # cellϵs[i] = SMatrix{3,3,Float64}(ϵ(v=crossVec, θ=-asin(norm(crossVec)/(norm(matrices.cellPerpAxes[i])))))
        # rotatedPoints = [Point{2,Float64}((cellϵs[i]*(pt.-cellPositions[i]))[2:end]) for pt in R[cellVertexOrders[i]]]
        # rotatedPoints = [Point{2,Float64}((cellϵs[i]*pt)[2:end]) for pt in R[cellVertexOrders[i]]]
        # cellAreas[i] = abs(area(rotatedPoints))
        cellAreas[i] = abs(area(Point{3,Float64}.(R[cellVertexOrders[i]])))
        # cellϵs[i] = SMatrix{3,3,Float64}(ϵ(v=cellPerpAxes[i]))
    end

    perpAxis = zeros(3)
    for j = 1:nEdges
        perpAxis .= params.surfaceCentre.-edgeMidpoints[j]
        matrices.edgeϵs[j] = MMatrix{3,3,Float64}(ϵ(v=perpAxis))
        # ϵ!(matrices.edgeϵs[j], v=perpAxis)
    end

    # Calculate cell pressures and tensions according to energy model choice 
    if energyModel == "log"
        # Model per Cowley et al. 2024 Section 2a
        # Calculate cell boundary tensions
        @.. thread = false cellTensions .= μ .* Γ .* cellL₀s .* log.(cellPerimeters ./ cellL₀s)
        # Calculate cell internal pressures
        @.. thread = false cellPressures .= μ .* cellA₀s .* log.(cellAreas ./ cellA₀s)
    else
        # Quadratic energy model
        # Calculate cell boundary tensions
        @.. thread = false cellTensions .= μ .* Γ .*(cellPerimeters - cellL₀s)
        # Calculate cell internal pressures
        @.. thread = false cellPressures .= μ .*(cellAreas - cellA₀s)
    end

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

    # for i=1:nCells 
    #     cellPerpAxes[i] = cellPositions[i].-params.surfaceCentre
    # end

    # Clockwise ordering of edges and vertices around cell face used in B 
    # means that the cross product of adjacent edge tangents defines a perpendicular 
    # vector into the cell face, with edge ordering then following the right hand rule 
    # around this perpendicular vector. 
    # for i=1:nCells 
    #     j = cellEdgeOrders[i][1]
    #     jj = cellEdgeOrders[i][2]
    #     cellPerpAxes[i] = (B[i,j].*edgeTangents[j])×(B[i,jj].*edgeTangents[jj]) # Don't need to normalize() this vector now because that is done later in the calculation of the rotation matrix
    #     # cellPerpAxes[i]⋅cellPositions[i] < 0 ? error("Flipped cell") : nothing
    #     j = cellEdgeOrders[i][2]
    #     jj = cellEdgeOrders[i][3]
    #     cellPerpAxes[i] += (B[i,j].*edgeTangents[j])×(B[i,jj].*edgeTangents[jj]) # Don't need to normalize() this vector now because that is done later in the calculation of the rotation matrix
    #     # cellPerpAxes[i]⋅cellPositions[i] < 0 ? error("Flipped cell") : nothing

    #     # Note: doing this cross product with 2 pairs of edges ensures that the process still works even if the edges in one pair are parallel
    # end
