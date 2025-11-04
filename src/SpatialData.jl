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
@from "AnalysisFunctions.jl" using AnalysisFunctions

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
        boundaryCells,
        cellPerimeters,
        cellOrientedAreas,
        cellShapeTensor,
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

    cellPolygons = makeCellPolygons(R, params, matrices)
    # cellPositions= C*R ./ cellEdgeCount
    for i in 1:nCells
        # Check whether the cell is on the boundary 
        # if boundaryCells[i]==1
            L_x,L_y = 10,10
            # get vertices; 
            # verts is a vector of points
            verts = [SVector(v[1], v[2]) for v in cellPolygons[i]]  # make immutable copies

            for k in 2:length(verts)
                dx = verts[k][1] - verts[1][1]
                if dx > L_x/2
                    verts[k] = SVector(verts[k][1]-L_x, verts[k][2])
                elseif dx < -L_x/2
                    verts[k] = SVector(verts[k][1]+L_x, verts[k][2])
                end
                dy = verts[k][2] - verts[1][2]
                if dy > L_y/2
                    verts[k] = SVector(verts[k][1], verts[k][2]-L_y)
                elseif dy < -L_y/2
                    verts[k] = SVector(verts[k][1], verts[k][2]+L_y)
                end
            end

            # now average
            mean_x = mean(v[1] for v in verts)
            mean_y = mean(v[2] for v in verts)
            cellPositions[i] = SVector(mod(mean_x, L_x), mod(mean_y, L_y))


        # end
    end
    
    # Computing edge data: 
    for j in 1:nEdges
        L_x,L_y = 10,10
        # Get the vertices of edge j: 
        verts = findall(x -> x!=0, A[j,:])
        v1, v2 = R[verts[1]], R[verts[2]]

        # unwrap for periodic boundaries:
        dx = v2[1] - v1[1]
        if dx >  L_x/2
            v2 = SVector(v2[1]-L_x, v2[2])
        elseif dx < -L_x/2
            v2 = SVector(v2[1]+L_x, v2[2])
        end

        dy = v2[2] - v1[2]
        if dy >  L_y/2
            v2 = SVector(v2[1], v2[2]-L_y)
        elseif dy < -L_y/2
            v2 = SVector(v2[1], v2[2]+L_y)
        end

        # tangent and length
        tangent = v2 - v1
        length  = norm(tangent)
        # midpoint (wrap back into domain if needed)
        midpoint = SVector(mod(v1[1] + 0.5*tangent[1], L_x), mod(v1[2] + 0.5*tangent[2], L_y))

        edgeTangents[j] = tangent
        edgeLengths[j]  = length
        edgeMidpoints[j] = midpoint
    end
    
    # edgeTangents   .= A*R
    
    # edgeLengths .= norm.(edgeTangents)

    # edgeMidpoints  .= 0.5.*Ā*R
    
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

    # Find cell areas and shape tensors 
    # for i = 1:nCells
    #     cellAreas[i] = abs(area(Point{2,Float64}.(R[cellVertexOrders[i]])))

    #     Rα = [R[kk].-matrices.cellPositions[i] for kk in cellVertexOrders[i]]
    #     cellShapeTensor[i] = sum(Rα.*transpose.(Rα))./cellEdgeCount[i]
    # end

    for i = 1:nCells
        L_x,L_y = 10,10
        verts = R[cellVertexOrders[i]]
        # unwrap all vertices relative to the first vertex
        x0, y0 = verts[1]
        unwrapped = [SVector(x0, y0)]
        for v in verts[2:end]
            dx = v[1] - x0
            dy = v[2] - y0
            if dx >  L_x/2; dx -= L_x
            elseif dx < -L_x/2; dx += L_x
            end
            if dy >  L_y/2; dy -= L_y
            elseif dy < -L_y/2; dy += L_y
            end
            push!(unwrapped, SVector(x0 + dx, y0 + dy))
        end
    
        # compute cell center, area, shape tensor
        cellPositions[i] = SVector(mod(mean(v[1] for v in unwrapped), L_x),
                                   mod(mean(v[2] for v in unwrapped), L_y))
        cellAreas[i] = abs(area(Point{2,Float64}.(unwrapped)))
        
        Rα = [v - cellPositions[i] for v in unwrapped]
        cellShapeTensor[i] = sum(Rα .* transpose.(Rα)) / length(Rα)
    end

    # Calculate cell pressures and tensions according to energy model choice 
    if energyModel == "log"
        # Model per Cowley et al. 2024 Section 2a
        # Calculate cell boundary tensions
        cellTensions .= μ .* Γ .* cellL₀s .* log.(cellPerimeters ./ cellL₀s)
        # Calculate cell internal pressures
        cellPressures .= μ .* cellA₀s .* log.(cellAreas ./ cellA₀s)
    else
        # Quadratic energy model
        # Calculate cell boundary tensions
        cellTensions .= μ .* Γ .*(cellPerimeters - cellL₀s)
        # Calculate cell internal pressures
        cellPressures .= μ .*(cellAreas - cellA₀s)
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
