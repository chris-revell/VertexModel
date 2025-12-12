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

function unwrapVertex(vert,v0,L_x,L_y)
    # Function to unwrap vert relative to vertex v0 in periodic box of size L_x, L_y
    dx = vert[1] - v0[1]
    dy = vert[2] - v0[2]

    if dx > L_x/2
        dx -= L_x
    elseif dx < -L_x/2
        dx += L_x
    end

    if dy > L_y/2
        dy -= L_y
    elseif dy < -L_y/2
        dy += L_y
    end
    unwrapped = SVector(v0[1] + dx, v0[2] + dy)

    return unwrapped
end

function unwrapCellVertices(R,cellVertsIndex,L_x,L_y)
    # Function to unwrap all vertices of a cell relative to the first vertex
    verts = R[cellVertsIndex]
    v0 = verts[1]
    unwrappedVerts = SVector{2,Float64}[]

    push!(unwrappedVerts, v0)
    for vert in verts[2:end]
        unwrappedVert = unwrapVertex(vert,v0,L_x,L_y)
        push!(unwrappedVerts, unwrappedVert)
    end

    return unwrappedVerts
end

function unwrapEdge(v1,v2,L_x,L_y)
    # Function to unwrap edge defined by vertices v1 and v2 in periodic box of size L_x, L_y

    dx = v2[1] - v1[1]
    dy = v2[2] - v1[2]

    if dx > L_x/2
        dx -= L_x
    elseif dx < -L_x/2
        dx += L_x
    end
    if dy > L_y/2
        dy -= L_y
    elseif dy < -L_y/2
        dy += L_y
    end

    # Unwrap v2 relative to v1
    unwrappedEdge = SVector(v1[1] + dx, v1[2] + dy)
    
    return unwrappedEdge
end

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
        cellPositions,
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
        Γ,
        Λs = matrices
    @unpack initialSystem,
        nCells,
        nEdges,
        nVerts,
        γ,
        L₀,
        A₀,
        energyModel,
        L_x,
        L_y = params

    # cellPolygons = makeCellPolygons(R, params, matrices)
    
    if initialSystem == "new"

        cellPositions = C*R ./ cellEdgeCount
        edgeTangents   .= A*R
        edgeLengths .= norm.(edgeTangents)
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

        # Find cell areas and shape tensors 
        for i = 1:nCells

            cellAreas[i] = abs(area(Point{2,Float64}.(R[cellVertexOrders[i]])))

            Rα = [R[kk].-matrices.cellPositions[i] for kk in cellVertexOrders[i]]
            cellShapeTensor[i] = sum(Rα.*transpose.(Rα))./cellEdgeCount[i]
        end

    elseif initialSystem == "periodic"
        
        # Computing edge data: 

        for j in 1:nEdges
            
            # Get the vertices of edge j: 
            
            a = findall(x -> x > 0, @view A[j, :])[1]
            b = findall(x -> x < 0, @view A[j, :])[1]
            
            v1 = R[a]
            v2 = R[b]
            # unwrap v2 relative to v1
            v2unwr  = unwrapVertex(v2,v1,L_x,L_y)
    
            # tangent and length
            tangent = v1 - v2unwr
            
            edgeTangents[j] = tangent
            edgeLengths[j]  = norm(tangent)
            edgeMidpoints[j] = SVector(mod(v1[1] + 0.5*tangent[1], L_x), mod(v1[2] + 0.5*tangent[2], L_y))
        end

        # Edge midpoint links:

        fill!(edgeMidpointLinks, SVector{2,Float64}(zeros(2)))
        dropzeros!(edgeMidpointLinks)
        nzC = findnz(C)
        ikPairs = tuple.(nzC[1], nzC[2])

        for (i, k) in ikPairs
            for j in cellEdgeOrders[i]
                edgeMidpointLinks[i, k] = edgeMidpointLinks[i, k] .+ 0.5 .* B[i, j] .* edgeTangents[j] .* Ā[j, k]
            end
        end

        # Cell perimeters

        cellPerimeters .= B̄ * edgeLengths

        # Compute cell positions, areas and shape tensors:

        for i = 1:nCells
            verts_unwrapped = unwrapCellVertices(R, cellVertexOrders[i], L_x, L_y)

            # Wrapped cell center
            cx = mod(mean(v[1] for v in verts_unwrapped), L_x)
            cy = mod(mean(v[2] for v in verts_unwrapped), L_y)
            cellPositions[i] = SVector(cx, cy)
            
            # area
            cellAreas[i] = abs(area(Point{2,Float64}.(verts_unwrapped)))
            
            # shape tensor
            Rα = [v - cellPositions[i] for v in verts_unwrapped]
            cellShapeTensor[i] = sum(Rα .* transpose.(Rα)) / length(Rα)
    
            if cellAreas[i] < 0
                println("Negative area for cell", i)
            end
        end

        if sum(cellAreas) > L_x*L_y + 1e-8 || sum(cellAreas) < L_x*L_y - 1e-8

            error("⚠️ Error: Total cell area exceeds box area in periodic system")
            println(sum(cellAreas) , " vs ", L_x*L_y)
        end

        # Calculate oriented cell areas
        # fill!(cellOrientedAreas,SMatrix{2,2}(zeros(2,2)))
        # for i=1:nCells
        #     for j in nzrange(Bᵀ,i)
        #         cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
        #     end
        #     cellAreas[i] = cellOrientedAreas[i][1,2]
        # end

        # Recalculate cells at the periodic boundary
        fill!(boundaryCells, 0)
        for i in 1:nCells
    
            # verts_unwrapped = unwrapCellVertices(R, cellVertexOrders[i], L_x, L_y)
            verts = R[cellVertexOrders[i]]

            for v in verts
                dx = abs(v[1] - cellPositions[i][1])
                dy = abs(v[2] - cellPositions[i][2])
                if dx > L_x/2 || dy > L_y/2
                    boundaryCells[i] = 1
                    break
                end
            end
        end

        # println("boundary Cells:", boundaryCells)

        # Compute vertex areas: 

        for k=1:nVerts
            # cells incident to vertex k: 
            k_is = findall(x->x !=0, @view C[:,k])

            if length(k_is) == 1
                error("Vertex with only one incident cell found in periodic system")
            elseif length(k_is) == 2

                #edges incident to vertex k: 
                k_js = findall(x->x !=0, @view A[:,k])

                # First adjacent cell
                edgesSharedBy_i1_And_k = findall(x -> x != 0, B[k_is[1], :] .* A[:, k])
                # Second adjacent cell
                edgesSharedBy_i2_And_k = findall(x -> x != 0, B[k_is[2], :] .* A[:, k])

                vertexAreas[k] = 0.5^3 *norm([edgeTangents[edgesSharedBy_i1_And_k[1]]..., 0.0] ×
                         [edgeTangents[edgesSharedBy_i1_And_k[2]]..., 0.0]) +
                    0.5^3 *
                    norm([edgeTangents[edgesSharedBy_i2_And_k[1]]..., 0.0] ×
                         [edgeTangents[edgesSharedBy_i2_And_k[2]]..., 0.0])

                         
            else
                # 3 or more adjacent cells → same midpoint-link formula
                vertexAreas[k] =
                    0.5 *norm([edgeMidpointLinks[k_is[1], k]..., 0.0] ×[edgeMidpointLinks[k_is[2], k]..., 0.0])

            end
            if !isfinite(vertexAreas[k]) || vertexAreas[k] <= 0
                println("⚠️ Warning: Vertex $k has suspicious area = ", vertexAreas[k])
            end



        end

    end
    # println("cell edge orders:",cellEdgeOrders[1])
    

    # Calculate cell pressures and tensions according to energy model choice 
    if energyModel == "log"
        # Model per Cowley et al. 2024 Section 2a
        # Calculate cell boundary tensions
        cellTensions .= μ .* Γ .* cellL₀s .* log.(cellPerimeters ./ cellL₀s)
        # Calculate cell internal pressures
        cellPressures .= μ .* cellA₀s .* log.(cellAreas ./ cellA₀s)
    elseif energyModel == "quadratic"
        # Quadratic energy model
        # Calculate cell boundary tensions
        cellTensions .= μ .* Γ .*(cellPerimeters - cellL₀s)
        # Calculate cell internal pressures
        cellPressures .= μ .*(cellAreas - cellA₀s)
    elseif energyModel == "quadratic2pops"
        # Quadratic energy model with two cell populations 
        cellTensions .= Γ .* cellPerimeters 
        cellPressures .= cellAreas .- 1.0
        
    end

    return nothing

end

export spatialData!

end



