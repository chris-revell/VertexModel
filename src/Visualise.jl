#
#  Visualise.jl
#  VertexModel
#
#  Created by Christopher Revell on 16/02/2021.
#
#
#

module Visualise

# Julia packages
using Printf
using LinearAlgebra
using ColorSchemes
using Colors
using UnPack
using GeometryBasics
using Random
using Makie
using CairoMakie
using StaticArrays
using SparseArrays
using CircularArrays
using FromFile

@from "OrderAroundCell.jl" using OrderAroundCell

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

function visualise(t, fig, ax1, ax2, mov, params, matrices)

    plotCells = 1
    plotEdges = 0
    scatterEdges = 0
    scatterVertices = 0
    scatterCells = 0
    plotForces = 0

    @unpack boundaryVertices, R, A, B, Bᵀ, C, cellPositions, edgeTangents, edgeMidpoints, F, ϵ = matrices
    @unpack nEdges, nVerts, nCells = params

    empty!(ax1)

    ax1.title = "t = $(@sprintf("%.2f", t))"
    xlims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    ylims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))

    # Plot cells
    if plotCells == 1
        for i = 1:nCells
            cellVertices = findall(x -> x != 0, C[i, :])
            vertexAngles = zeros(size(cellVertices))
            for (k, v) in enumerate(cellVertices)
                vertexAngles[k] = atan((R[v] .- cellPositions[i])...)
            end
            cellVertices .= cellVertices[sortperm(vertexAngles)]
            poly!(ax1, Point2f.(R[cellVertices]), color=(getRandomColor(i), 0.5))
        end
    end

    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax1, Point2f.(R), color=:green)
        annotations!(ax1, string.(collect(1:length(R))), Point2f.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax1, Point2f.(edgeMidpoints), color=:blue)
        annotations!(ax1, string.(collect(1:length(edgeMidpoints))), Point2f.(edgeMidpoints), color=:blue)
    end

    # Scatter cell positions
    if scatterCells == 1
        scatter!(ax1, Point2f.(cellPositions), color=:red)
        annotations!(ax1, string.(collect(1:length(cellPositions))), Point2f.(cellPositions), color=:red)
    end

    # Plot edges
    # For each edge, use A incidence matrix to find corresponding vertices x, and plot line between x[1] and x[2]
    if plotEdges == 1
        xs = Point2f[]
        us = Vec2f[]
        colours = Tuple{Symbol,Float64}[]
        for c = 1:nCells
            es = findall(x -> x != 0, B[c, :])
            for i in es
                vs = findall(x -> x != 0, A[i, :])
                colour = :black
                # Use B to set colour of edge depending on whether it runs with or against the orientation of the cell face
                if B[c, i] < 0
                    colour = (:red, 0.25)
                else
                    colour = (:blue, 0.25)
                end
                # Use A to set direction of arrow along edge
                if A[i, vs[1]] < 0
                    push!(xs, Point2f(R[vs[1]]))
                    push!(us, Vec2f(edgeTangents[i]))
                    push!(colours, colour)
                else
                    push!(xs, Point2f(R[vs[2]]))
                    push!(us, Vec2f(edgeTangents[i]))
                    push!(colours, colour)
                end
            end
        end
        arrows!(ax1, xs, us, color=colours, arrowsize=25, linewidth=5)
    end

    if plotForces == 1
        # Plot resultant forces on vertices (excluding external pressure)
        arrows!(ax1, Point2f.(R), Vec2f.(sum(F, dims=2)), color=:green)
        # Plot resultant forces on cells
        arrows!(ax1, Point2f.(cellPositions), Vec2f.(sum(F, dims=1)), color=:red)
    end

    # This routine will break if applied to a peripheral cell 
    empty!(ax2)
    centralCell = 1
    ax2.title = "Cell $centralCell force space"

    cellNeighbourMatrix = B * Bᵀ
    dropzeros!(cellNeighbourMatrix)
    neighbouringCells = findall(!iszero, cellNeighbourMatrix[centralCell, :]) # Note this includes centralCell itself 

    # Find and sort all vertices and edges around cells neighbouring centralCell, including centralCell itself
    cellVerticesDict = Dict()
    cellEdgesDict = Dict()
    for c in neighbouringCells
        cellVerticesDict[c], cellEdgesDict[c] = orderAroundCell(matrices,c)
    end

    # Sort neighbouringCells using shared edges 
    # This block will break if applied to a peripheral cell 
    orderedNeighbours = Int64[]
    for j in cellEdgesDict[centralCell]
        push!(orderedNeighbours, setdiff(findall(x->x!=0, B[:,j]), [centralCell])[1])
    end
    orderedNeighbours = CircularArray(orderedNeighbours)
    
    # Draw force network
    h = @SVector [0.0, 0.0]
    for (i, k) in enumerate(cellVerticesDict[centralCell])

        # For vertex k, the ith vertex in clockwise ordered list of vertices around centralCell, 
        # the corresponding edge and neighbour cell to loop around in force space are clockwise of k, 
        # and therefore found from index i+1: orderedNeighbours[i+1] and cellEdgesDict[centralCell][i+1]
        # Remember that these are circular arrays, so adding 1 to the last component takes you back to the start
        kNeighbourCell = orderedNeighbours[i+1]
        
        # Draw an arrow corresponding to the rotated force of centralCell on vertex k, starting from the last h space point
        arrows!(ax2, Point2f.([h]), Vec2f.([ϵ * F[k, centralCell]]), linewidth=4, arrowsize=16, color=(getRandomColor(centralCell), 0.75))
        # Update the last h space vector accordingly by adding the rotated force applied to vertex k by centralCell
        h = h + ϵ*F[k, centralCell]
        
        # List of static vectors to store h space vertices for corresponding neighbour cell 
        # hList = Array{SVector{2,Float64}}(undef, length(cellVerticesDict[kNeighbourCell]) + 1)
        hList = SVector{2,Float64}[]
        # List of rotated forces applied by neighbour cell to its vertices 
        rotatedForces = SVector{2,Float64}[]
        
        # Find the index of vertex k within ordered list of vertices around neighbouring cell kNeighbourCell
        index = findall(x->x==k, cellVerticesDict[kNeighbourCell])[1]
        
        # Store list of h space locations corresponding to rotated forces applied to vertices of cell kNeighbourCell by cell kNeighbourCell
        # Starting from last h space location, itself derived from the rotated force applied to vertex k by centralCell
        push!(hList,h)
        for m=1:length(cellVerticesDict[kNeighbourCell])
            push!(rotatedForces, ϵ*F[cellVerticesDict[kNeighbourCell][index+m-1], kNeighbourCell])
            push!(hList,hList[end]+rotatedForces[end])
        end
        
        # Plot rotated forces as arrows connecting points in h space
        arrows!(ax2, Point2f.(hList), Vec2f.(rotatedForces), color=(getRandomColor(kNeighbourCell), 0.75), linewidth=4, arrowsize=16)
    end

    # Need to set xlims and ylims in ax2 correctly

    recordframe!(mov)
    
    return nothing

end

export visualise

end
