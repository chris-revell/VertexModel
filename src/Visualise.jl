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
using DrWatson

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "AnalysisFunctions.jl" using AnalysisFunctions

function visualise(R, t, fig, ax1, mov, params, matrices, plotCells,scatterEdges,scatterVertices,scatterCells,plotForces,plotEdgeMidpointLinks)

    @unpack boundaryVertices, A, Ā, B, B̄, Bᵀ, C, cellPressures, cellTensions, cellPositions, edgeTangents, edgeLengths, edgeMidpoints, F, ϵ, edgeMidpointLinks, μ = matrices
    @unpack nEdges, nVerts, nCells = params

    empty!(ax1)

    ax1.title = "t = $(@sprintf("%.2f", t))"

    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R,params,matrices)
        for i=1:nCells
            # poly!(ax1,cellPolygons[i],color=(getRandomColor(i), 0.5),strokecolor=(:black,1.0),strokewidth=2)
            if matrices.μ[i] < 2.0
                poly!(ax1,cellPolygons[i],color=:green,strokecolor=(:black,1.0),strokewidth=2)
            else
                poly!(ax1,cellPolygons[i],color=:purple,strokecolor=(:black,1.0),strokewidth=2)
            end
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

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    if plotForces == 1
        arrows!(ax1, Point2f.(R), Vec2f.(sum(F, dims=2)), color=:green)
    end

    if plotEdgeMidpointLinks == 1
        for i=1:nCells
            orderedVertices, orderedEdges = orderAroundCell(matrices, i)
            for kk=1:length(orderedVertices)
                lines!(ax1, [Point2(edgeMidpoints[orderedEdges[kk]]...), Point2( (edgeMidpoints[orderedEdges[kk]].+edgeMidpointLinks[i,orderedVertices[kk]])... )], linestyle=:dot, color=:black)
            end
        end
    end

    # Set limits
    reset_limits!(ax1)

    recordframe!(mov)
    
    return nothing

end

export visualise

end
