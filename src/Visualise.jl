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

function visualise(R, t, fig, ax, mov, params, matrices, plotCells, scatterEdges, scatterVertices, scatterCells, plotForces, plotEdgeMidpointLinks)

    @unpack cellEdgeCount,
        cellVertexOrders,
        cellEdgeOrders,
        cellPositions,
        edgeMidpoints,
        F,
        edgeMidpointLinks,
        μ = matrices
    @unpack nEdges, 
    nVerts, 
    nCells, 
    cellsTypeA, 
    cellsTypeB = params


    empty!(ax)

    ax.title = "t = $(@sprintf("%.3f", t))"

    # Plot cells
    if plotCells == 1
        cellPolygons = makeCellPolygons(R, params, matrices)
        for i = 1:nCells

            if i in cellsTypeA

                poly!(ax, cellPolygons[i], color=RGB(102/255,178/255,255/255), strokecolor=(:black, 1.0), strokewidth=2)

            else

                poly!(ax, cellPolygons[i], color=RGB(255/255,178/255,102/255), strokecolor=(:black, 1.0), strokewidth=2)

            end
            
        end
    end

    # Scatter vertices
    if scatterVertices == 1
        scatter!(ax, Point{2,Float64}.(R), color=:green)
        annotations!(ax, string.(collect(1:length(R))), Point{2,Float64}.(R), color=:green)
    end

    # Scatter edge midpoints
    if scatterEdges == 1
        scatter!(ax, Point{2,Float64}.(edgeMidpoints), color=:blue)
        annotations!(ax, string.(collect(1:length(edgeMidpoints))), Point{2,Float64}.(edgeMidpoints), color=:blue)
    end

    # Scatter cell positions
    if scatterCells == 1
        scatter!(ax, Point{2,Float64}.(cellPositions), color=:red)
        annotations!(ax, string.(collect(1:length(cellPositions))), Point{2,Float64}.(cellPositions), color=:red)
    end

    # Plot resultant forces on vertices (excluding external pressure)
    # NB these forces will be those calculated in the previous integration step and thus will not be exactly up to date for the current vertex positions
    if plotForces == 1
        arrows!(ax, Point{2,Float64}.(R), Vec2f.(sum(F, dims=2)), color=:green)
    end

    if plotEdgeMidpointLinks == 1
        for i = 1:nCells
            for j = 1:cellEdgeCount[i]
                lines!(ax,
                    Point{2,Float64}.([edgeMidpoints[cellEdgeOrders[i][j]],(edgeMidpoints[cellEdgeOrders[i][j]] .+ edgeMidpointLinks[i, cellVertexOrders[i][j]])]),
                    linestyle=:dot,
                    color=:black)
            end
        end
    end

    # Set limits
    reset_limits!(ax)

    # Add frame to movie 
    recordframe!(mov)

    return nothing

end

export visualise

end
