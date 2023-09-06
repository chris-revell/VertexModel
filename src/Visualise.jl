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

function visualise(R, t, fig, ax1, mov, params, matrices)

    plotCells       = 1
    scatterEdges    = 0
    scatterVertices = 0
    scatterCells    = 0
    plotForces      = 1

    @unpack boundaryVertices, A, Ā, B, B̄, Bᵀ, C, cellPressures, cellTensions, cellPositions, edgeTangents, edgeLengths, edgeMidpoints, F, ϵ = matrices
    @unpack nEdges, nVerts, nCells = params

    empty!(ax1)

    ax1.title = "t = $(@sprintf("%.2f", t))"
    # xlims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    # ylims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))

    # Plot cells
    if plotCells == 1
        for i = 1:nCells
            orderedVertices, orderedEdges = orderAroundCell(matrices, i)
            poly!(ax1, Point2f.(R[orderedVertices]), color=(getRandomColor(i), 0.5))
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

    recordframe!(mov)
    
    return nothing

end

export visualise

end
