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
@from "$(projectdir("src","OrderAroundCell.jl"))" using OrderAroundCell

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

function visualise(R, t, fig, ax1, mov, params, matrices)

    plotCells = 1
    scatterEdges = 1
    scatterVertices = 1
    scatterCells = 1

    @unpack boundaryVertices, A, B, Bᵀ, C, cellPositions, edgeTangents, edgeMidpoints, F, ϵ = matrices
    @unpack nEdges, nVerts, nCells = params

    empty!(ax1)

    ax1.title = "t = $(@sprintf("%.2f", t))"
    xlims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    ylims!(ax1, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))

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
    arrows!(ax1, Point2f.(R), Vec2f.(sum(F, dims=2)), color=:green)
    # Plot resultant forces on cells
    # arrows!(ax1, Point2f.(cellPositions), Vec2f.(sum(F, dims=1)), color=:red)

    recordframe!(mov)
    
    return nothing

end

export visualise

end
