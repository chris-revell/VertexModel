# Script to produce a movie of defect evolution over cell network for a given system

using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using DataFrames
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/OrderAroundCell.jl" using OrderAroundCell

folderName = "newlongTest/L₀=0.75_realTimetMax=86400.0_t1Threshold=0.01_γ=0.2_23-03-08-20-49-23"

function neighbourColours(x)
    if x == 6
        return (:white, 0.0)
    elseif x == 5
        return (:red, 1.0)
    elseif x == 7
        return (:blue, 1.0)
    else
        return (:grey, 1.0)
    end
end

fig = CairoMakie.Figure(resolution=(1000,1000))
ax = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)

for t=0:100
    @unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", t)).jld2"))
    @unpack B, Bᵀ, C, cellPositions = matrices
    @unpack nCells, nVerts = params

    cellNeighbourMatrix = B*Bᵀ

    neighbourCounts = [cellNeighbourMatrix[i,i] for i in 1:nCells]

    empty!(ax)

    ax.title = "t = $(@sprintf("%.2f", t))"
    xlims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    ylims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))

    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices, i)
        poly!(ax, Point2f.(R[orderedVertices]), color=[neighbourCounts[i]-6] , colorrange=(-4, 4), colormap=colormap("RdBu", 5), strokecolor=(:black,1.0), strokewidth=1)
    end
    recordframe!(mov)
end

save(datadir(folderName,"defects.mp4"),mov)
