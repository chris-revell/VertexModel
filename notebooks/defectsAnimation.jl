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

# folderName = "annealing/L₀=3.0_γ=0.15_23-02-02-17-13-02"
folderName = "annealing/L₀=0.75_γ=0.2_23-02-02-17-39-54"
# folderName = "annealing/L₀=0.5_γ=0.1_23-02-02-18-06-16"

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
    @unpack matrices = load(datadir(folderName,"frames","matrices$(@sprintf("%03d", t)).jld2"))
    @unpack params = load(datadir(folderName,"frames","params$(@sprintf("%03d", t)).jld2"))

    @unpack B, Bᵀ, C, R, cellPositions = matrices
    @unpack nCells = params

    cellNeighbourMatrix = B*Bᵀ

    neighbourCounts = [cellNeighbourMatrix[i,i] for i in 1:nCells]

    empty!(ax)

    ax.title = "t = $(@sprintf("%.2f", t))"
    xlims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))
    ylims!(ax, 1.1*min(minimum(first.(R)), minimum(last.(R))), 1.1*max(maximum(first.(R)), maximum(last.(R))))

    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(matrices, i)
        poly!(ax, Point2f.(R[orderedVertices]), strokecolor=(:black,0.5), strokewidth=1, colorrange=(-4, 4),colormap=:bwr, fillcolor=[neighbourCounts[i]-6])#, color=neighbourColours(neighbourCounts[i])
        poly!(ax, Point2f.(R[orderedVertices]), color=[neighbourCounts[i]-6] , colorrange=(-4, 4), colormap=colormap("RdBu", 5))
    end
    recordframe!(mov)
end

save(datadir(folderName,"defects.mp4"),mov)
