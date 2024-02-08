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

function neighbourColours(x)

    if x==2
        return (:red, 0.8)
    elseif x==4
        return (:red,0.6)
    elseif x==4
        return (:red,0.4)
    elseif x == 5
        return (:red, 0.2)
    elseif x == 6
        return (:white, 0.0)
    elseif x == 7
        return (:blue, 0.2)
    elseif x == 8
        return (:blue, 0.4)
    elseif x == 9
        return (:red, 0.6)
    elseif x == 10
        return (:red, 0.8)
    elseif x == 11
        return (:red, 1.00)
    else
        return (:grey, 0.75)
    end
end


for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

    folderName = "sims/examples/$f"

    fig = CairoMakie.Figure(size=(1000,1000))
    ax = Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    mov = VideoStream(fig, framerate=5)

    for t=0:100
        @unpack R, matrices, params = load(datadir(folderName,"frameData","systemData$(@sprintf("%03d", t)).jld2"))
        @unpack B, Bᵀ, C, cellPositions = matrices
        @unpack nCells, nVerts = params

        cellNeighbourMatrix = B*Bᵀ
        neighbourCounts = [cellNeighbourMatrix[i,i] for i in 1:nCells]

        empty!(ax)
        # ax.title = "t = $(@sprintf("%.2f", t))"
        for i = 1:nCells
            orderedVertices, orderedEdges = orderAroundCell(matrices, i)
            poly!(ax, Point{2,Float64}.(R[orderedVertices]), color=neighbourColours(neighbourCounts[i]) , strokecolor=(:black,1.0), strokewidth=1)
        end
        reset_limits!(ax)
        recordframe!(mov)
    end

    save(datadir(folderName,"movieDefects.mp4"),mov)

end