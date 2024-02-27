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
using InvertedIndices

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

function neighbourColours(x)
    if x<=6
        return (:red, (6-x)/10.0)
    else
        return (:blue, (x-6)/10.0)
    end
end

# for f in [f for f in readdir(datadir("sims/examples")) if occursin("γ",f)]

folderName = "/Users/christopher/Postdoc/Code/VertexModel/data/sims/L₀=0.75_nCells=800_realTimetMax=17300.0_γ=0.2_24-02-27-08-48-35"

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
    neighbourCounts = [nnz(cellNeighbourMatrix[i,Not(i)]) for i in 1:nCells]

    empty!(ax)
    cellPolygons = makeCellPolygonsOld(R,params,matrices)
    for i = 1:nCells
        if matrices.μ[i]>1.5
            poly!(ax, cellPolygons[i], color=neighbourColours(neighbourCounts[i]) , strokecolor=(:black,1.0), strokewidth=2)
        else 
            poly!(ax, cellPolygons[i], color=neighbourColours(neighbourCounts[i]) , strokecolor=(:black,0.2), strokewidth=2)
        end
    end
    # annotations!(ax, string.(collect(1:nCells)), Point{2,Float64}.(cellPositions))
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir(folderName,"movieDefects.mp4"),mov)

# end