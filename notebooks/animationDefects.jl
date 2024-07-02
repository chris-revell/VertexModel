# Script to produce a movie of defect evolution over cell network for a given system

using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors
using InvertedIndices

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions

function neighbourColours(x)
    if x <= 6
        return (:red, (6 - x) / 10.0)
    else
        return (:blue, (x - 6) / 10.0)
    end
end

# folderName = "MCCComparison/pressureExternal=0.5_stiffnessFactor=10.0_γ=0.2_24-06-24-21-29-00"

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
mov = VideoStream(fig, framerate=5)
files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2",f)]
for t = 5:length(files)
    @show t
    @unpack R, matrices, params = load(files[t]; 
        typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer"=>ParametersContainer, 
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer"=>MatricesContainer))
    @unpack B, Bᵀ, C, cellPositions = matrices
    @unpack nCells, nVerts = params
    cellNeighbourMatrix = B * Bᵀ
    neighbourCounts = [nnz(cellNeighbourMatrix[i, Not(i)]) for i in 1:nCells]
    empty!(ax)
    cellPolygons = makeCellPolygonsOld(R, params, matrices)
    for i = 1:nCells
        if MCCs[t][i] == 0
            poly!(ax,
                cellPolygons[i],
                color=neighbourColours(neighbourCounts[i]),
                strokecolor=(:black,0.5),
                strokewidth=1,
            )
        end
    end
    for i = 1:nCells
        if MCCs[t][i] != 0
            poly!(ax,
                cellPolygons[i],
                color=neighbourColours(neighbourCounts[i]),
                strokecolor=(:black,1.0),
                strokewidth=3,
            )
        end
    end
    reset_limits!(ax)
    recordframe!(mov)
end

save(datadir("sims", folderName, "movieDefects.mp4"), mov)

