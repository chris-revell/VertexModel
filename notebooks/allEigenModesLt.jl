using JLD2
using SparseArrays
using LinearAlgebra
using DrWatson
using FromFile
using UnPack
using CairoMakie
using Printf
using Colors

@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/Potentials.jl" using Potentials
@from "$(projectdir())/src/Eigenmodes.jl" using Eigenmodes

frame = 100

folderName = "pressureExternal=0.5_stiffnessFactor=2.0_γ=0.2_24-06-18-17-39-57"

files = [datadir("sims", folderName, "frameData", f) for f in readdir(datadir("sims", folderName, "frameData")) if occursin(".jld2", f)]

@unpack R, matrices, params = load(files[end];
    typemap=Dict("VertexModel.../VertexModelContainers.jl.VertexModelContainers.ParametersContainer" => ParametersContainer,
        "VertexModel.../VertexModelContainers.jl.VertexModelContainers.MatricesContainer" => MatricesContainer))
@unpack B, Bᵀ, C, cellPositions = matrices
@unpack nCells, nVerts = params

fig = CairoMakie.Figure(size=(1000, 1000))
ax = Axis(fig[1, 1][1, 1], aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

mkpath(datadir("sims", folderName, "eigenmodesLt", "frame$(@sprintf("%03d", frame))"))

cellPolygons = makeCellPolygonsOld(R, params, matrices)
linkTriangles = makeLinkTriangles(R, params, matrices)

decomposition = eigenmodesLt(R, matrices, params)

for mode = 1:nVerts
    empty!(ax)
    lims = (-maximum(abs.(decomposition[:, mode])), maximum(abs.(decomposition[:, mode])))
    for k = 1:nVerts
        poly!(ax,
            linkTriangles[k],
            color=decomposition[k, mode],
            colorrange=lims,
            colormap=:bwr,
            strokewidth=1,
            strokecolor=(:black, 0.25),
        )
    end
    for i = 1:nCells
        poly!(ax,
            cellPolygons[i],
            color=(:white, 0.0),
            strokecolor=(:black, 1.0),
            strokewidth=1,
        )
    end
    save(datadir("sims", folderName, "eigenmodesLt", "frame$(@sprintf("%03d", frame))", "mode$(@sprintf("%03d", mode)).png"), fig)
end
