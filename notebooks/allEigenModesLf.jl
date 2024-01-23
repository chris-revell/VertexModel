# Script to produce a movie of Airy Stress evolution over cell network for a given system

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

folderName = "newlongTest/L₀=0.75_realTimetMax=86400.0_t1Threshold=0.01_γ=0.2_23-03-08-20-49-23"

@unpack R, matrices, params = load(datadir(folderName,"frames","systemData$(@sprintf("%03d", frame)).jld2"))
@unpack B, Bᵀ, C, cellPositions = matrices
@unpack nCells,nVerts = params

fig = CairoMakie.Figure(size=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

mkpath(datadir(folderName,"eigenmodesLf","frame$(@sprintf("%03d", frame))"))

cellPolygons = makeCellPolygons(R,params,matrices)
# linkTriangles = makeLinkTriangles(R,params,matrices)

decomposition = eigenmodesLf(R,matrices,params)

for mode=1:nCells
    empty!(ax)
    lims = (-maximum(abs.(decomposition[:,mode])),maximum(abs.(decomposition[:,mode])))
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=[decomposition[i,mode]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    save(datadir(folderName,"eigenmodesLf","frame$(@sprintf("%03d", frame))","mode$(@sprintf("%03d", mode)).png"),fig)
end
