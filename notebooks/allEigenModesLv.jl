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

fig = CairoMakie.Figure(resolution=(1000,1000))
ax = Axis(fig[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

mkpath(datadir(folderName,"eigenmodesLv","frame$(@sprintf("%03d", frame))"))

cellPolygons = makeCellPolygons(R,params,matrices)
linkTriangles = makeLinkTriangles(R,params,matrices)

decomposition = eigenmodesLv(R,matrices,params)

for mode=1:nVerts
    empty!(ax)
    lims = (-maximum(abs.(decomposition[:,mode])),maximum(abs.(decomposition[:,mode])))
    for k=1:nVerts
        poly!(ax,linkTriangles[k],color=[decomposition[k,mode]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
    end
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
    end
    save(datadir(folderName,"eigenmodesLv","frame$(@sprintf("%03d", frame))","mode$(@sprintf("%03d", mode)).png"),fig)
end
