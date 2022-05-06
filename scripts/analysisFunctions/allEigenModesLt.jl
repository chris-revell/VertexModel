# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using Random
using Colors
using JLD2
using Printf

# Local modules
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

function allEigenModesLt(dataDirectory)
    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
    matricesDict = load("$dataDirectory/matricesFinal.jld2")
    @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

    T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

    edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
    trapeziumAreas = abs.(area.(edgeTrapezia))

    linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
    linkTriangleAreas = abs.(area.(linkTriangles))

    cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)

    decomposition = (eigen(Matrix(Lₜ))).vectors

    isdir("$dataDirectory/eigenmodesLt") ? nothing : mkpath("$dataDirectory/eigenmodesLt")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/eigenmodesLt") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/eigenmodesLt")

    # Set up figure canvas
    fig = Figure(resolution=(1000,1000))
    grid = fig[1,1] = GridLayout()
    ax = Axis(grid[1,1],aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    for eigenvectorIndex=2:nVerts
        empty!(ax)
        lims = (-maximum(abs.(decomposition[:,eigenvectorIndex])),maximum(abs.(decomposition[:,eigenvectorIndex])))
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        save("$dataDirectory/eigenmodesLt/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
        save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/eigenmodesLt/eigenmode$(@sprintf("%03d", eigenvectorIndex)).png",fig)
    end
end
