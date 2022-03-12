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
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
includet("$(projectdir())/notebooks/functions.jl")

dataDirectory = "data/sims/2022-02-28-19-30-22"

isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg")

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

Lᵥ = makeLv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)

decomposition = (eigen(Matrix(Lᵥ))).vectors

# Set up figure canvas
fig = Figure(resolution=(600,2000))
grid = fig[1,1] = GridLayout()
for x=1:4
    for y=1:5
        eigenvectorIndex = ((y-1)*4 + x)+1
        lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
        ax = Axis(grid[y,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid[y,x,Bottom()],
                L"k=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end
for x=1:4
    for y=1:5
        eigenvectorIndex = ((y-1)*4 + x)+(nVerts-20)
        lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
        ax = Axis(grid[y+5,x],aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        for k=1:nVerts
            poly!(ax,linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
        end
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid[y+5,x,Bottom()],
                L"i=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end

display(fig)
save("$dataDirectory/eigenvectorTableauLv.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/eigenvectorTableauLv.svg",fig)
save("$dataDirectory/eigenvectorTableauLv.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/eigenvectorTableauLv.pdf",fig)
save("$dataDirectory/eigenvectorTableauLv.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/eigenvectorTableauLv.png",fig)
