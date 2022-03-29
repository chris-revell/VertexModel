# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using SparseArrays
using StaticArrays
using CairoMakie
using UnPack
using GeometryBasics
using DelimitedFiles
using Random
using Colors
using JLD2

# Local modules
#includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
includet("$(projectdir())/notebooks/functions.jl")

#dataDirectory = "data/sims/2022-02-28-19-30-22"

isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeLengths,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

centralCell=1 #14

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
ax1 = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)

# Find all cells neighbouring original cell
cellNeighbourMatrix = B*Bᵀ
dropzeros!(cellNeighbourMatrix)
neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])

# Find and sort all vertices around cells neighbouring centralCell
cellVerticesDict = makeCellVerticesDict(conditionsDict["params"],matricesDict["matrices"])

# Draw all cell and vertex positions with #annotations
cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

edgeMidpointPolygons = makeEdgeMidpointPolygons(conditionsDict["params"],matricesDict["matrices"])

for c in neighbouringCells
    poly!(ax1,cellPolygons[c],color=(getRandomColor(c),0.50),strokecolor=:black,strokewidth=8)
    poly!(ax1,edgeMidpointPolygons[c],color=(:white,0.0),strokecolor=(getRandomColor(c),1.0),strokewidth=8)
    # Plot all vertex positions
    scatter!(ax1,Point2f.(R[cellVerticesDict[c]]),color=:blue,markersize=16)
    scatter!(ax1,Point2f.(edgeMidpoints[findall(!iszero,B[c,:])]),color=:green,markersize=16)
end
scatter!(ax1,Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red,markersize=16)

display(fig)

save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cell$(centralCell)ForceNeighbourhood.pdf",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cell$(centralCell)ForceNeighbourhood.svg",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cell$(centralCell)ForceNeighbourhood.png",fig)
