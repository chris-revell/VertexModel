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
@unpack A,B,Bᵀ,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ = matricesDict["matrices"]

cellNeighbourMatrix = B*Bᵀ
dropzeros!(cellNeighbourMatrix)

onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints

# Find cell midpoint links T
T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

# Create vector of polygons for each cell
cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

centralCell=1 #14

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
ax = Axis(fig[1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    if cellNeighbourMatrix[centralCell,i] == 0
        poly!(ax,cellPolygons[i],color=(getRandomColor(i),0.25),strokecolor=(:black,0.5),strokewidth=1)
    else
        poly!(ax,cellPolygons[i],color=(getRandomColor(i),1.0),strokecolor=(:black,1.0),strokewidth=2)
    end
end

for j=1:nEdges
    edgeCells = findall(!iszero,B[:,j])
    if boundaryEdges[j] == 0
        lines!(ax,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    else
        lines!(ax,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    end
end

# # Scatter vertex locations
scatter!(ax,Point2f.(R),alpha=0.5,color=:blue)
# annotations!(ax,string.(collect(1:nVerts)),Point2f.(R),color=:blue)

# Edge labels
# annotations!(ax,string.(collect(1:nEdges)),Point2f.(edgeMidpoints),color=:green)

# Scatter cell centroid locations
scatter!(ax,Point2f.(cellPositions),color=:red)
# annotations!(ax,string.(collect(1:nCells)),Point2f.(cellPositions),color=:red)

display(fig)
save("$dataDirectory/fullSystem.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/fullSystem.pdf",fig)
save("$dataDirectory/fullSystem.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/fullSystem.svg",fig)
save("$dataDirectory/fullSystem.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/fullSystem.png",fig)
