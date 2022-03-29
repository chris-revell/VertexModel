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

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ = matricesDict["matrices"]


# Set up figure canvas
allForceFig = Figure(resolution=(1000,1000))
allForceGrid = allForceFig[1,1] = GridLayout()
allForceAx = Axis(allForceGrid[1,1],aspect=DataAspect())
allForceAx.title="All resultant forces in system"
hidedecorations!(allForceAx)
hidespines!(allForceAx)

# Plot cell polygons
for i=1:nCells
  cellVertices = findall(x->x!=0,C[i,:])
  vertexAngles = zeros(size(cellVertices))
  for (k,v) in enumerate(cellVertices)
     vertexAngles[k] = atan((R[v].-cellPositions[i])...)
  end
  cellVertices .= cellVertices[sortperm(vertexAngles)]
  poly!(allForceAx,Point2f.(R[cellVertices]),color=(getRandomColor(i),0.5))
end

# # Scatter vertex locations
scatter!(allForceAx,Point2f.(R),alpha=0.5,color=:blue)
annotations!(allForceAx,string.(collect(1:nVerts)),Point2f.(R),color=:blue)

# Edge labels
annotations!(allForceAx,string.(collect(1:nEdges)),Point2f.(edgeMidpoints),color=:green)

# Scatter cell centroid locations
scatter!(allForceAx,Point2f.(cellPositions),color=:red)
annotations!(allForceAx,string.(collect(14:14)),Point2f.(cellPositions[14:14]),color=:red)

# Plot resultant forces on vertices (excluding external pressure)
arrows!(allForceAx,Point2f.(R),Vec2f.(sum(F,dims=2)),color=:blue)

# Plot resultant forces on cells
arrows!(allForceAx,Point2f.(cellPositions),Vec2f.(sum(F,dims=1)),color=:red)

display(allForceFig)
save("$dataDirectory/allForces.png",allForceFig)
