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
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ = matricesDict["matrices"]

cellNeighbourMatrix = B*Bᵀ
dropzeros!(cellNeighbourMatrix)

# Find vector of cell-cell links
onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints
T = SVector{2,Float64}[]
for j=1:nEdges
    Tⱼ = @SVector zeros(2)
    for i=1:nCells
        Tⱼ = Tⱼ + B[i,j]*(cellPositions[i].-cᵖ[j])
    end
    push!(T,Tⱼ)
end

centralCell=14

# Set up figure canvas
allForceFig = Figure(resolution=(1000,1000))
allForceGrid = allForceFig[1,1] = GridLayout()
allForceAx = Axis(allForceGrid[1,1],aspect=DataAspect())
# allForceAx.title="All resultant forces in system"
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
    if cellNeighbourMatrix[centralCell,i] == 0
        poly!(allForceAx,Point2f.(R[cellVertices]),color=(getRandomColor(i),0.25),strokecolor=(:black,0.25),strokewidth=1)
    else
        poly!(allForceAx,Point2f.(R[cellVertices]),color=(getRandomColor(i),1.0),strokecolor=(:black,1.0),strokewidth=1)
    end
end

for j=1:nEdges
    edgeCells = findall(!iszero,B[:,j])
    if boundaryEdges[j] == 0
        lines!(allForceAx,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    else
        lines!(allForceAx,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    end
end



# # Scatter vertex locations
scatter!(allForceAx,Point2f.(R),alpha=0.5,color=:blue)
# annotations!(allForceAx,string.(collect(1:nVerts)),Point2f.(R),color=:blue)

# Edge labels
# annotations!(allForceAx,string.(collect(1:nEdges)),Point2f.(edgeMidpoints),color=:green)

# Scatter cell centroid locations
scatter!(allForceAx,Point2f.(cellPositions),color=:red)
# annotations!(allForceAx,string.(collect(1:nCells)),Point2f.(cellPositions),color=:red)

display(allForceFig)
save("$dataDirectory/fullSystem.png",allForceFig)
