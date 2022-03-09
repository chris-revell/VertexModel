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
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeLengths,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

centralCell=14

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
cellVerticesDict = Dict()
for c in neighbouringCells
    # Find vertices around cell
    cellVertices = findall(x->x!=0,C[c,:])
    # Find angles of vertices around cell
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[c])...)
    end
    # Sort vertices around cell by polar angle
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    # Store sorted cell vertices for this cell
    cellVerticesDict[c] = cellVertices
end

# Draw all cell and vertex positions with #annotations
cellPolygons = Vector{Point2f}[]
for i=1:nCells
    cellVertices = findall(x->x!=0,C[i,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[i])...)
    end
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    push!(cellPolygons,Point2f.(R[cellVertices]))
end

edgeMidpointPolygons = Vector{Point2f}[]
for i=1:nCells
    cellEdges = findall(!iszero,B[i,:])
    edgeAngles = zeros(size(cellEdges))
    for (k,v) in enumerate(cellEdges)
        edgeAngles[k] = atan((edgeMidpoints[v].-cellPositions[i])...)
    end
    cellEdges .= cellEdges[sortperm(edgeAngles)]
    push!(edgeMidpointPolygons,Point2f.(edgeMidpoints[cellEdges]))
end

# poly!(ax1,cellPolygons[centralCell],color=(getRandomColor(centralCell),1.0))
for c in neighbouringCells
    poly!(ax1,cellPolygons[c],color=(getRandomColor(c),0.50),strokecolor=:black,strokewidth=8)
    poly!(ax1,edgeMidpointPolygons[c],color=(:white,0.0),strokecolor=(getRandomColor(c),1.0),strokewidth=8)

    # linked = findall(!iszero,cellNeighbourMatrix[c,:])
    # for l in linked
    #     if l!=c
    #         lines!(ax1,Point2f.(cellPositions[[l,c]]),color=:white,linewidth=4)
    #     end
    # end

    # Plot all vertex positions
    scatter!(ax1,Point2f.(R[cellVerticesDict[c]]),color=:blue,markersize=16)
    scatter!(ax1,Point2f.(edgeMidpoints[findall(!iszero,B[c,:])]),color=:green,markersize=16)
end
scatter!(ax1,Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red,markersize=16)

display(fig)

save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/cell$(centralCell)ForceNeighbourhood.pdf",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/cell$(centralCell)ForceNeighbourhood.svg",fig)
save("$dataDirectory/cell$(centralCell)ForceNeighbourhood.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/cell$(centralCell)ForceNeighbourhood.png",fig)
