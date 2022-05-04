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
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

dataDirectory = "data/old/2022-02-28-19-30-22"

centralCell=14


# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,L₀,A₀,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeLengths,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

# Set up figure canvas
fig = Figure(resolution=(1000,500))
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
    poly!(ax1,cellPolygons[c],color=(getRandomColor(c),0.50),strokecolor=(:black,1.0),strokewidth=5)
    poly!(ax1,edgeMidpointPolygons[c],color=(:white,0.0),strokecolor=(getRandomColor(c),1.0),strokewidth=5)
    # Plot all vertex positions
    scatter!(ax1,Point2f.(R[cellVerticesDict[c]]),color=:blue,markersize=15)
    scatter!(ax1,Point2f.(edgeMidpoints[findall(!iszero,B[c,:])]),color=:green,markersize=15)
end
scatter!(ax1,Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red,markersize=15)

Label(fig[1,1,Bottom()],L"(a)",textsize = 32)

ax2 = Axis(fig[1,2],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)

# Find all cells neighbouring original cell
cellNeighbourMatrix = B*Bᵀ
dropzeros!(cellNeighbourMatrix)
neighbouringCells = findall(!iszero,cellNeighbourMatrix[centralCell,:])

# Find and sort all vertices around cells neighbouring centralCell
cellVerticesDict = makeCellVerticesDict(conditionsDict["params"],matricesDict["matrices"])

centralCellVertices = findall(x->x!=0,C[centralCell,:])
centralVertexAngles = zeros(size(centralCellVertices))
for (k,v) in enumerate(centralCellVertices)
   centralVertexAngles[k] = atan((R[v].-cellPositions[centralCell])...)
end
m = minimum(centralVertexAngles)
centralVertexAngles .-= m
centralCellVertices .= centralCellVertices[sortperm(centralVertexAngles)]

# Sort cells neighbouring centralCell by angle
setdiff!(neighbouringCells,[centralCell]) # Remove centralCell from neighbours list
neighbourAngles = zeros(length(neighbouringCells))
for (i,c) in enumerate(neighbouringCells)
    neighbourAngles[i] = atan((cellPositions[c].-cellPositions[centralCell])...)
end
neighbourAngles .+= (2π-m)
neighbourAngles = neighbourAngles.%(2π)
neighbouringCells .= neighbouringCells[sortperm(neighbourAngles)]

# Draw force network
startPosition = [SVector{2,Float64}([0.0,0.0])] #Make this a single component array of SVectors to avoid scope issues with immutable objects
for (i,v) in enumerate(cellVerticesDict[centralCell])
    arrows!(ax2,Point2f.([startPosition[1]]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=5,arrowsize=25,color=(getRandomColor(centralCell),0.9))
    startPosition[1] = startPosition[1] + ϵ*F[v,centralCell]
    H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
    cellForces = SVector{2, Float64}[]
    # Circular permutation of vertices to ensure vertex v is the first index
    # in the ordered cellVertices list around cell neighbouringCells[i]
    index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
    cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
    H[1] = startPosition[1]
    for (j,cv) in enumerate(cellVertices)
        push!(cellForces,ϵ*F[cv,neighbouringCells[i]])
        H[j+1] = H[j]+cellForces[end]
    end
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.9),linewidth=5,arrowsize=25)
end

Label(fig[1,2,Bottom()],L"(b)",textsize = 32)

# boxes = [Box(g,color=(:white,0.0)) for g in [fig[1,1],fig[1,2]]]

colgap!(fig.layout,1,Relative(0.0))

resize_to_layout!(fig)


display(fig)

save("$dataDirectory/png/figure4.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure4.eps",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure4.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/figure4.png",fig)
