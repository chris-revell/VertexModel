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
includet("$(projectdir())/src/TopologyChange.jl"); using .TopologyChange
includet("$(projectdir())/src/Initialise.jl"); using .Initialise
includet("$(projectdir())/src/SpatialData.jl"); using .SpatialData
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
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeLengths,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

centralCell=14

# Set up figure canvas
fig1 = Figure(resolution=(1000,1000))
grid1 = fig1[1,1] = GridLayout()
ax1 = Axis(grid1[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
# ax1.title = "Cell $centralCell Neighbourhood"
# New axis for force space
ax2 = Axis(grid1[1,2],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)
# ax2.title = "Cell $centralCell Force Network"

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
startPosition = @SVector [0.0,0.0]
for (i,v) in enumerate(cellVerticesDict[centralCell])
    arrows!(ax2,Point2f.([startPosition]),Vec2f.([ϵ*F[v,centralCell]]),linewidth=4,arrowsize=16,color=(getRandomColor(centralCell),0.75))
    annotations!(ax2,string.([v]),Point2f.([startPosition.+ϵ*F[v,centralCell]./2.0]),color=(getRandomColor(centralCell),0.75))
    startPosition = startPosition + ϵ*F[v,centralCell]
    H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[neighbouringCells[i]])+1)
    cellForces = SVector{2, Float64}[]
    # Circular permutation of vertices to ensure vertex v is the first index
    # in the ordered cellVertices list around cell neighbouringCells[i]
    index = findall(x->x==v, cellVerticesDict[neighbouringCells[i]])
    cellVertices = circshift(cellVerticesDict[neighbouringCells[i]],1-index[1])
    H[1] = startPosition
    for (j,cv) in enumerate(cellVertices)
        push!(cellForces,+ϵ*F[cv,neighbouringCells[i]])
        H[j+1] = H[j]+cellForces[end]
    end
    annotations!(ax2,string.(cellVertices),(Point2f.(H[1:end-1])+Vec2f.(cellForces)./2.0),color=(getRandomColor(neighbouringCells[i]),0.75))
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),1.0),linewidth=4,arrowsize=16)
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
poly!(ax1,cellPolygons[centralCell],color=(getRandomColor(centralCell),1.0))
for c in neighbouringCells
    poly!(ax1,cellPolygons[c],color=(getRandomColor(c),1.00))
    cellVerticesDict
    # Plot all vertex positions
    scatter!(ax1,Point2f.(R[cellVerticesDict[c]]),color=:blue)
    annotations!(ax1,string.(cellVerticesDict[c]),Point2f.(R[cellVerticesDict[c]]),color=:blue)
    # Plot resultant forces on vertices (excluding external pressure)
    arrows!(ax1,Point2f.(R[cellVerticesDict[c]]),Vec2f.(sum(F[cellVerticesDict[c],:],dims=2)),color=:blue)
end
# Scatter plot cell positions with ##annotations
scatter!(ax1,Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red)
annotations!(ax1,string.(neighbouringCells),Point2f.(cellPositions[push!(neighbouringCells,centralCell)]),color=:red)


# arrows!(ax2,[Point2f(0.1.*[-1,1])],[Vec2f(0.01.*[-0.001,-1])],color=:red)
# #annotations!(ax2,["-π"],[Point2f(0.1.*[-1,1].+0.02.*[-0.001,-1])])
# arrows!(ax2,[Point2f(0.1.*[-1,1])],[Vec2f(0.01.*[-1,0])],color=:green)
# #annotations!(ax2,["-π/2"],[Point2f(0.1.*[-1,1].+0.03.*[-1,0])])
# arrows!(ax2,[Point2f(0.1.*[-1,1])],[Vec2f(0.01.*[0,1])],color=:blue)
# #annotations!(ax2,["0"],[Point2f(0.1.*[-1,1].+0.01.*[0,1])])
# arrows!(ax2,[Point2f(0.1.*[-1,1])],[Vec2f(0.01.*[1,0])],color=:black)
# #annotations!(ax2,["π/2"],[Point2f(0.1.*[-1,1].+0.01.*[1,0])])

display(fig1)

save("$dataDirectory/png/cell$(centralCell)ForceSpace2.png",fig1)
