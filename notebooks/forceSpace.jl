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
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

initialSystem = "data/sims/2022-02-07-13-30-05"

# Import system data
conditionsDict    = load("$initialSystem/params.jld2")
@unpack γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]

importedArrays = load("$initialSystem/matricesFinal.jld2")
@unpack A,B,R,F = importedArrays["matrices"]

# Initialise system matrices
params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)
topologyChange!(matrices)
spatialData!(R,params,matrices)

@unpack ϵ,cellPositions = matrices
@unpack nVerts,nCells,nEdges,pressureExternal = params


# Set up figure canvas
fig1 = Figure(resolution=(1000,1000))
grid1 = fig1[1,1] = GridLayout()
ax1 = Axis(grid1[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)


# Plot for one set of 3 cells
centralCell=4

hidedecorations!(ax1)
hidespines!(ax1)

ax1.title = "Cell $centralCell Neighbourhood"

# Find all cells neighbouring original cell
cellNeighbourMatrix = matrices.B*matrices.Bᵀ
dropzeros!(cellNeighbourMatrix)
neighbouringCells = findnz(cellNeighbourMatrix[centralCell,:])[1]

# Draw all cell and vertex positions with #annotations
cellVerticesDict = Dict()
for c in neighbouringCells
    # Find vertices around cell
    cellVertices = findall(x->x!=0,matrices.C[c,:])
    # Find angles of vertices around cell
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-matrices.cellPositions[c])...)
    end
    # Sort vertices around cell by polar angle
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    # Store sorted cell vertices for this cell
    cellVerticesDict[c] = cellVertices
    # Draw cell as polygon using vertices
    poly!(ax1,Point2f.(R[cellVertices]),color=(getRandomColor(c),0.25))
    # Plot all vertex positions
    scatter!(ax1,Point2f.(R[cellVertices]),color=:black)
    annotations!(ax1,string.(cellVertices),Point2f.(R[cellVertices]),color=:blue)
    # Plot resultant forces on vertices (excluding external pressure)
    arrows!(ax1,Point2f.(R[cellVertices]),Vec2f.(sum(F[cellVertices,:],dims=2)),color=:blue,alpha=0.3)
end
# Scatter plot cell positions with #annotations
scatter!(ax1,Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)
annotations!(ax1,string.(neighbouringCells),Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)


# New axis for force space
ax2 = Axis(grid1[2,1],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)

ax2.title = "Cell $centralCell Force Network"

setdiff!(neighbouringCells,[centralCell])
neighbourAngles = zeros(length(neighbouringCells))
for (i,c) in enumerate(neighbouringCells)
    neighbourAngles[i] = atan((cellPositions[c].-cellPositions[centralCell])...)
end
neighbouringCells .= neighbouringCells[sortperm(neighbourAngles)]

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
        push!(cellForces,ϵ*F[cv,neighbouringCells[i]])
        H[j+1] = H[j]+cellForces[end]
    end

    annotations!(ax2,string.(cellVertices),(Point2f.(H[1:end-1])+Vec2f.(cellForces)./2.0),color=(getRandomColor(neighbouringCells[i]),0.75))
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(getRandomColor(neighbouringCells[i]),0.75),linewidth=4,arrowsize=16)

end

display(fig1)
