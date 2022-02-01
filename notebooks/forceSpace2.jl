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
using DelimitedFiles
using Random
using Colors
using JLD2

# Local modules
includet("$(projectdir())/src/TopologyChange.jl"); using .TopologyChange
includet("$(projectdir())/src/CreateRunDirectory.jl"); using .CreateRunDirectory
includet("$(projectdir())/src/Visualise.jl"); using .Visualise
includet("$(projectdir())/src/Initialise.jl"); using .Initialise
includet("$(projectdir())/src/Iterate.jl"); using .Iterate
includet("$(projectdir())/src/SpatialData.jl"); using .SpatialData

function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

# initialSystem = "data/sims/2022-01-24-20-29-36"
# initialSystem = "/Users/christopher/Dropbox (The University of Manchester)/Other/VertexModelStuff/2022-01-28-13-30-34"
# initialSystem = "$(projectdir())/data/sims/2022-01-20-14-52-51"
# initialSystem = "/Users/christopher/Dropbox (The University of Manchester)/Other/VertexModelStuff/2022-01-31-18-25-23"
# initialSystem = "data/sims/2022-01-31-13-34-41"
# initialSystem = "data/sims/2022-02-01-10-20-47"
initialSystem = "data/sims/2022-02-01-10-53-17"

# Import system data
conditionsDict    = load("$initialSystem/params.jld2")
@unpack γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]

# Initialise system matrices
params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)
topologyChange!(matrices)
spatialData!(matrices.R,params,matrices)

@unpack R,A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,F,ϵ,boundaryVertices = matrices
@unpack nVerts,nCells,nEdges,pressureExternal = params

# Calculate forces on vertices and cells
matrixF = Array{SVector{2,Float64}}(undef,params.nVerts,params.nCells)
fill!(matrixF,@SVector zeros(2))
for k=1:nVerts
    for i=1:nCells
        for j=1:nEdges
            matrixF[k,i] -= 0.5*cellPressures[i]*B[i,j]*Ā[j,k]*(ϵ*edgeTangents[j]) - cellTensions[i]*B̄[i,j]*A[j,k]*edgeTangents[j]/edgeLengths[j]
        end
    end
end

# Plot for one set of 3 cells
cell=63

fig1 = Figure(resolution=(1000,1000))
grid1 = fig1[1,:] = GridLayout()
ax1 = Axis(grid1[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)

# Find all cells neighbouring original cell
cellNeighbourMatrix = matrices.B*matrices.Bᵀ
dropzeros!(cellNeighbourMatrix)
neighbouringCells = findnz(cellNeighbourMatrix[cell,:])[1]

# Draw all cell and vertex positions with annotations
cellVerticesDict = Dict()
for c in neighbouringCells
    # Find vertices around cell
    cellVertices = findall(x->x!=0,matrices.C[c,:])
    # Find angles of vertices around cell
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((matrices.R[v].-matrices.cellPositions[c])...)
    end
    # Sort vertices around cell by polar angle
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    # Store sorted cell vertices for this cell
    cellVerticesDict[c] = cellVertices
    # Draw cell as polygon using vertices
    poly!(ax1,Point2f.(matrices.R[cellVertices]),color=(getRandomColor(c),0.25))
    # Plot all vertex positions
    scatter!(ax1,Point2f.(matrices.R[cellVertices]),color=:black)
    annotations!(ax1,string.(cellVertices),Point2f.(matrices.R[cellVertices]),color=:blue)
    # Plot resultant forces on vertices (excluding external pressure)
    arrows!(ax1,Point2f.(R[cellVertices]),Vec2f.(sum(matrixF[cellVertices,:],dims=2)),color=:blue,alpha=0.3)
end
scatter!(ax1,Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)
annotations!(ax1,string.(neighbouringCells),Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)






# New axis for force space
ax2 = Axis(grid1[2,1],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)

# firstVertex = cellVerticesDict[cell][1]
firstVertex = 128
#threeNeighbours = findnz(matrices.C[:,firstVertex])[1]
threeNeighbours = [cell,59,15]
startPosition = @SVector [0.0,0.0]

colors = [:red,:green,:blue]

for (i,cell) in enumerate(threeNeighbours)

    # Force space vectors H
    H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[cell])+1)
    cellForces = SVector{2, Float64}[]

    # Ensure first vertex is the first index in the ordered cellVertices list
    firstVertexIndex = findall(x -> x==firstVertex, cellVerticesDict[cell])
    cellVertices = circshift(cellVerticesDict[cell],1-firstVertexIndex[1])

    H[1] = startPosition#s[i]

    for (i,v) in enumerate(cellVertices)
        push!(cellForces,-ϵ*matrixF[v,cell])
        H[i+1] = H[i]+cellForces[end]
    end

    annotations!(ax2,string.(cellVertices),(Point2f.(H[1:end-1])+Vec2f.(cellForces)./2.0))
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(colors[i],0.5),linewidth=5,arrowsize=25)

    startPosition = startPosition -ϵ*matrixF[cellVertices[1],cell]
end

# resultantForce = ϵ*sum(matrixF[128,:])
# arrows!(ax2,Point2f.([startPosition]),Vec2f.([resultantForce]),color=(:black,1.0),linewidth=5,arrowsize=25)

display(fig1)
