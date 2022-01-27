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

# Import local modules
includet("$(projectdir())/src/TopologyChange.jl"); using .TopologyChange
includet("$(projectdir())/src/CreateRunDirectory.jl"); using .CreateRunDirectory
includet("$(projectdir())/src/Visualise.jl"); using .Visualise
includet("$(projectdir())/src/Initialise.jl"); using .Initialise
includet("$(projectdir())/src/Iterate.jl"); using .Iterate
includet("$(projectdir())/src/SpatialData.jl"); using .SpatialData
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers

# Function for colouring cells
function getRandomColor(seed)
    Random.seed!(seed)
    rand(RGB{})
end

# Path to simulation data
initialSystem = "data/sims/2022-01-26-20-38-54/"


# Import parameters from simulation data
conditionsDict    = load("$initialSystem/params.jld2")
@unpack γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]


# Run vertex model functions to obtain required data for force calculations
params,matrices = initialise(initialSystem,realTimetMax,γ,λ,preferredArea,pressureExternal,dt,viscousTimeScale,outputTotal,t1Threshold,realCycleTime)
topologyChange!(matrices)
spatialData!(matrices.R,params,matrices)
@unpack R,A,B,Ā,B̄,cellTensions,cellPressures,edgeLengths,edgeTangents,ϵ = matrices
@unpack nVerts,nCells,nEdges = params


# matrixF[v,c] gives force applied to vertex v by cell c and vice versa
matrixF = Array{SVector{2,Float64}}(undef,params.nVerts,params.nCells)
fill!(matrixF,@SVector zeros(2))
# Fill matrixF by calculating force
for v=1:nVerts
    for c=1:nCells
        for j=1:nEdges
            matrixF[v,c] -= 0.5*cellPressures[c]*B[c,j]*Ā[j,v]*(ϵ*edgeTangents[j]) - cellTensions[c]*B̄[c,j]*A[j,v]*edgeTangents[j]/edgeLengths[j]
        end
    end
end


# Work with cell 31
cell=62

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
    annotations!(ax1,string.(cellVertices),Point2f.(matrices.R[cellVertices]),color=:black)
end
scatter!(ax1,Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)
annotations!(ax1,string.(neighbouringCells),Point2f.(matrices.cellPositions[neighbouringCells]),color=:red)


# New axis for force space
ax2 = Axis(grid1[2,1],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)

firstVertex = cellVerticesDict[cell][1]
#threeNeighbours = findnz(matrices.C[:,firstVertex])[1]
threeNeighbours = [62,30,79]
startPosition = @SVector [0.0,0.0]
#     [0.0,0.0],
#     -ϵ*matrixF[cellVerticesDict[threeNeighbours[1]][1],cell],
#     ϵ*matrixF[cellVerticesDict[threeNeighbours[2]][1],cell]-ϵ*matrixF[cellVerticesDict[threeNeighbours[1]][1],cell]
# ])

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
    arrows!(ax2,Point2f.(H),Vec2f.(cellForces),color=(:red,0.5),linewidth=5,arrowsize=25)

    startPosition = startPosition -ϵ*matrixF[cellVertices[1],cell]
end

display(fig1)


# ax3 = Axis(grid1[2,2],aspect=DataAspect())
# hidedecorations!(ax3)
# hidespines!(ax3)
# cell = threeNeighbours[2]
# H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[cell])+1)
# cellForces = SVector{2, Float64}[]
# for (i,v) in enumerate(cellVerticesDict[cell])
#     push!(cellForces,-ϵ*matrixF[v,cell])
#     H[i+1] = H[i]+cellForces[end]
# end
# annotations!(ax3,string.(cellVerticesDict[cell]),(Point2f.(H[1:end-1]).+Vec2f.(cellForces)./2.0))
# ax3.title = "Cell $cell"
# arrows!(ax3,Point2f.(H),Vec2f.(cellForces),color=(:blue,0.5),linewidth=5,arrowsize=25)
#
#
#
#
# # for cell 61
# cell = threeNeighbours[3]
# ax4 = Axis(grid1[2,3],aspect=DataAspect(),)
# hidedecorations!(ax4)
# hidespines!(ax4)
# H = Array{SVector{2,Float64}}(undef,length(cellVerticesDict[cell])+1)
# cellForces = SVector{2, Float64}[]
# for (i,v) in enumerate(cellVerticesDict[cell])
#     push!(cellForces,-ϵ*matrixF[v,cell])
#     H[i+1] = H[i]+cellForces[end]
# end
# annotations!(ax4,string.(cellVerticesDict[cell]),(Point2f.(H[1:end-1]).+Vec2f.(cellForces)./2.0))
# ax4.title = "Cell $cell"
# arrows!(ax4,Point2f.(H),Vec2f.(cellForces),color=(:green,0.5),linewidth=5,arrowsize=25)
#
#
#
#
# display(fig1)
