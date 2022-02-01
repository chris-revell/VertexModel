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

# Set up figure canvas
allForceFig = Figure(resolution=(1000,1000))
allForceGrid = allForceFig[1,1] = GridLayout()
allForceAx = Axis(allForceGrid[1,1],aspect=DataAspect(),)
hidedecorations!(allForceAx)
hidespines!(allForceAx)

# Plot cell polygons
for i=1:params.nCells
  cellVertices = findall(x->x!=0,matrices.C[i,:])
  vertexAngles = zeros(size(cellVertices))
  for (k,v) in enumerate(cellVertices)
     vertexAngles[k] = atan((matrices.R[v].-matrices.cellPositions[i])...)
  end
  cellVertices .= cellVertices[sortperm(vertexAngles)]
  poly!(allForceAx,Point2f.(matrices.R[cellVertices]),color=(getRandomColor(i),0.25))
end

# Scatter vertex locations
scatter!(allForceAx,Point2f.(matrices.R),alpha=0.5,color=:blue)
annotations!(allForceAx,string.(collect(1:params.nVerts)),Point2f.(matrices.R),color=:red)

# Scatter cell centroid locations
scatter!(allForceAx,Point2f.(matrices.cellPositions),color=:red)
annotations!(allForceAx,string.(collect(1:params.nCells)),Point2f.(matrices.cellPositions),color=:red)

# Plot resultant forces on vertices (excluding external pressure)
arrows!(allForceAx,Point2f.(R),Vec2f.(sum(matrixF,dims=2)),linecolor=:blue,alpha=0.3)

# Plot resultant forces on cells
arrows!(allForceAx,Point2f.(matrices.cellPositions),Vec2f.(sum(matrixF,dims=1)),linecolor=:red)

display(allForceFig)
