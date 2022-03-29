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
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,cellPressures,cellTensions,cellPerimeters = matricesDict["matrices"]


# Create vector of polygons for each cell
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


pressuresEffective = cellPressures .- cellTensions.*cellPerimeters./(2.0*cellAreas)


# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
# ax1.title = "Effective pressures"
clims = (-maximum(abs.(-2.0.*pressuresEffective)),maximum(abs.(-2.0.*pressuresEffective)))
# Plot cell polygons
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[-2.0*pressuresEffective[i]],colormap=:bwr,colorrange=clims,strokecolor=(:black,1.0),strokewidth=5)
end
# Colorbar(fig[1, 2],limits=clims,colormap=:bwr)

display(fig)

save("$dataDirectory/effectivePressures.png",fig)
