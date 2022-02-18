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

initialSystem = "data/sims/2022-02-18-11-40-49"

# Import system data
conditionsDict    = load("$initialSystem/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$initialSystem/matricesFinal.jld2")
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,cellPressures,cellTensions,cellPerimeters = matricesDict["matrices"]


pressuresEffective = cellPressures .- cellTensions.*cellPerimeters./(2.0*cellAreas)


# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)

ax1.title = "Effective pressures"

clims = (-maximum(abs.(-2.0.*pressuresEffective)),maximum(abs.(-2.0.*pressuresEffective)))

# Plot cell polygons
for i=1:nCells
    cellVertices = findall(x->x!=0,C[i,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[i])...)
    end
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    poly!(ax1,Point2f.(R[cellVertices]),color=[-2.0*pressuresEffective[i]],colormap=:cork,colorrange=clims,strokecolor=(:black,1.0),strokewidth=5) #:bwr
end

Colorbar(fig[1, 2],limits=clims,colormap=:cork,flipaxis=false) #:bwr

display(fig)

save("$(datadir())/plots/effectivePressures.png",fig)
