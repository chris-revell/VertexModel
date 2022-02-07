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

initialSystem = "data/sims/2022-02-07-13-30-05"

# Import system data
conditionsDict    = load("$initialSystem/params.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$initialSystem/matricesFinal.jld2")
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas = matricesDict["matrices"]


cellCurls = Float64[]

# Plot for one set of 3 cells
for c=1:nCells

    cellVertices = findall(x->x!=0,C[c,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[c])...)
    end
    m = minimum(vertexAngles)
    vertexAngles .-= m
    cellVertices .= cellVertices[sortperm(vertexAngles)]

    cellEdges = findall(x->x!=0,B[c,:])
    edgeAngles = zeros(size(cellEdges))
    for (k,e) in enumerate(cellEdges)
        edgeAngles[k] = atan((edgeMidpoints[e].-cellPositions[c])...)
    end
    edgeAngles .+= (2π-m)
    edgeAngles .= edgeAngles.%(2π)
    cellEdges .= cellEdges[sortperm(edgeAngles)]


    h = @SVector [0.0,0.0]
    curlSum = 0
    for (i,e) in enumerate(cellEdges)
        h = h + ϵ*F[cellVertices[i],c]
        curlSum += B[c,e]*(h⋅edgeTangents[e])/cellAreas[c]
    end

    push!(cellCurls,curlSum)
end



# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1],aspect=DataAspect(),)
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    cellVertices = findall(x->x!=0,C[i,:])
    vertexAngles = zeros(size(cellVertices))
    for (k,v) in enumerate(cellVertices)
        vertexAngles[k] = atan((R[v].-cellPositions[i])...)
    end
    cellVertices .= cellVertices[sortperm(vertexAngles)]
    # poly!(ax,Point2f.(R[cellVertices]),color=(:black,1.0),strokecolor=(:black,1.0),strokewidth=5)
    poly!(ax,Point2f.(R[cellVertices]),color=(:green,cellCurls[i]),strokecolor=(:black,1.0),strokewidth=5)
end

display(fig)
