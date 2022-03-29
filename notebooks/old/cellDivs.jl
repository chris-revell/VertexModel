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
@unpack A,B,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas = matricesDict["matrices"]


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


cellDivs = Float64[]
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
    divSum = 0
    for (i,e) in enumerate(cellEdges)
        h = h + ϵ*F[cellVertices[i],c]
        divSum -= B[c,e]*(h⋅(ϵ*edgeTangents[e]))/cellAreas[c]
    end
    divSum *= (-0.5)
    push!(cellDivs,divSum)
end



# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
# ax1.title = "Cell divs"
clims = (-maximum(abs.(cellDivs)),maximum(abs.(cellDivs)))
# Plot cell polygons
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[cellDivs[i]],colormap=:bwr,colorrange=clims, strokecolor=(:black,1.0),strokewidth=5)
end

Colorbar(fig[1, 2],limits=clims,colormap=:bwr)

display(fig)
save("$dataDirectory/cellDivs.png",fig)
