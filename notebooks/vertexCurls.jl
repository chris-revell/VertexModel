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

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,B,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices = matricesDict["matrices"]

onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints
T = SVector{2,Float64}[]
for j=1:nEdges
    Tⱼ = @SVector zeros(2)
    for i=1:nCells
        Tⱼ = Tⱼ + B[i,j]*(cellPositions[i].-cᵖ[j])
    end
    push!(T,Tⱼ)
end

linkTriangles = Vector{Point2f}[]
for k=1:nVerts
    if boundaryVertices[k] == 0
        vertexCells = findall(x->x!=0,C[:,k])
        push!(linkTriangles, Point2f.(cellPositions[vertexCells]))
    else
        vertexCells = findall(x->x!=0,C[:,k])
        vertexEdges = findall(x->x!=0,A[:,k])
        boundaryVertexEdges = intersect(vertexEdges,findall(x->x!=0,boundaryEdges[1,:]))
        kiteVertices = [edgeMidpoints[boundaryVertexEdges]; cellPositions[vertexCells]]
        push!(kiteVertices,R[k])
        com = sum(kiteVertices)./length(kiteVertices)
        angles = Float64[]
        for p=1:length(kiteVertices)
            angle = atan((kiteVertices[p].-com)...)
            push!(angles,angle)
        end
        kiteVertices .= kiteVertices[sortperm(angles)]
        push!(linkTriangles,Point2f.(kiteVertices))
    end
end
linkTriangleAreas = abs.(area.(linkTriangles))

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


E = linkTriangleAreas

vertexCurls = Float64[]
# Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
for k=1:nVerts

    vertexEdges = findall(x->x!=0,A[:,k])
    edgeAngles = zeros(length(vertexEdges))
    for (i,e) in enumerate(vertexEdges)
        edgeAngles[i] = atan((-A[e,k].*edgeTangents[e])...)
    end
    m = minimum(edgeAngles)
    edgeAngles .-= m
    vertexEdges .= vertexEdges[sortperm(edgeAngles,rev=true)]

    vertexCells = findall(x->x!=0,C[:,k])
    cellAngles = zeros(length(vertexCells))
    for i=1:length(cellAngles)
        cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
    end
    cellAngles .+= 2π-m
    cellAngles .= cellAngles.%2π
    vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]


    h = @SVector [0.0,0.0]
    curlSum = 0
    if boundaryVertices[k] == 0
        for (i,j) in enumerate(vertexEdges)
            h = h + ϵ*F[k,vertexCells[i]]
            curlSum += A[j,k]*(T[j]⋅h)/E[k]
        end
    end
    push!(vertexCurls,curlSum)
end



# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)

ax1.title = "Vertex curls"

lims = (minimum(vertexCurls),maximum(vertexCurls))

for k=1:nVerts
    poly!(ax1,linkTriangles[k],color=[vertexCurls[k]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.0)) #:bwr
end

# Plot cell polygons
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=5) #:bwr
end

Colorbar(fig[:,2],limits=lims,colormap=:bwr,flipaxis=false) #:bwr

display(fig)
save("$dataDirectory/vertexCurls.png",fig)
