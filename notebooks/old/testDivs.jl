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
includet("$(projectdir())/notebooks/functions.jl")

dataDirectory = "data/sims/2022-02-28-19-30-22"

isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack R,A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellPositions,cellPerimeters,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,F,externalF,ϵ = matricesDict["matrices"]
@unpack initialSystem,nVerts,nCells,nEdges,γ,λ,preferredPerimeter,preferredArea,pressureExternal,dt,outputTotal,outputInterval,viscousTimeScale,realTimetMax,tMax,realCycleTime,nonDimCycleTime,t1Threshold = conditionsDict["params"]

# Create vector of polygons for each cell
cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

# Find cell midpoint links T
T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

# Create vector of triangles from midpoint links
linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

intersections = edgeLinkMidpoints(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)

q = calculateSpokes(conditionsDict["params"],matricesDict["matrices"])

onesVec = ones(1,nCells)
boundaryEdges = abs.(onesVec*B)
cᵖ = boundaryEdges'.*edgeMidpoints


vertexCurls = Float64[]
# Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
for k=1:nVerts
    curlSum = 0
    vertexCells = findall(x->x!=0,C[:,k])
    cellAngles = zeros(length(vertexCells))
    for i=1:length(cellAngles)
        cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
    end
    vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
    for i in vertexCells
        curlSum += (q[i,k]⋅(ϵ*F[k,i]))/linkTriangleAreas[k]
    end
    push!(vertexCurls,curlSum)
end

vertexCurlLims = (-maximum(abs.(vertexCurls)),maximum(abs.(vertexCurls)))


vertexDivs = Float64[]
# Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
for k=1:nVerts
    divSum = 0
    vertexCells = findall(x->x!=0,C[:,k])
    cellAngles = zeros(length(vertexCells))
    for i=1:length(cellAngles)
        cellAngles[i] = atan((cellPositions[vertexCells[i]].-R[k])...)
    end
    vertexCells .= vertexCells[sortperm(cellAngles,rev=true)]
    for i in vertexCells
        divSum += ((ϵ*q[i,k])⋅(ϵ*F[k,i]))/linkTriangleAreas[k]
    end
    push!(vertexDivs,divSum)
end

vertexDivLims = (-maximum(abs.(vertexDivs)),maximum(abs.(vertexDivs)))




# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=(:white,0.0),strokewidth=1)
end

for k=1:nVerts
    poly!(ax,linkTriangles[k],color=[vertexCurls[k]],colorrange=vertexCurlLims,colormap=:bwr,strokecolor=(:black,0.5),strokewidth=1)
end

for j=1:nEdges
    edgeCells = findall(!iszero,B[:,j])
    if boundaryEdges[j] == 0
        lines!(ax,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    else
        lines!(ax,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    end
end

scatter!(ax,Point2f.(R),alpha=0.5,color=:blue)
scatter!(ax,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
scatter!(ax,Point2f.(cellPositions),color=:red)
# scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=vertexCurlLims,colormap=:bwr,flipaxis=false)

display(fig)
save("$dataDirectory/pdf/vertexCurls.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/vertexCurls.pdf",fig)
save("$dataDirectory/svg/vertexCurls.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/vertexCurls.svg",fig)
save("$dataDirectory/png/vertexCurls.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/vertexCurls.png",fig)

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=(:white,0.0),strokewidth=1)
end

for k=1:nVerts
    poly!(ax,linkTriangles[k],color=[vertexDivs[k]],colorrange=vertexDivLims,colormap=:bwr,strokecolor=(:black,0.5),strokewidth=1)
end

for j=1:nEdges
    edgeCells = findall(!iszero,B[:,j])
    if boundaryEdges[j] == 0
        lines!(ax,Point2f.(cellPositions[edgeCells]),linewidth=1,color=(:white,1.0))
    else
        lines!(ax,Point2f.([cellPositions[edgeCells[1]],cᵖ[j]]),linewidth=1,color=(:white,1.0))
    end
end

scatter!(ax,Point2f.(R),alpha=0.5,color=:blue)
scatter!(ax,Point2f.(edgeMidpoints),alpha=0.5,color=:green)
scatter!(ax,Point2f.(cellPositions),color=:red)
# scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=vertexDivLims,colormap=:bwr,flipaxis=false)

display(fig)
save("$dataDirectory/pdf/vertexDivs.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/vertexDivs.pdf",fig)
save("$dataDirectory/svg/vertexDivs.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/vertexDivs.svg",fig)
save("$dataDirectory/png/vertexDivs.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/vertexDivs.png",fig)
