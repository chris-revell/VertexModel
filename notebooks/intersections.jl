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

dataDirectory = "data/sims/2022-03-16-16-02-03"

isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack initialSystem,nVerts,nCells,nEdges,γ,λ,preferredPerimeter,preferredArea,pressureExternal,dt,outputTotal,outputInterval,viscousTimeScale,realTimetMax,tMax,realCycleTime,nonDimCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack R,A,B,Aᵀ,Ā,Āᵀ,Bᵀ,B̄,B̄ᵀ,C,cellEdgeCount,boundaryVertices,cellPositions,cellPerimeters,cellAreas,cellTensions,cellPressures,edgeLengths,edgeTangents,edgeMidpoints,F,externalF,ϵ = matricesDict["matrices"]

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

vertexMidpointCurls = calculateVertexMidpointCurls(conditionsDict["params"],matricesDict["matrices"],intersections,linkTriangleAreas,q)
vertexMidpointCurlLims = (-maximum(abs.(vertexMidpointCurls)),maximum(abs.(vertexMidpointCurls)))

vertexMidpointDivs = calculateVertexMidpointDivs(conditionsDict["params"],matricesDict["matrices"],intersections,linkTriangleAreas,q)
vertexMidpointDivLims = (-maximum(abs.(vertexMidpointDivs)),maximum(abs.(vertexMidpointDivs)))

cellMidpointDivs = calculateCellMidpointDivs(conditionsDict["params"],matricesDict["matrices"],intersections,q)
cellMidpointDivLims = (-maximum(abs.(cellMidpointDivs)),maximum(abs.(cellMidpointDivs)))

cellMidpointCurls = calculateCellMidpointCurls(conditionsDict["params"],matricesDict["matrices"],intersections,q)
cellMidpointCurlLims = (-maximum(abs.(cellMidpointCurls)),maximum(abs.(cellMidpointCurls)))

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
    poly!(ax,linkTriangles[k],color=[vertexMidpointCurls[k]],colorrange=vertexMidpointCurlLims,colormap=:bwr,strokecolor=(:black,0.5),strokewidth=1)
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
scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=vertexMidpointCurlLims,colormap=:bwr,flipaxis=false)

save("$dataDirectory/pdf/vertexMidpointCurls.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/vertexMidpointCurls.pdf",fig)
save("$dataDirectory/svg/vertexMidpointCurls.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/vertexMidpointCurls.svg",fig)
save("$dataDirectory/png/vertexMidpointCurls.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/vertexMidpointCurls.png",fig)

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
    poly!(ax,linkTriangles[k],color=[vertexMidpointDivs[k]],colorrange=vertexMidpointDivLims,colormap=:bwr,strokecolor=(:black,0.5),strokewidth=1)
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
scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=vertexMidpointDivLims,colormap=:bwr,flipaxis=false)

save("$dataDirectory/pdf/vertexMidpointDivs.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/vertexMidpointDivs.pdf",fig)
save("$dataDirectory/svg/vertexMidpointDivs.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/vertexMidpointDivs.svg",fig)
save("$dataDirectory/png/vertexMidpointDivs.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/vertexMidpointDivs.png",fig)

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[cellMidpointDivs[i]],colorrange=cellMidpointCurlLims,colormap=:bwr,strokewidth=1)
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
scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=cellMidpointDivLims,colormap=:bwr,flipaxis=false)

save("$dataDirectory/pdf/cellMidpointDivs.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cellMidpointDivs.pdf",fig)
save("$dataDirectory/svg/cellMidpointDivs.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cellMidpointDivs.svg",fig)
save("$dataDirectory/png/cellMidpointDivs.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cellMidpointDivs.png",fig)

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()
ax = Axis(grid[1,1][1,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[cellMidpointCurls[i]],colorrange=cellMidpointCurlLims,colormap=:bwr,strokewidth=1)
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
scatter!(ax,Point2f.(intersections),color=:orange)
Colorbar(grid[1,1][1,2],limits=cellMidpointCurlLims,colormap=:bwr,flipaxis=false)

save("$dataDirectory/pdf/cellMidpointCurls.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/cellMidpointCurls.pdf",fig)
save("$dataDirectory/svg/cellMidpointCurls.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/cellMidpointCurls.svg",fig)
save("$dataDirectory/png/cellMidpointCurls.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/cellMidpointCurls.png",fig)
