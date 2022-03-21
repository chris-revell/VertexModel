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
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

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

edgeTrapezia = Vector{Point2f}[]
for j=1:nEdges
    edgeCells = findall(x->x!=0,B[:,j])
    edgeVertices = findall(x->x!=0,A[j,:])
    trapeziumVertices = [R[edgeVertices]; cellPositions[edgeCells]]
    com = sum(trapeziumVertices)./length(trapeziumVertices)
    angles = Float64[]
    for p=1:length(trapeziumVertices)
        angle = atan((trapeziumVertices[p].-com)...)
        push!(angles,angle)
    end
    trapeziumVertices .= trapeziumVertices[sortperm(angles)]
    push!(edgeTrapezia,Point2f.(trapeziumVertices))
end
trapeziumAreas = abs.(area.(edgeTrapezia))
#F = 2.0.*trapeziumAreas


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

H = Diagonal(cellAreas)
Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
boundaryEdgesFactor = abs.(boundaryEdges.-1)    # =1 for internal vertices, =0 for boundary vertices
invTₗ = inv(Tₗ)
boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1,:])
Lc = (H\B)*boundaryEdgesFactorMat*invTₗ*Bᵀ
dropzeros!(Lc)

decomposition = (eigen(Matrix(Lc))).vectors

# Set up figure canvas
fig = Figure(resolution=(550,1000))
grid = fig[1,1] = GridLayout()
for x=1:4
    for y=1:5
        eigenvectorIndex = ((y-1)*4 + x)+1
        lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
        ax = Axis(grid[y,x],aspect=DataAspect(),padding=0)
        hidedecorations!(ax)
        hidespines!(ax)
        # Plot cell polygons
        for i=1:nCells
            poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
        end
        Label(grid[y,x,Bottom()],
                L"i=%$eigenvectorIndex",
                textsize = 16,
        )
    end
end

eigenvectorIndex = 21+20*1
lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
ax = Axis(grid[6,1],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(grid[6,1,Bottom()],
        L"i=%$eigenvectorIndex",
        textsize = 16,
        #halign = :right
)

eigenvectorIndex = 21+20*2
lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
ax = Axis(grid[6,2],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
#ax.title = "Eigenmode $eigenvectorIndex"
# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(grid[6,2,Bottom()],
        L"i=%$eigenvectorIndex",
        textsize = 16,
        #halign = :right
)

eigenvectorIndex = 21+20*3
lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
ax = Axis(grid[6,3],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
#ax.title = "Eigenmode $eigenvectorIndex"
# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(grid[6,3,Bottom()],
        L"i=%$eigenvectorIndex",
        textsize = 16,
        #halign = :right
)

eigenvectorIndex = 21+20*4
lims = (minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex]))
ax = Axis(grid[6,4],aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)
#ax.title = "Eigenmode $eigenvectorIndex"
# Plot cell polygons
for i=1:nCells
    poly!(ax,cellPolygons[i],color=[decomposition[i,eigenvectorIndex]],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Label(grid[6,4,Bottom()],
        L"i=%$eigenvectorIndex",
        textsize = 16,
        #halign = :right
)

colgap!(grid, -5)
# rowgap!(grid, -2)

display(fig)
save("$dataDirectory/eigenvectorTableauLc.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/eigenvectorTableauLc.pdf",fig)
save("$dataDirectory/eigenvectorTableauLc.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/eigenvectorTableauLc.svg",fig)
save("$dataDirectory/eigenvectorTableauLc.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/eigenvectorTableauLc.png",fig)
