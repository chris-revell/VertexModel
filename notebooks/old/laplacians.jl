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

# Specify data folder
dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]


# Find vector of cell-cell links
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


# Create vector of edge trapezium polygons
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
# Calculate trapezium areas and F values
trapeziumAreas = abs.(area.(edgeTrapezia))
F = 2.0.*trapeziumAreas


# Create vector of polygons for triangles defined by cell-cell links around each vertex (note special consideration for boundary vertices)
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



# Define lagrangian matrices
H = Diagonal(cellAreas)
E = Diagonal(linkTriangleAreas)
Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
Lᵥ = (E\Aᵀ)*(Tₑ\A)
dropzeros!(Lᵥ)
Lₜ = (E\Aᵀ)*Tₗ*A
dropzeros!(Lₜ)
Lc = (H\B)*(Tₗ\Bᵀ)
dropzeros!(Lc)
Lf = (H\B)*Tₑ*Bᵀ
dropzeros!(Lf)

laplacians = Dict("Lv"=>Lᵥ,"Lt"=>Lₜ,"Lc"=>Lc,"Lf"=>Lf)

# Select which eigenvector to plot
column = 10

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()

ax1 = Axis(grid[1,1],aspect=DataAspect())
hidedecorations!(ax1)
hidespines!(ax1)
ax1.title = "Lᵥ Eigenvector $column"
decomposition = (eigen(Matrix(Lᵥ))).vectors
lims = (minimum(decomposition[:,column]),maximum(decomposition[:,column]))
for k=1:nVerts
    poly!(ax1,linkTriangles[k],color=[decomposition[k,column]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.5)) #:bwr
end
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,0.25),strokewidth=2) #:bwr
end

ax2 = Axis(grid[1,2],aspect=DataAspect())
hidedecorations!(ax2)
hidespines!(ax2)
ax2.title = "Lₜ Eigenvector $column"
decomposition = (eigen(Matrix(Lₜ))).vectors
lims = (minimum(decomposition[:,column]),maximum(decomposition[:,column]))
for k=1:nVerts
    poly!(ax2,linkTriangles[k],color=[decomposition[k,column]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.5)) #:bwr
end
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,0.25),strokewidth=2) #:bwr
end

ax3 = Axis(grid[2,1],aspect=DataAspect())
hidedecorations!(ax3)
hidespines!(ax3)
decomposition = (eigen(Matrix(Lc))).vectors
lims = (minimum(decomposition[:,column]),maximum(decomposition[:,column]))
ax3.title = "Lc Eigenvector $column"
for i=1:nCells
    poly!(ax3,cellPolygons[i],color=[decomposition[i,column]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.5)) #:bwr
end

ax4 = Axis(grid[2,2],aspect=DataAspect())
hidedecorations!(ax4)
hidespines!(ax4)
decomposition = (eigen(Matrix(Lf))).vectors
lims = (minimum(decomposition[:,column]),maximum(decomposition[:,column]))
ax4.title = "Lf Eigenvector $column"
for i=1:nCells
    poly!(ax4,cellPolygons[i],color=[decomposition[i,column]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.5)) #:bwr
end

display(fig)
save("$dataDirectory/eigenvectors$column.png",fig)
