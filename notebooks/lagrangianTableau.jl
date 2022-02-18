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

initialSystem = "data/sims/2022-02-18-11-59-24"

# Import system data
conditionsDict    = load("$initialSystem/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$initialSystem/matricesFinal.jld2")
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
F = 2.0.*trapeziumAreas


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



column = 10
decomposition = (eigen(Matrix(Lᵥ))).vectors

# Set up figure canvas
fig = Figure(resolution=(1000,1000))
grid = fig[1,1] = GridLayout()

axes = Axis[]
for x=1:4
    for y=1:4
        eigenvectorIndex = (y-1)*4 + x
        push!(axes,Axis(grid[y,x],aspect=DataAspect()))
        hidedecorations!(axes[end])
        hidespines!(axes[end])
        axes[end].title = "Eigenvector $eigenvectorIndex"
        for k=1:nVerts
            poly!(axes[end],linkTriangles[k],color=[decomposition[k,eigenvectorIndex]],colorrange=(minimum(decomposition[:,eigenvectorIndex]),maximum(decomposition[:,eigenvectorIndex])),colormap=:bwr,strokewidth=2,strokecolor=(:black,0.5)) #:bwr
        end
        # Plot cell polygons
        for i=1:nCells
            poly!(axes[end],cellPolygons[i],color=(:white,0.0),strokecolor=(:white,0.5),strokewidth=2) #:bwr
        end
    end
end


display(fig)
save("$(datadir())/plots/eigenvectorTableau.png",fig)
