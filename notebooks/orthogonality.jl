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

boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
Tₑ2 = Diagonal(diagonalComponent)
Lf = (H\B)*Tₑ2*Bᵀ
dropzeros!(Lf)

invTₗ = inv(Tₗ)
boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1,:])
Lc = (H\B)*boundaryEdgesFactorMat*invTₗ*Bᵀ
dropzeros!(Lc)

cellAreasMat = Diagonal(cellAreas)

decomposition = (eigen(Matrix(Lf))).vectors
dotProducts = Float64[]
for i=1:nCells-1
    for j=i+1:nCells
        push!(dotProducts,decomposition[:,i]'*cellAreasMat*decomposition[:,j])
    end
end
println("Lf max dot product = $(maximum(dotProducts))")

decomposition = (eigen(Matrix(Lc))).vectors
dotProducts = Float64[]
for i=1:nCells-1
    for j=i+1:nCells
        push!(dotProducts,decomposition[:,i]'*cellAreasMat*decomposition[:,j])
    end
end
println("Lc max dot product = $(maximum(dotProducts))")

triangleAreasMat = Diagonal(linkTriangleAreas)

decomposition = (eigen(Matrix(Lₜ))).vectors
dotProducts = Float64[]
for i=1:nVerts-1
    for j=i+1:nVerts
        push!(dotProducts,decomposition[:,i]'*triangleAreasMat*decomposition[:,j])
    end
end
println("Lₜ max dot product = $(maximum(dotProducts))")

decomposition = (eigen(Matrix(Lᵥ))).vectors
dotProducts = Float64[]
for i=1:nVerts-1
    for j=i+1:nVerts
        push!(dotProducts,decomposition[:,i]'*triangleAreasMat*decomposition[:,j])
    end
end
println("Lᵥ max dot product = $(maximum(dotProducts))")
