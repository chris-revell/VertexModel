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
includet("$(projectdir())/notebooks/functions.jl")

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

# Find cell midpoint links T
T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

Lᵥ = makeLv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)
Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)
Lc = makeLc(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)

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
