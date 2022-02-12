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
@unpack A,B,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas = matricesDict["matrices"]


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

F = Float64[]
for i=1:nEdges
    push!(F,T[i]⋅(ϵ*edgeTangents[i]))
end
F .= abs.(F)

trapeziumAreas = 0.5.*F

display(sum(cellAreas))
display(sum(trapeziumAreas))


linkTriangleAreas = Float64[]
