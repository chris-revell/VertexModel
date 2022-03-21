# # Import Julia packages
# using DrWatson
# @quickactivate
# using Revise
# using LinearAlgebra
# using DelimitedFiles
# using SparseArrays
# using StaticArrays
# using CairoMakie
# using UnPack
# using GeometryBasics
# using Random
# using Colors
# using JLD2
#
# # Local modules
# includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
# includet("$(projectdir())/notebooks/functions.jl")
#
# #dataDirectory = "data/sims/2022-02-28-19-30-22"
#
# # Import system data
# conditionsDict    = load("$dataDirectory/dataFinal.jld2")
# @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
# matricesDict = load("$dataDirectory/matricesFinal.jld2")
# @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]
#
# onesVec = ones(1,nCells)
# boundaryEdges = abs.(onesVec*B)
# cᵖ = boundaryEdges'.*edgeMidpoints
#
# s = Matrix{SVector{2,Float64}}(undef,nCells,nVerts)
# fill!(s,@SVector zeros(2))
# for i=1:nCells
#     for k=1:nVerts
#         for j=1:nEdges
#             s[i,k] += B[i,j].*cᵖ[j].*A[j,k]
#         end
#     end
# end
#
# for k=1:nVerts
#     if boundaryVertices[k] !=0
#         vertexEdges = intersect(findall(!iszero,A[:,k]),findall(!iszero,boundaryEdges[1,:]))
#         sTest = edgeMidpoints[vertexEdges[1]].-edgeMidpoints[vertexEdges[2]]
#         display(sTest*ϵ*s[i,k])
#     end
# end
