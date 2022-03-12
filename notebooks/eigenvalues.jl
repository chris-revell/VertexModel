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
using Printf

# Local modules
includet("$(projectdir())/src/VertexModelContainers.jl"); using .VertexModelContainers
includet("$(projectdir())/notebooks/functions.jl")

dataDirectory = "data/sims/2022-02-28-19-30-22"

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

Lᵥ = makeLv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)
eigenvaluesLᵥ = (eigen(Matrix(Lᵥ))).values
Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)
eigenvaluesLf = (eigen(Matrix(Lf))).values
Lc = makeLc(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)
eigenvaluesLc = (eigen(Matrix(Lc))).values
Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
eigenvaluesLₜ = (eigen(Matrix(Lₜ))).values

fig = Figure(resolution=(1000,700),fontsize = 24)
ax1 = Axis(fig[1,1], xlabel="Eigenvalue number", ylabel="Eigenvalue",fontsize=32)
lineLf = lines!(ax1,collect(1:nCells),eigenvaluesLf,linewidth=5)
lineLc = lines!(ax1,collect(1:nCells),eigenvaluesLc,linewidth=5)
lineLₜ = lines!(ax1,collect(1:nVerts),eigenvaluesLₜ,linewidth=5)
lineLᵥ = lines!(ax1,collect(1:nVerts),eigenvaluesLᵥ,linewidth=5)
Legend(fig[1, 2],[lineLf,lineLc,lineLₜ,lineLᵥ],["Lf","Lc","Lₜ","Lᵥ"])


display(fig)
save("$dataDirectory/eigenvalueSpectrum.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/eigenvalueSpectrum.svg",fig)
save("$dataDirectory/eigenvalueSpectrum.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/eigenvalueSpectrum.pdf",fig)
save("$dataDirectory/eigenvalueSpectrum.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/eigenvalueSpectrum.png",fig)
