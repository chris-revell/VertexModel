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

isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf")
isdir("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg")

# Import system data
conditionsDict    = load("$dataDirectory/dataFinal.jld2")
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
matricesDict = load("$dataDirectory/matricesFinal.jld2")
@unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
trapeziumAreas = abs.(area.(edgeTrapezia))

linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
linkTriangleAreas = abs.(area.(linkTriangles))

cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

Lf = makeLf(conditionsDict["params"],matricesDict["matrices"],trapeziumAreas)

# Calculate div on each cell
cellDivs = calculateCellDivs(conditionsDict["params"],matricesDict["matrices"])

onesVec = ones(nCells)
H = Diagonal(cellAreas)

eigenvectors = (eigen(Matrix(Lf))).vectors
eigenvalues = (eigen(Matrix(Lf))).values

ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(nCells))).*onesVec
ğ = cellDivs.-ḡ
ψ̆ = zeros(nCells)
eigenmodeAmplitudes = Float64[]
for k=2:nCells
    numerator = eigenvectors[:,k]'*H*ğ
    denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
    ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
    push!(eigenmodeAmplitudes,(numerator/denominator))
end

fig = Figure(resolution=(1000,700),fontsize = 24)
ax2 = Axis(fig[1,1], xlabel="Eigenmode number, i", ylabel="Amplitude",fontsize=32)
xlims!(ax2,0,nCells)
ylims!(ax2,0,1.1*maximum(abs.(eigenmodeAmplitudes)))
barplot!(ax2,collect(2:nCells),abs.(eigenmodeAmplitudes),width=1.0,color=:blue,strokecolor=:blue)

display(fig)
save("$dataDirectory/psicSpectrum.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/psicSpectrum.pdf",fig)
save("$dataDirectory/psicSpectrum.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/psicSpectrum.svg",fig)
save("$dataDirectory/psicSpectrum.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/psicSpectrum.png",fig)