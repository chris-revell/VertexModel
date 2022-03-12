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
@unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,pressureExternal,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
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
ϕ̆ = zeros(nCells)
eigenmodeAmplitudes = Float64[]
for k=2:nCells
    numerator = eigenvectors[:,k]'*H*ğ
    denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
    ϕ̆ .+= (numerator/denominator).*eigenvectors[:,k]
    push!(eigenmodeAmplitudes,(numerator/denominator))
end

ϕLims = (-maximum(abs.(ϕ̆)),maximum(abs.(ϕ̆)))

fig = Figure(resolution=(1000,1000),fontsize = 24)
ax1 = Axis(fig[1,1][1,1],aspect=DataAspect(),fontsize=32)
hidedecorations!(ax1)
hidespines!(ax1)
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[ϕ̆[i]],colormap=:bwr,colorrange=ϕLims, strokecolor=(:black,1.0),strokewidth=5)
end
Colorbar(fig[1,1][1,2],limits=ϕLims,colormap=:bwr,flipaxis=false,align=:left)

# ax2 = Axis(fig[2,:],title = LaTeXString("Eigenmode amplitudes"), xlabel=L"Eigenmode number, $i$", ylabel=LaTeXString("Amplitude"),fontsize=32)
# lines!(ax2,collect(2:nCells),abs.(eigenmodeAmplitudes),linewidth=3)

display(fig)
save("$dataDirectory/phicPotential.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/phicPotential.pdf",fig)
save("$dataDirectory/phicPotential.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/phicPotential.svg",fig)
save("$dataDirectory/phicPotential.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/phicPotential.png",fig)
