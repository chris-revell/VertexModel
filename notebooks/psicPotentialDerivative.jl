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


derivative = Lf*ψ̆

derivativeLims = (-maximum(abs.(derivative)),maximum(abs.(derivative)))

fig = Figure(resolution=(1000,1000),fontsize = 24)
ax1 = Axis(fig[1,1][1,1],aspect=DataAspect(),fontsize=32)
# ax1.title = LaTeXString("L_f\phi\breve")
hidedecorations!(ax1)
hidespines!(ax1)
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=[derivative[i]],colormap=:bwr,colorrange=derivativeLims, strokecolor=(:black,1.0),strokewidth=5)
end
Colorbar(fig[1,1][1,2],limits=derivativeLims,colormap=:bwr,flipaxis=false,align=:left)

ğLims = (-maximum(abs.(ğ)),maximum(abs.(ğ)))
fig2 = Figure(resolution=(1000,1000),fontsize = 24)
ax2 = Axis(fig2[1,1][1,1],aspect=DataAspect(),fontsize=32)
# ax1.title = LaTeXString("L_f\phi\breve")
hidedecorations!(ax2)
hidespines!(ax2)
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=[ğ[i]],colormap=:bwr,colorrange=ğLims, strokecolor=(:black,1.0),strokewidth=5)
end
Colorbar(fig2[1,1][1,2],limits=ğLims,colormap=:bwr,flipaxis=false,align=:left)

display(fig)
save("$dataDirectory/LfpsicDerivative.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/LfpsicDerivative.pdf",fig)
save("$dataDirectory/LfpsicDerivative.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/LfpsicDerivative.svg",fig)
save("$dataDirectory/LfpsicDerivative.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/LfpsicDerivative.png",fig)

save("$dataDirectory/gbreve.pdf",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/gbreve.pdf",fig2)
save("$dataDirectory/gbreve.svg",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/gbreve.svg",fig2)
save("$dataDirectory/gbreve.png",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/gbreve.png",fig2)