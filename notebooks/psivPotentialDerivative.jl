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

Lᵥ = makeLv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)
eigenvectors = (eigen(Matrix(Lᵥ))).vectors
eigenvalues = (eigen(Matrix(Lᵥ))).values


wideTildeVertexDivs = wideTildeVertexDiv(conditionsDict["params"],matricesDict["matrices"],linkTriangleAreas,trapeziumAreas)

onesVec = ones(nVerts)
E = Diagonal(linkTriangleAreas)

ḡ = ((onesVec'*E*wideTildeVertexDivs)/(onesVec'*E*ones(nVerts))).*onesVec
ğ = wideTildeVertexDivs.-ḡ
ψ̆ = zeros(nVerts)
eigenmodeAmplitudes = Float64[]
for k=2:nVerts
    numerator = -eigenvectors[:,k]'*E*ğ
    denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
    ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
    push!(eigenmodeAmplitudes,(numerator/denominator))
end

derivative = Lᵥ*ψ̆

derivativeLims = (-maximum(abs.(derivative)),maximum(abs.(derivative)))

fig = Figure(resolution=(1000,1000),fontsize = 24)
ax1 = Axis(fig[1,1][1,1],aspect=DataAspect(),fontsize=32)
hidedecorations!(ax1)
hidespines!(ax1)
for k=1:nVerts
    poly!(ax1,linkTriangles[k],color=[derivative[k]],colorrange=derivativeLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax1,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Colorbar(fig[1,1][1,2],limits=derivativeLims,colormap=:bwr,flipaxis=false,align=:left)

ğLims = (-maximum(abs.(ğ)),maximum(abs.(ğ)))
fig2 = Figure(resolution=(1000,1000),fontsize = 24)
ax2 = Axis(fig2[1,1][1,1],aspect=DataAspect(),fontsize=32)
ax2.title = LaTeXString("g\breve")
hidedecorations!(ax2)
hidespines!(ax2)
for k=1:nVerts
    poly!(ax2,linkTriangles[k],color=[-ğ[k]],colorrange=derivativeLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
end
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
end
Colorbar(fig2[1,1][1,2],limits=ğLims,colormap=:bwr,flipaxis=false,align=:left)


display(fig)
save("$dataDirectory/LvpsivDerivative.pdf",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/LvpsivDerivative.pdf",fig)
save("$dataDirectory/LvpsivDerivative.svg",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/LvpsivDerivative.svg",fig)
save("$dataDirectory/LvpsivDerivative.png",fig)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/LvpsivDerivative.png",fig)


save("$dataDirectory/gbreve_psiv.pdf",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/pdf/gbreve_psiv.pdf",fig2)
save("$dataDirectory/gbreve_psiv.svg",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/svg/gbreve_psiv.svg",fig2)
save("$dataDirectory/gbreve_psiv.png",fig2)
save("/Users/christopher/Dropbox (The University of Manchester)/VertexModelFigures/png/gbreve_psiv.png",fig2)
