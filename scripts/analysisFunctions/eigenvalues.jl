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
includet("$(projectdir())/scripts/analysisFunctions/functions.jl")

function eigenvalues(dataDirectory, show)

    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf")
    isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg") ? nothing : mkpath("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg")
    isdir("$dataDirectory/png") ? nothing : mkpath("$dataDirectory/png")
    isdir("$dataDirectory/pdf") ? nothing : mkpath("$dataDirectory/pdf")
    isdir("$dataDirectory/svg") ? nothing : mkpath("$dataDirectory/svg")

    # Import system data
    conditionsDict    = load("$dataDirectory/dataFinal.jld2")
    @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
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
    Lc = makeLc(conditionsDict["params"],matricesDict["matrices"],T,trapeziumAreas)
    eigenvaluesLc = (eigen(Matrix(Lc))).values
    Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)
    eigenvaluesLₜ = (eigen(Matrix(Lₜ))).values

    fig = Figure(resolution=(1000,700),fontsize = 32)
    ax1 = Axis(fig[1,1], xlabel="Eigenvalue number", ylabel="Eigenvalue",)
    lineLf = lines!(ax1,collect(1:nCells),eigenvaluesLf,linewidth=5)
    lineLc = lines!(ax1,collect(1:nCells),eigenvaluesLc,linewidth=5)
    lineLₜ = lines!(ax1,collect(1:nVerts),eigenvaluesLₜ,linewidth=5)
    lineLᵥ = lines!(ax1,collect(1:nVerts),eigenvaluesLᵥ,linewidth=5)
    Legend(fig[1, 2],[lineLf,lineLc,lineLₜ,lineLᵥ],["Lf","Lc","Lₜ","Lᵥ"])

    combinedEigenvalues = [eigenvaluesLᵥ[1:20] eigenvaluesLf[1:20] eigenvaluesLc[1:20] eigenvaluesLₜ[1:20]]

    writedlm("$dataDirectory/eigenvalues.txt",combinedEigenvalues,',')

    show==1 ? display(fig) : nothing
    save("$dataDirectory/svg/eigenvalueSpectrum.svg",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/svg/eigenvalueSpectrum.svg",fig)
    save("$dataDirectory/pdf/eigenvalueSpectrum.pdf",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/pdf/eigenvalueSpectrum.pdf",fig)
    save("$dataDirectory/png/eigenvalueSpectrum.png",fig)
    #save("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/$(splitdir(dataDirectory)[end])/png/eigenvalueSpectrum.png",fig)
end
